// Mm3Paper — reproduces every MM3-side measurement reported in the paper, so the
// published tables can be re-derived on any machine with `dotnet run -c Release`.
// Covers: mesh construction throughput and memory, query throughput (cached vs
// uncached), cache synchronization cost, symmetry canonicalization overhead,
// set-union kernels, and transpose thread scaling.
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime;
using System.Runtime.InteropServices;
using Numerical;

public readonly struct Vertex;
public readonly struct Cell;

internal static class Program
{
    // ------------------------------------------------------------------ utils

    /// <summary>
    ///     Runs <paramref name="f"/> a few times to warm the JIT, caches and branch
    ///     predictors, then reports the MINIMUM of <paramref name="reps"/> timed runs.
    ///     The minimum is the estimate least polluted by scheduling and frequency noise,
    ///     which matters on a laptop part under the powersave governor.
    /// </summary>
    private static double Best(Func<double> f, int reps, int warmup = 2)
    {
        for (var i = 0; i < warmup; i++) f();
        var best = double.MaxValue;
        for (var i = 0; i < reps; i++) best = Math.Min(best, f());
        return best;
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.NoInlining)]
    private static double TimeMs(Action a)
    {
        var sw = Stopwatch.StartNew();
        a();
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }

    /// <summary>
    ///     Times <paramref name="a"/> with collections suppressed, so a GC pause cannot
    ///     land inside the measured region. Falls back to a plain timing when the runtime
    ///     cannot satisfy the requested budget.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.NoInlining)]
    private static double TimeMsNoGC(Action a, long budgetBytes = 200L << 20)
    {
        GC.Collect(2, GCCollectionMode.Forced, true, true);
        GC.WaitForPendingFinalizers();
        var started = false;
        try { started = GC.TryStartNoGCRegion(budgetBytes, disallowFullBlockingGC: false); }
        catch (Exception) { started = false; }
        try
        {
            return TimeMs(a);
        }
        finally
        {
            if (started && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion)
                try { GC.EndNoGCRegion(); } catch (InvalidOperationException) { }
        }
    }

    private static long LiveBytes()
    {
        GC.Collect(2, GCCollectionMode.Forced, true, true);
        GC.WaitForPendingFinalizers();
        GC.Collect(2, GCCollectionMode.Forced, true, true);
        return GC.GetTotalMemory(true);
    }

    /// <summary>Structured quad mesh connectivity: n x n cells, (n+1)^2 vertices.</summary>
    private static int[][] QuadGrid(int n)
    {
        var cells = new int[n * n][];
        var t = 0;
        for (var r = 0; r < n; r++)
        for (var c = 0; c < n; c++)
        {
            int v00 = r * (n + 1) + c, v10 = v00 + 1;
            int v01 = (r + 1) * (n + 1) + c, v11 = v01 + 1;
            cells[t++] = new[] { v00, v10, v11, v01 };
        }
        return cells;
    }

    private static int SideFor(int targetElements) => (int)Math.Round(Math.Sqrt(targetElements));

    private static O2M BuildO2M(int[][] cells)
    {
        var o = new O2M(cells.Length);
        foreach (var c in cells) o.AppendElementCopy(c);
        return o;
    }

    // ------------------------------------------------------- 1. construction

    private static void Construction()
    {
        Console.WriteLine("\n## 1. Mesh construction (paper Table 10)\n");
        Console.WriteLine($"{"elements",12} {"time (ms)",11} {"elem/sec",13} {"memory (MB)",12}");
        foreach (var target in new[] { 10_000, 100_000, 1_000_000, 10_000_000 })
        {
            var n = SideFor(target);
            var cells = QuadGrid(n);
            var before = LiveBytes();
            O2M? keep = null;
            var ms = Best(() => TimeMs(() => keep = BuildO2M(cells)), 3);
            var after = LiveBytes();
            GC.KeepAlive(keep);
            var mb = (after - before) / (1024.0 * 1024.0);
            Console.WriteLine($"{cells.Length,12:N0} {ms,11:F1} {cells.Length / (ms / 1000.0),13:N0} {mb,12:F1}");
            keep = null;
        }
    }

    // ------------------------------------------------------ 2. query throughput

    private static void QueryThroughput()
    {
        Console.WriteLine("\n## 2. Query throughput on a 1M-element mesh (paper Table 11)\n");
        var n = SideFor(1_000_000);
        var cells = QuadGrid(n);
        var primary = BuildO2M(cells);
        using var m = new M2M(primary);
        m.EnsurePositionCaches();

        var elems = primary.Count;
        var nodes = primary.GetMaxNode() + 1;
        const int probes = 200_000;
        var rng = new Random(17);
        var eProbe = new int[probes];
        var nProbe = new int[probes];
        for (var i = 0; i < probes; i++) { eProbe[i] = rng.Next(elems); nProbe[i] = rng.Next(nodes); }

        // --- forward lookup, zero-copy span over the primary row
        long accS = 0;
        var msFwdSpan = Best(() => TimeMsNoGC(() =>
        {
            for (var i = 0; i < probes; i++) accS += primary[eProbe[i]].Length;
        }), 5);
        GC.KeepAlive(accS);
        Report("NodesOf<Element>  (span, zero-copy)", probes, msFwdSpan);

        // --- forward lookup, copying accessor
        var msFwd = Best(() => TimeMs(() =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += m.GetNodesForElement(eProbe[i]).Count;
            GC.KeepAlive(acc);
        }), 5);
        Report("NodesOf<Element>  (copying accessor)", probes, msFwd);

        // --- reverse lookup, zero-copy span over the cached transpose.
        // The delegate is allocated ONCE, outside the timed loop: allocating a
        // closure per probe would measure the allocator, not the lookup.
        long accR = 0;
        ReadOnlySpanAction<int> tally = s => accR += s.Length;
        var msRevSpan = Best(() => TimeMsNoGC(() =>
        {
            for (var i = 0; i < probes; i++)
                m.WithElementsForNodeSpan(nProbe[i], tally);
        }), 5);
        GC.KeepAlive(accR);
        Report("ElementsAt<Node>  (span, zero-copy)", probes, msRevSpan);

        // --- reverse lookup, copying accessor
        var msRev = Best(() => TimeMs(() =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += m.GetElementsForNode(nProbe[i]).Count;
            GC.KeepAlive(acc);
        }), 5);
        Report("ElementsAt<Node>  (copying accessor)", probes, msRev);

        // --- uncached baseline: linear scan for the same reverse lookup
        const int naiveProbes = 200;
        var msNaive = Best(() => TimeMs(() =>
        {
            long acc = 0;
            for (var i = 0; i < naiveProbes; i++)
            {
                var node = nProbe[i];
                for (var e = 0; e < primary.Count; e++)
                    foreach (var v in primary[e])
                        if (v == node) { acc++; break; }
            }
            GC.KeepAlive(acc);
        }), 1);
        Report("ElementsAt<Node> (uncached scan)", naiveProbes, msNaive);

        // --- neighbours
        const int nbProbes = 20_000;
        var msNb = Best(() => TimeMs(() =>
        {
            long acc = 0;
            for (var i = 0; i < nbProbes; i++) acc += m.GetElementNeighbors(eProbe[i]).Count;
            GC.KeepAlive(acc);
        }), 5);
        Report("GetDirectNeighbors", nbProbes, msNb);

        // --- k-hop (k = 2) through the typed API
        var topo = new Topology<TypeMap<Vertex, Cell>>();
        for (var i = 0; i < nodes; i++) topo.Add<Vertex>();
        topo.AddRangeParallel<Cell, Vertex>(cells, 4096);
        const int khopProbes = 2_000;
        var msK = Best(() => TimeMs(() =>
        {
            long acc = 0;
            for (var i = 0; i < khopProbes; i++) acc += topo.GetKHopNeighborhood<Cell, Vertex>(eProbe[i], 2).Count;
            GC.KeepAlive(acc);
        }), 3);
        Report("GetKHopNeighborhood(k=2)", khopProbes, msK);

        // --- connected components (whole mesh per call)
        var msCC = Best(() => TimeMs(() =>
        {
            var comps = topo.FindConnectedComponents<Cell, Vertex>();
            GC.KeepAlive(comps);
        }), 3);
        Report("FindConnectedComponents", 1, msCC);

        static void Report(string name, int ops, double ms)
            => Console.WriteLine($"{name,-36} {ops / (ms / 1000.0),15:N0} ops/sec   ({ms,9:F2} ms / {ops:N0} ops)");
    }

    // ---------------------------------------------- 3. cache synchronization

    private static void CacheSync()
    {
        Console.WriteLine("\n## 3. Cache synchronization after modification (paper Table 12)\n");
        Console.WriteLine($"{"elements",12} {"sync (ms)",11} {"query (us)",12} {"break-even",12}");
        foreach (var target in new[] { 10_000, 100_000, 1_000_000 })
        {
            var cells = QuadGrid(SideFor(target));
            var primary = BuildO2M(cells);
            using var m = new M2M(primary);
            var nodes = primary.GetMaxNode() + 1;

            // Force a rebuild: append one element, then time the first reverse query.
            var syncMs = Best(() =>
            {
                m.AppendElement(new List<int> { 0, 1, 2, 3 });
                return TimeMs(() => { var r = m.GetElementsForNode(1); GC.KeepAlive(r); });
            }, 3);

            // Steady-state reverse query once synchronized.
            const int probes = 100_000;
            var rng = new Random(5);
            var nProbe = new int[probes];
            for (var i = 0; i < probes; i++) nProbe[i] = rng.Next(nodes);
            var qMs = Best(() => TimeMs(() =>
            {
                long acc = 0;
                for (var i = 0; i < probes; i++) acc += m.GetElementsForNode(nProbe[i]).Count;
                GC.KeepAlive(acc);
            }), 3);
            var qUs = qMs * 1000.0 / probes;
            Console.WriteLine($"{primary.Count,12:N0} {syncMs,11:F2} {qUs,12:F3} {syncMs * 1000.0 / qUs,12:N0}");
        }
    }

    // ------------------------------------------- 4. symmetry canonicalization

    private static void SymmetryOverhead()
    {
        Console.WriteLine("\n## 4. Element addition with/without canonicalization (paper Table 13)\n");
        Console.WriteLine($"{"configuration",26} {"us/element",12} {"M elem/sec",12} {"overhead",10}");
        var quads = QuadGrid(SideFor(200_000));
        var tets = new int[quads.Length][];
        for (var i = 0; i < quads.Length; i++) tets[i] = new[] { quads[i][0], quads[i][1], quads[i][2], quads[i][3] };

        var baseUs = Measure(null);
        Row("No symmetry", baseUs, baseUs);
        Row("Dihedral D4 (quads)", Measure(Symmetry.Dihedral(4)), baseUs);
        Row("Full S4 (tetrahedra)", Measure(Symmetry.Full(4)), baseUs);

        double Measure(Symmetry? sym)
        {
            return Best(() =>
            {
                var topo = new Topology<TypeMap<Vertex, Cell>>();
                if (sym is not null) topo = topo.WithSymmetry<Cell>(sym);
                var nodeCount = quads.Max(q => q.Max()) + 1;
                for (var i = 0; i < nodeCount; i++) topo.Add<Vertex>();
                var ms = TimeMs(() =>
                {
                    foreach (var q in quads) topo.Add<Cell, Vertex>(q);
                });
                return ms * 1000.0 / quads.Length;   // microseconds per element
            }, 3);
        }

        static void Row(string name, double us, double baseUs)
            => Console.WriteLine($"{name,26} {us,12:F3} {1.0 / us,12:F2} {(us / baseUs - 1) * 100,9:F0}%");
    }

    // ---------------------------------------------------------- 5. set union

    private static void SetUnion()
    {
        Console.WriteLine("\n## 5. Set union: marker array vs hash set (paper Figure 14)\n");
        Console.WriteLine($"{"max node",10} {"marker (us)",13} {"hash (us)",12} {"speedup",9}");
        const int rows = 10_000, perRow = 4;
        foreach (var maxNode in new[] { 100, 1_000, 4_096, 8_192, 16_384 })
        {
            var a = RandomRows(rows, perRow, maxNode, 42);
            var b = RandomRows(rows, perRow, maxNode, 43);
            a.ParallelizationThreshold = int.MaxValue;   // compare kernels, not parallelism
            b.ParallelizationThreshold = int.MaxValue;

            var marker = Best(() => TimeMsNoGC(() => { var c = a | b; GC.KeepAlive(c); }), 9) * 1000.0;
            var hash = Best(() => TimeMsNoGC(() => { var c = HashUnion(a, b); GC.KeepAlive(c); }), 9) * 1000.0;
            Console.WriteLine($"{maxNode,10} {marker,13:F1} {hash,12:F1} {hash / marker,8:F2}x");
        }

        static O2M HashUnion(O2M x, O2M y)
        {
            var n = Math.Max(x.Count, y.Count);
            var res = new O2M(n);
            var seen = new HashSet<int>();
            var row = new List<int>(8);
            for (var i = 0; i < n; i++)
            {
                seen.Clear(); row.Clear();
                if (i < x.Count) foreach (var v in x[i]) if (seen.Add(v)) row.Add(v);
                if (i < y.Count) foreach (var v in y[i]) if (seen.Add(v)) row.Add(v);
                res.AppendElementCopy(CollectionsMarshal.AsSpan(row));
            }
            return res;
        }
    }

    private static O2M RandomRows(int rows, int perRow, int maxNode, int seed)
    {
        var rng = new Random(seed);
        var o = new O2M(rows);
        var seen = new HashSet<int>();
        var buf = new List<int>(perRow);
        for (var i = 0; i < rows; i++)
        {
            seen.Clear(); buf.Clear();
            while (buf.Count < perRow) { var v = rng.Next(maxNode); if (seen.Add(v)) buf.Add(v); }
            o.AppendElementCopy(CollectionsMarshal.AsSpan(buf));
        }
        return o;
    }

    // --------------------------------------------------- 6. transpose scaling

    private static void TransposeScaling()
    {
        Console.WriteLine("\n## 6. Transpose thread scaling (paper Section 9)\n");
        foreach (var target in new[] { 1_000_000, 4_000_000 })
        {
            var a = BuildO2M(QuadGrid(SideFor(target)));
            Console.WriteLine($"\n{a.Count:N0} elements ({a.Count * 4:N0} pairs)");
            Console.WriteLine($"{"threads",9} {"time (ms)",11} {"speedup",9} {"efficiency",11}");

            a.ParallelizationThreshold = int.MaxValue;
            var seq = Best(() => TimeMs(() => { var t = a.Transpose(); GC.KeepAlive(t); }), 5);
            Console.WriteLine($"{"1 (seq)",9} {seq,11:F1} {1.00,9:F2} {"--",11}");

            a.ParallelizationThreshold = 4096;
            foreach (var threads in ThreadCounts())
            {
                ParallelConfig.MaxDegreeOfParallelism = threads;
                var t = Best(() => TimeMs(() => { var x = a.Transpose(); GC.KeepAlive(x); }), 5);
                Console.WriteLine($"{threads,9} {t,11:F1} {seq / t,9:F2} {seq / t / threads * 100,10:F0}%");
            }
            ParallelConfig.MaxDegreeOfParallelism = Environment.ProcessorCount;
        }
    }

    /// <summary>Powers of two up to the number of cores this process may actually use.</summary>
    private static int[] ThreadCounts()
    {
        var max = Environment.ProcessorCount;
        var list = new List<int>();
        for (var t = 2; t <= max; t *= 2) list.Add(t);
        if (list.Count == 0 || list[^1] != max) list.Add(max);
        return list.ToArray();
    }

    private static void Main()
    {
        // Ask for the highest scheduling priority the runtime will grant, and keep the
        // GC in a low-latency mode outside the explicitly suppressed regions.
        try { Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.High; } catch (Exception) { /* needs privileges */ }
        try { Thread.CurrentThread.Priority = ThreadPriority.Highest; } catch (Exception) { }
        GCSettings.LatencyMode = GCLatencyMode.SustainedLowLatency;

        Console.WriteLine("MM3 paper measurements");
        Console.WriteLine($"machine: {Environment.ProcessorCount} usable cores, .NET {Environment.Version}");
        Console.WriteLine($"config : server GC={GCSettings.IsServerGC}, latency={GCSettings.LatencyMode}, "
                          + $"tiered={Environment.GetEnvironmentVariable("DOTNET_TieredCompilation") ?? "default"}");
        Construction();
        QueryThroughput();
        CacheSync();
        SymmetryOverhead();
        SetUnion();
        TransposeScaling();
        Console.WriteLine("\ndone.");
    }
}

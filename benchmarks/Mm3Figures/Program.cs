// Mm3Figures — regenerates the measured data behind Figures 8-13 of the paper.
// Emits one whitespace-separated .txt per figure (first column = abscissa, further
// columns = curves), which the plotgnu .data specs consume directly.
//
// Run through run-maxperf.sh so the measurements are taken pinned to the
// performance cores with tiering disabled; see BENCHMARK_COMPARISON.md.
using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Numerical;

public readonly struct Vertex;
public readonly struct Cell;

internal static class Program
{
    private static string _outDir = ".";

    // ------------------------------------------------------------------ utils

    [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.NoInlining)]
    private static double TimeMs(Action a)
    {
        var sw = Stopwatch.StartNew();
        a();
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }

    /// <summary>Warm up, then take the minimum of <paramref name="reps"/> timings.</summary>
    private static double Best(Action a, int reps = 3, int warmup = 1)
    {
        for (var i = 0; i < warmup; i++) a();
        var best = double.MaxValue;
        for (var i = 0; i < reps; i++) best = Math.Min(best, TimeMs(a));
        return best;
    }

    private static void Emit(string name, string header, IEnumerable<double[]> rows)
    {
        var path = Path.Combine(_outDir, name);
        using var w = new StreamWriter(path);
        w.WriteLine("# " + header);
        foreach (var r in rows)
            w.WriteLine(string.Join("  ", r.Select(v => v.ToString("G6", CultureInfo.InvariantCulture))));
        Console.WriteLine($"  wrote {path}");
    }

    /// <summary>Structured quad-grid connectivity for ~<paramref name="target"/> elements.</summary>
    private static int[][] Grid(int target)
    {
        var n = (int)Math.Round(Math.Sqrt(target));
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

    private static O2M Build(int[][] cells)
    {
        var o = new O2M(cells.Length);
        foreach (var c in cells) o.AppendElementCopy(c);
        return o;
    }

    private static void Free() { GC.Collect(2, GCCollectionMode.Forced, true, true); GC.WaitForPendingFinalizers(); }

    // ------------------------------------- Fig 8: assembling vs element count

    private static readonly int[] Sizes = { 1_000_000, 2_000_000, 4_000_000, 8_000_000, 16_000_000 };

    private static void FigAssembly()
    {
        Console.WriteLine("\n[Fig 8] assembling vs number of elements");
        var rows = new List<double[]>();
        foreach (var target in Sizes)
        {
            var cells = Grid(target);
            var n = cells.Length;

            var flex = Best(() => { var o = Build(cells); GC.KeepAlive(o); Free(); }, reps: 4);

            // CSR: the same connectivity assembled into the flat compressed form.
            var csr = Best(() =>
            {
                var o = Build(cells);
                var (rp, ci) = o.ToCsr();
                GC.KeepAlive(rp); GC.KeepAlive(ci);
                Free();
            }, reps: 4);

            rows.Add(new double[] { n / 1e6, flex, csr });
            Console.WriteLine($"    {n,12:N0}  flexible {flex,9:F1} ms   CSR {csr,9:F1} ms");
            cells = null; Free();
        }
        Emit("fig08_assembly.txt", "n_el[1e6]  flexible[ms]  csr[ms]", rows);
    }

    // ------------------------------------------ Fig 9: assembling vs threads

    private static void FigAssemblyThreads()
    {
        Console.WriteLine("\n[Fig 9] assembling vs thread count (bulk parallel path)");
        const int target = 4_000_000;
        var cells = Grid(target);
        var nodes = cells.Max(c => c.Max()) + 1;
        var rows = new List<double[]>();

        foreach (var threads in ThreadCounts())
        {
            ParallelConfig.MaxDegreeOfParallelism = threads;

            // (a) bulk add with a registered symmetry group: canonicalization is
            //     precomputed in parallel, insertion is sequential.
            var sym = Best(() =>
            {
                var t = new Topology<TypeMap<Vertex, Cell>>().WithSymmetry<Cell>(Symmetry.Dihedral(4));
                t.AddRange<Vertex, int>(new ReadOnlySpan<int>(new int[nodes]));   // bulk node add
                t.AddRangeParallel<Cell, Vertex>(cells, 4096);
                GC.KeepAlive(t); Free();
            }, reps: 2);

            // (b) bulk add without symmetry: no parallel stage worth the name.
            var plain = Best(() =>
            {
                var t = new Topology<TypeMap<Vertex, Cell>>();
                t.AddRange<Vertex, int>(new ReadOnlySpan<int>(new int[nodes]));   // bulk node add
                t.AddRangeParallel<Cell, Vertex>(cells, 4096);
                GC.KeepAlive(t); Free();
            }, reps: 2);

            rows.Add(new double[] { threads, sym, plain });
            Console.WriteLine($"    {threads,3} threads   with symmetry {sym,9:F1} ms   plain {plain,9:F1} ms");
        }
        ParallelConfig.MaxDegreeOfParallelism = Environment.ProcessorCount;
        Emit("fig09_assembly_threads.txt", "threads  with_symmetry[ms]  plain[ms]", rows);
    }

    // ------------------------------------------------- Fig 10: orderings

    private static void FigOrderings()
    {
        Console.WriteLine("\n[Fig 10] lexicographic ordering and compression/renumbering");
        var rows = new List<double[]>();
        foreach (var target in Sizes)
        {
            var cells = Grid(target);
            var o = Build(cells);
            var n = o.Count;
            var maxNode = o.GetMaxNode();

            var lex = Best(() => { var s = o.GetSortOrder(); GC.KeepAlive(s); }, reps: 2);

            // compression + renumbering: keep every element, identity node permutation
            var keep = new List<int>(n);
            for (var i = 0; i < n; i++) keep.Add(i);
            var perm = new List<int>(maxNode + 1);
            for (var i = 0; i <= maxNode; i++) perm.Add(i);
            var comp = Best(() =>
            {
                var c = (O2M)o.Clone();
                c.RearrangeAfterRenumbering(keep, perm);
                GC.KeepAlive(c);
            }, reps: 2);

            rows.Add(new double[] { n / 1e6, lex, comp });
            Console.WriteLine($"    {n,12:N0}  lexicographic {lex,9:F1} ms   compress {comp,9:F1} ms");
            cells = null; Free();
        }
        Emit("fig10_orderings.txt", "n_el[1e6]  lexicographic[ms]  compress_renumber[ms]", rows);
    }

    // --------------------------------- Fig 11: set operations, O2M and CSR

    /// <summary>
    ///     Figure 11 delegates to the sampling implementation in <see cref="OperationsFigure"/>.
    /// </summary>
    /// <remarks>
    ///     This used to be an inline min-of-two measurement. That is not enough to reject a
    ///     collection landing inside a timed region, and it produced a 2x outlier at 4M
    ///     elements on the intersection curve. The robust path is the only path now, so a
    ///     full-suite run cannot regenerate the figure with the weaker protocol.
    /// </remarks>
    private static void FigOperations() => OperationsFigure.Run(_outDir);

    // ------------------------------- Fig 12: set operations vs thread count

    private static void FigOperationsThreads()
    {
        Console.WriteLine("\n[Fig 12] set operations vs thread count");
        int[] densities = { 4_000_000, 8_000_000 };
        var all = new Dictionary<int, List<double[]>>();
        foreach (var target in densities)
        {
            var a = Build(Grid(target));
            var b = Build(Grid(target));
            a.ParallelizationThreshold = 4096;
            b.ParallelizationThreshold = 4096;
            var list = new List<double[]>();
            foreach (var threads in ThreadCounts())
            {
                ParallelConfig.MaxDegreeOfParallelism = threads;
                var u = Best(() => { var c = a | b; GC.KeepAlive(c); }, reps: 2);
                var s = Best(() => { var c = a - b; GC.KeepAlive(c); }, reps: 2);
                var i = Best(() => { var c = a & b; GC.KeepAlive(c); }, reps: 2);
                list.Add(new double[] { threads, u, s, i });
                Console.WriteLine($"    {target,10:N0}  {threads,3}t  + {u,8:F1}  - {s,8:F1}  & {i,8:F1} ms");
            }
            all[target] = list;
            Free();
        }
        ParallelConfig.MaxDegreeOfParallelism = Environment.ProcessorCount;
        var rows = all[densities[0]].Select((r, i) => new[]
        {
            r[0], r[1], r[2], r[3],
            all[densities[1]][i][1], all[densities[1]][i][2], all[densities[1]][i][3],
        });
        Emit("fig12_operations_threads.txt",
            "threads  u4M  d4M  i4M  u8M  d8M  i8M  [ms]", rows);
    }

    // --------------------------------- Fig 13: positions and clique formation

    private static void FigCliques()
    {
        Console.WriteLine("\n[Fig 13] local positions and clique formation");
        var rows = new List<double[]>();
        foreach (var target in new[] { 250_000, 500_000, 1_000_000, 2_000_000, 4_000_000 })
        {
            var primary = Build(Grid(target));
            var n = primary.Count;
            var transpose = primary.Transpose();

            var nodePos = Best(() =>
            {
                var p = O2M.GetNodePositions(primary, transpose);
                GC.KeepAlive(p);
            }, reps: 2);

            var elemPos = Best(() =>
            {
                var p = O2M.GetElementPositions(primary, transpose);
                GC.KeepAlive(p);
            }, reps: 2);

            var cliques = Best(() =>
            {
                var c = O2M.GetCliques(primary, transpose);
                GC.KeepAlive(c);
            }, reps: 2);

            rows.Add(new double[] { n / 1e6, elemPos, nodePos, cliques });
            Console.WriteLine($"    {n,12:N0}  elem-pos {elemPos,8:F1}  node-pos {nodePos,8:F1}  cliques {cliques,9:F1} ms");
            Free();
        }
        Emit("fig13_cliques.txt", "n_el[1e6]  element_pos[ms]  node_pos[ms]  cliques[ms]", rows);
    }

    private static int[] ThreadCounts()
    {
        var max = Environment.ProcessorCount;
        var l = new List<int> { 1 };
        for (var t = 2; t <= max; t *= 2) l.Add(t);
        if (l[^1] != max) l.Add(max);
        return l.ToArray();
    }

    private static void Main(string[] args)
    {
        try { Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.High; } catch (Exception) { }
        // `Mm3Figures ops <outdir>` re-measures only the set-operation figures.
        if (args.Length > 0 && args[0] == "ops")
        {
            var dir = args.Length > 1 ? args[1] : ".";
            Directory.CreateDirectory(dir);
            OperationsFigure.Run(dir);
            return;
        }

        _outDir = args.Length > 0 ? args[0] : ".";
        Directory.CreateDirectory(_outDir);
        Console.WriteLine($"MM3 figure data -> {Path.GetFullPath(_outDir)}");
        Console.WriteLine($"machine: {Environment.ProcessorCount} usable cores, .NET {Environment.Version}");

        FigAssembly();
        FigAssemblyThreads();
        FigOrderings();
        FigOperations();
        FigOperationsThreads();
        FigCliques();
        Console.WriteLine("\ndone.");
    }
}

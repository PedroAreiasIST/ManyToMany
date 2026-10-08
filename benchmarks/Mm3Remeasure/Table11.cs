// Full re-measurement of paper Table 11: every row, both columns.
//
// "Our Method" is measured on two access shapes where the framework offers both:
//   span  = zero-copy row read, read lock taken once for the batch
//   api   = public typed Topology call, lock and allocation per call
// "Naive" is an implementation holding ONLY the primary element->node O2M, with no
// transpose and no position caches, exactly as the paper's caption defines it. Every
// reverse or transitive query therefore scans.
using System.Diagnostics;
using Numerical;

public readonly struct V2;

public readonly struct C2;

/// <summary>Reverse and transitive queries with no transpose and no position caches.</summary>
internal sealed class Naive(int[][] cells, int nodeCount)
{
    public List<int> ElementsAt(int node)
    {
        var r = new List<int>();
        for (var e = 0; e < cells.Length; e++)
            foreach (var v in cells[e])
                if (v == node)
                {
                    r.Add(e);
                    break;
                }

        return r;
    }

    public IReadOnlyList<int> NodesOf(int element) => cells[element];

    /// <summary>Copying forward accessor, the shape a cache-free API would expose.</summary>
    public List<int> NodesOfCopy(int element) => [..cells[element]];

    public List<int> Neighbors(int element)
    {
        var mine = cells[element];
        var r = new List<int>();
        for (var e = 0; e < cells.Length; e++)
        {
            if (e == element) continue;
            var hit = false;
            foreach (var v in cells[e])
            {
                foreach (var w in mine)
                    if (v == w)
                    {
                        hit = true;
                        break;
                    }

                if (hit) break;
            }

            if (hit) r.Add(e);
        }

        return r;
    }

    public Dictionary<int, int> KHop(int seed, int k)
    {
        var dist = new Dictionary<int, int> { [seed] = 0 };
        var frontier = new List<int> { seed };
        for (var level = 1; level <= k; level++)
        {
            var next = new List<int>();
            foreach (var cur in frontier)
            foreach (var nb in Neighbors(cur))
                if (dist.TryAdd(nb, level))
                    next.Add(nb);
            frontier = next;
        }

        dist.Remove(seed);
        return dist;
    }

    /// <summary>Full connected-components pass with the uncached neighbour query.</summary>
    public int Components()
    {
        var visited = new bool[cells.Length];
        var count = 0;
        for (var i = 0; i < cells.Length; i++)
        {
            if (visited[i]) continue;
            count++;
            var q = new Queue<int>();
            q.Enqueue(i);
            visited[i] = true;
            while (q.Count > 0)
            foreach (var nb in Neighbors(q.Dequeue()))
                if (!visited[nb])
                {
                    visited[nb] = true;
                    q.Enqueue(nb);
                }
        }

        return count;
    }

    public int NodeCount => nodeCount;
}

internal static class Table11
{
    private static double TimeMs(Action a)
    {
        var sw = Stopwatch.StartNew();
        a();
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }

    private static double Best(Func<double> f, int reps, string tag)
    {
        var s = new double[reps];
        for (var i = 0; i < reps; i++)
        {
            GC.Collect(2, GCCollectionMode.Forced, true, true);
            GC.WaitForPendingFinalizers();
            GC.Collect(2, GCCollectionMode.Forced, true, true);
            s[i] = f();
        }

        Array.Sort(s);
        Console.WriteLine($"        [{tag}] {string.Join(" ", s.Select(x => x.ToString("F2")))}");
        return s[0];
    }

    private static int[][] QuadGrid(int n)
    {
        var cells = new int[n * n][];
        var t = 0;
        for (var r = 0; r < n; r++)
        for (var c = 0; c < n; c++)
        {
            int v00 = r * (n + 1) + c, v10 = v00 + 1;
            int v01 = (r + 1) * (n + 1) + c, v11 = v01 + 1;
            cells[t++] = [v00, v10, v11, v01];
        }

        return cells;
    }

    private static void Row(string name, long ops, double ms)
        => Console.WriteLine($"  {name,-46} {ops / (ms / 1000.0),16:N0} ops/sec  ({ms,10:F2} ms)");

    public static void Run()
    {
        Console.WriteLine("\n\n================ TABLE 11, FULL RE-MEASUREMENT ================");
        Console.WriteLine($".NET {Environment.Version}, serverGC={System.Runtime.GCSettings.IsServerGC}, " +
                          $"procs={Environment.ProcessorCount}\n");

        var cells = QuadGrid(1000);
        var primary = new O2M(cells.Length);
        foreach (var c in cells) primary.AppendElementCopy(c);
        using var m = new M2M(primary);
        m.EnsurePositionCaches();

        var elems = primary.Count;
        var nodes = primary.GetMaxNode() + 1;
        var naive = new Naive(cells, nodes);

        var topo = new Topology<TypeMap<V2, C2>>();
        for (var i = 0; i < nodes; i++) topo.Add<V2>();
        topo.AddRangeParallel<C2, V2>(cells, 4096);

        Console.WriteLine($"mesh: {elems:N0} quad elements, {nodes:N0} nodes\n");

        var rng = new Random(17);
        var eP = new int[200_000];
        var nP = new int[200_000];
        for (var i = 0; i < eP.Length; i++)
        {
            eP[i] = rng.Next(elems);
            nP[i] = rng.Next(nodes);
        }

        // ---- warm-up of every path
        long w = 0;
        for (var i = 0; i < 2000; i++) w += primary[eP[i]].Length + naive.NodesOf(eP[i]).Count;
        m.WithElementsFromNode(t => { for (var i = 0; i < 2000; i++) w += t[nP[i]].Length; });
        for (var i = 0; i < 200; i++) w += topo.GetDirectNeighbors<C2, V2>(eP[i]).Count;
        for (var i = 0; i < 50; i++) w += topo.GetKHopNeighborhood<C2, V2>(eP[i], 2).Count;
        for (var i = 0; i < 3; i++) w += naive.Neighbors(eP[i]).Count + naive.ElementsAt(nP[i]).Count;
        for (var i = 0; i < 2; i++) w += naive.KHop(eP[i], 2).Count;
        GC.KeepAlive(w);

        // ================= row 1: NodesOf<Element> =================
        Console.WriteLine("NodesOf<Element>");
        const int fwd = 200_000;
        long a1 = 0;
        Row("  our: primary span, no lock", fwd,
            Best(() => TimeMs(() => { for (var i = 0; i < fwd; i++) a1 += primary[eP[i]].Length; }), 7, "1a"));
        long a2 = 0;
        Row("  our: GetNodesForElement (copying API)", fwd,
            Best(() => TimeMs(() => { for (var i = 0; i < fwd; i++) a2 += m.GetNodesForElement(eP[i]).Count; }), 7,
                "1b"));
        long a3 = 0;
        Row("  naive: row read, no caches involved", fwd,
            Best(() => TimeMs(() => { for (var i = 0; i < fwd; i++) a3 += naive.NodesOf(eP[i]).Count; }), 7, "1c"));
        long a4 = 0;
        Row("  naive: copying accessor", fwd,
            Best(() => TimeMs(() => { for (var i = 0; i < fwd; i++) a4 += naive.NodesOfCopy(eP[i]).Count; }), 7, "1d"));
        GC.KeepAlive(a1 + a2 + a3 + a4);

        // ================= row 2: ElementsAt<Node> =================
        Console.WriteLine("\nElementsAt<Node>");
        long b1 = 0;
        Row("  our: transpose span, lock hoisted", fwd,
            Best(() => TimeMs(() => m.WithElementsFromNode(t =>
            {
                for (var i = 0; i < fwd; i++) b1 += t[nP[i]].Length;
            })), 7, "2a"));
        long b2 = 0;
        Row("  our: GetElementsForNode (copying API)", fwd,
            Best(() => TimeMs(() => { for (var i = 0; i < fwd; i++) b2 += m.GetElementsForNode(nP[i]).Count; }), 7,
                "2b"));
        const int scanProbes = 200;
        long b3 = 0;
        Row("  naive: full scan, no transpose", scanProbes,
            Best(() => TimeMs(() => { for (var i = 0; i < scanProbes; i++) b3 += naive.ElementsAt(nP[i]).Count; }), 3,
                "2c"));
        GC.KeepAlive(b1 + b2 + b3);

        // ================= row 3: GetDirectNeighbors =================
        Console.WriteLine("\nGetDirectNeighbors");
        const int nb = 20_000;
        long c1 = 0;
        Row("  our: Topology.GetDirectNeighbors", nb,
            Best(() => TimeMs(() => { for (var i = 0; i < nb; i++) c1 += topo.GetDirectNeighbors<C2, V2>(eP[i]).Count; }),
                5, "3a"));
        long c2 = 0;
        Row("  our: M2M.GetElementNeighbors (unsorted)", nb,
            Best(() => TimeMs(() => { for (var i = 0; i < nb; i++) c2 += m.GetElementNeighbors(eP[i], false).Count; }),
                5, "3b"));
        const int nbNaive = 100;
        long c3 = 0;
        Row("  naive: scan all elements per query", nbNaive,
            Best(() => TimeMs(() => { for (var i = 0; i < nbNaive; i++) c3 += naive.Neighbors(eP[i]).Count; }), 3,
                "3c"));
        GC.KeepAlive(c1 + c2 + c3);

        // ================= row 4: GetKHopNeighborhood(k=2) =================
        Console.WriteLine("\nGetKHopNeighborhood(k=2)");
        const int kh = 2_000;
        long d1 = 0;
        Row("  our: Topology.GetKHopNeighborhood", kh,
            Best(() => TimeMs(() => { for (var i = 0; i < kh; i++) d1 += topo.GetKHopNeighborhood<C2, V2>(eP[i], 2).Count; }),
                5, "4a"));
        const int khNaive = 10;
        long d2 = 0;
        Row("  naive: BFS over uncached neighbours", khNaive,
            Best(() => TimeMs(() => { for (var i = 0; i < khNaive; i++) d2 += naive.KHop(eP[i], 2).Count; }), 3, "4b"));
        GC.KeepAlive(d1 + d2);

        // ================= row 5: FindConnectedComponents =================
        Console.WriteLine("\nFindConnectedComponents");
        var msCc = Best(() => TimeMs(() =>
        {
            var comps = topo.FindConnectedComponents<C2, V2>();
            GC.KeepAlive(comps);
        }), 5, "5a");
        Row("  our: Topology.FindConnectedComponents", 1, msCc);

        // A naive components pass is n applications of the uncached neighbour query, so
        // it is quadratic in n. Measure it on three small meshes, check the quadratic
        // model holds, then extrapolate to 10^6 -- a direct run would take hours.
        Console.WriteLine("  naive: measured on small meshes, quadratic model, extrapolated to 10^6");
        double? prevMs = null;
        double prevN = 0;
        foreach (var (side, reps) in new[] { (80, 3), (160, 3), (320, 2) })
        {
            var sc = QuadGrid(side);
            var sn = new Naive(sc, (side + 1) * (side + 1));
            var ms = Best(() => TimeMs(() =>
            {
                var k = sn.Components();
                GC.KeepAlive(k);
            }), reps, $"5b/n={side * side}");
            double n = sc.Length;
            var model = prevMs is null ? double.NaN : ms / prevMs.Value / Math.Pow(n / prevN, 2);
            var extrapSec = ms / 1000.0 * Math.Pow(1e6 / n, 2);
            Console.WriteLine($"        n={n,7:N0}  {ms,9:F2} ms" +
                              (double.IsNaN(model) ? "        " : $"   n^2 model = {model,4:F2}") +
                              $"   -> 10^6 pass = {extrapSec,12:N0} s = {1.0 / extrapSec,10:E2} ops/sec");
            prevMs = ms;
            prevN = n;
        }
    }
}

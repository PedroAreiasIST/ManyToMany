// Paper Tables 12 and 13, plus the fast-path share of Table 4.
//
// Table 12: cost of the first reverse query after a topological modification (the
// transpose rebuild), the steady-state reverse query that follows, and the number of
// queries at which the amortized rebuild drops below the steady-state query cost.
// The transpose and the position caches carry separate validity flags, so both the
// transpose-only rebuild and the rebuild including the position caches are reported.
//
// Table 13: element insertion with and without a registered symmetry group.
//
// Fast-path share: the fraction of queries answered by a direct indexed read of a
// synchronized structure, for one completely specified workload. The workload is the
// query sequence a finite-element pass issues over a quadrilateral mesh: for every
// element, read its node row; for every node of that element, read the transpose row
// and the element-position entry; and, for one element in sixteen, take the
// direct-neighbor query, which is computed on demand rather than read.
using System.Diagnostics;
using Numerical;

public readonly struct V4;

public readonly struct C4;

internal static class Tables1213
{
    private static double TimeMs(Action a)
    {
        var sw = Stopwatch.StartNew();
        a();
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }

    private static double[] Samples(Func<double> f, int reps)
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
        return s;
    }

    private static string Show(double[] s) =>
        $"min {s[0]:F3}  median {s[s.Length / 2]:F3}  max {s[^1]:F3}";

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

    public static void Run()
    {
        Console.WriteLine("\n\n========= TABLES 12 & 13 + FAST-PATH SHARE =========");
        Console.WriteLine($".NET {Environment.Version}, serverGC={System.Runtime.GCSettings.IsServerGC}\n");

        // ---------------------------------------------------------------- Table 12
        Console.WriteLine("## Table 12 -- cache synchronization after one modification\n");
        Console.WriteLine($"{"mesh",10} {"transpose (ms)",16} {"+positions (ms)",17} " +
                          $"{"query (us)",12} {"break-even",12}");

        foreach (var target in new[] { 10_000, 100_000, 1_000_000 })
        {
            var side = (int)Math.Round(Math.Sqrt(target));
            var cells = QuadGrid(side);
            var n = cells.Length;
            var extra = new List<int> { 0, 1, 2, 3 };

            // transpose rebuild alone: append one element, then force the reverse path
            var tS = Samples(() =>
            {
                var p = new O2M(n);
                foreach (var c in cells) p.AppendElementCopy(c);
                using var m = new M2M(p);
                m.GetElementsForNode(0);            // first sync, not timed
                m.AppendElements(extra);             // invalidate
                return TimeMs(() => m.GetElementsForNode(0));
            }, 7);

            // rebuild including both position caches
            var tP = Samples(() =>
            {
                var p = new O2M(n);
                foreach (var c in cells) p.AppendElementCopy(c);
                using var m = new M2M(p);
                m.EnsurePositionCaches();
                m.AppendElements(extra);
                return TimeMs(() =>
                {
                    m.GetElementsForNode(0);
                    m.EnsurePositionCaches();
                });
            }, 7);

            // steady-state reverse query, zero-copy row read under one lock
            var prim = new O2M(n);
            foreach (var c in cells) prim.AppendElementCopy(c);
            using var mm = new M2M(prim);
            mm.EnsurePositionCaches();
            var nodes = prim.GetMaxNode() + 1;
            var rng = new Random(17);
            var probe = new int[200_000];
            for (var i = 0; i < probe.Length; i++) probe[i] = rng.Next(nodes);
            long acc = 0;
            mm.WithElementsFromNode(t => { for (var i = 0; i < 20_000; i++) acc += t[probe[i]].Length; });
            var q = Samples(() => TimeMs(() =>
                mm.WithElementsFromNode(t =>
                {
                    for (var i = 0; i < probe.Length; i++) acc += t[probe[i]].Length;
                })), 7);
            GC.KeepAlive(acc);
            var qUs = q[0] * 1000.0 / probe.Length;

            Console.WriteLine($"{target,10:N0} {tS[0],16:F2} {tP[0],17:F2} {qUs,12:F4} " +
                              $"{tP[0] * 1000.0 / qUs,12:N0}");
            Console.WriteLine($"           transpose  [{Show(tS)}]");
            Console.WriteLine($"           +positions [{Show(tP)}]");
            Console.WriteLine($"           query      [{Show(q)}] over {probe.Length:N0} probes");
        }

        // ---------------------------------------------------------------- Table 13
        Console.WriteLine("\n## Table 13 -- element addition with and without canonicalization\n");
        const int addN = 200_000;
        var quad = QuadGrid((int)Math.Round(Math.Sqrt(addN)));

        var maxNodeW = 0;
        foreach (var c in quad)
        foreach (var v in c)
            if (v > maxNodeW) maxNodeW = v;

        // Only the element-addition loop is timed; the nodes are created first, outside
        // the timed region, so the figure is the per-element insertion cost and nothing
        // else. Timing node creation as well dilutes the canonicalization overhead.
        double PerElement(Symmetry? sym)
        {
            var s = Samples(() =>
            {
                var topo = new Topology<TypeMap<V4, C4>>();
                if (sym is not null) topo.WithSymmetry<C4>(sym);
                for (var i = 0; i <= maxNodeW; i++) topo.Add<V4>();
                var ms = TimeMs(() =>
                {
                    foreach (var c in quad) topo.Add<C4, V4>(c);
                });
                GC.KeepAlive(topo);
                return ms;
            }, 5);
            return s[0] * 1000.0 / quad.Length; // microseconds per element
        }

        var plain = PerElement(null);
        var d4 = PerElement(Symmetry.Dihedral(4));
        var s4 = PerElement(Symmetry.Full(4));
        Console.WriteLine($"  {"configuration",-26} {"us/element",12} {"M/sec",10} {"overhead",10}");
        Console.WriteLine($"  {"No symmetry",-26} {plain,12:F3} {1 / plain,10:F1} {"---",10}");
        Console.WriteLine($"  {"Dihedral D4 (quads)",-26} {d4,12:F3} {1 / d4,10:F1} {(d4 / plain - 1) * 100,9:F0}%");
        Console.WriteLine($"  {"Full S4 (tetrahedra)",-26} {s4,12:F3} {1 / s4,10:F1} {(s4 / plain - 1) * 100,9:F0}%");

        // -------------------------------------------------------- fast-path share
        Console.WriteLine("\n## Fast-path share for the finite-element workload W\n");
        var sideW = (int)Math.Round(Math.Sqrt(1_000_000));
        var cellsW = QuadGrid(sideW);
        var pw = new O2M(cellsW.Length);
        foreach (var c in cellsW) pw.AppendElementCopy(c);
        using var mw = new M2M(pw);
        mw.EnsurePositionCaches();

        long fast = 0, slow = 0;
        var locs = mw.ElementLocations;
        mw.WithElementsFromNode(t =>
        {
            for (var e = 0; e < pw.Count; e++)
            {
                var row = pw[e];
                fast++;                                  // primary row read
                for (var p = 0; p < row.Length; p++)
                {
                    fast++;                              // transpose row read
                    var _ = t[row[p]].Length;
                    fast++;                              // element-position entry
                    var __ = locs[e][p];
                }

                if ((e & 15) == 0) slow++;               // neighbor query, computed
            }
        });
        Console.WriteLine($"  workload W over {pw.Count:N0} elements of 4 nodes:");
        Console.WriteLine($"    direct indexed reads (fast path) : {fast:N0}");
        Console.WriteLine($"    computed on demand (slow path)   : {slow:N0}");
        Console.WriteLine($"    fast-path share                  : {100.0 * fast / (fast + slow):F2}%");
        Console.WriteLine("  The share is a property of this call mix, not of the framework:");
        for (var every = 1; every <= 4; every++)
        {
            long s2 = pw.Count / (1 << (2 * every));
            Console.WriteLine($"    one neighbor query per {1 << (2 * every),3} elements -> " +
                              $"{100.0 * fast / (fast + s2):F2}%");
        }
    }
}

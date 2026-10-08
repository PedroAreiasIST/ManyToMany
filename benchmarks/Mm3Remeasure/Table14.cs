// Paper Table 14 rows for MM3 on the 100K-element quadrilateral mesh: build time,
// NodesOf latency, neighbour latency, cached and uncached. The Boost and NetworkX
// rows are third-party figures and are not re-measured here.
using System.Diagnostics;
using Numerical;

public readonly struct V3;

public readonly struct C3;

internal static class Table14
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
        Console.WriteLine($"        [{tag}] {string.Join(" ", s.Select(x => x.ToString("F3")))}");
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

    public static void Run()
    {
        Console.WriteLine("\n\n============ TABLE 14 -- MM3 rows, 100K quad mesh ============\n");

        var side = (int)Math.Round(Math.Sqrt(100_000));
        var cells = QuadGrid(side);
        var elems = cells.Length;

        // --- build: primary relation + the caches the cached row depends on
        var msBuildCached = Best(() => TimeMs(() =>
        {
            var p = new O2M(cells.Length);
            foreach (var c in cells) p.AppendElementCopy(c);
            using var mm = new M2M(p);
            mm.EnsurePositionCaches();
            GC.KeepAlive(mm);
        }), 7, "build cached");

        var msBuildPlain = Best(() => TimeMs(() =>
        {
            var p = new O2M(cells.Length);
            foreach (var c in cells) p.AppendElementCopy(c);
            GC.KeepAlive(p);
        }), 7, "build uncached");

        var primary = new O2M(cells.Length);
        foreach (var c in cells) primary.AppendElementCopy(c);
        using var m = new M2M(primary);
        m.EnsurePositionCaches();
        var nodes = primary.GetMaxNode() + 1;

        var topo = new Topology<TypeMap<V3, C3>>();
        for (var i = 0; i < nodes; i++) topo.Add<V3>();
        topo.AddRangeParallel<C3, V3>(cells, 4096);

        var naive = new Naive(cells, nodes);

        var rng = new Random(17);
        var eP = new int[200_000];
        for (var i = 0; i < eP.Length; i++) eP[i] = rng.Next(elems);

        long w = 0;
        for (var i = 0; i < 2000; i++) w += primary[eP[i]].Length + naive.NodesOf(eP[i]).Count;
        for (var i = 0; i < 200; i++) w += topo.GetDirectNeighbors<C3, V3>(eP[i]).Count;
        for (var i = 0; i < 20; i++) w += naive.Neighbors(eP[i]).Count;
        GC.KeepAlive(w);

        const int fwd = 200_000;
        long a = 0;
        var msFwd = Best(() => TimeMs(() => { for (var i = 0; i < fwd; i++) a += primary[eP[i]].Length; }), 7, "NodesOf cached");
        long b = 0;
        var msFwdN = Best(() => TimeMs(() => { for (var i = 0; i < fwd; i++) b += naive.NodesOf(eP[i]).Count; }), 7, "NodesOf uncached");

        const int nb = 20_000;
        long c1 = 0;
        var msNb = Best(() => TimeMs(() => { for (var i = 0; i < nb; i++) c1 += topo.GetDirectNeighbors<C3, V3>(eP[i]).Count; }), 5, "Neighbors cached");
        const int nbN = 300;
        long c2 = 0;
        var msNbN = Best(() => TimeMs(() => { for (var i = 0; i < nbN; i++) c2 += naive.Neighbors(eP[i]).Count; }), 3, "Neighbors uncached");
        GC.KeepAlive(a + b + c1 + c2);

        Console.WriteLine($"\n  {elems:N0} elements, {nodes:N0} nodes\n");
        Console.WriteLine($"  {"row",-24} {"Build (ms)",12} {"NodesOf (ns)",14} {"Neighbors (us)",16}");
        Console.WriteLine($"  {"Ours (cached)",-24} {msBuildCached,12:F1} {msFwd * 1e6 / fwd,14:F1} " +
                          $"{msNb * 1000.0 / nb,16:F2}");
        Console.WriteLine($"  {"Ours (uncached)",-24} {msBuildPlain,12:F1} {msFwdN * 1e6 / fwd,14:F1} " +
                          $"{msNbN * 1000.0 / nbN,16:F2}");
        Console.WriteLine($"\n  neighbour speedup cached/uncached = {(nb / msNb) / (nbN / msNbN),10:N0}x");
        Console.WriteLine($"  vs published Boost 12 us          = {12.0 / (msNb * 1000.0 / nb),10:F1}x faster");
    }
}

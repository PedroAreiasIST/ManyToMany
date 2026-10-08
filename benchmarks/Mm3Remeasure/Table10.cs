// Paper Table 10: mesh construction time, throughput and memory.
//
// "Construction" is the build of the primary element->node relation, which is what the
// O(m)-per-element append claim describes and what Table 14's build column reports.
// Memory is measured two ways, because the manuscript uses two accountings: the live
// bytes actually retained by the built structure (the default List<List<int>>
// representation), and the flat index count of Remark 6.2 (4 bytes per stored index).
using System.Diagnostics;
using Numerical;

internal static class Table10
{
    private static double TimeMs(Action a)
    {
        var sw = Stopwatch.StartNew();
        a();
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }

    private static long LiveBytes()
    {
        for (var i = 0; i < 3; i++)
        {
            GC.Collect(2, GCCollectionMode.Forced, true, true);
            GC.WaitForPendingFinalizers();
        }

        return GC.GetTotalMemory(true);
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

    private static O2M Build(int[][] cells)
    {
        var o = new O2M(cells.Length);
        foreach (var c in cells) o.AppendElementCopy(c);
        return o;
    }

    public static void Run()
    {
        Console.WriteLine("\n\n================ TABLE 10 -- construction ================");
        Console.WriteLine($".NET {Environment.Version}, serverGC={System.Runtime.GCSettings.IsServerGC}\n");
        Console.WriteLine($"{"target",10} {"elements",12} {"time (ms)",11} {"elem/sec",14} " +
                          $"{"live MB",10} {"flat MB",9} {"B/elem",8}");

        foreach (var target in new[] { 10_000, 100_000, 1_000_000, 10_000_000 })
        {
            var side = (int)Math.Round(Math.Sqrt(target));
            var cells = QuadGrid(side);
            var n = cells.Length;

            // warm the path at this size before timing
            var warm = Build(cells);
            GC.KeepAlive(warm);
            warm = null;

            var reps = n >= 5_000_000 ? 3 : 7;
            var samples = new double[reps];
            for (var i = 0; i < reps; i++)
            {
                GC.Collect(2, GCCollectionMode.Forced, true, true);
                GC.WaitForPendingFinalizers();
                GC.Collect(2, GCCollectionMode.Forced, true, true);
                O2M? keep = null;
                samples[i] = TimeMs(() => keep = Build(cells));
                GC.KeepAlive(keep);
                keep = null;
            }

            Array.Sort(samples);
            var ms = samples[0];

            // retained bytes of the built structure alone: `cells` is already live
            var before = LiveBytes();
            var built = Build(cells);
            var after = LiveBytes();
            var liveMb = (after - before) / (1024.0 * 1024.0);
            GC.KeepAlive(built);
            built = null;

            var flatMb = n * 4L * 4L / (1024.0 * 1024.0); // 4 nodes x 4 bytes per element

            // full topology as a user holds it: primary + cached transpose + both
            // position caches -- the quantity Figure 15 and Remark 6.2 report
            var beforeAll = LiveBytes();
            var prim = Build(cells);
            var m2m = new M2M(prim);
            m2m.EnsurePositionCaches();
            var afterAll = LiveBytes();
            var cachedMb = (afterAll - beforeAll) / (1024.0 * 1024.0);
            var cachedPerElem = (afterAll - beforeAll) / (double)n;
            GC.KeepAlive(m2m);
            m2m.Dispose();

            Console.WriteLine($"{target,10:N0} {n,12:N0} {ms,11:F1} {n / (ms / 1000.0),14:N0} " +
                              $"{liveMb,10:F1} {flatMb,9:F1} {(after - before) / (double)n,8:F1}");
            Console.WriteLine($"           samples [{string.Join(" ", samples.Select(s => s.ToString("F1")))}]");
            Console.WriteLine($"           with caches: {cachedMb,8:F1} MB live "
                              + $"({cachedPerElem:F1} B/elem), flat index count "
                              + $"{n * 16L * 4L / (1024.0 * 1024.0):F1} MB");
        }
    }
}

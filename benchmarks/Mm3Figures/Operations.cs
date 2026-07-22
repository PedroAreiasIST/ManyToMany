// Focused, statistically defensible re-measurement of the set-operation figures
// (Figure 11a/b). Invoked as `Mm3Figures ops <outdir>`.
//
// The first version of this measurement took the minimum of two timed repetitions,
// which is not enough to reject a collection landing inside a timed region: it
// produced a single 2x outlier at 4M elements on the intersection curve while the
// neighbouring sizes stayed on trend. Here every point is sampled REPS times with
// the heap settled before each sample, and the full sample set is printed so an
// outlier is visible in the log rather than silently averaged into the figure.
//
// The reported statistic is the MINIMUM of the samples. Each repetition allocates a
// complete result structure, so at the larger sizes the later samples carry growing
// allocation pressure -- the spread reaches ~2x at 16M elements while the minima stay
// on a straight line. The minimum is the sample least contaminated by that pressure
// and by scheduler interference, and is the conventional estimator for this kind of
// measurement; the median would import the allocator's behaviour into a figure that
// is meant to show the cost of the operation.
using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Numerical;

internal static class OperationsFigure
{
    private const int Warmup = 2;
    private const int Reps = 7;

    /// <summary>
    ///     Times one sample and reports whether a generation-2 collection ran inside the
    ///     timed region, so contaminated samples can be identified rather than guessed at.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization | MethodImplOptions.NoInlining)]
    private static (double Ms, bool Contaminated) TimeOnce(Action a)
    {
        // Settle the heap so a pending collection cannot be charged to this sample.
        GC.Collect(2, GCCollectionMode.Forced, true, true);
        GC.WaitForPendingFinalizers();
        GC.Collect(2, GCCollectionMode.Forced, true, true);
        var gen2Before = GC.CollectionCount(2);
        var sw = Stopwatch.StartNew();
        a();
        sw.Stop();
        return (sw.Elapsed.TotalMilliseconds, GC.CollectionCount(2) != gen2Before);
    }

    /// <summary>
    ///     Warm up, then return every timed sample sorted ascending, together with the
    ///     number of samples during which a gen-2 collection ran. Uncontaminated samples
    ///     are preferred; if every sample was contaminated they are all kept, so the
    ///     measurement degrades to "noisy" rather than to "missing".
    /// </summary>
    private static (double[] Sorted, int Contaminated) Samples(Action a)
    {
        for (var i = 0; i < Warmup; i++) a();
        var clean = new List<double>(Reps);
        var dirty = new List<double>(Reps);
        for (var i = 0; i < Reps; i++)
        {
            var (ms, bad) = TimeOnce(a);
            (bad ? dirty : clean).Add(ms);
        }
        var chosen = clean.Count > 0 ? clean : dirty;
        var s = chosen.ToArray();
        Array.Sort(s);
        return (s, dirty.Count);
    }

    private static double Median(double[] sorted) =>
        sorted.Length % 2 == 1
            ? sorted[sorted.Length / 2]
            : 0.5 * (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]);

    private static void Report(string name, double n, (double[] Sorted, int Contaminated) r)
    {
        var s = r.Sorted;
        var spread = s[^1] / s[0];
        var flag = r.Contaminated > 0 ? $"  gen2-hit {r.Contaminated}/{Reps}" : "";
        Console.WriteLine($"    {name,-14} n={n,5:F1}M  min {s[0],8:F1}  median {Median(s),8:F1}  max {s[^1],8:F1}  " +
                          $"spread {spread,4:F2}x  [{string.Join(" ", s.Select(v => v.ToString("F0", CultureInfo.InvariantCulture)))}]{flag}");
    }

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

    public static void Run(string outDir)
    {
        Console.WriteLine("\n[Fig 11] set operations, re-measured "
                          + $"({Warmup} warm-up + {Reps} timed samples per point, heap settled before each)\n");
        int[] sizes = { 1_000_000, 2_000_000, 4_000_000, 8_000_000, 16_000_000 };
        var o2m = new List<double[]>();
        var csr = new List<double[]>();

        foreach (var target in sizes)
        {
            var a = Build(Grid(target));
            var b = Build(Grid(target));
            a.ParallelizationThreshold = 4096;
            b.ParallelizationThreshold = 4096;
            var n = a.Count / 1e6;

            var su = Samples(() => { var c = a | b; GC.KeepAlive(c); });
            var sd = Samples(() => { var c = a - b; GC.KeepAlive(c); });
            var si = Samples(() => { var c = a & b; GC.KeepAlive(c); });
            Report("O2M union", n, su);
            Report("O2M diff", n, sd);
            Report("O2M inter", n, si);
            o2m.Add(new[] { n, su.Sorted[0], sd.Sorted[0], si.Sorted[0] });

            var (rpA, ciA) = a.ToCsr();
            var (rpB, ciB) = b.ToCsr();
            var nCols = Math.Max(a.GetMaxNode(), b.GetMaxNode()) + 1;
            var ca = new CSR(rpA, ciA, nCols, null, skipValidation: true);
            var cb = new CSR(rpB, ciB, nCols, null, skipValidation: true);

            var cu = Samples(() => { var c = ca + cb; GC.KeepAlive(c); });
            var cd = Samples(() => { var c = ca - cb; GC.KeepAlive(c); });
            var ci = Samples(() => { var c = ca & cb; GC.KeepAlive(c); });
            Report("CSR union", n, cu);
            Report("CSR diff", n, cd);
            Report("CSR inter", n, ci);
            csr.Add(new[] { n, cu.Sorted[0], cd.Sorted[0], ci.Sorted[0] });

            ca.Dispose(); cb.Dispose();
            Console.WriteLine();
            GC.Collect(2, GCCollectionMode.Forced, true, true);
            GC.WaitForPendingFinalizers();
        }

        Emit(Path.Combine(outDir, "fig11a_operations_o2m.txt"),
            "n_el[1e6]  union[ms]  difference[ms]  intersection[ms]   (minimum of 7)", o2m);
        Emit(Path.Combine(outDir, "fig11b_operations_csr.txt"),
            "n_el[1e6]  union[ms]  difference[ms]  intersection[ms]   (minimum of 7)", csr);
    }

    private static void Emit(string path, string header, List<double[]> rows)
    {
        using var w = new StreamWriter(path);
        w.WriteLine("# " + header);
        foreach (var r in rows)
            w.WriteLine(string.Join("  ", r.Select(v => v.ToString("G6", CultureInfo.InvariantCulture))));
        Console.WriteLine($"  wrote {path}");
    }
}

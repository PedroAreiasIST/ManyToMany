// Measures the thread scaling of O2M.Transpose() on this machine.
// The transpose assigns write positions in one sequential pass and then fills
// in parallel, so this quantifies the speedup actually attainable.
using System.Diagnostics;
using System.Runtime.InteropServices;
using Numerical;

static O2M BuildQuadMesh(int cellsPerSide)
{
    var n = cellsPerSide;
    var o = new O2M(n * n);
    var row = new int[4];
    for (var r = 0; r < n; r++)
    for (var c = 0; c < n; c++)
    {
        int v00 = r * (n + 1) + c, v10 = v00 + 1;
        int v01 = (r + 1) * (n + 1) + c, v11 = v01 + 1;
        row[0] = v00; row[1] = v10; row[2] = v11; row[3] = v01;
        o.AppendElementCopy(row);
    }
    return o;
}

static double TimeTranspose(O2M a, int reps)
{
    var sw = new Stopwatch();
    // warm-up
    for (var i = 0; i < 2; i++) { var w = a.Transpose(); GC.KeepAlive(w); }
    var best = double.MaxValue;
    for (var r = 0; r < reps; r++)
    {
        sw.Restart();
        var t = a.Transpose();
        sw.Stop();
        GC.KeepAlive(t);
        best = Math.Min(best, sw.Elapsed.TotalMilliseconds);
    }
    return best;
}

foreach (var cells in new[] { 1000, 2000 })      // 1M and 4M elements
{
    var a = BuildQuadMesh(cells);
    Console.WriteLine($"\n=== {a.Count:N0} elements ({a.Count * 4:N0} element-node pairs) ===");
    Console.WriteLine($"{"threads",8} {"time (ms)",12} {"speedup",9} {"efficiency",11}");

    // Sequential reference: threshold above the element count disables the parallel fill.
    a.ParallelizationThreshold = int.MaxValue;
    var seq = TimeTranspose(a, 5);
    Console.WriteLine($"{"1 (seq)",8} {seq,12:F1} {1.0,9:F2} {"--",11}");

    a.ParallelizationThreshold = 4096;   // library default: parallel fill enabled
    foreach (var threads in new[] { 2, 4, 8, 16, 24, 32 })
    {
        ParallelConfig.MaxDegreeOfParallelism = threads;
        var t = TimeTranspose(a, 5);
        Console.WriteLine($"{threads,8} {t,12:F1} {seq / t,9:F2} {seq / t / threads * 100,10:F0}%");
    }
}

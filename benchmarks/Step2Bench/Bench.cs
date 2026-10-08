// Bench.cs — measurement modes for the step2 paper.
//
// Every timing follows the protocol of the released MM3 harness: Release build,
// tiered compilation off, server GC, threads pinned to the performance cores by
// run-maxperf.sh, a full blocking collection before each timed region, and the
// minimum of several samples where a measurement can be repeated. Incremental
// steps are repeated exactly through the persistence layer itself: the step is
// applied, timed, rolled back to the version before it, and applied again.
//
//   mm3ref                         MM3 Table 12 protocol (baseline reproduction)
//   work  kind n steps every pers  per-step maintenance vs rebuild (crack|bisect|erode)
//   crossover n                    batch size at which a rebuild becomes cheaper
//   interleaved n1,n2,...          recursive bisection on MM3 (lazy rebuild) vs incremental
//   reads n                        single-thread read throughput
//   concurrent n seconds readers   readers during writes: snapshots vs reader-writer lock
//   persist n                      snapshot, rollback and retained-version costs
//   memory n                       live bytes of each representation
//   endtoend n                     maintenance vs assembly vs solve per crack step
//   cutpar n g threads reps        concurrent cutting of g^3 independent edge groups

using System.Diagnostics;
using System.Globalization;
using Numerical;

namespace Step2.Bench;

internal static partial class Bench
{
    private static readonly ElementSpec Tet = ElementSpec.Tetrahedron;

    // Output directory: $STEP2_DATA, else ./data under the working directory (the paper's
    // own copy resolved step2/benchmarks/data relative to its bin folder, which does not
    // exist in this repository).
    internal static string DataDir =>
        Directory.CreateDirectory(Environment.GetEnvironmentVariable("STEP2_DATA") ??
                                  Path.Combine(Environment.CurrentDirectory, "data")).FullName;

    public static void Run(string mode, string[] a)
    {
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
        Directory.CreateDirectory(DataDir);
        Console.WriteLine($"# step2 {mode} {string.Join(' ', a)} | .NET {Environment.Version} | " +
                          $"serverGC={System.Runtime.GCSettings.IsServerGC} | cpus={Environment.ProcessorCount}");
        int I(int i, int d) => a.Length > i ? int.Parse(a[i]) : d;
        if (mode != "mm3ref") WarmUp();
        switch (mode)
        {
            case "mm3ref": Mm3Reference(); break;
            case "work": Work(a.Length > 0 ? a[0] : "crack", I(1, 26), I(2, 0), I(3, 10), I(4, 1) == 1); break;
            case "crossover": Crossover(I(0, 55)); break;
            case "baseline-order": BaselineOrder(I(0, 55)); break;
            case "interleaved": Interleaved(a.Length > 0 ? a[0] : "4,6,8,10,12,14"); break;
            case "reads": Reads(I(0, 55)); break;
            case "concurrent": Concurrent(I(0, 55), I(1, 3), a.Length > 2 ? a[2] : "1,2,4,8,15"); break;
            case "persist": Persist(I(0, 55)); break;
            case "memory": Memory(I(0, 55)); break;
            case "endtoend": EndToEnd(I(0, 55), I(1, 100)); break;
            case "cutpar": CutParallelBench(I(0, 55), I(1, 4), a.Length > 2 ? a[2] : "1,2,4,8,16", I(3, 5)); break;
            case "batchcut": BatchCutBench(I(0, 55), I(1, 4), I(2, 3)); break;
            case "profile": ProfileRun(a.Length > 0 ? a[0] : "crack", I(1, 55), I(2, 1) == 1); break;
            default: throw new ArgumentException($"unknown mode {mode}");
        }
    }

    // ================================================================ helpers

    /// <summary>
    ///     Exercises every code path once on a small mesh so that just-in-time compilation
    ///     (tiering is off) does not land in a timed region.
    /// </summary>
    private static void WarmUp()
    {
        var saved = Console.Out;
        Console.SetOut(TextWriter.Null);
        try
        {
            for (var rep = 0; rep < 2; rep++)
            {
                var st = Build(5, persistent: true, mirror: true);
                var crack = new CrackPlane(5);
                while (crack.Row < crack.Steps) crack.Step(st);
                st.Mesh.Commit();
                var bis = new Bisection();
                bis.Refine(st, Bisection.InBall(st, 0.5, 0.5, 0.5, 0.3), 3);
                Erosion.Step(st, Bisection.InBall(st, 0.3, 0.3, 0.3, 0.2));
                var gen = new RandomEdits(1);
                for (var i = 0; i < 200; i++) gen.Duplicate(st);
                var snap = st.Mesh.Commit();
                st.Mesh.Rollback(snap);
                st.FlushMirror();
                GC.KeepAlive(st.Mirror!.Rebuild());
                GC.KeepAlive(FlatBuild.Build(Tet, st.Mesh.NodeCount, CurrentRows(st.Mesh), false));
                GC.KeepAlive(FlatBuild.Build(Tet, st.Mesh.NodeCount, CurrentRows(st.Mesh), true));
                var vals = new double[snap.TotalSlots];
                AssembleGapped(snap, st, vals);
            }
        }
        finally
        {
            Console.SetOut(saved);
        }
    }

    private static double TimeMs(Action a)
    {
        var sw = Stopwatch.StartNew();
        a();
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds;
    }

    private static void FullGC()
    {
        GC.Collect(2, GCCollectionMode.Forced, true, true);
        GC.WaitForPendingFinalizers();
        GC.Collect(2, GCCollectionMode.Forced, true, true);
    }

    private static double[] Samples(Func<double> f, int reps)
    {
        var s = new double[reps];
        for (var i = 0; i < reps; i++)
        {
            FullGC();
            s[i] = f();
        }

        Array.Sort(s);
        return s;
    }

    private static double Median(IList<double> x)
    {
        if (x.Count == 0) return double.NaN;
        var s = x.OrderBy(v => v).ToArray();
        return s.Length % 2 == 1 ? s[s.Length / 2] : 0.5 * (s[s.Length / 2 - 1] + s[s.Length / 2]);
    }

    private static double Quantile(IList<double> x, double q)
    {
        if (x.Count == 0) return double.NaN;
        var s = x.OrderBy(v => v).ToArray();
        var pos = q * (s.Length - 1);
        var lo = (int)Math.Floor(pos);
        var hi = Math.Min(s.Length - 1, lo + 1);
        return s[lo] + (pos - lo) * (s[hi] - s[lo]);
    }

    private static int[][] CurrentRows(IncrementalMesh m)
    {
        var rows = new int[m.ElementCount][];
        for (var e = 0; e < rows.Length; e++) rows[e] = m.Row(e).ToArray();
        return rows;
    }

    private static MeshState Build(int n, bool persistent, bool mirror)
    {
        var (rows, coords, nN) = Kuhn.Box(n);
        var flat = FlatBuild.Build(Tet, nN, rows, parallel: true);
        var mesh = new IncrementalMesh(flat, Tet, persistent);
        var st = new MeshState(mesh, coords, mirror ? new Mm3Mirror(rows) : null);
        FullGC();
        return st;
    }

    private static StreamWriter Csv(string name, string header)
    {
        var w = new StreamWriter(Path.Combine(DataDir, name));
        w.WriteLine(header);
        return w;
    }

    // ================================================================ E0 baseline reproduction

    private readonly struct V4;
    private readonly struct C4;

    /// <summary>
    ///     The MM3 paper's Table 12 protocol, verbatim in substance (figure_sources/Tables1213.cs):
    ///     quadrilateral grids, one appended element invalidates the caches, the first reverse
    ///     query (transpose) and the position-cache build are timed, minimum of seven samples.
    /// </summary>
    private static void Mm3Reference()
    {
        using var w = Csv("mm3ref.csv", "elements,transpose_ms_min,transpose_ms_median,positions_ms_min,positions_ms_median");
        Console.WriteLine($"{"elements",10} {"transpose (ms)",15} {"+positions (ms)",16}");
        foreach (var target in new[] { 10_000, 100_000, 1_000_000 })
        {
            var side = (int)Math.Round(Math.Sqrt(target));
            var cells = new int[side * side][];
            var t = 0;
            for (var r = 0; r < side; r++)
            for (var c = 0; c < side; c++)
            {
                int v00 = r * (side + 1) + c, v10 = v00 + 1, v01 = (r + 1) * (side + 1) + c, v11 = v01 + 1;
                cells[t++] = [v00, v10, v11, v01];
            }

            var extra = new List<int> { 0, 1, 2, 3 };
            var tS = Samples(() =>
            {
                var p = new O2M(cells.Length);
                foreach (var cell in cells) p.AppendElementCopy(cell);
                using var m = new M2M(p);
                m.GetElementsForNode(0);
                m.AppendElements(extra);
                return TimeMs(() => m.GetElementsForNode(0));
            }, 7);
            var tP = Samples(() =>
            {
                var p = new O2M(cells.Length);
                foreach (var cell in cells) p.AppendElementCopy(cell);
                using var m = new M2M(p);
                m.EnsurePositionCaches();
                m.AppendElements(extra);
                return TimeMs(() =>
                {
                    m.GetElementsForNode(0);
                    m.EnsurePositionCaches();
                });
            }, 7);
            Console.WriteLine($"{target,10:N0} {tS[0],15:F2} {tP[0],16:F2}   (medians {tS[3]:F2} / {tP[3]:F2})");
            w.WriteLine($"{target},{tS[0]:F4},{tS[3]:F4},{tP[0]:F4},{tP[3]:F4}");
        }
    }

    // ================================================================ phase profile

    private static void ProfileRun(string kind, int n, bool persistent)
    {
        if (!IncrementalMesh.Profile) throw new InvalidOperationException("set STEP2_PROFILE=1");
        var drv = new Driver { Kind = kind, N = n };
        drv.Steps = drv.DefaultSteps;
        drv.Crack = new CrackPlane(n);
        var st = Build(n, persistent, mirror: false);
        Array.Clear(IncrementalMesh.PhaseTicks);
        var e0 = st.Mesh.Edits;
        double total = 0;
        for (var s = 0; s < drv.Steps; s++)
        {
            var work = drv.Prepare(st, s);
            total += TimeMs(() =>
            {
                drv.Apply(st, work);
                st.Mesh.Commit();
            });
        }

        var edits = st.Mesh.Edits - e0;
        var ticks = IncrementalMesh.PhaseTicks.Sum();
        var phaseMs = ticks * 1000.0 / Stopwatch.Frequency;
        Console.WriteLine($"profile {kind} n={n} persistent={persistent}: {edits} edits, step total {total:F1} ms, " +
                          $"inside Apply {phaseMs:F1} ms ({1000 * phaseMs / edits:F2} us/edit), outside Apply " +
                          $"(workload logic, node creation, commit) {total - phaseMs:F1} ms");
        using var w = Csv($"profile_{kind}_{n}_{(persistent ? "p" : "t")}.csv", "phase,ms,us_per_edit,share");
        for (var i = 0; i < 8; i++)
        {
            var ms = IncrementalMesh.PhaseTicks[i] * 1000.0 / Stopwatch.Frequency;
            Console.WriteLine($"  {IncrementalMesh.PhaseNames[i],-40} {ms,9:F2} ms {1000 * ms / edits,7:F3} us/edit {100 * ms / phaseMs,6:F1} %");
            w.WriteLine($"\"{IncrementalMesh.PhaseNames[i]}\",{ms:F4},{1000 * ms / edits:F4},{ms / phaseMs:F4}");
        }
    }

    // ================================================================ E1 per-step maintenance

    private sealed class Driver
    {
        public string Kind = "crack";
        public int N;
        public CrackPlane? Crack;
        public Bisection Bis = new();
        public int Steps;

        public int DefaultSteps => Kind switch { "crack" => N + 1, _ => 40 };

        /// <summary>Untimed preparation (marking scans), returns the step's work list.</summary>
        public List<int>? Prepare(MeshState st, int s)
        {
            var r = 2.5 / N;
            var c = 0.15 + 0.7 * s / Math.Max(1, Steps - 1);
            return Kind switch
            {
                "bisect" => Bisection.InBall(st, c, 0.5, 0.5, r),
                "erode" => Bisection.InBall(st, c, 0.5, 0.5, r),
                _ => null
            };
        }

        public void Apply(MeshState st, List<int>? work)
        {
            switch (Kind)
            {
                case "crack": Crack!.Step(st); break;
                case "bisect": Bis.Refine(st, work!, maxLevel: 3); break;
                default: Erosion.Step(st, work!); break;
            }
        }
    }

    private static void Work(string kind, int n, int steps, int every, bool persistent)
    {
        var drv = new Driver { Kind = kind, N = n };
        drv.Steps = steps > 0 ? steps : drv.DefaultSteps;
        drv.Crack = new CrackPlane(n);
        var st = Build(n, persistent, mirror: true);
        var mesh = st.Mesh;
        var initial = mesh.Current;
        var workLists = new List<List<int>?>();
        Console.WriteLine($"{kind}: n={n}, {mesh.ElementCount:N0} elements, {mesh.NodeCount:N0} nodes, " +
                          $"{drv.Steps} steps, persistent={persistent}");

        var tag = persistent ? "p" : "t";
        using var w = Csv($"work_{kind}_{n}_{tag}.csv",
            "step,edits,new_nodes,elements,t_inc_ms,t_inc_min_ms,chunk_copies,dir_copies,relocations," +
            "t_mm3_ms,t_flat_core_par_ms,t_flat_all_seq_ms,t_flat_all_par_ms");
        var incTimes = new List<double>();
        var incFirst = new List<double>();
        var flatCore = new List<double>();
        var editCounts = new List<double>();
        var perEdit = new List<double>();
        var mm3 = new List<double>();
        var flatSeq = new List<double>();
        var flatPar = new List<double>();
        object? keep = null;

        for (var s = 0; s < drv.Steps; s++)
        {
            var work = drv.Prepare(st, s);
            workLists.Add(work);
            FullGC();
            var e0 = mesh.Edits;
            var n0 = mesh.NodeCount;
            var (c0, d0) = mesh.CopyCounters();
            var r0 = mesh.Relocations;
            var crackRow = drv.Crack.Row;

            // timed step (single shot), then exact repetitions through rollback (persistent mode)
            var before = mesh.Current;
            var t = TimeMs(() =>
            {
                drv.Apply(st, work);
                mesh.Commit();
            });
            var edits = mesh.Edits - e0;
            var nodes = mesh.NodeCount - n0;
            var (c1, d1) = mesh.CopyCounters();
            var r1 = mesh.Relocations;
            var tMin = t;
            if (persistent)
            {
                var after = mesh.Current;
                for (var rep = 0; rep < 4; rep++)
                {
                    mesh.Rollback(before);
                    st.DiscardMirrorLog();
                    drv.Crack.Row = crackRow;
                    tMin = Math.Min(tMin, TimeMs(() =>
                    {
                        drv.Apply(st, work);
                        mesh.Commit();
                    }));
                }

                GC.KeepAlive(after);
            }

            incTimes.Add(tMin);
            incFirst.Add(t);
            editCounts.Add(edits);
            if (edits > 0) perEdit.Add(1000.0 * tMin / edits);
            st.FlushMirror();

            double tm = double.NaN, tfc = double.NaN, tfs = double.NaN, tfp = double.NaN;
            if (s % every == 0 || s == drv.Steps - 1)
            {
                tm = Samples(() => TimeMs(() => keep = st.Mirror!.Rebuild()), 3)[0];
                keep = null;
                var rowsNow = CurrentRows(mesh);
                var nN = mesh.NodeCount;
                tfc = Samples(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN, rowsNow, true, false)), 3)[0];
                tfs = Samples(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN, rowsNow, false, true)), 3)[0];
                tfp = Samples(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN, rowsNow, true, true)), 3)[0];
                keep = null;
                mm3.Add(tm);
                flatCore.Add(tfc);
                flatSeq.Add(tfs);
                flatPar.Add(tfp);
            }

            w.WriteLine($"{s},{edits},{nodes},{mesh.ElementCount},{t:F4},{tMin:F4},{c1 - c0},{d1 - d0},{r1 - r0}," +
                        $"{tm:F3},{tfc:F3},{tfs:F3},{tfp:F3}");
        }

        GC.KeepAlive(keep);

        // pass B: the whole sequence back to back from the initial version, collections included
        st.Mirror = null;
        if (persistent)
        {
            mesh.Rollback(initial);
        }
        else
        {
            st = Build(n, persistent: false, mirror: false);
            mesh = st.Mesh;
        }

        drv.Crack = new CrackPlane(n);
        var gcBefore = GC.CollectionCount(0); // every collection includes generation 0
        FullGC();
        var hot = TimeMs(() =>
        {
            for (var s = 0; s < drv.Steps; s++)
            {
                drv.Apply(st, workLists[s]);
                mesh.Commit();
            }
        });
        var gcs = GC.CollectionCount(0) - gcBefore - 2; // minus the two forced by FullGC
        File.AppendAllText(Path.Combine(DataDir, $"work_{kind}_{n}_{tag}_hot.txt"),
            $"steps={drv.Steps} total_ms={hot:F3} per_step_ms={hot / drv.Steps:F4} " +
            $"edits={editCounts.Sum():F0} us_per_edit={1000 * hot / Math.Max(1, editCounts.Sum()):F4} gcs={gcs}" +
            Environment.NewLine);

        var line = $"{kind} n={n} ({mesh.ElementCount:N0} el.) {tag}: hot sequence {hot / drv.Steps:F3} ms/step " +
                   $"({1000 * hot / Math.Max(1, editCounts.Sum()):F2} us/edit, {gcs} collections); " +
                   $"cold median step {Median(incTimes):F3} ms " +
                   $"[p10 {Quantile(incTimes, 0.1):F3}, p90 {Quantile(incTimes, 0.9):F3}], single-shot median " +
                   $"{Median(incFirst):F3} ms, sequence total single-shot {incFirst.Sum():F1} ms; " +
                   $"median edits/step {Median(editCounts):F0}, median {Median(perEdit):F2} us/edit; " +
                   $"rebuild median: MM3 {Median(mm3):F1} ms, flat core par {Median(flatCore):F1} ms, " +
                   $"flat all seq {Median(flatSeq):F1} ms, flat all par {Median(flatPar):F1} ms; " +
                   $"speedup vs MM3 {Median(mm3) / Median(incTimes):F0}x, vs flat core par " +
                   $"{Median(flatCore) / Median(incTimes):F0}x, vs flat all par {Median(flatPar) / Median(incTimes):F0}x";
        Console.WriteLine(line);
        File.AppendAllText(Path.Combine(DataDir, "work_summary.txt"), line + Environment.NewLine);
    }
}

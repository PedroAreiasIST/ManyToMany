// Bench2.cs — crossover, interleaved, reads, concurrency, persistence, memory, end-to-end.

using System.Diagnostics;
using Numerical;

namespace Step2.Bench;

internal static partial class Bench
{
    // ================================================================ E2 crossover

    /// <summary>Order experiment for the crossover baseline: every sample of the tuned rebuild of the
    /// MM3 set, before and after the MM3 rebuild, each preceded by a full collection.</summary>
    private static void BaselineOrder(int n)
    {
        var st = Build(n, persistent: true, mirror: true);
        var rows0 = CurrentRows(st.Mesh);
        var nN0 = st.Mesh.NodeCount;
        object? keep = null;
        string S(Func<double> f) => string.Join(" ", Enumerable.Range(0, 5).Select(_ => { FullGC(); return f().ToString("F1"); }));
        Console.WriteLine("flat core par, first:   " + S(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN0, rows0, true, false))));
        Console.WriteLine("MM3 rebuild:            " + S(() => TimeMs(() => keep = st.Mirror!.Rebuild())));
        Console.WriteLine("flat core par, after:   " + S(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN0, rows0, true, false))));
        Console.WriteLine("flat all par:           " + S(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN0, rows0, true, true))));
        Console.WriteLine("flat core par, last:    " + S(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN0, rows0, true, false))));
        GC.KeepAlive(keep);
    }

    private static void Crossover(int n)
    {
        var st = Build(n, persistent: true, mirror: true);
        var mesh = st.Mesh;
        var nE = mesh.ElementCount;
        var baseSnap = mesh.Current;
        object? keep = null;

        var rows0 = CurrentRows(mesh);
        var nN0 = mesh.NodeCount;
        // warm-up: the first calls of a rebuild in the process run 2.4 times slower (baseline-order); one
        // untimed pass of every baseline puts all four in the state the workload measurements see
        for (var i = 0; i < 3; i++)
        {
            keep = st.Mirror!.Rebuild();
            keep = FlatBuild.Build(Tet, nN0, rows0, true, false);
            keep = FlatBuild.Build(Tet, nN0, rows0, false, true);
            keep = FlatBuild.Build(Tet, nN0, rows0, true, true);
        }

        var tFc = Samples(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN0, rows0, true, false)), 3)[0];
        var tFs = Samples(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN0, rows0, false, true)), 3)[0];
        var tFp = Samples(() => TimeMs(() => keep = FlatBuild.Build(Tet, nN0, rows0, true, true)), 3)[0];
        var tM = Samples(() => TimeMs(() => keep = st.Mirror!.Rebuild()), 3)[0];
        keep = null;
        Console.WriteLine($"crossover n={n}: {nE:N0} elements; rebuild MM3 {tM:F1} ms, flat core par {tFc:F1} ms, " +
                          $"flat all seq {tFs:F1} ms, flat all par {tFp:F1} ms");

        using var w = Csv($"crossover_{n}.csv",
            "batch,fraction,t_inc_ms,us_per_edit,t_mm3_ms,t_flat_core_par_ms,t_flat_all_seq_ms,t_flat_all_par_ms");
        var sizes = new List<int>();
        for (var b = 1.0; b <= nE; b *= Math.Sqrt(10)) sizes.Add((int)Math.Round(b));
        sizes.Add(nE);
        foreach (var b in sizes.Distinct())
        {
            var best = double.MaxValue;
            for (var trial = 0; trial < 3; trial++)
            {
                mesh.Rollback(baseSnap);
                st.DiscardMirrorLog();
                var gen = new RandomEdits(777);
                FullGC();
                best = Math.Min(best, TimeMs(() =>
                {
                    for (var i = 0; i < b; i++) gen.Duplicate(st);
                    mesh.Commit();
                }));
            }

            st.DiscardMirrorLog();
            Console.WriteLine($"  batch {b,9:N0} ({100.0 * b / nE,8:F4} %): {best,10:F3} ms  ({1000 * best / b,7:F2} us/edit)");
            w.WriteLine($"{b},{(double)b / nE:E4},{best:F4},{1000 * best / b:F4},{tM:F3},{tFc:F3},{tFs:F3},{tFp:F3}");
        }

        GC.KeepAlive(keep);
    }

    // ================================================================ E3 interleaved bisection on MM3

    /// <summary>Maubach bisection on an MM3 M2M: every edge query after an edit rebuilds the transpose.</summary>
    private sealed class M2MBisection(M2M m, List<double> x, List<int> tag, List<int> level)
    {
        public long Bisections;
        private int _nodes = x.Count / 3;

        private (int, int) RefEdge(int e)
        {
            var r = m[e];
            return (r[0], r[tag[e]]);
        }

        public void Bisect(int e, int depth = 0)
        {
            if (depth > 256) throw new InvalidOperationException("recursion");
            var (a, b) = RefEdge(e);
            List<int> patch;
            while (true)
            {
                patch = m.GetElementsWithNodes([a, b]);
                var inc = -1;
                foreach (var g in patch)
                {
                    var (ga, gb) = RefEdge(g);
                    if (!((ga == a && gb == b) || (ga == b && gb == a)))
                    {
                        inc = g;
                        break;
                    }
                }

                if (inc < 0) break;
                Bisect(inc, depth + 1);
            }

            var z = _nodes++;
            for (var c = 0; c < 3; c++) x.Add(0.5 * (x[3 * a + c] + x[3 * b + c]));
            foreach (var g in patch)
            {
                var row = m[g].ToArray();
                var d = tag[g];
                var child1 = (int[])row.Clone();
                child1[d] = z;
                var child2 = new int[4];
                for (var p = 0; p < d; p++) child2[p] = row[p + 1];
                child2[d] = z;
                for (var p = d + 1; p < 4; p++) child2[p] = row[p];
                var nd = d > 1 ? d - 1 : 3;
                var lv = level[g] + 1;
                m.ReplaceElement(g, child1.ToList());
                tag[g] = nd;
                level[g] = lv;
                m.AppendElement(child2.ToList());
                tag.Add(nd);
                level.Add(lv);
                Bisections++;
            }
        }
    }

    private static void Interleaved(string sizes)
    {
        {
            // warm the MM3 path (the shared warm-up covers the incremental path)
            var (wr, wc, _) = Kuhn.Box(3);
            using var wm = new M2M(wr.Select(r => r.ToList()).ToList());
            var wb = new M2MBisection(wm, wc, Enumerable.Repeat(3, wr.Length).ToList(),
                Enumerable.Repeat(0, wr.Length).ToList());
            for (var e = 0; e < 20; e++) wb.Bisect(e);
        }

        using var w = Csv("interleaved.csv",
            "n,elements,marked,bisections,t_inc_ms,t_mm3_ms,us_per_bisection_inc,us_per_bisection_mm3,identical");
        foreach (var n in sizes.Split(',').Select(int.Parse))
        {
            // incremental
            var st = Build(n, persistent: true, mirror: false);
            var marked = Bisection.InBall(st, 0.5, 0.5, 0.5, 0.25);
            var bis = new Bisection();
            FullGC();
            var tInc = TimeMs(() =>
            {
                foreach (var e in marked) bis.Bisect(st, e);
                st.Mesh.Commit();
            });

            // MM3 with lazy transpose rebuild
            var (rows, coords, _) = Kuhn.Box(n);
            using var m = new M2M(rows.Select(r => r.ToList()).ToList());
            var tag = Enumerable.Repeat(3, rows.Length).ToList();
            var level = Enumerable.Repeat(0, rows.Length).ToList();
            var mb = new M2MBisection(m, coords, tag, level);
            FullGC();
            var tMm3 = TimeMs(() =>
            {
                foreach (var e in marked) mb.Bisect(e);
            });

            // both runs take identical decisions: compare the final element rows
            var identical = m.Count == st.Mesh.ElementCount;
            for (var e = 0; identical && e < m.Count; e++)
                identical = m[e].SequenceEqual(st.Mesh.Row(e).ToArray());

            Console.WriteLine($"interleaved n={n}: {rows.Length:N0} el., {bis.Bisections:N0} bisections; " +
                              $"incremental {tInc:F1} ms ({1000 * tInc / bis.Bisections:F2} us/bisection), " +
                              $"MM3 {tMm3:F1} ms ({1000 * tMm3 / mb.Bisections:F1} us/bisection); identical={identical}");
            w.WriteLine($"{n},{rows.Length},{marked.Count},{bis.Bisections},{tInc:F3},{tMm3:F3}," +
                        $"{1000 * tInc / bis.Bisections:F4},{1000 * tMm3 / mb.Bisections:F4},{identical}");
        }
    }

    // ================================================================ E4a single-thread reads

    private readonly struct KNode;
    private readonly struct KTet;

    private static void Reads(int n)
    {
        var st = Build(n, persistent: true, mirror: false);
        var snap = st.Mesh.Current;
        var nE = snap.ElementCount;
        var nN = snap.NodeCount;
        const int probes = 2_000_000;
        var rng = new Random(17);
        var eP = new int[probes];
        var vP = new int[probes];
        for (var i = 0; i < probes; i++)
        {
            eP[i] = rng.Next(nE);
            vP[i] = rng.Next(nN);
        }

        var results = new List<(string Name, double MopsPerS)>();
        void Measure(string name, int ops, Func<long> body)
        {
            body();
            var ms = Samples(() => TimeMs(() => GC.KeepAlive(body())), 7)[0];
            results.Add((name, ops / ms / 1000.0));
        }

        Measure("step2 snapshot: element row", probes, () =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += snap.Row(eP[i])[3];
            return acc;
        });
        Measure("step2 snapshot: transpose row", probes, () =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += snap.TransposePairs(vP[i]).Length;
            return acc;
        });
        Measure("step2 snapshot: facet neighbours (4/probe)", 4 * probes, () =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++)
                for (var f = 0; f < 4; f++)
                    acc += snap.FacetNeighbour(eP[i], f);
            return acc;
        });
        Measure("step2 snapshot: destination block (16/probe)", probes, () =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += snap.Destinations(eP[i])[15];
            return acc;
        });

        // MM3 on the same mesh
        var rows = CurrentRows(st.Mesh);
        var o2m = new O2M(rows.Length);
        foreach (var r in rows) o2m.AppendElementCopy(r);
        var tr = o2m.Transpose();
        using var m2m = new M2M(o2m);
        m2m.EnsurePositionCaches();
        Measure("MM3 O2M: element row (no lock)", probes, () =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += o2m[eP[i]][3];
            return acc;
        });
        Measure("MM3 O2M: transpose row (no lock)", probes, () =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += tr[vP[i]].Length;
            return acc;
        });
        Measure("MM3 M2M: transpose row, per-call read lock", probes, () =>
        {
            long acc = 0;
            ReadOnlySpanAction<int> f = s => acc += s.Length; // one delegate, no per-call allocation
            for (var i = 0; i < probes; i++) m2m.WithElementsForNodeSpan(vP[i], f);
            return acc;
        });
        Measure("MM3 M2M: GetElementsForNode (per-call lock + copy)", probes, () =>
        {
            long acc = 0;
            for (var i = 0; i < probes; i++) acc += m2m.GetElementsForNode(vP[i]).Count;
            return acc;
        });

        // MM3 typed facet neighbours (shared count 3), computed per query
        var topo = new Topology<TypeMap<KNode, KTet>>();
        for (var v = 0; v < nN; v++) topo.Add<KNode>();
        foreach (var r in rows) topo.Add<KTet, KNode>(r);
        const int fprobes = 200_000;
        Measure("MM3 Topology: facet neighbours (shared count 3)", fprobes, () =>
        {
            long acc = 0;
            for (var i = 0; i < fprobes; i++) acc += topo.GetEntitiesWithSharedCount<KTet, KNode>(eP[i], 3).Count;
            return acc;
        });
        topo.Dispose();

        using var w = Csv($"reads_{n}.csv", "operation,Mops_per_s");
        foreach (var (name, mops) in results)
        {
            Console.WriteLine($"  {name,-52} {mops,9:F2} M ops/s");
            w.WriteLine($"\"{name}\",{mops:F4}");
        }
    }

    // ================================================================ E4b concurrent readers during writes

    private static void Concurrent(int n, int seconds, string readerList)
    {
        using var w = Csv($"concurrent_{n}.csv",
            "impl,readers,reads_per_s,window_s,writer_steps,writer_step_median_ms,writer_step_p90_ms");
        foreach (var readers in readerList.Split(',').Select(int.Parse))
        {
            ConcurrentStep2(n, seconds, readers, w);
            ConcurrentMm3(n, seconds, readers, w);
        }
    }

    private static void ConcurrentStep2(int n, int seconds, int readers, StreamWriter w)
    {
        var st = Build(n, persistent: true, mirror: false);
        var mesh = st.Mesh;
        var baseSnap = mesh.Current;
        var crack = new CrackPlane(n);
        var nN0 = mesh.NodeCount;
        var nE0 = mesh.ElementCount;
        var stop = 0;
        var counts = new long[readers];
        var threads = new Thread[readers];
        for (var r = 0; r < readers; r++)
        {
            var id = r;
            threads[r] = new Thread(() =>
            {
                var rng = new Random(100 + id);
                long acc = 0, ops = 0;
                while (Volatile.Read(ref stop) == 0)
                {
                    var s = mesh.Current;
                    for (var i = 0; i < 1000; i++)
                    {
                        acc += s.TransposePairs(rng.Next(nN0)).Length;
                        acc += s.FacetNeighbour(rng.Next(nE0), i & 3);
                    }

                    ops += 2000;
                }

                counts[id] = ops;
                GC.KeepAlive(acc);
            }) { IsBackground = true };
        }

        FullGC();
        var steps = new List<double>();
        var sw = Stopwatch.StartNew();
        foreach (var t in threads) t.Start();
        while (sw.Elapsed.TotalSeconds < seconds)
        {
            if (crack.Row >= crack.Steps)
            {
                mesh.Rollback(baseSnap);
                crack.Row = 0;
            }

            var t0 = sw.Elapsed.TotalMilliseconds;
            crack.Step(st);
            mesh.Commit();
            steps.Add(sw.Elapsed.TotalMilliseconds - t0);
        }

        Volatile.Write(ref stop, 1);
        foreach (var t in threads) t.Join();
        var elapsed = sw.Elapsed.TotalSeconds;
        var rate = counts.Sum() / elapsed;
        Console.WriteLine($"  step2   readers={readers,2}: {rate / 1e6,9:F1} M reads/s over {elapsed:F2} s; " +
                          $"writer {steps.Count} steps, median {Median(steps):F2} ms, p90 {Quantile(steps, 0.9):F2} ms");
        w.WriteLine($"step2,{readers},{rate:F0},{elapsed:F3},{steps.Count},{Median(steps):F4},{Quantile(steps, 0.9):F4}");
    }

    private static void ConcurrentMm3(int n, int seconds, int readers, StreamWriter w)
    {
        var (rows, coords, nN0) = Kuhn.Box(n);
        var nE0 = rows.Length;
        using var m = new M2M(rows.Select(r => r.ToList()).ToList());
        m.EnsurePositionCaches();
        var x = coords;
        var stop = 0;
        var counts = new long[readers];
        var threads = new Thread[readers];
        for (var r = 0; r < readers; r++)
        {
            var id = r;
            threads[r] = new Thread(() =>
            {
                var rng = new Random(100 + id);
                long acc = 0, ops = 0;
                ReadOnlySpanAction<int> f = s => acc += s.Length; // one delegate, no per-call allocation
                while (Volatile.Read(ref stop) == 0)
                {
                    for (var i = 0; i < 1000; i++)
                        m.WithElementsForNodeSpan(rng.Next(nN0), f);
                    ops += 1000;
                }

                counts[id] = ops;
                GC.KeepAlive(acc);
            }) { IsBackground = true };
        }

        FullGC();
        var steps = new List<double>();
        var ic = n / 2;
        var xPlane = (double)ic / n;
        var row = 0;
        var nodes = nN0;
        var sw = Stopwatch.StartNew();
        foreach (var t in threads) t.Start();
        while (sw.Elapsed.TotalSeconds < seconds && row <= n)
        {
            var t0 = sw.Elapsed.TotalMilliseconds;
            var edits = new SortedDictionary<int, List<int>>();
            for (var k = 0; k <= n; k++)
            {
                var v = Kuhn.Node(n, ic, row, k);
                var dup = -1;
                foreach (var e in m.GetElementsForNode(v))
                {
                    var r = edits.TryGetValue(e, out var pending) ? pending : m[e].ToList();
                    double cx = 0;
                    foreach (var u in r) cx += x[3 * u];
                    if (cx / 4 <= xPlane) continue;
                    if (dup < 0)
                    {
                        dup = nodes++;
                        x.Add(x[3 * v]);
                        x.Add(x[3 * v + 1]);
                        x.Add(x[3 * v + 2]);
                    }

                    r[r.IndexOf(v)] = dup;
                    edits[e] = r;
                }
            }

            foreach (var (e, r) in edits) m.ReplaceElement(e, r);
            m.Synchronize();
            steps.Add(sw.Elapsed.TotalMilliseconds - t0);
            row++;
        }

        Volatile.Write(ref stop, 1);
        foreach (var t in threads) t.Join();
        var elapsed = sw.Elapsed.TotalSeconds;
        var rate = counts.Sum() / elapsed;
        Console.WriteLine($"  MM3 M2M readers={readers,2}: {rate / 1e6,9:F1} M reads/s over {elapsed:F2} s; " +
                          $"writer {steps.Count} steps, median {Median(steps):F2} ms, p90 {Quantile(steps, 0.9):F2} ms");
        w.WriteLine($"mm3,{readers},{rate:F0},{elapsed:F3},{steps.Count},{Median(steps):F4},{Quantile(steps, 0.9):F4}");
    }

    // ================================================================ E5 persistence costs

    private static void Persist(int n)
    {
        using var w = Csv($"persist_{n}.csv", "quantity,value,unit");
        void Out(string q, double v, string unit)
        {
            Console.WriteLine($"  {q,-60} {v,14:F3} {unit}");
            w.WriteLine($"\"{q}\",{v:F6},{unit}");
        }

        var st = Build(n, persistent: true, mirror: false);
        var mesh = st.Mesh;
        Console.WriteLine($"persist n={n}: {mesh.ElementCount:N0} elements");

        // publication and rollback in isolation
        var commitUs = Samples(() => TimeMs(() =>
        {
            for (var i = 0; i < 10_000; i++) mesh.Commit();
        }), 5)[0] * 1000.0 / 10_000;
        Out("commit (publication, no edits)", commitUs, "us");
        var baseSnap = mesh.Current;
        var rollbackUs = Samples(() => TimeMs(() =>
        {
            for (var i = 0; i < 10_000; i++) mesh.Rollback(baseSnap);
        }), 5)[0] * 1000.0 / 10_000;
        Out("rollback to a published version", rollbackUs, "us");

        // first edit after a publication (directory and chunk copies) vs a transient mesh
        var gen = new RandomEdits(5);
        var firstEdit = Samples(() =>
        {
            mesh.Rollback(baseSnap);
            return TimeMs(() => gen.Duplicate(st));
        }, 31)[0];
        Out("first edit after publication (persistent)", firstEdit * 1000, "us");
        var trans = Build(n, persistent: false, mirror: false);
        var genT = new RandomEdits(5);
        var firstEditT = Samples(() => TimeMs(() => genT.Duplicate(trans)), 31)[0];
        Out("one edit (transient mesh)", firstEditT * 1000, "us");

        // retained versions: a complete crack sequence keeping every version alive
        mesh.Rollback(baseSnap);
        FullGC();
        var mem0 = GC.GetTotalMemory(true);
        var kept = new List<Snapshot>();
        var crack = new CrackPlane(n);
        long edits = 0;
        for (var s = 0; s < crack.Steps; s++)
        {
            var e0 = mesh.Edits;
            crack.Step(st);
            kept.Add(mesh.Commit());
            edits += mesh.Edits - e0;
        }

        FullGC();
        var memKept = GC.GetTotalMemory(true);
        var latest = mesh.Current;
        kept.Clear();
        FullGC();
        var memLatest = GC.GetTotalMemory(true);
        GC.KeepAlive(latest);
        var perVersion = (memKept - memLatest) / (double)(crack.Steps - 1);
        Out($"versions retained ({crack.Steps} crack steps, {edits} edits)", crack.Steps, "versions");
        Out("live memory, latest version only", memLatest / 1e6, "MB");
        Out("live memory, all versions retained", memKept / 1e6, "MB");
        Out("memory per retained version", perVersion / 1e6, "MB");
        Out("memory per retained edit", (memKept - memLatest) / (double)edits, "bytes");
        Out("growth of the latest version over the start", (memLatest - mem0) / 1e6, "MB");
    }

    // ================================================================ E6 memory footprint

    private static void Memory(int n)
    {
        using var w = Csv($"memory_{n}.csv", "representation,bytes,bytes_per_element");
        void Out(string name, long bytes, int nE)
        {
            Console.WriteLine($"  {name,-58} {bytes / 1e6,10:F1} MB  {bytes / (double)nE,8:F1} B/element");
            w.WriteLine($"\"{name}\",{bytes},{bytes / (double)nE:F2}");
        }

        FullGC();
        var m0 = GC.GetTotalMemory(true);
        var (rows, coords, nN) = Kuhn.Box(n);
        FullGC();
        var mRows = GC.GetTotalMemory(true);
        var nE = rows.Length;
        Console.WriteLine($"memory n={n}: {nE:N0} elements, {nN:N0} nodes");
        Out("element rows (int[][]), input", mRows - m0, nE);

        var flat = FlatBuild.Build(Tet, nN, rows, true);
        FullGC();
        var mFlat = GC.GetTotalMemory(true);
        Out("flat derived structures (rebuild output)", mFlat - mRows, nE);

        var mesh = new IncrementalMesh(flat, Tet, persistent: true);
        flat = null;
        FullGC();
        var mMesh = GC.GetTotalMemory(true);
        Out("incremental mesh (all structures, chunked, with rows)", mMesh - mRows, nE);
        GC.KeepAlive(mesh);
        mesh = null;
        FullGC();
        var mAfter = GC.GetTotalMemory(true);

        var o2m = new O2M(rows.Length);
        foreach (var r in rows) o2m.AppendElementCopy(r);
        FullGC();
        var mO2M = GC.GetTotalMemory(true);
        Out("MM3 O2M primary relation", mO2M - mAfter, nE);
        var t = o2m.Transpose();
        var ep = O2M.GetElementPositions(o2m, t);
        var cl = O2M.GetCliques(o2m, t);
        var nn = t * o2m;
        FullGC();
        var mDerived = GC.GetTotalMemory(true);
        Out("MM3 derived: transpose + positions + cliques + node-node", mDerived - mO2M, nE);
        Out("MM3 total (primary + derived)", mDerived - mAfter, nE);
        GC.KeepAlive((o2m, t, ep, cl, nn, coords));
    }
}

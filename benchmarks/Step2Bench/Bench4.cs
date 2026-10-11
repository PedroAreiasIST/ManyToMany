// Bench4.cs — concurrent cutting of independent edge groups.
//
// Copyright (C) 2026 Pedro Miguel de Almeida Areias
//
// This file is part of ManyToMany.
//
// ManyToMany is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ManyToMany is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ManyToMany. If not, see <https://www.gnu.org/licenses/>.
//
//
//   cutpar n g threads reps   Kuhn box of 6 n^3 tetrahedra; g^3 groups, group c marks the
//                             edges crossing a sphere of radius 0.35/g centred in grid cell c
//                             (disjoint closed shells); sequential Cut of the concatenated
//                             groups vs CutParallel with each thread count. Every sample
//                             starts from the same committed version (rollback), so the
//                             copy-on-write of the first edits after publication is included.
//   batchcut n g reps         the same marks: sequential local Cut against the batch cutter
//                             (Numerical.MeshRefinement.Refine, the C# translation of the published
//                             Fortran templates, which sorts all edges and builds a new mesh), for all
//                             g^3 groups and for the first group alone. The SimplexMesh input of the
//                             batch cutter is built once, outside the timing.
//   cutgc n g reps            the sequential cut of the cutpar marks in one process, with and without a second mesh
//                             alive and in two mark orders: the time, the collections and their pauses (the pause
//                             depends on what the process holds, the time without the pauses does not).
//
// Every sample records the stop-the-world pause of the collector inside the timed region (GcSamples): the logs and the
// CSV files give the time, the pause at the minimum, and the time without the pauses.

using Numerical;
using Step2;

namespace Step2.Bench;

internal static partial class Bench
{
    private static void CutParallelBench(int n, int g, string threadList, int reps)
    {
        var (rows, coords, nN) = Kuhn.Box(n);
        var flat = FlatBuild.Initial(Tet, nN, rows);
        var m = new IncrementalMesh(flat, Tet, persistent: true);
        var baseSnap = m.Current;
        var edges = CutVerify.Edges(m);
        var groups = new List<(int, int)>[g * g * g];
        for (var c = 0; c < groups.Length; c++) groups[c] = new List<(int, int)>();
        var rad = 0.35 / g;
        double Phi(int v, int c)
        {
            int i = c / (g * g), j = c / g % g, k = c % g;
            double dx = coords[3 * v] - (i + 0.5) / g, dy = coords[3 * v + 1] - (j + 0.5) / g,
                dz = coords[3 * v + 2] - (k + 0.5) / g;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz) - rad;
        }

        int CellOf(int v)
        {
            int Ix(int d) => Math.Min(g - 1, (int)(coords[3 * v + d] * g));
            return (Ix(0) * g + Ix(1)) * g + Ix(2);
        }

        foreach (var (a, b) in edges)
        {
            var c = CellOf(a);
            if (Phi(a, c) * Phi(b, c) < 0) groups[c].Add((a, b));
        }

        var list = groups.Select(l => (IReadOnlyList<(int, int)>)l).ToList();
        var concat = groups.SelectMany(l => l).ToList();
        Console.WriteLine($"  n={n}: {rows.Length} tetrahedra, {groups.Length} groups, {concat.Count} marked edges");

        using var w = Csv($"cutpar_{n}_{g}.csv",
            "impl,threads,ms_min,ms_median,plan_ms,reserve_ms,apply_ms,edges,shell_elements,groups_applied,merged,reserved_slots," +
            GcSamples.CsvHead);

        // sequential reference
        var seq = new GcSamples();
        CutRecord? rec = null;
        for (var r = 0; r < reps; r++)
        {
            m.Rollback(baseSnap);
            FullGC();
            seq.Add(TimeMs(() => rec = Cutter.Cut(m, concat)));
        }

        var nodesSeq = m.NodeCount;
        var elementsSeq = m.Current.ElementCount;
        m.Commit();
        var elementsAfter = m.ElementCount;
        Console.WriteLine($"  sequential: min {seq.Min:F2} ms, median {seq.Median:F2} ms, {seq.Text} " +
                          $"({rec!.Edges.Length} edges, {rec.ShellElements} shell elements)");
        w.WriteLine($"sequential,1,{seq.Min:F4},{seq.Median:F4},,,,{rec.Edges.Length},{rec.ShellElements},1,0,,{seq.CsvTail}");

        foreach (var t in threadList.Split(',').Select(int.Parse))
        {
            var s = new GcSamples();
            var best = (Plan: double.MaxValue, Reserve: 0.0, Apply: 0.0);
            ParallelCutResult? res = null;
            for (var r = 0; r < reps; r++)
            {
                m.Rollback(baseSnap);
                FullGC();
                var wall = TimeMs(() => res = Cutter.CutParallel(m, list, null, t));
                s.Add(wall);
                if (wall <= s.Min) best = (res!.PlanMs, res.ReserveMs, res.ApplyMs);
            }

            if (m.NodeCount != nodesSeq || m.ElementCount != elementsAfter)
                throw new Exception($"parallel cut: {m.NodeCount} nodes, {m.ElementCount} elements vs " +
                                    $"{nodesSeq}, {elementsAfter}");
            var shell = res!.Groups.Sum(x => x.ShellElements);
            Console.WriteLine($"  parallel {t,2} threads: min {s.Min:F2} ms (plan {best.Plan:F2}, reserve " +
                              $"{best.Reserve:F2}, apply {best.Apply:F2}), median {s.Median:F2} ms, speed-up " +
                              $"{seq.Min / s.Min:F2}, {s.Text}, speed-up without pauses {seq.FreeMin / s.FreeMin:F2}; " +
                              $"{res.Groups.Length} groups, {res.Merged} merged");
            w.WriteLine($"parallel,{t},{s.Min:F4},{s.Median:F4},{best.Plan:F4},{best.Reserve:F4},{best.Apply:F4}," +
                        $"{rec.Edges.Length},{shell},{res.Groups.Length},{res.Merged},{res.ReservedSlots},{s.CsvTail}");
        }

        GC.KeepAlive(elementsSeq);
    }

    private static void BatchCutBench(int n, int g, int reps)
    {
        var (rows, coords, nN) = Kuhn.Box(n);
        var flat = FlatBuild.Initial(Tet, nN, rows);
        var m = new IncrementalMesh(flat, Tet, persistent: true);
        var baseSnap = m.Current;
        var edges = CutVerify.Edges(m);
        var rad = 0.35 / g;
        var all = new List<(int, int)>();
        var first = new List<(int, int)>();
        foreach (var (a, b) in edges)
        {
            int Ix(int v, int d) => Math.Min(g - 1, (int)(coords[3 * v + d] * g));
            var c = (Ix(a, 0) * g + Ix(a, 1)) * g + Ix(a, 2);
            int i = c / (g * g), j = c / g % g, k = c % g;
            double Phi(int v)
            {
                double dx = coords[3 * v] - (i + 0.5) / g, dy = coords[3 * v + 1] - (j + 0.5) / g,
                    dz = coords[3 * v + 2] - (k + 0.5) / g;
                return Math.Sqrt(dx * dx + dy * dy + dz * dz) - rad;
            }

            if (Phi(a) * Phi(b) >= 0) continue;
            all.Add((a, b));
            if (c == 0) first.Add((a, b));
        }

        var sm = new SimplexMesh();
        sm.WithBatch(() =>
        {
            for (var v = 0; v < nN; v++) sm.AddNode(v);
            foreach (var r in rows) sm.Add<Tet4, Node>(r);
        });
        Console.WriteLine($"  n={n}: {rows.Length} tetrahedra, {g * g * g} groups, {all.Count} marked edges " +
                          $"({first.Count} in the first group)");
        using var w = Csv($"batchcut_{n}_{g}.csv", "impl,marks,edges,ms_min,ms_median,elements_out," + GcSamples.CsvHead);
        // the batch cutter returns element rows only: the structures the local cut maintains are then
        // rebuilt from scratch (all structures, parallel), timed separately on the extracted rows
        var stdout = Console.Out;
        foreach (var (name, marks) in new[] { ("all", all), ("one", first) })
        {
            var loc = new GcSamples();
            var bat = new GcSamples();
            var reb = new GcSamples();
            int outLocal = 0, outBatch = 0;
            for (var r = 0; r < reps; r++)
            {
                m.Rollback(baseSnap);
                FullGC();
                loc.Add(TimeMs(() => Cutter.Cut(m, marks)));
                outLocal = m.ElementCount;
                FullGC();
                Console.SetOut(TextWriter.Null);
                SimplexMesh? o = null;
                bat.Add(TimeMs(() => (o, _) = MeshRefinement.Refine(sm, marks)));
                Console.SetOut(stdout);
                outBatch = o!.Count<Tet4>();
                var rowsOut = new int[outBatch][];
                for (var e = 0; e < outBatch; e++) rowsOut[e] = o.NodesOf<Tet4, Node>(e).ToArray();
                var nOut = o.Count<Node>();
                FullGC();
                reb.Add(TimeMs(() => GC.KeepAlive(FlatBuild.Build(Tet, nOut, rowsOut, true, true))));
            }

            if (outLocal != outBatch)
                throw new Exception($"batchcut {name}: local {outLocal} elements, batch {outBatch}");
            Console.WriteLine($"  {name}: {marks.Count} edges, local min {loc.Min:F2} ms (median {loc.Median:F2}; {loc.Text}), " +
                              $"batch min {bat.Min:F2} ms (median {bat.Median:F2}; {bat.Text}), " +
                              $"rebuild min {reb.Min:F2} ms ({reb.Text}), " +
                              $"(batch + rebuild) / local {(bat.Min + reb.Min) / loc.Min:F1}, without pauses " +
                              $"{(bat.FreeMin + reb.FreeMin) / loc.FreeMin:F1}; {outLocal} elements");
            w.WriteLine($"local,{name},{marks.Count},{loc.Min:F4},{loc.Median:F4},{outLocal},{loc.CsvTail}");
            w.WriteLine($"batch,{name},{marks.Count},{bat.Min:F4},{bat.Median:F4},{outBatch},{bat.CsvTail}");
            w.WriteLine($"rebuild,{name},{marks.Count},{reb.Min:F4},{reb.Median:F4},{outBatch},{reb.CsvTail}");
        }
    }

    private static void CutGcBench(int n, int g, int reps)
    {
        var (rows, coords, nN) = Kuhn.Box(n);
        var flat = FlatBuild.Initial(Tet, nN, rows);
        var m = new IncrementalMesh(flat, Tet, persistent: true);
        var baseSnap = m.Current;
        var edges = CutVerify.Edges(m);
        var rad = 0.35 / g;
        var groups = new List<(int, int)>[g * g * g];
        for (var c = 0; c < groups.Length; c++) groups[c] = new List<(int, int)>();
        var edgeOrder = new List<(int, int)>();
        foreach (var (a, b) in edges)
        {
            int Ix(int v, int d) => Math.Min(g - 1, (int)(coords[3 * v + d] * g));
            var c = (Ix(a, 0) * g + Ix(a, 1)) * g + Ix(a, 2);
            int i = c / (g * g), j = c / g % g, k = c % g;
            double Phi(int v)
            {
                double dx = coords[3 * v] - (i + 0.5) / g, dy = coords[3 * v + 1] - (j + 0.5) / g,
                    dz = coords[3 * v + 2] - (k + 0.5) / g;
                return Math.Sqrt(dx * dx + dy * dy + dz * dz) - rad;
            }

            if (Phi(a) * Phi(b) >= 0) continue;
            groups[c].Add((a, b));
            edgeOrder.Add((a, b));
        }

        var groupOrder = groups.SelectMany(l => l).ToList();
        Console.WriteLine($"  n={n}: {rows.Length} tetrahedra, {edgeOrder.Count} marked edges");
        using var w = Csv($"cutgc_{n}_{g}.csv", "variant,ms_min,ms_median," + GcSamples.CsvHead + ",heap_mb");

        void Run(string name, List<(int, int)> marks)
        {
            var t = new GcSamples();
            for (var r = 0; r < reps; r++)
            {
                m.Rollback(baseSnap);
                FullGC();
                t.Add(TimeMs(() => Cutter.Cut(m, marks)));
            }

            var heap = GC.GetTotalMemory(false) / (1 << 20);
            Console.WriteLine($"  {name,-36} min {t.Min:F1} ms, median {t.Median:F1} ms, {t.Text}; heap {heap} MB");
            w.WriteLine($"{name},{t.Min:F4},{t.Median:F4},{t.CsvTail},{heap}");
        }

        Run("group order, one mesh", groupOrder);
        Run("edge order, one mesh", edgeOrder);
        var second = new SimplexMesh();
        second.WithBatch(() =>
        {
            for (var v = 0; v < nN; v++) second.AddNode(v);
            foreach (var r in rows) second.Add<Tet4, Node>(r);
        });
        Run("group order, second mesh alive", groupOrder);
        Run("edge order, second mesh alive", edgeOrder);
        GC.KeepAlive(second);
    }
}

// ParallelCutVerify.cs — the concurrent cutter against the sequential cutter and the oracle.
//
//   verify-cutpar n rounds seed threads
//
// Two copies of the mixed box mesh are cut in successive rounds with the same groups of
// marked edges: one by Cutter.CutParallel, one by Cutter.Cut of the concatenated groups.
// Groups are random edge subsets of cells of a grid over the box; the cells are
// separated by a gap of four initial element sizes (independent groups), by one
// (some groups merged through shared nodes), or touch while a sparse group spanning the
// whole box is added (everything merged). Each
// round requires: the same node count and node parents, the same multiset of rows
// together with the ancestor element recorded through ElementSplit, every derived
// structure equal to the oracle, conformity, and unchanged published versions (content
// hash of the previous version after the concurrent copy-on-write). Every third round
// reverts the cut of one group with Uncut on its own record and the same edges on the
// sequential record, and compares again; every fourth round rolls the parallel copy
// back, repeats the cut with another thread count and compares again.

using Step2;

namespace Step2.Bench;

public static class ParallelCutVerify
{
    private sealed class OriginObserver(Attribute<int> origin) : ICutObserver
    {
        public void NodeCreated(int node, int a, int b) { }
        public void ElementSplit(int parent, int child) => origin.Set(child, origin[parent]);
    }

    private static (IncrementalMesh M, Attribute<int> Origin, ICutObserver Obs) NewMesh(List<int[]> rows, int nN)
    {
        var flat = FlatBuild.Build(ElementSpec.MixedSimplex, nN, rows, parallel: false);
        var m = new IncrementalMesh(flat, ElementSpec.MixedSimplex, persistent: true);
        var origin = m.AddAttribute(AttributeKind.Element, -1);
        for (var e = 0; e < m.ElementCount; e++) origin.Set(e, e);
        m.Commit();
        return (m, origin, new OriginObserver(origin));
    }

    /// <summary>Live rows with the ancestor element appended (multiset comparison).</summary>
    private static List<int[]> Tagged(Snapshot s, Attribute<int> origin)
    {
        var o = s.View(origin);
        var list = new List<int[]>(s.ElementCount);
        for (var e = 0; e < s.ElementCount; e++)
        {
            var r = s.Row(e);
            if (r.Length == 0) continue;
            var t = new int[r.Length + 1];
            r.CopyTo(t);
            t[r.Length] = o[e];
            list.Add(t);
        }

        return list;
    }

    public static void Run(int n, int rounds, int seed, int threads)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var (rows, coords, nN) = CutVerify.MixedBox(n);
        var (mp, op, obsP) = NewMesh(rows, nN);
        var (ms, os, obsS) = NewMesh(rows, nN);
        var x = new List<double>(coords);
        var rnd = new Random(seed);
        long edgesCut = 0, groupsApplied = 0, merges = 0;
        var h = 1.0 / n;

        for (var r = 0; r < rounds; r++)
        {
            var before = mp.Current;
            var beforeHash = Oracle.Hash(before);
            var edges = CutVerify.Edges(mp);

            // ---- groups: cells of a g x g x g grid; margin > 0 separates them
            // a shell reaches one element diagonal (< 1.8 h) beyond its cell: a gap of 4 h
            // separates the groups, a gap of h lets some of them share nodes
            var g = 2 + rnd.Next(Math.Max(1, n / 8));
            var kind = r % 5 == 4 ? "spanning" : r % 3 == 1 ? "touching" : "separated";
            var margin = kind == "separated" ? 4 * h : kind == "touching" ? h : 0.0;
            var frac = new[] { 0.05, 0.2, 0.5 }[rnd.Next(3)];
            var groups = new List<List<(int, int)>>();
            for (var c = 0; c < g * g * g; c++) groups.Add(new List<(int, int)>());
            int Cell(int v)
            {
                var id = 0;
                for (var d = 0; d < 3; d++)
                {
                    var t = x[3 * v + d] * g;
                    var k = Math.Min(g - 1, (int)Math.Floor(t));
                    var off = (t - k) / g; // distance from the lower cell face
                    if (off < margin / 2 || 1.0 / g - off < margin / 2) return -1;
                    id = id * g + k;
                }

                return id;
            }

            // a pattern row holds at most 1023 neighbours (one chunk); repeated cuts of the same
            // region raise the valence of every node of the shells (each new node is adjacent to the
            // nodes of the shell of its edge), so an edge is marked only if no node of its shell has
            // more than 200 neighbours
            var shellScratch = new List<int>();
            bool Moderate(int a, int b)
            {
                shellScratch.Clear();
                Cutter.Shell(mp, a, b, shellScratch);
                foreach (var e in shellScratch)
                foreach (var v in mp.Row(e))
                    if (mp.PatternRow(v).Length > 200)
                        return false;
                return true;
            }

            foreach (var (a, b) in edges)
            {
                if (!Moderate(a, b)) continue;
                int ca = Cell(a), cb = Cell(b);
                if (ca < 0 || ca != cb || rnd.NextDouble() >= frac) continue;
                groups[ca].Add(rnd.Next(2) == 0 ? (a, b) : (b, a));
            }

            if (kind == "spanning")
                groups.Add(edges.Where(_ => rnd.NextDouble() < 0.01).ToList());
            foreach (var grp in groups)
                for (var i = grp.Count - 1; i > 0; i--)
                {
                    var j = rnd.Next(i + 1);
                    (grp[i], grp[j]) = (grp[j], grp[i]);
                }

            // ---- parallel and sequential cuts; every pattern row sees the same operations in the
            // same order, so the edit, slot-reuse and relocation counts must agree
            long eP = mp.Edits, rP = mp.SlotReuses, lP = mp.Relocations;
            long eS = ms.Edits, rS = ms.SlotReuses, lS = ms.Relocations;
            var res = Cutter.CutParallel(mp, groups.Select(l => (IReadOnlyList<(int, int)>)l).ToList(), obsP, threads);
            var concat = groups.SelectMany(l => l).ToList();
            var recS = Cutter.Cut(ms, concat, obsS);
            var sp = mp.Commit();
            var ss = ms.Commit();
            if (mp.Edits - eP != ms.Edits - eS || mp.SlotReuses - rP != ms.SlotReuses - rS ||
                mp.Relocations - lP != ms.Relocations - lS)
                throw new Exception($"round {r}: counters (edits, reuses, relocations) " +
                                    $"{mp.Edits - eP},{mp.SlotReuses - rP},{mp.Relocations - lP} vs " +
                                    $"{ms.Edits - eS},{ms.SlotReuses - rS},{ms.Relocations - lS}");
            foreach (var (a, b, _) in recS.Edges)
                for (var d = 0; d < 3; d++)
                    x.Add(0.5 * (x[3 * a + d] + x[3 * b + d]));

            if (Oracle.Hash(before) != beforeHash) throw new Exception($"round {r}: published version modified");
            Compare(sp, ss, op, os, res, recS, $"round {r}");
            edgesCut += recS.Edges.Length;
            groupsApplied += res.Groups.Length;
            merges += res.Merged;
            Console.WriteLine($"  round {r,2} ({kind,-9} g={g}, {frac,4:P0}): {res.Groups.Length,3} groups applied " +
                              $"({res.Merged,2} merged), {recS.Edges.Length,6} edges, {recS.ShellElements,7} shell " +
                              $"elements, {sp.ElementCount,8} elements: equal to sequential, oracle, conformity, " +
                              "immutability pass");

            // ---- rollback and repeat with another thread count: same rows
            if (r % 4 == 3)
            {
                mp.Rollback(before);
                var t2 = Math.Max(1, threads / 2 + rnd.Next(threads + 1));
                var res2 = Cutter.CutParallel(mp, groups.Select(l => (IReadOnlyList<(int, int)>)l).ToList(), obsP, t2);
                var sp2 = mp.Commit();
                Compare(sp2, ss, op, os, res2, recS, $"round {r} repeat ({t2} threads)");
                res = res2;
                Console.WriteLine($"           rollback and repeat with {t2} threads: equal");
            }

            // ---- revert one group on its own record, the same edges on the sequential record
            if (r % 3 == 2 && res.Groups.Length > 0)
            {
                var k = rnd.Next(res.Groups.Length);
                var keys = res.Groups[k].Edges.Select(t => Cutter.Key(t.A, t.B)).ToHashSet(EdgeKeyComparer.Instance);
                Cutter.Uncut(mp, res.Groups[k], keys, obsP);
                Cutter.Uncut(ms, recS, keys, obsS);
                var up = mp.Commit();
                var us = ms.Commit();
                Oracle.Verify(up);
                CutVerify.Conformity(up, $"round {r} revert");
                CutVerify.SameRows(Tagged(up, op), Tagged(us, os), $"round {r} revert of group {k}");
                Console.WriteLine($"           revert of group {k} ({keys.Count} edges): equal to the sequential revert");
            }
        }

        Console.WriteLine($"verify-cutpar n={n} threads={threads}: {rounds} rounds, {groupsApplied} groups applied, " +
                          $"{merges} merges, {edgesCut} cut edges: equal to the sequential cut, oracle, conformity " +
                          $"and immutability pass ({sw.Elapsed.TotalSeconds:F1} s)");
    }

    private static void Compare(Snapshot sp, Snapshot ss, Attribute<int> op, Attribute<int> os, ParallelCutResult res,
        CutRecord recS, string what)
    {
        if (sp.NodeCount != ss.NodeCount) throw new Exception($"{what}: {sp.NodeCount} nodes vs {ss.NodeCount}");
        var parP = res.Groups.SelectMany(g => g.Edges).OrderBy(t => t.Mid).ToArray();
        if (!parP.SequenceEqual(recS.Edges)) throw new Exception($"{what}: new nodes or parents differ");
        CutVerify.SameRows(Tagged(sp, op), Tagged(ss, os), what);
        Oracle.Verify(sp);
        CutVerify.Conformity(sp, what);
    }
}

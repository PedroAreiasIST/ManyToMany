// CutVerify.cs — the local cutter against the batch cutter and the oracle.
//
//   verify-cut n rounds seed
//
// A mixed mesh (Kuhn tetrahedra, the boundary triangles, bars on the box edges and
// point elements on boundary nodes) is cut in successive rounds. Each round marks a
// set of edges (a random fraction of all edges, or the edges crossing a sphere),
// cuts it locally on the incremental mesh and with the batch cutter
// (Numerical.MeshRefinement.Refine, the C# translation of remeshsimplex.f90) applied
// to the previous batch output, and requires: the same node count and node parents,
// the same multiset of rows for every element type, every derived structure equal
// to the oracle, and conformity (every triangular face of a tetrahedron shared by at
// most two tetrahedra, and the boundary faces of the tetrahedra equal, as node sets,
// to the triangle elements).

using Numerical;
using Step2;

namespace Step2.Bench;

public static class CutVerify
{
    public static (List<int[]> Rows, List<double> Coords, int NodeCount) MixedBox(int n, int pointStride = 5)
    {
        var (tets, coords, nN) = Kuhn.Box(n);
        var rows = new List<int[]>(tets);
        var flat = FlatBuild.Build(ElementSpec.Tetrahedron, nN, tets, parallel: false);
        var boundaryNodes = new SortedSet<int>();
        foreach (var id in flat.Boundary)
        {
            int e = id / 4, f = id % 4;
            var tri = ElementSpec.Tetrahedron.FacetLocal[f].Select(l => tets[e][l]).ToArray();
            rows.Add(tri);
            foreach (var v in tri) boundaryNodes.Add(v);
        }

        // bars along the twelve edges of the box
        for (var axis = 0; axis < 3; axis++)
        for (var a = 0; a <= n; a += n)
        for (var b = 0; b <= n; b += n)
        for (var t = 0; t < n; t++)
        {
            int[] p = axis switch { 0 => [t, a, b], 1 => [a, t, b], _ => [a, b, t] };
            int[] q = axis switch { 0 => [t + 1, a, b], 1 => [a, t + 1, b], _ => [a, b, t + 1] };
            rows.Add([Kuhn.Node(n, p[0], p[1], p[2]), Kuhn.Node(n, q[0], q[1], q[2])]);
        }

        // point elements (boundary conditions) on every pointStride-th boundary node
        var i = 0;
        foreach (var v in boundaryNodes)
            if (i++ % pointStride == 0)
                rows.Add([v]);
        return (rows, coords, nN);
    }

    private static SimplexMesh ToBatch(IReadOnlyList<int[]> rows, int nN)
    {
        var sm = new SimplexMesh();
        sm.WithBatch(() =>
        {
            for (var v = 0; v < nN; v++) sm.AddNode(v);
            foreach (var r in rows)
                switch (r.Length)
                {
                    case 1: sm.Add<Point, Node>(r); break;
                    case 2: sm.Add<Bar2, Node>(r); break;
                    case 3: sm.Add<Tri3, Node>(r); break;
                    case 4: sm.Add<Tet4, Node>(r); break;
                }
        });
        return sm;
    }

    private static List<int[]> BatchRows(SimplexMesh sm)
    {
        var rows = new List<int[]>();
        for (var e = 0; e < sm.Count<Point>(); e++) rows.Add(sm.NodesOf<Point, Node>(e).ToArray());
        for (var e = 0; e < sm.Count<Bar2>(); e++) rows.Add(sm.NodesOf<Bar2, Node>(e).ToArray());
        for (var e = 0; e < sm.Count<Tri3>(); e++) rows.Add(sm.NodesOf<Tri3, Node>(e).ToArray());
        for (var e = 0; e < sm.Count<Tet4>(); e++) rows.Add(sm.NodesOf<Tet4, Node>(e).ToArray());
        return rows;
    }

    private static string RowKey(ReadOnlySpan<int> r) => string.Join(',', r.ToArray());

    /// <summary>Multiset equality of live rows (row order and orientation included).</summary>
    public static void SameRows(IEnumerable<int[]> a, IEnumerable<int[]> b, string what)
    {
        var count = new Dictionary<string, int>();
        var na = 0;
        foreach (var r in a)
        {
            if (r.Length == 0) continue;
            var k = RowKey(r);
            count[k] = count.GetValueOrDefault(k) + 1;
            na++;
        }

        var nb = 0;
        foreach (var r in b)
        {
            if (r.Length == 0) continue;
            nb++;
            var k = RowKey(r);
            if (!count.TryGetValue(k, out var c) || c == 0)
            {
                if (Environment.GetEnvironmentVariable("STEP2_DEBUG") == "1")
                {
                    var ca = new Dictionary<string, int>();
                    foreach (var x in a) if (x.Length > 0) ca[RowKey(x)] = ca.GetValueOrDefault(RowKey(x)) + 1;
                    var cb = new Dictionary<string, int>();
                    foreach (var x in b) if (x.Length > 0) cb[RowKey(x)] = cb.GetValueOrDefault(RowKey(x)) + 1;
                    foreach (var key in ca.Keys.Union(cb.Keys))
                        if (ca.GetValueOrDefault(key) != cb.GetValueOrDefault(key))
                            Console.WriteLine($"    ({key}): first {ca.GetValueOrDefault(key)} second {cb.GetValueOrDefault(key)}");
                }

                throw new Exception($"{what}: row ({k}) of the second mesh missing from the first");
            }
            count[k] = c - 1;
        }

        if (na != nb) throw new Exception($"{what}: {na} live rows vs {nb}");
    }

    /// <summary>Conformity of tetrahedra and compatibility of the triangle elements with their faces.</summary>
    public static void Conformity(Snapshot s, string what)
    {
        var faces = new Dictionary<(int, int, int), int>();
        var tris = new Dictionary<(int, int, int), int>();
        for (var e = 0; e < s.ElementCount; e++)
        {
            var r = s.Row(e).ToArray();
            if (r.Length == 4)
                for (var f = 0; f < 4; f++)
                {
                    var key = Sorted(ElementSpec.Tetrahedron.FacetLocal[f].Select(l => r[l]).ToArray());
                    faces[key] = faces.GetValueOrDefault(key) + 1;
                }
            else if (r.Length == 3)
            {
                var key = Sorted(r);
                tris[key] = tris.GetValueOrDefault(key) + 1;
            }
        }

        var boundary = 0;
        foreach (var (key, c) in faces)
        {
            if (c > 2) throw new Exception($"{what}: face {key} shared by {c} tetrahedra");
            if (c == 1)
            {
                boundary++;
                if (tris.GetValueOrDefault(key) != 1)
                    throw new Exception($"{what}: boundary face {key} not covered by exactly one triangle");
            }
        }

        foreach (var (key, c) in tris)
            if (faces.GetValueOrDefault(key) != 1 || c != 1)
                throw new Exception($"{what}: triangle {key} is not exactly one boundary face");
        if (boundary != tris.Count) throw new Exception($"{what}: {boundary} boundary faces vs {tris.Count} triangles");
    }

    private static (int, int, int) Sorted(int[] t)
    {
        Array.Sort(t);
        return (t[0], t[1], t[2]);
    }

    /// <summary>Distinct edges (a &lt; b) of the live elements of length >= 2.</summary>
    public static List<(int, int)> Edges(IncrementalMesh m)
    {
        var list = new List<(int, int)>();
        for (var v = 0; v < m.NodeCount; v++)
        {
            var cols = m.PatternColumns(v);
            var cnts = m.PatternCounts(v);
            for (var t = 0; t < cols.Length; t++)
                if (cnts[t] > 0 && cols[t] > v)
                    list.Add((v, cols[t]));
        }

        return list;
    }

    public static void Run(int n, int rounds, int seed)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var (rows, coords, nN) = MixedBox(n);
        var flat = FlatBuild.Build(ElementSpec.MixedSimplex, nN, rows, parallel: false);
        var m = new IncrementalMesh(flat, ElementSpec.MixedSimplex, persistent: true);
        Oracle.Verify(m.Current);
        Conformity(m.Current, "initial mesh");
        var batch = ToBatch(rows, nN);
        var x = new List<double>(coords);
        var rnd = new Random(seed);
        var caseCount = new long[7];
        long totalMarked = 0, totalShell = 0;
        var quiet = TextWriter.Null;
        var stdout = Console.Out;

        for (var r = 0; r < rounds; r++)
        {
            // ---- marked edges: a random fraction, or the edges crossing a sphere
            var edges = Edges(m);
            List<(int, int)> marked;
            string how;
            if (r % 3 == 2)
            {
                double cx = rnd.NextDouble(), cy = rnd.NextDouble(), cz = rnd.NextDouble();
                var rad = 0.15 + 0.35 * rnd.NextDouble();
                double Phi(int v) =>
                    Math.Sqrt(Sq(x[3 * v] - cx) + Sq(x[3 * v + 1] - cy) + Sq(x[3 * v + 2] - cz)) - rad;
                marked = edges.Where(e => Phi(e.Item1) * Phi(e.Item2) < 0).ToList();
                how = $"sphere r={rad:F2}";
            }
            else
            {
                // dense fractions are restricted to a random sub-box so that the mesh stays small
                var frac = new[] { 0.005, 0.05, 0.3, 0.7, 0.97 }[rnd.Next(5)];
                double bx = rnd.NextDouble() * 0.7, by = rnd.NextDouble() * 0.7, bz = rnd.NextDouble() * 0.7;
                var w = frac > 0.5 ? 0.3 : 1.5;
                bool In(int v) => x[3 * v] >= bx && x[3 * v] <= bx + w && x[3 * v + 1] >= by &&
                                  x[3 * v + 1] <= by + w && x[3 * v + 2] >= bz && x[3 * v + 2] <= bz + w;
                marked = edges.Where(e => In(e.Item1) && In(e.Item2) && rnd.NextDouble() < frac).ToList();
                how = $"random {frac:P1}";
            }

            // random order and random orientation of each pair
            for (var i = marked.Count - 1; i > 0; i--)
            {
                var j = rnd.Next(i + 1);
                (marked[i], marked[j]) = (marked[j], marked[i]);
            }

            for (var i = 0; i < marked.Count; i++)
                if (rnd.Next(2) == 0)
                    marked[i] = (marked[i].Item2, marked[i].Item1);

            // ---- tetrahedron cases of this round (marked edges per tetrahedron)
            var markSet = new HashSet<long>(marked.Select(e => Cutter.Key(e.Item1, e.Item2)), EdgeKeyComparer.Instance);
            for (var e = 0; e < m.ElementCount; e++)
            {
                var row = m.Row(e);
                if (row.Length != 4) continue;
                var c = 0;
                foreach (var en in SubdivisionTemplates.EdgeNodes)
                    if (markSet.Contains(Cutter.Key(row[en[0]], row[en[1]])))
                        c++;
                caseCount[c]++;
            }

            // ---- local cut
            var rec = Cutter.Cut(m, marked);
            foreach (var (_, a, b) in rec.Edges.Select(t => (t.Mid, t.A, t.B)))
                for (var d = 0; d < 3; d++)
                    x.Add(0.5 * (x[3 * a + d] + x[3 * b + d]));
            var snap = m.Commit();
            totalMarked += rec.Edges.Length;
            totalShell += rec.ShellElements;

            // ---- batch cut of the previous batch output
            Console.SetOut(quiet);
            var (next, _) = MeshRefinement.Refine(batch, marked);
            Console.SetOut(stdout);
            batch = next;

            // ---- comparisons
            if (batch.Count<Node>() != m.NodeCount)
                throw new Exception($"round {r}: {batch.Count<Node>()} batch nodes vs {m.NodeCount}");
            foreach (var (a, b, c) in rec.Edges)
            {
                var p = batch.Get<Node, ParentNodes>(c);
                if (Math.Min(p.Parent1, p.Parent2) != a || Math.Max(p.Parent1, p.Parent2) != b)
                    throw new Exception($"round {r}: node {c} parents ({p.Parent1},{p.Parent2}) vs ({a},{b})");
            }

            var mine = new List<int[]>();
            for (var e = 0; e < snap.ElementCount; e++) mine.Add(snap.Row(e).ToArray());
            SameRows(mine, BatchRows(batch), $"round {r}");
            Oracle.Verify(snap);
            Conformity(snap, $"round {r}");
            Console.WriteLine($"  round {r,2} ({how,-14}): {rec.Edges.Length,6} edges, {rec.ShellElements,7} shell " +
                              $"elements, {snap.ElementCount,8} elements, {snap.NodeCount,7} nodes: equal to batch, " +
                              "oracle and conformity pass");
        }

        Console.WriteLine($"verify-cut n={n}: {rounds} rounds, {totalMarked} cut edges, {totalShell} split elements; " +
                          $"tetrahedra by marked edges 0..6: {string.Join(' ', caseCount)} ({sw.Elapsed.TotalSeconds:F1} s)");
    }

    private static double Sq(double v) => v * v;
}

public static class UncutVerify
{
    /// <summary>
    ///     verify-uncut n trials seed: cut a marked set M, revert a random subset S locally, and
    ///     compare with the batch cut of the uncut mesh with M \ S (labels mapped by rank); then
    ///     revert the rest and compare with the mesh before the cut. Two-level variant: a second
    ///     cut on top is reverted first.
    /// </summary>
    public static void Run(int n, int trials, int seed)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var (rows, coords, nN) = CutVerify.MixedBox(n);
        var flat = FlatBuild.Build(ElementSpec.MixedSimplex, nN, rows, parallel: false);
        var m = new IncrementalMesh(flat, ElementSpec.MixedSimplex, persistent: true);
        var rnd = new Random(seed);
        var quiet = TextWriter.Null;
        var stdout = Console.Out;
        long reverted = 0, kept = 0;

        for (var t = 0; t < trials; t++)
        {
            var before = m.Commit();
            var beforeRows = Rows(before);
            var first = before.NodeCount;
            var edges = CutVerify.Edges(m);
            var frac = new[] { 0.02, 0.1, 0.4 }[rnd.Next(3)];
            var M = edges.Where(_ => rnd.NextDouble() < frac).Select(e => rnd.Next(2) == 0 ? e : (e.Item2, e.Item1))
                .ToList();
            var rec = Cutter.Cut(m, M);

            // optional second level, reverted completely before the first is touched
            if (t % 2 == 1)
            {
                var e2 = CutVerify.Edges(m).Where(_ => rnd.NextDouble() < 0.05).ToList();
                var rec2 = Cutter.Cut(m, e2);
                Oracle.Verify(m.Commit());
                Cutter.Uncut(m, rec2, rec2.Edges.Select(x => Cutter.Key(x.A, x.B)).ToHashSet(EdgeKeyComparer.Instance));
                var s2 = m.Commit();
                Oracle.Verify(s2);
                var lvl = Rows(s2).Where(r => r.Length > 0).ToList();
                CutVerify.SameRows(lvl, AfterCut(beforeRows, M, first, quiet, stdout, nN0: first), $"trial {t} level 2 revert");
            }

            // revert a random subset S (sometimes all of M)
            var sFrac = t % 5 == 4 ? 1.0 : rnd.NextDouble();
            var S = rec.Edges.Where(_ => rnd.NextDouble() < sFrac).Select(x => Cutter.Key(x.A, x.B)).ToHashSet(EdgeKeyComparer.Instance);
            var rest = Cutter.Uncut(m, rec, S);
            var snap = m.Commit();
            Oracle.Verify(snap);
            CutVerify.Conformity(snap, $"trial {t}");

            // batch reference: cut of the mesh before the cut with M \ S, in M order
            var keepList = M.Where(e => !S.Contains(Cutter.Key(e.Item1, e.Item2))).ToList();
            var refRows = AfterCut(beforeRows, keepList, first, quiet, stdout, nN0: first);
            // our labels -> batch labels: old nodes unchanged, surviving new nodes by rank
            var map = new Dictionary<int, int>();
            var rank = 0;
            foreach (var (_, _, c) in rest.Edges) map[c] = first + rank++;
            var mapped = Rows(snap).Where(r => r.Length > 0)
                .Select(r => r.Select(v => v < first ? v : map.TryGetValue(v, out var w) ? w
                    : throw new Exception($"trial {t}: reverted node {v} still used")).ToArray()).ToList();
            CutVerify.SameRows(mapped, refRows, $"trial {t} partial revert");
            reverted += S.Count;
            kept += rest.Edges.Length;

            // revert the rest: the mesh before the cut, as a multiset of rows
            Cutter.Uncut(m, rest, rest.Edges.Select(x => Cutter.Key(x.A, x.B)).ToHashSet(EdgeKeyComparer.Instance));
            var back = m.Commit();
            Oracle.Verify(back);
            CutVerify.SameRows(Rows(back), beforeRows, $"trial {t} full revert");

            // continue from a cut mesh so that later trials revert inside refined regions
            Cutter.Cut(m, M);
        }

        Console.WriteLine($"verify-uncut n={n}: {trials} trials, {reverted} edges reverted, {kept} kept: partial " +
                          "reverts equal the batch cut with the remaining marks, full reverts restore the mesh, " +
                          $"oracle and conformity pass ({sw.Elapsed.TotalSeconds:F1} s)");
    }

    private static List<int[]> Rows(Snapshot s)
    {
        var r = new List<int[]>(s.ElementCount);
        for (var e = 0; e < s.ElementCount; e++) r.Add(s.Row(e).ToArray());
        return r;
    }

    /// <summary>Batch cut of a row set (nodes 0..nN0-1) with the given marks.</summary>
    private static List<int[]> AfterCut(List<int[]> rows, List<(int, int)> marks, int first, TextWriter quiet,
        TextWriter stdout, int nN0)
    {
        var sm = new SimplexMesh();
        sm.WithBatch(() =>
        {
            for (var v = 0; v < nN0; v++) sm.AddNode(v);
            foreach (var r in rows)
                switch (r.Length)
                {
                    case 1: sm.Add<Point, Node>(r); break;
                    case 2: sm.Add<Bar2, Node>(r); break;
                    case 3: sm.Add<Tri3, Node>(r); break;
                    case 4: sm.Add<Tet4, Node>(r); break;
                }
        });
        Console.SetOut(quiet);
        var (o, _) = MeshRefinement.Refine(sm, marks);
        Console.SetOut(stdout);
        var res = new List<int[]>();
        for (var e = 0; e < o.Count<Point>(); e++) res.Add(o.NodesOf<Point, Node>(e).ToArray());
        for (var e = 0; e < o.Count<Bar2>(); e++) res.Add(o.NodesOf<Bar2, Node>(e).ToArray());
        for (var e = 0; e < o.Count<Tri3>(); e++) res.Add(o.NodesOf<Tri3, Node>(e).ToArray());
        for (var e = 0; e < o.Count<Tet4>(); e++) res.Add(o.NodesOf<Tet4, Node>(e).ToArray());
        return res;
    }
}

// CrackVerify.cs — local crack insertion and growth against a global evaluation.
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
//   verify-crack n steps case
//
// The global evaluation applies the same definition (CrackInserter) to the whole mesh:
// crack edges by a scan of all edges, the batch cutter (MeshRefinement.Refine) for the
// cut, and the split of every surface node from a node-to-element map rebuilt from the
// rows. After every increment the incremental mesh must have the same rows (as a
// multiset), the same node count and coordinates, and the walk must find the same crack
// edges as a scan; the oracle and a conformity check that admits the crack faces (pairs
// of coincident boundary faces whose nodes lie on the surface) complete the test.

using Numerical;
using Step2;

namespace Step2.Bench;

public static class CrackVerify
{
    public static (CrackInserter.Field S, CrackInserter.Field R, string Name) Case(int c, int k, int steps)
    {
        var f = (k + 1.0) / steps;
        switch (c)
        {
            case 0: // inclined penny crack growing from the centre
            {
                CrackInserter.Field s = (x, y, z) => z - 0.5 - 0.17 * (x - 0.5) - 0.05 * (y - 0.5);
                var R = 0.08 + 0.3 * f;
                CrackInserter.Field r = (x, y, z) =>
                    Math.Sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) - R;
                return (s, r, $"inclined penny R={R:F3}");
            }
            case 1: // edge crack entering from the face x = 0
            {
                CrackInserter.Field s = (x, y, z) => z - 0.47 - 0.11 * y;
                var a = 0.12 + 0.6 * f;
                CrackInserter.Field r = (x, y, z) => x - a;
                return (s, r, $"edge crack a={a:F3}");
            }
            default: // curved (spherical cap) crack growing in angle
            {
                CrackInserter.Field s = (x, y, z) =>
                    Math.Sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z + 0.3) * (z + 0.3)) - 0.93;
                var R = 0.1 + 0.28 * f;
                CrackInserter.Field r = (x, y, z) =>
                    Math.Sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) - R;
                return (s, r, $"spherical cap R={R:F3}");
            }
        }
    }

    public static void Run(int n, int steps, int c)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var (rows, coords, nN) = CutVerify.MixedBox(n);
        var flat = FlatBuild.Build(ElementSpec.MixedSimplex, nN, rows, parallel: false);
        var m = new IncrementalMesh(flat, ElementSpec.MixedSimplex, persistent: true);
        var ins = new CrackInserter(m, coords.ToArray());
        var g = new GlobalCrack(rows, coords, nN);
        // a nodal field that is linear on every element of the initial mesh (random nodal values)
        var u = m.AddAttribute(AttributeKind.Node, 0.0);
        var rnd = new Random(17);
        var u0 = new double[nN];
        for (var v = 0; v < nN; v++) u.Set(v, u0[v] = rnd.NextDouble());
        ins.Follow(u);
        var rows0 = rows.Select(r => r.ToArray()).ToArray();
        double worst = 0;

        for (var k = 0; k < steps; k++)
        {
            var (phiS, phiR, name) = Case(c, k, steps);

            // the walk must find exactly the crack edges of a global scan
            if (ins.Front.Length > 0)
            {
                var scan = ins.ScanAll(phiS, phiR);
                var walk = ins.Walk(ins.Front, phiS, phiR, new Dictionary<long, double>(EdgeKeyComparer.Instance), out _);
                if (!scan.SequenceEqual(walk))
                    throw new Exception($"step {k}: walk found {walk.Count} crack edges, scan {scan.Count}");
            }

            var st = ins.Advance(phiS, phiR);
            var snap = m.Commit();
            g.Advance(phiS, phiR);

            if (Environment.GetEnvironmentVariable("STEP2_DEBUG") == "1")
            {
                Console.WriteLine($"  local splits:  {string.Join(' ', st.Splits.Select(p => $"{p.Node}->{p.Copy}"))}");
                Console.WriteLine($"  global splits: {string.Join(' ', g.SplitLog.Select(p => $"{p.Item1}->{p.Item2}"))}");
                Console.WriteLine($"  cut edges local {st.Edges.Length}, global marked {g.LastMarked}");
                var lc = new int[5]; var gc = new int[5];
                for (var e = 0; e < snap.ElementCount; e++) lc[snap.Row(e).Length]++;
                foreach (var r in g.Rows) gc[r.Length]++;
                Console.WriteLine($"  rows by length local {string.Join(' ', lc)} global {string.Join(' ', gc)}");
                g.SplitLog.Clear();
            }

            if (g.NodeCount != snap.NodeCount)
                throw new Exception($"step {k}: {snap.NodeCount} nodes vs {g.NodeCount} (global)");
            var mine = new List<int[]>();
            for (var e = 0; e < snap.ElementCount; e++) mine.Add(snap.Row(e).ToArray());
            CutVerify.SameRows(mine, g.Rows, $"step {k}");
            var xv = snap.View(ins.X);
            var yv = snap.View(ins.Y);
            var zv = snap.View(ins.Z);
            for (var v = 0; v < snap.NodeCount; v++)
                if (xv[v] != g.X[3 * v] || yv[v] != g.X[3 * v + 1] || zv[v] != g.X[3 * v + 2])
                    throw new Exception($"step {k}: coordinates of node {v} differ");
            Oracle.Verify(snap);
            CrackConformity(snap, ins, $"step {k}");
            worst = Math.Max(worst, TransferError(snap, ins, u, rows0, u0, coords));
            Console.WriteLine($"  step {k} ({name}): {st.Edges.Length,5} crack edges, {st.Cut.ShellElements,5} split " +
                              $"elements, {st.Splits.Length,5} split nodes, front {st.Front.Length,4}, walk visited " +
                              $"{st.VisitedTets,6} tets; equal to the global evaluation, oracle and conformity pass");
        }

        if (worst > 1e-12) throw new Exception($"transfer error {worst:E2}");
        Console.WriteLine($"verify-crack n={n} case {c}: {steps} increments verified, field transfer equal to the " +
                          $"interpolant of the initial field (max error {worst:E1}) ({sw.Elapsed.TotalSeconds:F1} s)");
    }

    /// <summary>
    ///     Largest difference between the transferred field and the P1 interpolant of the initial
    ///     field, evaluated at every node of every tetrahedron through barycentric coordinates in
    ///     its ancestor tetrahedron.
    /// </summary>
    private static double TransferError(Snapshot s, CrackInserter ins, Attribute<double> u, int[][] rows0,
        double[] u0, List<double> x0)
    {
        var uv = s.View(u);
        var origin = s.View(ins.Origin);
        var xv = s.View(ins.X);
        var yv = s.View(ins.Y);
        var zv = s.View(ins.Z);
        double worst = 0;
        Span<double> lam = stackalloc double[4];
        for (var e = 0; e < s.ElementCount; e++)
        {
            var r = s.Row(e);
            if (r.Length != 4) continue;
            var anc = rows0[origin[e]];
            foreach (var v in r)
            {
                Barycentric(x0, anc, xv[v], yv[v], zv[v], lam);
                var ref0 = 0.0;
                for (var k = 0; k < 4; k++) ref0 += lam[k] * u0[anc[k]];
                worst = Math.Max(worst, Math.Abs(ref0 - uv[v]));
                if (lam[0] < -1e-12 || lam[1] < -1e-12 || lam[2] < -1e-12 || lam[3] < -1e-12)
                    throw new Exception($"node {v} of element {e} lies outside its ancestor");
            }
        }

        return worst;
    }

    private static void Barycentric(List<double> x, int[] t, double px, double py, double pz, Span<double> lam)
    {
        double X(int v, int d) => x[3 * v + d];
        double ax = X(t[0], 0), ay = X(t[0], 1), az = X(t[0], 2);
        double m00 = X(t[1], 0) - ax, m01 = X(t[2], 0) - ax, m02 = X(t[3], 0) - ax;
        double m10 = X(t[1], 1) - ay, m11 = X(t[2], 1) - ay, m12 = X(t[3], 1) - ay;
        double m20 = X(t[1], 2) - az, m21 = X(t[2], 2) - az, m22 = X(t[3], 2) - az;
        double bx = px - ax, by = py - ay, bz = pz - az;
        var det = m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20);
        var l1 = (bx * (m11 * m22 - m12 * m21) - m01 * (by * m22 - m12 * bz) + m02 * (by * m21 - m11 * bz)) / det;
        var l2 = (m00 * (by * m22 - m12 * bz) - bx * (m10 * m22 - m12 * m20) + m02 * (m10 * bz - by * m20)) / det;
        var l3 = (m00 * (m11 * bz - by * m21) - m01 * (m10 * bz - by * m20) + bx * (m10 * m21 - m11 * m20)) / det;
        lam[0] = 1 - l1 - l2 - l3;
        lam[1] = l1;
        lam[2] = l2;
        lam[3] = l3;
    }

    /// <summary>Conformity with crack faces: unpaired tetrahedron faces are triangles or crack faces.</summary>
    public static void CrackConformity(Snapshot s, CrackInserter ins, string what)
    {
        var faces = new Dictionary<(int, int, int), int>();
        var tris = new HashSet<(int, int, int)>();
        var on = s.View(ins.OnSurface);
        var xv = s.View(ins.X);
        var yv = s.View(ins.Y);
        var zv = s.View(ins.Z);
        var crackFaces = new Dictionary<string, int>();
        for (var e = 0; e < s.ElementCount; e++)
        {
            var r = s.Row(e).ToArray();
            if (r.Length == 3)
            {
                var t = (int[])r.Clone();
                Array.Sort(t);
                tris.Add((t[0], t[1], t[2]));
            }

            if (r.Length != 4) continue;
            for (var f = 0; f < 4; f++)
            {
                var t = ElementSpec.Tetrahedron.FacetLocal[f].Select(l => r[l]).ToArray();
                Array.Sort(t);
                var key = (t[0], t[1], t[2]);
                faces[key] = faces.GetValueOrDefault(key) + 1;
            }
        }

        foreach (var (key, cnt) in faces)
        {
            if (cnt > 2) throw new Exception($"{what}: face {key} shared by {cnt} tetrahedra");
            if (cnt == 2 || tris.Contains(key)) continue;
            if (on[key.Item1] == 0 || on[key.Item2] == 0 || on[key.Item3] == 0)
                throw new Exception($"{what}: unpaired face {key} is neither a triangle nor a crack face");
            // crack faces come in pairs of coincident faces
            var pts = new[] { key.Item1, key.Item2, key.Item3 }
                .Select(v => (xv[v], yv[v], zv[v])).OrderBy(p => p).Select(p => $"{p}");
            var geo = string.Join("|", pts);
            crackFaces[geo] = crackFaces.GetValueOrDefault(geo) + 1;
        }

        foreach (var (geo, cnt) in crackFaces)
            if (cnt != 2) throw new Exception($"{what}: crack face at {geo} appears {cnt} times");
    }
}

/// <summary>The crack definition evaluated on the whole mesh (reference).</summary>
public sealed class GlobalCrack
{
    public List<int[]> Rows;
    public List<double> X;
    public List<byte> On;
    public List<int> P1, P2;
    public int NodeCount => On.Count;
    public double Tau = 0.05;
    public List<(int, int)> SplitLog = new();
    public int LastMarked;

    public GlobalCrack(List<int[]> rows, List<double> coords, int nN)
    {
        Rows = rows.Select(r => (int[])r.Clone()).ToList();
        X = new List<double>(coords);
        On = Enumerable.Repeat((byte)0, nN).ToList();
        P1 = Enumerable.Range(0, nN).ToList();
        P2 = Enumerable.Range(0, nN).ToList();
    }

    private int Sign(int v, CrackInserter.Field s) => On[v] != 0 ? 0 : s(X[3 * v], X[3 * v + 1], X[3 * v + 2]) >= 0 ? 1 : -1;

    public void Advance(CrackInserter.Field phiS, CrackInserter.Field phiR)
    {
        // ---- crack edges by a scan of all edges of all elements
        var edgeSet = new SortedSet<(int, int)>();
        foreach (var r in Rows)
            for (var p = 0; p < r.Length; p++)
            for (var q = p + 1; q < r.Length; q++)
                edgeSet.Add((Math.Min(r[p], r[q]), Math.Max(r[p], r[q])));
        var crossing = new List<(int, int)>();
        var rootOf = new Dictionary<long, double>(EdgeKeyComparer.Instance);
        foreach (var (a, b) in edgeSet)
        {
            if (Sign(a, phiS) * Sign(b, phiS) != -1) continue;
            var t = CrackInserter.RootOf(X[3 * a], X[3 * a + 1], X[3 * a + 2], X[3 * b], X[3 * b + 1], X[3 * b + 2], phiS);
            double x = X[3 * a] + t * (X[3 * b] - X[3 * a]), y = X[3 * a + 1] + t * (X[3 * b + 1] - X[3 * a + 1]),
                z = X[3 * a + 2] + t * (X[3 * b + 2] - X[3 * a + 2]);
            if (phiR(x, y, z) >= 0) continue;
            crossing.Add((a, b));
            rootOf[Cutter.Key(a, b)] = t;
        }

        // promotable: no tetrahedron of the star of v would have four surface nodes (node-to-element map
        // rebuilt from the rows)
        var tetsOf = new List<int>[NodeCount];
        for (var v = 0; v < NodeCount; v++) tetsOf[v] = new List<int>();
        for (var e = 0; e < Rows.Count; e++)
            if (Rows[e].Length == 4)
                foreach (var v in Rows[e])
                    tetsOf[v].Add(e);
        bool Promotable(int v, IReadOnlySet<int> done)
        {
            foreach (var e in tetsOf[v])
            {
                var on = 0;
                foreach (var w in Rows[e])
                    if (w != v && (On[w] != 0 || done.Contains(w)))
                        on++;
                if (on == 3) return false;
            }

            return true;
        }

        var (marked, promoted) = CrackInserter.Promote(crossing, rootOf, Tau, Promotable);
        foreach (var v in promoted) On[v] = 1;
        var roots = marked.Select(e => rootOf[Cutter.Key(e.Item1, e.Item2)]).ToList();
        LastMarked = marked.Count;

        // ---- batch cut
        var nN = NodeCount;
        var sm = new SimplexMesh();
        sm.WithBatch(() =>
        {
            for (var v = 0; v < nN; v++) sm.AddNode(v);
            foreach (var r in Rows)
                switch (r.Length)
                {
                    case 1: sm.Add<Point, Node>(r); break;
                    case 2: sm.Add<Bar2, Node>(r); break;
                    case 3: sm.Add<Tri3, Node>(r); break;
                    case 4: sm.Add<Tet4, Node>(r); break;
                }
        });
        var stdout = Console.Out;
        Console.SetOut(TextWriter.Null);
        var (o, _) = MeshRefinement.Refine(sm, marked);
        Console.SetOut(stdout);
        Rows = new List<int[]>();
        for (var e = 0; e < o.Count<Point>(); e++) Rows.Add(o.NodesOf<Point, Node>(e).ToArray());
        for (var e = 0; e < o.Count<Bar2>(); e++) Rows.Add(o.NodesOf<Bar2, Node>(e).ToArray());
        for (var e = 0; e < o.Count<Tri3>(); e++) Rows.Add(o.NodesOf<Tri3, Node>(e).ToArray());
        for (var e = 0; e < o.Count<Tet4>(); e++) Rows.Add(o.NodesOf<Tet4, Node>(e).ToArray());
        for (var i = 0; i < marked.Count; i++)
        {
            var (a, b) = marked[i];
            var t = roots[i];
            for (var d = 0; d < 3; d++) X.Add(X[3 * a + d] + t * (X[3 * b + d] - X[3 * a + d]));
            On.Add(1);
            P1.Add(a);
            P2.Add(b);
        }

        if (o.Count<Node>() != NodeCount) throw new Exception("global crack: node count mismatch after the cut");

        // ---- split every surface node whose star falls apart across crack faces
        var star = new Dictionary<int, List<int>>();
        for (var e = 0; e < Rows.Count; e++)
            foreach (var v in Rows[e])
            {
                if (!star.TryGetValue(v, out var l)) star[v] = l = new List<int>();
                l.Add(e);
            }

        var surface = Enumerable.Range(0, NodeCount).Where(v => On[v] != 0).ToList();
        foreach (var v in surface)
        {
            if (!star.TryGetValue(v, out var all)) continue;
            var tets = all.Where(e => Rows[e].Length == 4).ToList();
            if (tets.Count == 0) continue;
            // faces containing v -> tetrahedra
            var faceOwners = new Dictionary<(int, int), List<int>>();
            foreach (var e in tets)
            {
                var others = Rows[e].Where(w => w != v).ToArray();
                for (var i = 0; i < 3; i++)
                for (var j = i + 1; j < 3; j++)
                {
                    var key = (Math.Min(others[i], others[j]), Math.Max(others[i], others[j]));
                    if (!faceOwners.TryGetValue(key, out var l)) faceOwners[key] = l = new List<int>();
                    l.Add(e);
                }
            }

            var parent = tets.ToDictionary(e => e, e => e);
            int Find(int e)
            {
                while (parent[e] != e) e = parent[e] = parent[parent[e]];
                return e;
            }

            foreach (var ((a, b), owners) in faceOwners)
            {
                if (owners.Count != 2) continue;
                if (On[a] != 0 && On[b] != 0) continue; // crack face (v has sign 0)
                int ra = Find(owners[0]), rb = Find(owners[1]);
                if (ra != rb) parent[Math.Max(ra, rb)] = Math.Min(ra, rb);
            }

            var comps = tets.Select(Find).Distinct().ToList();
            if (comps.Count == 1) continue;
            if (comps.Count > 2) throw new Exception($"global crack: node {v} star splits into {comps.Count}");
            int plus = -1, minus = 0;
            var mixed = false;
            foreach (var r in comps)
            {
                int pn = 0, mn = 0;
                foreach (var e in tets.Where(e => Find(e) == r))
                    foreach (var w in Rows[e])
                    {
                        var sg = Sign(w, phiS);
                        if (sg > 0) pn++;
                        else if (sg < 0) mn++;
                    }

                if (pn > 0 && mn > 0) mixed = true;
                else if (pn > 0) plus = r;
                else minus++;
            }

            if (mixed) continue;
            if (plus < 0 || minus != 1) throw new Exception($"global crack: node {v} components on one side");

            var c = NodeCount;
            SplitLog.Add((v, c));
            for (var d = 0; d < 3; d++) X.Add(X[3 * v + d]);
            On.Add(1);
            P1.Add(v);
            P2.Add(v);
            star[c] = new List<int>();
            var move = tets.Where(e => Find(e) == plus).ToList();
            foreach (var e in all.Where(e => Rows[e].Length < 4).ToList())
            {
                if (Rows[e].Length == 1)
                {
                    Rows.Add([c]);
                    star[c].Add(Rows.Count - 1);
                    continue;
                }

                var side = 0;
                foreach (var w in Rows[e])
                    if (w != v && (side = Sign(w, phiS)) != 0)
                        break;
                if (side > 0) move.Add(e);
            }

            foreach (var e in move)
            {
                Rows[e][Array.IndexOf(Rows[e], v)] = c;
                star[v].Remove(e);
                star[c].Add(e);
            }
        }
    }
}

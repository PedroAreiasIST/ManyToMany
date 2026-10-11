// CrackInserter.cs — local insertion and growth of a crack given by two level sets.
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
// A crack increment is defined by a surface function phi_s and a region function
// phi_r (the crack occupies phi_s = 0, phi_r < 0). Nodes carry a sign: 0 for nodes on
// the surface (created by earlier increments, or joined to it), otherwise the sign of
// phi_s (phi_s >= 0 counts as +). One increment
//   1. finds the crack edges: edges joining a - node and a + node whose root on the
//      segment lies in the region; the root parameter is computed from the endpoint with
//      the smaller label, so it does not depend on how the edge was found;
//   2. lets an endpoint within Tau (edge parameter) of the root of one of its crack edges
//      join the surface instead of cutting, provided no surface node was created on one
//      of its edges before (which would leave a tetrahedron with four surface nodes);
//   3. cuts the remaining crack edges (Cutter) in ascending order of (a, b), placing each
//      new node at its root with sign 0;
//   4. splits every surface node whose tetrahedral star falls into two components once the
//      adjacency across crack faces (interior faces whose three nodes have sign 0) is
//      removed, one on each side of the surface: the + component receives a new node,
//      labelled in ascending order of the split node; a component holding nodes of both
//      signs marks a front node on the domain boundary, which is not split. Elements of
//      lower dimension at the node follow the side of their other nodes; point elements
//      are copied to the new node.
// Step 1 is local: starting from the tetrahedra around the previous front (the surface
// nodes that were not split), it walks through tetrahedra that contain an uncut crack
// edge, so its cost follows the advance of the crack, not the crack or the mesh. A
// global scan is used for the first increment, or on request for verification.

using Numerical;

namespace Step2;

public sealed class CrackStep
{
    public required (int A, int B)[] Edges;   // crack edges, ascending
    public required CutRecord Cut;
    public required (int Node, int Copy)[] Splits; // split surface nodes and their new copies
    public required int[] Front;              // surface nodes left unsplit (seeds of the next increment)
    public int VisitedTets;
}

public sealed class CrackInserter : ICutObserver
{
    public delegate double Field(double x, double y, double z);

    private readonly IncrementalMesh _m;
    public readonly Attribute<double> X, Y, Z;
    public readonly Attribute<byte> OnSurface;
    public readonly Attribute<int> Parent1, Parent2;
    public readonly Attribute<double> Weight;

    /// <summary>Ancestor of every element in the mesh the inserter was created on.</summary>
    public readonly Attribute<int> Origin;

    private readonly List<Attribute<double>> _fields = new();

    /// <summary>
    ///     Registers a nodal field that follows the mesh: a node created on edge (a, b) at
    ///     parameter t takes (1 - t) u_a + t u_b, a split copy takes the value of its node. For a
    ///     field that is linear on every element this is the interpolant on the new mesh,
    ///     because every new node lies in an element of the mesh it is created on.
    /// </summary>
    public void Follow(Attribute<double> field) => _fields.Add(field);

    private Dictionary<long, double> _roots = new(EdgeKeyComparer.Instance);
    private int[] _front = [];

    /// <summary>
    ///     A crack edge whose root lies within Tau (edge parameter) of an endpoint is not cut:
    ///     the endpoint joins the surface instead (sign 0), which avoids coincident nodes and
    ///     slivers when the surface passes through or next to a node.
    /// </summary>
    public double Tau = 0.05;

    public CrackInserter(IncrementalMesh m, ReadOnlySpan<double> xyz)
    {
        _m = m;
        X = m.AddAttribute(AttributeKind.Node, 0.0);
        Y = m.AddAttribute(AttributeKind.Node, 0.0);
        Z = m.AddAttribute(AttributeKind.Node, 0.0);
        OnSurface = m.AddAttribute(AttributeKind.Node, (byte)0);
        Parent1 = m.AddAttribute(AttributeKind.Node, -1);
        Parent2 = m.AddAttribute(AttributeKind.Node, -1);
        Weight = m.AddAttribute(AttributeKind.Node, 0.0);
        Origin = m.AddAttribute(AttributeKind.Element, -1);
        for (var e = 0; e < m.ElementCount; e++) Origin.Set(e, e);
        for (var v = 0; v < m.NodeCount; v++)
        {
            X.Set(v, xyz[3 * v]);
            Y.Set(v, xyz[3 * v + 1]);
            Z.Set(v, xyz[3 * v + 2]);
            Parent1.Set(v, v);
            Parent2.Set(v, v);
        }
    }

    public int[] Front => _front;

    /// <summary>
    ///     Restores the front after a rollback of the mesh to the version it was taken from (the front is the
    ///     inserter's only state outside the versioned attributes).
    /// </summary>
    public void RestoreFront(int[] front) => _front = (int[])front.Clone();

    private int Sign(int v, Field phiS) =>
        OnSurface[v] != 0 ? 0 : phiS(X[v], Y[v], Z[v]) >= 0 ? 1 : -1;

    private double Root(int a, int b, Field phiS) =>
        RootOf(X[a], Y[a], Z[a], X[b], Y[b], Z[b], phiS);

    /// <summary>
    ///     Root parameter on the segment from a to b (a the smaller label): secant start, then
    ///     bisection on the bracket; exact for affine fields up to round-off.
    /// </summary>
    public static double RootOf(double xa, double ya, double za, double xb, double yb, double zb, Field phiS)
    {
        double dx = xb - xa, dy = yb - ya, dz = zb - za;
        double fa = phiS(xa, ya, za), fb = phiS(xb, yb, zb);
        double lo = 0, hi = 1, flo = fa;
        var t = fa / (fa - fb);
        for (var it = 0; it < 80; it++)
        {
            var ft = phiS(xa + t * dx, ya + t * dy, za + t * dz);
            if (ft == 0) return t;
            if ((ft < 0) == (flo < 0))
            {
                lo = t;
                flo = ft;
            }
            else hi = t;

            if (hi - lo < 1e-15) break;
            t = 0.5 * (lo + hi);
        }

        return 0.5 * (lo + hi);
    }

    private bool CrackEdge(int a, int b, Field phiS, Field phiR, out double t)
    {
        t = 0;
        if (a > b) (a, b) = (b, a);
        if (Sign(a, phiS) * Sign(b, phiS) != -1) return false;
        t = Root(a, b, phiS);
        double x = X[a] + t * (X[b] - X[a]), y = Y[a] + t * (Y[b] - Y[a]), z = Z[a] + t * (Z[b] - Z[a]);
        return phiR(x, y, z) < 0;
    }

    /// <summary>Crack edges by a scan of all edges (first increment, verification).</summary>
    public List<(int, int)> ScanAll(Field phiS, Field phiR, Dictionary<long, double>? roots = null)
    {
        var list = new List<(int, int)>();
        for (var v = 0; v < _m.NodeCount; v++)
        {
            var cols = _m.PatternColumns(v);
            var cnts = _m.PatternCounts(v);
            for (var q = 0; q < cols.Length; q++)
                if (cnts[q] > 0 && cols[q] > v && CrackEdge(v, cols[q], phiS, phiR, out var t))
                {
                    list.Add((v, cols[q]));
                    if (roots != null) roots[Cutter.Key(v, cols[q])] = t;
                }
        }

        list.Sort(); // pattern rows are in slot order; the cut order (and the new labels) is ascending
        return list;
    }

    /// <summary>Crack edges by a walk from the tetrahedra around the seed nodes.</summary>
    public List<(int, int)> Walk(IEnumerable<int> seeds, Field phiS, Field phiR, Dictionary<long, double> roots,
        out int visited)
    {
        var seen = new HashSet<int>();
        var stack = new Stack<int>();
        foreach (var v in seeds)
        {
            var tr = _m.TransposePairs(v);
            for (var q = 0; q < tr.Length; q += 2)
                if (_m.Row(tr[q]).Length == 4 && seen.Add(tr[q]))
                    stack.Push(tr[q]);
        }

        var found = new HashSet<long>(EdgeKeyComparer.Instance);
        var list = new List<(int, int)>();
        visited = 0;
        while (stack.Count > 0)
        {
            var g = stack.Pop();
            visited++;
            var row = _m.Row(g).ToArray();
            var active = false;
            foreach (var en in SubdivisionTemplates.EdgeNodes)
            {
                int a = Math.Min(row[en[0]], row[en[1]]), b = Math.Max(row[en[0]], row[en[1]]);
                if (!CrackEdge(a, b, phiS, phiR, out var t)) continue;
                active = true;
                if (found.Add(Cutter.Key(a, b)))
                {
                    list.Add((a, b));
                    roots[Cutter.Key(a, b)] = t;
                }
            }

            if (!active) continue;
            foreach (var v in row)
            {
                var tr = _m.TransposePairs(v);
                for (var q = 0; q < tr.Length; q += 2)
                    if (_m.Row(tr[q]).Length == 4 && seen.Add(tr[q]))
                        stack.Push(tr[q]);
            }
        }

        list.Sort();
        return list;
    }

    /// <summary>One increment. With <paramref name="global"/> the crack edges come from a scan of all edges.</summary>
    public CrackStep Advance(Field phiS, Field phiR, bool global = false)
    {
        _roots = new Dictionary<long, double>(EdgeKeyComparer.Instance);
        var visited = 0;
        var crossing = global || _front.Length == 0
            ? ScanAll(phiS, phiR, _roots)
            : Walk(_front, phiS, phiR, _roots, out visited);
        var (edges, promoted) = Promote(crossing, _roots, Tau, Promotable);
        foreach (var v in promoted) OnSurface.Set(v, 1);
        var cut = Cutter.Cut(_m, edges, this);

        // ---- candidates: surface nodes of tetrahedra around the new nodes and the previous front
        var cand = new SortedSet<int>();
        void Around(int v)
        {
            var tr = _m.TransposePairs(v);
            for (var q = 0; q < tr.Length; q += 2)
                foreach (var w in _m.Row(tr[q]))
                    if (OnSurface[w] != 0)
                        cand.Add(w);
        }

        foreach (var (_, _, c) in cut.Edges) Around(c);
        foreach (var v in promoted) Around(v);
        foreach (var v in _front) Around(v);

        var splits = new List<(int, int)>();
        var front = new List<int>();
        foreach (var v in cand)
        {
            var copy = SplitNode(v, phiS);
            if (copy >= 0) splits.Add((v, copy));
            else if (copy == FrontNode) front.Add(v);
        }

        _front = front.ToArray();
        return new CrackStep
        {
            Edges = edges.ToArray(), Cut = cut, Splits = splits.ToArray(), Front = _front, VisitedTets = visited
        };
    }

    /// <summary>
    ///     Splits the crack edges into edges to cut and endpoints that join the surface: the endpoints
    ///     within tau of the root of one of their crack edges are candidates, decided in ascending label
    ///     order (each decision sees the earlier promotions), and edges with a promoted endpoint are not
    ///     cut. Depends only on the set of crack edges, their roots and the mesh.
    /// </summary>
    public static (List<(int, int)> Cut, List<int> Promoted) Promote(List<(int, int)> crossing,
        Dictionary<long, double> roots, double tau, Func<int, IReadOnlySet<int>, bool> promotable)
    {
        var candidates = new SortedSet<int>();
        foreach (var (a, b) in crossing)
        {
            var t = roots[Cutter.Key(a, b)];
            if (t <= tau) candidates.Add(a);
            if (t >= 1 - tau) candidates.Add(b);
        }

        var promoted = new HashSet<int>();
        var list = new List<int>();
        foreach (var v in candidates)
            if (promotable(v, promoted))
            {
                promoted.Add(v);
                list.Add(v);
            }

        var cut = crossing.Where(e => !promoted.Contains(e.Item1) && !promoted.Contains(e.Item2)).ToList();
        return (cut, list);
    }

    /// <summary>
    ///     A node can join the surface unless a tetrahedron of its star would then have four surface
    ///     nodes (its four faces would all be crack faces and cut it loose). A node on the surface itself
    ///     always qualifies, since four surface nodes on a plane cannot span a tetrahedron; a refused
    ///     candidate would otherwise be cut at a root next to it.
    /// </summary>
    private bool Promotable(int v, IReadOnlySet<int> promoted)
    {
        var tr = _m.TransposePairs(v);
        for (var q = 0; q < tr.Length; q += 2)
        {
            var row = _m.Row(tr[q]);
            if (row.Length != 4) continue;
            var on = 0;
            foreach (var w in row)
                if (w != v && (OnSurface[w] != 0 || promoted.Contains(w)))
                    on++;
            if (on == 3) return false;
        }

        return true;
    }

    private const int FrontNode = -1, Interior = -2;
    private readonly List<int> _star = new();
    private readonly Dictionary<int, int> _comp = new();

    /// <summary>
    ///     Splits v when its tetrahedral star falls apart across crack faces; returns the new node, or
    ///     FrontNode when the crack ends at v, or Interior when v needs no split (its star lies on one side).
    ///     A crack face is an interior face (two tetrahedra) whose nodes all have sign 0.
    /// </summary>
    private int SplitNode(int v, Field phiS)
    {
        _star.Clear();
        var tr = _m.TransposePairs(v);
        for (var q = 0; q < tr.Length; q += 2)
            if (_m.Row(tr[q]).Length == 4)
                _star.Add(tr[q]);
        if (_star.Count == 0) return Interior;

        // union-find over the star through non-crack faces containing v
        _comp.Clear();
        foreach (var g in _star) _comp[g] = g;
        int Find(int g)
        {
            while (_comp[g] != g) g = _comp[g] = _comp[_comp[g]];
            return g;
        }

        foreach (var g in _star)
        {
            var row = _m.Row(g);
            var pv = row.IndexOf(v);
            for (var f = 0; f < 4; f++)
            {
                if (f == pv) continue; // the face opposite v does not contain v
                var nb = _m.FacetNeighbour(g, f);
                if (nb < 0 || !_comp.ContainsKey(nb)) continue;
                var crackFace = true;
                for (var p = 0; p < 4 && crackFace; p++)
                    if (p != f && OnSurface[row[p]] == 0)
                        crackFace = false;
                if (crackFace) continue;
                int ra = Find(g), rb = Find(nb);
                if (ra != rb) _comp[Math.Max(ra, rb)] = Math.Min(ra, rb);
            }
        }

        var roots = new HashSet<int>();
        foreach (var g in _star) roots.Add(Find(g));
        if (roots.Count == 1)
        {
            // one component: a front node when it holds nodes of both signs (the crack ends at v), otherwise
            // a node whose other side was split off earlier (its star lies on one side)
            bool hasPlus = false, hasMinus = false;
            foreach (var g in _star)
            foreach (var w in _m.Row(g))
            {
                var sg = Sign(w, phiS);
                hasPlus |= sg > 0;
                hasMinus |= sg < 0;
            }

            return hasPlus && hasMinus ? FrontNode : Interior;
        }
        if (roots.Count > 8)
        {
            var sb = new System.Text.StringBuilder();
            foreach (var g in _star)
            {
                sb.Append($"\n  tet {g} comp {Find(g)}:");
                foreach (var w in _m.Row(g)) sb.Append($" {w}[s{Sign(w, phiS)},on{OnSurface[w]},p{Parent1[w]}/{Parent2[w]}]");
            }

            throw new InvalidOperationException($"node {v}: star splits into {roots.Count} parts{sb}");
        }

        // sides: every node of nonzero sign in a component must lie on one side; a component
        // holding both signs wraps around the extension of the surface beyond the front (a front
        // node on the domain boundary), and the node is not split
        // A node whose star falls into more than two components is split into one node per component: the first component
        // on the minus side keeps v, every other component receives a new node (with two components this is the
        // + component, as before). Such a star arises where an earlier increment left an unsplit node between split ones.
        var plusRoots = new List<int>();
        var minusList = new List<int>();
        foreach (var r in roots)
        {
            int plusN = 0, minusN = 0;
            foreach (var g in _star)
            {
                if (Find(g) != r) continue;
                foreach (var w in _m.Row(g))
                {
                    var sg = Sign(w, phiS);
                    if (sg > 0) plusN++;
                    else if (sg < 0) minusN++;
                }
            }

            if (plusN > 0 && minusN > 0) return FrontNode;
            if (plusN == 0 && minusN == 0) throw new InvalidOperationException($"node {v}: component without a signed node");
            if (plusN > 0) plusRoots.Add(r);
            else minusList.Add(r);
        }

        if (plusRoots.Count == 0 || minusList.Count == 0)
            throw new InvalidOperationException($"node {v}: components on the same side of the surface");
        plusRoots.Sort();
        minusList.Sort();
        var others = plusRoots.Concat(minusList.Skip(1)).ToList();
        var copies = new List<int>();
        foreach (var root in others)
        {
            var c = _m.AddNode();
            X.Set(c, X[v]);
            Y.Set(c, Y[v]);
            Z.Set(c, Z[v]);
            OnSurface.Set(c, 1);
            Parent1.Set(c, v);
            Parent2.Set(c, v);
            Weight.Set(c, 0);
            foreach (var u in _fields) u.Set(c, u[v]);
            foreach (var g in _star.Where(g => Find(g) == root).ToArray()) Relabel(g, v, c);
            copies.Add(c);
        }

        var cPlus = copies[0]; // the copy of the first + component
        // lower-dimensional elements at v: by the sign of their other nodes; points are copied
        var lower = new List<int>();
        for (var q = 0; q < tr.Length; q += 2)
            if (_m.Row(tr[q]).Length < 4)
                lower.Add(tr[q]);
        foreach (var g in lower)
        {
            var row = _m.Row(g);
            if (row.Length == 1)
            {
                foreach (var c in copies)
                {
                    var pe = _m.AddElement([c]);
                    Origin.Set(pe, Origin[g]);
                }

                continue;
            }

            var side = 0;
            foreach (var w in row)
                if (w != v && (side = Sign(w, phiS)) != 0)
                    break;
            if (side > 0) Relabel(g, v, cPlus);
        }

        return cPlus;
    }

    private void Relabel(int g, int v, int c)
    {
        var row = _m.Row(g).ToArray();
        row[Array.IndexOf(row, v)] = c;
        _m.ReplaceElement(g, row);
    }

    // ---------------------------------------------------------------- ICutObserver

    void ICutObserver.NodeCreated(int node, int a, int b)
    {
        var t = _roots[Cutter.Key(a, b)];
        X.Set(node, X[a] + t * (X[b] - X[a]));
        Y.Set(node, Y[a] + t * (Y[b] - Y[a]));
        Z.Set(node, Z[a] + t * (Z[b] - Z[a]));
        OnSurface.Set(node, 1);
        Parent1.Set(node, a);
        Parent2.Set(node, b);
        Weight.Set(node, t);
        foreach (var u in _fields) u.Set(node, (1 - t) * u[a] + t * u[b]);
    }

    void ICutObserver.ElementSplit(int parent, int child) => Origin.Set(child, Origin[parent]);
}

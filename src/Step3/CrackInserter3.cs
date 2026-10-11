// CrackInserter3.cs — crack increments for several cracks: level-set edges, promotion, cut, star-component splits.
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
// The increment of crack c is defined by a surface function phi_s and a region function
// phi_r, as in Step2 (the crack occupies phi_s = 0, phi_r < 0). The sign of a node is 0
// when bit c of its crack mask is set, otherwise the sign of phi_s; nodes of other cracks
// are ordinary material nodes for this increment, and the faces of other cracks are
// already boundary facets, so an increment that reaches another crack ends on it as a
// front ends on a free surface. The four stages are those of Step2 — crack edges by a
// walk from the seeds, promotion of endpoints within Tau of a root, the cut, and the
// split of surface nodes — with two generalizations: the root of an edge is computed
// from the endpoint with the smaller canonical label, so that the two copies of an edge
// cut by an earlier crack receive coincident nodes with one canonical label; and a node
// whose star falls into k >= 2 components across the crack faces of c receives one node
// per component, the component holding the smallest element keeping the original label.
// A component holding nodes of both signs marks a front node (the crack ends at it);
// a component with no signed node, or two components on the same side (unless allowed),
// reject the increment, which the driver rolls back.

using Numerical;
using Step2;

namespace Step3;

public delegate double Field3(V3 x);

public sealed class IncrementResult
{
    public required (int A, int B)[] Edges;
    public required CutRecord Cut;
    public required (int Node, int Copy)[] Splits;
    public required int[] FrontNodes;
    public required int[] Promoted;
    public int VisitedTets;
    public int SameSideSplits;
}

public sealed class CrackInserter3
{
    private readonly CrackMesh _cm;
    private IncrementalMesh M => _cm.M;

    /// <summary>Promotion tolerance (edge parameter).</summary>
    public double Tau = 0.05;

    /// <summary>Accept splits whose components lie on the same side of the surface (folds, junction wrap-around).</summary>
    public bool AllowSameSide;

    public int MaxComponents = 8;

    /// <summary>
    ///     Nodes the front fitter placed on the wall of the new front. A crack edge joining two of them has its root on
    ///     the wall and belongs to the increment regardless of the region function, which cannot resolve the wall to
    ///     better than the fit itself; the region margin can then be small.
    /// </summary>
    public HashSet<int>? OnWall;

    /// <summary>
    ///     A node of another crack within this distance of the surface phi_s = 0 (and inside the region) joins crack c:
    ///     the surfaces merge there (coalescence of coplanar cracks, junction lines). Set per increment by the driver.
    /// </summary>
    public double CoincidenceTolerance;
    private readonly HashSet<int> _coincident = new();

    /// <summary>
    ///     Optional geometric confirmation of a crack face: its centroid must lie within this distance of phi_s = 0
    ///     (infinity, the default, disables it). The topological rule below makes it unnecessary in the tested cases.
    /// </summary>
    public double FaceTolerance = double.PositiveInfinity;

    // A crack face of increment c is an interior face whose three nodes carry bit c and none of which lies BEHIND the
    // previous front: an old node of crack c that is neither a seed (a node of the front chain, including its mouth and
    // junction ends), nor promoted, nor created now. A face across the kink between the old surface and the new strip
    // always contains such a node, while every face of the new strip has only new nodes, promoted nodes and front nodes.
    private HashSet<int> _seeds = new();
    private HashSet<int> _promoted = new();
    private int _firstNew;

    public CrackInserter3(CrackMesh cm) => _cm = cm;

    private int Sign(int v, int bit, Field3 phiS)
    {
        var mask = _cm.Mask[v];
        if ((mask & bit) != 0) return 0;
        var f = phiS(_cm.P(v));
        if (mask != 0 && CoincidenceTolerance > 0 && Math.Abs(f) < CoincidenceTolerance)
        {
            _coincident.Add(v);
            return 0;
        }

        return f >= 0 ? 1 : -1;
    }

    /// <summary>
    ///     Root parameter on the segment from a to b, bracketed from the endpoint with the smaller canonical
    ///     label: secant start, then bisection; exact for affine fields up to round-off.
    /// </summary>
    public double Root(int a, int b, Field3 phiS)
    {
        var fromA = _cm.Canon[a] <= _cm.Canon[b];
        var p = fromA ? _cm.P(a) : _cm.P(b);
        var q = fromA ? _cm.P(b) : _cm.P(a);
        var t = RootOf(p, q, phiS);
        return fromA ? t : 1 - t;
    }

    public static double RootOf(V3 p, V3 q, Field3 phiS)
    {
        var d = q - p;
        double fa = phiS(p), fb = phiS(q);
        double lo = 0, hi = 1, flo = fa;
        var t = fa / (fa - fb);
        if (!(t > 0 && t < 1)) t = 0.5;
        for (var it = 0; it < 80; it++)
        {
            var ft = phiS(p + t * d);
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

    private bool CrackEdge(int a, int b, int bit, Field3 phiS, Field3 phiR, out double t)
    {
        t = 0;
        if (a > b) (a, b) = (b, a);
        if (Sign(a, bit, phiS) * Sign(b, bit, phiS) != -1) return false;
        t = Root(a, b, phiS);
        if (OnWall != null && OnWall.Contains(a) && OnWall.Contains(b)) return true;
        return phiR(V3.Lerp(_cm.P(a), _cm.P(b), t)) < 0;
    }

    /// <summary>Crack edges by a scan of all edges (first increment of a crack, verification).</summary>
    public List<(int, int)> ScanAll(int bit, Field3 phiS, Field3 phiR, Dictionary<long, double> roots)
    {
        var list = new List<(int, int)>();
        for (var v = 0; v < M.NodeCount; v++)
        {
            var cols = M.PatternColumns(v);
            var cnts = M.PatternCounts(v);
            for (var q = 0; q < cols.Length; q++)
                if (cnts[q] > 0 && cols[q] > v && CrackEdge(v, cols[q], bit, phiS, phiR, out var t))
                {
                    list.Add((v, cols[q]));
                    roots[Cutter.Key(v, cols[q])] = t;
                }
        }

        list.Sort();
        return list;
    }

    /// <summary>Crack edges by a walk from the tetrahedra around the seed nodes through tetrahedra with a crack edge.</summary>
    public List<(int, int)> Walk(IEnumerable<int> seeds, int bit, Field3 phiS, Field3 phiR, Dictionary<long, double> roots,
        out int visited)
    {
        var seen = new HashSet<int>();
        var stack = new Stack<int>();
        foreach (var v in seeds)
        {
            var tr = M.TransposePairs(v);
            for (var q = 0; q < tr.Length; q += 2)
                if (M.Row(tr[q]).Length == 4 && seen.Add(tr[q]))
                    stack.Push(tr[q]);
        }

        var found = new HashSet<long>(EdgeKeyComparer.Instance);
        var list = new List<(int, int)>();
        visited = 0;
        var row = new int[4];
        while (stack.Count > 0)
        {
            var g = stack.Pop();
            visited++;
            M.Row(g).CopyTo(row);
            var active = false;
            foreach (var en in SubdivisionTemplates.EdgeNodes)
            {
                int a = Math.Min(row[en[0]], row[en[1]]), b = Math.Max(row[en[0]], row[en[1]]);
                if (!CrackEdge(a, b, bit, phiS, phiR, out var t)) continue;
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
                var tr = M.TransposePairs(v);
                for (var q = 0; q < tr.Length; q += 2)
                    if (M.Row(tr[q]).Length == 4 && seen.Add(tr[q]))
                        stack.Push(tr[q]);
            }
        }

        list.Sort();
        return list;
    }

    /// <summary>
    ///     One increment of crack <paramref name="bit"/> (a single bit). The crack edges come from a walk from the
    ///     seeds, or from a scan of all edges when <paramref name="global"/> is set or no seed is given.
    /// </summary>
    public IncrementResult Advance(int bit, Field3 phiS, Field3 phiR, IReadOnlyCollection<int> seeds, bool global = false)
    {
        if (bit == 0 || (bit & (bit - 1)) != 0) throw new ArgumentException("bit must be a single crack bit", nameof(bit));
        var roots = new Dictionary<long, double>(EdgeKeyComparer.Instance);
        var visited = 0;
        _coincident.Clear();
        var crossing = global || seeds.Count == 0
            ? ScanAll(bit, phiS, phiR, roots)
            : Walk(seeds, bit, phiS, phiR, roots, out visited);
        var (edges, promoted) = Promote(crossing, roots, Tau, (v, set) => Promotable(v, bit, set));
        // nodes of other cracks lying on this surface inside the region join it
        foreach (var v in _coincident.OrderBy(x => x))
            if (phiR(_cm.P(v)) < 0 && (_cm.Mask[v] & bit) == 0 && Promotable(v, bit, new HashSet<int>(promoted)))
                promoted.Add(v);
        foreach (var v in promoted) _cm.Mask.Set(v, _cm.Mask[v] | bit);
        _seeds = new HashSet<int>(seeds);
        _promoted = new HashSet<int>(promoted);
        _firstNew = M.NodeCount;
        _cm.PendingRoots = roots;
        _cm.PendingMask = bit;
        var cut = Cutter.Cut(M, edges, _cm);

        var cand = new SortedSet<int>();
        void Around(int v)
        {
            var tr = M.TransposePairs(v);
            for (var q = 0; q < tr.Length; q += 2)
                foreach (var w in M.Row(tr[q]))
                    if ((_cm.Mask[w] & bit) != 0)
                        cand.Add(w);
        }

        foreach (var (_, _, c) in cut.Edges) Around(c);
        foreach (var v in promoted) Around(v);
        foreach (var v in seeds) Around(v);

        var splits = new List<(int, int)>();
        var front = new List<int>();
        var sameSide = 0;
        foreach (var v in cand)
        {
            var r = SplitNode(v, bit, phiS, splits, ref sameSide);
            if (r == FrontNode) front.Add(v);
        }

        return new IncrementResult
        {
            Edges = edges.ToArray(), Cut = cut, Splits = splits.ToArray(), FrontNodes = front.ToArray(),
            Promoted = promoted.ToArray(), VisitedTets = visited, SameSideSplits = sameSide
        };
    }

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

    /// <summary>A node joins the surface of crack c unless a tetrahedron of its star would then have four c-nodes.</summary>
    private bool Promotable(int v, int bit, IReadOnlySet<int> promoted)
    {
        var tr = M.TransposePairs(v);
        for (var q = 0; q < tr.Length; q += 2)
        {
            var row = M.Row(tr[q]);
            if (row.Length != 4) continue;
            var on = 0;
            foreach (var w in row)
                if (w != v && ((_cm.Mask[w] & bit) != 0 || promoted.Contains(w)))
                    on++;
            if (on == 3) return false;
        }

        return true;
    }

    public const int FrontNode = -1, Interior = -2, Split = -3;
    private readonly List<int> _star = new();
    private readonly Dictionary<int, int> _comp = new();

    /// <summary>
    ///     Splits v when its tetrahedral star falls apart across the crack faces of c (interior faces whose three nodes
    ///     carry bit c): one node per component, the component with the smallest element keeping v. Returns Split,
    ///     FrontNode (a component holds both signs: the crack ends at v) or Interior.
    /// </summary>
    private int SplitNode(int v, int bit, Field3 phiS, List<(int, int)> splits, ref int sameSide)
    {
        _star.Clear();
        var tr = M.TransposePairs(v);
        for (var q = 0; q < tr.Length; q += 2)
            if (M.Row(tr[q]).Length == 4)
                _star.Add(tr[q]);
        if (_star.Count == 0) return Interior;

        _comp.Clear();
        foreach (var g in _star) _comp[g] = g;
        int Find(int g)
        {
            while (_comp[g] != g) g = _comp[g] = _comp[_comp[g]];
            return g;
        }

        foreach (var g in _star)
        {
            var row = M.Row(g);
            var pv = row.IndexOf(v);
            for (var f = 0; f < 4; f++)
            {
                if (f == pv) continue;
                var nb = M.FacetNeighbour(g, f);
                if (nb < 0 || !_comp.ContainsKey(nb)) continue;
                var crackFace = true;
                for (var p = 0; p < 4 && crackFace; p++)
                {
                    if (p == f) continue;
                    var w = row[p];
                    if ((_cm.Mask[w] & bit) == 0) crackFace = false;
                    else if (w < _firstNew && !_promoted.Contains(w) && _seeds.Count > 0 && !_seeds.Contains(w)) crackFace = false;
                }

                if (crackFace && !double.IsPositiveInfinity(FaceTolerance))
                {
                    var cen = V3.Zero;
                    for (var p = 0; p < 4; p++)
                        if (p != f)
                            cen += _cm.P(row[p]) / 3;
                    if (Math.Abs(phiS(cen)) > FaceTolerance) crackFace = false;
                }

                if (crackFace) continue;
                int ra = Find(g), rb = Find(nb);
                if (ra != rb) _comp[Math.Max(ra, rb)] = Math.Min(ra, rb);
            }
        }

        var roots = new SortedSet<int>();
        foreach (var g in _star) roots.Add(Find(g));
        if (roots.Count == 1)
        {
            bool hasPlus = false, hasMinus = false;
            foreach (var g in _star)
            foreach (var w in M.Row(g))
            {
                var sg = Sign(w, bit, phiS);
                hasPlus |= sg > 0;
                hasMinus |= sg < 0;
            }

            return hasPlus && hasMinus ? FrontNode : Interior;
        }

        if (roots.Count > MaxComponents)
            throw new InvalidOperationException($"node {v}: star splits into {roots.Count} parts");

        // sides of the components
        var plusCount = new Dictionary<int, int>();
        var minusCount = new Dictionary<int, int>();
        foreach (var r in roots)
        {
            plusCount[r] = 0;
            minusCount[r] = 0;
        }

        foreach (var g in _star)
        {
            var r = Find(g);
            foreach (var w in M.Row(g))
            {
                var sg = Sign(w, bit, phiS);
                if (sg > 0) plusCount[r]++;
                else if (sg < 0) minusCount[r]++;
            }
        }

        int nPlus = 0, nMinus = 0;
        foreach (var r in roots)
        {
            if (plusCount[r] > 0 && minusCount[r] > 0) return FrontNode;
            if (plusCount[r] == 0 && minusCount[r] == 0)
                throw new InvalidOperationException($"node {v}: component without a signed node");
            if (plusCount[r] > 0) nPlus++;
            else nMinus++;
        }

        if (nPlus == 0 || nMinus == 0)
        {
            if (!AllowSameSide)
                throw new InvalidOperationException($"node {v}: {roots.Count} components on the same side of the surface");
            sameSide++;
        }
        else if (roots.Count > 2) sameSide += roots.Count - 2;

        // the component holding the smallest element keeps v (roots are the smallest elements of their components)
        var keep = roots.Min;
        var copies = new List<(int Root, int Node)>();
        foreach (var r in roots)
        {
            if (r == keep) continue;
            var c = _cm.Copy(v);
            foreach (var g in _star.Where(g => Find(g) == r).ToArray()) Relabel(g, v, c);
            copies.Add((r, c));
            splits.Add((v, c));
        }

        // lower-dimensional elements at v: a point is copied to every copy; a bar or triangle follows the
        // component of a tetrahedron that contains its other nodes
        var lower = new List<int>();
        for (var q = 0; q < tr.Length; q += 2)
            if (M.Row(tr[q]).Length < 4)
                lower.Add(tr[q]);
        foreach (var g in lower)
        {
            var row = M.Row(g).ToArray();
            if (row.Length == 1)
            {
                foreach (var (_, c) in copies)
                {
                    var pe = M.AddElement([c]);
                    _cm.Origin.Set(pe, _cm.Origin[g]);
                }

                continue;
            }

            foreach (var (r, c) in copies)
            {
                var found = false;
                foreach (var tet in _star)
                {
                    if (Find(tet) != r) continue;
                    var trow = M.Row(tet);
                    var all = true;
                    foreach (var w in row)
                        if (w != v && trow.IndexOf(w) < 0)
                        {
                            all = false;
                            break;
                        }

                    if (all)
                    {
                        found = true;
                        break;
                    }
                }

                if (found)
                {
                    Relabel(g, v, c);
                    break;
                }
            }
        }

        return Split;
    }

    private void Relabel(int g, int v, int c)
    {
        var row = M.Row(g).ToArray();
        row[Array.IndexOf(row, v)] = c;
        M.ReplaceElement(g, row);
    }
}

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
    public int Pushed;
}

public sealed class CrackInserter3
{
    private readonly CrackMesh _cm;
    private IncrementalMesh M => _cm.M;

    /// <summary>Promotion tolerance (edge parameter).</summary>
    public double Tau = 0.05;

    /// <summary>
    ///     A node within Tau of a root that cannot join the surface (a tetrahedron of its star would have its four nodes on
    ///     the crack) is pushed off the surface by this fraction of its edge length, away from the surface and within the
    ///     outer boundary, so that the cuts on its edges are not slivers (0 = no push).
    /// </summary>
    public double PushElements = 0.1;

    /// <summary>Forbids the promotion of a node when a tetrahedron of its star would get four crack nodes (the earlier rule).</summary>
    public bool ForbidFourOnCrack = true;

    /// <summary>
    ///     Components of a star all on one side of the surface reject the increment (the driver retries with the schedule);
    ///     false leaves the node unsplit instead.
    /// </summary>
    public bool RejectSameSide = true;
    public readonly HashSet<int> LastPushed = new();

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
    ///     Extra region margin for crack edges that touch no wall node: where the fit could not place a node on the wall,
    ///     the root nearest to the wall lies up to this distance beyond it and is still cut (an outward bulge of the
    ///     front, as in Step2, instead of a notch). Edges with a wall node keep the small margin of the region function.
    /// </summary>
    public double NonWallMargin;

    /// <summary>Candidates the fit could not place on the wall; the extra margin applies to their crack edges only.</summary>
    public HashSet<int>? NearWall;

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
    // previous front. A node of crack c is behind the front when it was split by an earlier increment (its canonical
    // label has several nodes) and it is neither a seed (a node of the front chain, including its mouth and junction
    // ends, which may be split by other cracks), nor promoted, nor created now. An unsplit old node of crack c lies at
    // the front (a chain node or a pinch node the previous increment left unsplit next to its wall). A face across the
    // kink between the old surface and the new strip always contains a split node.
    private HashSet<int> _seeds = new();
    private HashSet<int> _promoted = new();
    private HashSet<int> _oldSplit = new();
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
        if (OnWall != null)
        {
            // nodes on the wall and candidates within the move tolerance of it that could not be moved count alike: an
            // edge joining two of them has its root within a fraction of the move tolerance of the wall and is cut (the
            // front bulges by that fraction there instead of lagging by an element)
            bool wa = OnWall.Contains(a) || (NearWall?.Contains(a) ?? false), wb = OnWall.Contains(b) || (NearWall?.Contains(b) ?? false);
            if (wa && wb) return true;
            return phiR(V3.Lerp(_cm.P(a), _cm.P(b), t)) - NonWallMargin < 0;
        }

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
        // each edge is tested once: the test is a deterministic function of the two nodes (no node moves during the walk),
        // and an edge is reached from every tetrahedron around it (2026-10-12: the walk took seconds on 4e5 tetrahedra)
        var tested = new Dictionary<long, double>(EdgeKeyComparer.Instance);
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
                var key = Cutter.Key(a, b);
                if (tested.TryGetValue(key, out var seenT))
                {
                    if (double.IsNaN(seenT)) continue;
                    active = true;
                    continue;
                }

                if (!CrackEdge(a, b, bit, phiS, phiR, out var t))
                {
                    tested[key] = double.NaN;
                    continue;
                }

                tested[key] = t;
                active = true;
                if (found.Add(key))
                {
                    list.Add((a, b));
                    roots[key] = t;
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
        var dbg = Environment.GetEnvironmentVariable("STEP3_TIMING_DEBUG") == "1";
        var pw = System.Diagnostics.Stopwatch.StartNew();
        void Mark(string what) { if (dbg) { Console.WriteLine($"    inserter {what}: {pw.Elapsed.TotalMilliseconds:F1} ms"); pw.Restart(); } }
        var crossing = global || seeds.Count == 0
            ? ScanAll(bit, phiS, phiR, roots)
            : Walk(seeds, bit, phiS, phiR, roots, out visited);
        Mark($"crossing ({crossing.Count} edges, {visited} visited)");
        var (edges, promoted) = Promote(crossing, roots, Tau, (v, set) => Promotable(v, bit, set));
        Mark("promote");
        LastPushed.Clear();
        var pushed = PushElements > 0 ? PushBlocked(crossing, roots, promoted, bit, phiS, edges) : 0;
        // nodes of other cracks lying on this surface inside the region join it
        foreach (var v in _coincident.OrderBy(x => x))
            if (phiR(_cm.P(v)) < 0 && (_cm.Mask[v] & bit) == 0 && Promotable(v, bit, new HashSet<int>(promoted)))
                promoted.Add(v);
        Mark("push + coincident");
        foreach (var v in promoted) _cm.Mask.Set(v, _cm.Mask[v] | bit);
        _seeds = new HashSet<int>(seeds);
        _promoted = new HashSet<int>(promoted);
        _oldSplit = new HashSet<int>();
        var seenLabel = new HashSet<int>();
        for (var v = 0; v < M.NodeCount; v++)
            if ((_cm.Mask[v] & bit) != 0 && !seenLabel.Add(_cm.Canon[v]))
                _oldSplit.Add(_cm.Canon[v]);
        Mark("old-split scan");
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
        Mark("edge cut + candidates");
        foreach (var v in cand)
        {
            var r = SplitNode(v, bit, phiS, splits, ref sameSide);
            if (r == FrontNode) front.Add(v);
        }

        Mark($"split ({splits.Count} splits)");
        return new IncrementResult
        {
            Edges = edges.ToArray(), Cut = cut, Splits = splits.ToArray(), FrontNodes = front.ToArray(),
            Promoted = promoted.ToArray(), VisitedTets = visited, SameSideSplits = sameSide, Pushed = pushed };
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

    /// <summary>
    ///     Pushes the nodes within Tau of a root that were not promoted (see <see cref="PushElements"/>) and re-roots their
    ///     crossing edges; an edge that no longer crosses is dropped. Returns the number of nodes moved.
    /// </summary>
    private int PushBlocked(List<(int, int)> crossing, Dictionary<long, double> roots, List<int> promoted, int bit, Field3 phiS, List<(int, int)> edges)
    {
        var prom = new HashSet<int>(promoted);
        var blocked = new SortedSet<int>();
        foreach (var (a, b) in crossing)
        {
            var t = roots[Cutter.Key(a, b)];
            if (t <= Tau && !prom.Contains(a) && _cm.Mask[a] == 0) blocked.Add(a);
            if (t >= 1 - Tau && !prom.Contains(b) && _cm.Mask[b] == 0) blocked.Add(b);
        }

        var pushed = 0;
        foreach (var v in blocked)
        {
            var p = _cm.P(v);
            var tr = M.TransposePairs(v);
            double lsum = 0;
            var lcnt = 0;
            var bn = new List<V3>();
            for (var q = 0; q < tr.Length; q += 2)
            {
                var e = tr[q];
                var row = M.Row(e);
                if (row.Length != 4) continue;
                var pv = row.IndexOf(v);
                for (var k = 0; k < 4; k++)
                    if (k != pv)
                    {
                        lsum += (_cm.P(row[k]) - p).Norm;
                        lcnt++;
                    }

                for (var f = 0; f < 4; f++)
                {
                    if (f == pv || M.FacetNeighbour(e, f) >= 0) continue;
                    var o = _cm.P(row[f]);
                    var fp = new V3[3];
                    var c = 0;
                    for (var k = 0; k < 4; k++)
                        if (k != f)
                            fp[c++] = _cm.P(row[k]);
                    var nrm = V3.Cross(fp[1] - fp[0], fp[2] - fp[0]);
                    if (nrm.Norm < 1e-300) continue;
                    nrm = nrm.Unit();
                    if (V3.Dot(nrm, (fp[0] + fp[1] + fp[2]) / 3 - o) < 0) nrm = -nrm;
                    bn.Add(nrm);
                }
            }

            if (lcnt == 0) continue;
            var h = lsum / lcnt;
            var f0 = phiS(p);
            var d = 1e-3 * h;
            var g = new V3((phiS(p + new V3(d, 0, 0)) - phiS(p - new V3(d, 0, 0))) / (2 * d),
                (phiS(p + new V3(0, d, 0)) - phiS(p - new V3(0, d, 0))) / (2 * d),
                (phiS(p + new V3(0, 0, d)) - phiS(p - new V3(0, 0, d))) / (2 * d));
            if (g.Norm < 0.5) continue;
            var dir = (f0 >= 0 ? 1 : -1) * g.Unit();
            V3? n1 = null;
            foreach (var n in bn)
            {
                if (n1 == null)
                {
                    n1 = n;
                    dir -= V3.Dot(dir, n) * n;
                    continue;
                }

                if (Math.Abs(V3.Dot(n, n1.Value)) > 0.99) continue;
                var n2 = (n - V3.Dot(n, n1.Value) * n1.Value).Unit();
                dir -= V3.Dot(dir, n2) * n2;
            }

            if (dir.Norm < 0.3) continue;
            dir = dir.Unit();
            var moved = false;
            for (var push = PushElements * h; push > 0.02 * h && !moved; push *= 0.5)
            {
                var target = p + push * dir;
                if (Math.Abs(phiS(target)) <= Math.Abs(f0)) break;
                var ok = true;
                for (var q = 0; q < tr.Length && ok; q += 2)
                {
                    var row = M.Row(tr[q]);
                    if (row.Length != 4) continue;
                    V3[] pts = [_cm.P(row[0]), _cm.P(row[1]), _cm.P(row[2]), _cm.P(row[3])];
                    var before = Geo.Volume6(pts[0], pts[1], pts[2], pts[3]);
                    pts[row.IndexOf(v)] = target;
                    var after = Geo.Volume6(pts[0], pts[1], pts[2], pts[3]);
                    if (Math.Sign(after) != Math.Sign(before) || Math.Abs(after) < 0.2 * Math.Abs(before)) ok = false;
                }

                if (!ok) continue;
                _cm.X.Set(v, target.X);
                _cm.Y.Set(v, target.Y);
                _cm.Z.Set(v, target.Z);
                moved = true;
            }

            if (!moved) continue;
            pushed++;
            if (Environment.GetEnvironmentVariable("STEP3_DEBUG") == "3")
                Console.WriteLine($"    pushed node {v} from {p} to {_cm.P(v)} (|phi| {Math.Abs(f0):G3} -> {Math.Abs(phiS(_cm.P(v))):G3}, h {h:G3}, boundary normals {bn.Count})");
            LastPushed.Add(v);
            for (var i = edges.Count - 1; i >= 0; i--)
            {
                var (a, b) = edges[i];
                if (a != v && b != v) continue;
                if (Sign(a, bit, phiS) * Sign(b, bit, phiS) != -1)
                {
                    edges.RemoveAt(i);
                    continue;
                }

                roots[Cutter.Key(a, b)] = Root(a, b, phiS);
            }
        }

        return pushed;
    }

    /// <summary>
    ///     A node within Tau of a root joins the surface of crack c. A tetrahedron of its star may then have its four
    ///     nodes on the crack (it lies in a fold of the surface, or flat on it); the split classifies the faces of such a
    ///     tetrahedron by the sides of the two tetrahedra at the face, see <see cref="SideOf"/>; a component left without a signed node
    ///     rejects the increment, which the driver retries with a shifted surface.
    /// </summary>
    private bool Promotable(int v, int bit, IReadOnlySet<int> promoted)
    {
        if (!ForbidFourOnCrack) return true;
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

    private bool FourOnCrack(ReadOnlySpan<int> row, int bit)
    {
        if (row.Length != 4) return false;
        foreach (var w in row)
            if ((_cm.Mask[w] & bit) == 0)
                return false;
        return true;
    }

    /// <summary>Side of a tetrahedron: that of its first signed node, or of its centroid when all its nodes are on the crack.</summary>
    private int SideOf(int t, int bit, Field3 phiS)
    {
        var row = M.Row(t);
        foreach (var w in row)
        {
            var sg = Sign(w, bit, phiS);
            if (sg != 0) return sg;
        }

        var cen = V3.Zero;
        foreach (var w in row) cen += _cm.P(w) / row.Length;
        return phiS(cen) >= 0 ? 1 : -1;
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
                    else if (w < _firstNew && !_promoted.Contains(w) && _seeds.Count > 0 && !_seeds.Contains(w) && _oldSplit.Contains(_cm.Canon[w])) crackFace = false;
                }

                if (crackFace && !ForbidFourOnCrack && (FourOnCrack(row, bit) || FourOnCrack(M.Row(nb), bit)))
                {
                    // a tetrahedron with its four nodes on the crack lies in a fold of the surface (two nodes on each
                    // sheet): it belongs to the side of its centroid, and its face is a crack face only toward a
                    // tetrahedron of the other side
                    crackFace = SideOf(g, bit, phiS) != SideOf(nb, bit, phiS);
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

        // a component without a signed node (tetrahedra with their four nodes on the crack, flat on it or in a fold)
        // joins the neighbouring component across its largest shared face; a component that has no neighbour rejects
        for (var guard = 0; guard < 8; guard++)
        {
            var unsigned = roots.FirstOrDefault(r => plusCount[r] == 0 && minusCount[r] == 0, -1);
            if (unsigned < 0) break;
            var bestArea = -1.0;
            var bestOther = -1;
            foreach (var g in _star)
            {
                if (Find(g) != unsigned) continue;
                var row = M.Row(g);
                var pv = row.IndexOf(v);
                for (var f = 0; f < 4; f++)
                {
                    if (f == pv) continue;
                    var nb = M.FacetNeighbour(g, f);
                    if (nb < 0 || !_comp.ContainsKey(nb) || Find(nb) == unsigned) continue;
                    V3 a0 = _cm.P(row[(f + 1) % 4]), a1 = _cm.P(row[(f + 2) % 4]), a2 = _cm.P(row[(f + 3) % 4]);
                    var area = V3.Cross(a1 - a0, a2 - a0).Norm;
                    if (area > bestArea)
                    {
                        bestArea = area;
                        bestOther = Find(nb);
                    }
                }
            }

            if (bestOther < 0)
            {
                // no neighbour in the star (flat tetrahedra in the crack plane, bounded by crack facets on both faces):
                // the side is that of the tetrahedra across the faces opposite v, and the component joins the star
                // component of that side
                var sd = 0;
                foreach (var g in _star)
                {
                    if (Find(g) != unsigned || sd != 0) continue;
                    var row = M.Row(g);
                    var nb = M.FacetNeighbour(g, row.IndexOf(v));
                    if (nb >= 0) sd = SideOf(nb, bit, phiS);
                }

                if (sd != 0)
                    foreach (var r in roots)
                        if (r != unsigned && (sd > 0 ? plusCount[r] > 0 : minusCount[r] > 0))
                        {
                            bestOther = r;
                            break;
                        }
            }

            if (bestOther < 0)
            {
                if (Environment.GetEnvironmentVariable("STEP3_DEBUG") == "3")
                {
                    Console.WriteLine($"    unsigned component at node {v} (mask {_cm.Mask[v]}, canon {_cm.Canon[v]}, at {_cm.P(v)}, phi {phiS(_cm.P(v)):G3}), {roots.Count} components:");
                    foreach (var g in _star)
                    {
                        if (Find(g) != unsigned) continue;
                        var row = M.Row(g);
                        var nodes = string.Join(" ", row.ToArray().Select(w => $"{w}[m{_cm.Mask[w]} s{Sign(w, bit, phiS)} phi {phiS(_cm.P(w)):G2}{(w >= _firstNew ? " new" : "")}{(_promoted.Contains(w) ? " prom" : "")}]"));
                        var faces = string.Join(" ", Enumerable.Range(0, 4).Select(f => { var nb = M.FacetNeighbour(g, f); return nb < 0 ? "B" : !_comp.ContainsKey(nb) ? "out" : $"c{Find(nb)}"; }));
                        Console.WriteLine($"      tet {g}: {nodes} faces {faces}");
                    }
                }

                throw new InvalidOperationException($"node {v}: component without a signed node");
            }
            _comp[Math.Max(unsigned, bestOther)] = Math.Min(unsigned, bestOther);
            var keepR = Math.Min(unsigned, bestOther);
            var dropR = Math.Max(unsigned, bestOther);
            plusCount[keepR] += plusCount[dropR];
            minusCount[keepR] += minusCount[dropR];
            roots.Remove(dropR);
            if (roots.Count == 1) return plusCount[keepR] > 0 && minusCount[keepR] > 0 ? FrontNode : Interior;
        }

        int nPlus = 0, nMinus = 0;
        foreach (var r in roots)
        {
            if (plusCount[r] > 0 && minusCount[r] > 0) return FrontNode;
            if (plusCount[r] > 0) nPlus++;
            else nMinus++;
        }

        if (nPlus == 0 || nMinus == 0)
        {
            // every component lies on one side: the faces between them are not a crack through v (a chord across a
            // kink, or the surface touching the star from one side); v is not split
            if (RejectSameSide && !AllowSameSide)
                throw new InvalidOperationException($"node {v}: {roots.Count} components on the same side of the surface");
            if (!AllowSameSide) return Interior;
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

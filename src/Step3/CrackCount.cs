// CrackCount.cs — crack facets, fronts, mouths, junctions and sheets recovered from the mesh by counting.
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
// A boundary facet of the cut mesh is a crack facet when another boundary facet has the
// same canonical triangle (its twin, on the other side of the crack); a canonical
// triangle with one facet is an outer facet. Each canonical edge counts the canonical
// crack triangles that contain it: one makes a front edge, or a mouth edge when the edge
// also belongs to an outer facet; two make an interior crack edge; three or more make a
// junction edge, where crack surfaces meet. Sheets are the components of the crack
// triangles connected through interior edges, so a T-junction separates the sheets it
// joins. Front chains are the paths of front edges, oriented by the frame (N, B, T) of
// their triangles: B the oriented normal (positive side), N the in-plane direction away
// from the triangle, T = N x B. Orientations are propagated across interior edges from
// the triangles oriented in earlier scans, or seeded from the reference normal of the
// crack whose bit the triangles carry. The scan rebuilds everything from the boundary
// facets of a committed version, so a rollback needs no bookkeeping; events (merged or
// split sheets, new junctions, breakthroughs, vanished fronts) are differences of scans.

using Step2;

namespace Step3;

public enum EndKind : byte
{
    Closed,
    Mouth,
    Junction,
    Free
}

public sealed class FrontChain
{
    public int Sheet;
    public int Bits;
    public bool Closed;
    public required int[] Nodes; // canonical labels in order (a closed chain does not repeat its first node)
    public required V3[] P, N, B, T;
    public required double[] S;  // arc length at the nodes (closed: S[^1] + last segment = length)
    public EndKind Start = EndKind.Closed, End = EndKind.Closed;
    public double Length;

    public int Count => Nodes.Length;

    /// <summary>Segments: pairs of consecutive node indices (closing segment included for closed chains).</summary>
    public IEnumerable<(int, int)> Segments()
    {
        for (var j = 0; j + 1 < Nodes.Length; j++) yield return (j, j + 1);
        if (Closed && Nodes.Length > 2) yield return (Nodes.Length - 1, 0);
    }
}

public sealed class CrackCount
{
    private readonly CrackMesh _cm;
    private readonly Dictionary<(int, int, int), (int, int, int)> _orient = new();

    /// <summary>Reference normal (and a point) per crack bit, used to seed the orientation of a new sheet.</summary>
    public readonly Dictionary<int, (V3 Point, V3 Normal)> Reference = new();

    public Snapshot? Snap;
    public int PairCount;
    public int[] FacetPlus = [], FacetMinus = [];          // facet ids (e * 4 + f) of the twins on the + and - side
    public (int, int, int)[] Tri = [];                     // oriented canonical triangles (+ side = normal side)
    public int[] TriBits = [];
    public int[] SheetOfPair = [];
    public int SheetCount;
    public int[] SheetBits = [];
    public readonly Dictionary<int, int> PairOfFacet = new();
    public readonly Dictionary<long, int> EdgeCount = new();
    public readonly HashSet<long> OuterEdges = new();
    public readonly HashSet<int> OuterNodes = new();
    public readonly List<(int, int)> FrontEdges = new(), MouthEdges = new(), JunctionEdges = new();
    public readonly List<FrontChain> Chains = new();
    public readonly Dictionary<int, int> Rep = new();      // canonical label -> a mesh node with that label
    public readonly Dictionary<int, int> PairOfFrontEdge = new(); // (edge key hash) handled through _frontTri
    private readonly Dictionary<long, int> _frontTri = new();
    public int OrientationConflicts, SideConflicts, UnseededSheets;
    public double ScanMs;

    public CrackCount(CrackMesh cm) => _cm = cm;

    public static long EdgeKey(int a, int b) => a < b ? ((long)a << 32) | (uint)b : ((long)b << 32) | (uint)a;
    private static (int, int, int) Sorted(int a, int b, int c)
    {
        if (a > b) (a, b) = (b, a);
        if (b > c) (b, c) = (c, b);
        if (a > b) (a, b) = (b, a);
        return (a, b, c);
    }

    public V3 Point(int canonical) => _cm.P(Rep[canonical]);

    /// <summary>Recounts everything from the boundary facets of a committed version.</summary>
    public void Rescan(Snapshot s)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        Snap = s;
        var canon = s.View(_cm.Canon);
        var mask = s.View(_cm.Mask);
        PairOfFacet.Clear();
        EdgeCount.Clear();
        OuterEdges.Clear();
        OuterNodes.Clear();
        FrontEdges.Clear();
        MouthEdges.Clear();
        JunctionEdges.Clear();
        Chains.Clear();
        Rep.Clear();
        _frontTri.Clear();
        OrientationConflicts = SideConflicts = UnseededSheets = 0;

        // ---- group the boundary facets by canonical triangle
        var byTri = new Dictionary<(int, int, int), (int F1, int F2)>(s.BoundaryCount);
        var tri = new int[3];
        for (var i = 0; i < s.BoundaryCount; i++)
        {
            var id = s.BoundaryFacet(i);
            int e = id / 4, fl = id % 4;
            var row = s.Row(e);
            if (row.Length != 4) continue;
            var k = 0;
            for (var p = 0; p < 4; p++)
                if (p != fl)
                    tri[k++] = row[p];
            foreach (var v in tri) Rep.TryAdd(canon[v], v);
            var key = Sorted(canon[tri[0]], canon[tri[1]], canon[tri[2]]);
            if (byTri.TryGetValue(key, out var fs))
            {
                if (fs.F2 >= 0) throw new InvalidOperationException($"canonical triangle {key} has more than two boundary facets");
                byTri[key] = (fs.F1, id);
            }
            else byTri[key] = (id, -1);
        }

        // ---- pairs (crack facets) and outer facets
        var pairs = new List<((int, int, int) Key, int F1, int F2)>();
        foreach (var kv in byTri)
        {
            if (kv.Value.F2 < 0)
            {
                var (a, b, c) = kv.Key;
                OuterEdges.Add(EdgeKey(a, b));
                OuterEdges.Add(EdgeKey(b, c));
                OuterEdges.Add(EdgeKey(a, c));
                OuterNodes.Add(a);
                OuterNodes.Add(b);
                OuterNodes.Add(c);
            }
            else pairs.Add((kv.Key, kv.Value.F1, kv.Value.F2));
        }

        pairs.Sort((x, y) => x.Key.CompareTo(y.Key)); // deterministic order
        PairCount = pairs.Count;
        FacetPlus = new int[PairCount];
        FacetMinus = new int[PairCount];
        Tri = new (int, int, int)[PairCount];
        TriBits = new int[PairCount];
        var trisOfEdge = new Dictionary<long, List<int>>(3 * PairCount);
        for (var i = 0; i < PairCount; i++)
        {
            var (key, f1, f2) = pairs[i];
            PairOfFacet[f1] = i;
            PairOfFacet[f2] = i;
            TriBits[i] = FacetBits(s, mask, f1) | FacetBits(s, mask, f2);
            foreach (var ek in new[] { EdgeKey(key.Item1, key.Item2), EdgeKey(key.Item2, key.Item3), EdgeKey(key.Item1, key.Item3) })
            {
                if (!trisOfEdge.TryGetValue(ek, out var l)) trisOfEdge[ek] = l = new List<int>(2);
                l.Add(i);
            }
        }

        foreach (var kv in trisOfEdge) EdgeCount[kv.Key] = kv.Value.Count;

        // ---- sheets: components through interior (count 2) edges
        var parent = new int[PairCount];
        for (var i = 0; i < PairCount; i++) parent[i] = i;
        int Find(int a)
        {
            while (parent[a] != a) a = parent[a] = parent[parent[a]];
            return a;
        }

        foreach (var kv in trisOfEdge)
            if (kv.Value.Count == 2)
            {
                int ra = Find(kv.Value[0]), rb = Find(kv.Value[1]);
                if (ra != rb) parent[Math.Max(ra, rb)] = Math.Min(ra, rb);
            }

        SheetOfPair = new int[PairCount];
        var sheetId = new Dictionary<int, int>();
        for (var i = 0; i < PairCount; i++)
        {
            var r = Find(i);
            if (!sheetId.TryGetValue(r, out var sid)) sheetId[r] = sid = sheetId.Count;
            SheetOfPair[i] = sid;
        }

        SheetCount = sheetId.Count;
        SheetBits = new int[SheetCount];
        for (var i = 0; i < PairCount; i++) SheetBits[SheetOfPair[i]] |= TriBits[i];

        // ---- orientation: keep remembered orientations, propagate through interior edges, seed new sheets
        var oriented = new bool[PairCount];
        var queue = new Queue<int>();
        for (var i = 0; i < PairCount; i++)
            if (_orient.TryGetValue(pairs[i].Key, out var o))
            {
                Tri[i] = o;
                oriented[i] = true;
                queue.Enqueue(i);
            }

        void Propagate()
        {
            while (queue.Count > 0)
            {
                var i = queue.Dequeue();
                var (a, b, c) = Tri[i];
                foreach (var (p, q) in new[] { (a, b), (b, c), (c, a) })
                {
                    var l = trisOfEdge[EdgeKey(p, q)];
                    if (l.Count != 2) continue;
                    var j = l[0] == i ? l[1] : l[0];
                    var want = OrientWith(pairs[j].Key, q, p); // the neighbour traverses the edge the other way
                    if (oriented[j])
                    {
                        if (Tri[j] != want && !SameCyclic(Tri[j], want)) OrientationConflicts++;
                        continue;
                    }

                    Tri[j] = want;
                    oriented[j] = true;
                    queue.Enqueue(j);
                }
            }
        }

        Propagate();
        for (var i = 0; i < PairCount; i++)
        {
            if (oriented[i]) continue;
            // a sheet without remembered orientation: seed from the reference normal of its crack bit
            var (a, b, c) = pairs[i].Key;
            V3 pa = Point(a), pb = Point(b), pc = Point(c);
            var n = V3.Cross(pb - pa, pc - pa);
            var bits = SheetBits[SheetOfPair[i]];
            V3? refN = null;
            for (var bit = 1; bit != 0 && bit <= bits; bit <<= 1)
                if ((bits & bit) != 0 && Reference.TryGetValue(bit, out var r))
                {
                    refN = r.Normal;
                    break;
                }

            if (refN == null)
            {
                UnseededSheets++;
                refN = n.Z != 0 ? new V3(0, 0, Math.Sign(n.Z)) : n.Y != 0 ? new V3(0, Math.Sign(n.Y), 0) : new V3(1, 0, 0);
            }

            Tri[i] = V3.Dot(n, refN.Value) >= 0 ? (a, b, c) : (a, c, b);
            oriented[i] = true;
            queue.Enqueue(i);
            Propagate();
        }

        _orient.Clear();
        for (var i = 0; i < PairCount; i++) _orient[pairs[i].Key] = Tri[i];

        // ---- twins by side of the oriented normal
        for (var i = 0; i < PairCount; i++)
        {
            var (f1, f2) = (pairs[i].F1, pairs[i].F2);
            var (a, b, c) = Tri[i];
            V3 pa = Point(a), pb = Point(b), pc = Point(c);
            var nrm = V3.Cross(pb - pa, pc - pa);
            var cen = (pa + pb + pc) / 3;
            double s1 = V3.Dot(_cm.P(s.Row(f1 / 4)[f1 % 4]) - cen, nrm), s2 = V3.Dot(_cm.P(s.Row(f2 / 4)[f2 % 4]) - cen, nrm);
            if (s1 * s2 >= 0) SideConflicts++;
            if (s1 >= s2) { FacetPlus[i] = f1; FacetMinus[i] = f2; }
            else { FacetPlus[i] = f2; FacetMinus[i] = f1; }
        }

        // ---- edge classification
        foreach (var kv in trisOfEdge)
        {
            int a = (int)(kv.Key >> 32), b = (int)(kv.Key & 0xffffffff);
            if (kv.Value.Count == 1)
            {
                if (OuterEdges.Contains(kv.Key)) MouthEdges.Add((a, b));
                else
                {
                    FrontEdges.Add((a, b));
                    _frontTri[kv.Key] = kv.Value[0];
                }
            }
            else if (kv.Value.Count >= 3) JunctionEdges.Add((a, b));
        }

        FrontEdges.Sort();
        MouthEdges.Sort();
        JunctionEdges.Sort();
        BuildChains();
        ScanMs = sw.Elapsed.TotalMilliseconds;
    }

    private int FacetBits(Snapshot s, PView<int> mask, int facet)
    {
        int e = facet / 4, fl = facet % 4;
        var row = s.Row(e);
        var bits = -1;
        for (var p = 0; p < 4; p++)
            if (p != fl)
                bits &= mask[row[p]];
        return bits;
    }

    private static (int, int, int) OrientWith((int, int, int) key, int p, int q)
    {
        // oriented triple traversing p -> q
        var (a, b, c) = key;
        var r = a != p && a != q ? a : b != p && b != q ? b : c;
        return (p, q, r);
    }

    private static bool SameCyclic((int, int, int) x, (int, int, int) y) =>
        x == y || x == (y.Item2, y.Item3, y.Item1) || x == (y.Item3, y.Item1, y.Item2);

    /// <summary>Pair index of the triangle of a front edge, or -1.</summary>
    public int FrontTriangle(int a, int b) => _frontTri.TryGetValue(EdgeKey(a, b), out var i) ? i : -1;

    /// <summary>Frame (N, B, T) of a front edge from its single triangle: B the oriented normal, N away from the triangle.</summary>
    public (V3 N, V3 B, V3 T) EdgeFrame(int a, int b)
    {
        var i = FrontTriangle(a, b);
        var (c0, c1, c2) = Tri[i];
        V3 p0 = Point(c0), p1 = Point(c1), p2 = Point(c2);
        var bN = V3.Cross(p1 - p0, p2 - p0).Unit();
        var r = c0 != a && c0 != b ? c0 : c1 != a && c1 != b ? c1 : c2;
        V3 pa = Point(a), pb = Point(b), pr = Point(r);
        var e = (pb - pa).Unit();
        var w = 0.5 * (pa + pb) - pr;
        w = w.Perp(e).Perp(bN);
        var n = w.Unit();
        return (n, bN, V3.Cross(n, bN));
    }

    private void BuildChains()
    {
        var adj = new Dictionary<int, List<int>>();
        foreach (var (a, b) in FrontEdges)
        {
            if (!adj.TryGetValue(a, out var la)) adj[a] = la = new List<int>(2);
            la.Add(b);
            if (!adj.TryGetValue(b, out var lb)) adj[b] = lb = new List<int>(2);
            lb.Add(a);
        }

        var used = new HashSet<long>();
        var starts = adj.Where(kv => kv.Value.Count != 2).Select(kv => kv.Key).OrderBy(k => k).ToList();
        void WalkFrom(int start)
        {
            // open chains start at nodes whose front degree is not 2 (ends and branching nodes)
            foreach (var first in adj[start].OrderBy(k => k).ToArray())
            {
                if (used.Contains(EdgeKey(start, first))) continue;
                var nodes = new List<int> { start };
                int prev = start, cur = first;
                used.Add(EdgeKey(prev, cur));
                while (true)
                {
                    nodes.Add(cur);
                    var nb = adj[cur];
                    if (nb.Count != 2) break;
                    var next = nb[0] == prev ? nb[1] : nb[0];
                    if (used.Contains(EdgeKey(cur, next))) break;
                    used.Add(EdgeKey(cur, next));
                    prev = cur;
                    cur = next;
                }

                AddChain(nodes, false);
            }
        }

        foreach (var s0 in starts) WalkFrom(s0);
        // remaining edges belong to closed loops of degree-2 nodes
        foreach (var (a, b) in FrontEdges)
        {
            if (used.Contains(EdgeKey(a, b))) continue;
            var nodes = new List<int> { a };
            int prev = a, cur = b;
            used.Add(EdgeKey(a, b));
            while (cur != a)
            {
                nodes.Add(cur);
                var nb = adj[cur];
                var next = nb[0] == prev ? nb[1] : nb[0];
                if (used.Contains(EdgeKey(cur, next))) break;
                used.Add(EdgeKey(cur, next));
                prev = cur;
                cur = next;
            }

            AddChain(nodes, cur == a);
        }

        Chains.Sort((x, y) => x.Sheet != y.Sheet ? x.Sheet.CompareTo(y.Sheet) : x.Nodes[0].CompareTo(y.Nodes[0]));
    }

    private void AddChain(List<int> nodes, bool closed)
    {
        var m = nodes.Count;
        if (m < 2) return;
        var nseg = closed ? m : m - 1;
        var eN = new V3[nseg];
        var eB = new V3[nseg];
        var eT = new V3[nseg];
        var pairIdx = FrontTriangle(nodes[0], nodes[1 % m]);
        for (var j = 0; j < nseg; j++)
        {
            var (a, b) = (nodes[j], nodes[(j + 1) % m]);
            (eN[j], eB[j], eT[j]) = EdgeFrame(a, b);
        }

        // traversal direction: along T of the edges; reverse the chain when the majority disagrees
        var agree = 0.0;
        for (var j = 0; j < nseg; j++)
        {
            var d = Point(nodes[(j + 1) % m]) - Point(nodes[j]);
            agree += V3.Dot(d, eT[j]);
        }

        if (agree < 0)
        {
            nodes.Reverse();
            Array.Reverse(eN);
            Array.Reverse(eB);
            Array.Reverse(eT);
            if (closed)
            {
                // after reversing a closed loop the segment j joins nodes[j] and nodes[j+1] again only after a rotation
                var r = new V3[nseg];
                Array.Copy(eN, 1, r, 0, nseg - 1); r[nseg - 1] = eN[0]; eN = r;
                r = new V3[nseg];
                Array.Copy(eB, 1, r, 0, nseg - 1); r[nseg - 1] = eB[0]; eB = r;
                r = new V3[nseg];
                Array.Copy(eT, 1, r, 0, nseg - 1); r[nseg - 1] = eT[0]; eT = r;
            }
        }

        var P = new V3[m];
        var N = new V3[m];
        var B = new V3[m];
        var T = new V3[m];
        var S = new double[m];
        for (var j = 0; j < m; j++)
        {
            P[j] = Point(nodes[j]);
            // node frame: average of the adjacent segments
            var nAcc = V3.Zero;
            var bAcc = V3.Zero;
            var cnt = 0;
            if (j > 0 || closed)
            {
                var k = closed ? (j - 1 + nseg) % nseg : j - 1;
                nAcc += eN[k];
                bAcc += eB[k];
                cnt++;
            }

            if (j < nseg)
            {
                nAcc += eN[j];
                bAcc += eB[j];
                cnt++;
            }

            B[j] = bAcc.Unit();
            N[j] = nAcc.Perp(B[j]).Unit();
            T[j] = V3.Cross(N[j], B[j]);
            if (j > 0) S[j] = S[j - 1] + V3.Dist(P[j], P[j - 1]);
        }

        var length = S[m - 1] + (closed ? V3.Dist(P[m - 1], P[0]) : 0);
        var chain = new FrontChain
        {
            Nodes = nodes.ToArray(), P = P, N = N, B = B, T = T, S = S, Closed = closed, Length = length,
            Sheet = SheetOfPair[pairIdx], Bits = SheetBits[SheetOfPair[pairIdx]]
        };
        if (!closed)
        {
            chain.Start = EndKindOf(nodes[0]);
            chain.End = EndKindOf(nodes[^1]);
        }

        Chains.Add(chain);
    }

    private readonly HashSet<int> _junctionNodes = new();

    private EndKind EndKindOf(int c)
    {
        if (OuterNodes.Contains(c)) return EndKind.Mouth;
        if (_junctionNodes.Count == 0 && JunctionEdges.Count > 0)
            foreach (var (a, b) in JunctionEdges)
            {
                _junctionNodes.Add(a);
                _junctionNodes.Add(b);
            }

        return _junctionNodes.Contains(c) ? EndKind.Junction : EndKind.Free;
    }

    /// <summary>Mesh nodes (one per canonical label) of a chain.</summary>
    public int[] MeshNodes(FrontChain c) => c.Nodes.Select(k => Rep[k]).ToArray();

    /// <summary>All mesh nodes whose canonical label lies on a sheet's triangles (every copy).</summary>
    public HashSet<int> SheetNodes(int sheet)
    {
        var canon = new HashSet<int>();
        for (var i = 0; i < PairCount; i++)
            if (SheetOfPair[i] == sheet)
            {
                canon.Add(Tri[i].Item1);
                canon.Add(Tri[i].Item2);
                canon.Add(Tri[i].Item3);
            }

        var nodes = new HashSet<int>();
        var s = Snap!;
        for (var v = 0; v < s.NodeCount; v++)
            if (canon.Contains(_cm.Canon[v]))
                nodes.Add(v);
        return nodes;
    }

    /// <summary>Oriented triangles of one sheet (or all sheets with -1), optionally within a distance of a point set.</summary>
    public TriangleSoup Soup(int sheet, Func<V3, bool>? keep = null)
    {
        var soup = new TriangleSoup();
        for (var i = 0; i < PairCount; i++)
        {
            if (sheet >= 0 && SheetOfPair[i] != sheet) continue;
            var (a, b, c) = Tri[i];
            V3 pa = Point(a), pb = Point(b), pc = Point(c);
            if (keep != null && !keep((pa + pb + pc) / 3)) continue;
            soup.Add(pa, pb, pc, i);
        }

        soup.Build();
        return soup;
    }

    /// <summary>Soup of the outer boundary facets (oriented outward).</summary>
    public TriangleSoup OuterSoup(Snapshot s)
    {
        var soup = new TriangleSoup();
        var tri = new int[3];
        for (var i = 0; i < s.BoundaryCount; i++)
        {
            var id = s.BoundaryFacet(i);
            if (PairOfFacet.ContainsKey(id)) continue;
            int e = id / 4, fl = id % 4;
            var row = s.Row(e);
            if (row.Length != 4) continue;
            var k = 0;
            for (var p = 0; p < 4; p++)
                if (p != fl)
                    tri[k++] = row[p];
            V3 a = _cm.P(tri[0]), b = _cm.P(tri[1]), c = _cm.P(tri[2]);
            var d = _cm.P(row[fl]);
            soup.AddOriented(a, b, c, (a + b + c) / 3 - d, id);
        }

        soup.Build();
        return soup;
    }

    /// <summary>Area of a sheet's triangles.</summary>
    public double SheetArea(int sheet)
    {
        var area = 0.0;
        for (var i = 0; i < PairCount; i++)
            if (sheet < 0 || SheetOfPair[i] == sheet)
            {
                var (a, b, c) = Tri[i];
                area += 0.5 * Geo.AreaVector(Point(a), Point(b), Point(c)).Norm;
            }

        return area;
    }

    public string Summary()
    {
        var sb = new System.Text.StringBuilder();
        sb.Append($"crack facets {PairCount} sheets {SheetCount} front edges {FrontEdges.Count} mouth {MouthEdges.Count} junction {JunctionEdges.Count} chains {Chains.Count}");
        foreach (var c in Chains)
            sb.Append($"\n  chain sheet {c.Sheet} bits {c.Bits} nodes {c.Count} {(c.Closed ? "closed" : $"open {c.Start}-{c.End}")} length {c.Length:G5}");
        if (OrientationConflicts + SideConflicts + UnseededSheets > 0)
            sb.Append($"\n  conflicts orientation {OrientationConflicts} side {SideConflicts} unseeded {UnseededSheets}");
        return sb.ToString();
    }

    /// <summary>Events between a previous scan and this one, as text.</summary>
    public List<string> Events(CrackCount prev)
    {
        var ev = new List<string>();
        // sheet mapping through canonical triangles
        var prevIndex = new Dictionary<(int, int, int), int>();
        for (var i = 0; i < prev.PairCount; i++) prevIndex[Sorted(prev.Tri[i].Item1, prev.Tri[i].Item2, prev.Tri[i].Item3)] = prev.SheetOfPair[i];
        var curToPrev = new Dictionary<int, HashSet<int>>();
        var prevToCur = new Dictionary<int, HashSet<int>>();
        for (var i = 0; i < PairCount; i++)
        {
            var key = Sorted(Tri[i].Item1, Tri[i].Item2, Tri[i].Item3);
            if (!prevIndex.TryGetValue(key, out var ps)) continue;
            if (!curToPrev.TryGetValue(SheetOfPair[i], out var a)) curToPrev[SheetOfPair[i]] = a = new HashSet<int>();
            a.Add(ps);
            if (!prevToCur.TryGetValue(ps, out var b)) prevToCur[ps] = b = new HashSet<int>();
            b.Add(SheetOfPair[i]);
        }

        foreach (var kv in curToPrev)
            if (kv.Value.Count > 1)
                ev.Add($"coalescence: sheets {string.Join(",", kv.Value.OrderBy(x => x))} merged into sheet {kv.Key}");
        foreach (var kv in prevToCur)
            if (kv.Value.Count > 1)
                ev.Add($"split: sheet {kv.Key} became sheets {string.Join(",", kv.Value.OrderBy(x => x))}");
        if (SheetCount > prev.SheetCount + (curToPrev.Count(kv => kv.Value.Count > 1) == 0 ? 0 : 0))
        {
            var newSheets = Enumerable.Range(0, SheetCount).Where(sid => !curToPrev.ContainsKey(sid)).ToList();
            if (newSheets.Count > 0) ev.Add($"new sheets {string.Join(",", newSheets)}");
        }

        var newJunction = JunctionEdges.Count - prev.JunctionEdges.Count;
        if (prev.JunctionEdges.Count == 0 && JunctionEdges.Count > 0) ev.Add($"junction formed: {JunctionEdges.Count} junction edges");
        else if (newJunction > 0) ev.Add($"junction grew by {newJunction} edges");
        var prevMouth = prev.MouthEdges.Count;
        if (prevMouth == 0 && MouthEdges.Count > 0) ev.Add($"breakthrough: {MouthEdges.Count} mouth edges");
        else if (MouthEdges.Count > prevMouth) ev.Add($"mouth grew by {MouthEdges.Count - prevMouth} edges");
        int prevClosed = prev.Chains.Count(c => c.Closed), curClosed = Chains.Count(c => c.Closed);
        if (curClosed < prevClosed) ev.Add($"closed fronts {prevClosed} -> {curClosed}");
        if (Chains.Count < prev.Chains.Count) ev.Add($"fronts {prev.Chains.Count} -> {Chains.Count}");
        if (Chains.Count > prev.Chains.Count) ev.Add($"fronts {prev.Chains.Count} -> {Chains.Count}");
        return ev;
    }
}

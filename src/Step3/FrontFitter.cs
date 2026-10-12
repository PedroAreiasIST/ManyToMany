// FrontFitter.cs — fits the mesh to a front: nodes near the front wall move onto it, edges crossing it are cut.
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
// The wall of a front is the ruled surface through the new front curve along the crack
// normal (a strip of half-width a few element sizes). As in the planar examples of Step2,
// a node within a fraction of the element size of the wall moves onto it (along the wall
// normal) when no neighbour has moved and no tetrahedron of its star changes orientation
// or loses more than 90 % of its volume; the remaining edges that cross the wall more than
// a tenth of their length from either endpoint are cut at the crossing. The crack increment
// that follows then ends exactly on the wall, and the recovered front is a sharp curve
// rather than a path along template faces. Nodes that lie on a crack are never moved.

using Numerical;
using Step2;

namespace Step3;

public static class FrontFitter
{
    public sealed class Result
    {
        public int Candidates, Moved, EdgesCut, Rejected;

        /// <summary>Nodes placed on the wall: moved nodes and the nodes created on the crossing edges.</summary>
        public HashSet<int> WallNodes = new();

        /// <summary>Candidates within the move distance of the wall that could not be moved (rejected or blocked).</summary>
        public HashSet<int> Unmoved = new();
    }

    /// <summary>
    ///     Wall soup of a front: quads between p_i - w b_i and p_i + w b_i for consecutive knots (b the surface normal at
    ///     the new front, w the half-width). Closed fronts wrap.
    /// </summary>
    public static TriangleSoup Wall(V3[] p, V3[] b, double w, bool closed)
    {
        var soup = new TriangleSoup();
        var n = p.Length;
        var nq = closed ? n : n - 1;
        for (var i = 0; i < nq; i++)
        {
            var i1 = (i + 1) % n;
            V3 a0 = p[i] - w * b[i], a1 = p[i] + w * b[i], c0 = p[i1] - w * b[i1], c1 = p[i1] + w * b[i1];
            soup.Add(a0, c0, c1, i);
            soup.Add(a0, c1, a1, i);
        }

        soup.Build();
        return soup;
    }

    /// <summary>
    ///     Fits the mesh to the wall: <paramref name="h"/> the element size, <paramref name="moveTol"/> the move distance
    ///     in element sizes (0.2 in Step2), <paramref name="band"/> the distance from the wall within which edges are
    ///     examined (in element sizes). Only nodes without crack bits move.
    /// </summary>
    public static Result Fit(CrackMesh cm, TriangleSoup wall, double h, double moveTol = 0.2, double band = 1.5)
    {
        var m = cm.M;
        var res = new Result();
        var nN = m.NodeCount;
        var phi = new double[nN];
        var near = new bool[nN];
        var w = wall.MeanEdge; // the wall extent along the normal is 2 w_half; points farther than that are outside the strip
        Parallel.For(0, nN, v =>
        {
            if (m.TransposePairs(v).Length == 0) return;
            var p = cm.P(v);
            var (tri, cp, d2, feat) = wall.Nearest(p, band * h + 1e-9);
            if (tri < 0) { phi[v] = double.NaN; return; }
            // the signed distance to the wall plane of the nearest triangle; points beyond the strip's rim count as far
            var nrm = wall.Normal(tri);
            var d = V3.Dot(p - cp, nrm);
            var inPlane = (p - cp).Perp(nrm).Norm; // distance along the strip from the closest point (nonzero only beyond the rim)
            if (inPlane > 0.05 * h) { phi[v] = double.NaN; return; }
            phi[v] = d;
            near[v] = true;
        });

        // ---- nodes already on the wall (to a tolerance) count as wall nodes and are not moved or cut around
        var moved = new bool[nN];
        var onTol = 0.02 * h;
        for (var v = 0; v < nN; v++)
            if (near[v] && Math.Abs(phi[v]) <= onTol)
            {
                phi[v] = 0;
                res.WallNodes.Add(v);
            }

        // ---- moves: ascending label, no moved neighbour, no inverted or collapsed tetrahedron
        Span<int> row4 = stackalloc int[4];
        for (var v = 0; v < nN; v++)
        {
            if (!near[v] || phi[v] == 0 || Math.Abs(phi[v]) > moveTol * h || cm.Mask[v] != 0) continue;
            res.Candidates++;
            var cols = m.PatternColumns(v);
            var cnts = m.PatternCounts(v);
            var blocked = false;
            for (var q = 0; q < cols.Length && !blocked; q++)
                if (cnts[q] > 0 && moved[cols[q]])
                    blocked = true;
            if (blocked) continue;
            var p = cm.P(v);
            var (tri, cp, _, _) = wall.Nearest(p);
            var target = p - V3.Dot(p - cp, wall.Normal(tri)) * wall.Normal(tri);
            // check the star
            var ok = true;
            var tr = m.TransposePairs(v);
            for (var q = 0; q < tr.Length && ok; q += 2)
            {
                var row = m.Row(tr[q]);
                if (row.Length != 4) continue;
                V3[] pts = [cm.P(row[0]), cm.P(row[1]), cm.P(row[2]), cm.P(row[3])];
                var before = Geo.Volume6(pts[0], pts[1], pts[2], pts[3]);
                pts[row.IndexOf(v)] = target;
                var after = Geo.Volume6(pts[0], pts[1], pts[2], pts[3]);
                if (Math.Sign(after) != Math.Sign(before) || Math.Abs(after) < 0.1 * Math.Abs(before)) ok = false;
            }

            if (!ok)
            {
                res.Rejected++;
                continue;
            }

            cm.X.Set(v, target.X);
            cm.Y.Set(v, target.Y);
            cm.Z.Set(v, target.Z);
            moved[v] = true;
            phi[v] = 0;
            res.Moved++;
            res.WallNodes.Add(v);
        }

        // ---- edges crossing the wall: both endpoints near, opposite signs, crossing parameter in [0.1, 0.9] (a smaller
        // window leaves slivers that spoil the interaction integral near the front; the front then deviates from the
        // wall by at most a tenth of an edge at a candidate that could not move)
        var unmoved = new bool[nN];
        for (var v = 0; v < nN; v++)
            if (near[v] && !moved[v] && phi[v] != 0 && Math.Abs(phi[v]) <= moveTol * h && cm.Mask[v] == 0)
            {
                unmoved[v] = true;
                res.Unmoved.Add(v);
            }
        var marked = new List<(int, int)>();
        var roots = new Dictionary<long, double>(EdgeKeyComparer.Instance);
        for (var v = 0; v < nN; v++)
        {
            if (!near[v] || phi[v] == 0) continue;
            var cols = m.PatternColumns(v);
            var cnts = m.PatternCounts(v);
            for (var q = 0; q < cols.Length; q++)
            {
                var u = cols[q];
                if (cnts[q] == 0 || u <= v || !near[u] || phi[u] == 0) continue;
                if (phi[v] * phi[u] >= 0) continue;
                var window = 0.1;
                // root along the edge of the signed distance to the wall (bisection on the plane of the nearest triangle)
                V3 pv = cm.P(v), pu = cm.P(u);
                var t = CrackInserter3.RootOf(pv, pu, x =>
                {
                    var (tri, cp, _, _) = wall.Nearest(x);
                    return V3.Dot(x - cp, wall.Normal(tri));
                });
                if (t < window || t > 1 - window) continue;
                marked.Add((v, u));
                roots[Cutter.Key(v, u)] = t;
            }
        }

        if (marked.Count > 0)
        {
            marked.Sort();
            cm.PendingRoots = roots;
            cm.PendingMask = 0;
            var rec = Cutter.Cut(m, marked, cm);
            foreach (var (_, _, mid) in rec.Edges) res.WallNodes.Add(mid);
            res.EdgesCut = marked.Count;
        }

        return res;
    }
}

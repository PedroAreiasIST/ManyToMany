// TriangleSoup.cs — oriented triangle sets with a uniform grid: nearest point, signed distance, ray hits.
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
// A soup holds oriented triangles (the normal points to the positive side) with an
// integer tag and an optional scalar per vertex (used for the advance coordinate of
// increment strips). Queries run through a uniform grid over the bounding box, built
// once after the triangles are added. The signed distance takes the sign from the
// normals of all triangles at the minimum distance (within a relative tolerance), so a
// point near a fold or a crease is classified by the average normal of the fold.

namespace Step3;

public sealed class TriangleSoup
{
    private readonly List<double> _p = new();   // 9 per triangle
    private readonly List<double> _n = new();   // 3 per triangle, unit
    private readonly List<int> _tag = new();
    private readonly List<double> _val = new(); // 3 per triangle (vertex scalars), optional
    private double[] _P = [], _N = [], _Val = [];
    private int[] _Tag = [];
    private int[] _cellStart = [], _cellTris = [];
    private int _nx, _ny, _nz;
    private double _x0, _y0, _z0, _h;
    private bool _built;

    public int Count => _built ? _P.Length / 9 : _p.Count / 9;
    public bool IsEmpty => Count == 0;
    public ReadOnlySpan<double> Points => _P;
    public ReadOnlySpan<double> Normals => _N;
    public ReadOnlySpan<int> Tags => _Tag;
    public double MeanEdge { get; private set; }

    /// <summary>Adds the triangle (a, b, c) with the normal (b - a) x (c - a); a degenerate triangle is skipped.</summary>
    public int Add(V3 a, V3 b, V3 c, int tag = 0, double va = 0, double vb = 0, double vc = 0)
    {
        if (_built) throw new InvalidOperationException("soup already built");
        var n = V3.Cross(b - a, c - a);
        var l = n.Norm;
        if (l <= 1e-300) return -1;
        n /= l;
        _p.AddRange([a.X, a.Y, a.Z, b.X, b.Y, b.Z, c.X, c.Y, c.Z]);
        _n.AddRange([n.X, n.Y, n.Z]);
        _tag.Add(tag);
        _val.AddRange([va, vb, vc]);
        return _tag.Count - 1;
    }

    /// <summary>Adds the triangle with a prescribed orientation: the vertex order is flipped when needed.</summary>
    public int AddOriented(V3 a, V3 b, V3 c, V3 towards, int tag = 0, double va = 0, double vb = 0, double vc = 0)
    {
        var n = V3.Cross(b - a, c - a);
        return V3.Dot(n, towards) >= 0 ? Add(a, b, c, tag, va, vb, vc) : Add(a, c, b, tag, va, vc, vb);
    }

    public V3 Vertex(int t, int k) => new(_P[9 * t + 3 * k], _P[9 * t + 3 * k + 1], _P[9 * t + 3 * k + 2]);
    public V3 Normal(int t) => new(_N[3 * t], _N[3 * t + 1], _N[3 * t + 2]);
    public int Tag(int t) => _Tag[t];
    public double Value(int t, int k) => _Val[3 * t + k];

    /// <summary>Builds the grid; the cell size is <paramref name="cellFactor"/> times the mean edge length.</summary>
    public void Build(double cellFactor = 2.0)
    {
        _P = _p.ToArray();
        _N = _n.ToArray();
        _Tag = _tag.ToArray();
        _Val = _val.ToArray();
        _built = true;
        var nt = Count;
        if (nt == 0) { _nx = _ny = _nz = 1; _cellStart = new int[2]; _cellTris = []; _h = 1; return; }
        double x0 = double.MaxValue, y0 = x0, z0 = x0, x1 = -x0, y1 = x1, z1 = x1, esum = 0;
        for (var t = 0; t < nt; t++)
        {
            for (var k = 0; k < 3; k++)
            {
                var v = Vertex(t, k);
                x0 = Math.Min(x0, v.X); x1 = Math.Max(x1, v.X);
                y0 = Math.Min(y0, v.Y); y1 = Math.Max(y1, v.Y);
                z0 = Math.Min(z0, v.Z); z1 = Math.Max(z1, v.Z);
            }

            esum += V3.Dist(Vertex(t, 0), Vertex(t, 1)) + V3.Dist(Vertex(t, 1), Vertex(t, 2)) + V3.Dist(Vertex(t, 2), Vertex(t, 0));
        }

        MeanEdge = esum / (3 * nt);
        _h = Math.Max(cellFactor * MeanEdge, 1e-12);
        const int maxCells = 2_000_000;
        while ((long)Math.Ceiling((x1 - x0) / _h + 1) * (long)Math.Ceiling((y1 - y0) / _h + 1) * (long)Math.Ceiling((z1 - z0) / _h + 1) > maxCells)
            _h *= 1.5;
        _x0 = x0 - 0.5 * _h; _y0 = y0 - 0.5 * _h; _z0 = z0 - 0.5 * _h;
        _nx = (int)Math.Ceiling((x1 - _x0) / _h) + 1;
        _ny = (int)Math.Ceiling((y1 - _y0) / _h) + 1;
        _nz = (int)Math.Ceiling((z1 - _z0) / _h) + 1;
        var counts = new int[_nx * _ny * _nz + 1];
        var ranges = new int[6 * nt];
        for (var t = 0; t < nt; t++)
        {
            Range(t, ranges.AsSpan(6 * t, 6));
            for (var k = ranges[6 * t + 4]; k <= ranges[6 * t + 5]; k++)
            for (var j = ranges[6 * t + 2]; j <= ranges[6 * t + 3]; j++)
            for (var i = ranges[6 * t]; i <= ranges[6 * t + 1]; i++)
                counts[Cell(i, j, k) + 1]++;
        }

        for (var c = 0; c < counts.Length - 1; c++) counts[c + 1] += counts[c];
        _cellStart = counts;
        _cellTris = new int[counts[^1]];
        var fill = (int[])counts.Clone();
        for (var t = 0; t < nt; t++)
            for (var k = ranges[6 * t + 4]; k <= ranges[6 * t + 5]; k++)
            for (var j = ranges[6 * t + 2]; j <= ranges[6 * t + 3]; j++)
            for (var i = ranges[6 * t]; i <= ranges[6 * t + 1]; i++)
                _cellTris[fill[Cell(i, j, k)]++] = t;
    }

    private int Cell(int i, int j, int k) => i + _nx * (j + _ny * k);

    private void Range(int t, Span<int> r)
    {
        double x0 = double.MaxValue, y0 = x0, z0 = x0, x1 = -x0, y1 = x1, z1 = x1;
        for (var k = 0; k < 3; k++)
        {
            var v = Vertex(t, k);
            x0 = Math.Min(x0, v.X); x1 = Math.Max(x1, v.X);
            y0 = Math.Min(y0, v.Y); y1 = Math.Max(y1, v.Y);
            z0 = Math.Min(z0, v.Z); z1 = Math.Max(z1, v.Z);
        }

        r[0] = Ix(x0, _x0, _nx); r[1] = Ix(x1, _x0, _nx);
        r[2] = Ix(y0, _y0, _ny); r[3] = Ix(y1, _y0, _ny);
        r[4] = Ix(z0, _z0, _nz); r[5] = Ix(z1, _z0, _nz);
    }

    private int Ix(double x, double o, int n) => Math.Clamp((int)Math.Floor((x - o) / _h), 0, n - 1);

    /// <summary>Nearest triangle to p: index, closest point, squared distance and feature code (see Geo.ClosestOnTriangle).</summary>
    public (int Tri, V3 Point, double D2, int Feature) Nearest(V3 p, double maxDist = double.PositiveInfinity)
    {
        if (!_built) Build();
        if (Count == 0) return (-1, p, double.PositiveInfinity, -1);
        int ci = Ix(p.X, _x0, _nx), cj = Ix(p.Y, _y0, _ny), ck = Ix(p.Z, _z0, _nz);
        var best = (Tri: -1, Point: p, D2: maxDist < double.PositiveInfinity ? maxDist * maxDist : double.MaxValue, Feature: -1);
        var maxRing = Math.Max(_nx, Math.Max(_ny, _nz));
        for (var ring = 0; ring <= maxRing; ring++)
        {
            // the ring can only contain points closer than (ring - 1) h + distance of p to its own cell boundary; stop when
            // the best distance found is below the lower bound of this ring
            if (ring > 0)
            {
                var lower = (ring - 1) * _h; // conservative lower bound of the distance to any point in ring cells
                if (lower * lower > best.D2) break;
            }

            for (var k = ck - ring; k <= ck + ring; k++)
            {
                if (k < 0 || k >= _nz) continue;
                for (var j = cj - ring; j <= cj + ring; j++)
                {
                    if (j < 0 || j >= _ny) continue;
                    for (var i = ci - ring; i <= ci + ring; i++)
                    {
                        if (i < 0 || i >= _nx) continue;
                        if (Math.Max(Math.Abs(i - ci), Math.Max(Math.Abs(j - cj), Math.Abs(k - ck))) != ring) continue;
                        var c = Cell(i, j, k);
                        for (var q = _cellStart[c]; q < _cellStart[c + 1]; q++)
                        {
                            var t = _cellTris[q];
                            var cp = Geo.ClosestOnTriangle(p, Vertex(t, 0), Vertex(t, 1), Vertex(t, 2), out var feat);
                            var d2 = (cp - p).Norm2;
                            if (d2 < best.D2) best = (t, cp, d2, feat);
                        }
                    }
                }
            }
        }

        return best;
    }

    /// <summary>
    ///     Signed distance to the soup: positive on the side of the normals. The sign comes from the sum of the normals
    ///     of all triangles whose distance is within <paramref name="relTol"/> of the minimum (one triangle in the
    ///     interior of a face; several at a shared edge or vertex). The triangle index returned is the nearest one.
    /// </summary>
    public (double Phi, int Tri) SignedDistance(V3 p, double relTol = 1e-9)
    {
        var (t, cp, d2, _) = Nearest(p);
        if (t < 0) return (double.PositiveInfinity, -1);
        var d = Math.Sqrt(d2);
        var dir = p - cp;
        // gather the normals of all triangles at (almost) the same distance: the cells around cp
        var tol = Math.Max(relTol * Math.Max(d, MeanEdge), 1e-14);
        var nsum = V3.Zero;
        int ci = Ix(cp.X, _x0, _nx), cj = Ix(cp.Y, _y0, _ny), ck = Ix(cp.Z, _z0, _nz);
        var seen = 0;
        for (var k = Math.Max(0, ck - 1); k <= Math.Min(_nz - 1, ck + 1); k++)
        for (var j = Math.Max(0, cj - 1); j <= Math.Min(_ny - 1, cj + 1); j++)
        for (var i = Math.Max(0, ci - 1); i <= Math.Min(_nx - 1, ci + 1); i++)
        {
            var c = Cell(i, j, k);
            for (var q = _cellStart[c]; q < _cellStart[c + 1]; q++)
            {
                var u = _cellTris[q];
                var cq = Geo.ClosestOnTriangle(p, Vertex(u, 0), Vertex(u, 1), Vertex(u, 2), out _);
                if (Math.Abs((cq - p).Norm - d) <= tol)
                {
                    nsum += Normal(u);
                    seen++;
                }
            }
        }

        if (seen == 0) nsum = Normal(t);
        var s = V3.Dot(dir, nsum);
        if (d <= 1e-300) return (0, t);
        return (s >= 0 ? d : -d, t);
    }

    /// <summary>Interpolated vertex scalar at the closest point of triangle t to p.</summary>
    public double ValueAt(int t, V3 p)
    {
        var (l0, l1, l2) = Geo.Barycentric(p, Vertex(t, 0), Vertex(t, 1), Vertex(t, 2));
        l0 = Math.Clamp(l0, 0, 1); l1 = Math.Clamp(l1, 0, 1); l2 = Math.Clamp(l2, 0, 1);
        var s = l0 + l1 + l2;
        return (l0 * Value(t, 0) + l1 * Value(t, 1) + l2 * Value(t, 2)) / s;
    }

    /// <summary>First hit of the ray o + t d (unit d) with the soup for 0 &lt; t &lt;= tmax; triangle index -1 when none.</summary>
    public (double T, int Tri) FirstHit(V3 o, V3 d, double tmax)
    {
        if (!_built) Build();
        var best = (T: double.PositiveInfinity, Tri: -1);
        if (Count == 0) return best;
        // walk the cells along the ray (3D DDA)
        var p = o;
        var t0 = 0.0;
        // clip the ray to the grid box
        double bx1 = _x0 + _nx * _h, by1 = _y0 + _ny * _h, bz1 = _z0 + _nz * _h;
        double tEnter = 0, tExit = tmax;
        Clip(o.X, d.X, _x0, bx1, ref tEnter, ref tExit);
        Clip(o.Y, d.Y, _y0, by1, ref tEnter, ref tExit);
        Clip(o.Z, d.Z, _z0, bz1, ref tEnter, ref tExit);
        if (tEnter > tExit) return best;
        t0 = Math.Max(0, tEnter);
        p = o + t0 * d;
        int i = Ix(p.X, _x0, _nx), j = Ix(p.Y, _y0, _ny), k = Ix(p.Z, _z0, _nz);
        int si = d.X > 0 ? 1 : -1, sj = d.Y > 0 ? 1 : -1, sk = d.Z > 0 ? 1 : -1;
        double NextT(double x, double dx, double o0, int idx, int s) =>
            dx == 0 ? double.PositiveInfinity : (o0 + (idx + (s > 0 ? 1 : 0)) * _h - x) / dx;
        var tx = t0 + NextT(p.X, d.X, _x0, i, si);
        var ty = t0 + NextT(p.Y, d.Y, _y0, j, sj);
        var tz = t0 + NextT(p.Z, d.Z, _z0, k, sk);
        double dtx = d.X != 0 ? _h / Math.Abs(d.X) : double.PositiveInfinity;
        double dty = d.Y != 0 ? _h / Math.Abs(d.Y) : double.PositiveInfinity;
        double dtz = d.Z != 0 ? _h / Math.Abs(d.Z) : double.PositiveInfinity;
        var tcur = t0;
        var guard = 0;
        while (tcur <= tExit && guard++ < 10_000_000)
        {
            var c = Cell(i, j, k);
            for (var q = _cellStart[c]; q < _cellStart[c + 1]; q++)
            {
                var u = _cellTris[q];
                var th = Geo.RayTriangle(o, d, Vertex(u, 0), Vertex(u, 1), Vertex(u, 2));
                if (!double.IsNaN(th) && th > 1e-12 && th <= tmax && th < best.T) best = (th, u);
            }

            var tnext = Math.Min(tx, Math.Min(ty, tz));
            if (best.Tri >= 0 && best.T <= tnext) break;
            if (tnext == tx) { i += si; tx += dtx; }
            else if (tnext == ty) { j += sj; ty += dty; }
            else { k += sk; tz += dtz; }
            if (i < 0 || i >= _nx || j < 0 || j >= _ny || k < 0 || k >= _nz) break;
            tcur = tnext;
        }

        return best.Tri >= 0 ? best : (double.PositiveInfinity, -1);
    }

    private static void Clip(double o, double d, double lo, double hi, ref double tEnter, ref double tExit)
    {
        if (d == 0)
        {
            if (o < lo || o > hi) { tEnter = double.PositiveInfinity; tExit = double.NegativeInfinity; }
            return;
        }

        double t1 = (lo - o) / d, t2 = (hi - o) / d;
        if (t1 > t2) (t1, t2) = (t2, t1);
        tEnter = Math.Max(tEnter, t1);
        tExit = Math.Min(tExit, t2);
    }

    /// <summary>Distance from p to the soup (unsigned).</summary>
    public double Distance(V3 p) => Math.Sqrt(Nearest(p).D2);
}

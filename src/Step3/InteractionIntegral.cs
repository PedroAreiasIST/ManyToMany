// InteractionIntegral.cs — domain interaction integral along a front with free-surface and crack-face terms.
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
// The virtual extension is V = q N, q the product of a tube factor (1 up to the plateau
// gamma r_q, then linear to 0 at r_q) and a hat of half-width delta in arc length. The
// domain integral (Step2 Section 8.1) is evaluated with four points per tetrahedron in
// fixed blocks, so the result does not depend on the thread count. Two surface terms are
// subtracted where the tube meets a boundary facet: on free outer facets (the actual
// traction vanishes but the auxiliary one does not), and on crack facets when they are
// not aligned with the local crack plane (curved cracks), on both twins with their own
// displacements. Loaded, constrained and symmetry facets are excluded. The J-integral of
// the same domain, which needs no auxiliary field, is returned with its own free-surface
// term, and compared with the energy release rate of the factors by the driver.

using System.Collections.Concurrent;
using Step2;

namespace Step3;

public sealed class SifOptions
{
    /// <summary>Plateau fraction gamma of the tube factor (0 gives the linear factor of the penny-crack tables).</summary>
    public double Plateau = 0.5;

    public bool FreeSurfaceTerm = true;
    public bool CrackFaceTerm = true;
    public bool JSurfaceTerm = true;

    /// <summary>Hat half-width: max(HatFactor * knot spacing, HatRq * r_q).</summary>
    public double HatFactor = 1.0, HatRq = 1.0;

    public int Block = 64;

    /// <summary>
    ///     Surface-following auxiliary coordinates: r2 is the signed distance to the actual crack surface (extended ahead
    ///     of the front along N) and r1 the distance of its foot point from the front along the surface, so the branch
    ///     cut of the auxiliary field lies on the actual faces of a kinked or curved crack, as in the curvilinear
    ///     formulation of Gosz and Moran. The planar frame of the nearest front point is used when false.
    /// </summary>
    public bool SurfaceFollowing = true;

    /// <summary>Length of the surface extension ahead of the front, in tube radii.</summary>
    public double ExtensionRq = 1.5;
}

public sealed class SifResult
{
    public double KI, KII, KIII, J, HatLength, S;
    public double[] FreeTerm = new double[4], CrackTerm = new double[4];
    public int Elements;

    /// <summary>Energy release rate of the factors.</summary>
    public double G(Material m) => m.EnergyReleaseRate(KI, KII, KIII);
}

public sealed class InteractionIntegral
{
    private readonly Snapshot _s;
    private readonly double[] _X, _u;
    private readonly Material _m;
    private readonly FrontCurve3 _fr;
    private readonly double _rq;
    private readonly SifOptions _opt;

    // tube
    private readonly int[] _nodes, _nodeSeg;
    private readonly double[] _nodeS, _nodeQ;
    private readonly V3[] _nodeDir;
    private readonly int[] _elements, _elemHint;
    private readonly double[] _elemSMin, _elemSMax;
    private readonly List<(int Facet, FacetKind Kind)> _facets = new(); // boundary facets touching the tube
    private readonly bool[] _inTube;
    private readonly TriangleSoup? _surf; // crack surface near the front plus its extension ahead (surface-following map)
    private const int ExtTag = -7;

    public int TubeNodes => _nodes.Length;
    public int TubeElements => _elements.Length;
    public double Delta { get; }

    /// <summary>
    ///     Builds the tube of radius <paramref name="rq"/> around the front; <paramref name="facetKind"/> classifies the
    ///     boundary facet ids (e * 4 + f) that touch the tube.
    /// </summary>
    public InteractionIntegral(Snapshot s, double[] X, double[] u, Material m, FrontCurve3 fr, double rq, SifOptions opt,
        Func<int, FacetKind> facetKind, TriangleSoup? crackSurface = null)
    {
        _s = s;
        _X = X;
        _u = u;
        _m = m;
        _fr = fr;
        _rq = rq;
        _opt = opt;
        if (opt.SurfaceFollowing && crackSurface != null)
        {
            // the surface near the front (within 2 r_q) and the extension ahead along N
            _surf = new TriangleSoup();
            var keepR = 2.5 * rq;
            for (var t = 0; t < crackSurface.Count; t++)
            {
                V3 a = crackSurface.Vertex(t, 0), b = crackSurface.Vertex(t, 1), c = crackSurface.Vertex(t, 2);
                var cen = (a + b + c) / 3;
                if (fr.Nearest(cen).D2 > keepR * keepR) continue;
                _surf.Add(a, b, c, crackSurface.Tag(t));
            }

            var ext = opt.ExtensionRq * rq;
            var nseg = fr.SegmentCount;
            for (var j = 0; j < nseg; j++)
            {
                var j1 = (j + 1) % fr.M;
                V3 p0 = fr.P[j], p1 = fr.P[j1], q0 = p0 + ext * fr.N[j], q1 = p1 + ext * fr.N[j1];
                var towards = fr.B[j] + fr.B[j1];
                _surf.AddOriented(p0, p1, q1, towards, ExtTag);
                _surf.AddOriented(p0, q1, q0, towards, ExtTag);
            }

            _surf.Build();
        }
        Delta = Math.Max(opt.HatFactor * fr.Spacing, opt.HatRq * rq);
        var nN = s.NodeCount;
        _nodeSeg = new int[nN];
        _nodeS = new double[nN];
        _nodeQ = new double[nN];
        _nodeDir = new V3[nN];
        Array.Fill(_nodeSeg, -1);
        double x0 = double.MaxValue, y0 = x0, z0 = x0, x1 = -x0, y1 = x1, z1 = x1;
        foreach (var p in fr.P)
        {
            x0 = Math.Min(x0, p.X); x1 = Math.Max(x1, p.X);
            y0 = Math.Min(y0, p.Y); y1 = Math.Max(y1, p.Y);
            z0 = Math.Min(z0, p.Z); z1 = Math.Max(z1, p.Z);
        }

        x0 -= rq; y0 -= rq; z0 -= rq; x1 += rq; y1 += rq; z1 += rq;
        var inTube = new bool[nN];
        var gamma = opt.Plateau;
        Parallel.For(0, nN, v =>
        {
            if (s.TransposePairs(v).Length == 0) return;
            var p = V3.Of(X, v);
            if (p.X < x0 || p.X > x1 || p.Y < y0 || p.Y > y1 || p.Z < z0 || p.Z > z1) return;
            var (k, t, d2) = fr.Nearest(p);
            if (d2 >= rq * rq) return;
            var ft = fr.FootAt(p, k, t, d2);
            var d = ft.D;
            var tube = gamma <= 0 ? Math.Max(0, 1 - d / rq)
                : d <= gamma * rq ? 1.0 : Math.Max(0, 1 - (d - gamma * rq) / (rq - gamma * rq));
            if (tube <= 0) return;
            _nodeSeg[v] = k;
            _nodeS[v] = ft.S;
            _nodeQ[v] = tube;
            _nodeDir[v] = ft.N;
            inTube[v] = true;
        });
        _inTube = inTube;
        var nodes = new List<int>();
        for (var v = 0; v < nN; v++)
            if (inTube[v])
                nodes.Add(v);
        _nodes = nodes.ToArray();
        var els = new List<int>();
        var hint = new List<int>();
        var smin = new List<double>();
        var smax = new List<double>();
        for (var e = 0; e < s.ElementCount; e++)
        {
            var row = s.Row(e);
            if (row.Length != 4) continue;
            var any = false;
            var h = -1;
            double lo = double.MaxValue, hi = -1;
            foreach (var v in row)
                if (inTube[v])
                {
                    any = true;
                    if (h < 0) h = _nodeSeg[v];
                    lo = Math.Min(lo, _nodeS[v]);
                    hi = Math.Max(hi, _nodeS[v]);
                }

            if (!any) continue;
            els.Add(e);
            hint.Add(h);
            smin.Add(lo);
            smax.Add(hi);
        }

        _elements = els.ToArray();
        _elemHint = hint.ToArray();
        _elemSMin = smin.ToArray();
        _elemSMax = smax.ToArray();
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
            if (!(inTube[tri[0]] || inTube[tri[1]] || inTube[tri[2]])) continue;
            var kind = facetKind(id);
            if (kind == FacetKind.Prescribed || kind == FacetKind.Symmetry) continue;
            _facets.Add((id, kind));
        }
    }

    /// <summary>Local coordinates of x: planar frame of the nearest front point, or surface-following when a surface was given.</summary>
    public LocalFrame Map(V3 x, int hint)
    {
        var ft = _fr.Foot(x, hint);
        if (_surf == null) return new LocalFrame(ft.R1, ft.R2, ft.N, ft.B, ft.Tg, ft.K);
        var (tri, cp, _, _) = _surf.Nearest(x);
        if (tri < 0) return new LocalFrame(ft.R1, ft.R2, ft.N, ft.B, ft.Tg, ft.K);
        var nrm = _surf.Normal(tri);
        var r2 = V3.Dot(x - cp, nrm);
        var ahead = _surf.Tag(tri) == ExtTag;
        // foot of the surface point on the front: distance along the surface (chord) and the in-surface direction
        var fc = _fr.Foot(cp, ft.K);
        var d = (cp - fc.F).Norm;
        if (d < 1e-9)
        {
            // on the front: the front frame with the surface normal
            var n0 = ft.N.Perp(nrm).Unit();
            return new LocalFrame(ahead ? d : -d, r2, n0, nrm, V3.Cross(n0, nrm), fc.K);
        }

        var nDir = (ahead ? cp - fc.F : fc.F - cp).Perp(nrm).Unit();
        var bDir = nrm;
        if (V3.Dot(bDir, fc.B) < 0) { bDir = -bDir; r2 = -r2; }
        return new LocalFrame(ahead ? d : -d, r2, nDir, bDir, V3.Cross(nDir, bDir), fc.K);
    }

    /// <summary>Factors at every knot of the front.</summary>
    public SifResult[] EvaluateAll()
    {
        var r = new SifResult[_fr.M];
        for (var j = 0; j < _fr.M; j++) r[j] = Evaluate(j);
        return r;
    }

    /// <summary>Factors at knot j (hat of half-width Delta centered there).</summary>
    public SifResult Evaluate(int j) => EvaluateAt(_fr.S[j], Delta);

    /// <summary>Factors of a straight-through front: one hat equal to 1 along the whole curve.</summary>
    public SifResult EvaluateThrough() => EvaluateAt(0, double.PositiveInfinity);

    public SifResult EvaluateAt(double s0, double delta)
    {
        var nN = _s.NodeCount;
        var vn = new V3[nN];
        var any = new bool[nN];
        foreach (var v in _nodes)
        {
            var q = _nodeQ[v] * _fr.Hat(_nodeS[v], s0, delta);
            if (q <= 0) continue;
            any[v] = true;
            vn[v] = q * _nodeDir[v];
        }

        var sel = new List<int>();
        var wrap = _fr.Closed;
        for (var i = 0; i < _elements.Length; i++)
        {
            if (double.IsPositiveInfinity(delta)) { sel.Add(i); continue; }
            var lo = s0 - delta;
            var hi = s0 + delta;
            var ok = _elemSMax[i] >= lo && _elemSMin[i] <= hi;
            if (!ok && wrap)
                ok = _elemSMax[i] >= lo + _fr.Length && _elemSMin[i] <= hi + _fr.Length ||
                     _elemSMax[i] >= lo - _fr.Length && _elemSMin[i] <= hi - _fr.Length;
            if (ok) sel.Add(i);
        }

        var block = _opt.Block;
        var nb = (sel.Count + block - 1) / block;
        var partial = new double[nb][];
        var counted = 0;
        Parallel.For(0, nb, b =>
        {
            var part = new double[4];
            var cnt = 0;
            for (int k = b * block, end = Math.Min(sel.Count, k + block); k < end; k++)
                if (Term(sel[k], part, vn, any)) cnt++;
            partial[b] = part;
            Interlocked.Add(ref counted, cnt);
        });
        var sum = new double[4];
        foreach (var part in partial)
            for (var k = 0; k < 4; k++)
                sum[k] += part[k];

        var res = new SifResult { S = s0, Elements = counted };
        if (_opt.FreeSurfaceTerm || _opt.CrackFaceTerm) SurfaceTerms(vn, any, res);
        for (var k = 0; k < 3; k++) sum[k] -= res.FreeTerm[k] + res.CrackTerm[k];
        if (_opt.JSurfaceTerm) sum[3] -= res.FreeTerm[3] + res.CrackTerm[3];
        var l = _fr.HatLength(s0, delta);
        res.HatLength = l;
        var ep = _m.EPrime;
        res.KI = ep / 2 * sum[0] / l;
        res.KII = ep / 2 * sum[1] / l;
        res.KIII = _m.Mu * sum[2] / l;
        res.J = sum[3] / l;
        return res;
    }

    private bool Term(int idx, double[] part, V3[] vn, bool[] any)
    {
        int e = _elements[idx], hint = _elemHint[idx];
        var row = _s.Row(e);
        var touched = false;
        foreach (var v in row) touched |= any[v];
        if (!touched) return false;
        Span<double> g = stackalloc double[12];
        Span<double> G = stackalloc double[9];
        Span<double> S = stackalloc double[9];
        Span<double> GV = stackalloc double[9];
        Span<double> H = stackalloc double[9];
        Span<double> dv = stackalloc double[3];
        Span<double> dH = stackalloc double[27];
        var vol = P1Elasticity.Gradients(_X, row, g);
        var cen = Geo.Centroid(_X, row);
        var side = Map(cen, hint).R2 >= 0 ? 1 : -1;
        G.Clear();
        GV.Clear();
        for (var a = 0; a < 4; a++)
        {
            var vv = any[row[a]] ? vn[row[a]] : V3.Zero;
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
            {
                G[3 * i + j] += _u[3 * row[a] + i] * g[3 * a + j];
                GV[3 * i + j] += vv[i] * g[3 * a + j];
            }
        }

        P1Elasticity.Stress(G, _m, S);
        double W = 0, divV = GV[0] + GV[4] + GV[8];
        for (var i = 0; i < 3; i++)
        for (var j = 0; j < 3; j++)
            W += 0.5 * S[3 * i + j] * 0.5 * (G[3 * i + j] + G[3 * j + i]);
        var jv = -W * divV;
        for (var i = 0; i < 3; i++)
        for (var j = 0; j < 3; j++)
        for (var k = 0; k < 3; k++)
            jv += S[3 * i + j] * G[3 * i + k] * GV[3 * k + j];
        part[3] += jv * vol;

        const double qa = 0.5854101966249685, qb = 0.1381966011250105;
        for (var qp = 0; qp < 4; qp++)
        {
            var xq = V3.Zero;
            var Vq = V3.Zero;
            for (var a = 0; a < 4; a++)
            {
                var wa = a == qp ? qa : qb;
                xq += wa * V3.Of(_X, row[a]);
                if (any[row[a]]) Vq += wa * vn[row[a]];
            }

            for (var mode = 0; mode < 3; mode++)
            {
                Asymptotic.GradientGlobal(mode, xq, Map, hint, side, _m, H);
                var htr = H[0] + H[4] + H[8];
                double acc = 0, wint = 0;
                for (var i = 0; i < 3; i++)
                for (var j = 0; j < 3; j++)
                {
                    var sa = _m.Mu * (H[3 * i + j] + H[3 * j + i]) + (i == j ? _m.Lambda * htr : 0);
                    wint += S[3 * i + j] * 0.5 * (H[3 * i + j] + H[3 * j + i]);
                    for (var k = 0; k < 3; k++)
                        acc += (S[3 * i + j] * H[3 * i + k] + sa * G[3 * i + k]) * GV[3 * k + j];
                }

                acc -= wint * divV;
                Asymptotic.Derivatives(mode, xq, Map, hint, side, _m, dH, dv);
                for (var k = 0; k < 3; k++)
                {
                    double ck = 0;
                    for (var i = 0; i < 3; i++)
                    {
                        ck += dv[i] * G[3 * i + k];
                        for (var j = 0; j < 3; j++)
                            ck += S[3 * i + j] * (dH[9 * j + 3 * i + k] - dH[9 * k + 3 * i + j]);
                    }

                    acc += ck * Vq[k];
                }

                part[mode] += acc * vol / 4;
            }
        }

        return true;
    }

    /// <summary>Surface terms on the free outer facets and on the crack facets touching the hat.</summary>
    private void SurfaceTerms(V3[] vn, bool[] any, SifResult res)
    {
        var dbg = Environment.GetEnvironmentVariable("STEP3_DEBUG") == "4" ? new System.Text.StringBuilder() : null;
        Span<double> g = stackalloc double[12];
        Span<double> G = stackalloc double[9];
        Span<double> S = stackalloc double[9];
        Span<double> H = stackalloc double[9];
        var tri = new int[3];
        var free = new double[4];
        var crack = new double[4];
        foreach (var (id, kind) in _facets)
        {
            if (kind == FacetKind.Free && !_opt.FreeSurfaceTerm) continue;
            if (kind == FacetKind.Crack && !_opt.CrackFaceTerm) continue;
            int e = id / 4, fl = id % 4;
            var row = _s.Row(e);
            var k0 = 0;
            for (var p = 0; p < 4; p++)
                if (p != fl)
                    tri[k0++] = row[p];
            if (!(any[tri[0]] || any[tri[1]] || any[tri[2]])) continue;
            V3 a = V3.Of(_X, tri[0]), b = V3.Of(_X, tri[1]), c = V3.Of(_X, tri[2]);
            var av = Geo.AreaVector(a, b, c);
            var area = 0.5 * av.Norm;
            if (area <= 0) continue;
            var mN = av / (2 * area);
            var d = V3.Of(_X, row[fl]);
            if (V3.Dot(mN, (a + b + c) / 3 - d) < 0) mN = -mN; // outward from the element
            P1Elasticity.Gradients(_X, row, g);
            P1Elasticity.Gradient(_u, row, g, G);
            P1Elasticity.Stress(G, _m, S);
            var cen = Geo.Centroid(_X, row);
            var hint = _nodeSeg[tri[0]] >= 0 ? _nodeSeg[tri[0]] : _nodeSeg[tri[1]] >= 0 ? _nodeSeg[tri[1]] : _nodeSeg[tri[2]];
            var side = Map(cen, hint).R2 >= 0 ? 1 : -1;
            double W = 0;
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
                W += 0.5 * S[3 * i + j] * 0.5 * (G[3 * i + j] + G[3 * j + i]);
            var acc = kind == FacetKind.Crack ? crack : free;
            for (var qp = 0; qp < 3; qp++)
            {
                var xq = V3.Zero;
                var V = V3.Zero;
                for (var q = 0; q < 3; q++)
                {
                    var wq = q == qp ? 2.0 / 3 : 1.0 / 6;
                    xq += wq * V3.Of(_X, tri[q]);
                    if (any[tri[q]]) V += wq * vn[tri[q]];
                }

                var mv = V3.Dot(mN, V);
                // a crack facet point at or ahead of the polyline front (inside the chord sagitta of a curved front) sees
                // the auxiliary stress of the region ahead of the tip, which is not traction-free: skip it
                if (kind == FacetKind.Crack && Map(xq, hint).R1 > -0.02 * _rq) continue;
                // J: the actual traction vanishes on free and crack facets, so only -W m . V remains
                acc[3] += -W * mv * area / 3;
                for (var mode = 0; mode < 3; mode++)
                {
                    Asymptotic.GradientGlobal(mode, xq, Map, hint, side, _m, H);
                    var htr = H[0] + H[4] + H[8];
                    double wint = 0, bterm = 0;
                    for (var i = 0; i < 3; i++)
                    {
                        double sam = 0;
                        for (var j = 0; j < 3; j++)
                        {
                            var sa = _m.Mu * (H[3 * i + j] + H[3 * j + i]) + (i == j ? _m.Lambda * htr : 0);
                            sam += sa * mN[j];
                            wint += S[3 * i + j] * 0.5 * (H[3 * i + j] + H[3 * j + i]);
                        }

                        for (var k = 0; k < 3; k++) bterm += sam * G[3 * i + k] * V[k];
                    }

                    acc[mode] += (bterm - wint * mv) * area / 3;
                    if (dbg != null && mode == 0 && kind == FacetKind.Crack && qp == 0)
                    {
                        double sam0 = 0, sam1 = 0, sam2 = 0;
                        for (var j = 0; j < 3; j++)
                        {
                            sam0 += (_m.Mu * (H[0 + j] + H[3 * j + 0]) + (j == 0 ? _m.Lambda * htr : 0)) * mN[j];
                            sam1 += (_m.Mu * (H[3 + j] + H[3 * j + 1]) + (j == 1 ? _m.Lambda * htr : 0)) * mN[j];
                            sam2 += (_m.Mu * (H[6 + j] + H[3 * j + 2]) + (j == 2 ? _m.Lambda * htr : 0)) * mN[j];
                        }

                        var ft = Map(xq, hint);
                        dbg.AppendLine($"    facet {id} side {side} m {mN} V {V} mv {mv:G3} r1 {ft.R1:G4} r2 {ft.R2:G3} tracAux ({sam0:G3},{sam1:G3},{sam2:G3}) bterm {bterm:G3} wint {wint:G3} N {ft.N} B {ft.B}");
                    }
                }
            }
        }

        if (dbg != null)
        {
            Console.Write(dbg.ToString().Split('\n').Take(12).Aggregate("", (a, b) => a + b + "\n"));
            Console.WriteLine($"    crack-face term {crack[0]:G4} free term {free[0]:G4} (mode I sums)");
        }

        res.FreeTerm = free;
        res.CrackTerm = crack;
    }
}

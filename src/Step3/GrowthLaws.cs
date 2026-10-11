// GrowthLaws.cs — direction criteria, Paris law, front advance kinematics and the geometry of an increment.
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
// Direction: the maximum tangential stress criterion of Erdogan and Sih (kink from N
// toward B), the sigma_1' criterion of Schoellmann et al. (kink and twist with K_III) and
// the empirical criterion of Richard et al. Advance: a Paris law relative to the largest
// equivalent factor of the step, with cycle counting when C is given. Kinematics: the
// factors are extrapolated to the ends of open fronts, kink and advance are smoothed along
// the front, the kink is limited by an angle and by a minimum radius of curvature in
// element sizes, the end knots of a front that meets the outer boundary stay on it, and a
// knot whose advance reaches the outer boundary, another crack or another front within a
// snap distance is carried half an element beyond it (breakthrough, junction, coalescence).
// The increment geometry is the Hermite strip between the old front (the mesh nodes of the
// chain) and the new one, extended beyond the new front for the sign of the nodes ahead,
// with the old triangles near the front for sign continuity behind; the region function is
// the signed advance beyond the new front, minus a margin, and is positive behind the old
// front, on the old triangles, and on the far side of a crack the increment must stop at.

namespace Step3;

public enum DirectionCriterion
{
    MaximumTangentialStress,
    Schoellmann,
    Richard,
    Planar
}

public static class Criteria
{
    /// <summary>Kink angle of the maximum tangential stress criterion (Erdogan and Sih 1963), from N toward B.</summary>
    public static double MtsKink(double kI, double kII)
    {
        if (Math.Abs(kII) < 1e-14 * Math.Max(1e-300, Math.Abs(kI))) return 0;
        return 2 * Math.Atan((kI - Math.Sqrt(kI * kI + 8 * kII * kII)) / (4 * kII));
    }

    /// <summary>Asymptotic stresses times sqrt(2 pi r) at angle phi: (s_rr, s_pp, t_rp, t_pz) with s_zz = nu (s_rr + s_pp).</summary>
    public static (double Srr, double Spp, double Trp, double Tpz) NearTip(double kI, double kII, double kIII, double phi)
    {
        double c = Math.Cos(phi / 2), s = Math.Sin(phi / 2);
        var srr = kI * c * (1 + s * s) + kII * s * (1 - 3 * s * s);
        var spp = kI * c * c * c - 3 * kII * s * c * c;
        var trp = kI * s * c * c + kII * c * (1 - 3 * s * s);
        var tpz = kIII * c;
        return (srr, spp, trp, tpz);
    }

    /// <summary>Schoellmann et al. (2002): kink maximizing sigma_1' and the twist angle there; equivalent factor K_v.</summary>
    public static (double Kink, double Twist, double Kv) Schoellmann(double kI, double kII, double kIII, double nu)
    {
        double Sigma1(double phi)
        {
            var (srr, spp, _, tpz) = NearTip(kI, kII, kIII, phi);
            var szz = nu * (srr + spp);
            return 0.5 * (spp + szz) + 0.5 * Math.Sqrt((spp - szz) * (spp - szz) + 4 * tpz * tpz);
        }

        // coarse scan then golden-section refinement
        double best = double.NegativeInfinity, bphi = 0;
        const int n = 720;
        for (var i = 0; i <= n; i++)
        {
            var phi = -Math.PI * 0.95 + 1.9 * Math.PI * i / n;
            var v = Sigma1(phi);
            if (v > best)
            {
                best = v;
                bphi = phi;
            }
        }

        double a = bphi - 1.9 * Math.PI / n, b = bphi + 1.9 * Math.PI / n;
        const double gr = 0.6180339887498949;
        double x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1 = Sigma1(x1), f2 = Sigma1(x2);
        for (var it = 0; it < 60; it++)
        {
            if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = Sigma1(x1); }
            else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = Sigma1(x2); }
        }

        var phi0 = 0.5 * (a + b);
        var (srr0, spp0, _, tpz0) = NearTip(kI, kII, kIII, phi0);
        var szz0 = nu * (srr0 + spp0);
        var twist = Math.Abs(spp0 - szz0) + Math.Abs(tpz0) < 1e-300 ? 0 : 0.5 * Math.Atan2(2 * tpz0, spp0 - szz0);
        var kv = Sigma1(phi0); // sigma_1' times sqrt(2 pi r) at the kink equals K_v
        return (phi0, twist, kv);
    }

    /// <summary>Richard, Fulland and Sander (2005): kink, twist and equivalent factor (alpha1 = 1.155, alpha2 = 1).</summary>
    public static (double Kink, double Twist, double Kv) Richard(double kI, double kII, double kIII)
    {
        var den = Math.Max(kI, 0) + Math.Abs(kII) + Math.Abs(kIII);
        double kink = 0, twist = 0;
        if (den > 0)
        {
            var r2 = Math.Abs(kII) / den;
            var r3 = Math.Abs(kIII) / den;
            kink = -Math.Sign(kII) * (140 * r2 - 70 * r2 * r2) * Math.PI / 180;
            twist = -Math.Sign(kIII) * (78 * r3 - 33 * r3 * r3) * Math.PI / 180;
        }

        const double a1 = 1.155, a2 = 1.0;
        var kv = 0.5 * kI + 0.5 * Math.Sqrt(kI * kI + 4 * a1 * a1 * kII * kII + 4 * a2 * a2 * kIII * kIII);
        return (kink, twist, kv);
    }

    /// <summary>Kink, twist and equivalent factor for a criterion.</summary>
    public static (double Kink, double Twist, double Kv) Direction(DirectionCriterion c, double kI, double kII, double kIII, double nu)
    {
        switch (c)
        {
            case DirectionCriterion.Planar:
                return (0, 0, Math.Sqrt(kI * kI + kII * kII + kIII * kIII / (1 - nu)));
            case DirectionCriterion.MaximumTangentialStress:
            {
                var kink = MtsKink(kI, kII);
                // K_v of the MTS criterion: sigma_pp sqrt(2 pi r) at the kink
                double cc = Math.Cos(kink / 2);
                var kv = cc * (kI * cc * cc - 1.5 * kII * Math.Sin(kink));
                return (kink, 0, Math.Max(kv, Math.Sqrt(kI * kI + kII * kII)) == kv ? kv : Math.Sqrt(kI * kI + kII * kII));
            }
            case DirectionCriterion.Schoellmann:
                return Schoellmann(kI, kII, kIII, nu);
            default:
                return Richard(kI, kII, kIII);
        }
    }
}

public sealed class ParisLaw
{
    /// <summary>da/dN = C (Delta K)^m; C = 0 disables cycle counting (relative advances only).</summary>
    public double C, M = 3;

    /// <summary>Advance of a knot relative to the largest equivalent factor: da_max (K / K_max)^m.</summary>
    public double Advance(double daMax, double k, double kMax) => kMax <= 0 ? 0 : daMax * Math.Pow(Math.Max(k, 0) / kMax, M);

    /// <summary>Cycles for the advance da_max at the largest factor.</summary>
    public double Cycles(double daMax, double kMax) => C > 0 && kMax > 0 ? daMax / (C * Math.Pow(kMax, M)) : 0;
}

/// <summary>Factors and derived quantities at the knots of a resampled front.</summary>
public sealed class FrontFactors
{
    public required FrontChain Chain;
    public required FrontCurve3 Curve;
    public required SifResult[] K;        // per knot (interaction integral)
    public double[]? DctKI, DctKII, DctKIII;
    public required double[] Kink, Twist, Kv, G;
    public double Rq, H;
    public int Bit;

    public int M => Curve.M;
    public double KvMax => Kv.Length == 0 ? 0 : Kv.Max();
}

public sealed class AdvanceOptions
{
    public DirectionCriterion Direction = DirectionCriterion.MaximumTangentialStress;
    public double KinkMax = 45 * Math.PI / 180;
    public double CurvatureElements = 12;   // minimum radius of curvature of the path in local element sizes
    public int SmoothKink = 3, SmoothAdvance = 6, SmoothFront = 3;
    public int StripRows = 4, ExtensionRows = 2;
    public double ExtensionElements = 3;    // extension of the surface beyond the new front, in element sizes
    public double MarginElements = 0.25;    // region margin beyond the new front
    public double BehindElements = 0.25;    // old triangles kept in the surface soup within this many r_q of the front
    public double SnapElements = 1.0;       // a knot this close to a target is carried through it
    public double OvershootElements = 0.5;
    public bool Through;                    // straight-through front: one kink and advance for the whole chain
}

/// <summary>A surface the advance must not cross without snapping (outer boundary or another crack sheet).</summary>
public sealed class AdvanceTarget
{
    public required TriangleSoup Soup;
    public required string Name;

    /// <summary>When set, the far side of this surface (sign opposite to the front) is excluded from the increment.</summary>
    public bool ExcludeFarSide;
}

/// <summary>Geometry of one increment: the level sets of the cutter and the new front positions.</summary>
public sealed class IncrementGeometry
{
    public required TriangleSoup Surface;      // oriented: old triangles near the front (value -big), strip and extension (value = advance beyond the new front)
    public required V3[] NewDense;             // new positions of the chain's mesh nodes
    public required double[] DenseAdvance, DenseKink;
    public required bool[] DenseSnapped;
    public double Margin, Behind;
    public readonly List<(TriangleSoup Soup, int Side)> Excluded = new();
    public int Snapped;

    public double PhiS(V3 x) => Surface.SignedDistance(x).Phi;

    public double PhiR(V3 x)
    {
        var (tri, cp, _, _) = Surface.Nearest(x);
        if (tri < 0 || Surface.Tag(tri) < 0) return 1; // behind the old front
        var a = Surface.ValueAt(tri, cp);
        var phi = a - Margin;
        foreach (var (soup, side) in Excluded)
        {
            // the far side exists only where x projects onto the interior of a triangle of the other sheet
            var (t, cq, _, feat) = soup.Nearest(x);
            if (t < 0 || feat != 0) continue;
            var d = V3.Dot(x - cq, soup.Normal(t));
            if (side * d < -Behind) return 1;
        }

        return phi;
    }
}

public static class FrontAdvance
{
    private static void Smooth(double[] v, int passes, bool closed, bool[]? keep = null)
    {
        var n = v.Length;
        if (n < 3) return;
        var w = new double[n];
        for (var p = 0; p < passes; p++)
        {
            for (var j = 0; j < n; j++)
            {
                if (keep != null && keep[j]) { w[j] = v[j]; continue; }
                var jm = j == 0 ? (closed ? n - 1 : 0) : j - 1;
                var jp = j == n - 1 ? (closed ? 0 : n - 1) : j + 1;
                w[j] = 0.25 * v[jm] + 0.5 * v[j] + 0.25 * v[jp];
            }

            Array.Copy(w, v, n);
        }
    }

    private static void SmoothPositions(V3[] p, int passes, bool closed, bool[]? keep = null)
    {
        var n = p.Length;
        if (n < 3) return;
        var w = new V3[n];
        for (var q = 0; q < passes; q++)
        {
            for (var j = 0; j < n; j++)
            {
                if (!closed && (j == 0 || j == n - 1) || keep != null && keep[j]) { w[j] = p[j]; continue; }
                var jm = j == 0 ? n - 1 : j - 1;
                var jp = j == n - 1 ? 0 : j + 1;
                w[j] = 0.25 * p[jm] + 0.5 * p[j] + 0.25 * p[jp];
            }

            Array.Copy(w, p, n);
        }
    }

    /// <summary>Quadratic extrapolation of a knot quantity to the ends of an open front from the interior knots within range.</summary>
    public static void ExtrapolateEnds(double[] v, double[] s, double skip, double range, bool closed)
    {
        if (closed || v.Length < 4) return;
        void End(bool first)
        {
            var n = v.Length;
            var sEnd = first ? s[0] : s[n - 1];
            var xs = new List<double>();
            var ys = new List<double>();
            for (var j = 0; j < n; j++)
            {
                var d = Math.Abs(s[j] - sEnd);
                if (d >= skip && d <= range)
                {
                    xs.Add(d);
                    ys.Add(v[j]);
                }
            }

            if (xs.Count < 2) return;
            // least squares a + b d (+ c d^2 with 4 or more points)
            var deg = xs.Count >= 4 ? 2 : 1;
            var val = Fit(xs, ys, deg);
            var inner = new List<int>();
            for (var j = 0; j < n; j++)
                if (Math.Abs(s[j] - sEnd) < skip)
                    inner.Add(j);
            foreach (var j in inner)
            {
                var d = Math.Abs(s[j] - sEnd);
                v[j] = val[0] + val[1] * d + (deg == 2 ? val[2] * d * d : 0);
            }
        }

        End(true);
        End(false);
    }

    private static double[] Fit(List<double> x, List<double> y, int deg)
    {
        var n = deg + 1;
        var a = new double[n, n];
        var b = new double[n];
        for (var k = 0; k < x.Count; k++)
        {
            var pw = new double[2 * n];
            pw[0] = 1;
            for (var i = 1; i < 2 * n; i++) pw[i] = pw[i - 1] * x[k];
            for (var i = 0; i < n; i++)
            {
                b[i] += pw[i] * y[k];
                for (var j = 0; j < n; j++) a[i, j] += pw[i + j];
            }
        }

        // Gaussian elimination
        for (var i = 0; i < n; i++)
        {
            var piv = i;
            for (var r = i + 1; r < n; r++)
                if (Math.Abs(a[r, i]) > Math.Abs(a[piv, i]))
                    piv = r;
            if (piv != i)
            {
                for (var c = 0; c < n; c++) (a[i, c], a[piv, c]) = (a[piv, c], a[i, c]);
                (b[i], b[piv]) = (b[piv], b[i]);
            }

            if (Math.Abs(a[i, i]) < 1e-300) return new double[n];
            for (var r = i + 1; r < n; r++)
            {
                var f = a[r, i] / a[i, i];
                for (var c = i; c < n; c++) a[r, c] -= f * a[i, c];
                b[r] -= f * b[i];
            }
        }

        var sol = new double[n];
        for (var i = n - 1; i >= 0; i--)
        {
            var acc = b[i];
            for (var c = i + 1; c < n; c++) acc -= a[i, c] * sol[c];
            sol[i] = acc / a[i, i];
        }

        return sol;
    }

    /// <summary>Kink, twist and equivalent factors at the knots, with end extrapolation of the factors on open fronts.</summary>
    public static void Evaluate(FrontFactors f, AdvanceOptions opt, double nu, double hatDelta)
    {
        var m = f.M;
        var kI = f.K.Select(k => k.KI).ToArray();
        var kII = f.K.Select(k => k.KII).ToArray();
        var kIII = f.K.Select(k => k.KIII).ToArray();
        if (!f.Curve.Closed)
        {
            // knots within one hat of a free end are unreliable (vertex singularity): extrapolate from the interior
            ExtrapolateEnds(kI, f.Curve.S, hatDelta, 3 * hatDelta, false);
            ExtrapolateEnds(kII, f.Curve.S, hatDelta, 3 * hatDelta, false);
            ExtrapolateEnds(kIII, f.Curve.S, hatDelta, 3 * hatDelta, false);
        }

        for (var j = 0; j < m; j++)
        {
            var (kink, twist, kv) = Criteria.Direction(opt.Direction, kI[j], kII[j], kIII[j], nu);
            f.Kink[j] = kink;
            f.Twist[j] = twist;
            f.Kv[j] = kv;
            f.K[j].KI = kI[j];
            f.K[j].KII = kII[j];
            f.K[j].KIII = kIII[j];
        }
    }

    /// <summary>
    ///     New positions of the dense chain and the increment geometry. <paramref name="advanceOf"/> gives the advance of a
    ///     knot from its equivalent factor (the Paris law relative to the step's largest factor).
    /// </summary>
    public static IncrementGeometry Build(FrontFactors f, Func<double, double> advanceOf, AdvanceOptions opt, double h,
        IReadOnlyList<AdvanceTarget> targets, IReadOnlyList<FrontChain> otherFronts, TriangleSoup oldTriangles, double otherAdvance = 0)
    {
        var c = f.Curve;
        var m = c.M;
        var kink = (double[])f.Kink.Clone();
        var da = new double[m];
        for (var j = 0; j < m; j++) da[j] = advanceOf(f.Kv[j]);
        if (opt.Through)
        {
            var kMean = kink.Average();
            var dMean = da.Average();
            for (var j = 0; j < m; j++)
            {
                kink[j] = kMean;
                da[j] = dMean;
            }
        }
        else
        {
            Smooth(kink, opt.SmoothKink, c.Closed);
            Smooth(da, opt.SmoothAdvance, c.Closed);
        }

        // limits: angle and curvature radius (kink <= da / (curvature radius))
        for (var j = 0; j < m; j++)
        {
            var lim = Math.Min(opt.KinkMax, da[j] / (opt.CurvatureElements * h));
            kink[j] = Math.Clamp(kink[j], -lim, lim);
        }

        // interpolate kink and advance to the dense chain nodes by arc length
        var ch = f.Chain;
        var n = ch.Count;
        var dKink = new double[n];
        var dDa = new double[n];
        for (var i = 0; i < n; i++)
        {
            var s = ch.S[i];
            var (k, t) = c.Locate(s);
            var k1 = (k + 1) % m;
            dKink[i] = (1 - t) * kink[k] + t * kink[k1];
            dDa[i] = (1 - t) * da[k] + t * da[k1];
        }

        // directions and snapping
        var dir = new V3[n];
        var snapped = new bool[n];
        var nsnap = 0;
        var snap = opt.SnapElements * h;
        var over = opt.OvershootElements * h;
        for (var i = 0; i < n; i++)
        {
            dir[i] = (Math.Cos(dKink[i]) * ch.N[i] + Math.Sin(dKink[i]) * ch.B[i]).Unit();
            var reach = dDa[i] + snap;
            var best = double.PositiveInfinity;
            foreach (var tg in targets)
            {
                var (t, tri) = tg.Soup.FirstHit(ch.P[i], dir[i], reach);
                if (tri >= 0 && t < best) best = t;
            }

            // another front, which advances too in this step: snap when the two advances close the gap to below snap
            var reach2 = reach + otherAdvance;
            foreach (var of in otherFronts)
            {
                if (ReferenceEquals(of, ch)) continue;
                var q1 = ch.P[i] + reach2 * dir[i];
                foreach (var (a, b) in of.Segments())
                {
                    var d = Geo.SegmentSegment(ch.P[i], q1, of.P[a], of.P[b], out var s, out _);
                    if (d < snap && s * reach2 < best) best = s * reach2;
                }
            }

            if (best < reach)
            {
                dDa[i] = best + over;
                snapped[i] = true;
                nsnap++;
            }
        }

        var newP = new V3[n];
        for (var i = 0; i < n; i++) newP[i] = ch.P[i] + dDa[i] * dir[i];
        // ends of a front on the outer boundary stay on the boundary surface (projection onto the nearest outer facet plane)
        var outer = targets.FirstOrDefault(t => t.Name == "outer");
        if (!ch.Closed && outer != null)
        {
            foreach (var i in new[] { 0, n - 1 })
            {
                var kindEnd = i == 0 ? ch.Start : ch.End;
                if (kindEnd != EndKind.Mouth) continue;
                var (tri, _, _, _) = outer.Soup.Nearest(ch.P[i]);
                if (tri < 0) continue;
                var nrm = outer.Soup.Normal(tri);
                var p0 = outer.Soup.Vertex(tri, 0);
                newP[i] -= V3.Dot(newP[i] - p0, nrm) * nrm;
                snapped[i] = false;
            }
        }

        if (!opt.Through) SmoothPositions(newP, opt.SmoothFront, ch.Closed, snapped);

        // ---- surface soup: old triangles near the front, Hermite strip, extension
        var soup = new TriangleSoup();
        var nOld = oldTriangles.Count;
        for (var t = 0; t < nOld; t++)
            soup.Add(oldTriangles.Vertex(t, 0), oldTriangles.Vertex(t, 1), oldTriangles.Vertex(t, 2), -1, -1e9, -1e9, -1e9);
        var rows = opt.StripRows;
        var ext = opt.ExtensionElements * h;
        var nrow = rows + 1 + opt.ExtensionRows;
        var pts = new V3[nrow, n];
        var val = new double[nrow, n];
        for (var i = 0; i < n; i++)
        {
            var p0 = ch.P[i];
            var p1 = newP[i];
            var m0 = dDa[i] * ch.N[i];
            var m1 = dDa[i] * dir[i];
            for (var r = 0; r <= rows; r++)
            {
                var tau = (double)r / rows;
                double h00 = 2 * tau * tau * tau - 3 * tau * tau + 1, h10 = tau * tau * tau - 2 * tau * tau + tau,
                    h01 = -2 * tau * tau * tau + 3 * tau * tau, h11 = tau * tau * tau - tau * tau;
                pts[r, i] = h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1;
                val[r, i] = (tau - 1) * dDa[i];
            }

            for (var r = 1; r <= opt.ExtensionRows; r++)
            {
                var e = snapped[i] ? 0.0 : ext * r / opt.ExtensionRows;
                pts[rows + r, i] = p1 + e * dir[i];
                val[rows + r, i] = e;
            }
        }

        var nq = ch.Closed ? n : n - 1;
        for (var r = 0; r + 1 < nrow; r++)
        for (var i = 0; i < nq; i++)
        {
            var i1 = (i + 1) % n;
            var towards = ch.B[i] + ch.B[i1];
            V3 a = pts[r, i], b = pts[r, i1], cc = pts[r + 1, i1], d = pts[r + 1, i];
            soup.AddOriented(a, b, cc, towards, 1, val[r, i], val[r, i1], val[r + 1, i1]);
            soup.AddOriented(a, cc, d, towards, 1, val[r, i], val[r + 1, i1], val[r + 1, i]);
        }

        soup.Build();
        var geo = new IncrementGeometry
        {
            Surface = soup, NewDense = newP, DenseAdvance = dDa, DenseKink = dKink, DenseSnapped = snapped,
            Margin = opt.MarginElements * h, Behind = opt.BehindElements * h, Snapped = nsnap
        };
        return geo;
    }
}

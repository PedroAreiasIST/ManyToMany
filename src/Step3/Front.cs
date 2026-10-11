// Front.cs — resampled front polylines with frames, foot points, arc-length hats.
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
// A front curve is a polyline of knots with the frame (N, B, T) at every knot: N the
// direction of growth in the crack surface, B the normal of the surface (positive side)
// and T = N x B along the curve. Closed curves wrap. The foot point of a point x gives
// the local coordinates r1 = (x - F) . N and r2 = (x - F) . B used by the asymptotic
// fields, with the frame re-orthogonalized against the segment tangent as in Step2.

namespace Step3;

public sealed class FrontCurve3
{
    public readonly int M;
    public readonly bool Closed;
    public readonly V3[] P, N, B, T;
    public readonly double[] S;
    public readonly double Length;
    public EndKind Start, End;

    public FrontCurve3(V3[] p, V3[] n, V3[] b, bool closed, EndKind start = EndKind.Closed, EndKind end = EndKind.Closed)
    {
        M = p.Length;
        if (M < 2) throw new ArgumentException("a front needs two knots");
        Closed = closed;
        P = p;
        B = new V3[M];
        N = new V3[M];
        T = new V3[M];
        for (var j = 0; j < M; j++)
        {
            B[j] = b[j].Unit();
            N[j] = n[j].Perp(B[j]).Unit();
            T[j] = V3.Cross(N[j], B[j]);
        }

        S = new double[M];
        for (var j = 1; j < M; j++) S[j] = S[j - 1] + V3.Dist(P[j], P[j - 1]);
        Length = S[M - 1] + (closed ? V3.Dist(P[M - 1], P[0]) : 0);
        Start = closed ? EndKind.Closed : start;
        End = closed ? EndKind.Closed : end;
    }

    public int SegmentCount => Closed ? M : M - 1;
    private int Next(int j) => (j + 1) % M;

    /// <summary>Mean knot spacing.</summary>
    public double Spacing => Length / SegmentCount;

    /// <summary>Equal-arc-length resampling of a chain into knots about <paramref name="spacing"/> apart (at least minKnots).</summary>
    public static FrontCurve3 Resample(FrontChain c, double spacing, int minKnots = 8)
    {
        var m = Math.Max(minKnots, (int)Math.Round(c.Length / spacing) + (c.Closed ? 0 : 1));
        var src = new FrontCurve3(c.P, c.N, c.B, c.Closed, c.Start, c.End);
        return Resample(src, m);
    }

    /// <summary>Equal-arc-length resampling into m knots (for a closed curve, m segments).</summary>
    public static FrontCurve3 Resample(FrontCurve3 c, int m)
    {
        var nseg = c.Closed ? m : m - 1;
        var p = new V3[m];
        var n = new V3[m];
        var b = new V3[m];
        for (var j = 0; j < m; j++)
        {
            var s = c.Length * j / nseg;
            if (!c.Closed && j == m - 1) s = c.Length;
            var (k, t) = c.Locate(s);
            var k1 = c.Next(k);
            p[j] = V3.Lerp(c.P[k], c.P[k1], t);
            n[j] = V3.Lerp(c.N[k], c.N[k1], t);
            b[j] = V3.Lerp(c.B[k], c.B[k1], t);
        }

        if (!c.Closed)
        {
            p[0] = c.P[0];
            p[m - 1] = c.P[c.M - 1];
        }

        return new FrontCurve3(p, n, b, c.Closed, c.Start, c.End);
    }

    /// <summary>Segment and parameter of the arc length s (wrapped for closed curves).</summary>
    public (int K, double T) Locate(double s)
    {
        if (Closed)
        {
            s %= Length;
            if (s < 0) s += Length;
        }
        else s = Math.Clamp(s, 0, Length);

        var nseg = SegmentCount;
        var k = 0;
        while (k < nseg - 1 && SegmentEnd(k) < s) k++;
        var l = SegmentLength(k);
        return (k, l > 0 ? Math.Clamp((s - S[k]) / l, 0, 1) : 0);
    }

    private double SegmentEnd(int k) => k + 1 < M ? S[k + 1] : Length;
    public double SegmentLength(int k) => V3.Dist(P[k], P[Next(k)]);

    public V3 Position(double s)
    {
        var (k, t) = Locate(s);
        return V3.Lerp(P[k], P[Next(k)], t);
    }

    /// <summary>Signed arc distance from s0 to s along the curve (shortest way around for closed curves).</summary>
    public double ArcDistance(double s, double s0)
    {
        var d = s - s0;
        if (!Closed) return d;
        d %= Length;
        if (d > Length / 2) d -= Length;
        if (d < -Length / 2) d += Length;
        return d;
    }

    public struct FootInfo
    {
        public int K;
        public double T, R1, R2, S, D;
        public V3 N, B, Tg, F;
    }

    private (int K, double T, double D2) Search(V3 x, int lo, int hi)
    {
        var bk = lo;
        double bt = 0, bd = double.MaxValue;
        for (var k = lo; k <= hi; k++)
        {
            var a = P[k];
            var d = P[Next(k)] - a;
            var dd = d.Norm2;
            var t = dd > 0 ? V3.Dot(x - a, d) / dd : 0;
            t = Math.Clamp(t, 0, 1);
            var d2 = (x - a - t * d).Norm2;
            if (d2 < bd)
            {
                bd = d2;
                bk = k;
                bt = t;
            }
        }

        return (bk, bt, bd);
    }

    /// <summary>Nearest segment (index, parameter, squared distance), searched around a hint segment when given.</summary>
    public (int K, double T, double D2) Nearest(V3 x, int hint = -1, int win = 4)
    {
        var last = SegmentCount - 1;
        if (hint < 0 || Closed) return Search(x, 0, last);
        int lo = Math.Max(0, hint - win), hi = Math.Min(last, hint + win);
        var r = Search(x, lo, hi);
        if ((r.K == lo && lo > 0) || (r.K == hi && hi < last)) r = Search(x, 0, last);
        return r;
    }

    public FootInfo Foot(V3 x, int hint = -1)
    {
        var (k, t, d2) = Nearest(x, hint);
        return FootAt(x, k, t, d2);
    }

    public FootInfo FootAt(V3 x, int k, double t, double d2)
    {
        var k1 = Next(k);
        var f = new FootInfo { K = k, T = t, D = Math.Sqrt(d2) };
        var seg = P[k1] - P[k];
        var sl = seg.Norm;
        f.S = S[k] + t * sl;
        f.F = P[k] + t * seg;
        var tg = sl > 0 ? seg / sl : T[k];
        var n = V3.Lerp(N[k], N[k1], t);
        var b = V3.Lerp(B[k], B[k1], t);
        n = n.Perp(tg).Unit();
        var bb = V3.Cross(tg, n);
        if (V3.Dot(bb, b) < 0) bb = -bb;
        f.N = n;
        f.B = bb;
        f.Tg = V3.Cross(n, bb);
        var e = x - f.F;
        f.R1 = V3.Dot(e, n);
        f.R2 = V3.Dot(e, bb);
        return f;
    }

    /// <summary>Integral along the curve of the hat of half-width delta centered at s0 (clipped to open curves, wrapped on closed ones).</summary>
    public double HatLength(double s0, double delta)
    {
        if (double.IsPositiveInfinity(delta)) return Length;
        if (Closed)
        {
            var d = Math.Min(delta, Length / 2);
            return d; // integral of a unit hat of half-width d: d (two triangles of area d/2)
        }

        double lo = Math.Max(0, s0 - delta), hi = Math.Min(Length, s0 + delta);
        double Int(double a, double b) // integral of 1 - |s - s0|/delta from a to b (a <= b, within one side)
        {
            double Fa(double s) => s - Math.Sign(s - s0) * (s - s0) * (s - s0) / (2 * delta);
            return Fa(b) - Fa(a);
        }

        var left = lo < s0 ? Int(lo, Math.Min(s0, hi)) : 0;
        var right = hi > s0 ? Int(Math.Max(s0, lo), hi) : 0;
        return left + right;
    }

    /// <summary>Hat value at arc length s for the hat centered at s0.</summary>
    public double Hat(double s, double s0, double delta)
    {
        if (double.IsPositiveInfinity(delta)) return 1;
        var ds = Math.Abs(ArcDistance(s, s0));
        return ds >= delta ? 0 : 1 - ds / delta;
    }

    /// <summary>Copy with new positions (same frames); frames are re-orthonormalized against the new tangents.</summary>
    public FrontCurve3 WithPositions(V3[] p, V3[]? n = null, V3[]? b = null) =>
        new(p, n ?? N, b ?? B, Closed, Start, End);
}

// Vec.cs — small three-vector helpers for the geometric code of Step3.
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

namespace Step3;

/// <summary>Immutable three-vector with the usual operators; used by the geometric code, not by the assembly loops.</summary>
public readonly record struct V3(double X, double Y, double Z)
{
    public static readonly V3 Zero = new(0, 0, 0);
    public static V3 operator +(V3 a, V3 b) => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
    public static V3 operator -(V3 a, V3 b) => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
    public static V3 operator -(V3 a) => new(-a.X, -a.Y, -a.Z);
    public static V3 operator *(double s, V3 a) => new(s * a.X, s * a.Y, s * a.Z);
    public static V3 operator *(V3 a, double s) => new(s * a.X, s * a.Y, s * a.Z);
    public static V3 operator /(V3 a, double s) => new(a.X / s, a.Y / s, a.Z / s);
    public static double Dot(V3 a, V3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;
    public static V3 Cross(V3 a, V3 b) => new(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
    public double Norm => Math.Sqrt(X * X + Y * Y + Z * Z);
    public double Norm2 => X * X + Y * Y + Z * Z;
    public V3 Unit() { var n = Norm; return n > 0 ? this / n : this; }
    public static double Dist(V3 a, V3 b) => (a - b).Norm;
    public double this[int c] => c == 0 ? X : c == 1 ? Y : Z;
    public static V3 Of(ReadOnlySpan<double> x, int v) => new(x[3 * v], x[3 * v + 1], x[3 * v + 2]);
    public void To(Span<double> x, int v) { x[3 * v] = X; x[3 * v + 1] = Y; x[3 * v + 2] = Z; }
    public static V3 Lerp(V3 a, V3 b, double t) => new(a.X + t * (b.X - a.X), a.Y + t * (b.Y - a.Y), a.Z + t * (b.Z - a.Z));

    /// <summary>Component of this vector orthogonal to the unit vector <paramref name="n"/>.</summary>
    public V3 Perp(V3 n) => this - Dot(this, n) * n;

    /// <summary>Rotation of this vector by <paramref name="angle"/> about the unit axis <paramref name="k"/> (Rodrigues).</summary>
    public V3 Rotate(V3 k, double angle)
    {
        double c = Math.Cos(angle), s = Math.Sin(angle);
        return c * this + s * Cross(k, this) + (1 - c) * Dot(k, this) * k;
    }

    public override string ToString() => $"({X:G6}, {Y:G6}, {Z:G6})";
}

/// <summary>Geometric helpers on flat coordinate arrays (3 doubles per node).</summary>
public static class Geo
{
    /// <summary>Centroid of the nodes of a row.</summary>
    public static V3 Centroid(ReadOnlySpan<double> x, ReadOnlySpan<int> row)
    {
        double cx = 0, cy = 0, cz = 0;
        foreach (var v in row)
        {
            cx += x[3 * v];
            cy += x[3 * v + 1];
            cz += x[3 * v + 2];
        }

        return new V3(cx / row.Length, cy / row.Length, cz / row.Length);
    }

    /// <summary>Twice the area vector of a triangle (a, b, c): (b - a) x (c - a).</summary>
    public static V3 AreaVector(V3 a, V3 b, V3 c) => V3.Cross(b - a, c - a);

    /// <summary>Signed volume times six of the tetrahedron (a, b, c, d).</summary>
    public static double Volume6(V3 a, V3 b, V3 c, V3 d) => V3.Dot(V3.Cross(b - a, c - a), d - a);

    /// <summary>
    ///     Closest point of the triangle (a, b, c) to p (Ericson, Real-Time Collision Detection, 5.1.5) and its
    ///     feature: 0 = interior, 1..3 = edge opposite to vertex (feature - 1), 4..6 = vertex (feature - 4).
    /// </summary>
    public static V3 ClosestOnTriangle(V3 p, V3 a, V3 b, V3 c, out int feature)
    {
        V3 ab = b - a, ac = c - a, ap = p - a;
        double d1 = V3.Dot(ab, ap), d2 = V3.Dot(ac, ap);
        if (d1 <= 0 && d2 <= 0) { feature = 4; return a; }
        V3 bp = p - b;
        double d3 = V3.Dot(ab, bp), d4 = V3.Dot(ac, bp);
        if (d3 >= 0 && d4 <= d3) { feature = 5; return b; }
        var vc = d1 * d4 - d3 * d2;
        if (vc <= 0 && d1 >= 0 && d3 <= 0)
        {
            feature = 3; // edge ab, opposite c
            var v = d1 / (d1 - d3);
            return a + v * ab;
        }

        V3 cp = p - c;
        double d5 = V3.Dot(ab, cp), d6 = V3.Dot(ac, cp);
        if (d6 >= 0 && d5 <= d6) { feature = 6; return c; }
        var vb = d5 * d2 - d1 * d6;
        if (vb <= 0 && d2 >= 0 && d6 <= 0)
        {
            feature = 2; // edge ac, opposite b
            var w = d2 / (d2 - d6);
            return a + w * ac;
        }

        var va = d3 * d6 - d5 * d4;
        if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
        {
            feature = 1; // edge bc, opposite a
            var w = (d4 - d3) / (d4 - d3 + (d5 - d6));
            return b + w * (c - b);
        }

        var denom = 1.0 / (va + vb + vc);
        double v2 = vb * denom, w2 = vc * denom;
        feature = 0;
        return a + v2 * ab + w2 * ac;
    }

    /// <summary>Barycentric coordinates of p with respect to the triangle (a, b, c) projected onto its plane.</summary>
    public static (double L0, double L1, double L2) Barycentric(V3 p, V3 a, V3 b, V3 c)
    {
        V3 v0 = b - a, v1 = c - a, v2 = p - a;
        double d00 = V3.Dot(v0, v0), d01 = V3.Dot(v0, v1), d11 = V3.Dot(v1, v1), d20 = V3.Dot(v2, v0), d21 = V3.Dot(v2, v1);
        var den = d00 * d11 - d01 * d01;
        var l1 = (d11 * d20 - d01 * d21) / den;
        var l2 = (d00 * d21 - d01 * d20) / den;
        return (1 - l1 - l2, l1, l2);
    }

    /// <summary>Ray–triangle intersection (Möller–Trumbore); returns the ray parameter or NaN.</summary>
    public static double RayTriangle(V3 o, V3 d, V3 a, V3 b, V3 c, double eps = 1e-12)
    {
        V3 e1 = b - a, e2 = c - a;
        var h = V3.Cross(d, e2);
        var det = V3.Dot(e1, h);
        if (Math.Abs(det) < eps) return double.NaN;
        var inv = 1 / det;
        var s = o - a;
        var u = inv * V3.Dot(s, h);
        if (u < -1e-9 || u > 1 + 1e-9) return double.NaN;
        var q = V3.Cross(s, e1);
        var v = inv * V3.Dot(d, q);
        if (v < -1e-9 || u + v > 1 + 1e-9) return double.NaN;
        return inv * V3.Dot(e2, q);
    }

    /// <summary>Distance between the segments (p0, p1) and (q0, q1) and the parameters of the closest points.</summary>
    public static double SegmentSegment(V3 p0, V3 p1, V3 q0, V3 q1, out double s, out double t)
    {
        V3 d1 = p1 - p0, d2 = q1 - q0, r = p0 - q0;
        double a = V3.Dot(d1, d1), e = V3.Dot(d2, d2), f = V3.Dot(d2, r);
        const double eps = 1e-30;
        if (a <= eps && e <= eps) { s = t = 0; return r.Norm; }
        if (a <= eps) { s = 0; t = Math.Clamp(f / e, 0, 1); }
        else
        {
            var c = V3.Dot(d1, r);
            if (e <= eps) { t = 0; s = Math.Clamp(-c / a, 0, 1); }
            else
            {
                var b = V3.Dot(d1, d2);
                var den = a * e - b * b;
                s = den != 0 ? Math.Clamp((b * f - c * e) / den, 0, 1) : 0;
                t = (b * s + f) / e;
                if (t < 0) { t = 0; s = Math.Clamp(-c / a, 0, 1); }
                else if (t > 1) { t = 1; s = Math.Clamp((b - c) / a, 0, 1); }
            }
        }

        return (p0 + s * d1 - (q0 + t * d2)).Norm;
    }
}

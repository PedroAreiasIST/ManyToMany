// FrontCurve.cs — a star-shaped planar crack front rho = R(psi), periodic cubic spline.
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
// The front lies in the crack plane (in-plane coordinates p, q; polar rho, psi). R is
// interpolated through M equally spaced angles by a periodic cubic spline, so the
// curve, its unit tangent and outward normal and its curvature are continuous. A point
// of the plane is projected on the curve by Newton's method on the parameter psi
// (foot point, signed normal distance, positive outside). A circle and an ellipse are
// the two constructors the examples use; growth replaces the radii.

namespace Step2.Bench;

public sealed class FrontCurve
{
    public readonly int M;
    public readonly double[] Radii;
    private readonly double[] _m; // second derivatives of R at the knots
    private readonly double _h;

    public FrontCurve(double[] radii)
    {
        M = radii.Length;
        Radii = (double[])radii.Clone();
        _h = 2 * Math.PI / M;
        _m = new double[M];
        // periodic spline: m_{j-1} + 4 m_j + m_{j+1} = 6 (r_{j+1} - 2 r_j + r_{j-1}) / h^2 (Gauss-Seidel, dominant)
        for (var it = 0; it < 500; it++)
        {
            double change = 0;
            for (var j = 0; j < M; j++)
            {
                int jm = (j + M - 1) % M, jp = (j + 1) % M;
                var rhs = 6 * (Radii[jp] - 2 * Radii[j] + Radii[jm]) / (_h * _h);
                var mj = (rhs - _m[jm] - _m[jp]) / 4;
                change = Math.Max(change, Math.Abs(mj - _m[j]));
                _m[j] = mj;
            }

            if (change < 1e-15 * (1 + Radii.Max() / (_h * _h))) break;
        }
    }

    public static FrontCurve Circle(double a, int m) => new(Enumerable.Repeat(a, m).ToArray());

    /// <summary>Ellipse with semi-axes a (along p) and b (along q) in polar form.</summary>
    public static FrontCurve Ellipse(double a, double b, int m) =>
        new(Enumerable.Range(0, m).Select(j =>
        {
            var psi = 2 * Math.PI * j / m;
            return a * b / Math.Sqrt(Math.Pow(b * Math.Cos(psi), 2) + Math.Pow(a * Math.Sin(psi), 2));
        }).ToArray());

    private (int J, double T) Locate(double psi)
    {
        var x = psi / _h;
        var j = (int)Math.Floor(x);
        var t = x - j;
        j = ((j % M) + M) % M;
        return (j, t);
    }

    /// <summary>R, dR/dpsi, d2R/dpsi2.</summary>
    public (double R, double D1, double D2) Eval(double psi)
    {
        var (j, t) = Locate(psi);
        var jp = (j + 1) % M;
        double r0 = Radii[j], r1 = Radii[jp], m0 = _m[j], m1 = _m[jp], h = _h, s = 1 - t;
        var r = s * r0 + t * r1 + h * h / 6 * ((s * s * s - s) * m0 + (t * t * t - t) * m1);
        var d1 = (r1 - r0) / h + h / 6 * (-(3 * s * s - 1) * m0 + (3 * t * t - 1) * m1);
        var d2 = s * m0 + t * m1;
        return (r, d1, d2);
    }

    public double R(double psi) => Eval(psi).R;

    /// <summary>Point, derivative and second derivative of the curve with respect to psi (in-plane).</summary>
    public (double X, double Y, double Dx, double Dy, double Ddx, double Ddy) Point(double psi)
    {
        var (r, d1, d2) = Eval(psi);
        double c = Math.Cos(psi), s = Math.Sin(psi);
        return (r * c, r * s, d1 * c - r * s, d1 * s + r * c, d2 * c - 2 * d1 * s - r * c, d2 * s + 2 * d1 * c - r * s);
    }

    /// <summary>Outward unit normal and |dC/dpsi| at psi.</summary>
    public (double Nx, double Ny, double Speed) Normal(double psi)
    {
        var (_, _, dx, dy, _, _) = Point(psi);
        var sp = Math.Sqrt(dx * dx + dy * dy);
        return (dy / sp, -dx / sp, sp);
    }

    /// <summary>Foot point parameter and signed normal distance (positive outside) of the in-plane point (p, q).</summary>
    public (double Psi, double D) Project(double p, double q)
    {
        var psi = Math.Atan2(q, p);
        for (var it = 0; it < 40; it++)
        {
            var (x, y, dx, dy, ddx, ddy) = Point(psi);
            double ex = x - p, ey = y - q;
            var f = ex * dx + ey * dy;
            var fp = dx * dx + dy * dy + ex * ddx + ey * ddy;
            if (fp <= 0) fp = dx * dx + dy * dy;
            var step = f / fp;
            if (step > 0.2) step = 0.2;
            else if (step < -0.2) step = -0.2;
            psi -= step;
            if (Math.Abs(step) < 1e-15) break;
        }

        var (cx, cy, _, _, _, _) = Point(psi);
        var (nx, ny, _) = Normal(psi);
        return (psi, (p - cx) * nx + (q - cy) * ny);
    }

    /// <summary>Integral of a hat of half-width delta centered at psi0 along the curve (arc length).</summary>
    public double HatLength(double psi0, double delta, int n = 400, double psiLo = double.NegativeInfinity,
        double psiHi = double.PositiveInfinity)
    {
        // the part of the hat inside [psiLo, psiHi] (a front that ends on a free surface); unclipped: the whole hat
        double xa = Math.Max(-delta, psiLo - psi0), xb = Math.Min(delta, psiHi - psi0);
        if (xb <= xa) return 0;
        double sum = 0;
        for (var i = 0; i <= n; i++)
        {
            var x = xa + (xb - xa) * i / n;
            var w = i == 0 || i == n ? 0.5 : 1.0;
            sum += w * (1 - Math.Abs(x) / delta) * Normal(psi0 + x).Speed;
        }

        return sum * (xb - xa) / n;
    }
}

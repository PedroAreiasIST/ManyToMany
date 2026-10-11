// Reference.cs — closed-form and handbook stress intensity factors of the classical examples.
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
// Plane strain: the centre crack in a strip (Tada, Paris and Irwin, 0.1 %), the edge crack in
// a strip with free ends (Tada, 0.5 %), the inclined crack in an infinite plate, two collinear
// cracks of equal length (Erdogan 1962; Westergaard function with the single-valuedness
// constant c^2 = b^2 E(k)/K(k), k^2 = 1 - a^2/b^2, cracks on a <= |x| <= b), the kink angle
// of Erdogan and Sih. Three dimensions: Sneddon's penny crack, Kassir and Sih's inclined penny
// crack, Irwin's elliptical crack, Newman and Raju's surface crack (1981).

namespace Step3.Bench;

public static class Reference
{
    /// <summary>Complete elliptic integrals K(k) and E(k) of the modulus k (AGM; k in [0, 1)).</summary>
    public static (double K, double E) Elliptic(double k)
    {
        if (k >= 1) return (double.PositiveInfinity, 1);
        double a = 1, b = Math.Sqrt(1 - k * k), c = k, sum = c * c, pow = 1;
        for (var i = 0; i < 60 && Math.Abs(c) > 1e-17; i++)
        {
            var an = 0.5 * (a + b);
            var bn = Math.Sqrt(a * b);
            c = 0.5 * (a - b);
            a = an;
            b = bn;
            pow *= 2;
            sum += pow * c * c;
        }

        var kk = Math.PI / (2 * a);
        var ee = kk * (1 - 0.5 * sum);
        return (kk, ee);
    }

    /// <summary>Centre crack of half-length a in a strip of width W under tension: F = K / (sigma sqrt(pi a)).</summary>
    public static double CentreCrackStrip(double a, double W)
    {
        var al = 2 * a / W;
        return (1 - 0.025 * al * al + 0.06 * al * al * al * al) * Math.Sqrt(1 / Math.Cos(Math.PI * al / 2));
    }

    /// <summary>Edge crack of length a in a strip of width W under tension (ends free to rotate): F = K / (sigma sqrt(pi a)).</summary>
    public static double EdgeCrackStrip(double a, double W)
    {
        var al = a / W;
        var t = Math.Tan(Math.PI * al / 2);
        var s = Math.Sin(Math.PI * al / 2);
        return Math.Sqrt(2 / (Math.PI * al) * t) * (0.752 + 2.02 * al + 0.37 * Math.Pow(1 - s, 3)) / Math.Cos(Math.PI * al / 2);
    }

    /// <summary>Two collinear cracks occupying a &lt;= |x| &lt;= b under remote tension sigma: (K at the inner tips, K at the outer tips).</summary>
    public static (double Inner, double Outer) CollinearCracks(double a, double b, double sigma)
    {
        var k = Math.Sqrt(1 - a * a / (b * b));
        var (kk, ee) = Elliptic(k);
        var c2 = b * b * ee / kk;
        var inner = sigma * Math.Sqrt(Math.PI * a) * (c2 / (a * a) - 1) / Math.Sqrt(b * b / (a * a) - 1);
        var outer = sigma * Math.Sqrt(Math.PI * b) * (1 - c2 / (b * b)) / Math.Sqrt(1 - a * a / (b * b));
        return (inner, outer);
    }

    /// <summary>Inclined crack of half-length a in an infinite plate, beta the angle between the crack and the load axis.</summary>
    public static (double KI, double KII) InclinedCrack(double a, double sigma, double beta)
    {
        var s = Math.Sin(beta);
        return (sigma * Math.Sqrt(Math.PI * a) * s * s, sigma * Math.Sqrt(Math.PI * a) * s * Math.Cos(beta));
    }

    /// <summary>Kink angle of the maximum tangential stress criterion (magnitude, radians).</summary>
    public static double MtsKink(double kI, double kII) => Math.Abs(Criteria.MtsKink(kI, kII));

    /// <summary>Sneddon: penny crack of radius a under tension sigma.</summary>
    public static double Sneddon(double a, double sigma) => 2 * sigma * Math.Sqrt(a / Math.PI);

    /// <summary>Kassir and Sih: inclined penny crack (beta between the normal and the load) at the angle theta from the shear direction.</summary>
    public static (double KI, double KII, double KIII) KassirSih(double a, double sigma, double beta, double theta, double nu)
    {
        double sn = sigma * Math.Cos(beta) * Math.Cos(beta), tau = sigma * Math.Sin(beta) * Math.Cos(beta);
        var kI = 2 * sn * Math.Sqrt(a / Math.PI);
        var kII = 4 * tau * Math.Sqrt(a / Math.PI) * Math.Cos(theta) / (2 - nu);
        var kIII = 4 * (1 - nu) * tau * Math.Sqrt(a / Math.PI) * Math.Sin(theta) / (2 - nu);
        return (kI, kII, kIII);
    }

    /// <summary>Irwin: elliptical crack with semi-axes a (along x) and b at the parametric angle phi.</summary>
    public static double Irwin(double a, double b, double sigma, double phi)
    {
        var k = Math.Sqrt(1 - b * b / (a * a));
        var (_, ee) = Elliptic(k);
        var s = Math.Sin(phi);
        var c = Math.Cos(phi);
        return sigma * Math.Sqrt(Math.PI * b) * Math.Pow(s * s + b * b / (a * a) * c * c, 0.25) / ee;
    }

    /// <summary>
    ///     Newman and Raju (1981): K_I of a semi-elliptical surface crack (depth a, half-length c) in a plate of thickness t and
    ///     half-width b under tension sigma, at the parametric angle phi (a/c &lt;= 1).
    /// </summary>
    public static double NewmanRaju(double a, double c, double t, double b, double sigma, double phi)
    {
        var ac = a / c;
        var at = a / t;
        var q = 1 + 1.464 * Math.Pow(ac, 1.65);
        var m1 = 1.13 - 0.09 * ac;
        var m2 = -0.54 + 0.89 / (0.2 + ac);
        var m3 = 0.5 - 1 / (0.65 + ac) + 14 * Math.Pow(1 - ac, 24);
        var g = 1 + (0.1 + 0.35 * at * at) * Math.Pow(1 - Math.Sin(phi), 2);
        var fphi = Math.Pow(ac * ac * Math.Cos(phi) * Math.Cos(phi) + Math.Sin(phi) * Math.Sin(phi), 0.25);
        var fw = Math.Sqrt(1 / Math.Cos(Math.PI * c / (2 * b) * Math.Sqrt(at)));
        var fs = (m1 + m2 * at * at + m3 * Math.Pow(at, 4)) * g * fphi * fw;
        return sigma * Math.Sqrt(Math.PI * a / q) * fs;
    }
}

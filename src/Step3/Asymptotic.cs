// Asymptotic.cs — near-tip auxiliary fields of unit K in the frame of the nearest front point.
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
// Plane-strain (modes I, II) and antiplane (mode III) displacements of unit K in the
// local frame (N, B, T) of the foot point, their gradient rotated to global axes without
// the derivative along the front (Gosz and Moran 2002), and the derivatives of that
// gradient by central differences, which feed the non-equilibrium and incompatibility
// terms of the interaction integral. Behind the front the angle is continued across the
// cut on the side of the element, as in Step2.

namespace Step3;

/// <summary>Local coordinates of a point relative to the crack: r1 along N (positive ahead), r2 along B, and the frame there.</summary>
public readonly record struct LocalFrame(double R1, double R2, V3 N, V3 B, V3 T, int K);

/// <summary>Mapping of a point to its local coordinates (planar frame of the nearest front point, or surface-following).</summary>
public delegate LocalFrame LocalMap(V3 x, int hint);

public static class Asymptotic
{
    /// <summary>Displacement (local frame) of mode 0, 1 or 2 at (y1, y2); <paramref name="side"/> continues the angle across the cut.</summary>
    public static void Displacement(int mode, double y1, double y2, int side, Material m, Span<double> uOut)
    {
        var r = Math.Sqrt(y1 * y1 + y2 * y2);
        var phi = Math.Atan2(y2, y1);
        if (side > 0 && phi < -Math.PI / 2) phi += 2 * Math.PI;
        else if (side < 0 && phi > Math.PI / 2) phi -= 2 * Math.PI;
        var kappa = 3 - 4 * m.Nu;
        var sr = Math.Sqrt(r / (2 * Math.PI));
        double c = Math.Cos(phi / 2), s = Math.Sin(phi / 2);
        var mu = m.Mu;
        uOut.Clear();
        switch (mode)
        {
            case 0:
                uOut[0] = sr / (2 * mu) * c * (kappa - 1 + 2 * s * s);
                uOut[1] = sr / (2 * mu) * s * (kappa + 1 - 2 * c * c);
                break;
            case 1:
                uOut[0] = sr / (2 * mu) * s * (kappa + 1 + 2 * c * c);
                uOut[1] = -sr / (2 * mu) * c * (kappa - 1 - 2 * s * s);
                break;
            default:
                uOut[2] = 2 / mu * sr * s;
                break;
        }
    }

    /// <summary>Planar mapping: the frame of the nearest front point.</summary>
    public static LocalFrame PlanarMap(FrontCurve3 fr, V3 x, int hint)
    {
        var ft = fr.Foot(x, hint);
        return new LocalFrame(ft.R1, ft.R2, ft.N, ft.B, ft.Tg, ft.K);
    }

    /// <summary>Auxiliary gradient H[3 i + k] (global axes) at x for the frame of the foot point on the front.</summary>
    public static void GradientGlobal(int mode, V3 x, FrontCurve3 fr, int hint, int side, Material m, Span<double> Hg) =>
        GradientGlobal(mode, x, (p, h) => PlanarMap(fr, p, h), hint, side, m, Hg);

    /// <summary>Auxiliary gradient H[3 i + k] (global axes) at x with the local coordinates and frame of <paramref name="map"/>.</summary>
    public static void GradientGlobal(int mode, V3 x, LocalMap map, int hint, int side, Material m, Span<double> Hg)
    {
        var ft = map(x, hint);
        double r1 = ft.R1, r2 = ft.R2;
        var d = 1e-6 * Math.Max(Math.Sqrt(r1 * r1 + r2 * r2), 1e-12);
        Span<double> up = stackalloc double[3], um = stackalloc double[3];
        Span<double> Hl = stackalloc double[9];
        Hl.Clear();
        for (var j = 0; j < 2; j++)
        {
            Displacement(mode, r1 + (j == 0 ? d : 0), r2 + (j == 1 ? d : 0), side, m, up);
            Displacement(mode, r1 - (j == 0 ? d : 0), r2 - (j == 1 ? d : 0), side, m, um);
            for (var i = 0; i < 3; i++) Hl[3 * i + j] = (up[i] - um[i]) / (2 * d);
        }

        // basis columns N, B, T: Hg = Q Hl Q^T restricted to the two in-plane derivatives
        Span<double> Q = stackalloc double[9];
        Q[0] = ft.N.X; Q[3] = ft.N.Y; Q[6] = ft.N.Z;
        Q[1] = ft.B.X; Q[4] = ft.B.Y; Q[7] = ft.B.Z;
        Q[2] = ft.T.X; Q[5] = ft.T.Y; Q[8] = ft.T.Z;
        for (var i = 0; i < 3; i++)
        for (var k = 0; k < 3; k++)
        {
            double acc = 0;
            for (var aa = 0; aa < 3; aa++)
            for (var bb = 0; bb < 2; bb++)
                acc += Q[3 * i + aa] * Hl[3 * aa + bb] * Q[3 * k + bb];
            Hg[3 * i + k] = acc;
        }
    }

    /// <summary>dH[9 j + 3 i + k] = d H_ik / d x_j and the divergence of the auxiliary stress, by central differences.</summary>
    public static void Derivatives(int mode, V3 x, FrontCurve3 fr, int hint, int side, Material m, Span<double> dH, Span<double> div) =>
        Derivatives(mode, x, (p, h) => PlanarMap(fr, p, h), hint, side, m, dH, div);

    public static void Derivatives(int mode, V3 x, LocalMap map, int hint, int side, Material m, Span<double> dH, Span<double> div)
    {
        var ft = map(x, hint);
        var d = 1e-3 * Math.Max(Math.Sqrt(ft.R1 * ft.R1 + ft.R2 * ft.R2), 1e-9);
        Span<double> Hp = stackalloc double[9], Hm = stackalloc double[9];
        div.Clear();
        for (var j = 0; j < 3; j++)
        {
            var e = new V3(j == 0 ? d : 0, j == 1 ? d : 0, j == 2 ? d : 0);
            GradientGlobal(mode, x + e, map, ft.K, side, m, Hp);
            GradientGlobal(mode, x - e, map, ft.K, side, m, Hm);
            for (var t = 0; t < 9; t++) dH[9 * j + t] = (Hp[t] - Hm[t]) / (2 * d);
            double trp = Hp[0] + Hp[4] + Hp[8], trm = Hm[0] + Hm[4] + Hm[8];
            for (var i = 0; i < 3; i++)
            {
                var sp = m.Mu * (Hp[3 * i + j] + Hp[3 * j + i]) + (i == j ? m.Lambda * trp : 0);
                var sm = m.Mu * (Hm[3 * i + j] + Hm[3 * j + i]) + (i == j ? m.Lambda * trm : 0);
                div[i] += (sp - sm) / (2 * d);
            }
        }
    }
}

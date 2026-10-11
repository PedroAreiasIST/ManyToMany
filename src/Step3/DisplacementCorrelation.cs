// DisplacementCorrelation.cs — factors from the displacement jump across the crack faces behind the front.
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
// At a knot, points at distances r behind the front along -N are projected onto the
// crack facets (through the oriented triangle soup of the sheet, which carries the pair
// index), the displacements of the two twins are interpolated and the jump is resolved
// in the frame (N, B, T): K_I = E'/8 sqrt(2 pi / r) [u]_B, K_II = E'/8 sqrt(2 pi / r) [u]_N,
// K_III = mu/4 sqrt(2 pi / r) [u]_T (plane strain). The value at the largest sample radius
// is reported: linear tetrahedra underestimate the jump near the front, so the estimate is a
// check of the sign and the order of magnitude of the interaction-integral factors, not of
// their value.

using Step2;

namespace Step3;

public static class DisplacementCorrelation
{
    /// <summary>Displacement interpolated on a boundary facet (e * 4 + f) at point p (barycentric on the facet).</summary>
    private static V3 OnFacet(Snapshot s, double[] X, double[] u, int facet, V3 p)
    {
        int e = facet / 4, fl = facet % 4;
        var row = s.Row(e);
        Span<int> tri = stackalloc int[3];
        var k = 0;
        for (var q = 0; q < 4; q++)
            if (q != fl)
                tri[k++] = row[q];
        var (l0, l1, l2) = Geo.Barycentric(p, V3.Of(X, tri[0]), V3.Of(X, tri[1]), V3.Of(X, tri[2]));
        l0 = Math.Clamp(l0, 0, 1); l1 = Math.Clamp(l1, 0, 1); l2 = Math.Clamp(l2, 0, 1);
        var sum = l0 + l1 + l2;
        return (l0 * V3.Of(u, tri[0]) + l1 * V3.Of(u, tri[1]) + l2 * V3.Of(u, tri[2])) / sum;
    }

    /// <summary>
    ///     Factors at knot j from the jumps at the radii <paramref name="radii"/> behind the front; the soup must be the
    ///     oriented soup of the sheet with pair indices as tags (CrackCount.Soup). Returns NaN when no facet is found.
    /// </summary>
    public static (double KI, double KII, double KIII, int Samples) At(FrontCurve3 fr, int j, CrackCount cc, TriangleSoup soup,
        Snapshot s, double[] X, double[] u, Material m, double[] radii)
    {
        var kI = new List<(double R, double K)>();
        var kII = new List<(double R, double K)>();
        var kIII = new List<(double R, double K)>();
        foreach (var r in radii)
        {
            var p = fr.P[j] - r * fr.N[j];
            var (tri, cp, d2, _) = soup.Nearest(p, 2 * r);
            if (tri < 0 || Math.Sqrt(d2) > 0.5 * r) continue;
            var pair = soup.Tag(tri);
            var up = OnFacet(s, X, u, cc.FacetPlus[pair], cp);
            var um = OnFacet(s, X, u, cc.FacetMinus[pair], cp);
            var jump = up - um;
            var f = Math.Sqrt(2 * Math.PI / r);
            // the frame at the foot point of cp
            var ft = fr.Foot(cp, j);
            kI.Add((r, m.EPrime / 8 * f * V3.Dot(jump, ft.B)));
            kII.Add((r, m.EPrime / 8 * f * V3.Dot(jump, ft.N)));
            kIII.Add((r, m.Mu / 4 * f * V3.Dot(jump, ft.Tg)));
        }

        if (kI.Count == 0) return (double.NaN, double.NaN, double.NaN, 0);
        // linear tetrahedra underestimate the jump within the first elements, so an extrapolation to r = 0 amplifies the
        // error; the value at the largest sample radius is reported (a check of sign and order of magnitude)
        return (kI[^1].K, kII[^1].K, kIII[^1].K, kI.Count);
    }

    /// <summary>Least-squares line through (r, K) evaluated at r = 0 (the mean for one sample).</summary>
    private static double Extrapolate(List<(double R, double K)> pts)
    {
        if (pts.Count == 1) return pts[0].K;
        double sr = 0, sk = 0, srr = 0, srk = 0;
        foreach (var (r, k) in pts)
        {
            sr += r;
            sk += k;
            srr += r * r;
            srk += r * k;
        }

        var n = pts.Count;
        var den = n * srr - sr * sr;
        if (Math.Abs(den) < 1e-300) return sk / n;
        var slope = (n * srk - sr * sk) / den;
        return (sk - slope * sr) / n;
    }
}

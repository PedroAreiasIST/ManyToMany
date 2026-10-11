// Elasticity.cs — linear tetrahedra: stiffness into the destination table, block-Jacobi CG, loads and fixities.
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
// The stiffness is assembled in 3 x 3 blocks per pattern slot, addressed by the
// destination table of the snapshot and gathered row by row on all threads: row v sums
// the blocks of the pairs of its transpose row in element order, so the operator equals
// that of a sequential element loop bit for bit. The solver is the block-Jacobi
// preconditioned conjugate gradient of the step2 examples on the flat copy of the
// pattern, with prescribed unknowns decoupled and a warm start.

using System.Collections.Concurrent;
using Step2;

namespace Step3;

/// <summary>Isotropic linear elastic material.</summary>
public sealed class Material
{
    public readonly double E, Nu, Lambda, Mu;

    public Material(double e, double nu)
    {
        E = e;
        Nu = nu;
        Lambda = e * nu / ((1 + nu) * (1 - 2 * nu));
        Mu = e / (2 * (1 + nu));
    }

    /// <summary>E / (1 - nu^2), the plane-strain modulus of the factor formulas.</summary>
    public double EPrime => E / (1 - Nu * Nu);

    /// <summary>Energy release rate of the three factors (plane strain for I and II).</summary>
    public double EnergyReleaseRate(double kI, double kII, double kIII) => (kI * kI + kII * kII) / EPrime + kIII * kIII / (2 * Mu);
}

public static class P1Elasticity
{
    /// <summary>Gradients of the barycentric coordinates (g[3a + c]) and the volume of a tetrahedron.</summary>
    public static double Gradients(ReadOnlySpan<double> X, ReadOnlySpan<int> row, Span<double> g)
    {
        Span<double> j = stackalloc double[9];
        for (var c = 0; c < 3; c++)
        {
            j[c * 3 + 0] = X[3 * row[1] + c] - X[3 * row[0] + c];
            j[c * 3 + 1] = X[3 * row[2] + c] - X[3 * row[0] + c];
            j[c * 3 + 2] = X[3 * row[3] + c] - X[3 * row[0] + c];
        }

        var det = j[0] * (j[4] * j[8] - j[5] * j[7]) - j[1] * (j[3] * j[8] - j[5] * j[6]) +
                  j[2] * (j[3] * j[7] - j[4] * j[6]);
        Span<double> inv = stackalloc double[9];
        inv[0] = (j[4] * j[8] - j[5] * j[7]) / det;
        inv[1] = (j[2] * j[7] - j[1] * j[8]) / det;
        inv[2] = (j[1] * j[5] - j[2] * j[4]) / det;
        inv[3] = (j[5] * j[6] - j[3] * j[8]) / det;
        inv[4] = (j[0] * j[8] - j[2] * j[6]) / det;
        inv[5] = (j[2] * j[3] - j[0] * j[5]) / det;
        inv[6] = (j[3] * j[7] - j[4] * j[6]) / det;
        inv[7] = (j[1] * j[6] - j[0] * j[7]) / det;
        inv[8] = (j[0] * j[4] - j[1] * j[3]) / det;
        for (var c = 0; c < 3; c++)
        {
            g[3 + c] = inv[0 * 3 + c];
            g[6 + c] = inv[1 * 3 + c];
            g[9 + c] = inv[2 * 3 + c];
            g[c] = -(g[3 + c] + g[6 + c] + g[9 + c]);
        }

        return Math.Abs(det) / 6.0;
    }

    /// <summary>Signed volume times six of a tetrahedron row.</summary>
    public static double SignedVolume6(ReadOnlySpan<double> X, ReadOnlySpan<int> row) =>
        Geo.Volume6(V3.Of(X, row[0]), V3.Of(X, row[1]), V3.Of(X, row[2]), V3.Of(X, row[3]));

    /// <summary>Displacement gradient G[3 i + j] = u_i,j of a P1 tetrahedron from its barycentric gradients.</summary>
    public static void Gradient(ReadOnlySpan<double> u, ReadOnlySpan<int> row, ReadOnlySpan<double> g, Span<double> G)
    {
        G.Clear();
        for (var a = 0; a < 4; a++)
        for (var i = 0; i < 3; i++)
        for (var j = 0; j < 3; j++)
            G[3 * i + j] += u[3 * row[a] + i] * g[3 * a + j];
    }

    /// <summary>Stress of a displacement gradient: mu (G + G^T) + lambda tr(G) I.</summary>
    public static void Stress(ReadOnlySpan<double> G, Material m, Span<double> S)
    {
        var tr = G[0] + G[4] + G[8];
        for (var i = 0; i < 3; i++)
        for (var j = 0; j < 3; j++)
            S[3 * i + j] = m.Mu * (G[3 * i + j] + G[3 * j + i]) + (i == j ? m.Lambda * tr : 0);
    }

    public static double VonMises(ReadOnlySpan<double> S)
    {
        double p = (S[0] + S[4] + S[8]) / 3;
        double a = S[0] - p, b = S[4] - p, c = S[8] - p;
        return Math.Sqrt(1.5 * (a * a + b * b + c * c + 2 * (S[1] * S[1] + S[2] * S[2] + S[5] * S[5])));
    }

    /// <summary>Stiffness blocks (9 per pattern slot) of the live tetrahedra of a snapshot.</summary>
    public static double[] Assemble(Snapshot s, FlatPattern fp, double[] X, Material mat)
    {
        if (9L * s.TotalSlots > Array.MaxLength)
            throw new InvalidOperationException($"{s.TotalSlots} pattern slots: 3 x 3 blocks exceed the array size");
        var vals = new double[9 * s.TotalSlots];
        var rowStart = fp.RowStart;
        double lambda = mat.Lambda, mu = mat.Mu;
        Parallel.ForEach(Partitioner.Create(0, s.NodeCount, 1024), r =>
        {
            Span<double> g = stackalloc double[12];
            for (var v = r.Item1; v < r.Item2; v++)
            {
                var tr = s.TransposePairs(v);
                var rs = rowStart[v];
                for (var q = 0; q < tr.Length; q += 2)
                {
                    var e = tr[q];
                    var row = s.Row(e);
                    if (row.Length != 4) continue;
                    var a = tr[q + 1];
                    var vol = Gradients(X, row, g);
                    var dest = s.Destinations(e);
                    for (var b = 0; b < 4; b++)
                    {
                        var at = 9 * (rs + dest[4 * a + b]);
                        var gg = g[3 * a] * g[3 * b] + g[3 * a + 1] * g[3 * b + 1] + g[3 * a + 2] * g[3 * b + 2];
                        for (var i = 0; i < 3; i++)
                        for (var jj = 0; jj < 3; jj++)
                            vals[at + 3 * i + jj] += vol * (lambda * g[3 * a + i] * g[3 * b + jj] +
                                                            mu * g[3 * a + jj] * g[3 * b + i] + (i == jj ? mu * gg : 0));
                    }
                }
            }
        });

        return vals;
    }

    /// <summary>y = K x on the free unknowns (rows of fixed unknowns are zero).</summary>
    public static void Apply(FlatPattern fp, double[] vals, double[] x, double[] y, bool[] fixedDof)
    {
        int[] rowStart = fp.RowStart, rowLen = fp.RowLen, col = fp.Col;
        Parallel.ForEach(Partitioner.Create(0, fp.NodeCount, 1024), r =>
        {
            for (var v = r.Item1; v < r.Item2; v++)
            {
                double y0 = 0, y1 = 0, y2 = 0;
                var rs = rowStart[v];
                if (rs >= 0)
                    for (int t = rs, end = rs + rowLen[v]; t < end; t++)
                    {
                        var at = 9 * t;
                        var c = 3 * col[t];
                        double x0 = x[c], x1 = x[c + 1], x2 = x[c + 2];
                        y0 += vals[at] * x0 + vals[at + 1] * x1 + vals[at + 2] * x2;
                        y1 += vals[at + 3] * x0 + vals[at + 4] * x1 + vals[at + 5] * x2;
                        y2 += vals[at + 6] * x0 + vals[at + 7] * x1 + vals[at + 8] * x2;
                    }

                y[3 * v] = fixedDof[3 * v] ? 0 : y0;
                y[3 * v + 1] = fixedDof[3 * v + 1] ? 0 : y1;
                y[3 * v + 2] = fixedDof[3 * v + 2] ? 0 : y2;
            }
        });
    }

    /// <summary>Result of a solve.</summary>
    public sealed class SolveResult
    {
        public required double[] U;
        public int Iterations;
        public double Relative;
        public bool Converged;
    }

    /// <summary>
    ///     Block-Jacobi preconditioned CG on the free unknowns (prescribed values zero; nodes without elements are
    ///     fixed). <paramref name="fixedDof"/> is updated with the auto-fixed unknowns. Reductions are sequential, so the
    ///     result does not depend on the thread count.
    /// </summary>
    public static SolveResult Pcg(FlatPattern fp, double[] vals, double[] f, bool[] fixedDof, double tol, int maxIt, double[]? x0)
    {
        var nN = fp.NodeCount;
        var nd = 3 * nN;
        var minv = new double[9 * nN];
        for (var v = 0; v < nN; v++)
        {
            var d = new double[9];
            var rs = fp.RowStart[v];
            if (rs >= 0)
                for (int t = rs, end = rs + fp.RowLen[v]; t < end; t++)
                    if (fp.Col[t] == v)
                        for (var q = 0; q < 9; q++)
                            d[q] = vals[9 * t + q];
            for (var i = 0; i < 3; i++)
            {
                if (!fixedDof[3 * v + i] && d[4 * i] != 0) continue;
                fixedDof[3 * v + i] = true;
                for (var q = 0; q < 3; q++) d[3 * i + q] = d[3 * q + i] = 0;
                d[4 * i] = 1;
            }

            Invert3(d, minv.AsSpan(9 * v, 9));
        }

        var u = new double[nd];
        var r = (double[])f.Clone();
        var r0 = 0.0;
        for (var i = 0; i < nd; i++)
        {
            if (fixedDof[i]) r[i] = 0;
            r0 += r[i] * r[i];
        }

        r0 = Math.Sqrt(r0);
        if (r0 == 0) return new SolveResult { U = u, Iterations = 0, Relative = 0, Converged = true };
        if (x0 != null && x0.Length == nd)
        {
            for (var i = 0; i < nd; i++) u[i] = fixedDof[i] ? 0 : x0[i];
            var ax = new double[nd];
            Apply(fp, vals, u, ax, fixedDof);
            for (var i = 0; i < nd; i++) r[i] -= ax[i];
        }

        var z = new double[nd];
        var p = new double[nd];
        var q2 = new double[nd];
        void Precond(double[] src, double[] dst)
        {
            Parallel.For(0, nN, v =>
            {
                var m = minv.AsSpan(9 * v, 9);
                double a0 = src[3 * v], a1 = src[3 * v + 1], a2 = src[3 * v + 2];
                dst[3 * v] = fixedDof[3 * v] ? 0 : m[0] * a0 + m[1] * a1 + m[2] * a2;
                dst[3 * v + 1] = fixedDof[3 * v + 1] ? 0 : m[3] * a0 + m[4] * a1 + m[5] * a2;
                dst[3 * v + 2] = fixedDof[3 * v + 2] ? 0 : m[6] * a0 + m[7] * a1 + m[8] * a2;
            });
        }

        double DotV(double[] a, double[] b)
        {
            double acc = 0;
            for (var i = 0; i < nd; i++) acc += a[i] * b[i];
            return acc;
        }

        Precond(r, z);
        Array.Copy(z, p, nd);
        var rz = DotV(r, z);
        var it = 0;
        var rel = Math.Sqrt(DotV(r, r)) / r0;
        if (rel < tol) return new SolveResult { U = u, Iterations = 0, Relative = rel, Converged = true };
        for (; it < maxIt; it++)
        {
            Apply(fp, vals, p, q2, fixedDof);
            var alpha = rz / DotV(p, q2);
            for (var i = 0; i < nd; i++)
            {
                u[i] += alpha * p[i];
                r[i] -= alpha * q2[i];
            }

            rel = Math.Sqrt(DotV(r, r)) / r0;
            if (rel < tol) break;
            Precond(r, z);
            var rz2 = DotV(r, z);
            var beta = rz2 / rz;
            rz = rz2;
            for (var i = 0; i < nd; i++) p[i] = z[i] + beta * p[i];
        }

        return new SolveResult { U = u, Iterations = it + 1, Relative = rel, Converged = rel < tol };
    }

    public static void Invert3(ReadOnlySpan<double> a, Span<double> o)
    {
        var det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6]) +
                  a[2] * (a[3] * a[7] - a[4] * a[6]);
        o[0] = (a[4] * a[8] - a[5] * a[7]) / det;
        o[1] = (a[2] * a[7] - a[1] * a[8]) / det;
        o[2] = (a[1] * a[5] - a[2] * a[4]) / det;
        o[3] = (a[5] * a[6] - a[3] * a[8]) / det;
        o[4] = (a[0] * a[8] - a[2] * a[6]) / det;
        o[5] = (a[2] * a[3] - a[0] * a[5]) / det;
        o[6] = (a[3] * a[7] - a[4] * a[6]) / det;
        o[7] = (a[1] * a[6] - a[0] * a[7]) / det;
        o[8] = (a[0] * a[4] - a[1] * a[3]) / det;
    }

    /// <summary>Element von Mises stress of the live tetrahedra (zero for other rows).</summary>
    public static double[] ElementVonMises(Snapshot s, double[] X, double[] u, Material m)
    {
        var vm = new double[s.ElementCount];
        Parallel.For(0, s.ElementCount, e =>
        {
            var row = s.Row(e);
            if (row.Length != 4) return;
            Span<double> g = stackalloc double[12];
            Span<double> G = stackalloc double[9];
            Span<double> S = stackalloc double[9];
            Gradients(X, row, g);
            Gradient(u, row, g, G);
            Stress(G, m, S);
            vm[e] = VonMises(S);
        });
        return vm;
    }
}

/// <summary>Classification of a boundary facet for the surface terms of the interaction integral.</summary>
public enum FacetKind : byte
{
    /// <summary>Traction-free outer boundary: the free-surface term applies.</summary>
    Free = 0,

    /// <summary>Loaded or constrained outer boundary: excluded from the surface terms.</summary>
    Prescribed = 1,

    /// <summary>Symmetry or plane-strain face: excluded (its term vanishes for in-plane extension fields).</summary>
    Symmetry = 2,

    /// <summary>Crack face (a boundary facet with a twin): the crack-face term applies when enabled.</summary>
    Crack = 3
}

/// <summary>Loads, fixities and facet kinds of a problem on the current mesh.</summary>
public sealed class BoundaryConditions
{
    public required double[] Force;     // 3 per node
    public required bool[] FixedDof;    // 3 per node

    /// <summary>Kind of each boundary facet id (e * 4 + f) that is not free; absent ids are free.</summary>
    public required Dictionary<int, FacetKind> FacetKinds;

    /// <summary>Lumps the traction t (per unit area) on the boundary facets accepted by <paramref name="select"/>.</summary>
    public static void Traction(Snapshot s, double[] X, Func<V3, V3, bool> select, V3 t, double[] f, Dictionary<int, FacetKind> kinds)
    {
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
            V3 a = V3.Of(X, tri[0]), b = V3.Of(X, tri[1]), c = V3.Of(X, tri[2]);
            var av = Geo.AreaVector(a, b, c);
            var area = 0.5 * av.Norm;
            // outward normal: away from the fourth node
            var d = V3.Of(X, row[fl]);
            var n = av.Unit();
            if (V3.Dot(n, (a + b + c) / 3 - d) < 0) n = -n;
            if (!select((a + b + c) / 3, n)) continue;
            foreach (var v in tri)
            {
                f[3 * v] += t.X * area / 3;
                f[3 * v + 1] += t.Y * area / 3;
                f[3 * v + 2] += t.Z * area / 3;
            }

            kinds[id] = FacetKind.Prescribed;
        }
    }

    /// <summary>Marks the boundary facets accepted by <paramref name="select"/> (centroid, outward normal) with a kind.</summary>
    public static void Mark(Snapshot s, double[] X, Func<V3, V3, bool> select, FacetKind kind, Dictionary<int, FacetKind> kinds)
    {
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
            V3 a = V3.Of(X, tri[0]), b = V3.Of(X, tri[1]), c = V3.Of(X, tri[2]);
            var n = Geo.AreaVector(a, b, c).Unit();
            var d = V3.Of(X, row[fl]);
            if (V3.Dot(n, (a + b + c) / 3 - d) < 0) n = -n;
            if (select((a + b + c) / 3, n)) kinds[id] = kind;
        }
    }

    /// <summary>Potential energy of the converged solution under dead load: -f.u/2.</summary>
    public static double PotentialEnergy(double[] f, double[] u)
    {
        double acc = 0;
        for (var i = 0; i < f.Length; i++) acc += f[i] * u[i];
        return -0.5 * acc;
    }
}

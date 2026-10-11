// Meshes.cs — Kuhn boxes in Maubach order, recursive bisection to a size function, radius-ratio quality.
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
// Each cell of an nx x ny x nz grid is split into six tetrahedra sharing the main
// diagonal, the vertices listed along the monotone lattice path from the low to the high
// corner; with tag 3 these are the tagged initial simplices of Maubach's bisection. The
// refiner bisects every element longer than a factor times the size function at its
// centroid, with the compatibility closure, until nothing is marked; template children
// of a cut are not bisected (they carry no tag), and a closure that reaches one stops.

using Step2;

namespace Step3;

public static class KuhnMesh
{
    private static readonly int[][] Perms =
    [
        [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]
    ];

    /// <summary>Rows in Maubach order and coordinates of the box [lo, hi] with nx x ny x nz cells.</summary>
    public static (int[][] Rows, double[] Coords, int NodeCount) Box(int nx, int ny, int nz, V3 lo, V3 hi)
    {
        int Node(int i, int j, int k) => i + (nx + 1) * (j + (ny + 1) * k);
        var nN = (nx + 1) * (ny + 1) * (nz + 1);
        var coords = new double[3 * nN];
        for (var k = 0; k <= nz; k++)
        for (var j = 0; j <= ny; j++)
        for (var i = 0; i <= nx; i++)
        {
            var v = Node(i, j, k);
            coords[3 * v] = lo.X + (hi.X - lo.X) * i / nx;
            coords[3 * v + 1] = lo.Y + (hi.Y - lo.Y) * j / ny;
            coords[3 * v + 2] = lo.Z + (hi.Z - lo.Z) * k / nz;
        }

        var rows = new int[6 * nx * ny * nz][];
        var t = 0;
        Span<int> c = stackalloc int[3];
        for (var k = 0; k < nz; k++)
        for (var j = 0; j < ny; j++)
        for (var i = 0; i < nx; i++)
            foreach (var perm in Perms)
            {
                c[0] = i;
                c[1] = j;
                c[2] = k;
                var row = new int[4];
                row[0] = Node(c[0], c[1], c[2]);
                for (var s = 0; s < 3; s++)
                {
                    c[perm[s]]++;
                    row[s + 1] = Node(c[0], c[1], c[2]);
                }

                rows[t++] = row;
            }

        return (rows, coords, nN);
    }

    /// <summary>Builds the CrackMesh of a box (persistent mesh, first version committed).</summary>
    public static CrackMesh Build(int nx, int ny, int nz, V3 lo, V3 hi)
    {
        var (rows, coords, nN) = Box(nx, ny, nz, lo, hi);
        var flat = FlatBuild.Initial(ElementSpec.Tetrahedron, nN, rows);
        var cm = new CrackMesh(flat, coords);
        cm.Commit();
        return cm;
    }
}

public sealed class MaubachRefiner
{
    private readonly CrackMesh _cm;
    public long Bisections;
    public long Skipped; // elements whose closure reached a template child

    public MaubachRefiner(CrackMesh cm) => _cm = cm;

    private (int A, int B) RefEdge(int e)
    {
        var row = _cm.M.Row(e);
        return (row[0], row[_cm.Tag[e]]);
    }

    private static bool Same((int A, int B) x, int a, int b) => (x.A == a && x.B == b) || (x.A == b && x.B == a);

    private void ElementsWithEdge(int a, int b, List<int> result)
    {
        result.Clear();
        var ra = _cm.M.TransposePairs(a);
        var rb = _cm.M.TransposePairs(b);
        if (rb.Length < ra.Length)
        {
            var tmp = ra;
            ra = rb;
            rb = tmp;
        }
        for (var q = 0; q < ra.Length; q += 2)
            if (_cm.M.Row(ra[q]).Length == 4 && PairRow.IndexOf(rb, ra[q]) >= 0)
                result.Add(ra[q]);
    }

    private sealed class ClosureException : Exception;

    /// <summary>Bisects element e and, recursively, whatever its refinement edge requires; false if a template child blocks it.</summary>
    public bool Bisect(int e)
    {
        try
        {
            BisectCore(e, 0);
            return true;
        }
        catch (ClosureException)
        {
            Skipped++;
            return false;
        }
    }

    private void BisectCore(int e, int depth)
    {
        if (depth > 256) throw new InvalidOperationException("bisection recursion did not terminate");
        if (_cm.Bisectable[e] == 0) throw new ClosureException();
        var (a, b) = RefEdge(e);
        var patch = new List<int>();
        while (true)
        {
            ElementsWithEdge(a, b, patch);
            var incompatible = -1;
            foreach (var g in patch)
            {
                if (_cm.Bisectable[g] == 0) throw new ClosureException();
                if (!Same(RefEdge(g), a, b))
                {
                    incompatible = g;
                    break;
                }
            }

            if (incompatible < 0) break;
            BisectCore(incompatible, depth + 1);
        }

        var z = _cm.NodeOnEdge(a, b, 0.5, 0);
        foreach (var g in patch) Split(g, z);
    }

    private void Split(int g, int z)
    {
        var row = _cm.M.Row(g).ToArray();
        var d = _cm.Tag[g];
        var child1 = (int[])row.Clone();
        child1[d] = z;
        var child2 = new int[4];
        for (var p = 0; p < d; p++) child2[p] = row[p + 1];
        child2[d] = z;
        for (var p = d + 1; p < 4; p++) child2[p] = row[p];
        var nd = d > 1 ? d - 1 : 3;
        var level = _cm.Level[g] + 1;
        var origin = _cm.Origin[g];
        _cm.M.ReplaceElement(g, child1);
        _cm.Tag.Set(g, nd);
        _cm.Level.Set(g, level);
        var e2 = _cm.M.AddElement(child2);
        _cm.Tag.Set(e2, nd);
        _cm.Level.Set(e2, level);
        _cm.Origin.Set(e2, origin);
        _cm.Bisectable.Set(e2, 1);
        Bisections++;
    }

    /// <summary>Longest edge of a row.</summary>
    public double LongestEdge(int e)
    {
        var row = _cm.M.Row(e);
        var l2 = 0.0;
        for (var i = 0; i < 4; i++)
        for (var j = i + 1; j < 4; j++)
        {
            var d = (_cm.P(row[i]) - _cm.P(row[j])).Norm2;
            if (d > l2) l2 = d;
        }

        return Math.Sqrt(l2);
    }

    /// <summary>
    ///     Refines until no element has an edge longer than <paramref name="factor"/> times h at its centroid, or has
    ///     reached <paramref name="maxLevel"/>; <paramref name="filter"/> restricts the candidates. Returns the rounds.
    /// </summary>
    public int RefineToSize(Func<V3, double> h, double factor = 1.5, int maxLevel = 60, Func<int, bool>? filter = null)
    {
        var rounds = 0;
        while (true)
        {
            var marked = new List<int>();
            var n = _cm.M.ElementCount;
            for (var e = 0; e < n; e++)
            {
                var row = _cm.M.Row(e);
                if (row.Length != 4 || _cm.Bisectable[e] == 0 || _cm.Level[e] >= maxLevel) continue;
                if (filter != null && !filter(e)) continue;
                var c = Geo.Centroid(_cm.Coordinates4(row), [0, 1, 2, 3]);
                if (LongestEdge(e) > factor * h(c)) marked.Add(e);
            }

            if (marked.Count == 0) return rounds;
            foreach (var e in marked)
                if (_cm.M.Row(e).Length == 4 && _cm.Level[e] < maxLevel)
                    Bisect(e);
            rounds++;
            if (rounds > 200) throw new InvalidOperationException("refinement did not converge");
        }
    }
}

public static class CrackMeshExtensions
{
    /// <summary>Coordinates of the four nodes of a row as a flat array (local labels 0..3).</summary>
    public static double[] Coordinates4(this CrackMesh cm, ReadOnlySpan<int> row)
    {
        var x = new double[12];
        for (var k = 0; k < 4; k++)
        {
            x[3 * k] = cm.X[row[k]];
            x[3 * k + 1] = cm.Y[row[k]];
            x[3 * k + 2] = cm.Z[row[k]];
        }

        return x;
    }
}

public static class MeshQuality
{
    /// <summary>Radius ratio 3 r_in / r_circ (1 for the regular tetrahedron).</summary>
    public static double RadiusRatio(V3 a, V3 b, V3 c, V3 d)
    {
        var vol = Math.Abs(Geo.Volume6(a, b, c, d)) / 6;
        var area = 0.5 * (Geo.AreaVector(a, b, c).Norm + Geo.AreaVector(a, b, d).Norm + Geo.AreaVector(a, c, d).Norm +
                          Geo.AreaVector(b, c, d).Norm);
        var rin = 3 * vol / area;
        // circumradius from the Cayley-Menger-like formula: R = |ab x ac . ad| ... use the standard expression
        V3 ab = b - a, ac = c - a, ad = d - a;
        var num = (ab.Norm2 * V3.Cross(ac, ad) + ac.Norm2 * V3.Cross(ad, ab) + ad.Norm2 * V3.Cross(ab, ac)).Norm;
        var rc = num / (12 * vol);
        return 3 * rin / rc;
    }

    /// <summary>Minimum and mean radius ratio over the live tetrahedra of a snapshot.</summary>
    public static (double Min, double Mean, int Count) Stats(Snapshot s, double[] X)
    {
        double min = 1, sum = 0;
        var cnt = 0;
        for (var e = 0; e < s.ElementCount; e++)
        {
            var row = s.Row(e);
            if (row.Length != 4) continue;
            var q = RadiusRatio(V3.Of(X, row[0]), V3.Of(X, row[1]), V3.Of(X, row[2]), V3.Of(X, row[3]));
            min = Math.Min(min, q);
            sum += q;
            cnt++;
        }

        return (min, cnt > 0 ? sum / cnt : 0, cnt);
    }
}

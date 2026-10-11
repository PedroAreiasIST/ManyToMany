// Interfaces.cs — inclusion interfaces inserted by concurrent cutting.
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
//   interfaces n g threads out
//
// Kuhn box [-1,1]^3 of 6 n^3 tetrahedra; g^3 spheres of radius 0.35 (2/g) centered in the
// cells of a g x g x g grid. Every edge crossing a sphere is cut at its root (one group per
// sphere, the node placed by the observer from the endpoint with the smaller label), by
// CutParallel on one copy of the mesh and by the sequential Cut of the concatenated groups
// on another; the two must give the same rows. The interface of a sphere is the set of
// faces shared by a child inside and a child outside it, all three nodes on the sphere: a
// closed triangulated surface. The interfaces are written to <out>/interfaces.surf (nodes,
// triangles oriented outward, sphere of each triangle) for the vector figure.

using System.Diagnostics;
using Step2;

namespace Step2.Bench;

public static class Interfaces
{
    private sealed class RootPlacer(Attribute<double> x, Attribute<double> y, Attribute<double> z,
        Dictionary<long, double> roots) : ICutObserver
    {
        public void NodeCreated(int node, int a, int b)
        {
            int lo = Math.Min(a, b), hi = Math.Max(a, b);
            var t = roots[Cutter.Key(lo, hi)];
            x.Set(node, x[lo] + t * (x[hi] - x[lo]));
            y.Set(node, y[lo] + t * (y[hi] - y[lo]));
            z.Set(node, z[lo] + t * (z[hi] - z[lo]));
        }

        public void ElementSplit(int parent, int child) { }
    }

    public static void Run(string[] a)
    {
        int n = a.Length > 0 ? int.Parse(a[0]) : 26, g = a.Length > 1 ? int.Parse(a[1]) : 2;
        var threads = a.Length > 2 ? int.Parse(a[2]) : Environment.ProcessorCount;
        var outDir = a.Length > 3 ? a[3] : ".";
        var rad = 0.35 * 2.0 / g;
        var (rows, coords01, nN) = Kuhn.Box(n);
        var coords = coords01.Select(c => 2 * c - 1).ToArray();
        (double, double, double) Center(int c)
        {
            int i = c / (g * g), j = c / g % g, k = c % g;
            return (-1 + (i + 0.5) * 2.0 / g, -1 + (j + 0.5) * 2.0 / g, -1 + (k + 0.5) * 2.0 / g);
        }

        double Phi(double x, double y, double z, int c)
        {
            var (cx, cy, cz) = Center(c);
            return Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz)) - rad;
        }

        int CellOf(double x, double y, double z)
        {
            int Ix(double t) => Math.Clamp((int)((t + 1) / 2 * g), 0, g - 1);
            return (Ix(x) * g + Ix(y)) * g + Ix(z);
        }

        IncrementalMesh NewMesh(out Attribute<double> X, out Attribute<double> Y, out Attribute<double> Z)
        {
            var flat = FlatBuild.Initial(ElementSpec.Tetrahedron, nN, rows);
            var m = new IncrementalMesh(flat, ElementSpec.Tetrahedron, persistent: true);
            X = m.AddAttribute(AttributeKind.Node, 0.0);
            Y = m.AddAttribute(AttributeKind.Node, 0.0);
            Z = m.AddAttribute(AttributeKind.Node, 0.0);
            for (var v = 0; v < nN; v++)
            {
                X.Set(v, coords[3 * v]);
                Y.Set(v, coords[3 * v + 1]);
                Z.Set(v, coords[3 * v + 2]);
            }

            m.Commit();
            return m;
        }

        var mp = NewMesh(out var xp, out var yp, out var zp);
        var ms = NewMesh(out var xs, out var ys, out var zs);

        // marked edges: crossings of the sphere of the cell of the lower endpoint, with their roots
        var groups = new List<(int, int)>[g * g * g];
        for (var c = 0; c < groups.Length; c++) groups[c] = new List<(int, int)>();
        var roots = new Dictionary<long, double>(EdgeKeyComparer.Instance);
        foreach (var (p, q) in CutVerify.Edges(mp))
        {
            var c = CellOf(coords[3 * p], coords[3 * p + 1], coords[3 * p + 2]);
            double fp = Phi(coords[3 * p], coords[3 * p + 1], coords[3 * p + 2], c),
                fq = Phi(coords[3 * q], coords[3 * q + 1], coords[3 * q + 2], c);
            if (!(fp * fq < 0)) continue;
            int lo = Math.Min(p, q), hi = Math.Max(p, q);
            var cc = c;
            roots[Cutter.Key(lo, hi)] = CrackInserter.RootOf(coords[3 * lo], coords[3 * lo + 1], coords[3 * lo + 2],
                coords[3 * hi], coords[3 * hi + 1], coords[3 * hi + 2], (x, y, z) => Phi(x, y, z, cc));
            groups[c].Add((lo, hi));
        }

        foreach (var l in groups) l.Sort();
        var marked = groups.Sum(l => l.Count);
        var sw = Stopwatch.StartNew();
        var res = Cutter.CutParallel(mp, groups.Select(l => (IReadOnlyList<(int, int)>)l).ToList(),
            new RootPlacer(xp, yp, zp, roots), threads);
        var tPar = sw.Elapsed.TotalMilliseconds;
        sw.Restart();
        Cutter.Cut(ms, groups.SelectMany(l => l).ToList(), new RootPlacer(xs, ys, zs, roots));
        var tSeq = sw.Elapsed.TotalMilliseconds;
        var sp = mp.Commit();
        var ss = ms.Commit();
        CutVerify.SameRows(Enumerable.Range(0, sp.ElementCount).Select(e => sp.Row(e).ToArray()),
            Enumerable.Range(0, ss.ElementCount).Select(e => ss.Row(e).ToArray()), "interfaces");
        Oracle.Verify(sp, checkAssembly: false);
        Console.WriteLine($"  interfaces n={n}, {groups.Length} spheres: {marked} cut edges, {res.Groups.Sum(r => r.ShellElements)} " +
                          $"split elements, {res.Merged} merged; concurrent {tPar:F1} ms ({threads} threads), sequential " +
                          $"{tSeq:F1} ms; rows equal, oracle passes");

        // interface faces: shared by a child inside and a child outside, three nodes on the sphere
        var xv = sp.View(xp);
        var yv = sp.View(yp);
        var zv = sp.View(zp);
        var first = res.Groups.Length > 0 ? res.Groups[0].FirstNewNode : sp.NodeCount;
        bool On(int v) => v >= first;
        // a child lies on the side of its nodes off the sphere (every child keeps one, all on one side); the
        // centroid would not do: the sphere bulges beyond the flat interface faces
        int Side(ReadOnlySpan<int> r, out int cell)
        {
            cell = -1;
            foreach (var v in r)
            {
                if (On(v)) continue;
                cell = CellOf(xv[v], yv[v], zv[v]);
                return Phi(xv[v], yv[v], zv[v], cell) < 0 ? -1 : 1;
            }

            return 0;
        }

        var tris = new List<(int A, int B, int C, int Sphere)>();
        for (var e = 0; e < sp.ElementCount; e++)
        {
            var row = sp.Row(e);
            if (row.Length != 4) continue;
            if (Side(row, out var cell) >= 0) continue; // faces taken from the inside child
            for (var f = 0; f < 4; f++)
            {
                var nb = sp.FacetNeighbour(e, f);
                if (nb < 0) continue;
                int i0 = row[(f + 1) % 4], i1 = row[(f + 2) % 4], i2 = row[(f + 3) % 4];
                if (!On(i0) || !On(i1) || !On(i2)) continue;
                if (Side(sp.Row(nb), out _) <= 0) continue;
                // outward: away from the sphere center
                var (sx, sy, sz) = Center(cell);
                double ux = xv[i1] - xv[i0], uy = yv[i1] - yv[i0], uz = zv[i1] - zv[i0];
                double wx = xv[i2] - xv[i0], wy = yv[i2] - yv[i0], wz = zv[i2] - zv[i0];
                double nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
                var outward = nx * (xv[i0] - sx) + ny * (yv[i0] - sy) + nz * (zv[i0] - sz) > 0;
                tris.Add(outward ? (i0, i1, i2, cell) : (i0, i2, i1, cell));
            }
        }

        // closedness: every interface edge shared by two triangles of the same sphere
        var edgeUse = new Dictionary<long, int>(EdgeKeyComparer.Instance);
        foreach (var (p, q, r, _) in tris)
            foreach (var (u, w) in new[] { (p, q), (q, r), (r, p) })
                edgeUse[Cutter.Key(u, w)] = edgeUse.GetValueOrDefault(Cutter.Key(u, w)) + 1;
        var open = edgeUse.Count(kv => kv.Value != 2);
        double maxDev = 0;
        foreach (var (p, q, r, c) in tris)
            foreach (var v in new[] { p, q, r })
                maxDev = Math.Max(maxDev, Math.Abs(Phi(xv[v], yv[v], zv[v], c)));
        Console.WriteLine($"  interface: {tris.Count} triangles, {edgeUse.Count} edges, {open} not shared by two " +
                          $"triangles; max |phi| at the nodes {maxDev:E1}");

        Directory.CreateDirectory(outDir);
        var used = tris.SelectMany(t => new[] { t.A, t.B, t.C }).Distinct().OrderBy(v => v).ToArray();
        var map = new Dictionary<int, int>();
        for (var i = 0; i < used.Length; i++) map[used[i]] = i;
        using (var bw = new BinaryWriter(File.Create(Path.Combine(outDir, "interfaces.surf"))))
        {
            bw.Write(used.Length);
            bw.Write(tris.Count);
            foreach (var v in used)
            {
                bw.Write(xv[v]);
                bw.Write(yv[v]);
                bw.Write(zv[v]);
            }

            foreach (var t in tris)
            {
                bw.Write(map[t.A]);
                bw.Write(map[t.B]);
                bw.Write(map[t.C]);
            }

            foreach (var t in tris) bw.Write(t.Sphere);
        }

        Console.WriteLine($"interfaces n={n} g={g}: written {Path.Combine(outDir, "interfaces.surf")}");
    }
}

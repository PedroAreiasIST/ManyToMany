// EnsightGold.cs — EnSight Gold ASCII output: body, crack faces per sheet, fronts, junctions, outline.
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
// Files prefix.geoNNNN and prefix.<var>NNNN in the conventions of ist_ensight_gold.py
// (node id assign, one part per entity); the case file lists the steps. Parts: body
// (tetra4, optional), crack_face_k (tria3, the positive-side twins of sheet k), front_k
// (bar2, every chain), junction (bar2) and outline (bar2, the 12 edges of a box).
// Node variables: displacement, umag, mises (body), opening and sliding (crack faces).

using System.Globalization;
using System.Text;
using Step2;

namespace Step3;

public static class EnsightGold
{
    private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

    private sealed class Part
    {
        public int Id;
        public string Name = "";
        public int[] Nodes = [];
        public double[] Xyz = [];
        public string Type = "";
        public int[] Conn = [];
        public int NPerElem;
        public bool IsBody, IsCrack;
        public double[]? Opening, Sliding;
    }

    private static string F12(double v) => v.ToString("0.00000e+00", Ci).PadLeft(12);

    /// <summary>Writes one step; returns the number of parts written.</summary>
    public static int WriteStep(string dir, string prefix, int stepNo, Snapshot s, CrackMesh cm, CrackCount cc, double[] X,
        double[] u, Material mat, bool withBody, (V3 Lo, V3 Hi)? outline, string title = "Step3 crack growth")
    {
        Directory.CreateDirectory(dir);
        var tag = stepNo.ToString("D4", Ci);
        var nNodes = s.NodeCount;
        var parts = new List<Part>();
        if (withBody)
        {
            var live = new List<int>();
            for (var e = 0; e < s.ElementCount; e++)
                if (s.Row(e).Length == 4)
                    live.Add(e);
            var used = new bool[nNodes];
            foreach (var e in live)
                foreach (var v in s.Row(e))
                    used[v] = true;
            var newId = new int[nNodes];
            var order = new List<int>();
            for (var v = 0; v < nNodes; v++)
                if (used[v])
                {
                    order.Add(v);
                    newId[v] = order.Count;
                }

            var conn = new int[4 * live.Count];
            for (var k = 0; k < live.Count; k++)
            {
                var r = s.Row(live[k]);
                for (var c = 0; c < 4; c++) conn[4 * k + c] = newId[r[c]];
            }

            parts.Add(new Part { Name = "body", Nodes = order.ToArray(), Type = "tetra4", Conn = conn, NPerElem = 4, IsBody = true });
        }

        // crack faces per sheet: the positive-side twins, with opening and sliding from the jump to the negative twin
        for (var sheet = 0; sheet < cc.SheetCount; sheet++)
        {
            var tris = new List<int[]>();
            var jumpB = new List<double>();
            var jumpS = new List<double>();
            var cmap = new Dictionary<int, int>();
            var corder = new List<int>();
            var opening = new List<double>();
            var sliding = new List<double>();
            var tri = new int[3];
            for (var i = 0; i < cc.PairCount; i++)
            {
                if (cc.SheetOfPair[i] != sheet) continue;
                var fp = cc.FacetPlus[i];
                var fm = cc.FacetMinus[i];
                var row = s.Row(fp / 4);
                var k = 0;
                for (var p = 0; p < 4; p++)
                    if (p != fp % 4)
                        tri[k++] = row[p];
                var rowM = s.Row(fm / 4);
                var triM = new int[3];
                k = 0;
                for (var p = 0; p < 4; p++)
                    if (p != fm % 4)
                        triM[k++] = rowM[p];
                var (a, b, c) = cc.Tri[i];
                var nrm = V3.Cross(cc.Point(b) - cc.Point(a), cc.Point(c) - cc.Point(a)).Unit();
                var t = (int[])tri.Clone();
                tris.Add(t);
                foreach (var v in t)
                {
                    if (cmap.ContainsKey(v)) continue;
                    // the twin node: same canonical label on the minus facet
                    var twin = -1;
                    foreach (var w in triM)
                        if (cm.Canon[w] == cm.Canon[v])
                            twin = w;
                    var jump = twin >= 0 ? V3.Of(u, v) - V3.Of(u, twin) : V3.Zero;
                    var op = V3.Dot(jump, nrm);
                    cmap[v] = corder.Count + 1;
                    corder.Add(v);
                    opening.Add(op);
                    sliding.Add((jump - op * nrm).Norm);
                }
            }

            if (tris.Count == 0) continue;
            var cconn = new int[3 * tris.Count];
            for (var k = 0; k < tris.Count; k++)
                for (var c = 0; c < 3; c++)
                    cconn[3 * k + c] = cmap[tris[k][c]];
            parts.Add(new Part
            {
                Name = $"crack_face_{sheet}", Nodes = corder.ToArray(), Type = "tria3", Conn = cconn, NPerElem = 3, IsCrack = true,
                Opening = opening.ToArray(), Sliding = sliding.ToArray()
            });
        }

        for (var ci = 0; ci < cc.Chains.Count; ci++)
        {
            var ch = cc.Chains[ci];
            var nseg = ch.Closed ? ch.Count : ch.Count - 1;
            var fconn = new int[2 * nseg];
            for (var j = 0; j < nseg; j++)
            {
                fconn[2 * j] = j + 1;
                fconn[2 * j + 1] = (j + 1) % ch.Count + 1;
            }

            var xyz = new double[3 * ch.Count];
            for (var j = 0; j < ch.Count; j++) ch.P[j].To(xyz, j);
            parts.Add(new Part { Name = $"front_{ci}", Nodes = new int[ch.Count], Xyz = xyz, Type = "bar2", Conn = fconn, NPerElem = 2 });
        }

        if (cc.JunctionEdges.Count > 0)
        {
            var jn = new Dictionary<int, int>();
            var jx = new List<double>();
            var jc = new List<int>();
            foreach (var (a, b) in cc.JunctionEdges)
            {
                foreach (var v in new[] { a, b })
                    if (!jn.ContainsKey(v))
                    {
                        jn[v] = jn.Count + 1;
                        var p = cc.Point(v);
                        jx.AddRange([p.X, p.Y, p.Z]);
                    }

                jc.Add(jn[a]);
                jc.Add(jn[b]);
            }

            parts.Add(new Part { Name = "junction", Nodes = new int[jn.Count], Xyz = jx.ToArray(), Type = "bar2", Conn = jc.ToArray(), NPerElem = 2 });
        }

        if (outline != null)
        {
            var (lo, hi) = outline.Value;
            var corners = new double[24];
            for (var c = 0; c < 8; c++)
            {
                corners[3 * c] = (c & 1) != 0 ? hi.X : lo.X;
                corners[3 * c + 1] = (c & 2) != 0 ? hi.Y : lo.Y;
                corners[3 * c + 2] = (c & 4) != 0 ? hi.Z : lo.Z;
            }

            int[] edges = [1, 2, 1, 3, 2, 4, 3, 4, 5, 6, 5, 7, 6, 8, 7, 8, 1, 5, 2, 6, 3, 7, 4, 8];
            parts.Add(new Part { Name = "outline", Nodes = new int[8], Xyz = corners, Type = "bar2", Conn = edges, NPerElem = 2 });
        }

        for (var k = 0; k < parts.Count; k++)
        {
            parts[k].Id = k + 1;
            if (parts[k].Xyz.Length == 0)
            {
                var xs = new double[3 * parts[k].Nodes.Length];
                for (var q = 0; q < parts[k].Nodes.Length; q++)
                for (var c = 0; c < 3; c++)
                    xs[3 * q + c] = X[3 * parts[k].Nodes[q] + c];
                parts[k].Xyz = xs;
            }
        }

        var enc = new UTF8Encoding(false);
        using (var w = new StreamWriter(Path.Combine(dir, $"{prefix}.geo{tag}"), false, enc, 1 << 20))
        {
            w.WriteLine($"{title}, step {stepNo}");
            w.WriteLine("Step3: crack growth and intersection on the maintained structures");
            w.WriteLine("node id assign");
            w.WriteLine("element id assign");
            foreach (var part in parts)
            {
                var n = part.Xyz.Length / 3;
                w.WriteLine("part");
                w.WriteLine(part.Id.ToString(Ci).PadLeft(10));
                w.WriteLine(part.Name);
                w.WriteLine("coordinates");
                w.WriteLine(n.ToString(Ci).PadLeft(10));
                for (var c = 0; c < 3; c++)
                    for (var q = 0; q < n; q++)
                        w.WriteLine(F12(part.Xyz[3 * q + c]));
                w.WriteLine(part.Type);
                w.WriteLine((part.Conn.Length / part.NPerElem).ToString(Ci).PadLeft(10));
                var sb = new StringBuilder();
                for (var k = 0; k < part.Conn.Length / part.NPerElem; k++)
                {
                    sb.Clear();
                    for (var c = 0; c < part.NPerElem; c++) sb.Append(part.Conn[part.NPerElem * k + c].ToString(Ci).PadLeft(10));
                    w.WriteLine(sb.ToString());
                }
            }
        }

        var vmNode = new double[nNodes];
        if (withBody)
        {
            var vmElem = P1Elasticity.ElementVonMises(s, X, u, mat);
            var cnt = new int[nNodes];
            for (var e = 0; e < s.ElementCount; e++)
            {
                var row = s.Row(e);
                if (row.Length != 4) continue;
                foreach (var v in row)
                {
                    vmNode[v] += vmElem[e];
                    cnt[v]++;
                }
            }

            for (var v = 0; v < nNodes; v++)
                if (cnt[v] > 0)
                    vmNode[v] /= cnt[v];
        }

        void Variable(string file, string label, int ncomp, Func<Part, int, int, double> val)
        {
            using var w = new StreamWriter(Path.Combine(dir, $"{prefix}.{file}{tag}"), false, enc, 1 << 20);
            w.WriteLine(label);
            foreach (var part in parts)
            {
                w.WriteLine("part");
                w.WriteLine(part.Id.ToString(Ci).PadLeft(10));
                w.WriteLine("coordinates");
                for (var c = 0; c < ncomp; c++)
                    for (var q = 0; q < part.Nodes.Length; q++)
                        w.WriteLine(F12(val(part, q, c)));
            }
        }

        Variable("umag", "umag", 1, (part, q, _) => part.IsBody || part.IsCrack ? V3.Of(u, part.Nodes[q]).Norm : 0);
        Variable("open", "opening", 1, (part, q, _) => part.IsCrack ? part.Opening![q] : 0);
        Variable("slide", "sliding", 1, (part, q, _) => part.IsCrack ? part.Sliding![q] : 0);
        Variable("disp", "displacement", 3, (part, q, c) => part.IsBody || part.IsCrack ? u[3 * part.Nodes[q] + c] : 0);
        if (withBody) Variable("mises", "mises", 1, (part, q, _) => part.IsBody ? vmNode[part.Nodes[q]] : 0);
        return parts.Count;
    }

    public static void WriteCase(string dir, string prefix, bool withBody, IReadOnlyList<int> stepNos, IReadOnlyList<double> times)
    {
        using var w = new StreamWriter(Path.Combine(dir, $"{prefix}.case"));
        w.WriteLine("FORMAT");
        w.WriteLine("type: ensight gold");
        w.WriteLine();
        w.WriteLine("GEOMETRY");
        w.WriteLine($"model: 1 {prefix}.geo****");
        w.WriteLine();
        w.WriteLine("VARIABLE");
        w.WriteLine($"vector per node: 1 displacement {prefix}.disp****");
        w.WriteLine($"scalar per node: 1 umag {prefix}.umag****");
        w.WriteLine($"scalar per node: 1 opening {prefix}.open****");
        w.WriteLine($"scalar per node: 1 sliding {prefix}.slide****");
        if (withBody) w.WriteLine($"scalar per node: 1 mises {prefix}.mises****");
        w.WriteLine();
        w.WriteLine("TIME");
        w.WriteLine("time set: 1");
        w.WriteLine($"number of steps: {stepNos.Count}");
        var inc = stepNos.Count > 1 ? stepNos[1] - stepNos[0] : 1;
        var uniform = inc > 0 && Enumerable.Range(0, stepNos.Count - 1).All(i => stepNos[i + 1] - stepNos[i] == inc);
        if (uniform)
        {
            w.WriteLine($"filename start number: {stepNos[0]}");
            w.WriteLine($"filename increment: {inc}");
        }
        else
        {
            w.WriteLine("filename numbers:");
            foreach (var n in stepNos) w.WriteLine(n);
        }

        w.WriteLine("time values:");
        foreach (var t in times) w.WriteLine(t.ToString("0.00000000e+00", Ci));
    }
}

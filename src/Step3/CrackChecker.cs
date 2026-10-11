// CrackChecker.cs — independent checks of the crack facets, twins, fronts and chains of a version.
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
// The checker identifies crack facets geometrically, by coincident centroids of boundary
// facets, without the canonical labels CrackCount relies on, and compares: the same set of
// twin pairs, twins with their tetrahedra on opposite sides, a common crack bit on the
// three nodes of every crack facet, the front, mouth and junction edges recounted from the
// pairs, chains that cover every front edge once, and twins with identical coordinates.
// A mismatch throws.

using Step2;

namespace Step3;

public static class CrackChecker
{
    public sealed class Mismatch(string message) : Exception(message);

    private static long Hash(V3 p, double tol)
    {
        long a = (long)Math.Round(p.X / tol), b = (long)Math.Round(p.Y / tol), c = (long)Math.Round(p.Z / tol);
        return a * 73856093L ^ b * 19349663L ^ c * 83492791L;
    }

    public static (int Pairs, int Front, int Mouth, int Junction) Verify(Snapshot s, CrackMesh cm, CrackCount cc, double tol = 1e-9)
    {
        var tri = new int[3];
        var facets = new List<(int Id, V3 Cen, int[] Tri)>();
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
            var cen = (cm.P(tri[0]) + cm.P(tri[1]) + cm.P(tri[2])) / 3;
            facets.Add((id, cen, (int[])tri.Clone()));
        }

        // geometric twins: coincident centroids (within tol, through a hash with neighbour probing)
        var byHash = new Dictionary<long, List<int>>();
        for (var i = 0; i < facets.Count; i++)
        {
            var h = Hash(facets[i].Cen, tol * 100);
            if (!byHash.TryGetValue(h, out var l)) byHash[h] = l = new List<int>();
            l.Add(i);
        }

        var pairs = new HashSet<(int, int)>();
        var paired = new int[facets.Count];
        Array.Fill(paired, -1);
        foreach (var l in byHash.Values)
            for (var a = 0; a < l.Count; a++)
            for (var b = a + 1; b < l.Count; b++)
            {
                int i = l[a], j = l[b];
                if (V3.Dist(facets[i].Cen, facets[j].Cen) > tol) continue;
                // the three vertices must coincide too
                var ok = facets[i].Tri.All(v => facets[j].Tri.Any(w => V3.Dist(cm.P(v), cm.P(w)) <= tol));
                if (!ok) continue;
                if (paired[i] >= 0 || paired[j] >= 0) throw new Mismatch($"facet {facets[i].Id} or {facets[j].Id} has more than one twin");
                paired[i] = j;
                paired[j] = i;
                pairs.Add((Math.Min(facets[i].Id, facets[j].Id), Math.Max(facets[i].Id, facets[j].Id)));
            }

        var ccPairs = new HashSet<(int, int)>();
        for (var i = 0; i < cc.PairCount; i++)
            ccPairs.Add((Math.Min(cc.FacetPlus[i], cc.FacetMinus[i]), Math.Max(cc.FacetPlus[i], cc.FacetMinus[i])));
        if (!pairs.SetEquals(ccPairs))
            throw new Mismatch($"crack facet pairs: geometric {pairs.Count} vs counted {ccPairs.Count} (common {pairs.Intersect(ccPairs).Count()})");

        // sides, bits, coordinates of the twins
        for (var i = 0; i < cc.PairCount; i++)
        {
            int fp = cc.FacetPlus[i], fm = cc.FacetMinus[i];
            var rp = s.Row(fp / 4);
            var rm = s.Row(fm / 4);
            var k = 0;
            for (var p = 0; p < 4; p++)
                if (p != fp % 4)
                    tri[k++] = rp[p];
            V3 a = cm.P(tri[0]), b = cm.P(tri[1]), c = cm.P(tri[2]);
            var n = Geo.AreaVector(a, b, c);
            var cen = (a + b + c) / 3;
            double sp = V3.Dot(cm.P(rp[fp % 4]) - cen, n), sm = V3.Dot(cm.P(rm[fm % 4]) - cen, n);
            if (sp * sm >= 0) throw new Mismatch($"pair {i}: twins on the same side");
            var bits = -1;
            foreach (var v in tri) bits &= cm.Mask[v];
            var bitsM = -1;
            for (var p = 0; p < 4; p++)
                if (p != fm % 4)
                    bitsM &= cm.Mask[rm[p]];
            if (bits == 0 && bitsM == 0) throw new Mismatch($"pair {i}: no common crack bit on either twin");
        }

        // edge recount from the pairs by canonical labels (a different path: via the pair list, not the facet scan)
        var count = new Dictionary<long, int>();
        var outer = new HashSet<long>();
        var pairedIds = new HashSet<int>(ccPairs.SelectMany(p => new[] { p.Item1, p.Item2 }));
        foreach (var (id, _, t) in facets)
        {
            int c0 = cm.Canon[t[0]], c1 = cm.Canon[t[1]], c2 = cm.Canon[t[2]];
            var keys = new[] { CrackCount.EdgeKey(c0, c1), CrackCount.EdgeKey(c1, c2), CrackCount.EdgeKey(c0, c2) };
            if (!pairedIds.Contains(id))
            {
                foreach (var kk in keys) outer.Add(kk);
                continue;
            }

            if (id != Math.Min(id, PairTwin(cc, id))) continue; // count each pair once
            foreach (var kk in keys) count[kk] = count.GetValueOrDefault(kk) + 1;
        }

        int front = 0, mouth = 0, junction = 0;
        var frontSet = new HashSet<long>();
        foreach (var kv in count)
        {
            if (kv.Value == 1)
            {
                if (outer.Contains(kv.Key)) mouth++;
                else
                {
                    front++;
                    frontSet.Add(kv.Key);
                }
            }
            else if (kv.Value >= 3) junction++;
        }

        if (front != cc.FrontEdges.Count || mouth != cc.MouthEdges.Count || junction != cc.JunctionEdges.Count)
            throw new Mismatch($"edge counts: front {front}/{cc.FrontEdges.Count} mouth {mouth}/{cc.MouthEdges.Count} junction {junction}/{cc.JunctionEdges.Count}");
        foreach (var (a, b) in cc.FrontEdges)
            if (!frontSet.Contains(CrackCount.EdgeKey(a, b)))
                throw new Mismatch($"front edge ({a},{b}) not found by the recount");

        // chains cover every front edge once
        var covered = new HashSet<long>();
        foreach (var ch in cc.Chains)
            foreach (var (i, j) in ch.Segments())
            {
                var key = CrackCount.EdgeKey(ch.Nodes[i], ch.Nodes[j]);
                if (!frontSet.Contains(key)) throw new Mismatch($"chain segment ({ch.Nodes[i]},{ch.Nodes[j]}) is not a front edge");
                if (!covered.Add(key)) throw new Mismatch($"front edge ({ch.Nodes[i]},{ch.Nodes[j]}) covered twice");
            }

        if (covered.Count != frontSet.Count) throw new Mismatch($"chains cover {covered.Count} of {frontSet.Count} front edges");
        return (pairs.Count, front, mouth, junction);
    }

    private static int PairTwin(CrackCount cc, int facet)
    {
        var i = cc.PairOfFacet[facet];
        return cc.FacetPlus[i] == facet ? cc.FacetMinus[i] : cc.FacetPlus[i];
    }

    /// <summary>Largest error of a nodal field that is linear in the coordinates (transfer check).</summary>
    public static double LinearFieldError(Snapshot s, CrackMesh cm, Attribute<double> field, Func<V3, double> exact)
    {
        var view = s.View(field);
        var err = 0.0;
        for (var v = 0; v < s.NodeCount; v++)
            if (s.TransposePairs(v).Length > 0)
                err = Math.Max(err, Math.Abs(view[v] - exact(cm.P(v))));
        return err;
    }
}

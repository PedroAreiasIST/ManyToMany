// Workloads.cs — topology-changing edit sequences on Kuhn meshes.
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
//   CrackPlane     a planar crack x = 1/2 advancing in +y: at each step the nodes of
//                  one front row are duplicated and the elements on the + side take
//                  the duplicates (node duplication, as in crack insertion).
//   Bisection      Maubach newest-vertex bisection with recursive completion: to
//                  bisect an element, every element sharing its refinement edge is
//                  first made compatible by recursive bisection; queries (elements
//                  containing an edge) are interleaved with edits.
//   Erosion        elements whose centroid enters a moving ball are deleted.
//   RandomEdits    randomized edits for differential testing and crossover studies.
//
// All edits go through MeshState, which forwards them to the IncrementalMesh and,
// when present, to an MM3 O2M mirror so that the baseline sees the same rows.

using Numerical;

namespace Step2.Bench;

public sealed class Mm3Mirror
{
    public readonly O2M Rows;

    public Mm3Mirror(int[][] rows)
    {
        Rows = new O2M(rows.Length);
        foreach (var r in rows) Rows.AppendElementCopy(r);
    }

    public void Replace(int e, int[] row) => Rows.ReplaceElementCopy(e, row);
    public void Append(int[] row) => Rows.AppendElementCopy(row);
    public void Clear(int e) => Rows.ClearElement(e);

    /// <summary>
    ///     MM3 from-scratch derivation of the structures the incremental layer maintains:
    ///     transpose, element positions, node positions and clique destinations
    ///     (GetCliques), and the node-node pattern (transpose times primary).
    /// </summary>
    public object Rebuild()
    {
        var t = Rows.Transpose();
        var ep = O2M.GetElementPositions(Rows, t);
        var cl = O2M.GetCliques(Rows, t);
        var nn = t * Rows;
        return (t, ep, cl, nn);
    }
}

public sealed class MeshState
{
    public IncrementalMesh Mesh;
    public readonly Attribute<double> Xa, Ya, Za; // node coordinates, versioned with the mesh
    public readonly Attribute<int> Tag;            // Maubach tag per element
    public readonly Attribute<int> Level;          // refinement level per element
    public Mm3Mirror? Mirror;

    public MeshState(IncrementalMesh mesh, List<double> coords, Mm3Mirror? mirror)
    {
        Mesh = mesh;
        Mirror = mirror;
        Xa = mesh.AddAttribute(AttributeKind.Node, 0.0);
        Ya = mesh.AddAttribute(AttributeKind.Node, 0.0);
        Za = mesh.AddAttribute(AttributeKind.Node, 0.0);
        for (var v = 0; v < mesh.NodeCount; v++)
        {
            Xa.Set(v, coords[3 * v]);
            Ya.Set(v, coords[3 * v + 1]);
            Za.Set(v, coords[3 * v + 2]);
        }

        Tag = mesh.AddAttribute(AttributeKind.Element, 3);
        Level = mesh.AddAttribute(AttributeKind.Element, 0);
        mesh.Commit(); // the first committed version carries the attributes
    }

    public int AddNode(double x, double y, double z)
    {
        var v = Mesh.AddNode();
        Xa.Set(v, x);
        Ya.Set(v, y);
        Za.Set(v, z);
        return v;
    }

    // Edits destined for the MM3 mirror are logged and replayed outside timed regions.
    private readonly List<(int Kind, int E, int[]? Row)> _mirrorLog = new();

    public void Replace(int e, int[] row)
    {
        Mesh.ReplaceElement(e, row);
        if (Mirror != null) _mirrorLog.Add((0, e, row));
    }

    public int Add(int[] row, int tag, int level)
    {
        var e = Mesh.AddElement(row);
        Tag.Set(e, tag);
        Level.Set(e, level);
        if (Mirror != null) _mirrorLog.Add((1, e, row));
        return e;
    }

    public void Remove(int e)
    {
        Mesh.RemoveElement(e);
        if (Mirror != null) _mirrorLog.Add((2, e, null));
    }

    /// <summary>Applies the logged edits to the MM3 mirror (untimed).</summary>
    public void FlushMirror()
    {
        if (Mirror == null) return;
        foreach (var (kind, e, row) in _mirrorLog)
            switch (kind)
            {
                case 0: Mirror.Replace(e, row!); break;
                case 1: Mirror.Append(row!); break;
                default: Mirror.Clear(e); break;
            }

        _mirrorLog.Clear();
    }

    /// <summary>Discards logged edits (after a rollback of the incremental mesh).</summary>
    public void DiscardMirrorLog() => _mirrorLog.Clear();

    public double X(int v, int c) => c == 0 ? Xa[v] : c == 1 ? Ya[v] : Za[v];

    public (double X, double Y, double Z) Centroid(ReadOnlySpan<int> row)
    {
        double x = 0, y = 0, z = 0;
        foreach (var v in row)
        {
            x += Xa[v];
            y += Ya[v];
            z += Za[v];
        }

        return (x / row.Length, y / row.Length, z / row.Length);
    }
}

public sealed class CrackPlane
{
    private readonly int _n;
    private readonly int _ic;
    public int Row { get; set; }

    public CrackPlane(int n)
    {
        _n = n;
        _ic = n / 2;
    }

    public int Steps => _n + 1;

    /// <summary>Duplicates the nodes of the current front row; returns the number of edited elements.</summary>
    public int Step(MeshState st)
    {
        var mesh = st.Mesh;
        var xPlane = (double)_ic / _n;
        var edits = new SortedDictionary<int, int[]>();
        for (var k = 0; k <= _n; k++)
        {
            var v = Kuhn.Node(_n, _ic, Row, k);
            var pairs = mesh.TransposePairs(v).ToArray();
            int dup = -1;
            for (var q = 0; q < pairs.Length; q += 2)
            {
                var e = pairs[q];
                var p = pairs[q + 1];
                var row = edits.TryGetValue(e, out var pending) ? pending : mesh.Row(e).ToArray();
                if (st.Centroid(row).X <= xPlane) continue;
                if (dup < 0) dup = st.AddNode(st.X(v, 0), st.X(v, 1), st.X(v, 2));
                row[p] = dup;
                edits[e] = row;
            }
        }

        foreach (var (e, row) in edits) st.Replace(e, row);
        Row++;
        return edits.Count;
    }
}

public sealed class Bisection
{
    private readonly List<int> _patch = new();
    public long Bisections;

    private static (int A, int B) RefEdge(MeshState st, int e)
    {
        var row = st.Mesh.Row(e);
        return (row[0], row[st.Tag[e]]);
    }

    private static bool Same((int A, int B) x, int a, int b) => (x.A == a && x.B == b) || (x.A == b && x.B == a);

    private void ElementsWithEdge(MeshState st, int a, int b, List<int> result)
    {
        result.Clear();
        var ra = st.Mesh.TransposePairs(a);
        var rb = st.Mesh.TransposePairs(b);
        if (rb.Length < ra.Length)
        {
            var t = ra;
            ra = rb;
            rb = t;
        }

        for (var q = 0; q < ra.Length; q += 2)
            if (PairRow.IndexOf(rb, ra[q]) >= 0)
                result.Add(ra[q]);
    }

    /// <summary>Bisects element e and, recursively, whatever its refinement edge requires.</summary>
    public void Bisect(MeshState st, int e, int depth = 0)
    {
        if (depth > 256) throw new InvalidOperationException("bisection recursion did not terminate");
        var (a, b) = RefEdge(st, e);
        var patch = new List<int>();
        while (true)
        {
            ElementsWithEdge(st, a, b, patch);
            var incompatible = -1;
            foreach (var g in patch)
                if (!Same(RefEdge(st, g), a, b))
                {
                    incompatible = g;
                    break;
                }

            if (incompatible < 0) break;
            Bisect(st, incompatible, depth + 1);
        }

        var z = st.AddNode(0.5 * (st.X(a, 0) + st.X(b, 0)), 0.5 * (st.X(a, 1) + st.X(b, 1)),
            0.5 * (st.X(a, 2) + st.X(b, 2)));
        foreach (var g in patch) Split(st, g, z);
    }

    private void Split(MeshState st, int g, int z)
    {
        var row = st.Mesh.Row(g).ToArray();
        var d = st.Tag[g];
        var child1 = (int[])row.Clone();
        child1[d] = z;
        var child2 = new int[4];
        for (var p = 0; p < d; p++) child2[p] = row[p + 1];
        child2[d] = z;
        for (var p = d + 1; p < 4; p++) child2[p] = row[p];
        var nd = d > 1 ? d - 1 : 3;
        var level = st.Level[g] + 1;
        st.Replace(g, child1);
        st.Tag.Set(g, nd);
        st.Level.Set(g, level);
        st.Add(child2, nd, level);
        Bisections++;
    }

    /// <summary>Bisects every listed element whose level is still below maxLevel; returns elements bisected.</summary>
    public long Refine(MeshState st, IReadOnlyList<int> marked, int maxLevel)
    {
        var before = Bisections;
        foreach (var e in marked)
            if (st.Mesh.Row(e).Length > 0 && st.Level[e] < maxLevel)
                Bisect(st, e);
        return Bisections - before;
    }

    /// <summary>Elements with centroid inside the ball (scan; not part of the timed region).</summary>
    public static List<int> InBall(MeshState st, double cx, double cy, double cz, double r)
    {
        var list = new List<int>();
        var r2 = r * r;
        for (var e = 0; e < st.Mesh.ElementCount; e++)
        {
            var row = st.Mesh.Row(e);
            if (row.Length == 0) continue;
            var (x, y, z) = st.Centroid(row);
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz) <= r2) list.Add(e);
        }

        return list;
    }
}

public static class Erosion
{
    /// <summary>Deletes the listed elements; returns the number deleted.</summary>
    public static int Step(MeshState st, IReadOnlyList<int> victims)
    {
        var n = 0;
        foreach (var e in victims)
            if (st.Mesh.Row(e).Length > 0)
            {
                st.Remove(e);
                n++;
            }

        return n;
    }
}

public sealed class RandomEdits
{
    private readonly Random _rng;

    public RandomEdits(int seed) => _rng = new Random(seed);

    /// <summary>
    ///     One randomized edit: relabel a node to a fresh node (duplication), relabel to an
    ///     existing node, permute a row, delete an element, or create an element.
    /// </summary>
    public void One(MeshState st, bool allowNonManifold)
    {
        var mesh = st.Mesh;
        var kind = _rng.Next(allowNonManifold ? 5 : 2);
        var e = PickAlive(mesh);
        if (e < 0) kind = 4;
        switch (kind)
        {
            case 0:
            {
                // node duplication for one element (a hole opens at that node)
                var row = mesh.Row(e).ToArray();
                var p = _rng.Next(row.Length);
                var v = row[p];
                row[p] = st.AddNode(st.X(v, 0), st.X(v, 1), st.X(v, 2));
                st.Replace(e, row);
                break;
            }
            case 1:
            {
                // permutation of the row (orientation change, positions move)
                var row = mesh.Row(e).ToArray();
                var i = _rng.Next(row.Length);
                var j = _rng.Next(row.Length);
                (row[i], row[j]) = (row[j], row[i]);
                st.Replace(e, row);
                break;
            }
            case 2:
            {
                // relabel to an existing node not in the row (non-manifold configurations)
                var row = mesh.Row(e).ToArray();
                var p = _rng.Next(row.Length);
                var v = _rng.Next(mesh.NodeCount);
                if (Array.IndexOf(row, v) >= 0) break;
                row[p] = v;
                st.Replace(e, row);
                break;
            }
            case 3:
                st.Remove(e);
                break;
            default:
            {
                // new element on distinct random nodes
                var nN = mesh.NodeCount;
                if (nN < 4) break;
                var row = new int[4];
                for (var p = 0; p < 4; p++)
                {
                    int v;
                    do v = _rng.Next(nN);
                    while (Array.IndexOf(row, v, 0, p) >= 0);
                    row[p] = v;
                }

                st.Add(row, 3, 0);
                break;
            }
        }
    }

    /// <summary>Node duplication for one random element (the crack-insertion edit).</summary>
    public void Duplicate(MeshState st)
    {
        var mesh = st.Mesh;
        var e = PickAlive(mesh);
        if (e < 0) return;
        var row = mesh.Row(e).ToArray();
        var p = _rng.Next(row.Length);
        var v = row[p];
        row[p] = st.AddNode(st.X(v, 0), st.X(v, 1), st.X(v, 2));
        st.Replace(e, row);
    }

    private int PickAlive(IncrementalMesh mesh)
    {
        for (var attempt = 0; attempt < 64; attempt++)
        {
            var e = _rng.Next(mesh.ElementCount);
            if (mesh.Row(e).Length > 0) return e;
        }

        return -1;
    }
}

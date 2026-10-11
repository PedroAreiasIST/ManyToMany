// CrackMesh.cs — the mesh with its versioned node and element data: coordinates, lineage, canonical labels, crack masks.
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
// Every node carries, as versioned attributes of the IncrementalMesh, its coordinates,
// its parents (the edge it was created on, or itself twice for lattice nodes and split
// copies), the edge parameter, a crack mask (bit c set when the node lies on the surface
// of crack c) and a canonical label. The canonical label identifies the geometric point:
// a split copy inherits the label of the node it was split from, and a node created on
// an edge takes the label of the first node created on the same canonical edge at the
// same parameter, so that two copies of a cut edge (one on each side of an earlier crack)
// give coincident nodes with one label. Crack facets are then boundary facets whose
// canonical triangles occur twice (CrackCount), without reference to any level set.

using Step2;

namespace Step3;

public sealed class CrackMesh : ICutObserver
{
    public readonly IncrementalMesh M;
    public readonly Attribute<double> X, Y, Z, Weight;
    public readonly Attribute<int> Parent1, Parent2, Canon, Mask, Origin, Tag, Level;
    public readonly Attribute<byte> Bisectable;
    private readonly List<Attribute<double>> _fields = new();
    private readonly Dictionary<(int, int, long), int> _canonOfEdgeNode = new();

    /// <summary>Roots of the edges being cut (filled by the inserter before Cutter.Cut; midpoints when absent).</summary>
    public Dictionary<long, double> PendingRoots = new(EdgeKeyComparer.Instance);

    /// <summary>Crack mask given to the nodes created by the cut in progress (0 for refinement).</summary>
    public int PendingMask;

    public CrackMesh(FlatStructures flat, ReadOnlySpan<double> xyz, bool persistent = true, int[]? tags = null)
    {
        M = new IncrementalMesh(flat, ElementSpec.Tetrahedron, persistent);
        X = M.AddAttribute(AttributeKind.Node, 0.0);
        Y = M.AddAttribute(AttributeKind.Node, 0.0);
        Z = M.AddAttribute(AttributeKind.Node, 0.0);
        Weight = M.AddAttribute(AttributeKind.Node, 0.0);
        Parent1 = M.AddAttribute(AttributeKind.Node, -1);
        Parent2 = M.AddAttribute(AttributeKind.Node, -1);
        Canon = M.AddAttribute(AttributeKind.Node, -1);
        Mask = M.AddAttribute(AttributeKind.Node, 0);
        Origin = M.AddAttribute(AttributeKind.Element, -1);
        Tag = M.AddAttribute(AttributeKind.Element, 3);
        Level = M.AddAttribute(AttributeKind.Element, 0);
        Bisectable = M.AddAttribute(AttributeKind.Element, (byte)1);
        for (var e = 0; e < M.ElementCount; e++)
        {
            Origin.Set(e, e);
            if (tags != null) Tag.Set(e, tags[e]);
        }

        for (var v = 0; v < M.NodeCount; v++)
        {
            X.Set(v, xyz[3 * v]);
            Y.Set(v, xyz[3 * v + 1]);
            Z.Set(v, xyz[3 * v + 2]);
            Parent1.Set(v, v);
            Parent2.Set(v, v);
            Canon.Set(v, v);
        }
    }

    public int NodeCount => M.NodeCount;
    public int ElementCount => M.ElementCount;
    public V3 P(int v) => new(X[v], Y[v], Z[v]);

    /// <summary>Registers a nodal field that follows the mesh (interpolated on cut edges, copied on splits).</summary>
    public void Follow(Attribute<double> field) => _fields.Add(field);

    public IReadOnlyList<Attribute<double>> Fields => _fields;

    /// <summary>Flat coordinates (3 per node) of the current state.</summary>
    public double[] Coordinates()
    {
        var x = new double[3 * M.NodeCount];
        for (var v = 0; v < M.NodeCount; v++)
        {
            x[3 * v] = X[v];
            x[3 * v + 1] = Y[v];
            x[3 * v + 2] = Z[v];
        }

        return x;
    }

    /// <summary>Flat coordinates of a committed version.</summary>
    public static double[] Coordinates(Snapshot s, CrackMesh cm)
    {
        var x = new double[3 * s.NodeCount];
        var vx = s.View(cm.X);
        var vy = s.View(cm.Y);
        var vz = s.View(cm.Z);
        for (var v = 0; v < s.NodeCount; v++)
        {
            x[3 * v] = vx[v];
            x[3 * v + 1] = vy[v];
            x[3 * v + 2] = vz[v];
        }

        return x;
    }

    private static long Quantize(double t) => (long)Math.Round(t * 1e12);

    /// <summary>Creates a node on the edge (a, b) at parameter t, with the canonical label of that geometric point.</summary>
    public int NodeOnEdge(int a, int b, double t, int mask)
    {
        var v = M.AddNode();
        X.Set(v, X[a] + t * (X[b] - X[a]));
        Y.Set(v, Y[a] + t * (Y[b] - Y[a]));
        Z.Set(v, Z[a] + t * (Z[b] - Z[a]));
        Parent1.Set(v, a);
        Parent2.Set(v, b);
        Weight.Set(v, t);
        Mask.Set(v, mask | (Mask[a] & Mask[b])); // a node on an edge of a crack face lies on that crack
        int ca = Canon[a], cb = Canon[b];
        var key = ca < cb ? (ca, cb, Quantize(t)) : (cb, ca, Quantize(1 - t));
        if (_canonOfEdgeNode.TryGetValue(key, out var c) && c < v && Canon[c] == c && c < M.NodeCount)
            Canon.Set(v, c);
        else
        {
            Canon.Set(v, v);
            _canonOfEdgeNode[key] = v;
        }

        foreach (var u in _fields) u.Set(v, (1 - t) * u[a] + t * u[b]);
        return v;
    }

    /// <summary>Creates a copy of node v (same point, same canonical label and mask).</summary>
    public int Copy(int v)
    {
        var c = M.AddNode();
        X.Set(c, X[v]);
        Y.Set(c, Y[v]);
        Z.Set(c, Z[v]);
        Parent1.Set(c, v);
        Parent2.Set(c, v);
        Weight.Set(c, 0);
        Mask.Set(c, Mask[v]);
        Canon.Set(c, Canon[v]);
        foreach (var u in _fields) u.Set(c, u[v]);
        return c;
    }

    public Snapshot Commit() => M.Commit();

    /// <summary>Rolls back to a committed version and forgets canonical edge labels of nodes that no longer exist.</summary>
    public void Rollback(Snapshot s)
    {
        M.Rollback(s);
        var n = M.NodeCount;
        var stale = new List<(int, int, long)>();
        foreach (var kv in _canonOfEdgeNode)
            if (kv.Value >= n)
                stale.Add(kv.Key);
        foreach (var k in stale) _canonOfEdgeNode.Remove(k);
    }

    // ---------------------------------------------------------------- ICutObserver

    void ICutObserver.NodeCreated(int node, int a, int b)
    {
        // Cutter has already called AddNode; set the data of the new node here
        var t = PendingRoots.TryGetValue(Cutter.Key(a, b), out var r) ? r : 0.5;
        X.Set(node, X[a] + t * (X[b] - X[a]));
        Y.Set(node, Y[a] + t * (Y[b] - Y[a]));
        Z.Set(node, Z[a] + t * (Z[b] - Z[a]));
        Parent1.Set(node, a);
        Parent2.Set(node, b);
        Weight.Set(node, t);
        Mask.Set(node, PendingMask | (Mask[a] & Mask[b]));
        int ca = Canon[a], cb = Canon[b];
        var key = ca < cb ? (ca, cb, Quantize(t)) : (cb, ca, Quantize(1 - t));
        if (_canonOfEdgeNode.TryGetValue(key, out var c) && c < node && c < M.NodeCount && Canon[c] == c)
            Canon.Set(node, c);
        else
        {
            Canon.Set(node, node);
            _canonOfEdgeNode[key] = node;
        }

        foreach (var u in _fields) u.Set(node, (1 - t) * u[a] + t * u[b]);
    }

    void ICutObserver.ElementSplit(int parent, int child)
    {
        Origin.Set(child, Origin[parent]);
        Bisectable.Set(child, 0); // template children are not Maubach-tagged
        Level.Set(child, Level[parent]);
    }
}

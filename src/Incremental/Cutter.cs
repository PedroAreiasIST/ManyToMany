// Cutter.cs — local Steiner-point-free edge cutting on an IncrementalMesh.
//
// A cut takes a list of marked node pairs. Each distinct pair receives one new node,
// labelled above every existing node in the order of first appearance in the list
// (the order of the batch cutter). The elements containing a marked edge are the
// shell of that edge: the common elements of the transpose rows of its two nodes.
// Only these elements change: each is replaced by the children of its template
// (Templates), the first child keeping the parent's index and the others appended,
// and point elements at an endpoint of a marked edge are copied to its new node.
// Every derived structure follows through the delta rules of IncrementalMesh, so a
// cut costs time proportional to the shells, not to the mesh. Because a template
// depends only on the row, the midpoint labels and their relative order, the cut
// produces the same rows as the batch cutter applied to the whole mesh.

using Numerical;

namespace Step2;

/// <summary>Callbacks for data that follow the cut (coordinates, fields, history).</summary>
public interface ICutObserver
{
    /// <summary>A node was created on edge (a, b).</summary>
    void NodeCreated(int node, int a, int b);

    /// <summary>Element <paramref name="child"/> was produced by splitting <paramref name="parent"/>
    ///     (the first child is the parent index itself).</summary>
    void ElementSplit(int parent, int child);
}

/// <summary>What a cut did; enough to revert it locally (Uncut) or to transfer data.</summary>
public sealed class CutRecord
{
    /// <summary>Distinct marked edges (a &lt; b) in label order of their new nodes.</summary>
    public required (int A, int B, int Mid)[] Edges;

    /// <summary>Split elements: parent index, parent row, children (Children[0] = parent index).</summary>
    public required (int Parent, int[] Row, int[] Children)[] Splits;

    /// <summary>Point elements copied to new nodes: (new element, source point element).</summary>
    public required (int Element, int Source)[] Points;

    /// <summary>Node count before the cut (every new node is labelled at or above it).</summary>
    public required int FirstNewNode;

    public int ElementsBefore;
    public int ShellElements;
}

/// <summary>
///     Hash of an edge key (a &lt;&lt; 32 | b). The default hash of a long, the exclusive or of its
///     halves, takes few distinct values on mesh labels (1205 for the 256 000 edges of a Kuhn
///     box with n = 40), so lookups degrade to linear scans; Fibonacci hashing spreads them.
/// </summary>
public sealed class EdgeKeyComparer : IEqualityComparer<long>
{
    public static readonly EdgeKeyComparer Instance = new();
    public bool Equals(long x, long y) => x == y;
    public int GetHashCode(long x) => (int)(((ulong)x * 0x9E3779B97F4A7C15UL) >> 32);
}

public static partial class Cutter
{
    /// <summary>Elements (any length >= 2) containing both a and b, ascending.</summary>
    public static void Shell(IncrementalMesh m, int a, int b, List<int> result)
    {
        ReadOnlySpan<int> r0 = m.TransposePairs(a), r1 = m.TransposePairs(b);
        int i = 0, j = 0;
        while (i < r0.Length && j < r1.Length)
        {
            int x = r0[i], y = r1[j];
            if (x == y)
            {
                result.Add(x);
                i += 2;
                j += 2;
            }
            else if (x < y) i += 2;
            else j += 2;
        }
    }

    /// <summary>
    ///     Cuts the marked edges. Pairs are unordered; duplicates are ignored. A pair that is
    ///     not an edge of the mesh still receives a node, as in the batch cutter.
    /// </summary>
    public static CutRecord Cut(IncrementalMesh m, IReadOnlyList<(int, int)> marked, ICutObserver? obs = null)
    {
        var first = m.NodeCount;
        var e0 = m.ElementCount;
        var mid = new Dictionary<long, int>(marked.Count, EdgeKeyComparer.Instance);
        var edges = new List<(int, int, int)>(marked.Count);
        foreach (var (p, q) in marked)
        {
            if (p == q) throw new ArgumentException($"degenerate pair ({p},{p})");
            int a = Math.Min(p, q), b = Math.Max(p, q);
            var key = Key(a, b);
            if (mid.ContainsKey(key)) continue;
            var c = m.AddNode();
            mid[key] = c;
            edges.Add((a, b, c));
            obs?.NodeCreated(c, a, b);
        }

        // ---- point elements at the endpoints, copied to the new node (batch order: a's, then b's)
        var points = new List<(int, int)>();
        var scratch = new List<int>();
        foreach (var (a, b, c) in edges)
        {
            foreach (var v in (ReadOnlySpan<int>)[a, b])
            {
                scratch.Clear();
                var tr = m.TransposePairs(v);
                for (var q = 0; q < tr.Length; q += 2)
                    if (m.Row(tr[q]).Length == 1)
                        scratch.Add(tr[q]);
                foreach (var pe in scratch)
                {
                    var ne = m.AddElement([c]);
                    points.Add((ne, pe));
                    obs?.ElementSplit(pe, ne);
                }
            }
        }

        // ---- shells of the marked edges (elements of length >= 2), each element once
        var affected = new HashSet<int>();
        var shell = new List<int>();
        foreach (var (a, b, _) in edges)
        {
            shell.Clear();
            Shell(m, a, b, shell);
            foreach (var g in shell)
                if (m.Row(g).Length >= 2)
                    affected.Add(g);
        }

        var order = affected.ToArray();
        Array.Sort(order);
        var splits = new List<(int, int[], int[])>(order.Length);
        var children = new List<int[]>(8);
        Span<int> mid6 = stackalloc int[6];
        foreach (var e in order)
        {
            var row = m.Row(e).ToArray();
            children.Clear();
            Children(row, mid, mid6, children);
            var idx = new int[children.Count];
            idx[0] = e;
            m.ReplaceElement(e, children[0]);
            obs?.ElementSplit(e, e);
            for (var t = 1; t < children.Count; t++)
            {
                idx[t] = m.AddElement(children[t]);
                obs?.ElementSplit(e, idx[t]);
            }

            splits.Add((e, row, idx));
        }

        return new CutRecord
        {
            Edges = edges.ToArray(), Splits = splits.ToArray(), Points = points.ToArray(), FirstNewNode = first,
            ElementsBefore = e0, ShellElements = order.Length
        };
    }

    /// <summary>
    ///     Reverts the cut of the edges in <paramref name="revert"/> (keys from <see cref="Key"/>)
    ///     locally: every parent of <paramref name="rec"/> that contains a reverted edge is split
    ///     again with the remaining marks, reusing its children's indices; children no longer
    ///     needed and the point elements copied to reverted nodes are deleted, and the reverted
    ///     nodes are left without elements (labels are never reused). The surviving new nodes keep
    ///     their labels, whose order is that of the original cut, so the result equals the batch
    ///     cut of the parents with the remaining marks under an order-preserving relabeling.
    ///     The children must not have been edited since the cut; this is checked. Returns the
    ///     record of the remaining cut, which can itself be reverted.
    /// </summary>
    public static CutRecord Uncut(IncrementalMesh m, CutRecord rec, IReadOnlySet<long> revert,
        ICutObserver? obs = null)
    {
        var all = new Dictionary<long, int>(rec.Edges.Length, EdgeKeyComparer.Instance);
        var keep = new Dictionary<long, int>(rec.Edges.Length, EdgeKeyComparer.Instance);
        foreach (var (a, b, c) in rec.Edges)
        {
            all[Key(a, b)] = c;
            if (!revert.Contains(Key(a, b))) keep[Key(a, b)] = c;
        }

        var revertedNodes = new HashSet<int>(rec.Edges.Where(t => !keep.ContainsKey(Key(t.A, t.B))).Select(t => t.Mid));
        Span<int> mid6 = stackalloc int[6];
        var oldChildren = new List<int[]>(8);
        var newChildren = new List<int[]>(8);
        var splits = new List<(int, int[], int[])>(rec.Splits.Length);
        foreach (var (parent, row, idx) in rec.Splits)
        {
            // children as produced by the cut: they must still be in place
            oldChildren.Clear();
            Children(row, all, mid6, oldChildren);
            if (oldChildren.Count != idx.Length)
                throw new InvalidOperationException($"split of {parent}: record does not match the templates");
            for (var t = 0; t < idx.Length; t++)
                if (!m.Row(idx[t]).SequenceEqual(oldChildren[t]))
                    throw new InvalidOperationException($"child {idx[t]} of {parent} was edited after the cut");

            var touched = false;
            for (var p = 0; p < row.Length && !touched; p++)
            for (var q = p + 1; q < row.Length && !touched; q++)
                touched = revert.Contains(Key(row[p], row[q])) && all.ContainsKey(Key(row[p], row[q]));
            if (!touched)
            {
                splits.Add((parent, row, idx));
                continue;
            }

            newChildren.Clear();
            Children(row, keep, mid6, newChildren);
            var nidx = new int[newChildren.Count];
            for (var t = 0; t < newChildren.Count; t++)
            {
                if (t < idx.Length)
                {
                    nidx[t] = idx[t];
                    m.ReplaceElement(idx[t], newChildren[t]);
                }
                else nidx[t] = m.AddElement(newChildren[t]);

                obs?.ElementSplit(parent, nidx[t]);
            }

            for (var t = newChildren.Count; t < idx.Length; t++) m.RemoveElement(idx[t]);
            if (newChildren.Count > 1 || !newChildren[0].SequenceEqual(row)) splits.Add((parent, row, nidx));
        }

        var points = new List<(int, int)>(rec.Points.Length);
        foreach (var (pe, src) in rec.Points)
        {
            var r = m.Row(pe);
            if (r.Length == 1 && revertedNodes.Contains(r[0])) m.RemoveElement(pe);
            else points.Add((pe, src));
        }

        foreach (var v in revertedNodes)
            if (m.TransposePairs(v).Length != 0)
                throw new InvalidOperationException($"reverted node {v} still has elements");

        return new CutRecord
        {
            Edges = rec.Edges.Where(t => keep.ContainsKey(Key(t.A, t.B))).ToArray(), Splits = splits.ToArray(),
            Points = points.ToArray(), FirstNewNode = rec.FirstNewNode, ElementsBefore = rec.ElementsBefore,
            ShellElements = splits.Count
        };
    }

    /// <summary>Children of one row under the marked-edge map (new node labels).</summary>
    public static void Children(int[] row, Dictionary<long, int> mid, Span<int> mid6, List<int[]> children)
    {
        switch (row.Length)
        {
            case 2:
                SubdivisionTemplates.SplitBar(row[0], row[1], Lookup(mid, row[0], row[1]), children);
                break;
            case 3:
                SubdivisionTemplates.SplitTriangle(row[0], row[1], row[2], Lookup(mid, row[0], row[1]),
                    Lookup(mid, row[1], row[2]), Lookup(mid, row[2], row[0]), children);
                break;
            case 4:
                for (var ie = 0; ie < 6; ie++)
                {
                    var x = Lookup(mid, row[SubdivisionTemplates.EdgeNodes[ie][0]], row[SubdivisionTemplates.EdgeNodes[ie][1]]);
                    mid6[ie] = x >= 0 ? x : 0;
                }

                SubdivisionTemplates.SplitTet(row, mid6, children);
                break;
            default:
                throw new InvalidOperationException($"row of length {row.Length} cannot be split");
        }
    }

    public static long Key(int a, int b) => a < b ? ((long)a << 32) | (uint)b : ((long)b << 32) | (uint)a;

    private static int Lookup(Dictionary<long, int> mid, int a, int b) =>
        mid.TryGetValue(Key(a, b), out var c) ? c : -1;
}

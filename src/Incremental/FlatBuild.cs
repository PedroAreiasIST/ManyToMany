// FlatBuild.cs — from-scratch evaluation of every derived structure into flat arrays.
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
// This is the tuned rebuild: counting-sort transpose with node positions, element
// positions, the node-node pattern with pair counts, assembly destinations, facet
// adjacency, the boundary-facet list and a greedy colouring, each in one or two
// passes over the incidence. It serves three purposes: the initial state of an
// IncrementalMesh, the strongest rebuild baseline in the measurements, and (through
// a separate naive oracle in the test harness) the reference for verification.

using System.Runtime.CompilerServices;

namespace Step2;

/// <summary>Derived structures in flat (non-persistent) arrays.</summary>
public sealed class FlatStructures
{
    public int K, KK, Nf, NE, NN;
    public int[][] A = [];
    public int[] TPtr = [];      // NN + 1
    public int[] TPairs = [];    // 2 * nnz: (element, local position) sorted by element per node
    public int[] EP = [];        // NE * K: index of e in T[A[e][p]]
    public int[] PPtr = [];      // NN + 1
    public int[] PCol = [];      // pattern columns, discovery order per row
    public int[] PCnt = [];      // pair counts
    public int[] Dest = [];      // NE * K * K: slot of A[e][j] in pattern row A[e][i]
    public int[] FA = [];        // NE * Nf: smallest other element sharing facet f, or -1
    public int[] Boundary = [];  // facet ids e * Nf + f with FA = -1 (alive elements)
    public int[] Color = [];     // NE
    public int ColorCount;
}

/// <summary>Search helpers over transpose rows stored as (element, position) pairs.</summary>
public static class PairRow
{
    /// <summary>Pair index of element <paramref name="e"/>, or bitwise complement of the insertion point.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int IndexOf(ReadOnlySpan<int> pairs, int e)
    {
        int lo = 0, hi = (pairs.Length >> 1) - 1;
        while (lo <= hi)
        {
            var mid = (lo + hi) >> 1;
            var x = pairs[mid << 1];
            if (x == e) return mid;
            if (x < e) lo = mid + 1;
            else hi = mid - 1;
        }

        return ~lo;
    }

    /// <summary>
    ///     Smallest element other than <paramref name="self"/> present in all rows; -1 if none.
    ///     Rows are pair rows sorted by element.
    /// </summary>
    public static int MinOtherCommon(ReadOnlySpan<int> r0, ReadOnlySpan<int> r1, ReadOnlySpan<int> r2, int self)
    {
        // three-way merge of rows sorted by element (stride 2); first common element != self
        int i = 0, j = 0, k = 0;
        while (i < r0.Length && j < r1.Length && k < r2.Length)
        {
            int a = r0[i], b = r1[j], c = r2[k];
            if (a == b && b == c)
            {
                if (a != self) return a;
                i += 2;
                j += 2;
                k += 2;
                continue;
            }

            var m = Math.Max(a, Math.Max(b, c));
            if (a < m) i += 2;
            if (b < m) j += 2;
            if (c < m) k += 2;
        }

        return -1;
    }

    /// <summary>All elements other than self present in the three rows (three-way merge).</summary>
    public static void CommonOthers(ReadOnlySpan<int> r0, ReadOnlySpan<int> r1, ReadOnlySpan<int> r2, int self,
        List<int> result)
    {
        int i = 0, j = 0, k = 0;
        while (i < r0.Length && j < r1.Length && k < r2.Length)
        {
            int a = r0[i], b = r1[j], c = r2[k];
            if (a == b && b == c)
            {
                if (a != self) result.Add(a);
                i += 2;
                j += 2;
                k += 2;
                continue;
            }

            var m = Math.Max(a, Math.Max(b, c));
            if (a < m) i += 2;
            if (b < m) j += 2;
            if (c < m) k += 2;
        }
    }

    /// <summary>Two-row variant (triangles: facets are edges).</summary>
    public static int MinOtherCommon(ReadOnlySpan<int> r0, ReadOnlySpan<int> r1, int self)
    {
        if (r1.Length < r0.Length) { var t = r0; r0 = r1; r1 = t; }
        for (var q = 0; q < r0.Length; q += 2)
        {
            var x = r0[q];
            if (x == self) continue;
            if (IndexOf(r1, x) >= 0) return x;
        }

        return -1;
    }
}

public static class FlatBuild
{
    private const int ParallelThreshold = 4096;

    public static FlatStructures Build(ElementSpec spec, int nNodes, IReadOnlyList<int[]> rows, bool parallel) =>
        Build(spec, nNodes, rows, parallel, facetsAndColours: true);

    /// <summary>
    ///     The initial state of an incremental mesh: every phase threaded except the colouring, which is
    ///     the greedy colouring in element order, the rule the edits maintain (the threaded block colouring
    ///     uses more colours, and the colour repairs of later edits cost more on it).
    /// </summary>
    public static FlatStructures Initial(ElementSpec spec, int nNodes, IReadOnlyList<int[]> rows) =>
        Build(spec, nNodes, rows, parallel: true, facetsAndColours: true, greedyColouring: true);

    /// <summary>
    ///     From-scratch evaluation. With <paramref name="facetsAndColours"/> false only the
    ///     structures MM3 also derives are built: transpose with node positions, element
    ///     positions, the node-node pattern and the destinations.
    /// </summary>
    public static FlatStructures Build(ElementSpec spec, int nNodes, IReadOnlyList<int[]> rows, bool parallel,
        bool facetsAndColours, bool greedyColouring = false)
    {
        int K = spec.K, KK = K * K, Nf = spec.Nf;
        var nE = rows.Count;
        var A = rows as int[][] ?? rows.ToArray();
        var s = new FlatStructures { K = K, KK = KK, Nf = Nf, NE = nE, NN = nNodes, A = A };
        var par = parallel && nE >= ParallelThreshold;

        // ---- 1. transpose by counting sort; ascending element order gives sorted rows
        int[] tptr, tpairs;
        if (par) (tptr, tpairs) = TransposeChunked(A, nE, nNodes);
        else
        {
            tptr = new int[nNodes + 1];
            for (var e = 0; e < nE; e++)
                foreach (var v in A[e])
                    tptr[v + 1]++;
            for (var v = 0; v < nNodes; v++) tptr[v + 1] += tptr[v];
            tpairs = new int[2 * tptr[nNodes]];
            var fill = new int[nNodes];
            Array.Copy(tptr, fill, nNodes);
            for (var e = 0; e < nE; e++)
            {
                var row = A[e];
                for (var p = 0; p < row.Length; p++)
                {
                    var q = fill[row[p]]++;
                    tpairs[2 * q] = e;
                    tpairs[2 * q + 1] = p;
                }
            }
        }

        s.TPtr = tptr;
        s.TPairs = tpairs;

        // ---- 2. element positions
        var ep = GC.AllocateUninitializedArray<int>(nE * K);
        Array.Fill(ep, -1);
        void Ep(int v)
        {
            int b = tptr[v], n = tptr[v + 1] - b;
            for (var q = 0; q < n; q++)
                ep[tpairs[2 * (b + q)] * K + tpairs[2 * (b + q) + 1]] = q;
        }

        if (par) Parallel.For(0, nNodes, Ep);
        else
            for (var v = 0; v < nNodes; v++)
                Ep(v);
        s.EP = ep;

        // ---- 3. node-node pattern with counts, and destinations (two passes)
        var plen = new int[nNodes];
        if (par)
            Parallel.For(0, nNodes, () => new int[nNodes], (v, _, stamp) =>
            {
                plen[v] = CountRow(v, tptr, tpairs, A, stamp);
                return stamp;
            }, _ => { });
        else
        {
            var stamp = new int[nNodes];
            for (var v = 0; v < nNodes; v++) plen[v] = CountRow(v, tptr, tpairs, A, stamp);
        }

        var pptr = new int[nNodes + 1];
        for (var v = 0; v < nNodes; v++) pptr[v + 1] = pptr[v] + plen[v];
        var pcol = new int[pptr[nNodes]];
        var pcnt = new int[pptr[nNodes]];
        var dest = GC.AllocateUninitializedArray<int>(nE * KK);
        Array.Fill(dest, -1);
        if (par)
            Parallel.For(0, nNodes, () => (new int[nNodes], new int[nNodes]), (v, _, st) =>
            {
                FillRow(v, K, KK, tptr, tpairs, A, pptr, pcol, pcnt, dest, st.Item1, st.Item2);
                return st;
            }, _ => { });
        else
        {
            var stamp = new int[nNodes];
            var slot = new int[nNodes];
            for (var v = 0; v < nNodes; v++)
                FillRow(v, K, KK, tptr, tpairs, A, pptr, pcol, pcnt, dest, stamp, slot);
        }

        s.PPtr = pptr;
        s.PCol = pcol;
        s.PCnt = pcnt;
        s.Dest = dest;

        if (!facetsAndColours) return s;

        // ---- 4. facet adjacency
        var fa = GC.AllocateUninitializedArray<int>(nE * Nf);
        void Fa(int e)
        {
            var row = A[e];
            if (row.Length == 0)
            {
                for (var f = 0; f < Nf; f++) fa[e * Nf + f] = -1;
                return;
            }

            if (spec.Mixed)
            {
                var nf = spec.FacetCount(row.Length);
                for (var f = 0; f < Nf; f++)
                    fa[e * Nf + f] = f < nf ? MinOtherSameLength(spec.FacetLocalByLength[row.Length][f], row, e) : -1;
                return;
            }

            for (var f = 0; f < Nf; f++)
            {
                var loc = spec.FacetLocal[f];
                var r0 = Row(tptr, tpairs, row[loc[0]]);
                var r1 = Row(tptr, tpairs, row[loc[1]]);
                fa[e * Nf + f] = loc.Length == 3
                    ? PairRow.MinOtherCommon(r0, r1, Row(tptr, tpairs, row[loc[2]]), e)
                    : PairRow.MinOtherCommon(r0, r1, e);
            }
        }

        // smallest other element of the same row length containing every facet node
        int MinOtherSameLength(int[] loc, int[] row, int e)
        {
            var r0 = Row(tptr, tpairs, row[loc[0]]);
            for (var q = 0; q < r0.Length; q += 2)
            {
                var g = r0[q];
                if (g == e || A[g].Length != row.Length) continue;
                var all = true;
                for (var t = 1; t < loc.Length && all; t++)
                    all = PairRow.IndexOf(Row(tptr, tpairs, row[loc[t]]), g) >= 0;
                if (all) return g;
            }

            return -1;
        }

        if (par) Parallel.For(0, nE, Fa);
        else
            for (var e = 0; e < nE; e++)
                Fa(e);
        s.FA = fa;

        // ---- 5. boundary facets, in element order (per block, concatenated)
        void Boundary(int lo, int hi, List<int> list)
        {
            for (var e = lo; e < hi; e++)
            {
                if (A[e].Length == 0) continue;
                var nfe = spec.Mixed ? spec.FacetCount(A[e].Length) : Nf;
                for (var f = 0; f < nfe; f++)
                    if (fa[e * Nf + f] < 0)
                        list.Add(e * Nf + f);
            }
        }

        if (par)
        {
            var nb = Math.Max(1, Math.Min(64, nE / 4096));
            var parts = new List<int>[nb];
            Parallel.For(0, nb, b =>
            {
                parts[b] = new List<int>();
                Boundary((int)((long)nE * b / nb), (int)((long)nE * (b + 1) / nb), parts[b]);
            });
            var boundary = new int[parts.Sum(l => l.Count)];
            var at = 0;
            foreach (var l in parts)
            {
                l.CopyTo(boundary, at);
                at += l.Count;
            }

            s.Boundary = boundary;
        }
        else
        {
            var boundary = new List<int>();
            Boundary(0, nE, boundary);
            s.Boundary = boundary.ToArray();
        }

        // ---- 6. colouring over the vertex-adjacency graph
        if (par && !greedyColouring)
        {
            (s.Color, s.ColorCount) = ColourParallel(A, nE, tptr, tpairs);
            return s;
        }

        var color = new int[nE];
        Array.Fill(color, -1);
        var used = new int[64];
        var colorCount = 0;
        for (var e = 0; e < nE; e++)
        {
            var row = A[e];
            if (row.Length == 0) continue;
            var c = SmallestFree(e, row, tptr, tpairs, color, ref used);
            color[e] = c;
            if (c + 1 > colorCount) colorCount = c + 1;
        }

        s.Color = color;
        s.ColorCount = colorCount;
        return s;
    }

    /// <summary>
    ///     Transpose by a chunked counting sort: each of p contiguous element chunks counts its
    ///     incidences per node and gives them chunk-local ranks; an exclusive scan over the chunks
    ///     turns the ranks into positions. Chunk c precedes chunk c + 1 in every row, so the rows
    ///     are sorted by element, as in the sequential sort.
    /// </summary>
    private static (int[] Ptr, int[] Pairs) TransposeChunked(int[][] A, int nE, int nNodes)
    {
        var chunks = Math.Max(1, Math.Min(Environment.ProcessorCount, nE / 4096));
        chunks = (int)Math.Max(1, Math.Min(chunks, (512L << 20) / (4L * Math.Max(1, nNodes))));
        var size = (nE + chunks - 1) / chunks;
        var counts = new int[chunks][];
        Parallel.For(0, chunks, c =>
        {
            var local = new int[nNodes];
            for (int e = c * size, hi = Math.Min(nE, e + size); e < hi; e++)
                foreach (var v in A[e])
                    local[v]++;
            counts[c] = local;
        });

        var tptr = new int[nNodes + 1];
        Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(0, nNodes, 1 << 14), r =>
        {
            for (var v = r.Item1; v < r.Item2; v++)
            {
                var run = 0;
                for (var c = 0; c < chunks; c++)
                {
                    var here = counts[c][v];
                    counts[c][v] = run;
                    run += here;
                }

                tptr[v + 1] = run;
            }
        });
        for (var v = 0; v < nNodes; v++) tptr[v + 1] += tptr[v];

        var tpairs = new int[2 * tptr[nNodes]];
        Parallel.For(0, chunks, c =>
        {
            var next = counts[c];
            for (int e = c * size, hi = Math.Min(nE, e + size); e < hi; e++)
            {
                var row = A[e];
                for (var p = 0; p < row.Length; p++)
                {
                    var v = row[p];
                    var q = tptr[v] + next[v]++;
                    tpairs[2 * q] = e;
                    tpairs[2 * q + 1] = p;
                }
            }
        });
        return (tptr, tpairs);
    }

    /// <summary>
    ///     Deterministic parallel greedy colouring. The elements are split into contiguous blocks
    ///     (their number depends on the element count only); each block is coloured greedily in
    ///     element order against its own elements. An element that then has the colour of a smaller
    ///     neighbour in another block is uncoloured, and the uncoloured elements are recoloured in
    ///     rounds: an element is coloured once its uncoloured smaller neighbours are, with the
    ///     smallest colour unused by its coloured neighbours. Every colour depends on the colours
    ///     of earlier elements only, so the result does not depend on the thread schedule.
    /// </summary>
    private static (int[] Color, int Count) ColourParallel(int[][] A, int nE, int[] tptr, int[] tpairs)
    {
        var color = new int[nE];
        Array.Fill(color, -1);
        var nb = Math.Max(1, Math.Min(16, nE / 32768));
        int Lo(int b) => (int)((long)nE * b / nb);
        var block = new int[nE];
        Parallel.For(0, nb, b =>
        {
            for (int e = Lo(b), hi = Lo(b + 1); e < hi; e++) block[e] = b;
        });

        // ---- 1. greedy per block; cross[e]: e has a neighbour in another block
        var cross = new bool[nE];
        Parallel.For(0, nb, b =>
        {
            var used = new int[64];
            var stamp = 0;
            for (int e = Lo(b), hi = Lo(b + 1); e < hi; e++)
            {
                var row = A[e];
                if (row.Length == 0) continue;
                stamp++;
                var x = false;
                foreach (var v in row)
                    for (var q = tptr[v]; q < tptr[v + 1]; q++)
                    {
                        var g = tpairs[2 * q];
                        if (g == e) continue;
                        if (block[g] != b)
                        {
                            x = true;
                            continue;
                        }

                        var c = color[g];
                        if (c < 0) continue;
                        if (c >= used.Length) Array.Resize(ref used, Math.Max(c + 1, 2 * used.Length));
                        used[c] = stamp;
                    }

                var k = 0;
                while (k < used.Length && used[k] == stamp) k++;
                if (k == used.Length) Array.Resize(ref used, 2 * used.Length);
                color[e] = k;
                cross[e] = x;
            }
        });

        // ---- 2. conflicts across blocks: the larger element of a pair with equal colours
        var redo = new bool[nE];
        Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(0, nE, 1 << 14), r =>
        {
            for (var e = r.Item1; e < r.Item2; e++)
            {
                if (!cross[e]) continue;
                var ce = color[e];
                var conflict = false;
                foreach (var v in A[e])
                {
                    for (var q = tptr[v]; q < tptr[v + 1]; q++)
                    {
                        var g = tpairs[2 * q];
                        if (g < e && block[g] != block[e] && color[g] == ce)
                        {
                            conflict = true;
                            break;
                        }
                    }

                    if (conflict) break;
                }

                redo[e] = conflict;
            }
        });
        var pending = new List<int>();
        for (var e = 0; e < nE; e++)
            if (redo[e])
            {
                color[e] = -1;
                pending.Add(e);
            }

        // ---- 3. recolouring in rounds (the smallest pending element is always ready)
        var done = new bool[nE];
        while (pending.Count > 0)
        {
            var list = pending.ToArray();
            Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(0, list.Length, 256),
                () => (Used: new int[64], Stamp: 0), (r, _, st) =>
                {
                    for (var i = r.Item1; i < r.Item2; i++)
                    {
                        var e = list[i];
                        st.Stamp++;
                        var ready = true;
                        foreach (var v in A[e])
                        {
                            for (var q = tptr[v]; q < tptr[v + 1]; q++)
                            {
                                var g = tpairs[2 * q];
                                if (g == e) continue;
                                var c = Volatile.Read(ref color[g]);
                                if (c < 0)
                                {
                                    if (redo[g] && g < e)
                                    {
                                        ready = false;
                                        break;
                                    }

                                    continue;
                                }

                                if (c >= st.Used.Length) Array.Resize(ref st.Used, Math.Max(c + 1, 2 * st.Used.Length));
                                st.Used[c] = st.Stamp;
                            }

                            if (!ready) break;
                        }

                        if (!ready) continue;
                        var k = 0;
                        while (k < st.Used.Length && st.Used[k] == st.Stamp) k++;
                        if (k == st.Used.Length) Array.Resize(ref st.Used, 2 * st.Used.Length);
                        Volatile.Write(ref color[e], k);
                        done[e] = true;
                    }

                    return st;
                }, _ => { });
            pending.Clear();
            foreach (var e in list)
                if (!done[e])
                    pending.Add(e);
        }

        var count = 0;
        foreach (var c in color)
            if (c + 1 > count)
                count = c + 1;
        return (color, count);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ReadOnlySpan<int> Row(int[] tptr, int[] tpairs, int v) =>
        new(tpairs, 2 * tptr[v], 2 * (tptr[v + 1] - tptr[v]));

    private static int CountRow(int v, int[] tptr, int[] tpairs, int[][] A, int[] stamp)
    {
        var mark = v + 1;
        var len = 0;
        for (var q = tptr[v]; q < tptr[v + 1]; q++)
            foreach (var m in A[tpairs[2 * q]])
                if (stamp[m] != mark)
                {
                    stamp[m] = mark;
                    len++;
                }

        return len;
    }

    private static void FillRow(int v, int K, int KK, int[] tptr, int[] tpairs, int[][] A,
        int[] pptr, int[] pcol, int[] pcnt, int[] dest, int[] stamp, int[] slot)
    {
        var mark = v + 1;
        var b = pptr[v];
        var len = 0;
        for (var q = tptr[v]; q < tptr[v + 1]; q++)
        {
            var e = tpairs[2 * q];
            var i = tpairs[2 * q + 1];
            var row = A[e];
            for (var j = 0; j < row.Length; j++)
            {
                var m = row[j];
                if (stamp[m] != mark)
                {
                    stamp[m] = mark;
                    slot[m] = len;
                    pcol[b + len] = m;
                    len++;
                }

                var sl = slot[m];
                pcnt[b + sl]++;
                dest[e * KK + i * K + j] = sl;
            }
        }
    }

    internal static int SmallestFree(int e, int[] row, int[] tptr, int[] tpairs, int[] color, ref int[] used)
    {
        var mark = e + 1;
        foreach (var v in row)
            for (var q = tptr[v]; q < tptr[v + 1]; q++)
            {
                var g = tpairs[2 * q];
                if (g == e) continue;
                var c = color[g];
                if (c < 0) continue;
                if (c >= used.Length) Array.Resize(ref used, Math.Max(c + 1, used.Length * 2));
                used[c] = mark;
            }

        var k = 0;
        while (k < used.Length && used[k] == mark) k++;
        if (k == used.Length) Array.Resize(ref used, used.Length * 2);
        return k;
    }
}

// SubdivisionTemplates.cs
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

// The Steiner-point-free subdivision templates of edge refinement: the single C#
// implementation of the case logic of remeshsimplex.f90 (splittetwork1, splittria31,
// splitbar21), converted from 1-based to 0-based indexing. Used by MeshRefinement.Refine
// (which appends the children to a SimplexMesh) and by the incremental edge cutter of
// Step2.Core (which writes them into its own arrays); children are emitted as rows.
//
// Every decision is either the pattern of marked edges (a midpoint label or the
// sentinel) or an order comparison of global labels, Max(a, b) > Max(c, d) and
// mid1 > mid2: the subdivision of an element is a function of its row, its midpoint
// labels and the relative order of these labels, so two elements sharing a face choose
// the same diagonals and the refinement is conforming.

namespace Numerical;

/// <summary>
///     Subdivision templates of bars, triangles and tetrahedra with marked edges
///     (no Steiner points). Children are appended to an output list as node rows.
/// </summary>
public static class SubdivisionTemplates
{
    // local edges of a tetrahedron (0-based), node -> incident edges, opposite-edge nodes
    public static readonly int[][] EdgeNodes =
        [[0, 1], [1, 2], [0, 2], [2, 3], [3, 0], [1, 3]];

    private static readonly int[][] NodeEdges =
        [[0, 2, 4], [0, 1, 5], [1, 2, 3], [3, 4, 5]];

    private static readonly int[] Prev = [2, 0, 3, 0, 2, 2];
    private static readonly int[] Post = [3, 3, 1, 1, 1, 0];

    /// <summary>Number of tetrahedra produced for a given marked-edge pattern (for sizing).</summary>
    public static int TetChildren(int nmarked, ReadOnlySpan<int> midnodes)
    {
        switch (nmarked)
        {
            case 0: return 1;
            case 1: return 2;
            case 2: return MaxMarksAtNode(midnodes) == 2 ? 3 : 4;
            case 3:
            {
                var (n0, n3) = (0, 0);
                for (var i = 0; i < 4; i++)
                {
                    var c = MarksAt(i, midnodes);
                    if (c == 0) n0++;
                    if (c == 3) n3++;
                }

                return n3 > 0 ? 4 : n0 > 0 ? 4 : 5;
            }
            case 4: return 6;
            case 5: return 7;
            default: return 8;
        }
    }

    private static int MarksAt(int i, ReadOnlySpan<int> midnodes)
    {
        var c = 0;
        foreach (var ie in NodeEdges[i])
            if (midnodes[ie] > 0)
                c++;
        return c;
    }

    private static int MaxMarksAtNode(ReadOnlySpan<int> midnodes)
    {
        var m = 0;
        for (var i = 0; i < 4; i++) m = Math.Max(m, MarksAt(i, midnodes));
        return m;
    }

    // ------------------------------------------------------------------ bars and triangles

    public static void SplitBar(int v0, int v1, int mid, List<int[]> outRows)
    {
        if (mid < 0)
        {
            outRows.Add([v0, v1]);
            return;
        }

        outRows.Add([v0, mid]);
        outRows.Add([mid, v1]);
    }

    /// <summary>Triangle cases; a midpoint is -1 when its edge is not marked.</summary>
    public static void SplitTriangle(int v0, int v1, int v2, int m01, int m12, int m20, List<int[]> o)
    {
        if (m01 < 0 && m12 < 0 && m20 < 0)
            o.Add([v0, v1, v2]);
        else if (m01 >= 0 && m12 < 0 && m20 < 0)
        {
            o.Add([v0, m01, v2]);
            o.Add([m01, v1, v2]);
        }
        else if (m01 < 0 && m12 >= 0 && m20 < 0)
        {
            o.Add([v1, m12, v0]);
            o.Add([m12, v2, v0]);
        }
        else if (m01 < 0 && m12 < 0 && m20 >= 0)
        {
            o.Add([v2, m20, v1]);
            o.Add([m20, v0, v1]);
        }
        else if (m01 >= 0 && m12 >= 0 && m20 < 0)
            SplitTri2Edges(v0, m01, v1, m12, v2, o);
        else if (m01 >= 0 && m12 < 0 && m20 >= 0)
            SplitTri2Edges(v2, m20, v0, m01, v1, o);
        else if (m01 < 0 && m12 >= 0 && m20 >= 0)
            SplitTri2Edges(v1, m12, v2, m20, v0, o);
        else
        {
            o.Add([v0, m01, m20]);
            o.Add([m01, v1, m12]);
            o.Add([m12, v2, m20]);
            o.Add([m01, m12, m20]);
        }
    }

    private static void SplitTri2Edges(int nn1, int mid1, int nnMeet, int mid2, int nn5, List<int[]> o)
    {
        o.Add([mid1, nnMeet, mid2]);
        if (mid1 > mid2)
        {
            o.Add([nn1, mid1, nn5]);
            o.Add([mid1, mid2, nn5]);
        }
        else
        {
            o.Add([nn1, mid1, mid2]);
            o.Add([nn1, mid2, nn5]);
        }
    }

    // ------------------------------------------------------------------ tetrahedra

    /// <summary>
    ///     Subdivision of a tetrahedron; <paramref name="midnodes"/>[ie] is the midpoint
    ///     label of local edge ie (<see cref="EdgeNodes"/>) or 0 when the edge is not marked.
    /// </summary>
    public static void SplitTet(ReadOnlySpan<int> nodes, ReadOnlySpan<int> midnodes, List<int[]> o)
    {
        Span<int> marked = stackalloc int[6];
        var nmarked = 0;
        for (var i = 0; i < 6; i++)
            if (midnodes[i] > 0)
                marked[nmarked++] = i;

        Span<int> imarks = stackalloc int[5];
        Span<int> jmarks = stackalloc int[12];
        imarks[0] = 0;
        for (var i = 0; i < 4; i++)
        {
            imarks[i + 1] = imarks[i];
            for (var j = 0; j < 3; j++)
            {
                var ie = NodeEdges[i][j];
                if (midnodes[ie] > 0)
                {
                    jmarks[imarks[i + 1]] = ie;
                    imarks[i + 1]++;
                }
            }
        }

        switch (nmarked)
        {
            case 0:
                o.Add([nodes[0], nodes[1], nodes[2], nodes[3]]);
                break;
            case 1:
                Case1(nodes, midnodes, marked, o);
                break;
            case 2:
                Case2(nodes, midnodes, marked, imarks, o);
                break;
            case 3:
                Case3(nodes, midnodes, imarks, jmarks, o);
                break;
            case 4:
                Case4(nodes, midnodes, imarks, o);
                break;
            case 5:
                Case5(nodes, midnodes, o);
                break;
            default:
                Case6(nodes, midnodes, o);
                break;
        }
    }

    private static void Case1(ReadOnlySpan<int> nodes, ReadOnlySpan<int> mid, ReadOnlySpan<int> marked, List<int[]> o)
    {
        var n1 = EdgeNodes[marked[0]][0];
        var n2 = EdgeNodes[marked[0]][1];
        var (n3, n4) = From2Nodes(n1, n2);
        var n5 = mid[marked[0]];
        n1 = nodes[n1];
        n2 = nodes[n2];
        n3 = nodes[n3];
        n4 = nodes[n4];
        o.Add([n3, n4, n5, n2]);
        o.Add([n4, n3, n5, n1]);
    }

    private static void Case2(ReadOnlySpan<int> nodes, ReadOnlySpan<int> mid, ReadOnlySpan<int> marked,
        ReadOnlySpan<int> imarks, List<int[]> o)
    {
        var mm = 0;
        for (var i = 0; i < 4; i++) mm = Math.Max(mm, imarks[i + 1] - imarks[i]);
        if (mm == 2)
        {
            int n3 = -1, n4 = -1;
            for (var i = 0; i < 4; i++)
            {
                if (imarks[i + 1] - imarks[i] == 0) n3 = i;
                if (imarks[i + 1] - imarks[i] == 2) n4 = i;
            }

            var (n1, n2) = From2Nodes(n3, n4);
            var n5 = Split(n2, n4, mid);
            var n6 = Split(n1, n4, mid);
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            Pyramid(n1, n2, n5, n6, n3, o);
            o.Add([n6, n5, n3, n4]);
        }
        else
        {
            var ie1 = marked[0];
            var n1 = EdgeNodes[ie1][0];
            var n2 = EdgeNodes[ie1][1];
            var (n3, n4) = From2Nodes(n1, n2);
            var n5 = Split(n1, n2, mid);
            var n6 = Split(n3, n4, mid);
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            o.Add([n2, n3, n5, n6]);
            o.Add([n2, n6, n5, n4]);
            o.Add([n1, n5, n6, n4]);
            o.Add([n1, n5, n3, n6]);
        }
    }

    private static void Case3(ReadOnlySpan<int> nodes, ReadOnlySpan<int> mid, ReadOnlySpan<int> imarks,
        ReadOnlySpan<int> jmarks, List<int[]> o)
    {
        int n3 = -1, n4 = -1;
        for (var i = 0; i < 4; i++)
        {
            if (imarks[i + 1] - imarks[i] == 3) n4 = i;
            if (imarks[i + 1] - imarks[i] == 0) n3 = i;
        }

        if (n4 >= 0)
        {
            var n1 = (n4 + 1) % 4;
            int n2;
            (n2, n3) = From2Nodes(n1, n4);
            var n5 = Split(n3, n4, mid);
            var n6 = Split(n1, n4, mid);
            var n7 = Split(n2, n4, mid);
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            o.Add([n4, n7, n6, n5]);
            Prism(n2, n7, n6, n1, n3, n5, o);
        }
        else if (n3 >= 0)
        {
            var n4l = (n3 + 1) % 4;
            var (n1, n2) = From2Nodes(n3, n4l);
            var n5 = Split(n1, n4l, mid);
            var n6 = Split(n1, n2, mid);
            var n7 = Split(n2, n4l, mid);
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4l = nodes[n4l];
            o.Add([n1, n6, n3, n5]);
            o.Add([n2, n3, n6, n7]);
            o.Add([n6, n3, n5, n7]);
            o.Add([n5, n7, n3, n4l]);
        }
        else
        {
            int n1 = -1, n2 = -1;
            n3 = -1;
            n4 = -1;
            for (var i = 0; i < 4; i++)
            {
                if (imarks[i + 1] - imarks[i] != 1) continue;
                var ie = jmarks[imarks[i]];
                if (n1 < 0)
                {
                    n1 = i;
                    n4 = EdgeNodes[ie][0] == n1 ? EdgeNodes[ie][1] : EdgeNodes[ie][0];
                }
                else
                {
                    n2 = i;
                    n3 = EdgeNodes[ie][0] == n2 ? EdgeNodes[ie][1] : EdgeNodes[ie][0];
                }
            }

            if (n1 < 0 || n2 < 0) throw new InvalidOperationException("three marked edges: path case not found");
            var n5 = Split(n1, n4, mid);
            var n6 = Split(n2, n3, mid);
            var n7 = Split(n3, n4, mid);
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            Pyramid(n6, n7, n4, n2, n5, o);
            Pyramid(n3, n1, n5, n7, n6, o);
            o.Add([n1, n2, n6, n5]);
        }
    }

    private static void Case4(ReadOnlySpan<int> nodes, ReadOnlySpan<int> mid, ReadOnlySpan<int> imarks, List<int[]> o)
    {
        var mm = 0;
        for (var i = 0; i < 4; i++) mm = Math.Max(mm, imarks[i + 1] - imarks[i]);
        if (mm == 3)
        {
            int n3 = -1, n4 = -1;
            for (var i = 0; i < 4; i++)
            {
                if (imarks[i + 1] - imarks[i] == 1) n3 = i;
                if (imarks[i + 1] - imarks[i] == 3) n4 = i;
            }

            var (n1, n2) = From2Nodes(n3, n4);
            var n5 = Split(n1, n2, mid);
            var n6 = Split(n2, n4, mid);
            var n7 = Split(n3, n4, mid);
            var n8 = Split(n1, n4, mid);
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            o.Add([n8, n6, n7, n4]);
            o.Add([n7, n8, n5, n6]);
            Pyramid(n2, n3, n7, n6, n5, o);
            Pyramid(n1, n8, n7, n3, n5, o);
        }
        else
        {
            var ie = UnmarkedEdge(mid);
            var n1 = EdgeNodes[ie][0];
            var n2 = EdgeNodes[ie][1];
            var (n3, n4) = From2Nodes(n1, n2);
            var n5 = Split(n1, n4, mid);
            var n6 = Split(n1, n3, mid);
            var n7 = Split(n2, n3, mid);
            var n8 = Split(n2, n4, mid);
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            Prism(n8, n7, n3, n4, n5, n6, o);
            Prism(n1, n2, n8, n5, n6, n7, o);
        }
    }

    private static void Case5(ReadOnlySpan<int> nodes, ReadOnlySpan<int> mid, List<int[]> o)
    {
        var ie = UnmarkedEdge(mid);
        var n1 = EdgeNodes[ie][0];
        var n2 = EdgeNodes[ie][1];
        var (n3, n4) = From2Nodes(n1, n2);
        var n5 = Split(n1, n4, mid);
        var n6 = Split(n1, n3, mid);
        var n7 = Split(n2, n3, mid);
        var n8 = Split(n3, n4, mid);
        var n9 = Split(n2, n4, mid);
        n1 = nodes[n1];
        n2 = nodes[n2];
        n3 = nodes[n3];
        n4 = nodes[n4];
        o.Add([n9, n8, n5, n4]);
        o.Add([n3, n6, n7, n8]);
        Pyramid(n6, n7, n9, n5, n8, o);
        Prism(n1, n2, n9, n5, n6, n7, o);
    }

    private static void Case6(ReadOnlySpan<int> nodes, ReadOnlySpan<int> mid, List<int[]> o)
    {
        var n5 = Split(0, 3, mid);
        var n6 = Split(1, 3, mid);
        var n7 = Split(2, 3, mid);
        var n8 = Split(0, 2, mid);
        var n9 = Split(0, 1, mid);
        var n10 = Split(1, 2, mid);
        int n1 = nodes[0], n2 = nodes[1], n3 = nodes[2], n4 = nodes[3];
        o.Add([n1, n9, n8, n5]);
        o.Add([n2, n10, n9, n6]);
        o.Add([n3, n8, n10, n7]);
        o.Add([n5, n6, n7, n4]);
        Octahedron(n8, n6, n7, n5, n9, n10, o);
    }

    private static int UnmarkedEdge(ReadOnlySpan<int> mid)
    {
        for (var i = 0; i < 6; i++)
            if (mid[i] == 0)
                return i;
        throw new InvalidOperationException("no unmarked edge");
    }

    private static (int, int) From2Nodes(int i2, int i3)
    {
        for (var ie = 0; ie < 6; ie++)
        {
            if (EdgeNodes[ie][0] == i2 && EdgeNodes[ie][1] == i3) return (Prev[ie], Post[ie]);
            if (EdgeNodes[ie][0] == i3 && EdgeNodes[ie][1] == i2) return (Post[ie], Prev[ie]);
        }

        throw new InvalidOperationException($"edge ({i2},{i3}) not found");
    }

    private static int Split(int n1, int n2, ReadOnlySpan<int> mid)
    {
        for (var ie = 0; ie < 6; ie++)
            if ((EdgeNodes[ie][0] == n1 && EdgeNodes[ie][1] == n2) ||
                (EdgeNodes[ie][0] == n2 && EdgeNodes[ie][1] == n1))
            {
                if (mid[ie] == 0) throw new InvalidOperationException($"edge {ie} has no midpoint");
                return mid[ie];
            }

        throw new InvalidOperationException($"edge ({n1},{n2}) not found");
    }

    // ------------------------------------------------------------------ polyhedra (global-label rule)

    private static void Pyramid(int i1, int i2, int i3, int i4, int i5, List<int[]> o)
    {
        if (Math.Max(i1, i3) > Math.Max(i2, i4))
        {
            o.Add([i1, i2, i5, i3]);
            o.Add([i1, i5, i4, i3]);
        }
        else
        {
            o.Add([i1, i2, i5, i4]);
            o.Add([i2, i3, i5, i4]);
        }
    }

    private static void Prism(int i1, int i2, int i3, int i4, int i5, int i6, List<int[]> o)
    {
        var swap46 = Math.Max(i4, i6) > Math.Max(i3, i5);
        var swap13 = Math.Max(i1, i3) > Math.Max(i2, i4);
        var swap16 = Math.Max(i1, i6) > Math.Max(i2, i5);
        if (swap46)
        {
            if (swap13)
            {
                o.Add([i1, i2, i6, i3]);
                o.Add([i1, i6, i4, i3]);
                o.Add([i1, i6, i5, i4]);
            }
            else if (swap16)
            {
                o.Add([i1, i2, i6, i4]);
                o.Add([i1, i6, i5, i4]);
                o.Add([i2, i6, i4, i3]);
            }
            else
            {
                o.Add([i1, i2, i5, i4]);
                o.Add([i2, i3, i6, i4]);
                o.Add([i2, i6, i5, i4]);
            }
        }
        else if (swap13)
        {
            if (swap16)
            {
                o.Add([i1, i2, i6, i3]);
                o.Add([i1, i3, i5, i4]);
                o.Add([i1, i6, i5, i3]);
            }
            else
            {
                o.Add([i1, i3, i5, i4]);
                o.Add([i2, i5, i1, i3]);
                o.Add([i2, i6, i5, i3]);
            }
        }
        else
        {
            o.Add([i1, i2, i5, i4]);
            o.Add([i2, i3, i5, i4]);
            o.Add([i2, i6, i5, i3]);
        }
    }

    private static void Octahedron(int i1, int i2, int i3, int i4, int i5, int i6, List<int[]> o)
    {
        if (Math.Max(i3, i5) > Math.Max(i4, i6))
        {
            o.Add([i1, i3, i4, i5]);
            o.Add([i1, i3, i5, i6]);
            o.Add([i3, i4, i5, i2]);
            o.Add([i3, i5, i6, i2]);
        }
        else
        {
            o.Add([i1, i4, i6, i3]);
            o.Add([i1, i6, i4, i5]);
            o.Add([i2, i4, i6, i5]);
            o.Add([i2, i6, i4, i3]);
        }
    }
}

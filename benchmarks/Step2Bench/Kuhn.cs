// Kuhn.cs — structured tetrahedral box meshes (Kuhn / Freudenthal triangulation).
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
// Each cube of an n x n x n grid is split into six tetrahedra sharing the main
// diagonal; the vertices of every tetrahedron are listed along the monotone lattice
// path from the cube's low corner to its high corner. With this ordering and tag 3,
// the tetrahedra are the tagged initial simplices of Maubach's bisection, and the
// triangulation is conforming.

namespace Step2.Bench;

public static class Kuhn
{
    private static readonly int[][] Perms =
    [
        [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]
    ];

    public static int Node(int n, int i, int j, int k) => i + (n + 1) * (j + (n + 1) * k);

    /// <summary>Rows in Maubach order and node coordinates (unit cube, spacing 1/n).</summary>
    public static (int[][] Rows, List<double> Coords, int NodeCount) Box(int n)
    {
        var nN = (n + 1) * (n + 1) * (n + 1);
        var coords = new List<double>(3 * nN);
        for (var k = 0; k <= n; k++)
        for (var j = 0; j <= n; j++)
        for (var i = 0; i <= n; i++)
        {
            coords.Add((double)i / n);
            coords.Add((double)j / n);
            coords.Add((double)k / n);
        }

        var rows = new int[6 * n * n * n][];
        var t = 0;
        Span<int> c = stackalloc int[3];
        for (var k = 0; k < n; k++)
        for (var j = 0; j < n; j++)
        for (var i = 0; i < n; i++)
            foreach (var perm in Perms)
            {
                c[0] = i;
                c[1] = j;
                c[2] = k;
                var row = new int[4];
                row[0] = Node(n, c[0], c[1], c[2]);
                for (var s = 0; s < 3; s++)
                {
                    c[perm[s]]++;
                    row[s + 1] = Node(n, c[0], c[1], c[2]);
                }

                rows[t++] = row;
            }

        return (rows, coords, nN);
    }

    /// <summary>Cube count per side for a target element count (6 n^3 elements).</summary>
    public static int SideFor(long targetElements) => (int)Math.Max(1, Math.Round(Math.Cbrt(targetElements / 6.0)));
}

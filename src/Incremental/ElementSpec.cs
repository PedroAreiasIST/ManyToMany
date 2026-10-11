// ElementSpec.cs — local topology of a simplicial element type.
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

namespace Step2;

/// <summary>
///     Local topology of a simplex: K nodes, K facets, facet f opposite local node f.
///     The facet rule (facet opposite a node) is what makes the facet index of a
///     neighbour recoverable from the neighbour's row alone.
///     A mixed specification admits rows of 1 to K nodes in one mesh (point, bar,
///     triangle, tetrahedron); K is then the stride of the per-element tables, the
///     row length is the element type, and facets are matched only between elements
///     of the same length.
/// </summary>
public sealed class ElementSpec
{
    public readonly int K;            // nodes per element (maximum, and stride, when mixed)
    public readonly int Nf;           // facets per element (stride)
    public readonly int Fs;           // nodes per facet
    public readonly int[][] FacetLocal; // FacetLocal[f] = local nodes of facet f (all but f)
    public readonly bool Mixed;

    /// <summary>FacetLocalByLength[k][f]: local nodes of facet f of a k-node row (no facets for k &lt; 2).</summary>
    public readonly int[][][] FacetLocalByLength;

    private ElementSpec(int k, bool mixed)
    {
        K = k;
        Nf = k;
        Fs = k - 1;
        Mixed = mixed;
        FacetLocalByLength = new int[k + 1][][];
        for (var len = 0; len <= k; len++) FacetLocalByLength[len] = Facets(len);
        FacetLocal = FacetLocalByLength[k];
    }

    private static int[][] Facets(int k)
    {
        if (k < 2) return [];
        var t = new int[k][];
        for (var f = 0; f < k; f++)
        {
            var list = new int[k - 1];
            var n = 0;
            for (var p = 0; p < k; p++)
                if (p != f) list[n++] = p;
            t[f] = list;
        }

        return t;
    }

    /// <summary>Number of facets of a row of <paramref name="length"/> nodes.</summary>
    public int FacetCount(int length) => length < 2 ? 0 : length;

    public static readonly ElementSpec Triangle = new(3, false);
    public static readonly ElementSpec Tetrahedron = new(4, false);

    /// <summary>Points, bars, triangles and tetrahedra in one mesh (generalized elements).</summary>
    public static readonly ElementSpec MixedSimplex = new(4, true);
}

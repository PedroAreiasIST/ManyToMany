namespace Numerical;

/// <summary>
///     Predefined sub-entity definitions for standard finite element types.
/// </summary>
/// <remarks>
///     <para>
///         Provides ready-made edge and face definitions for common element topologies:
///         Bar2, Tri3, Quad4, Tet4, Hex8, Wedge6, and Pyramid5.
///     </para>
///     <para>
///         These are convenience definitions built on <see cref="SubEntityDefinition" />.
///         For custom or higher-order elements, construct a <see cref="SubEntityDefinition" /> directly.
///     </para>
/// </remarks>
public static class FiniteElementTopologies
{
    /// <summary>Bar2 (2-node line) edge: the element itself.</summary>
    public static SubEntityDefinition Bar2Edge => SubEntityDefinition.FromEdges((0, 1));

    /// <summary>Tri3 (3-node triangle) edges.</summary>
    public static SubEntityDefinition Tri3Edges => SubEntityDefinition.FromEdges((0, 1), (1, 2), (2, 0));

    /// <summary>Quad4 (4-node quadrilateral) edges.</summary>
    public static SubEntityDefinition Quad4Edges => SubEntityDefinition.FromEdges((0, 1), (1, 2), (2, 3), (3, 0));

    /// <summary>Tet4 (4-node tetrahedron) edges.</summary>
    public static SubEntityDefinition Tet4Edges => SubEntityDefinition.FromEdges(
        (0, 1), (1, 2), (0, 2), (2, 3), (0, 3), (1, 3));

    /// <summary>Tet4 (4-node tetrahedron) faces (outward normals with right-hand rule).</summary>
    public static SubEntityDefinition Tet4Faces => SubEntityDefinition.FromFaces(
        (0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3));

    /// <summary>Hex8 (8-node hexahedron) edges.</summary>
    public static SubEntityDefinition Hex8Edges => SubEntityDefinition.FromEdges(
        (0, 1), (1, 2), (2, 3), (3, 0), // Bottom face
        (4, 5), (5, 6), (6, 7), (7, 4), // Top face
        (0, 4), (1, 5), (2, 6), (3, 7)); // Vertical edges

    /// <summary>Hex8 (8-node hexahedron) faces.</summary>
    public static SubEntityDefinition Hex8Faces => SubEntityDefinition.FromQuadFaces(
        (0, 3, 2, 1), (4, 5, 6, 7), // Bottom, Top
        (0, 1, 5, 4), (2, 3, 7, 6), // Front, Back
        (0, 4, 7, 3), (1, 2, 6, 5)); // Left, Right

    /// <summary>Wedge6 (6-node prism) edges.</summary>
    public static SubEntityDefinition Wedge6Edges => SubEntityDefinition.FromEdges(
        (0, 1), (1, 2), (2, 0), // Bottom triangle
        (3, 4), (4, 5), (5, 3), // Top triangle
        (0, 3), (1, 4), (2, 5)); // Vertical edges

    /// <summary>Wedge6 (6-node prism) triangular faces.</summary>
    public static SubEntityDefinition Wedge6TriFaces => SubEntityDefinition.FromFaces(
        (0, 2, 1), (3, 4, 5)); // Bottom, Top triangles

    /// <summary>Pyramid5 (5-node pyramid) edges.</summary>
    public static SubEntityDefinition Pyramid5Edges => SubEntityDefinition.FromEdges(
        (0, 1), (1, 2), (2, 3), (3, 0), // Base quad
        (0, 4), (1, 4), (2, 4), (3, 4)); // Edges to apex

    /// <summary>Pyramid5 (5-node pyramid) triangular faces.</summary>
    public static SubEntityDefinition Pyramid5TriFaces => SubEntityDefinition.FromFaces(
        (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)); // Four triangular side faces

    /// <summary>
    ///     Selects the predefined edge definition for a given node count.
    /// </summary>
    /// <param name="nodesPerElement">Number of nodes per element (2=Bar2, 3=Tri3, 4=Tet4/Quad4, etc.)</param>
    /// <returns>The matching <see cref="SubEntityDefinition" />.</returns>
    /// <exception cref="ArgumentException">Thrown for unsupported element node counts.</exception>
    /// <remarks>
    ///     For Quad4 vs Tet4 disambiguation (both 4 nodes), use the specific property directly.
    ///     This method assumes Tet4 for 4-node elements.
    /// </remarks>
    public static SubEntityDefinition EdgeDefinition(int nodesPerElement)
    {
        return nodesPerElement switch
        {
            2 => Bar2Edge,
            3 => Tri3Edges,
            4 => Tet4Edges, // Assumes Tet4, not Quad4
            6 => Wedge6Edges,
            8 => Hex8Edges,
            _ => throw new ArgumentException(
                $"No predefined edge topology for {nodesPerElement}-node elements. " +
                $"Use DiscoverSubEntities with a custom SubEntityDefinition.",
                nameof(nodesPerElement))
        };
    }
}

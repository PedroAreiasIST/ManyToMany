// SimplexMesh.cs - Core mesh types, constants, and topology definitions
// Consolidates: SimplexMesh.cs + MeshConstants.cs + FiniteElementTopologies.cs
// License: GPLv3

namespace Numerical;

// ═══════════════════════════════════════════════════════════════════════════
// Shared Constants (from MeshConstants.cs)
// ═══════════════════════════════════════════════════════════════════════════

/// <summary>
///     Shared constants and tolerance values for mesh operations.
///     All mesh library classes should use these values for consistency.
/// </summary>
public static class MeshConstants
{
    /// <summary>
    ///     Default numerical tolerance for geometric comparisons.
    ///     Used for: degenerate element detection, point coincidence, etc.
    /// </summary>
    public const double Epsilon = 1e-10;

    /// <summary>
    ///     Tolerance for detecting degenerate triangle areas.
    /// </summary>
    public const double DegenerateAreaTolerance = 1e-14;

    /// <summary>
    ///     Tolerance for detecting degenerate tetrahedron volumes.
    /// </summary>
    public const double DegenerateVolumeTolerance = 1e-15;

    /// <summary>
    ///     Default tolerance for node merging operations.
    /// </summary>
    public const double NodeMergeTolerance = 1e-12;

    /// <summary>
    ///     Small perturbation factor for grid point generation.
    /// </summary>
    public const double GridPerturbationFactor = 0.15;

    /// <summary>
    ///     Hexagonal grid row spacing factor (sqrt(3)/2).
    /// </summary>
    public const double HexRowSpacing = 0.86602540378;
}

// ═══════════════════════════════════════════════════════════════════════════
// Element Type Markers and Data Structures
// ═══════════════════════════════════════════════════════════════════════════

public readonly struct Node;

public readonly struct Edge;

public readonly struct Point;

public readonly struct Bar2; // 2-node line element

public readonly struct Tri3; // 3-node triangle

public readonly struct Quad4; // 4-node quadrilateral

public readonly struct Tet4; // 4-node tetrahedron

// Data structures attached to elements
public readonly record struct ParentNodes(int Parent1, int Parent2);

public readonly record struct OriginalElement(int Index);

// ═══════════════════════════════════════════════════════════════════════════
// SimplexMesh Topology Container
// ═══════════════════════════════════════════════════════════════════════════

/// <summary>
///     Topology container for simplex meshes (nodes, edges, triangles, tetrahedra).
///     Provides mesh connectivity, element addition, and data attachment capabilities.
///     Uses the Topology library for efficient connectivity management.
/// </summary>
public sealed class SimplexMesh : Topology<TypeMap<Node, Edge, Point, Bar2, Tri3, Quad4, Tet4>>
{
    /// <summary>
    ///     Initialize mesh with symmetric edge connectivity.
    /// </summary>
    public SimplexMesh()
    {
        WithSymmetry<Edge>(Symmetry.Full(2));
    }

    #region Node Operations

    /// <summary>
    ///     Adds a node and sets its ParentNodes data.
    ///     For original nodes, parent index equals node index.
    /// </summary>
    public int AddNode(int parentIndex)
    {
        var idx = Add<Node>();
        Set<Node, ParentNodes>(idx, new ParentNodes(parentIndex, parentIndex));
        return idx;
    }

    /// <summary>
    ///     Adds a node as a midpoint between two parent nodes.
    ///     Used during mesh refinement to track node lineage.
    /// </summary>
    public int AddMidpointNode(int parent1, int parent2)
    {
        var idx = Add<Node>();
        Set<Node, ParentNodes>(idx, new ParentNodes(parent1, parent2));
        return idx;
    }

    #endregion

    #region 2D Element Operations

    /// <summary>
    ///     Adds a triangle element (3 nodes).
    /// </summary>
    public int AddTriangle(int n0, int n1, int n2)
    {
        var idx = Add<Tri3, Node>(n0, n1, n2);
        Set<Tri3, OriginalElement>(idx, new OriginalElement(idx));
        return idx;
    }

    /// <summary>
    ///     Adds a quadrilateral element (4 nodes).
    /// </summary>
    public int AddQuad(int n0, int n1, int n2, int n3)
    {
        var idx = Add<Quad4, Node>(n0, n1, n2, n3);
        Set<Quad4, OriginalElement>(idx, new OriginalElement(idx));
        return idx;
    }

    #endregion

    #region 3D Element Operations

    /// <summary>
    ///     Adds a tetrahedron element (4 nodes).
    /// </summary>
    public int AddTetrahedron(int n0, int n1, int n2, int n3)
    {
        var idx = Add<Tet4, Node>(n0, n1, n2, n3);
        Set<Tet4, OriginalElement>(idx, new OriginalElement(idx));
        return idx;
    }

    /// <summary>
    ///     Adds a tetrahedron element (alias for AddTetrahedron).
    /// </summary>
    public int AddTet(int n0, int n1, int n2, int n3)
    {
        return AddTetrahedron(n0, n1, n2, n3);
    }

    #endregion

    #region 1D and 0D Element Operations

    /// <summary>
    ///     Adds a bar (line) element (2 nodes).
    /// </summary>
    public int AddBar(int n0, int n1)
    {
        var idx = Add<Bar2, Node>(n0, n1);
        Set<Bar2, OriginalElement>(idx, new OriginalElement(idx));
        return idx;
    }

    /// <summary>
    ///     Adds a point element (1 node).
    /// </summary>
    public int AddPoint(int n0)
    {
        var idx = Add<Point, Node>(n0);
        Set<Point, OriginalElement>(idx, new OriginalElement(idx));
        return idx;
    }

    #endregion
}

// ═══════════════════════════════════════════════════════════════════════════
// Finite Element Topologies (from FiniteElementTopologies.cs)
// ═══════════════════════════════════════════════════════════════════════════

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

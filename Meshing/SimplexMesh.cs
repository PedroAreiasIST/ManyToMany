// SimplexMesh.cs - Core mesh topology for simplex elements
// Extracted from SimplexRemesher.cs as part of library refactoring
// License: GPLv3

namespace Numerical;

// Element type markers
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

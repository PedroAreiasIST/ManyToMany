// SimplexRemesher.cs - Mesh I/O, utilities, and backward-compatible wrappers
// Mesh file I/O (GiD, Gmsh, ASCII) and utility functions
// 
// This file provides:
//   ✓ Backward-compatible wrappers for refactored functions:
//     - CreateRectangularMesh() → calls MeshGeneration.CreateRectangularMesh()
//     - CreateBoxMesh()         → calls MeshGeneration.CreateBoxMesh()
//     - Refine()                → calls MeshRefinement.Refine()
//     - InterpolateCoordinates() → calls MeshRefinement.InterpolateCoordinates()
//   ✓ GiD I/O (SaveGiD, LoadGiD)
//   ✓ Gmsh I/O (SaveMSH, LoadMSH)
//   ✓ ASCII I/O (SaveASCII)
//   ✓ Utilities (PrintStats, DiscoverEdges)
//   ✓ Advanced (CreateCrackFromSignedField, CreateCrackFromSignedField3D)
//
// License: GPLv3

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using static Numerical.MeshGeometry;

namespace Numerical;


// NOTE: SimplexMesh class definition is in SimplexMesh.cs
// This file contains only the SimplexRemesher static class with I/O functions

/// <summary>
/// Simplex mesh refinement by edge bisection with conforming subdivision.
/// </summary>
public static class SimplexRemesher
{
    // ═══════════════════════════════════════════════════════════════════════════
    // Mesh Generation Functions (wrappers to MeshGeneration.cs)
    // ═══════════════════════════════════════════════════════════════════════════
    
    /// <summary>
    /// Creates a rectangular triangular mesh.
    /// </summary>
    /// <param name="nx">Number of divisions in x direction</param>
    /// <param name="ny">Number of divisions in y direction</param>
    /// <param name="xMin">Minimum x coordinate</param>
    /// <param name="xMax">Maximum x coordinate</param>
    /// <param name="yMin">Minimum y coordinate</param>
    /// <param name="yMax">Maximum y coordinate</param>
    /// <returns>Mesh and coordinates array</returns>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateRectangularMesh(
        int nx, int ny, double xMin, double xMax, double yMin, double yMax)
        => MeshGeneration.CreateRectangularMesh(nx, ny, xMin, xMax, yMin, yMax);

    /// <summary>
    /// Creates a rectangular triangular mesh with default bounds [0,1] x [0,1].
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateRectangularMesh(int nx, int ny)
        => MeshGeneration.CreateRectangularMesh(nx, ny, 0, 1, 0, 1);

    /// <summary>
    /// Creates a 3D box mesh of tetrahedra.
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateBoxMesh(
        int nx, int ny, int nz,
        double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
        => MeshGeneration.CreateBoxMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);

    /// <summary>
    /// Creates a 3D box mesh with default bounds [0,1] x [0,1] x [0,1].
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateBoxMesh(int nx, int ny, int nz)
        => MeshGeneration.CreateBoxMesh(nx, ny, nz, 0, 1, 0, 1, 0, 1);

    // ═══════════════════════════════════════════════════════════════════════════
    // Edge Refinement Functions (wrappers to MeshRefinement.cs)
    // ═══════════════════════════════════════════════════════════════════════════
    
    /// <summary>
    /// Refines a mesh by bisecting marked edges. Wrapper for MeshRefinement.Refine.
    /// </summary>
    /// <param name="mesh">Input mesh</param>
    /// <param name="markedEdges">Edges to refine (node pairs)</param>
    /// <returns>Refined mesh and node remapping array</returns>
    public static (SimplexMesh Mesh, int[] NodeRemap) Refine(
        SimplexMesh mesh, List<(int, int)> markedEdges)
        => MeshRefinement.Refine(mesh, markedEdges);

    /// <summary>
    /// Refines a mesh by bisecting marked edges. Wrapper for MeshRefinement.Refine.
    /// </summary>
    public static (SimplexMesh Mesh, int[] NodeRemap) Refine(
        SimplexMesh mesh, IEnumerable<(int, int)> markedEdges)
        => MeshRefinement.Refine(mesh, markedEdges.ToList());

    /// <summary>
    /// Interpolates coordinates for refined mesh nodes. Wrapper for MeshRefinement.InterpolateCoordinates.
    /// </summary>
    public static double[,] InterpolateCoordinates(SimplexMesh refinedMesh, double[,] originalCoords)
        => MeshRefinement.InterpolateCoordinates(refinedMesh, originalCoords);

    // ═══════════════════════════════════════════════════════════════════════════
    // Edge Refinement Internal Structures
    // ═══════════════════════════════════════════════════════════════════════════

    // Edge definitions: TetEdges[e] = (nodeA, nodeB)
    static readonly (int A, int B)[] TetEdges = new (int A, int B)[] { (0,1), (1,2), (0,2), (2,3), (0,3), (1,3) };
    
    // O(1) edge lookup: EdgeIdx[n1,n2] = edge index (-1 if invalid)
    static readonly int[,] EdgeIdx = BuildEdgeIndex();
    static int[,] BuildEdgeIndex()
    {
        var t = new int[4, 4];
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) t[i, j] = -1;
        for (int e = 0; e < 6; e++) { t[TetEdges[e].A, TetEdges[e].B] = e; t[TetEdges[e].B, TetEdges[e].A] = e; }
        return t;
    }

    public static void DiscoverEdges(SimplexMesh m)
    {
        if (m.Count<Edge>() > 0) return;
        if (m.Count<Bar2>() > 0) m.DiscoverSubEntities<Bar2, Edge, Node>(FiniteElementTopologies.Bar2Edge);
        if (m.Count<Tri3>() > 0) m.DiscoverSubEntities<Tri3, Edge, Node>(FiniteElementTopologies.Tri3Edges);
        if (m.Count<Tet4>() > 0) m.DiscoverSubEntities<Tet4, Edge, Node>(FiniteElementTopologies.Tet4Edges);
    }


    // ═══════════════════════════════════════════════════════════════════════════
    // Mesh I/O
    // ═══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Save mesh to Gmsh .msh (MSH 2.2 ASCII) format.
    /// </summary>
    public static void SaveMSH(SimplexMesh mesh, double[,] coordinates, string path)
    {
        ArgumentNullException.ThrowIfNull(mesh);
        ArgumentNullException.ThrowIfNull(coordinates);
        ArgumentNullException.ThrowIfNull(path);

        int numNodes = mesh.Count<Node>();
        int numPoints = mesh.Count<Point>();
        int numBars = mesh.Count<Bar2>();
        int numTris = mesh.Count<Tri3>();
        int numTets = mesh.Count<Tet4>();
        int numElements = numPoints + numBars + numTris + numTets;

        using var writer = new StreamWriter(path);
        
        // Header
        writer.WriteLine("$MeshFormat");
        writer.WriteLine("2.2 0 8");
        writer.WriteLine("$EndMeshFormat");

        // Nodes
        writer.WriteLine("$Nodes");
        writer.WriteLine(numNodes);
        for (int i = 0; i < numNodes; i++)
            writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                "{0} {1} {2} {3}", i + 1, coordinates[i, 0], coordinates[i, 1], coordinates[i, 2]));
        writer.WriteLine("$EndNodes");

        // Elements
        writer.WriteLine("$Elements");
        writer.WriteLine(numElements);
        
        int elemId = 1;
        // Gmsh element types: 15=point, 1=line, 2=triangle, 4=tetrahedron
        
        for (int i = 0; i < numPoints; i++)
        {
            var n = mesh.NodesOf<Point, Node>(i);
            writer.WriteLine($"{elemId++} 15 2 0 0 {n[0] + 1}");
        }
        for (int i = 0; i < numBars; i++)
        {
            var n = mesh.NodesOf<Bar2, Node>(i);
            writer.WriteLine($"{elemId++} 1 2 0 0 {n[0] + 1} {n[1] + 1}");
        }
        for (int i = 0; i < numTris; i++)
        {
            var n = mesh.NodesOf<Tri3, Node>(i);
            writer.WriteLine($"{elemId++} 2 2 0 0 {n[0] + 1} {n[1] + 1} {n[2] + 1}");
        }
        for (int i = 0; i < numTets; i++)
        {
            var n = mesh.NodesOf<Tet4, Node>(i);
            writer.WriteLine($"{elemId++} 4 2 0 0 {n[0] + 1} {n[1] + 1} {n[2] + 1} {n[3] + 1}");
        }
        
        writer.WriteLine("$EndElements");
    }

    /// <summary>
    /// Saves a SimplexMesh to Gmsh MSH format (v2.2 ASCII) with element tags.
    /// Each element has exactly 2 tags: physical tag and elementary tag.
    /// </summary>
    /// <param name="mesh">Mesh to save</param>
    /// <param name="coordinates">Node coordinates</param>
    /// <param name="path">Output file path</param>
    /// <param name="physicalNames">Optional dictionary mapping physical tag to (dimension, name)</param>
    /// <param name="pointTags">Optional physical tags for Point elements (length must match Point count)</param>
    /// <param name="barTags">Optional physical tags for Bar2 elements (length must match Bar2 count)</param>
    /// <param name="triTags">Optional physical tags for Tri3 elements (length must match Tri3 count)</param>
    /// <param name="tetTags">Optional physical tags for Tet4 elements (length must match Tet4 count)</param>
    public static void SaveMSH(SimplexMesh mesh, double[,] coordinates, string path,
        Dictionary<int, (int Dimension, string Name)>? physicalNames = null,
        int[]? pointTags = null,
        int[]? barTags = null,
        int[]? triTags = null,
        int[]? tetTags = null)
    {
        ArgumentNullException.ThrowIfNull(mesh);
        ArgumentNullException.ThrowIfNull(coordinates);
        ArgumentNullException.ThrowIfNull(path);

        int numNodes = mesh.Count<Node>();
        int numPoints = mesh.Count<Point>();
        int numBars = mesh.Count<Bar2>();
        int numTris = mesh.Count<Tri3>();
        int numTets = mesh.Count<Tet4>();
        int numElements = numPoints + numBars + numTris + numTets;

        // Validate tag array lengths
        if (pointTags != null && pointTags.Length != numPoints)
            throw new ArgumentException($"pointTags length ({pointTags.Length}) must match Point count ({numPoints})");
        if (barTags != null && barTags.Length != numBars)
            throw new ArgumentException($"barTags length ({barTags.Length}) must match Bar2 count ({numBars})");
        if (triTags != null && triTags.Length != numTris)
            throw new ArgumentException($"triTags length ({triTags.Length}) must match Tri3 count ({numTris})");
        if (tetTags != null && tetTags.Length != numTets)
            throw new ArgumentException($"tetTags length ({tetTags.Length}) must match Tet4 count ({numTets})");

        using var writer = new StreamWriter(path);
        
        // Header
        writer.WriteLine("$MeshFormat");
        writer.WriteLine("2.2 0 8");
        writer.WriteLine("$EndMeshFormat");

        // Physical names section
        if (physicalNames != null && physicalNames.Count > 0)
        {
            writer.WriteLine("$PhysicalNames");
            writer.WriteLine(physicalNames.Count);
            foreach (var (id, (dim, name)) in physicalNames.OrderBy(kv => kv.Key))
                writer.WriteLine($"{dim} {id} \"{name}\"");
            writer.WriteLine("$EndPhysicalNames");
        }

        // Nodes
        writer.WriteLine("$Nodes");
        writer.WriteLine(numNodes);
        for (int i = 0; i < numNodes; i++)
            writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                "{0} {1} {2} {3}", i + 1, coordinates[i, 0], coordinates[i, 1], coordinates[i, 2]));
        writer.WriteLine("$EndNodes");

        // Elements with physical tags
        // Format: elem-id type num-tags physical-tag elementary-tag node1 node2 ...
        writer.WriteLine("$Elements");
        writer.WriteLine(numElements);
        
        int elemId = 1;
        // Gmsh element types: 15=point, 1=line, 2=triangle, 4=tetrahedron
        
        for (int i = 0; i < numPoints; i++)
        {
            var n = mesh.NodesOf<Point, Node>(i);
            int physTag = pointTags?[i] ?? 0;
            writer.WriteLine($"{elemId++} 15 2 {physTag} 0 {n[0] + 1}");
        }
        for (int i = 0; i < numBars; i++)
        {
            var n = mesh.NodesOf<Bar2, Node>(i);
            int physTag = barTags?[i] ?? 0;
            writer.WriteLine($"{elemId++} 1 2 {physTag} 0 {n[0] + 1} {n[1] + 1}");
        }
        for (int i = 0; i < numTris; i++)
        {
            var n = mesh.NodesOf<Tri3, Node>(i);
            int physTag = triTags?[i] ?? 0;
            writer.WriteLine($"{elemId++} 2 2 {physTag} 0 {n[0] + 1} {n[1] + 1} {n[2] + 1}");
        }
        for (int i = 0; i < numTets; i++)
        {
            var n = mesh.NodesOf<Tet4, Node>(i);
            int physTag = tetTags?[i] ?? 0;
            writer.WriteLine($"{elemId++} 4 2 {physTag} 0 {n[0] + 1} {n[1] + 1} {n[2] + 1} {n[3] + 1}");
        }
        
        writer.WriteLine("$EndElements");
    }

    /// <summary>
    /// Saves a SimplexMesh to Gmsh MSH format with cracked/non-cracked element classification.
    /// Convenience overload that creates physical groups based on a set of cracked element indices.
    /// </summary>
    /// <param name="mesh">Mesh to save</param>
    /// <param name="coordinates">Node coordinates</param>
    /// <param name="path">Output file path</param>
    /// <param name="crackedTriangles">Set of cracked Tri3 element indices (0-based)</param>
    /// <param name="crackedTetrahedra">Set of cracked Tet4 element indices (0-based)</param>
    /// <param name="intactGroupId">Tag for intact elements (default: 1)</param>
    /// <param name="crackedGroupId">Tag for cracked elements (default: 2)</param>
    public static void SaveMSHWithCrackGroups(SimplexMesh mesh, double[,] coordinates, string path,
        HashSet<int>? crackedTriangles = null,
        HashSet<int>? crackedTetrahedra = null,
        int intactGroupId = 1,
        int crackedGroupId = 2)
    {
        ArgumentNullException.ThrowIfNull(mesh);
        ArgumentNullException.ThrowIfNull(coordinates);
        ArgumentNullException.ThrowIfNull(path);

        int numTris = mesh.Count<Tri3>();
        int numTets = mesh.Count<Tet4>();

        // Build tag arrays
        int[]? triTags = null;
        int[]? tetTags = null;

        if (crackedTriangles != null || numTris > 0)
        {
            triTags = new int[numTris];
            for (int i = 0; i < numTris; i++)
                triTags[i] = (crackedTriangles?.Contains(i) == true) ? crackedGroupId : intactGroupId;
        }

        if (crackedTetrahedra != null || numTets > 0)
        {
            tetTags = new int[numTets];
            for (int i = 0; i < numTets; i++)
                tetTags[i] = (crackedTetrahedra?.Contains(i) == true) ? crackedGroupId : intactGroupId;
        }

        // Determine dimensionality for physical names
        int dim = numTets > 0 ? 3 : 2;
        
        var physicalNames = new Dictionary<int, (int Dimension, string Name)>
        {
            { intactGroupId, (dim, "NonCracked") },
            { crackedGroupId, (dim, "Cracked") }
        };

        SaveMSH(mesh, coordinates, path, physicalNames,
            pointTags: null, barTags: null, triTags: triTags, tetTags: tetTags);
    }

    /// <summary>
    /// Saves mesh to simple ASCII format for debugging.
    /// </summary>
    public static void SaveASCII(SimplexMesh mesh, double[,] coordinates, string path)
    {
        ArgumentNullException.ThrowIfNull(mesh);
        ArgumentNullException.ThrowIfNull(coordinates);
        ArgumentNullException.ThrowIfNull(path);

        int numNodes = mesh.Count<Node>();
        int numTris = mesh.Count<Tri3>();
        int numTets = mesh.Count<Tet4>();

        using var writer = new StreamWriter(path);
        
        writer.WriteLine($"# SimplexMesh ASCII format");
        writer.WriteLine($"# Nodes: {numNodes}");
        writer.WriteLine($"# Tri3: {numTris}");
        writer.WriteLine($"# Tet4: {numTets}");
        writer.WriteLine();

        writer.WriteLine("COORDINATES");
        writer.WriteLine(numNodes);
        for (int i = 0; i < numNodes; i++)
            writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                "{0} {1} {2}", coordinates[i, 0], coordinates[i, 1], coordinates[i, 2]));
        writer.WriteLine();

        if (numTris > 0)
        {
            writer.WriteLine("TRIANGLES");
            writer.WriteLine(numTris);
            for (int i = 0; i < numTris; i++)
            {
                var n = mesh.NodesOf<Tri3, Node>(i);
                writer.WriteLine($"{n[0]} {n[1]} {n[2]}");
            }
            writer.WriteLine();
        }

        if (numTets > 0)
        {
            writer.WriteLine("TETRAHEDRA");
            writer.WriteLine(numTets);
            for (int i = 0; i < numTets; i++)
            {
                var n = mesh.NodesOf<Tet4, Node>(i);
                writer.WriteLine($"{n[0]} {n[1]} {n[2]} {n[3]}");
            }
        }
    }

    /// <summary>
    /// Prints mesh statistics to console for debugging.
    /// </summary>
    public static void PrintStats(SimplexMesh mesh, string label = "Mesh")
    {
        Console.WriteLine($"=== {label} ===");
        Console.WriteLine($"  Nodes: {mesh.Count<Node>()}");
        Console.WriteLine($"  Edges: {mesh.Count<Edge>()}");
        Console.WriteLine($"  Points: {mesh.Count<Point>()}");
        Console.WriteLine($"  Bar2: {mesh.Count<Bar2>()}");
        Console.WriteLine($"  Tri3: {mesh.Count<Tri3>()}");
        Console.WriteLine($"  Tet4: {mesh.Count<Tet4>()}");
    }
    // ═══════════════════════════════════════════════════════════════════════════
    // REMOVED: InterpolateCoordinates() Function (now in MeshRefinement.cs)
    // ═══════════════════════════════════════════════════════════════════════════
    // The InterpolateCoordinates() function has been REMOVED to eliminate duplication.
    //
    // Use instead: MeshRefinement.InterpolateCoordinates()
    //
    // Example usage:
    //   using static Numerical.MeshRefinement;
    //   var newCoords = InterpolateCoordinates(refinedMesh, oldCoords);
    //
    // The refactored version is identical in functionality.
    // ═══════════════════════════════════════════════════════════════════════════
    
    public static (SimplexMesh Mesh, double[,] Coordinates) LoadMSH(string path)
    {
        ArgumentNullException.ThrowIfNull(path);
        if (!File.Exists(path))
            throw new FileNotFoundException($"MSH file not found: {path}");

        using var reader = new StreamReader(path);
        return ParseMSH(reader);
    }

    /// <summary>
    /// Loads a SimplexMesh from a Gmsh MSH file, also returning physical tags per element type.
    /// </summary>
    /// <param name="path">Path to the .msh file</param>
    /// <returns>Loaded mesh, coordinates, and physical tags for each element type</returns>
    public static (SimplexMesh Mesh, double[,] Coordinates, 
        Dictionary<int, (int Dimension, string Name)> PhysicalNames,
        int[] PointTags, int[] BarTags, int[] TriTags, int[] TetTags) LoadMSHWithTags(string path)
    {
        ArgumentNullException.ThrowIfNull(path);
        if (!File.Exists(path))
            throw new FileNotFoundException($"MSH file not found: {path}");

        using var reader = new StreamReader(path);
        return ParseMSHWithTags(reader);
    }

    private static (SimplexMesh Mesh, double[,] Coordinates) ParseMSH(StreamReader reader)
    {
        var (mesh, coords, _, _, _, _, _) = ParseMSHWithTags(reader);
        return (mesh, coords);
    }

    private static (SimplexMesh Mesh, double[,] Coordinates,
        Dictionary<int, (int Dimension, string Name)> PhysicalNames,
        int[] PointTags, int[] BarTags, int[] TriTags, int[] TetTags) ParseMSHWithTags(StreamReader reader)
    {
        var physicalNames = new Dictionary<int, (int Dimension, string Name)>();
        double[,]? coordinates = null;
        var elements = new List<(int ElemType, int PhysTag, int[] Nodes)>();

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            line = line.Trim();
            
            if (line == "$MeshFormat")
            {
                line = reader.ReadLine()?.Trim();
                if (line != null)
                {
                    var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    double version = double.Parse(parts[0], CultureInfo.InvariantCulture);
                    if (version < 2.0 || version >= 3.0)
                        Console.WriteLine($"Warning: Expected Gmsh 2.x format, got {version}");
                }
                SkipToEndSection(reader, "$EndMeshFormat");
            }
            else if (line == "$PhysicalNames")
            {
                line = reader.ReadLine()?.Trim();
                int numNames = int.Parse(line!, CultureInfo.InvariantCulture);
                for (int i = 0; i < numNames; i++)
                {
                    line = reader.ReadLine()?.Trim();
                    if (line == null) break;
                    var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    int dim = int.Parse(parts[0], CultureInfo.InvariantCulture);
                    int id = int.Parse(parts[1], CultureInfo.InvariantCulture);
                    string name = parts[2].Trim('"');
                    physicalNames[id] = (dim, name);
                }
                SkipToEndSection(reader, "$EndPhysicalNames");
            }
            else if (line == "$Nodes")
            {
                line = reader.ReadLine()?.Trim();
                int nodeCount = int.Parse(line!, CultureInfo.InvariantCulture);
                coordinates = new double[nodeCount, 3];
                
                for (int i = 0; i < nodeCount; i++)
                {
                    line = reader.ReadLine()?.Trim();
                    if (line == null) break;
                    var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    int nodeId = int.Parse(parts[0], CultureInfo.InvariantCulture) - 1; // 0-based
                    coordinates[nodeId, 0] = double.Parse(parts[1], CultureInfo.InvariantCulture);
                    coordinates[nodeId, 1] = double.Parse(parts[2], CultureInfo.InvariantCulture);
                    coordinates[nodeId, 2] = double.Parse(parts[3], CultureInfo.InvariantCulture);
                }
                SkipToEndSection(reader, "$EndNodes");
            }
            else if (line == "$Elements")
            {
                line = reader.ReadLine()?.Trim();
                int numElements = int.Parse(line!, CultureInfo.InvariantCulture);
                
                for (int i = 0; i < numElements; i++)
                {
                    line = reader.ReadLine()?.Trim();
                    if (line == null) break;
                    var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    // Format: elem-id type num-tags [tags...] node1 node2 ...
                    int elemType = int.Parse(parts[1], CultureInfo.InvariantCulture);
                    int numTags = int.Parse(parts[2], CultureInfo.InvariantCulture);
                    int physTag = numTags > 0 ? int.Parse(parts[3], CultureInfo.InvariantCulture) : 0;
                    
                    int nodeStart = 3 + numTags;
                    var nodes = new int[parts.Length - nodeStart];
                    for (int j = 0; j < nodes.Length; j++)
                        nodes[j] = int.Parse(parts[nodeStart + j], CultureInfo.InvariantCulture) - 1; // 0-based
                    
                    elements.Add((elemType, physTag, nodes));
                }
                SkipToEndSection(reader, "$EndElements");
            }
            else if (line.StartsWith("$"))
            {
                // Skip unknown sections
                string endMarker = "$End" + line.Substring(1);
                SkipToEndSection(reader, endMarker);
            }
        }

        if (coordinates == null)
            throw new FormatException("MSH file missing $Nodes section");

        // Build mesh and collect tags by element type
        // Gmsh types: 15=point, 1=line, 2=triangle, 4=tetrahedron
        var pointList = new List<(int PhysTag, int[] Nodes)>();
        var barList = new List<(int PhysTag, int[] Nodes)>();
        var triList = new List<(int PhysTag, int[] Nodes)>();
        var tetList = new List<(int PhysTag, int[] Nodes)>();

        foreach (var (elemType, physTag, nodes) in elements)
        {
            switch (elemType)
            {
                case 15: // Point
                    pointList.Add((physTag, nodes));
                    break;
                case 1: // 2-node line
                case 8: // 3-node line (quadratic) - treat as linear
                    barList.Add((physTag, new[] { nodes[0], nodes[1] }));
                    break;
                case 2: // 3-node triangle
                case 9: // 6-node triangle (quadratic) - treat as linear
                    triList.Add((physTag, new[] { nodes[0], nodes[1], nodes[2] }));
                    break;
                case 4: // 4-node tetrahedron
                case 11: // 10-node tetrahedron (quadratic) - treat as linear
                    tetList.Add((physTag, new[] { nodes[0], nodes[1], nodes[2], nodes[3] }));
                    break;
            }
        }

        var mesh = new SimplexMesh();
        int numNodes = coordinates.GetLength(0);

        mesh.WithBatch(() =>
        {
            for (int i = 0; i < numNodes; i++)
            {
                mesh.Add<Node>();
                mesh.Set<Node, ParentNodes>(i, new ParentNodes(i, i));
            }

            foreach (var (_, nodes) in pointList)
            {
                int pi = mesh.Add<Point, Node>(nodes[0]);
                mesh.Set<Point, OriginalElement>(pi, new OriginalElement(pi));
            }
            foreach (var (_, nodes) in barList)
            {
                int bi = mesh.Add<Bar2, Node>(nodes[0], nodes[1]);
                mesh.Set<Bar2, OriginalElement>(bi, new OriginalElement(bi));
            }
            foreach (var (_, nodes) in triList)
            {
                int ti = mesh.Add<Tri3, Node>(nodes[0], nodes[1], nodes[2]);
                mesh.Set<Tri3, OriginalElement>(ti, new OriginalElement(ti));
            }
            foreach (var (_, nodes) in tetList)
            {
                int ei = mesh.Add<Tet4, Node>(nodes[0], nodes[1], nodes[2], nodes[3]);
                mesh.Set<Tet4, OriginalElement>(ei, new OriginalElement(ei));
            }
        });

        return (mesh, coordinates, physicalNames,
            pointList.Select(x => x.PhysTag).ToArray(),
            barList.Select(x => x.PhysTag).ToArray(),
            triList.Select(x => x.PhysTag).ToArray(),
            tetList.Select(x => x.PhysTag).ToArray());
    }

    private static void SkipToEndSection(StreamReader reader, string endMarker)
    {
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            if (line.Trim() == endMarker) break;
        }
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // GiD/CIMNE Mesh I/O
    // ═══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Loads a SimplexMesh from a GiD/CIMNE mesh file (.msh).
    /// </summary>
    /// <param name="path">Path to the .msh file</param>
    /// <returns>Loaded mesh and coordinates</returns>
    public static (SimplexMesh Mesh, double[,] Coordinates) LoadGiD(string path)
    {
        ArgumentNullException.ThrowIfNull(path);
        if (!File.Exists(path))
            throw new FileNotFoundException($"GiD mesh file not found: {path}");

        using var reader = new StreamReader(path);
        return ParseGiD(reader);
    }

    /// <summary>
    /// Loads a SimplexMesh from a GiD/CIMNE mesh file, also returning material IDs.
    /// </summary>
    /// <param name="path">Path to the .msh file</param>
    /// <returns>Loaded mesh, coordinates, and material IDs for each element type</returns>
    public static (SimplexMesh Mesh, double[,] Coordinates,
        int[] PointMaterials, int[] BarMaterials, int[] TriMaterials, int[] TetMaterials) LoadGiDWithMaterials(string path)
    {
        ArgumentNullException.ThrowIfNull(path);
        if (!File.Exists(path))
            throw new FileNotFoundException($"GiD mesh file not found: {path}");

        using var reader = new StreamReader(path);
        return ParseGiDWithMaterials(reader);
    }

    private static (SimplexMesh Mesh, double[,] Coordinates) ParseGiD(StreamReader reader)
    {
        var (mesh, coords, _, _, _, _) = ParseGiDWithMaterials(reader);
        return (mesh, coords);
    }

    private static (SimplexMesh Mesh, double[,] Coordinates,
        int[] PointMaterials, int[] BarMaterials, int[] TriMaterials, int[] TetMaterials) ParseGiDWithMaterials(StreamReader reader)
    {
        var nodeCoords = new Dictionary<int, (double X, double Y, double Z)>();
        var pointList = new List<(int Material, int[] Nodes)>();
        var barList = new List<(int Material, int[] Nodes)>();
        var triList = new List<(int Material, int[] Nodes)>();
        var tetList = new List<(int Material, int[] Nodes)>();

        int dimension = 3;
        string currentElemType = "";
        int currentNnode = 0;
        bool inCoordinates = false;
        bool inElements = false;

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            line = line.Trim();
            if (string.IsNullOrEmpty(line) || line.StartsWith("#")) continue;

            var upperLine = line.ToUpperInvariant();

            // Parse MESH header: MESH dimension D ElemType TypeName Nnode N
            if (upperLine.StartsWith("MESH"))
            {
                var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                for (int i = 0; i < parts.Length - 1; i++)
                {
                    var key = parts[i].ToUpperInvariant();
                    if (key == "DIMENSION" && i + 1 < parts.Length)
                        dimension = int.Parse(parts[i + 1], CultureInfo.InvariantCulture);
                    else if (key == "ELEMTYPE" && i + 1 < parts.Length)
                        currentElemType = parts[i + 1].ToUpperInvariant();
                    else if (key == "NNODE" && i + 1 < parts.Length)
                        currentNnode = int.Parse(parts[i + 1], CultureInfo.InvariantCulture);
                }
                continue;
            }

            if (upperLine == "COORDINATES")
            {
                inCoordinates = true;
                inElements = false;
                continue;
            }
            if (upperLine == "END COORDINATES")
            {
                inCoordinates = false;
                continue;
            }
            if (upperLine == "ELEMENTS")
            {
                inElements = true;
                inCoordinates = false;
                continue;
            }
            if (upperLine == "END ELEMENTS")
            {
                inElements = false;
                continue;
            }

            if (inCoordinates)
            {
                var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length >= 3)
                {
                    int nodeId = int.Parse(parts[0], CultureInfo.InvariantCulture);
                    double x = double.Parse(parts[1], CultureInfo.InvariantCulture);
                    double y = double.Parse(parts[2], CultureInfo.InvariantCulture);
                    double z = parts.Length > 3 ? double.Parse(parts[3], CultureInfo.InvariantCulture) : 0.0;
                    nodeCoords[nodeId] = (x, y, z);
                }
            }
            else if (inElements)
            {
                var parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length < 2) continue;

                // Format: elem-id node1 node2 ... [material]
                // Node count determined by currentNnode or element type
                int nNodes = currentNnode > 0 ? currentNnode : GuessNodeCount(currentElemType);
                if (parts.Length < 1 + nNodes) continue;

                var nodes = new int[nNodes];
                for (int j = 0; j < nNodes; j++)
                    nodes[j] = int.Parse(parts[1 + j], CultureInfo.InvariantCulture);
                
                // Material is the last column (optional)
                int material = parts.Length > 1 + nNodes 
                    ? int.Parse(parts[1 + nNodes], CultureInfo.InvariantCulture) 
                    : 0;

                // Classify by element type
                switch (currentElemType)
                {
                    case "POINT":
                        pointList.Add((material, nodes));
                        break;
                    case "LINEAR":
                    case "LINE":
                        barList.Add((material, new[] { nodes[0], nodes[1] }));
                        break;
                    case "TRIANGLE":
                        triList.Add((material, new[] { nodes[0], nodes[1], nodes[2] }));
                        break;
                    case "TETRAHEDRA":
                    case "TETRAHEDRON":
                        tetList.Add((material, new[] { nodes[0], nodes[1], nodes[2], nodes[3] }));
                        break;
                    default:
                        // Guess from node count
                        if (nNodes == 1) pointList.Add((material, nodes));
                        else if (nNodes == 2) barList.Add((material, nodes));
                        else if (nNodes == 3) triList.Add((material, nodes));
                        else if (nNodes == 4) tetList.Add((material, nodes));
                        break;
                }
            }
        }

        if (nodeCoords.Count == 0)
            throw new FormatException("GiD mesh file has no coordinates");

        // Build node array (GiD uses 1-based, need to map to 0-based)
        int maxNodeId = nodeCoords.Keys.Max();
        int minNodeId = nodeCoords.Keys.Min();
        int numNodes = nodeCoords.Count;
        var coordinates = new double[numNodes, 3];
        var nodeIdMap = new Dictionary<int, int>(); // old id -> new 0-based id

        int idx = 0;
        foreach (var nodeId in nodeCoords.Keys.OrderBy(x => x))
        {
            var (x, y, z) = nodeCoords[nodeId];
            coordinates[idx, 0] = x;
            coordinates[idx, 1] = y;
            coordinates[idx, 2] = z;
            nodeIdMap[nodeId] = idx;
            idx++;
        }

        // Remap element node IDs
        void RemapNodes(int[] nodes)
        {
            for (int i = 0; i < nodes.Length; i++)
                nodes[i] = nodeIdMap[nodes[i]];
        }

        foreach (var (_, nodes) in pointList) RemapNodes(nodes);
        foreach (var (_, nodes) in barList) RemapNodes(nodes);
        foreach (var (_, nodes) in triList) RemapNodes(nodes);
        foreach (var (_, nodes) in tetList) RemapNodes(nodes);

        // Build mesh
        var mesh = new SimplexMesh();
        mesh.WithBatch(() =>
        {
            for (int i = 0; i < numNodes; i++)
            {
                mesh.Add<Node>();
                mesh.Set<Node, ParentNodes>(i, new ParentNodes(i, i));
            }

            foreach (var (_, nodes) in pointList)
            {
                int pi = mesh.Add<Point, Node>(nodes[0]);
                mesh.Set<Point, OriginalElement>(pi, new OriginalElement(pi));
            }
            foreach (var (_, nodes) in barList)
            {
                int bi = mesh.Add<Bar2, Node>(nodes[0], nodes[1]);
                mesh.Set<Bar2, OriginalElement>(bi, new OriginalElement(bi));
            }
            foreach (var (_, nodes) in triList)
            {
                int ti = mesh.Add<Tri3, Node>(nodes[0], nodes[1], nodes[2]);
                mesh.Set<Tri3, OriginalElement>(ti, new OriginalElement(ti));
            }
            foreach (var (_, nodes) in tetList)
            {
                int ei = mesh.Add<Tet4, Node>(nodes[0], nodes[1], nodes[2], nodes[3]);
                mesh.Set<Tet4, OriginalElement>(ei, new OriginalElement(ei));
            }
        });

        return (mesh, coordinates,
            pointList.Select(x => x.Material).ToArray(),
            barList.Select(x => x.Material).ToArray(),
            triList.Select(x => x.Material).ToArray(),
            tetList.Select(x => x.Material).ToArray());
    }

    private static int GuessNodeCount(string elemType) => elemType switch
    {
        "POINT" => 1,
        "LINEAR" or "LINE" => 2,
        "TRIANGLE" => 3,
        "QUADRILATERAL" => 4,
        "TETRAHEDRA" or "TETRAHEDRON" => 4,
        "HEXAHEDRA" or "HEXAHEDRON" => 8,
        "PRISM" => 6,
        "PYRAMID" => 5,
        _ => 0
    };

    /// <summary>
    /// Saves a SimplexMesh to GiD/CIMNE mesh format (.msh).
    /// </summary>
    /// <param name="mesh">Mesh to save</param>
    /// <param name="coordinates">Node coordinates</param>
    /// <param name="path">Output file path</param>
    public static void SaveGiD(SimplexMesh mesh, double[,] coordinates, string path)
    {
        SaveGiD(mesh, coordinates, path, null, null, null, null, null);
    }

    /// <summary>
    /// Saves a SimplexMesh to GiD/CIMNE mesh format (.msh) with material IDs.
    /// GUARANTEES continuous node and element numbering with no gaps.
    /// </summary>
    /// <param name="mesh">Mesh to save</param>
    /// <param name="coordinates">Node coordinates</param>
    /// <param name="path">Output file path</param>
    /// <param name="pointMaterials">Optional material IDs for Point elements</param>
    /// <param name="barMaterials">Optional material IDs for Bar2 elements</param>
    /// <param name="triMaterials">Optional material IDs for Tri3 elements</param>
    /// <param name="tetMaterials">Optional material IDs for Tet4 elements</param>
    /// <param name="meshName">Optional mesh name for GiD header</param>
    public static void SaveGiD(SimplexMesh mesh, double[,] coordinates, string path,
        int[]? pointMaterials,
        int[]? barMaterials,
        int[]? triMaterials,
        int[]? quadMaterials,
        int[]? tetMaterials,
        string? meshName = null)
    {
        ArgumentNullException.ThrowIfNull(mesh);
        ArgumentNullException.ThrowIfNull(coordinates);
        ArgumentNullException.ThrowIfNull(path);

        int numPoints = mesh.Count<Point>();
        int numBars = mesh.Count<Bar2>();
        int numTris = mesh.Count<Tri3>();
        int numQuads = mesh.Count<Quad4>();
        int numTets = mesh.Count<Tet4>();
        int numElements = numPoints + numBars + numTris + numQuads + numTets;

        // Validate material array lengths
        if (pointMaterials != null && pointMaterials.Length != numPoints)
            throw new ArgumentException($"pointMaterials length ({pointMaterials.Length}) must match Point count ({numPoints})");
        if (barMaterials != null && barMaterials.Length != numBars)
            throw new ArgumentException($"barMaterials length ({barMaterials.Length}) must match Bar2 count ({numBars})");
        if (triMaterials != null && triMaterials.Length != numTris)
            throw new ArgumentException($"triMaterials length ({triMaterials.Length}) must match Tri3 count ({numTris})");
        if (quadMaterials != null && quadMaterials.Length != numQuads)
            throw new ArgumentException($"quadMaterials length ({quadMaterials.Length}) must match Quad4 count ({numQuads})");
        if (tetMaterials != null && tetMaterials.Length != numTets)
            throw new ArgumentException($"tetMaterials length ({tetMaterials.Length}) must match Tet4 count ({numTets})");

        // STEP 1: Collect all active nodes (nodes referenced by elements)
        var activeNodes = new HashSet<int>();
        
        for (int i = 0; i < numPoints; i++)
        {
            var n = mesh.NodesOf<Point, Node>(i);
            activeNodes.Add(n[0]);
        }
        
        for (int i = 0; i < numBars; i++)
        {
            var n = mesh.NodesOf<Bar2, Node>(i);
            activeNodes.Add(n[0]);
            activeNodes.Add(n[1]);
        }
        
        for (int i = 0; i < numTris; i++)
        {
            var n = mesh.NodesOf<Tri3, Node>(i);
            activeNodes.Add(n[0]);
            activeNodes.Add(n[1]);
            activeNodes.Add(n[2]);
        }
        
        for (int i = 0; i < numQuads; i++)
        {
            var n = mesh.NodesOf<Quad4, Node>(i);
            activeNodes.Add(n[0]);
            activeNodes.Add(n[1]);
            activeNodes.Add(n[2]);
            activeNodes.Add(n[3]);
        }
        
        for (int i = 0; i < numTets; i++)
        {
            var n = mesh.NodesOf<Tet4, Node>(i);
            activeNodes.Add(n[0]);
            activeNodes.Add(n[1]);
            activeNodes.Add(n[2]);
            activeNodes.Add(n[3]);
        }
        
        // STEP 2: Create continuous node numbering: oldID -> newID (1, 2, 3, ...)
        var nodeMap = new Dictionary<int, int>();
        var sortedNodes = activeNodes.OrderBy(n => n).ToArray();
        
        for (int i = 0; i < sortedNodes.Length; i++)
        {
            nodeMap[sortedNodes[i]] = i + 1;  // GiD uses 1-based indexing
        }
        
        int numActiveNodes = sortedNodes.Length;
        
        Console.WriteLine($"[SaveGiD] Mesh has {mesh.Count<Node>()} node IDs, {numActiveNodes} active nodes");
        if (numActiveNodes < mesh.Count<Node>())
        {
            Console.WriteLine($"[SaveGiD] Removed {mesh.Count<Node>() - numActiveNodes} unused nodes (gaps eliminated)");
        }

        // Determine mesh dimension
        int dimension = numTets > 0 ? 3 : (numTris > 0 ? 2 : 1);

        // Ensure path has .msh extension
        if (!path.EndsWith(".msh", StringComparison.OrdinalIgnoreCase))
        {
            path = path + ".msh";
        }

        using var writer = new StreamWriter(path);
        
        // Write mesh information header
        writer.WriteLine("# GiD/CIMNE Mesh File");
        writer.WriteLine($"# Generated by SimplexRemesher with continuous numbering");
        if (!string.IsNullOrEmpty(meshName))
            writer.WriteLine($"# Mesh name: {meshName}");
        writer.WriteLine($"# Dimension: {dimension}");
        writer.WriteLine($"# Nodes: {numActiveNodes} (continuous 1-based)");
        writer.WriteLine($"# Elements: {numElements} total");
        writer.WriteLine("#   Tets: {0}, Quads: {1}, Tris: {2}, Bars: {3}, Points: {4}",
            numTets, numQuads, numTris, numBars, numPoints);
        writer.WriteLine();
        
        int currentElementId = 1; // Track continuous element numbering across types

        // Write bulk elements first (highest dimension), then lower-dimensional
        // GiD groups elements by type in separate MESH blocks
        // Each MESH block MUST have its own Coordinates section

        // Tetrahedra
        if (numTets > 0)
        {
            writer.WriteLine($"MESH dimension {dimension} ElemType Tetrahedra Nnode 4");
            WriteGiDCoordinatesWithMapping(writer, coordinates, sortedNodes, dimension);
            writer.WriteLine("Elements");
            for (int i = 0; i < numTets; i++)
            {
                var n = mesh.NodesOf<Tet4, Node>(i);
                int mat = tetMaterials?[i] ?? 0;
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "  {0}  {1}  {2}  {3}  {4}  {5}", 
                    currentElementId++,  // Continuous numbering
                    nodeMap[n[0]], nodeMap[n[1]], nodeMap[n[2]], nodeMap[n[3]], 
                    mat));
            }
            writer.WriteLine("End Elements");
            writer.WriteLine();
        }

        // Quadrilaterals
        if (numQuads > 0)
        {
            writer.WriteLine($"MESH dimension {dimension} ElemType Quadrilateral Nnode 4");
            WriteGiDCoordinatesWithMapping(writer, coordinates, sortedNodes, dimension);
            writer.WriteLine("Elements");
            for (int i = 0; i < numQuads; i++)
            {
                var n = mesh.NodesOf<Quad4, Node>(i);
                int mat = quadMaterials?[i] ?? 0;
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "  {0}  {1}  {2}  {3}  {4}  {5}", 
                    currentElementId++,  // Continuous numbering
                    nodeMap[n[0]], nodeMap[n[1]], nodeMap[n[2]], nodeMap[n[3]], 
                    mat));
            }
            writer.WriteLine("End Elements");
            writer.WriteLine();
        }

        // Triangles
        if (numTris > 0)
        {
            writer.WriteLine($"MESH dimension {dimension} ElemType Triangle Nnode 3");
            WriteGiDCoordinatesWithMapping(writer, coordinates, sortedNodes, dimension);
            writer.WriteLine("Elements");
            for (int i = 0; i < numTris; i++)
            {
                var n = mesh.NodesOf<Tri3, Node>(i);
                int mat = triMaterials?[i] ?? 0;
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "  {0}  {1}  {2}  {3}  {4}", 
                    currentElementId++,  // Continuous numbering
                    nodeMap[n[0]], nodeMap[n[1]], nodeMap[n[2]], 
                    mat));
            }
            writer.WriteLine("End Elements");
            writer.WriteLine();
        }

        // Lines/Bars
        if (numBars > 0)
        {
            writer.WriteLine($"MESH dimension {dimension} ElemType Linear Nnode 2");
            WriteGiDCoordinatesWithMapping(writer, coordinates, sortedNodes, dimension);
            writer.WriteLine("Elements");
            for (int i = 0; i < numBars; i++)
            {
                var n = mesh.NodesOf<Bar2, Node>(i);
                int mat = barMaterials?[i] ?? 0;
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "  {0}  {1}  {2}  {3}", 
                    currentElementId++,  // Continuous numbering
                    nodeMap[n[0]], nodeMap[n[1]], 
                    mat));
            }
            writer.WriteLine("End Elements");
            writer.WriteLine();
        }

        // Points
        if (numPoints > 0)
        {
            writer.WriteLine($"MESH dimension {dimension} ElemType Point Nnode 1");
            WriteGiDCoordinatesWithMapping(writer, coordinates, sortedNodes, dimension);
            writer.WriteLine("Elements");
            for (int i = 0; i < numPoints; i++)
            {
                var n = mesh.NodesOf<Point, Node>(i);
                int mat = pointMaterials?[i] ?? 0;
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "  {0}  {1}  {2}", 
                    currentElementId++,  // Continuous numbering
                    nodeMap[n[0]], 
                    mat));
            }
            writer.WriteLine("End Elements");
            writer.WriteLine();
        }
        
        Console.WriteLine($"[SaveGiD] ✓ Saved to {path} with continuous numbering (no gaps)");
    }

    private static void WriteGiDCoordinatesWithMapping(StreamWriter writer, double[,] coordinates, int[] sortedNodes, int dimension)
    {
        writer.WriteLine("Coordinates");
        for (int i = 0; i < sortedNodes.Length; i++)
        {
            int oldNodeId = sortedNodes[i];
            int newNodeId = i + 1;  // GiD uses 1-based indexing
            
            if (dimension == 2)
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "  {0}  {1}  {2}", newNodeId, coordinates[oldNodeId, 0], coordinates[oldNodeId, 1]));
            else
                writer.WriteLine(string.Format(CultureInfo.InvariantCulture,
                    "  {0}  {1}  {2}  {3}", newNodeId, coordinates[oldNodeId, 0], coordinates[oldNodeId, 1], coordinates[oldNodeId, 2]));
        }
        writer.WriteLine("End Coordinates");
        writer.WriteLine();
    }
    
    // ========================================================================
    // SIGNED FIELD CRACK INSERTION
    // ========================================================================
    
    /// <summary>
    /// Signed field function delegate for level set crack insertion.
    /// Returns: positive on one side, negative on other side, zero on crack.
    /// </summary>
    /// <example>
    /// SINGLE LEVEL SET EXAMPLES:
    /// Horizontal crack at y=0.5: (x,y,z) => y - 0.5
    /// Circular crack: (x,y,z) => Math.Sqrt((x-cx)^2 + (y-cy)^2) - radius
    /// Planar 3D crack: (x,y,z) => z - 0.5
    /// 
    /// TWO LEVEL SET EXAMPLES:
    /// Level set 1 defines crack surface, level set 2 defines active region.
    /// Crack forms where surface crosses zero AND region is ≤ 0.
    /// 
    /// 2D horizontal crack from x=0.2 to x=1.5:
    ///   crackSurface = (x,y,z) => y - 0.5
    ///   activeRegion = (x,y,z) => (x < 0.2) ? (x - 0.2) : (x > 1.5) ? (x - 1.5) : -0.1
    /// 
    /// 3D circular crack (plane z=0.5, radius 0.3 at center):
    ///   crackSurface = (x,y,z) => z - 0.5
    ///   activeRegion = (x,y,z) => Math.Sqrt((x-0.5)^2 + (y-0.5)^2) - 0.3
    /// </example>
    public delegate double SignedFieldFunction(double x, double y, double z);
    
    /// <summary>
    /// Create crack in mesh using signed field (level set method).
    /// 
    /// Algorithm:
    /// 1. Find edges where field changes sign
    /// 2. Refine those edges, positioning new nodes at exact zero-crossing
    /// 3. Classify original nodes by field sign (positive/negative)
    /// 4. Duplicate all crack nodes except tips
    /// 5. Elements with all positive original nodes use duplicates
    /// 
    /// Works for all crack types:
    /// - Boundary-to-boundary
    /// - Edge-to-interior (with singular tip)
    /// - Closed loops
    /// - Any shape defined by implicit function
    /// </summary>
    /// <param name="mesh">Input mesh</param>
    /// <param name="coordinates">Node coordinates</param>
    /// <param name="signedField">Signed distance/level set function (defines crack surface)</param>
    /// <param name="regionField">Optional second level set to define crack active region (crack exists where regionField ≤ 0). If null, crack extends along entire zero level set.</param>
    /// <param name="enableSmoothing">Enable CVT smoothing after refinement (disable for multi-crack patterns)</param>
    /// <returns>Cracked mesh and coordinates</returns>
    
    /// <summary>
    /// Helper to safely get z-coordinate from 2D or 3D coordinate array
    /// </summary>
    private static double GetZ(double[,] coords, int nodeId)
    {
        return coords.GetLength(1) == 3 ? coords[nodeId, 2] : 0.0;
    }
    
    public static (SimplexMesh, double[,]) CreateCrackFromSignedField(
        SimplexMesh mesh,
        double[,] coordinates,
        SignedFieldFunction signedField,
        SignedFieldFunction? regionField = null,
        bool enableMeshPerturbation = false,
        bool enableSmoothing = false,
        double visualizationOffset = 0.0)
    {
        int originalNodeCount = mesh.Count<Node>();
        int coordDim = coordinates.GetLength(1);  // 2 for 2D, 3 for 3D
        
        // Clone coordinates to avoid modifying input
        coordinates = (double[,])coordinates.Clone();
        
        // OPTIONAL: MESH PERTURBATION to avoid exact zeros in signed field
        if (enableMeshPerturbation)
        {
            Console.WriteLine($"  → Applying mesh perturbation...");
            
            // Calculate average edge length
            DiscoverEdges(mesh);
            double avgEdgeLength = 0.0;
            int numEdges = 0;
            for (int i = 0; i < mesh.Count<Edge>(); i++)
            {
                var nodes = mesh.NodesOf<Edge, Node>(i);
                int n1 = nodes[0], n2 = nodes[1];
                double dx = coordinates[n2, 0] - coordinates[n1, 0];
                double dy = coordinates[n2, 1] - coordinates[n1, 1];
                double dz = coordDim == 3 ? coordinates[n2, 2] - coordinates[n1, 2] : 0.0;
                avgEdgeLength += Math.Sqrt(dx * dx + dy * dy + dz * dz);
                numEdges++;
            }
            avgEdgeLength /= numEdges;
            
            // Perturb each node by a tiny random amount (1e-6 * avgEdgeLength)
            double perturbationScale = 1e-6 * avgEdgeLength;
            var random = new Random(12345);  // Fixed seed for reproducibility
            
            for (int i = 0; i < mesh.Count<Node>(); i++)
            {
                coordinates[i, 0] += (random.NextDouble() * 2.0 - 1.0) * perturbationScale;
                coordinates[i, 1] += (random.NextDouble() * 2.0 - 1.0) * perturbationScale;
                if (coordDim == 3)
                    coordinates[i, 2] += (random.NextDouble() * 2.0 - 1.0) * perturbationScale;
            }
            
            Console.WriteLine($"  → Mesh perturbed by ±{perturbationScale:E3}");
        }
        else
        {
            DiscoverEdges(mesh);
        }
        
        // Step 1: Find ALL edges that cross the signed field
        var edgesToRefine = new List<(int, int)>();
        
        Console.WriteLine($"  → Total edges in mesh: {mesh.Count<Edge>()}");
        
        int crossingEdges = 0;
        int regionFiltered = 0;
        
        for (int i = 0; i < mesh.Count<Edge>(); i++)
        {
            var nodes = mesh.NodesOf<Edge, Node>(i);
            int n1 = nodes[0], n2 = nodes[1];
            
            double f1 = signedField(coordinates[n1, 0], coordinates[n1, 1], GetZ(coordinates, n1));
            double f2 = signedField(coordinates[n2, 0], coordinates[n2, 1], GetZ(coordinates, n2));
            
            // Check if crack surface crosses this edge (sign change)
            // With perturbation, f1*f2 == 0 should be extremely rare
            bool crossesCrack = f1 * f2 < 0;
            
            if (crossesCrack)
            {
                crossingEdges++;
                
                // If region field is provided, check if edge is in active region
                if (regionField != null)
                {
                    double r1 = regionField(coordinates[n1, 0], coordinates[n1, 1], GetZ(coordinates, n1));
                    double r2 = regionField(coordinates[n2, 0], coordinates[n2, 1], GetZ(coordinates, n2));
                    
                    // Check if crossing point is in active region
                    // Approximate crossing point using linear interpolation
                    double t = Math.Abs(f2 - f1) > 1e-14 ? -f1 / (f2 - f1) : 0.5;
                    t = Math.Max(0.0, Math.Min(1.0, t));
                    
                    double xCross = coordinates[n1, 0] + t * (coordinates[n2, 0] - coordinates[n1, 0]);
                    double yCross = coordinates[n1, 1] + t * (coordinates[n2, 1] - coordinates[n1, 1]);
                    double zCross = GetZ(coordinates, n1) + t * (GetZ(coordinates, n2) - GetZ(coordinates, n1));
                    double rCross = regionField(xCross, yCross, zCross);
                    
                    // Only refine edge if crossing point is in active region (r ≤ 0)
                    if (rCross <= 0)
                    {
                        edgesToRefine.Add((n1, n2));
                    }
                    else
                    {
                        regionFiltered++;
                    }
                }
                else
                {
                    // No region constraint - refine all crossing edges
                    edgesToRefine.Add((n1, n2));
                }
            }
        }
        
        Console.WriteLine($"  → Edges with sign change: {crossingEdges}");
        Console.WriteLine($"  → Filtered by region: {regionFiltered}");
        Console.WriteLine($"  → Edges to refine: {edgesToRefine.Count}");
        
        if (edgesToRefine.Count == 0)
        {
            Console.WriteLine($"  ⚠️  No edges cross the crack surface in active region");
            return (mesh, coordinates);
        }
        
        Console.WriteLine($"  → Refining {edgesToRefine.Count} edges that cross crack");
        
        // Step 2: Refine those edges once
        var refinedMesh = mesh;
        var refinedCoords = coordinates;
        
        {
            // Cache signed field values for all nodes (avoid recomputation)
            var cachedFieldValues = new double[refinedMesh.Count<Node>()];
            var cachedRegionValues = new double[refinedMesh.Count<Node>()];

            for (int i = 0; i < refinedMesh.Count<Node>(); i++)
            {
                cachedFieldValues[i] = signedField(refinedCoords[i, 0], refinedCoords[i, 1], GetZ(refinedCoords, i));
                if (regionField != null)
                    cachedRegionValues[i] = regionField(refinedCoords[i, 0], refinedCoords[i, 1], GetZ(refinedCoords, i));
            }

            // Find edges to refine in current mesh
            var currentEdgesToRefine = new HashSet<(int, int)>();

            for (int i = 0; i < refinedMesh.Count<Tri3>(); i++)
            {
                var n = refinedMesh.NodesOf<Tri3, Node>(i);
                int n1 = n[0], n2 = n[1], n3 = n[2];

                double f1 = cachedFieldValues[n1];
                double f2 = cachedFieldValues[n2];
                double f3 = cachedFieldValues[n3];

                // Check region constraint
                bool checkRegion = (regionField != null);
                bool refineEdge12 = false, refineEdge23 = false, refineEdge31 = false;

                if (f1 * f2 <= 0)
                {
                    if (checkRegion)
                    {
                        double r1 = cachedRegionValues[n1];
                        double r2 = cachedRegionValues[n2];
                        refineEdge12 = (r1 <= 0 || r2 <= 0);
                    }
                    else refineEdge12 = true;
                }

                if (f2 * f3 <= 0)
                {
                    if (checkRegion)
                    {
                        double r2 = cachedRegionValues[n2];
                        double r3 = cachedRegionValues[n3];
                        refineEdge23 = (r2 <= 0 || r3 <= 0);
                    }
                    else refineEdge23 = true;
                }

                if (f3 * f1 <= 0)
                {
                    if (checkRegion)
                    {
                        double r3 = cachedRegionValues[n3];
                        double r1 = cachedRegionValues[n1];
                        refineEdge31 = (r3 <= 0 || r1 <= 0);
                    }
                    else refineEdge31 = true;
                }

                if (refineEdge12) currentEdgesToRefine.Add(n1 < n2 ? (n1, n2) : (n2, n1));
                if (refineEdge23) currentEdgesToRefine.Add(n2 < n3 ? (n2, n3) : (n3, n2));
                if (refineEdge31) currentEdgesToRefine.Add(n1 < n3 ? (n1, n3) : (n3, n1));
            }

            if (currentEdgesToRefine.Count > 0)
            {
                Console.WriteLine($"  → Refining {currentEdgesToRefine.Count} edges");
                (refinedMesh, _) = MeshRefinement.Refine(refinedMesh, currentEdgesToRefine.ToList());
                refinedCoords = MeshRefinement.InterpolateCoordinates(refinedMesh, refinedCoords);
            }
        }
        
        // Step 3: Snap new crack nodes to EXACT zero-crossing (surface = 0)
        int snappedCount = 0;
        double maxSnapError = 0.0;
        int skippedTooClose = 0;
        int skippedOutsideRegion = 0;
        int skippedWouldInvert = 0;
        
        for (int i = originalNodeCount; i < refinedMesh.Count<Node>(); i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            
            if (parents.Parent1 != parents.Parent2)
            {
                int p1 = parents.Parent1, p2 = parents.Parent2;
                
                double x1 = refinedCoords[p1, 0], y1 = refinedCoords[p1, 1], z1 = GetZ(refinedCoords, p1);
                double x2 = refinedCoords[p2, 0], y2 = refinedCoords[p2, 1], z2 = GetZ(refinedCoords, p2);
                
                double f1 = signedField(x1, y1, z1);
                double f2 = signedField(x2, y2, z2);
                
                                // Snap to zero-crossing if edge crosses crack

                // Robust edge–surface intersection (curved implicit surfaces may intersect with same-sign endpoints)
                if (TryFindEdgeRootOnSegment(signedField,
                        x1, y1, z1,
                        x2, y2, z2,
                        out double tRoot,
                        out double xRoot, out double yRoot, out double zRoot,
                        out double phiRoot))
                {
                    // Current interpolated position
                    double xCurr = refinedCoords[i, 0];
                    double yCurr = refinedCoords[i, 1];
                    double zCurr = refinedCoords[i, 2];
                    double fCurr = signedField(xCurr, yCurr, zCurr);

                    // Don't snap if already very close to surface
                    double edgeLength = Math.Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
                    double snapTolerance = Math.Max(1e-6, 1e-5 * edgeLength); // 0.001% of edge

                    if (Math.Abs(fCurr) < snapTolerance)
                    {
                        skippedTooClose++;
                        continue;
                    }

                    // CRITICAL: Only snap if the new position is INSIDE the crack region
                    if (regionField != null)
                    {
                        double regionValue = regionField(xRoot, yRoot, zRoot);
                        if (regionValue > 0)
                        {
                            skippedOutsideRegion++;
                            continue;
                        }
                    }

                    double snapDist = Math.Sqrt((xRoot-xCurr)*(xRoot-xCurr) + (yRoot-yCurr)*(yRoot-yCurr) + (zRoot-zCurr)*(zRoot-zCurr));

                    if (snapDist > snapTolerance)
                    {
                        bool wouldInvert = false; // (kept disabled as before)
                        if (!wouldInvert)
                        {
                            refinedCoords[i, 0] = xRoot;
                            refinedCoords[i, 1] = yRoot;
                            refinedCoords[i, 2] = zRoot;

                            double finalError = Math.Abs(signedField(refinedCoords[i, 0], refinedCoords[i, 1], refinedCoords[i, 2]));
                            maxSnapError = Math.Max(maxSnapError, finalError);
                            snappedCount++;
                        }
                        else
                        {
                            skippedWouldInvert++;
                        }
                    }
                    else
                    {
                        skippedTooClose++;
                    }
                }

            }
        }
        
        Console.WriteLine($"  → Snapped {snappedCount} nodes to surface=0 (max error: {maxSnapError:E3})");
        
        // Step 4: RECALCULATE nodal values of surface and region on refined mesh
        Console.WriteLine($"  → Recalculating surface and region values on ALL refined nodes...");
        
        var surfaceValues = new double[refinedMesh.Count<Node>()];
        var regionValues = new double[refinedMesh.Count<Node>()];
        
        for (int i = 0; i < refinedMesh.Count<Node>(); i++)
        {
            surfaceValues[i] = signedField(refinedCoords[i, 0], refinedCoords[i, 1], GetZ(refinedCoords, i));
            
            if (regionField != null)
            {
                regionValues[i] = regionField(refinedCoords[i, 0], refinedCoords[i, 1], GetZ(refinedCoords, i));
            }
            else
            {
                regionValues[i] = -1.0; // All nodes in active region if no region field
            }
        }
        
        // Discover edges before computing average length
        DiscoverEdges(refinedMesh);
        
        // Compute average edge length for tolerance
        double avgEdgeLen = 0;
        int numEdges2 = 0;
        for (int i = 0; i < refinedMesh.Count<Edge>(); i++)
        {
            var n = refinedMesh.NodesOf<Edge, Node>(i);
            double dx = refinedCoords[n[1],0] - refinedCoords[n[0],0];
            double dy = refinedCoords[n[1],1] - refinedCoords[n[0],1];
            double dz = coordDim == 3 ? refinedCoords[n[1],2] - refinedCoords[n[0],2] : 0;
            avgEdgeLen += Math.Sqrt(dx*dx + dy*dy + dz*dz);
            numEdges2++;
        }
        avgEdgeLen /= numEdges2;
        
        // Step 5: Identify crack nodes (nodes where |surface| ≈ 0 AND region < -tol)
        // Use 0.1% of avgEdgeLen to exclude crack tip boundary
        const double zeroTolerance = 1e-5;
        double regionTolerance = 0.001 * avgEdgeLen;  // 0.1% of average edge
        var crackNodes = new HashSet<int>();
        int rejectedByRegion = 0;
        
        for (int i = 0; i < refinedMesh.Count<Node>(); i++)
        {
            if (Math.Abs(surfaceValues[i]) < zeroTolerance)
            {
                if (regionValues[i] < -regionTolerance)
                {
                    crackNodes.Add(i);
                }
                else
                {
                    rejectedByRegion++;
                }
            }
        }
        
        if (rejectedByRegion > 0)
        {
            Console.WriteLine($"  → Rejected {rejectedByRegion} nodes at crack tip (r ≥ {-regionTolerance:E2})");
        }
        
        Console.WriteLine($"  → Found {crackNodes.Count} crack nodes (|surface| < {zeroTolerance:E1}, region < {-regionTolerance:E2})");
        
        if (crackNodes.Count == 0)
        {
            Console.WriteLine($"  ⚠️  No crack nodes found!");
            return (mesh, coordinates);
        }
        
        // Step 6: Identify ORIGINAL mesh boundary (before refinement)
        Console.WriteLine($"  → Identifying original mesh boundary...");
        var originalBoundaryNodes = FindBoundaryNodes(mesh);
        Console.WriteLine($"  → Found {originalBoundaryNodes.Count} boundary nodes in original mesh");
        
        // Map refined nodes back to original boundary
        // A refined node is on original boundary if its parents are both on original boundary
        var refinedBoundaryNodes = new HashSet<int>();
        for (int i = 0; i < refinedMesh.Count<Node>(); i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            if (parents.Parent1 == parents.Parent2)
            {
                // Original node (not created by refinement)
                if (originalBoundaryNodes.Contains(parents.Parent1))
                    refinedBoundaryNodes.Add(i);
            }
            else
            {
                // Refined node (created on edge between Parent1 and Parent2)
                if (originalBoundaryNodes.Contains(parents.Parent1) && 
                    originalBoundaryNodes.Contains(parents.Parent2))
                    refinedBoundaryNodes.Add(i);
            }
        }
        
        Console.WriteLine($"  → Mapped to {refinedBoundaryNodes.Count} boundary nodes in refined mesh");
        
        // Step 7: Identify tip nodes (crack nodes at crack boundary/termination)
        // A tip node is a crack node that has neighbors OUTSIDE the crack region
        var tipNodes = new HashSet<int>();
        
        if (regionField != null)
        {
            // Build node-to-node connectivity
            var nodeNeighbors = BuildNodeNeighborsTri(refinedMesh);
            
            // Check each crack node for neighbors outside active region
            // Interior tips (crack terminates inside mesh): DON'T duplicate
            // Boundary tips (crack reaches mesh edge): DO duplicate
            foreach (int crackNode in crackNodes)
            {
                // Check if any neighbor is outside the crack region (r > 0)
                bool hasFarNeighbor = false;
                foreach (int neighbor in nodeNeighbors[crackNode])
                {
                    if (regionValues[neighbor] > 0)
                    {
                        hasFarNeighbor = true;
                        break;
                    }
                }
                
                // If this is a tip AND not on mesh boundary, don't duplicate it
                if (hasFarNeighbor && !refinedBoundaryNodes.Contains(crackNode))
                {
                    tipNodes.Add(crackNode);
                }
            }
        }
        
        Console.WriteLine($"  → Found {tipNodes.Count} tip nodes (will NOT be duplicated)");
        
        // Step 8: Nodes to duplicate = crack nodes - tip nodes
        var nodesToDuplicate = new HashSet<int>(crackNodes);
        nodesToDuplicate.ExceptWith(tipNodes);
        
        Console.WriteLine($"  → Will duplicate {nodesToDuplicate.Count} nodes");
        
        if (nodesToDuplicate.Count == 0)
        {
            Console.WriteLine($"  → No nodes to duplicate - returning refined mesh");
            return (refinedMesh, refinedCoords);
        }
        
        // Step 9: Duplicate nodes CAREFULLY
        return DuplicateNodesCarefully(
            refinedMesh,
            refinedCoords,
            nodesToDuplicate,
            surfaceValues,
            visualizationOffset,
            enableSmoothing);
    }
    
    /// <summary>
    /// Duplicate crack nodes to create discontinuity.
    /// Element-based approach: for each element, determine which side based on non-crack nodes.
    /// </summary>
    private static (SimplexMesh, double[,]) DuplicateNodesCarefully(
        SimplexMesh mesh,
        double[,] coords,
        HashSet<int> nodesToDuplicate,
        double[] surfaceValues,
        double visualizationOffset = 0.0,
        bool enableSmoothing = false)
    {
        Console.WriteLine($"  → Starting element-based crack duplication...");
        
        int totalNodes = mesh.Count<Node>();
        int coordDim = coords.GetLength(1);
        
        // CREATE ALL NODES (originals + duplicates)
        var crackedMesh = new SimplexMesh();
        int finalNodeCount = totalNodes + nodesToDuplicate.Count;
        var crackedCoords = new double[finalNodeCount, coordDim];
        
        var originalMap = new Dictionary<int, int>();
        var duplicateMap = new Dictionary<int, int>();
        
        // Copy all original nodes
        for (int i = 0; i < totalNodes; i++)
        {
            int newId = crackedMesh.Add<Node>();
            originalMap[i] = newId;
            
            for (int d = 0; d < coordDim; d++)
                crackedCoords[newId, d] = coords[i, d];
            
            var parents = mesh.Get<Node, ParentNodes>(i);
            crackedMesh.Set<Node, ParentNodes>(newId, parents);
        }
        
        // Create duplicates
        foreach (int nodeId in nodesToDuplicate)
        {
            int dupId = crackedMesh.Add<Node>();
            duplicateMap[nodeId] = dupId;
            
            for (int d = 0; d < coordDim; d++)
                crackedCoords[dupId, d] = coords[nodeId, d];
            
            var parents = mesh.Get<Node, ParentNodes>(nodeId);
            crackedMesh.Set<Node, ParentNodes>(dupId, parents);
        }
        
        Console.WriteLine($"  → Created {totalNodes} original + {nodesToDuplicate.Count} duplicate nodes");
        
        // RECREATE ELEMENTS - determine side based on non-crack nodes
        int elementsCreated = 0;
        int positiveSideElements = 0;
        int negativeSideElements = 0;
        int duplicatesUsed = 0;
        
        for (int i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            var newNodes = new int[3];
            
            // Determine element side from first non-crack node
            double sideValue = 0.0;
            bool hasNonCrackNode = false;
            
            for (int j = 0; j < 3; j++)
            {
                if (!nodesToDuplicate.Contains(nodes[j]))
                {
                    sideValue = surfaceValues[nodes[j]];
                    hasNonCrackNode = true;
                    break;
                }
            }
            
            // If element has no non-crack nodes, default to negative side
            bool usePositiveSide = hasNonCrackNode && (sideValue > 0);
            
            if (usePositiveSide)
                positiveSideElements++;
            else
                negativeSideElements++;
            
            // Assign nodes
            for (int j = 0; j < 3; j++)
            {
                if (!nodesToDuplicate.Contains(nodes[j]))
                {
                    // Non-crack node
                    newNodes[j] = originalMap[nodes[j]];
                }
                else
                {
                    // Crack node - use duplicate for positive, original for negative
                    if (usePositiveSide)
                    {
                        newNodes[j] = duplicateMap[nodes[j]];
                        duplicatesUsed++;
                    }
                    else
                    {
                        newNodes[j] = originalMap[nodes[j]];
                    }
                }
            }
            
            crackedMesh.Add<Tri3, Node>(newNodes[0], newNodes[1], newNodes[2]);
            elementsCreated++;
        }
        
        Console.WriteLine($"  → Created {elementsCreated} elements");
        Console.WriteLine($"  → Positive side: {positiveSideElements} elements, Negative side: {negativeSideElements} elements");
        Console.WriteLine($"  → Duplicate nodes used {duplicatesUsed} times in elements");
        Console.WriteLine($"  → Final mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tri3>()} triangles");
        
        // OPTIONAL: VISUALIZATION OFFSET (separate crack surfaces)
        if (visualizationOffset > 0)
        {
            Console.WriteLine($"  → Applying visualization offset: {visualizationOffset}mm");
            
            // Simple approach: offset along surface normal approximation
            foreach (var crackNode in nodesToDuplicate)
            {
                int origId = originalMap[crackNode];
                int dupId = duplicateMap[crackNode];
                
                // Approximate normal from surface gradient (simplified)
                double nx = 0, ny = 0;
                
                // Use surface values to estimate gradient direction
                double surfVal = surfaceValues[crackNode];
                if (Math.Abs(surfVal) < 1e-6)
                {
                    // On crack - use simple approximation
                    nx = 0;
                    ny = 1.0;  // Default to y-direction
                }
                
                // Normalize
                double len = Math.Sqrt(nx * nx + ny * ny + 1e-10);
                nx /= len;
                ny /= len;
                
                // Offset nodes in opposite directions
                crackedCoords[origId, 0] -= (visualizationOffset / 2.0) * nx;
                crackedCoords[origId, 1] -= (visualizationOffset / 2.0) * ny;
                
                crackedCoords[dupId, 0] += (visualizationOffset / 2.0) * nx;
                crackedCoords[dupId, 1] += (visualizationOffset / 2.0) * ny;
            }
        }
        
        // OPTIONAL: SMOOTHING (with crack nodes fixed as boundaries)
        if (enableSmoothing)
        {
            Console.WriteLine($"  → Smoothing mesh (all boundaries fixed)...");
            
            var boundaryNodes = new HashSet<int>();
            
            // Fix ALL boundary nodes in the cracked mesh (includes crack surfaces)
            var crackedBoundary = FindBoundaryNodes(crackedMesh);
            foreach (var node in crackedBoundary)
                boundaryNodes.Add(node);
            
            Console.WriteLine($"  → Fixed {boundaryNodes.Count} boundary nodes (mesh edges + crack surfaces)");
            
            // Fast conservative Laplacian smoothing
            var smoothedCoords = (double[,])crackedCoords.Clone();
            int iterations = 3;  // Reduced from 5
            double relaxation = 0.25;  // More conservative (was 0.5)
            
            // Build neighbor map once (FAST)
            var nodeNeighbors = BuildNodeNeighborsTri(crackedMesh);
            
            for (int iter = 0; iter < iterations; iter++)
            {
                var newCoords = (double[,])smoothedCoords.Clone();
                int rejected = 0;
                
                for (int nodeId = 0; nodeId < crackedMesh.Count<Node>(); nodeId++)
                {
                    if (boundaryNodes.Contains(nodeId)) continue;
                    
                    var neighbors = nodeNeighbors[nodeId];
                    if (neighbors.Count == 0) continue;
                    
                    // Compute centroid
                    double[] centroid = new double[coordDim];
                    foreach (var nbr in neighbors)
                        for (int d = 0; d < coordDim; d++)
                            centroid[d] += smoothedCoords[nbr, d];
                    
                    for (int d = 0; d < coordDim; d++)
                        centroid[d] /= neighbors.Count;
                    
                    // Conservative move
                    double[] newPos = new double[coordDim];
                    for (int d = 0; d < coordDim; d++)
                        newPos[d] = smoothedCoords[nodeId, d] + relaxation * (centroid[d] - smoothedCoords[nodeId, d]);
                    
                    // Quality check: reject if move inverts any adjacent triangle
                    bool valid = true;
                    var incidentTris = crackedMesh.ElementsAt<Tri3, Node>(nodeId);
                    foreach (var elemId in incidentTris)
                    {
                        var n = crackedMesh.NodesOf<Tri3, Node>(elemId);

                        // Compute area with newPos substituted for nodeId (avoids full array clone)
                        double x0 = n[0] == nodeId ? newPos[0] : smoothedCoords[n[0], 0];
                        double y0 = n[0] == nodeId ? newPos[1] : smoothedCoords[n[0], 1];
                        double x1 = n[1] == nodeId ? newPos[0] : smoothedCoords[n[1], 0];
                        double y1 = n[1] == nodeId ? newPos[1] : smoothedCoords[n[1], 1];
                        double x2 = n[2] == nodeId ? newPos[0] : smoothedCoords[n[2], 0];
                        double y2 = n[2] == nodeId ? newPos[1] : smoothedCoords[n[2], 1];

                        double area2 = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);

                        if (area2 <= 0)
                        {
                            valid = false;
                            rejected++;
                            break;
                        }
                    }
                    
                    if (valid)
                    {
                        for (int d = 0; d < coordDim; d++)
                            newCoords[nodeId, d] = newPos[d];
                    }
                }
                
                smoothedCoords = newCoords;
            }
            
            CorrectTriangleOrientations(crackedMesh, smoothedCoords);
            return (crackedMesh, smoothedCoords);
        }
        
        CorrectTriangleOrientations(crackedMesh, crackedCoords);
        return (crackedMesh, crackedCoords);
    }
    
    /// <summary>
    /// Find boundary nodes (nodes that share boundary edges).
    /// </summary>
    public static HashSet<int> FindBoundaryNodes(SimplexMesh mesh)
    {
        var boundaryNodes = new HashSet<int>();
        
        DiscoverEdges(mesh);
        
        // Boundary edges are shared by exactly 1 triangle
        for (int e = 0; e < mesh.Count<Edge>(); e++)
        {
            if (mesh.CountElementsSharingSubEntity<Tri3, Edge, Node>(e) == 1)
            {
                var edgeNodes = mesh.NodesOf<Edge, Node>(e);
                boundaryNodes.Add(edgeNodes[0]);
                boundaryNodes.Add(edgeNodes[1]);
            }
        }
        
        return boundaryNodes;
    }
    
    #region Crack Insertion - 3D (Tetrahedral Meshes)
    
    /// <summary>
/// Robustly finds an intersection (root) of an implicit signed field on a line segment.
/// This is needed for curved fields (e.g. cylinders) where the segment may intersect even if
/// both endpoints have the same sign (two roots are possible in general).
/// </summary>
private static bool TryFindEdgeRootOnSegment(
    SignedFieldFunction phi,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    out double tRoot,
    out double xRoot, out double yRoot, out double zRoot,
    out double phiAtRoot)
{
    // Conservative sampling to detect sign changes / near-zero points.
    // Keep this small for performance; 9 points catches most curved intersections.
    const int S = 9; // includes endpoints
    const double nearZero = 1e-12;

    double bestAbs = double.PositiveInfinity;
    double bestT = 0.5;
    double bestPhi = double.PositiveInfinity;

    double prevT = 0.0;
    double prevPhi = phi(x1, y1, z1);

    bestAbs = Math.Abs(prevPhi);
    bestT = 0.0;
    bestPhi = prevPhi;

    // Track best bracketing interval (closest to mid-segment)
    bool haveBracket = false;
    double aT = 0.0, bT = 1.0;
    double aPhi = prevPhi, bPhi = double.NaN;
    double bestBracketScore = double.PositiveInfinity;

    for (int k = 1; k < S; k++)
    {
        double t = (double)k / (S - 1);
        double x = x1 + t * (x2 - x1);
        double y = y1 + t * (y2 - y1);
        double z = z1 + t * (z2 - z1);
        double v = phi(x, y, z);

        double av = Math.Abs(v);
        if (av < bestAbs)
        {
            bestAbs = av;
            bestT = t;
            bestPhi = v;
        }

        // Bracket if sign changes or either endpoint is near zero
        if ((prevPhi == 0.0) || (v == 0.0) || (prevPhi * v < 0.0) || (Math.Abs(prevPhi) < nearZero) || (Math.Abs(v) < nearZero))
        {
            double mid = 0.5 * (prevT + t);
            double score = Math.Abs(mid - 0.5); // prefer closest to midpoint (stable, symmetric)
            if (score < bestBracketScore)
            {
                bestBracketScore = score;
                haveBracket = true;
                aT = prevT; bT = t;
                aPhi = prevPhi; bPhi = v;
            }
        }

        prevT = t;
        prevPhi = v;
    }

    // If we didn't find a bracket but we got extremely close somewhere, treat it as an intersection.
    if (!haveBracket)
    {
        if (bestAbs > 1e-9) // not close enough -> no reliable intersection
        {
            tRoot = 0; xRoot = yRoot = zRoot = 0; phiAtRoot = 0;
            return false;
        }

        tRoot = bestT;
        xRoot = x1 + tRoot * (x2 - x1);
        yRoot = y1 + tRoot * (y2 - y1);
        zRoot = z1 + tRoot * (z2 - z1);
        phiAtRoot = bestPhi;
        return true;
    }

    // Bisection on the selected bracket (works even if phi is nonlinear, as long as continuous)
    double loT = aT, hiT = bT;
    double loPhi = aPhi, hiPhi = bPhi;

    // If the bracket is degenerate (exact zero at one endpoint), return that
    if (Math.Abs(loPhi) < nearZero)
    {
        tRoot = loT;
        xRoot = x1 + tRoot * (x2 - x1);
        yRoot = y1 + tRoot * (y2 - y1);
        zRoot = z1 + tRoot * (z2 - z1);
        phiAtRoot = loPhi;
        return true;
    }
    if (Math.Abs(hiPhi) < nearZero)
    {
        tRoot = hiT;
        xRoot = x1 + tRoot * (x2 - x1);
        yRoot = y1 + tRoot * (y2 - y1);
        zRoot = z1 + tRoot * (z2 - z1);
        phiAtRoot = hiPhi;
        return true;
    }

    // Ensure opposite signs; if not, still bisect toward minimum |phi|
    bool opposite = loPhi * hiPhi < 0.0;

    double tMid = 0.5 * (loT + hiT);
    double phiMid = 0.0;

    for (int iter = 0; iter < 40; iter++)
    {
        tMid = 0.5 * (loT + hiT);
        double xm = x1 + tMid * (x2 - x1);
        double ym = y1 + tMid * (y2 - y1);
        double zm = z1 + tMid * (z2 - z1);
        phiMid = phi(xm, ym, zm);

        if (Math.Abs(phiMid) < 1e-12) break;

        if (opposite)
        {
            // standard sign bisection
            if (loPhi * phiMid < 0.0)
            {
                hiT = tMid; hiPhi = phiMid;
            }
            else
            {
                loT = tMid; loPhi = phiMid;
            }
        }
        else
        {
            // fallback: shrink interval toward smaller |phi|
            if (Math.Abs(loPhi) < Math.Abs(hiPhi))
            {
                hiT = tMid; hiPhi = phiMid;
            }
            else
            {
                loT = tMid; loPhi = phiMid;
            }
        }
    }

    tRoot = tMid;
    xRoot = x1 + tRoot * (x2 - x1);
    yRoot = y1 + tRoot * (y2 - y1);
    zRoot = z1 + tRoot * (z2 - z1);
    phiAtRoot = phiMid;
    return true;
}

public static (SimplexMesh, double[,]) CreateCrackFromSignedField3D(
        SimplexMesh mesh,
        double[,] coordinates,
        SignedFieldFunction signedField,
        SignedFieldFunction? regionField = null,
        bool enableMeshPerturbation = false,
        bool enableSmoothing = false)
    {
        Console.WriteLine($"  → Creating 3D crack with mesh refinement and node duplication...");
        Console.WriteLine($"  → Initial mesh: {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tetrahedra");
        
        // Find edges that cross the crack surface
        DiscoverEdges(mesh);
        var edgesToRefine = new List<(int, int)>();
        
        int totalCrossingEdges = 0;
        int regionFilteredEdges = 0;
        
        for (int i = 0; i < mesh.Count<Edge>(); i++)
        {
            var nodes = mesh.NodesOf<Edge, Node>(i);
            int n1 = nodes[0], n2 = nodes[1];
            
            // Robust edge–surface intersection (curved implicit surfaces may intersect with same-sign endpoints)
if (TryFindEdgeRootOnSegment(signedField,
        coordinates[n1, 0], coordinates[n1, 1], coordinates[n1, 2],
        coordinates[n2, 0], coordinates[n2, 1], coordinates[n2, 2],
        out double tRoot,
        out double xRoot, out double yRoot, out double zRoot,
        out double phiRoot))
{
    totalCrossingEdges++;

    if (regionField != null)
    {
        // Filter based on the intersection point (NOT only endpoints)
        double r = regionField(xRoot, yRoot, zRoot);
        if (r <= 0)
            edgesToRefine.Add((n1, n2));
        else
            regionFilteredEdges++;
    }
    else
    {
        edgesToRefine.Add((n1, n2));
    }
}
        }
        
        Console.WriteLine($"  → Edges crossing crack surface: {totalCrossingEdges}");
        if (regionField != null)
            Console.WriteLine($"  → Filtered by region: {regionFilteredEdges}");
        Console.WriteLine($"  → Edges to refine: {edgesToRefine.Count}");
        
        if (edgesToRefine.Count == 0)
        {
            Console.WriteLine($"  ⚠️  No edges to refine - returning original mesh");
            return (mesh, coordinates);
        }
        
        // CHECK INITIAL MESH JACOBIANS (before refinement)
        Console.WriteLine($"  → Checking initial mesh quality...");
        int initialNegative = 0;
        for (int i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            double jac = ComputeTetrahedronJacobian(coordinates, nodes[0], nodes[1], nodes[2], nodes[3]);
            if (jac <= 0) initialNegative++;
        }
        if (initialNegative > 0)
            Console.WriteLine($"  ⚠️  Initial mesh has {initialNegative}/{mesh.Count<Tet4>()} elements with non-positive Jacobian!");
        else
            Console.WriteLine($"  ✓ Initial mesh: all elements have positive Jacobians");
        
        // Store original node count before refinement
        int originalNodeCount = mesh.Count<Node>();
        
        // Use MeshRefinement.Refine to split edges AND subdivide tetrahedra
        Console.WriteLine($"  → Calling MeshRefinement.Refine with closure enforcement...");
        var (refinedMesh, _) = MeshRefinement.Refine(mesh, edgesToRefine, 
            enforceClosureForTets: false,
            validateTopology: true,
            inputCoordinates: coordinates);
        
        Console.WriteLine($"  → Refined mesh: {refinedMesh.Count<Node>()} nodes, {refinedMesh.Count<Tet4>()} tetrahedra");
        
        // Use MeshRefinement's own coordinate interpolation
        Console.WriteLine($"  → Using MeshRefinement.InterpolateCoordinates...");
        var refinedCoords = MeshRefinement.InterpolateCoordinates(refinedMesh, coordinates);
        
        // Save intermediate refined mesh for visualization
        Console.WriteLine($"  → Saving refined mesh before snapping...");
        SimplexRemesher.SaveGiD(refinedMesh, refinedCoords, "refined_crack_mesh.post.msh");
        
        // Snap new nodes to exact crack surface (surface = 0) ONLY IF INSIDE REGION
        Console.WriteLine($"  → Snapping crack nodes to exact zero crossing...");
        int snappedCount = 0;
        int skippedTooClose = 0;
        int skippedWouldInvert = 0;
        int skippedOutsideRegion = 0;
        double maxSnapError = 0.0;
        
        for (int i = originalNodeCount; i < refinedMesh.Count<Node>(); i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            
            if (parents.Parent1 != parents.Parent2)
            {
                int p1 = parents.Parent1, p2 = parents.Parent2;
                
                double x1 = refinedCoords[p1, 0], y1 = refinedCoords[p1, 1], z1 = refinedCoords[p1, 2];
                double x2 = refinedCoords[p2, 0], y2 = refinedCoords[p2, 1], z2 = refinedCoords[p2, 2];
                
                double f1 = signedField(x1, y1, z1);
                double f2 = signedField(x2, y2, z2);
                
                // Snap to zero-crossing if edge crosses crack
                if (Math.Abs(f2 - f1) > 1e-10 && f1 * f2 <= 0)
                {
                    // Current interpolated position
                    double xCurr = refinedCoords[i, 0];
                    double yCurr = refinedCoords[i, 1];
                    double zCurr = refinedCoords[i, 2];
                    double fCurr = signedField(xCurr, yCurr, zCurr);
                    
                    // Don't snap if already very close to surface
                    double edgeLength = Math.Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
                    double snapTolerance = Math.Max(1e-6, 1e-5 * edgeLength); // Much tighter: 0.001% of edge
                    
                    if (Math.Abs(fCurr) < snapTolerance)
                    {
                        skippedTooClose++;
                        continue;
                    }
                    
                    // Use bisection for nonlinear fields
                    double tMin = 0.0, tMax = 1.0;
                    double t = 0.5;
                    
                    for (int iter = 0; iter < 30; iter++)  // More iterations for curved surfaces
                    {
                        t = (tMin + tMax) / 2.0;
                        double xm = x1 + t * (x2 - x1);
                        double ym = y1 + t * (y2 - y1);
                        double zm = z1 + t * (z2 - z1);
                        double fm = signedField(xm, ym, zm);
                        
                        if (Math.Abs(fm) < 1e-12) break;  // Tighter convergence
                        
                        if (fm * f1 < 0)
                        {
                            tMax = t;
                            f2 = fm;
                        }
                        else
                        {
                            tMin = t;
                            f1 = fm;
                        }
                    }
                    
                    // Check if snap distance is significant
                    double xNew = x1 + t * (x2 - x1);
                    double yNew = y1 + t * (y2 - y1);
                    double zNew = z1 + t * (z2 - z1);
                    
                    // CRITICAL FIX: Only snap if the new position is INSIDE the crack region
                    if (regionField != null)
                    {
                        double regionValue = regionField(xNew, yNew, zNew);
                        if (regionValue > 0)
                        {
                            // Outside the crack region (e.g., outside the ellipse)
                            // Do NOT snap - leave the node at its interpolated position
                            skippedOutsideRegion++;
                            continue;
                        }
                    }
                    
                    double snapDist = Math.Sqrt((xNew-xCurr)*(xNew-xCurr) + (yNew-yCurr)*(yNew-yCurr) + (zNew-zCurr)*(zNew-zCurr));
                    
                    if (snapDist > snapTolerance)
                    {
                        // Validate: check if snapping would invert any adjacent tet
                        bool wouldInvert = false;
                        
                        // TEMPORARILY DISABLED - diagnose snapping issues
                        /*
                        var incidentTets = refinedMesh.ElementsAt<Tet4, Node>(i);
                        {
                            var tempCoords = (double[,])refinedCoords.Clone();
                            tempCoords[i, 0] = xNew;
                            tempCoords[i, 1] = yNew;
                            tempCoords[i, 2] = zNew;
                            
                            foreach (var tetIdx in incidentTets)
                            {
                                var tetNodes = refinedMesh.NodesOf<Tet4, Node>(tetIdx);
                                double jac = ComputeTetrahedronJacobian(tempCoords, tetNodes[0], tetNodes[1], tetNodes[2], tetNodes[3]);
                                if (jac <= 0)
                                {
                                    wouldInvert = true;
                                    break;
                                }
                            }
                        }
                        */
                        
                        if (!wouldInvert)
                        {
                            refinedCoords[i, 0] = xNew;
                            refinedCoords[i, 1] = yNew;
                            refinedCoords[i, 2] = zNew;
                            
                            double finalError = Math.Abs(signedField(refinedCoords[i, 0], refinedCoords[i, 1], refinedCoords[i, 2]));
                            maxSnapError = Math.Max(maxSnapError, finalError);
                            snappedCount++;
                        }
                        else
                        {
                            skippedWouldInvert++;
                        }
                    }
                    else
                    {
                        skippedTooClose++;
                    }
                }
            }
        }
        
        Console.WriteLine($"  → Snapped {snappedCount} nodes to surface=0 (max error: {maxSnapError:E3})");
        if (skippedTooClose > 0)
            Console.WriteLine($"  → Skipped {skippedTooClose} nodes already close enough to surface");
        if (skippedWouldInvert > 0)
            Console.WriteLine($"  → Skipped {skippedWouldInvert} nodes that would invert tets");
        if (skippedOutsideRegion > 0)
            Console.WriteLine($"  → Skipped {skippedOutsideRegion} nodes outside crack region (tip area)");

        MeshRefinement.CheckJacobians(refinedMesh, refinedCoords, "After snapping");
        MeshRefinement.FixNegativeJacobians(refinedMesh, refinedCoords);
        MeshRefinement.CheckJacobians(refinedMesh, refinedCoords, "Final mesh");

        // ============================================================
        // NODE DUPLICATION - Element-based approach (matching 2D)
        // ============================================================
        Console.WriteLine($"  → Starting element-based crack node duplication...");
        
        // Step 1: Pre-compute signed field values for ALL nodes
        var surfaceValues = new double[refinedMesh.Count<Node>()];
        var regionValues = new double[refinedMesh.Count<Node>()];
        
        for (int i = 0; i < refinedMesh.Count<Node>(); i++)
        {
            surfaceValues[i] = signedField(refinedCoords[i, 0], refinedCoords[i, 1], refinedCoords[i, 2]);
            regionValues[i] = regionField != null 
                ? regionField(refinedCoords[i, 0], refinedCoords[i, 1], refinedCoords[i, 2]) 
                : -1.0; // Default: all nodes in active region
        }
        
        // Compute average edge length for adaptive tolerances
        double avgEdgeLen = 0;
        int numEdgesCount = 0;
        for (int t = 0; t < refinedMesh.Count<Tet4>(); t++)
        {
            var n = refinedMesh.NodesOf<Tet4, Node>(t);
            for (int i = 0; i < 4; i++)
            {
                for (int j = i + 1; j < 4; j++)
                {
                    double dx = refinedCoords[n[j], 0] - refinedCoords[n[i], 0];
                    double dy = refinedCoords[n[j], 1] - refinedCoords[n[i], 1];
                    double dz = refinedCoords[n[j], 2] - refinedCoords[n[i], 2];
                    avgEdgeLen += Math.Sqrt(dx*dx + dy*dy + dz*dz);
                    numEdgesCount++;
                }
            }
        }
        avgEdgeLen /= Math.Max(1, numEdgesCount);
        Console.WriteLine($"  → Average edge length: {avgEdgeLen:F4}");
        
        // Step 2: Identify crack nodes by SIGN CHANGE (not tolerance)
        // A node is a crack node if it was created by edge refinement (Parent1 != Parent2)
        // AND both parents have opposite signs
        
        var crackNodes = new HashSet<int>();
        
        for (int i = originalNodeCount; i < refinedMesh.Count<Node>(); i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            
            if (parents.Parent1 != parents.Parent2)
            {
                int p1 = parents.Parent1, p2 = parents.Parent2;
                
                double f1 = signedField(coordinates[p1, 0], coordinates[p1, 1], coordinates[p1, 2]);
                double f2 = signedField(coordinates[p2, 0], coordinates[p2, 1], coordinates[p2, 2]);
                
                // Sign change = edge crosses surface
                if (f1 * f2 < 0)
                {
                    // Check if inside region
                    if (regionField == null || regionField(refinedCoords[i, 0], refinedCoords[i, 1], refinedCoords[i, 2]) <= 0)
                    {
                        crackNodes.Add(i);
                    }
                }
            }
        }
        
        Console.WriteLine($"  → Found {crackNodes.Count} crack nodes (by sign change)");
        
        if (crackNodes.Count == 0)
        {
            Console.WriteLine($"  ⚠️  No crack nodes found - returning refined mesh");
            return (refinedMesh, refinedCoords);
        }
        
        // Step 3: Identify boundary nodes (to exclude from tip detection)
        var boundaryNodes3D = FindBoundaryNodes3D(refinedMesh);
        Console.WriteLine($"  → Found {boundaryNodes3D.Count} boundary nodes");
        
        // Step 4: Identify tip nodes (crack nodes with neighbors outside crack region)
        var tipNodes = new HashSet<int>();
        
        if (regionField != null)
        {
            // Build node-to-node connectivity from tets
            var nodeNeighbors = BuildNodeNeighborsTet(refinedMesh);
            
            // Threshold for "near region boundary" - should be wider than regionTolerance
            double regionBoundaryThreshold = 0.1 * avgEdgeLen;  // 10% of edge length
            
            // Check each crack node for:
            // 1. Neighbors outside active region (classic tip detection)
            // 2. The node itself being near the region boundary
            int tipByNeighbor = 0;
            int tipByBoundary = 0;
            
            foreach (int crackNode in crackNodes)
            {
                // Check 1: Is this node itself near the region boundary?
                if (regionValues[crackNode] > -regionBoundaryThreshold)
                {
                    // Node is near the crack front (boundary of ellipse, etc.)
                    if (!boundaryNodes3D.Contains(crackNode))
                    {
                        tipNodes.Add(crackNode);
                        tipByBoundary++;
                        continue;
                    }
                }
                
                // Check 2: Does this node have neighbors outside the region?
                bool hasFarNeighbor = false;
                foreach (int neighbor in nodeNeighbors[crackNode])
                {
                    if (regionValues[neighbor] > 0)
                    {
                        hasFarNeighbor = true;
                        break;
                    }
                }
                
                // If tip AND not on mesh boundary, don't duplicate
                if (hasFarNeighbor && !boundaryNodes3D.Contains(crackNode))
                {
                    tipNodes.Add(crackNode);
                    tipByNeighbor++;
                }
            }
            
            Console.WriteLine($"  → Tip detection: {tipByBoundary} by region boundary, {tipByNeighbor} by neighbor check");
        }
        
        Console.WriteLine($"  → Found {tipNodes.Count} tip nodes (will NOT be duplicated)");
        
        // Step 5: Nodes to duplicate = crack nodes - tip nodes
        var nodesToDuplicate = new HashSet<int>(crackNodes);
        nodesToDuplicate.ExceptWith(tipNodes);
        
        Console.WriteLine($"  → Will duplicate {nodesToDuplicate.Count} nodes");
        
        // DIAGNOSTIC: Print bounding box of crack nodes and nodes to duplicate
        if (nodesToDuplicate.Count > 0)
        {
            double minX = double.MaxValue, maxX = double.MinValue;
            double minY = double.MaxValue, maxY = double.MinValue;
            double minZ = double.MaxValue, maxZ = double.MinValue;
            double minRegion = double.MaxValue, maxRegion = double.MinValue;
            
            foreach (int nodeId in nodesToDuplicate)
            {
                minX = Math.Min(minX, refinedCoords[nodeId, 0]);
                maxX = Math.Max(maxX, refinedCoords[nodeId, 0]);
                minY = Math.Min(minY, refinedCoords[nodeId, 1]);
                maxY = Math.Max(maxY, refinedCoords[nodeId, 1]);
                minZ = Math.Min(minZ, refinedCoords[nodeId, 2]);
                maxZ = Math.Max(maxZ, refinedCoords[nodeId, 2]);
                minRegion = Math.Min(minRegion, regionValues[nodeId]);
                maxRegion = Math.Max(maxRegion, regionValues[nodeId]);
            }
            
            Console.WriteLine($"  → Crack bbox: X=[{minX:F2},{maxX:F2}], Y=[{minY:F2},{maxY:F2}], Z=[{minZ:F2},{maxZ:F2}]");
            Console.WriteLine($"  → Region values of crack nodes: [{minRegion:F3},{maxRegion:F3}]");
        }
        
        if (nodesToDuplicate.Count == 0)
        {
            Console.WriteLine($"  → No nodes to duplicate - returning refined mesh");
            return (refinedMesh, refinedCoords);
        }
        
        // Step 6: Create new mesh with duplicated nodes
        int totalNodes = refinedMesh.Count<Node>();
        int finalNodeCount = totalNodes + nodesToDuplicate.Count;
        var crackedMesh = new SimplexMesh();
        var crackedCoords = new double[finalNodeCount, 3];
        
        var originalMap = new Dictionary<int, int>();
        var duplicateMap = new Dictionary<int, int>();
        
        // Copy all original nodes
        for (int i = 0; i < totalNodes; i++)
        {
            int newId = crackedMesh.Add<Node>();
            originalMap[i] = newId;
            
            crackedCoords[newId, 0] = refinedCoords[i, 0];
            crackedCoords[newId, 1] = refinedCoords[i, 1];
            crackedCoords[newId, 2] = refinedCoords[i, 2];
            
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            crackedMesh.Set<Node, ParentNodes>(newId, parents);
        }

        // Create duplicates
        foreach (int nodeId in nodesToDuplicate)
        {
            int dupId = crackedMesh.Add<Node>();
            duplicateMap[nodeId] = dupId;

            crackedCoords[dupId, 0] = refinedCoords[nodeId, 0];
            crackedCoords[dupId, 1] = refinedCoords[nodeId, 1];
            crackedCoords[dupId, 2] = refinedCoords[nodeId, 2];

            var parents = refinedMesh.Get<Node, ParentNodes>(nodeId);
            crackedMesh.Set<Node, ParentNodes>(dupId, parents);
        }
        
        Console.WriteLine($"  → Created {totalNodes} original + {nodesToDuplicate.Count} duplicate nodes");
        
        // Step 7: Recreate elements - determine side based on nodes CLEARLY off the crack surface
        // CRITICAL: Avoid using nodes near the crack plane for side classification
        int elementsCreated = 0;
        int positiveSideElements = 0;
        int negativeSideElements = 0;
        int duplicatesUsed = 0;
        int classifiedByClearNode = 0;
        int classifiedByCentroid = 0;
        
        // Threshold for "clearly off the crack surface" - use 10% of average edge length
        double clearThreshold = 0.1 * avgEdgeLen;
        
        for (int t = 0; t < refinedMesh.Count<Tet4>(); t++)
        {
            var nodes = refinedMesh.NodesOf<Tet4, Node>(t);
            var newNodes = new int[4];
            
            // Count how many nodes are crack nodes to be duplicated
            int crackNodeCount = 0;
            for (int j = 0; j < 4; j++)
            {
                if (nodesToDuplicate.Contains(nodes[j]))
                    crackNodeCount++;
            }
            
            // If element has NO crack nodes, just use originals for all
            if (crackNodeCount == 0)
            {
                for (int j = 0; j < 4; j++)
                    newNodes[j] = originalMap[nodes[j]];
                crackedMesh.Add<Tet4, Node>(newNodes[0], newNodes[1], newNodes[2], newNodes[3]);
                elementsCreated++;
                continue;
            }
            
            // VOTING-BASED ELEMENT ASSIGNMENT
            // Count positive vs negative non-crack nodes
            int positiveVotes = 0;
            int negativeVotes = 0;
            
            for (int j = 0; j < 4; j++)
            {
                int nodeId = nodes[j];
                if (!nodesToDuplicate.Contains(nodeId))
                {
                    // Non-crack node - check which side
                    if (surfaceValues[nodeId] > 0)
                        positiveVotes++;
                    else
                        negativeVotes++;
                }
            }
            
            // Assign element to majority side
            bool usePositiveSide;
            if (positiveVotes != negativeVotes)
            {
                usePositiveSide = positiveVotes > negativeVotes;
                classifiedByClearNode++;
            }
            else
            {
                // Tie or all nodes are crack nodes - use centroid
                double cx = 0, cy = 0, cz = 0;
                for (int j = 0; j < 4; j++)
                {
                    cx += refinedCoords[nodes[j], 0];
                    cy += refinedCoords[nodes[j], 1];
                    cz += refinedCoords[nodes[j], 2];
                }
                cx /= 4; cy /= 4; cz /= 4;
                usePositiveSide = signedField(cx, cy, cz) > 0;
                classifiedByCentroid++;
            }
            
            if (usePositiveSide)
                positiveSideElements++;
            else
                negativeSideElements++;
            
            // Assign nodes
            for (int j = 0; j < 4; j++)
            {
                if (!nodesToDuplicate.Contains(nodes[j]))
                {
                    // Non-crack node - use original
                    newNodes[j] = originalMap[nodes[j]];
                }
                else
                {
                    // Crack node - use duplicate for positive, original for negative
                    if (usePositiveSide)
                    {
                        newNodes[j] = duplicateMap[nodes[j]];
                        duplicatesUsed++;
                    }
                    else
                    {
                        newNodes[j] = originalMap[nodes[j]];
                    }
                }
            }
            
            crackedMesh.Add<Tet4, Node>(newNodes[0], newNodes[1], newNodes[2], newNodes[3]);
            elementsCreated++;
        }
        
        Console.WriteLine($"  → Created {elementsCreated} tetrahedra");
        Console.WriteLine($"  → Positive side: {positiveSideElements}, Negative side: {negativeSideElements}");
        Console.WriteLine($"  → Classification: {classifiedByClearNode} by clear node, {classifiedByCentroid} by centroid");
        Console.WriteLine($"  → Duplicate nodes used {duplicatesUsed} times");
        
        refinedMesh = crackedMesh;
        refinedCoords = crackedCoords;
        
        // Step 8: Validate mesh quality
        Console.WriteLine($"  → Validating cracked mesh...");
        int negativeCount = 0;
        for (int i = 0; i < refinedMesh.Count<Tet4>(); i++)
        {
            var nodes = refinedMesh.NodesOf<Tet4, Node>(i);
            double jac = ComputeTetrahedronJacobian(refinedCoords, nodes[0], nodes[1], nodes[2], nodes[3]);
            if (jac <= 0) negativeCount++;
        }
        
        if (negativeCount > 0)
            Console.WriteLine($"  ⚠️  WARNING: {negativeCount} elements with non-positive Jacobian");
        else
            Console.WriteLine($"  ✓ All elements have positive Jacobians");
        
        // OPTIONAL: 3D SMOOTHING (with boundary and crack nodes fixed)
        if (enableSmoothing)
        {
            Console.WriteLine($"  → Applying 3D smoothing (all boundaries fixed)...");
            
            var boundaryNodes = FindBoundaryNodes3D(refinedMesh);
            Console.WriteLine($"  → Fixed {boundaryNodes.Count} boundary nodes");
            
            // Build neighbor map from tetrahedra
            var nodeNeighbors = BuildNodeNeighborsTet(refinedMesh);
            
            // Conservative Laplacian smoothing
            var smoothedCoords = (double[,])refinedCoords.Clone();
            int iterations = 3;
            double relaxation = 0.25;
            
            for (int iter = 0; iter < iterations; iter++)
            {
                var newCoords = (double[,])smoothedCoords.Clone();
                int smoothedCount = 0;
                int rejectedCount = 0;
                
                for (int nodeId = 0; nodeId < refinedMesh.Count<Node>(); nodeId++)
                {
                    if (boundaryNodes.Contains(nodeId)) continue;
                    
                    var neighbors = nodeNeighbors[nodeId];
                    if (neighbors.Count == 0) continue;
                    
                    // Compute centroid
                    double cx = 0, cy = 0, cz = 0;
                    foreach (var nbr in neighbors)
                    {
                        cx += smoothedCoords[nbr, 0];
                        cy += smoothedCoords[nbr, 1];
                        cz += smoothedCoords[nbr, 2];
                    }
                    cx /= neighbors.Count;
                    cy /= neighbors.Count;
                    cz /= neighbors.Count;
                    
                    // Conservative move
                    double newX = smoothedCoords[nodeId, 0] + relaxation * (cx - smoothedCoords[nodeId, 0]);
                    double newY = smoothedCoords[nodeId, 1] + relaxation * (cy - smoothedCoords[nodeId, 1]);
                    double newZ = smoothedCoords[nodeId, 2] + relaxation * (cz - smoothedCoords[nodeId, 2]);
                    
                    // Quality check: reject if move inverts any adjacent tet
                    bool valid = true;
                    foreach (var tetId in refinedMesh.ElementsAt<Tet4, Node>(nodeId))
                    {
                        var tn = refinedMesh.NodesOf<Tet4, Node>(tetId);
                        
                        // Temporarily update coords for this node
                        double oldX = smoothedCoords[nodeId, 0];
                        double oldY = smoothedCoords[nodeId, 1];
                        double oldZ = smoothedCoords[nodeId, 2];
                        
                        smoothedCoords[nodeId, 0] = newX;
                        smoothedCoords[nodeId, 1] = newY;
                        smoothedCoords[nodeId, 2] = newZ;
                        
                        double jac = ComputeTetrahedronJacobian(smoothedCoords, tn[0], tn[1], tn[2], tn[3]);
                        
                        // Restore
                        smoothedCoords[nodeId, 0] = oldX;
                        smoothedCoords[nodeId, 1] = oldY;
                        smoothedCoords[nodeId, 2] = oldZ;
                        
                        if (jac <= 0)
                        {
                            valid = false;
                            rejectedCount++;
                            break;
                        }
                    }
                    
                    if (valid)
                    {
                        newCoords[nodeId, 0] = newX;
                        newCoords[nodeId, 1] = newY;
                        newCoords[nodeId, 2] = newZ;
                        smoothedCount++;
                    }
                }
                
                smoothedCoords = newCoords;
                Console.WriteLine($"     Iteration {iter + 1}/{iterations}: smoothed {smoothedCount} nodes, rejected {rejectedCount}");
            }
            
            refinedCoords = smoothedCoords;
            Console.WriteLine($"  ✓ 3D smoothing complete");
        }
        
        return (refinedMesh, refinedCoords);
    }
    public static HashSet<int> FindBoundaryNodes3D(SimplexMesh mesh)
    {
        var boundaryNodes = new HashSet<int>();
        
        // Count how many tetrahedra share each face
        var faceCount = new Dictionary<(int, int, int), int>();
        
        for (int i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            
            // Four faces per tetrahedron
            var faces = new[]
            {
                (nodes[0], nodes[1], nodes[2]),
                (nodes[0], nodes[1], nodes[3]),
                (nodes[0], nodes[2], nodes[3]),
                (nodes[1], nodes[2], nodes[3])
            };
            
            foreach (var face in faces)
            {
                // Sort to get canonical form
                var sorted = new[] { face.Item1, face.Item2, face.Item3 };
                Array.Sort(sorted);
                var key = (sorted[0], sorted[1], sorted[2]);
                
                if (!faceCount.ContainsKey(key))
                    faceCount[key] = 0;
                faceCount[key]++;
            }
        }
        
        // Boundary faces have count == 1
        foreach (var (face, count) in faceCount)
        {
            if (count == 1)
            {
                boundaryNodes.Add(face.Item1);
                boundaryNodes.Add(face.Item2);
                boundaryNodes.Add(face.Item3);
            }
        }
        
        return boundaryNodes;
    }
    
    #endregion
    
    #region Node Connectivity Helpers
    
    /// <summary>
    /// Build node-to-node connectivity from triangles.
    /// Returns a dictionary mapping each node to its set of topological neighbors.
    /// </summary>
    private static Dictionary<int, HashSet<int>> BuildNodeNeighborsTri(SimplexMesh mesh)
    {
        var neighbors = new Dictionary<int, HashSet<int>>();
        for (int i = 0; i < mesh.Count<Node>(); i++)
            neighbors[i] = new HashSet<int>();
        
        for (int elemId = 0; elemId < mesh.Count<Tri3>(); elemId++)
        {
            var n = mesh.NodesOf<Tri3, Node>(elemId);
            neighbors[n[0]].Add(n[1]); neighbors[n[0]].Add(n[2]);
            neighbors[n[1]].Add(n[0]); neighbors[n[1]].Add(n[2]);
            neighbors[n[2]].Add(n[0]); neighbors[n[2]].Add(n[1]);
        }
        return neighbors;
    }
    
    /// <summary>
    /// Build node-to-node connectivity from tetrahedra.
    /// Returns a dictionary mapping each node to its set of topological neighbors.
    /// </summary>
    private static Dictionary<int, HashSet<int>> BuildNodeNeighborsTet(SimplexMesh mesh)
    {
        var neighbors = new Dictionary<int, HashSet<int>>();
        for (int i = 0; i < mesh.Count<Node>(); i++)
            neighbors[i] = new HashSet<int>();
        
        for (int elemId = 0; elemId < mesh.Count<Tet4>(); elemId++)
        {
            var n = mesh.NodesOf<Tet4, Node>(elemId);
            for (int i = 0; i < 4; i++)
                for (int j = i + 1; j < 4; j++)
                {
                    neighbors[n[i]].Add(n[j]);
                    neighbors[n[j]].Add(n[i]);
                }
        }
        return neighbors;
    }
    
    #endregion

    #region Orientation Correction

    /// <summary>Correct triangle orientations in-place to ensure CCW (positive signed area)</summary>
    private static void CorrectTriangleOrientations(SimplexMesh mesh, double[,] coords)
    {
        int flipped = 0;
        
        for (int i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var n = mesh.NodesOf<Tri3, Node>(i);
            double area2 = (coords[n[1],0] - coords[n[0],0]) * (coords[n[2],1] - coords[n[0],1]) -
                          (coords[n[1],1] - coords[n[0],1]) * (coords[n[2],0] - coords[n[0],0]);
            
            const double tol = 1e-14;
            if (area2 < -tol)
            {
                // CW → CCW: swap nodes 1 and 2 (in-place)
                mesh.ReplaceElementNodes<Tri3, Node>(i, n[0], n[2], n[1]);
                flipped++;
            }
        }
        
        if (flipped > 0)
            Console.WriteLine($"  → Corrected {flipped} CW→CCW triangles");
    }

    /// <summary>Correct quad orientations in-place to ensure positive Jacobians</summary>
    private static void CorrectQuadOrientations(SimplexMesh mesh, double[,] coords)
    {
        int flipped = 0;
        
        for (int i = 0; i < mesh.Count<Quad4>(); i++)
        {
            var n = mesh.NodesOf<Quad4, Node>(i);
            
            double x1 = coords[n[0],0], y1 = coords[n[0],1];
            double x2 = coords[n[1],0], y2 = coords[n[1],1];
            double x3 = coords[n[2],0], y3 = coords[n[2],1];
            double x4 = coords[n[3],0], y4 = coords[n[3],1];
            
            double dxdxi = 0.25 * (-x1 + x2 + x3 - x4);
            double dydxi = 0.25 * (-y1 + y2 + y3 - y4);
            double dxdeta = 0.25 * (-x1 - x2 + x3 + x4);
            double dydeta = 0.25 * (-y1 - y2 + y3 + y4);
            
            double jac = dxdxi * dydeta - dydxi * dxdeta;
            
            if (jac < 0)
            {
                mesh.ReplaceElementNodes<Quad4, Node>(i, n[0], n[3], n[2], n[1]);
                flipped++;
            }
        }
        
        if (flipped > 0)
            Console.WriteLine($"  → Corrected {flipped} inverted quads");
    }

    /// <summary>Correct tetrahedron orientations in-place to ensure positive volumes</summary>
    private static void CorrectTetOrientations(SimplexMesh mesh, double[,] coords)
    {
        int flipped = 0;
        
        for (int i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var n = mesh.NodesOf<Tet4, Node>(i);
            
            double ax = coords[n[1],0] - coords[n[0],0];
            double ay = coords[n[1],1] - coords[n[0],1];
            double az = coords[n[1],2] - coords[n[0],2];
            double bx = coords[n[2],0] - coords[n[0],0];
            double by = coords[n[2],1] - coords[n[0],1];
            double bz = coords[n[2],2] - coords[n[0],2];
            double cx = coords[n[3],0] - coords[n[0],0];
            double cy = coords[n[3],1] - coords[n[0],1];
            double cz = coords[n[3],2] - coords[n[0],2];
            
            double vol6 = ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx);
            
            if (vol6 < 0)
            {
                // Swap nodes 0 and 1 to fix orientation (in-place)
                mesh.ReplaceElementNodes<Tet4, Node>(i, n[1], n[0], n[2], n[3]);
                flipped++;
            }
        }
        
        if (flipped > 0)
            Console.WriteLine($"  → Corrected {flipped} inverted tetrahedra");
    }
    
    private static (int, int, int) CanonicalFace(int n1, int n2, int n3)
    {
        var sorted = new[] { n1, n2, n3 };
        Array.Sort(sorted);
        return (sorted[0], sorted[1], sorted[2]);
    }
    
    #endregion
}
// MeshOperations.cs - Mesh manipulation operations
// Consolidates: MeshOptimization.cs, MeshRefinement.cs, QuadConversion.cs
// FIXED: RemoveDegenerateTetrahedra now fully implemented
// FIXED: Uses shared Jacobian functions instead of duplicates
// License: GPLv3

using static Numerical.MeshConstants;
using static Numerical.MeshGeometry;

namespace Numerical;

#region Mesh Smoothing (from MeshOptimization.cs)

/// <summary>
///     Mesh optimization algorithms: Laplacian smoothing, CVT smoothing, element cleanup.
/// </summary>
public static class MeshOptimization
{
    #region Laplacian Smoothing

    /// <summary>
    ///     Laplacian mesh smoothing: moves each node to the average of its neighbors.
    /// </summary>
    public static double[,] LaplacianSmoothing(
        SimplexMesh mesh, double[,] coords,
        int iterations = 5, HashSet<int>? fixedNodes = null,
        double relaxation = 1.0, bool keepBoundaryFixed = true)
    {
        var nNodes = coords.GetLength(0);
        var newCoords = (double[,])coords.Clone();

        var boundaryNodes = keepBoundaryFixed ? IdentifyBoundaryNodes(mesh) : new HashSet<int>();
        var allFixedNodes = new HashSet<int>(boundaryNodes);
        if (fixedNodes != null) allFixedNodes.UnionWith(fixedNodes);

        Console.WriteLine($"[LaplacianSmoothing] Starting {iterations} iterations...");
        Console.WriteLine($"[LaplacianSmoothing] Fixed nodes: {allFixedNodes.Count}/{nNodes}");

        for (var iter = 0; iter < iterations; iter++)
        {
            var neighbors = BuildNodeNeighbors(mesh, nNodes);
            var smoothedCount = 0;

            for (var nodeId = 0; nodeId < nNodes; nodeId++)
            {
                if (allFixedNodes.Contains(nodeId)) continue;

                var nodeNeighbors = neighbors[nodeId];
                if (nodeNeighbors.Count == 0) continue;

                double avgX = 0, avgY = 0, avgZ = 0;
                foreach (var neighborId in nodeNeighbors)
                {
                    avgX += newCoords[neighborId, 0];
                    avgY += newCoords[neighborId, 1];
                    avgZ += newCoords[neighborId, 2];
                }

                var count = nodeNeighbors.Count;
                avgX /= count;
                avgY /= count;
                avgZ /= count;

                newCoords[nodeId, 0] += relaxation * (avgX - newCoords[nodeId, 0]);
                newCoords[nodeId, 1] += relaxation * (avgY - newCoords[nodeId, 1]);
                newCoords[nodeId, 2] += relaxation * (avgZ - newCoords[nodeId, 2]);

                smoothedCount++;
            }

            if ((iter + 1) % Math.Max(1, iterations / 5) == 0 || iter == iterations - 1)
                Console.WriteLine($"[LaplacianSmoothing] Iteration {iter + 1}/{iterations}: smoothed {smoothedCount} nodes");
        }

        Console.WriteLine("[LaplacianSmoothing] ✓ Complete");
        return newCoords;
    }

    #endregion

    #region CVT Smoothing

    /// <summary>
    ///     CVT smoothing: moves nodes toward centroids of their dual Voronoi cells.
    /// </summary>
    public static double[,] CVTSmoothing(
        SimplexMesh mesh, double[,] coords,
        int iterations = 5, HashSet<int>? fixedNodes = null, double relaxation = 1.0)
    {
        var nNodes = coords.GetLength(0);
        var newCoords = (double[,])coords.Clone();

        var boundaryNodes = IdentifyBoundaryNodes(mesh);
        var allFixedNodes = new HashSet<int>(boundaryNodes);
        if (fixedNodes != null) allFixedNodes.UnionWith(fixedNodes);

        Console.WriteLine($"[CVTSmoothing] Starting {iterations} iterations...");

        for (var iter = 0; iter < iterations; iter++)
        {
            var nodeTris = BuildNodeToTriangles(mesh, nNodes);
            var smoothedCount = 0;

            for (var nodeId = 0; nodeId < nNodes; nodeId++)
            {
                if (allFixedNodes.Contains(nodeId)) continue;

                var triangles = nodeTris[nodeId];
                if (triangles.Count == 0) continue;

                double totalArea = 0;
                double centroidX = 0, centroidY = 0;

                foreach (var triId in triangles)
                {
                    var nodes = mesh.NodesOf<Tri3, Node>(triId);

                    var cx = (newCoords[nodes[0], 0] + newCoords[nodes[1], 0] + newCoords[nodes[2], 0]) / 3.0;
                    var cy = (newCoords[nodes[0], 1] + newCoords[nodes[1], 1] + newCoords[nodes[2], 1]) / 3.0;

                    var area = ComputeTriangleArea(newCoords, nodes[0], nodes[1], nodes[2]);

                    totalArea += area;
                    centroidX += area * cx;
                    centroidY += area * cy;
                }

                if (totalArea > Epsilon)
                {
                    centroidX /= totalArea;
                    centroidY /= totalArea;

                    newCoords[nodeId, 0] += relaxation * (centroidX - newCoords[nodeId, 0]);
                    newCoords[nodeId, 1] += relaxation * (centroidY - newCoords[nodeId, 1]);

                    smoothedCount++;
                }
            }

            if ((iter + 1) % Math.Max(1, iterations / 5) == 0 || iter == iterations - 1)
                Console.WriteLine($"[CVTSmoothing] Iteration {iter + 1}/{iterations}: smoothed {smoothedCount} nodes");
        }

        Console.WriteLine("[CVTSmoothing] ✓ Complete");
        return newCoords;
    }

    #endregion

    #region Element Cleanup

    /// <summary>
    ///     Remove degenerate triangles (zero or near-zero area).
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) RemoveDegenerateTriangles(
        SimplexMesh mesh, double[,] coords, double tolerance = Epsilon)
    {
        Console.WriteLine($"[RemoveDegenerateTriangles] Checking {mesh.Count<Tri3>()} triangles...");

        var degenerateTris = new HashSet<int>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            if (IsTriangleDegenerate(coords, nodes[0], nodes[1], nodes[2], tolerance))
                degenerateTris.Add(i);
        }

        Console.WriteLine($"[RemoveDegenerateTriangles] Found {degenerateTris.Count} degenerate triangles");

        if (degenerateTris.Count == 0)
        {
            Console.WriteLine("[RemoveDegenerateTriangles] ✓ No cleanup needed");
            return (mesh, coords);
        }

        return RebuildMeshWithoutElements(mesh, coords, degenerateTris, new HashSet<int>());
    }

    /// <summary>
    ///     Remove degenerate tetrahedra (zero or near-zero volume).
    ///     FIXED: Now fully implemented instead of returning input unchanged.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) RemoveDegenerateTetrahedra(
        SimplexMesh mesh, double[,] coords, double tolerance = Epsilon)
    {
        Console.WriteLine($"[RemoveDegenerateTetrahedra] Checking {mesh.Count<Tet4>()} tetrahedra...");

        var degenerateTets = new HashSet<int>();

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            if (IsTetrahedronDegenerate(coords, nodes[0], nodes[1], nodes[2], nodes[3], tolerance))
                degenerateTets.Add(i);
        }

        Console.WriteLine($"[RemoveDegenerateTetrahedra] Found {degenerateTets.Count} degenerate tetrahedra");

        if (degenerateTets.Count == 0)
        {
            Console.WriteLine("[RemoveDegenerateTetrahedra] ✓ No cleanup needed");
            return (mesh, coords);
        }

        return RebuildMeshWithoutElements(mesh, coords, new HashSet<int>(), degenerateTets);
    }

    /// <summary>
    ///     Rebuild mesh excluding specified degenerate elements.
    /// </summary>
    private static (SimplexMesh mesh, double[,] coords) RebuildMeshWithoutElements(
        SimplexMesh mesh, double[,] coords,
        HashSet<int> excludeTriangles, HashSet<int> excludeTetrahedra)
    {
        // Collect nodes that are still used
        var usedNodes = new HashSet<int>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            if (excludeTriangles.Contains(i)) continue;
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            usedNodes.Add(nodes[0]);
            usedNodes.Add(nodes[1]);
            usedNodes.Add(nodes[2]);
        }

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            if (excludeTetrahedra.Contains(i)) continue;
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            usedNodes.Add(nodes[0]);
            usedNodes.Add(nodes[1]);
            usedNodes.Add(nodes[2]);
            usedNodes.Add(nodes[3]);
        }

        for (var i = 0; i < mesh.Count<Bar2>(); i++)
        {
            var nodes = mesh.NodesOf<Bar2, Node>(i);
            usedNodes.Add(nodes[0]);
            usedNodes.Add(nodes[1]);
        }

        for (var i = 0; i < mesh.Count<Point>(); i++)
        {
            var nodes = mesh.NodesOf<Point, Node>(i);
            usedNodes.Add(nodes[0]);
        }

        // Create node mapping
        var nodeMap = new Dictionary<int, int>();
        var newCoordsList = new List<double[]>();

        foreach (var oldId in usedNodes.OrderBy(x => x))
        {
            nodeMap[oldId] = newCoordsList.Count;
            newCoordsList.Add(new[] { coords[oldId, 0], coords[oldId, 1], coords[oldId, 2] });
        }

        var newCoords = new double[newCoordsList.Count, 3];
        for (var i = 0; i < newCoordsList.Count; i++)
        {
            newCoords[i, 0] = newCoordsList[i][0];
            newCoords[i, 1] = newCoordsList[i][1];
            newCoords[i, 2] = newCoordsList[i][2];
        }

        // Build new mesh
        var newMesh = new SimplexMesh();

        newMesh.WithBatch(() =>
        {
            for (var i = 0; i < newCoordsList.Count; i++) newMesh.AddNode(i);

            for (var i = 0; i < mesh.Count<Tri3>(); i++)
            {
                if (excludeTriangles.Contains(i)) continue;
                var nodes = mesh.NodesOf<Tri3, Node>(i);
                var idx = newMesh.AddTriangle(nodeMap[nodes[0]], nodeMap[nodes[1]], nodeMap[nodes[2]]);
                var orig = mesh.Get<Tri3, OriginalElement>(i);
                newMesh.Set<Tri3, OriginalElement>(idx, orig);
            }

            for (var i = 0; i < mesh.Count<Tet4>(); i++)
            {
                if (excludeTetrahedra.Contains(i)) continue;
                var nodes = mesh.NodesOf<Tet4, Node>(i);
                var idx = newMesh.AddTetrahedron(nodeMap[nodes[0]], nodeMap[nodes[1]], nodeMap[nodes[2]], nodeMap[nodes[3]]);
                var orig = mesh.Get<Tet4, OriginalElement>(i);
                newMesh.Set<Tet4, OriginalElement>(idx, orig);
            }

            for (var i = 0; i < mesh.Count<Bar2>(); i++)
            {
                var nodes = mesh.NodesOf<Bar2, Node>(i);
                var idx = newMesh.AddBar(nodeMap[nodes[0]], nodeMap[nodes[1]]);
                var orig = mesh.Get<Bar2, OriginalElement>(i);
                newMesh.Set<Bar2, OriginalElement>(idx, orig);
            }

            for (var i = 0; i < mesh.Count<Point>(); i++)
            {
                var nodes = mesh.NodesOf<Point, Node>(i);
                var idx = newMesh.AddPoint(nodeMap[nodes[0]]);
                var orig = mesh.Get<Point, OriginalElement>(i);
                newMesh.Set<Point, OriginalElement>(idx, orig);
            }
        });

        var removedCount = excludeTriangles.Count + excludeTetrahedra.Count;
        Console.WriteLine($"[RebuildMesh] ✓ Removed {removedCount} degenerate elements");
        Console.WriteLine($"[RebuildMesh] Result: {newMesh.Count<Node>()} nodes, {newMesh.Count<Tri3>()} tris, {newMesh.Count<Tet4>()} tets");

        return (newMesh, newCoords);
    }

    #endregion

    #region Helper Functions

    /// <summary>
    ///     Identify boundary nodes (nodes on mesh boundary).
    /// </summary>
    public static HashSet<int> IdentifyBoundaryNodes(SimplexMesh mesh)
    {
        var boundaryNodes = new HashSet<int>();

        if (mesh.Count<Tri3>() > 0)
        {
            var edgeCount = new Dictionary<(int, int), int>();

            for (var i = 0; i < mesh.Count<Tri3>(); i++)
            {
                var nodes = mesh.NodesOf<Tri3, Node>(i);
                IncrementEdgeCount(edgeCount, nodes[0], nodes[1]);
                IncrementEdgeCount(edgeCount, nodes[1], nodes[2]);
                IncrementEdgeCount(edgeCount, nodes[2], nodes[0]);
            }

            foreach (var (edge, count) in edgeCount)
            {
                if (count == 1)
                {
                    boundaryNodes.Add(edge.Item1);
                    boundaryNodes.Add(edge.Item2);
                }
            }
        }

        if (mesh.Count<Tet4>() > 0)
        {
            var faceCount = new Dictionary<(int, int, int), int>();

            for (var i = 0; i < mesh.Count<Tet4>(); i++)
            {
                var nodes = mesh.NodesOf<Tet4, Node>(i);
                IncrementFaceCount(faceCount, nodes[0], nodes[1], nodes[2]);
                IncrementFaceCount(faceCount, nodes[0], nodes[1], nodes[3]);
                IncrementFaceCount(faceCount, nodes[0], nodes[2], nodes[3]);
                IncrementFaceCount(faceCount, nodes[1], nodes[2], nodes[3]);
            }

            foreach (var (face, count) in faceCount)
            {
                if (count == 1)
                {
                    boundaryNodes.Add(face.Item1);
                    boundaryNodes.Add(face.Item2);
                    boundaryNodes.Add(face.Item3);
                }
            }
        }

        return boundaryNodes;
    }

    private static List<HashSet<int>> BuildNodeNeighbors(SimplexMesh mesh, int nNodes)
    {
        var neighbors = new List<HashSet<int>>(nNodes);
        for (var i = 0; i < nNodes; i++) neighbors.Add(new HashSet<int>());

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            neighbors[nodes[0]].Add(nodes[1]);
            neighbors[nodes[0]].Add(nodes[2]);
            neighbors[nodes[1]].Add(nodes[0]);
            neighbors[nodes[1]].Add(nodes[2]);
            neighbors[nodes[2]].Add(nodes[0]);
            neighbors[nodes[2]].Add(nodes[1]);
        }

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            for (var j = 0; j < 4; j++)
                for (var k = 0; k < 4; k++)
                    if (j != k) neighbors[nodes[j]].Add(nodes[k]);
        }

        return neighbors;
    }

    private static List<HashSet<int>> BuildNodeToTriangles(SimplexMesh mesh, int nNodes)
    {
        var nodeTris = new List<HashSet<int>>(nNodes);
        for (var i = 0; i < nNodes; i++) nodeTris.Add(new HashSet<int>());

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            nodeTris[nodes[0]].Add(i);
            nodeTris[nodes[1]].Add(i);
            nodeTris[nodes[2]].Add(i);
        }

        return nodeTris;
    }

    private static void IncrementEdgeCount(Dictionary<(int, int), int> dict, int n0, int n1)
    {
        var edge = n0 < n1 ? (n0, n1) : (n1, n0);
        dict[edge] = dict.GetValueOrDefault(edge, 0) + 1;
    }

    private static void IncrementFaceCount(Dictionary<(int, int, int), int> dict, int n0, int n1, int n2)
    {
        var sorted = new[] { n0, n1, n2 }.OrderBy(x => x).ToArray();
        var face = (sorted[0], sorted[1], sorted[2]);
        dict[face] = dict.GetValueOrDefault(face, 0) + 1;
    }

    #endregion
}

#endregion

#region Quad Conversion (from QuadConversion.cs)

/// <summary>
///     Convert triangular meshes to quad-dominant meshes.
/// </summary>
public static class QuadConversion
{
    /// <summary>
    ///     Convert triangular mesh to quad-dominant mesh.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) ConvertToQuads(
        SimplexMesh mesh, double[,] coords,
        int passes = 2, double minQualityThreshold = 0.3)
    {
        Console.WriteLine("[QuadConversion] Starting triangle-to-quad conversion...");
        Console.WriteLine($"  Input: {mesh.Count<Tri3>()} triangles");

        var adjacency = BuildTriangleAdjacency(mesh);
        var converted = new HashSet<int>();
        var quads = new List<(int n0, int n1, int n2, int n3)>();

        for (var pass = 0; pass < passes; pass++)
        {
            Console.WriteLine($"  Pass {pass + 1}/{passes}...");
            var quadsThisPass = 0;

            quadsThisPass += PairByQuality(mesh, coords, adjacency, converted, quads, minQualityThreshold);
            quadsThisPass += PairByValence(mesh, coords, adjacency, converted, quads, minQualityThreshold);
            quadsThisPass += PairByGeometry(mesh, coords, adjacency, converted, quads, minQualityThreshold);

            Console.WriteLine($"    → Created {quadsThisPass} quads in pass {pass + 1}");
            if (quadsThisPass == 0) break;
        }

        Console.WriteLine($"  Total quads created: {quads.Count}");

        var newMesh = BuildQuadMesh(mesh, quads, converted);
        return (newMesh, coords);
    }

    #region Triangle Adjacency

    private static Dictionary<int, List<int>> BuildTriangleAdjacency(SimplexMesh mesh)
    {
        var edgeToTriangles = new Dictionary<(int, int), List<int>>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            var edges = new[]
            {
                MakeEdge(nodes[0], nodes[1]),
                MakeEdge(nodes[1], nodes[2]),
                MakeEdge(nodes[2], nodes[0])
            };

            foreach (var edge in edges)
            {
                if (!edgeToTriangles.ContainsKey(edge))
                    edgeToTriangles[edge] = new List<int>();
                edgeToTriangles[edge].Add(i);
            }
        }

        var adjacency = new Dictionary<int, List<int>>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            adjacency[i] = new List<int>();
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            var edges = new[]
            {
                MakeEdge(nodes[0], nodes[1]),
                MakeEdge(nodes[1], nodes[2]),
                MakeEdge(nodes[2], nodes[0])
            };

            foreach (var edge in edges)
                foreach (var neighbor in edgeToTriangles[edge])
                    if (neighbor != i && !adjacency[i].Contains(neighbor))
                        adjacency[i].Add(neighbor);
        }

        return adjacency;
    }

    private static (int, int) MakeEdge(int v0, int v1) => v0 < v1 ? (v0, v1) : (v1, v0);

    #endregion

    #region Pairing Strategies

    private static int PairByQuality(SimplexMesh mesh, double[,] coords,
        Dictionary<int, List<int>> adjacency, HashSet<int> converted,
        List<(int, int, int, int)> quads, double minQuality)
    {
        var count = 0;

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            if (converted.Contains(i)) continue;

            var nodesI = mesh.NodesOf<Tri3, Node>(i);
            var bestQuality = minQuality;
            var bestNeighbor = -1;
            var bestQuad = (0, 0, 0, 0);

            foreach (var j in adjacency[i])
            {
                if (converted.Contains(j)) continue;

                var nodesJ = mesh.NodesOf<Tri3, Node>(j);
                var quadNodes = TryFormQuad(nodesI, nodesJ);
                if (!quadNodes.HasValue) continue;

                var quality = ComputeQuadQuality(coords,
                    quadNodes.Value.n0, quadNodes.Value.n1,
                    quadNodes.Value.n2, quadNodes.Value.n3);

                if (quality > bestQuality)
                {
                    bestQuality = quality;
                    bestNeighbor = j;
                    bestQuad = quadNodes.Value;
                }
            }

            if (bestNeighbor >= 0)
            {
                quads.Add(bestQuad);
                converted.Add(i);
                converted.Add(bestNeighbor);
                count++;
            }
        }

        return count;
    }

    private static int PairByValence(SimplexMesh mesh, double[,] coords,
        Dictionary<int, List<int>> adjacency, HashSet<int> converted,
        List<(int, int, int, int)> quads, double minQuality)
    {
        var valence = ComputeNodeValence(mesh, converted);
        var count = 0;

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            if (converted.Contains(i)) continue;

            var nodesI = mesh.NodesOf<Tri3, Node>(i);
            var bestNeighbor = -1;
            var bestScore = double.MinValue;
            var bestQuad = (0, 0, 0, 0);

            foreach (var j in adjacency[i])
            {
                if (converted.Contains(j)) continue;

                var nodesJ = mesh.NodesOf<Tri3, Node>(j);
                var quadNodes = TryFormQuad(nodesI, nodesJ);
                if (!quadNodes.HasValue) continue;

                var quality = ComputeQuadQuality(coords,
                    quadNodes.Value.n0, quadNodes.Value.n1,
                    quadNodes.Value.n2, quadNodes.Value.n3);

                if (quality < minQuality) continue;

                double valenceScore = 0;
                var nodes = new[] { quadNodes.Value.n0, quadNodes.Value.n1, quadNodes.Value.n2, quadNodes.Value.n3 };
                foreach (var node in nodes)
                {
                    var v = valence.GetValueOrDefault(node, 0);
                    if (v == 3 || v == 5) valenceScore += 1.0;
                }

                var score = quality + 0.2 * valenceScore;
                if (score > bestScore)
                {
                    bestScore = score;
                    bestNeighbor = j;
                    bestQuad = quadNodes.Value;
                }
            }

            if (bestNeighbor >= 0)
            {
                quads.Add(bestQuad);
                converted.Add(i);
                converted.Add(bestNeighbor);
                count++;
            }
        }

        return count;
    }

    private static int PairByGeometry(SimplexMesh mesh, double[,] coords,
        Dictionary<int, List<int>> adjacency, HashSet<int> converted,
        List<(int, int, int, int)> quads, double minQuality)
    {
        var count = 0;

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            if (converted.Contains(i)) continue;

            var nodesI = mesh.NodesOf<Tri3, Node>(i);

            foreach (var j in adjacency[i])
            {
                if (converted.Contains(j)) continue;

                var nodesJ = mesh.NodesOf<Tri3, Node>(j);
                var quadNodes = TryFormQuad(nodesI, nodesJ);
                if (!quadNodes.HasValue) continue;

                var quality = ComputeQuadQuality(coords,
                    quadNodes.Value.n0, quadNodes.Value.n1,
                    quadNodes.Value.n2, quadNodes.Value.n3);

                if (quality >= minQuality)
                {
                    quads.Add(quadNodes.Value);
                    converted.Add(i);
                    converted.Add(j);
                    count++;
                    break;
                }
            }
        }

        return count;
    }

    #endregion

    #region Quad Formation and Quality

    private static (int n0, int n1, int n2, int n3)? TryFormQuad(IReadOnlyList<int> tri1, IReadOnlyList<int> tri2)
    {
        var shared = FindSharedEdge(tri1, tri2);
        if (!shared.HasValue) return null;

        var (s0, s1) = shared.Value;

        int unique1 = -1, unique2 = -1;

        foreach (var v in tri1)
            if (v != s0 && v != s1) { unique1 = v; break; }

        foreach (var v in tri2)
            if (v != s0 && v != s1) { unique2 = v; break; }

        if (unique1 < 0 || unique2 < 0) return null;

        return (s0, unique1, s1, unique2);
    }

    private static (int, int)? FindSharedEdge(IReadOnlyList<int> tri1, IReadOnlyList<int> tri2)
    {
        var edges1 = new[] { MakeEdge(tri1[0], tri1[1]), MakeEdge(tri1[1], tri1[2]), MakeEdge(tri1[2], tri1[0]) };
        var edges2 = new[] { MakeEdge(tri2[0], tri2[1]), MakeEdge(tri2[1], tri2[2]), MakeEdge(tri2[2], tri2[0]) };

        foreach (var e1 in edges1)
            foreach (var e2 in edges2)
                if (e1 == e2) return e1;

        return null;
    }

    private static double ComputeQuadQuality(double[,] coords, int n0, int n1, int n2, int n3)
    {
        if (!IsQuadConvex(coords, n0, n1, n2, n3)) return 0.0;

        // Use shared EdgeLength2D from MeshGeometry
        var e0 = EdgeLength2D(coords, n0, n1);
        var e1 = EdgeLength2D(coords, n1, n2);
        var e2 = EdgeLength2D(coords, n2, n3);
        var e3 = EdgeLength2D(coords, n3, n0);

        var minEdge = Math.Min(Math.Min(e0, e1), Math.Min(e2, e3));
        var maxEdge = Math.Max(Math.Max(e0, e1), Math.Max(e2, e3));

        if (maxEdge < Epsilon) return 0.0;

        var aspectRatio = minEdge / maxEdge;

        var angles = new double[4];
        angles[0] = ComputeAngle(coords, n3, n0, n1);
        angles[1] = ComputeAngle(coords, n0, n1, n2);
        angles[2] = ComputeAngle(coords, n1, n2, n3);
        angles[3] = ComputeAngle(coords, n2, n3, n0);

        double angleQuality = 0;
        foreach (var angle in angles)
        {
            var dev = Math.Abs(angle - 90.0);
            angleQuality += Math.Max(0, 1.0 - dev / 90.0);
        }
        angleQuality /= 4.0;

        return 0.5 * aspectRatio + 0.5 * angleQuality;
    }

    private static double ComputeAngle(double[,] coords, int a, int center, int b)
    {
        var ax = coords[a, 0] - coords[center, 0];
        var ay = coords[a, 1] - coords[center, 1];
        var bx = coords[b, 0] - coords[center, 0];
        var by = coords[b, 1] - coords[center, 1];

        var dot = ax * bx + ay * by;
        var lenA = Math.Sqrt(ax * ax + ay * ay);
        var lenB = Math.Sqrt(bx * bx + by * by);

        if (lenA < Epsilon || lenB < Epsilon) return 90.0;

        var cosAngle = Math.Clamp(dot / (lenA * lenB), -1.0, 1.0);
        return Math.Acos(cosAngle) * 180.0 / Math.PI;
    }

    private static Dictionary<int, int> ComputeNodeValence(SimplexMesh mesh, HashSet<int> excludeTriangles)
    {
        var valence = new Dictionary<int, int>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            if (excludeTriangles.Contains(i)) continue;
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            foreach (var node in nodes)
                valence[node] = valence.GetValueOrDefault(node, 0) + 1;
        }

        return valence;
    }

    #endregion

    #region Mesh Building

    private static SimplexMesh BuildQuadMesh(SimplexMesh oldMesh, List<(int n0, int n1, int n2, int n3)> quads, HashSet<int> convertedTriangles)
    {
        var newMesh = new SimplexMesh();

        newMesh.WithBatch(() =>
        {
            for (var i = 0; i < oldMesh.Count<Node>(); i++) newMesh.AddNode(i);

            foreach (var (n0, n1, n2, n3) in quads)
                newMesh.AddQuad(n0, n1, n2, n3);

            for (var i = 0; i < oldMesh.Count<Tri3>(); i++)
            {
                if (convertedTriangles.Contains(i)) continue;
                var nodes = oldMesh.NodesOf<Tri3, Node>(i);
                newMesh.AddTriangle(nodes[0], nodes[1], nodes[2]);
            }
        });

        Console.WriteLine($"[QuadConversion] Built mesh: {newMesh.Count<Quad4>()} quads, {newMesh.Count<Tri3>()} triangles");

        return newMesh;
    }

    #endregion
}

#endregion

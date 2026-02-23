// MeshOperations.cs - Mesh manipulation operations
// Consolidates: MeshOptimization.cs, QuadConversion.cs, CrackOperations.cs
// FIXED: RemoveDegenerateTetrahedra now fully implemented
// FIXED: Uses shared Jacobian functions instead of duplicates
// License: GPLv3

using static Numerical.MeshConstants;
using static Numerical.MeshGeometry;
using static Numerical.MeshRefinement;

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

#region Crack Operations (from CrackOperations.cs)

/// <summary>
///     Crack node duplication for forming cracks after edge refinement.
///     Creates physical cracks in meshes via node duplication and level-set based side assignment.
/// </summary>
public static class CrackDuplication
{
    #region Mesh Renumbering

    /// <summary>
    ///     Renumber mesh to ensure continuous node and element numbering (no gaps).
    /// </summary>
    private static (SimplexMesh mesh, double[,] coords) RenumberMesh(SimplexMesh mesh, double[,] coords)
    {
        Console.WriteLine("[RenumberMesh] Ensuring continuous numbering...");

        var nNodes = mesh.Count<Node>();

        // Collect all used nodes from elements
        var usedNodes = new HashSet<int>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            usedNodes.Add(nodes[0]);
            usedNodes.Add(nodes[1]);
            usedNodes.Add(nodes[2]);
        }

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
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

        // Check if renumbering is needed
        var sortedNodes = usedNodes.OrderBy(x => x).ToList();
        var needsRenumbering = sortedNodes.Count != nNodes;

        if (!needsRenumbering)
        {
            for (int i = 0; i < sortedNodes.Count; i++)
            {
                if (sortedNodes[i] != i)
                {
                    needsRenumbering = true;
                    break;
                }
            }
        }

        if (!needsRenumbering)
        {
            Console.WriteLine($"[RenumberMesh] ✓ Nodes: {nNodes} (0 to {nNodes - 1}, continuous)");
            return (mesh, coords);
        }

        for (int i = mesh.Count<Edge>() - 1; i >= 0; i--)
            mesh.Remove<Edge>(i);

        for (int i = 0; i < nNodes; i++)
            if (!usedNodes.Contains(i))
                mesh.Remove<Node>(i);

        mesh.Compress();

        var newCoords = new double[sortedNodes.Count, 3];
        for (int i = 0; i < sortedNodes.Count; i++)
        {
            var oldId = sortedNodes[i];
            newCoords[i, 0] = coords[oldId, 0];
            newCoords[i, 1] = coords[oldId, 1];
            newCoords[i, 2] = coords[oldId, 2];
        }

        Console.WriteLine($"[RenumberMesh] ✓ Renumbered: {nNodes} → {sortedNodes.Count} nodes (removed {nNodes - sortedNodes.Count} unused)");
        Console.WriteLine($"[RenumberMesh] ✓ Elements: {mesh.Count<Tri3>()} tris + {mesh.Count<Tet4>()} tets");

        return (mesh, newCoords);
    }

    #endregion

    #region Public API

    /// <summary>
    ///     Creates a crack by refining edges, smoothing, and duplicating nodes.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) CreateCrack(
        SimplexMesh mesh,
        double[,] coords,
        List<(int, int)> crackEdges,
        Func<double, double, double> levelSetFunction,
        int smoothingIterations = 5)
    {
        var originalNodeCount = mesh.Count<Node>();

        Console.WriteLine($"[CreateCrack] Refining {crackEdges.Count} crack edges...");
        var (refinedMesh, _) = Refine(mesh, crackEdges);
        var refinedCoords = InterpolateCoordinates(refinedMesh, coords);

        Console.WriteLine($"[CreateCrack] After refinement: {refinedMesh.Count<Node>()} nodes");

        var newNodes = IdentifyMidpointNodes(refinedMesh, originalNodeCount);
        Console.WriteLine($"[CreateCrack] Found {newNodes.Count} new midpoint nodes");

        if (smoothingIterations > 0)
        {
            Console.WriteLine($"[CreateCrack] Smoothing mesh ({smoothingIterations} iterations)...");
            refinedCoords = MeshOptimization.LaplacianSmoothing(
                refinedMesh, refinedCoords, smoothingIterations, newNodes);
            Console.WriteLine("[CreateCrack] ✓ Smoothing complete");
        }

        var tipNodes = IdentifyTipNodes(refinedMesh, newNodes, crackEdges);
        Console.WriteLine($"[CreateCrack] Identified {tipNodes.Count} tip nodes (will NOT be duplicated)");

        var nodesToDuplicate = new HashSet<int>(newNodes);
        nodesToDuplicate.ExceptWith(tipNodes);

        Console.WriteLine($"[CreateCrack] Duplicating {nodesToDuplicate.Count} interior crack nodes...");

        var (crackedMesh, crackedCoords) = DuplicateNodesAndAssignSides(
            refinedMesh, refinedCoords, nodesToDuplicate, newNodes, levelSetFunction);

        return RenumberMesh(crackedMesh, crackedCoords);
    }

    /// <summary>
    ///     Creates a crack from an ALREADY REFINED mesh (for exact geometry).
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) CreateCrackFromRefinedMesh(
        SimplexMesh refinedMesh,
        double[,] refinedCoords,
        int originalNodeCount,
        List<(int, int)> crackEdges,
        Func<double, double, double> levelSetFunction,
        int smoothingIterations = 5)
    {
        Console.WriteLine("[CreateCrackFromRefinedMesh] Using pre-refined mesh with exact geometry");
        Console.WriteLine($"[CreateCrackFromRefinedMesh] Refined mesh: {refinedMesh.Count<Node>()} nodes");

        var newNodes = IdentifyMidpointNodes(refinedMesh, originalNodeCount);
        Console.WriteLine($"[CreateCrackFromRefinedMesh] Found {newNodes.Count} crack nodes");

        DiagnoseZeroAreaTriangles(refinedMesh, refinedCoords);

        var nodeMapping = new Dictionary<int, int>();
        (refinedMesh, refinedCoords, nodeMapping) = MergeDuplicateNodesWithMapping(
            refinedMesh, refinedCoords, newNodes);

        if (nodeMapping.Count > 0)
        {
            var updatedNewNodes = new HashSet<int>();
            foreach (var node in newNodes)
            {
                var canonical = node;
                while (nodeMapping.ContainsKey(canonical))
                    canonical = nodeMapping[canonical];
                updatedNewNodes.Add(canonical);
            }
            newNodes = updatedNewNodes;
            Console.WriteLine($"[CreateCrackFromRefinedMesh] After merging: {newNodes.Count} unique crack nodes");
        }

        (refinedMesh, refinedCoords) = MeshOptimization.RemoveDegenerateTriangles(refinedMesh, refinedCoords);

        if (smoothingIterations > 0)
        {
            Console.WriteLine($"[CreateCrackFromRefinedMesh] Smoothing mesh ({smoothingIterations} iterations)...");
            refinedCoords = MeshOptimization.LaplacianSmoothing(
                refinedMesh, refinedCoords, smoothingIterations, newNodes);
            Console.WriteLine("[CreateCrackFromRefinedMesh] ✓ Smoothing complete");
        }

        var tipNodes = IdentifyTipNodes(refinedMesh, newNodes, crackEdges);
        Console.WriteLine($"[CreateCrackFromRefinedMesh] Identified {tipNodes.Count} TRUE tip nodes");

        var nodesToDuplicate = new HashSet<int>(newNodes);
        nodesToDuplicate.ExceptWith(tipNodes);

        Console.WriteLine($"[CreateCrackFromRefinedMesh] Duplicating {nodesToDuplicate.Count} crack nodes...");

        var (crackedMesh, crackedCoords) = DuplicateNodesAndAssignSides(
            refinedMesh, refinedCoords, nodesToDuplicate, newNodes, levelSetFunction);

        return RenumberMesh(crackedMesh, crackedCoords);
    }

    #endregion

    #region Crack-Specific Core Logic

    private static HashSet<int> IdentifyMidpointNodes(SimplexMesh mesh, int originalNodeCount)
    {
        var midpointNodes = new HashSet<int>();
        for (var i = originalNodeCount; i < mesh.Count<Node>(); i++)
        {
            var parents = mesh.Get<Node, ParentNodes>(i);
            if (parents.Parent1 != parents.Parent2) midpointNodes.Add(i);
        }
        return midpointNodes;
    }

    private static HashSet<int> IdentifyTipNodes(
        SimplexMesh mesh, HashSet<int> newNodes, List<(int, int)> originalCrackEdges)
    {
        var tipNodes = new HashSet<int>();
        var refinedCrackEdges = new HashSet<(int, int)>();

        foreach (var (n0, n1) in originalCrackEdges)
        {
            var nodesOnEdge = new List<int>();
            foreach (var nodeId in newNodes)
            {
                var parents = mesh.Get<Node, ParentNodes>(nodeId);
                if ((parents.Parent1 == n0 && parents.Parent2 == n1) ||
                    (parents.Parent1 == n1 && parents.Parent2 == n0))
                    nodesOnEdge.Add(nodeId);
            }

            if (nodesOnEdge.Count > 0)
            {
                nodesOnEdge.Sort();
                for (var i = 0; i < nodesOnEdge.Count - 1; i++)
                {
                    var edge = (Math.Min(nodesOnEdge[i], nodesOnEdge[i + 1]),
                        Math.Max(nodesOnEdge[i], nodesOnEdge[i + 1]));
                    refinedCrackEdges.Add(edge);
                }
            }
        }

        var edgeCount = new Dictionary<int, int>();
        foreach (var (n0, n1) in refinedCrackEdges)
        {
            if (newNodes.Contains(n0))
                edgeCount[n0] = edgeCount.GetValueOrDefault(n0, 0) + 1;
            if (newNodes.Contains(n1))
                edgeCount[n1] = edgeCount.GetValueOrDefault(n1, 0) + 1;
        }

        foreach (var (nodeId, count) in edgeCount)
            if (count == 1)
                tipNodes.Add(nodeId);

        Console.WriteLine($"[IdentifyTipNodes] Analyzed {refinedCrackEdges.Count} refined crack edges");
        Console.WriteLine($"[IdentifyTipNodes] Found {tipNodes.Count} tip nodes: {string.Join(", ", tipNodes.Take(10))}");

        return tipNodes;
    }

    private static (SimplexMesh mesh, double[,] coords) DuplicateNodesAndAssignSides(
        SimplexMesh mesh, double[,] coords,
        HashSet<int> nodesToDuplicate, HashSet<int> allCrackNodes,
        Func<double, double, double> levelSetFunction)
    {
        var nNodes = mesh.Count<Node>();
        var nNewNodes = nNodes + nodesToDuplicate.Count;
        var newCoords = new double[nNewNodes, 3];

        for (var i = 0; i < nNodes; i++)
        {
            newCoords[i, 0] = coords[i, 0];
            newCoords[i, 1] = coords[i, 1];
            newCoords[i, 2] = coords[i, 2];
        }

        var nodeDuplicates = new Dictionary<int, int>();
        var nextId = nNodes;

        foreach (var nodeId in nodesToDuplicate)
        {
            var dupId = nextId++;
            nodeDuplicates[nodeId] = dupId;
            newCoords[dupId, 0] = coords[nodeId, 0];
            newCoords[dupId, 1] = coords[nodeId, 1];
            newCoords[dupId, 2] = coords[nodeId, 2];
        }

        Console.WriteLine($"[DuplicateNodes] Created {nodeDuplicates.Count} duplicate nodes");
        Console.WriteLine($"[DuplicateNodes] Total nodes: {nNodes} → {nNewNodes}");

        var newMesh = new SimplexMesh();
        newMesh.WithBatch(() =>
        {
            for (var i = 0; i < nNewNodes; i++) newMesh.AddNode(i);
            CopyElementsWithSideAssignment(mesh, newMesh, allCrackNodes, nodeDuplicates, newCoords, levelSetFunction);
        });

        return (newMesh, newCoords);
    }

    private static void CopyElementsWithSideAssignment(
        SimplexMesh mesh, SimplexMesh newMesh,
        HashSet<int> crackNodes, Dictionary<int, int> nodeDuplicates,
        double[,] coords, Func<double, double, double> levelSet)
    {
        var useOriginalCount = 0;
        var useDuplicateCount = 0;
        var levelSetStats = new List<(int elemId, bool useDup, double[] nodeValues)>();

        bool ShouldUseDuplicates(IReadOnlyList<int> nodes, out double[] levelSetValues)
        {
            levelSetValues = new double[nodes.Count];
            for (var i = 0; i < nodes.Count; i++)
            {
                var nodeId = nodes[i];
                if (!crackNodes.Contains(nodeId))
                {
                    var phi = levelSet(coords[nodeId, 0], coords[nodeId, 1]);
                    levelSetValues[i] = phi;
                    if (phi <= 0) return false;
                }
                else
                {
                    levelSetValues[i] = double.NaN;
                }
            }
            return true;
        }

        int RemapNode(int nodeId, bool useDup)
        {
            if (useDup && nodeDuplicates.ContainsKey(nodeId))
                return nodeDuplicates[nodeId];
            return nodeId;
        }

        for (var i = 0; i < mesh.Count<Point>(); i++)
        {
            var nodes = mesh.NodesOf<Point, Node>(i);
            var idx = newMesh.AddPoint(nodes[0]);
            newMesh.Set<Point, OriginalElement>(idx, mesh.Get<Point, OriginalElement>(i));
        }

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            var useDup = ShouldUseDuplicates(nodes, out var phiValues);
            if (useDup) useDuplicateCount++; else useOriginalCount++;
            if (levelSetStats.Count < 10) levelSetStats.Add((i, useDup, phiValues));

            var idx = newMesh.AddTriangle(RemapNode(nodes[0], useDup), RemapNode(nodes[1], useDup), RemapNode(nodes[2], useDup));
            newMesh.Set<Tri3, OriginalElement>(idx, mesh.Get<Tri3, OriginalElement>(i));
        }

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            var useDup = ShouldUseDuplicates(nodes, out _);
            var idx = newMesh.AddTetrahedron(RemapNode(nodes[0], useDup), RemapNode(nodes[1], useDup),
                RemapNode(nodes[2], useDup), RemapNode(nodes[3], useDup));
            newMesh.Set<Tet4, OriginalElement>(idx, mesh.Get<Tet4, OriginalElement>(i));
        }

        Console.WriteLine($"[CopyElements] Created {newMesh.Count<Tri3>()} triangles, {newMesh.Count<Tet4>()} tetrahedra");
        Console.WriteLine($"[CopyElements] Side assignment: {useOriginalCount} original, {useDuplicateCount} duplicate");

        if (levelSetStats.Count > 0)
        {
            Console.WriteLine("[CopyElements] Sample element classifications:");
            foreach (var (elemId, useDup, phiValues) in levelSetStats)
            {
                var side = useDup ? "DUPLICATE" : "ORIGINAL";
                var values = string.Join(", ", phiValues.Select(v => double.IsNaN(v) ? "crack" : $"{v:F6}"));
                Console.WriteLine($"  Tri {elemId}: {side} (φ values: [{values}])");
            }
        }

        var totalElements = useOriginalCount + useDuplicateCount;
        if (totalElements > 0)
        {
            var ratio = (double)useDuplicateCount / totalElements;
            if (ratio < 0.1 || ratio > 0.9)
            {
                Console.WriteLine("[CopyElements] ⚠ WARNING: Imbalanced side assignment!");
                Console.WriteLine($"[CopyElements]   Duplicate side: {useDuplicateCount}/{totalElements} ({ratio * 100:F1}%)");
            }
        }
    }

    #endregion

    #region Node Merging and Cleanup

    private static (SimplexMesh mesh, double[,] coords, Dictionary<int, int> mapping)
        MergeDuplicateNodesWithMapping(
            SimplexMesh mesh, double[,] coords, HashSet<int> crackNodes, double tolerance = 1e-12)
    {
        Console.WriteLine($"[MergeDuplicateNodes] Checking for duplicate nodes (tolerance: {tolerance:E2})...");

        var nNodes = mesh.Count<Node>();
        var nodeMapping = new Dictionary<int, int>();
        var spatialBuckets = new Dictionary<(int, int, int), List<int>>();
        var bucketSize = Math.Max(tolerance * 100, 1e-6);

        for (var i = 0; i < nNodes; i++)
        {
            var key = ((int)Math.Floor(coords[i, 0] / bucketSize),
                       (int)Math.Floor(coords[i, 1] / bucketSize),
                       (int)Math.Floor(coords[i, 2] / bucketSize));
            if (!spatialBuckets.ContainsKey(key))
                spatialBuckets[key] = new List<int>();
            spatialBuckets[key].Add(i);
        }

        var merged = new HashSet<int>();
        var mergeCount = 0;

        foreach (var bucket in spatialBuckets.Values)
        {
            if (bucket.Count < 2) continue;
            for (var i = 0; i < bucket.Count; i++)
            {
                var id1 = bucket[i];
                if (merged.Contains(id1)) continue;
                for (var j = i + 1; j < bucket.Count; j++)
                {
                    var id2 = bucket[j];
                    if (merged.Contains(id2)) continue;

                    var dx = coords[id1, 0] - coords[id2, 0];
                    var dy = coords[id1, 1] - coords[id2, 1];
                    var dz = coords[id1, 2] - coords[id2, 2];
                    var dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);

                    if (dist < tolerance)
                    {
                        if (!(crackNodes.Contains(id1) && crackNodes.Contains(id2)))
                        {
                            var canonicalId = Math.Min(id1, id2);
                            var mergedId = Math.Max(id1, id2);
                            nodeMapping[mergedId] = canonicalId;
                            merged.Add(mergedId);
                            mergeCount++;
                            Console.WriteLine($"[MergeDuplicateNodes]   Merging {mergedId} → {canonicalId} (dist: {dist:E3})");
                        }
                    }
                }
            }
        }

        Console.WriteLine($"[MergeDuplicateNodes] Found {mergeCount} duplicate nodes to merge");
        if (mergeCount == 0) return (mesh, coords, nodeMapping);

        var GetCanonicalNode = (int nodeId) =>
        {
            while (nodeMapping.ContainsKey(nodeId)) nodeId = nodeMapping[nodeId];
            return nodeId;
        };

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            var c0 = GetCanonicalNode(nodes[0]);
            var c1 = GetCanonicalNode(nodes[1]);
            var c2 = GetCanonicalNode(nodes[2]);
            if (c0 != nodes[0] || c1 != nodes[1] || c2 != nodes[2])
                mesh.ReplaceElementNodes<Tri3, Node>(i, c0, c1, c2);
        }

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            var c0 = GetCanonicalNode(nodes[0]);
            var c1 = GetCanonicalNode(nodes[1]);
            var c2 = GetCanonicalNode(nodes[2]);
            var c3 = GetCanonicalNode(nodes[3]);
            if (c0 != nodes[0] || c1 != nodes[1] || c2 != nodes[2] || c3 != nodes[3])
                mesh.ReplaceElementNodes<Tet4, Node>(i, c0, c1, c2, c3);
        }

        var usedNodes = new HashSet<int>();
        for (var i = 0; i < mesh.Count<Tri3>(); i++)
            foreach (var n in mesh.NodesOf<Tri3, Node>(i)) usedNodes.Add(n);
        for (var i = 0; i < mesh.Count<Tet4>(); i++)
            foreach (var n in mesh.NodesOf<Tet4, Node>(i)) usedNodes.Add(n);

        for (int i = mesh.Count<Edge>() - 1; i >= 0; i--) mesh.Remove<Edge>(i);
        for (int i = 0; i < nNodes; i++)
            if (!usedNodes.Contains(i)) mesh.Remove<Node>(i);

        var sortedUsed = usedNodes.OrderBy(x => x).ToList();
        mesh.Compress();

        var newCoords = new double[sortedUsed.Count, 3];
        for (int i = 0; i < sortedUsed.Count; i++)
        {
            var oldId = sortedUsed[i];
            newCoords[i, 0] = coords[oldId, 0];
            newCoords[i, 1] = coords[oldId, 1];
            newCoords[i, 2] = coords[oldId, 2];
        }

        Console.WriteLine($"[MergeDuplicateNodes] ✓ Merged {mergeCount} nodes: {nNodes} → {mesh.Count<Node>()} nodes");

        var compactMapping = new Dictionary<int, int>();
        for (int i = 0; i < sortedUsed.Count; i++) compactMapping[sortedUsed[i]] = i;

        var finalMapping = new Dictionary<int, int>();
        foreach (var (oldId, _) in nodeMapping)
        {
            var finalCanonical = GetCanonicalNode(oldId);
            if (compactMapping.ContainsKey(finalCanonical))
                finalMapping[oldId] = compactMapping[finalCanonical];
        }

        return (mesh, newCoords, finalMapping);
    }

    private static void DiagnoseZeroAreaTriangles(SimplexMesh mesh, double[,] coords, double tolerance = 1e-10)
    {
        var zeroAreaCount = 0;
        var badTris = new List<int>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            if (IsTriangleDegenerate(coords, nodes[0], nodes[1], nodes[2], tolerance))
            {
                zeroAreaCount++;
                if (badTris.Count < 5) badTris.Add(i);
            }
        }

        if (zeroAreaCount > 0)
        {
            Console.WriteLine($"[Diagnostic] ⚠ Found {zeroAreaCount} zero-area triangles!");
            Console.WriteLine($"[Diagnostic] First few: {string.Join(", ", badTris)}");
        }
        else
        {
            Console.WriteLine("[Diagnostic] ✓ No zero-area triangles found");
        }
    }

    #endregion
}

#endregion

#region Topology-Based Node Duplication

/// <summary>
///     Topology-based node duplication: splits a mesh along a set of interface nodes
///     using only graph connectivity (no level-set or geometry required).
/// </summary>
/// <remarks>
///     <para>
///         <b>Algorithm:</b> Given a set of nodes to duplicate ("interface nodes"),
///         builds an element adjacency graph where two elements are connected only
///         if they share a non-interface node. Connected components of this graph
///         identify the independent "sides" of the interface. Component 0 keeps the
///         original interface nodes; each subsequent component receives fresh copies.
///     </para>
///     <para>
///         <b>Advantages over level-set approach:</b>
///         <list type="bullet">
///             <item>Works with arbitrary interfaces (cracks, contacts, material boundaries)</item>
///             <item>Correctly handles branching cracks (3+ components)</item>
///             <item>No signed-field evaluation or snapping required</item>
///             <item>Pure topological operation — geometry enters only via coordinate copying</item>
///         </list>
///     </para>
/// </remarks>
public static class NodeDuplication
{
    /// <summary>
    ///     Result of a topology-based node duplication operation.
    /// </summary>
    /// <param name="Mesh">The new mesh with duplicated connectivity.</param>
    /// <param name="Coords">Coordinate array for all nodes (originals + duplicates).</param>
    /// <param name="ComponentCount">Number of connected components found (≥ 1).</param>
    /// <param name="DuplicateMap">
    ///     Maps (original interface node, component index) → new node index.
    ///     Component 0 maps to the original node indices.
    /// </param>
    public readonly record struct DuplicationResult(
        SimplexMesh Mesh,
        double[,] Coords,
        int ComponentCount,
        Dictionary<(int OriginalNode, int Component), int> DuplicateMap);

    /// <summary>
    ///     Duplicates a set of interface nodes using topology-based side assignment.
    /// </summary>
    /// <param name="mesh">Input mesh.</param>
    /// <param name="coords">Node coordinate array (nNodes × dim).</param>
    /// <param name="interfaceNodes">Set of node indices to duplicate.</param>
    /// <returns>A <see cref="DuplicationResult"/> with the new mesh, coordinates, and mapping.</returns>
    /// <exception cref="ArgumentNullException">Thrown when any argument is null.</exception>
    /// <exception cref="ArgumentException">Thrown when interfaceNodes is empty.</exception>
    public static DuplicationResult DuplicateNodes(
        SimplexMesh mesh,
        double[,] coords,
        IReadOnlySet<int> interfaceNodes)
    {
        ArgumentNullException.ThrowIfNull(mesh);
        ArgumentNullException.ThrowIfNull(coords);
        ArgumentNullException.ThrowIfNull(interfaceNodes);
        if (interfaceNodes.Count == 0)
            throw new ArgumentException("Interface node set must not be empty.", nameof(interfaceNodes));

        int nNodes = mesh.Count<Node>();
        int coordDim = coords.GetLength(1);
        int nTri = mesh.Count<Tri3>();
        int nTet = mesh.Count<Tet4>();

        // ── Step 1: Find connected components along the interface ──

        var (componentOf, componentCount) = FindComponentsAlongInterface(mesh, interfaceNodes);

        // ── Step 2: Create duplicate nodes for components > 0 ──

        var duplicateMap = new Dictionary<(int, int), int>();

        // Component 0 keeps original node indices for interface nodes
        foreach (int n in interfaceNodes)
            duplicateMap[(n, 0)] = n;

        // Count extra nodes: (componentCount - 1) copies per interface node
        int extraNodes = interfaceNodes.Count * (componentCount - 1);
        int finalNodeCount = nNodes + extraNodes;
        var newCoords = new double[finalNodeCount, coordDim];

        // Copy original coordinates
        for (int i = 0; i < nNodes; i++)
            for (int d = 0; d < coordDim; d++)
                newCoords[i, d] = coords[i, d];

        // Create duplicates for components 1..componentCount-1
        int nextId = nNodes;
        for (int comp = 1; comp < componentCount; comp++)
        {
            foreach (int n in interfaceNodes)
            {
                int dupId = nextId++;
                duplicateMap[(n, comp)] = dupId;

                for (int d = 0; d < coordDim; d++)
                    newCoords[dupId, d] = coords[n, d];
            }
        }

        // ── Step 3: Build new mesh with remapped connectivity ──

        var newMesh = new SimplexMesh();
        newMesh.WithBatch(() =>
        {
            for (int i = 0; i < finalNodeCount; i++)
            {
                int idx = newMesh.Add<Node>();
                if (i < nNodes)
                {
                    var parents = mesh.Get<Node, ParentNodes>(i);
                    newMesh.Set<Node, ParentNodes>(idx, parents);
                }
                else
                {
                    newMesh.Set<Node, ParentNodes>(idx, new ParentNodes(i, i));
                }
            }

            for (int e = 0; e < nTri; e++)
            {
                var nodes = mesh.NodesOf<Tri3, Node>(e);
                int comp = componentOf[e];

                int n0 = RemapNode(nodes[0], comp, interfaceNodes, duplicateMap);
                int n1 = RemapNode(nodes[1], comp, interfaceNodes, duplicateMap);
                int n2 = RemapNode(nodes[2], comp, interfaceNodes, duplicateMap);

                int idx = newMesh.AddTriangle(n0, n1, n2);
                newMesh.Set<Tri3, OriginalElement>(idx, mesh.Get<Tri3, OriginalElement>(e));
            }

            for (int e = 0; e < nTet; e++)
            {
                var nodes = mesh.NodesOf<Tet4, Node>(e);
                int comp = componentOf[nTri + e];

                int n0 = RemapNode(nodes[0], comp, interfaceNodes, duplicateMap);
                int n1 = RemapNode(nodes[1], comp, interfaceNodes, duplicateMap);
                int n2 = RemapNode(nodes[2], comp, interfaceNodes, duplicateMap);
                int n3 = RemapNode(nodes[3], comp, interfaceNodes, duplicateMap);

                int idx = newMesh.AddTetrahedron(n0, n1, n2, n3);
                newMesh.Set<Tet4, OriginalElement>(idx, mesh.Get<Tet4, OriginalElement>(e));
            }
        });

        return new DuplicationResult(newMesh, newCoords, componentCount, duplicateMap);
    }

    /// <summary>
    ///     Remaps a node index for a given component.
    ///     Non-interface nodes keep their original index; interface nodes use the duplicate map.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    private static int RemapNode(
        int nodeId, int component,
        IReadOnlySet<int> interfaceNodes,
        Dictionary<(int, int), int> duplicateMap)
    {
        if (interfaceNodes.Contains(nodeId))
            return duplicateMap[(nodeId, component)];
        return nodeId;
    }

    /// <summary>
    ///     Finds connected components of elements when cutting along interface nodes.
    ///     Useful for previewing the effect of node duplication without modifying the mesh.
    /// </summary>
    /// <param name="mesh">Input mesh.</param>
    /// <param name="interfaceNodes">Nodes that define the cut interface.</param>
    /// <returns>
    ///     componentOf: array where result[i] is the component index for element i
    ///     (0..nTri-1 = triangles, nTri..nTri+nTet-1 = tetrahedra).
    ///     componentCount: total number of components found.
    /// </returns>
    public static (int[] componentOf, int componentCount) FindComponentsAlongInterface(
        SimplexMesh mesh,
        IReadOnlySet<int> interfaceNodes)
    {
        ArgumentNullException.ThrowIfNull(mesh);
        ArgumentNullException.ThrowIfNull(interfaceNodes);

        int nTri = mesh.Count<Tri3>();
        int nTet = mesh.Count<Tet4>();
        int totalElements = nTri + nTet;

        if (totalElements == 0)
            return (Array.Empty<int>(), 0);

        // Build node→element map for non-interface nodes only
        var nodeToElements = new Dictionary<int, List<int>>();

        for (int e = 0; e < nTri; e++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(e);
            for (int j = 0; j < nodes.Count; j++)
            {
                int n = nodes[j];
                if (!interfaceNodes.Contains(n))
                {
                    if (!nodeToElements.TryGetValue(n, out var list))
                    {
                        list = new List<int>(4);
                        nodeToElements[n] = list;
                    }
                    list.Add(e);
                }
            }
        }

        for (int e = 0; e < nTet; e++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(e);
            for (int j = 0; j < nodes.Count; j++)
            {
                int n = nodes[j];
                if (!interfaceNodes.Contains(n))
                {
                    if (!nodeToElements.TryGetValue(n, out var list))
                    {
                        list = new List<int>(4);
                        nodeToElements[n] = list;
                    }
                    list.Add(nTri + e);
                }
            }
        }

        // Build element adjacency (shared non-interface node)
        var adj = new List<HashSet<int>>(totalElements);
        for (int i = 0; i < totalElements; i++)
            adj.Add(new HashSet<int>());

        foreach (var (_, elements) in nodeToElements)
        {
            for (int i = 0; i < elements.Count; i++)
                for (int j = i + 1; j < elements.Count; j++)
                {
                    adj[elements[i]].Add(elements[j]);
                    adj[elements[j]].Add(elements[i]);
                }
        }

        // BFS connected components
        var componentOf = new int[totalElements];
        Array.Fill(componentOf, -1);
        int componentCount = 0;

        for (int seed = 0; seed < totalElements; seed++)
        {
            if (componentOf[seed] >= 0) continue;

            int comp = componentCount++;
            var queue = new Queue<int>();
            queue.Enqueue(seed);
            componentOf[seed] = comp;

            while (queue.Count > 0)
            {
                int cur = queue.Dequeue();
                foreach (int nb in adj[cur])
                {
                    if (componentOf[nb] < 0)
                    {
                        componentOf[nb] = comp;
                        queue.Enqueue(nb);
                    }
                }
            }
        }

        return (componentOf, componentCount);
    }
}

#endregion

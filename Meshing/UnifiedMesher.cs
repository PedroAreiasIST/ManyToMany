// UnifiedMesher.cs - High-level meshing API
// REFACTORED: Thin wrapper using shared utilities
// Original: ~3,800 lines → Refactored: ~300 lines (92% reduction)
// License: GPLv3

using static Numerical.MeshConstants;
using static Numerical.MeshGeometry;
using static Numerical.GeometryUtilities;
using static Numerical.CurveUtilities;
using static Numerical.MeshGeneration;
using static Numerical.MeshOptimization;
using static Numerical.QuadConversion;

namespace Numerical;

/// <summary>
///     Unified meshing engine combining:
///     - Bowyer-Watson Delaunay triangulation (via MeshGeneration)
///     - CVT smoothing (via MeshOptimization)
///     - Multi-strategy quad conversion (via QuadConversion)
///     - MATLAB-compatible quality metrics (via MeshGeometry)
/// </summary>
public static class UnifiedMesher
{
    #region Public API

    /// <summary>
    ///     Generate Delaunay triangulation for a 2D polygon.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) Triangulate(
        double[,] boundaryCoords,
        double[,]? interiorCoords = null,
        bool refine = true,
        double maxArea = 0.0,
        double sizeGradation = 1.3,
        bool convertToQuads = true,
        bool enableSmoothing = true)
    {
        return TriangulateWithHoles(boundaryCoords, null, interiorCoords, refine, maxArea, sizeGradation,
            convertToQuads, enableSmoothing);
    }

    /// <summary>
    ///     Generate Delaunay triangulation for a 2D polygon with holes/voids.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) TriangulateWithHoles(
        double[,] outerBoundary,
        List<double[,]>? holes = null,
        double[,]? interiorCoords = null,
        bool refine = true,
        double maxArea = 0.0,
        double sizeGradation = 1.3,
        bool convertToQuads = true,
        bool enableSmoothing = true)
    {
        Console.WriteLine("\n=== TRIANGULATION START ===");
        Console.WriteLine($"Input: {outerBoundary.GetLength(0)} boundary points");
        Console.WriteLine($"Refine: {refine}, SizeGradation: {sizeGradation}");

        // Ensure CCW orientation using shared utility
        if (!IsBoundaryCCW(outerBoundary))
        {
            Console.WriteLine("  ⚠️  Outer boundary is CW, reversing to CCW");
            outerBoundary = ReverseBoundary(outerBoundary);
        }

        // Compute target edge length
        var avgEdgeLength = ComputeAverageEdgeLength(outerBoundary);
        double targetEdgeLength;

        if (sizeGradation > 0)
        {
            targetEdgeLength = avgEdgeLength * sizeGradation;
            Console.WriteLine($"  Adaptive sizing: target edge = {targetEdgeLength:F4}");
        }
        else if (maxArea > 0)
        {
            targetEdgeLength = Math.Sqrt(maxArea * 2.0); // Approximate from area
            Console.WriteLine($"  Fixed sizing: target edge = {targetEdgeLength:F4}");
        }
        else
        {
            targetEdgeLength = avgEdgeLength * 1.5;
            Console.WriteLine($"  Auto sizing: target edge = {targetEdgeLength:F4}");
        }

        // Use shared Delaunay triangulation
        var (points, triangles) = DelaunayTriangulate(outerBoundary, holes, interiorCoords, targetEdgeLength);

        // Convert to SimplexMesh
        var (mesh, coords) = DelaunayToSimplexMesh(points, triangles);

        Console.WriteLine($"  Initial mesh: {mesh.Count<Tri3>()} triangles");

        // Apply CVT smoothing if enabled
        if (enableSmoothing && mesh.Count<Tri3>() > 0)
        {
            Console.WriteLine("  Applying CVT smoothing...");
            coords = CVTSmoothing(mesh, coords, iterations: 5);
        }

        // Convert to quads if requested
        if (convertToQuads && mesh.Count<Tri3>() > 0)
        {
            Console.WriteLine("  Converting to quad-dominant mesh...");
            (mesh, coords) = ConvertToQuads(mesh, coords, passes: 2);
        }

        // Validate and report final mesh
        var stats = ComputeQualityStatistics(mesh, coords);
        Console.WriteLine("\n=== MESH COMPLETE ===");
        Console.WriteLine(stats);

        return (mesh, coords);
    }

    #endregion

    #region Convenience Methods

    /// <summary>
    ///     Create a simple rectangular mesh.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) CreateRectangle(
        double xMin, double yMin, double xMax, double yMax,
        double targetEdgeLength = 0.0,
        bool convertToQuads = true)
    {
        // Create boundary
        var boundary = new double[4, 3]
        {
            { xMin, yMin, 0 },
            { xMax, yMin, 0 },
            { xMax, yMax, 0 },
            { xMin, yMax, 0 }
        };

        // If targetEdgeLength is provided, convert it to maxArea so TriangulateWithHoles
        // picks fixed sizing and reproduces the requested edge size
        if (targetEdgeLength > 0)
        {
            var maxArea = 0.5 * targetEdgeLength * targetEdgeLength;
            return Triangulate(boundary, null, true, maxArea, 0, convertToQuads);
        }

        return Triangulate(boundary, null, true, 0, 1.3, convertToQuads);
    }

    /// <summary>
    ///     Create a mesh from a polygon defined by points.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) CreateFromPolygon(
        List<(double x, double y)> polygon,
        bool convertToQuads = true,
        double sizeGradation = 1.3)
    {
        var boundary = ConvertToArray(polygon, 3);
        return Triangulate(boundary, null, true, 0, sizeGradation, convertToQuads);
    }

    /// <summary>
    ///     Create a mesh with a hole.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) CreateWithHole(
        List<(double x, double y)> outerPolygon,
        List<(double x, double y)> holePolygon,
        bool convertToQuads = true,
        double sizeGradation = 1.3)
    {
        var outer = ConvertToArray(outerPolygon, 3);
        var hole = ConvertToArray(holePolygon, 3);
        var holes = new List<double[,]> { hole };

        return TriangulateWithHoles(outer, holes, null, true, 0, sizeGradation, convertToQuads);
    }

    #endregion

    #region Quality Metrics

    /// <summary>
    ///     Compute and print mesh quality metrics.
    /// </summary>
    public static void PrintQualityReport(SimplexMesh mesh, double[,] coords, string title = "Mesh Quality Report")
    {
        Console.WriteLine($"\n=== {title} ===");

        var stats = ComputeQualityStatistics(mesh, coords);
        Console.WriteLine(stats);

        // Check for inverted elements
        var inverted = ValidateMeshOrientation(mesh, coords, out int invertedTris, out int invertedTets);
        if (inverted > 0)
        {
            Console.WriteLine($"⚠️  WARNING: {invertedTris} inverted triangles, {invertedTets} inverted tetrahedra");
        }
        else
        {
            Console.WriteLine("✓ All elements have correct orientation");
        }

        // Check for degenerate elements
        var (degTris, degTets) = FindDegenerateElements(mesh, coords);
        if (degTris.Length > 0 || degTets.Length > 0)
        {
            Console.WriteLine($"⚠️  WARNING: {degTris.Length} degenerate triangles, {degTets.Length} degenerate tetrahedra");
        }
        else
        {
            Console.WriteLine("✓ No degenerate elements");
        }
    }

    /// <summary>
    ///     Get MATLAB-compatible quality metrics.
    /// </summary>
    public static (double minAngle, double maxAspectRatio, double avgAspectRatio) GetQualityMetrics(
        SimplexMesh mesh, double[,] coords)
    {
        var stats = ComputeQualityStatistics(mesh, coords);

        if (stats.TriangleCount > 0)
        {
            return (stats.MinTriangleAngleDegrees, stats.MaxTriangleAspectRatio, stats.AvgTriangleAspectRatio);
        }
        else if (stats.TetrahedronCount > 0)
        {
            return (0, stats.MaxTetrahedronAspectRatio, stats.AvgTetrahedronAspectRatio);
        }

        return (0, 0, 0);
    }

    #endregion

    #region Mesh Export

    /// <summary>
    ///     Export mesh to simple format (for debugging/visualization).
    /// </summary>
    public static void ExportToConsole(SimplexMesh mesh, double[,] coords)
    {
        Console.WriteLine("\n--- Mesh Export ---");
        Console.WriteLine($"Nodes: {mesh.Count<Node>()}");

        for (int i = 0; i < Math.Min(mesh.Count<Node>(), 10); i++)
        {
            Console.WriteLine($"  Node {i}: ({coords[i, 0]:F4}, {coords[i, 1]:F4}, {coords[i, 2]:F4})");
        }

        if (mesh.Count<Node>() > 10)
            Console.WriteLine($"  ... and {mesh.Count<Node>() - 10} more nodes");

        Console.WriteLine($"Triangles: {mesh.Count<Tri3>()}");
        for (int i = 0; i < Math.Min(mesh.Count<Tri3>(), 5); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            Console.WriteLine($"  Tri {i}: ({nodes[0]}, {nodes[1]}, {nodes[2]})");
        }

        Console.WriteLine($"Quads: {mesh.Count<Quad4>()}");
        for (int i = 0; i < Math.Min(mesh.Count<Quad4>(), 5); i++)
        {
            var nodes = mesh.NodesOf<Quad4, Node>(i);
            Console.WriteLine($"  Quad {i}: ({nodes[0]}, {nodes[1]}, {nodes[2]}, {nodes[3]})");
        }

        Console.WriteLine($"Tetrahedra: {mesh.Count<Tet4>()}");
    }

    #endregion
}

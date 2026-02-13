// MeshGeneration.cs - Mesh generation algorithms and high-level meshing API
// Consolidates: MeshGeneration.cs + DelaunayTriangulation.cs + UnifiedMesher.cs
// Provides structured meshes, Delaunay triangulation, and unified meshing API
// License: GPLv3

using static Numerical.MeshConstants;
using static Numerical.MeshGeometry;
using static Numerical.GeometryUtilities;
using static Numerical.CurveUtilities;
using static Numerical.MeshOptimization;
using static Numerical.QuadConversion;

namespace Numerical;

/// <summary>
///     Mesh generation algorithms: structured meshes and Delaunay triangulation.
/// </summary>
public static class MeshGeneration
{
    #region 2D Rectangular Meshes

    /// <summary>
    ///     Creates a rectangular triangular mesh with 2 triangles per quad cell.
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateRectangularMesh(
        int nx, int ny, double xMin, double xMax, double yMin, double yMax)
    {
        var mesh = new SimplexMesh();
        var coords = new double[(nx + 1) * (ny + 1), 3];

        double dx = (xMax - xMin) / nx;
        double dy = (yMax - yMin) / ny;

        mesh.WithBatch(() =>
        {
            for (int j = 0; j <= ny; j++)
            {
                for (int i = 0; i <= nx; i++)
                {
                    int idx = j * (nx + 1) + i;
                    coords[idx, 0] = xMin + i * dx;
                    coords[idx, 1] = yMin + j * dy;
                    coords[idx, 2] = 0.0;
                    mesh.AddNode(idx);
                }
            }

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    int n0 = j * (nx + 1) + i;
                    int n1 = n0 + 1;
                    int n2 = n0 + (nx + 1);
                    int n3 = n2 + 1;

                    mesh.AddTriangle(n0, n1, n2);
                    mesh.AddTriangle(n1, n3, n2);
                }
            }
        });

        return (mesh, coords);
    }

    /// <summary>
    ///     Creates a rectangular quadrilateral mesh.
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateRectangularQuadMesh(
        int nx, int ny, double xMin, double xMax, double yMin, double yMax)
    {
        var mesh = new SimplexMesh();
        var coords = new double[(nx + 1) * (ny + 1), 3];

        double dx = (xMax - xMin) / nx;
        double dy = (yMax - yMin) / ny;

        mesh.WithBatch(() =>
        {
            for (int j = 0; j <= ny; j++)
            {
                for (int i = 0; i <= nx; i++)
                {
                    int idx = j * (nx + 1) + i;
                    coords[idx, 0] = xMin + i * dx;
                    coords[idx, 1] = yMin + j * dy;
                    coords[idx, 2] = 0.0;
                    mesh.AddNode(idx);
                }
            }

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    int n0 = j * (nx + 1) + i;
                    int n1 = n0 + 1;
                    int n2 = n0 + (nx + 1);
                    int n3 = n2 + 1;

                    mesh.AddQuad(n0, n1, n3, n2);
                }
            }
        });

        return (mesh, coords);
    }

    /// <summary>
    ///     Create a unit square mesh [0,1] × [0,1].
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateUnitSquareMesh(int n)
        => CreateRectangularMesh(n, n, 0.0, 1.0, 0.0, 1.0);

    #endregion

    #region 3D Box Meshes

    /// <summary>
    ///     Creates a 3D box mesh with tetrahedral elements using body-diagonal Kuhn subdivision.
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateBoxMesh(
        int nx, int ny, int nz,
        double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
    {
        var mesh = new SimplexMesh();
        int nodeCount = (nx + 1) * (ny + 1) * (nz + 1);
        var coords = new double[nodeCount, 3];

        double dx = (xMax - xMin) / nx;
        double dy = (yMax - yMin) / ny;
        double dz = (zMax - zMin) / nz;

        int NodeIndex(int i, int j, int k) => k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;

        mesh.WithBatch(() =>
        {
            for (int k = 0; k <= nz; k++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    for (int i = 0; i <= nx; i++)
                    {
                        int idx = NodeIndex(i, j, k);
                        coords[idx, 0] = xMin + i * dx;
                        coords[idx, 1] = yMin + j * dy;
                        coords[idx, 2] = zMin + k * dz;
                        mesh.AddNode(idx);
                    }
                }
            }

            for (int k = 0; k < nz; k++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        int n0 = NodeIndex(i, j, k);
                        int n1 = NodeIndex(i + 1, j, k);
                        int n2 = NodeIndex(i, j + 1, k);
                        int n3 = NodeIndex(i + 1, j + 1, k);
                        int n4 = NodeIndex(i, j, k + 1);
                        int n5 = NodeIndex(i + 1, j, k + 1);
                        int n6 = NodeIndex(i, j + 1, k + 1);
                        int n7 = NodeIndex(i + 1, j + 1, k + 1);

                        // Body-diagonal subdivision (6 tets per cube)
                        mesh.AddTetrahedron(n0, n1, n3, n7);
                        mesh.AddTetrahedron(n0, n3, n2, n7);
                        mesh.AddTetrahedron(n0, n2, n6, n7);
                        mesh.AddTetrahedron(n0, n6, n4, n7);
                        mesh.AddTetrahedron(n0, n4, n5, n7);
                        mesh.AddTetrahedron(n0, n5, n1, n7);
                    }
                }
            }
        });

        FixInvertedTetrahedra(mesh, coords);

        Console.WriteLine($"[CreateBoxMesh] Created {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tetrahedra");

        return (mesh, coords);
    }

    /// <summary>
    ///     Create a unit cube mesh [0,1]³.
    /// </summary>
    public static (SimplexMesh Mesh, double[,] Coordinates) CreateUnitCubeMesh(int n)
        => CreateBoxMesh(n, n, n, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    private static void FixInvertedTetrahedra(SimplexMesh mesh, double[,] coords)
    {
        var tetsToFix = new List<(int index, int[] nodes)>();

        for (int i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            // Use shared MeshGeometry.ComputeTetrahedronJacobian instead of duplicate
            double jac = ComputeTetrahedronJacobian(coords, nodes[0], nodes[1], nodes[2], nodes[3]);

            if (jac <= 0)
                tetsToFix.Add((i, new[] { nodes[1], nodes[0], nodes[2], nodes[3] }));
        }

        tetsToFix.Sort((a, b) => b.index.CompareTo(a.index));

        foreach (var (index, swappedNodes) in tetsToFix)
        {
            mesh.Remove<Tet4>(index);
            mesh.AddTetrahedron(swappedNodes[0], swappedNodes[1], swappedNodes[2], swappedNodes[3]);
        }

        if (tetsToFix.Count > 0)
            Console.WriteLine($"[CreateBoxMesh] Fixed {tetsToFix.Count} inverted tetrahedra");
    }

    #endregion

    #region Delaunay Triangulation

    /// <summary>
    ///     Create Delaunay triangulation for a polygon with optional holes using Bowyer-Watson algorithm.
    /// </summary>
    public static (List<(double x, double y)> points, List<(int v0, int v1, int v2)> triangles)
        DelaunayTriangulate(
            double[,] outerBoundary,
            List<double[,]>? holes = null,
            double[,]? interiorPoints = null,
            double targetEdgeLength = 0.0)
    {
        Console.WriteLine("[DelaunayTriangulation] Starting...");
        Console.WriteLine($"  Outer boundary: {outerBoundary.GetLength(0)} points");

        var outerPts = ExtractCurve2D(outerBoundary);
        var holePts = new List<List<(double x, double y)>>();

        if (holes != null)
        {
            foreach (var hole in holes) holePts.Add(ExtractCurve2D(hole));
            Console.WriteLine($"  Holes: {holePts.Count}");
        }

        var points = new List<(double x, double y)>();
        points.AddRange(outerPts);

        foreach (var hole in holePts)
            points.AddRange(hole);

        // Generate interior points if needed
        if (interiorPoints == null && targetEdgeLength > 0)
        {
            var interiorList = GenerateInteriorGrid(outerBoundary, holes, targetEdgeLength);
            points.AddRange(interiorList);
            Console.WriteLine($"  Generated {interiorList.Count} interior points");
        }
        else if (interiorPoints != null)
        {
            for (var i = 0; i < interiorPoints.GetLength(0); i++)
                points.Add((interiorPoints[i, 0], interiorPoints[i, 1]));
        }
        else
        {
            var avgEdgeLength = ComputeAverageEdgeLength(outerPts);
            var autoSpacing = avgEdgeLength * 1.5;
            var interiorList = GenerateInteriorGrid(outerBoundary, holes, autoSpacing);
            points.AddRange(interiorList);
            Console.WriteLine($"  Generated {interiorList.Count} interior points");
        }

        var totalPoints = points.Count;

        // Compute bounding box and create super-triangle
        var (xMin, xMax, yMin, yMax) = ComputeBoundingBox2D(points);
        var margin = Math.Max(xMax - xMin, yMax - yMin) * 0.1;

        var superTriangle = CreateSuperTriangle(xMin - margin, xMax + margin, yMin - margin, yMax + margin);

        var st0 = points.Count;
        var st1 = st0 + 1;
        var st2 = st0 + 2;
        points.Add(superTriangle.v0);
        points.Add(superTriangle.v1);
        points.Add(superTriangle.v2);

        var triangles = new List<(int v0, int v1, int v2)> { (st0, st1, st2) };
        var circumcenters = new List<(double x, double y, double rSq)>();
        circumcenters.Add(ComputeCircumcircle(superTriangle.v0, superTriangle.v1, superTriangle.v2));

        Console.WriteLine($"  Inserting {totalPoints} points via Bowyer-Watson...");

        for (var i = 0; i < totalPoints; i++)
            AddPointBW(i, points, triangles, circumcenters);

        Console.WriteLine($"  → Triangulation complete: {triangles.Count} triangles");

        // Remove triangles connected to super-triangle
        var validTriangles = triangles
            .Where(t => t.v0 < totalPoints && t.v1 < totalPoints && t.v2 < totalPoints)
            .ToList();

        // Filter triangles outside polygon or inside holes
        validTriangles = FilterTriangles(validTriangles, points, outerPts, holePts);

        Console.WriteLine($"  → After filtering: {validTriangles.Count} triangles");

        points.RemoveRange(totalPoints, 3);

        return (points, validTriangles);
    }

    /// <summary>
    ///     Convert Delaunay result to SimplexMesh format.
    /// </summary>
    public static (SimplexMesh mesh, double[,] coords) DelaunayToSimplexMesh(
        List<(double x, double y)> points, List<(int v0, int v1, int v2)> triangles)
    {
        var mesh = new SimplexMesh();
        var coords = new double[points.Count, 3];

        mesh.WithBatch(() =>
        {
            for (int i = 0; i < points.Count; i++)
            {
                coords[i, 0] = points[i].x;
                coords[i, 1] = points[i].y;
                coords[i, 2] = 0;
                mesh.AddNode(i);
            }

            foreach (var (v0, v1, v2) in triangles)
                mesh.AddTriangle(v0, v1, v2);
        });

        return (mesh, coords);
    }

    #region Bowyer-Watson Implementation

    private static ((double x, double y) v0, (double x, double y) v1, (double x, double y) v2)
        CreateSuperTriangle(double xMin, double xMax, double yMin, double yMax)
    {
        var dx = xMax - xMin;
        var dy = yMax - yMin;
        var dMax = Math.Max(dx, dy);

        var cx = (xMin + xMax) / 2;
        var cy = (yMin + yMax) / 2;

        var scale = 3.0;
        return (
            (cx - scale * dMax, cy - scale * dMax),
            (cx + scale * dMax, cy - scale * dMax),
            (cx, cy + scale * dMax)
        );
    }

    private static (double x, double y, double rSq) ComputeCircumcircle(
        (double x, double y) p0, (double x, double y) p1, (double x, double y) p2)
    {
        var ax = p1.x - p0.x;
        var ay = p1.y - p0.y;
        var bx = p2.x - p0.x;
        var by = p2.y - p0.y;

        var d = 2 * (ax * by - ay * bx);

        if (Math.Abs(d) < Epsilon)
            return (0, 0, double.MaxValue);

        var aSq = ax * ax + ay * ay;
        var bSq = bx * bx + by * by;

        var ux = (by * aSq - ay * bSq) / d;
        var uy = (ax * bSq - bx * aSq) / d;

        return (p0.x + ux, p0.y + uy, ux * ux + uy * uy);
    }

    private static void AddPointBW(
        int pointIdx,
        List<(double x, double y)> points,
        List<(int v0, int v1, int v2)> triangles,
        List<(double x, double y, double rSq)> circumcenters)
    {
        var point = points[pointIdx];
        var badTriangles = new List<int>();

        for (var i = 0; i < triangles.Count; i++)
        {
            var (cx, cy, rSq) = circumcenters[i];
            var dx = point.x - cx;
            var dy = point.y - cy;
            var distSq = dx * dx + dy * dy;

            if (distSq < rSq - Epsilon) badTriangles.Add(i);
        }

        if (badTriangles.Count == 0) return;

        var polygon = new List<(int v0, int v1)>();

        foreach (var ti in badTriangles)
        {
            var (v0, v1, v2) = triangles[ti];
            var edges = new[] { (v0, v1), (v1, v2), (v2, v0) };

            foreach (var edge in edges)
            {
                var isShared = false;

                foreach (var tj in badTriangles)
                {
                    if (ti == tj) continue;

                    var (u0, u1, u2) = triangles[tj];

                    if (SharesEdge(edge, (u0, u1)) ||
                        SharesEdge(edge, (u1, u2)) ||
                        SharesEdge(edge, (u2, u0)))
                    {
                        isShared = true;
                        break;
                    }
                }

                if (!isShared) polygon.Add(edge);
            }
        }

        badTriangles.Sort();
        for (var i = badTriangles.Count - 1; i >= 0; i--)
        {
            var idx = badTriangles[i];
            triangles.RemoveAt(idx);
            circumcenters.RemoveAt(idx);
        }

        foreach (var (v0, v1) in polygon)
        {
            triangles.Add((pointIdx, v0, v1));
            circumcenters.Add(ComputeCircumcircle(points[pointIdx], points[v0], points[v1]));
        }
    }

    private static bool SharesEdge((int a, int b) edge1, (int a, int b) edge2)
        => (edge1.a == edge2.a && edge1.b == edge2.b) || (edge1.a == edge2.b && edge1.b == edge2.a);

    private static List<(int v0, int v1, int v2)> FilterTriangles(
        List<(int v0, int v1, int v2)> triangles,
        List<(double x, double y)> points,
        List<(double x, double y)> outerBoundary,
        List<List<(double x, double y)>> holes)
    {
        var filtered = new List<(int v0, int v1, int v2)>();

        var outerArray = ConvertToArray(outerBoundary);
        var holeArrays = holes.Select(h => ConvertToArray(h)).ToList();

        foreach (var tri in triangles)
        {
            var p0 = points[tri.v0];
            var p1 = points[tri.v1];
            var p2 = points[tri.v2];

            var centroid = ((p0.x + p1.x + p2.x) / 3.0, (p0.y + p1.y + p2.y) / 3.0);

            if (!IsPointInPolygon(centroid, outerArray))
                continue;

            var inHole = holeArrays.Any(hole => IsPointInPolygon(centroid, hole));
            if (inHole) continue;

            filtered.Add(tri);
        }

        return filtered;
    }

    #endregion

    #region Interior Point Generation

    /// <summary>
    ///     Generate interior grid points with hexagonal packing.
    /// </summary>
    public static List<(double x, double y)> GenerateInteriorGrid(
        double[,] outerBoundary, List<double[,]>? holes, double spacing)
    {
        var interiorPoints = new List<(double x, double y)>();

        var n = outerBoundary.GetLength(0);
        double minX = double.MaxValue, maxX = double.MinValue;
        double minY = double.MaxValue, maxY = double.MinValue;

        for (var i = 0; i < n; i++)
        {
            var x = outerBoundary[i, 0];
            var y = outerBoundary[i, 1];
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }

        var random = new Random(12345);
        var perturbation = spacing * GridPerturbationFactor;

        var rowIndex = 0;
        for (var y = minY + spacing; y < maxY; y += spacing * HexRowSpacing)
        {
            var xOffset = rowIndex % 2 == 1 ? spacing * 0.5 : 0.0;
            for (var x = minX + spacing + xOffset; x < maxX; x += spacing)
            {
                var px = x + (random.NextDouble() - 0.5) * 2 * perturbation;
                var py = y + (random.NextDouble() - 0.5) * 2 * perturbation;

                var point = (px, py);
                if (!IsPointInPolygon(point, outerBoundary))
                    continue;

                var inHole = false;
                if (holes != null)
                {
                    foreach (var hole in holes)
                    {
                        if (IsPointInPolygon(point, hole))
                        {
                            inHole = true;
                            break;
                        }
                    }
                }

                if (!inHole)
                    interiorPoints.Add(point);
            }

            rowIndex++;
        }

        return interiorPoints;
    }

    #endregion

    #endregion
}

// ═══════════════════════════════════════════════════════════════════════════
// Unified Mesher - High-Level API (from UnifiedMesher.cs)
// ═══════════════════════════════════════════════════════════════════════════

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
        var (points, triangles) = MeshGeneration.DelaunayTriangulate(outerBoundary, holes, interiorCoords, targetEdgeLength);

        // Convert to SimplexMesh
        var (mesh, coords) = MeshGeneration.DelaunayToSimplexMesh(points, triangles);

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

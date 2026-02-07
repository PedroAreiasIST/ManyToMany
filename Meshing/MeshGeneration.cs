// MeshGeneration.cs - Mesh generation algorithms
// Consolidates: Original MeshGeneration.cs + DelaunayTriangulation.cs
// Provides structured meshes (rectangles, boxes) and Delaunay triangulation
// License: GPLv3

using static Numerical.MeshConstants;
using static Numerical.MeshGeometry;
using static Numerical.GeometryUtilities;
using static Numerical.CurveUtilities;

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

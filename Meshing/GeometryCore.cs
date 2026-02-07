// GeometryCore.cs - Unified geometric operations for mesh library
// Consolidates: MeshGeometry.cs, GeometryUtilities.cs, CurveUtilities.cs
// Eliminates all duplicate geometric functions across the library
// License: GPLv3

using System.Text;
using static Numerical.MeshConstants;

namespace Numerical;

#region Element Geometry (from MeshGeometry.cs)

/// <summary>
///     Geometric operations for mesh elements: Jacobians, areas, volumes, quality metrics.
/// </summary>
public static class MeshGeometry
{
    #region 2D Triangle Operations

    /// <summary>
    ///     Compute 2D Jacobian determinant for a triangle.
    ///     Returns twice the signed area: positive = CCW, negative = CW.
    /// </summary>
    public static double ComputeTriangleJacobian(double[,] coords, int n0, int n1, int n2)
    {
        double x0 = coords[n0, 0], y0 = coords[n0, 1];
        double x1 = coords[n1, 0], y1 = coords[n1, 1];
        double x2 = coords[n2, 0], y2 = coords[n2, 1];
        return (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    }

    /// <summary>
    ///     Compute signed area of a triangle (half of Jacobian).
    /// </summary>
    public static double ComputeTriangleArea(double[,] coords, int n0, int n1, int n2)
        => 0.5 * Math.Abs(ComputeTriangleJacobian(coords, n0, n1, n2));

    /// <summary>
    ///     Check if triangle has counter-clockwise (CCW) orientation.
    /// </summary>
    public static bool IsTriangleCCW(double[,] coords, int n0, int n1, int n2)
        => ComputeTriangleJacobian(coords, n0, n1, n2) > Epsilon;

    /// <summary>
    ///     Check if triangle is degenerate (zero or near-zero area).
    /// </summary>
    public static bool IsTriangleDegenerate(double[,] coords, int n0, int n1, int n2, double tolerance = Epsilon)
        => Math.Abs(ComputeTriangleJacobian(coords, n0, n1, n2)) < tolerance;

    /// <summary>
    ///     Compute aspect ratio of a triangle.
    /// </summary>
    public static double ComputeTriangleAspectRatio(double[,] coords, int n0, int n1, int n2)
    {
        var area = ComputeTriangleArea(coords, n0, n1, n2);
        if (area < Epsilon) return double.PositiveInfinity;

        var e1 = EdgeLength2D(coords, n0, n1);
        var e2 = EdgeLength2D(coords, n1, n2);
        var e3 = EdgeLength2D(coords, n2, n0);

        var perimeter = e1 + e2 + e3;
        var inradius = 2.0 * area / perimeter;
        var maxEdge = Math.Max(e1, Math.Max(e2, e3));

        return maxEdge / (2.0 * inradius);
    }

    /// <summary>
    ///     Compute minimum interior angle of a triangle in radians.
    /// </summary>
    public static double ComputeTriangleMinAngle(double[,] coords, int n0, int n1, int n2)
    {
        var e1 = EdgeLength2D(coords, n0, n1);
        var e2 = EdgeLength2D(coords, n1, n2);
        var e3 = EdgeLength2D(coords, n2, n0);

        var angle0 = Math.Acos(Math.Clamp((e1 * e1 + e3 * e3 - e2 * e2) / (2 * e1 * e3), -1, 1));
        var angle1 = Math.Acos(Math.Clamp((e1 * e1 + e2 * e2 - e3 * e3) / (2 * e1 * e2), -1, 1));
        var angle2 = Math.Acos(Math.Clamp((e2 * e2 + e3 * e3 - e1 * e1) / (2 * e2 * e3), -1, 1));

        return Math.Min(angle0, Math.Min(angle1, angle2));
    }

    #endregion

    #region 3D Tetrahedron Operations

    /// <summary>
    ///     Compute 3D Jacobian determinant for a tetrahedron.
    ///     Returns 6 times the signed volume.
    /// </summary>
    public static double ComputeTetrahedronJacobian(double[,] coords, int n0, int n1, int n2, int n3)
    {
        var v1x = coords[n1, 0] - coords[n0, 0];
        var v1y = coords[n1, 1] - coords[n0, 1];
        var v1z = coords[n1, 2] - coords[n0, 2];

        var v2x = coords[n2, 0] - coords[n0, 0];
        var v2y = coords[n2, 1] - coords[n0, 1];
        var v2z = coords[n2, 2] - coords[n0, 2];

        var v3x = coords[n3, 0] - coords[n0, 0];
        var v3y = coords[n3, 1] - coords[n0, 1];
        var v3z = coords[n3, 2] - coords[n0, 2];

        return v1x * (v2y * v3z - v2z * v3y) -
               v1y * (v2x * v3z - v2z * v3x) +
               v1z * (v2x * v3y - v2y * v3x);
    }

    /// <summary>
    ///     Compute volume of a tetrahedron (1/6 of Jacobian absolute value).
    /// </summary>
    public static double ComputeTetrahedronVolume(double[,] coords, int n0, int n1, int n2, int n3)
        => Math.Abs(ComputeTetrahedronJacobian(coords, n0, n1, n2, n3)) / 6.0;

    /// <summary>
    ///     Check if tetrahedron has correct orientation (positive Jacobian).
    /// </summary>
    public static bool IsTetrahedronCorrectOrientation(double[,] coords, int n0, int n1, int n2, int n3)
        => ComputeTetrahedronJacobian(coords, n0, n1, n2, n3) > Epsilon;

    /// <summary>
    ///     Check if tetrahedron is degenerate (zero or near-zero volume).
    /// </summary>
    public static bool IsTetrahedronDegenerate(double[,] coords, int n0, int n1, int n2, int n3, double tolerance = Epsilon)
        => Math.Abs(ComputeTetrahedronJacobian(coords, n0, n1, n2, n3)) < tolerance;

    /// <summary>
    ///     Compute aspect ratio of a tetrahedron.
    /// </summary>
    public static double ComputeTetrahedronAspectRatio(double[,] coords, int n0, int n1, int n2, int n3)
    {
        var volume = ComputeTetrahedronVolume(coords, n0, n1, n2, n3);
        if (volume < Epsilon) return double.PositiveInfinity;

        var e01 = EdgeLength3D(coords, n0, n1);
        var e02 = EdgeLength3D(coords, n0, n2);
        var e03 = EdgeLength3D(coords, n0, n3);
        var e12 = EdgeLength3D(coords, n1, n2);
        var e13 = EdgeLength3D(coords, n1, n3);
        var e23 = EdgeLength3D(coords, n2, n3);

        var maxEdge = Math.Max(e01, Math.Max(e02, Math.Max(e03, Math.Max(e12, Math.Max(e13, e23)))));
        return maxEdge * maxEdge * maxEdge / (8.48 * volume);
    }

    #endregion

    #region 2D Quadrilateral Operations

    /// <summary>
    ///     Check if quadrilateral has counter-clockwise (CCW) orientation.
    /// </summary>
    public static bool IsQuadCCW(double[,] coords, int n0, int n1, int n2, int n3)
    {
        var ax = coords[n1, 0] - coords[n0, 0];
        var ay = coords[n1, 1] - coords[n0, 1];
        var bx = coords[n3, 0] - coords[n0, 0];
        var by = coords[n3, 1] - coords[n0, 1];
        return ax * by - ay * bx > Epsilon;
    }

    /// <summary>
    ///     Compute area of a quadrilateral (sum of two triangles).
    /// </summary>
    public static double ComputeQuadArea(double[,] coords, int n0, int n1, int n2, int n3)
        => ComputeTriangleArea(coords, n0, n1, n2) + ComputeTriangleArea(coords, n0, n2, n3);

    /// <summary>
    ///     Check if quadrilateral is convex.
    /// </summary>
    public static bool IsQuadConvex(double[,] coords, int n0, int n1, int n2, int n3)
    {
        var cross0 = CrossProductAtVertex(coords, n0, n1, n2);
        var cross1 = CrossProductAtVertex(coords, n1, n2, n3);
        var cross2 = CrossProductAtVertex(coords, n2, n3, n0);
        var cross3 = CrossProductAtVertex(coords, n3, n0, n1);

        return (cross0 > 0 && cross1 > 0 && cross2 > 0 && cross3 > 0) ||
               (cross0 < 0 && cross1 < 0 && cross2 < 0 && cross3 < 0);
    }

    private static double CrossProductAtVertex(double[,] coords, int n0, int n1, int n2)
    {
        var ax = coords[n1, 0] - coords[n0, 0];
        var ay = coords[n1, 1] - coords[n0, 1];
        var bx = coords[n2, 0] - coords[n1, 0];
        var by = coords[n2, 1] - coords[n1, 1];
        return ax * by - ay * bx;
    }

    #endregion

    #region Edge Length Functions (PUBLIC - used by many modules)

    /// <summary>
    ///     Compute 2D edge length between two nodes.
    /// </summary>
    public static double EdgeLength2D(double[,] coords, int n0, int n1)
    {
        var dx = coords[n1, 0] - coords[n0, 0];
        var dy = coords[n1, 1] - coords[n0, 1];
        return Math.Sqrt(dx * dx + dy * dy);
    }

    /// <summary>
    ///     Compute 3D edge length between two nodes.
    /// </summary>
    public static double EdgeLength3D(double[,] coords, int n0, int n1)
    {
        var dx = coords[n1, 0] - coords[n0, 0];
        var dy = coords[n1, 1] - coords[n0, 1];
        var dz = coords[n1, 2] - coords[n0, 2];
        return Math.Sqrt(dx * dx + dy * dy + dz * dz);
    }

    #endregion

    #region Mesh-wide Operations

    /// <summary>
    ///     Validate mesh quality: check all elements have positive Jacobians.
    /// </summary>
    public static int ValidateMeshOrientation(SimplexMesh mesh, double[,] coords,
        out int invertedTris, out int invertedTets)
    {
        invertedTris = 0;
        invertedTets = 0;

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            if (!IsTriangleCCW(coords, nodes[0], nodes[1], nodes[2]))
                invertedTris++;
        }

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            if (!IsTetrahedronCorrectOrientation(coords, nodes[0], nodes[1], nodes[2], nodes[3]))
                invertedTets++;
        }

        return invertedTris + invertedTets;
    }

    /// <summary>
    ///     Find degenerate elements.
    /// </summary>
    public static (int[] degenerateTris, int[] degenerateTets) FindDegenerateElements(
        SimplexMesh mesh, double[,] coords, double tolerance = Epsilon)
    {
        var badTris = new List<int>();
        var badTets = new List<int>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            if (IsTriangleDegenerate(coords, nodes[0], nodes[1], nodes[2], tolerance))
                badTris.Add(i);
        }

        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            if (IsTetrahedronDegenerate(coords, nodes[0], nodes[1], nodes[2], nodes[3], tolerance))
                badTets.Add(i);
        }

        return (badTris.ToArray(), badTets.ToArray());
    }

    /// <summary>
    ///     Compute mesh quality statistics.
    /// </summary>
    public static MeshQualityStats ComputeQualityStatistics(SimplexMesh mesh, double[,] coords)
    {
        var stats = new MeshQualityStats();

        if (mesh.Count<Tri3>() > 0)
        {
            var minAspect = double.PositiveInfinity;
            double maxAspect = 0, sumAspect = 0;
            var minAngle = double.PositiveInfinity;
            double maxAngle = 0;

            for (var i = 0; i < mesh.Count<Tri3>(); i++)
            {
                var nodes = mesh.NodesOf<Tri3, Node>(i);
                var aspect = ComputeTriangleAspectRatio(coords, nodes[0], nodes[1], nodes[2]);
                minAspect = Math.Min(minAspect, aspect);
                maxAspect = Math.Max(maxAspect, aspect);
                sumAspect += aspect;

                var angle = ComputeTriangleMinAngle(coords, nodes[0], nodes[1], nodes[2]);
                minAngle = Math.Min(minAngle, angle);
                maxAngle = Math.Max(maxAngle, angle);
            }

            stats.TriangleCount = mesh.Count<Tri3>();
            stats.MinTriangleAspectRatio = minAspect;
            stats.MaxTriangleAspectRatio = maxAspect;
            stats.AvgTriangleAspectRatio = sumAspect / mesh.Count<Tri3>();
            stats.MinTriangleAngleDegrees = minAngle * 180.0 / Math.PI;
        }

        if (mesh.Count<Tet4>() > 0)
        {
            var minAspect = double.PositiveInfinity;
            double maxAspect = 0, sumAspect = 0;

            for (var i = 0; i < mesh.Count<Tet4>(); i++)
            {
                var nodes = mesh.NodesOf<Tet4, Node>(i);
                var aspect = ComputeTetrahedronAspectRatio(coords, nodes[0], nodes[1], nodes[2], nodes[3]);
                minAspect = Math.Min(minAspect, aspect);
                maxAspect = Math.Max(maxAspect, aspect);
                sumAspect += aspect;
            }

            stats.TetrahedronCount = mesh.Count<Tet4>();
            stats.MinTetrahedronAspectRatio = minAspect;
            stats.MaxTetrahedronAspectRatio = maxAspect;
            stats.AvgTetrahedronAspectRatio = sumAspect / mesh.Count<Tet4>();
        }

        return stats;
    }

    #endregion
}

/// <summary>
///     Container for mesh quality statistics.
/// </summary>
public class MeshQualityStats
{
    public int TriangleCount { get; set; }
    public int TetrahedronCount { get; set; }
    public double MinTriangleAspectRatio { get; set; }
    public double MaxTriangleAspectRatio { get; set; }
    public double AvgTriangleAspectRatio { get; set; }
    public double MinTetrahedronAspectRatio { get; set; }
    public double MaxTetrahedronAspectRatio { get; set; }
    public double AvgTetrahedronAspectRatio { get; set; }
    public double MinTriangleAngleDegrees { get; set; }

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine("=== Mesh Quality Statistics ===");
        if (TriangleCount > 0)
        {
            sb.AppendLine($"Triangles: {TriangleCount}");
            sb.AppendLine($"  Aspect Ratio: {MinTriangleAspectRatio:F2} - {MaxTriangleAspectRatio:F2} (avg: {AvgTriangleAspectRatio:F2})");
            sb.AppendLine($"  Min Angle: {MinTriangleAngleDegrees:F1}°");
        }
        if (TetrahedronCount > 0)
        {
            sb.AppendLine($"Tetrahedra: {TetrahedronCount}");
            sb.AppendLine($"  Aspect Ratio: {MinTetrahedronAspectRatio:F2} - {MaxTetrahedronAspectRatio:F2} (avg: {AvgTetrahedronAspectRatio:F2})");
        }
        return sb.ToString();
    }
}

#endregion

#region Point/Polygon Geometry (from GeometryUtilities.cs)

/// <summary>
///     Fundamental 2D/3D geometric operations: point-in-polygon, distances, projections.
/// </summary>
public static class GeometryUtilities
{
    #region Point-in-Polygon Tests

    /// <summary>
    ///     Check if point is inside polygon using ray casting algorithm.
    /// </summary>
    public static bool IsPointInPolygon((double x, double y) point, double[,] polygon)
    {
        var n = polygon.GetLength(0);
        var inside = false;

        for (int i = 0, j = n - 1; i < n; j = i++)
        {
            double xi = polygon[i, 0], yi = polygon[i, 1];
            double xj = polygon[j, 0], yj = polygon[j, 1];

            if (yi > point.y != yj > point.y &&
                point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi)
                inside = !inside;
        }

        return inside;
    }

    /// <summary>
    ///     Check if point is inside polygon (list format).
    /// </summary>
    public static bool IsPointInPolygon((double x, double y) point, List<(double x, double y)> polygon)
    {
        var n = polygon.Count;
        var inside = false;

        for (int i = 0, j = n - 1; i < n; j = i++)
        {
            double xi = polygon[i].x, yi = polygon[i].y;
            double xj = polygon[j].x, yj = polygon[j].y;

            if (yi > point.y != yj > point.y &&
                point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi)
                inside = !inside;
        }

        return inside;
    }

    /// <summary>
    ///     Check if point is on polygon boundary.
    /// </summary>
    public static bool IsPointOnPolygonBoundary((double x, double y) point, double[,] polygon, double tolerance = Epsilon)
    {
        var n = polygon.GetLength(0);
        for (var i = 0; i < n; i++)
        {
            var j = (i + 1) % n;
            var dist = DistancePointToSegment(point, (polygon[i, 0], polygon[i, 1]), (polygon[j, 0], polygon[j, 1]));
            if (dist < tolerance) return true;
        }
        return false;
    }

    #endregion

    #region Distance Calculations

    /// <summary>
    ///     Compute distance from point to infinite line.
    /// </summary>
    public static double DistancePointToLine((double x, double y) point, (double x, double y) lineStart, (double x, double y) lineEnd)
    {
        var dx = lineEnd.x - lineStart.x;
        var dy = lineEnd.y - lineStart.y;
        var lineLength = Math.Sqrt(dx * dx + dy * dy);

        if (lineLength < Epsilon) return Distance2D(point, lineStart);

        var cross = Math.Abs((lineEnd.x - lineStart.x) * (lineStart.y - point.y) -
                             (lineStart.x - point.x) * (lineEnd.y - lineStart.y));
        return cross / lineLength;
    }

    /// <summary>
    ///     Compute distance from point to line segment.
    /// </summary>
    public static double DistancePointToSegment((double x, double y) point, (double x, double y) segStart, (double x, double y) segEnd)
    {
        var dx = segEnd.x - segStart.x;
        var dy = segEnd.y - segStart.y;
        var lengthSquared = dx * dx + dy * dy;

        if (lengthSquared < Epsilon) return Distance2D(point, segStart);

        var t = Math.Clamp(((point.x - segStart.x) * dx + (point.y - segStart.y) * dy) / lengthSquared, 0, 1);
        return Distance2D(point, (segStart.x + t * dx, segStart.y + t * dy));
    }

    /// <summary>
    ///     Compute 2D Euclidean distance between two points.
    /// </summary>
    public static double Distance2D((double x, double y) p1, (double x, double y) p2)
    {
        var dx = p2.x - p1.x;
        var dy = p2.y - p1.y;
        return Math.Sqrt(dx * dx + dy * dy);
    }

    /// <summary>
    ///     Compute 3D Euclidean distance between two points.
    /// </summary>
    public static double Distance3D((double x, double y, double z) p1, (double x, double y, double z) p2)
    {
        var dx = p2.x - p1.x;
        var dy = p2.y - p1.y;
        var dz = p2.z - p1.z;
        return Math.Sqrt(dx * dx + dy * dy + dz * dz);
    }

    #endregion

    #region Point Projections

    /// <summary>
    ///     Project point onto infinite line.
    /// </summary>
    public static (double x, double y) ProjectPointToLine((double x, double y) point, (double x, double y) lineStart, (double x, double y) lineEnd)
    {
        var dx = lineEnd.x - lineStart.x;
        var dy = lineEnd.y - lineStart.y;
        var lengthSquared = dx * dx + dy * dy;

        if (lengthSquared < Epsilon) return lineStart;

        var t = ((point.x - lineStart.x) * dx + (point.y - lineStart.y) * dy) / lengthSquared;
        return (lineStart.x + t * dx, lineStart.y + t * dy);
    }

    /// <summary>
    ///     Project point onto line segment.
    /// </summary>
    public static (double x, double y) ProjectPointToSegment((double x, double y) point, (double x, double y) segStart, (double x, double y) segEnd)
    {
        var dx = segEnd.x - segStart.x;
        var dy = segEnd.y - segStart.y;
        var lengthSquared = dx * dx + dy * dy;

        if (lengthSquared < Epsilon) return segStart;

        var t = Math.Clamp(((point.x - segStart.x) * dx + (point.y - segStart.y) * dy) / lengthSquared, 0, 1);
        return (segStart.x + t * dx, segStart.y + t * dy);
    }

    #endregion

    #region Line Intersections

    /// <summary>
    ///     Find intersection point of two infinite lines.
    /// </summary>
    public static ((double x, double y) point, bool intersects) IntersectLines(
        (double x, double y) line1Start, (double x, double y) line1End,
        (double x, double y) line2Start, (double x, double y) line2End)
    {
        double x1 = line1Start.x, y1 = line1Start.y, x2 = line1End.x, y2 = line1End.y;
        double x3 = line2Start.x, y3 = line2Start.y, x4 = line2End.x, y4 = line2End.y;

        var denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (Math.Abs(denom) < Epsilon) return ((0, 0), false);

        var t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
        return ((x1 + t * (x2 - x1), y1 + t * (y2 - y1)), true);
    }

    /// <summary>
    ///     Find intersection point of two line segments.
    /// </summary>
    public static ((double x, double y) point, bool intersects) IntersectSegments(
        (double x, double y) seg1Start, (double x, double y) seg1End,
        (double x, double y) seg2Start, (double x, double y) seg2End)
    {
        double x1 = seg1Start.x, y1 = seg1Start.y, x2 = seg1End.x, y2 = seg1End.y;
        double x3 = seg2Start.x, y3 = seg2Start.y, x4 = seg2End.x, y4 = seg2End.y;

        var denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (Math.Abs(denom) < Epsilon) return ((0, 0), false);

        var t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
        var u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;

        if (t >= 0 && t <= 1 && u >= 0 && u <= 1)
            return ((x1 + t * (x2 - x1), y1 + t * (y2 - y1)), true);

        return ((0, 0), false);
    }

    #endregion

    #region Bounding Box Operations

    /// <summary>
    ///     Compute axis-aligned bounding box.
    /// </summary>
    public static (double xMin, double xMax, double yMin, double yMax) ComputeBoundingBox2D(double[,] points)
    {
        var n = points.GetLength(0);
        if (n == 0) return (0, 0, 0, 0);

        double xMin = points[0, 0], xMax = points[0, 0];
        double yMin = points[0, 1], yMax = points[0, 1];

        for (var i = 1; i < n; i++)
        {
            xMin = Math.Min(xMin, points[i, 0]);
            xMax = Math.Max(xMax, points[i, 0]);
            yMin = Math.Min(yMin, points[i, 1]);
            yMax = Math.Max(yMax, points[i, 1]);
        }

        return (xMin, xMax, yMin, yMax);
    }

    /// <summary>
    ///     Compute axis-aligned bounding box (list format).
    /// </summary>
    public static (double xMin, double xMax, double yMin, double yMax) ComputeBoundingBox2D(List<(double x, double y)> points)
    {
        if (points.Count == 0) return (0, 0, 0, 0);

        double xMin = points[0].x, xMax = points[0].x;
        double yMin = points[0].y, yMax = points[0].y;

        for (var i = 1; i < points.Count; i++)
        {
            xMin = Math.Min(xMin, points[i].x);
            xMax = Math.Max(xMax, points[i].x);
            yMin = Math.Min(yMin, points[i].y);
            yMax = Math.Max(yMax, points[i].y);
        }

        return (xMin, xMax, yMin, yMax);
    }

    /// <summary>
    ///     Check if point is inside bounding box.
    /// </summary>
    public static bool IsPointInBoundingBox((double x, double y) point, (double xMin, double xMax, double yMin, double yMax) bbox)
        => point.x >= bbox.xMin && point.x <= bbox.xMax && point.y >= bbox.yMin && point.y <= bbox.yMax;

    #endregion

    #region Angular and Vector Operations

    /// <summary>
    ///     Compute angle between two 2D vectors in radians [0, π].
    /// </summary>
    public static double AngleBetweenVectors2D((double x, double y) v1, (double x, double y) v2)
    {
        var dot = v1.x * v2.x + v1.y * v2.y;
        var mag1 = Math.Sqrt(v1.x * v1.x + v1.y * v1.y);
        var mag2 = Math.Sqrt(v2.x * v2.x + v2.y * v2.y);

        if (mag1 < Epsilon || mag2 < Epsilon) return 0;

        return Math.Acos(Math.Clamp(dot / (mag1 * mag2), -1.0, 1.0));
    }

    /// <summary>
    ///     Compute signed angle from v1 to v2 in radians [-π, π].
    /// </summary>
    public static double SignedAngleBetweenVectors2D((double x, double y) v1, (double x, double y) v2)
    {
        var angle = AngleBetweenVectors2D(v1, v2);
        var cross = v1.x * v2.y - v1.y * v2.x;
        return cross < 0 ? -angle : angle;
    }

    /// <summary>
    ///     Compute 2D cross product (scalar).
    /// </summary>
    public static double CrossProduct2D((double x, double y) v1, (double x, double y) v2)
        => v1.x * v2.y - v1.y * v2.x;

    /// <summary>
    ///     Compute 2D dot product.
    /// </summary>
    public static double DotProduct2D((double x, double y) v1, (double x, double y) v2)
        => v1.x * v2.x + v1.y * v2.y;

    /// <summary>
    ///     Normalize 2D vector to unit length.
    /// </summary>
    public static (double x, double y) Normalize2D((double x, double y) v)
    {
        var length = Math.Sqrt(v.x * v.x + v.y * v.y);
        if (length < Epsilon) return (0, 0);
        return (v.x / length, v.y / length);
    }

    /// <summary>
    ///     Rotate 2D vector by angle (radians, counter-clockwise).
    /// </summary>
    public static (double x, double y) Rotate2D((double x, double y) v, double angle)
    {
        var cos = Math.Cos(angle);
        var sin = Math.Sin(angle);
        return (v.x * cos - v.y * sin, v.x * sin + v.y * cos);
    }

    #endregion

    #region Interpolation

    /// <summary>
    ///     Linear interpolation between two points.
    /// </summary>
    public static (double x, double y) Lerp2D((double x, double y) p1, (double x, double y) p2, double t)
        => (p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y));

    /// <summary>
    ///     Barycentric interpolation in a triangle.
    /// </summary>
    public static (double x, double y) BarycentricInterpolation(
        (double x, double y) p0, (double x, double y) p1, (double x, double y) p2,
        double w0, double w1, double w2)
        => (w0 * p0.x + w1 * p1.x + w2 * p2.x, w0 * p0.y + w1 * p1.y + w2 * p2.y);

    /// <summary>
    ///     Compute barycentric coordinates of point in triangle.
    /// </summary>
    public static (double w0, double w1, double w2) ComputeBarycentricCoordinates(
        (double x, double y) point, (double x, double y) p0, (double x, double y) p1, (double x, double y) p2)
    {
        var v0x = p1.x - p0.x;
        var v0y = p1.y - p0.y;
        var v1x = p2.x - p0.x;
        var v1y = p2.y - p0.y;
        var v2x = point.x - p0.x;
        var v2y = point.y - p0.y;

        var d00 = v0x * v0x + v0y * v0y;
        var d01 = v0x * v1x + v0y * v1y;
        var d11 = v1x * v1x + v1y * v1y;
        var d20 = v2x * v0x + v2y * v0y;
        var d21 = v2x * v1x + v2y * v1y;

        var denom = d00 * d11 - d01 * d01;
        if (Math.Abs(denom) < Epsilon) return (0, 0, 0);

        var w1 = (d11 * d20 - d01 * d21) / denom;
        var w2 = (d00 * d21 - d01 * d20) / denom;
        return (1.0 - w1 - w2, w1, w2);
    }

    #endregion
}

#endregion

#region Curve/Boundary Utilities (from CurveUtilities.cs)

/// <summary>
///     Curve and boundary manipulation: refinement, resampling, orientation.
/// </summary>
public static class CurveUtilities
{
    #region Boundary Refinement

    /// <summary>
    ///     Refines boundary edges to ensure good mesh quality.
    /// </summary>
    public static double[,] RefineBoundary(double[,] boundary, double maxEdgeRatio = 2.0, int minEdgeCount = 200)
    {
        var n = boundary.GetLength(0);
        var dim = boundary.GetLength(1);

        if (n < 3) return boundary;

        var edgeLengths = new double[n];
        var totalLength = 0.0;
        var minLength = double.MaxValue;

        for (var i = 0; i < n; i++)
        {
            var j = (i + 1) % n;
            var length = Distance2DArray(boundary, i, j);
            edgeLengths[i] = length;
            totalLength += length;
            if (length > Epsilon && length < minLength) minLength = length;
        }

        if (totalLength < Epsilon) return boundary;

        var targetFromCount = totalLength / minEdgeCount;
        var targetFromRatio = maxEdgeRatio * minLength;
        var targetLength = Math.Min(targetFromCount, targetFromRatio);

        var refinedPoints = new List<double[]>();

        for (var i = 0; i < n; i++)
        {
            var j = (i + 1) % n;

            var currentPoint = new double[dim];
            for (var d = 0; d < dim; d++) currentPoint[d] = boundary[i, d];
            refinedPoints.Add(currentPoint);

            var edgeLength = edgeLengths[i];
            if (edgeLength < Epsilon) continue;

            if (edgeLength > targetLength * 1.01)
            {
                var nSegments = Math.Max(1, (int)Math.Round(edgeLength / targetLength));

                for (var k = 1; k < nSegments; k++)
                {
                    var t = (double)k / nSegments;
                    var interpPoint = new double[dim];
                    for (var d = 0; d < dim; d++)
                        interpPoint[d] = boundary[i, d] + t * (boundary[j, d] - boundary[i, d]);
                    refinedPoints.Add(interpPoint);
                }
            }
        }

        var nRefined = refinedPoints.Count;
        var result = new double[nRefined, dim];
        for (var i = 0; i < nRefined; i++)
            for (var d = 0; d < dim; d++)
                result[i, d] = refinedPoints[i][d];

        return result;
    }

    private static double Distance2DArray(double[,] coords, int i, int j)
    {
        var dx = coords[j, 0] - coords[i, 0];
        var dy = coords[j, 1] - coords[i, 1];
        return Math.Sqrt(dx * dx + dy * dy);
    }

    #endregion

    #region Curve Resampling

    /// <summary>
    ///     Resample curve to have specified number of points uniformly distributed by arc length.
    /// </summary>
    public static List<(double x, double y)> ResampleCurve(List<(double x, double y)> curve, int targetPoints, bool closedCurve = false)
    {
        if (curve.Count >= targetPoints) return new List<(double x, double y)>(curve);

        double totalLength = 0;
        for (var i = 0; i < curve.Count - 1; i++)
        {
            var dx = curve[i + 1].x - curve[i].x;
            var dy = curve[i + 1].y - curve[i].y;
            totalLength += Math.Sqrt(dx * dx + dy * dy);
        }

        if (closedCurve)
        {
            var dx = curve[0].x - curve[^1].x;
            var dy = curve[0].y - curve[^1].y;
            totalLength += Math.Sqrt(dx * dx + dy * dy);
        }

        if (totalLength < Epsilon) return new List<(double x, double y)>(curve);

        var targetSpacing = totalLength / (closedCurve ? targetPoints : targetPoints - 1);
        var resampled = new List<(double x, double y)> { curve[0] };

        double accumulatedLength = 0;
        var nextTarget = targetSpacing;

        for (var i = 0; i < curve.Count - 1; i++)
        {
            var dx = curve[i + 1].x - curve[i].x;
            var dy = curve[i + 1].y - curve[i].y;
            var segLength = Math.Sqrt(dx * dx + dy * dy);

            while (accumulatedLength + segLength >= nextTarget && resampled.Count < targetPoints - 1)
            {
                var t = (nextTarget - accumulatedLength) / segLength;
                resampled.Add((curve[i].x + t * dx, curve[i].y + t * dy));
                nextTarget += targetSpacing;
            }

            accumulatedLength += segLength;
        }

        if (!closedCurve || resampled.Count < targetPoints) resampled.Add(curve[^1]);

        return resampled;
    }

    /// <summary>
    ///     Resample curve from array format.
    /// </summary>
    public static double[,] ResampleCurveArray(double[,] curve, int targetPoints, bool closedCurve = false)
    {
        var curveList = ExtractCurve2D(curve);
        var resampled = ResampleCurve(curveList, targetPoints, closedCurve);

        var result = new double[resampled.Count, 2];
        for (var i = 0; i < resampled.Count; i++)
        {
            result[i, 0] = resampled[i].x;
            result[i, 1] = resampled[i].y;
        }

        return result;
    }

    #endregion

    #region Arc Length and Geometry

    /// <summary>
    ///     Compute total arc length of a curve.
    /// </summary>
    public static double ComputeArcLength(List<(double x, double y)> curve, bool closedCurve = false)
    {
        double length = 0;

        for (var i = 0; i < curve.Count - 1; i++)
        {
            var dx = curve[i + 1].x - curve[i].x;
            var dy = curve[i + 1].y - curve[i].y;
            length += Math.Sqrt(dx * dx + dy * dy);
        }

        if (closedCurve && curve.Count > 0)
        {
            var dx = curve[0].x - curve[^1].x;
            var dy = curve[0].y - curve[^1].y;
            length += Math.Sqrt(dx * dx + dy * dy);
        }

        return length;
    }

    /// <summary>
    ///     Compute average edge length of a boundary.
    /// </summary>
    public static double ComputeAverageEdgeLength(double[,] boundary)
    {
        var n = boundary.GetLength(0);
        if (n < 2) return 0;

        double totalLength = 0;
        for (var i = 0; i < n; i++)
        {
            var j = (i + 1) % n;
            totalLength += Distance2DArray(boundary, i, j);
        }

        return totalLength / n;
    }

    /// <summary>
    ///     Compute average edge length of a boundary (list format).
    /// </summary>
    public static double ComputeAverageEdgeLength(List<(double x, double y)> points)
    {
        var n = points.Count;
        if (n < 2) return 0;

        double totalLength = 0;
        for (var i = 0; i < n; i++)
        {
            var p1 = points[i];
            var p2 = points[(i + 1) % n];
            var dx = p2.x - p1.x;
            var dy = p2.y - p1.y;
            totalLength += Math.Sqrt(dx * dx + dy * dy);
        }

        return totalLength / n;
    }

    /// <summary>
    ///     Compute centroid of a boundary.
    /// </summary>
    public static (double x, double y) ComputeCentroid(double[,] boundary)
    {
        var n = boundary.GetLength(0);
        if (n == 0) return (0, 0);

        double sumX = 0, sumY = 0;
        for (var i = 0; i < n; i++)
        {
            sumX += boundary[i, 0];
            sumY += boundary[i, 1];
        }

        return (sumX / n, sumY / n);
    }

    /// <summary>
    ///     Compute signed area of a boundary (2D). Positive = CCW, negative = CW.
    /// </summary>
    public static double ComputeSignedArea(double[,] boundary)
    {
        var n = boundary.GetLength(0);
        if (n < 3) return 0;

        double area = 0;
        for (var i = 0; i < n; i++)
        {
            var j = (i + 1) % n;
            area += boundary[i, 0] * boundary[j, 1] - boundary[j, 0] * boundary[i, 1];
        }

        return area / 2.0;
    }

    #endregion

    #region Orientation

    /// <summary>
    ///     Check if boundary has counter-clockwise (CCW) orientation.
    /// </summary>
    public static bool IsBoundaryCCW(double[,] boundary)
        => ComputeSignedArea(boundary) > 0;

    /// <summary>
    ///     Reverse boundary orientation.
    /// </summary>
    public static double[,] ReverseBoundary(double[,] boundary)
    {
        var n = boundary.GetLength(0);
        var dim = boundary.GetLength(1);
        var reversed = new double[n, dim];

        for (var i = 0; i < n; i++)
            for (var d = 0; d < dim; d++)
                reversed[i, d] = boundary[n - 1 - i, d];

        return reversed;
    }

    /// <summary>
    ///     Ensure boundary has CCW orientation.
    /// </summary>
    public static double[,] EnsureCCW(double[,] boundary)
        => IsBoundaryCCW(boundary) ? boundary : ReverseBoundary(boundary);

    /// <summary>
    ///     Ensure boundary has CW orientation.
    /// </summary>
    public static double[,] EnsureCW(double[,] boundary)
        => IsBoundaryCCW(boundary) ? ReverseBoundary(boundary) : boundary;

    #endregion

    #region Conversion Utilities

    /// <summary>
    ///     Extract 2D curve from array format to list of tuples.
    /// </summary>
    public static List<(double x, double y)> ExtractCurve2D(double[,] curve)
    {
        var points = new List<(double x, double y)>();
        var n = curve.GetLength(0);
        for (var i = 0; i < n; i++) points.Add((curve[i, 0], curve[i, 1]));
        return points;
    }

    /// <summary>
    ///     Convert list of tuples to array format.
    /// </summary>
    public static double[,] ConvertToArray(List<(double x, double y)> curve, int dim = 2)
    {
        var n = curve.Count;
        var result = new double[n, dim];

        for (var i = 0; i < n; i++)
        {
            result[i, 0] = curve[i].x;
            result[i, 1] = curve[i].y;
            for (var d = 2; d < dim; d++) result[i, d] = 0;
        }

        return result;
    }

    #endregion
}

#endregion

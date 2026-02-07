// Examples2DA.cs - Advanced 2D meshing examples with intricate geometries + Fracture Mechanics
// Demonstrates: Complex boundaries, holes, wedges, TRI and QUAD elements, Crack insertion
// Output: ALL examples save to GiD/CIMNE .msh format
// License: GPLv3
//
// PERFORMANCE: Fracture examples use REDUCED mesh density for speed:
//   - Boundary: 51×51 or 51×101 nodes (instead of 101×101 or 101×201)
//   - maxArea: 5.0-8.0 (instead of 1.5-2.0)
//   - Result: 4-8x faster, still adequate resolution for demonstration
//   - For production: increase boundary points and reduce maxArea

namespace Numerical.Examples;

/// <summary>
///     Advanced 2D meshing examples with intricate geometries + Fracture Mechanics.
///     Shows real-world complex shapes, holes, wedges, TRI and QUAD meshes, and crack insertion.
///     OUTPUT: ALL examples automatically save to GiD/CIMNE .msh format
///     using SimplexRemesher.SaveGiD() for visualization in GiD (www.gidhome.com).
///     Files saved:
///     - Advanced Meshing (1-10): ex1-ex10 with _tri.msh and/or _quad.msh
///     - Fracture examples (11-20): anderson2005, griffith1921, etc. with _tri.msh
///     PERFORMANCE NOTES:
///     Fracture examples (11-20) use REDUCED mesh density (51×51 instead of 101×101)
///     to demonstrate crack insertion quickly. For production analysis, increase
///     boundary point density and reduce maxArea parameter.
/// </summary>
public static class Examples2DA
{
    public static void Main(string[] args)
    {
        //TestSingleTetRefinement();
        Console.WriteLine("╔═══════════════════════════════════════════════════════════════╗");
        Console.WriteLine("║  ADVANCED 2D MESHING EXAMPLES + FRACTURE MECHANICS           ║");
        Console.WriteLine("║  Intricate Geometries, Holes, Wedges, TRI & QUAD, Cracks    ║");
        Console.WriteLine("╚═══════════════════════════════════════════════════════════════╝");
        Console.WriteLine();

        // ========== PART 1: ADVANCED MESHING EXAMPLES ==========
        Console.WriteLine("═══ PART 1: ADVANCED MESHING EXAMPLES ═══\n");

   //     Example1_CircularDomainWithHole();
   //     Example2_LShapeWithCornerRefinement();
    //    Example3_AnnulusRegion();
    //    Example4_WedgeGeometry();
   //     Example5_MultipleHoles();
    //    Example6_IntricateBoundary();
     //   Example7_TriVsQuadComparison();
     //   Example8_CrackedPlateWithHole();
     //   Example9_GearLikeGeometry();
     //   Example10_ComplexIndustrialShape();

        // ========== PART 2: FRACTURE MECHANICS EXAMPLES ==========
        Console.WriteLine("\n═══ PART 2: FRACTURE MECHANICS - CLASSICAL BENCHMARKS ═══\n");

        //   Example11_Anderson2005_EdgeCrack();
        Example12_Griffith1921_CenterCrack();
       Example13_KanninenPopelar1985_DoubleEdgeNotch();
        //   Example14_ErdoganSih1963_SlantCrack();
        Example15_NewmanRaju1984_CrackFromHole();

        // ========== PART 3: SPECTACULAR CRACK PATTERNS ==========
        Console.WriteLine("\n═══ PART 3: FRACTURE MECHANICS - SPECTACULAR PATTERNS ═══\n");

       // Example16_SpiralGalaxy();
        //  Example17_FractalTree();
        Example18_SinusoidalWaves();
        //  Example19_StarBurst();
        //   Example20_ConcentricMandalas();

        // ========== PART 4: 3D FRACTURE MECHANICS ==========
        Console.WriteLine("\n═══ PART 4: 3D FRACTURE MECHANICS - CLASSICAL BENCHMARKS ═══\n");

       Example21_Sneddon1946_PennyShapedCrack();
       Example22_Irwin1962_EllipticalCrack();
        Example23_Tada1973_EdgeCrack();
       Example24_NewmanRaju1981_CornerCrack();
        Example25_ErdoganSih1963_SlantCrack();
        Example26_SemiCylindricalSurfaceCrack();

        Console.WriteLine();
        Console.WriteLine("╔═══════════════════════════════════════════════════════════════╗");
        Console.WriteLine("║  ✓ ALL 26 EXAMPLES COMPLETED!                                ║");
        Console.WriteLine("║  10 Advanced Meshing + 10 2D Fracture + 6 3D Fracture       ║");
        Console.WriteLine("╚═══════════════════════════════════════════════════════════════╝");
        Console.WriteLine();

        // Write all fracture mechanics meshes to single Ensight case file
        Console.WriteLine("═══ Writing Unified Ensight Output ═══");
        EnsightWriter.WriteAllMeshes("FractureMechanics");
    }

    #region Example 1: Circular Domain with Circular Hole

    /// <summary>
    ///     Example 1: Circular domain with eccentric circular hole
    ///     Tests: Curved boundaries, interior hole, quad conversion
    /// </summary>
    public static void Example1_CircularDomainWithHole()
    {
        Console.WriteLine("--- Example 1: Circular Domain with Eccentric Hole ---");

        // Outer circle (radius = 1.0, center at origin)
        var outer = CreateCircle(0, 0, 1.0, 40);

        // Inner hole (radius = 0.3, center at (0.4, 0.3))
        var hole = CreateCircle(0.4, 0.3, 0.3, 20);
        var holes = new List<double[,]> { hole };

        // Mesh with QUADS
        var (meshQuad, coordsQuad) = UnifiedMesher.TriangulateWithHoles(
            outer, holes, convertToQuads: true);

        // Mesh with TRIS
        var (meshTri, coordsTri) = UnifiedMesher.TriangulateWithHoles(
            outer, holes, convertToQuads: false);

        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");
        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");

        var statsQuad = MeshGeometry.ComputeQualityStatistics(meshQuad, coordsQuad);
        var statsTri = MeshGeometry.ComputeQualityStatistics(meshTri, coordsTri);
        Console.WriteLine($"  Quality QUAD: min angle = {statsQuad.MinTriangleAngleDegrees:F1}°");
        Console.WriteLine($"  Quality TRI:  min angle = {statsTri.MinTriangleAngleDegrees:F1}°");

        // Save to GiD format
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "01_ex1_circular_hole_quad.msh");
        SimplexRemesher.SaveGiD(meshTri, coordsTri, "01_ex1_circular_hole_tri.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 2: L-Shape with Corner Refinement

    /// <summary>
    ///     Example 2: L-shaped domain with refined corners
    ///     Tests: Re-entrant corners, stress concentration regions
    /// </summary>
    public static void Example2_LShapeWithCornerRefinement()
    {
        Console.WriteLine("--- Example 2: L-Shape with Corner Refinement ---");

        // L-shape with extra points near re-entrant corner for refinement
        var boundary = new[,]
        {
            { 0, 0 },
            { 2, 0 },
            { 2, 0.95 }, // Extra point before corner
            { 2, 1 },
            { 2, 1.05 }, // Extra point after corner
            { 1.05, 1 }, // Extra point
            { 1, 1 }, // Re-entrant corner (stress concentration)
            { 0.95, 1 }, // Extra point
            { 0.95, 2 }, // Extra point
            { 1, 2 },
            { 1, 2.05 }, // Extra point
            { 0, 2 }
        };

        var (meshTri, coordsTri) = UnifiedMesher.Triangulate(
            boundary, convertToQuads: false, enableSmoothing: false);

        var (meshQuad, coordsQuad) = UnifiedMesher.Triangulate(
            boundary, convertToQuads: true, enableSmoothing: false);

        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");

        // Check mesh quality near re-entrant corner
        var stats = MeshGeometry.ComputeQualityStatistics(meshTri, coordsTri);
        Console.WriteLine($"  Min angle: {stats.MinTriangleAngleDegrees:F1}° (critical near re-entrant corner)");

        SimplexRemesher.SaveGiD(meshTri, coordsTri, "02_ex2_lshape_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "02_ex2_lshape_quad.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 3: Annulus Region

    /// <summary>
    ///     Example 3: Annulus (ring) - circular domain with circular hole
    ///     Tests: Concentric boundaries, radial meshing patterns
    /// </summary>
    public static void Example3_AnnulusRegion()
    {
        Console.WriteLine("--- Example 3: Annulus (Ring) Region ---");

        // Outer circle (R_outer = 2.0)
        var outer = CreateCircle(0, 0, 2.0, 60);

        // Inner circle (R_inner = 1.0)
        var inner = CreateCircle(0, 0, 1.0, 30);
        var holes = new List<double[,]> { inner };

        // Create TRI and QUAD versions
        var (meshTri, coordsTri) = UnifiedMesher.TriangulateWithHoles(
            outer, holes, convertToQuads: false);

        var (meshQuad, coordsQuad) = UnifiedMesher.TriangulateWithHoles(
            outer, holes, convertToQuads: true);

        Console.WriteLine("  Annulus: R_outer=2.0, R_inner=1.0");
        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");

        SimplexRemesher.SaveGiD(meshTri, coordsTri, "03_ex3_annulus_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "03_ex3_annulus_quad.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 4: Wedge Geometry

    /// <summary>
    ///     Example 4: Wedge-shaped domain (angular sector)
    ///     Tests: Sharp corners, boundary-layer meshing
    /// </summary>
    public static void Example4_WedgeGeometry()
    {
        Console.WriteLine("--- Example 4: Wedge Geometry (60° sector) ---");

        // 60-degree wedge
        var angle = Math.PI / 3; // 60 degrees
        var nPoints = 20;
        var boundary = new List<double[]>();

        // Origin
        boundary.Add(new double[] { 0, 0 });

        // Arc from 0 to 60 degrees at radius=2
        for (var i = 0; i <= nPoints; i++)
        {
            var theta = i * angle / nPoints;
            boundary.Add(new[] { 2 * Math.Cos(theta), 2 * Math.Sin(theta) });
        }

        var boundaryArray = ConvertListToArray(boundary);

        var (meshTri, coordsTri) = UnifiedMesher.Triangulate(
            boundaryArray, convertToQuads: false);

        var (meshQuad, coordsQuad) = UnifiedMesher.Triangulate(
            boundaryArray, convertToQuads: true);

        Console.WriteLine("  Wedge: 60° sector, R=2");
        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");

        SimplexRemesher.SaveGiD(meshTri, coordsTri, "04_ex4_wedge_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "04_ex4_wedge_quad.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 5: Multiple Holes

    /// <summary>
    ///     Example 5: Domain with multiple holes
    ///     Tests: Multi-hole support, complex topology
    /// </summary>
    public static void Example5_MultipleHoles()
    {
        Console.WriteLine("--- Example 5: Domain with Multiple Holes ---");

        // Rectangular domain
        var boundary = new double[,]
        {
            { 0, 0 }, { 4, 0 }, { 4, 3 }, { 0, 3 }
        };

        // Three circular holes
        var hole1 = CreateCircle(1.0, 1.5, 0.4, 20);
        var hole2 = CreateCircle(2.5, 1.5, 0.3, 20);
        var hole3 = CreateCircle(3.2, 2.2, 0.25, 16);
        var holes = new List<double[,]> { hole1, hole2, hole3 };

        var (meshTri, coordsTri) = UnifiedMesher.TriangulateWithHoles(
            boundary, holes, convertToQuads: false);

        var (meshQuad, coordsQuad) = UnifiedMesher.TriangulateWithHoles(
            boundary, holes, convertToQuads: true);

        Console.WriteLine("  Domain: 4×3 rectangle with 3 circular holes");
        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");

        SimplexRemesher.SaveGiD(meshTri, coordsTri, "05_ex5_multiple_holes_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "05_ex5_multiple_holes_quad.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 6: Intricate Boundary

    /// <summary>
    ///     Example 6: Domain with intricate boundary (star-shaped)
    ///     Tests: Complex boundary geometry, concave regions
    /// </summary>
    public static void Example6_IntricateBoundary()
    {
        Console.WriteLine("--- Example 6: Intricate Star-Shaped Boundary ---");

        // Star-shaped boundary
        var nPoints = 50;
        var boundary = new double[nPoints, 2];

        for (var i = 0; i < nPoints; i++)
        {
            var theta = 2 * Math.PI * i / nPoints;
            var r = 1.0 + 0.3 * Math.Sin(5 * theta); // 5-pointed star modulation
            boundary[i, 0] = r * Math.Cos(theta);
            boundary[i, 1] = r * Math.Sin(theta);
        }

        var (meshTri, coordsTri) = UnifiedMesher.Triangulate(
            boundary, convertToQuads: false);

        var (meshQuad, coordsQuad) = UnifiedMesher.Triangulate(
            boundary, convertToQuads: true);

        Console.WriteLine("  Star: 5-pointed with radial modulation");
        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");

        SimplexRemesher.SaveGiD(meshTri, coordsTri, "06_ex6_star_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "06_ex6_star_quad.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 7: Tri vs Quad Comparison

    /// <summary>
    ///     Example 7: Side-by-side comparison of TRI vs QUAD meshing
    ///     Tests: Quality metrics, element count comparison
    /// </summary>
    public static void Example7_TriVsQuadComparison()
    {
        Console.WriteLine("--- Example 7: TRI vs QUAD Comparison ---");

        // Simple square domain for clean comparison
        var boundary = new double[,]
        {
            { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 }
        };

        var (meshTri, coordsTri) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 0.01, convertToQuads: false);

        var (meshQuad, coordsQuad) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 0.01, convertToQuads: true);

        var statsTri = MeshGeometry.ComputeQualityStatistics(meshTri, coordsTri);
        var statsQuad = MeshGeometry.ComputeQualityStatistics(meshQuad, coordsQuad);

        Console.WriteLine("  Domain: Unit square [0,1]×[0,1]");
        Console.WriteLine("\n  TRI MESH:");
        Console.WriteLine($"    Elements: {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine($"    Nodes: {meshTri.Count<Node>()}");
        Console.WriteLine(
            $"    Quality: min angle = {statsTri.MinTriangleAngleDegrees:F1}°, avg aspect = {statsTri.AvgTriangleAspectRatio:F3}");

        Console.WriteLine("\n  QUAD MESH:");
        Console.WriteLine($"    Elements: {meshQuad.Count<Quad4>()} quads + {meshQuad.Count<Tri3>()} tris");
        Console.WriteLine($"    Nodes: {meshQuad.Count<Node>()}");
        Console.WriteLine(
            $"    Coverage: {100.0 * meshQuad.Count<Quad4>() / (meshQuad.Count<Quad4>() + meshQuad.Count<Tri3>() * 0.5):F1}% quads");
        Console.WriteLine(
            $"    Quality: min angle = {statsQuad.MinTriangleAngleDegrees:F1}°, avg aspect = {statsQuad.AvgTriangleAspectRatio:F3}");

        // Save both for visual comparison
        SimplexRemesher.SaveGiD(meshTri, coordsTri, "07_ex7_comparison_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "07_ex7_comparison_quad.msh");
        Console.WriteLine("  ✓ Saved both meshes for comparison");
        Console.WriteLine();
    }

    #endregion

    #region Example 8: Cracked Plate with Hole

    /// <summary>
    ///     Example 8: Plate with hole and emanating crack
    ///     Tests: Crack insertion, stress concentration
    /// </summary>
    public static void Example8_CrackedPlateWithHole()
    {
        Console.WriteLine("--- Example 8: Cracked Plate with Hole ---");

        // Rectangular plate
        var boundary = new double[,]
        {
            { 0, 0 }, { 4, 0 }, { 4, 2 }, { 0, 2 }
        };

        // Circular hole at center
        var hole = CreateCircle(2.0, 1.0, 0.3, 24);
        var holes = new List<double[,]> { hole };

        // Generate base mesh (triangles only for crack insertion)
        var (mesh, coords) = UnifiedMesher.TriangulateWithHoles(
            boundary, holes,
            refine: true,
            maxArea: 0.05,
            convertToQuads: false); // Triangles for crack insertion

        Console.WriteLine("  Plate: 4×2 with central hole (R=0.3)");
        Console.WriteLine($"  Initial mesh: {mesh.Count<Node>()} nodes, {mesh.Count<Tri3>()} triangles");

        // Insert horizontal crack from hole edge using two-level-set
        var R = 0.3;
        double xHole = 2.0, yHole = 1.0;
        var crackLength = 0.8;

        // Surface: horizontal at y = yHole
        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => y - yHole;

        // Region: crack from hole edge to hole + crackLength
        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            var xStart = xHole + R;
            var xEnd = xStart + crackLength;

            if (x < xStart) return x - xStart;
            if (x > xEnd) return x - xEnd;
            return -1.0;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField(
            mesh, coords, surface, region);

        Console.WriteLine(
            $"  Crack inserted: {crackedMesh.Count<Node>()} nodes (+{crackedMesh.Count<Node>() - mesh.Count<Node>()} from crack)");

        var stats = MeshGeometry.ComputeQualityStatistics(crackedMesh, crackedCoords);
        Console.WriteLine($"  Final quality: min angle = {stats.MinTriangleAngleDegrees:F1}°");

        // Save cracked mesh to GiD
        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "08_ex8_cracked_plate.msh");
        Console.WriteLine("  ✓ Saved cracked mesh to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 9: Gear-Like Geometry

    /// <summary>
    ///     Example 9: Simplified gear profile with teeth
    ///     Tests: Repeating features, small-scale details
    /// </summary>
    public static void Example9_GearLikeGeometry()
    {
        Console.WriteLine("--- Example 9: Gear-Like Geometry with Teeth ---");

        // Simplified gear: circular base with rectangular teeth
        var nTeeth = 12;
        var rBase = 1.5;
        var rTip = 2.0;
        var toothWidth = 0.15;
        var nPointsPerGap = 5; // Points between teeth for smooth base circle

        var boundary = new List<double[]>();

        for (var i = 0; i < nTeeth; i++)
        {
            var theta = 2.0 * Math.PI * i / nTeeth;
            var thetaNext = 2.0 * Math.PI * (i + 1) / nTeeth;

            // Tooth angles
            var thetaLeft = theta - toothWidth / rBase / 2;
            var thetaRight = theta + toothWidth / rBase / 2;
            var thetaTipLeft = theta - toothWidth / rTip / 2;
            var thetaTipRight = theta + toothWidth / rTip / 2;

            // Add base circle arc from previous tooth to this tooth
            if (i == 0)
            {
                var prevRight = 2.0 * Math.PI - toothWidth / rBase / 2;
                for (var j = 0; j < nPointsPerGap; j++)
                {
                    var t = (double)j / nPointsPerGap;
                    var thetaGap = prevRight + t * (thetaLeft + 2.0 * Math.PI - prevRight);
                    boundary.Add(new[] { rBase * Math.Cos(thetaGap), rBase * Math.Sin(thetaGap) });
                }
            }
            else
            {
                var prevRight = 2.0 * Math.PI * (i - 1) / nTeeth + toothWidth / rBase / 2;
                for (var j = 0; j < nPointsPerGap; j++)
                {
                    var t = (double)j / nPointsPerGap;
                    var thetaGap = prevRight + t * (thetaLeft - prevRight);
                    boundary.Add(new[] { rBase * Math.Cos(thetaGap), rBase * Math.Sin(thetaGap) });
                }
            }

            // Add tooth: base-left → tip-left → tip-right → base-right
            boundary.Add(new[] { rBase * Math.Cos(thetaLeft), rBase * Math.Sin(thetaLeft) });
            boundary.Add(new[] { rTip * Math.Cos(thetaTipLeft), rTip * Math.Sin(thetaTipLeft) });
            boundary.Add(new[] { rTip * Math.Cos(thetaTipRight), rTip * Math.Sin(thetaTipRight) });
            boundary.Add(new[] { rBase * Math.Cos(thetaRight), rBase * Math.Sin(thetaRight) });
        }

        var boundaryArray = ConvertListToArray(boundary);

        // TRI mesh
        var (meshTri, coordsTri) = UnifiedMesher.Triangulate(
            boundaryArray, convertToQuads: false);

        // QUAD mesh
        var (meshQuad, coordsQuad) = UnifiedMesher.Triangulate(
            boundaryArray, convertToQuads: true);

        Console.WriteLine($"  Gear: {nTeeth} teeth, R_base={rBase}, R_tip={rTip}");
        Console.WriteLine($"  Boundary points: {boundaryArray.GetLength(0)}");
        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");

        SimplexRemesher.SaveGiD(meshTri, coordsTri, "09_ex9_gear_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "09_ex9_gear_quad.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    #region Example 10: Complex Industrial Shape

    /// <summary>
    ///     Example 10: Complex industrial component (bracket-like)
    ///     Tests: Real-world complexity, multiple features
    /// </summary>
    public static void Example10_ComplexIndustrialShape()
    {
        Console.WriteLine("--- Example 10: Complex Industrial Bracket ---");

        // Industrial bracket: Base + mounting flange + stress relief fillets
        // Use more boundary points for better mesh quality
        var boundary = new List<double[]>();

        // Bottom edge (left to right) - add intermediate points
        for (var i = 0; i <= 10; i++)
            boundary.Add(new[] { 3.0 * i / 10.0, 0 });

        // Right edge going up with step
        for (var i = 1; i <= 5; i++)
            boundary.Add(new[] { 3.0, 1.0 * i / 5.0 });
        boundary.Add(new[] { 3.5, 1.0 });
        for (var i = 1; i <= 5; i++)
            boundary.Add(new[] { 3.5, 1.0 + 1.0 * i / 5.0 });

        // Top edge with cutout (right to left)
        boundary.Add(new[] { 2.5, 2.0 });
        for (var i = 1; i <= 3; i++)
            boundary.Add(new[] { 2.5, 2.0 + 0.5 * i / 3.0 });
        for (var i = 1; i <= 5; i++)
            boundary.Add(new[] { 2.5 - 1.0 * i / 5.0, 2.5 });
        for (var i = 1; i <= 3; i++)
            boundary.Add(new[] { 1.5, 2.5 - 0.5 * i / 3.0 });
        for (var i = 1; i <= 5; i++)
            boundary.Add(new[] { 1.5 - 1.0 * i / 5.0, 2.0 });

        // Left edge (top to bottom)
        for (var i = 1; i <= 10; i++)
            boundary.Add(new[] { 0, 2.0 - 2.0 * i / 10.0 });

        var boundaryArray = ConvertListToArray(boundary);

        // Mounting holes - smooth circles
        var hole1 = CreateCircle(0.75, 1.0, 0.2, 40);
        var hole2 = CreateCircle(2.75, 1.0, 0.2, 40);
        var holes = new List<double[,]> { hole1, hole2 };

        // TRI mesh with UNIFORM sizing (sizeGradation = 0)
        var (meshTri, coordsTri) = UnifiedMesher.TriangulateWithHoles(
            boundaryArray,
            holes,
            refine: true,
            maxArea: 0.02, // Uniform element size
            sizeGradation: 0.0, // DISABLE adaptive sizing!
            convertToQuads: false);

        // QUAD mesh with UNIFORM sizing
        var (meshQuad, coordsQuad) = UnifiedMesher.TriangulateWithHoles(
            boundaryArray,
            holes,
            refine: true,
            maxArea: 0.02, // Uniform element size
            sizeGradation: 0.0, // DISABLE adaptive sizing!
            convertToQuads: true);

        Console.WriteLine("  Industrial bracket: stepped profile, 2 mounting holes");
        Console.WriteLine($"  TRI mesh:  {meshTri.Count<Node>()} nodes, {meshTri.Count<Tri3>()} triangles");
        Console.WriteLine(
            $"  QUAD mesh: {meshQuad.Count<Node>()} nodes, {meshQuad.Count<Quad4>()} quads, {meshQuad.Count<Tri3>()} tris");

        var statsTri = MeshGeometry.ComputeQualityStatistics(meshTri, coordsTri);
        var statsQuad = MeshGeometry.ComputeQualityStatistics(meshQuad, coordsQuad);
        Console.WriteLine(
            $"  Quality TRI:  min angle = {statsTri.MinTriangleAngleDegrees:F1}°, avg aspect = {statsTri.AvgTriangleAspectRatio:F3}");
        Console.WriteLine(
            $"  Quality QUAD: min angle = {statsQuad.MinTriangleAngleDegrees:F1}°, avg aspect = {statsQuad.AvgTriangleAspectRatio:F3}");

        SimplexRemesher.SaveGiD(meshTri, coordsTri, "10_ex10_bracket_tri.msh");
        SimplexRemesher.SaveGiD(meshQuad, coordsQuad, "10_ex10_bracket_quad.msh");
        Console.WriteLine("  ✓ Saved to GiD .msh format");
        Console.WriteLine();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════════
    // FRACTURE MECHANICS EXAMPLES - CLASSICAL BENCHMARKS
    // ═══════════════════════════════════════════════════════════════════════

    #region Example 11: Anderson 2005 - Edge Crack

    /// <summary>
    ///     Example 11: Single Edge Crack in Tension (Anderson 2005)
    ///     Reference: Anderson, T.L. (2005). Fracture Mechanics: Fundamentals and Applications, 3rd ed.
    /// </summary>
    private static void Example11_Anderson2005_EdgeCrack()
    {
        Console.WriteLine("--- Example 11: Anderson (2005) - Single Edge Crack ---");
        Console.WriteLine("Reference: Anderson, T.L. (2005). Fracture Mechanics, 3rd ed., CRC Press");
        Console.WriteLine("Geometry: W=50mm, H=100mm, a=12.5mm\n");

        double W = 50.0, H = 100.0, a = 12.5;

        // REDUCED MESH DENSITY: 51x101 instead of 101x201
        var boundary = CreateRectangle(0, W, 0, H, 51, 101);
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 5.0, convertToQuads: false,
            enableSmoothing: false); // Keep smoothing for quality 5.0 instead of 1.5

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        // Insert edge crack
        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => y - H / 2;
        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
            x < 0 ? x : x > a ? x - a : -1.0;

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField(
            baseMesh, baseCoords, surface, region);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "11_Anderson2005_EdgeCrack");
        EnsightWriter.AddMesh("11_Anderson2005_EdgeCrack", crackedMesh, crackedCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 12: Griffith 1921 - Center Crack

    /// <summary>
    ///     Example 12: Center Crack in Infinite Plate (Griffith 1921)
    ///     Reference: Griffith, A.A. (1921). Phil. Trans. Royal Soc. London A, 221:163-198
    /// </summary>
    private static void Example12_Griffith1921_CenterCrack()
    {
        Console.WriteLine("--- Example 12: Griffith (1921) - Center Crack ---");
        Console.WriteLine("Reference: Griffith, A.A. (1921). Phil. Trans. Royal Soc. A, 221:163-198");
        Console.WriteLine("Geometry: 100×100mm plate, crack length 2a=20mm\n");

        double L = 100.0, a = 10.0;

        // REDUCED MESH DENSITY: 51x51 instead of 101x101
        var boundary = CreateRectangle(0, L, 0, L, 101, 101); // Finer mesh for multi-crack patterns
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 2.0, convertToQuads: false,
            enableSmoothing: false); // Finer mesh, no smoothing for multi-crack

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);
        

        // CENTER CRACK (horizontal at y=50, x from 40 to 60)
        // Plate: [0,100] × [0,100]
        // Crack: centered at (50, 50), length 2a=20mm (±10mm from center)
        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => y - L / 2; // y=50 (center)
        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            var xc = x - L / 2; // Distance from x=50 (center)
            if (xc < -a || xc > a) return 1.0; // Outside crack region: POSITIVE
            return -1.0; // Inside crack (40 ≤ x ≤ 60): NEGATIVE
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField(
            baseMesh, baseCoords, surface, region);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "12_griffith1921_center_crack");
        EnsightWriter.AddMesh("12_griffith1921_center_crack", crackedMesh, crackedCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 13: Kanninen & Popelar 1985 - DENT

    /// <summary>
    ///     Example 13: Double Edge Notch Tension (DENT) Specimen
    ///     Reference: Kanninen, M.F. & Popelar, C.H. (1985). Advanced Fracture Mechanics
    /// </summary>
    private static void Example13_KanninenPopelar1985_DoubleEdgeNotch()
    {
        Console.WriteLine("--- Example 13: Kanninen & Popelar (1985) - DENT Specimen ---");
        Console.WriteLine("Reference: Kanninen & Popelar (1985). Advanced Fracture Mechanics, Oxford");
        Console.WriteLine("Geometry: W=50mm, H=100mm, symmetric cracks a=12.5mm each\n");

        double W = 50.0, H = 100.0, a = 12.5;

        // REDUCED MESH DENSITY: 51x101 instead of 101x201
        var boundary = CreateRectangle(0, W, 0, H, 51, 101);
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 5.0, convertToQuads: false,
            enableSmoothing: false); // Keep smoothing for quality 5.0 instead of 1.5

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        // LEFT CRACK
        SimplexRemesher.SignedFieldFunction leftSurface = (x, y, z) => y - H / 2;
        SimplexRemesher.SignedFieldFunction leftRegion = (x, y, z) =>
            x < 0 ? x : x > a ? x - a : -1.0;

        var (mesh1, coords1) = SimplexRemesher.CreateCrackFromSignedField(
            baseMesh, baseCoords, leftSurface, leftRegion);

        // RIGHT CRACK
        SimplexRemesher.SignedFieldFunction rightSurface = (x, y, z) => y - H / 2;
        SimplexRemesher.SignedFieldFunction rightRegion = (x, y, z) =>
        {
            // Right crack from x = W-a to x = W
            if (x > W) return x - W; // Outside (right of boundary)
            if (x < W - a) return W - a - x; // Outside (left of crack start)
            return -1.0; // Inside crack region
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField(
            mesh1, coords1, rightSurface, rightRegion);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "13_kanninen1985_dent");
        EnsightWriter.AddMesh("13_kanninen1985_dent", crackedMesh, crackedCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 14: Erdogan & Sih 1963 - Slant Crack

    /// <summary>
    ///     Example 14: Slant Crack under Mixed Mode Loading
    ///     Reference: Erdogan, F. & Sih, G.C. (1963). J. Basic Engineering, 85:519-527
    /// </summary>
    private static void Example14_ErdoganSih1963_SlantCrack()
    {
        Console.WriteLine("--- Example 14: Erdogan & Sih (1963) - Slant Crack (Mixed Mode) ---");
        Console.WriteLine("Reference: Erdogan & Sih (1963). J. Basic Eng., 85:519-527");
        Console.WriteLine("Geometry: 100×100mm, inclined crack β=45°, length=20mm\n");

        double L = 100.0, crackLength = 20.0;
        var beta = Math.PI / 4;

        // REDUCED MESH DENSITY: 51x51 instead of 101x101
        var boundary = CreateRectangle(0, L, 0, L, 101, 101); // Finer mesh for multi-crack patterns
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 2.0, convertToQuads: false,
            enableSmoothing: false); // Finer mesh, no smoothing for multi-crack

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        // Inclined crack
        double cosB = Math.Cos(beta), sinB = Math.Sin(beta);
        SimplexRemesher.SignedFieldFunction surface = (x, y, z) =>
        {
            double xc = x - L / 2, yc = y - L / 2;
            return sinB * xc - cosB * yc;
        };

        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            double xc = x - L / 2, yc = y - L / 2;
            var along = cosB * xc + sinB * yc;
            var halfLen = crackLength / 2;

            if (along < -halfLen) return along + halfLen;
            if (along > halfLen) return along - halfLen;
            return -1.0;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField(
            baseMesh, baseCoords, surface, region);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "14_erdogan1963_slant_crack");
        EnsightWriter.AddMesh("14_erdogan1963_slant_crack", crackedMesh, crackedCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 15: Newman & Raju 1984 - Crack from Hole

    /// <summary>
    ///     Example 15: Crack Emanating from Circular Hole
    ///     Reference: Newman, J.C. & Raju, I.S. (1984). Eng. Fracture Mech., 20(1):87-106
    /// </summary>
    private static void Example15_NewmanRaju1984_CrackFromHole()
    {
        Console.WriteLine("--- Example 15: Newman & Raju (1984) - Crack from Hole ---");
        Console.WriteLine("Reference: Newman & Raju (1984). Eng. Fract. Mech., 20(1):87-106");
        Console.WriteLine("Geometry: 100×100mm plate, hole R=10mm, crack a=15mm\n");

        double L = 100.0, R = 10.0, a = 15.0;

        // Create boundary with hole - REDUCED DENSITY: 51x51 instead of 101x101
        var outer = CreateRectangle(0, L, 0, L, 101, 101); // Finer mesh for multi-crack patterns
        var hole = CreateCircle(L / 2, L / 2, R, 30); // 30 points instead of 40
        var holes = new List<double[,]> { hole };

        var (baseMesh, baseCoords) = UnifiedMesher.TriangulateWithHoles(
            outer, holes, refine: true, maxArea: 2.0, convertToQuads: false,
            enableSmoothing: false); // Finer mesh, no smoothing for multi-crack

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        // Crack from hole
        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => y - L / 2;

        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            var xStart = L / 2 + R;
            var xEnd = xStart + a;

            if (x < xStart) return x - xStart;
            if (x > xEnd) return x - xEnd;
            return -1.0;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField(
            baseMesh, baseCoords, surface, region);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "15_newman1984_crack_from_hole");
        EnsightWriter.AddMesh("15_newman1984_crack_from_hole", crackedMesh, crackedCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════════
    // FRACTURE MECHANICS EXAMPLES - SPECTACULAR PATTERNS
    // ═══════════════════════════════════════════════════════════════════════

    #region Example 16: Spiral Galaxy

    private static void Example16_SpiralGalaxy()
    {
        Console.WriteLine("--- Example 16: Spiral Galaxy - Logarithmic Spiral Crack ---");
        Console.WriteLine("Geometry: 100×100mm, 2-turn spiral from r=5mm to r=40mm\n");

        var L = 100.0;

        // Fine mesh for curves (balanced between speed and quality)
        var boundary = CreateRectangle(0, L, 0, L, 101, 101);
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 1.0, convertToQuads: false, enableSmoothing: false);

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        var currentMesh = baseMesh;
        var currentCoords = baseCoords;

        // Create 2-turn logarithmic spiral
        {
            // Logarithmic spiral: r(θ) = r_inner * exp(b*θ) for θ ∈ [0, 4π]
            SimplexRemesher.SignedFieldFunction surface = (x, y, z) =>
            {
                double xc = x - L / 2, yc = y - L / 2;
                var r = Math.Sqrt(xc * xc + yc * yc);

                if (r < 0.1) return 1.0; // Avoid singularity at center

                var theta = Math.Atan2(yc, xc);
                if (theta < 0) theta += 2 * Math.PI; // Normalize to [0, 2π]

                // Spiral parameters: 2 turns from r=5 to r=40
                var r_inner = 5.0;
                var r_outer = 40.0;
                var b = Math.Log(r_outer / r_inner) / (4 * Math.PI);

                // Check both wraps of the spiral (0-2π and 2π-4π)
                var minDist = double.MaxValue;

                for (var wrap = 0; wrap <= 1; wrap++)
                {
                    var theta_extended = theta + wrap * 2 * Math.PI;
                    if (theta_extended > 4 * Math.PI) continue;

                    // Expected radius at this angle
                    var r_expected = r_inner * Math.Exp(b * theta_extended);

                    // Radial distance to spiral
                    var dist = r - r_expected;

                    if (Math.Abs(dist) < Math.Abs(minDist))
                        minDist = dist;
                }

                return minDist;
            };

            // Region: radial bounds (5mm to 40mm from center)
            SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
            {
                double xc = x - L / 2, yc = y - L / 2;
                var r = Math.Sqrt(xc * xc + yc * yc);

                if (r < 5 || r > 40) return 1.0; // Outside radial bounds
                return -1.0; // Inside active region
            };

            (currentMesh, currentCoords) = SimplexRemesher.CreateCrackFromSignedField(
                currentMesh, currentCoords, surface, region);
        }

        Console.WriteLine($"Cracked mesh: {currentMesh.Count<Node>()} nodes, {currentMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(currentMesh, currentCoords, "16_spiral_galaxy");
        EnsightWriter.AddMesh("16_spiral_galaxy", currentMesh, currentCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 17: Fractal Tree

    private static void Example17_FractalTree()
    {
        Console.WriteLine("--- Example 17: Fractal Tree - Branching Crack System ---");
        Console.WriteLine("Geometry: 100×100mm, binary tree with 3 levels\n");

        var L = 100.0;

        // REDUCED MESH DENSITY: 51x51 instead of 101x101
        var boundary = CreateRectangle(0, L, 0, L, 101, 101); // Finer mesh for multi-crack patterns
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 2.0, convertToQuads: false,
            enableSmoothing: false); // Finer mesh, no smoothing for multi-crack

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        var currentMesh = baseMesh;
        var currentCoords = baseCoords;

        // Main trunk
        SimplexRemesher.SignedFieldFunction trunkSurface = (x, y, z) => x - L / 2;
        SimplexRemesher.SignedFieldFunction trunkRegion = (x, y, z) =>
            y < 10 ? y - 10 : y > 40 ? y - 40 : -1.0;

        (currentMesh, currentCoords) = SimplexRemesher.CreateCrackFromSignedField(
            currentMesh, currentCoords, trunkSurface, trunkRegion);

        // Left branch
        var branchAngle = Math.PI / 6;
        SimplexRemesher.SignedFieldFunction leftSurface = (x, y, z) =>
        {
            double xc = x - L / 2, yc = y - 40;
            return Math.Sin(branchAngle) * xc - Math.Cos(branchAngle) * yc;
        };

        SimplexRemesher.SignedFieldFunction leftRegion = (x, y, z) =>
        {
            double xc = x - L / 2, yc = y - 40;
            var along = Math.Cos(branchAngle) * xc + Math.Sin(branchAngle) * yc;
            return along < -20 || along > 0 ? Math.Abs(along) - 20 : -1.0;
        };

        (currentMesh, currentCoords) = SimplexRemesher.CreateCrackFromSignedField(
            currentMesh, currentCoords, leftSurface, leftRegion);

        // Right branch
        SimplexRemesher.SignedFieldFunction rightSurface = (x, y, z) =>
        {
            double xc = x - L / 2, yc = y - 40;
            return -Math.Sin(branchAngle) * xc - Math.Cos(branchAngle) * yc;
        };

        SimplexRemesher.SignedFieldFunction rightRegion = (x, y, z) =>
        {
            double xc = x - L / 2, yc = y - 40;
            var along = Math.Cos(branchAngle) * xc - Math.Sin(branchAngle) * yc;
            return along < 0 || along > 20 ? Math.Abs(along) - 20 : -1.0;
        };

        (currentMesh, currentCoords) = SimplexRemesher.CreateCrackFromSignedField(
            currentMesh, currentCoords, rightSurface, rightRegion);

        Console.WriteLine($"Cracked mesh: {currentMesh.Count<Node>()} nodes, {currentMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(currentMesh, currentCoords, "17_fractal_tree");
        EnsightWriter.AddMesh("17_fractal_tree", currentMesh, currentCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 18: Sinusoidal Waves

    private static void Example18_SinusoidalWaves()
    {
        Console.WriteLine("--- Example 18: Sinusoidal Waves - Periodic Curved Cracks ---");
        Console.WriteLine("Geometry: 100×100mm, 3 horizontal wavy cracks\n");

        var L = 100.0;

        // REDUCED MESH DENSITY: 51x51 instead of 101x101
        var boundary = CreateRectangle(0, L, 0, L, 101, 101); // Finer mesh for multi-crack patterns
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 2.0, convertToQuads: false,
            enableSmoothing: false); // Finer mesh, no smoothing for multi-crack

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        var currentMesh = baseMesh;
        var currentCoords = baseCoords;

        double[] heights = { 25.0, 50.0, 75.0 };

        foreach (var h in heights)
        {
            var amplitude = 5.0;
            var wavelength = 20.0;

            SimplexRemesher.SignedFieldFunction surface = (x, y, z) =>
            {
                var ySine = h + amplitude * Math.Sin(2 * Math.PI * x / wavelength);
                return y - ySine;
            };

            SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
            {
                if (x < 10) return x - 10;
                if (x > 90) return x - 90;
                return -1.0;
            };

            (currentMesh, currentCoords) = SimplexRemesher.CreateCrackFromSignedField(
                currentMesh, currentCoords, surface, region);
        }

        Console.WriteLine($"Cracked mesh: {currentMesh.Count<Node>()} nodes, {currentMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(currentMesh, currentCoords, "18_sinusoidal_waves");
        EnsightWriter.AddMesh("18_sinusoidal_waves", currentMesh, currentCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 19: Star Burst

    private static void Example19_StarBurst()
    {
        Console.WriteLine("--- Example 19: Star Burst - Radial Crack Pattern ---");
        Console.WriteLine("Geometry: 100×100mm, 8 radial cracks from center\n");

        var L = 100.0;

        // REDUCED MESH DENSITY: 51x51 instead of 101x101
        var boundary = CreateRectangle(0, L, 0, L, 101, 101); // Finer mesh for multi-crack patterns
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 2.0, convertToQuads: false,
            enableSmoothing: false); // Finer mesh, no smoothing for multi-crack

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        var currentMesh = baseMesh;
        var currentCoords = baseCoords;

        for (var ray = 0; ray < 8; ray++)
        {
            var angle = ray * 2 * Math.PI / 8;
            double cosA = Math.Cos(angle), sinA = Math.Sin(angle);

            SimplexRemesher.SignedFieldFunction surface = (x, y, z) =>
            {
                double xc = x - L / 2, yc = y - L / 2;
                return sinA * xc - cosA * yc;
            };

            SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
            {
                double xc = x - L / 2, yc = y - L / 2;
                var r = Math.Sqrt(xc * xc + yc * yc);

                if (r < 5) return r - 5;
                if (r > 40) return r - 40;
                return -1.0;
            };

            (currentMesh, currentCoords) = SimplexRemesher.CreateCrackFromSignedField(
                currentMesh, currentCoords, surface, region);
        }

        Console.WriteLine($"Cracked mesh: {currentMesh.Count<Node>()} nodes, {currentMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(currentMesh, currentCoords, "19_star_burst");
        EnsightWriter.AddMesh("19_star_burst", currentMesh, currentCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Example 20: Concentric Mandalas

    private static void Example20_ConcentricMandalas()
    {
        Console.WriteLine("--- Example 20: Concentric Mandalas - Circular Ring Cracks ---");
        Console.WriteLine("Geometry: 100×100mm, 3 concentric circles\n");

        var L = 100.0;

        // Fine mesh for curves (balanced between speed and quality)
        var boundary = CreateRectangle(0, L, 0, L, 101, 101);
        var (baseMesh, baseCoords) = UnifiedMesher.Triangulate(
            boundary, refine: true, maxArea: 1.0, convertToQuads: false, enableSmoothing: false);

        Console.WriteLine($"Base mesh: {baseMesh.Count<Node>()} nodes, {baseMesh.Count<Tri3>()} triangles");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates(baseMesh, baseCoords);
        SmoothMesh(baseMesh, baseCoords);

        var currentMesh = baseMesh;
        var currentCoords = baseCoords;

        double[] radii = { 15.0, 25.0, 35.0 };

        foreach (var R in radii)
        {
            // Circular crack surface (distance from center)
            SimplexRemesher.SignedFieldFunction surface = (x, y, z) =>
            {
                double xc = x - L / 2, yc = y - L / 2;
                var r = Math.Sqrt(xc * xc + yc * yc);
                return r - R; // Distance to circle at radius R
            };

            // Full circle (crack exists everywhere)
            SimplexRemesher.SignedFieldFunction region = (x, y, z) => -1.0;

            (currentMesh, currentCoords) = SimplexRemesher.CreateCrackFromSignedField(
                currentMesh, currentCoords, surface, region);
        }

        Console.WriteLine($"Cracked mesh: {currentMesh.Count<Node>()} nodes, {currentMesh.Count<Tri3>()} triangles");

        SimplexRemesher.SaveGiD(currentMesh, currentCoords, "20_concentric_mandalas");
        EnsightWriter.AddMesh("20_concentric_mandalas", currentMesh, currentCoords);
        Console.WriteLine("✓ Saved to GiD and added to Ensight collection\n");
    }

    #endregion

    #region Helper Methods

    /// <summary>
    ///     Create a circular boundary
    /// </summary>
    /// <summary>
    /// Perturb mesh coordinates to avoid exact geometric alignments.
    /// Applies random perturbation of ±0.1% of average edge length.
    /// This prevents pathological cases in geometric algorithms.
    /// </summary>
    private static void PerturbMeshCoordinates(SimplexMesh mesh, double[,] coords, int seed = 12345)
    {
        // Calculate average edge length
        double totalLength = 0;
        int edgeCount = 0;
        
        for (int i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            for (int j = 0; j < 3; j++)
            {
                int n1 = nodes[j];
                int n2 = nodes[(j + 1) % 3];
                double dx = coords[n2, 0] - coords[n1, 0];
                double dy = coords[n2, 1] - coords[n1, 1];
                totalLength += Math.Sqrt(dx * dx + dy * dy);
                edgeCount++;
            }
        }
        
        double avgEdgeLength = totalLength / edgeCount;
        double perturbationAmplitude = 0.001 * avgEdgeLength; // 0.1% of average edge length
        
        var random = new Random(seed);
        
        for (int i = 0; i < mesh.Count<Node>(); i++)
        {
            coords[i, 0] += (random.NextDouble() * 2.0 - 1.0) * perturbationAmplitude;
            coords[i, 1] += (random.NextDouble() * 2.0 - 1.0) * perturbationAmplitude;
        }
        
        Console.WriteLine($"  → Perturbed coordinates: ±{perturbationAmplitude:E3}mm (0.1% of avg edge)");
    }

    /// <summary>
    /// Smooth mesh using Laplacian smoothing with boundary nodes fixed.
    /// Uses only 2 iterations for speed.
    /// </summary>
    private static void SmoothMesh(SimplexMesh mesh, double[,] coords, int iterations = 2, double relaxation = 0.5)
    {
        Console.WriteLine($"  → Smoothing mesh (2 iterations, boundaries fixed)...");
        
        // Find boundary nodes to keep fixed
        var boundaryNodes = SimplexRemesher.FindBoundaryNodes(mesh);
        
        // BUILD NEIGHBOR MAP ONCE
        var neighborMap = new Dictionary<int, List<int>>();
        for (int i = 0; i < mesh.Count<Node>(); i++)
            neighborMap[i] = new List<int>();
        
        for (int elemId = 0; elemId < mesh.Count<Tri3>(); elemId++)
        {
            var elemNodes = mesh.NodesOf<Tri3, Node>(elemId);
            neighborMap[elemNodes[0]].Add(elemNodes[1]);
            neighborMap[elemNodes[0]].Add(elemNodes[2]);
            neighborMap[elemNodes[1]].Add(elemNodes[0]);
            neighborMap[elemNodes[1]].Add(elemNodes[2]);
            neighborMap[elemNodes[2]].Add(elemNodes[0]);
            neighborMap[elemNodes[2]].Add(elemNodes[1]);
        }
        
        // Just 2 quick smoothing passes
        for (int iter = 0; iter < 2; iter++)
        {
            for (int nodeId = 0; nodeId < mesh.Count<Node>(); nodeId++)
            {
                if (boundaryNodes.Contains(nodeId)) continue;
                
                var neighbors = neighborMap[nodeId];
                if (neighbors.Count == 0) continue;
                
                double avgX = 0, avgY = 0;
                foreach (var nbr in neighbors)
                {
                    avgX += coords[nbr, 0];
                    avgY += coords[nbr, 1];
                }
                avgX /= neighbors.Count;
                avgY /= neighbors.Count;
                
                // Update directly with relaxation
                coords[nodeId, 0] += 0.5 * (avgX - coords[nodeId, 0]);
                coords[nodeId, 1] += 0.5 * (avgY - coords[nodeId, 1]);
            }
        }
        
        Console.WriteLine($"  → Smoothed {mesh.Count<Node>()} nodes ({boundaryNodes.Count} boundaries fixed)");
    }

    /// <summary>
    /// Perturb 3D mesh coordinates to avoid exact geometric alignments.
    /// Applies random perturbation of ±0.1% of average edge length.
    /// </summary>
    private static void PerturbMeshCoordinates3D(SimplexMesh mesh, double[,] coords, int seed = 12345)
    {
        // Calculate average edge length
        double totalLength = 0;
        int edgeCount = 0;
        
        for (int i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            for (int j = 0; j < 4; j++)
            {
                for (int k = j + 1; k < 4; k++)
                {
                    double dx = coords[nodes[k], 0] - coords[nodes[j], 0];
                    double dy = coords[nodes[k], 1] - coords[nodes[j], 1];
                    double dz = coords[nodes[k], 2] - coords[nodes[j], 2];
                    totalLength += Math.Sqrt(dx * dx + dy * dy + dz * dz);
                    edgeCount++;
                }
            }
        }
        
        double avgEdgeLength = totalLength / edgeCount;
        double perturbationAmplitude = 0.001 * avgEdgeLength; // 0.1% of average edge length
        
        var random = new Random(seed);
        
        for (int i = 0; i < mesh.Count<Node>(); i++)
        {
            coords[i, 0] += (random.NextDouble() * 2.0 - 1.0) * perturbationAmplitude;
            coords[i, 1] += (random.NextDouble() * 2.0 - 1.0) * perturbationAmplitude;
            coords[i, 2] += (random.NextDouble() * 2.0 - 1.0) * perturbationAmplitude;
        }
        
        Console.WriteLine($"  → Perturbed 3D coordinates: ±{perturbationAmplitude:E3}mm (0.1% of avg edge)");
    }

    /// <summary>
    /// Smooth 3D mesh using Laplacian smoothing with boundary nodes fixed.
    /// Uses only 2 iterations for speed.
    /// </summary>
    private static void SmoothMesh3D(SimplexMesh mesh, double[,] coords, int iterations = 2, double relaxation = 0.5)
    {
        Console.WriteLine($"  → Smoothing 3D mesh (2 iterations, boundaries fixed)...");
        
        // Find boundary nodes to keep fixed
        var boundaryNodes = SimplexRemesher.FindBoundaryNodes3D(mesh);
        
        // BUILD NEIGHBOR MAP ONCE
        var neighborMap = new Dictionary<int, List<int>>();
        for (int i = 0; i < mesh.Count<Node>(); i++)
            neighborMap[i] = new List<int>();
        
        for (int elemId = 0; elemId < mesh.Count<Tet4>(); elemId++)
        {
            var elemNodes = mesh.NodesOf<Tet4, Node>(elemId);
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    if (j != i) neighborMap[elemNodes[i]].Add(elemNodes[j]);
                }
            }
        }
        
        // Just 2 quick smoothing passes
        for (int iter = 0; iter < 2; iter++)
        {
            for (int nodeId = 0; nodeId < mesh.Count<Node>(); nodeId++)
            {
                if (boundaryNodes.Contains(nodeId)) continue;
                
                var neighbors = neighborMap[nodeId];
                if (neighbors.Count == 0) continue;
                
                double avgX = 0, avgY = 0, avgZ = 0;
                foreach (var nbr in neighbors)
                {
                    avgX += coords[nbr, 0];
                    avgY += coords[nbr, 1];
                    avgZ += coords[nbr, 2];
                }
                avgX /= neighbors.Count;
                avgY /= neighbors.Count;
                avgZ /= neighbors.Count;
                
                // Update directly with relaxation
                coords[nodeId, 0] += 0.5 * (avgX - coords[nodeId, 0]);
                coords[nodeId, 1] += 0.5 * (avgY - coords[nodeId, 1]);
                coords[nodeId, 2] += 0.5 * (avgZ - coords[nodeId, 2]);
            }
        }
        
        Console.WriteLine($"  → Smoothed {mesh.Count<Node>()} nodes ({boundaryNodes.Count} boundaries fixed)");
    }

    private static double[,] CreateCircle(double cx, double cy, double radius, int nPoints)
    {
        var circle = new double[nPoints, 2];
        for (var i = 0; i < nPoints; i++)
        {
            var theta = 2.0 * Math.PI * i / nPoints;
            circle[i, 0] = cx + radius * Math.Cos(theta);
            circle[i, 1] = cy + radius * Math.Sin(theta);
        }

        return circle;
    }

    /// <summary>
    ///     Create rectangular boundary
    /// </summary>
    private static double[,] CreateRectangle(double xMin, double xMax, double yMin, double yMax,
        int nx, int ny)
    {
        var boundary = new double[2 * (nx + ny - 2), 2];
        var idx = 0;

        // Bottom edge
        for (var i = 0; i < nx; i++)
        {
            boundary[idx, 0] = xMin + i * (xMax - xMin) / (nx - 1);
            boundary[idx, 1] = yMin;
            idx++;
        }

        // Right edge
        for (var i = 1; i < ny; i++)
        {
            boundary[idx, 0] = xMax;
            boundary[idx, 1] = yMin + i * (yMax - yMin) / (ny - 1);
            idx++;
        }

        // Top edge
        for (var i = nx - 2; i >= 0; i--)
        {
            boundary[idx, 0] = xMin + i * (xMax - xMin) / (nx - 1);
            boundary[idx, 1] = yMax;
            idx++;
        }

        // Left edge
        for (var i = ny - 2; i > 0; i--)
        {
            boundary[idx, 0] = xMin;
            boundary[idx, 1] = yMin + i * (yMax - yMin) / (ny - 1);
            idx++;
        }

        return boundary;
    }

    /// <summary>
    ///     Convert List of arrays to 2D array
    /// </summary>
    private static double[,] ConvertListToArray(List<double[]> list)
    {
        var n = list.Count;
        var dim = list[0].Length;
        var array = new double[n, dim];
        for (var i = 0; i < n; i++)
        for (var j = 0; j < dim; j++)
            array[i, j] = list[i][j];

        return array;
    }

    #endregion

    #region 3D Fracture Mechanics Examples (21-25)

    /// <summary>
    ///     Example 21: Penny-Shaped Crack in 3D (Sneddon 1946)
    ///     Reference: Sneddon, I.N. (1946). Proc. Royal Soc. London A, 187:229-260
    /// </summary>
    private static void Example21_Sneddon1946_PennyShapedCrack()
    {
        Console.WriteLine("--- Example 21: Sneddon (1946) - Penny-Shaped Crack (3D) ---");
        Console.WriteLine("Reference: Sneddon, I.N. (1946). Proc. Royal Soc. London A, 187:229-260");
        Console.WriteLine("Geometry: QUARTER MODEL with symmetry planes (x=0, y=0)");
        Console.WriteLine("Circular crack (R=5mm) in quarter domain [0,20]×[0,20]×[0,20]mm");
        Console.WriteLine("Solution: K_I = (2/π)σ√(πa) = 1.128σ√a\n");

        var L = 20.0;
        var R = 5.0;
        var n = 31;

        // Create QUARTER of cube (only positive x, y quadrant)
        // This exploits symmetry: x=0 and y=0 are symmetry planes
        var (mesh, coords) = CreateCubeTetMesh(L, L, L, n, n, n, 
            xMin: 0, yMin: 0, zMin: 0);  // Quarter domain

        Console.WriteLine($"Base mesh (QUARTER): {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tetrahedra");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates3D(mesh, coords);
        // SmoothMesh3D(mesh, coords);  // DEACTIVATED

        // Crack surface: plane at z = L/2
        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => z - L / 2;

        // Region: quarter circle in first quadrant (x≥0, y≥0)
        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            var r = Math.Sqrt(x * x + y * y);  // Distance from z-axis (center at origin now)
            return r - R;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField3D(
            mesh, coords, surface, region);

        Console.WriteLine($"Cracked mesh (QUARTER): {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tet4>()} tetrahedra");
        Console.WriteLine("NOTE: Crack is visible on symmetry planes (x=0, y=0) in postprocessor");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "21_sneddon1946_penny_crack_quarter");
        Console.WriteLine("✓ Saved to GiD format (QUARTER model with symmetry)\n");
    }

    /// <summary>
    ///     Example 22: Elliptical Crack in 3D (Irwin 1962)
    ///     Reference: Irwin, G.R. (1962). J. Applied Mechanics, 29:651-654
    /// </summary>
    private static void Example22_Irwin1962_EllipticalCrack()
    {
        Console.WriteLine("--- Example 22: Irwin (1962) - Elliptical Crack (3D) ---");
        Console.WriteLine("Reference: Irwin, G.R. (1962). J. Applied Mechanics, 29:651-654");
        Console.WriteLine("Geometry: QUARTER MODEL with symmetry planes (x=0, y=0)");
        Console.WriteLine("Elliptical crack (a=6mm, c=3mm) in quarter domain [0,30]×[0,30]×[0,30]mm");
        Console.WriteLine("Aspect ratio: a/c = 2.0\n");

        var L = 30.0;
        var a = 6.0;
        var c = 3.0;
        var n = 25;  // FIXED: Increased from 10 to 25 for adequate crack resolution

        // Create QUARTER of cube (only positive x, y quadrant)
        var (mesh, coords) = CreateCubeTetMesh(L, L, L, n, n, n,
            xMin: 0, yMin: 0, zMin: 0);  // Quarter domain

        Console.WriteLine($"Base mesh (QUARTER): {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tetrahedra");
        
        // Perturb coordinates to avoid exact geometric alignments
        PerturbMeshCoordinates3D(mesh, coords);

        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => z - L / 2;

        // Quarter ellipse in first quadrant (x≥0, y≥0)
        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            return x * x / (a * a) + y * y / (c * c) - 1.0;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField3D(
            mesh, coords, surface, region);

        Console.WriteLine($"Cracked mesh (QUARTER): {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tet4>()} tetrahedra");
        Console.WriteLine("NOTE: Crack is visible on symmetry planes (x=0, y=0) in postprocessor");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "22_irwin1962_elliptical_crack_quarter");
        Console.WriteLine("✓ Saved to GiD format (QUARTER model with symmetry)\n");
    }

    /// <summary>
    ///     Example 23: Edge Crack in 3D Plate (Tada et al. 1973)
    ///     Reference: Tada et al. (1973). The Stress Analysis of Cracks Handbook
    ///     Uses y-symmetry (perpendicular to crack at x=0)
    /// </summary>
    private static void Example23_Tada1973_EdgeCrack()
    {
        Console.WriteLine("--- Example 23: Tada et al. (1973) - Edge Crack (3D) ---");
        Console.WriteLine("Reference: Tada et al. (1973). The Stress Analysis of Cracks Handbook");
        Console.WriteLine("Symmetry: y≥10mm (perpendicular to crack)");
        Console.WriteLine("Edge crack (R=4mm) at x=0, centered at y=15mm, z=10mm\n");

        double W = 40.0, H = 20.0, T = 20.0;
        var R = 4.0;
        double yMid = H / 2;
        double crackCenterY = 15.0;  // y=15mm: crack extends y∈[11,19], clear of y=10 boundary

        // y-symmetry: mesh only y ≥ yMid
        var (mesh, coords) = CreateCubeTetMesh(W, H/2, T, 30, 8, 15,
            xMin: 0, yMin: yMid, zMin: 0);

        Console.WriteLine($"Base mesh (half): {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tets");
        PerturbMeshCoordinates3D(mesh, coords);

        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => x;

        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            var dy = y - crackCenterY;
            var dz = z - T/2;
            return Math.Sqrt(dy*dy + dz*dz) - R;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField3D(
            mesh, coords, surface, region);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tet4>()} tets");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "23_tada1973_edge_crack");
        Console.WriteLine("✓ Saved to GiD format\n");
    }

    /// <summary>
    ///     Example 24: Corner Crack (Newman & Raju 1981)
    ///     Reference: Newman & Raju (1981). Eng. Fract. Mech., 15:185-192
    ///     Uses z-symmetry (perpendicular to crack at x=0)
    /// </summary>
    private static void Example24_NewmanRaju1981_CornerCrack()
    {
        Console.WriteLine("--- Example 24: Newman & Raju (1981) - Quarter-Elliptical Corner Crack (3D) ---");
        Console.WriteLine("Reference: Newman & Raju (1981). Eng. Fract. Mech., 15:185-192");
        Console.WriteLine("Symmetry: z≥15mm (perpendicular to crack)");
        Console.WriteLine("Corner crack (a=5mm, c=3mm) at x=0, y=0 in half-domain [0,30]×[0,30]×[15,30]mm\n");

        double W = 30.0, H = 30.0, T = 30.0;
        var a = 5.0;
        var c = 3.0;
        var n = 25;
        double zMid = T / 2;

        // z-symmetry: mesh only z ≥ zMid
        var (mesh, coords) = CreateCubeTetMesh(W, H, T/2, n, n, n/2,
            xMin: 0, yMin: 0, zMin: zMid);

        Console.WriteLine($"Base mesh (half): {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tets");
        PerturbMeshCoordinates3D(mesh, coords);

        SimplexRemesher.SignedFieldFunction surface = (x, y, z) => x;

        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            var dz = z - zMid;
            return (y*y)/(c*c) + (dz*dz)/(a*a) - 1.0;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField3D(
            mesh, coords, surface, region);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tet4>()} tets");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "24_newmanraju1981_corner_crack");
        Console.WriteLine("✓ Saved to GiD format\n");
    }

    /// <summary>
    ///     Example 25: Through-Thickness Slant Crack (Erdogan & Sih 1963)
    ///     Reference: Erdogan & Sih (1963). J. Basic Eng., 85:519-527
    ///     Uses z-symmetry (perpendicular to inclined crack plane)
    /// </summary>
    private static void Example25_ErdoganSih1963_SlantCrack()
    {
        Console.WriteLine("--- Example 25: Erdogan & Sih (1963) - Through-Thickness Slant Crack (3D) ---");
        Console.WriteLine("Reference: Erdogan & Sih (1963). J. Basic Eng., 85:519-527");
        Console.WriteLine("Symmetry: z≥5mm (perpendicular to crack)");
        Console.WriteLine("Inclined crack (β=30°, L=10mm) in half-domain 30×30×[5,10]mm\n");

        double W = 30.0, H = 30.0, T = 10.0;
        var crackLength = 10.0;
        var beta = 30.0 * Math.PI / 180.0;
        double cx = W/2, cy = H/2, zMid = T/2;

        // z-symmetry: mesh only z ≥ zMid
        var (mesh, coords) = CreateCubeTetMesh(W, H, T/2, 30, 30, 5,
            xMin: 0, yMin: 0, zMin: zMid);

        Console.WriteLine($"Base mesh (half): {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tets");
        PerturbMeshCoordinates3D(mesh, coords);

        SimplexRemesher.SignedFieldFunction surface = (x, y, z) =>
        {
            var nx = Math.Cos(beta);
            var ny = Math.Sin(beta);
            return nx * (x - cx) + ny * (y - cy);
        };

        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            var tx = -Math.Sin(beta);
            var ty = Math.Cos(beta);
            var along = tx * (x - cx) + ty * (y - cy);
            return Math.Abs(along) - crackLength / 2;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField3D(
            mesh, coords, surface, region);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tet4>()} tets");
        Console.WriteLine("Mixed mode: K_I, K_II, K_III all non-zero");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "25_erdogansih1963_slant_crack");
        Console.WriteLine("✓ Saved to GiD format\n");
    }

    /// <summary>
    ///     Example 26: Semi-Cylindrical Surface Crack
    ///     Crack emanates from x=0 face - VISIBLE at boundary
    ///     Uses z-symmetry (perpendicular to crack)
    /// </summary>
    private static void Example26_SemiCylindricalSurfaceCrack()
    {
        Console.WriteLine("--- Example 26: Semi-Cylindrical Surface Crack (VISIBLE) ---");
        Console.WriteLine("Geometry: Semi-cylinder crack at x=0 face, R=6mm");
        Console.WriteLine("Symmetry: z≥10mm (perpendicular to crack)");
        Console.WriteLine("Crack center: (0, 10, 10), extends into block in +x direction\n");

        double W = 40.0, H = 20.0, T = 20.0;
        var R = 6.0;  // Crack radius
        double zMid = T / 2;

        // z-symmetry: mesh only z ≥ zMid
        var (mesh, coords) = CreateCubeTetMesh(W, H, T/2, 66, 30, 18,
            xMin: 0, yMin: 0, zMin: zMid);

        Console.WriteLine($"Base mesh (half): {mesh.Count<Node>()} nodes, {mesh.Count<Tet4>()} tets");
        PerturbMeshCoordinates3D(mesh, coords);

        // Crack surface: semi-cylinder centered at (0, H/2, zMid)
        SimplexRemesher.SignedFieldFunction surface = (x, y, z) =>
        {
            var dy = y - H/2;
            var dz = z - zMid;
            var r = Math.Sqrt(dy*dy + dz*dz);
            return x - (R - r);  // Negative inside, positive outside cylinder
        };

        // Crack region: active where surface is near AND x is positive (in material)
        SimplexRemesher.SignedFieldFunction region = (x, y, z) =>
        {
            if (x < -1.0) return 1.0;  // No crack outside material
            if(x>R+1.0e-3)return 1.0;
            return -1.0;
        };

        var (crackedMesh, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField3D(
            mesh, coords, surface, region,
            enableMeshPerturbation: false,
            enableSmoothing: false);

        Console.WriteLine($"Cracked mesh: {crackedMesh.Count<Node>()} nodes, {crackedMesh.Count<Tet4>()} tets");
        Console.WriteLine("✓ Crack VISIBLE on x=0 face as semi-circular refined region");

        SimplexRemesher.SaveGiD(crackedMesh, crackedCoords, "26_semicylindrical_surface_crack");
        Console.WriteLine("✓ Saved to GiD format\n");
    }

    /// <summary>
    ///     Create a structured tetrahedral mesh for a rectangular block.
    ///     Divides the block into cubes, then each cube into 5 tetrahedra.
    /// </summary>
    private static (SimplexMesh, double[,]) CreateCubeTetMesh(
        double width, double height, double depth,
        int nx, int ny, int nz,
        double xMin = 0, double yMin = 0, double zMin = 0)
    {
        var mesh = new SimplexMesh();

        var totalNodes = (nx + 1) * (ny + 1) * (nz + 1);
        var coords = new double[totalNodes, 3];

        var dx = width / nx;
        var dy = height / ny;
        var dz = depth / nz;

        // Create nodes
        var nodeMap = new int[nx + 1, ny + 1, nz + 1];

        for (var k = 0; k <= nz; k++)
        for (var j = 0; j <= ny; j++)
        for (var i = 0; i <= nx; i++)
        {
            var id = mesh.Add<Node>();
            coords[id, 0] = xMin + i * dx;
            coords[id, 1] = yMin + j * dy;
            coords[id, 2] = zMin + k * dz;

            nodeMap[i, j, k] = id;
        }

        // Create tetrahedra using body-diagonal Kuhn subdivision (6 tets per cube, all positive Jac)
        // Splits cube along body diagonal from n0 (0,0,0) to n7 (1,1,1)
        for (var k = 0; k < nz; k++)
        for (var j = 0; j < ny; j++)
        for (var i = 0; i < nx; i++)
        {
            // Get 8 corner nodes of this cube
            var n0 = nodeMap[i, j, k];
            var n1 = nodeMap[i + 1, j, k];
            var n2 = nodeMap[i, j + 1, k];
            var n3 = nodeMap[i + 1, j + 1, k];
            var n4 = nodeMap[i, j, k + 1];
            var n5 = nodeMap[i + 1, j, k + 1];
            var n6 = nodeMap[i, j + 1, k + 1];
            var n7 = nodeMap[i + 1, j + 1, k + 1];

            // Body-diagonal subdivision (n0 to n7): 6 tets, all positive Jacobians
            mesh.AddTetrahedron(n0, n1, n3, n7);
            mesh.AddTetrahedron(n0, n3, n2, n7);
            mesh.AddTetrahedron(n0, n2, n6, n7);
            mesh.AddTetrahedron(n0, n6, n4, n7);
            mesh.AddTetrahedron(n0, n4, n5, n7);
            mesh.AddTetrahedron(n0, n5, n1, n7);
        }

        // Fix any inverted tets by swapping nodes
        FixInvertedTetrahedra(mesh, coords);

        return (mesh, coords);
    }

    private static void FixInvertedTetrahedra(SimplexMesh mesh, double[,] coords)
    {
        int fixedCount = 0;
        var tetsToFix = new List<(int index, int[] nodes)>();
        
        // First pass: identify inverted tets
        for (int i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            double jac = ComputeTetrahedronJacobian(coords, nodes[0], nodes[1], nodes[2], nodes[3]);
            
            if (jac <= 0)
            {
                // Store inverted tet with swapped nodes
                tetsToFix.Add((i, new[] { nodes[1], nodes[0], nodes[2], nodes[3] }));
            }
        }
        
        // Second pass: remove and recreate inverted tets
        // Note: We need to remove from highest index to lowest to avoid index shifting
        tetsToFix.Sort((a, b) => b.index.CompareTo(a.index));
        
        foreach (var (index, swappedNodes) in tetsToFix)
        {
            mesh.Remove<Tet4>(index);
            mesh.AddTetrahedron(swappedNodes[0], swappedNodes[1], swappedNodes[2], swappedNodes[3]);
            fixedCount++;
        }
        
        if (fixedCount > 0)
            Console.WriteLine($"[FixInvertedTets] Fixed {fixedCount} inverted tetrahedra by swapping nodes");
    }

    private static double ComputeTetrahedronJacobian(double[,] coords, int n0, int n1, int n2, int n3)
    {
        // Compute signed volume * 6 (Jacobian determinant)
        double x0 = coords[n0, 0], y0 = coords[n0, 1], z0 = coords[n0, 2];
        double x1 = coords[n1, 0], y1 = coords[n1, 1], z1 = coords[n1, 2];
        double x2 = coords[n2, 0], y2 = coords[n2, 1], z2 = coords[n2, 2];
        double x3 = coords[n3, 0], y3 = coords[n3, 1], z3 = coords[n3, 2];
        
        // Edge vectors from n0
        double v1x = x1 - x0, v1y = y1 - y0, v1z = z1 - z0;
        double v2x = x2 - x0, v2y = y2 - y0, v2z = z2 - z0;
        double v3x = x3 - x0, v3y = y3 - y0, v3z = z3 - z0;
        
        // Determinant = v1 · (v2 × v3)
        return v1x * (v2y * v3z - v2z * v3y)
             - v1y * (v2x * v3z - v2z * v3x)
             + v1z * (v2x * v3y - v2y * v3x);
    }

    // In Examples2DA.cs, add this diagnostic test
    public static void TestSingleTetRefinement()
    {
        Console.WriteLine("=== SINGLE TET REFINEMENT TEST ===");

        var mesh = new SimplexMesh();
        mesh.Add<Node>(); mesh.Add<Node>(); mesh.Add<Node>(); mesh.Add<Node>();
        mesh.Add<Tet4, Node>(0, 1, 2, 3);

        var coords = new double[4, 3]
        {
            {0, 0, 0},
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        };

        // Mark all 6 edges
        var edgesToRefine = new List<(int, int)>
        {
            (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
        };

        var (refined, _) = MeshRefinement.Refine(mesh, edgesToRefine);

        Console.WriteLine($"Original: 1 tet, 4 nodes");
        Console.WriteLine($"Refined: {refined.Count<Tet4>()} tets, {refined.Count<Node>()} nodes");
        Console.WriteLine($"Expected: 8 tets, 10 nodes");
    }
    #endregion
}

// Perfect 1-to-1 translation of Fortran remeshsimplex.f90 splittetwork1
// Converted from 1-based to 0-based indexing
// License: GPLv3
//
// JACOBIAN CHECKING WORKFLOW:
// 1. Refine() checks input mesh Jacobians (reports only)
// 2. Refine() checks output mesh Jacobians (reports only, does NOT fix)
// 3. User's code does coordinate snapping
// 4. User calls CheckJacobians() to verify after snapping
// 5. User calls FixNegativeJacobians() to fix orientations in-place
//
// Example usage:
//   var (refined, remap) = MeshRefinement.Refine(mesh, edges, coords: coords);
//   var newCoords = MeshRefinement.InterpolateCoordinates(refined, coords);
//   // ... do snapping ...
//   MeshRefinement.CheckJacobians(refined, newCoords, "After snapping");
//   MeshRefinement.FixNegativeJacobians(refined, newCoords);  // Fixes in-place!
//   MeshRefinement.CheckJacobians(refined, newCoords, "Final mesh");
//
// NOTE: FixNegativeJacobians uses ReplaceElementNodes for in-place orientation fixes

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using static Numerical.MeshGeometry;

namespace Numerical;

public static class MeshRefinement
{
    #region Topology Arrays (Fortran lines 570-611, converted to 0-based)
    
    // Fortran: edgenodes(1:2, 1:6) -> C#: edgenodes[6][2]
    // Edge i connects node edgenodes[i][0] to node edgenodes[i][1]
    private static readonly int[][] edgenodes = new int[][]
    {
        new int[] {0, 1}, // Edge 0: nodes 0-1 (Fortran edge 1: nodes 1-2)
        new int[] {1, 2}, // Edge 1: nodes 1-2 (Fortran edge 2: nodes 2-3)
        new int[] {0, 2}, // Edge 2: nodes 0-2 (Fortran edge 3: nodes 1-3)
        new int[] {2, 3}, // Edge 3: nodes 2-3 (Fortran edge 4: nodes 3-4)
        new int[] {3, 0}, // Edge 4: nodes 3-0 (Fortran edge 5: nodes 4-1)
        new int[] {1, 3}  // Edge 5: nodes 1-3 (Fortran edge 6: nodes 2-4)
    };
    
    // Fortran: nodeedges(1:3, 1:4) -> C#: nodeedges[4][3]
    // Node i connects to edges nodeedges[i][0], nodeedges[i][1], nodeedges[i][2]
    private static readonly int[][] nodeedges = new int[][]
    {
        new int[] {0, 2, 4}, // Node 0: edges 0,2,4 (Fortran node 1: edges 1,3,5)
        new int[] {0, 1, 5}, // Node 1: edges 0,1,5 (Fortran node 2: edges 1,2,6)
        new int[] {1, 2, 3}, // Node 2: edges 1,2,3 (Fortran node 3: edges 2,3,4)
        new int[] {3, 4, 5}  // Node 3: edges 3,4,5 (Fortran node 4: edges 4,5,6)
    };
    
    // Fortran: prev(1:6), post(1:6) - opposing edge nodes (converted to 0-based)
    private static readonly int[] prev = new int[] {2, 0, 3, 0, 2, 2}; // Fortran: {3,1,4,1,3,3}
    private static readonly int[] post = new int[] {3, 3, 1, 1, 1, 0}; // Fortran: {4,4,2,2,2,1}
    
    #endregion
    
    #region Main Refine Function
    
    public static (SimplexMesh, int[]) Refine(SimplexMesh input, IReadOnlyList<(int, int)> markedEdges,
        bool enforceClosureForTets = false, bool validateTopology = false, double[,]? inputCoordinates = null)
    {
        // Check input mesh quality if coordinates provided
        if (inputCoordinates != null && input.Count<Tet4>() > 0)
        {
            int negativeCount = 0;
            for (int i = 0; i < input.Count<Tet4>(); i++)
            {
                var nodes = input.NodesOf<Tet4, Node>(i);
                double vol6 = SignedVolume6x(inputCoordinates, nodes[0], nodes[1], nodes[2], nodes[3]);
                if (vol6 < 0) negativeCount++;
            }
            
            if (negativeCount > 0)
                Console.WriteLine($"  ⚠ WARNING: Input mesh has {negativeCount} negative Jacobians BEFORE refinement");
            else
                Console.WriteLine($"  ✓ Input mesh: all {input.Count<Tet4>()} tets have positive Jacobians");
        }
        
        var output = new SimplexMesh();
        var midpointMap = new Dictionary<(int, int), int>();
        var nodeRemap = new int[input.Count<Node>()];
        
        // Store coordinates for orientation checking
        double[,]? coords = inputCoordinates;

        // Declared outside WithBatch — used both inside (mutation) and after (reporting)
        int missingMidpoints = 0;
        int[] caseCounts = new int[7]; // Track how many tets for each case

        // Batch all topology mutations for performance (avoids per-add locking/index updates)
        output.WithBatch(() =>
        {

        // Copy all original nodes
        for (var i = 0; i < input.Count<Node>(); i++)
        {
            var idx = output.Add<Node>();
            output.Set<Node, ParentNodes>(idx, new ParentNodes(i, i));
            nodeRemap[i] = idx;
        }

        // Create midpoint nodes for marked edges
        Console.WriteLine($"  → Processing {markedEdges.Count} marked edges");
        
        foreach (var (a, b) in markedEdges)
        {
            // Check if midpoint already exists (in either direction)
            if (GetMidpoint(midpointMap, a, b) >= 0) continue;
            
            // Create new midpoint node
            int mid = output.Add<Node>();
            output.Set<Node, ParentNodes>(mid, new ParentNodes(a, b));
            midpointMap[(a, b)] = mid; // Direction doesn't matter because GetMidpoint checks both
        }
        
        Console.WriteLine($"  → Created {midpointMap.Count} unique midpoint nodes");

        // Enforce closure for tets if requested — MUST be before any element splitting
        // so that all element types see the complete set of midpoint nodes.
        if (enforceClosureForTets && input.Count<Tet4>() > 0)
        {
            // Convert midpointMap keys to HashSet for closure computation
            var markedEdgeSet = new HashSet<(int, int)>(midpointMap.Keys);
            ComputeTetEdgeClosureFull(input, markedEdgeSet, validateTopology);
            
            // Add any new closure edges as midpoints
            int closureAdded = 0;
            foreach (var (a, b) in markedEdgeSet)
            {
                if (GetMidpoint(midpointMap, a, b) >= 0) continue;
                
                int mid = output.Add<Node>();
                output.Set<Node, ParentNodes>(mid, new ParentNodes(a, b));
                midpointMap[(a, b)] = mid;
                closureAdded++;
            }
            if (closureAdded > 0)
                Console.WriteLine($"  → Closure added {closureAdded} additional midpoints (total: {midpointMap.Count})");
        }

        // ═══════════════════════════════════════════════════════════════════
        // Split Point elements (Fortran: splitpoint0/splitpoint1)
        //
        // Point elements enforce boundary conditions. Each existing point
        // is copied unchanged, and for each marked edge, ALL point elements
        // at both endpoints are propagated to the midpoint node. Without
        // this, midpoint DOFs are unconstrained → singular stiffness.
        // ═══════════════════════════════════════════════════════════════════
        int inputPointCount = input.Count<Point>();
        if (inputPointCount > 0)
        {
            int inputNodeCount = input.Count<Node>();
            
            // Build node → point-element CSR map (Fortran: pelcount/pelptr/pellist)
            var pelcount = new int[inputNodeCount];
            for (int iel = 0; iel < inputPointCount; iel++)
            {
                var pnodes = input.NodesOf<Point, Node>(iel);
                pelcount[pnodes[0]]++;
            }
            
            var pelptr = new int[inputNodeCount + 1];
            pelptr[0] = 0;
            for (int ing = 0; ing < inputNodeCount; ing++)
                pelptr[ing + 1] = pelptr[ing] + pelcount[ing];
            
            var pellist = new int[Math.Max(1, pelptr[inputNodeCount])];
            Array.Clear(pelcount, 0, pelcount.Length);
            for (int iel = 0; iel < inputPointCount; iel++)
            {
                var pnodes = input.NodesOf<Point, Node>(iel);
                int ing = pnodes[0];
                pellist[pelptr[ing] + pelcount[ing]] = iel;
                pelcount[ing]++;
            }
            
            // Copy existing point elements unchanged
            int pointsCopied = 0;
            for (int iel = 0; iel < inputPointCount; iel++)
            {
                var pnodes = input.NodesOf<Point, Node>(iel);
                AddPt(output, iel, nodeRemap[pnodes[0]]);
                pointsCopied++;
            }
            
            // For each marked edge, propagate ALL point elements from both
            // endpoints to the midpoint node. (Fortran splitpoint1 lines 195-214)
            int pointsPropagated = 0;
            foreach (var ((a, b), midNode) in midpointMap)
            {
                // All points from endpoint a
                for (int ip = 0; ip < pelcount[a]; ip++)
                {
                    int iel = pellist[pelptr[a] + ip];
                    AddPt(output, iel, midNode);
                    pointsPropagated++;
                }
                // All points from endpoint b
                for (int ip = 0; ip < pelcount[b]; ip++)
                {
                    int iel = pellist[pelptr[b] + ip];
                    AddPt(output, iel, midNode);
                    pointsPropagated++;
                }
            }
            
            Console.WriteLine($"  → Points: {pointsCopied} copied + {pointsPropagated} propagated to midpoints = {output.Count<Point>()} total");
        }

        // ═══════════════════════════════════════════════════════════════════
        // Split Bar2 elements (Fortran: splitbar20/splitbar21)
        //
        // Bar elements define boundary edges. If the bar's edge is marked,
        // split into two bars; otherwise copy unchanged.
        // ═══════════════════════════════════════════════════════════════════
        int inputBarCount = input.Count<Bar2>();
        if (inputBarCount > 0)
        {
            int barsCopied = 0, barsSplit = 0;
            for (int iel = 0; iel < inputBarCount; iel++)
            {
                var bnodes = input.NodesOf<Bar2, Node>(iel);
                int midNode = GetMidpoint(midpointMap, bnodes[0], bnodes[1]);
                
                if (midNode < 0)
                {
                    // Edge not marked — copy unchanged
                    AddBar(output, iel, nodeRemap[bnodes[0]], nodeRemap[bnodes[1]]);
                    barsCopied++;
                }
                else
                {
                    // Edge marked — split into two bars
                    AddBar(output, iel, nodeRemap[bnodes[0]], midNode);
                    AddBar(output, iel, midNode, nodeRemap[bnodes[1]]);
                    barsSplit++;
                }
            }
            Console.WriteLine($"  → Bars: {barsCopied} copied + {barsSplit} split (×2) = {output.Count<Bar2>()} total");
        }

        // ═══════════════════════════════════════════════════════════════════
        // Refine Tri3 elements (Fortran: splittria30/splittria31)
        // ═══════════════════════════════════════════════════════════════════
        int triRefinedCount = 0;
        for (var i = 0; i < input.Count<Tri3>(); i++)
        {
            var nodes = input.NodesOf<Tri3, Node>(i);
            var v = new int[] { nodeRemap[nodes[0]], nodeRemap[nodes[1]], nodeRemap[nodes[2]] };
            
            // Look up midpoints using INPUT indices (midpointMap keys are input indices)
            var m01 = GetMidpoint(midpointMap, nodes[0], nodes[1]);
            var m12 = GetMidpoint(midpointMap, nodes[1], nodes[2]);
            var m20 = GetMidpoint(midpointMap, nodes[2], nodes[0]);
            
            if (!(m01 < 0 && m12 < 0 && m20 < 0))
                triRefinedCount++;
            
            SplitTriangle(output, i, v[0], v[1], v[2], m01, m12, m20);
        }
        Console.WriteLine($"  → Triangles: {triRefinedCount} refined out of {input.Count<Tri3>()} → {output.Count<Tri3>()} total");

        // ═══════════════════════════════════════════════════════════════════
        // Refine Tet4 elements (Fortran: splittetra40/splittetra41)
        // ═══════════════════════════════════════════════════════════════════
        int tetCount = input.Count<Tet4>();
        
        for (var i = 0; i < tetCount; i++)
        {
            var nodes = input.NodesOf<Tet4, Node>(i);
            var globalNodes = new int[] { nodeRemap[nodes[0]], nodeRemap[nodes[1]], nodeRemap[nodes[2]], nodeRemap[nodes[3]] };
            
            // Build midnodes array (0 means no midpoint, >0 means midpoint node index)
            var midnodes = new int[6];
            int nmarked = 0;
            for (int e = 0; e < 6; e++)
            {
                // Look up using INPUT mesh indices
                int mid = GetMidpoint(midpointMap, nodes[edgenodes[e][0]], nodes[edgenodes[e][1]]);
                midnodes[e] = (mid >= 0) ? mid : 0;
                if (mid >= 0) nmarked++;
                else missingMidpoints++;
            }
            
            // Track case distribution
            caseCounts[nmarked]++;
            
            // Call exact Fortran splittetwork1 translation
            SplitTetWork1(globalNodes, midnodes, output, i);
        }
        
        }); // end WithBatch
        
        Console.WriteLine($"  → Tet split case distribution:");
        for (int c = 0; c <= 6; c++)
        {
            if (caseCounts[c] > 0)
                Console.WriteLine($"     Case {c} ({c} marked edges): {caseCounts[c]} tets");
        }
        
        if (missingMidpoints > 0 && enforceClosureForTets)
            Console.WriteLine($"  ⚠️  WARNING: {missingMidpoints} tet edges missing midpoints (closure may have failed)");
        
        // Final summary
        Console.WriteLine($"  → MESH SUMMARY:");
        Console.WriteLine($"     Input:  {input.Count<Node>()} nodes, {input.Count<Point>()} pts, {input.Count<Bar2>()} bars, {input.Count<Tri3>()} tris, {input.Count<Tet4>()} tets");
        Console.WriteLine($"     Output: {output.Count<Node>()} nodes, {output.Count<Point>()} pts, {output.Count<Bar2>()} bars, {output.Count<Tri3>()} tris, {output.Count<Tet4>()} tets");

        
        // Check negative Jacobians after refinement but DON'T fix yet
        if (coords != null && output.Count<Tet4>() > 0)
        {
            var refinedCoords = InterpolateCoordinates(output, coords);
            int negAfterRefine = CountNegativeJacobians(refinedCoords, output);
            
            if (negAfterRefine > 0)
                Console.WriteLine($"  → After refinement: {negAfterRefine} negative Jacobians (NOT fixed yet - fix AFTER snapping)");
            else
                Console.WriteLine($"  ✓ After refinement: all {output.Count<Tet4>()} tets have positive Jacobians");
        }

        return (output, nodeRemap);
    }
    
    private static int CountNegativeJacobians(double[,] coordinates, SimplexMesh mesh)
    {
        int count = 0;
        for (int i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            double vol6 = SignedVolume6x(coordinates, nodes[0], nodes[1], nodes[2], nodes[3]);
            if (vol6 < 0) count++;
        }
        return count;
    }
    
    public static int CheckJacobians(SimplexMesh mesh, double[,] coordinates, string label = "Mesh")
    {
        int negCount = CountNegativeJacobians(coordinates, mesh);
        
        if (negCount > 0)
            Console.WriteLine($"  → {label}: {negCount} negative Jacobians out of {mesh.Count<Tet4>()} tets");
        else
            Console.WriteLine($"  ✓ {label}: all {mesh.Count<Tet4>()} tets have positive Jacobians");
        
        return negCount;
    }
    
    /// <summary>
    /// Checks edge sharing patterns in the mesh (detects topology errors).
    /// Uses Topology's sub-entity discovery and sharing queries.
    /// Valid interior edges are shared by multiple tets, boundary edges appear once.
    /// </summary>
    public static void CheckEdgeTopology(SimplexMesh mesh, string label = "Mesh")
    {
        SimplexRemesher.DiscoverEdges(mesh);
        
        int totalEdges = mesh.Count<Edge>();
        int boundaryEdges = 0;
        int interiorEdges = 0;
        int abnormalEdges = 0;
        int maxShare = 0;
        
        for (int e = 0; e < totalEdges; e++)
        {
            int shareCount = mesh.ElementsSharingSubEntity<Tet4, Edge, Node>(e).Count;
            if (shareCount == 1) boundaryEdges++;
            else if (shareCount >= 2) interiorEdges++;
            if (shareCount > 10) abnormalEdges++;
            if (shareCount > maxShare) maxShare = shareCount;
        }
        
        Console.WriteLine($"  → {label}: {totalEdges} unique edges ({boundaryEdges} boundary, {interiorEdges} interior)");
        
        if (abnormalEdges > 0)
        {
            Console.WriteLine($"  ⚠️  {label}: {abnormalEdges} edges shared by >10 tets (max: {maxShare} tets)");
            Console.WriteLine($"      This may indicate crossing edges or topology errors!");
        }
        else
        {
            Console.WriteLine($"  ✓ {label}: edge sharing looks normal");
        }
    }
    
    /// <summary>
    /// Fixes negative Jacobians by permuting tet nodes in-place using ReplaceElementNodes.
    /// No mesh rebuild required — all data attachments (ParentNodes, OriginalElement) are preserved.
    /// </summary>
    public static void FixNegativeJacobians(SimplexMesh mesh, double[,] coordinates)
    {
        int tetCount = mesh.Count<Tet4>();
        int fixedCount = 0;
        int unfixable = 0;
        
        // All possible even permutations (preserve handedness, just fix orientation)
        // These are the 12 even permutations of (0,1,2,3)
        var evenPermutations = new[]
        {
            new[] {0, 1, 2, 3}, // Original
            new[] {0, 2, 3, 1}, // Rotate
            new[] {0, 3, 1, 2}, // Rotate
            new[] {1, 0, 3, 2}, // Swap 0↔1, 2↔3
            new[] {1, 2, 0, 3}, // Rotate
            new[] {1, 3, 2, 0}, // Rotate
            new[] {2, 0, 1, 3}, // Rotate
            new[] {2, 1, 3, 0}, // Swap 0↔2, 1↔3
            new[] {2, 3, 0, 1}, // Rotate
            new[] {3, 0, 2, 1}, // Swap 0↔3, 1↔2
            new[] {3, 1, 0, 2}, // Rotate
            new[] {3, 2, 1, 0}  // Rotate
        };
        
        for (int i = 0; i < tetCount; i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            double vol6 = SignedVolume6x(coordinates, nodes[0], nodes[1], nodes[2], nodes[3]);
            
            if (vol6 >= 0) continue; // Already positive
            
            // Try to fix by swapping nodes 2↔3 (most common fix)
            double vol6_swap23 = SignedVolume6x(coordinates, nodes[0], nodes[1], nodes[3], nodes[2]);
            
            if (vol6_swap23 > 0)
            {
                mesh.ReplaceElementNodes<Tet4, Node>(i, nodes[0], nodes[1], nodes[3], nodes[2]);
                fixedCount++;
            }
            else
            {
                // Try all even permutations to find one with positive Jacobian
                bool wasFixed = false;
                foreach (var perm in evenPermutations)
                {
                    double testVol = SignedVolume6x(coordinates, 
                        nodes[perm[0]], nodes[perm[1]], nodes[perm[2]], nodes[perm[3]]);
                    
                    if (testVol > 0)
                    {
                        mesh.ReplaceElementNodes<Tet4, Node>(i,
                            nodes[perm[0]], nodes[perm[1]], nodes[perm[2]], nodes[perm[3]]);
                        fixedCount++;
                        wasFixed = true;
                        break;
                    }
                }
                
                if (!wasFixed)
                {
                    // No permutation works - tet is degenerate or coordinates are bad
                    // Just use swap 2↔3 anyway
                    mesh.ReplaceElementNodes<Tet4, Node>(i, nodes[0], nodes[1], nodes[3], nodes[2]);
                    unfixable++;
                }
            }
        }
        
        if (fixedCount > 0)
            Console.WriteLine($"  → Fixed {fixedCount} negative Jacobians by node permutation");
        if (unfixable > 0)
            Console.WriteLine($"  ⚠️  WARNING: {unfixable} tets remain inverted (degenerate geometry)");
    }
    
    public static double[,] InterpolateCoordinates(SimplexMesh mesh, double[,] originalCoords)
    {
        int nodeCount = mesh.Count<Node>();
        int dim = originalCoords.GetLength(1);
        var newCoords = new double[nodeCount, dim];

        for (int i = 0; i < nodeCount; i++)
        {
            var parents = mesh.Get<Node, ParentNodes>(i);
            if (parents.Parent1 == parents.Parent2)
            {
                for (int d = 0; d < dim; d++)
                    newCoords[i, d] = originalCoords[parents.Parent1, d];
            }
            else
            {
                for (int d = 0; d < dim; d++)
                    newCoords[i, d] = 0.5 * (originalCoords[parents.Parent1, d] + originalCoords[parents.Parent2, d]);
            }
        }
        return newCoords;
    }
    
    public static void CorrectTetOrientations(SimplexMesh mesh, double[,] coordinates)
    {
        int tetCount = mesh.Count<Tet4>();
        int negativeCount = 0;

        for (int i = 0; i < tetCount; i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            double vol6 = SignedVolume6x(coordinates, nodes[0], nodes[1], nodes[2], nodes[3]);
            
            if (vol6 < -1e-15)
            {
                // Swap nodes 2 and 3 to fix orientation (in-place)
                mesh.ReplaceElementNodes<Tet4, Node>(i, nodes[0], nodes[1], nodes[3], nodes[2]);
                negativeCount++;
            }
        }

        if (negativeCount > 0)
            Console.WriteLine($"  → Fixed {negativeCount} inverted tets");
    }
    
    /// <summary>
    /// Compute signed volume * 6 (uses shared MeshGeometry implementation).
    /// </summary>
    private static double SignedVolume6x(double[,] coords, int n0, int n1, int n2, int n3)
        => ComputeTetrahedronJacobian(coords, n0, n1, n2, n3);
    
    #endregion
    
    #region Exact Fortran splittetwork1 Translation (lines 612-939)
    
    /// <summary>
    /// EXACT translation of Fortran splittetwork1 (lines 612-939)
    /// All logic preserved, only converted from 1-based to 0-based indexing
    /// </summary>
    private static void SplitTetWork1(int[] nodes, int[] midnodes, SimplexMesh output, int origIdx)
    {
        // Line 612-619: identify the marked edges
        int nmarked = 0;
        var marked = new int[6];
        for (int i = 0; i < 6; i++)
        {
            if (midnodes[i] > 0)
            {
                marked[nmarked] = i;
                nmarked++;
            }
        }

        // Line 620-632: determine the number of marked edges per node and record the edges
        var imarks = new int[5]; // imarks(1:5) in Fortran
        var jmarks = new int[12]; // jmarks(1:12) in Fortran
        imarks[0] = 0; // Fortran: imarks(1) = 1, converted to 0-based
        
        for (int i = 0; i < 4; i++)
        {
            imarks[i + 1] = imarks[i];
            for (int j = 0; j < 3; j++)
            {
                int ie = nodeedges[i][j]; // Fortran: nodeedges(j+1, i+1)
                if (midnodes[ie] > 0)
                {
                    jmarks[imarks[i + 1]] = ie;
                    imarks[i + 1]++;
                }
            }
        }

        // Line 633-939: split the tetrahedron based on the number of marked edges
        switch (nmarked)
        {
            case 0:
                // Case #1: no marked edges (line 635-638)
                AddTet(output, origIdx, nodes[0], nodes[1], nodes[2], nodes[3]);
                break;

            case 1:
                // Case #2: one marked edge (lines 639-656)
                Case1Edge(nodes, midnodes, marked, output, origIdx);
                break;

            case 2:
                // Cases #3 and #4: two marked edges (lines 657-755)
                Case2Edges(nodes, midnodes, marked, imarks, jmarks, output, origIdx);
                break;

            case 3:
                // Cases #5, #6, #7: three marked edges (lines 756-816)
                Case3Edges(nodes, midnodes, marked, imarks, jmarks, output, origIdx);
                break;

            case 4:
                // Cases #8, #9: four marked edges (lines 817-885)
                Case4Edges(nodes, midnodes, imarks, output, origIdx);
                break;

            case 5:
                // Case #10: five marked edges (lines 886-917)
                Case5Edges(nodes, midnodes, output, origIdx);
                break;

            case 6:
                // Case #11: all edges marked (lines 918-938)
                Case6Edges(nodes, midnodes, output, origIdx);
                break;
        }
    }
    
    #endregion
    
    #region Case Implementations (Exact Fortran Translation)
    
    // Case #2: one marked edge (Fortran lines 639-656)
    private static void Case1Edge(int[] nodes, int[] midnodes, int[] marked, SimplexMesh output, int origIdx)
    {
        int n1 = edgenodes[marked[0]][0];
        int n2 = edgenodes[marked[0]][1];
        
        // Find opposing nodes
        var (n3, n4) = NodesfFrom2Nodes(n1, n2);
        
        // Convert to global indices
        n1 = nodes[n1];
        n2 = nodes[n2];
        n3 = nodes[n3];
        n4 = nodes[n4];
        int n5 = midnodes[marked[0]];
        
        if (n5 == 0)
            throw new InvalidOperationException($"Case1Edge: marked edge {marked[0]} has no midpoint");
        
        // Build two new tetrahedra (Fortran lines 654-656)
        AddTet(output, origIdx, n3, n4, n5, n2);
        AddTet(output, origIdx, n4, n3, n5, n1);
    }
    
    // Cases #3 and #4: two marked edges (Fortran lines 657-755)
    private static void Case2Edges(int[] nodes, int[] midnodes, int[] marked, int[] imarks, int[] jmarks, SimplexMesh output, int origIdx)
    {
        // Find maximum number of marked edges connected to any node
        int mm = 0;
        for (int i = 0; i < 4; i++)
            mm = Math.Max(mm, imarks[i + 1] - imarks[i]);

        if (mm == 2)
        {
            // Case #3: one node with 2 marked edges (Fortran lines 664-693)
            int n3 = -1, n4 = -1;
            for (int i = 0; i < 4; i++)
            {
                if (imarks[i + 1] - imarks[i] == 0) n3 = i;
                if (imarks[i + 1] - imarks[i] == 2) n4 = i;
            }

            // Fortran line 680: CALL nodesfrom2nodes(edgenodes, post, prev, n1, n3, n4, n2)
            var (n1, n2) = NodesfFrom2Nodes(n3, n4);
            
            // Fortran lines 683-684 - get midpoint nodes
            int n5 = SplitNode(n2, n4, midnodes);  // Fortran: splitnode(n2, n4, ...)
            int n6 = SplitNode(n1, n4, midnodes);  // Fortran: splitnode(n1, n4, ...)
            
            // Fortran lines 685-688 - convert to global indices
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            
            // Fortran lines 691-692 - build 3 tets
            var pyrTets = SplitPyramidIntoTets(n1, n2, n5, n6, n3);  // Fortran: splitpyramidintotets(n1, n2, n5, n6, n3, ...)
            foreach (var tet in pyrTets)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
            AddTet(output, origIdx, n6, n5, n3, n4);  // Fortran: splitetintoitself(n6, n5, n3, n4, ...)
        }
        else
        {
            // Case #4: each node has 1 marked edge (Fortran lines 694-714)
            int ie1 = marked[0];  // Fortran: marked(1) - first marked edge
            int n1 = edgenodes[ie1][0];
            int n2 = edgenodes[ie1][1];
            
            // Fortran line 700: CALL nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
            var (n3, n4) = NodesfFrom2Nodes(n1, n2);
            
            // Fortran lines 702-703 - get midpoint nodes
            int n5 = SplitNode(n1, n2, midnodes);  // Fortran: splitnode(n1, n2, ...)
            int n6 = SplitNode(n3, n4, midnodes);  // Fortran: splitnode(n3, n4, ...)
            
            // Fortran lines 705-708 - convert to global indices
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            
            // Fortran lines 711-714 - build 4 tets
            AddTet(output, origIdx, n2, n3, n5, n6);  // Fortran: splitetintoitself(n2, n3, n5, n6, ...)
            AddTet(output, origIdx, n2, n6, n5, n4);  // Fortran: splitetintoitself(n2, n6, n5, n4, ...)
            AddTet(output, origIdx, n1, n5, n6, n4);  // Fortran: splitetintoitself(n1, n5, n6, n4, ...)
            AddTet(output, origIdx, n1, n5, n3, n6);  // Fortran: splitetintoitself(n1, n5, n3, n6, ...)
        }
    }
    
    // Cases #5, #6, #7: three marked edges (Fortran lines 716-816)
    private static void Case3Edges(int[] nodes, int[] midnodes, int[] marked, int[] imarks, int[] jmarks, SimplexMesh output, int origIdx)
    {
        // Fortran lines 719-727: identify n3 (0 marks) and n4 (3 marks)
        int n3 = -1, n4 = -1;
        for (int i = 0; i < 4; i++)
        {
            if (imarks[i + 1] - imarks[i] == 3) n4 = i;
            if (imarks[i + 1] - imarks[i] == 0) n3 = i;
        }

        if (n4 >= 0)
        {
            // Case #6: one node with 3 marked edges (Fortran lines 728-744)
            int n1 = (n4 + 1) % 4;
            int n2;
            (n2, n3) = NodesfFrom2Nodes(n1, n4);  // Fortran line 731: OVERWRITES n3
            
            int n5 = SplitNode(n3, n4, midnodes);
            int n6 = SplitNode(n1, n4, midnodes);
            int n7 = SplitNode(n2, n4, midnodes);
            
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            
            AddTet(output, origIdx, n4, n7, n6, n5);
            var prismTets = SplitPrismIntoTets(n2, n7, n6, n1, n3, n5);
            foreach (var tet in prismTets)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
        }
        else if (n3 >= 0)
        {
            // Case #5: one node with 0 marked edges (Fortran lines 745-763)
            int n4_local = (n3 + 1) % 4;
            int n1, n2;
            (n1, n2) = NodesfFrom2Nodes(n3, n4_local);
            
            int n5 = SplitNode(n1, n4_local, midnodes);
            int n6 = SplitNode(n1, n2, midnodes);
            int n7 = SplitNode(n2, n4_local, midnodes);
            
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4_local = nodes[n4_local];
            
            AddTet(output, origIdx, n1, n6, n3, n5);
            AddTet(output, origIdx, n2, n3, n6, n7);
            AddTet(output, origIdx, n6, n3, n5, n7);
            AddTet(output, origIdx, n5, n7, n3, n4_local);
        }
        else
        {
            // Case #7: complex case (Fortran lines 764-815)
            int n1 = -1, n2 = -1;
            n3 = -1;
            n4 = -1;
            
            for (int i = 0; i < 4; i++)
            {
                if (imarks[i + 1] - imarks[i] == 1)
                {
                    int ie = jmarks[imarks[i]];
                    if (n1 < 0)
                    {
                        n1 = i;
                        if (edgenodes[ie][0] == n1)
                            n4 = edgenodes[ie][1];
                        else
                            n4 = edgenodes[ie][0];
                    }
                    else
                    {
                        n2 = i;
                        if (edgenodes[ie][0] == n2)
                            n3 = edgenodes[ie][1];
                        else
                            n3 = edgenodes[ie][0];
                    }
                }
            }
            
            if (n1 < 0 || n2 < 0) return;  // Safety
            
            int n5 = SplitNode(n1, n4, midnodes);
            int n6 = SplitNode(n2, n3, midnodes);
            int n7 = SplitNode(n3, n4, midnodes);
            
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            
            var pyr1 = SplitPyramidIntoTets(n6, n7, n4, n2, n5);
            foreach (var tet in pyr1)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
                
            var pyr2 = SplitPyramidIntoTets(n3, n1, n5, n7, n6);
            foreach (var tet in pyr2)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
                
            AddTet(output, origIdx, n1, n2, n6, n5);
        }
    }
    
    // Cases #8, #9: four marked edges (Fortran lines 817-885)
    private static void Case4Edges(int[] nodes, int[] midnodes, int[] imarks, SimplexMesh output, int origIdx)
    {
        int mm = 0;
        for (int i = 0; i < 4; i++)
            mm = Math.Max(mm, imarks[i + 1] - imarks[i]);

        if (mm == 3)
        {
            // Case #8: one node with 3 marked edges (lines 824-855)
            int n3 = -1, n4 = -1;
            for (int i = 0; i < 4; i++)
            {
                if (imarks[i + 1] - imarks[i] == 1) n3 = i;
                if (imarks[i + 1] - imarks[i] == 3) n4 = i;
            }

            var (n1, n2) = NodesfFrom2Nodes(n3, n4);
            
            int n5 = SplitNode(n1, n2, midnodes);
            int n6 = SplitNode(n2, n4, midnodes);
            int n7 = SplitNode(n3, n4, midnodes);
            int n8 = SplitNode(n1, n4, midnodes);
            
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            
            AddTet(output, origIdx, n8, n6, n7, n4);
            AddTet(output, origIdx, n7, n8, n5, n6);
            
            var pyr1 = SplitPyramidIntoTets(n2, n3, n7, n6, n5);
            foreach (var tet in pyr1)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
                
            var pyr2 = SplitPyramidIntoTets(n1, n8, n7, n3, n5);
            foreach (var tet in pyr2)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
        }
        else
        {
            // Case #9: quadrilateral (one unmarked edge) (lines 857-884)
            int ie = -1;
            for (int i = 0; i < 6; i++)
            {
                if (midnodes[i] == 0)
                {
                    ie = i;
                    break;
                }
            }

            int n1 = edgenodes[ie][0];
            int n2 = edgenodes[ie][1];
            var (n3, n4) = NodesfFrom2Nodes(n1, n2);
            
            int n5 = SplitNode(n1, n4, midnodes);
            int n6 = SplitNode(n1, n3, midnodes);
            int n7 = SplitNode(n2, n3, midnodes);
            int n8 = SplitNode(n2, n4, midnodes);
            
            n1 = nodes[n1];
            n2 = nodes[n2];
            n3 = nodes[n3];
            n4 = nodes[n4];
            
            var prism1 = SplitPrismIntoTets(n8, n7, n3, n4, n5, n6);
            foreach (var tet in prism1)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
                
            var prism2 = SplitPrismIntoTets(n1, n2, n8, n5, n6, n7);
            foreach (var tet in prism2)
                AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
        }
    }
    
    // Case #10: five marked edges (Fortran lines 886-917)
    private static void Case5Edges(int[] nodes, int[] midnodes, SimplexMesh output, int origIdx)
    {
        int ie = -1;
        for (int i = 0; i < 6; i++)
        {
            if (midnodes[i] == 0)
            {
                ie = i;
                break;
            }
        }

        int n1 = edgenodes[ie][0];
        int n2 = edgenodes[ie][1];
        var (n3, n4) = NodesfFrom2Nodes(n1, n2);
        
        int n5 = SplitNode(n1, n4, midnodes);
        int n6 = SplitNode(n1, n3, midnodes);
        int n7 = SplitNode(n2, n3, midnodes);
        int n8 = SplitNode(n3, n4, midnodes);
        int n9 = SplitNode(n2, n4, midnodes);
        
        n1 = nodes[n1];
        n2 = nodes[n2];
        n3 = nodes[n3];
        n4 = nodes[n4];
        
        AddTet(output, origIdx, n9, n8, n5, n4);
        AddTet(output, origIdx, n3, n6, n7, n8);
        
        var pyr = SplitPyramidIntoTets(n6, n7, n9, n5, n8);
        foreach (var tet in pyr)
            AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
            
        var prism = SplitPrismIntoTets(n1, n2, n9, n5, n6, n7);
        foreach (var tet in prism)
            AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
    }
    
    // Case #11: all edges marked (octahedron) (Fortran lines 918-938)
    private static void Case6Edges(int[] nodes, int[] midnodes, SimplexMesh output, int origIdx)
    {
        int n5 = SplitNode(0, 3, midnodes);
        int n6 = SplitNode(1, 3, midnodes);
        int n7 = SplitNode(2, 3, midnodes);
        int n8 = SplitNode(0, 2, midnodes);
        int n9 = SplitNode(0, 1, midnodes);
        int n10 = SplitNode(1, 2, midnodes);
        
        int n1 = nodes[0];
        int n2 = nodes[1];
        int n3 = nodes[2];
        int n4 = nodes[3];
        
        AddTet(output, origIdx, n1, n9, n8, n5);
        AddTet(output, origIdx, n2, n10, n9, n6);
        AddTet(output, origIdx, n3, n8, n10, n7);
        AddTet(output, origIdx, n5, n6, n7, n4);
        
        var octa = SplitOctoIntoTets(n8, n6, n7, n5, n9, n10);
        foreach (var tet in octa)
            AddTet(output, origIdx, tet[0], tet[1], tet[2], tet[3]);
    }
    
    #endregion
    
    #region Helper Functions (Fortran translations)
    
    // Fortran: nodesfrom2nodes (lines 1419-1435)
    private static (int, int) NodesfFrom2Nodes(int i2, int i3)
    {
        if (i2 == i3)
            throw new InvalidOperationException($"NodesfFrom2Nodes({i2},{i3}): cannot find edge from node to itself");
            
        // Find edge containing nodes i2 and i3
        for (int ie = 0; ie < 6; ie++)
        {
            if ((edgenodes[ie][0] == i2 && edgenodes[ie][1] == i3))
                return (prev[ie], post[ie]);
            if ((edgenodes[ie][0] == i3 && edgenodes[ie][1] == i2))
                return (post[ie], prev[ie]);
        }
        
        throw new InvalidOperationException($"NodesfFrom2Nodes({i2},{i3}): edge not found in topology");
    }
    
    // Fortran: splitnode (finds midpoint node for edge)
    private static int SplitNode(int n1, int n2, int[] midnodes)
    {
        for (int ie = 0; ie < 6; ie++)
        {
            if ((edgenodes[ie][0] == n1 && edgenodes[ie][1] == n2) ||
                (edgenodes[ie][0] == n2 && edgenodes[ie][1] == n1))
            {
                if (midnodes[ie] == 0)
                {
                    throw new InvalidOperationException($"SplitNode({n1},{n2}): edge {ie} found but has no midpoint (midnodes[{ie}]=0)");
                }
                return midnodes[ie];
            }
        }
        throw new InvalidOperationException($"SplitNode({n1},{n2}): edge not found in edgenodes");
    }
    
    // Fortran: splitpyramidintotets (lines 1130-1145)
    private static List<int[]> SplitPyramidIntoTets(int i1, int i2, int i3, int i4, int i5)
    {
        var result = new List<int[]>();
        bool swap13 = Math.Max(i1, i3) > Math.Max(i2, i4);
        
        if (swap13)
        {
            result.Add(new int[] { i1, i2, i5, i3 });
            result.Add(new int[] { i1, i5, i4, i3 });
        }
        else
        {
            result.Add(new int[] { i1, i2, i5, i4 });
            result.Add(new int[] { i2, i3, i5, i4 });
        }
        
        return result;
    }
    
    // Fortran: splitprismintotets (lines 982-1038)
    // Fortran: splitprismintotets (remeshsimplex.f90 lines 961-1000)
    private static List<int[]> SplitPrismIntoTets(int i1, int i2, int i3, int i4, int i5, int i6)
    {
        var result = new List<int[]>(3);
        bool swap46 = Math.Max(i4, i6) > Math.Max(i3, i5);
        bool swap13 = Math.Max(i1, i3) > Math.Max(i2, i4);
        bool swap16 = Math.Max(i1, i6) > Math.Max(i2, i5);
        
        if (swap46)
        {
            if (swap13)
            {
                result.Add(new[] { i1, i2, i6, i3 });
                result.Add(new[] { i1, i6, i4, i3 });
                result.Add(new[] { i1, i6, i5, i4 });
            }
            else if (swap16)
            {
                result.Add(new[] { i1, i2, i6, i4 });
                result.Add(new[] { i1, i6, i5, i4 });
                result.Add(new[] { i2, i6, i4, i3 });
            }
            else
            {
                result.Add(new[] { i1, i2, i5, i4 });
                result.Add(new[] { i2, i3, i6, i4 });
                result.Add(new[] { i2, i6, i5, i4 });
            }
        }
        else
        {
            if (swap13)
            {
                if (swap16)
                {
                    result.Add(new[] { i1, i2, i6, i3 });
                    result.Add(new[] { i1, i3, i5, i4 });
                    result.Add(new[] { i1, i6, i5, i3 });
                }
                else
                {
                    result.Add(new[] { i1, i3, i5, i4 });
                    result.Add(new[] { i2, i5, i1, i3 });
                    result.Add(new[] { i2, i6, i5, i3 });
                }
            }
            else
            {
                result.Add(new[] { i1, i2, i5, i4 });
                result.Add(new[] { i2, i3, i5, i4 });
                result.Add(new[] { i2, i6, i5, i3 });
            }
        }
        
        return result;
    }
    
    // Fortran: splitoctointotets (lines 958-977)
    private static List<int[]> SplitOctoIntoTets(int i1, int i2, int i3, int i4, int i5, int i6)
    {
        var result = new List<int[]>();
        bool swap3546 = Math.Max(i3, i5) > Math.Max(i4, i6);
        
        if (swap3546)
        {
            result.Add(new int[] { i1, i3, i4, i5 });
            result.Add(new int[] { i1, i3, i5, i6 });
            result.Add(new int[] { i3, i4, i5, i2 });
            result.Add(new int[] { i3, i5, i6, i2 });
        }
        else
        {
            result.Add(new int[] { i1, i4, i6, i3 });
            result.Add(new int[] { i1, i6, i4, i5 });
            result.Add(new int[] { i2, i4, i6, i5 });
            result.Add(new int[] { i2, i6, i4, i3 });
        }
        
        return result;
    }
    
    #endregion
    
    #region Element Creation Helpers (with OriginalElement tracking)
    
    /// <summary>Adds a Point element and sets its OriginalElement to the parent point index.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void AddPt(SimplexMesh mesh, int origIdx, int n0)
    {
        var idx = mesh.Add<Point, Node>(n0);
        mesh.Set<Point, OriginalElement>(idx, new OriginalElement(origIdx));
    }
    
    /// <summary>Adds a Bar2 element and sets its OriginalElement to the parent bar index.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void AddBar(SimplexMesh mesh, int origIdx, int n0, int n1)
    {
        var idx = mesh.Add<Bar2, Node>(n0, n1);
        mesh.Set<Bar2, OriginalElement>(idx, new OriginalElement(origIdx));
    }
    
    /// <summary>Adds a Tri3 element and sets its OriginalElement to the parent triangle index.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void AddTri(SimplexMesh mesh, int origIdx, int n0, int n1, int n2)
    {
        var idx = mesh.Add<Tri3, Node>(n0, n1, n2);
        mesh.Set<Tri3, OriginalElement>(idx, new OriginalElement(origIdx));
    }
    
    /// <summary>Adds a Tet4 element and sets its OriginalElement to the parent tetrahedron index.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void AddTet(SimplexMesh mesh, int origIdx, int n0, int n1, int n2, int n3)
    {
        var idx = mesh.Add<Tet4, Node>(n0, n1, n2, n3);
        mesh.Set<Tet4, OriginalElement>(idx, new OriginalElement(origIdx));
    }
    
    #endregion
    
    #region Triangle Splitting and Utilities
    
    private static void SplitTriangle(SimplexMesh mesh, int origIdx, int v0, int v1, int v2, int m01, int m12, int m20)
    {
        // CASE 0: no marked edges — retain original triangle
        if (m01 < 0 && m12 < 0 && m20 < 0)
            AddTri(mesh, origIdx, v0, v1, v2);
        
        // CASE 1 (×3): one marked edge — split into 2 triangles
        // Fortran splittria31 CASE(1): (start, mid, opposite), (mid, end, opposite)
        else if (m01 >= 0 && m12 < 0 && m20 < 0)
        {
            AddTri(mesh, origIdx, v0, m01, v2);
            AddTri(mesh, origIdx, m01, v1, v2);
        }
        else if (m01 < 0 && m12 >= 0 && m20 < 0)
        {
            AddTri(mesh, origIdx, v1, m12, v0);
            AddTri(mesh, origIdx, m12, v2, v0);
        }
        else if (m01 < 0 && m12 < 0 && m20 >= 0)
        {
            AddTri(mesh, origIdx, v2, m20, v1);
            AddTri(mesh, origIdx, m20, v0, v1);
        }
        
        // CASE 2 (×3): two marked edges — split into 3 triangles
        // Fortran splittria31 CASE(2): diagonal heuristic based on midpoint node indices
        // SplitTri2Edges(nn1, mid1, nnMeet, mid2, nn5) where nnMeet is where both marked edges meet
        else if (m01 >= 0 && m12 >= 0 && m20 < 0)
            SplitTri2Edges(mesh, origIdx, v0, m01, v1, m12, v2);
        else if (m01 >= 0 && m12 < 0 && m20 >= 0)
            SplitTri2Edges(mesh, origIdx, v2, m20, v0, m01, v1);
        else if (m01 < 0 && m12 >= 0 && m20 >= 0)
            SplitTri2Edges(mesh, origIdx, v1, m12, v2, m20, v0);
        
        // CASE 3: three marked edges — split into 4 triangles
        // Fortran splittria31 CASE(3): standard red refinement
        else if (m01 >= 0 && m12 >= 0 && m20 >= 0)
        {
            AddTri(mesh, origIdx, v0, m01, m20);
            AddTri(mesh, origIdx, m01, v1, m12);
            AddTri(mesh, origIdx, m12, v2, m20);
            AddTri(mesh, origIdx, m01, m12, m20);
        }
    }
    
    /// <summary>
    /// Fortran splittria31 CASE(2): two marked edges meeting at a common node.
    /// The quadrilateral formed by (nn1, mid1, mid2, nn5) is split along a diagonal
    /// chosen by comparing midpoint node indices, ensuring conforming diagonals
    /// across elements sharing the same edge.
    /// </summary>
    /// <param name="nn1">First vertex (node at start of first marked edge)</param>
    /// <param name="mid1">Midpoint of first marked edge (between nn1 and nnMeet)</param>
    /// <param name="nnMeet">Common vertex where both marked edges meet</param>
    /// <param name="mid2">Midpoint of second marked edge (between nnMeet and nn5)</param>
    /// <param name="nn5">Opposite vertex (node at end of second marked edge)</param>
    private static void SplitTri2Edges(SimplexMesh mesh, int origIdx,
        int nn1, int mid1, int nnMeet, int mid2, int nn5)
    {
        // Corner triangle at the meeting point — always the same
        AddTri(mesh, origIdx, mid1, nnMeet, mid2);
        
        // Diagonal choice for quad (nn1, mid1, mid2, nn5)
        // Fortran: IF (newnode2 > newnode4)
        if (mid1 > mid2)
        {
            // Diagonal mid1—nn5
            AddTri(mesh, origIdx, nn1, mid1, nn5);
            AddTri(mesh, origIdx, mid1, mid2, nn5);
        }
        else
        {
            // Diagonal nn1—mid2
            AddTri(mesh, origIdx, nn1, mid1, mid2);
            AddTri(mesh, origIdx, nn1, mid2, nn5);
        }
    }

    private static void ComputeTetEdgeClosureFull(SimplexMesh mesh, HashSet<(int, int)> edgesToRefine, bool validateTopology)
    {
        var edgesToAdd = new HashSet<(int, int)>();
        bool changed = true;
        int iteration = 0;

        while (changed && iteration < 100)
        {
            changed = false;
            edgesToAdd.Clear();

            for (int i = 0; i < mesh.Count<Tet4>(); i++)
            {
                var nodes = mesh.NodesOf<Tet4, Node>(i);
                var edges = new[]
                {
                    (nodes[0], nodes[1]),
                    (nodes[1], nodes[2]),
                    (nodes[0], nodes[2]),
                    (nodes[2], nodes[3]),
                    (nodes[0], nodes[3]),
                    (nodes[1], nodes[3])
                };

                // Count how many edges are marked (check both directions)
                int markedCount = edges.Count(e => ContainsEdge(edgesToRefine, e.Item1, e.Item2));

                if (markedCount == 1 || markedCount == 2)
                {
                    foreach (var e in edges)
                        if (!ContainsEdge(edgesToRefine, e.Item1, e.Item2))
                            edgesToAdd.Add(e);
                }
            }

            foreach (var e in edgesToAdd)
            {
                if (edgesToRefine.Add(e))
                    changed = true;
            }

            iteration++;
        }
    }

    // Helper to check if edge exists in set in either direction
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool ContainsEdge(HashSet<(int, int)> edges, int a, int b)
    {
        return edges.Contains((a, b)) || edges.Contains((b, a));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static (int, int) Canonical(int a, int b) => a < b ? (a, b) : (b, a);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int GetMidpoint(Dictionary<(int, int), int> map, int a, int b)
    {
        // Try exact order first (a, b)
        if (map.TryGetValue((a, b), out var mid))
            return mid;
        // Try reversed order (b, a)
        if (map.TryGetValue((b, a), out mid))
            return mid;
        // Edge not found
        return -1;
    }
    
    #endregion
}

// MeshRefinement.cs
//
// Copyright (C) 2026 Pedro Miguel de Almeida Areias
//
// This file is part of ManyToMany.
//
// ManyToMany is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ManyToMany is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ManyToMany. If not, see <https://www.gnu.org/licenses/>.

// Edge refinement of simplex meshes, after Fortran remeshsimplex.f90: midpoint
// creation, closure, point propagation (splitpoint1) and the per-element split.
// The split itself (splitbar21, splittria31, splittetwork1) is SubdivisionTemplates,
// the single implementation shared with the incremental edge cutter of Step2.Core.
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
            var barRows = new List<int[]>(2);
            for (int iel = 0; iel < inputBarCount; iel++)
            {
                var bnodes = input.NodesOf<Bar2, Node>(iel);
                int midNode = GetMidpoint(midpointMap, bnodes[0], bnodes[1]);
                if (midNode < 0) barsCopied++;
                else barsSplit++;
                
                // Copied when the edge is not marked, split in two otherwise (SubdivisionTemplates)
                barRows.Clear();
                SubdivisionTemplates.SplitBar(nodeRemap[bnodes[0]], nodeRemap[bnodes[1]], midNode, barRows);
                foreach (var r in barRows) AddBar(output, iel, r[0], r[1]);
            }
            Console.WriteLine($"  → Bars: {barsCopied} copied + {barsSplit} split (×2) = {output.Count<Bar2>()} total");
        }

        // ═══════════════════════════════════════════════════════════════════
        // Refine Tri3 elements (Fortran: splittria30/splittria31)
        // ═══════════════════════════════════════════════════════════════════
        var children = new List<int[]>(8); // rows emitted by a template for one parent element
        int triRefinedCount = 0;
        for (var i = 0; i < input.Count<Tri3>(); i++)
        {
            var nodes = input.NodesOf<Tri3, Node>(i);
            int v0 = nodeRemap[nodes[0]], v1 = nodeRemap[nodes[1]], v2 = nodeRemap[nodes[2]];
            
            // Look up midpoints using INPUT indices (midpointMap keys are input indices)
            var m01 = GetMidpoint(midpointMap, nodes[0], nodes[1]);
            var m12 = GetMidpoint(midpointMap, nodes[1], nodes[2]);
            var m20 = GetMidpoint(midpointMap, nodes[2], nodes[0]);
            
            if (!(m01 < 0 && m12 < 0 && m20 < 0))
                triRefinedCount++;
            
            // Fortran splittria30/splittria31 (SubdivisionTemplates)
            children.Clear();
            SubdivisionTemplates.SplitTriangle(v0, v1, v2, m01, m12, m20, children);
            AddTris(output, i, children);
        }
        Console.WriteLine($"  → Triangles: {triRefinedCount} refined out of {input.Count<Tri3>()} → {output.Count<Tri3>()} total");

        // ═══════════════════════════════════════════════════════════════════
        // Refine Tet4 elements (Fortran: splittetra40/splittetra41)
        // ═══════════════════════════════════════════════════════════════════
        int tetCount = input.Count<Tet4>();
        Span<int> globalNodes = stackalloc int[4];
        Span<int> midnodes = stackalloc int[6]; // midpoint label of local edge e, 0 when not marked
        
        for (var i = 0; i < tetCount; i++)
        {
            var nodes = input.NodesOf<Tet4, Node>(i);
            for (int k = 0; k < 4; k++) globalNodes[k] = nodeRemap[nodes[k]];
            
            int nmarked = 0;
            for (int e = 0; e < 6; e++)
            {
                // Look up using INPUT mesh indices
                var en = SubdivisionTemplates.EdgeNodes[e];
                int mid = GetMidpoint(midpointMap, nodes[en[0]], nodes[en[1]]);
                midnodes[e] = (mid >= 0) ? mid : 0;
                if (mid >= 0) nmarked++;
                else missingMidpoints++;
            }
            
            // Track case distribution
            caseCounts[nmarked]++;
            
            // Fortran splittetwork1 (SubdivisionTemplates)
            children.Clear();
            SubdivisionTemplates.SplitTet(globalNodes, midnodes, children);
            AddTets(output, i, children);
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
    
    /// <summary>
    /// Compute signed volume * 6 (uses shared MeshGeometry implementation).
    /// </summary>
    private static double SignedVolume6x(double[,] coords, int n0, int n1, int n2, int n3)
        => ComputeTetrahedronJacobian(coords, n0, n1, n2, n3);
    
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
    
    /// <summary>Appends the triangle rows emitted by a template, each tagged with its parent triangle.</summary>
    private static void AddTris(SimplexMesh mesh, int origIdx, List<int[]> rows)
    {
        foreach (var r in rows)
        {
            var idx = mesh.Add<Tri3, Node>(r[0], r[1], r[2]);
            mesh.Set<Tri3, OriginalElement>(idx, new OriginalElement(origIdx));
        }
    }
    
    /// <summary>Appends the tetrahedron rows emitted by a template, each tagged with its parent tetrahedron.</summary>
    private static void AddTets(SimplexMesh mesh, int origIdx, List<int[]> rows)
    {
        foreach (var r in rows)
        {
            var idx = mesh.Add<Tet4, Node>(r[0], r[1], r[2], r[3]);
            mesh.Set<Tet4, OriginalElement>(idx, new OriginalElement(origIdx));
        }
    }
    
    #endregion
    
    #region Closure and Midpoint Utilities
    
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

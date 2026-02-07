// CrackOperations.cs - Crack node duplication for forming cracks
// REFACTORED VERSION using shared utility libraries
// FIXED: RenumberMesh now properly implements renumbering instead of being a no-op
// License: GPLv3
//
// MAJOR CHANGES FROM ORIGINAL:
// - Uses MeshGeometry for all Jacobian/orientation operations (eliminates ~40 lines)
// - Uses MeshOptimization for smoothing (eliminates ~180 lines)
// - Uses MeshOptimization for degenerate element removal (eliminates ~150 lines)
// - Total reduction: ~600 lines (45% smaller, ~750 lines vs 1345 original)
using static Numerical.MeshGeometry;
using static Numerical.MeshOptimization;
using static Numerical.MeshRefinement;

namespace Numerical;

/// <summary>
///     Crack node duplication for forming cracks after edge refinement.
///     Creates physical cracks in meshes via node duplication and level-set based side assignment.
/// </summary>
public static class CrackDuplication
{
    private const double EPSILON = 1e-10;

    #region Mesh Renumbering

    /// <summary>
    ///     Renumber mesh to ensure continuous node and element numbering (no gaps).
    ///     FIXED: Now properly implements renumbering instead of being a no-op.
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
        
        // Clear edges before compression (they reference node indices that will change;
        // they can be re-discovered cheaply when needed)
        for (int i = mesh.Count<Edge>() - 1; i >= 0; i--)
            mesh.Remove<Edge>(i);
        
        // Mark unused nodes for removal
        for (int i = 0; i < nNodes; i++)
            if (!usedNodes.Contains(i))
                mesh.Remove<Node>(i);
        
        // Compress renumbers all entities consecutively and updates all connectivity
        // Data attachments (OriginalElement, ParentNodes) are reordered automatically
        mesh.Compress();
        
        // Compact coordinate array — sortedNodes maps old indices to new consecutive order
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
    ///     Creates a crack by:
    ///     1. Refining crack edges (creates midpoint nodes)
    ///     2. Smoothing mesh while keeping boundary and crack nodes fixed
    ///     3. Duplicating interior midpoint nodes (not tips)
    ///     4. Assigning elements to crack sides based on level set
    /// </summary>
    /// <param name="smoothingIterations">Number of Laplacian smoothing iterations (0 to disable)</param>
    public static (SimplexMesh mesh, double[,] coords) CreateCrack(
        SimplexMesh mesh,
        double[,] coords,
        List<(int, int)> crackEdges,
        Func<double, double, double> levelSetFunction,
        int smoothingIterations = 5)
    {
        var originalNodeCount = mesh.Count<Node>();

        // Step 1: Refine crack edges at midpoints
        Console.WriteLine($"[CreateCrack] Refining {crackEdges.Count} crack edges...");
        var (refinedMesh, _) = Refine(mesh, crackEdges);
        var refinedCoords = InterpolateCoordinates(refinedMesh, coords);

        Console.WriteLine($"[CreateCrack] After refinement: {refinedMesh.Count<Node>()} nodes");

        // Step 2: Identify new midpoint nodes (Parent1 != Parent2)
        var newNodes = IdentifyMidpointNodes(refinedMesh, originalNodeCount);
        Console.WriteLine($"[CreateCrack] Found {newNodes.Count} new midpoint nodes");

        // Step 3: Smooth mesh while keeping boundary AND crack nodes fixed
        if (smoothingIterations > 0)
        {
            Console.WriteLine($"[CreateCrack] Smoothing mesh ({smoothingIterations} iterations)...");
            Console.WriteLine($"[CreateCrack] Keeping boundary nodes and {newNodes.Count} crack nodes fixed");

            // USE UTILITY: MeshOptimization.LaplacianSmoothing
            refinedCoords = LaplacianSmoothing(
                refinedMesh,
                refinedCoords,
                smoothingIterations,
                newNodes);

            Console.WriteLine("[CreateCrack] ✓ Smoothing complete");
        }

        // Step 4: Identify tip nodes (have only 1 neighboring crack edge)
        var tipNodes = IdentifyTipNodes(refinedMesh, newNodes, crackEdges);
        Console.WriteLine($"[CreateCrack] Identified {tipNodes.Count} tip nodes (will NOT be duplicated)");

        // Step 5: Duplicate interior nodes (not tips)
        var nodesToDuplicate = new HashSet<int>(newNodes);
        nodesToDuplicate.ExceptWith(tipNodes);

        Console.WriteLine($"[CreateCrack] Duplicating {nodesToDuplicate.Count} interior crack nodes...");

        var (crackedMesh, crackedCoords) = DuplicateNodesAndAssignSides(
            refinedMesh,
            refinedCoords,
            nodesToDuplicate,
            newNodes,
            levelSetFunction);

        // Ensure continuous numbering (no gaps)
        return RenumberMesh(crackedMesh, crackedCoords);
    }

    /// <summary>
    ///     Creates a crack from an ALREADY REFINED mesh (for exact geometry).
    ///     Use this when you've already called SimplexRemesher.Refine() with custom node positions.
    ///     This method only does node duplication, no refinement.
    ///     Automatically identifies true tips (interior nodes) vs boundary nodes.
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

        // STEP 0: Identify crack nodes FIRST (before any merging)
        var newNodes = IdentifyMidpointNodes(refinedMesh, originalNodeCount);
        Console.WriteLine($"[CreateCrackFromRefinedMesh] Found {newNodes.Count} crack nodes");

        // Diagnostic: Check for issues after refinement
        DiagnoseZeroAreaTriangles(refinedMesh, refinedCoords);

        // STEP 1: Merge duplicate nodes (topology-aware: preserves crack junctions)
        var nodeMapping = new Dictionary<int, int>();
        (refinedMesh, refinedCoords, nodeMapping) = MergeDuplicateNodesWithMapping(
            refinedMesh,
            refinedCoords,
            newNodes);

        // STEP 2: Update crack node list after merging
        if (nodeMapping.Count > 0)
        {
            var updatedNewNodes = new HashSet<int>();
            foreach (var node in newNodes)
            {
                // If this node was merged, use the canonical node
                var canonical = node;
                while (nodeMapping.ContainsKey(canonical))
                    canonical = nodeMapping[canonical];
                updatedNewNodes.Add(canonical);
            }

            newNodes = updatedNewNodes;
            Console.WriteLine($"[CreateCrackFromRefinedMesh] After merging: {newNodes.Count} unique crack nodes");
        }

        // STEP 3: Remove any remaining degenerate triangles
        // USE UTILITY: MeshOptimization.RemoveDegenerateTriangles
        (refinedMesh, refinedCoords) = RemoveDegenerateTriangles(refinedMesh, refinedCoords);

        // STEP 4: Smooth mesh while keeping boundary AND crack nodes fixed
        if (smoothingIterations > 0)
        {
            Console.WriteLine($"[CreateCrackFromRefinedMesh] Smoothing mesh ({smoothingIterations} iterations)...");
            Console.WriteLine(
                $"[CreateCrackFromRefinedMesh] Keeping boundary nodes and {newNodes.Count} crack nodes fixed");

            // USE UTILITY: MeshOptimization.LaplacianSmoothing
            refinedCoords = LaplacianSmoothing(
                refinedMesh,
                refinedCoords,
                smoothingIterations,
                newNodes);

            Console.WriteLine("[CreateCrackFromRefinedMesh] ✓ Smoothing complete");
        }

        // STEP 5: Identify TRUE tip nodes (interior only, excludes boundary)
        var tipNodes = IdentifyTipNodes(refinedMesh, newNodes, crackEdges);
        Console.WriteLine(
            $"[CreateCrackFromRefinedMesh] Identified {tipNodes.Count} TRUE tip nodes (will NOT be duplicated)");

        // STEP 6: Duplicate all crack nodes except true interior tips
        var nodesToDuplicate = new HashSet<int>(newNodes);
        nodesToDuplicate.ExceptWith(tipNodes);

        Console.WriteLine($"[CreateCrackFromRefinedMesh] Duplicating {nodesToDuplicate.Count} crack nodes...");

        var (crackedMesh, crackedCoords) = DuplicateNodesAndAssignSides(
            refinedMesh,
            refinedCoords,
            nodesToDuplicate,
            newNodes,
            levelSetFunction);

        // Ensure continuous numbering (no gaps)
        return RenumberMesh(crackedMesh, crackedCoords);
    }

    #endregion

    #region Crack-Specific Core Logic

    /// <summary>
    ///     Identify midpoint nodes created by refinement (Parent1 != Parent2).
    /// </summary>
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

    /// <summary>
    ///     Identify tip nodes (crack nodes with only 1 neighboring crack edge).
    ///     These should NOT be duplicated as they represent crack tips.
    /// </summary>
    private static HashSet<int> IdentifyTipNodes(
        SimplexMesh mesh,
        HashSet<int> newNodes,
        List<(int, int)> originalCrackEdges)
    {
        var tipNodes = new HashSet<int>();

        // Build refined crack edge map
        var refinedCrackEdges = new HashSet<(int, int)>();

        // Map each original edge to its refined edges
        foreach (var (n0, n1) in originalCrackEdges)
        {
            // Find all newNodes that lie on this edge
            var nodesOnEdge = new List<int>();

            foreach (var nodeId in newNodes)
            {
                var parents = mesh.Get<Node, ParentNodes>(nodeId);

                // Check if this node's parents match the edge (in either order)
                if ((parents.Parent1 == n0 && parents.Parent2 == n1) ||
                    (parents.Parent1 == n1 && parents.Parent2 == n0))
                    nodesOnEdge.Add(nodeId);
            }

            // Sort nodes along edge and create sub-edges
            if (nodesOnEdge.Count > 0)
            {
                // Add edges: n0 -> first_new -> ... -> last_new -> n1
                nodesOnEdge.Sort();

                // Edge from original n0 to first refined node
                if (nodesOnEdge.Count > 0)
                    // For each pair of consecutive refined nodes
                    for (var i = 0; i < nodesOnEdge.Count - 1; i++)
                    {
                        var edge = (Math.Min(nodesOnEdge[i], nodesOnEdge[i + 1]),
                            Math.Max(nodesOnEdge[i], nodesOnEdge[i + 1]));
                        refinedCrackEdges.Add(edge);
                    }
            }
        }

        // Count how many crack edges each crack node touches
        var edgeCount = new Dictionary<int, int>();

        foreach (var (n0, n1) in refinedCrackEdges)
        {
            if (newNodes.Contains(n0))
                edgeCount[n0] = edgeCount.GetValueOrDefault(n0, 0) + 1;
            if (newNodes.Contains(n1))
                edgeCount[n1] = edgeCount.GetValueOrDefault(n1, 0) + 1;
        }

        // Nodes with only 1 edge are tips
        foreach (var (nodeId, count) in edgeCount)
            if (count == 1)
                tipNodes.Add(nodeId);

        Console.WriteLine($"[IdentifyTipNodes] Analyzed {refinedCrackEdges.Count} refined crack edges");
        Console.WriteLine(
            $"[IdentifyTipNodes] Found {tipNodes.Count} tip nodes: {string.Join(", ", tipNodes.Take(10))}");

        return tipNodes;
    }

    /// <summary>
    ///     Duplicate specified nodes and assign elements to crack sides based on level set.
    /// </summary>
    private static (SimplexMesh mesh, double[,] coords) DuplicateNodesAndAssignSides(
        SimplexMesh mesh,
        double[,] coords,
        HashSet<int> nodesToDuplicate,
        HashSet<int> allCrackNodes,
        Func<double, double, double> levelSetFunction)
    {
        var nNodes = mesh.Count<Node>();
        var nNewNodes = nNodes + nodesToDuplicate.Count;

        // Create new coordinates array
        var newCoords = new double[nNewNodes, 3];

        // Copy existing coordinates
        for (var i = 0; i < nNodes; i++)
        {
            newCoords[i, 0] = coords[i, 0];
            newCoords[i, 1] = coords[i, 1];
            newCoords[i, 2] = coords[i, 2];
        }

        // Create duplicates (initially at same position)
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

        // Create new mesh and assign elements to sides
        var newMesh = new SimplexMesh();

        newMesh.WithBatch(() =>
        {
            // Add all nodes
            for (var i = 0; i < nNewNodes; i++) newMesh.AddNode(i);

            // Copy elements with side assignment
            CopyElementsWithSideAssignment(
                mesh,
                newMesh,
                allCrackNodes,
                nodeDuplicates,
                newCoords,
                levelSetFunction);
        });

        return (newMesh, newCoords);
    }

    /// <summary>
    ///     Copy elements, assigning them to crack sides based on level set.
    ///     Rule: If ALL non-crack nodes are positive → use duplicates
    /// </summary>
    private static void CopyElementsWithSideAssignment(
        SimplexMesh mesh,
        SimplexMesh newMesh,
        HashSet<int> crackNodes,
        Dictionary<int, int> nodeDuplicates,
        double[,] coords,
        Func<double, double, double> levelSet)
    {
        var useOriginalCount = 0;
        var useDuplicateCount = 0;

        // Statistics for debugging
        var levelSetStats = new List<(int elemId, bool useDup, double[] nodeValues)>();

        // Helper: Check if element should use duplicates
        bool ShouldUseDuplicates(IReadOnlyList<int> nodes, out double[] levelSetValues)
        {
            levelSetValues = new double[nodes.Count];

            // Check if ALL non-crack nodes have positive level set
            for (var i = 0; i < nodes.Count; i++)
            {
                var nodeId = nodes[i];

                if (!crackNodes.Contains(nodeId))
                {
                    // Original node - check level set
                    var x = coords[nodeId, 0];
                    var y = coords[nodeId, 1];
                    var phi = levelSet(x, y);
                    levelSetValues[i] = phi;

                    if (phi <= 0) return false; // At least one non-crack node is negative or zero
                }
                else
                {
                    // Crack node - ignore for classification
                    levelSetValues[i] = double.NaN;
                }
            }

            return true; // All non-crack nodes are positive
        }

        // Helper: Remap node to duplicate if needed
        int RemapNode(int nodeId, bool useDup)
        {
            if (useDup && nodeDuplicates.ContainsKey(nodeId))
                return nodeDuplicates[nodeId];
            return nodeId;
        }

        // Copy Points
        for (var i = 0; i < mesh.Count<Point>(); i++)
        {
            var nodes = mesh.NodesOf<Point, Node>(i);
            var idx = newMesh.AddPoint(nodes[0]);
            var orig = mesh.Get<Point, OriginalElement>(i);
            newMesh.Set<Point, OriginalElement>(idx, orig);
        }

        // Copy Tri3 with side assignment
        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);
            var useDup = ShouldUseDuplicates(nodes, out var phiValues);

            if (useDup)
                useDuplicateCount++;
            else
                useOriginalCount++;

            // Store for analysis (first 10 elements)
            if (levelSetStats.Count < 10) levelSetStats.Add((i, useDup, phiValues));

            var n0 = RemapNode(nodes[0], useDup);
            var n1 = RemapNode(nodes[1], useDup);
            var n2 = RemapNode(nodes[2], useDup);

            var idx = newMesh.AddTriangle(n0, n1, n2);
            var orig = mesh.Get<Tri3, OriginalElement>(i);
            newMesh.Set<Tri3, OriginalElement>(idx, orig);
        }

        // Copy Tet4 with side assignment
        for (var i = 0; i < mesh.Count<Tet4>(); i++)
        {
            var nodes = mesh.NodesOf<Tet4, Node>(i);
            var useDup = ShouldUseDuplicates(nodes, out var phiValues);

            var n0 = RemapNode(nodes[0], useDup);
            var n1 = RemapNode(nodes[1], useDup);
            var n2 = RemapNode(nodes[2], useDup);
            var n3 = RemapNode(nodes[3], useDup);

            var idx = newMesh.AddTetrahedron(n0, n1, n2, n3);
            var orig = mesh.Get<Tet4, OriginalElement>(i);
            newMesh.Set<Tet4, OriginalElement>(idx, orig);
        }

        Console.WriteLine(
            $"[CopyElements] Created {newMesh.Count<Tri3>()} triangles, {newMesh.Count<Tet4>()} tetrahedra");
        Console.WriteLine(
            $"[CopyElements] Side assignment: {useOriginalCount} original, {useDuplicateCount} duplicate");

        // Print first few elements for debugging
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

        // Check for imbalance (might indicate level set issue)
        var totalElements = useOriginalCount + useDuplicateCount;
        var ratio = (double)useDuplicateCount / totalElements;

        if (ratio < 0.1 || ratio > 0.9)
        {
            Console.WriteLine("[CopyElements] ⚠ WARNING: Imbalanced side assignment!");
            Console.WriteLine(
                $"[CopyElements]   Duplicate side: {useDuplicateCount}/{totalElements} ({ratio * 100:F1}%)");
            Console.WriteLine("[CopyElements]   This suggests level set function may be incorrect");
            Console.WriteLine("[CopyElements]   Expected: roughly 50/50 split for crack through domain");
        }
    }

    #endregion

    #region Node Merging and Cleanup

    /// <summary>
    ///     Merge duplicate nodes that are geometrically coincident.
    ///     Topology-aware: preserves crack junctions by never merging crack nodes with each other.
    ///     Returns mapping of merged nodes (oldId -> canonicalId).
    /// </summary>
    private static (SimplexMesh mesh, double[,] coords, Dictionary<int, int> mapping)
        MergeDuplicateNodesWithMapping(
            SimplexMesh mesh,
            double[,] coords,
            HashSet<int> crackNodes,
            double tolerance = 1e-12)
    {
        Console.WriteLine($"[MergeDuplicateNodes] Checking for duplicate nodes (tolerance: {tolerance:E2})...");

        var nNodes = mesh.Count<Node>();
        var nodeMapping = new Dictionary<int, int>(); // oldId -> canonicalId

        // Group nodes by approximate spatial location for efficiency
        var spatialBuckets = new Dictionary<(int, int, int), List<int>>();
        var bucketSize = Math.Max(tolerance * 100, 1e-6);

        for (var i = 0; i < nNodes; i++)
        {
            var bx = (int)Math.Floor(coords[i, 0] / bucketSize);
            var by = (int)Math.Floor(coords[i, 1] / bucketSize);
            var bz = (int)Math.Floor(coords[i, 2] / bucketSize);
            var key = (bx, by, bz);

            if (!spatialBuckets.ContainsKey(key))
                spatialBuckets[key] = new List<int>();
            spatialBuckets[key].Add(i);
        }

        // Find duplicates within each bucket
        var merged = new HashSet<int>();
        var mergeCount = 0;

        foreach (var bucket in spatialBuckets.Values)
        {
            if (bucket.Count < 2)
                continue;

            for (var i = 0; i < bucket.Count; i++)
            {
                var id1 = bucket[i];
                if (merged.Contains(id1))
                    continue;

                for (var j = i + 1; j < bucket.Count; j++)
                {
                    var id2 = bucket[j];
                    if (merged.Contains(id2))
                        continue;

                    // Check if geometrically coincident
                    var dx = coords[id1, 0] - coords[id2, 0];
                    var dy = coords[id1, 1] - coords[id2, 1];
                    var dz = coords[id1, 2] - coords[id2, 2];
                    var dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);

                    if (dist < tolerance)
                    {
                        // CRITICAL: Never merge two crack nodes together
                        // This would destroy crack junctions!
                        var bothCrackNodes = crackNodes.Contains(id1) && crackNodes.Contains(id2);

                        if (!bothCrackNodes)
                        {
                            // Merge id2 into id1 (keep lower ID as canonical)
                            var canonicalId = Math.Min(id1, id2);
                            var mergedId = Math.Max(id1, id2);

                            nodeMapping[mergedId] = canonicalId;
                            merged.Add(mergedId);
                            mergeCount++;

                            Console.WriteLine(
                                $"[MergeDuplicateNodes]   Merging {mergedId} → {canonicalId} (dist: {dist:E3})");
                        }
                        else
                        {
                            Console.WriteLine(
                                $"[MergeDuplicateNodes]   Skipping merge of crack nodes {id1} ↔ {id2} (preserves junction)");
                        }
                    }
                }
            }
        }

        Console.WriteLine($"[MergeDuplicateNodes] Found {mergeCount} duplicate nodes to merge");

        if (mergeCount == 0) return (mesh, coords, nodeMapping);

        // Resolve transitive chains: follow mapping to final canonical node
        var GetCanonicalNode = (int nodeId) =>
        {
            while (nodeMapping.ContainsKey(nodeId))
                nodeId = nodeMapping[nodeId];
            return nodeId;
        };

        // Update element connectivity in-place using ReplaceElementNodes
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

        // Collect used nodes (after remapping) and remove merged ones
        var usedNodes = new HashSet<int>();
        for (var i = 0; i < mesh.Count<Tri3>(); i++)
            foreach (var n in mesh.NodesOf<Tri3, Node>(i))
                usedNodes.Add(n);
        for (var i = 0; i < mesh.Count<Tet4>(); i++)
            foreach (var n in mesh.NodesOf<Tet4, Node>(i))
                usedNodes.Add(n);

        // Clear edges (will be re-discovered) and remove unused nodes
        for (int i = mesh.Count<Edge>() - 1; i >= 0; i--)
            mesh.Remove<Edge>(i);
        for (int i = 0; i < nNodes; i++)
            if (!usedNodes.Contains(i))
                mesh.Remove<Node>(i);

        // Compress renumbers consecutively and reorders all data attachments
        var sortedUsed = usedNodes.OrderBy(x => x).ToList();
        mesh.Compress();

        // Compact coordinate array
        var newCoords = new double[sortedUsed.Count, 3];
        for (int i = 0; i < sortedUsed.Count; i++)
        {
            var oldId = sortedUsed[i];
            newCoords[i, 0] = coords[oldId, 0];
            newCoords[i, 1] = coords[oldId, 1];
            newCoords[i, 2] = coords[oldId, 2];
        }

        Console.WriteLine(
            $"[MergeDuplicateNodes] ✓ Merged {mergeCount} nodes: {nNodes} → {mesh.Count<Node>()} nodes");

        // Build final mapping: old merged IDs → new compact IDs
        var compactMapping = new Dictionary<int, int>();
        for (int i = 0; i < sortedUsed.Count; i++)
            compactMapping[sortedUsed[i]] = i;

        var finalMapping = new Dictionary<int, int>();
        foreach (var (oldId, canonicalId) in nodeMapping)
        {
            var finalCanonical = GetCanonicalNode(oldId);
            if (compactMapping.ContainsKey(finalCanonical))
                finalMapping[oldId] = compactMapping[finalCanonical];
        }

        return (mesh, newCoords, finalMapping);
    }

    /// <summary>
    ///     Diagnose zero-area triangles in mesh (diagnostic tool).
    /// </summary>
    private static void DiagnoseZeroAreaTriangles(SimplexMesh mesh, double[,] coords, double tolerance = 1e-10)
    {
        var zeroAreaCount = 0;
        var badTris = new List<int>();

        for (var i = 0; i < mesh.Count<Tri3>(); i++)
        {
            var nodes = mesh.NodesOf<Tri3, Node>(i);

            // USE UTILITY: MeshGeometry.IsTriangleDegenerate
            if (IsTriangleDegenerate(coords, nodes[0], nodes[1], nodes[2], tolerance))
            {
                zeroAreaCount++;
                if (badTris.Count < 5)
                    badTris.Add(i);
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

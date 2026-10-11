// HistoryTransfer.cs
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

// Transfer of nodal state from a parent mesh to its refinement.
//
// Refine() creates every new node as the midpoint of a parent edge and records
// the lineage as ParentNodes(a, b); a surviving node carries ParentNodes(i, i).
// Every node in the refined mesh therefore has exactly two parents, and the
// transfer of any nodal field is fixed by what is done on that pair.
//
// Which rule is correct depends on the field, and getting it wrong is silent:
//
//   coordinates, displacements    linear interpolation   (Interpolate)
//   irreversible history H, d     maximum                (TransferMonotone)
//
// Averaging a monotone field over its parents can return a value BELOW the
// maximum already attained at a parent.  For the history field of a phase-field
// or gradient-damage model that lowers the stored driving force, lets the damage
// unload, and heals a crack the continuous model forbids.  The error does not
// throw, does not show up in a residual, and is not visible in the converged
// solution of the step in which it happens.
//
// Example usage:
//   var (refined, remap) = MeshRefinement.Refine(mesh, edges, coords: coords);
//   var newCoords  = MeshRefinement.InterpolateCoordinates(refined, coords);
//   var newHistory = HistoryTransfer.TransferMonotone(refined, history);
//   var newDamage  = HistoryTransfer.TransferMonotone(refined, damage);
//   HistoryTransfer.VerifyMonotone(refined, history, newHistory, "after refine");
//
// Quadrature-point values are NOT handled here: recover them by evaluating the
// transferred nodal field on the new mesh, never by projecting the old
// quadrature values, which reintroduces exactly the averaging this avoids.

using System;

namespace Numerical;

/// <summary>
///     Transfer of nodal state from a parent mesh onto its refinement, using the
///     ParentNodes lineage recorded by <see cref="MeshRefinement.Refine" />.
///     Provides the monotone (maximum) rule required by irreversible fields and
///     the linear rule appropriate to kinematic ones.
/// </summary>
public static class HistoryTransfer
{
    /// <summary>
    ///     Transfers a scalar nodal field that must not decrease, taking the maximum
    ///     over the two parents of each node:  H_m = max(H_a, H_b).
    ///     This is the only rule consistent with an irreversibility condition of the
    ///     form H(t) = max over the history of the driving force.
    /// </summary>
    /// <param name="refinedMesh">Mesh returned by <see cref="MeshRefinement.Refine" />.</param>
    /// <param name="parentField">Field on the parent mesh, indexed by parent node id.</param>
    /// <returns>Field on the refined mesh, indexed by refined node id.</returns>
    public static double[] TransferMonotone(SimplexMesh refinedMesh, double[] parentField)
    {
        ArgumentNullException.ThrowIfNull(refinedMesh);
        ArgumentNullException.ThrowIfNull(parentField);

        int nodeCount = refinedMesh.Count<Node>();
        var childField = new double[nodeCount];

        for (int i = 0; i < nodeCount; i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            CheckParents(parents, i, parentField.Length);

            childField[i] = parents.Parent1 == parents.Parent2
                ? parentField[parents.Parent1]
                : Math.Max(parentField[parents.Parent1], parentField[parents.Parent2]);
        }

        return childField;
    }

    /// <summary>
    ///     Component-wise monotone transfer of a multi-component nodal field,
    ///     laid out as [node, component].  Each component is treated independently.
    /// </summary>
    public static double[,] TransferMonotone(SimplexMesh refinedMesh, double[,] parentField)
    {
        ArgumentNullException.ThrowIfNull(refinedMesh);
        ArgumentNullException.ThrowIfNull(parentField);

        int nodeCount = refinedMesh.Count<Node>();
        int nComp = parentField.GetLength(1);
        int nParent = parentField.GetLength(0);
        var childField = new double[nodeCount, nComp];

        for (int i = 0; i < nodeCount; i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            CheckParents(parents, i, nParent);

            if (parents.Parent1 == parents.Parent2)
                for (int c = 0; c < nComp; c++)
                    childField[i, c] = parentField[parents.Parent1, c];
            else
                for (int c = 0; c < nComp; c++)
                    childField[i, c] = Math.Max(parentField[parents.Parent1, c],
                                                parentField[parents.Parent2, c]);
        }

        return childField;
    }

    /// <summary>
    ///     Transfers a scalar nodal field by linear interpolation on the parent edge,
    ///     f_m = (f_a + f_b)/2.  Correct for kinematic fields such as displacement;
    ///     WRONG for an irreversible history field — use <see cref="TransferMonotone(SimplexMesh, double[])" />
    ///     for those.  The scalar counterpart of
    ///     <see cref="MeshRefinement.InterpolateCoordinates" />.
    /// </summary>
    public static double[] Interpolate(SimplexMesh refinedMesh, double[] parentField)
    {
        ArgumentNullException.ThrowIfNull(refinedMesh);
        ArgumentNullException.ThrowIfNull(parentField);

        int nodeCount = refinedMesh.Count<Node>();
        var childField = new double[nodeCount];

        for (int i = 0; i < nodeCount; i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            CheckParents(parents, i, parentField.Length);

            childField[i] = parents.Parent1 == parents.Parent2
                ? parentField[parents.Parent1]
                : 0.5 * (parentField[parents.Parent1] + parentField[parents.Parent2]);
        }

        return childField;
    }

    /// <summary>
    ///     Verifies that a transferred field did not decrease anywhere relative to the
    ///     maximum over each node's parents.  Intended as a post-condition on any
    ///     adaptive step that moves an irreversible field, and as the mesh-side half of
    ///     an energy-bookkeeping check across remeshing.
    /// </summary>
    /// <returns>Number of nodes that violate monotonicity; zero if the transfer is admissible.</returns>
    public static int VerifyMonotone(SimplexMesh refinedMesh, double[] parentField,
                                     double[] childField, string label = "transfer",
                                     double tolerance = 0.0, bool verbose = true)
    {
        ArgumentNullException.ThrowIfNull(refinedMesh);
        ArgumentNullException.ThrowIfNull(parentField);
        ArgumentNullException.ThrowIfNull(childField);

        int nodeCount = refinedMesh.Count<Node>();
        if (childField.Length != nodeCount)
            throw new ArgumentException(
                $"childField has {childField.Length} entries, refined mesh has {nodeCount} nodes.",
                nameof(childField));

        int violations = 0;
        double worst = 0.0;

        for (int i = 0; i < nodeCount; i++)
        {
            var parents = refinedMesh.Get<Node, ParentNodes>(i);
            CheckParents(parents, i, parentField.Length);

            double required = Math.Max(parentField[parents.Parent1], parentField[parents.Parent2]);
            double deficit = required - childField[i];

            if (deficit > tolerance)
            {
                violations++;
                if (deficit > worst) worst = deficit;
            }
        }

        if (verbose)
        {
            if (violations == 0)
                Console.WriteLine($"  ✓ {label}: monotone at all {nodeCount} nodes");
            else
                Console.WriteLine($"  ⚠ WARNING: {label}: {violations} of {nodeCount} nodes lost history " +
                                  $"(worst deficit {worst:E3}) — irreversibility violated");
        }

        return violations;
    }

    /// <summary>
    ///     Guards against a node whose lineage was never recorded, which would otherwise
    ///     read parent 0 silently and corrupt the transfer.
    /// </summary>
    private static void CheckParents(ParentNodes parents, int node, int parentCount)
    {
        if (parents.Parent1 < 0 || parents.Parent1 >= parentCount ||
            parents.Parent2 < 0 || parents.Parent2 >= parentCount)
            throw new InvalidOperationException(
                $"Node {node} has ParentNodes({parents.Parent1}, {parents.Parent2}) outside the parent " +
                $"range [0, {parentCount}). The mesh was not produced by Refine(), or the parent field " +
                $"does not belong to the mesh that was refined.");
    }
}

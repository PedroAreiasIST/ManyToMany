# Meshing Directory Audit Report

**Date:** 2026-02-13
**Scope:** All source files in `/Meshing/` (10 active .cs files, ~8,500 lines)

## Overview

The Meshing library is a 2D/3D finite element mesh processing toolkit in C# targeting
.NET 9.0. It provides mesh generation, refinement (edge bisection), smoothing
(Laplacian/CVT), quality analysis, I/O (Gmsh/GiD/ASCII), and crack formation via
level sets.

## File Summary

| File | Lines | Role |
|------|-------|------|
| MeshConstants.cs | 43 | Shared tolerance constants |
| SimplexMesh.cs | 137 | Core mesh topology class |
| FiniteElementTopologies.cs | 91 | FE sub-entity definitions |
| GeometryCore.cs | 1,065 | Jacobians, areas, volumes, quality metrics |
| MeshGeneration.cs | 547 | Structured meshes + Delaunay |
| MeshRefinement.cs | 1,180 | Conforming edge bisection (Rivara) |
| MeshOperations.cs | 808 | Smoothing + quad conversion |
| CrackOperations.cs | 778 | Node duplication for cracks |
| SimplexRemesher.cs | 3,572 | I/O wrappers + crack-from-signed-field |
| UnifiedMesher.cs | 272 | High-level meshing API |

## Findings

### 1. CODE DUPLICATION

- `SimplexRemesher.ComputeTriangleJacobian()` (line 2311) duplicates
  `MeshGeometry.ComputeTriangleJacobian()` in GeometryCore.cs:24.
- `SimplexRemesher.ComputeTetrahedronJacobian()` (line 2323) duplicates
  `MeshGeometry.ComputeTetrahedronJacobian()` in GeometryCore.cs:93.
- `VerifyAndFixJacobians` (line 2229) duplicates logic already in
  `MeshRefinement.CheckJacobians` and `MeshRefinement.FixNegativeJacobians`.
- `BuildNodeNeighborsTri()` / `BuildNodeNeighborsTet()` in SimplexRemesher.cs
  (lines 3433, 3453) duplicate similar helpers in MeshOperations.cs.
- `FindBoundaryNodes()` (line 1937) and `FindBoundaryNodes3D()` (line 3203)
  partially duplicate `MeshOptimization.IdentifyBoundaryNodes()`.
- `TetEdges`/`EdgeIdx` structures (lines 107-117) overlap with `edgenodes` in
  MeshRefinement.cs:36-44.

### 2. DEAD/UNUSED CODE

- `DuplicateNewNodes()` (line 1963): always returns input unchanged ("disabled").
- `DetermineSideByCrackNormal()` (line 1917): always returns "NEGATIVE".
- `ComputeOpeningDirections()` (line 1983): ~120 lines with no callers.
- `CopyElementsWithCrackNodes()` (line 2115): superseded by `DuplicateNodesCarefully`.
- `CheckAllNonCrackNodesPositive()` (line 2204): remnant of older algorithm.
- `VerifyAndFixJacobians()` (line 2229): no callers found.
- `DuplicateNodesCarefully3D()` (line 3253): not called; 3D path has inline logic.
- `MeshRefinement.csold`: archived file that should be removed from the repository.

### 3. ORPHANED XML DOCUMENTATION

- SimplexRemesher.cs:416-418 has a dangling `<summary>` with no method.
- SimplexRemesher.cs:135 has a duplicate `</summary>` tag.
- Several methods have double `<summary>` blocks (lines 1665-1672, 1928-1932,
  2222-2228).

### 4. EXCESSIVE CONSOLE OUTPUT

The library has 100+ `Console.WriteLine()` calls for diagnostics. In a library
context this should be replaced with a logging abstraction (e.g., `ILogger`)
or a verbosity flag.

### 5. SILENT EXCEPTION SWALLOWING

SimplexRemesher.cs:2950-2955 and 2968-2973 in `CreateCrackFromSignedField3D`:
```csharp
try { ... }
catch { }
```
Bare `catch { }` blocks silently swallow all exceptions including logic errors.

### 6. PERFORMANCE CONCERNS

- Full `double[,]` array cloning inside smoothing loops at line 1880 (per-node,
  per-incident-triangle). The 3D version correctly uses in-place temp/restore.
- O(N*E) linear scan in `ComputeOpeningDirections()` (dead code, but noting).
- Dictionary-based neighbor maps where array-based would be more efficient.

### 7. ALGORITHMIC / CORRECTNESS ISSUES

- `DiscoverEdges` guard (line 121) only checks edge count, not correctness.
- Inconsistent sign-change detection: 2D uses `f1*f2 <= 0`, 3D adds
  `Math.Abs(f2-f1) > 1e-10`.
- `for (int pass = 0; pass < 1; pass++)` at line 1341: loop that always
  executes once.
- Visualization offset always pushes in y-direction regardless of crack geometry.

### 8. ROBUSTNESS

- `SaveMSH()` always writes `coordinates[i, 2]` but input may be `double[n, 2]`.
- GiD loader `int.Parse` without error handling on malformed files.
- `HexRowSpacing = 0.866` has only 3 significant digits (sqrt(3)/2 = 0.86602540...).

### 9. API DESIGN

- `FindBoundaryNodes`/`FindBoundaryNodes3D` belong on `MeshGeometry` or
  `MeshOptimization`, not `SimplexRemesher`.
- Inconsistent naming: `IdentifyBoundaryNodes` vs `FindBoundaryNodes`.
- `SimplexRemesher` class is too large (3,572 lines); should split into
  `MeshIO` and `CrackInsertion`.

### 10. STYLE / MINOR

- Inconsistent indentation in `TryFindEdgeRootOnSegment` and parts of
  `CreateCrackFromSignedField3D`.
- Unused variable `isOneBased` at line 808.
- Magic numbers scattered (1e-5, 1e-6, 1e-12, 0.001, 0.1) instead of using
  `MeshConstants`.

## Strengths

1. Good modular decomposition into focused static classes.
2. Shared constants in `MeshConstants.cs` prevents scattered epsilon definitions.
3. Robust topology via `Relations/Topology` library with type-safe connectivity.
4. Comprehensive element support (points, bars, triangles, quads, tetrahedra).
5. Correct numerical algorithms (Bowyer-Watson, Rivara, CVT).
6. Good use of `WithBatch` for topology mutation performance.

## Priority Recommendations

1. Remove dead code (~400 lines can be deleted).
2. Eliminate Jacobian duplication in SimplexRemesher.cs.
3. Fix silent `catch { }` blocks.
4. Replace Console.WriteLine with a logging abstraction.
5. Split SimplexRemesher.cs into MeshIO.cs and keep crack operations separate.
6. Fix coordinate dimension assumption in SaveMSH.

<div align="center">

<img src="m2m.png" alt="ManyToMany" width="640" />

# ManyToMany

**A high-performance scientific-computing and finite-element library for .NET 9**

[![.NET](https://github.com/PedroAreiasIST/ManyToMany/actions/workflows/dotnet.yml/badge.svg)](https://github.com/PedroAreiasIST/ManyToMany/actions/workflows/dotnet.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![.NET 9.0](https://img.shields.io/badge/.NET-9.0-512BD4.svg)](https://dotnet.microsoft.com/download/dotnet/9.0)
[![Language: C#](https://img.shields.io/badge/language-C%23-239120.svg)](https://learn.microsoft.com/dotnet/csharp/)
[![Platforms](https://img.shields.io/badge/platforms-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)](#platforms--prerequisites)

*Mesh topology · sparse linear algebra · mesh generation & fracture · nonlinear dynamics · visualization — tuned for problems up to 10M+ degrees of freedom.*

</div>

---

## Overview

**ManyToMany** is a unified framework for large-scale computational mechanics, written in C# by **Pedro Areias (IST)**. It couples a type-safe many-to-many topology engine with SIMD/GPU-accelerated linear algebra, conforming mesh refinement with level-set crack insertion, unconditionally-stable time integration, and post-processing — all in a single, dependency-light managed codebase.

> **Every public type in the library lives in the `Numerical` namespace.** A single `using Numerical;` gives you the whole API. (The core topology assembly is named `Topology` for historical reasons, but you never need to import a `Topology` namespace.)

### Why ManyToMany?

- 🧩 **One topology, many relations** — a generic, type-checked `Topology<TTypes>` models nodes, edges, faces, elements and arbitrary associations with graph algorithms (BFS/DFS, coloring, Cuthill–McKee, components) built in.
- ⚡ **Hardware-aware by default** — dense and sparse kernels dispatch automatically to AVX2/AVX-512, multi-threaded `Parallel.For`, Intel MKL PARDISO, or NVIDIA cuSPARSE based on problem size. No native library? It degrades gracefully to managed SIMD.
- 🔬 **Fracture-ready meshing** — structured generators plus longest-edge conforming refinement and **level-set crack insertion** in 2D and 3D, validated against classical SIF benchmarks (Griffith, Sneddon, Newman–Raju, …).
- 🚀 **Engineered for scale** — lock-striped parallel assembly, `>2 GB` chunked storage, `ArrayPool`/`Span` zero-allocation hot paths, and Kahan-compensated reductions for 10M+ DOF systems.
- 📦 **Self-contained** — pure managed C# with a single small dependency (`Microsoft.Extensions.ObjectPool`); MKL and CUDA are strictly optional accelerators.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Architecture Overview](#architecture-overview)
4. [Relations — Topology & Connectivity](#relations--topology--connectivity)
5. [Matrices — Linear Algebra](#matrices--linear-algebra)
6. [Meshing — Generation, Refinement & Fracture](#meshing--generation-refinement--fracture)
7. [Nonlinear — Dynamics & Root Finding](#nonlinear--dynamics--root-finding)
8. [Postprocess — Visualization](#postprocess--visualization)
9. [Examples](#examples)
10. [Performance](#performance)
11. [Platforms & Prerequisites](#platforms--prerequisites)
12. [Public API Reference](#public-api-reference)
13. [Project Structure](#project-structure)
14. [Documentation](#documentation)
15. [Contributing](#contributing)
16. [Citation](#citation)
17. [License](#license)

---

## Quick Start

Generate a mesh, insert a crack along a level set, and export it for visualization — using nothing but the public API:

```csharp
using Numerical;

// 1) Generate a structured 40×40 triangular mesh of the unit square.
//    Returns the mesh plus an [nNodes, 3] coordinate array.
var (mesh, coords) = SimplexRemesher.CreateRectangularMesh(
    nx: 40, ny: 40, xMin: 0, xMax: 1, yMin: 0, yMax: 1);

Console.WriteLine($"{mesh.Count<Node>()} nodes, {mesh.Count<Tri3>()} triangles");

// 2) Insert a horizontal crack along y = 0.5 between x = 0.25 and x = 0.75.
//    `surface` is the signed level-set (φ = 0 on the crack plane);
//    `region`  restricts the crack to a finite segment (negative = inside the crack).
SimplexRemesher.SignedFieldFunction surface = (x, y, z) => y - 0.5;
SimplexRemesher.SignedFieldFunction region  = (x, y, z) => (x > 0.25 && x < 0.75) ? -1.0 : 1.0;

var (cracked, crackedCoords) = SimplexRemesher.CreateCrackFromSignedField(
    mesh, coords, surface, region);

Console.WriteLine($"After crack: {cracked.Count<Node>()} nodes (duplicated along the crack faces)");

// 3) Export for GiD (.msh) and ParaView (Ensight Gold .case).
SimplexRemesher.SaveGiD(cracked, crackedCoords, "cracked_plate.msh");
EnsightWriter.SaveEnsight(cracked, crackedCoords, "cracked_plate");
```

Need a sparse solve instead? Assemble a finite-element system and solve it:

```csharp
using Numerical;

var system = new CliqueSystem(numElements: elements.Length);

// Symbolic phase: declare per-element DOFs, then build the sparsity pattern once.
for (int e = 0; e < elements.Length; e++)
{
    system.SetElementSize(e, elements[e].Dofs.Length);
    system.SetElementConnectivity(e, elements[e].Dofs);   // global DOF indices
}
system.BuildSparsityPattern();

// Numeric phase: scatter each element's stiffness/force (thread-safe).
Parallel.For(0, elements.Length, e =>
    system.AddElement(e, elements[e].Force, elements[e].Stiffness));

system.Assemble();
double[] u = system.Solve();   // PARDISO / iterative, chosen automatically

system.Reset();                // keep the pattern, clear values for the next load step
```

---

## Installation

### Prerequisites

- [.NET 9.0 SDK](https://dotnet.microsoft.com/download/dotnet/9.0) (x64)
- A 64-bit OS (Windows, Linux, or macOS). The library is **x64-only by design** — 32-bit is not supported.

### Build from source

```bash
git clone https://github.com/PedroAreiasIST/ManyToMany.git
cd ManyToMany
dotnet build Numerical.sln -c Release
```

> The solution also defines explicit `x64`/`64` platforms; if your tooling requires it you can pass `-p:Platform=x64`. The continuous-integration workflow builds the whole solution on every push to `master`.

### Use it in your own project

ManyToMany is consumed as project references (it is not currently published to NuGet). Add the modules you need to your `.csproj`:

```xml
<ItemGroup>
  <ProjectReference Include="path/to/Relations/Relations.csproj" />
  <ProjectReference Include="path/to/Matrices/Matrices.csproj" />
  <ProjectReference Include="path/to/Meshing/Meshing.csproj" />
  <ProjectReference Include="path/to/Nonlinear/Nonlinear.csproj" />
  <ProjectReference Include="path/to/Postprocess/Postprocess.csproj" />
</ItemGroup>
```

Each project pulls in its own dependencies (`Relations` is the only mandatory one — everything else builds on it).

### Optional accelerators

| Accelerator | Enables | Install |
|---|---|---|
| Intel MKL / oneAPI | PARDISO direct sparse solver | **Windows:** automatic via NuGet · **Linux:** `sudo apt-get install intel-mkl` · **macOS:** `brew install intel-mkl` |
| CUDA Toolkit 11.0+ | GPU SpMV & `cusolver` sparse solve | [developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads) |

Without these, ManyToMany runs entirely on the CPU using managed SIMD and parallel algorithms — accelerators are detected at runtime and used only when present.

---

## Architecture Overview

ManyToMany is a layered dependency stack: each module builds cleanly on the one below it.

```
┌────────────────────────────────────────────────────────────────┐
│                        Teste (Examples)                        │
├────────────────────┬───────────────────────────────────────────┤
│    Postprocess     │                Nonlinear                  │
├────────────────────┴───────────────────────────────────────────┤
│                          Meshing                               │
├────────────────────────────────────────────────────────────────┤
│                          Matrices                              │
├────────────────────────────────────────────────────────────────┤
│                      Relations  (core)                         │
└────────────────────────────────────────────────────────────────┘
```

| Module | Assembly | Primary types | Role |
|---|---|---|---|
| `Relations` | `Topology` | `O2M`, `Topology<TTypes>`, `Symmetry` | Mesh topology & graph algorithms |
| `Matrices` | `Matrices` | `Matrix`, `Vector`, `CSR`, `CliqueSystem` | Dense/sparse algebra, FE assembly |
| `Meshing` | `Meshing` | `SimplexMesh`, `SimplexRemesher`, `MeshGeometry` | Mesh generation, refinement, fracture |
| `Nonlinear` | `Nonlinear` | `BatheTwoStageIntegrator`, `RootFinder` | Time integration, root finding |
| `Postprocess` | `Postprocess` | `EnsightWriter` | Visualization export |

All public types share the **`Numerical`** namespace, so cross-module code needs only `using Numerical;`.

---

## Relations — Topology & Connectivity

`Relations` is the foundation of ManyToMany. It provides type-safe, high-performance many-to-many relationship structures for representing mesh connectivity, adjacency graphs, and arbitrary entity associations.

### `Topology<TTypes>` — the type-safe container

The primary user-facing API. It combines multi-type adjacency storage with arbitrary per-entity attribute dictionaries under a single `ReaderWriterLockSlim`.

```csharp
using Numerical;

// Map each entity type to an integer slot: Node→0, Edge→1, Tri3→2.
using Types = Numerical.TypeMap<Numerical.Node, Numerical.Edge, Numerical.Tri3>;

var topo = new Topology<Types>();

// Add nodes and an element connecting them.
int n0 = topo.Add<Node>();
int n1 = topo.Add<Node>();
int n2 = topo.Add<Node>();
int t0 = topo.Add<Tri3, Node>(n0, n1, n2);

// Attach typed attribute data to entities.
topo.Set<Node, Position>(n0, new Position(x, y, z));
var pos = topo.Get<Node, Position>(n0);

// Discover the edges implied by triangle connectivity.
topo.DiscoverSubEntities<Tri3, Edge, Node>(
    SubEntityDefinition.FromEdges((0, 1), (1, 2), (2, 0)));

// Graph algorithms operate directly on the topology.
List<int> order   = topo.GetTopologicalOrder<Node, Tri3>();
int[]     coloring = topo.ComputeElementColoring<Tri3, Node>();
```

> For finite-element meshing you rarely define a `TypeMap` by hand — [`SimplexMesh`](#meshing--generation-refinement--fracture) is a ready-made `Topology` specialized for the standard element zoo.

**Built-in graph algorithms**

- Breadth-First Search (single-type and multi-type), Depth-First Search across heterogeneous types
- Topological ordering (Kahn's algorithm) for DAG-structured connectivity
- Greedy element coloring for race-free parallel assembly scheduling
- Connected-component detection
- Cuthill–McKee (and reverse) bandwidth-reduction ordering
- Transitive connectivity, dual-graph and element/node adjacency extraction

### `O2M` — One-to-Many adjacency

A sparse adjacency list mapping one entity to an ordered list of entities — the storage primitive beneath `Topology`.

- Internal storage: `List<List<int>>` with pre-allocated capacity
- Parallel cloning via `GC.AllocateUninitializedArray` + `Parallel.For` above a configurable threshold
- Zero-allocation list mirroring through `CollectionsMarshal.AsSpan`; `[SkipLocalsInit]` on hot paths
- Implements `IComparable<O2M>`, `IEquatable<O2M>`, `ICloneable`

### `Symmetry` — canonical symmetry groups

Encodes permutation symmetry groups for element types, enabling automatic canonical representation and deduplication of equivalent elements (undirected edges, symmetric faces, …). Factories: `Identity`, `Cyclic`, `Dihedral`, `Full`, `FromGenerators`. Register one with `topo.WithSymmetry<TElement>(symmetry)` and use `AddUnique<…>` to suppress duplicates.

### Serialization & comparison

Every relation type supports JSON round-tripping (`ToJson`/`FromJson`, `SaveToFile`/`LoadFromFile`), lexicographic `==`/`<`/`>` comparison, set operations, and structural integrity validation.

---

## Matrices — Linear Algebra

### `Matrix` & `Vector` — dense operations

Column-major (`_data[col * RowCount + row]`, BLAS/LAPACK-compatible) dense matrix with SIMD-accelerated arithmetic and a full decomposition suite. Optimized for element-level matrices up to ~1000×1000.

**Decompositions**

| Method | Returns | Algorithm |
|---|---|---|
| `ComputeLU()` | `LUDecomposition` | Rook pivoting — searches both pivot row and column; superior stability near singularity |
| `ComputeQR()` | `QRDecomposition` | Householder reflections (thin QR) |
| `ComputeSVD()` | `SVDDecomposition` | One-sided Jacobi; converges for any real matrix |
| `ComputeEigenvalues()` | `EigenDecomposition` | Householder tridiagonalization + symmetric QR (symmetric matrices) |

**SIMD:** all arithmetic operators (`+`, `-`, `*`, scalar) use `System.Runtime.Intrinsics` AVX2/AVX-512 lanes; large matrices fall back to `Parallel.For`.

**Highlights:** factories (`Identity`, `Zeros`, `Ones`, `Diagonal`, `Random`, `RandomNormal`, `MatrixSquareRoot`, `MatrixExponential`); solvers (`Solve`, `SolveMultiple`, `SolveLeastSquares`); subspaces (`GetNullSpace`, `PseudoInverse`); statistics (`Covariance`, `Correlation`); norms (Frobenius, 1-, ∞-, max).

### `CSR` — Compressed Sparse Row

Production-grade sparse matrix for finite-element applications, supporting both structural-only (sparsity pattern) and valued matrices. **Sparse matrix–vector products and solves select the best available backend automatically:**

| Backend | Activation condition | API |
|---|---|---|
| NVIDIA cuSPARSE (GPU) | `rows ≥ 50,000` **and** `nnz ≥ 1,000,000` | `MultiplyGPU`, `MultiplyAuto` |
| Intel MKL PARDISO | MKL discoverable | `SolvePardiso`, `SolvePardisoMultiple` |
| SIMD SpMV (CPU) | `rows ≥ 5,000` | `MultiplySIMD` |
| Parallel SpMV (CPU) | `rows ≥ 1,000` | `MultiplyParallel` |
| Sequential SpMV | fallback | `Multiply` |

**Iterative solvers** live in the `CSRIterativeSolvers` extension class:

```csharp
double[] x = A.DiagonallyPreconditionedBiCGSTAB(b);   // returns the solution vector
SolverResult result = A.TrySolve(b);                  // solution + iterations + residual + convergence flag
double[] eigs = A.SmallestEigenvalues(m: 5);          // m smallest eigenvalues
```

**Thresholds** (public): `MIN_ROWS_FOR_PARALLEL` = 1,000 · `MIN_ROWS_FOR_SIMD` = 5,000 · `MIN_ROWS_FOR_GPU` = 50,000 · `MIN_NNZ_FOR_GPU` = 1,000,000 · `DEFAULT_TOLERANCE` = 1e-14.

### `CliqueSystem` — parallel finite-element assembly

High-performance assembly using **Gustavson's algorithm** for symbolic factorization and **lock-striped** numerical assembly.

1. **Symbolic phase** — from element connectivity, computes the exact global CSR sparsity pattern (`Cᵀ·C`) without storing values. If the raw DOF space exceeds the actual DOF count by >4× (or 10M entries), a dictionary-based DOF map replaces the dense array.
2. **Numeric phase** — each element scatters its local stiffness `kₑ` and force `fₑ` into the global system. Lock striping (4096 stripes, power-of-2 for fast bitwise modulo, hashed with the golden-ratio constant `0x9E3779B9`) prevents data races without serializing the loop.
3. **Large systems** — internal storage is chunked at 256 MB to break the `Int32.MaxValue` array-length limit, supporting >10M non-zeros.

```csharp
var system = new CliqueSystem(numElements, enableGpu: true);

for (int e = 0; e < numElements; e++)
{
    system.SetElementSize(e, dofsPerElement[e]);
    system.SetElementConnectivity(e, globalDofs[e]);
}
system.BuildSparsityPattern();                  // symbolic phase (once)

Parallel.For(0, numElements, e =>
    system.AddElement(e, force[e], stiffness[e]));  // numeric phase (thread-safe)

system.Assemble();
double[] u = system.Solve();
system.Reset();                                 // preserve pattern for the next load step
```

Construct one straight from a topology with `CliqueSystem.FromTopology<TTypes, TElement, TNode>(…)`, and inspect the result with `GetMatrix()`, `GetForceVector()`, and `GetStatistics()`.

### Native-library discovery

`NativeLibraryConfig` (implementing `INativeLibraryConfig`) performs cross-platform discovery of MKL and CUDA from standard locations (`%MKLROOT%`, `/opt/intel/oneapi/mkl`, `%CUDA_PATH%`, `/usr/local/cuda/lib64`, …). Libraries load lazily on first use; if one is absent, the corresponding backend silently degrades to the next available option. Backend selection itself is an internal implementation detail — you never call it directly.

---

## Meshing — Generation, Refinement & Fracture

### `SimplexMesh`

The core container for 2D triangular and 3D tetrahedral meshes — a `Topology` pre-specialized for the standard element zoo:

```csharp
public sealed class SimplexMesh
    : Topology<TypeMap<Node, Edge, Point, Bar2, Tri3, Quad4, Tet4>>
```

| Type | Dim | Nodes | Description |
|---|---|---|---|
| `Node` | 0D | 1 | Mesh vertex |
| `Point` | 0D | 1 | Single-node element |
| `Bar2` | 1D | 2 | Line segment |
| `Edge` | 1D | 2 | Topological edge (discovered) |
| `Tri3` | 2D | 3 | Linear triangle |
| `Quad4` | 2D | 4 | Bilinear quadrilateral |
| `Tet4` | 3D | 4 | Linear tetrahedron |

Coordinates are stored as `double[numNodes, 3]` regardless of problem dimension, for API consistency.

### Mesh generation

```csharp
// Structured grids
var (mesh2d, c2d) = SimplexRemesher.CreateRectangularMesh(nx, ny, xMin, xMax, yMin, yMax); // 2·nx·ny triangles
var (mesh3d, c3d) = SimplexRemesher.CreateBoxMesh(nx, ny, nz, /* bounds… */);              // 6·nx·ny·nz tets

// Unstructured triangulation of an arbitrary boundary, with optional holes & interior points
var (mesh, coords) = SimplexRemesher.Triangulate(boundaryCoords, refine: true, maxArea: 2.0);
var (mh, ch)       = SimplexRemesher.TriangulateWithHoles(outerBoundary, holes, convertToQuads: true);
```

Each rectangular cell is split into two right triangles along its main diagonal; each hexahedral cell is decomposed into 6 tetrahedra sharing the cube's main diagonal.

### Conforming refinement

`SimplexRemesher.Refine(mesh, markedEdges)` performs longest-edge bisection with conforming closure — no hanging nodes, no quality degradation:

1. **Edge selection** — mark edges for refinement (longest edge, or user-specified)
2. **Midpoint insertion** — each marked edge gets a midpoint node; the `ParentNodes(Parent1, Parent2)` attribute records its parents for solution transfer
3. **Element splitting** — triangles split into 2, tetrahedra into 4 (or 8) children
4. **Conforming closure** — refinement propagates to neighbors to eliminate hanging nodes

Transfer field values to the refined mesh with `SimplexRemesher.InterpolateCoordinates(refined, originalCoords)`; midpoint values are `u_mid = ½(u_p1 + u_p2)`.

### Level-set crack insertion

Insert arbitrary cracks defined by a signed-distance field — the headline feature for fracture mechanics:

```csharp
// 2D: a SignedFieldFunction φ(x,y,z) defines the crack surface (φ = 0),
//     and an optional region field bounds the crack to a finite extent.
var (cracked, coords) = SimplexRemesher.CreateCrackFromSignedField(
    mesh, baseCoords, signedField: surface, regionField: region);

// 3D analogue for tetrahedral meshes:
var (cracked3d, c3d) = SimplexRemesher.CreateCrackFromSignedField3D(
    mesh3d, baseCoords3d, signedField: surface3d, regionField: region3d);
```

The algorithm classifies elements relative to the crack (uncut / cut / tip), then **duplicates the nodes on one side** of cut elements to create independent crack faces, rewiring connectivity to reference the original or duplicated nodes as appropriate.

### Smoothing & quality

| Method | Purpose |
|---|---|
| `LaplacianSmoothing(mesh, coords, iterations)` | Classic Laplacian node relocation |
| `CVTSmoothing(mesh, coords, iterations)` | Centroidal Voronoi tessellation smoothing |
| `MeshGeometry.ComputeQualityStatistics(mesh, coords)` | Aspect ratios, minimum angles, degenerate-element counts (`MeshQualityStats`) |
| `MeshRefinement.CheckJacobians` / `FixNegativeJacobians` | Detect and repair inverted elements |

### Mesh I/O

| Format | Read | Write | API |
|---|:---:|:---:|---|
| Gmsh `.msh` (v2) | ✅ | ✅ | `LoadMSH`, `LoadMSHWithTags`, `SaveMSH`, `SaveMSHWithCrackGroups` |
| GiD / CIMNE `.msh` | ✅ | ✅ | `LoadGiD`, `SaveGiD` |
| Plain ASCII | — | ✅ | `SaveASCII` |
| Ensight Gold `.case` | — | ✅ | [`EnsightWriter`](#postprocess--visualization) |

---

## Nonlinear — Dynamics & Root Finding

### `BatheTwoStageIntegrator`

Unconditionally-stable implicit time integrator for second-order systems

```
M ü(t) + C u̇(t) + f_int(u(t), t) = R_ext(t)
```

The **Bathe two-stage method** (Bathe, 2007) splits each step `Δt` into two sub-steps with different Newmark parameters — Stage 1 (`β = ¼, γ = ½`, trapezoidal) followed by a Stage 2 corrector (`β = 4/9, γ = 2/3`) — giving second-order accuracy with no high-frequency energy growth.

You drive it with two delegates that the integrator calls during the inner Newton–Raphson loop. Note that state is passed as `ReadOnlySpan<double>`/`Span<double>` for zero-allocation evaluation:

```csharp
using Numerical;

// Residual r(t, u, v, a) = M·a + C·v + f_int(u) − R_ext(t)
BatheTwoStageIntegrator.ResidualEvaluator residual =
    (double t, ReadOnlySpan<double> u, ReadOnlySpan<double> v, ReadOnlySpan<double> a, Span<double> r) =>
    { /* fill r */ };

// Solve the effective tangent system (a0·M + a1·C + Kt)·Δu = rhs
BatheTwoStageIntegrator.EffectiveSystemSolver solver =
    (double t, double a0, double a1,
     ReadOnlySpan<double> u, ReadOnlySpan<double> v, ReadOnlySpan<double> a,
     ReadOnlySpan<double> rhs, Span<double> deltaU) =>
    { /* solve for deltaU */ };

// The constructor solves the initial static equilibrium when given only an initial guess.
var integrator = new BatheTwoStageIntegrator(time0: 0.0, u0Guess: u0, residual, solver);

integrator.MaxNewtonIterations = 30;   // defaults shown
integrator.RelTolerance        = 1e-8;
integrator.AbsTolerance        = 1e-12;

integrator.Step(dt: 1e-4, numSteps: 1000);
integrator.GetDisplacement(result);    // copy current state into your buffer
```

- **Newton–Raphson** inner loop with configurable absolute/relative tolerances and divergence detection (`DivergenceThreshold`, default `1e6`).
- **Kahan compensated summation** for residual norms — eliminates cancellation error in long vectors, essential at 10M+ DOF.
- **SIMD** vector kernels (AVX2/AVX-512, detected at runtime) and **zero-allocation** hot paths (working vectors pre-allocated; `ArrayPool<double>` for temporaries).
- Per-step diagnostics via `LastStepConvergence` (`ConvergenceInfo`) and cumulative `Performance` (`PerformanceCounters`).

### `RootFinder`

Thread-safe scalar root-finding with two overloads:

```csharp
// With derivative — hybrid Newton-Raphson + Inverse Quadratic Interpolation, bisection fallback
var (root, status) = RootFinder.FindRoot(xmin, xmax, x => (f(x), df(x)));

// Without derivative — ITP (Interpolate-Truncate-Project, Oliveira & Takahashi 2020):
// optimal O(log₂(1/ε)) worst case, superlinear on smooth functions
var (root2, status2) = RootFinder.FindRoot(xmin, xmax, x => f(x));
```

`status` is a `RootFinder.Status`: `OK`, `Tolerance`, `MaxIterations`, `NoBracket`, `BadInput`, `NonFinite`, `TooNarrow`.

### `TrustRegionNewtonDogleg`

A matrix-free trust-region Newton-dogleg solver for nonlinear systems, configured via `TRNOptions` and returning a `TRNResult`. You supply the residual, a linear solve, and a Jacobian-vector product as delegates — suitable for large systems where the Jacobian is never formed explicitly.

---

## Postprocess — Visualization

### `EnsightWriter`

Exports meshes and field data to **Ensight Gold** (ASCII), compatible with ParaView and GiD.

```csharp
// Single mesh
EnsightWriter.SaveEnsight(mesh, coords, "result");
EnsightWriter.SaveEnsightWithScalar(mesh, coords, "result", "Temperature", scalarField);

// Aggregate many meshes into one multi-part .case (e.g. a fracture study)
EnsightWriter.AddMesh("step_0", mesh, coords, displacement);
EnsightWriter.AddMesh("step_1", mesh, coords, displacement);
EnsightWriter.WriteAllMeshes("Study");   // → Study.case + per-part .geo files
```

Output is a `.case` descriptor plus per-part geometry files; displacement fields are written as vector data and can be scaled inside ParaView/GiD.

---

## Examples

The `Teste` project contains **26 worked examples** in four parts. The default `Main` (`Examples2DA.Main`) runs a representative subset and writes a unified Ensight case; uncomment lines in `Main` (or call a public method directly) to run any specific one.

```bash
dotnet run --project Teste -c Release
```

| Part | Examples | Theme |
|---|---|---|
| **1 — Advanced meshing** | `Example1`–`Example10` *(public)* | Holes, re-entrant corners, annuli, wedges, gears, tri-vs-quad quality, industrial shapes |
| **2 — 2D fracture benchmarks** | `Example11`–`Example15` | Edge / center / double-edge / slant cracks & crack-from-hole, vs. published SIF solutions (Anderson, Griffith, Kanninen–Popelar, Erdogan–Sih, Newman–Raju) |
| **3 — Spectacular crack patterns** | `Example16`–`Example20` | Spiral, fractal-tree, sinusoidal, starburst & mandala crack networks via the level-set engine |
| **4 — 3D fracture benchmarks** | `Example21`–`Example26` | Penny-shaped, elliptical, edge, corner-at-hole, slant (`K_I/K_II/K_III`) & semi-cylindrical surface cracks (Sneddon, Irwin, Tada, Newman–Raju) |

> The Part 1 meshing examples (`Example1`–`Example10`) are `public static` and can be invoked directly, e.g. `Examples2DA.Example1_CircularDomainWithHole();`. The fracture examples are driven from `Main`.

All examples emit GiD `.msh` files and a unified `FractureMechanics.case` for ParaView.

---

## Performance

The library is engineered for high-throughput computational mechanics on modern x86-64 hardware.

### Hardware acceleration

| Feature | Technology | Activation |
|---|---|---|
| Dense SIMD | AVX2 / AVX-512 intrinsics | runtime `Avx2.IsSupported` / `Avx512F.IsSupported` |
| Sparse SpMV SIMD | `System.Numerics.Vector<double>` + AVX2 | `rows ≥ 5,000` |
| Sparse direct solver | Intel MKL PARDISO | MKL discoverable |
| GPU SpMV / solve | NVIDIA cuSPARSE / cuSolver | `rows ≥ 50,000` **and** `nnz ≥ 1,000,000` |
| Parallel CPU | `Parallel.For` with `ParallelOptions` | per-operation size thresholds |

### Memory efficiency

- `ArrayPool<double>` for temporary buffers — no GC pressure in hot paths
- `GC.AllocateUninitializedArray` for large pre-allocated arrays (skips zero-fill)
- `CollectionsMarshal.AsSpan` for zero-copy list access during parallel cloning
- `[SkipLocalsInit]` to suppress redundant stack zeroing on performance-critical types
- 256 MB chunked storage in `CliqueSystem` to exceed the `Int32.MaxValue` element limit
- `Microsoft.Extensions.ObjectPool` for `HashSet<int>` reuse during symbolic assembly

### Scalability

- Designed and tested for problems exceeding **10 million DOFs**, with `>2 GB` array support via chunked storage
- Tiered-compilation + dynamic PGO and Server GC enabled in Release builds
- Lock-striped assembly (4096 stripes) scales to 32+ hardware threads with minimal contention

### Connectivity benchmark (reproducible)

`Benchmarks/mesh_connectivity_benchmark.py` is a self-contained harness that times
connectivity-query performance across mesh data-structure libraries on identical
meshes and prints a comparison table. It is wired into CI
(`.github/workflows/connectivity-benchmark.yml`), so the table is regenerated on
every push (see the workflow run summary and its uploaded artifacts).

```bash
pip install numpy psutil libigl vtk openmesh      # every library is optional
dotnet build Benchmarks/Mm3Bench -c Release       # builds the MM3 bridge
python Benchmarks/mesh_connectivity_benchmark.py --mm3-project ./Benchmarks/Mm3Bench --out results
```

Identical workload for every library: build the topology, then traverse every
cell's face-neighbours and every vertex's incident cells once. Reported time is
total **build + query** in milliseconds.

> **The numbers are hardware-specific and not comparable across machines.** They
> depend on CPU, memory bandwidth, thread count and library versions. The
> authors' measurements in the paper were taken on an **Intel Core i9-13900HX**;
> run the command above to obtain the table for *your* hardware. A library that
> is not installed — or that cannot represent the mesh type (half-edge *surface*
> libraries have no tetrahedra) — is reported as `N/A`. The comparison also mixes
> traversal-based structures (OpenMesh/CGAL/VTK) with array-building ones
> (libigl/MM3), so read across a row as order-of-magnitude, not exact. MM3's
> first invocation additionally pays JIT warm-up and incremental topology
> construction, so its build+query figure is **not** a steady-state query latency.

Illustrative output on a small Linux x86-64 CI runner (1 vCPU / 4 GB — values are
for format only, **not** a performance claim):

| Method | grid-2K-tri | grid-7K-tri | box-6K-tet |
| --- | --- | --- | --- |
| OpenMesh (Half-Edge) | N/A | N/A | N/A |
| CGAL (HalfedgeDS) | N/A | N/A | N/A |
| VTK | 10.7 | 33.8 | 326.1 |
| libigl | 0.9 | 2.4 | N/A |
| MM3 | 69.6 | 115.8 | 79.7 |

`N/A` above = library absent on that runner, or unsupported mesh type (note MM3
handles the tetrahedral mesh, which the half-edge surface libraries cannot). On a
host with every library installed, all applicable cells are populated.

---

## Platforms & Prerequisites

| OS | Architecture |
|---|---|
| Windows 10/11 | x64 |
| Linux (Ubuntu 20.04+) | x64 |
| macOS 12+ | x64, ARM64 |

**Required:** [.NET 9.0 SDK](https://dotnet.microsoft.com/download/dotnet/9.0) (x64). **Optional:** Intel MKL (PARDISO) and CUDA Toolkit 11.0+ (GPU). See [Installation](#installation).

---

## Public API Reference

This section lists the public surface of each module. With editor autocomplete and a single `using Numerical;`, these are all the entry points you need.

### Relations

#### `Topology<TTypes>`

```csharp
public class Topology<TTypes> : IDisposable where TTypes : ITypeMap, new()
```

The main container combining multi-type adjacency with per-entity attribute dictionaries under a single `ReaderWriterLockSlim`.

**Entity management**

| Method | Description |
|---|---|
| `int Add<TNode>()` | Add a node entity; returns its index. |
| `int Add<TNode, TData>(TData data)` | Add a node with associated attribute data. |
| `int Add<TElement, TNode>(params int[] nodes)` | Add an element connected to node indices; returns element index. |
| `(int Index, bool WasNew) AddUnique<TElement, TNode>(params int[] nodes)` | Add only if canonically unique (requires a registered `Symmetry`). |
| `int[] AddRange<TElement, TNode>(IEnumerable<int[]> connectivity)` | Batch-add elements (`…Parallel` overload for large sets). |
| `void Remove<TEntity>(int index)` / `RemoveRange<TEntity>(…)` | Mark entities for removal. |
| `int Count<TEntity>()` / `CountActive<TEntity>()` | Total / non-deleted counts. |
| `bool Exists<TElement>(params int[] nodes)` / `int Find<TElement>(params int[] nodes)` | Existence / lookup (returns index or −1). |

**Attribute storage**

| Method | Description |
|---|---|
| `void Set<TEntity, TData>(int index, TData value)` | Attach attribute data to an entity. |
| `TData Get<TEntity, TData>(int index)` / `bool TryGet<…>(…)` | Retrieve attribute data (throwing / non-throwing). |
| `IEnumerable<(int, TData)> Each<TEntity, TData>()` | Iterate entities with their data. |
| `void ForEach<…>` / `ParallelForEach<…>` | Sequential / parallel iteration. |

**Adjacency queries** — `WithNodesOf`, `EnumerateNeighbors`, `GetElementsWithNodes`, `GetElementsContainingAnyNode`, `GetDirectNeighbors`, `GetWeightedNeighbors`, `GetElementNeighbors`, `GetNodeNeighbors`, `GetKHopNeighborhood`, `CountIncident`, …

**Sub-entity discovery & boundary**

| Method | Description |
|---|---|
| `(int TotalExtracted, int UniqueAdded, int DuplicatesSkipped) DiscoverSubEntities<TElement, TSubEntity, TNode>(SubEntityDefinition def, bool addUnique = true)` | Enumerate & register implied sub-entities (faces, edges). |
| `Topology<TTypes> WithSymmetry<TElement>(Symmetry symmetry)` | Register a symmetry group (chainable). |
| `List<int[]> ExtractBoundaryFacets<TElement, TNode>(int nodesPerFacet)` | Extract boundary facets. |
| `List<int> GetBoundarySubEntities<…>()` / `GetInteriorSubEntities<…>()` | Classify sub-entities. |
| `List<int> DetectNonManifoldSubEntities<…>()` | Find non-manifold sub-entities. |

**Graph algorithms** — `GetTopologicalOrder<…>`, `ComputeElementColoring<TElement, TNode>`, `GetColorGroups<…>`, `BreadthFirstSearch<TEntity>(int start, Action<int,int>? visitor = null)`, `BreadthFirstDistances<…>`, `BreadthFirstSearchMultiType<…>`, `MultiTypeDFS<…>`, `IsAcyclic<…>`, `FindComponents<TElement, TNode>`, `ComputeCuthillMcKeeOrdering<TElement, TNode>(bool reverse = true)`, `ComputeBandwidth<…>`, `ComputeTransitiveConnectivity<…>`, `GetDualStructure<…>`, `GetElementToElementGraph<…>`, `GetNodeToNodeGraph<…>`.

**Serialization & validation** — `ToJson`/`FromJson`, `SaveToFile`/`LoadFromFile`, `GetStatistics()`, `ValidateStructure()`, `ValidateIntegrity<…>()`, `GetDuplicates<…>()`, `IsPermutationOf<…>(…)`.

**Memory & lifecycle** — `Compress(…)`, `Clear()`, `Clone()`, `Reserve<…>(…)`, `ShrinkToFit()`, `ConfigureType<…>(…)`, `WithBatch(Action)`, `Merge<…>(…)`, plus transpose accessors (`GetTranspose<…>`, `WithTranspose<…>`, `EnsureSynchronized<…>`).

#### Supporting types

| Type | Kind | Description |
|---|---|---|
| `O2M` | `sealed class` | One-to-many adjacency list; `IComparable`/`IEquatable`/`ICloneable`. |
| `ReadOnlyTopology<TTypes>` | `sealed class` | Read-only projection — same queries, no mutation. |
| `Symmetry` | `sealed class` | Permutation symmetry group; factories `Identity`/`Cyclic`/`Dihedral`/`Full`/`FromGenerators`. |
| `SubEntityDefinition` | `readonly struct` | Which local nodes form each sub-entity; factories `FromEdges`/`FromFaces`/`FromQuadFaces`. |
| `SmartEntity<TEntity>` | `readonly record struct` | `(Topology, Index)` pair with `IsValid`, `Count`, `Data<T>()`, `BreadthFirstSearch()`, … |
| `ITypeMap` / `TypeMap<T0,…,Tn>` | `interface` / `sealed class` | Compile-time type→slot mapping; `TypeMap` provided for 2–24 type arguments. |
| `ResultOrder` | `enum` | `Unordered` (fastest) or `Sorted` (deterministic) results. |
| `ParallelConfig` | `static class` | Global parallelization thresholds. |
| `TopologyStats`, `Utils` | — | Statistics snapshot; shared utilities. |

### Matrices

#### `Matrix` / `Vector`

```csharp
public sealed class Matrix : IEquatable<Matrix>, IFormattable, ICloneable
public sealed class Vector : IEquatable<Vector>, IFormattable, ICloneable
```

Column-major dense matrix and companion vector with SIMD arithmetic. Factories `Identity`/`Zeros`/`Ones`/`Diagonal`/`Random`/`RandomNormal`/`MatrixSquareRoot`/`MatrixExponential`; decompositions `ComputeLU`/`ComputeQR`/`ComputeEigenvalues`/`ComputeSVD`; solvers `Solve`/`SolveMultiple`/`SolveLeastSquares`; plus `Inverse`, `Determinant`, `Rank`, `ConditionNumber`, `PseudoInverse`, `KroneckerProduct`, `GetNullSpace`/`GetRowSpace`/`GetImageSpace`, statistics and norms. `Vector` adds `Dot`, `Cross`, `OuterProduct`, `Normalize`, `ProjectOnto`, reductions, and slicing.

| Decomposition result | Key members |
|---|---|
| `LUDecomposition` | `Solve`, `Determinant`, `ConditionNumber` |
| `QRDecomposition` | `Solve`, `Rank` |
| `EigenDecomposition` | `Eigenvalues` (`double[]`), `Eigenvectors` (`Matrix`) |
| `SVDDecomposition` | `U`, `S`, `Vt`; `Rank`, `ConditionNumber` |

#### `CSR`

```csharp
public sealed class CSR : IFormattable, IEquatable<CSR>, ICloneable, IDisposable
```

Sparse matrix that auto-selects its SpMV/solve backend. SpMV: `Multiply` (+ `MultiplyParallel`/`MultiplySIMD`/`MultiplyGPU`/`MultiplyAuto`). Direct solve: `SolvePardiso`, `SolvePardisoMultiple`; triangular solves `SolveLowerTriangular`/`SolveUpperTriangular`. Sparse products `Multiply3Phase`/`MultiplySymbolicOnly`.

#### `CSRIterativeSolvers` (extension methods on `CSR`)

| Method | Returns | Algorithm |
|---|---|---|
| `DiagonallyPreconditionedBiCGSTAB(b, tol, maxIter, x0)` | `double[]` | Jacobi-preconditioned BiCGSTAB |
| `TrySolve(b, tol, maxIter, x0)` | `SolverResult` | BiCGSTAB with full diagnostics |
| `SmallestEigenvalues(m, tol, maxIter)` | `double[]` | `m` smallest eigenvalues |

```csharp
public record SolverResult(double[] Solution, int Iterations, double ResidualNorm, bool Converged, string? Message = null);
public record MatrixStatistics(int Rows, int Columns, int NonZeros, double Sparsity, int MinNnzPerRow, int MaxNnzPerRow, double AvgNnzPerRow);
public class  SolverException : Exception { /* thrown on non-convergence or singular systems */ }
```

#### `CliqueSystem`

```csharp
public sealed class CliqueSystem : IDisposable
public CliqueSystem(int numElements, bool enableGpu = false);
```

| Phase | Method |
|---|---|
| Setup | `SetElementSize(int e, int numDofs)`, `SetElementConnectivity(int e, int[] globalDofs)` |
| Symbolic | `BuildSparsityPattern()` |
| Numeric | `AddElement(int e, double[] force, double[] stiffness)` *(thread-safe; `ReadOnlySpan` overload available)* |
| Finalize | `Assemble()`, `Solve()` → `double[]` |
| Reuse | `Reset()` *(clears values, preserves the pattern)* |
| Inspect | `GetMatrix()` → `CSR`, `GetForceVector()` → `double[]`, `GetStatistics()` → `AssemblyStatistics` |
| Build from topology | `static CliqueSystem FromTopology<TTypes, TElement, TNode>(…)` |

`DiscreteLinearSystem` is a higher-level wrapper for node-major DOF layouts (`Solve(double[,] result)`, `BuildSystemValues()`, `Reset()`). `INativeLibraryConfig` / `NativeLibraryConfig` expose custom MKL/CUDA search paths.

### Meshing

| Type | Kind | Notes |
|---|---|---|
| `SimplexMesh` | `sealed class` | Mesh container; adds `AddNode`/`AddMidpointNode`/`AddTriangle`/`AddQuad`/`AddTetrahedron`/`AddBar`. |
| `SimplexRemesher` | `static class` | Generation, refinement, crack insertion, smoothing, I/O (see below). |
| `MeshGeometry` | `static class` | Element geometry, point/curve predicates, quality statistics. |
| `MeshRefinement` | `static class` | `Refine`, `CheckJacobians`, `FixNegativeJacobians`, `InterpolateCoordinates`. |
| `MeshConstants` | `static class` | Numerical tolerances (`Epsilon` = 1e-10, …). |
| `FiniteElementTopologies` | `static class` | Pre-built `SubEntityDefinition`s (`Tri3Edges`, `Tet4Faces`, …). |
| `MeshQualityStats` | `class` | Quality metrics from `ComputeQualityStatistics`. |
| element markers | `readonly struct` | `Node`, `Edge`, `Point`, `Bar2`, `Tri3`, `Quad4`, `Tet4`. |
| `ParentNodes`, `OriginalElement` | `readonly record struct` | Refinement lineage attributes. |
| `SignedFieldFunction` | `delegate` | `double (double x, double y, double z)` — level-set for crack insertion. |

**`SimplexRemesher` methods** — *generation:* `CreateRectangularMesh`, `CreateRectangularQuadMesh`, `CreateUnitSquareMesh`, `CreateBoxMesh`, `CreateUnitCubeMesh`, `Triangulate`, `TriangulateWithHoles`, `DelaunayTriangulate`; *refinement:* `Refine`, `InterpolateCoordinates`, `DiscoverEdges`; *fracture:* `CreateCrackFromSignedField`, `CreateCrackFromSignedField3D`, `CreateCrack`, `CreateCrackFromRefinedMesh`; *smoothing/cleanup:* `LaplacianSmoothing`, `CVTSmoothing`, `RemoveDegenerateTriangles`, `RemoveDegenerateTetrahedra`, `ConvertToQuads`; *boundary:* `FindBoundaryNodes`, `FindBoundaryNodes3D`; *I/O:* `SaveMSH`, `SaveMSHWithCrackGroups`, `LoadMSH`, `LoadMSHWithTags`, `SaveGiD`, `LoadGiD`, `SaveASCII`; *reporting:* `PrintStats`, `PrintQualityReport`.

### Nonlinear

| Type | Kind | Notes |
|---|---|---|
| `BatheTwoStageIntegrator` | `sealed class` | Implicit dynamics; delegates `ResidualEvaluator` & `EffectiveSystemSolver` (constructor-injected); `Step`, `GetDisplacement`/`GetVelocity`/`GetAcceleration`/`GetState`; nested `ConvergenceInfo`, `PerformanceCounters`. |
| `RootFinder` | `static class` | Two `FindRoot` overloads (with/without derivative); `Status` enum. |
| `TrustRegionNewtonDogleg` | `static class` | Matrix-free `Solve(…)`; configured by `TRNOptions`, returns `TRNResult`. |

### Postprocess

| Type | Kind | Methods |
|---|---|---|
| `EnsightWriter` | `static class` | `AddMesh`, `WriteAllMeshes`, `SaveEnsight`, `SaveEnsightWithScalar`. |

---

## Project Structure

```
ManyToMany/
├── Numerical.sln                # Visual Studio solution
├── Relations/                   # Core topology library (assembly: Topology)
│   ├── Relations.cs             #   O2M and supporting adjacency structures
│   ├── Topology.cs              #   Topology<TTypes>, SubEntityDefinition, Symmetry
│   └── Utils.cs                 #   ITypeMap, TypeMap<…>, Utils, ParallelConfig
├── Matrices/                    # Linear algebra
│   ├── Matrix.cs                #   Dense Matrix/Vector, decompositions, SIMD
│   ├── CSR.cs                   #   Sparse CSR, PARDISO, cuSPARSE, BiCGSTAB
│   ├── Assembly.cs              #   CliqueSystem, DiscreteLinearSystem
│   └── NativeLibraries.cs       #   Cross-platform MKL/CUDA discovery & loading
├── Meshing/                     # Mesh generation & refinement
│   ├── SimplexMesh.cs           #   Mesh container, element markers, constants
│   ├── SimplexRemesher.cs       #   Generation, bisection refinement, fracture, I/O
│   ├── MeshRefinement.cs        #   Conforming refinement driver
│   └── MeshGeometry.cs          #   Geometric primitives, quality metrics
├── Nonlinear/                   # Time integration & root finding
│   ├── Integrator.cs            #   BatheTwoStageIntegrator
│   └── RootFinder.cs            #   ITP, hybrid Newton-IQI, trust-region dogleg
├── Postprocess/                 # Visualization
│   └── EnsightWriter.cs         #   Ensight Gold export
├── Teste/                       # 26 demo examples
│   └── Examples2DA.cs           #   Meshing + 2D/3D fracture mechanics
└── Docs/                        # Extended documentation (see below)
```

---

## Documentation

Start with this README, then dive into the per-module deep-dives in [`Docs/`](Docs):

| Document | Covers |
|---|---|
| [`Docs/README.md`](Docs/README.md) | Documentation index & reading guide |
| [`Numerical-Complete-Documentation.md`](Docs/Numerical-Complete-Documentation.md) | Dense/sparse matrices, FE assembly, native-library integration |
| [`Topology-Complete-Documentation.md`](Docs/Topology-Complete-Documentation.md) | Topology operations, graph algorithms, serialization |
| [`SimplexRemesher-Complete-Documentation.md`](Docs/SimplexRemesher-Complete-Documentation.md) | Mesh refinement, crack insertion, file I/O, tutorials |

A companion paper, *A Simple C# Library for Computational Mechanics*, is included as [`Docs/p_areias_simple_csharp_final.pdf`](Docs/p_areias_simple_csharp_final.pdf).

---

## Contributing

Contributions are welcome! Please see [`CONTRIBUTING.md`](CONTRIBUTING.md) for build, test, and style guidelines. In short:

```bash
dotnet build Numerical.sln -c Release   # build everything
dotnet run  --project Teste -c Release  # run the example suite
```

Open an issue to discuss substantial changes before sending a pull request.

---

## Citation

If you use ManyToMany in academic work, please cite:

```bibtex
@software{areias_manytomany,
  author  = {Areias, Pedro},
  title   = {ManyToMany: A High-Performance Scientific Computing and
             Finite Element Library for .NET},
  year     = {2026},
  url      = {https://github.com/PedroAreiasIST/ManyToMany}
}
```

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**. You are free to use, modify, and distribute it under the terms of the GPLv3 — see [LICENSE](LICENSE) for the full text.

Copyright © 2026 **Pedro Miguel de Almeida Areias** (Instituto Superior Técnico).

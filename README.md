![ManyToMany](m2m.png)

# ManyToMany

A high-performance scientific computing and finite element analysis library for .NET 9.0, written in C# by Pedro Areias (IST).

ManyToMany provides a unified framework for managing complex mesh topologies, sparse linear algebra, mesh generation with crack insertion, nonlinear time integration, and post-processing — all tuned for large-scale computational mechanics problems up to 10M+ DOF.

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Relations — Topology & Connectivity](#relations--topology--connectivity)
3. [Matrices — Linear Algebra](#matrices--linear-algebra)
4. [Meshing — Generation & Refinement](#meshing--generation--refinement)
5. [Nonlinear — Dynamics & Root Finding](#nonlinear--dynamics--root-finding)
6. [Postprocess — Visualization](#postprocess--visualization)
7. [Project Structure](#project-structure)
8. [Examples](#examples)
9. [Performance](#performance)
10. [Building & Running](#building--running)
11. [Platforms & Prerequisites](#platforms--prerequisites)
12. [Documentation](#documentation)
13. [License](#license)

---

## Architecture Overview

ManyToMany is organized as a layered dependency stack: each module builds cleanly on the one below it.

```
┌────────────────────────────────────────────────────────────────┐
│                        Teste (Examples)                        │
├────────────────────┬───────────────────────────────────────────┤
│    Postprocess     │         Nonlinear                         │
├────────────────────┴───────────────────────────────────────────┤
│                    Meshing                                     │
├────────────────────────────────────────────────────────────────┤
│                    Matrices                                    │
├────────────────────────────────────────────────────────────────┤
│                    Relations  (core)                           │
└────────────────────────────────────────────────────────────────┘
```

| Module | Primary Types | Role |
|---|---|---|
| `Relations` | `O2M`, `M2M`, `MM2M`, `Topology<TTypes>`, `Symmetry` | Mesh topology, graph algorithms |
| `Matrices` | `Matrix`, `CSR`, `CliqueSystem` | Dense/sparse algebra, FE assembly |
| `Meshing` | `SimplexMesh`, `SimplexRemesher`, `CrackOperations` | Mesh generation, refinement, fracture |
| `Nonlinear` | `BatheTwoStageIntegrator`, `RootFinder` | Time integration, root finding |
| `Postprocess` | `EnsightWriter` | Visualization export |

---

## Relations — Topology & Connectivity

The `Relations` library is the foundation of ManyToMany. It provides type-safe, high-performance many-to-many relationship structures for representing mesh connectivity, adjacency graphs, and arbitrary entity associations.

### Relationship Hierarchy

#### `O2M` — One-to-Many

A sparse adjacency list mapping one entity to an ordered list of entities. Single-threaded, not thread-safe by design (thread safety is added at the `M2M` level).

- Internal storage: `List<List<int>>` with pre-allocated capacity
- Parallel cloning via `GC.AllocateUninitializedArray` + `Parallel.For` above a configurable parallelization threshold
- Deep-copy uses `CollectionsMarshal.AsSpan` for zero-allocation list mirroring
- `[SkipLocalsInit]` attribute eliminates redundant zero-fills in hot paths
- Implements `IComparable<O2M>`, `IEquatable<O2M>`, `ICloneable`

#### `M2M` — Many-to-Many

Thread-safe sparse adjacency structure wrapping `O2M` with `ReaderWriterLockSlim` synchronization. Provides a cached transpose and O(1) position lookups.

- **Transpose caching**: transposed adjacency is computed once and cached; invalidated on mutation
- **Position lookup**: given `(entity_i, entity_j)`, returns the index of `entity_j` in the adjacency list of `entity_i` in O(n_neighbors) time
- **Set operations**: intersection, union, and difference on adjacency sets
- **Lexicographical ordering** of all adjacency relations for deterministic processing

#### `MM2M` — Multi-Many-to-Many

Generalizes `M2M` to heterogeneous entity types. Stores a 2D array of `M2M` blocks indexed by `(type_i, type_j)`, enabling traversal across entity kinds (e.g., element-to-node, node-to-edge, face-to-element).

- Block-structured layout: `_blocks[typeCount, typeCount]`
- Type indices are resolved at compile time via the generic `ITypeMap` constraint
- `[Obsolete]` indexer for internal use (only `Topology<TTypes>` holds the lifecycle lock to prevent stale references)

#### `Topology<TTypes>` — Type-Safe Topology Container

The primary user-facing API. Combines `MM2M` adjacency storage with arbitrary per-entity attribute dictionaries under a single `ReaderWriterLockSlim`.

```csharp
// Define entity types
struct MyTypes : ITypeMap { ... }

// Build topology
var topo = new Topology<MyTypes>();
int n0 = topo.Add<Node>();
int n1 = topo.Add<Node>();
int e0 = topo.Add<Bar2, Node>(n0, n1);

// Attribute storage
topo.Set<Node, Position>(n0, new Position(x, y, z));
var pos = topo.Get<Node, Position>(n0);

// Sub-entity discovery
topo.DiscoverSubEntities<Tet4, Face, Node>(
    SubEntityDefinition.FromFaces((0,1,2), (0,1,3), (0,2,3), (1,2,3)));

// Graph algorithms
var order = topo.TopologicalOrder<Node, Bar2>();
var coloring = topo.ColorEntities<Node>();
```

**Supported graph algorithms:**
- Depth-First Search (DFS) and Breadth-First Search (BFS) on entity adjacency
- Topological ordering (Kahn's algorithm) for DAG-structured connectivity
- Entity coloring (greedy graph coloring) for parallel assembly scheduling

#### `Symmetry` — Canonical Symmetry Groups

Encodes permutation symmetry groups for element types, enabling automatic canonical representation and deduplication of equivalent elements.

- Represents symmetry as a set of permutations acting on local node indices
- Canonical form: lexicographically smallest permutation of the node index sequence
- Eliminates duplicated entries arising from orientation-insensitive operations (e.g., undirected edges, symmetric faces)

### Serialization & Comparison

All relation types support:
- Binary serialization/deserialization for mesh checkpointing
- `==`, `<`, `>` operators via lexicographic comparison of adjacency lists
- Set operations: intersection, union, symmetric difference

---

## Matrices — Linear Algebra

### `Matrix` — Dense Operations

Column-major dense matrix with SIMD-accelerated arithmetic and full decomposition suite. Optimized for element-level matrices up to ~1000×1000.

**Storage:** `double[]` in column-major (Fortran) order for BLAS/LAPACK compatibility:
```
Element [row, col] → _data[col * RowCount + row]
```

**Decompositions:**

| Method | Algorithm | Notes |
|---|---|---|
| `ComputeLU()` | Rook pivoting | Searches both pivot row and column; better stability than partial pivoting for near-singular matrices |
| `ComputeQR()` | Householder reflections | Numerically stable; produces thin QR |
| `ComputeSVD()` | One-sided Jacobi | Iterative; converges for any real matrix |
| `ComputeEigenvalues()` | Symmetric QR + Householder tridiagonalization | For symmetric matrices only; eigenvalues sorted descending |

**SIMD acceleration:** All arithmetic operators (`+`, `-`, `*`, scalar multiply) use `System.Runtime.Intrinsics` AVX2/AVX-512 vector lanes. Matrices above a size threshold are processed with `Parallel.For`.

**Factory methods:** `Identity(n)`, `Zero(m,n)`, `Diagonal(v)`, `Random(m,n,seed)`, `RandomNormal(m,n,mean,stddev)`

**Statistical functions:** `ColumnMeans()`, `ColumnStdDev()`, `Covariance()`, `Correlation()`

**Norms:** Frobenius, 1-norm (max column sum), ∞-norm (max row sum), max-norm

**Performance (50×50 matrix, Release x64):**

| Operation | Approx. time |
|---|---|
| Matrix multiply | ~0.15 ms |
| LU decomposition | ~0.28 ms |
| QR decomposition | ~0.65 ms |
| SVD | ~3.5 ms |
| Inverse | ~0.42 ms |

---

### `CSR` — Compressed Sparse Row

Production-grade sparse matrix for finite element applications. Supports structural-only matrices (sparsity pattern only) as well as valued matrices.

**Format:** Three arrays defining the non-zero structure:
```
rowPointers[i..i+1]  →  column indices and values for row i
columnIndices[k]     →  column of the k-th non-zero
values[k]            →  value of the k-th non-zero
```

**Backends (in priority order):**

| Backend | Condition | API |
|---|---|---|
| NVIDIA cuSPARSE (GPU) | `rows >= 50,000` and `nnz >= 1,000,000` | `CuSparseBackend` |
| Intel MKL PARDISO | MKL available | `PardisoSolver` |
| SIMD SpMV (CPU) | `rows >= 5,000` | `OptimizedSpMV` |
| Parallel SpMV (CPU) | `rows >= 1,000` | `Parallel.For` over rows |
| Sequential SpMV | fallback | standard row loop |

**`OptimizedSpMV`:** Hand-written AVX2 kernel processes 4 doubles per cycle. Handles row-aligned and tail remainder cases. Falls back to `System.Numerics.Vector<double>` when AVX2 is unavailable.

**`PardisoSolver`:** Wraps Intel MKL PARDISO via P/Invoke. Supports symmetric positive definite (type 2), symmetric indefinite (type -2), and general unsymmetric (type 11) matrices. Reuse factorization across multiple RHS vectors.

**`CuSparseBackend`:** Wraps `cusparseSpMV` and `cusolverSpDcsrlsvqr`. Automatically migrates matrix to device memory, runs the kernel, and retrieves results. GPU path activated only when problem size justifies transfer overhead.

**Iterative solvers (`CSRIterativeSolvers`):**
- **BiCGSTAB** — Bi-Conjugate Gradient Stabilized for general unsymmetric systems
- **GMRES** — Generalized Minimum Residual with restarts

**Object pooling:** `HashSet<int>` instances used during symbolic assembly are pooled via `Microsoft.Extensions.ObjectPool` to eliminate allocation pressure.

**Thresholds:**

| Constant | Value | Meaning |
|---|---|---|
| `MIN_ROWS_FOR_PARALLEL` | 1,000 | Minimum rows for parallel SpMV |
| `MIN_ROWS_FOR_SIMD` | 5,000 | Minimum rows for AVX2 SpMV |
| `MIN_ROWS_FOR_GPU` | 50,000 | Minimum rows for GPU consideration |
| `MIN_NNZ_FOR_GPU` | 1,000,000 | Minimum nnz for GPU consideration |
| `DEFAULT_TOLERANCE` | 1e-14 | Near-zero comparison tolerance |

---

### `CliqueSystem` — Finite Element Assembly

High-performance parallel finite element assembly using Gustavson's algorithm for symbolic factorization and lock-striped numerical assembly.

**Algorithm:**

1. **Symbolic assembly (C^T × C pattern):** Given element connectivity arrays, computes the global sparsity pattern using Gustavson's algorithm — the same technique used inside sparse direct solvers. This produces the exact CSR non-zero structure without storing values.

2. **DOF compression:** If the raw DOF space is more than 4× the actual DOF count (or exceeds 10M entries), a dictionary-based mapping is used instead of a dense array to avoid excessive memory allocation.

3. **Numerical assembly:** Each element contributes its local stiffness matrix `k_e` and force vector `f_e` to the global system. Lock-striped parallelism with **4096 stripes** (power-of-2 for fast bitwise modulo) prevents data races without serializing the assembly loop.

4. **Lock stripe hash:** Global DOF index → stripe index via `(uint)(dof * 0x9E3779B9) & 0xFFF`. The Fibonacci/golden-ratio multiplier `0x9E3779B9` ensures excellent distribution across stripes even for sequentially numbered DOFs.

5. **Large matrix support:** Internal storage is chunked at 256 MB per chunk (`33,554,432` doubles) to break the `Int32.MaxValue` array size limit for problems with >10M non-zeros.

**Key constants:**

| Constant | Value | Meaning |
|---|---|---|
| `LOCK_STRIPE_COUNT` | 4096 | Number of assembly lock stripes |
| `MIN_ELEMENTS_FOR_PARALLEL` | 100 | Minimum elements for parallel assembly |
| `MIN_DOFS_FOR_PARALLEL` | 10,000 | Minimum DOFs for parallel assembly |
| `MIN_DOFS_FOR_UNROLLED` | 8 | DOFs/element for unrolled inner loop |
| `MAX_DENSE_DOF_ARRAY_SIZE` | 10,000,000 | Switch to dict-based DOF map above this |

**Lifecycle:**

```csharp
var sys = new CliqueSystem(numElements, enableGpu: true);

// 1. Register element DOFs
for (int e = 0; e < numElements; e++)
    sys.SetElementDofs(e, dofIndices[e]);

// 2. Build sparsity pattern (symbolic phase)
sys.BuildSparsityPattern();

// 3. Assemble (numeric phase)
Parallel.For(0, numElements, e => {
    sys.AddElementMatrix(e, k_e[e]);
    sys.AddElementVector(e, f_e[e]);
});
sys.Assemble();

// 4. Apply boundary conditions and solve
sys.ApplyDirichletBC(fixedDofs, values);
var u = sys.Solve();

// 5. Reset for next load step (preserves sparsity pattern)
sys.Reset();
```

---

## Meshing — Generation & Refinement

### `SimplexMesh`

Core mesh container for 2D triangular and 3D tetrahedral meshes. Entity types:

| Type | Dimension | Nodes | Description |
|---|---|---|---|
| `Node` | 0D | 1 | Mesh vertex |
| `Point` | 0D | 1 | Single-node element |
| `Bar2` | 1D | 2 | Line segment |
| `Tri3` | 2D | 3 | Linear triangle |
| `Quad4` | 2D | 4 | Bilinear quadrilateral |
| `Tet4` | 3D | 4 | Linear tetrahedron |
| `Edge` | 1D | 2 | Topological edge (discovered) |

Coordinates are stored as `double[numNodes, 3]` regardless of problem dimension for API consistency.

### Structured Mesh Generators

#### Rectangular (2D)

```
CreateRectangularMesh(nx, ny, xMin, xMax, yMin, yMax)
```

Produces `(nx+1)*(ny+1)` nodes and `2*nx*ny` triangles. Each rectangular cell is split into two right triangles sharing the main diagonal:

```
Cell (i,j):   n00──n10      Triangle 1: (n00, n10, n11)
               │  ╲  │      Triangle 2: (n00, n11, n01)
              n01──n11
```

#### Box (3D)

```
CreateBoxMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax)
```

Produces `(nx+1)*(ny+1)*(nz+1)` nodes and `6*nx*ny*nz` tetrahedra. Each hexahedral cell is decomposed into 6 tetrahedra sharing the main diagonal:

```
Tet decomposition of hex cell (6 tets sharing vertex n000 and n111):
  Tet 1: (n000, n100, n110, n111)
  Tet 2: (n000, n110, n010, n111)
  Tet 3: (n000, n010, n011, n111)
  Tet 4: (n000, n011, n001, n111)
  Tet 5: (n000, n001, n101, n111)
  Tet 6: (n000, n101, n100, n111)
```

### `SimplexRemesher` — Conforming Refinement

Longest-edge bisection refinement, ensuring no element quality degradation:

1. **Edge selection:** Mark edges for refinement (longest edge, or user-specified)
2. **Midpoint insertion:** Add midpoint node at each marked edge; record parent nodes for solution transfer via the `ParentNodes` struct
3. **Element splitting:** Each element containing a refined edge is split into 2 (triangle) or 4 (tetrahedron) child elements
4. **Conforming closure:** Propagate refinement to neighboring elements to eliminate hanging nodes
5. **Parent tracking:** `ParentNodes` struct stores `(Parent1, Parent2)` for each new node; midpoint solution values computed as `u_mid = 0.5*(u_parent1 + u_parent2)`

**Canonical edge representation:** All edges stored as `(min(a,b), max(a,b))` tuples for unique identification in hash sets.

### `CrackOperations` — Level-Set Crack Insertion

Level-set based crack insertion for arbitrary 2D and 3D crack geometries:

1. **Level-set evaluation:** A signed-distance function `φ(x)` defines the crack geometry. Nodes on opposite sides of the crack front satisfy `φ(n_i) * φ(n_j) < 0`.

2. **Element classification:** Elements are classified as:
   - **Uncut:** entirely on one side of the crack
   - **Cut:** the crack front passes through the element interior
   - **Tip elements:** contain the crack tip/front

3. **Node duplication:** For cut elements, nodes on the positive-`φ` side are duplicated to produce independent crack faces. Connectivity is updated to reference the appropriate original or duplicated nodes.

4. **Crack-tip enrichment (XFEM):** The crack-tip region supports level-set enrichment functions for the near-tip singular stress field:
   ```
   {√r sin(θ/2),  √r cos(θ/2),  √r sin(θ/2)sin(θ),  √r cos(θ/2)sin(θ)}
   ```
   where `(r, θ)` are polar coordinates relative to the crack tip.

### Multi-Format I/O

| Format | Read | Write | Notes |
|---|---|---|---|
| VTK Legacy | ✓ | ✓ | ASCII `.vtk`, supports Tri3, Tet4, Quad4 |
| MSH (Gmsh) | ✓ | ✓ | v2 format |
| GiD/CIMNE `.msh` | — | ✓ | CIMNE GiD post-processing |
| Ensight Gold | — | ✓ | Binary/ASCII, multi-part |

---

## Nonlinear — Dynamics & Root Finding

### `BatheTwoStageIntegrator`

Unconditionally stable implicit time integrator for second-order dynamical systems:

```
M ü(t) + C u̇(t) + f_int(u(t), t) = R_ext(t)
```

**Bathe two-stage method** (Bathe 2007) splits each time step `Δt` into two sub-steps using different Newmark parameters:

| Stage | Sub-step | β | γ | Scheme |
|---|---|---|---|---|
| 1 | `Δt/2` | `1/4` | `1/2` | Trapezoidal rule (unconditionally stable) |
| 2 | `Δt/2` | `4/9` | `2/3` | Bathe corrector stage |

**Stage 1 Newmark predictor–corrector:**
```
ũ^(n+1/2) = u^n + (Δt/2)v^n + (Δt/2)²(1/2 - β₁)a^n
ṽ^(n+1/2) = v^n + (Δt/2)(1 - γ₁)a^n

Solve: [a₀M + a₁C + K_t] Δu = R_ext - f_int - M(a₀(u-ũ)) - C(a₁(u-ṽ))
  where a₀ = 1/(β₁(Δt/2)²),  a₁ = γ₁/(β₁(Δt/2))
```

**Stage 2 is analogous** with `(β₂, γ₂) = (4/9, 2/3)`, starting from the Stage 1 result.

**Newton-Raphson inner loop** with:
- Absolute residual tolerance: user-configurable `AbsTolerance`
- Relative residual tolerance: user-configurable `RelTolerance`
- Divergence detection: residual growth by `DivergenceThreshold` factor (default 1000×) within first 5 iterations flags divergence
- **Kahan compensated summation** for residual norm computation: eliminates floating-point cancellation error in long vectors, crucial for 10M+ DOF systems

**SIMD vectorization:** All vector operations (`VectorAdd`, `VectorNegate`, `ComputePredictor`) dispatch to AVX2 (8 doubles/cycle) or AVX-512 (8 doubles/cycle with ZMM registers) when detected at runtime via `Avx2.IsSupported` / `Avx512F.IsSupported`.

**Memory:** Zero-allocation hot paths — all working vectors are pre-allocated at construction. `ArrayPool<double>` used for temporary buffers in sub-routines.

**Initial static equilibrium:** Before time-stepping begins, the integrator solves the static problem `f_int(u₀, t₀) = R_ext(t₀)` with `v = a = 0` to find a consistent initial state.

---

### `RootFinder`

Scalar root-finding with two distinct API overloads:

#### With derivative — Hybrid Newton-Raphson + IQI

For functions providing both `f(x)` and `f'(x)`:

```csharp
var (root, status) = RootFinder.FindRoot(xmin, xmax, x => (f(x), df(x)));
```

Algorithm: Newton-Raphson steps when the Newton update stays within the bracket; falls back to Inverse Quadratic Interpolation (IQI) otherwise. Bisection as ultimate fallback.

#### Without derivative — ITP Algorithm

For functions providing only `f(x)`:

```csharp
var (root, status) = RootFinder.FindRoot(xmin, xmax, x => f(x));
```

**ITP (Interpolate-Truncate-Project)** algorithm (Oliveira & Takahashi, 2020) achieves optimal worst-case convergence of `O(log₂(1/ε))` iterations — matching bisection — while being superlinearly fast for smooth functions. Parameters:

| Parameter | Value | Role |
|---|---|---|
| `n₀` | 1 | Extra iterations over pure bisection |
| `κ_tr` | 0.2 | Truncation factor `κ ∈ (0, ∞)` |
| `p` | 2.0 | Super-linear convergence exponent |

**Convergence criteria and tolerances:**

| Constant | Value | Meaning |
|---|---|---|
| `FTOL` | 1e-10 | Absolute function value tolerance |
| `RTOL` | 1e-8 | Relative bracket width tolerance |
| `ATOL` | 1e-12 | Absolute bracket width tolerance |
| `MAX_ITER` | 100 | Maximum iterations |
| `EPS_MACHINE` | 2.22e-16 | IEEE 754 double epsilon |
| `GOLDEN_RATIO_COMPLEMENT` | 0.3820 | `(3 − √5)/2` for golden-section fallback |

**Status codes:** `OK`, `Tolerance`, `MaxIterations`, `NoBracket`, `BadInput`, `NonFinite`, `TooNarrow`

---

## Postprocess — Visualization

### `EnsightWriter`

Exports mesh and field data to Ensight Gold format, compatible with GiD/CIMNE and ParaView.

- **Multi-part case files:** aggregate many mesh/field pairs under a single `.case` descriptor
- Supports nodal scalar, vector, and tensor fields
- `EnsightWriter.WriteAllMeshes(caseName)` consolidates all registered meshes into one Ensight case file

---

## Project Structure

```
ManyToMany/
├── Numerical.sln                # Visual Studio solution
├── Relations/                   # Core topology library
│   ├── Relations.cs             #   O2M, M2M, MM2M implementations
│   ├── Topology.cs              #   Topology<TTypes>, SubEntityDefinition, Symmetry
│   └── Utils.cs                 #   Shared utilities, ParallelConfig
├── Matrices/                    # Linear algebra
│   ├── Matrix.cs                #   Dense matrix, decompositions, SIMD
│   ├── CSR.cs                   #   Sparse CSR, PARDISO, cuSPARSE, BiCGSTAB, GMRES
│   ├── Assembly.cs              #   CliqueSystem (Gustavson + lock-striped assembly)
│   └── NativeLibraries.cs       #   Cross-platform MKL/CUDA discovery & loading
├── Meshing/                     # Mesh generation & refinement
│   ├── SimplexMesh.cs           #   Core mesh container
│   ├── SimplexRemesher.cs       #   Structured generation, bisection refinement, I/O
│   ├── MeshRefinement.cs        #   Adaptive refinement drivers
│   ├── MeshGeometry.cs          #   Geometric primitives (level sets, distances)
│   └── CrackOperations.cs       #   Level-set crack insertion
├── Nonlinear/                   # Time integration & root finding
│   ├── Integrator.cs            #   BatheTwoStageIntegrator
│   └── RootFinder.cs            #   ITP + hybrid Newton-IQI algorithms
├── Postprocess/                 # Visualization
│   └── EnsightWriter.cs         #   Ensight Gold export
├── Teste/                       # 26 demo examples
│   └── Examples2DA.cs           #   Meshing + 2D/3D fracture mechanics
└── Docs/                        # Extended documentation
    ├── Numerical-Complete-Documentation.md
    ├── Topology-Complete-Documentation.md
    └── SimplexRemesher-Complete-Documentation.md
```

---

## Examples

The `Teste` project contains 26 examples across four parts.

### Part 1 — Advanced Meshing (Examples 1–10)

| # | Example | Features |
|---|---|---|
| 1 | Circular domain with eccentric hole | Curved boundary, interior void, Tri3/Quad4 |
| 2 | L-shape with corner refinement | Re-entrant corner, local refinement |
| 3 | Annulus region | Concentric boundaries, structured-to-unstructured |
| 4 | Wedge geometry | Non-convex boundary, mixed element types |
| 5 | Multiple holes | Multiple interior voids, conforming connectivity |
| 6 | Intricate boundary | High-curvature boundary discretization |
| 7 | Tri vs. Quad comparison | Side-by-side element type quality analysis |
| 8 | Cracked plate with hole | Combined hole and pre-existing crack |
| 9 | Gear-like geometry | Periodic boundary features, sharp re-entrant angles |
| 10 | Complex industrial shape | Multi-feature domain with 200+ boundary segments |

### Part 2 — Classical 2D Fracture Mechanics Benchmarks (Examples 11–15)

Each example meshes a reference domain, inserts a crack via level sets, and exports results for stress intensity factor (SIF) validation against published analytical solutions.

| # | Reference | Crack type | Analytical SIF |
|---|---|---|---|
| 11 | Anderson (2005) | Edge crack in tension | `K_I = σ√(πa) F(a/W)` |
| 12 | Griffith (1921) | Center crack in infinite plate | `K_I = σ√(πa)` |
| 13 | Kanninen & Popelar (1985) | Double edge notch | `K_I = σ√(πa) · 1.12` |
| 14 | Erdogan & Sih (1963) | Slant crack, mixed-mode | `K_I, K_II` as function of inclination angle |
| 15 | Newman & Raju (1984) | Crack emanating from circular hole | `K_I = σ√(πa) F(a/r, a/t)` |

### Part 3 — Spectacular Crack Patterns (Examples 16–20)

Artistic/scientific demonstrations of arbitrary crack geometries using the level-set insertion engine:

| # | Pattern | Crack geometry |
|---|---|---|
| 16 | Spiral galaxy | Archimedean spiral crack network |
| 17 | Fractal tree | Self-similar branching crack pattern |
| 18 | Sinusoidal waves | Periodic sinusoidal crack paths |
| 19 | Starburst | Radially symmetric crack fan |
| 20 | Concentric mandalas | Nested closed-loop cracks |

### Part 4 — 3D Fracture Mechanics Benchmarks (Examples 21–26)

| # | Reference | Crack type | Domain |
|---|---|---|---|
| 21 | Sneddon (1946) | Penny-shaped crack | Infinite solid under remote tension |
| 22 | Irwin (1962) | Elliptical crack | Infinite solid, semi-axes `(a, c)` |
| 23 | Tada (1973) | Edge crack | Semi-infinite body under bending |
| 24 | Newman & Raju (1981) | Corner crack at hole | Plate with through-thickness hole |
| 25 | Erdogan & Sih (1963) | Slant crack (3D) | General mixed-mode `K_I, K_II, K_III` |
| 26 | — | Semi-cylindrical surface crack | Pressurized cylindrical geometry |

All examples output:
- GiD/CIMNE `.msh` files for visualization in GiD (`www.gidhome.com`)
- A unified Ensight case file (`FractureMechanics.case`) for ParaView

---

## Performance

The library is engineered for high-throughput computational mechanics on modern x86-64 hardware.

### Hardware Acceleration

| Feature | Technology | Activation condition |
|---|---|---|
| Dense SIMD | AVX2 / AVX-512 intrinsics | Runtime `Avx2.IsSupported` / `Avx512F.IsSupported` |
| Sparse SpMV SIMD | `System.Numerics.Vector<double>` + AVX2 | `rows >= 5,000` |
| Sparse direct solver | Intel MKL PARDISO | MKL library discoverable |
| GPU SpMV/solve | NVIDIA cuSPARSE / cuSolver | `rows >= 50,000` AND `nnz >= 1,000,000` |
| Parallel CPU | `Parallel.For` with `ParallelOptions` | Problem size above per-operation thresholds |

### Memory Efficiency

- **`ArrayPool<double>`** for temporary buffers in hot paths (no GC pressure)
- **`GC.AllocateUninitializedArray`** for large pre-allocated arrays (skips zero-fill)
- **`CollectionsMarshal.AsSpan`** for zero-copy list access in parallel cloning
- **`[SkipLocalsInit]`** on performance-critical types to suppress redundant stack zeroing
- **Chunked storage** (`256 MB / chunk`) in `CliqueSystem` to exceed `Int32.MaxValue` element limit
- **Object pooling** (`Microsoft.Extensions.ObjectPool`) for `HashSet<int>` reuse in symbolic assembly
- **Stack-allocated `Span<T>`** for small temporary buffers in leaf routines

### Scalability

- Tested and designed for problems exceeding **10 million DOFs**
- `>2 GB` array support via chunked storage
- **JIT profile-guided optimization (PGO)** via `[AggressiveOptimization]` method attribute on hottest loops
- Lock-striped parallel assembly with 4096 stripes scales to 32+ hardware threads with minimal contention

### Native Library Discovery (`NativeLibraries.cs`)

Cross-platform automatic discovery of native accelerators:

| Library | Windows path hints | Linux path hints | macOS path hints |
|---|---|---|---|
| Intel MKL | `%MKLROOT%\redist`, Program Files | `/opt/intel/oneapi/mkl`, `/usr/lib` | `/opt/intel/oneapi/mkl` |
| CUDA Runtime | `%CUDA_PATH%\bin` | `/usr/local/cuda/lib64` | — |
| cuSPARSE | bundled with CUDA | bundled with CUDA | — |

Libraries are loaded lazily on first use. If a library is absent, the corresponding backend silently degrades to the next available option.

---

## Building & Running

### Build

```bash
dotnet build Numerical.sln -c Release -p:Platform=64
```

### Run all 26 examples

```bash
dotnet run --project Teste -c Release -p:Platform=64
```

Output files are written to the working directory:
- `ex1_tri.msh` … `ex10_quad.msh` — GiD mesh files for Part 1
- `griffith1921_tri.msh`, `anderson2005_tri.msh`, … — GiD mesh files for Parts 2–4
- `FractureMechanics.case` + companion files — unified Ensight case for ParaView

### Run a single example

```csharp
// In Teste/Examples2DA.cs, call any static method directly:
Examples2DA.Example12_Griffith1921_CenterCrack();
```

---

## Platforms & Prerequisites

### Supported Platforms

| OS | Architecture |
|---|---|
| Windows 10/11 | x64 |
| Linux (Ubuntu 20.04+) | x64 |
| macOS 12+ | x64, ARM64 |

### Requirements

- [.NET 9.0 SDK](https://dotnet.microsoft.com/download/dotnet/9.0) (x64)

### Optional Accelerators

| Accelerator | Purpose | Installation |
|---|---|---|
| Intel MKL / oneAPI | PARDISO direct solver | Windows: auto via NuGet · Linux: `sudo apt-get install intel-mkl` · macOS: `brew install intel-mkl` |
| CUDA Toolkit 11.0+ | GPU SpMV, cuSolver | From [developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads) |

Without optional accelerators the library runs entirely on CPU using managed SIMD and parallel algorithms.

---

## Documentation

Extended documentation is available in the `Docs/` directory:

| File | Contents |
|---|---|
| `Numerical-Complete-Documentation.md` | Dense/sparse matrices, FE assembly, native library integration, full API reference |
| `Topology-Complete-Documentation.md` | Topology operations, graph algorithms, serialization, full API reference |
| `SimplexRemesher-Complete-Documentation.md` | Mesh refinement algorithms, crack insertion, file I/O, tutorial examples |

---

## License

GPLv3

## Author

Pedro Areias (IST)

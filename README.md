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
| `Relations` | `O2M`, `Topology<TTypes>`, `Symmetry` | Mesh topology, graph algorithms |
| `Matrices` | `Matrix`, `CSR`, `CliqueSystem` | Dense/sparse algebra, FE assembly |
| `Meshing` | `SimplexMesh`, `SimplexRemesher`, `MeshGeometry` | Mesh generation, refinement, fracture |
| `Nonlinear` | `BatheTwoStageIntegrator`, `RootFinder` | Time integration, root finding |
| `Postprocess` | `EnsightWriter` | Visualization export |

---

## Relations — Topology & Connectivity

The `Relations` library is the foundation of ManyToMany. It provides type-safe, high-performance many-to-many relationship structures for representing mesh connectivity, adjacency graphs, and arbitrary entity associations.

### Relationship Hierarchy

#### `O2M` — One-to-Many

A sparse adjacency list mapping one entity to an ordered list of entities. Single-threaded by design; thread safety is provided at the `Topology<TTypes>` level.

- Internal storage: `List<List<int>>` with pre-allocated capacity
- Parallel cloning via `GC.AllocateUninitializedArray` + `Parallel.For` above a configurable parallelization threshold
- Deep-copy uses `CollectionsMarshal.AsSpan` for zero-allocation list mirroring
- `[SkipLocalsInit]` attribute eliminates redundant zero-fills in hot paths
- Implements `IComparable<O2M>`, `IEquatable<O2M>`, `ICloneable`

#### `Topology<TTypes>` — Type-Safe Topology Container

The primary user-facing API. Combines multi-type adjacency storage with arbitrary per-entity attribute dictionaries under a single `ReaderWriterLockSlim`.

```csharp
// Define entity types via TypeMap
using Types = TypeMap<Node, Edge, Tri3>;

// Build topology
var topo = new Topology<Types>();
int n0 = topo.Add<Node>();
int n1 = topo.Add<Node>();
int n2 = topo.Add<Node>();
int e0 = topo.Add<Tri3, Node>(n0, n1, n2);

// Attribute storage
topo.Set<Node, Position>(n0, new Position(x, y, z));
var pos = topo.Get<Node, Position>(n0);

// Sub-entity discovery
topo.DiscoverSubEntities<Tet4, Face, Node>(
    SubEntityDefinition.FromFaces((0,1,2), (0,1,3), (0,2,3), (1,2,3)));

// Graph algorithms
var order = topo.GetTopologicalOrder<Node, Tri3>();
var coloring = topo.ComputeElementColoring<Tri3, Node>();
```

**Supported graph algorithms:**
- Breadth-First Search (BFS) on entity adjacency, including multi-type BFS
- Multi-type Depth-First Search (DFS) across heterogeneous entity types
- Topological ordering (Kahn's algorithm) for DAG-structured connectivity
- Element coloring (greedy graph coloring) for parallel assembly scheduling
- Connected component detection
- Cuthill-McKee bandwidth reduction ordering
- Transitive connectivity computation and dual structure extraction

#### `Symmetry` — Canonical Symmetry Groups

Encodes permutation symmetry groups for element types, enabling automatic canonical representation and deduplication of equivalent elements.

- Represents symmetry as a set of permutations acting on local node indices
- Canonical form: lexicographically smallest permutation of the node index sequence
- Eliminates duplicated entries arising from orientation-insensitive operations (e.g., undirected edges, symmetric faces)

### Serialization & Comparison

All relation types support:
- JSON serialization/deserialization via `ToJson()` / `FromJson()` and `SaveToFile()` / `LoadFromFile()`
- `==`, `<`, `>` operators via lexicographic comparison of adjacency lists (on `O2M`)
- Set operations: intersection, union, difference (via `Topology<TTypes>`)
- Structural validation and integrity checks

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
| NVIDIA cuSPARSE (GPU) | `rows >= 50,000` and `nnz >= 1,000,000` | GPU sparse backend |
| Intel MKL PARDISO | MKL available | Direct sparse solver |
| SIMD SpMV (CPU) | `rows >= 5,000` | AVX2 vectorized SpMV |
| Parallel SpMV (CPU) | `rows >= 1,000` | `Parallel.For` over rows |
| Sequential SpMV | fallback | standard row loop |

**SIMD SpMV:** Hand-written AVX2 kernel processes 4 doubles per cycle. Handles row-aligned and tail remainder cases. Falls back to `System.Numerics.Vector<double>` when AVX2 is unavailable.

**PARDISO backend:** Wraps Intel MKL PARDISO via P/Invoke. Supports symmetric positive definite (type 2), symmetric indefinite (type -2), and general unsymmetric (type 11) matrices. Reuse factorization across multiple RHS vectors.

**GPU backend:** Wraps `cusparseSpMV` and `cusolverSpDcsrlsvqr`. Automatically migrates matrix to device memory, runs the kernel, and retrieves results. GPU path activated only when problem size justifies transfer overhead.

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

### Crack Insertion (via `SimplexRemesher`)

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
│   ├── Relations.cs             #   O2M and supporting adjacency structures
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
│   └── MeshGeometry.cs          #   Geometric primitives (level sets, distances)
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

## Public API Reference

This section lists every public type in each module. With autocomplete, these are all the entry points a user needs.

---

### Relations

#### Core relationship types

| Type | Kind | Description |
|---|---|---|
| `O2M` | `sealed class` | One-to-many sparse adjacency list. Implements `IComparable<O2M>`, `IEquatable<O2M>`, `ICloneable`. Internal storage is `List<List<int>>`; parallel cloning via `GC.AllocateUninitializedArray` above a configurable threshold. |

#### `Topology<TTypes>` — primary user-facing API

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
| `int Add<TElement, TNode, TData>(TData data, params int[] nodes)` | Add an element with connectivity and attribute data. |
| `(int Index, bool WasNew) AddUnique<TElement, TNode>(params int[] nodes)` | Add an element only if canonically unique (requires symmetry). |
| `int[] AddRange<TNode, TData>(IEnumerable<TData> dataItems)` | Batch-add nodes with data. |
| `int[] AddRange<TElement, TNode>(IEnumerable<int[]> connectivityList)` | Batch-add elements. |
| `int[] AddRangeParallel<TElement, TNode>(int[][] connectivityList, ...)` | Parallel batch-add for large element sets. |
| `void Remove<TEntity>(int index)` | Mark an entity for removal. |
| `void RemoveRange<TEntity>(IEnumerable<int> indices)` | Mark multiple entities for removal. |
| `int Count<TEntity>()` | Total number of entities of the given type (including marked). |
| `int CountActive<TEntity>()` | Number of non-deleted entities. |
| `List<int> GetActive<TEntity>()` | Indices of all non-deleted entities. |
| `bool Exists<TElement>(params int[] nodes)` | Check whether an element with the given nodes exists. |
| `int Find<TElement>(params int[] nodes)` | Find an element by its nodes; returns index or -1. |

**Attribute storage**

| Method | Description |
|---|---|
| `void Set<TEntity, TData>(int index, TData value)` | Attach attribute data to an entity. |
| `void SetRange<TEntity, TData>(int startIndex, ReadOnlySpan<TData> values)` | Batch-set attributes for a contiguous range. |
| `TData Get<TEntity, TData>(int index)` | Retrieve attribute data from an entity. |
| `bool TryGet<TEntity, TData>(int index, out TData value)` | Non-throwing attribute retrieval. |
| `IEnumerable<(int Index, TData Data)> Each<TEntity, TData>()` | Iterate over all entities with their data. |
| `IEnumerable<int> Each<TEntity>()` | Iterate over all entity indices of a type. |
| `void ForEach<TEntity, TData>(Action<int, TData> action)` | Execute an action on each entity with data. |
| `void ParallelForEach<TEntity, TData>(Action<int, TData> action, ...)` | Parallel iteration over entities with data. |

**Adjacency queries**

| Method | Description |
|---|---|
| `void WithNodesOf<TElement, TNode>(int element, Action<ReadOnlySpan<int>> action)` | Access an element's nodes via a zero-allocation span callback. |
| `IEnumerable<int> EnumerateNeighbors<TElement, TNode>(int element)` | Enumerate nodes connected to an element. |
| `List<int> GetElementsWithNodes<TElement, TNode>(List<int> nodes)` | Get elements connected to all specified nodes (intersection). |
| `List<int> GetElementsContainingAnyNode<TElement, TNode>(List<int> nodes)` | Get elements connected to any specified node (union). |
| `List<int> GetElementsFromNodes<TElement, TNode>(List<int> nodes)` | Get elements whose node sets are subsets of the given nodes. |
| `List<int> ElementsContainingAllNodes<TElement, TNode>(params int[] nodes)` | Find elements containing all specified nodes. |
| `int CountRelated<TEntity, TRelated>(int entityIndex)` | Count of related entities. |
| `int CountIncident<TElement, TNode>(int nodeIndex)` | Count of elements incident to a node. |
| `List<int> GetDirectNeighbors<TEntity, TRelated>(int entityIndex, ...)` | Direct neighbors via shared related entities. |
| `List<(int EntityIndex, int SharedCount)> GetWeightedNeighbors<TEntity, TRelated>(int entityIndex, ...)` | Neighbors with shared-entity counts. |
| `List<int> GetElementNeighbors<TElement, TNode>(int element, ...)` | Element-to-element neighbors via shared nodes. |
| `List<int> GetNodeNeighbors<TElement, TNode>(int node, ...)` | Node-to-node neighbors via shared elements. |
| `Dictionary<int, int> GetKHopNeighborhood<TEntity, TRelated>(int entityIndex, int maxHops)` | K-hop neighborhood with distances. |

**Sub-entity discovery & boundary**

| Method | Description |
|---|---|
| `void DiscoverSubEntities<TElement, TSubEntity, TNode>(SubEntityDefinition def)` | Enumerate and register all sub-entities (faces, edges) implied by element connectivity. |
| `Topology<TTypes> WithSymmetry<TElement>(Symmetry symmetry)` | Register a symmetry group for canonical deduplication; returns `this` for chaining. |
| `List<int> ElementsSharingSubEntity<TParent, TSubEntity, TNode>(int subEntityIndex)` | Find parent elements sharing a sub-entity. |
| `int CountElementsSharingSubEntity<TParent, TSubEntity, TNode>(int subEntityIndex)` | Count parent elements sharing a sub-entity. |
| `int GetLocalNodeIndex<TElement, TNode>(int element, int node)` | Get the local position of a node within an element. |
| `List<int[]> ExtractBoundaryFacets<TElement, TNode>(int nodesPerBoundaryFacet)` | Extract boundary faces/facets. |
| `List<(int[] Nodes, int, int)> FindInternalFacets<TElement, TNode>(int nodesPerFacet)` | Find interior faces shared by two elements. |
| `List<int> GetBoundarySubEntities<TParent, TSubEntity, TNode>()` | Get boundary sub-entity indices. |
| `List<int> GetInteriorSubEntities<TParent, TSubEntity, TNode>()` | Get interior sub-entity indices. |
| `bool IsSubEntityOnBoundary<TParent, TSubEntity, TNode>(int subEntityIndex)` | Check if a sub-entity is on the boundary. |
| `List<int> DetectNonManifoldSubEntities<TParent, TSubEntity, TNode>()` | Find non-manifold sub-entities. |
| `SubEntityBoundaryResult DetectSubEntityBoundary<TParent, TSubEntity, TNode>()` | Full boundary detection returning classified sub-entities. |

**Graph algorithms**

| Method | Description |
|---|---|
| `List<int> GetTopologicalOrder<TEntity>()` | Topological sort within a single entity type. |
| `List<int> GetTopologicalOrder<TEntity, TRelated>()` | Topological sort over a bipartite entity–related DAG. |
| `List<int> GetSortOrder<TEntity>()` | Lexicographic sort order of entity adjacency. |
| `int[] ComputeElementColoring<TElement, TNode>()` | Greedy graph coloring; returns color index per element. |
| `List<List<int>> GetColorGroups<TElement, TNode>()` | Group elements by assigned color. |
| `List<int> BreadthFirstSearch<TEntity>(int startEntity, ...)` | BFS from a starting entity with optional visitor callback. |
| `Dictionary<int, int> BreadthFirstDistances<TEntity>(int startEntity)` | BFS distances from a starting entity. |
| `List<(int, int)> BreadthFirstSearchMultiType<TStartEntity>(int startEntity)` | BFS across heterogeneous entity types. |
| `List<(int, int)> MultiTypeDFS<TNode>(int nodeIndex)` | DFS across all types starting from a node. |
| `bool IsAcyclic<TEntity, TRelated>()` | Check if the entity–related graph is acyclic. |
| `IReadOnlyList<IReadOnlyList<int>> FindComponents<TElement, TNode>()` | Find connected components. |
| `List<List<int>> FindConnectedComponents<TEntity, TRelated>(int minShared)` | Connected components with minimum shared count. |
| `int[] ComputeCuthillMcKeeOrdering<TElement, TNode>(bool reverse = true)` | Cuthill-McKee bandwidth reduction ordering. |
| `int ComputeBandwidth<TElement, TNode>()` | Compute the bandwidth of the adjacency. |
| `O2M ComputeTransitiveConnectivity<TEntity, TRelated>()` | Transitive closure of the entity graph. |
| `O2M GetDualStructure<TEntity, TRelated>()` | Dual graph as O2M. |
| `O2M GetElementToElementGraph<TElement, TNode>()` | Element-to-element adjacency graph. |
| `O2M GetNodeToNodeGraph<TElement, TNode>()` | Node-to-node adjacency graph. |

**Serialization & validation**

| Method | Description |
|---|---|
| `string ToJson(JsonSerializerOptions? options)` | Serialize topology to JSON string. |
| `static Topology<TTypes> FromJson(string json, ...)` | Reconstruct from a JSON string. |
| `void SaveToFile(string path, ...)` | Save topology to a JSON file. |
| `static Topology<TTypes> LoadFromFile(string path, ...)` | Load topology from a JSON file. |
| `TopologyStats GetStatistics()` | Summary of entity counts and adjacency density. |
| `bool ValidateStructure()` | Structural consistency check. |
| `ValidationResult ValidateIntegrity<TElement, TNode>()` | Detailed integrity validation for a relationship. |
| `List<int> GetDuplicates<TEntity>()` | Find duplicate entities of a type. |
| `bool IsPermutationOf<TElement, TNode>(Topology<TTypes> other)` | Check if two topologies are permutations of each other. |

**Memory & lifecycle**

| Method | Description |
|---|---|
| `void Compress(bool removeDuplicates, bool shrinkMemory, bool validate)` | Remove deleted entities and optionally shrink memory. |
| `void Clear()` | Remove all entities and data. |
| `Topology<TTypes> Clone()` | Deep copy. |
| `void Reserve<TElement, TNode>(int capacity)` | Pre-allocate capacity for a relationship. |
| `void ShrinkToFit()` | Release excess memory. |
| `void ConfigureType<TEntity>(int parallelizationThreshold, int? reserveCapacity)` | Configure parallelization threshold for a type. |
| `void WithBatch(Action action)` | Execute multiple operations atomically under a single write lock. |
| `int Merge<TElement, TNode>(Topology<TTypes> other)` | Merge another topology's entities into this one. |

**Transpose & internal structure**

| Method | Description |
|---|---|
| `O2M GetTranspose<TElement, TNode>()` | Get the transposed adjacency (node-to-elements). |
| `O2M GetTranspose<TElement, TNode>(int maxNodeCap)` | Transposed adjacency with a maximum node capacity. |
| `O2M GetTransposeStrict<TElement, TNode>()` | Strict transpose that throws on invalid node indices. |
| `void WithTranspose<TElement, TNode>(Action<O2M> action)` | Access transpose via callback with automatic locking. |
| `void WithElementsForNodeSpan<TElement, TNode>(int nodeIndex, ReadOnlySpanAction<int> action)` | Zero-allocation access to elements incident to a node. |
| `void EnsureSynchronized<TElement, TNode>()` | Force transpose cache synchronization. |
| `void EnsurePositionCaches<TElement, TNode>()` | Force position cache rebuild. |

**Smart entity wrapper**

| Type | Description |
|---|---|
| `SmartEntity<TEntity>` | `readonly record struct` pairing a `Topology` reference with an entity index. Provides `IsValid`, `IsMarked`, `Count`, `Data<T>()`, `SetData<T>()`, `BreadthFirstSearch()`, and `BreadthFirstDistances()`. |

**Nested result types**

| Type | Kind | Description |
|---|---|---|
| `ConnectivityStatistics` | `sealed class` | Entity counts, average/max/min adjacency degree. |
| `ValidationResult` | `readonly struct` | Pass/fail + error description from `ValidateIntegrity()`. |
| `SubEntityBoundaryResult` | `readonly struct` | Output of boundary sub-entity detection. |
| `ColoringStatistics` | `readonly struct` | Number of colors used, distribution per color. |
| `DualGraph` | `sealed class` | Element-to-element adjacency graph derived from shared nodes, with BFS. |
| `ElementStatistics` | `readonly struct` | Per-element quality and connectivity metrics. |

#### `ReadOnlyTopology<TTypes>`

```csharp
public sealed class ReadOnlyTopology<TTypes> where TTypes : ITypeMap, new()
```

Read-only projection of a `Topology<TTypes>`. Exposes the same query methods but disallows mutation. Useful for passing topology to solvers that must not modify connectivity.

#### `Symmetry`

```csharp
public sealed class Symmetry
```

Encodes a permutation symmetry group for an element type. Used with `Topology.WithSymmetry<TEntity>()` to canonicalize and deduplicate entries (e.g., undirected edges, symmetric faces).

| Factory | Description |
|---|---|
| `Symmetry.Identity(int nodeCount)` | Identity symmetry (no permutations beyond identity). |
| `Symmetry.Cyclic(int n)` | Cyclic group C_n (n rotations). |
| `Symmetry.Dihedral(int n)` | Dihedral group D_n (n rotations + n reflections = 2n elements). |
| `Symmetry.Full(int n)` | Full symmetric group S_n (all n! permutations, max n=8). |
| `Symmetry.FromGenerators(List<List<int>> generators)` | Build group from generating permutations. |

#### `SubEntityDefinition`

```csharp
public readonly struct SubEntityDefinition
```

Describes which local node indices form each sub-entity (face, edge) of a parent element.

| Factory | Description |
|---|---|
| `SubEntityDefinition.FromEdges(params (int, int)[] edges)` | Edges (2-node sub-entities). |
| `SubEntityDefinition.FromFaces(params (int, int, int)[] faces)` | Triangular faces (3-node sub-entities). |
| `SubEntityDefinition.FromQuadFaces(params (int, int, int, int)[] faces)` | Quadrilateral faces (4-node sub-entities). |

#### `ResultOrder` enum

Controls the ordering of results returned by graph algorithms: `Unordered` (insertion order, fastest) and `Sorted` (deterministic ordering by type and entity index).

#### Type-mapping infrastructure

`ITypeMap` is the compile-time interface that maps C# types to integer indices within a `Topology`. Users implement it by using one of the pre-built `TypeMap<T0, …, Tn>` generic classes (provided for 2 up to 25 type arguments):

```csharp
// Example: 3-type map
TypeMap<Node, Edge, Tri3>  // Node→0, Edge→1, Tri3→2
```

#### Serialization DTOs

Plain data objects used for JSON export/import of topology state:

`TopologyDto`, `AdjacencyDto`, `DataListDto`, `SymmetryDto`, `CanonicalIndexDto`, `CanonicalEntryDto`.

#### Utility types

| Type | Kind | Description |
|---|---|---|
| `Topology` (non-generic) | `static class` | Factory and extension methods for building topologies from common patterns. |
| `TopologyStats` | `sealed class` | Aggregate statistics snapshot. |
| `Utils` | `static class` | Shared utility methods (sorting, span helpers, etc.). |
| `ParallelConfig` | `static class` | Global parallelization thresholds (`ProcessorCount`, `DefaultParallelThreshold`). |
| `ListComparer<T>` | `sealed class` | `IComparer<List<T>>` for lexicographic list ordering. |
| `ListEqualityComparer<T>` | `sealed class` | `IEqualityComparer<List<T>>` for set operations on adjacency lists. |

---

### Matrices

#### `Matrix` — dense matrix

```csharp
public sealed class Matrix : IEquatable<Matrix>, IFormattable, ICloneable
```

Column-major (`_data[col * RowCount + row]`) dense matrix with SIMD-accelerated arithmetic (AVX2/AVX-512) and a full decomposition suite.

**Factory methods**

| Method | Description |
|---|---|
| `Matrix.Identity(int n)` | n×n identity. |
| `Matrix.Zeros(int rows, int cols)` | Zero matrix. |
| `Matrix.Ones(int rows, int cols)` | All-ones matrix. |
| `Matrix.Diagonal(params double[] values)` | Diagonal matrix from values. |
| `Matrix.Random(int rows, int cols, int? seed)` | Uniform random entries. |
| `Matrix.RandomNormal(int rows, int cols, double mean, double stdDev, int? seed)` | Gaussian random entries. |
| `Matrix.MatrixSquareRoot(Matrix A)` | Principal square root via Schur decomposition. |
| `Matrix.MatrixExponential(Matrix A, int order)` | Matrix exponential via Padé approximant. |

**Decompositions**

| Method | Returns | Algorithm |
|---|---|---|
| `ComputeLU()` | `LUDecomposition` | Rook pivoting (better stability than partial for near-singular). |
| `ComputeQR()` | `QRDecomposition` | Householder reflections. |
| `ComputeEigenvalues()` | `EigenDecomposition` | Symmetric QR + Householder tridiagonalization (symmetric matrices). |
| `ComputeSVD(int? seed)` | `SVDDecomposition` | One-sided Jacobi (any real matrix). |

**Key instance methods**

| Method | Description |
|---|---|
| `Transpose()` / `TransposeInPlace()` | Matrix transpose. |
| `Inverse()` | Full matrix inverse via LU. |
| `Determinant()` | Scalar determinant. |
| `Solve(Vector b)` | Solve Ax = b (LU). |
| `SolveMultiple(Matrix B)` | Solve AX = B. |
| `SolveLeastSquares(Vector b)` | Minimum-norm least-squares via QR. |
| `Rank(double tol)` | Numerical rank. |
| `ConditionNumber()` | 2-norm condition number. |
| `FrobeniusNorm()` / `OneNorm()` / `InfinityNorm()` | Matrix norms. |
| `IsSymmetric()` / `IsDiagonal()` / `IsUpperTriangular()` / `IsLowerTriangular()` | Structure predicates. |
| `GetRow(int)` / `GetColumn(int)` / `GetSubMatrix(...)` | Slicing. |
| `HorizontalConcat(Matrix)` / `VerticalConcat(Matrix)` | Augmentation. |
| `KroneckerProduct(Matrix)` | Kronecker product. |
| `PseudoInverse(double tol)` | Moore-Penrose pseudoinverse via SVD. |
| `GetNullSpace()` / `GetRowSpace()` / `GetImageSpace()` | Fundamental subspaces. |
| `ColumnMeans()` / `RowMeans()` / `Covariance()` / `Correlation()` | Statistics. |
| `Hadamard(Matrix)` | Element-wise product. |
| `Apply(Func<double,double>)` / `Map(Func<int,int,double,double>)` | Element-wise transforms. |
| `static Multiply(A, B)` / `MultiplyAtB(A,B)` / `MultiplyABt(A,B)` / `MultiplyAtBt(A,B)` | Named multiply variants (SIMD-accelerated). |

**Operators:** `+`, `-`, `*` (matrix-matrix and scalar), unary `-`, `==`, `!=`, implicit/explicit casts to/from `Vector`.

#### `Vector`

```csharp
public sealed class Vector : IEquatable<Vector>, IFormattable, ICloneable
```

Dense vector companion to `Matrix`.

| Method | Description |
|---|---|
| `Vector.Zeros(int n)` / `Vector.Ones(int n)` / `Vector.Random(int n, int? seed)` | Factory methods. |
| `Dot(Vector)` | Inner product. |
| `Norm()` / `Norm1()` / `NormInf()` | Vector norms. |
| `Normalize()` | Unit vector. |
| `Cross(Vector)` | 3D cross product. |
| `OuterProduct(Vector)` | Returns `Matrix`. |
| `ProjectOnto(Vector)` | Orthogonal projection. |
| `AngleTo(Vector)` / `DistanceTo(Vector)` | Geometric helpers. |
| `ElementwiseMultiply(Vector)` / `ElementwiseDivide(Vector)` | Component-wise ops. |
| `Slice(int start, int length)` | Sub-vector. |
| `Map(Func<double,double>)` | Element-wise transform. |
| `Sum()` / `Mean()` / `Max()` / `Min()` / `ArgMax()` / `ArgMin()` | Reductions. |

#### Decomposition result types

| Type | Key members |
|---|---|
| `LUDecomposition` | `Solve(Vector b)`, `Determinant()`, `ConditionNumber()` |
| `QRDecomposition` | `Solve(Vector b)`, `Rank(double tol)` |
| `EigenDecomposition` | `Eigenvalues` (`double[]`), `Eigenvectors` (`Matrix`) |
| `SVDDecomposition` | `U`, `S`, `Vt` (`Matrix`); `Rank(double tol)`, `ConditionNumber()` |

#### `CSR` — compressed sparse row

```csharp
public sealed class CSR : IFormattable, IEquatable<CSR>, ICloneable, IDisposable
```

Production-grade sparse matrix. Automatically selects the best available backend: NVIDIA cuSPARSE (GPU), Intel MKL PARDISO, AVX2 SpMV, parallel SpMV, or sequential fallback.

Key methods: `Multiply(double[] x, double[] y)` (SpMV), `Solve(double[] rhs)`, `Clone()`, `ToArray()`, plus constants `MIN_ROWS_FOR_PARALLEL`, `MIN_ROWS_FOR_SIMD`, `MIN_ROWS_FOR_GPU`, `MIN_NNZ_FOR_GPU`, `DEFAULT_TOLERANCE`.

#### `CSRIterativeSolvers`

```csharp
public static class CSRIterativeSolvers
```

| Method | Algorithm |
|---|---|
| `BiCGSTAB(CSR A, double[] b, ...)` | Bi-Conjugate Gradient Stabilized — general unsymmetric systems. |
| `GMRES(CSR A, double[] b, ...)` | Generalized Minimum Residual with restarts. |

Both return a `SolverResult` record.

#### `SolverResult` / `MatrixStatistics`

```csharp
public record SolverResult(double[] Solution, int Iterations, double ResidualNorm, bool Converged);
public record MatrixStatistics(int Rows, int Columns, long NonZeroCount, double SparsityRatio, ...);
```

#### `SolverException`

```csharp
public class SolverException : Exception
```

Thrown when a solver fails to converge or encounters a singular system.

#### `CliqueSystem` — finite element assembly

```csharp
public sealed class CliqueSystem : IDisposable
```

High-performance parallel FE assembly using Gustavson's algorithm for symbolic factorization and lock-striped (4096 stripes) numerical assembly. Supports problems exceeding 10M DOFs via 256 MB chunked storage.

**Lifecycle**

| Step | Method | Description |
|---|---|---|
| 1 | `CliqueSystem(int numElements, bool enableGpu = false)` | Constructor. |
| 2a | `SetElementSize(int e, int numDofs)` | Declare per-element DOF count before building structure. |
| 2b | `SetElementConnectivity(int e, int[] globalDofs)` | Assign global DOF indices to element `e`. |
| 2c | `BuildSparsityPattern()` | Symbolic phase — computes the CSR sparsity pattern. |
| 3 | `AddElement(int e, double[] force, double[] stiffness)` | Accumulate element matrices (thread-safe, call from `Parallel.For`). |
| 4 | `Assemble()` | Finalize global matrix and force vector. |
| 5 | `Solve()` | Returns `double[]` displacement vector. |
| 6 | `Reset()` | Clear values, preserve sparsity pattern for next load step. |

**Other methods**

| Method | Description |
|---|---|
| `static CliqueSystem FromTopology<TTypes,TElement,TNode>(...)` | Construct directly from a `Topology` instance. |
| `GetMatrix()` | Returns assembled `CSR`. |
| `GetForceVector()` | Returns assembled `double[]` RHS. |
| `GetStatistics()` | Returns `AssemblyStatistics`. |
| `GetSystemInfo()` | Human-readable summary string. |

#### `AssemblyStatistics`

```csharp
public sealed class AssemblyStatistics
```

Carries `TotalDofs`, `NonZeroCount`, `SparsityRatio`, and timing fields populated after `Assemble()`.

#### `DiscreteLinearSystem`

```csharp
public class DiscreteLinearSystem : IDisposable
```

Higher-level wrapper that holds an assembled `CSR` system and exposes `Solve(double[,] result)`, `BuildSystemValues()`, `Reset()`, and boundary-condition application. Useful when DOF numbering follows node-major layout rather than pure DOF indices.

#### Native library helpers

| Type | Kind | Description |
|---|---|---|
| `INativeLibraryConfig` | `interface` | Contract for specifying custom search paths for MKL/CUDA. |
| `NativeLibraryConfig` | `sealed class` | Default implementation of `INativeLibraryConfig`. |
| `LibraryStatus` | `class` | Reports load status, version, and path for one native library. |
| `LibraryAvailability` | `static class` | `IsMklAvailable`, `IsCudaAvailable`, `GetMklStatus()`, etc. |
| `NativeLibraryStatus` | `static class` | Runtime summary of all discovered accelerators. |
| `NuGetLibraryChecker` | `static class` | Detects NuGet-distributed native packages (e.g., `MathNet.Numerics.MKL`). |
| `RobustNativeLibraryLoader` | `static class` | Cross-platform `LoadLibrary`/`dlopen` wrapper with fallback search. |
| `SparseBackendFactory` | `static class` | Selects the optimal sparse backend given problem size. |
| `HybridScheduler` | `sealed class` | Runtime dispatcher between CPU and GPU backends. Exposes `BackendType` and `OperationType` enums. |

---

### Meshing

#### Element-type marker structs

Zero-size `readonly struct` types used exclusively as generic type arguments in `Topology<TTypes>` and `SimplexMesh`:

| Type | Dimension | Nodes |
|---|---|---|
| `Node` | 0D | 1 — mesh vertex |
| `Edge` | 1D | 2 — topological edge (discovered automatically) |
| `Point` | 0D | 1 — single-node element |
| `Bar2` | 1D | 2 — line segment |
| `Tri3` | 2D | 3 — linear triangle |
| `Quad4` | 2D | 4 — bilinear quadrilateral |
| `Tet4` | 3D | 4 — linear tetrahedron |

#### `ParentNodes` / `OriginalElement`

```csharp
public readonly record struct ParentNodes(int Parent1, int Parent2);
public readonly record struct OriginalElement(int Index);
```

Attribute types attached to nodes/elements during refinement. `ParentNodes` stores the two parent node indices of a midpoint node; `OriginalElement` stores the pre-refinement element index.

#### `SimplexMesh`

```csharp
public sealed class SimplexMesh : Topology<TypeMap<Node, Edge, Point, Bar2, Tri3, Quad4, Tet4>>
```

The main mesh container. Inherits all `Topology<TTypes>` methods and adds mesh-specific helpers:

| Method | Description |
|---|---|
| `AddNode(int parentIndex)` | Add a vertex node; sets `ParentNodes` to itself. |
| `AddMidpointNode(int p1, int p2)` | Add a refinement midpoint; records both parents. |
| `AddTriangle(int n0, int n1, int n2)` | Add a `Tri3` element. |
| `AddQuad(int n0, int n1, int n2, int n3)` | Add a `Quad4` element. |
| `AddTetrahedron(int n0, int n1, int n2, int n3)` | Add a `Tet4` element. |
| `AddBar(int n0, int n1)` | Add a `Bar2` element. |

#### `MeshConstants`

```csharp
public static class MeshConstants
```

Shared numerical tolerances: `Epsilon` (1e-10), `DegenerateAreaTolerance` (1e-14), `DegenerateVolumeTolerance` (1e-15), `NodeMergeTolerance` (1e-12), `GridPerturbationFactor` (0.15), `HexRowSpacing` (√3/2).

#### `FiniteElementTopologies`

```csharp
public static class FiniteElementTopologies
```

Pre-built `SubEntityDefinition` instances for standard element types (e.g., `Tri3Edges`, `Tet4Faces`, `Tet4Edges`). Pass directly to `SimplexMesh.DiscoverSubEntities<...>()`.

#### `SimplexRemesher`

```csharp
public static class SimplexRemesher
```

Structured mesh generators, conforming bisection refinement, and multi-format I/O.

**Mesh generation**

| Method | Description |
|---|---|
| `CreateRectangularMesh(nx, ny, xMin, xMax, yMin, yMax)` | `2·nx·ny` triangles over a rectangle; returns `(SimplexMesh mesh, double[,] coords)`. |
| `CreateBoxMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax)` | `6·nx·ny·nz` tetrahedra over a box; returns `(SimplexMesh mesh, double[,] coords)`. |

**Refinement**

| Method | Description |
|---|---|
| `InterpolateCoordinates(SimplexMesh refined, double[,] originalCoords)` | Transfer coordinate array to a refined mesh using `ParentNodes`. |
| `DiscoverEdges(SimplexMesh mesh)` | Populate all `Edge` entities from existing elements. |

**File I/O**

| Method | Format |
|---|---|
| `SaveMSH(mesh, coords, path)` | Gmsh v2 `.msh` |
| `SaveMSHWithCrackGroups(mesh, coords, path, ...)` | Gmsh v2 with named crack groups |
| `SaveASCII(mesh, coords, path)` | Plain text |
| `SaveGiD(mesh, coords, path)` | GiD/CIMNE `.msh` |
| `FindBoundaryNodes(SimplexMesh)` | Returns `HashSet<int>` of boundary node indices (2D). |
| `FindBoundaryNodes3D(SimplexMesh)` | Returns `HashSet<int>` of boundary node indices (3D). |
| `PrintStats(SimplexMesh, string label)` | Print element/node counts to console. |

**Smoothing**

| Method | Description |
|---|---|
| `LaplacianSmoothing(mesh, coords, iterations)` | Classic Laplacian node relocation. |
| `CVTSmoothing(mesh, coords, iterations)` | Centroidal Voronoi tessellation smoothing. |

#### `MeshRefinement`

```csharp
public static class MeshRefinement
```

Adaptive refinement drivers that operate on an existing `SimplexMesh`:

| Method | Description |
|---|---|
| `CheckJacobians(mesh, coords, label)` | Count and report negative-Jacobian elements; returns count. |
| `CheckEdgeTopology(mesh, label)` | Validate edge-to-element connectivity. |
| `FixNegativeJacobians(mesh, coords)` | Flip element connectivity to correct negative Jacobians. |
| `InterpolateCoordinates(mesh, originalCoords)` | Coordinate transfer after refinement. |

#### `MeshGeometry`

```csharp
public static class MeshGeometry
```

Geometric primitives for mesh construction and quality analysis.

**Element geometry**

| Method | Description |
|---|---|
| `ComputeTriangleJacobian(coords, n0, n1, n2)` | Signed 2× triangle area. |
| `ComputeTriangleArea(coords, n0, n1, n2)` | Unsigned area. |
| `ComputeTriangleAspectRatio(coords, n0, n1, n2)` | Circumradius / inradius ratio. |
| `ComputeTriangleMinAngle(coords, n0, n1, n2)` | Minimum interior angle (radians). |
| `IsTriangleCCW / IsTriangleDegenerate` | Orientation and degeneracy tests. |
| `ComputeTetrahedronJacobian / Volume / AspectRatio` | 3D analogues. |
| `IsTetrahedronCorrectOrientation / IsTetrahedronDegenerate` | 3D tests. |
| `ComputeQuadArea / IsQuadCCW / IsQuadConvex` | Quadrilateral helpers. |
| `EdgeLength2D / EdgeLength3D` | Edge length from coordinate array. |

**Point/curve geometry**

| Method | Description |
|---|---|
| `IsPointInPolygon(point, polygon)` | Ray-casting containment test (array or list overload). |
| `IsPointOnPolygonBoundary(point, polygon, tol)` | Boundary proximity test. |
| `DistancePointToLine / DistancePointToSegment` | Point-to-line/segment distances. |
| `Distance2D / Distance3D` | Euclidean distances. |
| `AngleBetweenVectors2D / SignedAngleBetweenVectors2D` | Vector angles. |
| `CrossProduct2D / DotProduct2D` | 2D vector operations. |
| `RefineBoundary(boundary, maxEdgeRatio, minEdgeCount)` | Refine a boundary polygon. |
| `ResampleCurve / ResampleCurveArray` | Resample a curve to a target point count. |
| `ComputeArcLength / ComputeAverageEdgeLength / ComputeSignedArea` | Curve metrics. |
| `EnsureCCW / EnsureCW / ReverseBoundary / IsBoundaryCCW` | Boundary orientation. |

**Mesh-level utilities**

| Method | Description |
|---|---|
| `ValidateMeshOrientation(mesh, coords, ...)` | Count mis-oriented elements. |
| `ComputeQualityStatistics(mesh, coords)` | Returns `MeshQualityStats`. |
| `GenerateInteriorGrid(...)` | Generate a perturbed interior point cloud for mesh seeding. |
| `PrintQualityReport / ExportToConsole` | Console diagnostics. |
| `IdentifyBoundaryNodes(mesh)` | Alternative boundary-node finder. |

**Delegate type**

```csharp
public delegate double SignedFieldFunction(double x, double y, double z);
```

Used throughout `MeshGeometry` and `SimplexRemesher` to define level-set functions for crack insertion.

#### `MeshQualityStats`

```csharp
public class MeshQualityStats
```

Aggregated quality metrics returned by `MeshGeometry.ComputeQualityStatistics()`: min/max/mean aspect ratio, min/max/mean minimum angle, count of degenerate elements.

---

### Nonlinear

#### `BatheTwoStageIntegrator`

```csharp
public sealed class BatheTwoStageIntegrator
```

Unconditionally stable implicit time integrator for second-order dynamical systems `M ü + C u̇ + f_int(u, t) = R_ext(t)`. Uses the Bathe two-stage Newmark method (Stage 1: β=¼, γ=½; Stage 2: β=4/9, γ=2/3).

**Delegates (set before stepping)**

```csharp
public delegate void ResidualEvaluator(double t, double[] u, double[] v, double[] a, double[] residual);
public delegate void EffectiveSystemSolver(double[] effectiveForce, double[] displacement);
```

**Configuration properties**

| Property | Description |
|---|---|
| `int Dimension` | System size (number of DOFs). |
| `double AbsTolerance` | Absolute residual convergence tolerance. |
| `double RelTolerance` | Relative residual convergence tolerance. |
| `int MaxNewtonIterations` | Maximum Newton-Raphson iterations per step. |
| `double DivergenceThreshold` | Residual growth factor to flag divergence (default 1000). |
| `double Time` | Current simulation time. |

**Stepping**

| Method | Description |
|---|---|
| `Step(double dt)` | Advance one time step of size `dt`. Solves initial static equilibrium on first call. |
| `Step(double dt, int numSteps)` | Advance `numSteps` time steps. |

**State access**

| Method | Description |
|---|---|
| `GetDisplacement(Span<double> dest)` | Copy current displacement into `dest`. |
| `GetVelocity(Span<double> dest)` | Copy current velocity. |
| `GetAcceleration(Span<double> dest)` | Copy current acceleration. |
| `GetState(Span<double> u, Span<double> v, Span<double> a)` | Copy all three in one call. |
| `ResetPerformanceCounters()` | Zero cumulative counters. |

**Nested types**

| Type | Description |
|---|---|
| `ConvergenceInfo` | Per-step Newton data: `Iterations`, `InitialResidualNorm`, `FinalResidualNorm`, `MaxResidualNorm`, `Converged`. |
| `PerformanceCounters` | Cumulative totals: total Newton iterations, total steps, average iterations/step. |

#### `RootFinder`

```csharp
public static class RootFinder
```

Scalar root-finding with two overloads, both thread-safe (no static mutable state).

| Overload | Algorithm | Signature |
|---|---|---|
| With derivative | Hybrid Newton-Raphson + IQI, bisection fallback | `FindRoot(xmin, xmax, Func<double, (double f, double df)>)` |
| Without derivative | ITP (Interpolate-Truncate-Project) | `FindRoot(xmin, xmax, Func<double, double>)` |

Both return `(double root, RootFinder.Status status)`.

**`RootFinder.Status` enum**

`OK`, `Tolerance`, `MaxIterations`, `NoBracket`, `BadInput`, `NonFinite`, `TooNarrow`.

#### `TrustRegionNewtonDogleg`

```csharp
public static class TrustRegionNewtonDogleg
```

Trust-region Newton dogleg method for scalar or vector nonlinear systems. Configured via `TRNOptions` and returns `TRNResult`.

#### `TRNOptions` / `TRNResult`

```csharp
public readonly record struct TRNOptions(double InitialRadius, double MaxRadius, double Tolerance, int MaxIterations, ...);
public readonly record struct TRNResult(double Root, int Iterations, double Residual, RootFinder.Status Status);
```

---

### Postprocess

#### `EnsightWriter`

```csharp
public static class EnsightWriter
```

Exports mesh and field data to Ensight Gold format (ASCII), compatible with ParaView and GiD.

| Method | Description |
|---|---|
| `AddMesh(string name, SimplexMesh mesh, double[,] coords, double[,]? displacement)` | Register a mesh in the internal collection. |
| `WriteAllMeshes(string basename)` | Write all registered meshes as a multi-timestep Ensight case (`.case` + `_XXXX.geo` files). |
| `SaveEnsight(SimplexMesh mesh, double[,] coords, string basename)` | Write a single mesh directly to Ensight format. |
| `SaveEnsightWithScalar(SimplexMesh mesh, double[,] coords, string basename, ...)` | Write mesh plus a nodal scalar field. |

Output consists of a `.case` descriptor and per-timestep `.geo` geometry files. Displacement fields are written as `.CrackOpening` vector files and can be scaled inside ParaView/GiD.

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**.

You are free to use, modify, and distribute this software under the terms of the GPLv3. See the [LICENSE](LICENSE) file for the full license text.

Copyright (C) 2026 Pedro Miguel de Almeida Areias

## Author

Pedro Areias (IST)

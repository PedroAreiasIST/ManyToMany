![ManyToMany](m2m.png)

# ManyToMany

A high-performance scientific computing and finite element analysis library for .NET 9.0, written in C# by Pedro Areias (IST).

ManyToMany provides a unified framework for managing complex mesh topologies, sparse linear algebra, mesh generation with crack insertion, nonlinear time integration, and post-processing — all tuned for large-scale computational mechanics on x64 hardware.

---

## Table of Contents

- [Features](#features)
- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Building](#building)
- [Module Reference](#module-reference)
  - [Relations — Topology & Connectivity](#relations--topology--connectivity)
  - [Matrices — Linear Algebra](#matrices--linear-algebra)
  - [Meshing — Generation & Refinement](#meshing--generation--refinement)
  - [Nonlinear — Dynamics & Root Finding](#nonlinear--dynamics--root-finding)
  - [Postprocess — Visualization](#postprocess--visualization)
- [Examples](#examples)
- [Performance](#performance)
- [Architecture & Design](#architecture--design)
- [Platforms](#platforms)
- [Documentation](#documentation)
- [License](#license)

---

## Features

| Module | Highlights |
|--------|-----------|
| **Relations** | Type-safe many-to-many mesh topology; O2M/M2M/MM2M hierarchy; graph algorithms (BFS, DFS, Dijkstra); JSON serialization; symmetry groups |
| **Matrices** | Dense (BLAS/LAPACK-compatible) and CSR sparse matrices; SIMD (AVX2/AVX-512); GPU via CUDA/cuSPARSE; PARDISO direct solver; LU/QR/SVD/eigenvalues; FEM assembly engine |
| **Meshing** | 2D/3D simplex mesh generation; longest-edge bisection refinement; level-set crack insertion; GiD/CIMNE and Ensight I/O |
| **Nonlinear** | Bathe two-stage implicit integrator (unconditionally stable, 2nd-order); Newton-Raphson, Bisection, Secant, Brent root finders |
| **Postprocess** | Ensight 6.0 ASCII export for ParaView/GiD visualization with multi-step time series support |

---

## Project Structure

```
ManyToMany/
├── Numerical.sln            # Solution file
├── Relations/               # Topology & connectivity (core library)
│   ├── Topology.cs          # Topology<TTypes>, Symmetry (~10 000 lines)
│   ├── Relations.cs         # O2M, M2M, MM2M (~9 300 lines)
│   └── Utils.cs             # ITypeMap, TypeMap<T0..T5>, ParallelConfig (~3 600 lines)
├── Matrices/                # Dense/sparse linear algebra
│   ├── CSR.cs               # Compressed sparse row matrix
│   ├── Matrix.cs            # Dense matrix
│   ├── Assembly.cs          # CliqueSystem FEM assembler
│   └── Solvers.cs           # PARDISO, cuSPARSE, iterative solvers
├── Meshing/                 # Mesh generation, refinement, crack insertion
├── Nonlinear/               # Time integration & root finding
├── Postprocess/             # Visualization export
├── Teste/                   # 26 worked examples
└── Docs/                    # API reference documentation
```

**Codebase size:** ~49 600 lines across 15 source files.

---

## Prerequisites

- [.NET 9.0 SDK](https://dotnet.microsoft.com/download/dotnet/9.0), x64 only
- **Optional — Intel MKL** (for PARDISO direct solver):
  - Windows: resolved automatically via NuGet
  - Linux: `sudo apt-get install intel-mkl`
  - macOS: `brew install intel-mkl`
- **Optional — CUDA 11 or 12** (for cuSPARSE GPU acceleration, NVIDIA GPU required)

---

## Building

```bash
# Release build (x64 optimised)
dotnet build Numerical.sln -c Release -p:Platform=64

# Debug build
dotnet build Numerical.sln -c Debug -p:Platform=64
```

The project enables aggressive JIT settings in Release mode: tiered PGO, AVX2/AVX-512 instruction sets, no overflow checks, and ReadyToRun compilation.

---

## Module Reference

### Relations — Topology & Connectivity

The Relations module is the architectural core of the library. It provides type-safe, high-performance many-to-many relationship management for mesh entities, backed by a three-level hierarchy: **O2M → M2M → MM2M → Topology**.

#### Relationship Hierarchy

| Class | Thread-safe | Description |
|-------|-------------|-------------|
| `O2M` | No | One-to-Many — sparse adjacency lists (element → nodes) |
| `M2M` | Yes | Many-to-Many — thread-safe O2M with cached transpose and position lookups |
| `MM2M` | Yes | Multi-type Many-to-Many — 2D array of M2M structures indexed by `[elementType, nodeType]` |
| `Topology<TTypes>` | Yes | Generic typed container combining MM2M with per-entity data storage |

#### O2M — One-to-Many

O2M stores element-to-node connectivity as a list of adjacency lists. It is the low-level workhorse — single-threaded by design for maximum throughput.

```csharp
var o2m = new O2M();
o2m.Add(new[] { 0, 1, 2 });          // element 0 connects to nodes 0,1,2
o2m.Add(new[] { 1, 2, 3 });          // element 1 connects to nodes 1,2,3

O2M transpose = o2m.Transpose();      // node → elements (O(n), parallel-sort)
O2M strict    = o2m.TransposeStrict(); // throws if negative node indices found

// Graph queries
var bfsOrder  = o2m.BreadthFirstSearch(startElement: 0, transpose);
var distances = o2m.BreadthFirstDistances(startElement: 0, transpose);
var paths     = o2m.DijkstraShortestPaths(0, (from, to, node) => 1.0);

// Graph extraction
O2M elemGraph = o2m.GetElementToElementGraph(transpose); // adjacency via shared nodes
O2M nodeGraph = o2m.GetNodeToNodeGraph();                // adjacency via shared elements
```

Key details:
- `TransposeSkipsInvalidNodes = true` (default) — negative node indices are ignored in transpose (useful for sentinel values)
- `ParallelizationThreshold` — configures when parallel processing engages (default 1 000)
- Supports `IComparable<O2M>` and `IEquatable<O2M>` for set operations

#### M2M — Many-to-Many (thread-safe)

M2M wraps O2M with a `ReaderWriterLockSlim` and adds lazily-computed transpose caching and position caches. Batch operations defer synchronization until the batch exits.

```csharp
var m2m = new M2M();
m2m.Add(nodeIndex: 5, elementConnectivity: new[] { 0, 1, 5 });

// Zero-copy access to transpose (no allocation)
m2m.WithElementsFromNode(transpose => {
    var elements = transpose[nodeIndex: 5];
});

// Batch multiple writes — lock held for entire duration
m2m.WithBatch(() => {
    for (int i = 0; i < 1000; i++) m2m.Add(i, ...);
});

m2m.Compress(); // compact storage — WARNING: invalidates any cached M2M references
```

> **Warning (P0.4):** After `Compress()`, any previously stored M2M references become stale and will throw `ObjectDisposedException`. Use `WithBlock()` instead of the direct indexer `[]` for safe access.

#### MM2M — Multi-type Many-to-Many

MM2M manages a `[numberOfTypes × numberOfTypes]` grid of M2M structures for heterogeneous entity types.

```csharp
var mm2m = new MM2M(numberOfTypes: 3);

// Safe access via WithBlock (survives Compress)
mm2m.WithBlock(elementType: 0, nodeType: 1, m2m => {
    m2m.Add(nodeIndex, connectivity);
});

// Cross-type BFS
var visited = mm2m.BreadthFirstSearchMultiType(startElement, startType);
```

#### Topology\<TTypes\> — Generic Typed Container

`Topology<TTypes>` is the primary API surface. It combines MM2M with per-entity typed data storage and compile-time type safety via `ITypeMap`.

```csharp
// Define entity types at compile time
class MyMesh : Topology<TypeMap<Node, Edge, Triangle, Quad>> { }

var mesh = new MyMesh();

// Add entities
int tri  = mesh.Add<Triangle, Node>(new[] { 0, 1, 2 });
int quad = mesh.Add<Quad, Node>(new[] { 0, 1, 2, 3 });
int n    = mesh.Add<Node, NodeData>(new NodeData(x: 1.0, y: 2.0));

// Batch insert for performance
mesh.WithBatch(() => {
    for (int i = 0; i < 100_000; i++) mesh.Add<Triangle, Node>(...);
});

// Query connectivity
IReadOnlyList<int> nodes    = mesh.NodesOf<Triangle, Node>(tri);
IReadOnlyList<int> tris     = mesh.ElementsAt<Triangle, Node>(nodeIndex: 0);
IReadOnlyList<int> allTris  = mesh.ElementsAtAll<Triangle, Node>(new[] { 0, 1 }); // at ALL given nodes
IReadOnlyList<int> anyTris  = mesh.ElementsAtAny<Triangle, Node>(new[] { 0, 1 }); // at ANY given node
IReadOnlyList<int> neighbors = mesh.Neighbors<Triangle, Node>(tri);

// Iteration
mesh.ForEach<Triangle, TriData>((index, data) => { ... });
mesh.ParallelForEach<Node, NodeData>((index, data) => { ... });

// Structural
mesh.Remove<Triangle>(tri);
mesh.Compress(removeDuplicates: true, shrinkMemory: true, validate: true);
var copy = mesh.Clone();
var sub  = mesh.CloneWhere<Triangle, Node>(t => keepCondition(t));

// Symmetry-aware deduplication
var sym = Symmetry.Cyclic(3);           // triangle rotation symmetry
mesh.WithSymmetry<Triangle>(sym);       // deduplicate by canonical form
mesh.AddUnique<Triangle, Node>(...);    // inserts only if not symmetrically equivalent

// JSON persistence
string json = mesh.ToJson();
mesh.SaveToFile("mesh.json");
var loaded = MyMesh.LoadFromFile("mesh.json");

// Read-only view
IReadOnlyTopology<TypeMap<...>> ro = mesh.AsReadOnly();
```

**Graph algorithms on Topology:**

```csharp
var bfsOrder  = mesh.BreadthFirstSearch<Triangle, Node>(startElement: 0);
var distances = mesh.BreadthFirstDistances<Triangle, Node>(startElement: 0);
var paths     = mesh.DijkstraShortestPaths<Triangle, Node>(0,
                    (from, to, sharedNode) => edgeWeight(from, to));
```

#### ITypeMap and TypeMap Variants

`ITypeMap` provides the compile-time type-to-integer mapping that underpins `Topology<TTypes>`. It is resolved at JIT time with zero runtime overhead.

```csharp
public interface ITypeMap {
    int Count { get; }
    int IndexOf<T>();
    bool TryIndexOf<T>(out int index);
}

// Implementations provided for 2 to 25 type parameters:
TypeMap<T0, T1>
TypeMap<T0, T1, T2>
TypeMap<T0, T1, T2, T3>
// ... up to TypeMap<T0, ..., T24>
```

#### Symmetry

`Symmetry` represents a permutation group and computes canonical forms for automatic element deduplication.

```csharp
// Cyclic symmetry: [0,1,2], [1,2,0], [2,0,1] all equivalent
var cyclic  = Symmetry.Cyclic(3);

// Dihedral symmetry: rotations + reflections
var dihedral = Symmetry.Dihedral(4);

// Custom symmetry
var sym = new Symmetry(new List<List<int>> {
    new() { 0, 1, 2 },   // identity
    new() { 2, 0, 1 },   // rotation
    new() { 1, 2, 0 },   // rotation
});

// Canonical form — lexicographically smallest permutation
List<int> canon = sym.Canonical(new[] { 5, 3, 7 }); // → [3, 5, 7]

// Equivalence check
bool eq = sym.AreEquivalent(new[] { 0, 1, 2 }, new[] { 2, 0, 1 }); // true

// For large arrays: span-based zero-allocation variant
sym.CanonicalSpan(nodes.AsSpan(), destination.AsSpan());
```

#### Graph Algorithms

| Algorithm | Method | Complexity | Notes |
|-----------|--------|-----------|-------|
| BFS | `BreadthFirstSearch()` | O(V+E) with pre-computed transpose | Returns elements in discovery order; optional visitor callback |
| BFS distances | `BreadthFirstDistances()` | O(V+E) | Hop count to each element |
| DFS | `DepthFirstSearchFromANode()` | O(V+E) | Three-coloring (white/gray/black); cycle detection |
| Dijkstra | `DijkstraShortestPaths()` | O((V+E) log V) | Edge weight `Func<from, to, sharedNode, double>` |
| Element graph | `GetElementToElementGraph()` | O(V·E) | Elements as vertices, shared-node edges |
| Node graph | `GetNodeToNodeGraph()` | O(V·E) | Nodes as vertices, shared-element edges |
| Multi-type BFS | `BreadthFirstSearchMultiType()` | O(V+E) | Cross-entity-type traversal on MM2M |

---

### Matrices — Linear Algebra

#### CSR — Compressed Sparse Row

High-performance sparse matrix with automatic SIMD and GPU acceleration.

```csharp
// Construct from row/column/value arrays
var csr = new CSR(rowPointers, columnIndices, nCols, values);

// Element access
double v = csr[row, col];
csr.Set(row, col, 3.14);
csr.AddToElement(row, col, delta);   // thread-safe with object-level lock

// Sparse operations
double[] y = csr.SparseMV(x);        // y = A·x  (SIMD/GPU accelerated)
CSR      C = csr.SparseMMOptimized(B); // C = A·B
CSR      T = csr.Transpose();

// Direct solvers
double[] x = csr.SolvePardiso(b, matrixType: 11);  // Intel MKL PARDISO
double[] x = csr.SolvePardiso(B, nrhs: 5, ...);    // multiple RHS

// Iterative solvers
var result = CSRIterativeSolvers.DiagonallyPreconditionedBiCGSTAB(csr, b,
                 tolerance: 1e-10, maxIterations: 1000);

// Decompositions
var lu  = csr.LU();
var qr  = csr.QR();
var svd = csr.SVD();

// Diagnostics
double norm  = csr.Norm();      // Frobenius norm
double trace = csr.Trace();
```

**Automatic acceleration thresholds:**

| Condition | Acceleration |
|-----------|-------------|
| ≥ 1 000 rows | Parallel CPU processing |
| ≥ 5 000 rows | SIMD vectorization (AVX2/AVX-512) |
| ≥ 50 000 rows **and** ≥ 1 000 000 non-zeros | GPU via cuSPARSE |

**GPU acceleration** is resolved dynamically at runtime. The library probes for `cudart64_12.dll` / `libcudart.so.12` (CUDA 12) then `_11` variants on both Windows and Linux. If absent, it silently falls back to CPU.

#### Matrix — Dense

Column-major dense matrix compatible with BLAS/LAPACK conventions.

```csharp
var A = new Matrix(4, 4);           // zero-filled
var B = Matrix.Identity(4);
var C = Matrix.Random(100, 100);

// Decompositions
var lu   = A.LU();                  // Rook-pivoted LU
var qr   = A.QR();                  // Householder QR
var svd  = A.SVD();                 // Jacobi SVD
var eig  = A.Eigenvalues();         // QR algorithm (symmetric matrices)

// Solving
double[] x  = A.Solve(b);           // via LU
Matrix   X  = A.SolveMultiple(B);   // multiple RHS
double[] xl = A.SolveLeastSquares(b); // via QR

// Operations
Matrix D = A * B;
Matrix E = A + B;
double d = A.Determinant();
Matrix I = A.Inverse();
double t = A.Trace();
double n = A.Norm();
```

> **Performance note:** Iterate column-first (`for j, for i, A[i,j]`) to exploit column-major storage and avoid cache misses.

**Representative timings (50×50, AVX2+FMA, modern x64):**

| Operation | Time |
|-----------|------|
| Matrix multiply | ~0.15 ms |
| LU decomposition | ~0.28 ms |
| QR decomposition | ~0.65 ms |
| SVD | ~3.5 ms |
| Matrix inverse | ~0.42 ms |

#### CliqueSystem — FEM Assembly Engine

High-throughput finite element global system assembler. Operates directly on DOFs — no node abstraction.

```csharp
var sys = new CliqueSystem(numElements: 10_000, enableGpu: true);

// Symbolic phase: specify connectivity
for (int e = 0; e < numElements; e++)
    sys.SetElementDofs(e, localDofs[e]);   // int[]

sys.DetermineSystemSize();                 // computes total DOF count

// Numeric phase: fill local matrices / vectors
for (int e = 0; e < numElements; e++) {
    sys.SetElementMatrix(e, Ke[e]);        // double[ndof, ndof]
    sys.SetElementVector(e, fe[e]);        // double[ndof]
}

// Assemble and solve
sys.Assemble();                           // Gustavson C^T×C, lock-striped parallel
double[] u = sys.Solve();                 // PARDISO or cuSPARSE

// Incremental reassembly (same sparsity pattern, new values)
sys.Reset();
// ... refill element matrices ...
sys.Assemble();
double[] u2 = sys.Solve();

// Diagnostics
var stats = sys.Statistics;
Console.WriteLine($"DOFs: {stats.TotalDofs}, nnz: {stats.NonZeroCount}");
```

The assembler uses 4 096 lock stripes with golden-ratio hashing to minimise contention across threads during parallel numeric assembly.

---

### Meshing — Generation & Refinement

`SimplexMesh` inherits from `Topology<TypeMap<Node, Edge, Point, Bar2, Tri3, Quad4, Tet4>>`. All mesh operations go through `SimplexRemesher` (a static utility class).

#### Mesh Generation

```csharp
// 2D structured meshes
SimplexMesh tri  = SimplexRemesher.CreateRectangularMesh(nx: 10, ny: 10,
                       xMin: 0, xMax: 1, yMin: 0, yMax: 1);
SimplexMesh quad = SimplexRemesher.CreateRectangularQuadMesh(nx: 10, ny: 10, ...);

// 3D structured mesh
SimplexMesh box  = SimplexRemesher.CreateBoxMesh(nx: 5, ny: 5, nz: 5, ...);

// Constrained Delaunay triangulation (2D)
SimplexMesh mesh = SimplexRemesher.Triangulate(boundary: boundaryPoints,
                       convertToQuads: false, enableSmoothing: true);

// With interior holes
SimplexMesh mesh = SimplexRemesher.TriangulateWithHoles(
                       outer: outerBoundary,
                       holes: new[] { hole1, hole2 });
```

#### Mesh Refinement

```csharp
// Mark edges for bisection
int[] markedEdges = SelectEdgesNearCrack(mesh);

// Refine (longest-edge bisection; no hanging nodes)
SimplexMesh refined = SimplexRemesher.Refine(mesh, markedEdges);

// Reposition new midpoint nodes
double[,] newCoords = MeshRefinement.InterpolateCoordinates(refined, originalCoords);
```

#### Crack Insertion

```csharp
// Insert a level-set crack path into an existing mesh
// Two coincident nodes are placed at each crack location,
// enabling crack-opening visualization
SimplexMesh cracked = SimplexRemesher.InsertCrack(mesh, coordinates, crackPath);
```

#### Quad Conversion

```csharp
// Convert triangular mesh to mixed tri/quad mesh (greedy pairing)
SimplexMesh quadMesh = SimplexRemesher.ConvertToQuads(triMesh);
```

#### Mesh I/O

```csharp
// GiD/CIMNE format (.msh)
SimplexRemesher.SaveGiD(mesh, coordinates, "output.msh");
SimplexMesh loaded = SimplexRemesher.LoadGiD("input.msh");

// Ensight 6.0 format (single mesh)
EnsightWriter.SaveEnsight(mesh, coordinates, "output");
// → output.case + output.geo

// Ensight time series (multiple meshes / deformed states)
EnsightWriter.AddMesh("Step0", mesh, coords0);
EnsightWriter.AddMesh("Step1", mesh, coords1, displacement: delta1);
EnsightWriter.AddMesh("Step2", mesh, coords2, displacement: delta2);
EnsightWriter.WriteAllMeshes("result");
// → result.case + result_0.geo + result_1.geo + result_2.geo
```

---

### Nonlinear — Dynamics & Root Finding

#### BatheTwoStageIntegrator

Unconditionally stable implicit time integrator for second-order structural dynamics:

$$M\ddot{u}(t) + C\dot{u}(t) + f_\text{int}(u,t) = R_\text{ext}(t)$$

The method applies two half-steps per time increment:
- **Stage 1** (t → t + Δt/2): Trapezoidal rule (β = 1/4, γ = 1/2)
- **Stage 2** (t + Δt/2 → t + Δt): Bathe composite (β = 4/9, γ = 2/3)

```csharp
var integrator = new BatheTwoStageIntegrator {
    Dimension             = nDofs,
    MaxNewtonIterations   = 20,
    AbsTolerance          = 1e-10,
    RelTolerance          = 1e-8,
    ThrowOnDivergence     = true,
    DivergenceThreshold   = 1e6,
};

integrator.Initialize(
    u0: initialDisplacement,
    v0: initialVelocity,
    systemSolver: (t, beta, gamma, u, v, a, rhs) => {
        // Assemble and solve: [K + beta*Dt²*M] Δu = rhs
        return SolveTangentSystem(t, beta, gamma, u, v, a, rhs);
    },
    residualEvaluator: (t, u, v, a, residual) => {
        // Fill residual vector: r = M·a + C·v + f_int(u) - R_ext(t)
        ComputeResidual(t, u, v, a, residual);
    }
);

// Solve for static equilibrium before time marching
integrator.SolveInitialStaticEquilibrium();

// Time march
for (double t = 0; t < T_end; t += dt)
    integrator.Step(dt);
```

The integrator uses Kahan compensated summation throughout to avoid accumulation of floating-point errors over long time histories.

#### RootFinder

```csharp
double root = RootFinder.NewtonRaphson(f, df, x0: 1.0, tol: 1e-12);
double root = RootFinder.Bisection(f, a: 0.0, b: 2.0, tol: 1e-10);
double root = RootFinder.Secant(f, x0: 0.5, x1: 1.5, tol: 1e-10);
double root = RootFinder.Brent(f, a: 0.0, b: 2.0, tol: 1e-12);
```

---

### Postprocess — Visualization

```csharp
// Single-step export (Ensight 6.0 ASCII — open in ParaView or GiD)
EnsightWriter.SaveEnsight(mesh, coordinates, basename: "crack_result");
// Produces: crack_result.case, crack_result.geo

// Multi-step time series
for (int step = 0; step < nSteps; step++) {
    double[,] deformed = ComputeDeformedCoords(u[step]);
    double[,] disp     = u[step];
    EnsightWriter.AddMesh($"t={step * dt:F3}", mesh, deformed, displacement: disp);
}
EnsightWriter.WriteAllMeshes("simulation");
// Produces: simulation.case, simulation_0.geo ... simulation_N.geo
//           simulation_0.CrackOpening ... (if displacement provided)
```

Supported element types in Ensight output: `tria3` (Tri3), `quad4` (Quad4), `tetra4` (Tet4).

---

## Examples

The `Teste` project contains 26 progressively complex examples. Run all:

```bash
dotnet run --project Teste -c Release -p:Platform=64
```

### Part 1 — Advanced 2D Meshing (Examples 1–10)
Circular domains, L-shapes, annuli, wedges, multiple holes, gear-like geometries, and triangular-vs-quad comparison. Each example generates both GiD `.msh` and Ensight output.

### Part 2 — 2D Fracture Mechanics (Examples 11–15)
Classical benchmark problems:

| Example | Reference | Description |
|---------|-----------|-------------|
| 11 | Anderson (2005) | Edge crack in tension |
| 12 | Griffith (1921) | Centre crack — original fracture mechanics benchmark |
| 13 | Kanninen & Popelar (1985) | Double-edge notch |
| 14 | Erdogan & Sih (1963) | Slant crack, mixed-mode |
| 15 | Newman & Raju (1984) | Crack emanating from a circular hole |

### Part 3 — Crack Patterns (Examples 16–20)
Artistic/synthetic crack geometries for stress-testing the crack insertion algorithm: spiral galaxy, fractal tree, sinusoidal waves, starburst, concentric mandalas.

### Part 4 — 3D Fracture Mechanics (Examples 21–26)

| Example | Reference | Description |
|---------|-----------|-------------|
| 21 | Sneddon (1946) | Penny-shaped crack |
| 22 | Irwin (1962) | Elliptical crack |
| 23 | Tada (1973) | 3D edge crack |
| 24 | Newman & Raju (1981) | Corner crack |
| 25 | Erdogan & Sih (1963) | 3D slant crack |
| 26 | — | Semi-cylindrical surface crack |

All examples write GiD/CIMNE `.msh` files and a unified Ensight case file. Mesh quality statistics are printed to stdout.

---

## Performance

### Vectorisation

The library uses Roslyn intrinsics throughout hot paths:

| Instruction Set | Usage |
|-----------------|-------|
| AVX2 (`Vector256<double>`) | Sparse matrix-vector products, dense matrix ops |
| AVX-512 (`Vector512<double>`) | Dense matrix multiply (width-8 FMA) |
| FMA | Fused multiply-add in GEMM and integrator |

SIMD is selected at runtime via `Vector512.IsHardwareAccelerated` / `Vector256.IsHardwareAccelerated` with scalar fallback.

### Parallelism

```
ParallelConfig.MaxDegreeOfParallelism  — global CPU thread cap
ParallelConfig.IsGPUAvailable          — check GPU availability at runtime
ParallelConfig.SetMklNumThreads(n)     — limit Intel MKL thread count
```

| Operation | Parallel threshold |
|-----------|--------------------|
| O2M graph traversal | ≥ 1 000 elements |
| CSR row operations | ≥ 1 000 rows |
| CSR SIMD | ≥ 5 000 rows |
| CSR GPU (cuSPARSE) | ≥ 50 000 rows **and** ≥ 1 000 000 nnz |
| FEM assembly lock stripes | 4 096 stripes (golden-ratio hashed) |

### Memory

- `ArrayPool<T>` for large temporary buffers in CSR and Matrix operations
- `GC.AllocateUninitializedArray` to skip zero-fill on known-overwritten arrays
- Pinned arrays for large matrices (≥ several MB) to reduce GC pressure
- `stackalloc` for small per-call buffers inside hot loops
- `ObjectPool<HashSet<int>>` for duplicate detection in sparse assembly

### Scale

- Tested with systems exceeding **10 million DOF**
- `gcAllowVeryLargeObjects` enabled — arrays larger than 2 GB are supported
- `GC.Server = true` in the Release runtime configuration for throughput-oriented collection

---

## Architecture & Design

### Type Safety via ITypeMap

`Topology<TTypes>` is parameterised by a `TTypes : ITypeMap, new()`. The type map translates a C# type argument (`typeof(Triangle)`) to an integer index at JIT compile time, enabling zero-overhead generic dispatch without dictionaries or reflection at runtime:

```csharp
// Resolved at JIT time — no runtime Dictionary lookup
int idx = typeMap.IndexOf<Triangle>(); // inlined to a constant
```

### Thread-Safety Model

| Layer | Strategy |
|-------|----------|
| `O2M` | None — single-threaded |
| `M2M` | `ReaderWriterLockSlim` — concurrent reads, exclusive writes |
| `MM2M` | Per-block `ReaderWriterLockSlim` |
| `Topology<TTypes>` | Inherited from MM2M; `WithBatch()` holds write lock for duration |
| `CSR.AddToElement` | Per-object `lock` for fine-grained element updates |
| `CliqueSystem.Assemble` | Lock-striped (4 096 stripes) |

### Batch Operations

`WithBatch()` uses an RAII pattern: acquire the write lock, execute the action, defer all synchronisation (transpose recompute, position cache invalidation) until the batch exits. Batches can be nested.

### Stale References After Compress

After `M2M.Compress()` or `MM2M.Compress()`, previously stored M2M references are disposed. The `Version` property increments on each structural change. Always use `WithBlock()` rather than caching a reference from `mm2m[type1, type2]`:

```csharp
// UNSAFE — reference becomes stale after Compress()
M2M m = mm2m[0, 1];
mm2m.Compress();
m.Add(...);   // throws ObjectDisposedException

// SAFE
mm2m.WithBlock(0, 1, m => m.Add(...));  // reference valid only inside lambda
```

### JSON Serialization Format

The full topology — entity counts, adjacency lists, per-entity data, symmetry groups, and canonical index maps — is round-trip serialized to JSON via `System.Text.Json`. Custom DTOs (`TopologyDto`, `AdjacencyDto`, `SymmetryDto`, …) are used to control the wire format independently of internal representation changes.

---

## Platforms

| OS | Architecture | Notes |
|----|-------------|-------|
| Windows | x64 | Full support including PARDISO and CUDA |
| Linux | x64 | Full support; MKL via apt, CUDA via standard toolkit |
| macOS | x64, ARM64 | CPU-only (no PARDISO/CUDA on ARM64) |

The project targets `net9.0` with `LangVersion=latest` and requires a 64-bit process (`<PlatformTarget>x64</PlatformTarget>`).

---

## Documentation

Detailed API reference is in `Docs/`:

| File | Contents |
|------|----------|
| `Topology-Complete-Documentation.md` | Topology operations, graph algorithms, serialization, ITypeMap design |
| `Numerical-Complete-Documentation.md` | Dense/sparse matrices, FEM assembly, PARDISO, cuSPARSE integration |
| `SimplexRemesher-Complete-Documentation.md` | Mesh generation, refinement, crack insertion, I/O formats |

---

## License

GPLv3

## Author

Pedro Areias (IST — Instituto Superior Técnico, Lisbon)

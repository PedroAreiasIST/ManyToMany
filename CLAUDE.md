# CLAUDE.md

> Guide for AI assistants working on the ManyToMany (Numerical) codebase.

## Project Overview

ManyToMany is a high-performance scientific computing and finite element analysis (FEA) library for **.NET 9.0**, written in C# by Pedro Areias (IST). It provides a unified framework for mesh topology management, sparse linear algebra, mesh generation with crack insertion, nonlinear time integration, and post-processing for computational mechanics.

**License:** GPLv3

## Repository Structure

```
ManyToMany/
├── Numerical.sln             # Main solution (6 projects)
├── Relations/                # Core: topology & connectivity (namespace: Topology)
│   ├── Topology.cs           # Generic type-safe topology container, graph algorithms
│   ├── O2M.cs                # One-to-many sparse adjacency (single-threaded)
│   ├── M2M.cs                # Many-to-many relationships (thread-safe, cached)
│   ├── MM2M.cs               # Multi-type many-to-many traversal
│   ├── Utils.cs              # TypeMap generics, DFS/BFS/topological sort
│   ├── Symmetry.cs           # Symmetry group canonical representation
│   └── ParallelConfig.cs     # Global parallelism + GPU/MKL configuration
├── Matrices/                 # Dense/sparse linear algebra (namespace: Numerical)
│   ├── Matrix.cs             # Dense matrix ops with SIMD (LU/QR/SVD)
│   ├── CSR.cs                # Compressed Sparse Row, iterative/direct solvers
│   ├── Assembly.cs           # Lock-striped parallel FEM assembly
│   └── NativeLibraries.cs    # CUDA/cuSPARSE/MKL/PARDISO detection & loading
├── Meshing/                  # Mesh generation & refinement (namespace: Numerical)
│   ├── SimplexMesh.cs        # Core mesh topology for simplex elements
│   ├── SimplexRemesher.cs    # Mesh I/O (GiD, Gmsh, ASCII) + refinement wrappers
│   ├── MeshGeneration.cs     # Structured mesh + Delaunay triangulation
│   ├── MeshRefinement.cs     # Longest-edge bisection refinement
│   ├── CrackOperations.cs    # Level-set based crack insertion
│   ├── MeshOperations.cs     # Mesh utilities and optimizations
│   ├── GeometryCore.cs       # Geometry computation helpers
│   ├── FiniteElementTopologies.cs  # Element type definitions (Tri3, Quad4, Tet4)
│   └── MeshConstants.cs      # Tolerances and constants
├── Nonlinear/                # Time integration & root finding (namespace: Numerical)
│   ├── Integrator.cs         # Bathe two-stage implicit integrator
│   └── RootFinder.cs         # Newton-Raphson, IQI, ITP algorithms
├── Postprocess/              # Visualization export (namespace: Numerical)
│   └── EnsightWriter.cs      # Ensight 6.0 format export for GiD
├── Teste/                    # Examples & demos (26 examples)
│   └── Examples2DA.cs        # Advanced meshing + 2D/3D fracture mechanics
└── Docs/                     # Detailed documentation
    ├── Topology-Complete-Documentation.md
    ├── Numerical-Complete-Documentation.md
    └── SimplexRemesher-Complete-Documentation.md
```

## Build Commands

```bash
# Build (Release, 64-bit) -- the standard build command
dotnet build Numerical.sln -c Release -p:Platform=64

# Build (Debug, 64-bit)
dotnet build Numerical.sln -c Debug -p:Platform=64

# Run examples
dotnet run --project Teste -c Release -p:Platform=64
```

**Platform constraint:** x64 only. No 32-bit support. The `-p:Platform=64` flag is required.

## Project Dependencies

```
Teste (executable)
├── Matrices     → Relations
├── Meshing      → Relations
├── Nonlinear    (no internal deps)
├── Postprocess  → Meshing, Relations
└── Relations    (core, no deps)
```

**NuGet packages:**
- `Microsoft.Extensions.ObjectPool` 9.0.0 (used by Relations, Matrices)
- `Intel.oneAPI.MKL.redist.win` 2024.2.1 (Windows only, conditional)

**Optional native libraries** (auto-detected at runtime):
- Intel MKL / PARDISO (sparse direct solver)
- CUDA / cuSPARSE (GPU acceleration)

## Namespaces

Despite the project names, there are only two root namespaces:
- **`Topology`** -- used by the Relations project
- **`Numerical`** -- used by Matrices, Meshing, Nonlinear, and Postprocess

## Testing

There is no unit test framework (no xUnit/NUnit/MSTest). Validation is done through the **Teste** project, which contains 26 example-based demonstrations covering meshing and fracture mechanics. Run with:

```bash
dotnet run --project Teste -c Release -p:Platform=64
```

Examples output GiD `.msh` files and Ensight case files for visual verification.

## Key Conventions

### Naming
- **Classes/Methods:** PascalCase (`Matrix`, `GetNodesForElement`, `WithBatch`)
- **Private fields:** `_camelCase` (`_o2m`, `_rwLock`, `_adjacencies`)
- **Parameters:** camelCase (`xmin`, `numElements`, `func`)
- **Generic type parameters:** `T`, `TTypes`, `TElement`

### Code Patterns
- **Regions:** `#region`/`#endregion` used extensively for method grouping
- **XML docs:** Comprehensive `/// <summary>`, `/// <remarks>`, `/// <param>` on public APIs
- **Static factory methods:** preferred over constructors (e.g., `SubEntityDefinition.FromEdges()`)
- **Builder patterns:** `Topology.WithBatch()`, `Topology.WithSymmetry()`

### Performance Patterns (critical -- do not regress)
- **SIMD intrinsics:** Direct `Vector512<T>`, `Vector256<T>`, `Avx2`, `Sse41` usage in hot paths
- **Zero-allocation:** `Span<T>`, `stackalloc`, `GC.AllocateUninitializedArray<T>` on hot paths
- **Object pooling:** `ObjectPool<T>` for temporary allocations
- **ArrayPool:** `ArrayPool<T>.Shared` for temporary buffer management
- **Parallel.For:** All parallelism routed through `ParallelConfig.Options` with tunable thresholds
- **Cache-aware:** Column-major matrix storage, blocked GEMM operations

### Thread Safety
- **ReaderWriterLockSlim:** Used in `M2M` and `MM2M` for concurrent read access
- **Lock ordering:** `MM2M._rwLock` (outer) then `M2M._rwLock` (inner) -- never reverse
- **Double-checked locking:** For lazy initialization (e.g., MKL loading in `ParallelConfig`)
- **Volatile fields:** For thread-safe caching (`_cachedOptions`, `_isGPUAvailable`)
- **Immutable types:** Records and readonly structs for value types

### Error Handling
- Precondition validation at method entry (`ArgumentNullException`, `ArgumentOutOfRangeException`)
- `IDisposable` with `ThrowIfDisposed()` guard pattern
- Graceful fallback when native libraries are unavailable

## Architecture Notes

### Relations (Core)
The foundational layer. `Topology<TTypes>` provides compile-time type-safe mesh entity management using `TypeMap<>` for type-to-index mapping. `O2M` is the low-level single-threaded sparse adjacency; `M2M` wraps it with thread safety and caching; `MM2M` provides multi-type traversal.

### Matrices
Dense (`Matrix`) and sparse (`CSR`) linear algebra with hardware acceleration. `Assembly` implements lock-striped parallel FEM matrix assembly with a clique-based element connectivity system. `NativeLibraries` handles runtime detection and loading of MKL/CUDA via `NativeLibrary.TryLoad()` with platform-specific library names.

### Meshing
Simplex mesh generation and manipulation. `MeshGeneration` creates structured meshes and Delaunay triangulations. `MeshRefinement` implements longest-edge bisection. `CrackOperations` handles level-set crack insertion for fracture mechanics. `SimplexRemesher` provides multi-format I/O (GiD, Gmsh, VTK, ASCII).

### Nonlinear
`BatheTwoStageIntegrator` implements an unconditionally stable 2nd-order implicit time integrator for structural dynamics. `RootFinder` provides Newton-Raphson and bracketing methods (IQI, ITP) with compensated (Kahan) summation.

## Things to Watch Out For

1. **Always build with `-p:Platform=64`** -- the solution uses a custom platform name `64`, not `x64` or `AnyCPU`.
2. **Nullable warnings suppressed:** `CS8600`, `CS8602`, `CS8603`, `CS8604` are suppressed in Relations via `<NoWarn>`. Be aware of potential null reference issues.
3. **Unsafe code:** `AllowUnsafeBlocks` is enabled in Relations, Matrices, Nonlinear, and Teste for P/Invoke and SIMD intrinsics.
4. **No CI pipeline:** There are no GitHub Actions or other CI workflows. Build verification is manual.
5. **Large objects:** The runtime is configured to allow >2GB arrays (`AllowVeryLargeObjects`). Keep this in mind when working with memory-sensitive code.
6. **Globalization invariant:** Matrices Release builds use `InvariantGlobalization`, so locale-dependent string operations will not behave as expected.

## Commit Style

Based on repository history, commits use short imperative descriptions:
- `Fix thread-safety and consistency issues across O2M and Topology`
- `Fix off-by-one in MM2M.GetNodeRangeForTypeUnlocked skipping last node`
- `Improve README with comprehensive project documentation`

Prefix with the action: `Fix`, `Add`, `Update`, `Improve`, `Remove`.

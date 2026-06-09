# ManyToMany — Documentation

This directory holds the in-depth reference documentation for **ManyToMany**. For a high-level tour, the [project README](../README.md) is the best starting point — it includes a quick start, installation instructions, and a verified API reference for every module.

> **Namespace:** every public type in the library lives in the single `Numerical` namespace. All example code in these documents assumes `using Numerical;`.

## Reading guide

| If you want to… | Read |
|---|---|
| Get running in five minutes | [README → Quick Start](../README.md#quick-start) |
| Understand the module layout | [README → Architecture Overview](../README.md#architecture-overview) |
| Work with dense/sparse matrices & FE assembly | [Numerical-Complete-Documentation.md](Numerical-Complete-Documentation.md) |
| Model connectivity & run graph algorithms | [Topology-Complete-Documentation.md](Topology-Complete-Documentation.md) |
| Generate meshes, refine, and insert cracks | [SimplexRemesher-Complete-Documentation.md](SimplexRemesher-Complete-Documentation.md) |
| See the theory & motivation | [p_areias_simple_csharp_final.pdf](p_areias_simple_csharp_final.pdf) |

## Documents

### [Numerical-Complete-Documentation.md](Numerical-Complete-Documentation.md)
Dense matrices (`Matrix`, `Vector`) and their decompositions (LU with rook pivoting, QR, SVD, symmetric eigen), the sparse `CSR` type with its automatic backend selection (PARDISO / cuSPARSE / SIMD / parallel), the `CliqueSystem` finite-element assembler, and cross-platform native-library integration.

### [Topology-Complete-Documentation.md](Topology-Complete-Documentation.md)
The `Topology<TTypes>` container: entity management, typed attributes, adjacency queries, sub-entity discovery, the full suite of graph algorithms (BFS/DFS, coloring, topological order, Cuthill–McKee, components), symmetry/canonicalization, and JSON serialization.

### [SimplexRemesher-Complete-Documentation.md](SimplexRemesher-Complete-Documentation.md)
Structured and unstructured mesh generation, longest-edge conforming refinement with lineage tracking, level-set crack insertion in 2D and 3D, mesh smoothing and quality analysis, and file I/O.

## A note on accuracy

The [project README](../README.md) is the **canonical, code-verified API reference** — its method names and signatures have been checked against the source. The deep-dive documents in this folder are broader reference and tutorial material; where you spot any disagreement, trust the README and the source code. Found an error? Please [open an issue](https://github.com/PedroAreiasIST/ManyToMany/issues).

## Supported file formats

ManyToMany reads and writes the following mesh formats (see `SimplexRemesher` and `EnsightWriter`):

| Format | Read | Write |
|---|:---:|:---:|
| Gmsh `.msh` (v2) | ✅ | ✅ |
| GiD / CIMNE `.msh` | ✅ | ✅ |
| Plain ASCII | — | ✅ |
| Ensight Gold `.case` (ParaView) | — | ✅ |

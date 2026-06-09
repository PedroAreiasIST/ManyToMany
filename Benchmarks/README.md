# Benchmarks

Reproducible connectivity-query benchmark comparing mesh data-structure
libraries (OpenMesh, CGAL, VTK, libigl, and MM3) on identical meshes, emitted in
the comparison-table format used in the paper.

## Files

| File | Purpose |
| --- | --- |
| `mesh_connectivity_benchmark.py` | The cross-library harness. No required dependencies to run (NumPy only for synthetic meshes). Each benchmarked library is optional. |
| `Mm3Bench/` | A tiny .NET console project the harness shells out to so it can time MM3 alongside the other libraries (references `Relations` + `Meshing`). |
| `Mm3Bdn/` | A standalone **BenchmarkDotNet** project to independently reproduce and *assess* the MM3 numbers with full statistics (mean/error/stddev/median + allocations). Needs no Python or external library — see [below](#benchmarkdotnet-harness-mm3bdn). |

## Run

```bash
pip install numpy psutil libigl vtk openmesh    # each optional; CGAL via conda if wanted
dotnet build Mm3Bench -c Release
python mesh_connectivity_benchmark.py --mm3-project ./Mm3Bench --out results
# fast smoke test:
python mesh_connectivity_benchmark.py --quick --mm3-project ./Mm3Bench
# real Stanford Bunny/Dragon instead of synthetic meshes:
python mesh_connectivity_benchmark.py --download --mm3-project ./Mm3Bench
```

Outputs a table to stdout and, with `--out results`, `results.csv` + `results.md`
+ `results.tex` (a plain `\begin{tabular}` — no extra packages — that you can
`\input{}` straight into the paper).

## Installing the comparison libraries on Windows (verified 2026-06)

The `pip install …` line above is the Linux-friendly path. On Windows the four
comparison libraries split into two environments because **OpenMesh has no wheel
for CPython ≥ 3.10 and its binaries are built against the NumPy 1.x ABI**:

- **Modern stack — Python 3.11+, NumPy 2.x** (`pip install vtk cgal libigl`):
  gets **VTK**, **libigl**, and **CGAL**. Use the PyPI package **`cgal`** (the
  maintained CGAL SWIG bindings, ships a `cp311` Windows wheel) — *not* the
  abandoned `cgal-bindings`. OpenMesh is **not installable** here (no `cp311`
  wheel anywhere; a source build fails on MSVC with `error C2065: 'ssize_t'`).

- **Full stack incl. OpenMesh — Python 3.9 + NumPy 1.x** (single environment for
  all five, used for the table above):

  ```powershell
  # 3.9 interpreter (any per-user install works); then an isolated venv:
  py -3.9 -m venv .venv39 ; .\.venv39\Scripts\Activate.ps1
  python -m pip install -U pip
  python -m pip install "numpy<2"                       # OpenMesh needs the 1.x ABI
  python -m pip install --only-binary=:all: openmesh vtk cgal trimesh psutil "libigl==2.6.1"
  ```

  `--only-binary=:all:` makes pip fail fast instead of dropping into a long MSVC
  build; pin `libigl==2.6.1` because the newer 2.6.2 has no Windows wheel. Note
  the `import` names differ from the install names: `import igl` (libigl) and
  `import CGAL` (cgal). `libigl` does not expose `tet_tet_adjacency`, so it is
  `N/A` on the tet mesh regardless of platform.

Verify the environment before benchmarking:

```powershell
python -c "import openmesh, vtk, igl, numpy; from CGAL.CGAL_Polyhedron_3 import Polyhedron_3; print('all five OK', numpy.__version__)"
```

## Workload

For every library: build the topology/adjacency, then traverse every cell's
face-neighbours and every vertex's incident cells once. The reported number is
total **build + query** milliseconds (median of `--repeats` runs; build and query
are also reported separately in the CSV).

## Honest scope

- **Host-specific.** Numbers depend on CPU, memory bandwidth, thread count and
  library versions; they differ between machines and are not comparable across
  hosts. They are not meant to reproduce a specific published table — the paper's
  measurements were taken on an Intel Core i9-13900HX.
- **Mixed models & bindings.** OpenMesh/CGAL/VTK answer queries by per-element
  traversal of an explicit topology, driven here through their **Python** bindings,
  so their query time includes a Python call per element. libigl and MM3 build full
  adjacency arrays once and answer in compiled code. Read across a row as
  order-of-magnitude; the fairest single comparison is **MM3 vs libigl** — both are
  compiled libraries that build full adjacency arrays.
- **Surface vs volume.** Half-edge surface libraries (OpenMesh, CGAL) do not
  represent tetrahedra, so they are `N/A` on the tet mesh; this libigl build also
  lacks `tet_tet_adjacency`. MM3 and VTK handle both.
- **Steady state, not cold start.** MM3 is the only JIT-compiled entrant. The
  bridge therefore gives it one *discarded* warm-up and reports the median of
  `--repeats` iterations run in a single warm process — the same fully-optimized
  regime the AOT-compiled native libraries and the libigl/CGAL wheels are always
  in. It also builds the topology with the bulk `AddRangeParallel` API (how a large
  mesh is built in practice), not one locked `Add` per element. So the MM3 figure
  is steady-state throughput; earlier revisions that re-spawned the process per
  sample were measuring .NET start-up + cold JIT, which understated it.

## MM3 bridge (`Mm3Bench/Program.cs`)

Implements the contract the harness shells out to:

```
dotnet run -c Release --project Mm3Bench -- <meshfile> <tri|tet> <repeats>
# prints final line:  BUILD_MS QUERY_MS
```

It reads the mesh once (file I/O is never timed), builds a
`Topology<TypeMap<Node, Tri3>>` (or `Tet4`) with `AddRangeParallel`, then times
`GetElementToElementGraph<…>()` + `GetTranspose<…>()` plus the neighbour
traversal. It runs one discarded warm-up and `<repeats>` timed iterations in this
one warm process and prints the **median** build/query — so re-running it is not
needed; pass the harness's `--repeats` straight through. If you change those
public signatures, update the calls in `Program.cs`.

## BenchmarkDotNet harness (`Mm3Bdn`)

`Mm3Bench` exists to slot MM3 into the cross-library table. To **independently
assess** MM3's numbers with proper statistics, use the standalone BenchmarkDotNet
project `Mm3Bdn/`. It generates the same synthetic meshes in C# (no Python, no
external library) and benchmarks both *build + query* and *query only* (adjacency
on a pre-built topology) for triangle and tetrahedral meshes, reporting
mean/error/stddev/median and allocations.

```bash
dotnet run -c Release --project Benchmarks/Mm3Bdn                      # all
dotnet run -c Release --project Benchmarks/Mm3Bdn -- --filter '*Tet*'  # just tets
dotnet run -c Release --project Benchmarks/Mm3Bdn -- --list flat       # list cases
```

It is intentionally **not** part of `Numerical.sln` (so a solution build stays
lean) and must be run in `Release` — BenchmarkDotNet refuses a debug build.

## CI

`.github/workflows/connectivity-benchmark.yml` runs this on every push: it builds
the bridge, installs the Python libraries, runs the harness (`--quick`), writes
the table to the job summary, and uploads `results.md`/`results.csv`/`results.tex`
as artifacts.

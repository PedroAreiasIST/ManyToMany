# Benchmarks

Reproducible connectivity-query benchmark comparing mesh data-structure
libraries (OpenMesh, CGAL, VTK, libigl, and MM3) on identical meshes, emitted in
the comparison-table format used in the paper.

## Files

| File | Purpose |
| --- | --- |
| `mesh_connectivity_benchmark.py` | The harness. No required dependencies to run (NumPy only for synthetic meshes). Each benchmarked library is optional. |
| `Mm3Bench/` | A tiny .NET console project that lets the harness time MM3 (references `Relations` + `Meshing`). |

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

Outputs a table to stdout and, with `--out results`, `results.csv` + `results.md`.

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
- **Mixed models.** OpenMesh/CGAL/VTK answer queries by per-element traversal of
  an explicit topology; libigl and MM3 build full adjacency arrays once. Read
  across a row as order-of-magnitude.
- **Surface vs volume.** Half-edge surface libraries (OpenMesh, CGAL) do not
  represent tetrahedra, so they are `N/A` on the tet mesh.
- **JIT.** MM3's first invocation pays JIT warm-up; pass `--repeats >= 3` for a
  steadier figure.

## MM3 bridge (`Mm3Bench/Program.cs`)

Implements the contract the harness shells out to:

```
dotnet run -c Release --project Mm3Bench -- <meshfile> <tri|tet> <repeats>
# prints final line:  BUILD_MS QUERY_MS
```

It builds a `Topology<TypeMap<Node, Tri3>>` (or `Tet4`), then times
`GetElementToElementGraph<…>()` + `GetTranspose<…>()`. If you change those public
signatures, update the two calls in `Program.cs`.

## CI

`.github/workflows/connectivity-benchmark.yml` runs this on every push: it builds
the bridge, installs the Python libraries, runs the harness (`--quick`), writes
the table to the job summary, and uploads `results.md`/`results.csv` as artifacts.

# Step2: incremental, persistent and concurrent connectivity for mesh cutting and crack growth

Code of the paper *Incremental, persistent and concurrent finite element connectivity for mesh cutting and crack growth*
(P. Areias and N. Silvestre). It extends the MM3 connectivity kernel of this repository with counting relations,
chunk-shared persistent versions and a local Steiner-point-free cutter.

| Project | Content | Paper |
|---|---|---|
| `Step2.Core` | the library (3 998 lines of C#) | Sections 2 to 6 |
| `Step2.Bench` | verification, measurement and example driver (8 088 lines) | Sections 7 to 9 |

## Library (`Step2.Core`)

| File | Content | Section |
|---|---|---|
| `PArray.cs` | two-level persistent arrays, epoch tags, commit, rollback | 4 |
| `IncrementalMesh.cs` | rows, transpose with positions, gapped pattern with stable slots, destinations, facet adjacency, boundary set, coloring; Algorithm 1 | 2, 3 |
| `Snapshot.cs` | read-only handle on a committed version | 4 |
| `Attributes.cs` | versioned per-element and per-node attributes | 4 |
| `FlatBuild.cs`, `FlatPattern.cs` | tuned parallel rebuild; flat copy of the pattern for solvers | 6 |
| `Cutter.cs`, `ParallelCutter.cs` (templates: `Meshing/SubdivisionTemplates.cs`, shared with `MeshRefinement`) | templates of the batch cutter; local cut, revert and concurrent cut | 5 |
| `CrackInserter.cs` | crack increments from two level sets | 5.3 |
| `ElementSpec.cs` | element types | 2 |

## Build

.NET SDK 9.0.200 (`global.json`), x64. `Step2.Bench` references the `Relations` and `Meshing` projects of this repository
(MM3 kernel, batch cutter and the structured mesh generator) by relative paths; the library is `src/Incremental` (assembly `Step2.Core`, namespace `Step2`) and the driver `benchmarks/Step2Bench`.

```bash
dotnet build Step2.Bench -c Release
```

## Verification (Section 7)

```bash
./run-maxperf.sh verify-fuzz 3 3000 5 12345     # oracle after every batch, immutability, rollback, destination stability
./run-maxperf.sh verify-work 6                  # crack, bisection and erosion workloads with the oracle after every step
./run-maxperf.sh verify-cut 4 9 1               # local cut against the batch cutter
./run-maxperf.sh verify-uncut 4 10 1            # local revert
./run-maxperf.sh verify-crack 8 6 0             # level-set crack insertion
./run-maxperf.sh verify-cutpar 6 8 1 8          # concurrent cut against the sequential cut
./run-maxperf.sh verify-rowlimit 6 20           # a cut that would overflow a pattern row is rejected
./run-maxperf.sh verify-flatbuild 30            # tuned parallel rebuild against the sequential one
```

The fuzz gates of the paper are `verify-fuzz` with `3 3000 5 12345`, `2 5000 3 1`, `4 2000 20 7`, `5 1000 50 99` and `3 20000 2 2024`
(160 000 edits, 837 rollbacks, 31 000 immutability checks, 6 467 348 destination blocks of unedited elements compared).

## Measurements (Section 8)

`run-maxperf.sh` pins the process to the physical cores (`STEP2_CPUS`, default `0-15`), disables tiered compilation and ReadyToRun,
and runs the server collector with concurrent collection off and a 512 MB generation-0 budget (`STEP2_GEN0`, hexadecimal bytes).
Results are written to the directory `STEP2_DATA`. Every timed region records the collections and the pause time of the collector;
the cutting tables use the time without pauses. The final run is `benchmarks/run-fin4-timed-amd2.sh`
(one measurement at a time on an idle machine), the node runs are `benchmarks/run-cutpar7-deu.sbatch` and `benchmarks/run-pennyctl-deu.sbatch`.

| Mode | Table or figure |
|---|---|
| `mm3ref` | MM3 rebuild protocol (Section 8.1) |
| `work crack|bisect|erode n` | workloads, scaling (Section 8.2) |
| `crossover`, `interleaved` | batch size, interleaved queries and edits (8.3, 8.4) |
| `reads`, `concurrent`, `persist`, `memory` | reads, versioning, memory (8.5, 8.6) |
| `endtoend`, `profile` | complete step (8.7) |
| `batchcut`, `cutpar`, `cutgc`, `interfaces` | local and concurrent cutting (8.8) |
| `penny`, `growth`, `surface`, `corner` | Section 9 |

## License

GPL-3.0, as the repository.

## Scripts that regenerate the tables and figures

`benchmarks/` (copy next to the project): `run-*.sh` and `run-*.sbatch` are the drivers of the runs behind Sections 7 to 9
(`run-fin4-timed-amd2.sh` is the single timed run of Section 8; `run-fin5` to `run-fin8` are the follow-up controls;
`run-cutpar7-deu.sbatch` and `run-pennyctl-deu.sbatch` the Deucalion runs). `table_work.py` and `compare_logs.py` reduce the logs to
table rows; `tools/` holds the helpers that write the numbers into the LyX tables. `figs/` holds `make_data.py`, `make_examples.py`,
the gnuplot files, `rerender.sh`, the Figure 1 source `structures.tex` and `corner/` (the EnSight EPS processing of the corner-crack figure).
The data (CSV) and logs of the runs are not included; the corner-crack EnSight Gold files (thirteen steps, about 5 GB) are available on request.

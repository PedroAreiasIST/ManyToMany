# Mm3Remeasure: harness of the revised ADES Tables 10-14

Re-measurement harness for the MM3 paper (P. Areias, M. Sardinha, N. Silvestre,
*Advances in Engineering Software* 222 (2026) 104296), revision 2. It regenerates
the MM3 numbers of these tables of the revised manuscript:

| Command | File | Paper table |
| --- | --- | --- |
| `table10` | `Table10.cs` | **Table 10**: construction time, throughput and memory (10K-10M quads; live bytes and the flat index count) |
| `table11` | `Table11.cs` | **Table 11**: query throughput, every row and both columns; "Our Method" on the span (zero-copy, lock hoisted) and the public API access path, against the caption-exact **Naive** baseline (primary element-to-node relation only: no transpose, no position caches) |
| `tables1213` | `Tables1213.cs` | **Table 12** (cache synchronization after one modification, break-even query count), **Table 13** (element insertion with and without a symmetry group) and the fast-path share quoted in **Table 4** |
| `table14` | `Table14.cs` | **Table 14**: the MM3 rows of the framework comparison on the 100K quad mesh (build, `NodesOf`, neighbours; cached and uncached). The Boost/NetworkX rows are third-party figures and are not re-measured |
| `all` | | all of the above, in order |

`RESULTS.md` is the record of the measurement that went into the revision
(2026-07-22, i9-13900HX, .NET 9.0.2, taken against `src/Relations` at `6eeac3c`),
including the comparison against the originally published values.

## Run

```bash
dotnet build benchmarks/Mm3Remeasure -c Release
benchmarks/Mm3Remeasure/run-maxperf.sh table14     # pinned, ASLR off, tiering off
dotnet run -c Release --project benchmarks/Mm3Remeasure -- table11   # unpinned
```

`table14` is the default and the quickest (about ten seconds). `table10` goes up to
10M elements and `table11` runs the quadratic Naive connected-components pass on
growing meshes: both need several GB of memory and minutes to hours.

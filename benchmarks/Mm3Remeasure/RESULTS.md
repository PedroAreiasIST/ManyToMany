# Table 11 re-measurement — every row, both columns

Date 2026-07-22. Intel Core i9-13900HX, .NET 9.0.2, Linux, Release.
`taskset -c 0-15` (P-cores) + `setarch -R`; tiered compilation, TieredPGO and
ReadyToRun off; server GC, concurrent GC off; full blocking collection before every
sample; minimum of 5–7 samples; every sample printed. Mesh 1,000,000 quadrilateral
elements, 1,002,001 nodes, position caches built before timing.
Harness `scratchpad/Mm3Remeasure` against `src/Relations/Relations.csproj` at `6eeac3c`.
Logs: `table11_run1.log`, `table11_run2.log`.

`Naive` holds only the primary element→node relation — no transpose, no position
caches — as the published caption defines it. Every reverse or transitive query
therefore recomputes from the primary relation.

## Measured (run 1 / run 2)

| Operation | Our Method | Naive | Speedup |
|---|---:|---:|---:|
| `NodesOf<Element>` span | 52.5M / 50.5M | 59.9M / 57.2M | **1×** |
| `NodesOf<Element>` copying API | 11.0M / 10.7M | 12.1M / 8.1M | 1× |
| `ElementsAt<Node>` span, lock hoisted | 52.6M / 51.7M | 514 / 517 | **~1.0×10⁵** |
| `ElementsAt<Node>` copying API | 5.8M / 6.1M | 514 / 517 | 1.2×10⁴ |
| `GetDirectNeighbors` | 1.23M / 1.20M | 211 / 226 | **~5.3×10³** |
| `GetKHopNeighborhood(k=2)` | 215K / 200K | 23 / 25 | **~8.0×10³** |
| `FindConnectedComponents` | 3.8 / 3.3 | 2.3×10⁻⁴ (extrap.) | **~1.4×10⁴** |

Naive `FindConnectedComponents` is n applications of the uncached neighbour query and
is quadratic in n. Measured at n = 6,400 (176 ms), 25,600 (2,817 ms) and 102,400
(45,734 ms); the n² model holds to 1.00 and 1.01, and all three extrapolate to
4,298–4,361 s per pass at n = 10⁶, i.e. 2.3×10⁻⁴ ops/sec.

## Against the published table

| Operation | Published our / naive / speedup | Measured our / naive / speedup |
|---|---|---|
| `NodesOf<Element>` | 50M / 500K / 100× | 50M / 57M / **1×** |
| `ElementsAt<Node>` | 50M / 80K / 625× | 52M / 517 / **1.0×10⁵** |
| `GetDirectNeighbors` | 400K / 8K / 50× | 1.2M / 226 / **5.3×10³** |
| `GetKHopNeighborhood(k=2)` | 20K / 400 / 50× | 200K / 25 / **8.0×10³** |
| `FindConnectedComponents` | 800 / 20 / 40× | 3.3 / 2.3×10⁻⁴ / **1.4×10⁴** |

1. **`ElementsAt<Node>` at 50M is confirmed** (51.7–52.6M) on the zero-copy transpose
   row read with the read lock taken once for the batch — the same shape as the
   `NodesOf<Element>` row, also confirmed (50.5–52.5M). 1/52M = 19 ns, agreeing with
   the 20 ns of Table 14 and the 20–50 ns of Table 4.
2. **The per-call public API is 8× slower** on the reverse path: `GetElementsForNode`
   gives 5.8–6.1M because `EnterSynchronizedReadLock`/`ExitReadLock` costs ~125 ns
   against a 19 ns row read. Delegate allocation is not a factor — hoisting the lambda
   out of the probe loop changes nothing.
3. **`FindConnectedComponents` at 800 ops/sec is unattainable**: measured 3.3–3.8 ops/sec
   (264–303 ms per pass). 800 ops/sec needs 1.25 ms per pass, i.e. 0.8G entity-visits/s.
   `Topology.cs:4186` visits all 10⁶ entities, calling `GetElementNeighbors` per entity
   with a `List<int>` allocation each, over a 10⁶-entry `HashSet<int>`. The structured
   grid is a single component of 1,000,000 elements, so one BFS covers the whole mesh.
4. **`GetDirectNeighbors` and `GetKHopNeighborhood` are 3× and 10× better than
   published** (1.2M vs 400K; 200K vs 20K).
5. **The whole Naive column is too fast for its stated definition.** A cache-free
   reverse lookup on 10⁶ elements must scan ~4×10⁶ entries; 80,000 ops/sec would
   require 3×10¹¹ comparisons/sec. Measured 517.
6. **`NodesOf<Element>` gains nothing from the caches, correctly** — a forward lookup
   reads the stored primary row in both configurations. The published 100× is not
   reproducible under any reading of the caption.

## Consequence for the claims

The cache-dependent speedup range becomes **5.3×10³–1.0×10⁵×**, replacing 40–625×.
The forward-lookup row should be stated as 1×, which is the honest result and shows
the caches cost nothing on the path that does not use them.

The ratio is large because the baseline recomputes from the primary relation on every
query, so it grows linearly with mesh size. It measures the cost of not caching at
this mesh size; it is not a comparison against another library, and the paper must say
so where the number appears.

Query throughput now spans 50M down to 3.3 ops/sec — seven orders of magnitude, not
five (`mm3_revisednew.tex:2519`).

Table 14 (`tab:framework-comparison`, 100K mesh) carries the same over-fast uncached
row: 125 µs for neighbours, i.e. 8K ops/sec. Scaling the measured 226 ops/sec at 10⁶
to 10⁵ elements gives ~2.3K ops/sec (~440 µs). That row needs the same correction.


## Table 10 — construction, re-measured (same conditions)

| Mesh | Time (ms) | Elements/sec | Memory (MB) | With caches (MB) |
|---|---:|---:|---:|---:|
| 10K | 0.2 | 61M | 0.8 | 3.6 |
| 100K | 1.8 | 56M | 8.0 | 34 |
| 1M | 20 | 49M | 80 | 337 |
| 10M | 360 | 28M | 800 | 3,469 |

Published was 8.3/8.3/8.0/7.8M el/s and 0.6/6.4/64/640 MB. Throughput is 3.5–7×
higher and is **not** flat: it falls 2.2× from 10^4 to 10^7 elements, total time
growing as n^1.09, because one container is allocated per row. The "scales linearly,
within 6% of peak" claim was corrected.

The published memory column was the **index count** (4 B/index over primary +
transpose + both position caches = 64 B/element): 0.6/6.1/61/610 MiB reproduces
0.6/6.4/64/640 MB exactly. Measured resident footprint of the default
`List<List<int>>` is 80 B/element for the primary and ~340 B/element with caches, so
the structure counted as 64 MB occupies 337 MB. Both accountings are now stated, and
Figure 15's caption says which one it plots.

## Figure lightening — the two families need different levels

Tone was compared at page scale (100 dpi rasterisation of the compiled page, which is
how a reviewer reads it on screen; low-dpi rendering snaps hairlines to a full pixel):

| | original | shipped |
|---|---:|---:|
| Fig 18, 3D meshes | 169 | 219 |
| Fig 17, 2D meshes | 155 | 200 |

The 3D figures carry real gray *fills* (0.47–0.70) over large areas, so stroke 0.62 +
fill factor 0.45 lightens them on paper and on screen alike. The 2D figures are
near-white fills (0.96) plus 31.5K black hairlines at 0.048 pt, so their tone is
almost entirely stroke coverage: at 600 dpi they were already light (area-averaged
230) and stroke 0.62 washed them out (245, crack barely visible). Stroke **0.45** was
chosen instead — it removes the on-screen muddiness (155 → 200, matching Fig 18) while
keeping the crack legible in print. The dark 0.16 fills that mark the crack are
preserved by the `fill_min` threshold.

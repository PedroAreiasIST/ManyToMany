using Microsoft.Extensions.ObjectPool;
using System.Buffers;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.X86;
using System.Runtime.Intrinsics;
using System.Text.Json.Serialization;
using System.Text;
using static System.Runtime.CompilerServices.MethodImplOptions;

namespace Numerical;

/// <summary>
///     Represents a one-to-many relationship structure using sparse adjacency lists.
///     High-performance implementation with parallel processing support.
/// </summary>
/// <remarks>
///     <para>
///         <b>THREAD SAFETY:</b> This class is NOT thread-safe. For thread-safe access,
///         use <see cref="M2M" /> which wraps O2M with ReaderWriterLockSlim synchronization.
///     </para>
///     <para>
///         <b>P1-6 CLARIFICATION:</b> The internal <c>_maxNodeIndexCache</c> field is not
///         protected by any synchronization. Callers must ensure external synchronization
///         if accessing O2M from multiple threads. Single-threaded usage is the expected
///         pattern for direct O2M access.
///     </para>
/// </remarks>
[SkipLocalsInit]
public sealed class O2M : IComparable<O2M>, IEquatable<O2M>, ICloneable
{
    #region Cloning

    /// <summary>
    ///     Creates a deep copy of this O2M instance.
    ///     Uses parallel cloning with pre-allocated array.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public object Clone()
    {
        var count = _adjacencies.Count;
        var clonedO2m = new O2M(count)
        {
            ParallelizationThreshold = ParallelizationThreshold
        };

        if (count == 0)
            return clonedO2m;

        if (count >= ParallelizationThreshold)
        {
            var clonedRows = GC.AllocateUninitializedArray<List<int>>(count);

            Parallel.For(0, count, ParallelConfig.Options, i =>
            {
                var src = _adjacencies[i];
                var dst = new List<int>(src.Count);
                if (src.Count > 0)
                {
                    CollectionsMarshal.SetCount(dst, src.Count);
                    CollectionsMarshal.AsSpan(src).CopyTo(CollectionsMarshal.AsSpan(dst));
                }

                clonedRows[i] = dst;
            });

            clonedO2m._adjacencies.AddRange(clonedRows);
        }
        else
        {
            foreach (var row in _adjacencies)
            {
                var newRow = new List<int>(row.Count);
                if (row.Count > 0)
                {
                    CollectionsMarshal.SetCount(newRow, row.Count);
                    CollectionsMarshal.AsSpan(row).CopyTo(CollectionsMarshal.AsSpan(newRow));
                }

                clonedO2m._adjacencies.Add(newRow);
            }
        }

        // The source cache is valid for the cloned data — a deep copy has identical content,
        // so the max node index is the same. If the source cache is null (not yet computed),
        // the clone will also compute on demand.
        clonedO2m._maxNodeIndexCache = _maxNodeIndexCache;

        return clonedO2m;
    }

    #endregion

    #region Factory Methods

    /// <summary>
    ///     Generates a random O2M for testing and benchmarking.
    ///     Parallel generation for large structures.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static O2M GetRandomO2M(
        int elementCount,
        int nodeCount,
        double density,
        int? seed = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(elementCount);
        ArgumentOutOfRangeException.ThrowIfNegative(nodeCount);
        ArgumentOutOfRangeException.ThrowIfNegative(density);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(density, 1.0);

        var o = new O2M(elementCount);
        var expected = (int)Math.Round(nodeCount * density);

        if (elementCount >= DefaultParallelizationThreshold)
        {
            var rows = new List<int>[elementCount];
            var baseSeed = seed ?? Environment.TickCount;

            Parallel.For(0, elementCount, ParallelConfig.Options, i =>
            {
                var rnd = new Random(baseSeed + i);
                var row = new List<int>(Math.Max(expected, 4));

                for (var j = 0; j < nodeCount; j++)
                    if (rnd.NextDouble() < density)
                        row.Add(j);

                rows[i] = row;
            });

            o._adjacencies.AddRange(rows);
        }
        else
        {
            var rnd = seed.HasValue ? new Random(seed.Value) : Random.Shared;

            for (var i = 0; i < elementCount; i++)
            {
                var row = new List<int>(Math.Max(expected, 4));

                for (var j = 0; j < nodeCount; j++)
                    if (rnd.NextDouble() < density)
                        row.Add(j);

                o._adjacencies.Add(row);
            }
        }

        return o;
    }

    #endregion

    #region HashSet Pool Policy

    private sealed class HashSetPoolPolicy : IPooledObjectPolicy<HashSet<int>>
    {
        private const int MaxCapacity = 4096;
        private const int InitialCapacity = 256;

        /// <summary>
        ///     Maximum internal bucket array capacity we tolerate before rejecting the
        ///     set from the pool. Prevents a single large operation from permanently
        ///     bloating pooled instances.
        /// </summary>
        private const int MaxRetainedBucketCapacity = InitialCapacity * 4;

        public HashSet<int> Create()
        {
            return new HashSet<int>(InitialCapacity);
        }

        public bool Return(HashSet<int> obj)
        {
            // Reject sets that grew too large during use — they would waste memory
            // sitting in the pool. Let GC reclaim them; fresh ones will be created.
            if (obj.Count > MaxCapacity) return false;

            obj.Clear();

            // EnsureCapacity(0) with a cleared set returns the current internal bucket
            // array capacity without triggering any resize. Reject if the internal
            // storage grew beyond our retention threshold.
            var currentBucketCapacity = obj.EnsureCapacity(0);
            return currentBucketCapacity <= MaxRetainedBucketCapacity;
        }
    }

    #endregion

    #region Transpose Operation

    /// <summary>
    ///     Transposes the structure: swaps element and node roles.
    ///     Fully parallel implementation with guaranteed thread-safe access.
    /// </summary>
    /// <remarks>
    ///     <b>ORDERING GUARANTEE:</b> The transposed structure maintains sorted adjacency lists.
    ///     Each node's element list is sorted in ascending order.
    ///     This property is critical for correctness and is preserved in both parallel
    ///     and sequential execution paths.
    ///     <b>THREAD SAFETY:</b> Uses arrays for concurrent writes instead of List spans.
    ///     Atomic offset tracking ensures unique positions without race conditions.
    ///     <b>PERFORMANCE:</b> Three-phase algorithm:
    ///     1. Count (parallel with thread-local aggregation)
    ///     2. Allocate (sequential, pre-sized arrays)
    ///     3. Fill (parallel with atomic offsets) + Sort (parallel)
    ///     <b>NEGATIVE NODE HANDLING:</b> Nodes with negative indices (e.g., sentinel
    ///     values like -1) are silently skipped during transpose.
    ///     <b>MEMORY WARNING (P0-2):</b> This method allocates arrays sized by (maxNode + 1),
    ///     where maxNode is the largest non-negative node index in the structure.
    ///     A single element with node index 100,000,000 will cause a 400MB+ allocation.
    ///     For sparse structures with large node indices, use <see cref="Transpose(int)" />
    ///     with an explicit cap.
    ///     If you need strict validation of node indices, use <see cref="TransposeStrict" />.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    public O2M Transpose()
    {
        // Behavior controlled by TransposeSkipsInvalidNodes property
        return TransposeCore(!TransposeSkipsInvalidNodes);
    }

    /// <summary>
    ///     Transposes the structure with an explicit maximum node limit.
    /// </summary>
    /// <param name="maxNodeCap">
    ///     Maximum node index to include (inclusive).
    ///     Nodes > maxNodeCap are silently skipped.
    /// </param>
    /// <returns>Transposed O2M where nodes become elements.</returns>
    /// <remarks>
    ///     <b>USE CASE (P0-2):</b> When you have sparse large node indices and want to avoid
    ///     excessive memory allocation.
    ///     <b>EXAMPLE:</b> If your data has node indices [0, 5, 1000000] but you only
    ///     care about nodes ≤ 100, call <c>Transpose(100)</c> to allocate only 101 entries
    ///     instead of 1,000,001.
    ///     <b>SORTING:</b> Maintains sorted adjacency lists in the result.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public O2M Transpose(int maxNodeCap)
    {
        if (maxNodeCap < 0)
            throw new ArgumentOutOfRangeException(nameof(maxNodeCap),
                "Maximum node cap must be non-negative.");

        return TransposeCore(false, maxNodeCap);
    }

    /// <summary>
    ///     Transposes the structure with strict validation of node indices.
    /// </summary>
    /// <remarks>
    ///     Unlike <see cref="Transpose()" />, this method always throws an exception if any
    ///     node index is negative or out of range, regardless of <see cref="TransposeSkipsInvalidNodes" />.
    ///     Use this when you want to catch data integrity issues early rather than silently dropping connectivity.
    /// </remarks>
    /// <exception cref="InvalidOperationException">
    ///     Thrown when any element contains a negative node index.
    /// </exception>
    [MethodImpl(AggressiveOptimization)]
    public O2M TransposeStrict()
    {
        return TransposeCore(true);
    }

    /// <summary>
    ///     Core transpose implementation with optional strict validation.
    /// </summary>
    /// <param name="strict">If true, throw on negative node indices.</param>
    /// <param name="maxNodeCap">If set, skip nodes exceeding this value.</param>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    private O2M TransposeCore(bool strict, int? maxNodeCap = null)
    {
        var sourceCount = Count;
        if (sourceCount == 0)
            return new O2M { ParallelizationThreshold = ParallelizationThreshold };

        // P1.4 FIX: Validate for negative nodes if strict mode
        if (strict)
            for (var elementIdx = 0; elementIdx < sourceCount; elementIdx++)
            {
                var nodes = CollectionsMarshal.AsSpan(_adjacencies[elementIdx]);
                foreach (var node in nodes)
                    if (node < 0)
                        throw new InvalidOperationException(
                            $"Element {elementIdx} contains negative node index {node}. " +
                            $"Negative indices are not allowed in strict mode. " +
                            $"Use Transpose() instead if you want to silently skip invalid nodes.");
            }

        // P0.2 FIX: CRITICAL - Always recompute maxNode, ignore cache
        // External mutations via indexer or AdjacenciesMutable can invalidate cache
        // Stale cache causes silent loss of connectivity in transposed structure
        // Also respect maxNodeCap if specified to limit memory allocation
        var maxNode = -1;
        for (var i = 0; i < sourceCount; i++)
        {
            var span = CollectionsMarshal.AsSpan(_adjacencies[i]);
            foreach (var node in span)
            {
                if (node < 0) continue; // Skip negative nodes
                if (maxNodeCap.HasValue && node > maxNodeCap.Value) continue; // P0-2: Skip if exceeds cap
                if (node > maxNode)
                    maxNode = node;
            }
        }

        if (maxNode < 0)
            return new O2M { ParallelizationThreshold = ParallelizationThreshold };

        var numberOfNodes = maxNode + 1;
        var result = new O2M(numberOfNodes)
        {
            ParallelizationThreshold = ParallelizationThreshold
        };

        // ============================================================================
        // PHASE 1+2 COMBINED: Sequential counting with duplicate detection + position
        // assignment, followed by allocation.
        //
        // P0 FIX 1.1: Duplicate-node detection per element using cheap O(k²) check.
        // For FEM meshes with k ≤ 20 nodes/element this adds negligible overhead.
        // Duplicates within an element would produce repeated elementIdx entries in the
        // transpose while the old invariant check (offsets == counts) would still pass,
        // silently corrupting downstream algorithms that assume set semantics.
        //
        // P1 FIX 2.1: Flat pooled buffer for element write positions instead of one
        // int[] per element. Reduces O(#elements) small allocations to 1 large allocation,
        // dramatically reducing GC pressure and improving cache locality.
        //
        // Algorithm:
        //   Pass 0: Sequential scan to count unique (non-duplicate) node references per
        //           target node AND record per-element write positions into a flat buffer.
        //   Allocate: Size target arrays from the duplicate-free counts.
        //   Pass 1: Parallel fill using precomputed positions (embarrassingly parallel).
        //
        // This merges counting + duplicate detection + position assignment into one pass,
        // then allocates with exact sizes. The fill phase is fully parallel with no sync.
        // ============================================================================

        var counts = new int[numberOfNodes];
        var writePositions = new int[numberOfNodes];

        // P1 FIX 2.1: Compute element offsets for flat position buffer
        var elemDegrees = new int[sourceCount];
        var elemOffsets = new int[sourceCount + 1];
        for (var i = 0; i < sourceCount; i++)
            elemDegrees[i] = _adjacencies[i].Count;

        // Prefix sum for flat buffer offsets
        elemOffsets[0] = 0;
        for (var i = 0; i < sourceCount; i++)
            elemOffsets[i + 1] = elemOffsets[i] + elemDegrees[i];

        var totalSlots = elemOffsets[sourceCount];

        // Single flat buffer replaces sourceCount individual int[] allocations
        var flatPositions = totalSlots > 0
            ? ArrayPool<int>.Shared.Rent(totalSlots)
            : Array.Empty<int>();

        try
        {
        // Pass 0: Sequential counting + duplicate detection + position assignment
        for (var elementIdx = 0; elementIdx < sourceCount; elementIdx++)
        {
            var nodes = CollectionsMarshal.AsSpan(_adjacencies[elementIdx]);
            var posBase = elemOffsets[elementIdx];

            for (var i = 0; i < nodes.Length; i++)
            {
                var node = nodes[i];
                if ((uint)node >= (uint)numberOfNodes)
                {
                    flatPositions[posBase + i] = -1; // Mark as invalid/skipped
                    continue;
                }

                // P0 FIX 1.1: O(k²) duplicate check within this element (k ≤ ~20 for FEM)
                var isDuplicate = false;
                for (var j = 0; j < i; j++)
                {
                    if (nodes[j] == node)
                    {
                        isDuplicate = true;
                        break;
                    }
                }

                if (isDuplicate)
                {
                    if (strict)
                        throw new InvalidOperationException(
                            $"Element {elementIdx} contains duplicate node index {node}. " +
                            $"Duplicate nodes within an element are not allowed in strict mode. " +
                            $"Use Transpose() with TransposeSkipsInvalidNodes=true to skip duplicates.");

                    // Non-strict: skip duplicate to prevent repeated elementIdx in transpose
                    flatPositions[posBase + i] = -1;
                }
                else
                {
                    counts[node]++;
                    flatPositions[posBase + i] = writePositions[node]++;
                }
            }
        }

        // ============================================================================
        // PHASE 2: Allocate arrays with exact sizes (duplicate-free counts)
        // Uses ARRAYS instead of Lists for guaranteed thread-safe concurrent access
        // ============================================================================
        var resultArrays = new int[numberOfNodes][];
        for (var node = 0; node < numberOfNodes; node++)
            resultArrays[node] = counts[node] > 0
                ? new int[counts[node]]
                : Array.Empty<int>();

        // ============================================================================
        // PHASE 3: Parallel fill using precomputed positions
        // No atomics, no synchronization, no sorting - positions guarantee sorted output
        // ============================================================================
        if (sourceCount >= ParallelizationThreshold)
            Parallel.For(0, sourceCount, ParallelConfig.Options, elementIdx =>
            {
                var nodes = CollectionsMarshal.AsSpan(_adjacencies[elementIdx]);
                var posBase = elemOffsets[elementIdx];

                for (var i = 0; i < nodes.Length; i++)
                {
                    var node = nodes[i];
                    var pos = flatPositions[posBase + i];
                    if (pos >= 0) resultArrays[node][pos] = elementIdx;
                }
            });
        else
            // Sequential fill for small graphs
            for (var elementIdx = 0; elementIdx < sourceCount; elementIdx++)
            {
                var nodes = CollectionsMarshal.AsSpan(_adjacencies[elementIdx]);
                var posBase = elemOffsets[elementIdx];

                for (var i = 0; i < nodes.Length; i++)
                {
                    var node = nodes[i];
                    var pos = flatPositions[posBase + i];
                    if (pos >= 0) resultArrays[node][pos] = elementIdx;
                }
            }

        // Copy final positions for invariant check
        var offsets = writePositions;

        // ============================================================================
        // INVARIANT CHECK: Verify offsets match allocated counts
        // A mismatch indicates data corruption (concurrent mutation, or algorithm bug).
        // Duplicate nodes are already handled in Pass 0 above.
        // ============================================================================
        for (var node = 0; node < numberOfNodes; node++)
            if (offsets[node] != resultArrays[node].Length)
                throw new InvalidOperationException(
                    $"Transpose invariant violation at node {node}: expected {resultArrays[node].Length} elements, " +
                    $"but filled {offsets[node]}. This indicates either: " +
                    $"(1) concurrent modification of the O2M during transpose, or " +
                    $"(2) a bug in the transpose algorithm.");

        // ============================================================================
        // PHASE 4: Convert arrays to Lists
        // Now that arrays are filled and sorted, convert to O2M's internal format
        // ============================================================================
        for (var node = 0; node < numberOfNodes; node++)
        {
            var array = resultArrays[node];
            if (array.Length == 0)
            {
                result._adjacencies.Add(new List<int>());
            }
            else
            {
                var list = new List<int>(array.Length);
                CollectionsMarshal.SetCount(list, array.Length);
                var listSpan = CollectionsMarshal.AsSpan(list);
                array.AsSpan().CopyTo(listSpan);
                result._adjacencies.Add(list);
            }
        }

        } // end try
        finally
        {
            // P1 FIX 2.1: Return flat positions buffer to ArrayPool
            if (totalSlots > 0)
                ArrayPool<int>.Shared.Return(flatPositions);
        }

        // ============================================================================
        // Cache management
        // The transpose's "nodes" are element indices from the source.
        // Leave cache unknown (null) to be computed on demand.
        // ============================================================================
        result._maxNodeIndexCache = null;

        return result;
    }

    #endregion

    #region Fields and Pooling

    private static readonly ObjectPool<HashSet<int>> _hashSetPool =
        ObjectPool.Create(new HashSetPoolPolicy());

    [JsonInclude] private List<List<int>> _adjacencies;
    private int? _maxNodeIndexCache;

    public static int DefaultParallelizationThreshold { get; set; } = 4096;
    public int ParallelizationThreshold { get; set; }

    /// <summary>
    ///     Controls how Transpose() handles invalid node indices (negative or out of range).
    ///     If true (default), invalid indices are silently skipped (legacy behavior).
    ///     If false, invalid indices cause InvalidOperationException via TransposeStrict().
    /// </summary>
    public bool TransposeSkipsInvalidNodes { get; set; } = true;

    #endregion

    #region Constructor

    public O2M()
    {
        _adjacencies = [];
        _maxNodeIndexCache = null;
        ParallelizationThreshold = DefaultParallelizationThreshold;
    }

    public O2M(int reservedCapacity)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(reservedCapacity);
        _adjacencies = new List<List<int>>(reservedCapacity);
        _maxNodeIndexCache = null;
        ParallelizationThreshold = DefaultParallelizationThreshold;
    }

    public O2M(List<List<int>> adjacenciesList)
    {
        ArgumentNullException.ThrowIfNull(adjacenciesList);
        _adjacencies = adjacenciesList;
        _maxNodeIndexCache = null;
        ParallelizationThreshold = DefaultParallelizationThreshold;
    }

    #endregion

    #region Properties and Indexers

    public int Count => _adjacencies.Count;

    /// <summary>
    ///     Gets the nodes for the specified element as a read-only span.
    /// </summary>
    /// <param name="rowIndex">The element index.</param>
    /// <returns>Read-only span of node indices. Cannot be modified.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when rowIndex is out of range.</exception>
    /// <remarks>
    ///     P0.3 FIX: Returns ReadOnlySpan instead of IReadOnlyList to prevent
    ///     external mutation via downcast to List&lt;int&gt;. Previous implementation
    ///     exposed internal List which could be cast and mutated, breaking invariants
    ///     and invalidating caches.
    /// </remarks>
    public ReadOnlySpan<int> this[int rowIndex]
    {
        get
        {
            if ((uint)rowIndex >= (uint)_adjacencies.Count)
                throw new ArgumentOutOfRangeException(nameof(rowIndex));
            return CollectionsMarshal.AsSpan(_adjacencies[rowIndex]);
        }
    }

    public int this[int rowIndex, int columnIndex]
    {
        get
        {
            if ((uint)rowIndex >= (uint)_adjacencies.Count)
                throw new ArgumentOutOfRangeException(nameof(rowIndex));
            var row = _adjacencies[rowIndex];
            if ((uint)columnIndex >= (uint)row.Count)
                throw new ArgumentOutOfRangeException(nameof(columnIndex));
            return row[columnIndex];
        }
    }

    internal List<List<int>> AdjacenciesMutable => _adjacencies;

    #endregion

    #region Core Modification Methods

    public void ClearElement(int elementIndex)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(elementIndex);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(elementIndex, Count);

        if (_adjacencies[elementIndex].Count > 0)
            _maxNodeIndexCache = null;

        _adjacencies[elementIndex].Clear();
    }

    /// <summary>
    ///     Replaces nodes for an element. Takes ownership of the list — caller must not mutate it afterward.
    /// </summary>
    /// <remarks>
    ///     <b>P0.1 FIX:</b> Renamed from ReplaceElement to make ownership transfer explicit.
    ///     The passed list is stored directly (no copy). Caller must not modify it after this call.
    ///     For a safe alternative, use <see cref="ReplaceElementCopy"/>.
    /// </remarks>
    public void ReplaceElement(int elementIndex, List<int> newNodes)
    {
        ArgumentNullException.ThrowIfNull(newNodes);
        ArgumentOutOfRangeException.ThrowIfNegative(elementIndex);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(elementIndex, Count);

        _adjacencies[elementIndex] = newNodes;
        _maxNodeIndexCache = null;
    }

    /// <summary>
    ///     Replaces nodes for an element by copying from a span. Safe for external API boundaries.
    /// </summary>
    /// <remarks>
    ///     <b>P0.1 FIX:</b> Safe alternative that copies data, so the caller retains full ownership
    ///     of their buffer. Uses <see cref="CollectionsMarshal.SetCount"/> for zero-overhead allocation.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public void ReplaceElementCopy(int elementIndex, ReadOnlySpan<int> newNodes)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(elementIndex);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(elementIndex, Count);

        var list = new List<int>(newNodes.Length);
        if (newNodes.Length > 0)
        {
            CollectionsMarshal.SetCount(list, newNodes.Length);
            newNodes.CopyTo(CollectionsMarshal.AsSpan(list));
        }

        _adjacencies[elementIndex] = list;
        _maxNodeIndexCache = null;
    }

    /// <summary>
    ///     Appends a new element with the specified nodes.
    ///     Takes ownership of the list — caller must not mutate it afterward.
    ///     SIMD max tracking.
    /// </summary>
    /// <remarks>
    ///     <b>P0.1 FIX:</b> The passed list is stored directly (no copy). Caller must not modify it
    ///     after this call or caches and invariants will be silently corrupted.
    ///     For a safe alternative, use <see cref="AppendElementCopy"/>.
    /// </remarks>
    [MethodImpl(AggressiveOptimization | AggressiveInlining)]
    public int AppendElement(List<int> nodes)
    {
        ArgumentNullException.ThrowIfNull(nodes);

        _adjacencies.Add(nodes);

        if (_maxNodeIndexCache.HasValue)
        {
            var span = CollectionsMarshal.AsSpan(nodes);
            var length = span.Length;

            if (length > 0)
            {
                var currentMax = _maxNodeIndexCache.Value;
                currentMax = FindMaxInSpan(span, currentMax);
                _maxNodeIndexCache = currentMax;
            }
        }

        return Count - 1;
    }

    /// <summary>
    ///     Appends a new element by copying from a span. Safe for external API boundaries.
    /// </summary>
    /// <remarks>
    ///     <b>P0.1 FIX:</b> Safe alternative that copies data, so the caller retains full ownership
    ///     of their buffer.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public int AppendElementCopy(ReadOnlySpan<int> nodes)
    {
        var list = new List<int>(nodes.Length);
        if (nodes.Length > 0)
        {
            CollectionsMarshal.SetCount(list, nodes.Length);
            nodes.CopyTo(CollectionsMarshal.AsSpan(list));
        }

        return AppendElement(list);
    }

    /// <summary>
    ///     Appends multiple elements in a batch.
    ///     Parallel max computation for large batches.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public void AppendElements(params List<int>[] nodes)
    {
        ArgumentNullException.ThrowIfNull(nodes);
        if (nodes.Length == 0) return;

        foreach (var nodeList in nodes)
        {
            if (nodeList is null)
                throw new ArgumentNullException(nameof(nodes), "Cannot append null node list.");
            _adjacencies.Add(nodeList);
        }

        // Parallel max computation for large batches
        // P0.3 FIX: Use stable per-worker slot assignment instead of threadId % slots
        // to avoid false sharing and CAS contention from collisions
        if (_maxNodeIndexCache.HasValue && nodes.Length >= 64)
        {
            var currentMax = _maxNodeIndexCache.Value;
            var localMaxes = new int[ParallelConfig.MaxDegreeOfParallelism];
            Array.Fill(localMaxes, currentMax);
            var nextSlot = 0;

            Parallel.For(0, nodes.Length,
                ParallelConfig.Options,
                () =>
                {
                    // P0.3 FIX: Each worker gets a unique, stable slot
                    var slot = Interlocked.Increment(ref nextSlot) - 1;
                    slot %= localMaxes.Length;
                    return (slot, localMax: currentMax);
                },
                (i, _, state) =>
                {
                    var span = CollectionsMarshal.AsSpan(nodes[i]);
                    return (state.slot, localMax: FindMaxInSpan(span, state.localMax));
                },
                state =>
                {
                    int current;
                    do
                    {
                        current = Volatile.Read(ref localMaxes[state.slot]);
                        if (state.localMax <= current) break;
                    } while (Interlocked.CompareExchange(ref localMaxes[state.slot], state.localMax, current) != current);
                });

            _maxNodeIndexCache = localMaxes.Max();
        }
        else
        {
            _maxNodeIndexCache = null;
        }
    }

    /// <summary>
    ///     SIMD-accelerated max finding in span.
    /// </summary>
    [MethodImpl(AggressiveOptimization | AggressiveInlining)]
    private static unsafe int FindMaxInSpan(ReadOnlySpan<int> span, int currentMax)
    {
        if (span.Length == 0) return currentMax;

        fixed (int* ptr = span)
        {
            var length = span.Length;

            if (Avx2.IsSupported && length >= 8)
            {
                var maxVec = Vector256.Create(currentMax);
                var i = 0;

                for (; i + 7 < length; i += 8)
                {
                    var vec = Avx.LoadVector256(ptr + i);
                    maxVec = Avx2.Max(maxVec, vec);
                }

                // Horizontal max of vector
                var temp = stackalloc int[8];
                Avx.Store(temp, maxVec);
                for (var j = 0; j < 8; j++)
                    if (temp[j] > currentMax)
                        currentMax = temp[j];

                // Handle remainder
                for (; i < length; i++)
                    if (ptr[i] > currentMax)
                        currentMax = ptr[i];
            }
            else if (Vector.IsHardwareAccelerated && length >= Vector<int>.Count)
            {
                var maxVec = new Vector<int>(currentMax);
                var i = 0;

                for (; i + Vector<int>.Count - 1 < length; i += Vector<int>.Count)
                {
                    var vec = new Vector<int>(span.Slice(i, Vector<int>.Count));
                    maxVec = Vector.Max(maxVec, vec);
                }

                for (var j = 0; j < Vector<int>.Count; j++)
                    if (maxVec[j] > currentMax)
                        currentMax = maxVec[j];

                for (; i < length; i++)
                    if (ptr[i] > currentMax)
                        currentMax = ptr[i];
            }
            else
            {
                // Unrolled scalar
                var i = 0;
                for (; i + 3 < length; i += 4)
                {
                    if (ptr[i] > currentMax) currentMax = ptr[i];
                    if (ptr[i + 1] > currentMax) currentMax = ptr[i + 1];
                    if (ptr[i + 2] > currentMax) currentMax = ptr[i + 2];
                    if (ptr[i + 3] > currentMax) currentMax = ptr[i + 3];
                }

                for (; i < length; i++)
                    if (ptr[i] > currentMax)
                        currentMax = ptr[i];
            }
        }

        return currentMax;
    }

    public void AppendNodeToElement(int elementIndex, int nodeValue)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(elementIndex);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(elementIndex, Count);
        ArgumentOutOfRangeException.ThrowIfNegative(nodeValue);

        _adjacencies[elementIndex].Add(nodeValue);

        if (_maxNodeIndexCache.HasValue && nodeValue > _maxNodeIndexCache.Value)
            _maxNodeIndexCache = nodeValue;
    }

    public bool RemoveNodeFromElement(int elementIndex, int nodeValue)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(elementIndex);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(elementIndex, Count);

        var removed = _adjacencies[elementIndex].Remove(nodeValue);
        if (removed) _maxNodeIndexCache = null;
        return removed;
    }

    public void ClearAll()
    {
        _adjacencies.Clear();
        _maxNodeIndexCache = null;
    }

    public void Reserve(int reservedCapacity)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(reservedCapacity);
        if (reservedCapacity > _adjacencies.Capacity)
            _adjacencies.Capacity = reservedCapacity;
    }

    /// <summary>
    ///     Reduces memory usage by trimming excess capacity.
    ///     Parallel trimming for large structures.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public void ShrinkToFit()
    {
        _adjacencies.TrimExcess();

        var count = _adjacencies.Count;
        if (count >= ParallelizationThreshold)
            Parallel.For(0, count, ParallelConfig.Options, i => _adjacencies[i].TrimExcess());
        else
            for (var i = 0; i < count; i++)
                _adjacencies[i].TrimExcess();
    }

    #endregion

    #region Query Methods

    /// <summary>
    ///     Gets the maximum node index across all elements.
    ///     Parallel SIMD search.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    public int GetMaxNode()
    {
        if (_maxNodeIndexCache.HasValue) return _maxNodeIndexCache.Value;

        if (_adjacencies.Count == 0)
        {
            _maxNodeIndexCache = -1;
            return -1;
        }

        var count = _adjacencies.Count;
        int max;

        if (count >= ParallelizationThreshold)
        {
            // FIXED: Use proper partitioning instead of thread-ID bucketing
            // The old code used thread-ID % numThreads which causes collisions
            var partitioner = Partitioner.Create(0, count);
            var localMaxes = new ConcurrentBag<int>();

            Parallel.ForEach(partitioner, ParallelConfig.Options, (range, state) =>
            {
                var localMax = int.MinValue;
                for (var i = range.Item1; i < range.Item2; i++)
                {
                    var span = CollectionsMarshal.AsSpan(_adjacencies[i]);
                    localMax = FindMaxInSpan(span, localMax);
                }

                localMaxes.Add(localMax);
            });

            // Merge partition results
            max = int.MinValue;
            foreach (var localMax in localMaxes)
                if (localMax > max)
                    max = localMax;
        }
        else
        {
            max = int.MinValue;
            for (var i = 0; i < count; i++)
            {
                var span = CollectionsMarshal.AsSpan(_adjacencies[i]);
                max = FindMaxInSpan(span, max);
            }
        }

        _maxNodeIndexCache = max == int.MinValue ? -1 : max;
        return _maxNodeIndexCache.Value;
    }

    /// <summary>
    ///     Validates that all node indices are non-negative and unique within each element.
    ///     Parallel validation with early termination.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public bool IsValid()
    {
        var count = _adjacencies.Count;
        if (count == 0) return true;

        if (count >= ParallelizationThreshold)
        {
            var isValid = 1; // Use int for Interlocked

            Parallel.For(0, count, ParallelConfig.Options, (i, loopState) =>
            {
                if (Volatile.Read(ref isValid) == 0)
                {
                    loopState.Stop();
                    return;
                }

                var set = _hashSetPool.Get();
                try
                {
                    var span = CollectionsMarshal.AsSpan(_adjacencies[i]);
                    foreach (var value in span)
                        if (value < 0 || !set.Add(value))
                        {
                            Interlocked.Exchange(ref isValid, 0);
                            loopState.Stop();
                            return;
                        }
                }
                finally
                {
                    _hashSetPool.Return(set);
                }
            });

            return isValid == 1;
        }

        {
            var set = _hashSetPool.Get();
            try
            {
                foreach (var row in _adjacencies)
                {
                    set.Clear();
                    var span = CollectionsMarshal.AsSpan(row);
                    foreach (var value in span)
                        if (value < 0 || !set.Add(value))
                            return false;
                }

                return true;
            }
            finally
            {
                _hashSetPool.Return(set);
            }
        }
    }

    /// <summary>
    ///     Checks whether all adjacency lists are sorted in ascending order.
    /// </summary>
    /// <returns>True if all lists are sorted; false otherwise.</returns>
    /// <remarks>
    ///     <b>INVARIANT:</b> Many algorithms (binary search, position computation, clique detection)
    ///     assume sorted adjacency lists. This method verifies that invariant.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public bool IsSorted()
    {
        var count = _adjacencies.Count;
        if (count == 0) return true;

        if (count >= ParallelizationThreshold)
        {
            var isSorted = 1;

            Parallel.For(0, count, ParallelConfig.Options, (i, loopState) =>
            {
                if (Volatile.Read(ref isSorted) == 0)
                {
                    loopState.Stop();
                    return;
                }

                var span = CollectionsMarshal.AsSpan(_adjacencies[i]);
                for (var j = 1; j < span.Length; j++)
                    if (span[j - 1] >= span[j]) // Not strictly ascending
                    {
                        Interlocked.Exchange(ref isSorted, 0);
                        loopState.Stop();
                        return;
                    }
            });

            return isSorted == 1;
        }

        foreach (var row in _adjacencies)
        {
            var span = CollectionsMarshal.AsSpan(row);
            for (var j = 1; j < span.Length; j++)
                if (span[j - 1] >= span[j])
                    return false;
        }

        return true;
    }

    /// <summary>
    ///     Performs comprehensive validation and returns detailed error information.
    /// </summary>
    /// <returns>Null if valid; otherwise a description of the first error found.</returns>
    /// <remarks>
    ///     Checks:
    ///     <list type="bullet">
    ///         <item>All node indices are non-negative</item>
    ///         <item>No duplicate nodes within any element</item>
    ///         <item>All adjacency lists are sorted in ascending order</item>
    ///     </list>
    /// </remarks>
    public string? ValidateStrict()
    {
        var seen = new HashSet<int>();
        for (var i = 0; i < _adjacencies.Count; i++)
        {
            seen.Clear();
            var span = CollectionsMarshal.AsSpan(_adjacencies[i]);

            for (var j = 0; j < span.Length; j++)
            {
                var node = span[j];

                if (node < 0)
                    return $"Element {i} has negative node index {node} at position {j}.";

                if (!seen.Add(node))
                    return $"Element {i} has duplicate node {node} at position {j}.";

                if (j > 0 && span[j - 1] >= node)
                    return
                        $"Element {i} is not sorted: node {span[j - 1]} at position {j - 1} >= node {node} at position {j}.";
            }
        }

        return null; // Valid
    }

    /// <summary>
    ///     Counts total number of edges (node references) across all elements.
    ///     NEW METHOD: Useful for memory estimation and statistics.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public long GetTotalEdgeCount()
    {
        var count = _adjacencies.Count;
        if (count == 0) return 0;

        if (count >= ParallelizationThreshold)
        {
            long total = 0;
            Parallel.For(0, count,
                ParallelConfig.Options,
                () => 0L,
                (i, _, localSum) => localSum + _adjacencies[i].Count,
                localSum => Interlocked.Add(ref total, localSum));
            return total;
        }
        else
        {
            long total = 0;
            foreach (var row in _adjacencies)
                total += row.Count;
            return total;
        }
    }

    /// <summary>
    ///     Gets statistics about the structure.
    ///     NEW METHOD: Returns min, max, average degree and other metrics.
    /// </summary>
    public (int MinDegree, int MaxDegree, double AvgDegree, long TotalEdges) GetStatistics()
    {
        var count = _adjacencies.Count;
        if (count == 0) return (0, 0, 0, 0);

        var minDeg = int.MaxValue;
        var maxDeg = int.MinValue;
        long total = 0;

        if (count >= ParallelizationThreshold)
        {
            // P0.3 FIX: Removed Environment.CurrentManagedThreadId % numThreads bucketing.
            // Thread ID collisions cause false sharing and CAS contention (same bug fixed in
            // AppendElements). The finalize callback runs only once per worker thread, so
            // direct Interlocked merge has negligible contention.
            var sharedMin = int.MaxValue;
            var sharedMax = int.MinValue;
            long sharedSum = 0;

            Parallel.For(0, count, ParallelConfig.Options,
                () => (Min: int.MaxValue, Max: int.MinValue, Sum: 0L),
                (i, _, local) =>
                {
                    var deg = _adjacencies[i].Count;
                    return (Math.Min(local.Min, deg), Math.Max(local.Max, deg), local.Sum + deg);
                },
                local =>
                {
                    // Atomic min
                    int current;
                    do
                    {
                        current = Volatile.Read(ref sharedMin);
                    } while (local.Min < current &&
                             Interlocked.CompareExchange(ref sharedMin, local.Min, current) != current);

                    // Atomic max
                    do
                    {
                        current = Volatile.Read(ref sharedMax);
                    } while (local.Max > current &&
                             Interlocked.CompareExchange(ref sharedMax, local.Max, current) != current);

                    Interlocked.Add(ref sharedSum, local.Sum);
                });

            minDeg = sharedMin;
            maxDeg = sharedMax;
            total = sharedSum;
        }
        else
        {
            foreach (var row in _adjacencies)
            {
                var deg = row.Count;
                if (deg < minDeg) minDeg = deg;
                if (deg > maxDeg) maxDeg = deg;
                total += deg;
            }
        }

        return (minDeg, maxDeg, (double)total / count, total);
    }

    #endregion

    #region Graph Operations

    /// <summary>
    ///     Tests if the structure represents an acyclic graph.
    /// </summary>
    public bool IsAcyclic()
    {
        var maxNodeValue = GetMaxNode();
        if (maxNodeValue < 0) return true;

        var nodeCount = Math.Max(Count, maxNodeValue + 1);
        var state = new byte[nodeCount];
        var stack = new Stack<int>(Math.Min(nodeCount, 1024));

        for (var start = 0; start < nodeCount; start++)
        {
            if (state[start] != 0) continue;

            stack.Push(start);

            while (stack.Count > 0)
            {
                var u = stack.Peek();

                if (state[u] == 0)
                {
                    state[u] = 1;

                    if (u < Count)
                    {
                        var span = CollectionsMarshal.AsSpan(_adjacencies[u]);
                        foreach (var v in span)
                        {
                            if ((uint)v >= (uint)nodeCount) continue;
                            if (state[v] == 1) return false;
                            if (state[v] == 0) stack.Push(v);
                        }
                    }
                }
                else
                {
                    stack.Pop();
                    state[u] = 2;
                }
            }
        }

        return true;
    }

    /// <summary>
    ///     Computes a topological ordering of the graph nodes.
    ///     Parallel in-degree calculation.
    /// </summary>
    public List<int> GetTopOrder()
    {
        var order = new List<int>();
        var inDegree = new int[Count];

        if (Count >= ParallelizationThreshold)
            Parallel.For(0, Count, ParallelConfig.Options, i =>
            {
                var nodes = _adjacencies[i];
                foreach (var node in nodes)
                    if ((uint)node < (uint)Count)
                        Interlocked.Increment(ref inDegree[node]);
            });
        else
            for (var i = 0; i < Count; i++)
            {
                var nodes = _adjacencies[i];
                foreach (var node in nodes)
                    if ((uint)node < (uint)Count)
                        inDegree[node]++;
            }

        var queue = new Queue<int>();
        for (var i = 0; i < inDegree.Length; i++)
            if (inDegree[i] == 0)
                queue.Enqueue(i);

        while (queue.Count > 0)
        {
            var cur = queue.Dequeue();
            order.Add(cur);

            var nodes = _adjacencies[cur];
            foreach (var nbr in nodes)
                if ((uint)nbr < (uint)Count)
                    if (--inDegree[nbr] == 0)
                        queue.Enqueue(nbr);
        }

        if (order.Count != Count)
            throw new InvalidOperationException("The relation contains cycles, topological sort not possible.");

        return order;
    }

    #endregion

    #region Sorting and Ordering

    /// <summary>
    ///     Computes a sort order for elements based on lexicographic comparison.
    ///     Uses radix sort - O(n × k) instead of O(n × k × log n).
    ///     Parallelized for large datasets.
    /// </summary>
    /// <remarks>
    ///     Algorithm: Two-phase radix sort
    ///     1. Counting sort by row length (groups rows of same size)
    ///     2. LSD radix sort within each length group by column values
    ///     Complexity: O(n × k) where n = row count, k = max row length.
    ///     This is optimal for bounded node indices.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public List<int> GetSortOrder()
    {
        var n = _adjacencies.Count;
        if (n == 0) return [];
        if (n == 1) return [0];

        // Phase 0: Find max row length and max node value (parallel for large n)
        int maxLen = 0, maxNode = 0;

        if (n >= ParallelizationThreshold)
        {
            var localMaxLen = 0;
            var localMaxNode = 0;
            var lockObj = new object();

            Parallel.ForEach(
                Partitioner.Create(0, n),
                ParallelConfig.Options,
                () => (maxLen: 0, maxNode: 0),
                (range, _, local) =>
                {
                    for (var i = range.Item1; i < range.Item2; i++)
                    {
                        var row = _adjacencies[i];
                        if (row.Count > local.maxLen) local.maxLen = row.Count;
                        var span = CollectionsMarshal.AsSpan(row);
                        foreach (var v in span)
                            if (v > local.maxNode)
                                local.maxNode = v;
                    }

                    return local;
                },
                local =>
                {
                    lock (lockObj)
                    {
                        if (local.maxLen > localMaxLen) localMaxLen = local.maxLen;
                        if (local.maxNode > localMaxNode) localMaxNode = local.maxNode;
                    }
                });

            maxLen = localMaxLen;
            maxNode = localMaxNode;
        }
        else
        {
            for (var i = 0; i < n; i++)
            {
                var row = _adjacencies[i];
                if (row.Count > maxLen) maxLen = row.Count;
                var span = CollectionsMarshal.AsSpan(row);
                foreach (var v in span)
                    if (v > maxNode)
                        maxNode = v;
            }
        }

        // Handle edge case: all empty rows
        if (maxLen == 0)
        {
            var result = new List<int>(n);
            for (var i = 0; i < n; i++) result.Add(i);
            return result;
        }

        // Phase 1: Counting sort by length
        var indices = CountingSortByLength(n, maxLen);

        // Phase 2: Radix sort within each length group
        var start = 0;
        while (start < n)
        {
            var len = _adjacencies[indices[start]].Count;
            var end = start + 1;
            while (end < n && _adjacencies[indices[end]].Count == len)
                end++;

            // Radix sort this group by columns (LSD: right to left)
            if (end - start > 1 && len > 0)
                RadixSortByColumns(indices, start, end, len, maxNode);

            start = end;
        }

        return new List<int>(indices);
    }

    /// <summary>
    ///     Counting sort indices by row length. O(n + maxLen).
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    private int[] CountingSortByLength(int n, int maxLen)
    {
        var count = new int[maxLen + 2];

        // Count lengths (parallel for large n)
        if (n >= ParallelizationThreshold)
        {
            var localCounts = new ThreadLocal<int[]>(() => new int[maxLen + 2], true);

            Parallel.ForEach(Partitioner.Create(0, n), ParallelConfig.Options, range =>
            {
                var local = localCounts.Value!;
                for (var i = range.Item1; i < range.Item2; i++)
                    local[_adjacencies[i].Count + 1]++;
            });

            // Merge thread-local counts
            foreach (var local in localCounts.Values)
                for (var i = 0; i < count.Length; i++)
                    count[i] += local[i];

            localCounts.Dispose();
        }
        else
        {
            for (var i = 0; i < n; i++)
                count[_adjacencies[i].Count + 1]++;
        }

        // Prefix sum
        for (var i = 1; i < count.Length; i++)
            count[i] += count[i - 1];

        // Place elements (sequential - maintains stability)
        var output = GC.AllocateUninitializedArray<int>(n);
        for (var i = 0; i < n; i++)
        {
            var len = _adjacencies[i].Count;
            output[count[len]++] = i;
        }

        return output;
    }

    /// <summary>
    ///     LSD Radix sort a range of indices by column values.
    ///     Processes columns from right to left for stable sort.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    private void RadixSortByColumns(int[] indices, int start, int end, int numCols, int maxNode)
    {
        var n = end - start;
        var bucketCount = maxNode + 2; // +1 for zero-indexing, +1 for prefix sum space
        var temp = GC.AllocateUninitializedArray<int>(n);

        // For large sorts with many buckets, use parallel counting
        var useParallel = n >= ParallelizationThreshold && bucketCount <= 100000;

        // LSD radix sort: from last column to first (ensures stability)
        for (var col = numCols - 1; col >= 0; col--)
        {
            int[] count;

            if (useParallel)
            {
                // Parallel counting phase
                var localCounts = new ThreadLocal<int[]>(() => new int[bucketCount], true);

                Parallel.ForEach(Partitioner.Create(start, end), ParallelConfig.Options, range =>
                {
                    var local = localCounts.Value!;
                    for (var i = range.Item1; i < range.Item2; i++)
                    {
                        var val = _adjacencies[indices[i]][col];
                        local[val + 1]++;
                    }
                });

                // Merge counts
                count = new int[bucketCount];
                foreach (var local in localCounts.Values)
                    for (var i = 0; i < bucketCount; i++)
                        count[i] += local[i];

                localCounts.Dispose();
            }
            else
            {
                // Sequential counting
                count = new int[bucketCount];
                for (var i = start; i < end; i++)
                {
                    var val = _adjacencies[indices[i]][col];
                    count[val + 1]++;
                }
            }

            // Prefix sum (always sequential - small array)
            for (var i = 1; i < bucketCount; i++)
                count[i] += count[i - 1];

            // Place elements (sequential - stability requires order preservation)
            for (var i = start; i < end; i++)
            {
                var val = _adjacencies[indices[i]][col];
                temp[count[val]++] = indices[i];
            }

            // Copy back
            Array.Copy(temp, 0, indices, start, n);
        }
    }

    /// <summary>
    ///     Finds indices of duplicate elements.
    ///     Parallel duplicate detection after sorting.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public List<int> GetDuplicates()
    {
        var n = Count;
        if (n <= 1) return [];

        var sorted = GetSortOrder();

        if (n >= ParallelizationThreshold)
        {
            // Parallel detection with concurrent collection
            var duplicates = new ConcurrentBag<int>();

            Parallel.For(0, n - 1, ParallelConfig.Options, k =>
            {
                var row1 = _adjacencies[sorted[k]];
                var row2 = _adjacencies[sorted[k + 1]];
                if (CompareRows(row1, row2) == 0)
                    duplicates.Add(sorted[k + 1]);
            });

            // Review Issue #5: Sort for deterministic output
            var result = duplicates.ToList();
            result.Sort();
            return result;
        }
        else
        {
            var duplicates = new List<int>();
            for (var k = 0; k < n - 1; k++)
                if (CompareRows(_adjacencies[sorted[k]], _adjacencies[sorted[k + 1]]) == 0)
                    duplicates.Add(sorted[k + 1]);
            return duplicates;
        }
    }

    /// <summary>
    ///     Compares two rows lexicographically.
    ///     Uses SIMD comparison where possible.
    /// </summary>
    [MethodImpl(AggressiveOptimization | AggressiveInlining)]
    private static int CompareRows(List<int> left, List<int> right)
    {
        var comparison = left.Count.CompareTo(right.Count);
        if (comparison != 0) return comparison;

        return CollectionsMarshal.AsSpan(left)
            .SequenceCompareTo(CollectionsMarshal.AsSpan(right));
    }

    #endregion

    #region Permutation and Compression

    /// <summary>
    ///     Reorders elements according to a permutation mapping.
    ///     Parallel reordering.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public void PermuteElements(List<int> oldToNewElementMap)
    {
        ArgumentNullException.ThrowIfNull(oldToNewElementMap);

        var n = Count;
        if (oldToNewElementMap.Count != n || n == 0) return;

        // Validate permutation (parallel for large n)
        var seen = new int[n]; // Use int for Interlocked
        var isValid = 1;

        if (n >= ParallelizationThreshold)
            Parallel.For(0, n, ParallelConfig.Options, (oldIdx, loopState) =>
            {
                if (Volatile.Read(ref isValid) == 0)
                {
                    loopState.Stop();
                    return;
                }

                var newIdx = oldToNewElementMap[oldIdx];
                if ((uint)newIdx >= (uint)n || Interlocked.Exchange(ref seen[newIdx], 1) == 1)
                {
                    Interlocked.Exchange(ref isValid, 0);
                    loopState.Stop();
                }
            });
        else
            for (var oldIdx = 0; oldIdx < n; oldIdx++)
            {
                var newIdx = oldToNewElementMap[oldIdx];
                if ((uint)newIdx >= (uint)n || seen[newIdx] != 0)
                {
                    isValid = 0;
                    break;
                }

                seen[newIdx] = 1;
            }

        if (isValid == 0)
            throw new InvalidOperationException("PermuteElements requires a valid permutation mapping.");

        var reordered = new List<int>[n];

        if (n >= ParallelizationThreshold)
            Parallel.For(0, n, ParallelConfig.Options,
                oldIdx => { reordered[oldToNewElementMap[oldIdx]] = _adjacencies[oldIdx]; });
        else
            for (var oldIdx = 0; oldIdx < n; oldIdx++)
                reordered[oldToNewElementMap[oldIdx]] = _adjacencies[oldIdx];

        _adjacencies = new List<List<int>>(reordered);
        // _maxNodeIndexCache remains valid: node values are unchanged by element reordering
    }

    /// <summary>
    ///     Renumbers all nodes according to a mapping.
    ///     Parallel with SIMD where applicable.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public unsafe void PermuteNodes(List<int> oldToNewNodeMap)
    {
        ArgumentNullException.ThrowIfNull(oldToNewNodeMap);

        var count = _adjacencies.Count;
        var mapSpan = CollectionsMarshal.AsSpan(oldToNewNodeMap);
        var mapLength = mapSpan.Length;

        if (count >= ParallelizationThreshold)
            Parallel.For(0, count, ParallelConfig.Options, i =>
            {
                var row = _adjacencies[i];
                var rowSpan = CollectionsMarshal.AsSpan(row);
                var localMapSpan = CollectionsMarshal.AsSpan(oldToNewNodeMap);

                fixed (int* rowPtr = rowSpan)
                fixed (int* mapPtr = localMapSpan)
                {
                    var len = rowSpan.Length;
                    var j = 0;

                    // Unrolled permutation
                    for (; j + 3 < len; j += 4)
                    {
                        var idx0 = rowPtr[j];
                        var idx1 = rowPtr[j + 1];
                        var idx2 = rowPtr[j + 2];
                        var idx3 = rowPtr[j + 3];

                        rowPtr[j] = (uint)idx0 < (uint)mapLength ? mapPtr[idx0] : idx0;
                        rowPtr[j + 1] = (uint)idx1 < (uint)mapLength ? mapPtr[idx1] : idx1;
                        rowPtr[j + 2] = (uint)idx2 < (uint)mapLength ? mapPtr[idx2] : idx2;
                        rowPtr[j + 3] = (uint)idx3 < (uint)mapLength ? mapPtr[idx3] : idx3;
                    }

                    for (; j < len; j++)
                    {
                        var idx = rowPtr[j];
                        rowPtr[j] = (uint)idx < (uint)mapLength ? mapPtr[idx] : idx;
                    }
                }
            });
        else
            fixed (int* mapPtr = mapSpan)
            {
                for (var i = 0; i < count; i++)
                {
                    var rowSpan = CollectionsMarshal.AsSpan(_adjacencies[i]);
                    fixed (int* rowPtr = rowSpan)
                    {
                        for (var j = 0; j < rowSpan.Length; j++)
                        {
                            var idx = rowPtr[j];
                            rowPtr[j] = (uint)idx < (uint)mapLength ? mapPtr[idx] : idx;
                        }
                    }
                }
            }

        _maxNodeIndexCache = null;
    }

    /// <summary>
    ///     Compresses elements by selecting a subset according to a mapping.
    ///     Parallel copying with span operations.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public void CompressElements(List<int> newToOldElementMap)
    {
        ArgumentNullException.ThrowIfNull(newToOldElementMap);

        var newCount = newToOldElementMap.Count;
        if (newCount == 0)
        {
            _adjacencies.Clear();
            _maxNodeIndexCache = null;
            return;
        }

        var compressed = GC.AllocateUninitializedArray<List<int>>(newCount);
        var oldCount = _adjacencies.Count;

        if (newCount >= ParallelizationThreshold)
            Parallel.For(0, newCount, ParallelConfig.Options, i =>
            {
                var oldIdx = newToOldElementMap[i];
                if ((uint)oldIdx < (uint)oldCount)
                {
                    var src = _adjacencies[oldIdx];
                    var dst = new List<int>(src.Count);
                    if (src.Count > 0)
                    {
                        CollectionsMarshal.SetCount(dst, src.Count);
                        CollectionsMarshal.AsSpan(src).CopyTo(CollectionsMarshal.AsSpan(dst));
                    }

                    compressed[i] = dst;
                }
                else
                {
                    compressed[i] = new List<int>();
                }
            });
        else
            for (var i = 0; i < newCount; i++)
            {
                var oldIdx = newToOldElementMap[i];
                if ((uint)oldIdx < (uint)oldCount)
                {
                    var src = _adjacencies[oldIdx];
                    var dst = new List<int>(src.Count);
                    if (src.Count > 0)
                    {
                        CollectionsMarshal.SetCount(dst, src.Count);
                        CollectionsMarshal.AsSpan(src).CopyTo(CollectionsMarshal.AsSpan(dst));
                    }

                    compressed[i] = dst;
                }
                else
                {
                    compressed[i] = new List<int>();
                }
            }

        _adjacencies.Clear();
        _adjacencies.AddRange(compressed);
        _maxNodeIndexCache = null;
    }

    public void RearrangeAfterRenumbering(List<int> newToOldElementMap, List<int> oldToNewNodeMap)
    {
        ArgumentNullException.ThrowIfNull(newToOldElementMap);
        ArgumentNullException.ThrowIfNull(oldToNewNodeMap);

        CompressElements(newToOldElementMap);
        PermuteNodes(oldToNewNodeMap);
    }

    #endregion

    #region Set Operations (Operators)

    /// <summary>
    ///     Matrix multiplication: computes A × B.
    /// </summary>
    /// <summary>
    ///     Matrix multiplication with adaptive safe/unsafe strategy.
    ///     Uses fast path without bounds checking when all nodes in left are valid indices.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static O2M operator *(O2M left, O2M right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        var leftMaxNode = left.GetMaxNode();

        // Validate dimension compatibility: left's node indices must be valid row indices in right.
        // leftMaxNode is the largest column index in left; it must be < right.Count (the number of rows).
        if (leftMaxNode >= right.Count && left.Count > 0 && right.Count > 0)
            throw new InvalidOperationException(
                $"Matrix dimension mismatch: left references node index {leftMaxNode}, " +
                $"but right only has {right.Count} rows (valid range [0, {right.Count - 1}]).");

        // ADAPTIVE STRATEGY: If all nodes in left are < right.Count, we can skip bounds checks
        var product = leftMaxNode < right.Count
            ? PerformSymbolicMultiplicationUnsafe(left._adjacencies, right._adjacencies)
            : PerformSymbolicMultiplicationCpp(left._adjacencies, right._adjacencies);

        return new O2M(product)
        {
            ParallelizationThreshold = Math.Max(left.ParallelizationThreshold, right.ParallelizationThreshold)
        };
    }

    /// <summary>
    ///     Union operator: combines adjacency lists element-wise.
    ///     Parallel with thread-local markers.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static O2M operator |(O2M left, O2M right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        var maxElements = Math.Max(left.Count, right.Count);
        var result = new O2M(maxElements)
        {
            ParallelizationThreshold = Math.Max(left.ParallelizationThreshold, right.ParallelizationThreshold)
        };

        var leftMax = left.GetMaxNode();
        var rightMax = right.GetMaxNode();
        var maxNode = Math.Max(leftMax, rightMax);

        if (maxNode < 0)
        {
            for (var i = 0; i < maxElements; i++)
                result._adjacencies.Add(new List<int>());
            return result;
        }

        var resultRows = new List<int>[maxElements];

        // P1.3 FIX: Use consistent threshold from both operands
        var threshold = Math.Max(left.ParallelizationThreshold, right.ParallelizationThreshold);

        if (maxElements >= threshold)
        {
            // P0.4 FIX: Use touched-index clearing instead of full Array.Clear.
            // P2.1 FIX: Removed (element << 1) + 1 generation which could overflow for large element counts.
            // Instead, use a simple marker approach with per-element touched-index clearing.
            var markerSize = maxNode + 1;

            Parallel.For(0, maxElements, ParallelConfig.Options,
                () =>
                {
                    var marker = ArrayPool<int>.Shared.Rent(markerSize);
                    Array.Clear(marker, 0, markerSize); // Initial clear only once per worker
                    var touched = ArrayPool<int>.Shared.Rent(markerSize);
                    return (marker, touched);
                },
                (element, _, state) =>
                {
                    var marker = state.marker;
                    var touched = state.touched;
                    var count = 0;

                    // First pass: count unique nodes using marker=1 for "seen"
                    if (element < left.Count)
                    {
                        var elementA = left._adjacencies[element];
                        foreach (var node in elementA)
                            if ((uint)node <= (uint)maxNode && marker[node] == 0)
                            {
                                marker[node] = 1;
                                touched[count] = node;
                                count++;
                            }
                    }

                    if (element < right.Count)
                    {
                        var elementB = right._adjacencies[element];
                        foreach (var node in elementB)
                            if ((uint)node <= (uint)maxNode && marker[node] == 0)
                            {
                                marker[node] = 1;
                                touched[count] = node;
                                count++;
                            }
                    }

                    // Build result list
                    var temp = new List<int>(count);
                    if (count > 0)
                    {
                        CollectionsMarshal.SetCount(temp, count);
                        var tempSpan = CollectionsMarshal.AsSpan(temp);
                        for (var i = 0; i < count; i++)
                            tempSpan[i] = touched[i];
                    }

                    resultRows[element] = temp;

                    // P0.4 FIX: Clear only touched indices
                    for (var i = 0; i < count; i++)
                        marker[touched[i]] = 0;

                    return state;
                },
                state =>
                {
                    ArrayPool<int>.Shared.Return(state.marker);
                    ArrayPool<int>.Shared.Return(state.touched);
                });
        }
        else
        {
            // Sequential with single marker — P0.4 FIX: touched-index clearing, single-pass
            var marker = new int[maxNode + 1];
            var touched = new int[maxNode + 1];

            for (var element = 0; element < maxElements; element++)
            {
                var count = 0;

                if (element < left.Count)
                    foreach (var node in left._adjacencies[element])
                        if ((uint)node <= (uint)maxNode && marker[node] == 0)
                        {
                            marker[node] = 1;
                            touched[count] = node;
                            count++;
                        }

                if (element < right.Count)
                    foreach (var node in right._adjacencies[element])
                        if ((uint)node <= (uint)maxNode && marker[node] == 0)
                        {
                            marker[node] = 1;
                            touched[count] = node;
                            count++;
                        }

                var temp = new List<int>(count);
                if (count > 0)
                {
                    CollectionsMarshal.SetCount(temp, count);
                    var tempSpan = CollectionsMarshal.AsSpan(temp);
                    for (var i = 0; i < count; i++)
                        tempSpan[i] = touched[i];
                }

                resultRows[element] = temp;

                // P0.4 FIX: Clear only touched indices
                for (var i = 0; i < count; i++)
                    marker[touched[i]] = 0;
            }
        }

        result._adjacencies.AddRange(resultRows);
        return result;
    }

    public static O2M operator +(O2M left, O2M right)
    {
        return left | right;
    }

    /// <summary>
    ///     Intersection operator.
    ///     Parallel per-element with vectorized merge.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static O2M operator &(O2M left, O2M right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        var elementCount = Math.Min(left.Count, right.Count);
        var resultRows = new List<int>[elementCount];

        var threshold = Math.Max(left.ParallelizationThreshold, right.ParallelizationThreshold);

        if (elementCount >= threshold)
            Parallel.For(0, elementCount, ParallelConfig.Options, element =>
            {
                var elementA = left._adjacencies[element];
                var elementB = right._adjacencies[element];

                // Adaptive strategy: HashSet for small lists, sort+merge for large
                const int hashSetThreshold = 100;
                if (elementA.Count + elementB.Count < hashSetThreshold)
                {
                    resultRows[element] = GetIntersectionHashSet(elementA, elementB);
                }
                else
                {
                    var copyA = new List<int>(elementA);
                    var copyB = new List<int>(elementB);
                    copyA.Sort();
                    copyB.Sort();
                    resultRows[element] = GetIntersectionSorted(
                        CollectionsMarshal.AsSpan(copyA),
                        CollectionsMarshal.AsSpan(copyB));
                }
            });
        else
            for (var element = 0; element < elementCount; element++)
            {
                var elementA = left._adjacencies[element];
                var elementB = right._adjacencies[element];

                // Adaptive strategy: HashSet for small lists, sort+merge for large
                const int hashSetThreshold = 100;
                if (elementA.Count + elementB.Count < hashSetThreshold)
                {
                    resultRows[element] = GetIntersectionHashSet(elementA, elementB);
                }
                else
                {
                    var copyA = new List<int>(elementA);
                    var copyB = new List<int>(elementB);
                    copyA.Sort();
                    copyB.Sort();
                    resultRows[element] = GetIntersectionSorted(
                        CollectionsMarshal.AsSpan(copyA),
                        CollectionsMarshal.AsSpan(copyB));
                }
            }

        var result = new O2M(elementCount) { ParallelizationThreshold = threshold };
        result._adjacencies.AddRange(resultRows);
        return result;
    }

    /// <summary>
    ///     Difference operator.
    ///     Parallel per-element.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static O2M operator -(O2M left, O2M right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        var elementCount = left.Count;
        var resultRows = new List<int>[elementCount];
        var threshold = Math.Max(left.ParallelizationThreshold, right.ParallelizationThreshold);

        if (elementCount >= threshold)
            Parallel.For(0, elementCount, ParallelConfig.Options, element =>
            {
                if (element < right.Count)
                {
                    var elementA = new List<int>(left._adjacencies[element]);
                    var elementB = new List<int>(right._adjacencies[element]);

                    elementA.Sort();
                    elementB.Sort();

                    resultRows[element] = GetDifferenceSorted(
                        CollectionsMarshal.AsSpan(elementA),
                        CollectionsMarshal.AsSpan(elementB));
                }
                else
                {
                    var src = left._adjacencies[element];
                    var dst = new List<int>(src.Count);
                    if (src.Count > 0)
                    {
                        CollectionsMarshal.SetCount(dst, src.Count);
                        CollectionsMarshal.AsSpan(src).CopyTo(CollectionsMarshal.AsSpan(dst));
                    }

                    resultRows[element] = dst;
                }
            });
        else
            for (var element = 0; element < elementCount; element++)
                if (element < right.Count)
                {
                    var elementA = new List<int>(left._adjacencies[element]);
                    var elementB = new List<int>(right._adjacencies[element]);

                    elementA.Sort();
                    elementB.Sort();

                    resultRows[element] = GetDifferenceSorted(
                        CollectionsMarshal.AsSpan(elementA),
                        CollectionsMarshal.AsSpan(elementB));
                }
                else
                {
                    resultRows[element] = new List<int>(left._adjacencies[element]);
                }

        var result = new O2M(elementCount) { ParallelizationThreshold = threshold };
        result._adjacencies.AddRange(resultRows);
        return result;
    }

    /// <summary>
    ///     Symmetric difference operator: nodes that appear in exactly one operand per element.
    ///     Single-pass implementation with marker arrays.
    /// </summary>
    /// <remarks>
    ///     <b>P1 FIX:</b> Previous implementation <c>(left | right) - (left &amp; right)</c>
    ///     allocated 3 intermediate O2M instances and performed 3 full passes.
    ///     This single-pass version uses the same touched-index marker pattern as union,
    ///     with a two-state marker: 1 = seen once (include), 2 = seen in both (exclude).
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public static O2M operator ^(O2M left, O2M right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        var maxElements = Math.Max(left.Count, right.Count);
        var result = new O2M(maxElements)
        {
            ParallelizationThreshold = Math.Max(left.ParallelizationThreshold, right.ParallelizationThreshold)
        };

        var leftMax = left.GetMaxNode();
        var rightMax = right.GetMaxNode();
        var maxNode = Math.Max(leftMax, rightMax);

        if (maxNode < 0)
        {
            for (var i = 0; i < maxElements; i++)
                result._adjacencies.Add(new List<int>());
            return result;
        }

        var resultRows = new List<int>[maxElements];
        var threshold = Math.Max(left.ParallelizationThreshold, right.ParallelizationThreshold);

        if (maxElements >= threshold)
        {
            var markerSize = maxNode + 1;

            Parallel.For(0, maxElements, ParallelConfig.Options,
                () =>
                {
                    var marker = ArrayPool<int>.Shared.Rent(markerSize);
                    Array.Clear(marker, 0, markerSize);
                    var touched = ArrayPool<int>.Shared.Rent(markerSize);
                    return (marker, touched);
                },
                (element, _, state) =>
                {
                    var marker = state.marker;
                    var touched = state.touched;
                    var touchedCount = 0;

                    // Pass 1: mark left nodes as 1 (seen once)
                    if (element < left.Count)
                    {
                        var leftNodes = left._adjacencies[element];
                        foreach (var node in leftNodes)
                            if ((uint)node <= (uint)maxNode && marker[node] == 0)
                            {
                                marker[node] = 1;
                                touched[touchedCount++] = node;
                            }
                    }

                    // Pass 2: for right nodes, toggle: 1→2 (in both, exclude), 0→1 (right only, include)
                    if (element < right.Count)
                    {
                        var rightNodes = right._adjacencies[element];
                        foreach (var node in rightNodes)
                            if ((uint)node <= (uint)maxNode)
                            {
                                if (marker[node] == 1)
                                {
                                    marker[node] = 2; // In both → exclude
                                }
                                else if (marker[node] == 0)
                                {
                                    marker[node] = 1; // Right only → include
                                    touched[touchedCount++] = node;
                                }
                                // marker[node] == 2 means already excluded (duplicate in right), skip
                            }
                    }

                    // Collect nodes with marker == 1 (appeared in exactly one operand)
                    var count = 0;
                    for (var i = 0; i < touchedCount; i++)
                        if (marker[touched[i]] == 1)
                            count++;

                    var temp = new List<int>(count);
                    if (count > 0)
                    {
                        CollectionsMarshal.SetCount(temp, count);
                        var tempSpan = CollectionsMarshal.AsSpan(temp);
                        var idx = 0;
                        for (var i = 0; i < touchedCount; i++)
                            if (marker[touched[i]] == 1)
                                tempSpan[idx++] = touched[i];
                    }

                    resultRows[element] = temp;

                    // Clear only touched indices
                    for (var i = 0; i < touchedCount; i++)
                        marker[touched[i]] = 0;

                    return state;
                },
                state =>
                {
                    ArrayPool<int>.Shared.Return(state.marker);
                    ArrayPool<int>.Shared.Return(state.touched);
                });
        }
        else
        {
            // Sequential with single marker
            var marker = new int[maxNode + 1];
            var touched = new int[maxNode + 1];

            for (var element = 0; element < maxElements; element++)
            {
                var touchedCount = 0;

                if (element < left.Count)
                    foreach (var node in left._adjacencies[element])
                        if ((uint)node <= (uint)maxNode && marker[node] == 0)
                        {
                            marker[node] = 1;
                            touched[touchedCount++] = node;
                        }

                if (element < right.Count)
                    foreach (var node in right._adjacencies[element])
                        if ((uint)node <= (uint)maxNode)
                        {
                            if (marker[node] == 1)
                                marker[node] = 2;
                            else if (marker[node] == 0)
                            {
                                marker[node] = 1;
                                touched[touchedCount++] = node;
                            }
                        }

                var count = 0;
                for (var i = 0; i < touchedCount; i++)
                    if (marker[touched[i]] == 1)
                        count++;

                var temp = new List<int>(count);
                if (count > 0)
                {
                    CollectionsMarshal.SetCount(temp, count);
                    var tempSpan = CollectionsMarshal.AsSpan(temp);
                    var idx = 0;
                    for (var i = 0; i < touchedCount; i++)
                        if (marker[touched[i]] == 1)
                            tempSpan[idx++] = touched[i];
                }

                resultRows[element] = temp;

                for (var i = 0; i < touchedCount; i++)
                    marker[touched[i]] = 0;
            }
        }

        result._adjacencies.AddRange(resultRows);
        return result;
    }

    #endregion

    #region Matrix Multiplication Implementation

    /// <summary>
    ///     Fast multiplication path without bounds checking.
    ///     PRECONDITION: All nodes in leftAdj must be valid indices into rightAdj (leftMaxNode &lt; rightCount).
    ///     PERFORMANCE: 10-30% faster than safe version due to eliminated bounds checks.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    private static List<List<int>> PerformSymbolicMultiplicationUnsafe(
        List<List<int>> leftAdj,
        List<List<int>> rightAdj)
    {
        var leftCount = leftAdj.Count;
        var rightCount = rightAdj.Count;

        if (rightCount == 0)
        {
            var empty = new List<List<int>>(leftCount);
            for (var i = 0; i < leftCount; i++)
                empty.Add(new List<int>());
            return empty;
        }

        // Find max node in right adjacencies
        var maxNode = -1;
        for (var i = 0; i < rightCount; i++)
        {
            var span = CollectionsMarshal.AsSpan(rightAdj[i]);
            foreach (var node in span)
                if (node > maxNode)
                    maxNode = node;
        }

        if (maxNode < 0)
        {
            var empty = new List<List<int>>(leftCount);
            for (var i = 0; i < leftCount; i++)
                empty.Add(new List<int>());
            return empty;
        }

        var result = new List<List<int>>(leftCount);
        var rowSizes = new int[leftCount];

        var markerSize = maxNode + 1;

        // P0.4 FIX: Use touched-index clearing instead of full Array.Clear.
        // For FEM-typical meshes (element arity ~3–20), each row touches a small fraction
        // of the marker array. Clearing only touched indices is proportional to local
        // neighborhood size, not global nodeCount.

        // P1.3 FIX: Apply parallelization threshold consistently
        var threshold = DefaultParallelizationThreshold;

        // PASS 1: Calculate row sizes
        if (leftCount >= threshold)
        {
            Parallel.For(0, leftCount,
                ParallelConfig.Options,
                () =>
                {
                    var marker = ArrayPool<int>.Shared.Rent(markerSize);
                    Array.Clear(marker, 0, markerSize); // Initial clear only once per worker
                    var touched = ArrayPool<int>.Shared.Rent(markerSize);
                    return (marker, touched);
                },
                (ra, loopState, state) =>
                {
                    var marker = state.marker;
                    var touched = state.touched;
                    var generation = ra + 1;
                    var len = 0;

                    var raNodes = leftAdj[ra];
                    foreach (var ca in raNodes)
                    {
                        var rbNodes = rightAdj[ca];
                        foreach (var cb in rbNodes)
                            if (marker[cb] != generation)
                            {
                                marker[cb] = generation;
                                touched[len] = cb;
                                len++;
                            }
                    }

                    rowSizes[ra] = len;
                    // P0.4 FIX: Clear only touched indices instead of entire marker
                    for (var i = 0; i < len; i++)
                        marker[touched[i]] = 0;
                    return (marker, touched);
                },
                state =>
                {
                    ArrayPool<int>.Shared.Return(state.marker);
                    ArrayPool<int>.Shared.Return(state.touched);
                });
        }
        else
        {
            // Sequential path
            var marker = new int[markerSize];
            var touched = new int[markerSize];
            for (var ra = 0; ra < leftCount; ra++)
            {
                var generation = ra + 1;
                var len = 0;

                var raNodes = leftAdj[ra];
                foreach (var ca in raNodes)
                {
                    var rbNodes = rightAdj[ca];
                    foreach (var cb in rbNodes)
                        if (marker[cb] != generation)
                        {
                            marker[cb] = generation;
                            touched[len] = cb;
                            len++;
                        }
                }

                rowSizes[ra] = len;
                for (var i = 0; i < len; i++)
                    marker[touched[i]] = 0;
            }
        }

        // Allocate exact sizes
        for (var ra = 0; ra < leftCount; ra++)
        {
            var list = new List<int>(rowSizes[ra]);
            if (rowSizes[ra] > 0)
                CollectionsMarshal.SetCount(list, rowSizes[ra]);
            result.Add(list);
        }

        // PASS 2: Fill values
        if (leftCount >= threshold)
        {
            Parallel.For(0, leftCount,
                ParallelConfig.Options,
                () =>
                {
                    var marker = ArrayPool<int>.Shared.Rent(markerSize);
                    Array.Clear(marker, 0, markerSize); // Initial clear only once per worker
                    var touched = ArrayPool<int>.Shared.Rent(markerSize);
                    return (marker, touched);
                },
                (ra, loopState, state) =>
                {
                    var marker = state.marker;
                    var touched = state.touched;
                    var generation = ra + 1;
                    var len = 0;

                    var raNodes = leftAdj[ra];
                    var rcSpan = CollectionsMarshal.AsSpan(result[ra]);

                    foreach (var ca in raNodes)
                    {
                        var rbNodes = rightAdj[ca];
                        foreach (var cb in rbNodes)
                            if (marker[cb] != generation)
                            {
                                marker[cb] = generation;
                                touched[len] = cb;
                                rcSpan[len++] = cb;
                            }
                    }

                    // P0.4 FIX: Clear only touched indices
                    for (var i = 0; i < len; i++)
                        marker[touched[i]] = 0;
                    return (marker, touched);
                },
                state =>
                {
                    ArrayPool<int>.Shared.Return(state.marker);
                    ArrayPool<int>.Shared.Return(state.touched);
                });
        }
        else
        {
            // Sequential path
            var marker = new int[markerSize];
            var touched = new int[markerSize];
            for (var ra = 0; ra < leftCount; ra++)
            {
                var generation = ra + 1;
                var len = 0;

                var raNodes = leftAdj[ra];
                var rcSpan = CollectionsMarshal.AsSpan(result[ra]);

                foreach (var ca in raNodes)
                {
                    var rbNodes = rightAdj[ca];
                    foreach (var cb in rbNodes)
                        if (marker[cb] != generation)
                        {
                            marker[cb] = generation;
                            touched[len] = cb;
                            rcSpan[len++] = cb;
                        }
                }

                for (var i = 0; i < len; i++)
                    marker[touched[i]] = 0;
            }
        }

        return result;
    }

    /// <summary>
    ///     Safe multiplication path with bounds checking.
    ///     Used when leftMaxNode >= rightCount (nodes may be out of range).
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    private static List<List<int>> PerformSymbolicMultiplicationCpp(
        List<List<int>> leftAdj,
        List<List<int>> rightAdj)
    {
        var leftCount = leftAdj.Count;
        var rightCount = rightAdj.Count;

        if (rightCount == 0)
        {
            var empty = new List<List<int>>(leftCount);
            for (var i = 0; i < leftCount; i++)
                empty.Add(new List<int>());
            return empty;
        }

        // Find max node in right adjacencies
        var maxNode = -1;
        for (var i = 0; i < rightCount; i++)
        {
            var span = CollectionsMarshal.AsSpan(rightAdj[i]);
            foreach (var node in span)
                if (node > maxNode)
                    maxNode = node;
        }

        if (maxNode < 0)
        {
            var empty = new List<List<int>>(leftCount);
            for (var i = 0; i < leftCount; i++)
                empty.Add(new List<int>());
            return empty;
        }

        var result = new List<List<int>>(leftCount);
        var rowSizes = new int[leftCount];
        var camax = rightCount - 1;

        var markerSize = maxNode + 1;

        // P1.3 FIX: Apply parallelization threshold consistently
        var threshold = DefaultParallelizationThreshold;

        // PASS 1: Calculate row sizes — P0.4 FIX: touched-index clearing
        if (leftCount >= threshold)
        {
            Parallel.For(0, leftCount,
                ParallelConfig.Options,
                () =>
                {
                    var marker = ArrayPool<int>.Shared.Rent(markerSize);
                    Array.Clear(marker, 0, markerSize);
                    var touched = ArrayPool<int>.Shared.Rent(markerSize);
                    return (marker, touched);
                },
                (ra, loopState, state) =>
                {
                    var marker = state.marker;
                    var touched = state.touched;
                    var generation = ra + 1;
                    var len = 0;

                    var raNodes = leftAdj[ra];
                    foreach (var ca in raNodes)
                    {
                        if ((uint)ca > (uint)camax) continue;

                        var rbNodes = rightAdj[ca];
                        foreach (var cb in rbNodes)
                            if (marker[cb] != generation)
                            {
                                marker[cb] = generation;
                                touched[len] = cb;
                                len++;
                            }
                    }

                    rowSizes[ra] = len;
                    for (var i = 0; i < len; i++)
                        marker[touched[i]] = 0;
                    return state;
                },
                state =>
                {
                    ArrayPool<int>.Shared.Return(state.marker);
                    ArrayPool<int>.Shared.Return(state.touched);
                });
        }
        else
        {
            var marker = new int[markerSize];
            var touched = new int[markerSize];
            for (var ra = 0; ra < leftCount; ra++)
            {
                var generation = ra + 1;
                var len = 0;

                var raNodes = leftAdj[ra];
                foreach (var ca in raNodes)
                {
                    if ((uint)ca > (uint)camax) continue;
                    var rbNodes = rightAdj[ca];
                    foreach (var cb in rbNodes)
                        if (marker[cb] != generation)
                        {
                            marker[cb] = generation;
                            touched[len] = cb;
                            len++;
                        }
                }

                rowSizes[ra] = len;
                for (var i = 0; i < len; i++)
                    marker[touched[i]] = 0;
            }
        }

        // Allocate exact sizes
        for (var ra = 0; ra < leftCount; ra++)
        {
            var list = new List<int>(rowSizes[ra]);
            if (rowSizes[ra] > 0)
                CollectionsMarshal.SetCount(list, rowSizes[ra]);
            result.Add(list);
        }

        // PASS 2: Fill values — P0.4 FIX: touched-index clearing
        if (leftCount >= threshold)
        {
            Parallel.For(0, leftCount,
                ParallelConfig.Options,
                () =>
                {
                    var marker = ArrayPool<int>.Shared.Rent(markerSize);
                    Array.Clear(marker, 0, markerSize);
                    var touched = ArrayPool<int>.Shared.Rent(markerSize);
                    return (marker, touched);
                },
                (ra, loopState, state) =>
                {
                    var marker = state.marker;
                    var touched = state.touched;
                    var generation = ra + 1;
                    var len = 0;

                    var raNodes = leftAdj[ra];
                    var rcSpan = CollectionsMarshal.AsSpan(result[ra]);

                    foreach (var ca in raNodes)
                    {
                        if ((uint)ca > (uint)camax) continue;
                        var rbNodes = rightAdj[ca];
                        foreach (var cb in rbNodes)
                            if (marker[cb] != generation)
                            {
                                marker[cb] = generation;
                                touched[len] = cb;
                                rcSpan[len++] = cb;
                            }
                    }

                    for (var i = 0; i < len; i++)
                        marker[touched[i]] = 0;
                    return state;
                },
                state =>
                {
                    ArrayPool<int>.Shared.Return(state.marker);
                    ArrayPool<int>.Shared.Return(state.touched);
                });
        }
        else
        {
            var marker = new int[markerSize];
            var touched = new int[markerSize];
            for (var ra = 0; ra < leftCount; ra++)
            {
                var generation = ra + 1;
                var len = 0;

                var raNodes = leftAdj[ra];
                var rcSpan = CollectionsMarshal.AsSpan(result[ra]);

                foreach (var ca in raNodes)
                {
                    if ((uint)ca > (uint)camax) continue;
                    var rbNodes = rightAdj[ca];
                    foreach (var cb in rbNodes)
                        if (marker[cb] != generation)
                        {
                            marker[cb] = generation;
                            touched[len] = cb;
                            rcSpan[len++] = cb;
                        }
                }

                for (var i = 0; i < len; i++)
                    marker[touched[i]] = 0;
            }
        }

        return result;
    }

    /// <summary>
    ///     Gets intersection using HashSet strategy (optimal for small lists).
    ///     Builds HashSet from smaller list for better performance.
    /// </summary>
    [MethodImpl(AggressiveOptimization | AggressiveInlining)]
    private static List<int> GetIntersectionHashSet(List<int> leftRow, List<int> rightRow)
    {
        // Build HashSet from rightRow; iterate leftRow to preserve its ordering in the output
        var hashSet = new HashSet<int>(rightRow);
        var result = new List<int>();

        foreach (var node in leftRow)
            if (hashSet.Contains(node))
                result.Add(node);

        return result;
    }

    /// <summary>
    ///     Gets intersection of two sorted spans.
    ///     Uses galloping search for skewed distributions.
    /// </summary>
    [MethodImpl(AggressiveOptimization | AggressiveInlining)]
    private static List<int> GetIntersectionSorted(ReadOnlySpan<int> sortedA, ReadOnlySpan<int> sortedB)
    {
        if (sortedA.Length == 0 || sortedB.Length == 0)
            return new List<int>();

        // Estimate result size
        var result = new List<int>(Math.Min(sortedA.Length, sortedB.Length));

        int i = 0, j = 0;

        // Use galloping for very skewed sizes
        if (sortedA.Length > sortedB.Length * 8)
            // Gallop through A
            while (i < sortedA.Length && j < sortedB.Length)
            {
                var target = sortedB[j];

                // Galloping search in A
                var step = 1;
                while (i + step < sortedA.Length && sortedA[i + step] < target)
                    step *= 2;

                // Binary search in [i, i+step]
                int lo = i, hi = Math.Min(i + step, sortedA.Length - 1);
                while (lo < hi)
                {
                    var mid = lo + (hi - lo) / 2;
                    if (sortedA[mid] < target)
                        lo = mid + 1;
                    else
                        hi = mid;
                }

                i = lo;

                if (i < sortedA.Length && sortedA[i] == target)
                {
                    result.Add(target);
                    i++;
                }

                j++;
            }
        else if (sortedB.Length > sortedA.Length * 8)
            // Gallop through B (symmetric case)
            while (i < sortedA.Length && j < sortedB.Length)
            {
                var target = sortedA[i];

                var step = 1;
                while (j + step < sortedB.Length && sortedB[j + step] < target)
                    step *= 2;

                int lo = j, hi = Math.Min(j + step, sortedB.Length - 1);
                while (lo < hi)
                {
                    var mid = lo + (hi - lo) / 2;
                    if (sortedB[mid] < target)
                        lo = mid + 1;
                    else
                        hi = mid;
                }

                j = lo;

                if (j < sortedB.Length && sortedB[j] == target)
                {
                    result.Add(target);
                    j++;
                }

                i++;
            }
        else
            // Standard merge
            while (i < sortedA.Length && j < sortedB.Length)
                if (sortedA[i] < sortedB[j])
                {
                    i++;
                }
                else if (sortedA[i] > sortedB[j])
                {
                    j++;
                }
                else
                {
                    result.Add(sortedA[i]);
                    i++;
                    j++;
                }

        return result;
    }

    /// <summary>
    ///     Gets set difference of two sorted spans (A - B).
    ///     Uses galloping for skewed distributions.
    /// </summary>
    [MethodImpl(AggressiveOptimization | AggressiveInlining)]
    private static List<int> GetDifferenceSorted(ReadOnlySpan<int> sortedA, ReadOnlySpan<int> sortedB)
    {
        if (sortedA.Length == 0)
            return new List<int>();

        if (sortedB.Length == 0)
        {
            var result = new List<int>(sortedA.Length);
            CollectionsMarshal.SetCount(result, sortedA.Length);
            sortedA.CopyTo(CollectionsMarshal.AsSpan(result));
            return result;
        }

        var diff = new List<int>(sortedA.Length);
        int i = 0, j = 0;

        while (i < sortedA.Length && j < sortedB.Length)
            if (sortedA[i] < sortedB[j])
            {
                diff.Add(sortedA[i]);
                i++;
            }
            else if (sortedA[i] > sortedB[j])
            {
                j++;
            }
            else
            {
                i++;
                j++;
            }

        // Add remaining from A
        while (i < sortedA.Length)
        {
            diff.Add(sortedA[i]);
            i++;
        }

        return diff;
    }

    #endregion

    #region Comparison and Equality

    public int CompareTo(O2M? other)
    {
        if (other is null) return 1;
        if (Count != other.Count) return Count.CompareTo(other.Count);

        for (var i = 0; i < Count; i++)
        {
            var comparison = CompareRows(_adjacencies[i], other._adjacencies[i]);
            if (comparison != 0) return comparison;
        }

        return 0;
    }

    /// <summary>
    ///     Determines equality with another O2M.
    ///     Parallel comparison for large structures.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public bool Equals(O2M? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;

        var a = _adjacencies;
        var b = other._adjacencies;

        if (a.Count != b.Count) return false;

        var count = a.Count;

        if (count >= ParallelizationThreshold)
        {
            var isEqual = 1;

            Parallel.For(0, count, ParallelConfig.Options, (i, loopState) =>
            {
                if (Volatile.Read(ref isEqual) == 0)
                {
                    loopState.Stop();
                    return;
                }

                var ar = a[i];
                var br = b[i];

                if (ar.Count != br.Count ||
                    !CollectionsMarshal.AsSpan(ar).SequenceEqual(CollectionsMarshal.AsSpan(br)))
                {
                    Interlocked.Exchange(ref isEqual, 0);
                    loopState.Stop();
                }
            });

            return isEqual == 1;
        }

        for (var i = 0; i < count; i++)
        {
            var ar = a[i];
            var br = b[i];

            if (ar.Count != br.Count) return false;
            if (!CollectionsMarshal.AsSpan(ar).SequenceEqual(CollectionsMarshal.AsSpan(br))) return false;
        }

        return true;
    }

    public bool IsPermutationOf(O2M? other)
    {
        if (other is null) return false;
        if (Count != other.Count) return false;
        if (ReferenceEquals(this, other)) return true;
        if (Count == 0) return true;

        var orderA = GetSortOrder();
        var orderB = other.GetSortOrder();

        for (var i = 0; i < Count; i++)
            if (CompareRows(_adjacencies[orderA[i]], other._adjacencies[orderB[i]]) != 0)
                return false;

        return true;
    }

    public override bool Equals(object? obj)
    {
        return Equals(obj as O2M);
    }

    /// <summary>
    ///     Computes a hash code for this structure.
    ///     <b>PERFORMANCE TRADE-OFF:</b> This implementation uses a sampling strategy rather
    ///     than hashing all elements. It samples 5 strategic rows (first, last, middle, quarter,
    ///     three-quarter) plus the maximum node value. This is much faster than hashing every
    ///     element but has a higher chance of hash collision.
    ///     <b>RATIONALE:</b> For large matrices with thousands of rows, computing a full hash
    ///     would be prohibitively expensive (O(total edges)). The sampling approach provides
    ///     good distribution for typical use cases (e.g., as dictionary keys) while maintaining
    ///     O(1) expected time complexity.
    ///     <b>WARNING:</b> If using O2M as dictionary key with many topologically-similar
    ///     instances (same row count, similar boundary rows), consider using
    ///     <see cref="FullContentHashCode" /> or reference equality via
    ///     <c>RuntimeHelpers.GetHashCode(obj)</c> to avoid pathological collision chains.
    /// </summary>
    public override int GetHashCode()
    {
        var hashCode = new HashCode();
        var count = _adjacencies.Count;

        hashCode.Add(count);

        if (count > 0)
        {
            // Hash first row
            foreach (var item in _adjacencies[0])
                hashCode.Add(item);

            // Hash last row if different
            if (count > 1)
                foreach (var item in _adjacencies[count - 1])
                    hashCode.Add(item);

            // Hash middle row
            if (count > 2)
                foreach (var item in _adjacencies[count / 2])
                    hashCode.Add(item);

            // Hash quarter points for better distribution
            if (count > 8)
            {
                foreach (var item in _adjacencies[count / 4])
                    hashCode.Add(item);
                foreach (var item in _adjacencies[count * 3 / 4])
                    hashCode.Add(item);
            }
        }

        hashCode.Add(GetMaxNode());
        return hashCode.ToHashCode();
    }

    /// <summary>
    ///     Returns a hash code computed from ALL content (every node in every row).
    ///     Use this instead of <see cref="GetHashCode" /> when storing many topologically-similar
    ///     O2M instances in hash-based collections and collision avoidance is critical.
    /// </summary>
    /// <remarks>
    ///     <b>COMPLEXITY:</b> O(total edges) - significantly slower than <see cref="GetHashCode" />.
    ///     <b>USE CASE:</b> Dictionary keys where you have many meshes with identical row counts
    ///     and similar structure (e.g., uniform refinement levels, regular grids).
    /// </remarks>
    public int FullContentHashCode()
    {
        var hashCode = new HashCode();
        hashCode.Add(_adjacencies.Count);

        foreach (var row in _adjacencies)
        {
            hashCode.Add(row.Count);
            foreach (var node in row)
                hashCode.Add(node);
        }

        return hashCode.ToHashCode();
    }

    public static bool operator ==(O2M? left, O2M? right)
    {
        return ReferenceEquals(left, right) || left?.Equals(right) == true;
    }

    public static bool operator !=(O2M? left, O2M? right)
    {
        return !(left == right);
    }

    public static bool operator <(O2M? left, O2M? right)
    {
        return left is null ? right is not null : left.CompareTo(right) < 0;
    }

    public static bool operator >(O2M? left, O2M? right)
    {
        return right < left;
    }

    public static bool operator <=(O2M? left, O2M? right)
    {
        return !(left > right);
    }

    public static bool operator >=(O2M? left, O2M? right)
    {
        return !(left < right);
    }

    #endregion

    #region Conversion and Serialization

    /// <summary>
    ///     Converts to Compressed Sparse Row (CSR) format.
    ///     Parallel counting and copying.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public (int[] rowPtr, int[] columnIndices) ToCsr()
    {
        var m = _adjacencies.Count;
        var rowPtr = new int[m + 1];

        if (m == 0)
            return (rowPtr, Array.Empty<int>());

        // Parallel count
        if (m >= ParallelizationThreshold)
            Parallel.For(0, m, ParallelConfig.Options, i => rowPtr[i + 1] = _adjacencies[i].Count);
        else
            for (var i = 0; i < m; i++)
                rowPtr[i + 1] = _adjacencies[i].Count;

        // Sequential prefix sum
        for (var i = 1; i <= m; i++)
            rowPtr[i] += rowPtr[i - 1];

        var col = new int[rowPtr[m]];

        // Parallel copy
        if (m >= ParallelizationThreshold)
            Parallel.For(0, m, ParallelConfig.Options,
                i => { CollectionsMarshal.AsSpan(_adjacencies[i]).CopyTo(col.AsSpan(rowPtr[i])); });
        else
            for (var i = 0; i < m; i++)
                CollectionsMarshal.AsSpan(_adjacencies[i]).CopyTo(col.AsSpan(rowPtr[i]));

        return (rowPtr, col);
    }

    /// <summary>
    ///     Creates an O2M from Compressed Sparse Row (CSR) format.
    ///     Parallel row construction.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static O2M FromCsr(int[] rowPointers, int[] columnIndices)
    {
        ArgumentNullException.ThrowIfNull(rowPointers);
        ArgumentNullException.ThrowIfNull(columnIndices);

        if (rowPointers.Length == 0) return new O2M();

        var m = rowPointers.Length - 1;
        var rows = new List<int>[m];

        if (m >= DefaultParallelizationThreshold)
            Parallel.For(0, m, ParallelConfig.Options, i =>
            {
                var start = rowPointers[i];
                var end = rowPointers[i + 1];
                var row = new List<int>(end - start);
                if (end > start)
                {
                    CollectionsMarshal.SetCount(row, end - start);
                    columnIndices.AsSpan(start, end - start).CopyTo(CollectionsMarshal.AsSpan(row));
                }

                rows[i] = row;
            });
        else
            for (var i = 0; i < m; i++)
            {
                var start = rowPointers[i];
                var end = rowPointers[i + 1];
                var row = new List<int>(end - start);
                for (var k = start; k < end; k++)
                    row.Add(columnIndices[k]);
                rows[i] = row;
            }

        return new O2M(new List<List<int>>(rows));
    }

    /// <summary>
    ///     Converts to a dense boolean matrix.
    ///     Parallel row processing.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public bool[,] ToBooleanMatrix()
    {
        var maxNodeValue = GetMaxNode();
        var rowCount = Count;
        var colCount = maxNodeValue + 1;

        if (colCount <= 0) return new bool[rowCount, 0];

        var matrix = new bool[rowCount, colCount];

        if (rowCount >= ParallelizationThreshold)
            Parallel.For(0, rowCount, ParallelConfig.Options, i =>
            {
                var span = CollectionsMarshal.AsSpan(_adjacencies[i]);
                foreach (var v in span)
                    if ((uint)v < (uint)colCount)
                        matrix[i, v] = true;
            });
        else
            for (var i = 0; i < rowCount; i++)
            {
                var span = CollectionsMarshal.AsSpan(_adjacencies[i]);
                foreach (var v in span)
                    if ((uint)v < (uint)colCount)
                        matrix[i, v] = true;
            }

        return matrix;
    }

    /// <summary>
    ///     Creates an O2M from a dense boolean matrix.
    ///     Parallel row construction.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static O2M FromBooleanMatrix(bool[,] matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);

        var r = matrix.GetLength(0);
        var c = matrix.GetLength(1);
        var rows = new List<int>[r];

        if (r >= DefaultParallelizationThreshold)
            Parallel.For(0, r, ParallelConfig.Options, i =>
            {
                // Count non-zeros
                var nnz = 0;
                for (var j = 0; j < c; j++)
                    if (matrix[i, j])
                        nnz++;

                // Fill row
                var row = new List<int>(nnz);
                for (var j = 0; j < c; j++)
                    if (matrix[i, j])
                        row.Add(j);

                rows[i] = row;
            });
        else
            for (var i = 0; i < r; i++)
            {
                var nnz = 0;
                for (var j = 0; j < c; j++)
                    if (matrix[i, j])
                        nnz++;

                var row = new List<int>(nnz);
                for (var j = 0; j < c; j++)
                    if (matrix[i, j])
                        row.Add(j);

                rows[i] = row;
            }

        return new O2M(new List<List<int>>(rows));
    }

    public override string ToString()
    {
        if (Count == 0) return string.Empty;

        var estimatedCapacity = Count * 20;
        var sb = new StringBuilder(estimatedCapacity);

        for (var i = 0; i < Count; i++)
        {
            sb.Append('[').Append(i).Append("] -> ");

            var row = _adjacencies[i];
            if (row.Count > 0)
                for (var j = 0; j < row.Count; j++)
                {
                    if (j > 0) sb.Append(", ");
                    sb.Append(row[j]);
                }

            sb.AppendLine();
        }

        return sb.ToString();
    }

    #endregion

    #region Implicit Conversions

    public static implicit operator O2M(List<List<int>> nodes)
    {
        return new O2M(nodes);
    }

    public static implicit operator O2M(List<int> elements)
    {
        ArgumentNullException.ThrowIfNull(elements);

        var nodes = new List<List<int>>(elements.Count);
        for (var i = 0; i < elements.Count; i++)
            nodes.Add([elements[i]]);

        return new O2M(nodes);
    }

    /// <summary>
    ///     Returns a deep copy of the internal adjacency lists.
    ///     Callers may freely mutate the returned structure without affecting this O2M.
    /// </summary>
    /// <remarks>
    ///     <b>P0.3 FIX:</b> The previous implicit conversion to <c>List&lt;List&lt;int&gt;&gt;</c>
    ///     exposed the raw internal list, allowing external mutation that silently invalidated
    ///     <c>_maxNodeIndexCache</c> and broke safety invariants established by the
    ///     <see cref="ReadOnlySpan{T}"/>-returning indexer.
    ///     Replaced with an explicit method that returns a defensive deep copy.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public List<List<int>> ToAdjacencyLists()
    {
        var count = _adjacencies.Count;
        var copy = new List<List<int>>(count);

        for (var i = 0; i < count; i++)
        {
            var src = _adjacencies[i];
            var dst = new List<int>(src.Count);
            if (src.Count > 0)
            {
                CollectionsMarshal.SetCount(dst, src.Count);
                CollectionsMarshal.AsSpan(src).CopyTo(CollectionsMarshal.AsSpan(dst));
            }

            copy.Add(dst);
        }

        return copy;
    }

    #endregion

    #region Node Position Analysis

    /// <summary>
    ///     Computes node positions within elements.
    ///     OPTIMIZED: Node-centric parallel iteration (no sorting required).
    /// </summary>
    /// <remarks>
    ///     <b>KEY OPTIMIZATION:</b> Instead of iterating elements and collecting positions
    ///     into ConcurrentBags that require sorting, we iterate nodes in parallel.
    ///     For each node, we look up its elements (from elementsFromNode, already sorted)
    ///     and find the node's position within each element via linear scan.
    ///     <b>COMPLEXITY:</b> O(nodeCount × avgElementsPerNode × avgNodesPerElement)
    ///     Typically O(N × 6 × 4) for tetrahedral meshes ≈ O(24N).
    ///     <b>WHY THIS IS FASTER:</b>
    ///     - No ConcurrentBag allocation (nodeCount bags eliminated)
    ///     - No LINQ OrderBy().ToList() sorting overhead
    ///     - No tuple allocations (elem, loc)
    ///     - Better cache locality: each thread works on contiguous node data
    ///     - Embarrassingly parallel: no synchronization between nodes
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    public static List<List<int>> GetNodePositions(O2M nodesFromElement, O2M elementsFromNode)
    {
        ArgumentNullException.ThrowIfNull(nodesFromElement);
        ArgumentNullException.ThrowIfNull(elementsFromNode);

        var nodeCount = elementsFromNode.Count;
        var elementCount = nodesFromElement.Count;

        if (nodeCount == 0 || elementCount == 0)
            return new List<List<int>>(0);

        // Pre-allocate result with exact sizes (using SetCount for zero-copy)
        var nodePos = new List<List<int>>(nodeCount);
        for (var n = 0; n < nodeCount; n++)
        {
            var elemCount = elementsFromNode._adjacencies[n].Count;
            var list = new List<int>(elemCount);
            CollectionsMarshal.SetCount(list, elemCount);
            nodePos.Add(list);
        }

        // Process nodes in parallel - each node is completely independent
        if (nodeCount >= nodesFromElement.ParallelizationThreshold)
            Parallel.For(0, nodeCount,
                ParallelConfig.Options,
                node =>
                {
                    // Get all elements containing this node (already sorted from Transpose)
                    var elements = CollectionsMarshal.AsSpan(elementsFromNode._adjacencies[node]);
                    var positions = CollectionsMarshal.AsSpan(nodePos[node]);

                    for (var le = 0; le < elements.Length; le++)
                    {
                        var element = elements[le];
                        if ((uint)element >= (uint)nodesFromElement.Count)
                        {
                            positions[le] = -1;
                            continue;
                        }
                        var elementNodes = CollectionsMarshal.AsSpan(nodesFromElement._adjacencies[element]);

                        // Linear search for node's position within element
                        // Element size is typically 3-8 nodes, so linear search is optimal
                        // (branch prediction + cache locality beats binary search here)
                        var pos = -1;
                        for (var k = 0; k < elementNodes.Length; k++)
                            if (elementNodes[k] == node)
                            {
                                pos = k;
                                break;
                            }

                        positions[le] = pos;
                    }
                });
        else
            // Sequential version
            for (var node = 0; node < nodeCount; node++)
            {
                var elements = CollectionsMarshal.AsSpan(elementsFromNode._adjacencies[node]);
                var positions = CollectionsMarshal.AsSpan(nodePos[node]);

                for (var le = 0; le < elements.Length; le++)
                {
                    var element = elements[le];
                    if ((uint)element >= (uint)nodesFromElement.Count)
                    {
                        positions[le] = -1;
                        continue;
                    }
                    var elementNodes = CollectionsMarshal.AsSpan(nodesFromElement._adjacencies[element]);

                    var pos = -1;
                    for (var k = 0; k < elementNodes.Length; k++)
                        if (elementNodes[k] == node)
                        {
                            pos = k;
                            break;
                        }

                    positions[le] = pos;
                }
            }

        return nodePos;
    }

    /// <summary>
    ///     Computes element positions assuming sorted node-to-element lists.
    ///     Parallel binary search.
    /// </summary>
    [MethodImpl(AggressiveOptimization)]
    public static List<List<int>> GetElementPositions(O2M nodesFromElement, O2M elementsFromNode)
    {
        ArgumentNullException.ThrowIfNull(nodesFromElement);
        ArgumentNullException.ThrowIfNull(elementsFromNode);

        var elemCount = nodesFromElement.Count;
        var elemPos = new List<List<int>>(elemCount);

        for (var e = 0; e < elemCount; e++)
        {
            var list = new List<int>(nodesFromElement._adjacencies[e].Count);
            CollectionsMarshal.SetCount(list, nodesFromElement._adjacencies[e].Count);
            elemPos.Add(list);
        }

        void Process(int e)
        {
            var nodes = CollectionsMarshal.AsSpan(nodesFromElement._adjacencies[e]);
            var positions = CollectionsMarshal.AsSpan(elemPos[e]);

            for (var j = 0; j < nodes.Length; j++)
            {
                var n = nodes[j];
                if ((uint)n >= (uint)elementsFromNode.Count)
                {
                    positions[j] = -1;
                    continue;
                }

                var pos = elementsFromNode._adjacencies[n].BinarySearch(e);
                if (pos < 0)
                    throw new InvalidOperationException(
                        $"Graph inconsistency: Element {e} → Node {n}, but reverse mapping missing.");

                positions[j] = pos;
            }
        }

        if (elemCount >= nodesFromElement.ParallelizationThreshold)
            Parallel.For(0, elemCount, ParallelConfig.Options, Process);
        else
            for (var e = 0; e < elemCount; e++)
                Process(e);

        return elemPos;
    }

    #endregion

    #region EPS Visualization

    public string ToEpsString()
    {
        const int Margin = 40;
        const int ElemSpacing = 20;
        const int NodeSpacing = 20;
        const int ElemRadius = 4;
        const int NodeRadius = 4;
        const double LineWidth = 0.5;
        const string ElementColor = "0 0 0";
        const string NodeColor = "0 0 0";
        const string LineColor = "0.5 0.5 0.5";
        const string TextColor = "0 0 0";
        const int FontSize = 12;
        const bool DrawElementLabels = true;
        const bool DrawNodeLabels = true;

        var elementCount = Count;
        var maxNodeValue = GetMaxNode();

        var elementsAreaWidth = elementCount > 0 ? 2 * ElemRadius : 0;
        var nodesAreaWidth = maxNodeValue >= 0
            ? NodeSpacing + maxNodeValue * NodeSpacing + 2 * NodeRadius
            : 0;
        var contentWidth = elementsAreaWidth + (nodesAreaWidth > 0 ? nodesAreaWidth : 0);
        var elementsAreaHeight = elementCount > 0
            ? 2 * ElemRadius + (elementCount - 1) * ElemSpacing
            : 0;

        var finalWidth = 2 * Margin + contentWidth;
        var finalHeight = 2 * Margin + elementsAreaHeight;

        var sb = new StringBuilder();
        sb.AppendLine("%!PS-Adobe-3.0 EPSF-3.0");
        sb.AppendLine($"%%BoundingBox: 0 0 {finalWidth} {finalHeight}");
        sb.AppendLine($"%%Title: {EscapePS("Sparse relation")}");
        sb.AppendLine("%%Creator: O2M.ToEpsString");
        sb.AppendLine("%%EndComments\n");
        sb.AppendLine($"/Times-Roman findfont {FontSize} scalefont setfont\n");

        // Draw elements
        sb.AppendLine($"{ElementColor} setrgbcolor");
        for (var i = 0; i < elementCount; i++)
        {
            double x = Margin + ElemRadius;
            double y = finalHeight - Margin - ElemRadius - i * ElemSpacing;
            sb.AppendLine($"{x} {y} {ElemRadius} 0 360 arc fill");

            if (DrawElementLabels)
                sb.AppendLine(
                    $"{TextColor} setrgbcolor " +
                    $"{x + ElemRadius + FontSize / 3.0} {y - FontSize / 3.0} " +
                    $"moveto ({i}) show");
        }

        // Draw nodes
        if (maxNodeValue >= 0)
        {
            sb.AppendLine($"{NodeColor} setrgbcolor");
            for (var j = 0; j <= maxNodeValue; j++)
            {
                double x = Margin + elementsAreaWidth + NodeSpacing + j * NodeSpacing;
                double y = Margin + NodeRadius;
                sb.AppendLine($"{x} {y} {NodeRadius} 0 360 arc fill");

                if (DrawNodeLabels)
                    sb.AppendLine(
                        $"{TextColor} setrgbcolor " +
                        $"{x} {y + NodeRadius + FontSize * 0.5} " +
                        $"moveto ({j}) dup stringwidth pop 2 div neg 0 rmoveto show");
            }
        }

        // Draw edges
        if (elementCount > 0 && maxNodeValue >= 0)
        {
            sb.AppendLine($"{LineWidth} setlinewidth");
            sb.AppendLine($"{LineColor} setrgbcolor");

            for (var i = 0; i < elementCount; i++)
            {
                double elemX = Margin + ElemRadius;
                double elemY = finalHeight - Margin - ElemRadius - i * ElemSpacing;
                var startX = elemX + ElemRadius;
                var span = CollectionsMarshal.AsSpan(_adjacencies[i]);

                foreach (var node in span)
                {
                    if (node > maxNodeValue) continue;

                    double nodeX = Margin + elementsAreaWidth + NodeSpacing + node * NodeSpacing;
                    double nodeY = Margin + NodeRadius;
                    sb.AppendLine($"{startX} {elemY} moveto {nodeX} {nodeY} lineto stroke");
                }
            }
        }

        sb.AppendLine("\nshowpage\n%%EOF");
        return sb.ToString();
    }

    private static string EscapePS(string text)
    {
        if (string.IsNullOrEmpty(text)) return string.Empty;
        return text.Replace("\\", "\\\\").Replace("(", "\\(").Replace(")", "\\)");
    }

    #endregion

    #region Clique Computation

    /// <summary>
    ///     Computes clique indices for finite element matrix assembly.
    ///     OPTIMIZED: Uses optimized GetNodePositions and reduces allocation overhead.
    /// </summary>
    /// <param name="nodesFromElement">The element-to-node adjacency (forward mapping).</param>
    /// <param name="elementsFromNode">The node-to-element adjacency (transpose).</param>
    /// <returns>
    ///     A list where result[element] is a flattened ns×ns array.
    ///     For local nodes i,j: result[element][j + i*ns] gives the clique index.
    /// </returns>
    /// <remarks>
    ///     <b>OPTIMIZATIONS:</b>
    ///     1. Uses optimized GetNodePositions (no sorting overhead)
    ///     2. Avoids ArrayPool for small-medium meshes (direct allocation is cheaper)
    ///     3. Thread-local arrays are reused across iterations
    ///     <b>PURPOSE:</b> Essential for sparse matrix assembly in FEM.
    ///     Identifies which node pairs share the same "clique" across elements.
    ///     <b>ALGORITHM:</b> Uses a generation-based marker approach for O(1) lookup.
    ///     Each node1 processes all its elements, marking node2's with unique indices.
    ///     When the same node2 appears again (same generation), it reuses the index.
    ///     <b>THREAD SAFETY:</b> The method is fully parallelized with thread-local markers.
    ///     Each thread maintains its own marker arrays to avoid synchronization overhead.
    ///     <b>COMPLEXITY:</b> O(sum of nnodes²) time, O(nodeCount) space per thread.
    /// </remarks>
    /// <example>
    ///     <code>
    ///     var m2m = new M2M();
    ///     // ... populate m2m with element-node connectivity ...
    ///     var cliques = m2m.GetCliques();
    ///     // For element e with nodes [n0, n1, n2]:
    ///     // cliques[e] contains 9 values (3×3 matrix flattened)
    ///     // cliques[e][j + i*3] is the clique index for pair (nodes[i], nodes[j])
    ///     </code>
    /// </example>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    public static List<List<int>> GetCliques(O2M nodesFromElement, O2M elementsFromNode)
    {
        return GetCliquesCore(nodesFromElement, elementsFromNode, strict: false);
    }

    /// <summary>
    ///     Computes clique indices with strict validation that all node indices are in range.
    /// </summary>
    /// <remarks>
    ///     P0 FIX 1.2: Unlike <see cref="GetCliques"/>, this method validates that every node
    ///     index in nodesFromElement is within [0, elementsFromNode.Count). Invalid indices
    ///     typically indicate upstream bugs (stale transpose, broken topology).
    /// </remarks>
    /// <exception cref="InvalidOperationException">
    ///     Thrown when any element contains a node index outside [0, elementsFromNode.Count).
    /// </exception>
    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    public static List<List<int>> GetCliquesStrict(O2M nodesFromElement, O2M elementsFromNode)
    {
        return GetCliquesCore(nodesFromElement, elementsFromNode, strict: true);
    }

    [MethodImpl(AggressiveOptimization)]
    [SkipLocalsInit]
    private static List<List<int>> GetCliquesCore(O2M nodesFromElement, O2M elementsFromNode, bool strict)
    {
        ArgumentNullException.ThrowIfNull(nodesFromElement);
        ArgumentNullException.ThrowIfNull(elementsFromNode);

        var elementCount = nodesFromElement.Count;
        var nodeCount = elementsFromNode.Count;

        // P0 FIX 1.2: Strict validation - verify all node indices are within range
        if (strict)
        {
            for (var e = 0; e < elementCount; e++)
            {
                var nodes = nodesFromElement[e];
                for (var j = 0; j < nodes.Length; j++)
                {
                    var node = nodes[j];
                    if (node < 0 || node >= nodeCount)
                        throw new InvalidOperationException(
                            $"Element {e} contains node index {node} (at position {j}) which is outside " +
                            $"the valid range [0, {nodeCount}). This indicates an inconsistency between " +
                            $"nodesFromElement and elementsFromNode (stale transpose or corrupted topology). " +
                            $"Use GetCliques() instead if you want to skip validation.");
                }
            }
        }

        if (elementCount == 0 || nodeCount == 0)
            return new List<List<int>>(0);

        // Pre-compute node positions within elements (OPTIMIZED)
        var nodeLocation = GetNodePositions(nodesFromElement, elementsFromNode);

        // Allocate result structure with exact sizes
        var cliques = new List<List<int>>(elementCount);
        for (var e = 0; e < elementCount; e++)
        {
            var ns = nodesFromElement._adjacencies[e].Count;
            var list = new List<int>(ns * ns);
            CollectionsMarshal.SetCount(list, ns * ns);
            cliques.Add(list);
        }

        // P0.2 FIX: Estimate max neighborhood size for touched buffer sizing.
        // The marker/markerGen arrays are nodeCount-sized for O(1) lookup, and the
        // touched buffer must also be nodeCount-sized since a hub node's 1-ring neighborhood
        // can contain up to nodeCount unique node2 values.
        const int PoolThreshold = 50000;
        var maxTouchedEstimate = nodeCount;

        if (nodeCount >= nodesFromElement.ParallelizationThreshold)
        {
            if (nodeCount <= PoolThreshold)
                Parallel.For(0, nodeCount,
                    ParallelConfig.Options,
                    () =>
                    {
                        var marker = new int[nodeCount];
                        var markerGen = new int[nodeCount];
                        Array.Fill(markerGen, -1);
                        var touched = new int[maxTouchedEstimate];
                        return (marker, markerGen, touched);
                    },
                    (node1, _, local) =>
                    {
                        ProcessNode1ForCliques(
                            node1,
                            nodesFromElement,
                            elementsFromNode,
                            nodeLocation,
                            cliques,
                            local.marker,
                            local.markerGen,
                            local.touched,
                            out var touchedCount);
                        // P0.2 FIX: Clear only the indices we actually touched
                        for (var i = 0; i < touchedCount; i++)
                        {
                            local.markerGen[local.touched[i]] = -1;
                            local.marker[local.touched[i]] = 0;
                        }
                        return local;
                    },
                    _ => { });
            else
                Parallel.For(0, nodeCount,
                    ParallelConfig.Options,
                    () =>
                    {
                        var marker = ArrayPool<int>.Shared.Rent(nodeCount);
                        var markerGen = ArrayPool<int>.Shared.Rent(nodeCount);
                        Array.Clear(marker, 0, nodeCount);
                        Array.Fill(markerGen, -1, 0, nodeCount);
                        var touched = ArrayPool<int>.Shared.Rent(maxTouchedEstimate);
                        return (marker, markerGen, touched);
                    },
                    (node1, _, local) =>
                    {
                        ProcessNode1ForCliques(
                            node1,
                            nodesFromElement,
                            elementsFromNode,
                            nodeLocation,
                            cliques,
                            local.marker,
                            local.markerGen,
                            local.touched,
                            out var touchedCount);
                        for (var i = 0; i < touchedCount; i++)
                        {
                            local.markerGen[local.touched[i]] = -1;
                            local.marker[local.touched[i]] = 0;
                        }
                        return local;
                    },
                    local =>
                    {
                        ArrayPool<int>.Shared.Return(local.marker);
                        ArrayPool<int>.Shared.Return(local.markerGen);
                        ArrayPool<int>.Shared.Return(local.touched);
                    });
        }
        else
        {
            // Sequential version for small graphs
            var marker = new int[nodeCount];
            var markerGen = new int[nodeCount];
            Array.Fill(markerGen, -1);
            var touched = new int[maxTouchedEstimate];

            for (var node1 = 0; node1 < nodeCount; node1++)
            {
                ProcessNode1ForCliques(
                    node1,
                    nodesFromElement,
                    elementsFromNode,
                    nodeLocation,
                    cliques,
                    marker,
                    markerGen,
                    touched,
                    out var touchedCount);
                // P0.2 FIX: Clear only touched indices
                for (var i = 0; i < touchedCount; i++)
                {
                    markerGen[touched[i]] = -1;
                    marker[touched[i]] = 0;
                }
            }
        }

        return cliques;
    }

    /// <summary>
    ///     Processes a single node1 for clique computation.
    /// </summary>
    /// <remarks>
    ///     <b>P0.2 FIX:</b> For FEM-typical meshes (element arity ≤ ~20), each node1's neighborhood
    ///     is small. This method now uses the marker arrays only for the indices actually encountered
    ///     (tracked via a touched list), making it proportional to local neighborhood size rather
    ///     than global nodeCount. The marker arrays are still passed in for O(1) lookup, but only
    ///     touched indices are cleared after each node1.
    /// </remarks>
    [MethodImpl(AggressiveOptimization | AggressiveInlining)]
    private static void ProcessNode1ForCliques(
        int node1,
        O2M nodesFromElement,
        O2M elementsFromNode,
        List<List<int>> nodeLocation,
        List<List<int>> cliques,
        int[] marker,
        int[] markerGen,
        int[] touched,
        out int touchedCount)
    {
        var generation = node1;
        var nnz = 0;
        touchedCount = 0;

        // Get all elements containing node1
        var node1Elements = elementsFromNode[node1];
        var node1Locations = nodeLocation[node1];
        var numElements = node1Elements.Length;

        for (var le = 0; le < numElements; le++)
        {
            var lnode1 = node1Locations[le];
            var element = node1Elements[le];

            var elementNodes = nodesFromElement[element];
            var esize = elementNodes.Length;

            var cliqueSpan = CollectionsMarshal.AsSpan(cliques[element]);

            for (var lnode2 = 0; lnode2 < esize; lnode2++)
            {
                var node2 = elementNodes[lnode2];

                if (markerGen[node2] != generation)
                {
                    // First time seeing node2 for this node1 — track it
                    markerGen[node2] = generation;
                    marker[node2] = nnz;
                    cliqueSpan[lnode2 + lnode1 * esize] = nnz;
                    // P0.2 FIX: Track touched indices for efficient clearing
                    touched[touchedCount++] = node2;
                    nnz++;
                }
                else
                {
                    cliqueSpan[lnode2 + lnode1 * esize] = marker[node2];
                }
            }
        }
    }

    #endregion

    // ============================================================================
    // GRAPH ALGORITHMS - BFS, DIJKSTRA
    // Integrated: December 15, 2024
    // ============================================================================

    #region Graph Algorithms - BFS

    /// <summary>
    ///     Performs breadth-first search starting from a single element.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="visitor">Optional visitor callback invoked for each discovered element.</param>
    /// <returns>List of elements in BFS order.</returns>
    /// <remarks>
    ///     <b>ALGORITHM:</b> Standard BFS using queue. Visits elements level by level.
    ///     <b>COMPLEXITY:</b> O(V × E) without transpose, O(V + E) with transpose.
    ///     <b>THREAD SAFETY:</b> Safe for concurrent reads. Not safe during modifications.
    ///     <b>USAGE:</b> Flood fill, shortest unweighted paths, level sets.
    ///     <b>PERFORMANCE:</b> For better performance, use the overload with transpose parameter.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public List<int> BreadthFirstSearch(int startElement, Action<int, int>? visitor = null)
    {
        return BreadthFirstSearch(startElement, null, visitor);
    }

    /// <summary>
    ///     Performs breadth-first search using pre-computed transpose for O(1) neighbor lookup.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="transpose">Pre-computed transpose (node → elements). If null, uses O(n) fallback.</param>
    /// <param name="visitor">Optional visitor callback invoked for each discovered element.</param>
    /// <returns>List of elements in BFS order.</returns>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> When transpose is provided, neighbor lookup is O(1) instead of O(n).
    ///     Use <see cref="Transpose" /> to compute the transpose once, then pass it to multiple calls.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public List<int> BreadthFirstSearch(int startElement, O2M? transpose, Action<int, int>? visitor = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(startElement);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(startElement, Count);

        var visited = new bool[Count];
        var result = new List<int>();
        var queue = new Queue<int>();

        queue.Enqueue(startElement);
        visited[startElement] = true;
        result.Add(startElement);
        visitor?.Invoke(startElement, 0);

        var currentDepth = 0;
        var elementsAtCurrentDepth = 1;
        var elementsAtNextDepth = 0;

        while (queue.Count > 0)
        {
            var current = queue.Dequeue();
            var nodes = this[current];

            foreach (var node in nodes)
                // Use transpose if available for O(1) lookup, otherwise O(n) scan
                if (transpose != null && node >= 0 && node < transpose.Count)
                {
                    var connectedElements = transpose[node];
                    foreach (var neighbor in connectedElements)
                    {
                        if (visited[neighbor]) continue;

                        visited[neighbor] = true;
                        queue.Enqueue(neighbor);
                        result.Add(neighbor);
                        elementsAtNextDepth++;
                        visitor?.Invoke(neighbor, currentDepth + 1);
                    }
                }
                else
                {
                    // Fallback: O(n) scan
                    for (var elem = 0; elem < Count; elem++)
                    {
                        if (visited[elem]) continue;
                        if (!this[elem].Contains(node)) continue;

                        visited[elem] = true;
                        queue.Enqueue(elem);
                        result.Add(elem);
                        elementsAtNextDepth++;
                        visitor?.Invoke(elem, currentDepth + 1);
                    }
                }

            elementsAtCurrentDepth--;
            if (elementsAtCurrentDepth == 0)
            {
                currentDepth++;
                elementsAtCurrentDepth = elementsAtNextDepth;
                elementsAtNextDepth = 0;
            }
        }

        return result;
    }

    /// <summary>
    ///     Performs BFS and returns distances from start element.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <returns>Dictionary mapping element indices to their BFS distance (hop count).</returns>
    /// <remarks>
    ///     <b>DISTANCES:</b> Unweighted shortest path distances (hop count).
    ///     <b>UNREACHABLE:</b> Elements not in the dictionary are unreachable.
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public Dictionary<int, int> BreadthFirstDistances(int startElement)
    {
        return BreadthFirstDistances(startElement, null);
    }

    /// <summary>
    ///     Performs BFS and returns distances, using pre-computed transpose for efficiency.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="transpose">Pre-computed transpose (node → elements). If null, uses O(n) fallback.</param>
    /// <returns>Dictionary mapping element indices to their BFS distance (hop count).</returns>
    [MethodImpl(AggressiveOptimization)]
    public Dictionary<int, int> BreadthFirstDistances(int startElement, O2M? transpose)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(startElement);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(startElement, Count);

        var distances = new Dictionary<int, int> { [startElement] = 0 };
        var queue = new Queue<int>();
        queue.Enqueue(startElement);

        while (queue.Count > 0)
        {
            var current = queue.Dequeue();
            var currentDist = distances[current];
            var nodes = this[current];

            foreach (var node in nodes)
                if (transpose != null && node >= 0 && node < transpose.Count)
                {
                    var connectedElements = transpose[node];
                    foreach (var neighbor in connectedElements)
                    {
                        if (distances.ContainsKey(neighbor)) continue;

                        distances[neighbor] = currentDist + 1;
                        queue.Enqueue(neighbor);
                    }
                }
                else
                {
                    for (var elem = 0; elem < Count; elem++)
                    {
                        if (distances.ContainsKey(elem)) continue;
                        if (!this[elem].Contains(node)) continue;

                        distances[elem] = currentDist + 1;
                        queue.Enqueue(elem);
                    }
                }
        }

        return distances;
    }

    #endregion

    #region Graph Algorithms - Dijkstra

    /// <summary>
    ///     Computes shortest paths from a start element using Dijkstra's algorithm.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="edgeWeight">Function computing edge weight between two elements sharing a node.</param>
    /// <returns>Dictionary mapping element indices to (distance, predecessor) tuples.</returns>
    /// <remarks>
    ///     <b>ALGORITHM:</b> Dijkstra with binary heap priority queue.
    ///     <b>COMPLEXITY:</b> O((V + E) log V) with binary heap.
    ///     <b>WEIGHTS:</b> Edge weight function receives (fromElement, toElement, sharedNode).
    ///     <b>NEGATIVE WEIGHTS:</b> Not supported. Use Bellman-Ford if needed.
    ///     <b>RECONSTRUCTION:</b> Use predecessor to reconstruct shortest path.
    ///     <example>
    ///         // Geodesic distance on mesh
    ///         var distances = mesh.DijkstraShortestPaths(startElem, (from, to, node) =>
    ///         Vector3.Distance(positions[from], positions[to]));
    ///     </example>
    /// </remarks>
    [MethodImpl(AggressiveOptimization)]
    public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths(
        int startElement,
        Func<int, int, int, double> edgeWeight)
    {
        return DijkstraShortestPaths(startElement, null, edgeWeight);
    }

    /// <summary>
    ///     Computes shortest paths using pre-computed transpose for efficient neighbor lookup.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="transpose">Pre-computed transpose (node → elements). If null, uses O(n) fallback.</param>
    /// <param name="edgeWeight">Function computing edge weight between two elements sharing a node.</param>
    /// <returns>Dictionary mapping element indices to (distance, predecessor) tuples.</returns>
    [MethodImpl(AggressiveOptimization)]
    public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths(
        int startElement,
        O2M? transpose,
        Func<int, int, int, double> edgeWeight)
    {
        ArgumentNullException.ThrowIfNull(edgeWeight);
        ArgumentOutOfRangeException.ThrowIfNegative(startElement);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(startElement, Count);

        var distances = new Dictionary<int, (double Distance, int Predecessor)>
        {
            [startElement] = (0.0, -1)
        };

        var pq = new PriorityQueue<int, double>();
        pq.Enqueue(startElement, 0.0);

        var finalized = new HashSet<int>();

        while (pq.Count > 0)
        {
            var current = pq.Dequeue();

            if (finalized.Contains(current)) continue;
            finalized.Add(current);

            var currentDist = distances[current].Distance;
            var nodes = this[current];

            foreach (var node in nodes)
                if (transpose != null && node >= 0 && node < transpose.Count)
                {
                    var connectedElements = transpose[node];
                    foreach (var neighbor in connectedElements)
                    {
                        if (finalized.Contains(neighbor)) continue;
                        if (neighbor == current) continue;

                        var weight = edgeWeight(current, neighbor, node);
                        if (weight < 0)
                            throw new ArgumentException(
                                $"Negative edge weight {weight} detected between elements {current} and {neighbor}. " +
                                "Dijkstra's algorithm does not support negative weights.");

                        var newDist = currentDist + weight;

                        if (!distances.TryGetValue(neighbor, out var neighborData) ||
                            newDist < neighborData.Distance)
                        {
                            distances[neighbor] = (newDist, current);
                            pq.Enqueue(neighbor, newDist);
                        }
                    }
                }
                else
                {
                    for (var neighbor = 0; neighbor < Count; neighbor++)
                    {
                        if (finalized.Contains(neighbor)) continue;
                        if (neighbor == current) continue;
                        if (!this[neighbor].Contains(node)) continue;

                        var weight = edgeWeight(current, neighbor, node);
                        if (weight < 0)
                            throw new ArgumentException(
                                $"Negative edge weight {weight} detected between elements {current} and {neighbor}. " +
                                "Dijkstra's algorithm does not support negative weights.");

                        var newDist = currentDist + weight;

                        if (!distances.TryGetValue(neighbor, out var neighborData) ||
                            newDist < neighborData.Distance)
                        {
                            distances[neighbor] = (newDist, current);
                            pq.Enqueue(neighbor, newDist);
                        }
                    }
                }
        }

        return distances;
    }

    /// <summary>
    ///     Reconstructs shortest path from Dijkstra result.
    /// </summary>
    /// <param name="dijkstraResult">Result from DijkstraShortestPaths.</param>
    /// <param name="targetElement">Target element.</param>
    /// <returns>Path from start to target, or null if unreachable.</returns>
    [MethodImpl(AggressiveOptimization)]
    public static List<int>? ReconstructPath(
        Dictionary<int, (double Distance, int Predecessor)> dijkstraResult,
        int targetElement)
    {
        if (!dijkstraResult.ContainsKey(targetElement))
            return null;

        var path = new List<int>();
        var current = targetElement;

        while (current != -1)
        {
            path.Add(current);
            current = dijkstraResult[current].Predecessor;
        }

        path.Reverse();
        return path;
    }

    #endregion
}


/// <summary>
///     Represents a thread-safe many-to-many relationship with cached transpose and position lookups.
/// </summary>
/// <remarks>
///     <b>THREAD SAFETY:</b> All methods are thread-safe. Uses ReaderWriterLockSlim to allow
///     concurrent read operations while ensuring exclusive access for modifications.
///     Multiple threads can execute query operations simultaneously once the cache is synchronized.
///     <b>DESIGN:</b> Wraps an O2M instance (composition) and adds:
///     - Thread-safe access with reader-writer locking
///     - Automatic transpose caching (ElementsFromNode)
///     - Position caching (ElementLocations, NodeLocations) for spatial queries
///     - Smart synchronization with batch operation support
///     <b>PERFORMANCE:</b> Query operations after initial synchronization are highly concurrent.
///     Modifications invalidate cache, requiring resynchronization on next query.
///     Use batch operations to amortize synchronization cost across multiple modifications.
///     <b>DISPOSAL:</b> Implements IDisposable with finalizer support. Disposal is OPTIONAL -
///     the finalizer will automatically clean up resources when the object is garbage collected.
///     Explicit disposal via using statement or Dispose() is recommended for deterministic cleanup:
///     <code>
/// // Recommended (immediate cleanup):
/// using (var m2m = new M2M()) { }
/// // Also valid (finalizer cleans up eventually):
/// var m2m = new M2M();
/// </code>
///     <b>USAGE GUIDANCE:</b>
///     Use M2M when you need:
///     - Thread-safe multi-threaded access
///     - Frequent transpose queries (GetElementsWithNodes, GetNeighbors)
///     - Cached position lookups for spatial algorithms
///     Use plain O2M when you have:
///     - Single-threaded access
///     - No need for transpose caching
///     - Maximum performance requirements (no locking overhead)
/// </remarks>
public sealed class M2M : IComparable<M2M>, IEquatable<M2M>, IDisposable
{
    #region Clique Computation (Add to M2M class)

    /// <summary>
    ///     Computes clique indices for finite element matrix assembly.
    /// </summary>
    /// <returns>
    ///     A list where result[element] is a flattened ns×ns array mapping
    ///     local node pairs to unique clique indices.
    /// </returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>PURPOSE:</b> Essential for sparse matrix assembly in FEM applications.
    ///     Identifies which node pairs share connectivity patterns across elements.
    ///     <b>THREAD SAFETY:</b> Thread-safe with automatic synchronization.
    ///     <b>PERFORMANCE:</b> First call after modification triggers synchronization.
    ///     Subsequent calls are efficient.
    ///     <b>USAGE:</b>
    ///     <code>
    ///         using var m2m = new M2M();
    ///         // ... build mesh connectivity ...
    ///         var cliques = m2m.GetCliques();
    ///         // Use for matrix assembly:
    ///         for (int e = 0; e < m2m.Count; e++)
    ///                               {
    ///                               var elemNodes= m2m.GetNodesForElement( e);
    ///                               int ns= elemNodes.Count;
    ///                               for ( int i= 0; i < ns; i++)
    ///                                                     for ( int j= 0; j < ns; j++)
    ///                                                                           {
    ///                                                                           int cliqueIdx= cliques[ e][ j + i * ns];
    ///                                                                           // Accumulate element matrix contribution at
    ///                                                                           cliqueIdx
    ///                 }
    ///         }
    ///         </code>
    /// </remarks>
    public List<List<int>> GetCliques()
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return O2M.GetCliquesStrict(_o2m, _elementsFromNode!);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Fields and Properties

    /// <summary>
    ///     The wrapped O2M structure containing the element-to-node adjacency.
    /// </summary>
    private readonly O2M _o2m;

    /// <summary>
    ///     Reader-writer lock for thread-safe access.
    ///     Allows multiple concurrent readers or single writer.
    /// </summary>
    private readonly ReaderWriterLockSlim _rwLock = new(LockRecursionPolicy.SupportsRecursion);

    /// <summary>
    ///     Tracks whether resources have been disposed (0 = not disposed, 1 = disposed).
    ///     Uses int for Interlocked.CompareExchange compatibility.
    /// </summary>
    private int _disposed;

    /// <summary>
    ///     Monotonic instance ID for deterministic lock ordering.
    ///     Guarantees a unique, consistent ordering even when RuntimeHelpers.GetHashCode collides.
    /// </summary>
    private static long _nextInstanceId;
    private readonly long _instanceId = Interlocked.Increment(ref _nextInstanceId);

    /// <summary>
    ///     Tracks whether disposal completed unsuccessfully (lock timeout).
    ///     When true, the object is in an inconsistent state and should not be used.
    /// </summary>
    private volatile bool _disposalIncomplete;

    /// <summary>
    ///     Gets whether disposal was attempted but did not complete successfully.
    ///     When true, the lock could not be acquired within the timeout, indicating
    ///     that other threads may still be accessing the object.
    /// </summary>
    /// <remarks>
    ///     <b>RECOVERY:</b> If this returns true, ensure all other threads have completed
    ///     their operations and call Dispose() again to retry cleanup.
    /// </remarks>
    public bool IsDisposalIncomplete => _disposalIncomplete;

    /// <summary>
    ///     Cached transpose for fast node-to-element lookups.
    ///     Lazily computed on first query after structure modification.
    /// </summary>
    private O2M? _elementsFromNode;

    /// <summary>
    ///     Element position cache for spatial queries.
    ///     Lazily computed only when accessed via ElementLocations property.
    /// </summary>
    private IReadOnlyList<IReadOnlyList<int>>? _elemeloc;

    /// <summary>
    ///     Node position cache for spatial queries.
    ///     Lazily computed only when accessed via NodeLocations property.
    /// </summary>
    private IReadOnlyList<IReadOnlyList<int>>? _nodeloc;

    /// <summary>
    ///     Tracks whether transpose cache is synchronized with current structure.
    ///     MUST be volatile for double-check locking pattern in EnterSynchronizedReadLock().
    /// </summary>
    private volatile bool _isInSync;

    /// <summary>
    ///     Tracks whether position caches (element locations and node locations) have been computed.
    /// </summary>
    private bool _positionCachesComputed;

    /// <summary>
    ///     Nesting level for batch operations.
    ///     When > 0, synchronization is deferred until batch completes.
    /// </summary>
    private int _batchNesting;

    /// <summary>
    ///     Gets a clone of the cached transpose for node-to-element lookups.
    /// </summary>
    /// <remarks>
    ///     Thread Safety: Safe to call from any thread. Automatically synchronizes if needed.
    ///     Performance: O(n × m) for clone. If you need repeated access, cache the returned object.
    ///     Returns a defensive copy to prevent external mutation of internal state.
    ///     The returned O2M is independent and can be freely modified without affecting this M2M.
    /// </remarks>
    /// <exception cref="InvalidOperationException">Thrown when accessed during a batch update.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public O2M ElementsFromNode
    {
        get
        {
            ThrowIfDisposed();
            ThrowIfInBatch();
            EnterSynchronizedReadLock();
            try
            {
                // Return a clone to prevent external mutation of internal cache
                return (O2M)_elementsFromNode!.Clone();
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    /// <summary>
    ///     Executes an action with direct read access to the transpose cache (no cloning).
    /// </summary>
    /// <param name="action">Action to execute with the transpose. Do NOT store references.</param>
    /// <remarks>
    ///     <para>
    ///         <b>PERFORMANCE:</b> Unlike <see cref="ElementsFromNode" /> which clones O(n×m),
    ///         this provides O(1) access for read-only operations in hot paths.
    ///     </para>
    ///     <para>
    ///         <b>SAFETY:</b> The O2M reference is only valid during action execution.
    ///         Do NOT store references to the O2M or its internal lists.
    ///         Read lock is held during execution - do not call other M2M methods.
    ///     </para>
    /// </remarks>
    /// <exception cref="InvalidOperationException">Thrown when accessed during a batch update.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public void WithElementsFromNode(Action<O2M> action)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        ArgumentNullException.ThrowIfNull(action);
        EnterSynchronizedReadLock();
        try
        {
            action(_elementsFromNode!);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes a function with direct read access to the transpose cache (no cloning).
    /// </summary>
    /// <typeparam name="TResult">The result type.</typeparam>
    /// <param name="func">Function to execute with the transpose.</param>
    /// <returns>The function result.</returns>
    /// <remarks>
    ///     <para>
    ///         <b>PERFORMANCE:</b> Zero-copy access for hot-path queries.
    ///     </para>
    ///     <para>
    ///         <b>SAFETY:</b> Do NOT return or store references to the O2M or its contents.
    ///     </para>
    /// </remarks>
    public TResult WithElementsFromNode<TResult>(Func<O2M, TResult> func)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        ArgumentNullException.ThrowIfNull(func);
        EnterSynchronizedReadLock();
        try
        {
            return func(_elementsFromNode!);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets a read-only clone of the underlying element-to-node adjacency.
    /// </summary>
    /// <remarks>
    ///     Thread Safety: Safe to call from any thread.
    ///     Performance: O(n × m) for deep copy.
    /// </remarks>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public O2M NodesFromElement
    {
        get
        {
            ThrowIfDisposed();
            _rwLock.EnterReadLock();
            try
            {
                return (O2M)_o2m.Clone();
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    #region Async Operations (Priority 3)

    /// <summary>
    ///     Asynchronously gets the transpose (elements from node) for large datasets.
    /// </summary>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <returns>A task that represents the asynchronous operation, containing the transposed O2M.</returns>
    /// <remarks>
    ///     Async version that offloads computation to the thread pool.
    ///     <b>LOCKING BEHAVIOR:</b> If the transpose needs to be computed, this method
    ///     acquires a write lock that blocks all other operations until the transpose
    ///     completes. For very large meshes (&gt;1M elements), transpose computation can
    ///     take several seconds during which the structure is locked.
    ///     <b>PERFORMANCE TIP:</b> For maximum concurrency, pre-compute the transpose by
    ///     calling EnsureSynchronized() before making async calls. After initial synchronization,
    ///     this method only takes a brief read lock to clone the cached transpose.
    ///     <b>USAGE:</b>
    ///     <code>
    ///         // First call after modification: May hold write lock during transpose
    ///         var transpose = await m2m.ElementsFromNodeAsync(cancellationToken);
    ///         
    ///         // Alternative for better concurrency:
    ///         m2m.EnsureSynchronized();  // Sync under lock
    ///         var transpose = await m2m.ElementsFromNodeAsync();  // Fast clone only
    ///         </code>
    /// </remarks>
    /// <exception cref="InvalidOperationException">Thrown when called during a batch update.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <exception cref="OperationCanceledException">Thrown if cancellation is requested.</exception>
    public async Task<O2M> ElementsFromNodeAsync(CancellationToken cancellationToken = default)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();

        return await Task.Run(() =>
        {
            cancellationToken.ThrowIfCancellationRequested();

            // P0.3 FIX: Use EnterSynchronizedReadLock to close TOCTOU gap
            ThrowIfDisposed();
            EnterSynchronizedReadLock();

            O2M result;
            try
            {
                result = (O2M)_elementsFromNode!.Clone();
            }
            finally
            {
                _rwLock.ExitReadLock();
            }

            return result;
        }, cancellationToken).ConfigureAwait(false);
    }

    #endregion

    /// <summary>
    ///     Gets or sets the threshold for using parallel algorithms in the underlying O2M.
    /// </summary>
    /// <remarks>
    ///     Operations on structures with more elements than this threshold will use parallelization.
    ///     Default is 4096. Tune based on your hardware and typical problem sizes.
    ///     Thread Safety: Property is thread-safe.
    /// </remarks>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public int ParallelizationThreshold
    {
        get
        {
            ThrowIfDisposed();
            _rwLock.EnterReadLock();
            try
            {
                return _o2m.ParallelizationThreshold;
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
        set
        {
            ThrowIfDisposed();
            _rwLock.EnterWriteLock();
            try
            {
                _o2m.ParallelizationThreshold = value;
            }
            finally
            {
                _rwLock.ExitWriteLock();
            }
        }
    }

    /// <summary>
    ///     Gets or sets the default parallelization threshold for new O2M instances.
    /// </summary>
    /// <remarks>
    ///     This static property affects all newly created O2M instances.
    ///     Thread Safety: This property is thread-safe.
    /// </remarks>
    public static int DefaultParallelizationThreshold
    {
        get => O2M.DefaultParallelizationThreshold;
        set => O2M.DefaultParallelizationThreshold = value;
    }

    /// <summary>
    ///     Gets the element position cache for spatial queries.
    /// </summary>
    /// <remarks>
    ///     For each element e and position p in e's adjacency list, ElementLocations[e][p] gives
    ///     the position of element e within the adjacency list of node at position p.
    ///     Thread Safety: Safe to call from any thread.
    ///     Performance: O(1) if cached, O(n × m × log k) if synchronization required.
    /// </remarks>
    /// <exception cref="InvalidOperationException">Thrown when accessed during a batch update.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public IReadOnlyList<IReadOnlyList<int>> ElementLocations
    {
        get
        {
            ThrowIfDisposed();
            ThrowIfInBatch();
            EnterPositionCachedReadLock();
            try
            {
                return _elemeloc!;
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    /// <summary>
    ///     Gets the node position cache for spatial queries.
    /// </summary>
    /// <remarks>
    ///     For each node n, NodeLocations[n] gives the positions where node n appears
    ///     within each element's adjacency list.
    ///     Thread Safety: Safe to call from any thread.
    ///     Performance: O(1) if cached, O(n × m) if synchronization required.
    /// </remarks>
    /// <exception cref="InvalidOperationException">Thrown when accessed during a batch update.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public IReadOnlyList<IReadOnlyList<int>> NodeLocations
    {
        get
        {
            ThrowIfDisposed();
            ThrowIfInBatch();
            EnterPositionCachedReadLock();
            try
            {
                return _nodeloc!;
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    /// <summary>
    ///     Gets the number of elements in the structure.
    /// </summary>
    /// <remarks>
    ///     Thread Safety: Safe to call from any thread.
    ///     Time Complexity: O(1)
    /// </remarks>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public int Count
    {
        get
        {
            ThrowIfDisposed();
            _rwLock.EnterReadLock();
            try
            {
                return _o2m.Count;
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    /// <summary>
    ///     Gets a read-only view of the adjacency list for the specified element.
    /// </summary>
    /// <param name="rowIndex">The element index.</param>
    /// <returns>Read-only list of node indices connected to this element.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when rowIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Safe to call from any thread.
    ///     Time Complexity: O(m) where m = nodes in element (defensive copy via ToArray).
    /// </remarks>
    public IReadOnlyList<int> this[int rowIndex]
    {
        get
        {
            ThrowIfDisposed();
            _rwLock.EnterReadLock();
            try
            {
                return _o2m[rowIndex].ToArray();
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    /// <summary>
    ///     Gets the node value at the specified element and position.
    /// </summary>
    /// <param name="rowIndex">The element index.</param>
    /// <param name="columnIndex">The position within the element's adjacency list.</param>
    /// <returns>The node index at the specified position.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when indices are out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Safe to call from any thread.
    ///     Time Complexity: O(1)
    /// </remarks>
    public int this[int rowIndex, int columnIndex]
    {
        get
        {
            ThrowIfDisposed();
            _rwLock.EnterReadLock();
            try
            {
                return _o2m[rowIndex, columnIndex];
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    #endregion

    #region Constructors

    /// <summary>
    ///     Initializes an empty M2M instance.
    /// </summary>
    /// <remarks>
    ///     Time Complexity: O(1)
    /// </remarks>
    public M2M()
    {
        _o2m = new O2M();
        _isInSync = false;
        _positionCachesComputed = false;
    }

    /// <summary>
    ///     Initializes an M2M instance with pre-allocated capacity.
    /// </summary>
    /// <param name="reservedCapacity">Initial capacity for the adjacency list.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when reservedCapacity is negative.</exception>
    /// <remarks>
    ///     Time Complexity: O(1)
    ///     Use this constructor when you know the approximate number of elements
    ///     to avoid reallocation during construction.
    /// </remarks>
    public M2M(int reservedCapacity)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(reservedCapacity);
        _o2m = new O2M(reservedCapacity);
        _isInSync = false;
        _positionCachesComputed = false;
    }

    /// <summary>
    ///     Initializes an M2M instance from an existing adjacency list.
    /// </summary>
    /// <param name="adjacencies">The adjacency list to wrap.</param>
    /// <exception cref="ArgumentNullException">Thrown when adjacencies is null.</exception>
    /// <remarks>
    ///     Time Complexity: O(1) - the list is used directly by the wrapped O2M
    ///     WARNING: The provided list becomes owned by the M2M instance.
    /// </remarks>
    public M2M(List<List<int>> adjacencies)
    {
        ArgumentNullException.ThrowIfNull(adjacencies);
        _o2m = new O2M(adjacencies);
        _isInSync = false;
        _positionCachesComputed = false;
    }

    /// <summary>
    ///     Initializes an M2M instance by wrapping an existing O2M.
    /// </summary>
    /// <param name="o2m">The O2M instance to wrap.</param>
    /// <exception cref="ArgumentNullException">Thrown when o2m is null.</exception>
    /// <remarks>
    ///     Time Complexity: O(n × m) for deep copy
    ///     The O2M is cloned to ensure the M2M has exclusive ownership.
    /// </remarks>
    public M2M(O2M o2m)
    {
        ArgumentNullException.ThrowIfNull(o2m);
        _o2m = (O2M)o2m.Clone();
        _isInSync = false;
        _positionCachesComputed = false;
    }

    /// <summary>
    ///     Private constructor that takes ownership of an O2M without cloning.
    /// </summary>
    /// <remarks>
    ///     <b>P0.3 FIX:</b> Used by <see cref="Clone"/>, <see cref="GetElementToElementGraph"/>,
    ///     and <see cref="GetNodeToNodeGraph"/> to avoid double-cloning. The caller guarantees
    ///     the O2M is freshly created and not shared with any other owner.
    /// </remarks>
    private M2M(O2M o2m, bool takeOwnership)
    {
        Debug.Assert(takeOwnership, "This constructor takes ownership — pass true");
        _o2m = o2m;
        _isInSync = false;
        _positionCachesComputed = false;
    }

    #endregion

    #region Core Modification Methods

    /// <summary>
    ///     Appends a new element with the specified nodes.
    /// </summary>
    /// <param name="nodes">The list of nodes for the new element.</param>
    /// <returns>The index of the newly added element.</returns>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(1) amortized, plus synchronization cost if not batching.
    ///     Use batch operations when adding many elements to avoid repeated synchronization:
    ///     <code>
    /// using (m2m.BeginBatchUpdate())
    /// {
    ///     m2m.AppendElement([1, 2, 3]);
    ///     m2m.AppendElement([4, 5, 6]);
    /// } // Single synchronization on dispose
    /// </code>
    /// </remarks>
    public int AppendElement(List<int> nodes)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(nodes);
        _rwLock.EnterWriteLock();
        try
        {
            InvalidateCache();
            return _o2m.AppendElement(nodes);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Appends multiple elements in a batch.
    /// </summary>
    /// <param name="nodes">Array of node lists to append.</param>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null or contains null lists.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(k) where k is the number of elements to append.
    /// </remarks>
    public void AppendElements(params List<int>[] nodes)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(nodes);
        if (nodes.Length == 0) return;
        _rwLock.EnterWriteLock();
        try
        {
            InvalidateCache();
            _o2m.AppendElements(nodes);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Appends a single node to an existing element.
    /// </summary>
    /// <param name="elementIndex">The element index to modify.</param>
    /// <param name="nodeValue">The node value to append.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementIndex is out of range or nodeValue is negative.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(1) amortized.
    /// </remarks>
    public void AppendNodeToElement(int elementIndex, int nodeValue)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            InvalidateCache();
            _o2m.AppendNodeToElement(elementIndex, nodeValue);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Removes a specific node from an element.
    /// </summary>
    /// <param name="elementIndex">The element index to modify.</param>
    /// <param name="nodeValue">The node value to remove.</param>
    /// <returns>True if the node was found and removed; false otherwise.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(m) where m is the number of nodes in the element.
    /// </remarks>
    public bool RemoveNodeFromElement(int elementIndex, int nodeValue)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            var removed = _o2m.RemoveNodeFromElement(elementIndex, nodeValue);
            if (removed) InvalidateCache();
            return removed;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }


    /// <summary>
    ///     Clears all nodes from the specified element.
    /// </summary>
    /// <param name="elementIndex">The element index to clear.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(1) amortized.
    /// </remarks>
    public void ClearElement(int elementIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            _o2m.ClearElement(elementIndex);
            InvalidateCache();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Replaces the nodes for the specified element.
    /// </summary>
    /// <param name="elementIndex">The element index to replace.</param>
    /// <param name="newNodes">The new list of nodes.</param>
    /// <exception cref="ArgumentNullException">Thrown when newNodes is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(1).
    /// </remarks>
    public void ReplaceElement(int elementIndex, List<int> newNodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            _o2m.ReplaceElement(elementIndex, newNodes);
            InvalidateCache();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Compresses elements by selecting a subset according to a mapping.
    /// </summary>
    /// <param name="newToOldElementMap">Maps each new index to an old index to keep.</param>
    /// <exception cref="ArgumentNullException">Thrown when mapping is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(k) where k = newToOldElementMap.Count.
    /// </remarks>
    public void CompressElements(List<int> newToOldElementMap)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(newToOldElementMap);
        if (newToOldElementMap.Count == 0) return;
        _rwLock.EnterWriteLock();
        try
        {
            InvalidateCache();
            _o2m.CompressElements(newToOldElementMap);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Reorders elements according to a permutation mapping.
    /// </summary>
    /// <param name="oldToNewElementMap">Maps each old index to its new index.</param>
    /// <exception cref="ArgumentNullException">Thrown when mapping is null.</exception>
    /// <exception cref="InvalidOperationException">Thrown when mapping is not a valid permutation.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(n) where n = element count.
    /// </remarks>
    public void PermuteElements(List<int> oldToNewElementMap)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(oldToNewElementMap);
        if (oldToNewElementMap.Count == 0) return;
        _rwLock.EnterWriteLock();
        try
        {
            InvalidateCache();
            _o2m.PermuteElements(oldToNewElementMap);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Renumbers all nodes according to a mapping.
    /// </summary>
    /// <param name="oldToNewNodeMap">Maps each old node index to its new index.</param>
    /// <exception cref="ArgumentNullException">Thrown when mapping is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(n × m) where n = elements, m = average nodes per element.
    /// </remarks>
    public void PermuteNodes(List<int> oldToNewNodeMap)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(oldToNewNodeMap);
        if (oldToNewNodeMap.Count == 0) return;
        _rwLock.EnterWriteLock();
        try
        {
            InvalidateCache();
            _o2m.PermuteNodes(oldToNewNodeMap);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Combined operation: compress elements then renumber nodes.
    /// </summary>
    /// <param name="newToOldElementMap">Element compression mapping.</param>
    /// <param name="oldToNewNodeMap">Node renumbering mapping.</param>
    /// <exception cref="ArgumentNullException">Thrown when any mapping is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(k × m) where k = new element count, m = average nodes per element.
    /// </remarks>
    public void RearrangeAfterRenumbering(List<int> newToOldElementMap, List<int> oldToNewNodeMap)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            InvalidateCache();
            _o2m.RearrangeAfterRenumbering(newToOldElementMap, oldToNewNodeMap);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Clears all elements and resets caches.
    /// </summary>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with exclusive write lock.
    ///     Time Complexity: O(1).
    /// </remarks>
    public void ClearAll()
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            _o2m.ClearAll();
            _elementsFromNode = null;
            _elemeloc = null;
            _nodeloc = null;
            _isInSync = false;
            _positionCachesComputed = false;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Query Methods

    /// <summary>
    ///     Gets the nodes connected to the specified element.
    /// </summary>
    /// <param name="elementIndex">The element index.</param>
    /// <returns>Read-only list of node indices.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(m) where m = nodes in element (defensive copy via ToArray).
    /// </remarks>
    public IReadOnlyList<int> GetNodesForElement(int elementIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m[elementIndex].ToArray();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #region Zero-Copy Span Access (Thread-Safe)

    /// <summary>
    ///     Executes an action with zero-copy span access to element nodes.
    /// </summary>
    /// <param name="elementIndex">The element index.</param>
    /// <param name="action">Action to execute with the span. Invoked under read lock.</param>
    /// <exception cref="ArgumentNullException">Thrown when action is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe. The action is executed under read lock.
    ///     Time Complexity: O(1) + O(action execution time).
    ///     <para>
    ///         <b>IMPLEMENTATION NOTE:</b> Requires O2M indexer to return List&lt;int&gt;.
    ///         Falls back to array copy if implementation changes.
    ///     </para>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void WithNodesSpan(int elementIndex, Action<ReadOnlySpan<int>> action)
    {
        ArgumentNullException.ThrowIfNull(action);
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            var nodes = _o2m[elementIndex];
            // O2M now returns ReadOnlySpan<int> directly - just use it
            action(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes a function with zero-copy span access to element nodes and returns a result.
    /// </summary>
    /// <typeparam name="TResult">The type of result to return.</typeparam>
    /// <param name="elementIndex">The element index.</param>
    /// <param name="func">Function to execute with the span. Invoked under read lock.</param>
    /// <returns>The result of the function.</returns>
    /// <exception cref="ArgumentNullException">Thrown when func is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe. The function is executed under read lock.
    ///     Time Complexity: O(1) + O(function execution time).
    ///     <para>
    ///         <b>IMPLEMENTATION NOTE:</b> Requires O2M indexer to return List&lt;int&gt;.
    ///         Falls back to array copy if implementation changes.
    ///     </para>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public TResult WithNodesSpan<TResult>(int elementIndex, Func<ReadOnlySpan<int>, TResult> func)
    {
        ArgumentNullException.ThrowIfNull(func);
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            var nodes = _o2m[elementIndex];
            // O2M now returns ReadOnlySpan<int> directly - just use it
            return func(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    /// <summary>
    ///     Gets the number of nodes in the transpose (max node index + 1).
    /// </summary>
    /// <returns>The count of nodes in the transposed structure, or 0 if empty.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <exception cref="InvalidOperationException">Thrown when called during a batch update.</exception>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> O(sync) + O(1). Much more efficient than ElementsFromNode.Count
    ///     which clones O(n×m) just to get a count.
    /// </remarks>
    public int GetTransposeNodeCount()
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return _elementsFromNode?.Count ?? 0;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the count of elements containing a specific node.
    /// </summary>
    /// <param name="node">The node index to look up.</param>
    /// <returns>The count of elements containing the node, or 0 if node not found.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <exception cref="InvalidOperationException">Thrown when called during a batch update.</exception>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> O(sync) + O(1). Much more efficient than ElementsFromNode[node].Count
    ///     which clones O(n×m) just to get a count.
    /// </remarks>
    public int GetElementCountForNode(int node)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            if ((uint)node >= (uint)_elementsFromNode!.Count)
                return 0;
            return _elementsFromNode[node].Length;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the elements containing a specific node (returns a defensive copy).
    /// </summary>
    /// <param name="node">The node index to look up.</param>
    /// <returns>List of element indices containing the specified node, or empty list if node not found.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <exception cref="InvalidOperationException">Thrown when called during a batch update.</exception>
    /// <remarks>
    ///     <para>
    ///         <b>PERFORMANCE:</b> O(sync) + O(k) where k = elements containing node.
    ///         Much more efficient than <c>ElementsFromNode[node]</c> which clones O(n×m).
    ///     </para>
    ///     <para>
    ///         <b>SAFETY:</b> Returns a defensive copy. For zero-copy access, use
    ///         <see cref="WithElementsForNode(int, Action{IReadOnlyList{int}})" />.
    ///     </para>
    /// </remarks>
    public List<int> GetElementsForNode(int node)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            if ((uint)node >= (uint)_elementsFromNode!.Count)
                return [];

            var internal_list = _elementsFromNode[node];
            // Return defensive copy - convert span to list
            var copy = new List<int>(internal_list.Length);
            for (var i = 0; i < internal_list.Length; i++)
                copy.Add(internal_list[i]);
            return copy;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes an action with access to a node's element list.
    /// </summary>
    /// <param name="node">The node index to look up.</param>
    /// <param name="action">Action to execute with the element list. Do NOT store references.</param>
    /// <remarks>
    ///     <para>
    ///         <b>NOTE:</b> This method allocates a copy due to the IReadOnlyList signature.
    ///         For true zero-copy access in hot paths, use <see cref="WithElementsForNodeSpan" />.
    ///     </para>
    ///     <para>
    ///         <b>SAFETY:</b> Read lock is held during execution - do not call other M2M methods.
    ///     </para>
    /// </remarks>
    public void WithElementsForNode(int node, Action<IReadOnlyList<int>> action)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        ArgumentNullException.ThrowIfNull(action);
        EnterSynchronizedReadLock();
        try
        {
            if ((uint)node >= (uint)_elementsFromNode!.Count)
            {
                action(Array.Empty<int>());
                return;
            }

            // Note: ToArray() is required here because the action expects IReadOnlyList<int>
            // and the internal O2M indexer returns ReadOnlySpan<int>.
            // For true zero-copy access, use the public WithElementsForNodeSpan method.
            action(_elementsFromNode[node].ToArray());
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes an action with true zero-copy span access to a node's element list.
    ///     P1-1 FIX: This is the recommended method for hot paths.
    /// </summary>
    /// <param name="node">The node index to look up.</param>
    /// <param name="action">Action to execute with the element span. Reference is only valid during execution.</param>
    /// <remarks>
    ///     <para>
    ///         <b>PERFORMANCE:</b> True zero-copy access. No allocations.
    ///     </para>
    ///     <para>
    ///         <b>SAFETY:</b> The span is only valid during action execution.
    ///         Do NOT store the span. Read lock is held during execution.
    ///     </para>
    /// </remarks>
    public void WithElementsForNodeSpan(int node, ReadOnlySpanAction<int> action)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        ArgumentNullException.ThrowIfNull(action);
        EnterSynchronizedReadLock();
        try
        {
            if ((uint)node >= (uint)_elementsFromNode!.Count)
            {
                action(ReadOnlySpan<int>.Empty);
                return;
            }

            action(_elementsFromNode[node]);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all elements that contain all specified nodes.
    /// </summary>
    /// <param name="nodes">The nodes that must all be present.</param>
    /// <returns>List of element indices containing all specified nodes.</returns>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe. Automatically upgrades to write lock for synchronization if needed.
    ///     Time Complexity: O(k × m) where k = nodes.Count, m = avg elements per node (after synchronization).
    ///     First query after modification: O(n × m) for sync + O(k × m) for query.
    ///     Space Complexity: O(result size)
    ///     Uses intersection strategy: starts with smallest candidate set, then filters.
    /// </remarks>
    public List<int> GetElementsWithNodes(List<int> nodes)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(nodes);
        if (nodes.Count == 0) return [];

        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            // Find the smallest element list for optimal intersection
            ReadOnlySpan<int> smallestList = default;
            var smallestLength = int.MaxValue;
            var smallestIdx = -1;

            for (var i = 0; i < nodes.Count; i++)
            {
                var node = nodes[i];
                if ((uint)node >= (uint)_elementsFromNode!.Count)
                    return []; // Node out of range

                var currentList = _elementsFromNode[node];
                if (smallestIdx == -1 || currentList.Length < smallestLength)
                {
                    smallestList = currentList;
                    smallestLength = currentList.Length;
                    smallestIdx = i;
                }
            }

            if (smallestLength == 0)
                return [];

            // Single node case - fast path
            if (nodes.Count == 1)
                return new List<int>(smallestList.ToArray());

            // Use HashSet for efficient intersection
            var resultSet = new HashSet<int>(smallestList.ToArray());

            // Intersect with the other lists
            for (var i = 0; i < nodes.Count; i++)
            {
                if (i == smallestIdx) continue; // Skip the one we already added

                var nodeElements = _elementsFromNode![nodes[i]];
                resultSet.IntersectWith(nodeElements.ToArray());

                if (resultSet.Count == 0)
                    return []; // Early exit
            }

            var result = new List<int>(resultSet);
            result.Sort();
            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all elements that contain at least one of the specified nodes (union semantics).
    /// </summary>
    /// <param name="nodes">The nodes to search for (at least one must be present).</param>
    /// <returns>List of element indices containing at least one of the specified nodes.</returns>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe. Automatically upgrades to write lock for synchronization if needed.
    ///     Time Complexity: O(k × m) where k = nodes.Count, m = avg elements per node.
    ///     First query after modification: O(n × m) for sync + O(k × m) for query.
    ///     <b>DIFFERENCE FROM GetElementsWithNodes:</b>
    ///     - GetElementsWithNodes returns elements containing ALL specified nodes (intersection)
    ///     - GetElementsContainingAnyNode returns elements containing AT LEAST ONE (union)
    /// </remarks>
    public List<int> GetElementsContainingAnyNode(List<int> nodes)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(nodes);
        if (nodes.Count == 0) return [];

        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            // Use HashSet for automatic deduplication (union semantics)
            var resultSet = new HashSet<int>();

            foreach (var node in nodes)
            {
                // Skip out-of-range nodes
                if ((uint)node >= (uint)_elementsFromNode!.Count)
                    continue;

                var elementsForNode = _elementsFromNode[node];
                foreach (var elem in elementsForNode)
                    resultSet.Add(elem);
            }

            var result = new List<int>(resultSet);
            result.Sort();
            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all elements that contain exactly the specified nodes (no more, no less).
    /// </summary>
    /// <param name="nodes">The exact set of nodes to match.</param>
    /// <returns>List of element indices with exact node match.</returns>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(k × m + r) where k = nodes.Count, m = avg elements per node, r = result size.
    ///     First query after modification: O(n × m) for sync + O(k × m + r) for query.
    ///     First finds all elements containing all specified nodes (superset),
    ///     then filters to only those with exact cardinality match.
    /// </remarks>
    public List<int> GetElementsFromNodes(List<int> nodes)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(nodes);

        var nodeCount = nodes.Count;
        var candidates = GetElementsWithNodes(nodes); // Already synchronized

        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var result = new List<int>(candidates.Count);

            foreach (var e in candidates)
                if (e < _o2m.Count && _o2m[e].Length == nodeCount)
                    result.Add(e);

            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all elements that share at least one node with the specified element.
    /// </summary>
    /// <param name="element">The element index to find neighbors for.</param>
    /// <param name="sorted">Whether to return neighbors in sorted order (default: true).</param>
    /// <returns>List of neighboring element indices.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when element is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(m × k) where m = nodes in element, k = avg elements per node.
    ///     Add O(n log n) if sorted=true where n = number of neighbors.
    ///     First query after modification: O(n × m) for sync + query cost.
    ///     Space Complexity: O(neighbors)
    ///     Returns elements that share at least one node, excluding the query element itself.
    /// </remarks>
    public List<int> GetElementNeighbors(int element, bool sorted = true)
    {
        ThrowIfDisposed();
        ArgumentOutOfRangeException.ThrowIfNegative(element);

        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            if (element >= _o2m.Count)
                throw new ArgumentOutOfRangeException(
                    nameof(element),
                    $"Element {element} is out of range. Structure has {_o2m.Count} elements.");

            var neighbours = new HashSet<int>();

            foreach (var node in _o2m[element])
            {
                // P0 FIX: Skip out-of-range nodes (consistent with other methods at lines 1099, 1135, etc.)
                if ((uint)node >= (uint)_elementsFromNode!.Count) continue;

                foreach (var elem in _elementsFromNode[node])
                    if (elem != element)
                        neighbours.Add(elem);
            }

            var result = new List<int>(neighbours);
            if (sorted) result.Sort();
            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all nodes that share at least one element with the specified node.
    /// </summary>
    /// <param name="node">The node index to find neighbors for.</param>
    /// <param name="sorted">Whether to return neighbors in sorted order (default: true).</param>
    /// <returns>List of neighboring node indices.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when node is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(k × m) where k = elements containing node, m = avg nodes per element.
    ///     Add O(n log n) if sorted=true where n = number of neighbors.
    ///     First query after modification: O(n × m) for sync + query cost.
    ///     Space Complexity: O(neighbors)
    ///     Returns nodes that share at least one element, excluding the query node itself.
    /// </remarks>
    public List<int> GetNodeNeighbors(int node, bool sorted = true)
    {
        ThrowIfDisposed();
        ArgumentOutOfRangeException.ThrowIfNegative(node);

        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            if (node >= _elementsFromNode!.Count)
                throw new ArgumentOutOfRangeException(
                    nameof(node),
                    $"Node {node} is out of range. Transpose has {_elementsFromNode.Count} nodes.");

            var neighbours = new HashSet<int>();

            foreach (var elem in _elementsFromNode[node])
            foreach (var n in _o2m[elem])
                if (n != node)
                    neighbours.Add(n);

            var result = new List<int>(neighbours);
            if (sorted) result.Sort();
            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the maximum node index across all elements.
    /// </summary>
    /// <returns>The maximum node index, or -1 if structure is empty.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(1) if O2M has cached the value, O(n × m) otherwise.
    /// </remarks>
    public int GetMaxNode()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.GetMaxNode();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if the structure contains the specified element index.
    /// </summary>
    /// <param name="elementIndex">The element index to check.</param>
    /// <returns>True if element exists; false otherwise.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(1)
    /// </remarks>
    public bool HasElement(int elementIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return elementIndex >= 0 && elementIndex < _o2m.Count;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if the specified node is connected to any element.
    /// </summary>
    /// <param name="nodeIndex">The node index to check.</param>
    /// <returns>True if node is connected to at least one element; false otherwise.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × m) in worst case if not synchronized; O(1) if synchronized.
    /// </remarks>
    public bool HasNode(int nodeIndex)
    {
        ThrowIfDisposed();
        if (nodeIndex < 0) return false;

        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return nodeIndex < _elementsFromNode!.Count && _elementsFromNode[nodeIndex].Length > 0;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if the specified element contains the specified node.
    /// </summary>
    /// <param name="elementIndex">The element index.</param>
    /// <param name="nodeIndex">The node index.</param>
    /// <returns>True if element contains node; false otherwise.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(m) where m = nodes in element.
    /// </remarks>
    public bool ElementContainsNode(int elementIndex, int nodeIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            if (elementIndex < 0 || elementIndex >= _o2m.Count) return false;
            return _o2m[elementIndex].Contains(nodeIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Validates that all node indices are non-negative and unique within each element.
    /// </summary>
    /// <returns>True if structure is valid; false otherwise.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × m) where n = elements, m = average nodes per element.
    /// </remarks>
    public bool IsValid()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.IsValid();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the element-to-element connectivity matrix.
    /// </summary>
    /// <returns>O2M representing which elements connect to which elements through shared nodes.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × k × m) where n = elements, k = avg nodes per element,
    ///     m = avg elements per node. First query after modification adds O(n × m) for sync.
    ///     Computes A × A^T where A is the element-to-node adjacency matrix.
    ///     Two elements are connected if they share at least one common node.
    ///     Useful for finite element mesh connectivity analysis.
    /// </remarks>
    public O2M GetElementsToElements()
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return _o2m * _elementsFromNode!;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the element-to-element connectivity as an M2M with caching.
    /// </summary>
    /// <returns>M2M representing element-to-element adjacency.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Similar to GetElementsToElements() but returns M2M for additional query capabilities.
    ///     Useful for graph algorithms on the element connectivity graph.
    /// </remarks>
    public M2M GetElementToElementGraph()
    {
        // P0.3 FIX: Use ownership constructor — e2e is a fresh O2M from multiplication.
        var e2e = GetElementsToElements();
        return new M2M(e2e, takeOwnership: true);
    }

    /// <summary>
    ///     Gets the node-to-node connectivity matrix.
    /// </summary>
    /// <returns>O2M representing which nodes connect to which nodes through shared elements.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × k × m) where n = nodes, k = avg elements per node,
    ///     m = avg nodes per element. First query after modification adds O(n × m) for sync.
    ///     Computes A^T × A where A is the element-to-node adjacency matrix.
    ///     Two nodes are connected if they belong to at least one common element.
    ///     Useful for graph algorithms on the dual mesh.
    /// </remarks>
    public O2M GetNodesToNodes()
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return _elementsFromNode! * _o2m;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the node-to-node connectivity as an M2M with caching.
    /// </summary>
    /// <returns>M2M representing node-to-node adjacency.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Similar to GetNodesToNodes() but returns M2M for additional query capabilities.
    ///     Useful for graph algorithms on node connectivity.
    /// </remarks>
    public M2M GetNodeToNodeGraph()
    {
        // P0.3 FIX: Use ownership constructor — n2n is a fresh O2M from multiplication.
        var n2n = GetNodesToNodes();
        return new M2M(n2n, takeOwnership: true);
    }

    /// <summary>
    ///     Computes a topological ordering of elements based on their dependencies.
    /// </summary>
    /// <returns>List of element indices in topological order.</returns>
    /// <exception cref="InvalidOperationException">Thrown when the graph contains cycles.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(V + E) where V = elements, E = total connections.
    ///     Uses Kahn's algorithm with parallel atomic in-degree calculation.
    ///     Very useful for determining optimal processing order in FEM assembly
    ///     or any dependency-based computation.
    ///     Example: If element A depends on element B (A's nodes include B as an index),
    ///     then B will appear before A in the topological order.
    /// </remarks>
    public List<int> GetTopologicalOrder()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.GetTopOrder();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Tests if the structure represents an acyclic graph.
    /// </summary>
    /// <returns>True if the graph has no cycles; false otherwise.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(V + E) where V = elements, E = total connections.
    ///     Uses DFS with three-coloring (white/gray/black) for cycle detection.
    ///     Essential for validating dependency graphs before topological sorting.
    /// </remarks>
    public bool IsAcyclic()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.IsAcyclic();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Computes a sort order for elements based on lexicographic comparison of their adjacency lists.
    /// </summary>
    /// <returns>List of element indices in sorted order.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n log n × m) where n = elements, m = avg nodes per element.
    ///     Returns a permutation such that: this[order[i]] ≤ this[order[i+1]] lexicographically.
    ///     Useful for creating canonical representations or deterministic output ordering.
    /// </remarks>
    public List<int> GetSortOrder()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.GetSortOrder();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Finds indices of duplicate elements (elements with identical adjacency lists).
    /// </summary>
    /// <returns>List of element indices that are duplicates of earlier elements.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n log n × m) where n = elements, m = avg nodes per element.
    ///     Elements are first sorted, then consecutive equal elements are identified.
    ///     Useful for mesh cleanup after refinement or merging operations.
    ///     The returned indices can be used with CompressElements to remove duplicates.
    ///     Example usage:
    ///     <code>
    /// var duplicates = m2m.GetDuplicates();
    /// if (duplicates.Count > 0)
    /// {
    ///     var toKeep = Enumerable.Range(0, m2m.Count).Except(duplicates).ToList();
    ///     m2m.CompressElements(toKeep);
    /// }
    /// </code>
    /// </remarks>
    public List<int> GetDuplicates()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.GetDuplicates();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Tests if this M2M is a permutation of another M2M.
    /// </summary>
    /// <param name="other">The M2M to compare with.</param>
    /// <returns>True if the structures contain the same elements in possibly different order.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if either object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock on both instances.
    ///     Time Complexity: O(n log n × m) where n = elements, m = avg nodes per element.
    ///     Two structures are permutations if they have the same elements (same adjacency lists)
    ///     but possibly in different order. Useful for validating that reordering operations
    ///     didn't lose or corrupt data.
    /// </remarks>
    public bool IsPermutationOf(M2M? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;

        ThrowIfDisposed();
        other.ThrowIfDisposed();

        var (first, second) = GetLockOrder(this, other);

        first._rwLock.EnterReadLock();
        try
        {
            second._rwLock.EnterReadLock();
            try
            {
                return _o2m.IsPermutationOf(other._o2m);
            }
            finally
            {
                second._rwLock.ExitReadLock();
            }
        }
        finally
        {
            first._rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Reserves capacity for the specified number of elements to avoid reallocation.
    /// </summary>
    /// <param name="capacity">The minimum capacity to ensure.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when capacity is negative.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with write lock.
    ///     Time Complexity: O(n) in worst case if reallocation occurs.
    ///     Use this before bulk operations to avoid multiple reallocations.
    /// </remarks>
    public void Reserve(int capacity)
    {
        ThrowIfDisposed();
        ArgumentOutOfRangeException.ThrowIfNegative(capacity);
        _rwLock.EnterWriteLock();
        try
        {
            _o2m.Reserve(capacity);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Reduces memory usage by trimming excess capacity.
    /// </summary>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with write lock.
    ///     Time Complexity: O(n × m) where n = elements, m = avg nodes per element.
    ///     Call after bulk operations complete to reclaim unused memory.
    /// </remarks>
    public void ShrinkToFit()
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            _o2m.ShrinkToFit();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Batch Operations

    /// <summary>
    ///     Begins a batch update scope that defers cache synchronization.
    /// </summary>
    /// <returns>A disposable batch scope that synchronizes on disposal.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>IMPORTANT:</b> While inside a batch scope, queries that depend on cached data
    ///     (ElementsFromNode, ElementLocations, NodeLocations, HasNode, GetElementsWithNodes,
    ///     GetElementNeighbors, GetNodeNeighbors) will throw <see cref="InvalidOperationException" />.
    ///     This prevents reading stale data during modifications.
    ///     <b>Thread Safety:</b> Thread-safe, but batch state is SHARED across all threads.
    ///     When any thread enters a batch scope, cache-dependent queries will throw for
    ///     ALL threads until the outermost batch scope exits. Batch scopes can be nested,
    ///     but nesting is global to the M2M instance, not per-thread.
    ///     Use batch operations when performing multiple modifications to avoid
    ///     repeated cache invalidation and resynchronization:
    ///     <code>
    /// using (m2m.BeginBatchUpdate())
    /// {
    ///     for (int i = 0; i &lt; 1000; i++)
    ///     {
    ///         m2m.AppendElement([...]);
    ///     }
    ///     // DO NOT query ElementsFromNode or similar here - will throw!
    ///     // This applies to ALL threads, not just this thread.
    /// } // Cache synchronization occurs on next query after batch ends
    /// // Safe to query after batch scope exits (if no other thread holds a batch)
    /// var neighbors = m2m.GetElementNeighbors(0);
    ///     </code>
    ///     Batch scopes can be nested. Synchronization only occurs when the outermost scope exits
    ///     and a query requiring cached data is performed.
    /// </remarks>
    public IDisposable BeginBatchUpdate()
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            _batchNesting++;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }

        return new BatchUpdateScope(this);
    }

    /// <summary>
    ///     Forces synchronization of the internal caches.
    /// </summary>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with write lock.
    ///     Time Complexity: O(n × m) for transpose only.
    ///     Position caches are computed lazily when accessed.
    ///     Normally called automatically by query methods. You can call this manually
    ///     after batch updates if you want to control when synchronization occurs.
    /// </remarks>
    public void Synchronize()
    {
        ThrowIfDisposed();
        EnsureSynchronized();
    }

    /// <summary>
    ///     Forces synchronization of the transpose and position caches.
    /// </summary>
    /// <exception cref="InvalidOperationException">Thrown when called during a batch update.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with write lock.
    ///     Time Complexity: O(n × m) for transpose + position cache computation.
    ///     Use this method to pre-compute all caches before running hot-path queries
    ///     that rely on <see cref="ElementLocations" /> or <see cref="NodeLocations" />.
    /// </remarks>
    public void EnsurePositionCaches()
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        _rwLock.EnterWriteLock();
        try
        {
            if (!_isInSync)
                SynchronizeTranspose();
            if (!_positionCachesComputed)
                ComputePositionCaches();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Internal batch update scope that manages nesting counter.
    /// </summary>
    /// <remarks>
    ///     <b>THREAD SAFETY:</b> Must be disposed on the same thread that created it.
    ///     This is required because the underlying ReaderWriterLockSlim tracks thread affinity.
    /// </remarks>
    private sealed class BatchUpdateScope : IDisposable
    {
        private readonly M2M _owner;
        private readonly int _creatingThreadId;
        private bool _disposed;

        public BatchUpdateScope(M2M owner)
        {
            _owner = owner;
            _creatingThreadId = Environment.CurrentManagedThreadId;
        }

        public void Dispose()
        {
            if (_disposed) return;
            _disposed = true;

            // Don't throw if owner is disposed - just exit gracefully
            if (Volatile.Read(ref _owner._disposed) != 0) return;

            // Validate thread affinity - ReaderWriterLockSlim requires same-thread access
            var currentThreadId = Environment.CurrentManagedThreadId;
            if (currentThreadId != _creatingThreadId)
                throw new InvalidOperationException(
                    $"BatchUpdateScope must be disposed on the same thread that created it. " +
                    $"Created on thread {_creatingThreadId}, but Dispose called on thread {currentThreadId}. " +
                    $"Ensure 'using' statements do not cross async boundaries without ConfigureAwait(true).");

            _owner._rwLock.EnterWriteLock();
            try
            {
                _owner._batchNesting--;
                Debug.Assert(_owner._batchNesting >= 0, "Batch nesting should never be negative");
            }
            finally
            {
                _owner._rwLock.ExitWriteLock();
            }
        }
    }

    #endregion

    #region Synchronization Logic

    /// <summary>
    ///     Ensures the transpose cache is synchronized, upgrading from read lock to write lock if needed.
    /// </summary>
    /// <remarks>
    ///     This method implements the double-check locking pattern:
    ///     1. Check if sync needed under read lock (fast path)
    ///     2. If needed, upgrade to write lock
    ///     3. Check again (another thread might have synchronized)
    ///     4. Perform synchronization if still needed
    ///     Synchronization is deferred while inside a batch update scope (_batchNesting > 0).
    /// </remarks>
    private void EnsureSynchronized()
    {
        // Fast path: check under read lock
        _rwLock.EnterReadLock();
        try
        {
            if (_isInSync || _batchNesting > 0)
                return; // Already synchronized or batching
        }
        finally
        {
            _rwLock.ExitReadLock();
        }

        // Slow path: need to synchronize
        _rwLock.EnterWriteLock();
        try
        {
            // Double-check: another thread might have synchronized while we waited for write lock
            // Also check batch nesting again in case a batch started
            if (_isInSync || _batchNesting > 0)
                return;

            SynchronizeTranspose();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Enters a read lock with the transpose cache guaranteed to be synchronized.
    ///     Returns with the read lock HELD — caller MUST release via _rwLock.ExitReadLock().
    /// </summary>
    /// <remarks>
    ///     <b>P0.3 FIX:</b> Eliminates the TOCTOU race in the previous pattern where
    ///     <c>EnsureSynchronized()</c> released all locks, then callers re-acquired a read lock.
    ///     A writer could interpose between those two steps and null out <c>_elementsFromNode</c>.
    ///     <para/>
    ///     This method uses the RWLS downgrade pattern: if synchronization is needed, it acquires
    ///     the write lock, synchronizes, then acquires a read lock <i>while still holding the
    ///     write lock</i> (allowed by <see cref="LockRecursionPolicy.SupportsRecursion"/>),
    ///     and finally releases the write lock. This leaves only the read lock held, with no
    ///     window for a writer to invalidate the cache.
    /// </remarks>
    private void EnterSynchronizedReadLock()
    {
        // Fast path: cache is already valid
        _rwLock.EnterReadLock();
        if (_isInSync)
            return; // Read lock held, _elementsFromNode is guaranteed non-null
        _rwLock.ExitReadLock();

        // Slow path: synchronize under write lock, then downgrade to read lock
        _rwLock.EnterWriteLock();
        try
        {
            // Double-check: another thread may have synchronized while we waited
            if (!_isInSync && _batchNesting == 0)
                SynchronizeTranspose();
        }
        finally
        {
            // Downgrade: acquire read lock while holding write lock, then release write lock.
            // SupportsRecursion allows EnterReadLock while the write lock is held.
            // After ExitWriteLock, only the read lock remains — no gap for a writer.
            _rwLock.EnterReadLock();
            _rwLock.ExitWriteLock();
        }
        // Returns with read lock held
    }

    /// <summary>
    ///     Enters a read lock with both the transpose and position caches guaranteed valid.
    ///     Returns with the read lock HELD — caller MUST release via _rwLock.ExitReadLock().
    /// </summary>
    /// <remarks>
    ///     <b>P0.3 FIX:</b> Same TOCTOU fix as <see cref="EnterSynchronizedReadLock"/> but
    ///     additionally ensures position caches (<c>_elemeloc</c>, <c>_nodeloc</c>) are computed.
    /// </remarks>
    private void EnterPositionCachedReadLock()
    {
        // Fast path: everything already computed
        _rwLock.EnterReadLock();
        if (_isInSync && _positionCachesComputed)
            return; // Read lock held, all caches valid
        _rwLock.ExitReadLock();

        // Slow path: compute under write lock, then downgrade
        _rwLock.EnterWriteLock();
        try
        {
            if (!_isInSync && _batchNesting == 0)
                SynchronizeTranspose();
            // FIX: Guard against calling ComputePositionCaches when transpose is not synchronized.
            // ComputePositionCaches requires _isInSync (asserted via Debug.Assert).
            // If we're in a batch (_batchNesting > 0), we cannot synchronize, so we must skip
            // position cache computation. The caller will get the read lock but with stale caches,
            // which is acceptable since batch operations should not rely on position caches.
            if (_isInSync && !_positionCachesComputed)
                ComputePositionCaches();
        }
        finally
        {
            _rwLock.EnterReadLock();
            _rwLock.ExitWriteLock();
        }
        // Returns with read lock held
    }

    /// <summary>
    ///     Internal synchronization logic that rebuilds transpose cache only.
    ///     Must be called within write lock.
    /// </summary>
    private void SynchronizeTranspose()
    {
        Debug.Assert(_rwLock.IsWriteLockHeld, "SynchronizeTranspose must be called under write lock");

        // Compute transpose (output is sorted, required by position calculations if needed later)
        _elementsFromNode = _o2m.Transpose();

        // _isInSync is declared volatile, which provides the necessary release semantics:
        // the _elementsFromNode assignment above is guaranteed visible to other threads
        // before _isInSync is observed as true in the fast path of EnterSynchronizedReadLock.
        _isInSync = true;

        // Mark position caches as invalid since structure changed
        _positionCachesComputed = false;
        _elemeloc = null;
        _nodeloc = null;
    }

    /// <summary>
    ///     Computes position caches. Must be called within write lock with transpose already computed.
    /// </summary>
    private void ComputePositionCaches()
    {
        Debug.Assert(_rwLock.IsWriteLockHeld, "ComputePositionCaches must be called under write lock");
        Debug.Assert(_isInSync, "Transpose must be synchronized before computing position caches");

        // Compute position caches
        var elemPositions = O2M.GetElementPositions(_o2m, _elementsFromNode!);
        var nodePositions = O2M.GetNodePositions(_o2m, _elementsFromNode!);

        // Convert to read-only for public properties
        var elemLocList = new List<IReadOnlyList<int>>(elemPositions.Count);
        for (var i = 0; i < elemPositions.Count; i++)
            elemLocList.Add(elemPositions[i].AsReadOnly());
        _elemeloc = elemLocList.AsReadOnly();

        var nodeLocList = new List<IReadOnlyList<int>>(nodePositions.Count);
        for (var i = 0; i < nodePositions.Count; i++)
            nodeLocList.Add(nodePositions[i].AsReadOnly());
        _nodeloc = nodeLocList.AsReadOnly();

        _positionCachesComputed = true;
    }

    /// <summary>
    ///     Marks all caches as invalid and releases cached memory.
    ///     Must be called within write lock.
    /// </summary>
    /// <remarks>
    ///     <b>FIX (Priority 1):</b> Now sets cached structures to null to release memory.
    ///     Previously, invalidation only set flags but retained references to large
    ///     cached objects (_elementsFromNode, _elemeloc, _nodeloc), causing unnecessary
    ///     memory retention for large meshes. Setting these to null allows the GC to
    ///     reclaim the memory immediately when the cache is invalidated.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void InvalidateCache()
    {
        Debug.Assert(_rwLock.IsWriteLockHeld, "InvalidateCache must be called under write lock");
        _isInSync = false;
        _positionCachesComputed = false;

        // FIX: Release references to potentially large cached structures
        // This allows GC to reclaim memory for large meshes
        _elementsFromNode = null;
        _elemeloc = null;
        _nodeloc = null;
    }

    /// <summary>
    ///     Throws ObjectDisposedException if this instance has been disposed.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ThrowIfDisposed()
    {
        if (Volatile.Read(ref _disposed) != 0)
        {
            if (_disposalIncomplete)
                throw new ObjectDisposedException(nameof(M2M),
                    "Object disposal was attempted but did not complete successfully. " +
                    "The lock could not be acquired within timeout, indicating concurrent access. " +
                    "Check IsDisposalIncomplete property and retry Dispose() after ensuring all threads have finished.");
            throw new ObjectDisposedException(nameof(M2M));
        }
    }

    /// <summary>
    ///     Throws InvalidOperationException if currently inside a batch update.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ThrowIfInBatch()
    {
        _rwLock.EnterReadLock();
        try
        {
            if (_batchNesting > 0)
                throw new InvalidOperationException(
                    "Cannot read M2M properties while inside a BatchUpdate. " +
                    "Exit the batch scope first or call Synchronize() after the batch.");
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Cloning

    /// <summary>
    ///     Creates a deep copy of this M2M instance.
    /// </summary>
    /// <returns>New M2M instance with copied data and synchronized cache if available.</returns>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × m) for structure copy.
    ///     If the source is synchronized, the clone will also be synchronized with copied caches.
    ///     This avoids resynchronization cost on first query to the cloned instance.
    /// </remarks>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public object Clone()
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        _rwLock.EnterReadLock();
        try
        {
            // P0.3 FIX: Use ownership constructor to avoid double-clone.
            // M2M(O2M) would clone again; the O2M is already a fresh deep copy.
            var clonedO2m = (O2M)_o2m.Clone();
            var cloned = new M2M(clonedO2m, takeOwnership: true);

            // Copy transpose cache if synchronized
            if (_isInSync)
            {
                cloned._elementsFromNode = (O2M)_elementsFromNode!.Clone();
                cloned._isInSync = true;
            }

            // Copy position caches if computed
            if (_positionCachesComputed)
            {
                cloned._elemeloc = _elemeloc; // ReadOnlyCollection, safe to share
                cloned._nodeloc = _nodeloc; // ReadOnlyCollection, safe to share
                cloned._positionCachesComputed = true;
            }

            return cloned;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Creates a strongly-typed deep copy of this M2M instance.
    /// </summary>
    /// <returns>New M2M instance with copied data.</returns>
    public M2M CloneTyped()
    {
        return (M2M)Clone();
    }

    #endregion

    #region Comparison and Equality

    /// <summary>
    ///     Compares this instance to another M2M lexicographically.
    /// </summary>
    /// <param name="other">The M2M to compare with.</param>
    /// <returns>
    ///     - Negative if this &lt; other
    ///     - Zero if this == other
    ///     - Positive if this &gt; other
    /// </returns>
    /// <exception cref="ObjectDisposedException">Thrown if either object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × m) where n = elements, m = average nodes per element.
    ///     Lock ordering uses object identity hash codes to prevent deadlocks when
    ///     comparing two M2M instances from different threads.
    /// </remarks>
    public int CompareTo(M2M? other)
    {
        if (other is null) return 1;
        if (ReferenceEquals(this, other)) return 0;

        ThrowIfDisposed();
        other.ThrowIfDisposed();

        // Acquire locks in consistent order to prevent deadlocks
        var (first, second) = GetLockOrder(this, other);

        first._rwLock.EnterReadLock();
        try
        {
            second._rwLock.EnterReadLock();
            try
            {
                return _o2m.CompareTo(other._o2m);
            }
            finally
            {
                second._rwLock.ExitReadLock();
            }
        }
        finally
        {
            first._rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Determines equality with another M2M.
    /// </summary>
    /// <param name="other">The M2M to compare with.</param>
    /// <returns>True if structures are equal; false otherwise.</returns>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × m) where n = elements, m = average nodes per element.
    ///     Lock ordering uses object identity hash codes to prevent deadlocks when
    ///     comparing two M2M instances from different threads.
    /// </remarks>
    public bool Equals(M2M? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;

        // Don't throw on disposed - just return false for safety
        if (Volatile.Read(ref _disposed) != 0 || Volatile.Read(ref other._disposed) != 0) return false;

        // Acquire locks in consistent order to prevent deadlocks
        var (first, second) = GetLockOrder(this, other);

        first._rwLock.EnterReadLock();
        try
        {
            second._rwLock.EnterReadLock();
            try
            {
                return _o2m.Equals(other._o2m);
            }
            finally
            {
                second._rwLock.ExitReadLock();
            }
        }
        finally
        {
            first._rwLock.ExitReadLock();
        }
    }

    public override bool Equals(object? obj)
    {
        return Equals(obj as M2M);
    }

    public override int GetHashCode()
    {
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            return _o2m.GetHashCode();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Determines consistent lock ordering to prevent deadlocks.
    /// </summary>
    private static (M2M first, M2M second) GetLockOrder(M2M a, M2M b)
    {
        // Use RuntimeHelpers.GetHashCode for consistent ordering based on object identity
        var aId = RuntimeHelpers.GetHashCode(a);
        var bId = RuntimeHelpers.GetHashCode(b);

        if (aId != bId)
            return aId < bId ? (a, b) : (b, a);

        // Extremely rare: same identity hash for different objects
        // Use internal O2M hash as tiebreaker
        var aO2mId = RuntimeHelpers.GetHashCode(a._o2m);
        var bO2mId = RuntimeHelpers.GetHashCode(b._o2m);

        if (aO2mId != bO2mId)
            return aO2mId < bO2mId ? (a, b) : (b, a);

        // Still tied — use unique monotonic instance ID as final tiebreaker.
        // Unlike Unsafe.IsAddressLessThan on stack parameters (which compares local
        // variable addresses, not heap object addresses), _instanceId is guaranteed
        // unique per instance and consistent regardless of call-site parameter order.
        return a._instanceId < b._instanceId ? (a, b) : (b, a);
    }

    public static bool operator ==(M2M? left, M2M? right)
    {
        return ReferenceEquals(left, right) || left?.Equals(right) == true;
    }

    public static bool operator !=(M2M? left, M2M? right)
    {
        return !(left == right);
    }

    public static bool operator <(M2M? left, M2M? right)
    {
        return left is null ? right is not null : left.CompareTo(right) < 0;
    }

    public static bool operator >(M2M? left, M2M? right)
    {
        return right < left;
    }

    public static bool operator <=(M2M? left, M2M? right)
    {
        return !(left > right);
    }

    public static bool operator >=(M2M? left, M2M? right)
    {
        return !(left < right);
    }

    #endregion

    #region Static Factory Methods

    /// <summary>
    ///     Creates a new M2M instance from CSR (Compressed Sparse Row) format.
    /// </summary>
    /// <param name="rowPointers">Row pointer array (length = rows + 1).</param>
    /// <param name="columnIndices">Column indices array.</param>
    /// <returns>New M2M instance.</returns>
    /// <exception cref="ArgumentNullException">Thrown when arrays are null.</exception>
    /// <remarks>
    ///     Time Complexity: O(nnz) where nnz = number of non-zeros.
    /// </remarks>
    public static M2M FromCsr(int[] rowPointers, int[] columnIndices)
    {
        var o2m = O2M.FromCsr(rowPointers, columnIndices);
        return new M2M(o2m, takeOwnership: true);
    }

    /// <summary>
    ///     Creates a new M2M instance from a boolean matrix.
    /// </summary>
    /// <param name="matrix">Boolean matrix where true indicates a connection.</param>
    /// <returns>New M2M instance.</returns>
    /// <exception cref="ArgumentNullException">Thrown when matrix is null.</exception>
    /// <remarks>
    ///     Time Complexity: O(r × c) where r = rows, c = columns.
    /// </remarks>
    public static M2M FromBooleanMatrix(bool[,] matrix)
    {
        var o2m = O2M.FromBooleanMatrix(matrix);
        return new M2M(o2m, takeOwnership: true);
    }

    /// <summary>
    ///     Generates a random M2M for testing and benchmarking.
    /// </summary>
    /// <param name="elementCount">Number of elements to generate.</param>
    /// <param name="nodeCount">Range of node indices [0, nodeCount).</param>
    /// <param name="density">Probability of including each node (0.0 to 1.0).</param>
    /// <param name="seed">Optional random seed for reproducibility.</param>
    /// <returns>Randomly generated M2M instance.</returns>
    /// <remarks>
    ///     Time Complexity: O(elementCount × nodeCount × density).
    /// </remarks>
    public static M2M CreateRandom(int elementCount, int nodeCount, double density, int? seed = null)
    {
        var o2m = O2M.GetRandomO2M(elementCount, nodeCount, density, seed);
        return new M2M(o2m, takeOwnership: true);
    }

    #endregion

    #region Conversion Methods

    /// <summary>
    ///     Converts to Compressed Sparse Row (CSR) format.
    /// </summary>
    /// <returns>Tuple of (row pointers, column indices) arrays.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × m).
    /// </remarks>
    public (int[] rowPtr, int[] columnIndices) ToCsr()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.ToCsr();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Converts to a dense boolean matrix.
    /// </summary>
    /// <returns>Boolean matrix where true indicates a connection.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × (m + k)) where k = max node index.
    ///     WARNING: Can consume large amounts of memory for sparse structures.
    /// </remarks>
    public bool[,] ToBooleanMatrix()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _o2m.ToBooleanMatrix();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Generates a string representation.
    /// </summary>
    /// <returns>String showing element indices and their connected nodes.</returns>
    /// <remarks>
    ///     Thread Safety: Thread-safe with read lock.
    ///     Time Complexity: O(n × m).
    /// </remarks>
    public override string ToString()
    {
        // Return a safe string for disposed instances instead of throwing.
        // The read lock may already be disposed, so we cannot acquire it.
        if (Volatile.Read(ref _disposed) != 0) return "[Disposed M2M]";

        _rwLock.EnterReadLock();
        try
        {
            return _o2m.ToString();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Unsafe Method Protection (Conditional Verification)

    /// <summary>
    ///     Verifies that a lock is held. Compiled only in DEBUG builds.
    /// </summary>
    [Conditional("DEBUG")]
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void VerifyLockHeld()
    {
        if (_rwLock.CurrentReadCount == 0 && !_rwLock.IsWriteLockHeld)
            throw new InvalidOperationException(
                "Unsafe method called without holding lock! " +
                "Call from within EnterReadLock or EnterWriteLock block. " +
                "This indicates a serious bug in the calling code that could cause race conditions.");
    }

    /// <summary>
    ///     Verifies that the cache is synchronized. Compiled only in DEBUG builds.
    /// </summary>
    [Conditional("DEBUG")]
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void VerifySynchronized()
    {
        if (!_isInSync)
            throw new InvalidOperationException(
                "Cache not synchronized! " +
                "Call EnsureSynchronized() before accessing cached data. " +
                "This indicates incorrect usage that could return stale data.");
    }

    /// <summary>
    ///     Verifies that cache access is safe (not in batch mode with stale cache).
    ///     This check runs in ALL builds to prevent silent data corruption.
    /// </summary>
    /// <remarks>
    ///     <b>CRITICAL:</b> During batch operations, the transpose cache may be stale.
    ///     Accessing stale cache data can return incorrect adjacency information,
    ///     leading to silent corruption in downstream algorithms.
    ///     This check is unconditional because the consequences of stale data are severe.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ThrowIfCacheUnsafe()
    {
        // If we're in a batch (_batchNesting > 0) AND the cache is not in sync,
        // the caller is attempting to read potentially stale data.
        // This is always an error, even in Release builds.
        if (_batchNesting > 0 && !_isInSync)
            throw new InvalidOperationException(
                "Cannot access cached adjacency data during batch operation after structural mutation. " +
                "The transpose cache may contain stale data. Either: " +
                "(1) Complete the batch operation first, or " +
                "(2) Use the safe accessor methods that handle synchronization automatically.");
    }

    #endregion


    #region Internal Non-Cloning Accessors (Performance Critical)

    /// <summary>
    ///     Delegate for actions that receive a read-only span.
    /// </summary>
    public delegate void ReadOnlySpanAction<T>(ReadOnlySpan<T> span);

    /// <summary>
    ///     Gets elements for a node as a snapshot copy. INTERNAL USE ONLY.
    /// </summary>
    /// <param name="node">The node index.</param>
    /// <returns>Copy of internal list. Safe to use after lock release.</returns>
    /// <remarks>
    ///     <b>P1-1 FIX:</b> Renamed from GetElementsForNodeUnsafe. This method returns a COPY.
    ///     The "Unsafe" suffix refers to the requirement to hold a lock, not zero-copy semantics.
    ///     <b>PERFORMANCE:</b> O(n) where n = number of elements for this node.
    ///     For true zero-copy access, use <see cref="WithElementsForNodeSpan" />.
    ///     <b>REQUIREMENTS:</b>
    ///     - Caller MUST hold appropriate lock (_rwLock.EnterReadLock or EnterWriteLock)
    ///     - Caller MUST ensure cache is synchronized
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal int[] GetElementsForNodeSnapshot(int node)
    {
        VerifyLockHeld();
        VerifySynchronized();
        ThrowIfCacheUnsafe();

        if ((uint)node >= (uint)_elementsFromNode!.Count)
            return Array.Empty<int>();

        return _elementsFromNode[node].ToArray();
    }

    /// <summary>
    ///     Gets elements for a node without cloning. INTERNAL USE ONLY.
    /// </summary>
    [Obsolete("P1-1: Use GetElementsForNodeSnapshot (returns copy) or WithElementsForNodeSpan (true zero-copy).")]
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal IReadOnlyList<int> GetElementsForNodeUnsafe(int node)
    {
        return GetElementsForNodeSnapshot(node);
    }

    /// <summary>
    ///     Gets nodes for an element as a snapshot copy. INTERNAL USE ONLY.
    /// </summary>
    /// <param name="element">The element index.</param>
    /// <returns>Copy of internal list. Safe to use after lock release.</returns>
    /// <remarks>
    ///     <b>P1-1 FIX:</b> Renamed from GetNodesForElementUnsafe. This method returns a COPY.
    ///     <b>PERFORMANCE:</b> O(n) where n = number of nodes for this element.
    ///     For true zero-copy access, use <see cref="WithNodesForElementSpan" />.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal int[] GetNodesForElementSnapshot(int element)
    {
        VerifyLockHeld();

        if ((uint)element >= (uint)_o2m.Count)
            return Array.Empty<int>();

        return _o2m[element].ToArray();
    }

    /// <summary>
    ///     Gets nodes for an element without cloning. INTERNAL USE ONLY.
    /// </summary>
    [Obsolete("P1-1: Use GetNodesForElementSnapshot (returns copy) or WithNodesForElementSpan (true zero-copy).")]
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal IReadOnlyList<int> GetNodesForElementUnsafe(int element)
    {
        return GetNodesForElementSnapshot(element);
    }

    /// <summary>
    ///     Executes an action with true zero-copy span access to an element's node list.
    /// </summary>
    /// <param name="element">The element index to look up.</param>
    /// <param name="action">Action to execute with the node span.</param>
    /// <remarks>
    ///     <b>P1-1 FIX:</b> True zero-copy access. No allocations.
    ///     See <see cref="WithElementsForNodeSpan" /> for safety requirements.
    /// </remarks>
    internal void WithNodesForElementSpan(int element, ReadOnlySpanAction<int> action)
    {
        VerifyLockHeld();

        if ((uint)element >= (uint)_o2m.Count)
        {
            action(ReadOnlySpan<int>.Empty);
            return;
        }

        // O2M indexer already returns ReadOnlySpan<int> - true zero-copy
        action(_o2m[element]);
    }

    #endregion

    #region IDisposable

    /// <summary>
    ///     Finalizer for M2M.
    /// </summary>
    /// <remarks>
    ///     NOTE: This finalizer does NOT dispose the ReaderWriterLockSlim (managed resources
    ///     are not disposed from finalizers per standard dispose pattern). The finalizer only
    ///     sets the disposed flag. The ReaderWriterLockSlim will be cleaned up by its own
    ///     finalizer eventually, but for deterministic cleanup, call Dispose() explicitly.
    /// </remarks>
    ~M2M()
    {
        Dispose(false);
    }

    /// <summary>
    ///     Releases resources used by the M2M instance.
    /// </summary>
    /// <remarks>
    ///     <b>EXPLICIT DISPOSAL RECOMMENDED:</b> Call Dispose() or use a using statement
    ///     for deterministic cleanup of the ReaderWriterLockSlim.
    ///     If not explicitly disposed, the ReaderWriterLockSlim will eventually be
    ///     cleaned up by its own finalizer during garbage collection, but this is
    ///     non-deterministic and should not be relied upon.
    ///     <code>
    /// // Recommended (deterministic cleanup):
    /// using (var m2m = new M2M())
    /// {
    ///     // ... use m2m ...
    /// }
    ///     </code>
    ///     After disposal, the M2M instance should not be used.
    /// </remarks>
    public void Dispose()
    {
        Dispose(true);
        GC.SuppressFinalize(this);
    }

    /// <summary>
    ///     Protected implementation of Dispose pattern.
    /// </summary>
    /// <param name="disposing">True if called from Dispose(), false if called from finalizer.</param>
    private void Dispose(bool disposing)
    {
        // First-time disposal: atomically mark as disposed
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0)
        {
            // Already marked disposed. Allow retry only if previous disposal was incomplete.
            // This handles the case where the write lock timed out on the first attempt.
            if (!Volatile.Read(ref _disposalIncomplete))
                return;
        }

        if (disposing)
        {
            // CRITICAL: Acquire write lock to wait for all readers to finish.
            // This ensures no thread is holding the lock when we dispose it.
            var lockAcquired = false;
            try
            {
                // Try to acquire with timeout to avoid deadlock if something went wrong
                lockAcquired = _rwLock.TryEnterWriteLock(TimeSpan.FromSeconds(30));

                if (!lockAcquired)
                {
                    // Mark disposal as incomplete - object is disposed but not fully cleaned up
                    // The object remains marked as disposed to prevent further use.
                    // Use Volatile.Write for cross-thread visibility of the retry flag.
                    Volatile.Write(ref _disposalIncomplete, true);
                    
                    // Log warning for debugging - visible in debug output
                    Debug.WriteLine(
                        "WARNING: M2M.Dispose() could not acquire write lock within timeout. " +
                        "Object is marked disposed but cleanup is incomplete. " +
                        "Check IsDisposalIncomplete property and retry Dispose() after ensuring all threads have finished.");
                    
                    return;
                }
            }
            catch (ObjectDisposedException)
            {
                // Lock already disposed (shouldn't happen, but be defensive)
                return;
            }

            // Lock was acquired - cleanup succeeded, clear incomplete flag
            Volatile.Write(ref _disposalIncomplete, false);

            try
            {
                _rwLock.ExitWriteLock();
            }
            catch
            {
                /* Suppress */
            }

            try
            {
                _rwLock.Dispose();
            }
            catch
            {
                /* Suppress */
            }
        }

        // Note: No unmanaged resources to release
    }

    #endregion

    // ============================================================================
    // GRAPH ALGORITHMS - BFS, DIJKSTRA (THREAD-SAFE)
    // Integrated: December 15, 2024
    // ============================================================================

    #region Graph Algorithms - BFS

    /// <summary>
    ///     Performs thread-safe breadth-first search starting from a single element.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="visitor">
    ///     Optional visitor callback invoked for each discovered element.
    ///     <b>WARNING:</b> Callback is invoked while holding read lock.
    ///     Do NOT call any methods on this M2M instance from the callback.
    ///     Do NOT perform long-running operations in the callback.
    /// </param>
    /// <returns>List of elements in BFS order.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>ALGORITHM:</b> Delegates to O2M.BreadthFirstSearch with cached transpose.
    ///     <b>COMPLEXITY:</b> O(V + E) where V = element count, E = total node references.
    ///     <b>THREAD SAFETY:</b> Thread-safe with automatic synchronization.
    ///     Acquires read lock for duration of traversal.
    ///     <b>PERFORMANCE:</b> First call after modification triggers cache synchronization.
    ///     <b>USAGE:</b> Flood fill, shortest unweighted paths, connected regions.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public List<int> BreadthFirstSearch(int startElement, Action<int, int>? visitor = null)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return _o2m.BreadthFirstSearch(startElement, _elementsFromNode!, visitor);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Performs BFS and returns distances from start element.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <returns>Dictionary mapping element indices to their BFS distance (hop count).</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>DISTANCES:</b> Unweighted shortest path distances (hop count).
    ///     <b>UNREACHABLE:</b> Elements not in the dictionary are unreachable.
    ///     <b>THREAD SAFETY:</b> Thread-safe with read lock.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<int, int> BreadthFirstDistances(int startElement)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return _o2m.BreadthFirstDistances(startElement, _elementsFromNode!);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Graph Algorithms - Dijkstra

    /// <summary>
    ///     Computes shortest paths from a start element using Dijkstra's algorithm.
    /// </summary>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="edgeWeight">Function computing edge weight between two elements sharing a node.</param>
    /// <returns>Dictionary mapping element indices to (distance, predecessor) tuples.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <exception cref="ArgumentException">Thrown if negative edge weights detected.</exception>
    /// <remarks>
    ///     <b>ALGORITHM:</b> Delegates to O2M.DijkstraShortestPaths with cached transpose.
    ///     <b>COMPLEXITY:</b> O((V + E) log V) with binary heap.
    ///     <b>WEIGHTS:</b> Edge weight function receives (fromElement, toElement, sharedNode).
    ///     Must return non-negative values.
    ///     <b>THREAD SAFETY:</b> Thread-safe with read lock held during computation.
    ///     <b>WARNING:</b> edgeWeight callback is invoked while holding read lock.
    ///     Do NOT call any methods on this M2M instance from the callback.
    ///     <b>GEODESICS:</b> For mesh geodesics, use Euclidean distance between element centroids.
    ///     <example>
    ///         // Compute weighted shortest paths
    ///         var paths = m2m.DijkstraShortestPaths(0, (from, to, node) =>
    ///         ComputeDistance(from, to));
    ///         var shortestPath = M2M.ReconstructPath(paths, targetElement);
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths(
        int startElement,
        Func<int, int, int, double> edgeWeight)
    {
        ThrowIfDisposed();
        ThrowIfInBatch();
        EnterSynchronizedReadLock();
        try
        {
            return _o2m.DijkstraShortestPaths(startElement, _elementsFromNode!, edgeWeight);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Reconstructs shortest path from Dijkstra result.
    /// </summary>
    /// <param name="dijkstraResult">Result from DijkstraShortestPaths.</param>
    /// <param name="targetElement">Target element.</param>
    /// <returns>Path from start to target, or null if unreachable.</returns>
    /// <remarks>
    ///     <b>USAGE:</b> Static helper method to extract path from Dijkstra result.
    ///     <b>PATH FORMAT:</b> Returned list starts with start element and ends with target.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static List<int>? ReconstructPath(
        Dictionary<int, (double Distance, int Predecessor)> dijkstraResult,
        int targetElement)
    {
        return O2M.ReconstructPath(dijkstraResult, targetElement);
    }

    #endregion
}


/// <summary>
///     Represents a multi-type many-to-many relationship structure with advanced optimization capabilities.
/// </summary>
/// <remarks>
///     <para>
///         <b>THREAD SAFETY:</b> All public methods are thread-safe. Uses ReaderWriterLockSlim
///         for efficient concurrent read access while ensuring exclusive write access.
///     </para>
///     <para>
///         <b>LOCK ORDERING (P1-4):</b> This class uses a hierarchical locking strategy:
///         <list type="number">
///             <item>MM2M._rwLock (outer lock) - always acquired first</item>
///             <item>M2M._rwLock (inner lock) - acquired after MM2M lock if needed</item>
///         </list>
///         <b>IMPORTANT:</b> Never acquire MM2M lock while holding an M2M lock.
///         All public MM2M methods acquire the outer lock first, then delegate to M2M.
///         Internal code paths have been audited to ensure consistent lock ordering.
///     </para>
///     <para>
///         <b>DESIGN:</b> Manages relationships between multiple entity types using a matrix
///         of M2M structures. Supports marking elements for erasure, batch compression,
///         duplicate detection, structure validation, and performance tuning.
///     </para>
///     <para>
///         <b>PERFORMANCE:</b> Query operations use O(1) HashSet lookups where possible.
///         The Compress operation is O(n × m) where n is total elements and m is relationships.
///         Supports parallelization tuning per type and memory optimization.
///     </para>
///     <para>
///         <b>ENHANCEMENTS:</b> This version leverages advanced O2M and M2M features including:
///         - Duplicate detection before compression
///         - Structure validation with cycle detection
///         - Topological ordering per type
///         - Parallelization threshold configuration
///         - Memory management (Reserve/ShrinkToFit)
///         - Optional unsorted queries for performance
///     </para>
///     <para>
///         <b>STALE REFERENCE WARNING:</b> After <see cref="Compress" /> is called, any previously
///         obtained M2M references from the indexer become invalid (disposed). Use
///         <see cref="Version" /> to detect when the structure has been rebuilt. Increment
///         <see cref="Version" /> on structural changes for consumer reference validation.
///     </para>
/// </remarks>
public sealed class MM2M : IDisposable
{
    #region Constructor

    /// <summary>
    ///     Initializes a new MM2M instance with the specified number of entity types.
    /// </summary>
    /// <param name="numberOfTypes">The number of distinct entity types to manage.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when numberOfTypes is zero or negative.</exception>
    public MM2M(int numberOfTypes)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(numberOfTypes);

        _mat = new M2M[numberOfTypes, numberOfTypes];
        for (var i = 0; i < numberOfTypes; i++)
        for (var j = 0; j < numberOfTypes; j++)
            _mat[i, j] = new M2M();

        NumberOfTypes = numberOfTypes;

        _listOfMarked = new HashSet<int>[numberOfTypes];
        for (var i = 0; i < numberOfTypes; i++)
            _listOfMarked[i] = new HashSet<int>();
    }

    #endregion

    #region Fields

    private readonly HashSet<int>[] _listOfMarked;
    private readonly M2M[,] _mat;
    private readonly ReaderWriterLockSlim _rwLock = new();

    /// <summary>
    ///     Monotonically increasing version number. Incremented on structural changes
    ///     that invalidate previously obtained M2M references (e.g., Compress).
    /// </summary>
    private long _version;

    /// <summary>
    ///     Tracks whether resources have been disposed (0 = not disposed, 1 = disposed).
    ///     Uses int for Interlocked.CompareExchange compatibility.
    /// </summary>
    private int _disposed;

    #endregion

    #region Properties

    /// <summary>
    ///     Gets the number of entity types managed by this structure.
    /// </summary>
    public int NumberOfTypes { get; }

    /// <summary>
    ///     Gets the current structure version. Incremented when structural operations
    ///     invalidate previously obtained M2M references (e.g., after Compress).
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>⚠️ IMPORTANT (P0.4):</b> Prefer using <see cref="WithBlock" /> methods instead
    ///         of caching M2M references manually. WithBlock guarantees reference validity
    ///         and is less error-prone.
    ///     </para>
    ///     <para>
    ///         <b>USE CASE:</b> Consumers who must cache M2M references from the indexer should
    ///         store the Version at cache time and verify it hasn't changed before use.
    ///         If Version has changed, the cached M2M may be disposed.
    ///     </para>
    ///     <para>
    ///         <b>EXAMPLE:</b>
    ///         <code>
    ///         var cachedVersion = mm2m.Version;
    ///         var cachedM2M = mm2m[i, j];
    ///         // ... later ...
    ///         if (mm2m.Version != cachedVersion)
    ///             throw new InvalidOperationException("Cached M2M is stale");
    ///         </code>
    ///     </para>
    /// </remarks>
    public long Version
    {
        get
        {
            ThrowIfDisposed();
            return Interlocked.Read(ref _version);
        }
    }

    /// <summary>
    ///     Gets a read-only view of the elements marked for erasure by type.
    /// </summary>
    /// <remarks>
    ///     Returns a snapshot of the current marked elements. Changes to the underlying
    ///     structure after this call are not reflected in the returned dictionary.
    /// </remarks>
    public IReadOnlyDictionary<int, IReadOnlySet<int>> ListOfMarked
    {
        get
        {
            ThrowIfDisposed();
            _rwLock.EnterReadLock();
            try
            {
                var result = new Dictionary<int, IReadOnlySet<int>>(_listOfMarked.Length);
                for (var i = 0; i < _listOfMarked.Length; i++)
                    result[i] = new HashSet<int>(_listOfMarked[i]);
                return result;
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
    }

    /// <summary>
    ///     Gets or sets the M2M relationship between two entity types.
    /// </summary>
    /// <param name="elementType">The element (row) type index.</param>
    /// <param name="nodeType">The node (column) type index.</param>
    /// <returns>The M2M structure for the specified type pair.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ArgumentNullException">Thrown when setting a null value.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <para>
    ///         <b>⚠️ P0.4 WARNING - STALE REFERENCES:</b>
    ///     </para>
    ///     <para>
    ///         Cached M2M references obtained from this indexer become <b>INVALID</b> after
    ///         <see cref="Compress" /> is called. The old M2M instances are disposed and
    ///         replaced with new ones. Using a cached reference will throw ObjectDisposedException.
    ///     </para>
    ///     <para>
    ///         <b>RECOMMENDED:</b> Use <see cref="WithBlock(int, int, Action{M2M})" /> or
    ///         <see cref="WithBlock{TResult}(int, int, Func{M2M, TResult})" /> instead, which
    ///         guarantee reference validity.
    ///     </para>
    ///     <para>
    ///         <b>If you must cache:</b> Check <see cref="Version" /> before using cached references.
    ///         Version increments on every Compress() call.
    ///     </para>
    ///     <para>
    ///         <b>LOCKING DESIGN (Review Issue #4):</b>
    ///         Each M2M block has independent locking. MM2M's lock protects access to
    ///         which M2M block is retrieved, but does NOT synchronize operations across
    ///         multiple blocks. If you need to maintain invariants spanning multiple blocks
    ///         (e.g., transpose symmetry between _mat[i,j] and _mat[j,i]), you must
    ///         coordinate externally.
    ///     </para>
    ///     <para>
    ///         This design allows maximum concurrency for independent relationship types,
    ///         but requires care if cross-block invariants are needed.
    ///     </para>
    /// </remarks>
    [Obsolete("P0-4 FIX: Use WithBlock() to avoid stale references after Compress(). " +
              "Cached M2M references become disposed after Compress() is called. " +
              "This indexer may be removed in a future version.")]
    public M2M this[int elementType, int nodeType]
    {
        get
        {
            ThrowIfDisposed();
            ValidateTypeIndex(elementType);
            ValidateTypeIndex(nodeType);
            _rwLock.EnterReadLock();
            try
            {
                return _mat[elementType, nodeType];
            }
            finally
            {
                _rwLock.ExitReadLock();
            }
        }
        set
        {
            ThrowIfDisposed();
            ValidateTypeIndex(elementType);
            ValidateTypeIndex(nodeType);
            ArgumentNullException.ThrowIfNull(value);
            _rwLock.EnterWriteLock();
            try
            {
                _mat[elementType, nodeType] = value;
            }
            finally
            {
                _rwLock.ExitWriteLock();
            }
        }
    }

    #endregion

    #region Safe Block Access (P0.4 FIX)

    /// <summary>
    ///     Executes an action on the specified M2M block with automatic lock management.
    /// </summary>
    /// <param name="elementType">The element type index.</param>
    /// <param name="nodeType">The node type index.</param>
    /// <param name="action">Action to execute on the M2M block.</param>
    /// <exception cref="ArgumentNullException">Thrown when action is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>P0.4 FIX:</b> This is the safe way to access M2M blocks. Unlike the indexer,
    ///     this method guarantees the block reference is valid for the duration of the action.
    ///     The block reference becomes invalid after Compress() disposes and recreates blocks.
    ///     <para>
    ///         <b>UNSAFE PATTERN (DON'T DO THIS):</b>
    ///     </para>
    ///     <code>
    ///     var block = mm2m[3, 5];  // Get reference
    ///     mm2m.Compress();         // Block now disposed!
    ///     block.Add(1, 2);         // BOOM! ObjectDisposedException
    ///     </code>
    ///     <para>
    ///         <b>SAFE PATTERN (DO THIS):</b>
    ///     </para>
    ///     <code>
    ///     mm2m.WithBlock(3, 5, block => {
    ///         block.Add(1, 2);
    ///         block.Add(3, 4);
    ///     });
    ///     </code>
    /// </remarks>
    public void WithBlock(int elementType, int nodeType, Action<M2M> action)
    {
        ArgumentNullException.ThrowIfNull(action);
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);

        _rwLock.EnterReadLock();
        try
        {
            action(_mat[elementType, nodeType]);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes a function on the specified M2M block and returns its result.
    /// </summary>
    /// <typeparam name="TResult">The return type of the function.</typeparam>
    /// <param name="elementType">The element type index.</param>
    /// <param name="nodeType">The node type index.</param>
    /// <param name="func">Function to execute on the M2M block.</param>
    /// <returns>The result of the function.</returns>
    /// <exception cref="ArgumentNullException">Thrown when func is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>P0.4 FIX:</b> Safe pattern for querying M2M blocks. See <see cref="WithBlock(int, int, Action{M2M})" />
    ///     for details on why this is safer than using the indexer directly.
    ///     <para>
    ///         <b>EXAMPLE:</b>
    ///     </para>
    ///     <code>
    ///     int count = mm2m.WithBlock(3, 5, block => block.Count);
    ///     var nodes = mm2m.WithBlock(3, 5, block => block.GetNodes(element));
    ///     </code>
    /// </remarks>
    public TResult WithBlock<TResult>(int elementType, int nodeType, Func<M2M, TResult> func)
    {
        ArgumentNullException.ThrowIfNull(func);
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);

        _rwLock.EnterReadLock();
        try
        {
            return func(_mat[elementType, nodeType]);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Query Methods

    /// <summary>
    ///     Gets all elements of any type that reference the specified node.
    /// </summary>
    /// <param name="nodeType">The type of the node.</param>
    /// <param name="node">The node index.</param>
    /// <param name="sorted">Whether to return results in sorted order (default: true).</param>
    /// <returns>A list of (ElementType, ElementIndex) tuples.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when nodeType is out of range or node is negative.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Time Complexity: O(k × m) where k = number of types, m = avg elements per node.
    ///     Add O(n log n) if sorted=true where n = result count.
    ///     Set sorted=false for ~15-20% performance improvement when order doesn't matter.
    /// </remarks>
    public List<(int ElemType, int Elem)> GetAllElements(int nodeType, int node, bool sorted = true)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(nodeType);
        ArgumentOutOfRangeException.ThrowIfNegative(node);
        _rwLock.EnterReadLock();
        try
        {
            return GetAllElementsUnlocked(nodeType, node, sorted);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the number of nodes of a specific type connected to an element.
    /// </summary>
    /// <param name="elementType">The element type.</param>
    /// <param name="element">The element index.</param>
    /// <param name="nodeType">The node type to count.</param>
    /// <returns>The count of connected nodes.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when element index or type indices are out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public int GetNumberOfNodes(int elementType, int element, int nodeType)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);
        ArgumentOutOfRangeException.ThrowIfNegative(element);

        _rwLock.EnterReadLock();
        try
        {
            var block = _mat[elementType, nodeType];
            if (element >= block.Count)
                throw new ArgumentOutOfRangeException(nameof(element),
                    $"Element index {element} is out of range. Valid range: 0 to {block.Count - 1}.");
            return block[element].Count;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Tries to get the number of nodes of a specific type connected to an element.
    /// </summary>
    /// <param name="elementType">The element type.</param>
    /// <param name="element">The element index.</param>
    /// <param name="nodeType">The node type to count.</param>
    /// <param name="count">When this method returns, contains the node count if successful; otherwise, 0.</param>
    /// <returns>True if the element exists; otherwise, false.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public bool TryGetNumberOfNodes(int elementType, int element, int nodeType, out int count)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);

        if (element < 0)
        {
            count = 0;
            return false;
        }

        _rwLock.EnterReadLock();
        try
        {
            var block = _mat[elementType, nodeType];
            if (element >= block.Count)
            {
                count = 0;
                return false;
            }
            count = block[element].Count;
            return true;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the number of elements of a specific type connected to a node.
    /// </summary>
    /// <param name="nodeType">The node type.</param>
    /// <param name="node">The node index.</param>
    /// <param name="elementType">The element type to count.</param>
    /// <returns>The count of connected elements, or 0 if the node doesn't exist.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public int GetNumberOfElements(int nodeType, int node, int elementType)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(nodeType);
        ValidateTypeIndex(elementType);
        _rwLock.EnterReadLock();
        try
        {
            // Use GetElementCountForNode which is O(1) instead of ElementsFromNode[node].Count
            // which clones O(n×m) twice (once for bounds check, once for count)
            return _mat[elementType, nodeType].GetElementCountForNode(node);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the total number of elements of a specific type.
    /// </summary>
    /// <param name="elementType">The element type.</param>
    /// <returns>The total count of elements of the specified type.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementType is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>FIX:</b> This method now searches all nodeType columns for the maximum count.
    ///     The previous implementation only checked the diagonal _mat[type, type], which
    ///     failed for elements that don't have self-relationships (e.g., Tri3 elements
    ///     stored in _mat[Tri3, Node]).
    ///     For nodes with self-relationships, the diagonal still contains the correct count.
    ///     For elements stored with a different nodeType, the search finds them.
    /// </remarks>
    public int GetNumberOfElements(int elementType)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        _rwLock.EnterReadLock();
        try
        {
            // FIX: Search all nodeType columns for this elementType
            // Elements may be stored at _mat[elementType, nodeType] where nodeType != elementType
            var maxCount = 0;
            for (var nodeType = 0; nodeType < NumberOfTypes; nodeType++)
            {
                var count = _mat[elementType, nodeType].Count;
                if (count > maxCount)
                    maxCount = count;
            }

            return maxCount;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the number of active (not marked for erasure) elements of a specific type.
    /// </summary>
    /// <param name="elementType">The element type.</param>
    /// <returns>The count of active elements.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementType is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>FIX:</b> Now uses GetNumberOfElementsUnlocked to get correct total count
    ///     regardless of which nodeType the elements are stored with.
    /// </remarks>
    public int GetNumberOfActiveElements(int elementType)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        _rwLock.EnterReadLock();
        try
        {
            // FIX: Get total count from all nodeType columns
            var totalCount = GetNumberOfElementsUnlocked(elementType);
            var markedSet = _listOfMarked[elementType];
            var activeCount = 0;

            // Count elements that are not marked for erasure
            for (var i = 0; i < totalCount; i++)
                if (!markedSet.Contains(i))
                    activeCount++;

            return activeCount;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the maximum node index + 1 for a given node type across all element types.
    ///     This is the effective range of node indices [0, result) for that type.
    /// </summary>
    /// <remarks>
    ///     <b>FIX (Issue 3):</b> The previous code used <c>_mat[nodeType, nodeType].Count</c>
    ///     (the diagonal block's element count) as a proxy for node range. When the diagonal
    ///     block is empty (nodes typically don't reference other nodes), this returned 0 and
    ///     made <see cref="GetAllElements(int, bool)"/> return empty results.
    ///     <para/>
    ///     MUST be called under read or write lock.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int GetNodeRangeForTypeUnlocked(int nodeType)
    {
        // FIX: Initialize to -1 so we can distinguish "no nodes found" from "only node 0 exists".
        // Previously initialized to 0, which caused node index 0 to be indistinguishable
        // from the empty case, and the return value was maxNode instead of maxNode + 1,
        // causing GetAllElements(nodeType) to skip the highest-indexed node.
        var maxNode = -1;
        for (var et = 0; et < NumberOfTypes; et++)
        {
            var block = _mat[et, nodeType];
            if (block.Count > 0)
            {
                var blockMax = block.GetMaxNode();
                if (blockMax > maxNode)
                    maxNode = blockMax;
            }
        }

        // Return maxNode + 1 to match documented behavior: "the effective range [0, result)"
        // When no nodes exist, returns 0 (-1 + 1 = 0). Otherwise returns maxIndex + 1.
        return maxNode + 1;
    }

    /// <summary>
    ///     Internal helper to get element count without lock acquisition.
    ///     MUST be called under read or write lock.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int GetNumberOfElementsUnlocked(int elementType)
    {
        var maxCount = 0;
        for (var nodeType = 0; nodeType < NumberOfTypes; nodeType++)
        {
            var count = _mat[elementType, nodeType].Count;
            if (count > maxCount)
                maxCount = count;
        }

        return maxCount;
    }

    /// <summary>
    ///     Gets all elements that reference any node of the specified type.
    /// </summary>
    /// <param name="nodeType">The node type.</param>
    /// <param name="sorted">Whether to return results in sorted order (default: true).</param>
    /// <returns>A list of unique (Type, Element) tuples.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when nodeType is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public List<(int Type, int Node)> GetAllElements(int nodeType, bool sorted = true)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(nodeType);
        _rwLock.EnterReadLock();
        try
        {
            var resultSet = new HashSet<(int Type, int Node)>();
            // FIX (Issue 3): Use actual node range instead of diagonal block count.
            // _mat[nodeType, nodeType].Count was 0 when nodes don't self-reference.
            var nodeCount = GetNodeRangeForTypeUnlocked(nodeType);

            for (var n = 0; n < nodeCount; n++)
            {
                var elements = GetAllElementsUnlocked(nodeType, n, false);
                resultSet.UnionWith(elements);
            }

            var result = new List<(int Type, int Node)>(resultSet.Count);
            foreach (var item in resultSet)
                result.Add(item);

            if (sorted)
                result.Sort((x, y) => x.CompareTo(y));

            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all nodes of any type connected to a specific element.
    /// </summary>
    /// <param name="elementType">The element type.</param>
    /// <param name="elementNumber">The element index.</param>
    /// <param name="sorted">Whether to return results in sorted order (default: true).</param>
    /// <returns>A list of unique (Type, Node) tuples.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementType is out of range or elementNumber is negative.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public List<(int Type, int Node)> GetAllNodes(int elementType, int elementNumber, bool sorted = true)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ArgumentOutOfRangeException.ThrowIfNegative(elementNumber);
        _rwLock.EnterReadLock();
        try
        {
            return GetAllNodesUnlocked(elementType, elementNumber, sorted);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all nodes of any type connected to any element of the specified type.
    /// </summary>
    /// <param name="elementType">The element type.</param>
    /// <param name="sorted">Whether to return results in sorted order (default: true).</param>
    /// <returns>A list of unique (Type, Node) tuples.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementType is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public List<(int Type, int Node)> GetAllNodes(int elementType, bool sorted = true)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        _rwLock.EnterReadLock();
        try
        {
            var resultSet = new HashSet<(int Type, int Node)>();
            // FIX: Use GetNumberOfElementsUnlocked instead of diagonal count
            var elementCount = GetNumberOfElementsUnlocked(elementType);

            for (var e = 0; e < elementCount; e++)
            {
                var nodes = GetAllNodesUnlocked(elementType, e, false);
                resultSet.UnionWith(nodes);
            }

            var result = new List<(int Type, int Node)>(resultSet.Count);
            foreach (var item in resultSet)
                result.Add(item);

            if (sorted)
                result.Sort();

            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Performs a depth-first search starting from a node to find all reachable elements.
    /// </summary>
    /// <param name="nodeType">The starting node's type.</param>
    /// <param name="node">The starting node index.</param>
    /// <param name="sorted">Whether to return results in sorted order (default: true).</param>
    /// <returns>A list of all reachable (Type, Element) tuples including the start node.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when nodeType is out of range or node is negative.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public List<(int ElemType, int Elem)> DepthFirstSearchFromANode(int nodeType, int node, bool sorted = true)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(nodeType);
        ArgumentOutOfRangeException.ThrowIfNegative(node);
        _rwLock.EnterReadLock();
        try
        {
            return DepthFirstSearchFromANodeUnlocked(nodeType, node, sorted);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets elements of a specific type that contain all specified nodes.
    /// </summary>
    /// <param name="elementType">The element type to search.</param>
    /// <param name="nodeType">The type of the nodes.</param>
    /// <param name="nodes">The list of nodes that must all be present.</param>
    /// <returns>List of matching element indices.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public List<int> GetElementsFromNodes(int elementType, int nodeType, List<int> nodes)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);
        ArgumentNullException.ThrowIfNull(nodes);
        _rwLock.EnterReadLock();
        try
        {
            return _mat[elementType, nodeType].GetElementsFromNodes(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets elements of a specific type that contain at least one of the specified nodes.
    /// </summary>
    /// <param name="elementType">The element type to search.</param>
    /// <param name="nodeType">The type of the nodes.</param>
    /// <param name="nodes">The list of nodes to search for.</param>
    /// <returns>List of matching element indices.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public List<int> GetElementsWithNodes(int elementType, int nodeType, List<int> nodes)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);
        ArgumentNullException.ThrowIfNull(nodes);
        _rwLock.EnterReadLock();
        try
        {
            return _mat[elementType, nodeType].GetElementsWithNodes(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets elements of a specific type that contain at least one of the specified nodes (union semantics).
    /// </summary>
    /// <param name="elementType">The element type to search.</param>
    /// <param name="nodeType">The type of the nodes.</param>
    /// <param name="nodes">The list of nodes (at least one must be present).</param>
    /// <returns>List of matching element indices.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>DIFFERENCE FROM GetElementsWithNodes:</b>
    ///     - GetElementsFromNodes returns elements containing ALL specified nodes (intersection)
    ///     - GetElementsWithAnyNode returns elements containing AT LEAST ONE (union)
    /// </remarks>
    public List<int> GetElementsWithAnyNode(int elementType, int nodeType, List<int> nodes)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);
        ArgumentNullException.ThrowIfNull(nodes);
        _rwLock.EnterReadLock();
        try
        {
            return _mat[elementType, nodeType].GetElementsContainingAnyNode(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Modification Methods

    /// <summary>
    ///     Appends a new element with connections to nodes of a specific type.
    /// </summary>
    /// <param name="elementType">The type of the new element.</param>
    /// <param name="nodeType">The type of nodes being connected.</param>
    /// <param name="nodes">The list of node indices to connect.</param>
    /// <returns>The index of the newly created element.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range.</exception>
    /// <exception cref="ArgumentNullException">Thrown when nodes is null.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public int AppendElement(int elementType, int nodeType, List<int> nodes)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);
        ArgumentNullException.ThrowIfNull(nodes);
        _rwLock.EnterWriteLock();
        try
        {
            return _mat[elementType, nodeType].AppendElement(nodes);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Marks a node and all its dependent elements for erasure.
    /// </summary>
    /// <param name="nodeType">The type of the node to mark.</param>
    /// <param name="node">The node index to mark.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when nodeType is out of range or node is negative.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     This method performs a depth-first search to find all elements that depend on
    ///     the specified node and marks them all for erasure. The actual removal occurs
    ///     when Compress() is called.
    /// </remarks>
    public void MarkToErase(int nodeType, int node)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(nodeType);
        ArgumentOutOfRangeException.ThrowIfNegative(node);
        _rwLock.EnterWriteLock();
        try
        {
            // Check if already marked (inside lock to prevent race condition)
            if (_listOfMarked[nodeType].Contains(node))
                return;

            // Mark the initial node
            _listOfMarked[nodeType].Add(node);

            // Find and mark all dependent elements
            var dependentElements = DepthFirstSearchFromANodeUnlocked(nodeType, node, false);
            foreach (var pair in dependentElements)
                _listOfMarked[pair.ElemType].Add(pair.Elem);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Compresses the structure by removing all elements marked for erasure.
    /// </summary>
    /// <param name="removeDuplicates">
    ///     If true, automatically detects and marks duplicate elements before compression.
    ///     Duplicates are elements with identical adjacency lists.
    /// </param>
    /// <param name="shrinkMemory">
    ///     If true, reclaims excess memory after compression by calling ShrinkToFit on all M2M structures.
    /// </param>
    /// <returns>
    ///     A list of remapping tuples for each type, where each tuple contains:
    ///     - newNodesFromOld: Maps old indices to new indices (-1 if removed)
    ///     - oldNodesFromNew: Maps new indices back to original indices
    ///     Returns null if no elements were marked for erasure.
    /// </returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     This operation remaps all element indices to be contiguous after removing
    ///     marked elements. All M2M structures are updated to reflect the new indices.
    ///     <b>EXCEPTION SAFETY:</b> This method is transactional. If any error occurs during
    ///     compression, the structure is left unchanged (strong exception guarantee).
    ///     This is achieved using clone-and-swap: changes are applied to clones first,
    ///     then atomically swapped into place only on success.
    ///     <b>EXTERNAL DATA:</b> The returned remapping can be used by callers to update
    ///     external data structures (e.g., attribute lists) that are indexed in parallel
    ///     with the MM2M structure. Use oldNodesFromNew to reorder data lists:
    ///     <code>
    ///         var maps = mm2m.Compress(removeDuplicates: true, shrinkMemory: true);
    ///         if (maps != null)
    ///         {
    ///             var reordered = maps[typeIndex].oldNodesFromNew
    ///                 .Select(oldIdx => dataList[oldIdx]).ToList();
    ///         }
    ///         </code>
    ///     <b>DUPLICATE DETECTION:</b> When removeDuplicates=true, this method scans each type
    ///     for duplicate elements (5-10% typical reduction in refined meshes) and marks them
    ///     for removal before compression. This is particularly useful after mesh refinement
    ///     or merging operations.
    ///     <b>MEMORY OPTIMIZATION:</b> When shrinkMemory=true, excess capacity is reclaimed
    ///     after compression (typically 10-20% memory reduction). This is recommended after
    ///     large batch operations complete.
    /// </remarks>
    public List<(List<int> newNodesFromOld, List<int> oldNodesFromNew)>? Compress(
        bool removeDuplicates = false,
        bool shrinkMemory = false)
    {
        ThrowIfDisposed();

        _rwLock.EnterWriteLock();
        try
        {
            // FIX (Issue 2): Use composite duplicate detection across all node types.
            // Previous code only checked _mat[type, type] (diagonal), missing cross-type
            // connectivity (e.g. Tri3→Node) which is where duplicates typically occur.
            if (removeDuplicates)
                for (var type = 0; type < NumberOfTypes; type++)
                {
                    var duplicates = GetCompositeDuplicatesUnlocked(type);
                    foreach (var dup in duplicates)
                        _listOfMarked[type].Add(dup);
                }

            // Check if there's anything to compress
            var hasMarkedElements = false;
            for (var mi = 0; mi < _listOfMarked.Length; mi++)
                if (_listOfMarked[mi].Count > 0)
                {
                    hasMarkedElements = true;
                    break;
                }

            if (!hasMarkedElements)
                return null;

            // PHASE 1: Build remapping for each type (read-only, no state changes)
            var remap = new List<(List<int> newNodesFromOld, List<int> oldNodesFromNew)>(NumberOfTypes);
            for (var type = 0; type < NumberOfTypes; type++)
            {
                // HashSet is unordered, so use ToSortedList instead of ToList.
                var listOfElementsToKill = Utils.ToSortedList(_listOfMarked[type]);
                // FIX: Use GetNumberOfElementsUnlocked instead of diagonal count
                var numberOfNodesFromType = GetNumberOfElementsUnlocked(type);
                var maxNodeValue = numberOfNodesFromType > 0 ? numberOfNodesFromType - 1 : -1;
                remap.Add(Utils.GetNodeMapsFromKillList(maxNodeValue, listOfElementsToKill));
            }

            // PHASE 2: Validate remaps before modifying state (Issue #8 fix - enhanced validation)
            for (var type = 0; type < NumberOfTypes; type++)
            {
                var (newFromOld, oldFromNew) = remap[type];
                // FIX: Use GetNumberOfElementsUnlocked instead of diagonal count
                var matrixCount = GetNumberOfElementsUnlocked(type);

                // Validate old->new mapping indices are in valid range
                foreach (var oldIdx in oldFromNew)
                    if (oldIdx < 0 || oldIdx >= matrixCount)
                        throw new InvalidOperationException(
                            $"Invalid remap for type {type}: old index {oldIdx} out of range [0, {matrixCount})");

                // Validate new->old mapping 
                for (var oldIdx = 0; oldIdx < newFromOld.Count; oldIdx++)
                {
                    var newIdx = newFromOld[oldIdx];
                    if (newIdx != -1 && (newIdx < 0 || newIdx >= oldFromNew.Count))
                        throw new InvalidOperationException(
                            $"Invalid remap for type {type}: new index {newIdx} for old index {oldIdx} " +
                            $"out of range [0, {oldFromNew.Count})");
                }

                // Validate bijection: oldFromNew should have no duplicate values
                // Review Issue #7: Each old index should appear at most once
                var seenOldIndices = new HashSet<int>();
                foreach (var oldIndexValue in oldFromNew) // Iterate over VALUES (old indices)
                    if (!seenOldIndices.Add(oldIndexValue))
                        throw new InvalidOperationException(
                            $"Invalid remap for type {type}: old index {oldIndexValue} appears multiple times in mapping");
            }

            // PHASE 3: Clone-and-modify for exception safety
            // Clone all M2M structures, apply rearrangement to clones, then swap atomically
            var clonedMat = new M2M[NumberOfTypes, NumberOfTypes];

            try
            {
                // Clone all matrices
                for (var i = 0; i < NumberOfTypes; i++)
                for (var j = 0; j < NumberOfTypes; j++)
                    clonedMat[i, j] = _mat[i, j].CloneTyped();

                // Apply rearrangement to clones
                for (var elementType = 0; elementType < NumberOfTypes; ++elementType)
                {
                    var elementMap = remap[elementType];
                    var oldElementsFromNew = elementMap.oldNodesFromNew;

                    for (var nodeType = 0; nodeType < NumberOfTypes; ++nodeType)
                    {
                        var nodeMap = remap[nodeType];
                        var newNodesFromOld = nodeMap.newNodesFromOld;
                        clonedMat[elementType, nodeType].RearrangeAfterRenumbering(oldElementsFromNew, newNodesFromOld);
                    }
                }
            }
            catch (Exception ex)
            {
                // Clones failed - dispose them and rethrow
                // Original _mat is unchanged (strong exception guarantee)
                for (var i = 0; i < NumberOfTypes; i++)
                for (var j = 0; j < NumberOfTypes; j++)
                    clonedMat[i, j]?.Dispose();

                throw new InvalidOperationException(
                    "Compression failed during remap application. " +
                    "The structure remains unchanged (exception-safe). " +
                    "Marked elements were not cleared - fix the issue and retry.", ex);
            }

            // PHASE 4: Swap references (pure reference assignments, cannot throw),
            // then dispose old matrices separately for exception safety.
            var oldMat = new M2M[NumberOfTypes, NumberOfTypes];
            for (var i = 0; i < NumberOfTypes; i++)
            for (var j = 0; j < NumberOfTypes; j++)
            {
                oldMat[i, j] = _mat[i, j];
                _mat[i, j] = clonedMat[i, j];
            }

            for (var i = 0; i < NumberOfTypes; i++)
            for (var j = 0; j < NumberOfTypes; j++)
                oldMat[i, j].Dispose();

            // Increment version to signal stale reference invalidation
            // Any consumers with cached M2M references should check Version before use
            Interlocked.Increment(ref _version);

            // PHASE 5: Cleanup (only after successful compression)
            // Clear all marked lists
            for (var mi = 0; mi < _listOfMarked.Length; mi++)
                _listOfMarked[mi].Clear();

            // NEW: Optional memory shrinking after compression
            if (shrinkMemory)
                for (var i = 0; i < NumberOfTypes; i++)
                for (var j = 0; j < NumberOfTypes; j++)
                    _mat[i, j].ShrinkToFit();

            return remap;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Validation Methods

    /// <summary>
    ///     Validates the structural integrity of the entire MM2M.
    /// </summary>
    /// <returns>True if all M2M structures are valid and type dependencies are acyclic; false otherwise.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     This method performs comprehensive validation:
    ///     - All M2M structures have valid node indices (non-negative, no duplicates)
    ///     - Diagonal M2M structures (self-referential types) are acyclic
    ///     - Type-level dependencies form a valid DAG
    ///     Time Complexity: O(T² × n × m) where T = number of types, n = avg elements, m = avg nodes.
    ///     This is an expensive operation - use during debugging or after major structural changes.
    /// </remarks>
    public bool ValidateStructure()
    {
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            // Validate all M2M structures
            for (var i = 0; i < NumberOfTypes; i++)
            for (var j = 0; j < NumberOfTypes; j++)
            {
                if (!_mat[i, j].IsValid())
                    return false;

                // Check for cycles in self-referential types
                if (i == j && _mat[i, j].Count > 0)
                    if (!_mat[i, j].IsAcyclic())
                        return false;
            }

            // Validate type-level dependencies are acyclic
            return AreTypesAcyclicUnlocked();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets duplicate elements for a specific type.
    /// </summary>
    /// <param name="elementType">The element type to check for duplicates.</param>
    /// <returns>List of element indices that are duplicates of earlier elements.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementType is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Time Complexity: O(n log n × m) where n = elements, m = avg nodes per element.
    ///     Useful for identifying redundant elements before compression or during mesh quality checks.
    /// </remarks>
    public List<int> GetDuplicates(int elementType)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);

        _rwLock.EnterReadLock();
        try
        {
            // FIX (Issue 2): Use composite duplicate detection across all node types
            return GetCompositeDuplicatesUnlocked(elementType);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all duplicate elements across all types.
    /// </summary>
    /// <returns>Dictionary mapping type index to list of duplicate element indices.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Time Complexity: O(T × n log n × m) where T = types, n = avg elements, m = avg nodes.
    ///     Only returns types that have duplicates (empty types are omitted).
    /// </remarks>
    public Dictionary<int, List<int>> GetAllDuplicates()
    {
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            var result = new Dictionary<int, List<int>>();

            // FIX (Issue 2): Use composite duplicate detection across all node types
            for (var type = 0; type < NumberOfTypes; type++)
            {
                var duplicates = GetCompositeDuplicatesUnlocked(type);
                if (duplicates.Count > 0)
                    result[type] = duplicates;
            }

            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Performance Tuning Methods

    /// <summary>
    ///     Configures performance parameters for a specific type.
    /// </summary>
    /// <param name="typeIndex">The type to configure.</param>
    /// <param name="parallelizationThreshold">
    ///     Threshold for using parallel algorithms. Operations on structures with more elements
    ///     than this value will use parallelization. Default is 4096.
    /// </param>
    /// <param name="reserveCapacity">
    ///     Optional capacity to pre-allocate. Use when you know approximate type size in advance
    ///     to avoid reallocation overhead during construction.
    /// </param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when typeIndex is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     This method configures all M2M structures involving the specified type.
    ///     Use this to optimize for specific workload characteristics:
    ///     - Small types (vertices): Low threshold, no pre-allocation
    ///     - Large types (volume elements): High threshold, pre-allocate
    ///     Example usage for multi-material mesh:
    ///     <code>
    /// // Many small surface elements - enable parallelization early
    /// mm2m.ConfigureType(surfaceType, parallelizationThreshold: 1000, reserveCapacity: 10000);
    /// // Few large volume elements - disable parallelization
    /// mm2m.ConfigureType(volumeType, parallelizationThreshold: 100000, reserveCapacity: 100);
    /// </code>
    /// </remarks>
    public void ConfigureType(int typeIndex, int parallelizationThreshold, int? reserveCapacity = null)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(typeIndex);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(parallelizationThreshold);
        _rwLock.EnterWriteLock();
        try
        {
            // Configure all M2M structures involving this type
            for (var i = 0; i < NumberOfTypes; i++)
            {
                _mat[typeIndex, i].ParallelizationThreshold = parallelizationThreshold;
                _mat[i, typeIndex].ParallelizationThreshold = parallelizationThreshold;

                if (reserveCapacity.HasValue)
                {
                    _mat[typeIndex, i].Reserve(reserveCapacity.Value);
                    if (i != typeIndex)
                        _mat[i, typeIndex].Reserve(reserveCapacity.Value / 4); // Heuristic for cross-type
                }
            }
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Reserves capacity for a specific type pair relationship.
    /// </summary>
    /// <param name="elementType">The element type.</param>
    /// <param name="nodeType">The node type.</param>
    /// <param name="capacity">The capacity to reserve.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when type indices are out of range or capacity is negative.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Use this for fine-grained control when you know the relationship size in advance.
    ///     More efficient than ConfigureType when only specific relationships are large.
    /// </remarks>
    public void Reserve(int elementType, int nodeType, int capacity)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        ValidateTypeIndex(nodeType);
        ArgumentOutOfRangeException.ThrowIfNegative(capacity);
        _rwLock.EnterWriteLock();
        try
        {
            _mat[elementType, nodeType].Reserve(capacity);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Reclaims excess memory from all M2M structures.
    /// </summary>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Time Complexity: O(T² × n × m) where T = types, n = avg elements, m = avg nodes.
    ///     Typically reclaims 10-20% memory after bulk operations or compression.
    ///     Recommended to call after construction or compression completes.
    /// </remarks>
    public void ShrinkToFit()
    {
        ThrowIfDisposed();

        _rwLock.EnterWriteLock();
        try
        {
            for (var i = 0; i < NumberOfTypes; i++)
            for (var j = 0; j < NumberOfTypes; j++)
                _mat[i, j].ShrinkToFit();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Topology Methods

    /// <summary>
    ///     Gets a topological ordering of entity types based on their dependencies.
    /// </summary>
    /// <returns>A list of type indices in topological order.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     A type A depends on type B if there exists at least one element of type A
    ///     that references a node of type B. If there are no dependencies, returns
    ///     types in natural order [0, 1, 2, ...].
    /// </remarks>
    public List<int> GetTypeTopOrder()
    {
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            var (typeDeps, hasDeps) = BuildTypeDependencyGraphUnlocked();
            return hasDeps ? typeDeps.GetTopOrder() : CreateRangeList(NumberOfTypes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Determines whether the type dependency graph is acyclic.
    /// </summary>
    /// <returns>True if there are no circular dependencies between types; otherwise false.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    public bool AreTypesAcyclic()
    {
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            return AreTypesAcyclicUnlocked();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Builds the type-level dependency graph (caller must hold read lock).
    /// </summary>
    private (O2M TypeDeps, bool HasDeps) BuildTypeDependencyGraphUnlocked()
    {
        var typeDeps = new O2M(NumberOfTypes);
        for (var e = 0; e < NumberOfTypes; e++)
            typeDeps.AppendElement([]);

        var hasDeps = false;
        for (var e = 0; e < NumberOfTypes; e++)
        for (var n = 0; n < NumberOfTypes; n++)
            if (n != e && _mat[e, n].Count > 0)
            {
                typeDeps.AppendNodeToElement(e, n);
                hasDeps = true;
            }

        return (typeDeps, hasDeps);
    }

    /// <summary>
    ///     Checks acyclicity without acquiring lock (caller must hold read lock).
    /// </summary>
    private bool AreTypesAcyclicUnlocked()
    {
        var (typeDeps, hasDeps) = BuildTypeDependencyGraphUnlocked();
        return !hasDeps || typeDeps.IsAcyclic();
    }

    /// <summary>
    ///     Gets the topological ordering of elements within a specific type.
    /// </summary>
    /// <param name="elementType">The element type to order.</param>
    /// <returns>List of element indices in topological order.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementType is out of range.</exception>
    /// <exception cref="InvalidOperationException">Thrown when the element graph contains cycles.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Returns the optimal processing order for elements of this type based on their
    ///     dependency structure (self-referential elements). Very useful for:
    ///     - FEM assembly order optimization
    ///     - Dependency-aware processing
    ///     - Constraint ordering
    ///     For types without self-references, returns natural order [0, 1, 2, ...].
    ///     If cycles are detected, throws InvalidOperationException - use GetSortOrder() instead.
    /// </remarks>
    public List<int> GetElementTopologicalOrder(int elementType)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        _rwLock.EnterReadLock();
        try
        {
            // Check if type has self-references
            // FIX (Issue 4): Use unlocked helper to avoid recursive lock acquisition
            if (_mat[elementType, elementType].Count == 0)
                return CreateRangeList(GetNumberOfElementsUnlocked(elementType));

            // Get topological order (throws if cyclic)
            return _mat[elementType, elementType].GetTopologicalOrder();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets a canonical sort order for elements of a specific type.
    /// </summary>
    /// <param name="elementType">The element type to sort.</param>
    /// <returns>List of element indices in lexicographic sort order.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when elementType is out of range.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     Time Complexity: O(n log n × m) where n = elements, m = avg nodes per element.
    ///     Useful for deterministic output, testing, or when topological order fails due to cycles.
    /// </remarks>
    public List<int> GetElementSortOrder(int elementType)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(elementType);
        _rwLock.EnterReadLock();
        try
        {
            if (_mat[elementType, elementType].Count == 0)
                return CreateRangeList(GetNumberOfElementsUnlocked(elementType));
            return _mat[elementType, elementType].GetSortOrder();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Internal Methods (Unlocked)

    /// <summary>
    ///     Internal implementation of GetAllElements that assumes lock is already held.
    /// </summary>
    /// <remarks>
    ///     <b>FIX (Issue 1):</b> Removed <c>if (et == nodeType) continue</c> guard that silently
    ///     skipped self-type elements. This caused <see cref="MarkToErase"/> to miss elements of
    ///     the same type as the starting node.
    ///     <para/>
    ///     <b>FIX (Issue 5):</b> Replaced <c>ElementsFromNode</c> property (deep-clones the entire
    ///     transpose O2M per type) with <c>HasNode</c> + <c>GetElementsForNode</c> (copies only the
    ///     single node's element list). Avoids O(T × n × m) cloning overhead.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private List<(int ElemType, int Elem)> GetAllElementsUnlocked(int nodeType, int node, bool sorted)
    {
        if (sorted)
        {
            // Use SortedSet for automatic ordering
            var resultSet = new SortedSet<(int ElemType, int Elem)>();

            for (var et = 0; et < NumberOfTypes; et++)
            {
                var block = _mat[et, nodeType];
                if (block.Count == 0 || !block.HasNode(node)) continue;

                foreach (var elem in block.GetElementsForNode(node))
                    resultSet.Add((et, elem));
            }

            var result = new List<(int ElemType, int Elem)>(resultSet.Count);
            foreach (var item in resultSet)
                result.Add(item);
            return result;
        }
        else
        {
            // Use regular list for better performance when order doesn't matter
            var result = new List<(int ElemType, int Elem)>();

            for (var et = 0; et < NumberOfTypes; et++)
            {
                var block = _mat[et, nodeType];
                if (block.Count == 0 || !block.HasNode(node)) continue;

                foreach (var elem in block.GetElementsForNode(node))
                    result.Add((et, elem));
            }

            return result;
        }
    }

    /// <summary>
    ///     Internal implementation of GetAllNodes that assumes lock is already held.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private List<(int Type, int Node)> GetAllNodesUnlocked(int elementType, int elementNumber, bool sorted)
    {
        var resultSet = new HashSet<(int Type, int Node)>();

        for (var nt = 0; nt < NumberOfTypes; nt++)
            if (elementNumber < _mat[elementType, nt].Count)
                foreach (var node in _mat[elementType, nt][elementNumber])
                    resultSet.Add((nt, node));

        var result = new List<(int Type, int Node)>(resultSet.Count);
        foreach (var item in resultSet)
            result.Add(item);

        if (sorted)
            result.Sort();

        return result;
    }

    /// <summary>
    ///     Internal implementation of DFS that assumes lock is already held.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private List<(int ElemType, int Elem)> DepthFirstSearchFromANodeUnlocked(int nodeType, int node, bool sorted)
    {
        var visited = new HashSet<(int ElemType, int Elem)>();
        var stack = new Stack<(int ElemType, int Elem)>();
        stack.Push((nodeType, node));

        while (stack.Count > 0)
        {
            var curr = stack.Pop();
            if (!visited.Add(curr)) continue;

            foreach (var e in GetAllElementsUnlocked(curr.ElemType, curr.Elem, false))
                stack.Push(e);
        }

        var result = new List<(int ElemType, int Elem)>(visited.Count);
        foreach (var item in visited)
            result.Add(item);

        if (sorted)
            result.Sort();

        return result;
    }

    /// <summary>
    ///     Detects duplicate elements of the given type by comparing composite adjacency
    ///     signatures across ALL node types, not just the diagonal block.
    /// </summary>
    /// <remarks>
    ///     <b>FIX (Issue 2):</b> The previous implementation called
    ///     <c>_mat[type, type].GetDuplicates()</c>, which only inspected the self-referential
    ///     block. In typical FEM meshes the defining connectivity is cross-type (e.g. Tri3→Node),
    ///     so diagonal-only detection was effectively a no-op.
    ///     <para/>
    ///     Two elements of type T are duplicates iff their adjacency lists are identical
    ///     across every node type: <c>_mat[T, 0]</c>, <c>_mat[T, 1]</c>, …, <c>_mat[T, N-1]</c>.
    ///     <para/>
    ///     MUST be called under read or write lock.
    /// </remarks>
    private List<int> GetCompositeDuplicatesUnlocked(int elementType)
    {
        var elementCount = GetNumberOfElementsUnlocked(elementType);
        if (elementCount <= 1) return [];

        var firstSeen = new Dictionary<string, int>(elementCount);
        var duplicates = new List<int>();

        for (var elem = 0; elem < elementCount; elem++)
        {
            var sig = BuildElementSignatureUnlocked(elementType, elem);
            if (!firstSeen.TryAdd(sig, elem))
                duplicates.Add(elem);
        }

        return duplicates;
    }

    /// <summary>
    ///     Builds a composite string signature for an element by concatenating its adjacency
    ///     lists from all node types. Used for duplicate detection.
    /// </summary>
    /// <remarks>
    ///     MUST be called under read or write lock.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private string BuildElementSignatureUnlocked(int elementType, int element)
    {
        var sb = new System.Text.StringBuilder();
        for (var nt = 0; nt < NumberOfTypes; nt++)
        {
            var block = _mat[elementType, nt];
            if (element < block.Count)
            {
                var nodes = block[element];
                for (var i = 0; i < nodes.Count; i++)
                {
                    if (i > 0) sb.Append(',');
                    sb.Append(nodes[i]);
                }
            }

            sb.Append('|');
        }

        return sb.ToString();
    }

    /// <summary>
    ///     Creates a list containing integers from 0 to count-1.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static List<int> CreateRangeList(int count)
    {
        var result = new List<int>(count);
        for (var i = 0; i < count; i++)
            result.Add(i);
        return result;
    }

    #endregion

    #region Validation

    /// <summary>
    ///     Validates that a type index is within the valid range.
    /// </summary>
    /// <param name="typeIndex">The type index to validate.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when typeIndex is out of range.</exception>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ValidateTypeIndex(int typeIndex)
    {
        if (typeIndex < 0 || typeIndex >= NumberOfTypes)
            throw new ArgumentOutOfRangeException(
                nameof(typeIndex),
                $"Type index {typeIndex} is out of range. Valid range is [0, {NumberOfTypes - 1}].");
    }

    /// <summary>
    ///     Throws ObjectDisposedException if this instance has been disposed.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ThrowIfDisposed()
    {
        if (Volatile.Read(ref _disposed) != 0)
            throw new ObjectDisposedException(nameof(MM2M));
    }

    #endregion

    #region IDisposable

    /// <summary>
    ///     Finalizer for MM2M.
    /// </summary>
    /// <remarks>
    ///     NOTE: This finalizer does NOT dispose the ReaderWriterLockSlim or contained M2M
    ///     instances (managed resources are not disposed from finalizers per standard dispose
    ///     pattern). The finalizer only sets the disposed flag. For deterministic cleanup,
    ///     call Dispose() explicitly.
    /// </remarks>
    ~MM2M()
    {
        Dispose(false);
    }

    /// <summary>
    ///     Releases all resources used by the MM2M instance.
    /// </summary>
    /// <remarks>
    ///     <b>EXPLICIT DISPOSAL RECOMMENDED:</b> Call Dispose() or use a using statement
    ///     for deterministic cleanup of the ReaderWriterLockSlim and all contained M2M instances.
    /// </remarks>
    public void Dispose()
    {
        Dispose(true);
        GC.SuppressFinalize(this);
    }

    /// <summary>
    ///     Protected implementation of Dispose pattern.
    /// </summary>
    /// <param name="disposing">True if called from Dispose(), false if called from finalizer.</param>
    private void Dispose(bool disposing)
    {
        // Atomic check-and-set to ensure disposal happens only once
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0)
            return;

        if (disposing)
        {
            // STEP 1: Acquire write lock to quiesce all operations
            var lockAcquired = false;
            try
            {
                lockAcquired = _rwLock.TryEnterWriteLock(TimeSpan.FromSeconds(30));
                if (!lockAcquired)
                    Debug.WriteLine(
                        "WARNING: MM2M.Dispose() could not acquire write lock within timeout. " +
                        "Skipping child M2M disposal to avoid corrupting concurrent readers. " +
                        "Child M2M instances will be reclaimed by GC.");
            }
            catch (ObjectDisposedException)
            {
                // Lock already disposed (shouldn't happen, but be defensive)
            }

            try
            {
                // STEP 2: Only dispose child M2M instances if we have exclusive access.
                // Without the write lock, concurrent readers may still be using the M2M
                // instances; disposing them would cause ObjectDisposedException in those readers.
                if (lockAcquired)
                    for (var i = 0; i < NumberOfTypes; i++)
                    for (var j = 0; j < NumberOfTypes; j++)
                        try
                        {
                            _mat[i, j].Dispose();
                        }
                        catch
                        {
                            /* Continue disposing others */
                        }
            }
            finally
            {
                // STEP 3: Release and dispose our lock LAST
                if (lockAcquired)
                    try
                    {
                        _rwLock.ExitWriteLock();
                    }
                    catch
                    {
                        /* Suppress */
                    }

                try
                {
                    _rwLock.Dispose();
                }
                catch
                {
                    /* Suppress */
                }
            }
        }
    }

    #endregion

    // ============================================================================
    // MULTI-TYPE GRAPH ALGORITHMS - BFS, DIJKSTRA
    // Integrated: December 15, 2024
    // ============================================================================

    #region Graph Algorithms - Multi-Type BFS

    /// <summary>
    ///     Performs multi-type breadth-first search starting from an element of a specific type.
    /// </summary>
    /// <param name="startElementType">Starting element type index.</param>
    /// <param name="startElement">Starting element index within that type.</param>
    /// <param name="visitor">Optional visitor callback invoked for each discovered element.</param>
    /// <returns>List of (elementType, elementIndex) tuples in BFS order.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when indices are out of range.</exception>
    /// <remarks>
    ///     <b>ALGORITHM:</b> BFS that traverses across entity types through shared nodes.
    ///     <b>COMPLEXITY:</b> O(V + E) where V = total elements across all types, E = total relationships.
    ///     <b>THREAD SAFETY:</b> Thread-safe with read lock held during traversal.
    ///     <b>USAGE:</b> Discover all entities reachable from a starting entity, across type boundaries.
    ///     <example>
    ///         // Find all entities connected to a vertex
    ///         var reachable = mm2m.BreadthFirstSearchMultiType(vertexType, vertexId);
    ///         foreach (var (type, id) in reachable) {
    ///         Console.WriteLine($"Type {type}, ID {id}");
    ///         }
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public List<(int ElementType, int Element)> BreadthFirstSearchMultiType(
        int startElementType,
        int startElement,
        Action<int, int, int>? visitor = null)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(startElementType);
        ArgumentOutOfRangeException.ThrowIfNegative(startElement);
        _rwLock.EnterReadLock();
        try
        {
            // Track visited entities per type
            var visited = new Dictionary<int, HashSet<int>>();
            for (var t = 0; t < NumberOfTypes; t++)
                visited[t] = new HashSet<int>();

            var result = new List<(int, int)>();
            var queue = new Queue<(int ElemType, int Elem, int Depth)>();

            // Start BFS
            queue.Enqueue((startElementType, startElement, 0));
            visited[startElementType].Add(startElement);
            result.Add((startElementType, startElement));
            visitor?.Invoke(startElementType, startElement, 0);

            while (queue.Count > 0)
            {
                var (currentType, current, depth) = queue.Dequeue();

                // Explore all node types this element is connected to
                for (var nodeType = 0; nodeType < NumberOfTypes; nodeType++)
                {
                    var m2m = _mat[currentType, nodeType];
                    if (m2m.Count == 0) continue;
                    if (current >= m2m.Count) continue;

                    var nodes = m2m.GetNodesForElement(current);

                    // For each node, find elements of other types
                    foreach (var node in nodes)
                        // Check all element types that connect to this node type
                        for (var neighborType = 0; neighborType < NumberOfTypes; neighborType++)
                        {
                            var neighborM2M = _mat[neighborType, nodeType];
                            if (neighborM2M.Count == 0) continue;
                            if (node > neighborM2M.GetMaxNode()) continue;

                            // FIX (Issue 6): Use GetElementsForNode to avoid per-node list allocation
                            var connectedElems = neighborM2M.GetElementsForNode(node);

                            foreach (var neighborElem in connectedElems)
                            {
                                if (visited[neighborType].Contains(neighborElem)) continue;

                                visited[neighborType].Add(neighborElem);
                                queue.Enqueue((neighborType, neighborElem, depth + 1));
                                result.Add((neighborType, neighborElem));
                                visitor?.Invoke(neighborType, neighborElem, depth + 1);
                            }
                        }
                }
            }

            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Performs multi-type BFS and returns distances from start element.
    /// </summary>
    /// <param name="startElementType">Starting element type index.</param>
    /// <param name="startElement">Starting element index.</param>
    /// <returns>Dictionary mapping (elementType, element) to BFS distance (hop count).</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <remarks>
    ///     <b>DISTANCES:</b> Hop count across entity types.
    ///     <b>UNREACHABLE:</b> Entities not in dictionary are unreachable.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<(int ElementType, int Element), int> BreadthFirstDistancesMultiType(
        int startElementType,
        int startElement)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(startElementType);
        ArgumentOutOfRangeException.ThrowIfNegative(startElement);
        _rwLock.EnterReadLock();
        try
        {
            var distances = new Dictionary<(int, int), int>
            {
                [(startElementType, startElement)] = 0
            };

            var queue = new Queue<(int ElemType, int Elem)>();
            queue.Enqueue((startElementType, startElement));

            while (queue.Count > 0)
            {
                var (currentType, current) = queue.Dequeue();
                var currentDist = distances[(currentType, current)];

                for (var nodeType = 0; nodeType < NumberOfTypes; nodeType++)
                {
                    var m2m = _mat[currentType, nodeType];
                    if (m2m.Count == 0 || current >= m2m.Count) continue;

                    var nodes = m2m.GetNodesForElement(current);

                    foreach (var node in nodes)
                        for (var neighborType = 0; neighborType < NumberOfTypes; neighborType++)
                        {
                            var neighborM2M = _mat[neighborType, nodeType];
                            if (neighborM2M.Count == 0) continue;
                            if (node > neighborM2M.GetMaxNode()) continue;

                            // FIX (Issue 6): Use GetElementsForNode to avoid per-node list allocation
                            var connectedElems = neighborM2M.GetElementsForNode(node);

                            foreach (var neighborElem in connectedElems)
                            {
                                var key = (neighborType, neighborElem);
                                if (distances.ContainsKey(key)) continue;

                                distances[key] = currentDist + 1;
                                queue.Enqueue((neighborType, neighborElem));
                            }
                        }
                }
            }

            return distances;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Graph Algorithms - Multi-Type Dijkstra

    /// <summary>
    ///     Computes shortest paths across entity types using Dijkstra's algorithm.
    /// </summary>
    /// <param name="startElementType">Starting element type index.</param>
    /// <param name="startElement">Starting element index.</param>
    /// <param name="edgeWeight">
    ///     Function computing edge weight: (fromType, fromElem, toType, toElem, sharedNodeType,
    ///     sharedNode) -> weight.
    /// </param>
    /// <returns>Dictionary mapping (elementType, element) to (distance, predecessor) tuples.</returns>
    /// <exception cref="ObjectDisposedException">Thrown if the object has been disposed.</exception>
    /// <exception cref="ArgumentException">Thrown if negative weights detected.</exception>
    /// <remarks>
    ///     <b>ALGORITHM:</b> Dijkstra with priority queue, traversing across entity types.
    ///     <b>COMPLEXITY:</b> O((V + E) log V) where V = total entities, E = total relationships.
    ///     <b>WEIGHTS:</b> Edge weight function receives full context including types.
    ///     <b>THREAD SAFETY:</b> Thread-safe with read lock.
    ///     <b>RECONSTRUCTION:</b> Use ReconstructPathMultiType to extract path.
    ///     <example>
    ///         // Weighted shortest path across vertex-edge-face hierarchy
    ///         var paths = mm2m.DijkstraShortestPathsMultiType(
    ///         vertexType, startVertex,
    ///         (fromType, fromElem, toType, toElem, nodeType, node) =>
    ///         ComputeWeight(fromType, fromElem, toType, toElem));
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<(int ElementType, int Element), (double Distance, (int PredType, int PredElem))>
        DijkstraShortestPathsMultiType(
            int startElementType,
            int startElement,
            Func<int, int, int, int, int, int, double> edgeWeight)
    {
        ThrowIfDisposed();
        ValidateTypeIndex(startElementType);
        ArgumentOutOfRangeException.ThrowIfNegative(startElement);
        ArgumentNullException.ThrowIfNull(edgeWeight);
        _rwLock.EnterReadLock();
        try
        {
            var distances = new Dictionary<(int, int), (double, (int, int))>
            {
                [(startElementType, startElement)] = (0.0, (-1, -1))
            };

            var pq = new PriorityQueue<(int Type, int Elem), double>();
            pq.Enqueue((startElementType, startElement), 0.0);
            var finalized = new HashSet<(int, int)>();

            while (pq.Count > 0)
            {
                var (currentType, current) = pq.Dequeue();
                var key = (currentType, current);

                if (finalized.Contains(key)) continue;
                finalized.Add(key);

                var currentDist = distances[key].Item1;

                // Explore neighbors across types
                for (var nodeType = 0; nodeType < NumberOfTypes; nodeType++)
                {
                    var m2m = _mat[currentType, nodeType];
                    if (m2m.Count == 0 || current >= m2m.Count) continue;

                    var nodes = m2m.GetNodesForElement(current);

                    foreach (var node in nodes)
                        for (var neighborType = 0; neighborType < NumberOfTypes; neighborType++)
                        {
                            var neighborM2M = _mat[neighborType, nodeType];
                            if (neighborM2M.Count == 0) continue;
                            if (node > neighborM2M.GetMaxNode()) continue;

                            // FIX (Issue 6): Use GetElementsForNode to avoid per-node list allocation
                            var connectedElems = neighborM2M.GetElementsForNode(node);

                            foreach (var neighborElem in connectedElems)
                            {
                                var neighborKey = (neighborType, neighborElem);
                                if (finalized.Contains(neighborKey)) continue;
                                if (neighborKey == key) continue;

                                var weight = edgeWeight(
                                    currentType, current,
                                    neighborType, neighborElem,
                                    nodeType, node);

                                if (weight < 0)
                                    throw new ArgumentException(
                                        $"Negative weight {weight} detected between " +
                                        $"({currentType},{current}) and ({neighborType},{neighborElem}).");

                                var newDist = currentDist + weight;

                                if (!distances.TryGetValue(neighborKey, out var neighborData) ||
                                    newDist < neighborData.Item1)
                                {
                                    distances[neighborKey] = (newDist, key);
                                    pq.Enqueue((neighborType, neighborElem), newDist);
                                }
                            }
                        }
                }
            }

            return distances;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Reconstructs shortest path from multi-type Dijkstra result.
    /// </summary>
    /// <param name="dijkstraResult">Result from DijkstraShortestPathsMultiType.</param>
    /// <param name="targetElementType">Target element type.</param>
    /// <param name="targetElement">Target element index.</param>
    /// <returns>Path from start to target, or null if unreachable.</returns>
    /// <remarks>
    ///     <b>PATH FORMAT:</b> List of (elementType, elementIndex) from start to target.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static List<(int ElementType, int Element)>? ReconstructPathMultiType(
        Dictionary<(int ElementType, int Element), (double Distance, (int PredType, int PredElem))> dijkstraResult,
        int targetElementType,
        int targetElement)
    {
        var key = (targetElementType, targetElement);
        if (!dijkstraResult.ContainsKey(key))
            return null;

        var path = new List<(int, int)>();
        var current = key;

        while (current != (-1, -1))
        {
            path.Add(current);
            var pred = dijkstraResult[current].Item2;
            current = pred;
        }

        path.Reverse();
        return path;
    }

    #endregion
}
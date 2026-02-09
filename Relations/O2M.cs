using System.Buffers;
using System.Collections.Concurrent;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Text;
using System.Text.Json.Serialization;
using Microsoft.Extensions.ObjectPool;
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
                        if ((uint)element >= (uint)nodesFromElement.Count) continue;
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
                    if ((uint)element >= (uint)nodesFromElement.Count) continue;
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
                if ((uint)n >= (uint)elementsFromNode.Count) continue;

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

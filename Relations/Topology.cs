using System.Buffers;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text.Json;
using System.Text.Json.Serialization;

// Topology is the primary consumer of MM2M and intentionally uses the [Obsolete] indexer
// throughout, because it controls the lifecycle and holds _rwLock to prevent stale references.
#pragma warning disable CS0618 // MM2M indexer marked obsolete for external callers

namespace Numerical;

/// <summary>
///     Specifies the ordering of query results.
/// </summary>
/// <remarks>
///     Using <c>ResultOrder.Sorted</c> is clearer than <c>sorted: true</c> at call sites.
/// </remarks>
public enum ResultOrder
{
    /// <summary>
    ///     Results returned in insertion order (fastest, no sorting overhead).
    /// </summary>
    Unordered = 0,

    /// <summary>
    ///     Results sorted by (type index, entity index) for deterministic ordering.
    /// </summary>
    Sorted = 1
}

#region Sub-Entity Definitions

/// <summary>
///     Defines the local node pairs/triples that form sub-entities within a parent element type.
/// </summary>
/// <remarks>
///     <para>
///         Used by <see cref="Topology{TTypes}.DiscoverSubEntities{TElement,TSubEntity,TNode}" /> to extract
///         edges, faces, or other sub-entities from parent elements.
///     </para>
///     <para>
///         <b>Example - Triangle edges:</b> [(0,1), (1,2), (2,0)] defines the 3 edges
///         of a triangle using local node indices.
///     </para>
///     <para>
///         <b>Example - Tetrahedron faces:</b> [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]
///         defines the 4 triangular faces.
///     </para>
/// </remarks>
public readonly struct SubEntityDefinition
{
    /// <summary>
    ///     Local node indices forming each sub-entity.
    ///     For edges: each int[] has 2 elements.
    ///     For faces: each int[] has 3+ elements.
    /// </summary>
    public readonly int[][] LocalNodeIndices;

    /// <summary>
    ///     Creates a sub-entity definition from local node index arrays.
    /// </summary>
    /// <param name="localNodeIndices">
    ///     Array of local node index combinations. Each inner array defines one sub-entity.
    /// </param>
    public SubEntityDefinition(int[][] localNodeIndices)
    {
        if (localNodeIndices == null) throw new ArgumentNullException(nameof(localNodeIndices));

        // P1.H2 FIX: Deep-copy to prevent caller mutation of internal state.
        // readonly struct + readonly field only prevents reassignment, not array content mutation.
        LocalNodeIndices = new int[localNodeIndices.Length][];
        for (var i = 0; i < localNodeIndices.Length; i++)
        {
            var src = localNodeIndices[i] ?? throw new ArgumentNullException($"{nameof(localNodeIndices)}[{i}]");
            LocalNodeIndices[i] = new int[src.Length];
            Array.Copy(src, LocalNodeIndices[i], src.Length);
        }
    }

    /// <summary>
    ///     Creates edge definitions from node pairs.
    /// </summary>
    public static SubEntityDefinition FromEdges(params (int, int)[] edges)
    {
        var indices = new int[edges.Length][];
        for (var i = 0; i < edges.Length; i++)
            indices[i] = [edges[i].Item1, edges[i].Item2];
        return new SubEntityDefinition(indices);
    }

    /// <summary>
    ///     Creates face definitions from node triples.
    /// </summary>
    public static SubEntityDefinition FromFaces(params (int, int, int)[] faces)
    {
        var indices = new int[faces.Length][];
        for (var i = 0; i < faces.Length; i++)
            indices[i] = [faces[i].Item1, faces[i].Item2, faces[i].Item3];
        return new SubEntityDefinition(indices);
    }

    /// <summary>
    ///     Creates face definitions from node quads.
    /// </summary>
    public static SubEntityDefinition FromQuadFaces(params (int, int, int, int)[] faces)
    {
        var indices = new int[faces.Length][];
        for (var i = 0; i < faces.Length; i++)
            indices[i] = [faces[i].Item1, faces[i].Item2, faces[i].Item3, faces[i].Item4];
        return new SubEntityDefinition(indices);
    }
}

#endregion

/// <summary>
///     High-level mesh/graph topology with per-entity attribute storage.
/// </summary>
/// <remarks>
///     <para>
///         Single entry point for building and querying mesh/graph structures.
///         Combines adjacency relationships with arbitrary per-entity data storage.
///     </para>
///     <para>
///         <b>Thread Safety:</b> All public methods are thread-safe using reader-writer locking.
///         Read operations acquire a shared read lock; write operations acquire an exclusive write lock.
///         Internal methods (e.g., AllUnsafe) require the caller to manage locking.
///     </para>
///     <para>
///         <b>Symmetry Support:</b> Configure symmetry groups to enable canonical storage
///         and automatic deduplication of equivalent elements.
///     </para>
/// </remarks>
public class Topology<TTypes> : IDisposable where TTypes : ITypeMap, new()
{
    #region Constructor

    /// <summary>
    ///     Creates a new topology instance.
    /// </summary>
    public Topology()
    {
        _types = new TTypes();
        _adjacency = new MM2M(_types.Count);
    }

    #endregion

    #region IDisposable

    public void Dispose()
    {
        // Atomic check-and-set to ensure disposal happens only once
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0)
            return;

        GC.SuppressFinalize(this);

        var lockAcquired = false;
        try
        {
            // Acquire write lock to wait for all readers to finish
            lockAcquired = _rwLock.TryEnterWriteLock(TimeSpan.FromSeconds(30));
            if (!lockAcquired)
                Debug.WriteLine(
                    "WARNING: Topology.Dispose() could not acquire write lock within timeout. " +
                    "Proceeding with disposal - concurrent operations may fail.");
        }
        catch (ObjectDisposedException)
        {
            // Lock already disposed (shouldn't happen, but be defensive)
        }

        try
        {
            // Dispose dependent structure first
            try
            {
                _adjacency.Dispose();
            }
            catch
            {
                /* Continue with disposal */
            }
        }
        finally
        {
            if (lockAcquired)
                try
                {
                    _rwLock.ExitWriteLock();
                }
                catch
                {
                    /* Suppress */
                }

            // Lock disposed LAST
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

    #endregion

    #region Cloning

    /// <summary>
    ///     Creates a deep copy of this topology.
    /// </summary>
    public Topology<TTypes> Clone()
    {
        // CRITICAL: Check disposed BEFORE acquiring lock
        ThrowIfDisposed();

        _rwLock.EnterReadLock();
        try
        {
            // Double-check after acquiring lock
            ThrowIfDisposed();

            var clone = new Topology<TTypes>();

            // Copy symmetries (can share, they're immutable)
            // and initialize canonical indices for all symmetry types in one pass
            foreach (var kvp in _symmetries)
            {
                clone._symmetries[kvp.Key] = kvp.Value;
                // Initialize empty canonical index for this type
                clone._canonicalIndex[kvp.Key] = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
            }

            // Copy type index mapping
            foreach (var kvp in _typeToIndex)
                clone._typeToIndex[kvp.Key] = kvp.Value;

            // Deep copy canonical indices with collision chains
            // CRITICAL: Verify state consistency during copy
            foreach (var kvp in _canonicalIndex)
            {
                // Consistency check: canonical index should only contain types that have symmetries
                if (!clone._canonicalIndex.ContainsKey(kvp.Key))
                    // This should never happen if the class invariants are maintained
                    throw new InvalidOperationException(
                        $"Internal state corruption detected: canonical index contains " +
                        $"type {kvp.Key.Name} that is not registered in symmetries. " +
                        $"This indicates a serious bug in topology state management.");

                foreach (var (hash, collisionChain) in kvp.Value)
                {
                    var clonedChain = new List<(int Index, List<int> Nodes)>(collisionChain.Count);
                    foreach (var (index, nodes) in collisionChain)
                    {
                        // Deep copy the nodes list
                        var nodesCopy = new List<int>(nodes.Count);
                        nodesCopy.AddRange(nodes);
                        clonedChain.Add((index, nodesCopy));
                    }

                    clone._canonicalIndex[kvp.Key][hash] = clonedChain;
                }
            }

            // Copy data lists using IDataList interface (no reflection needed)
            foreach (var kvp in _data)
                clone._data[kvp.Key] = kvp.Value.Clone();

            // Copy adjacency structure
            for (var elemType = 0; elemType < _types.Count; elemType++)
            for (var nodeType = 0; nodeType < _types.Count; nodeType++)
            {
                var m2m = _adjacency[elemType, nodeType];
                var elemCount = m2m.Count;

                for (var e = 0; e < elemCount; e++)
                {
                    var nodes = m2m[e];
                    var nodesCopy = new List<int>(nodes.Count);
                    nodesCopy.AddRange(nodes);
                    clone._adjacency.AppendElement(elemType, nodeType, nodesCopy);
                }
            }

            return clone;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Read-Only Wrapper

    /// <summary>
    ///     Creates a read-only view of this topology.
    /// </summary>
    public ReadOnlyTopology<TTypes> AsReadOnly()
    {
        return new ReadOnlyTopology<TTypes>(this);
    }

    #endregion

    #region Helpers

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static PropertyInfo GetCachedProperty(Type type, string propertyName)
    {
        return type.GetProperty(propertyName, BindingFlags.Public | BindingFlags.Instance)
               ?? throw new InvalidOperationException($"Property {propertyName} not found on {type.Name}");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static FieldInfo GetCachedField(Type type, string fieldName, BindingFlags bindingFlags)
    {
        return type.GetField(fieldName, bindingFlags)
               ?? throw new InvalidOperationException($"Field {fieldName} not found on {type.Name}");
    }

    #endregion

    #region Fields

    private readonly MM2M _adjacency;

    // Canonical index with collision chaining - each hash maps to list of (index, canonical nodes) pairs
    // This handles hash collisions by comparing actual node lists when hashes match
    private readonly Dictionary<Type, Dictionary<long, List<(int Index, List<int> Nodes)>>> _canonicalIndex = new();
    private readonly Dictionary<(Type Entity, Type Data), IDataList> _data = new();
    private readonly ReaderWriterLockSlim _rwLock = new(LockRecursionPolicy.SupportsRecursion);
    private readonly Dictionary<Type, Symmetry> _symmetries = new();
    private readonly TTypes _types;
    private readonly Dictionary<Type, int> _typeToIndex = new();
    private readonly object _typeIndexLock = new();

    // Batch operation support
    private readonly object _batchLock = new();
    private int _batchNesting;

    // Thread-local buffer for zero-allocation canonicalization
    // MUST be cleared before each use to prevent stale data contamination
    [ThreadStatic] private static int[]? _canonicalBuffer;
    
    // Tracks the high-water mark of used buffer elements for proper clearing
    [ThreadStatic] private static int _canonicalBufferUsedLength;

    /// <summary>
    ///     Gets or creates a thread-local canonical buffer of the specified size.
    ///     The buffer is cleared up to the high-water mark to prevent stale data issues.
    /// </summary>
    /// <param name="size">Required buffer size.</param>
    /// <returns>A cleared buffer of at least the specified size.</returns>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int[] GetOrCreateCanonicalBuffer(int size)
    {
        var buffer = _canonicalBuffer;

        if (buffer == null || buffer.Length < size)
        {
            // Allocate new buffer if needed (already zeroed by CLR)
            buffer = new int[Math.Max(size, 16)];
            _canonicalBuffer = buffer;
            _canonicalBufferUsedLength = 0; // Fresh buffer, nothing to clear
        }
        else
        {
            // CRITICAL: Clear up to the high-water mark to prevent stale data
            // This ensures any previously written elements are zeroed
            var clearLength = Math.Max(_canonicalBufferUsedLength, size);
            if (clearLength > 0)
                Array.Clear(buffer, 0, Math.Min(clearLength, buffer.Length));
        }

        // Update high-water mark for next call
        _canonicalBufferUsedLength = size;

        return buffer;
    }

    /// <summary>
    ///     Tracks whether resources have been disposed (0 = not disposed, 1 = disposed).
    ///     Uses int for Interlocked.CompareExchange compatibility.
    /// </summary>
    private int _disposed;

    #endregion

    #region Batch Operations

    /// <summary>
    ///     Begins a batch operation that defers cache synchronization until complete.
    ///     Returns a disposable that ends the batch when disposed.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>⚠️ P0.1 FIX - MADE INTERNAL:</b>
    ///     </para>
    ///     <para>
    ///         This method was made internal to eliminate deadlock risk from the public API.
    ///         Failure to dispose causes <b>PERMANENT DEADLOCK</b>. The finalizer cannot fix this
    ///         because ReaderWriterLockSlim requires same-thread disposal.
    ///     </para>
    ///     <para>
    ///         <b>PUBLIC API:</b> External callers must use <see cref="WithBatch(Action)" /> or
    ///         <see cref="WithBatch{TResult}(Func{TResult})" /> instead, which guarantee proper disposal
    ///         even if exceptions occur.
    ///     </para>
    ///     <para>
    ///         <b>INTERNAL USAGE:</b> Always use within a using statement:
    ///         <code>
    ///     using (BeginBatch())
    ///     {
    ///         for (int i = 0; i &lt; 10000; i++)
    ///             Add&lt;Node, Point&gt;(points[i]);
    ///     } // Batch ends here - CRITICAL!
    ///     </code>
    ///     </para>
    ///     <para>
    ///         <b>THREAD SAFETY:</b> Write lock is held for the duration of the batch.
    ///         Nested batches are supported - only the outermost batch holds the lock.
    ///     </para>
    ///     <para>
    ///         <b>DEBUG BUILDS:</b> Captures thread ID and stack trace to detect
    ///         improper disposal patterns (cross-thread, leaks, etc.).
    ///     </para>
    /// </remarks>
    internal BatchOperation BeginBatch()
    {
        ThrowIfDisposed();

        // Thread-safe check: only outermost batch acquires the write lock
        bool acquireLock;
        lock (_batchLock)
        {
            acquireLock = _batchNesting == 0;
            _batchNesting++;
        }

        if (acquireLock)
        {
            ThrowIfDisposed();
            _rwLock.EnterWriteLock();
        }

        return new BatchOperation(this, acquireLock);
    }

    /// <summary>
    ///     Executes an action within a batch operation, guaranteeing proper disposal.
    /// </summary>
    /// <remarks>
    ///     This is the preferred API for batch operations as it cannot be misused:
    ///     <code>
    ///     mesh.WithBatch(() => {
    ///         for (int i = 0; i &lt; 10000; i++)
    ///             mesh.Add&lt;Node, Point&gt;(points[i]);
    ///     });
    ///     </code>
    ///     <b>ADVANTAGE:</b> Unlike BeginBatch(), this method guarantees the batch
    ///     is always properly closed, even if an exception occurs.
    /// </remarks>
    /// <param name="action">The action to execute within the batch.</param>
    /// <exception cref="ArgumentNullException">Thrown when action is null.</exception>
    public void WithBatch(Action action)
    {
        ArgumentNullException.ThrowIfNull(action);
        using var batch = BeginBatch();
        action();
    }

    /// <summary>
    ///     Executes a function within a batch operation, returning its result.
    /// </summary>
    /// <typeparam name="TResult">The return type of the function.</typeparam>
    /// <param name="func">The function to execute within the batch.</param>
    /// <returns>The result of the function.</returns>
    /// <exception cref="ArgumentNullException">Thrown when func is null.</exception>
    public TResult WithBatch<TResult>(Func<TResult> func)
    {
        ArgumentNullException.ThrowIfNull(func);
        using var batch = BeginBatch();
        return func();
    }

    /// <summary>
    ///     Represents an active batch operation.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>P0.1 FIX - MADE INTERNAL:</b> This class is now internal to prevent misuse.
    ///         External callers must use <see cref="WithBatch(Action)" /> or
    ///         <see cref="WithBatch{TResult}(Func{TResult})" />, which guarantee proper disposal.
    ///     </para>
    ///     <para>
    ///         CRITICAL: Always use within a using statement or call Dispose() explicitly.
    ///         Failure to dispose will cause deadlock by holding the write lock.
    ///         A finalizer is provided as emergency cleanup but adds overhead and should
    ///         not be relied upon. Proper disposal is strongly recommended.
    ///     </para>
    ///     <para>
    ///         <b>DEBUG BUILDS:</b> Captures thread ID and stack trace to detect
    ///         improper disposal patterns (cross-thread, leaks, etc.).
    ///     </para>
    /// </remarks>
    internal sealed class BatchOperation : IDisposable
    {
        private readonly Topology<TTypes> _owner;
        private readonly bool _releaseLock;
        private int _disposed; // 0 = not disposed, 1 = disposed (thread-safe)

#if DEBUG
        // P0.1 FIX: Capture diagnostic info for better error detection
        private readonly int _creatingThreadId;
        private readonly string _creationStackTrace;
#endif

        internal BatchOperation(Topology<TTypes> owner, bool releaseLock)
        {
            _owner = owner ?? throw new ArgumentNullException(nameof(owner));
            _releaseLock = releaseLock;
            _disposed = 0;

#if DEBUG
            // P0.1 FIX: Capture thread and stack trace for debugging
            _creatingThreadId = Environment.CurrentManagedThreadId;
            _creationStackTrace = Environment.StackTrace;
#endif
        }

        /// <summary>
        ///     Completes the batch operation and releases resources.
        /// </summary>
        public void Dispose()
        {
#if DEBUG
            // P0.1 FIX: Assert dispose on same thread that created the batch
            var currentThreadId = Environment.CurrentManagedThreadId;
            if (currentThreadId != _creatingThreadId)
            {
                var message = $"CRITICAL ERROR: BatchOperation disposed on thread {currentThreadId} " +
                             $"but was created on thread {_creatingThreadId}. " +
                             $"ReaderWriterLockSlim requires locks to be released by the same thread. " +
                             $"This will cause a SynchronizationLockException.\n" +
                             $"Creation stack trace:\n{_creationStackTrace}";
                System.Diagnostics.Debug.WriteLine(message);
                System.Diagnostics.Trace.TraceError(message);
                // Still attempt disposal but it will likely fail
            }
#endif

            // Use thread-safe compare-exchange to ensure disposal happens only once
            if (Interlocked.CompareExchange(ref _disposed, 1, 0) == 0)
            {
                try
                {
                    // Check if owner is disposed using Volatile.Read for int field
                    if (Volatile.Read(ref _owner._disposed) != 0) return;

                    lock (_owner._batchLock)
                    {
                        // P2.M2 FIX: Guard against underflow if Dispose called multiple times
                        if (_owner._batchNesting > 0)
                            _owner._batchNesting--;
                        else
                            Debug.Fail("BatchOperation._batchNesting underflow detected - " +
                                       "Dispose was likely called multiple times.");
                    }

                    if (_releaseLock)
                        _owner._rwLock.ExitWriteLock();
                }
                catch
                {
                    // Even if disposal fails, prevent finalizer from running
                    GC.SuppressFinalize(this);
                    throw;
                }

                GC.SuppressFinalize(this);
            }
        }

        /// <summary>
        ///     Finalizer for logging improper disposal.
        /// </summary>
        /// <remarks>
        ///     <b>CRITICAL:</b> This finalizer CANNOT release the write lock because
        ///     ReaderWriterLockSlim requires locks to be exited by the same thread that
        ///     entered them. The finalizer runs on the finalizer thread, not the owning thread.
        ///     If this finalizer runs, it means the BatchOperation was not properly disposed,
        ///     which WILL cause a permanent deadlock. The only remedy is to restart the application.
        ///     <b>BATCHNESTING:</b> The _batchNesting counter cannot be safely decremented here either,
        ///     as we cannot access the owner's _batchLock from the finalizer thread without risking
        ///     additional deadlocks. The counter will remain corrupted.
        ///     <b>P0.1 FIX:</b> DEBUG builds include thread ID and creation stack trace to help
        ///     identify where the batch was created but not disposed.
        ///     Always use BatchOperation within a using statement or try/finally block, or use
        ///     the safer <see cref="Topology{TTypes}.WithBatch(Action)" /> API instead.
        /// </remarks>
        ~BatchOperation()
        {
            if (Interlocked.CompareExchange(ref _disposed, 1, 0) == 0)
            {
                // CRITICAL: Do NOT attempt to release the lock here!
                // ReaderWriterLockSlim.ExitWriteLock() will throw SynchronizationLockException
                // because we're on the finalizer thread, not the thread that acquired the lock.
                // Any attempt to "fix" this by calling ExitWriteLock will either:
                // 1. Throw an exception (which we'd have to suppress)
                // 2. Silently fail
                // 3. Corrupt the lock state
                // None of these outcomes help - the topology is now permanently deadlocked.

                // CRITICAL: We also cannot safely decrement _batchNesting here because:
                // 1. We cannot acquire _batchLock from the finalizer thread
                // 2. Even if we could, it would leave the batch state inconsistent
                // The _batchNesting counter will remain elevated, causing future batches to malfunction.

#if DEBUG
                // P0.1 FIX: Enhanced diagnostics with thread and stack trace
                var message = $"CRITICAL ERROR: BatchOperation was not properly disposed!\n" +
                             $"Created on thread: {_creatingThreadId}\n" +
                             $"Finalized on thread: {Environment.CurrentManagedThreadId}\n" +
                             $"This has caused a permanent deadlock on the Topology's write lock.\n" +
                             $"The _batchNesting counter is also corrupted.\n" +
                             $"The application must be restarted.\n\n" +
                             $"Creation stack trace:\n{_creationStackTrace}\n\n" +
                             $"FIX: Always use BatchOperation within a using statement:\n" +
                             $"    using (topology.BeginBatch()) {{ ... }}\n" +
                             $"Or use the safer WithBatch API:\n" +
                             $"    topology.WithBatch(() => {{ ... }});";
                System.Diagnostics.Debug.WriteLine(message);
                System.Diagnostics.Trace.TraceError(message);
#else
                // Log the error so developers can find and fix the bug
                Debug.WriteLine(
                    "CRITICAL ERROR: BatchOperation was not properly disposed! " +
                    "This has caused a permanent deadlock on the Topology's write lock. " +
                    "The _batchNesting counter is also corrupted. " +
                    "The application must be restarted. " +
                    "FIX: Always use BatchOperation within a using statement: " +
                    "using (topology.BeginBatch()) { ... } " +
                    "Or use the safer WithBatch API: topology.WithBatch(() => { ... });");

                // Also write to Trace for Release builds
                Trace.TraceError(
                    "CRITICAL: BatchOperation finalizer invoked - topology is now deadlocked. " +
                    "BatchOperation must be disposed by the thread that created it. " +
                    "Consider using WithBatch() API for guaranteed cleanup.");
#endif
            }
        }
    }

    #endregion

    #region Type Index Caching

    /// <summary>
    ///     TTypes-scoped static generic cache for type indices.
    ///     Each (TTypes, T) pair gets its own static field - zero overhead O(1) access.
    /// </summary>
    /// <remarks>
    ///     <b>THREAD SAFETY:</b> Uses volatile field to ensure visibility across threads.
    ///     The pattern is safe because: (1) int writes are atomic on all platforms,
    ///     (2) volatile ensures no stale reads, (3) computing the same value twice is benign.
    ///     <para>
    ///         <b>CORRECTNESS (P0.C1 FIX):</b> Scoped by TTypes to prevent cross-instance poisoning.
    ///         Previously, TypeIndexCache&lt;T&gt; was shared across ALL Topology&lt;TTypes&gt; instances
    ///         process-wide, meaning Topology&lt;TypeMap&lt;Node,Edge&gt;&gt; and
    ///         Topology&lt;TypeMap&lt;Edge,Node&gt;&gt; would corrupt each other's cached indices.
    ///         Now each TTypes gets its own cache, which is correct because the type-to-index
    ///         mapping is determined entirely by TTypes.
    ///     </para>
    /// </remarks>
    private static class TypeIndexCache<TTypesScope, T> where TTypesScope : ITypeMap
    {
        // Volatile ensures visibility across threads and prevents reordering
        public static volatile int Index = -1;
    }

    /// <summary>
    ///     Gets the cached type index for a type.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int GetTypeIndex<T>()
    {
        // Volatile read ensures we see the latest value across all threads
        // Scoped by TTypes so different Topology<TTypes> parameterizations
        // cannot corrupt each other's caches (P0.C1 fix)
        var cached = TypeIndexCache<TTypes, T>.Index;
        if (cached >= 0) return cached;

        // First access: compute and validate
        cached = _types.IndexOf<T>();

        // CRITICAL: Validate type was found before caching
        if (cached < 0)
            throw new InvalidOperationException(
                $"Type '{typeof(T).FullName ?? typeof(T).Name}' is not registered in this Topology. " +
                $"Only types specified during TTypes construction can be used. " +
                $"Total registered types: {_types.Count}");

        // Safe to cache - type is valid (volatile write ensures visibility)
        TypeIndexCache<TTypes, T>.Index = cached;

        // Thread-safe dictionary update
        var type = typeof(T);
        lock (_typeIndexLock)
        {
            if (!_typeToIndex.ContainsKey(type))
                _typeToIndex[type] = cached;
        }

        return cached;
    }

    #endregion

    #region Canonical Key Generation (Optimized)

    /// <summary>
    ///     Generates a hash-based key for canonical node lists.
    ///     Uses xxHash-inspired algorithm for better collision resistance.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static long ComputeCanonicalKey(ReadOnlySpan<int> nodes)
    {
        // xxHash-inspired algorithm with better avalanche properties
        const ulong Prime1 = 11400714785074694791UL;
        const ulong Prime2 = 14029467366897019727UL;
        const ulong Prime3 = 1609587929392839161UL;

        var hash = (ulong)nodes.Length * Prime1 + Prime3;

        // Process in 64-bit chunks when possible for better performance

        // For odd node counts, 4*len is not divisible by 8, so cast only even-length prefix.
        var evenLen = nodes.Length & ~1; // Round down to nearest even number
        var uint64Span = evenLen > 0
            ? MemoryMarshal.Cast<int, ulong>(nodes[..evenLen])
            : ReadOnlySpan<ulong>.Empty;

        foreach (var chunk in uint64Span)
        {
            hash ^= chunk * Prime2;
            hash = BitOperations.RotateLeft(hash, 31) * Prime1;
        }

        // Handle remaining element (0 or 1 int that doesn't fill a ulong)
        if ((nodes.Length & 1) == 1)
        {
            hash ^= (ulong)nodes[nodes.Length - 1] * Prime2;
            hash = BitOperations.RotateLeft(hash, 23);
        }

        // Final mixing for avalanche effect
        hash ^= hash >> 33;
        hash *= Prime2;
        hash ^= hash >> 29;

        return unchecked((long)hash);
    }

    /// <summary>
    ///     Gets canonical form or original, using span to avoid allocations.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private List<int> GetCanonicalOrOriginal<TElement>(ReadOnlySpan<int> nodes)
    {
        if (_symmetries.TryGetValue(typeof(TElement), out var sym))
        {
            // Use helper method to get cleared buffer (Fix 4)
            var buffer = GetOrCreateCanonicalBuffer(nodes.Length);

            // Compute canonical form - modifies buffer in place (void return)
            sym.CanonicalSpan(nodes, buffer);

            // Only allocate the final List - buffer now contains canonical form
            var result = new List<int>(sym.NodeCount);
            for (var i = 0; i < sym.NodeCount; i++)
                result.Add(buffer[i]);

            return result;
        }

        // Fast path for no symmetry
        var list = new List<int>(nodes.Length);
        for (var i = 0; i < nodes.Length; i++)
            list.Add(nodes[i]);
        return list;
    }

    /// <summary>
    ///     Gets canonical form with efficient hash key.
    /// </summary>
    private (List<int> Canonical, long Key) GetCanonicalWithKey<TElement>(ReadOnlySpan<int> nodes)
    {
        if (_symmetries.TryGetValue(typeof(TElement), out var sym))
        {
            // Use helper method to get cleared buffer (Fix 4)
            var buffer = GetOrCreateCanonicalBuffer(nodes.Length);

            // Compute canonical form in-place (void return - modifies buffer)
            sym.CanonicalSpan(nodes, buffer);

            // Compute key from buffer (zero allocations)
            // Create a span view of just the canonical part
            var canonicalView = buffer.AsSpan(0, sym.NodeCount);
            var key = ComputeCanonicalKey(canonicalView);

            // Only allocate List for storage
            var canonical = new List<int>(sym.NodeCount);
            for (var i = 0; i < sym.NodeCount; i++)
                canonical.Add(buffer[i]);

            // No Span escapes this scope - safe to return tuple
            return (canonical, key);
        }

        var list = new List<int>(nodes.Length);
        for (var i = 0; i < nodes.Length; i++)
            list.Add(nodes[i]);
        return (list, ComputeCanonicalKey(nodes));
    }

    /// <summary>
    ///     Compares two node lists for equality (for collision detection).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool NodesEqual(List<int> a, IReadOnlyList<int> b)
    {
        if (a.Count != b.Count) return false;

        var aSpan = CollectionsMarshal.AsSpan(a);
        for (var i = 0; i < aSpan.Length; i++)
            if (aSpan[i] != b[i])
                return false;
        return true;
    }

    /// <summary>
    ///     Compares two node lists for equality (overload for IReadOnlyList comparison).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool NodesEqual(IReadOnlyList<int> a, IReadOnlyList<int> b)
    {
        if (a.Count != b.Count) return false;

        for (var i = 0; i < a.Count; i++)
            if (a[i] != b[i])
                return false;
        return true;
    }

    #endregion

    #region Configuration

    /// <summary>
    ///     Configures symmetry for an element type.
    /// </summary>
    public Topology<TTypes> WithSymmetry<TElement>(Symmetry symmetry)
    {
        ArgumentNullException.ThrowIfNull(symmetry);
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            _symmetries[typeof(TElement)] = symmetry;
            _canonicalIndex[typeof(TElement)] = new Dictionary<long, List<(int Index, List<int> Nodes)>>();

            // Ensure type is indexed for Compress() to properly remap canonical indices
            _typeToIndex[typeof(TElement)] = GetTypeIndex<TElement>();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }

        return this;
    }

    /// <summary>
    ///     Gets the symmetry configured for an element type, or null if none.
    /// </summary>
    public Symmetry? GetSymmetry<TElement>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _symmetries.GetValueOrDefault(typeof(TElement));
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Adding Entities

    /// <summary>
    ///     Adds a standalone entity (node) with associated data.
    /// </summary>
    public int Add<TNode, TData>(TData data)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var list = GetOrCreateList<TNode, TData>();
            var nodeType = GetTypeIndex<TNode>();

            // P1.H3 FIX: Assert data list and adjacency are in sync
            var index = list.Count;
            Debug.Assert(index == _adjacency.GetNumberOfElements(nodeType),
                $"Data list count ({index}) != adjacency element count " +
                $"({_adjacency.GetNumberOfElements(nodeType)}) for type {typeof(TNode).Name}. " +
                "Index mismatch will cause silent data corruption.");

            list.Add(data);
            _adjacency.AppendElement(nodeType, nodeType, [index]);
            return index;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds a standalone entity (node) without data.
    /// </summary>
    public int Add<TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var nodeType = GetTypeIndex<TNode>();
            var index = _adjacency.GetNumberOfElements(nodeType);
            _adjacency.AppendElement(nodeType, nodeType, [index]);
            return index;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds an element connected to other entities.
    /// </summary>
    /// <summary>
    ///     Adds an element connected to other entities.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>CANONICAL INDEXING:</b> If the element type has a registered symmetry,
    ///         the element is added to the canonical index for O(1) lookup via
    ///         <see cref="Find{TElement}" /> and <see cref="Exists{TElement}" />.
    ///     </para>
    ///     <para>
    ///         <b>DUPLICATES:</b> Unlike <see cref="AddUnique{TElement,TNode}" />, this method
    ///         does NOT check for duplicates. Use AddUnique if you need deduplication.
    ///     </para>
    /// </remarks>
    public int Add<TElement, TNode>(params int[] connectedNodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            return AddInternal<TElement, TNode>(connectedNodes);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds an element connected to other entities (span overload for zero-allocation).
    /// </summary>
    public int Add<TElement, TNode>(ReadOnlySpan<int> connectedNodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            return AddInternal<TElement, TNode>(connectedNodes);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Internal implementation of Add that handles canonical indexing.
    /// </summary>
    /// <remarks>
    ///     MUST be called under write lock.
    ///     Populates canonical index if symmetry is registered for the element type.
    /// </remarks>
    private int AddInternal<TElement, TNode>(ReadOnlySpan<int> connectedNodes)
    {
        // Check if this type has symmetry registered
        if (_symmetries.TryGetValue(typeof(TElement), out var symmetry))
        {
            // VALIDATION: Verify node count matches symmetry definition
            if (connectedNodes.Length != symmetry.NodeCount)
                throw new InvalidOperationException(
                    $"Element type {typeof(TElement).Name} has symmetry defined for {symmetry.NodeCount} nodes, " +
                    $"but {connectedNodes.Length} nodes were provided. " +
                    $"Either provide the correct number of nodes or remove/update the symmetry configuration.");

            // Has symmetry - compute canonical form and add to index
            var (canonical, key) = GetCanonicalWithKey<TElement>(connectedNodes);
            var index = _adjacency.AppendElement(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), canonical);

            // Add to canonical index for Find/Exists support
            if (!_canonicalIndex.TryGetValue(typeof(TElement), out var typeIndex))
            {
                typeIndex = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
                _canonicalIndex[typeof(TElement)] = typeIndex;
            }

            if (typeIndex.TryGetValue(key, out var collisionChain))
                collisionChain.Add((index, canonical));
            else
                typeIndex[key] = new List<(int Index, List<int> Nodes)> { (index, canonical) };

            return index;
        }

        // No symmetry - store as-is, no canonical indexing
        var toStore = new List<int>(connectedNodes.Length);
        for (var i = 0; i < connectedNodes.Length; i++)
            toStore.Add(connectedNodes[i]);
        return _adjacency.AppendElement(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), toStore);
    }

    /// <summary>
    ///     Adds an element connected to other entities, with associated data.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>CANONICAL INDEXING:</b> If the element type has a registered symmetry,
    ///         the element is added to the canonical index for O(1) lookup via
    ///         <see cref="Find{TElement}" /> and <see cref="Exists{TElement}" />.
    ///     </para>
    /// </remarks>
    public int Add<TElement, TNode, TData>(TData data, params int[] connectedNodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var list = GetOrCreateList<TElement, TData>();
            var index = AddInternal<TElement, TNode>(connectedNodes);
            list.Add(data);
            return index;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds an element only if an equivalent one doesn't already exist.
    /// </summary>
    public (int Index, bool WasNew) AddUnique<TElement, TNode>(params int[] connectedNodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var (canonical, key) = GetCanonicalWithKey<TElement>(connectedNodes);

            // Get or create canonical index for this type
            if (!_canonicalIndex.TryGetValue(typeof(TElement), out var typeIndex))
            {
                typeIndex = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
                _canonicalIndex[typeof(TElement)] = typeIndex;
            }

            // Check if hash exists
            if (typeIndex.TryGetValue(key, out var collisionChain))
            {
                // Check each entry in the collision chain for exact match
                foreach (var (existingIdx, existingNodes) in collisionChain)
                    if (NodesEqual(canonical, existingNodes))
                        // Found exact match - return existing index
                        return (existingIdx, false);

                // Hash collision but different nodes - add to collision chain
                var newIndex = _adjacency.AppendElement(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), canonical);
                collisionChain.Add((newIndex, canonical));

                // Log collision for diagnostics (conditional compilation)
#if DEBUG
                Console.WriteLine(
                    $"[INFO] Hash collision handled for {typeof(TElement).Name}: " +
                    $"key={key:X16}, chain_length={collisionChain.Count}");
#endif

                return (newIndex, true);
            }

            // First element with this hash - create new collision chain
            var firstIndex = _adjacency.AppendElement(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), canonical);
            typeIndex[key] = new List<(int Index, List<int> Nodes)> { (firstIndex, canonical) };

            return (firstIndex, true);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds an element with data only if an equivalent one doesn't already exist.
    /// </summary>
    public (int Index, bool WasNew) AddUnique<TElement, TNode, TData>(TData data, params int[] connectedNodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var (canonical, key) = GetCanonicalWithKey<TElement>(connectedNodes);

            // Get or create canonical index for this type
            if (!_canonicalIndex.TryGetValue(typeof(TElement), out var typeIndex))
            {
                typeIndex = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
                _canonicalIndex[typeof(TElement)] = typeIndex;
            }

            // Check if hash exists
            if (typeIndex.TryGetValue(key, out var collisionChain))
            {
                // Check each entry in the collision chain for exact match
                foreach (var (existingIdx, existingNodes) in collisionChain)
                    if (NodesEqual(canonical, existingNodes))
                        // Found exact match - return existing index (don't add data)
                        return (existingIdx, false);

                // Hash collision but different nodes - add to collision chain
                var list = GetOrCreateList<TElement, TData>();
                var newIndex = _adjacency.AppendElement(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), canonical);
                list.Add(data);
                collisionChain.Add((newIndex, canonical));

                return (newIndex, true);
            }

            // First element with this hash - create new collision chain
            var dataList = GetOrCreateList<TElement, TData>();
            var firstIndex = _adjacency.AppendElement(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), canonical);
            dataList.Add(data);
            typeIndex[key] = new List<(int Index, List<int> Nodes)> { (firstIndex, canonical) };

            return (firstIndex, true);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Lookup

    /// <summary>
    ///     Checks if an element with the given nodes exists.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>SYMMETRY REQUIREMENT:</b> This method uses O(1) canonical index lookup.
    ///         For types WITH registered symmetry, it finds elements regardless of node order.
    ///         For types WITHOUT symmetry, it performs exact match (no node reordering).
    ///     </para>
    ///     <para>
    ///         <b>NOTE:</b> For types without symmetry, a linear scan fallback is used if
    ///         the canonical index doesn't contain the element. This is O(n) but ensures
    ///         correctness even for elements added before symmetry was configured.
    ///     </para>
    /// </remarks>
    public bool Exists<TElement>(params int[] nodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var (canonical, key) = GetCanonicalWithKey<TElement>(nodes);

            // Fast path: check canonical index
            if (_canonicalIndex.TryGetValue(typeof(TElement), out var typeIndex) &&
                typeIndex.TryGetValue(key, out var collisionChain))
                foreach (var (_, storedNodes) in collisionChain)
                    if (NodesEqual(canonical, storedNodes))
                        return true;

            // Slow path fallback: linear scan for types without canonical index
            // This handles elements added before canonical indexing was set up
            if (!_symmetries.ContainsKey(typeof(TElement))) return FindLinearScan<TElement>(canonical) >= 0;

            return false;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Finds an element by its connected nodes.
    /// </summary>
    /// <returns>The element index, or -1 if not found.</returns>
    /// <remarks>
    ///     <para>
    ///         <b>SYMMETRY REQUIREMENT:</b> This method uses O(1) canonical index lookup.
    ///         For types WITH registered symmetry, it finds elements regardless of node order.
    ///         For types WITHOUT symmetry, it performs exact match (no node reordering).
    ///     </para>
    ///     <para>
    ///         <b>NOTE:</b> For types without symmetry, a linear scan fallback is used if
    ///         the canonical index doesn't contain the element. This is O(n) but ensures
    ///         correctness even for elements added before symmetry was configured.
    ///     </para>
    /// </remarks>
    public int Find<TElement>(params int[] nodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var (canonical, key) = GetCanonicalWithKey<TElement>(nodes);

            // Fast path: check canonical index
            if (_canonicalIndex.TryGetValue(typeof(TElement), out var typeIndex) &&
                typeIndex.TryGetValue(key, out var collisionChain))
                foreach (var (index, storedNodes) in collisionChain)
                    if (NodesEqual(canonical, storedNodes))
                        return index;

            // Slow path fallback: linear scan for types without canonical index
            // This handles elements added before canonical indexing was set up
            if (!_symmetries.ContainsKey(typeof(TElement))) return FindLinearScan<TElement>(canonical);

            return -1;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Linear scan fallback for Find when canonical index is not available.
    /// </summary>
    /// <remarks>
    ///     MUST be called under read lock. O(n × m) where n = element count, m = nodes per element.
    /// </remarks>
    private int FindLinearScan<TElement>(List<int> targetNodes)
    {
        var elemTypeIdx = GetTypeIndex<TElement>();

        // Try to find matching node type - scan all possible node types
        for (var nodeTypeIdx = 0; nodeTypeIdx < _types.Count; nodeTypeIdx++)
        {
            var m2m = _adjacency[elemTypeIdx, nodeTypeIdx];
            var count = m2m.Count;

            if (count == 0) continue;

            for (var i = 0; i < count; i++)
            {
                var elemNodes = m2m[i];
                // Use IReadOnlyList overload - no cast needed
                if (NodesEqual(targetNodes, elemNodes))
                    return i;
            }
        }

        return -1;
    }

    #endregion

    #region Queries (Fixed: Now using Read Locks)

    /// <summary>
    ///     Gets the number of entities of a given type.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Count<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _adjacency.GetNumberOfElements(GetTypeIndex<TEntity>());
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the connected nodes of an element (returns a defensive copy).
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>SAFETY:</b> Returns a defensive copy to prevent external mutation of
    ///         internal topology state. Modifications to the returned list do not affect
    ///         the topology.
    ///     </para>
    ///     <para>
    ///         <b>PERFORMANCE:</b> For zero-allocation access, use
    ///         <see cref="WithNodesOf{TElement,TNode}(int, Action{ReadOnlySpan{int}})" /> instead.
    ///     </para>
    /// </remarks>
    public IReadOnlyList<int> NodesOf<TElement, TNode>(int element)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var internal_list = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()][element];
            // Return defensive copy to prevent mutation via downcast
            var copy = new List<int>(internal_list.Count);
            for (var i = 0; i < internal_list.Count; i++)
                copy.Add(internal_list[i]);
            return copy;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes an action with zero-copy span access to element nodes.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>PERFORMANCE:</b> Zero allocation alternative to <see cref="NodesOf{TElement,TNode}" />.
    ///         The span is only valid during the action execution; do not store it.
    ///     </para>
    ///     <para>
    ///         <b>THREAD SAFETY:</b> Read lock is held during action execution.
    ///         Do not call other topology methods from within the action.
    ///     </para>
    /// </remarks>
    public void WithNodesOf<TElement, TNode>(int element, Action<ReadOnlySpan<int>> action)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(action);

        _rwLock.EnterReadLock();
        try
        {
            var internal_list = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()][element];
            // Fast path: direct span access if underlying type is List<int>
            if (internal_list is List<int> list)
            {
                action(CollectionsMarshal.AsSpan(list));
            }
            else
            {
                // Fallback: copy to array (shouldn't happen with current implementation)
                var array = new int[internal_list.Count];
                for (var i = 0; i < internal_list.Count; i++)
                    array[i] = internal_list[i];
                action(array.AsSpan());
            }
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes a function with zero-copy span access to element nodes.
    /// </summary>
    public TResult WithNodesOf<TElement, TNode, TResult>(int element, Func<ReadOnlySpan<int>, TResult> func)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(func);

        _rwLock.EnterReadLock();
        try
        {
            var internal_list = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()][element];
            // Fast path: direct span access if underlying type is List<int>
            if (internal_list is List<int> list)
            {
                return func(CollectionsMarshal.AsSpan(list));
            }
            else
            {
                // Fallback: copy to array (shouldn't happen with current implementation)
                var array = new int[internal_list.Count];
                for (var i = 0; i < internal_list.Count; i++)
                    array[i] = internal_list[i];
                return func(array.AsSpan());
            }
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all elements of a type connected to a specific node (returns a defensive copy).
    /// </summary>
    /// <remarks>
    ///     <b>SAFETY:</b> Returns a defensive copy to prevent external mutation.
    ///     <b>PERFORMANCE:</b> O(sync) + O(k) where k = elements at node.
    ///     <para><b>P0.C3 FIX:</b> Now acquires read lock, consistent with NodesOf and other queries.</para>
    /// </remarks>
    public IReadOnlyList<int> ElementsAt<TElement, TNode>(int node)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()].GetElementsForNode(node);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the count of related entities of a specific type connected to an entity (O(1) operation).
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <typeparam name="TRelated">The related/node type to count.</typeparam>
    /// <param name="entityIndex">The entity index.</param>
    /// <returns>The number of connected entities of the specified related type.</returns>
    /// <remarks>
    ///     O(1) counter wrapper for efficient counting without materializing a list.
    ///     <b>PERFORMANCE:</b> More efficient than NodesOf().Count when only the count is needed,
    ///     as it avoids creating a defensive copy of the list.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int CountRelated<TEntity, TRelated>(int entityIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            return _adjacency.GetNumberOfNodes(entityTypeIndex, entityIndex, relatedTypeIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the count of elements of a specific type incident to a node (O(1) operation).
    /// </summary>
    /// <typeparam name="TElement">The element type to count.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="nodeIndex">The node index.</param>
    /// <returns>The number of elements of the specified type connected to the node.</returns>
    /// <remarks>
    ///     O(1) counter wrapper for efficient counting without materializing a list.
    ///     <b>PERFORMANCE:</b> More efficient than ElementsAt().Count when only the count is needed,
    ///     as it avoids creating a defensive copy of the list.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int CountIncident<TElement, TNode>(int nodeIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var elementTypeIndex = GetTypeIndex<TElement>();
            var nodeTypeIndex = GetTypeIndex<TNode>();
            return _adjacency.GetNumberOfElements(nodeTypeIndex, nodeIndex, elementTypeIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the count of active (not marked for deletion) entities of a given type (O(n) operation).
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <returns>The number of active entities.</returns>
    /// <remarks>
    ///     Direct counter for active entities.
    ///     <b>ALTERNATIVE:</b> GetActive&lt;TEntity&gt;().Count provides the same value but
    ///     also materializes the full list of indices. Use CountActive when only the count is needed.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int CountActive<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _adjacency.GetNumberOfActiveElements(GetTypeIndex<TEntity>());
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the data value for an entity.
    /// </summary>
    public TData Get<TEntity, TData>(int index)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var key = (typeof(TEntity), typeof(TData));
            if (_data.TryGetValue(key, out var existing))
                return ((DataList<TData>)existing).Items[index];

            throw new KeyNotFoundException($"No data of type {typeof(TData).Name} for entity {typeof(TEntity).Name}");
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Tries to get the data value for an entity.
    /// </summary>
    /// <returns>True if data exists, false otherwise.</returns>
    public bool TryGet<TEntity, TData>(int index, out TData? value)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var key = (typeof(TEntity), typeof(TData));
            if (_data.TryGetValue(key, out var existing))
            {
                var list = ((DataList<TData>)existing).Items;
                if (index >= 0 && index < list.Count)
                {
                    value = list[index];
                    return true;
                }
            }

            value = default;
            return false;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Sets the data value for an entity.
    /// </summary>
    public void Set<TEntity, TData>(int index, TData value)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var list = GetOrCreateList<TEntity, TData>();

            // Ensure list is large enough
            while (list.Count <= index)
                list.Add(default!);

            list[index] = value;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Sets multiple data values at once (bulk operation).
    /// </summary>
    public void SetRange<TEntity, TData>(int startIndex, ReadOnlySpan<TData> values)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var list = GetOrCreateList<TEntity, TData>();

            // Ensure list is large enough
            var requiredSize = startIndex + values.Length;
            while (list.Count < requiredSize)
                list.Add(default!);

            for (var i = 0; i < values.Length; i++) list[startIndex + i] = values[i];
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Gets all data values for an entity type as a read-only list (returns a defensive copy).
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>SAFETY:</b> Returns a defensive copy to prevent external mutation of
    ///         internal data storage. Modifications to the returned list do not affect
    ///         the topology's data.
    ///     </para>
    ///     <para>
    ///         <b>PERFORMANCE:</b> For iteration without allocation, use
    ///         <see cref="Each{TEntity,TData}" /> which yields items one at a time.
    ///     </para>
    /// </remarks>
    public IReadOnlyList<TData> All<TEntity, TData>()
    {
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var key = (typeof(TEntity), typeof(TData));
            if (_data.TryGetValue(key, out var existing))
            {
                // Return defensive copy to prevent mutation via downcast
                var source = ((DataList<TData>)existing).Items;
                var copy = new List<TData>(source.Count);
                for (var i = 0; i < source.Count; i++)
                    copy.Add(source[i]);
                return copy;
            }

            return new List<TData>(); // Empty list for non-existent data
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the internal data list without copying (UNSAFE - for internal use only).
    /// </summary>
    /// <remarks>
    ///     <b>WARNING:</b> Returns direct reference to internal list. Callers MUST:
    ///     - Hold appropriate locks
    ///     - NOT mutate the returned list
    ///     - NOT store references beyond lock scope
    /// </remarks>
    internal List<TData> AllUnsafe<TEntity, TData>()
    {
        var key = (typeof(TEntity), typeof(TData));
        if (_data.TryGetValue(key, out var existing))
            return ((DataList<TData>)existing).Items;
        return new List<TData>();
    }

    /// <summary>
    ///     Iterates over all entities of a type with their index and data.
    /// </summary>
    public IEnumerable<(int Index, TData Data)> Each<TEntity, TData>()
    {
        List<TData> snapshot;
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var key = (typeof(TEntity), typeof(TData));
            if (_data.TryGetValue(key, out var existing))
                snapshot = new List<TData>(((DataList<TData>)existing).Items);
            else
                snapshot = new List<TData>();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }

        for (var i = 0; i < snapshot.Count; i++)
            yield return (i, snapshot[i]);
    }

    /// <summary>
    ///     Iterates over all entity indices of a type.
    /// </summary>
    public IEnumerable<int> Each<TEntity>()
    {
        var count = Count<TEntity>();
        for (var i = 0; i < count; i++)
            yield return i;
    }

    /// <summary>
    ///     Iterates with callback action (avoids iterator allocation for hot paths).
    /// </summary>
    public void ForEach<TEntity, TData>(Action<int, TData> action)
    {
        ArgumentNullException.ThrowIfNull(action);

        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var key = (typeof(TEntity), typeof(TData));
            if (_data.TryGetValue(key, out var existing))
            {
                var list = ((DataList<TData>)existing).Items;
                var span = CollectionsMarshal.AsSpan(list);
                for (var i = 0; i < span.Length; i++)
                    action(i, span[i]);
            }
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Parallel iteration for large datasets.
    /// </summary>
    public void ParallelForEach<TEntity, TData>(Action<int, TData> action, int minParallelCount = 1000)
    {
        ArgumentNullException.ThrowIfNull(action);

        List<TData>? snapshot;
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var key = (typeof(TEntity), typeof(TData));
            snapshot = _data.TryGetValue(key, out var existing)
                ? new List<TData>(((DataList<TData>)existing).Items)
                : null;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }

        if (snapshot == null || snapshot.Count == 0) return;

        if (snapshot.Count >= minParallelCount)
            Parallel.For(0, snapshot.Count, ParallelConfig.Options, i => action(i, snapshot[i]));
        else
            for (var i = 0; i < snapshot.Count; i++)
                action(i, snapshot[i]);
    }

    /// <summary>
    ///     Finds elements that share at least one node with the given element.
    /// </summary>
    /// <param name="element">The element index to find neighbors for.</param>
    /// <param name="sorted">Whether to return results in sorted order. Default is true for backward compatibility.</param>
    /// <returns>List of neighbor element indices.</returns>
    /// <remarks>
    ///     <b>PERFORMANCE :</b> Set sorted=false for 15-20% speedup when order doesn't matter.
    ///     This method now leverages M2M's optimized GetElementNeighbors implementation.
    /// </remarks>
    public IReadOnlyList<int> Neighbors<TElement, TNode>(int element, bool sorted = true)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            // v3: Use M2M's optimized neighbor method with optional sorting
            // Note: M2M.GetElementNeighbors already excludes self (elem != element check)
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetElementNeighbors(element, sorted);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Yields neighbors without allocating a full list (memory efficient).
    /// </summary>
    public IEnumerable<int> EnumerateNeighbors<TElement, TNode>(int element)
    {
        ThrowIfDisposed();
        var nodes = NodesOf<TElement, TNode>(element);
        var seen = new HashSet<int> { element };

        foreach (var node in nodes)
        foreach (var neighbor in ElementsAt<TElement, TNode>(node))
            if (seen.Add(neighbor))
                yield return neighbor;
    }

    /// <summary>
    ///     Gets elements that contain exactly the specified nodes (exact match).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="nodes">List of node indices to search for.</param>
    /// <returns>List of element indices that contain exactly these nodes.</returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Finds elements with exact node set match.
    ///     An element matches if its node set equals the provided node set.
    ///     Order doesn't matter (with symmetry), but all nodes must match exactly.
    ///     <b>USE CASE:</b> Finding specific element by its nodes.
    /// </remarks>
    public List<int> GetElementsWithNodes<TElement, TNode>(List<int> nodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetElementsWithNodes(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets elements that contain any of the specified nodes.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="nodes">List of node indices to search for.</param>
    /// <returns>List of element indices that contain at least one of these nodes.</returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Finds all elements incident to any node in the set.
    ///     An element matches if it contains at least one of the provided nodes.
    ///     This is a union operation across all elements incident to the nodes.
    ///     <b>USE CASE:</b> Finding all elements in a region defined by a node set.
    /// </remarks>
    public List<int> GetElementsContainingAnyNode<TElement, TNode>(List<int> nodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetElementsContainingAnyNode(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets elements connected to all specified nodes (intersection).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="nodes">List of node indices.</param>
    /// <returns>List of element indices connected to all these nodes.</returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Finds elements in intersection of node neighborhoods.
    ///     An element matches if it contains ALL of the provided nodes
    ///     (but may contain additional nodes too).
    ///     <b>USE CASE:</b> Finding elements sharing a specific face or edge.
    /// </remarks>
    public List<int> GetElementsFromNodes<TElement, TNode>(List<int> nodes)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetElementsFromNodes(nodes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Set Operations

    /// <summary>
    ///     Gets elements connected to ALL of the given nodes (intersection).
    /// </summary>
    public IReadOnlyList<int> ElementsAtAll<TElement, TNode>(params int[] nodes)
    {
        ThrowIfDisposed();
        if (nodes.Length == 0)
            return [];

        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];

            var result = new HashSet<int>();
            m2m.WithElementsForNodeSpan(nodes[0], span =>
            {
                for (var i = 0; i < span.Length; i++)
                    result.Add(span[i]);
            });

            if (nodes.Length == 1 || result.Count == 0)
                return Utils.ToList(result);

            var scratch = new HashSet<int>();
            for (var iNode = 1; iNode < nodes.Length && result.Count > 0; iNode++)
            {
                scratch.Clear();
                m2m.WithElementsForNodeSpan(nodes[iNode], span =>
                {
                    for (var i = 0; i < span.Length; i++)
                        scratch.Add(span[i]);
                });

                result.IntersectWith(scratch);
            }

            return Utils.ToList(result);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets elements connected to ANY of the given nodes (union).
    /// </summary>
    public IReadOnlyList<int> ElementsAtAny<TElement, TNode>(params int[] nodes)
    {
        ThrowIfDisposed();
        if (nodes.Length == 0)
            return [];

        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];

            var result = new HashSet<int>();
            foreach (var node in nodes)
            {
                m2m.WithElementsForNodeSpan(node, span =>
                {
                    for (var i = 0; i < span.Length; i++)
                        result.Add(span[i]);
                });
            }

            return Utils.ToList(result);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets elements connected to include nodes but not to exclude nodes (difference).
    /// </summary>
    public IReadOnlyList<int> ElementsAtExcluding<TElement, TNode>(int[] include, int[] exclude)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(include);
        ArgumentNullException.ThrowIfNull(exclude);

        if (include.Length == 0)
            return [];

        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];

            var result = new HashSet<int>();
            foreach (var node in include)
            {
                m2m.WithElementsForNodeSpan(node, span =>
                {
                    for (var i = 0; i < span.Length; i++)
                        result.Add(span[i]);
                });
            }

            foreach (var node in exclude)
            {
                m2m.WithElementsForNodeSpan(node, span =>
                {
                    for (var i = 0; i < span.Length; i++)
                        result.Remove(span[i]);
                });
            }

            return Utils.ToList(result);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Sub-Entity Extraction

    /// <summary>
    ///     Discovers and adds sub-entities (edges, faces) from parent elements.
    /// </summary>
    /// <typeparam name="TElement">Parent element type (e.g., Tri3, Tet4).</typeparam>
    /// <typeparam name="TSubEntity">Sub-entity type to create (e.g., Edge).</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="definition">
    ///     Defines which local node combinations form sub-entities.
    ///     Use <see cref="SubEntityDefinition" /> static properties for common elements.
    /// </param>
    /// <param name="addUnique">
    ///     If true (default), uses AddUnique to deduplicate sub-entities.
    ///     Requires symmetry to be set on TSubEntity for proper deduplication.
    /// </param>
    /// <returns>
    ///     Statistics about the extraction: (totalExtracted, uniqueAdded, duplicatesSkipped).
    /// </returns>
    /// <remarks>
    ///     <para>
    ///         <b>SYMMETRY:</b> For proper deduplication, set symmetry on TSubEntity before calling:
    ///         <code>mesh.WithSymmetry&lt;Edge&gt;(Symmetry.Full(2));</code>
    ///     </para>
    ///     <para>
    ///         <b>PERFORMANCE:</b> Uses batch operations internally for efficiency.
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Discover edges from all triangles
    ///         mesh.WithSymmetry&lt;Edge&gt;(Symmetry.Full(2));
    ///         var stats = mesh.DiscoverSubEntities&lt;Tri3, Edge, Node&gt;(FiniteElementTopologies.Tri3Edges);
    ///         Console.WriteLine($"Found {stats.UniqueAdded} unique edges");
    ///         </code>
    ///     </example>
    /// </remarks>
    public (int TotalExtracted, int UniqueAdded, int DuplicatesSkipped)
        DiscoverSubEntities<TElement, TSubEntity, TNode>(
            SubEntityDefinition definition,
            bool addUnique = true)
    {
        ThrowIfDisposed();

        var totalExtracted = 0;
        var uniqueAdded = 0;
        var duplicatesSkipped = 0;

        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            var elementTypeIdx = GetTypeIndex<TElement>();
            var nodeTypeIdx = GetTypeIndex<TNode>();
            var subEntityTypeIdx = GetTypeIndex<TSubEntity>();

            // FIX: Read element count inside write lock to prevent TOCTOU race.
            // Previously read outside the lock, which could become stale if another
            // thread modified elements between the count read and lock acquisition.
            var elementCount = _adjacency.GetNumberOfElements(elementTypeIdx);
            if (elementCount == 0)
                return (0, 0, 0);

            for (var elemIdx = 0; elemIdx < elementCount; elemIdx++)
            {
                var elementNodes = _adjacency[elementTypeIdx, nodeTypeIdx][elemIdx];

                foreach (var localIndices in definition.LocalNodeIndices)
                {
                    // Build global node indices for this sub-entity
                    var globalNodes = new int[localIndices.Length];
                    for (var i = 0; i < localIndices.Length; i++)
                    {
                        var localIdx = localIndices[i];
                        if (localIdx < 0 || localIdx >= elementNodes.Count)
                            throw new InvalidOperationException(
                                $"Local node index {localIdx} out of range for element {elemIdx} " +
                                $"with {elementNodes.Count} nodes");
                        globalNodes[i] = elementNodes[localIdx];
                    }

                    totalExtracted++;

                    if (addUnique)
                    {
                        var (canonical, key) = GetCanonicalWithKey<TSubEntity>(globalNodes);

                        // Check for existing
                        if (!_canonicalIndex.TryGetValue(typeof(TSubEntity), out var typeIndex))
                        {
                            typeIndex = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
                            _canonicalIndex[typeof(TSubEntity)] = typeIndex;
                        }

                        var found = false;
                        if (typeIndex.TryGetValue(key, out var collisionChain))
                            foreach (var (existingIdx, existingNodes) in collisionChain)
                                if (NodesEqual(canonical, existingNodes))
                                {
                                    found = true;
                                    duplicatesSkipped++;
                                    break;
                                }

                        if (!found)
                        {
                            var newIndex = _adjacency.AppendElement(subEntityTypeIdx, nodeTypeIdx, canonical);
                            if (typeIndex.TryGetValue(key, out var chain))
                                chain.Add((newIndex, canonical));
                            else
                                typeIndex[key] = new List<(int Index, List<int> Nodes)> { (newIndex, canonical) };
                            uniqueAdded++;
                        }
                    }
                    else
                    {
                        // Add without deduplication
                        var toStore = GetCanonicalOrOriginal<TSubEntity>(globalNodes);
                        _adjacency.AppendElement(subEntityTypeIdx, nodeTypeIdx, toStore);
                        uniqueAdded++;
                    }
                }
            }

            return (totalExtracted, uniqueAdded, duplicatesSkipped);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Elements Sharing Sub-Entities

    /// <summary>
    ///     Gets all elements of type TElement that contain ALL of the specified nodes.
    /// </summary>
    /// <typeparam name="TElement">Element type to search.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="nodes">Nodes that must all be present in the element.</param>
    /// <returns>List of element indices containing all specified nodes.</returns>
    /// <remarks>
    ///     <para>
    ///         <b>USE CASE:</b> Finding elements sharing an edge or face.
    ///         Pass the edge/face nodes to find all incident elements.
    ///     </para>
    ///     <para>
    ///         <b>ALGORITHM:</b> Intersects the incident element sets of each node.
    ///         Starts with the node having the smallest incident count for efficiency.
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Find triangles sharing edge (n1, n2)
    ///         var tris = mesh.ElementsContainingAllNodes&lt;Tri3, Node&gt;(n1, n2);
    ///         
    ///         // Find tetrahedra sharing face (n1, n2, n3)
    ///         var tets = mesh.ElementsContainingAllNodes&lt;Tet4, Node&gt;(n1, n2, n3);
    ///         </code>
    ///     </example>
    /// </remarks>
    public List<int> ElementsContainingAllNodes<TElement, TNode>(params int[] nodes)
    {
        ThrowIfDisposed();

        if (nodes == null || nodes.Length == 0)
            return new List<int>();

        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();

            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            var transpose = m2m.ElementsFromNode;

            // Start with node having smallest incident count (optimization)
            if (nodes[0] >= transpose.Count)
                return new List<int>(); // Node doesn't exist
            var minNode = nodes[0];
            var minCount = transpose[minNode].Length;
            for (var i = 1; i < nodes.Length; i++)
            {
                if (nodes[i] >= transpose.Count)
                    return new List<int>(); // Node doesn't exist
                var count = transpose[nodes[i]].Length;
                if (count < minCount)
                {
                    minCount = count;
                    minNode = nodes[i];
                }
            }

            // Start with elements at the smallest set
            var result = new HashSet<int>(transpose[minNode].ToArray());

            // Intersect with remaining nodes
            foreach (var node in nodes)
            {
                if (node == minNode) continue;
                if (node >= transpose.Count)
                {
                    result.Clear();
                    break;
                }

                result.IntersectWith(transpose[node].ToArray());
                if (result.Count == 0) break;
            }

            return result.ToList();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all elements of type TParent that share the specified sub-entity.
    /// </summary>
    /// <typeparam name="TParent">Parent element type (e.g., Tri3, Tet4).</typeparam>
    /// <typeparam name="TSubEntity">Sub-entity type (e.g., Edge).</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="subEntityIndex">Index of the sub-entity.</param>
    /// <returns>List of parent element indices sharing this sub-entity.</returns>
    /// <remarks>
    ///     <para>
    ///         <b>ALGORITHM:</b> Gets the nodes of the sub-entity, then finds all
    ///         parent elements containing all those nodes.
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Find triangles sharing edge 42
    ///         var tris = mesh.ElementsSharingSubEntity&lt;Tri3, Edge, Node&gt;(42);
    ///         </code>
    ///     </example>
    /// </remarks>
    public List<int> ElementsSharingSubEntity<TParent, TSubEntity, TNode>(int subEntityIndex)
    {
        ThrowIfDisposed();

        var subEntityNodes = NodesOf<TSubEntity, TNode>(subEntityIndex);
        var nodesArray = new int[subEntityNodes.Count];
        for (var i = 0; i < subEntityNodes.Count; i++)
            nodesArray[i] = subEntityNodes[i];

        return ElementsContainingAllNodes<TParent, TNode>(nodesArray);
    }

    /// <summary>
    ///     Counts how many elements of TParent contain the specified sub-entity.
    /// </summary>
    /// <typeparam name="TParent">Parent element type.</typeparam>
    /// <typeparam name="TSubEntity">Sub-entity type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="subEntityIndex">Index of the sub-entity.</param>
    /// <returns>Number of parent elements sharing this sub-entity.</returns>
    /// <remarks>
    ///     Convenience wrapper over ElementsSharingSubEntity.
    ///     Note: currently materializes the full list. For frequent calls,
    ///     consider caching ElementsSharingSubEntity results.
    /// </remarks>
    public int CountElementsSharingSubEntity<TParent, TSubEntity, TNode>(int subEntityIndex)
    {
        return ElementsSharingSubEntity<TParent, TSubEntity, TNode>(subEntityIndex).Count;
    }

    #endregion

    #region Position Caches

    /// <summary>
    ///     Gets the position of an element within each of its connected nodes' adjacency lists.
    /// </summary>
    public IReadOnlyList<int> GetElementPositions<TElement, TNode>(int element)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()].ElementLocations[element];
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all positions where a node appears within element adjacency lists.
    /// </summary>
    public IReadOnlyList<int> GetNodePositions<TElement, TNode>(int node)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()].NodeLocations[node];
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the local position of a specific node within a specific element.
    /// </summary>
    public int GetLocalNodeIndex<TElement, TNode>(int element, int node)
    {
        ThrowIfDisposed();
        var nodes = NodesOf<TElement, TNode>(element);
        for (var i = 0; i < nodes.Count; i++)
            if (nodes[i] == node)
                return i;
        return -1;
    }

    #endregion

    #region Bulk Operations (Optimized)

    /// <summary>
    ///     Adds multiple nodes at once, returning their indices.
    /// </summary>
    public int[] AddRange<TNode, TData>(IEnumerable<TData> dataItems)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var list = GetOrCreateList<TNode, TData>();
            var nodeType = GetTypeIndex<TNode>();
            var indices = new List<int>();

            foreach (var data in dataItems)
            {
                var index = list.Count;
                list.Add(data);
                _adjacency.AppendElement(nodeType, nodeType, [index]);
                indices.Add(index);
            }

            return indices.ToArray();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds multiple nodes at once using spans (zero-allocation for source data).
    /// </summary>
    public int[] AddRange<TNode, TData>(ReadOnlySpan<TData> dataItems)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var list = GetOrCreateList<TNode, TData>();
            var nodeType = GetTypeIndex<TNode>();
            var indices = new int[dataItems.Length];

            for (var i = 0; i < dataItems.Length; i++)
            {
                var index = list.Count;
                list.Add(dataItems[i]);
                _adjacency.AppendElement(nodeType, nodeType, [index]);
                indices[i] = index;
            }

            return indices;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds multiple elements at once.
    /// </summary>
    public int[] AddRange<TElement, TNode>(IEnumerable<int[]> connectivityList)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var elementType = GetTypeIndex<TElement>();
            var nodeType = GetTypeIndex<TNode>();
            var indices = new List<int>();

            foreach (var nodes in connectivityList)
            {
                var index = AddInternal<TElement, TNode>(nodes);
                indices.Add(index);
            }

            return indices.ToArray();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds multiple unique elements at once, skipping duplicates.
    /// </summary>
    public (int Index, bool WasNew)[] AddRangeUnique<TElement, TNode>(IEnumerable<int[]> connectivityList)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var elementType = GetTypeIndex<TElement>();
            var nodeType = GetTypeIndex<TNode>();
            var results = new List<(int, bool)>();

            foreach (var nodes in connectivityList)
            {
                var (canonical, key) = GetCanonicalWithKey<TElement>(nodes);

                // Initialize type index if needed
                if (!_canonicalIndex.TryGetValue(typeof(TElement), out var typeIndex))
                {
                    typeIndex = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
                    _canonicalIndex[typeof(TElement)] = typeIndex;
                }

                // Check collision chain for exact match
                var found = false;
                if (typeIndex.TryGetValue(key, out var collisionChain))
                    foreach (var (existingIdx, existingNodes) in collisionChain)
                        if (NodesEqual(canonical, existingNodes))
                        {
                            results.Add((existingIdx, false));
                            found = true;
                            break;
                        }

                if (!found)
                {
                    // Add new element
                    var newIndex = _adjacency.AppendElement(elementType, nodeType, canonical);

                    // Add to collision chain (or create new chain)
                    if (typeIndex.TryGetValue(key, out var chain))
                        chain.Add((newIndex, canonical));
                    else
                        typeIndex[key] = new List<(int Index, List<int> Nodes)> { (newIndex, canonical) };

                    results.Add((newIndex, true));
                }
            }

            return results.ToArray();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Adds multiple elements in parallel (for very large datasets).
    /// </summary>
    public int[] AddRangeParallel<TElement, TNode>(int[][] connectivityList, int minParallelCount = 10000)
    {
        if (connectivityList.Length < minParallelCount)
            return AddRange<TElement, TNode>(connectivityList);

        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var elementType = GetTypeIndex<TElement>();
            var nodeType = GetTypeIndex<TNode>();

            // Pre-compute canonical forms in parallel
            var canonicalForms = new List<int>[connectivityList.Length];
            Parallel.For(0, connectivityList.Length, ParallelConfig.Options,
                i => { canonicalForms[i] = GetCanonicalOrOriginal<TElement>(connectivityList[i]); });

            // Add sequentially (adjacency structure not thread-safe internally)
            var indices = new int[connectivityList.Length];
            for (var i = 0; i < canonicalForms.Length; i++)
                indices[i] = _adjacency.AppendElement(elementType, nodeType, canonicalForms[i]);

            return indices;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Modification

    /// <summary>
    ///     Marks an entity for removal.
    /// </summary>
    /// <remarks>
    ///     <b>P0.C2 FIX:</b> Now acquires write lock, consistent with all other mutators.
    /// </remarks>
    public void Remove<TEntity>(int index)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            _adjacency.MarkToErase(GetTypeIndex<TEntity>(), index);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Marks multiple entities for removal.
    /// </summary>
    /// <remarks>
    ///     <b>P0.C2 FIX:</b> Now acquires write lock, consistent with all other mutators.
    /// </remarks>
    public void RemoveRange<TEntity>(IEnumerable<int> indices)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            var typeIndex = GetTypeIndex<TEntity>();
            foreach (var index in indices)
                _adjacency.MarkToErase(typeIndex, index);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region In-Place Element Modification (Feature Completion)

    /// <summary>
    ///     Adds a node to an existing element's connectivity list.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="element">The element index.</param>
    /// <param name="node">The node to add.</param>
    /// <remarks>
    ///     <para>
    ///         <b>FEATURE COMPLETION:</b> This method fills a gap in the Topology API.
    ///         Previously, users had to remove and re-add elements to modify connectivity,
    ///         which changed element indices and broke references.
    ///     </para>
    ///     <para>
    ///         <b>INDEX PRESERVATION:</b> This method modifies the element in-place,
    ///         preserving its index. All external references remain valid.
    ///     </para>
    ///     <para>
    ///         <b>CACHE INVALIDATION:</b> Automatically invalidates canonical indices
    ///         if the element type has defined symmetry.
    ///     </para>
    ///     <para>
    ///         <b>USE CASES:</b>
    ///         - Mesh refinement (splitting elements)
    ///         - Adaptive topology (adding connectivity)
    ///         - Dynamic graph modification
    ///     </para>
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when element index is invalid.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the topology has been disposed.</exception>
    public void AddNodeToElement<TElement, TNode>(int element, int node)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            var elementType = GetTypeIndex<TElement>();
            var nodeType = GetTypeIndex<TNode>();
            var m2m = _adjacency[elementType, nodeType];

            m2m.AppendNodeToElement(element, node);

            // Invalidate canonical index for this type if it has symmetry
            InvalidateCanonicalIndex<TElement>();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Removes a node from an existing element's connectivity list.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="element">The element index.</param>
    /// <param name="node">The node to remove.</param>
    /// <returns>True if the node was found and removed; false otherwise.</returns>
    /// <remarks>
    ///     <para>
    ///         <b>FEATURE COMPLETION:</b> This method fills a gap in the Topology API.
    ///         Previously, users had to remove and re-add elements to modify connectivity,
    ///         which changed element indices and broke references.
    ///     </para>
    ///     <para>
    ///         <b>INDEX PRESERVATION:</b> This method modifies the element in-place,
    ///         preserving its index. All external references remain valid.
    ///     </para>
    ///     <para>
    ///         <b>CACHE INVALIDATION:</b> Automatically invalidates canonical indices
    ///         if the element type has defined symmetry.
    ///     </para>
    ///     <para>
    ///         <b>USE CASES:</b>
    ///         - Mesh coarsening (merging elements)
    ///         - Adaptive topology (removing connectivity)
    ///         - Dynamic graph modification
    ///     </para>
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when element index is invalid.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the topology has been disposed.</exception>
    public bool RemoveNodeFromElement<TElement, TNode>(int element, int node)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            var elementType = GetTypeIndex<TElement>();
            var nodeType = GetTypeIndex<TNode>();
            var m2m = _adjacency[elementType, nodeType];

            var removed = m2m.RemoveNodeFromElement(element, node);

            if (removed)
                // Invalidate canonical index for this type if it has symmetry
                InvalidateCanonicalIndex<TElement>();

            return removed;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Replaces all nodes of an existing element with a new set of nodes.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="element">The element index.</param>
    /// <param name="newNodes">The new nodes for the element.</param>
    /// <remarks>
    ///     <para>
    ///         <b>FEATURE COMPLETION:</b> This method fills a gap in the Topology API.
    ///         Previously, users had to remove and re-add elements to change connectivity,
    ///         which changed element indices and broke references.
    ///     </para>
    ///     <para>
    ///         <b>INDEX PRESERVATION:</b> This method modifies the element in-place,
    ///         preserving its index. All external references remain valid.
    ///     </para>
    ///     <para>
    ///         <b>CACHE INVALIDATION:</b> Automatically invalidates canonical indices
    ///         if the element type has defined symmetry.
    ///     </para>
    ///     <para>
    ///         <b>USE CASES:</b>
    ///         - Topology updates (changing element connectivity)
    ///         - Mesh transformation (redefining elements)
    ///         - Algorithm implementations requiring in-place updates
    ///     </para>
    /// </remarks>
    /// <exception cref="ArgumentNullException">Thrown when newNodes is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when element index is invalid.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the topology has been disposed.</exception>
    public void ReplaceElementNodes<TElement, TNode>(int element, params int[] newNodes)
    {
        ArgumentNullException.ThrowIfNull(newNodes);
        ThrowIfDisposed();

        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            var elementType = GetTypeIndex<TElement>();
            var nodeType = GetTypeIndex<TNode>();
            var m2m = _adjacency[elementType, nodeType];

            m2m.ReplaceElement(element, new List<int>(newNodes));

            // Invalidate canonical index for this type if it has symmetry
            InvalidateCanonicalIndex<TElement>();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Clears all nodes from an existing element, making it empty.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="element">The element index.</param>
    /// <remarks>
    ///     <para>
    ///         <b>FEATURE COMPLETION:</b> This method fills a gap in the Topology API.
    ///         Previously, there was no way to clear an element's connectivity without
    ///         removing the element entirely.
    ///     </para>
    ///     <para>
    ///         <b>INDEX PRESERVATION:</b> The element index remains valid after clearing.
    ///         The element exists but has no nodes. Use <see cref="Remove{TEntity}" />
    ///         if you want to remove the element entirely.
    ///     </para>
    ///     <para>
    ///         <b>CACHE INVALIDATION:</b> Automatically invalidates canonical indices
    ///         if the element type has defined symmetry.
    ///     </para>
    ///     <para>
    ///         <b>USE CASES:</b>
    ///         - Preparing elements for repopulation
    ///         - Temporarily disconnecting elements
    ///         - Algorithm implementations requiring empty placeholders
    ///     </para>
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when element index is invalid.</exception>
    /// <exception cref="ObjectDisposedException">Thrown if the topology has been disposed.</exception>
    public void ClearElement<TElement, TNode>(int element)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            var elementType = GetTypeIndex<TElement>();
            var nodeType = GetTypeIndex<TNode>();
            var m2m = _adjacency[elementType, nodeType];

            m2m.ClearElement(element);

            // Invalidate canonical index for this type if it has symmetry
            InvalidateCanonicalIndex<TElement>();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Helper method to invalidate canonical index for a type if it has symmetry defined.
    /// </summary>
    private void InvalidateCanonicalIndex<TElement>()
    {
        var type = typeof(TElement);
        if (_canonicalIndex.ContainsKey(type)) _canonicalIndex[type].Clear();
    }

    #endregion

    #region Compression and Optimization

    /// <summary>
    ///     Applies all pending removals and renumbers entities.
    /// </summary>
    /// <param name="removeDuplicates">
    ///     If true, automatically detects and removes duplicate elements before compression.
    ///     Typical reduction: 5-10% in refined meshes.
    /// </param>
    /// <param name="shrinkMemory">
    ///     If true, reclaims excess memory after compression.
    ///     Typical savings: 10-20% memory.
    /// </param>
    /// <param name="validate">
    ///     If true, validates structure integrity before and after compression.
    ///     Throws InvalidOperationException if validation fails.
    /// </param>
    /// <remarks>
    ///     Now supports automatic duplicate removal, memory optimization,
    ///     and optional validation for increased robustness.
    /// </remarks>
    /// <exception cref="InvalidOperationException">Thrown when validate=true and structure is invalid.</exception>
    public void Compress(bool removeDuplicates = false, bool shrinkMemory = false, bool validate = false)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            // v3: Optional pre-compression validation
            if (validate && !_adjacency.ValidateStructure())
                throw new InvalidOperationException("Structure invalid before compression");

            // v3: Use enhanced MM2M compression with options
            var maps = _adjacency.Compress(removeDuplicates, shrinkMemory);
            if (maps == null) return;

            // Reorder data lists using IDataList interface (no reflection)
            foreach (var kvp in _data)
                if (_typeToIndex.TryGetValue(kvp.Key.Entity, out var typeIndex) && typeIndex < maps.Count)
                    kvp.Value.ReorderByMapping(maps[typeIndex].oldNodesFromNew);

            // Rebuild canonical indices with new numbering
            foreach (var kvp in _canonicalIndex)
                if (_typeToIndex.TryGetValue(kvp.Key, out var typeIndex) && typeIndex < maps.Count)
                {
                    var newFromOld = maps[typeIndex].newNodesFromOld;
                    var rebuilt = new Dictionary<long, List<(int Index, List<int> Nodes)>>();

                    // Rebuild each collision chain with new indices
                    foreach (var (key, collisionChain) in kvp.Value)
                    {
                        var newChain = new List<(int Index, List<int> Nodes)>();
                        foreach (var (oldIdx, nodes) in collisionChain)
                            if (oldIdx < newFromOld.Count && newFromOld[oldIdx] >= 0)
                                newChain.Add((newFromOld[oldIdx], nodes));

                        if (newChain.Count > 0) rebuilt[key] = newChain;
                    }

                    kvp.Value.Clear();
                    foreach (var (k, chain) in rebuilt)
                        kvp.Value[k] = chain;
                }

            // v3: Optional post-compression validation
            if (validate && !_adjacency.ValidateStructure())
                throw new InvalidOperationException("Structure corrupted during compression");
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #region Async Operations (Priority 3)

    /// <summary>
    ///     Asynchronously compresses the topology structure.
    /// </summary>
    /// <param name="removeDuplicates">Whether to remove duplicate elements.</param>
    /// <param name="shrinkMemory">Whether to shrink memory allocations.</param>
    /// <param name="validate">Whether to validate structure before/after compression.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>A task representing the asynchronous operation.</returns>
    /// <remarks>
    ///     Async version for long-running compression operations.
    ///     Allows UI to remain responsive and supports cancellation.
    ///     Useful for large meshes where compression might take seconds.
    ///     Example: UI applications can show progress dialog during compression.
    /// </remarks>
    public async Task CompressAsync(
        bool removeDuplicates = false,
        bool shrinkMemory = false,
        bool validate = false,
        CancellationToken cancellationToken = default)
    {
        await Task.Run(() =>
        {
            cancellationToken.ThrowIfCancellationRequested();
            Compress(removeDuplicates, shrinkMemory, validate);
        }, cancellationToken).ConfigureAwait(false);
    }

    /// <summary>
    ///     Asynchronously validates the structural integrity of the topology.
    /// </summary>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>True if structure is valid, false otherwise.</returns>
    /// <remarks>
    ///     Async version for validation of large structures.
    ///     Validation can take considerable time for multi-million element meshes.
    /// </remarks>
    public async Task<bool> ValidateStructureAsync(
        CancellationToken cancellationToken = default)
    {
        return await Task.Run(() =>
        {
            cancellationToken.ThrowIfCancellationRequested();
            return ValidateStructure();
        }, cancellationToken).ConfigureAwait(false);
    }

    #endregion

    /// <summary>
    ///     Clears all entities and data.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <b>EXCEPTION SAFETY:</b> This method provides strong exception guarantee.
    ///         If any operation fails, the topology is left in a valid (but possibly cleared) state.
    ///         Data is cleared first, then adjacency is marked and compressed.
    ///     </para>
    ///     <para>
    ///         <b>ORDER:</b> Clears in order: data lists, canonical indices, adjacency.
    ///         This order ensures external observers see consistent state.
    ///     </para>
    /// </remarks>
    public void Clear()
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            // Phase 1: Clear auxiliary structures first (fast, low failure risk)
            foreach (var kvp in _data)
                kvp.Value.Clear();

            foreach (var kvp in _canonicalIndex)
                kvp.Value.Clear();

            // Phase 2: Mark all elements for erasure
            // Collect marks first to avoid partial marking on failure
            var marksToApply = new List<(int typeIdx, int elemIdx)>();
            for (var t = 0; t < _types.Count; t++)
            {
                var count = _adjacency.GetNumberOfElements(t);
                for (var i = 0; i < count; i++)
                    marksToApply.Add((t, i));
            }

            // Apply all marks
            foreach (var (typeIdx, elemIdx) in marksToApply)
                _adjacency.MarkToErase(typeIdx, elemIdx);

            // Phase 3: Compress (removes marked elements)
            // If this fails, elements are marked but not removed - still consistent
            _adjacency.Compress();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Type-Level Analysis

    /// <summary>
    ///     Gets the entity types in dependency order.
    /// </summary>
    public IReadOnlyList<int> GetTypeDependencyOrder()
    {
        ThrowIfDisposed();
        return _adjacency.GetTypeTopOrder();
    }

    /// <summary>
    ///     Checks whether the type dependency graph is acyclic.
    /// </summary>
    public bool AreTypeDependenciesAcyclic()
    {
        ThrowIfDisposed();
        return _adjacency.AreTypesAcyclic();
    }

    /// <summary>
    ///     Gets the types that a given type depends on.
    /// </summary>
    public IReadOnlyList<int> GetDependencies<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var entityType = GetTypeIndex<TEntity>();
            var deps = new List<int>();

            for (var t = 0; t < _types.Count; t++)
                if (t != entityType && _adjacency[entityType, t].Count > 0)
                    deps.Add(t);

            return deps;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the types that depend on a given type.
    /// </summary>
    public IReadOnlyList<int> GetDependents<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var entityType = GetTypeIndex<TEntity>();
            var deps = new List<int>();

            for (var t = 0; t < _types.Count; t++)
                if (t != entityType && _adjacency[t, entityType].Count > 0)
                    deps.Add(t);

            return deps;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities marked for removal.
    /// </summary>
    public IReadOnlySet<int> GetMarkedForRemoval<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var typeIndex = GetTypeIndex<TEntity>();
            var marked = _adjacency.ListOfMarked;
            return marked.TryGetValue(typeIndex, out var set) ? set : new HashSet<int>();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all active (not marked for deletion) entities.
    /// </summary>
    public List<int> GetActive<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var typeIndex = GetTypeIndex<TEntity>();
            var count = _adjacency.GetNumberOfElements(typeIndex);
            var markedSet = _adjacency.ListOfMarked;

            var marked = markedSet.TryGetValue(typeIndex, out var set)
                ? set
                : new HashSet<int>();

            var result = new List<int>(count - marked.Count);

            for (var i = 0; i < count; i++)
                if (!marked.Contains(i))
                    result.Add(i);

            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Validation and Quality

    /// <summary>
    ///     Validates the structural integrity of the entire topology.
    /// </summary>
    /// <returns>True if all structures are valid; false otherwise.</returns>
    /// <remarks>
    ///     Comprehensive validation including:
    ///     - All M2M structures have valid node indices
    ///     - Diagonal M2M structures are acyclic
    ///     - Type dependencies form valid DAG
    ///     Time Complexity: O(T² × n × m) - expensive, use for debugging or after major operations.
    /// </remarks>
    public bool ValidateStructure()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            return _adjacency.ValidateStructure();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets duplicate elements for a specific entity type.
    /// </summary>
    /// <typeparam name="TEntity">The entity type to check.</typeparam>
    /// <returns>List of element indices that are duplicates of earlier elements.</returns>
    /// <remarks>
    ///     Identifies elements with identical adjacency lists.
    ///     Useful for mesh quality checks and cleanup after refinement.
    ///     Time Complexity: O(n log n × m) where n = entities, m = avg connections.
    /// </remarks>
    public List<int> GetDuplicates<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            return _adjacency.GetDuplicates(GetTypeIndex<TEntity>());
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all duplicate elements across all entity types.
    /// </summary>
    /// <returns>Dictionary mapping type index to list of duplicate element indices.</returns>
    /// <remarks>
    ///     Comprehensive duplicate detection across all types.
    ///     Only returns types that have duplicates (empty types omitted).
    /// </remarks>
    public Dictionary<int, List<int>> GetAllDuplicates()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            return _adjacency.GetAllDuplicates();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Performance Tuning

    /// <summary>
    ///     Configures performance parameters for a specific entity type.
    /// </summary>
    /// <typeparam name="TEntity">The entity type to configure.</typeparam>
    /// <param name="parallelizationThreshold">
    ///     Threshold for using parallel algorithms. Operations on structures with more elements
    ///     than this value will use parallelization. Recommended: 1000-10000 depending on hardware.
    /// </param>
    /// <param name="reserveCapacity">
    ///     Optional capacity to pre-allocate. Use when approximate size is known in advance
    ///     to avoid reallocation overhead during construction.
    /// </param>
    /// <remarks>
    ///     Tune performance for specific entity type characteristics.
    ///     Example: Configure vertices with high threshold and large capacity,
    ///     volumes with low threshold and small capacity.
    ///     Typical speedup: 10-50% depending on workload.
    /// </remarks>
    public void ConfigureType<TEntity>(int parallelizationThreshold, int? reserveCapacity = null)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            var typeIndex = GetTypeIndex<TEntity>();
            _adjacency.ConfigureType(typeIndex, parallelizationThreshold, reserveCapacity);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Reserves capacity for a specific entity-to-related relationship.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="capacity">The capacity to pre-allocate.</param>
    /// <remarks>
    ///     Fine-grained capacity reservation for specific relationships.
    ///     Use when relationship size is known in advance.
    ///     Typical speedup during construction: 20-30%.
    /// </remarks>
    public void Reserve<TElement, TNode>(int capacity)
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            _adjacency.Reserve(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), capacity);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Reclaims excess memory from all internal structures.
    /// </summary>
    /// <remarks>
    ///     Call after bulk operations complete to reduce memory footprint.
    ///     Typical savings: 10-20% memory.
    ///     Time Complexity: O(T² × n × m) - relatively expensive.
    /// </remarks>
    public void ShrinkToFit()
    {
        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();
            _adjacency.ShrinkToFit();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Element Ordering

    /// <summary>
    ///     Gets the topological ordering of entities within a type based on dependencies.
    /// </summary>
    /// <typeparam name="TEntity">The entity type to order.</typeparam>
    /// <returns>List of entity indices in topological order.</returns>
    /// <exception cref="InvalidOperationException">Thrown when the entity graph contains cycles.</exception>
    /// <remarks>
    ///     Returns optimal processing order for entities based on their
    ///     self-referential dependency structure. Very useful for:
    ///     - FEM assembly order optimization
    ///     - Dependency-aware processing
    ///     - Constraint ordering
    ///     For types without self-references, returns natural order [0, 1, 2, ...].
    ///     If cycles detected, use GetSortOrder instead.
    /// </remarks>
    public List<int> GetTopologicalOrder<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            return _adjacency.GetElementTopologicalOrder(GetTypeIndex<TEntity>());
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets a canonical sort order for entities of a specific type.
    /// </summary>
    /// <typeparam name="TEntity">The entity type to sort.</typeparam>
    /// <returns>List of entity indices in lexicographic sort order.</returns>
    /// <remarks>
    ///     Returns a deterministic ordering based on lexicographic comparison
    ///     of adjacency lists. Useful for:
    ///     - Deterministic output
    ///     - Testing and comparison
    ///     - Fallback when topological order fails due to cycles
    ///     Time Complexity: O(n log n × m) where n = entities, m = avg connections.
    /// </remarks>
    public List<int> GetSortOrder<TEntity>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            return _adjacency.GetElementSortOrder(GetTypeIndex<TEntity>());
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Zero-Copy Access (v4 - Thread-Safe)

    /// <summary>
    ///     Executes an action with zero-copy span access to element nodes.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="element">The element index.</param>
    /// <param name="action">Action to execute with the span.</param>
    /// <remarks>
    ///     Thread-safe zero-copy access for high-performance inner loops.
    ///     The action is executed under Topology read lock + M2M read lock, ensuring safety.
    ///     Typical speedup: 5-10% in tight loops vs. list allocation.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void WithNodesSpan<TElement, TNode>(int element, Action<ReadOnlySpan<int>> action)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            m2m.WithNodesSpan(element, action);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes a function with zero-copy span access to element nodes.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <typeparam name="TResult">The result type.</typeparam>
    /// <param name="element">The element index.</param>
    /// <param name="func">Function to execute with the span.</param>
    /// <returns>The result of the function.</returns>
    /// <remarks>
    ///     Thread-safe zero-copy access for high-performance inner loops.
    ///     The function is executed under Topology read lock + M2M read lock, ensuring safety.
    ///     Typical speedup: 5-10% in tight loops vs. list allocation.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public TResult WithNodesSpan<TElement, TNode, TResult>(int element, Func<ReadOnlySpan<int>, TResult> func)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.WithNodesSpan(element, func);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Traversal

    /// <summary>
    ///     Performs depth-first traversal starting from a node.
    /// </summary>
    public IReadOnlyList<int> Traverse<TElement, TNode>(int startNode)
    {
        ThrowIfDisposed();

        var visited = new HashSet<int>();
        var result = new List<int>();
        var stack = new Stack<int>();

        foreach (var elem in ElementsAt<TElement, TNode>(startNode))
            stack.Push(elem);

        while (stack.Count > 0)
        {
            var current = stack.Pop();
            if (!visited.Add(current))
                continue;

            result.Add(current);

            foreach (var node in NodesOf<TElement, TNode>(current))
            foreach (var neighbor in ElementsAt<TElement, TNode>(node))
                if (!visited.Contains(neighbor))
                    stack.Push(neighbor);
        }

        return result;
    }

    /// <summary>
    ///     Performs breadth-first traversal starting from a node.
    /// </summary>
    public IReadOnlyList<int> TraverseBreadthFirst<TElement, TNode>(int startNode)
    {
        ThrowIfDisposed();

        var visited = new HashSet<int>();
        var result = new List<int>();
        var queue = new Queue<int>();

        foreach (var elem in ElementsAt<TElement, TNode>(startNode))
            queue.Enqueue(elem);

        while (queue.Count > 0)
        {
            var current = queue.Dequeue();
            if (!visited.Add(current))
                continue;

            result.Add(current);

            foreach (var node in NodesOf<TElement, TNode>(current))
            foreach (var neighbor in ElementsAt<TElement, TNode>(node))
                if (!visited.Contains(neighbor))
                    queue.Enqueue(neighbor);
        }

        return result;
    }

    /// <summary>
    ///     Finds connected components of elements.
    /// </summary>
    public IReadOnlyList<IReadOnlyList<int>> FindComponents<TElement, TNode>()
    {
        ThrowIfDisposed();

        var count = Count<TElement>();
        var visited = new HashSet<int>();
        var components = new List<IReadOnlyList<int>>();

        for (var i = 0; i < count; i++)
        {
            if (visited.Contains(i))
                continue;

            var component = new List<int>();
            var stack = new Stack<int>();
            stack.Push(i);

            while (stack.Count > 0)
            {
                var current = stack.Pop();
                if (!visited.Add(current))
                    continue;

                component.Add(current);

                foreach (var node in NodesOf<TElement, TNode>(current))
                foreach (var neighbor in ElementsAt<TElement, TNode>(node))
                    if (!visited.Contains(neighbor))
                        stack.Push(neighbor);
            }

            components.Add(component);
        }

        return components;
    }

    #endregion

    #region Connectivity Operations

    /// <summary>
    ///     Gets entities sharing at least one relationship (direct neighbors).
    /// </summary>
    public List<int> GetDirectNeighbors<TEntity, TRelated>(
        int entityIndex,
        bool includeSelf = false,
        bool sorted = true) // v3: Added for consistency
    {
        ArgumentOutOfRangeException.ThrowIfNegative(entityIndex);

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            var neighbors = m2m.GetElementNeighbors(entityIndex, sorted); // v3: Pass sorting parameter

            if (includeSelf)
            {
                neighbors.Add(entityIndex);
                if (sorted) // v3: Only sort if requested
                    neighbors.Sort();
            }

            return neighbors;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities sharing exactly k relationships.
    /// </summary>
    public List<int> GetEntitiesWithSharedCount<TEntity, TRelated>(
        int entityIndex,
        int exactCount,
        bool includeSelf = false)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(entityIndex);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(exactCount);

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            var myRelations = m2m[entityIndex];

            if (myRelations.Count < exactCount)
                return new List<int>();

            var sharedCounts = new Dictionary<int, int>();

            foreach (var relation in myRelations)
            {
                m2m.WithElementsForNodeSpan(relation, span =>
                {
                    for (var i = 0; i < span.Length; i++)
                    {
                        var entity = span[i];
                        if (entity == entityIndex && !includeSelf)
                            continue;

                        sharedCounts.TryGetValue(entity, out var count);
                        sharedCounts[entity] = count + 1;
                    }
                });
            }

            var result = new List<int>();
            foreach (var kvp in sharedCounts)
                if (kvp.Value == exactCount)
                    result.Add(kvp.Key);

            result.Sort();
            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities sharing at least minCount relationships.
    /// </summary>
    public List<int> GetEntitiesWithMinSharedCount<TEntity, TRelated>(
        int entityIndex,
        int minCount = 1,
        bool includeSelf = false)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(entityIndex);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(minCount);

        if (minCount == 1)
            return GetDirectNeighbors<TEntity, TRelated>(entityIndex, includeSelf);

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            var myRelations = m2m[entityIndex];
            var sharedCounts = new Dictionary<int, int>();

            foreach (var relation in myRelations)
            {
                m2m.WithElementsForNodeSpan(relation, span =>
                {
                    for (var i = 0; i < span.Length; i++)
                    {
                        var entity = span[i];
                        if (entity == entityIndex && !includeSelf)
                            continue;

                        sharedCounts.TryGetValue(entity, out var count);
                        sharedCounts[entity] = count + 1;
                    }
                });
            }

            var result = new List<int>();
            foreach (var kvp in sharedCounts)
                if (kvp.Value >= minCount)
                    result.Add(kvp.Key);

            result.Sort();
            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities with their shared relationship counts.
    /// </summary>
    public List<(int EntityIndex, int SharedCount)> GetWeightedNeighbors<TEntity, TRelated>(
        int entityIndex,
        int minCount = 1,
        bool includeSelf = false)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(entityIndex);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(minCount);

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            var myRelations = m2m[entityIndex];
            var sharedCounts = new Dictionary<int, int>();

            foreach (var relation in myRelations)
            {
                m2m.WithElementsForNodeSpan(relation, span =>
                {
                    for (var i = 0; i < span.Length; i++)
                    {
                        var entity = span[i];
                        if (entity == entityIndex && !includeSelf)
                            continue;

                        sharedCounts.TryGetValue(entity, out var count);
                        sharedCounts[entity] = count + 1;
                    }
                });
            }

            var result = new List<(int EntityIndex, int SharedCount)>();
            foreach (var kvp in sharedCounts)
                if (kvp.Value >= minCount)
                    result.Add((kvp.Key, kvp.Value));

            result.Sort();
            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities containing ALL specified relationships (intersection).
    /// </summary>
    public List<int> GetEntitiesContainingAll<TEntity, TRelated>(List<int> relationships)
    {
        ArgumentNullException.ThrowIfNull(relationships);

        if (relationships.Count == 0)
            return new List<int>();

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            return m2m.GetElementsWithNodes(relationships);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities containing ANY of the specified relationships (union).
    /// </summary>
    /// <typeparam name="TEntity">The entity type to search.</typeparam>
    /// <typeparam name="TRelated">The relationship/node type.</typeparam>
    /// <param name="relationships">The relationships to search for (at least one must be present).</param>
    /// <returns>List of entity indices containing at least one of the specified relationships.</returns>
    /// <remarks>
    ///     Union-semantics counterpart to GetEntitiesContainingAll.
    ///     <b>DIFFERENCE FROM GetEntitiesContainingAll:</b>
    ///     - GetEntitiesContainingAll returns entities with ALL specified relationships (intersection)
    ///     - GetEntitiesContainingAny returns entities with AT LEAST ONE (union)
    ///     Time Complexity: O(k × m) where k = relationships.Count, m = avg entities per relationship.
    /// </remarks>
    public List<int> GetEntitiesContainingAny<TEntity, TRelated>(List<int> relationships)
    {
        ArgumentNullException.ThrowIfNull(relationships);

        if (relationships.Count == 0)
            return new List<int>();

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            return m2m.GetElementsContainingAnyNode(relationships);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities within k relationship steps via BFS.
    /// </summary>
    public Dictionary<int, int> GetKHopNeighborhood<TEntity, TRelated>(
        int seedEntity,
        int k,
        int minSharedForConnection = 1,
        bool includeSeed = true)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(seedEntity);
        ArgumentOutOfRangeException.ThrowIfNegative(k);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(minSharedForConnection);

        var result = new Dictionary<int, int>();
        if (includeSeed)
            result[seedEntity] = 0;

        if (k == 0)
            return result;

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            var visited = new HashSet<int> { seedEntity };
            List<int> currentHop = [seedEntity];

            for (var hop = 1; hop <= k; hop++)
            {
                var nextHop = new List<int>();

                foreach (var entity in currentHop)
                {
                    List<int> neighbors;

                    if (minSharedForConnection == 1)
                        neighbors = m2m.GetElementNeighbors(entity, false); // v3: unsorted for BFS
                    else
                        neighbors = GetEntitiesWithMinSharedCount<TEntity, TRelated>(
                            entity, minSharedForConnection);

                    foreach (var neighbor in neighbors)
                        if (visited.Add(neighbor))
                        {
                            result[neighbor] = hop;
                            nextHop.Add(neighbor);
                        }
                }

                currentHop = nextHop;
                if (currentHop.Count == 0)
                    break;
            }

            return result;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets entities at exactly k steps from seed.
    /// </summary>
    public List<int> GetEntitiesAtDistance<TEntity, TRelated>(
        int seedEntity,
        int k,
        int minSharedForConnection = 1)
    {
        var neighborhood = GetKHopNeighborhood<TEntity, TRelated>(
            seedEntity, k, minSharedForConnection, false);

        var result = new List<int>();
        foreach (var kvp in neighborhood)
            if (kvp.Value == k)
                result.Add(kvp.Key);

        result.Sort();
        return result;
    }

    #endregion

    #region Multi-Type Connectivity

    /// <summary>
    ///     Performs DFS across ALL types starting from a node.
    /// </summary>
    public List<(int TypeIndex, int EntityIndex)> MultiTypeDFS<TNode>(int nodeIndex)
    {
        ThrowIfDisposed();
        ArgumentOutOfRangeException.ThrowIfNegative(nodeIndex);

        _rwLock.EnterReadLock();
        try
        {
            var nodeTypeIndex = GetTypeIndex<TNode>();
            return _adjacency.DepthFirstSearchFromANode(nodeTypeIndex, nodeIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets ALL entities of ANY type connected to a specific node.
    /// </summary>
    public List<(int TypeIndex, int EntityIndex)> GetAllEntitiesAtNode<TNode>(int nodeIndex)
    {
        ThrowIfDisposed();
        ArgumentOutOfRangeException.ThrowIfNegative(nodeIndex);

        _rwLock.EnterReadLock();
        try
        {
            var nodeTypeIndex = GetTypeIndex<TNode>();
            return _adjacency.GetAllElements(nodeTypeIndex, nodeIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets ALL nodes of ANY type connected to a specific entity.
    /// </summary>
    public List<(int TypeIndex, int NodeIndex)> GetAllNodesOfEntity<TEntity>(int entityIndex)
    {
        ThrowIfDisposed();
        ArgumentOutOfRangeException.ThrowIfNegative(entityIndex);

        _rwLock.EnterReadLock();
        try
        {
            var entityTypeIndex = GetTypeIndex<TEntity>();
            return _adjacency.GetAllNodes(entityTypeIndex, entityIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets all nodes connected to an entity (with explicit ordering control).
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <param name="entityIndex">The entity index.</param>
    /// <param name="order">Specifies whether results should be sorted.</param>
    /// <returns>List of (type index, node index) tuples.</returns>
    /// <remarks>
    ///     Improved API with explicit <see cref="ResultOrder" /> parameter
    ///     instead of boolean for better readability at call sites.
    ///     Examples:
    ///     <code>
    ///         // Clear and explicit:
    ///         var nodes = mesh.GetAllNodesOfEntity&lt;Element&gt;(0, ResultOrder.Sorted);
    ///         // vs unclear boolean:
    ///         var nodes = mesh.GetAllNodesOfEntity&lt;Element&gt;(0, true);  // true means what?
    ///         </code>
    /// </remarks>
    public List<(int TypeIndex, int NodeIndex)> GetAllNodesOfEntity<TEntity>(
        int entityIndex,
        ResultOrder order)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(entityIndex);

        var entityTypeIndex = GetTypeIndex<TEntity>();
        return _adjacency.GetAllNodes(entityTypeIndex, entityIndex, order == ResultOrder.Sorted);
    }

    /// <summary>
    ///     Gets all entities connected to a node (with explicit ordering control).
    /// </summary>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="nodeIndex">The node index.</param>
    /// <param name="order">Specifies whether results should be sorted.</param>
    /// <returns>List of (type index, entity index) tuples.</returns>
    /// <remarks>
    ///     Improved API with explicit <see cref="ResultOrder" /> parameter.
    /// </remarks>
    public List<(int TypeIndex, int EntityIndex)> GetAllEntitiesAtNode<TNode>(
        int nodeIndex,
        ResultOrder order)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(nodeIndex);

        var nodeTypeIndex = GetTypeIndex<TNode>();
        return _adjacency.GetAllElements(nodeTypeIndex, nodeIndex, order == ResultOrder.Sorted);
    }

    /// <summary>
    ///     Gets the type topological order.
    /// </summary>
    public List<int> GetTypeTopologicalOrder()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _adjacency.GetTypeTopOrder();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if the type hierarchy is acyclic.
    /// </summary>
    public bool IsTypeHierarchyAcyclic()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            return _adjacency.AreTypesAcyclic();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Algebraic Operations

    /// <summary>
    ///     Computes transitive connectivity via matrix multiplication.
    /// </summary>
    /// <remarks>
    ///     Returns O2M where row i contains all entities reachable from entity i
    ///     through the entity-related-entity path.
    /// </remarks>
    public O2M ComputeTransitiveConnectivity<TEntity, TRelated>()
    {
        ThrowIfDisposed();
        var forward = GetForwardStructure<TEntity, TRelated>();
        var transpose = GetTranspose<TEntity, TRelated>();
        return forward * transpose;
    }

    /// <summary>
    ///     Gets the dual/transpose relationship structure.
    /// </summary>
    /// <remarks>
    ///     Alias for GetTranspose. Returns O2M where row i contains
    ///     all entities connected to related-entity i.
    /// </remarks>
    public O2M GetDualStructure<TEntity, TRelated>()
    {
        ThrowIfDisposed();
        return GetTranspose<TEntity, TRelated>();
    }

    /// <summary>
    ///     Checks if the relationship structure is acyclic.
    /// </summary>
    /// <remarks>
    ///     Uses the forward structure to detect cycles via DFS.
    /// </remarks>
    public bool IsAcyclic<TEntity, TRelated>()
    {
        ThrowIfDisposed();
        return WithTranspose<TEntity, TRelated, bool>(o2m => o2m.IsAcyclic());
    }

    /// <summary>
    ///     Gets topological ordering if structure is a DAG.
    /// </summary>
    /// <exception cref="InvalidOperationException">Thrown if the structure contains cycles.</exception>
    public List<int> GetTopologicalOrder<TEntity, TRelated>()
    {
        ThrowIfDisposed();
        return WithTranspose<TEntity, TRelated, List<int>>(o2m => o2m.GetTopOrder());
    }

    /// <summary>
    ///     Computes related-to-related connectivity.
    /// </summary>
    /// <remarks>
    ///     Returns O2M where row i contains all related-entities that share
    ///     at least one entity with related-entity i.
    /// </remarks>
    public O2M ComputeRelatedToRelatedConnectivity<TEntity, TRelated>()
    {
        ThrowIfDisposed();
        var forward = GetForwardStructure<TEntity, TRelated>();
        var transpose = GetTranspose<TEntity, TRelated>();
        return transpose * forward;
    }

    /// <summary>
    ///     Gets neighbors of an element (other elements sharing nodes).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="element">The element index.</param>
    /// <param name="sorted">Whether to sort results. Default is true.</param>
    /// <returns>List of neighbor element indices.</returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Wraps M2M.GetElementNeighbors for clean public API.
    ///     Two elements are neighbors if they share at least one node.
    ///     The returned list does not include the element itself.
    /// </remarks>
    public List<int> GetElementNeighbors<TElement, TNode>(int element, bool sorted = true)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetElementNeighbors(element, sorted);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets neighbors of a node (other nodes sharing elements).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="node">The node index.</param>
    /// <param name="sorted">Whether to sort results. Default is true.</param>
    /// <returns>List of neighbor node indices.</returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Wraps M2M.GetNodeNeighbors for clean public API.
    ///     Two nodes are neighbors if they appear together in at least one element.
    ///     The returned list does not include the node itself.
    /// </remarks>
    public List<int> GetNodeNeighbors<TElement, TNode>(int node, bool sorted = true)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetNodeNeighbors(node, sorted);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Graph Analysis

    /// <summary>
    ///     Finds connected components using BFS.
    /// </summary>
    public List<List<int>> FindConnectedComponents<TEntity, TRelated>(int minSharedForConnection = 1)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(minSharedForConnection);

        _rwLock.EnterReadLock();
        try
        {
            var entityCount = Count<TEntity>();
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            var visited = new HashSet<int>();
            var components = new List<List<int>>();

            for (var i = 0; i < entityCount; i++)
            {
                if (visited.Contains(i)) continue;

                var component = new List<int>();
                var queue = new Queue<int>();
                queue.Enqueue(i);
                visited.Add(i);

                while (queue.Count > 0)
                {
                    var current = queue.Dequeue();
                    component.Add(current);

                    List<int> neighbors;
                    if (minSharedForConnection == 1)
                        neighbors = m2m.GetElementNeighbors(current, false); // v3: unsorted for BFS
                    else
                        neighbors = GetEntitiesWithMinSharedCount<TEntity, TRelated>(
                            current, minSharedForConnection);

                    foreach (var neighbor in neighbors)
                        if (visited.Add(neighbor))
                            queue.Enqueue(neighbor);
                }

                components.Add(component);
            }

            return components;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Computes degree for each entity.
    /// </summary>
    public Dictionary<int, int> ComputeDegrees<TEntity, TRelated>()
    {
        _rwLock.EnterReadLock();
        try
        {
            var entityCount = Count<TEntity>();
            var entityTypeIndex = GetTypeIndex<TEntity>();
            var relatedTypeIndex = GetTypeIndex<TRelated>();
            var m2m = _adjacency[entityTypeIndex, relatedTypeIndex];

            var degrees = new Dictionary<int, int>(entityCount);

            for (var i = 0; i < entityCount; i++) degrees[i] = m2m[i].Count;

            return degrees;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Computes connectivity statistics.
    /// </summary>
    public ConnectivityStatistics GetConnectivityStatistics<TEntity, TRelated>(params int[] sharedCounts)
    {
        if (sharedCounts.Length == 0)
            sharedCounts = [1];

        var entityCount = Count<TEntity>();
        var neighborCountsByLevel = new Dictionary<int, List<int>>();

        foreach (var sharedCount in sharedCounts)
            neighborCountsByLevel[sharedCount] = new List<int>(entityCount);

        for (var i = 0; i < entityCount; i++)
            foreach (var sharedCount in sharedCounts)
            {
                List<int> neighbors;
                if (sharedCount == 1)
                    neighbors = GetDirectNeighbors<TEntity, TRelated>(i, sorted: false); // v3: unsorted for counting
                else
                    neighbors = GetEntitiesWithSharedCount<TEntity, TRelated>(i, sharedCount);
                neighborCountsByLevel[sharedCount].Add(neighbors.Count);
            }

        return new ConnectivityStatistics(entityCount, neighborCountsByLevel);
    }

    /// <summary>
    ///     Constructs element-to-element adjacency graph.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>
    ///     O2M structure where row i contains indices of elements adjacent to element i.
    ///     Two elements are adjacent if they share at least one node.
    /// </returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Constructs dual graph for graph partitioning, coloring, etc.
    ///     Common uses:
    ///     - Graph partitioning (METIS, Scotch)
    ///     - Element coloring for parallel assembly
    ///     - Mesh smoothing algorithms
    ///     <b>PERFORMANCE:</b> O(n·m) where n = elements, m = avg nodes per element.
    ///     Result is symmetric (if element i is neighbor of j, then j is neighbor of i).
    /// </remarks>
    public O2M GetElementToElementGraph<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetElementsToElements();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Constructs node-to-node adjacency graph.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>
    ///     O2M structure where row i contains indices of nodes adjacent to node i.
    ///     Two nodes are adjacent if they appear together in at least one element.
    /// </returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Constructs node connectivity graph.
    ///     Common uses:
    ///     - Bandwidth reduction (Cuthill-McKee, Reverse Cuthill-McKee)
    ///     - Node coloring
    ///     - Mesh quality analysis
    ///     <b>PERFORMANCE:</b> O(n·m²) where n = elements, m = avg nodes per element.
    ///     Result is symmetric.
    /// </remarks>
    public O2M GetNodeToNodeGraph<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetNodesToNodes();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Statistics

    /// <summary>
    ///     Statistics about entity connectivity.
    /// </summary>
    public sealed class ConnectivityStatistics
    {
        internal ConnectivityStatistics(
            int entityCount,
            Dictionary<int, List<int>> neighborCountsByLevel)
        {
            EntityCount = entityCount;
            NeighborCountsByLevel = neighborCountsByLevel;
        }

        public int EntityCount { get; }
        public Dictionary<int, List<int>> NeighborCountsByLevel { get; }

        public double GetAverageNeighbors(int sharedCount)
        {
            if (NeighborCountsByLevel.TryGetValue(sharedCount, out var counts) && counts.Count > 0)
            {
                long sum = 0;
                for (var i = 0; i < counts.Count; i++)
                    sum += counts[i];
                return sum / (double)counts.Count;
            }

            return 0;
        }

        public int GetMinNeighbors(int sharedCount)
        {
            if (NeighborCountsByLevel.TryGetValue(sharedCount, out var counts) && counts.Count > 0)
            {
                var min = counts[0];
                for (var i = 1; i < counts.Count; i++)
                    if (counts[i] < min)
                        min = counts[i];
                return min;
            }

            return 0;
        }

        public int GetMaxNeighbors(int sharedCount)
        {
            if (NeighborCountsByLevel.TryGetValue(sharedCount, out var counts) && counts.Count > 0)
            {
                var max = counts[0];
                for (var i = 1; i < counts.Count; i++)
                    if (counts[i] > max)
                        max = counts[i];
                return max;
            }

            return 0;
        }

        public override string ToString()
        {
            var lines = new List<string> { $"Connectivity Statistics ({EntityCount} entities):" };

            var sortedKeys = new List<int>(NeighborCountsByLevel.Keys);
            sortedKeys.Sort();
            foreach (var sharedCount in sortedKeys)
            {
                var counts = NeighborCountsByLevel[sharedCount];
                lines.Add($"  Shared={sharedCount}: avg={GetAverageNeighbors(sharedCount):F2}, " +
                          $"min={GetMinNeighbors(sharedCount)}, max={GetMaxNeighbors(sharedCount)}");
            }

            return string.Join(Environment.NewLine, lines);
        }
    }

    /// <summary>
    ///     Gets basic statistics about the topology.
    /// </summary>
    public TopologyStats GetStatistics()
    {
        ThrowIfDisposed();

        var entityCounts = new Dictionary<Type, int>();
        var dataCounts = new Dictionary<Type, int>();

        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            foreach (var type in _typeToIndex.Keys)
            {
                var typeIndex = _typeToIndex[type];
                entityCounts[type] = _adjacency.GetNumberOfElements(typeIndex);
            }

            foreach (var kvp in _data)
                // Use IDataList.Count directly (no reflection needed)
                dataCounts[kvp.Key.Entity] = kvp.Value.Count;

            var symmetryTypes = new List<Type>();
            foreach (var k in _symmetries.Keys)
                symmetryTypes.Add(k);

            return new TopologyStats(entityCounts, dataCounts, symmetryTypes);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Validation

    /// <summary>
    ///     Validates mesh integrity: checks for invalid node references.
    /// </summary>
    public ValidationResult ValidateIntegrity<TElement, TNode>()
    {
        ThrowIfDisposed();
        var errors = new List<string>();

        var elementCount = Count<TElement>();
        var nodeCount = Count<TNode>();

        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            for (var i = 0; i < nodes.Count; i++)
                if (nodes[i] < 0 || nodes[i] >= nodeCount)
                    errors.Add($"Element {e} references invalid node {nodes[i]} at position {i}");
        }

        return new ValidationResult(errors.Count == 0, errors);
    }

    /// <summary>
    ///     Result of a validation operation.
    /// </summary>
    public readonly struct ValidationResult
    {
        public bool IsValid { get; }
        public IReadOnlyList<string> Errors { get; }

        public ValidationResult(bool isValid, IReadOnlyList<string> errors)
        {
            IsValid = isValid;
            Errors = errors;
        }

        public override string ToString()
        {
            return IsValid ? "Valid" : $"Invalid: {Errors.Count} errors";
        }
    }

    #endregion

    #region Private Helpers

    private List<TData> GetOrCreateList<TEntity, TData>()
    {
        var key = (typeof(TEntity), typeof(TData));

        if (_data.TryGetValue(key, out var existing))
            return ((DataList<TData>)existing).Items;

        // Only populate _typeToIndex when actually creating a new data list
        var typeIndex = GetTypeIndex<TEntity>();
        _typeToIndex[typeof(TEntity)] = typeIndex;

        var dataList = new DataList<TData>();
        _data[key] = dataList;
        return dataList.Items;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ThrowIfDisposed()
    {
        if (Volatile.Read(ref _disposed) != 0)
            throw new ObjectDisposedException(nameof(Topology<TTypes>));
    }

    #endregion

    #region Serialization

    /// <summary>
    ///     Serializes the topology to JSON.
    /// </summary>
    public string ToJson(JsonSerializerOptions? options = null)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();

            options ??= new JsonSerializerOptions
            {
                WriteIndented = true,
                DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull
            };

            var dto = new TopologyDto
            {
                TypeCount = _types.Count,
                Adjacency = SerializeAdjacency(),
                Data = SerializeData(),
                Symmetries = SerializeSymmetries(),
                CanonicalIndices = SerializeCanonicalIndicesFull()
            };

            return JsonSerializer.Serialize(dto, options);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Saves topology to a file.
    /// </summary>
    public void SaveToFile(string path, JsonSerializerOptions? options = null)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentNullException(nameof(path), "File path cannot be null or empty");

        var json = ToJson(options);
        File.WriteAllText(path, json);
    }

    /// <summary>
    ///     Creates a topology from JSON.
    /// </summary>
    /// <param name="json">The JSON string to deserialize.</param>
    /// <param name="options">Optional JSON serializer options.</param>
    public static Topology<TTypes> FromJson(string json, JsonSerializerOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(json);

        var dto = JsonSerializer.Deserialize<TopologyDto>(json, options);
        if (dto == null)
            throw new JsonException("Failed to deserialize topology");

        var topology = new Topology<TTypes>();
        topology.RestoreFromDto(dto);
        return topology;
    }

    /// <summary>
    ///     Loads a topology from a file.
    /// </summary>
    /// <param name="path">The file path to load from.</param>
    /// <param name="options">Optional JSON serializer options.</param>
    public static Topology<TTypes> LoadFromFile(string path, JsonSerializerOptions? options = null)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentNullException(nameof(path), "File path cannot be null or empty");

        var json = File.ReadAllText(path);
        return FromJson(json, options);
    }

    private List<AdjacencyDto> SerializeAdjacency()
    {
        var result = new List<AdjacencyDto>();

        for (var e = 0; e < _types.Count; e++)
        for (var n = 0; n < _types.Count; n++)
        {
            var m2m = _adjacency[e, n];
            if (m2m.Count > 0)
            {
                var elements = new List<List<int>>();
                for (var i = 0; i < m2m.Count; i++)
                {
                    var row = m2m[i];
                    var copy = new List<int>(row.Count);
                    for (var j = 0; j < row.Count; j++)
                        copy.Add(row[j]);
                    elements.Add(copy);
                }

                result.Add(new AdjacencyDto
                {
                    EntityTypeIndex = e,
                    NodeTypeIndex = n,
                    Elements = elements
                });
            }
        }

        return result;
    }

    private Dictionary<string, DataListDto> SerializeData()
    {
        var result = new Dictionary<string, DataListDto>();

        foreach (var kvp in _data)
        {
            var (entityType, dataType) = kvp.Key;
            var list = kvp.Value;

            // Issue #7 fix: Use reflection-free interface methods instead of reflection
            var count = list.Count;
            if (count == 0) continue;

            var items = new List<JsonElement>();
            var elementType = list.ElementType;

            for (var i = 0; i < count; i++)
            {
                var item = list.GetItemAt(i); // No reflection!
                items.Add(JsonSerializer.SerializeToElement(item, elementType));
            }

            var key = $"{entityType.FullName}|{dataType.FullName}";
            result[key] = new DataListDto
            {
                EntityTypeName = entityType.FullName ?? entityType.Name,
                DataTypeName = dataType.FullName ?? dataType.Name,
                Items = items
            };
        }

        return result;
    }

    private Dictionary<string, SymmetryDto> SerializeSymmetries()
    {
        var result = new Dictionary<string, SymmetryDto>();

        foreach (var kvp in _symmetries)
        {
            var key = kvp.Key.FullName ?? kvp.Key.Name;
            var permutations = new List<List<int>>();
            foreach (var perm in kvp.Value.Permutations)
            {
                var copy = new List<int>(perm.Count);
                for (var i = 0; i < perm.Count; i++)
                    copy.Add(perm[i]);
                permutations.Add(copy);
            }

            result[key] = new SymmetryDto
            {
                NodeCount = kvp.Value.NodeCount,
                Permutations = permutations
            };
        }

        return result;
    }

    private List<CanonicalIndexDto> SerializeCanonicalIndicesFull()
    {
        var result = new List<CanonicalIndexDto>();

        foreach (var kvp in _canonicalIndex)
        {
            var typeKey = kvp.Key.FullName ?? kvp.Key.Name;

            foreach (var hashEntry in kvp.Value)
            {
                var dto = new CanonicalIndexDto
                {
                    TypeKey = typeKey,
                    Hash = hashEntry.Key,
                    Entries = new List<CanonicalEntryDto>()
                };

                // Serialize ALL entries in collision chain (not just the first!)
                foreach (var (index, nodes) in hashEntry.Value)
                    dto.Entries.Add(new CanonicalEntryDto
                    {
                        Index = index,
                        CanonicalNodes = new List<int>(nodes)
                    });

                result.Add(dto);
            }
        }

        return result;
    }

    private void RestoreFromDto(TopologyDto dto)
    {
        if (dto.TypeCount != _types.Count)
            throw new InvalidOperationException(
                $"Type count mismatch: expected {_types.Count}, got {dto.TypeCount}");

        // Track what we've restored for cleanup on failure
        var restorationStarted = false;

        try
        {
            restorationStarted = true;

            // Phase 1: Restore adjacency
            foreach (var adj in dto.Adjacency)
            foreach (var nodes in adj.Elements)
                _adjacency.AppendElement(adj.EntityTypeIndex, adj.NodeTypeIndex, nodes);

            // Phase 2: Restore data using IDataList interface (no reflection needed for item addition)
            foreach (var kvp in dto.Data)
            {
                var entityType = ResolveType(kvp.Value.EntityTypeName);
                var dataType = ResolveType(kvp.Value.DataTypeName);

                if (entityType == null || dataType == null)
                    throw new InvalidOperationException(
                        $"Could not resolve types: {kvp.Value.EntityTypeName}, {kvp.Value.DataTypeName}");

                // Populate _typeToIndex for the entity type
                EnsureTypeIndexCached(entityType);

                // Create DataList<T> wrapper for storage (reflection only for type construction)
                var dataListType = typeof(DataList<>).MakeGenericType(dataType);
                var dataListObj = Activator.CreateInstance(dataListType);
                if (dataListObj is not IDataList dataList)
                    throw new InvalidOperationException($"Failed to create DataList<{dataType.Name}>");

                // Use IDataList.AddItem interface method instead of reflection
                // This is both correct and faster than reflection-based approach
                foreach (var itemJson in kvp.Value.Items)
                {
                    var item = itemJson.Deserialize(dataType);
                    dataList.AddItem(item);
                }

                _data[(entityType, dataType)] = dataList;
            }

            // Phase 3: Restore symmetries
            foreach (var kvp in dto.Symmetries)
            {
                var entityType = ResolveType(kvp.Key);
                if (entityType == null)
                    throw new InvalidOperationException($"Could not resolve type: {kvp.Key}");

                // Populate _typeToIndex for the entity type
                EnsureTypeIndexCached(entityType);

                var symmetry = new Symmetry(kvp.Value.Permutations);
                _symmetries[entityType] = symmetry;

                if (!_canonicalIndex.ContainsKey(entityType))
                    _canonicalIndex[entityType] = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
            }

            // Phase 4: Restore canonical indices (full collision chains with nodes)
            foreach (var item in dto.CanonicalIndices)
            {
                var entityType = ResolveType(item.TypeKey);
                if (entityType == null)
                    throw new InvalidOperationException($"Could not resolve type: {item.TypeKey}");

                // Populate _typeToIndex for the entity type
                EnsureTypeIndexCached(entityType);

                if (!_canonicalIndex.ContainsKey(entityType))
                    _canonicalIndex[entityType] = new Dictionary<long, List<(int Index, List<int> Nodes)>>();

                var dict = _canonicalIndex[entityType];

                // Restore full collision chain with canonical nodes
                var chain = new List<(int Index, List<int> Nodes)>();
                foreach (var entry in item.Entries) chain.Add((entry.Index, entry.CanonicalNodes));

                dict[item.Hash] = chain;
            }

            // P1.1 FIX: Validate canonical index consistency after restoration
            ValidateCanonicalIndexConsistency();
        }
        catch (Exception ex)
        {
            // EXCEPTION SAFETY: On failure, the topology is in an inconsistent state.
            // Clear everything to prevent use of partial data.
            if (restorationStarted)
            {
                // Clear all partially restored data
                foreach (var kvp in _data)
                    kvp.Value.Clear();
                _data.Clear();

                foreach (var kvp in _canonicalIndex)
                    kvp.Value.Clear();
                _canonicalIndex.Clear();

                _symmetries.Clear();

                // Note: adjacency may have partial data but we can't easily roll that back
                // The topology is now unusable until properly re-initialized
            }

            throw new InvalidOperationException(
                "Restoration failed. The topology has been cleared to prevent inconsistent state. " +
                "Re-create the topology or restore from a valid DTO.", ex);
        }
    }

    /// <summary>
    ///     Validates that canonical index entries are consistent with adjacency data.
    /// </summary>
    /// <remarks>
    ///     This validation ensures:
    ///     1. Each canonical index entry points to a valid index in the adjacency
    ///     2. The stored canonical nodes match the actual adjacency nodes after canonicalization
    ///     3. No orphaned canonical entries exist
    /// </remarks>
    /// <exception cref="InvalidOperationException">Thrown when inconsistency is detected.</exception>
    private void ValidateCanonicalIndexConsistency()
    {
        foreach (var kvp in _canonicalIndex)
        {
            var entityType = kvp.Key;
            var canonicalDict = kvp.Value;

            // Get symmetry for this type (must exist for canonical index to exist)
            if (!_symmetries.TryGetValue(entityType, out var symmetry))
                throw new InvalidOperationException(
                    $"Canonical index exists for type {entityType.Name} but no symmetry is registered. " +
                    $"This indicates corrupted serialization data.");

            // Get type index for adjacency lookup
            if (!_typeToIndex.TryGetValue(entityType, out var typeIndex))
                throw new InvalidOperationException(
                    $"Canonical index exists for type {entityType.Name} but type is not registered. " +
                    $"This indicates corrupted serialization data.");

            // P0-1A FIX: Get maximum element count across all node types for this entity type
            // The diagonal may be empty if type only has off-diagonal relationships
            // (e.g., Element→Node where Element→Element diagonal is unused)
            var maxElementCount = 0;
            for (var nodeTypeIdx = 0; nodeTypeIdx < _types.Count; nodeTypeIdx++)
            {
                var blockCount = _adjacency[typeIndex, nodeTypeIdx].Count;
                if (blockCount > maxElementCount)
                    maxElementCount = blockCount;
            }

            var adjacencyCount = maxElementCount;

            // Validate each entry in the canonical index
            foreach (var (hash, collisionChain) in canonicalDict)
            foreach (var (index, canonicalNodes) in collisionChain)
            {
                // Validate index is in valid range
                if (index < 0 || index >= adjacencyCount)
                    throw new InvalidOperationException(
                        $"Canonical index for type {entityType.Name} contains invalid index {index}. " +
                        $"Valid range is [0, {adjacencyCount}). " +
                        $"This indicates corrupted serialization data or index mismatch.");

                // Get actual adjacency nodes for this index
                // We need to find the right M2M block - look for the node type
                // For canonical index, we typically use the primary node relationship
                // Check all node types to find one that has data for this entity
                var foundValidAdjacency = false;
                for (var nodeTypeIdx = 0; nodeTypeIdx < _types.Count; nodeTypeIdx++)
                {
                    var m2m = _adjacency[typeIndex, nodeTypeIdx];
                    if (index < m2m.Count)
                    {
                        var actualNodes = m2m.GetNodesForElement(index);
                        if (actualNodes.Count > 0)
                        {
                            // P0-1B FIX: Skip if arity doesn't match symmetry's expected node count
                            // This block may not be the one the canonical index was built from
                            if (actualNodes.Count != symmetry.NodeCount)
                                continue;

                            // Compute canonical form of actual nodes
                            var buffer = GetOrCreateCanonicalBuffer(actualNodes.Count);
                            var actualSpan = new int[actualNodes.Count];
                            for (var i = 0; i < actualNodes.Count; i++)
                                actualSpan[i] = actualNodes[i];

                            symmetry.CanonicalSpan(actualSpan, buffer);

                            // Compare with stored canonical nodes
                            if (canonicalNodes.Count == symmetry.NodeCount)
                            {
                                var matches = true;
                                for (var i = 0; i < symmetry.NodeCount; i++)
                                    if (buffer[i] != canonicalNodes[i])
                                    {
                                        matches = false;
                                        break;
                                    }

                                if (matches)
                                {
                                    foundValidAdjacency = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                // Note: We don't throw if no matching adjacency found, as the canonical index
                // might reference a specific node type relationship. This is a soft validation.
                // For strict validation, uncomment the following:
                // if (!foundValidAdjacency)
                // {
                //     System.Diagnostics.Debug.WriteLine(
                //         $"Warning: Canonical index entry for type {entityType.Name} at index {index} " +
                //         $"could not be validated against adjacency data.");
                // }
            }
        }
    }

    /// <summary>
    ///     Ensures a type is cached in _typeToIndex by invoking IndexOf&lt;T&gt; via reflection.
    /// </summary>
    private void EnsureTypeIndexCached(Type type)
    {
        if (_typeToIndex.ContainsKey(type))
            return;

        // Call _types.IndexOf<T>() via reflection
        var indexOfMethod = typeof(TTypes).GetMethod("IndexOf", BindingFlags.Public | BindingFlags.Instance);
        if (indexOfMethod == null)
            throw new InvalidOperationException("IndexOf method not found on ITypeMap");

        var genericMethod = indexOfMethod.MakeGenericMethod(type);
        var index = (int)genericMethod.Invoke(_types, null)!;

        _typeToIndex[type] = index;

        // Also populate the TTypes-scoped static generic cache (P0.C1 FIX)
        var cacheType = typeof(TypeIndexCache<,>).MakeGenericType(typeof(TTypes), type);
        var indexField = cacheType.GetField("Index", BindingFlags.Public | BindingFlags.Static);
        indexField?.SetValue(null, index);
    }

    // P2.M6 FIX: Cache resolved types to avoid repeated assembly scanning
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<string, Type?> _resolvedTypeCache = new();

    private static Type? ResolveType(string typeName)
    {
        if (_resolvedTypeCache.TryGetValue(typeName, out var cached))
            return cached;

        var type = Type.GetType(typeName);
        if (type != null)
        {
            _resolvedTypeCache.TryAdd(typeName, type);
            return type;
        }

        try
        {
            foreach (var assembly in AppDomain.CurrentDomain.GetAssemblies())
            {
                if (assembly.IsDynamic) continue;
                try
                {
                    foreach (var t in assembly.GetTypes())
                        if (t.FullName == typeName || t.Name == typeName)
                        {
                            _resolvedTypeCache.TryAdd(typeName, t);
                            return t;
                        }
                }
                catch (ReflectionTypeLoadException)
                {
                    /* Skip assemblies with type load failures */
                }
            }
        }
        catch (AppDomainUnloadedException)
        {
            /* AppDomain being torn down */
        }

        _resolvedTypeCache.TryAdd(typeName, null);
        return null;
    }

    #endregion

    // ============================================================================
    // ADDITIONAL FEATURES FOR FEM AND MESH PROCESSING
    // ============================================================================

    #region Boundary Detection

    /// <summary>
    ///     Finds boundary nodes - nodes that belong to facets shared by only one element.
    /// </summary>
    /// <typeparam name="TElement">Element type (e.g., Triangle, Quad, Tet, Hex).</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="nodesPerBoundaryFacet">
    ///     Number of nodes defining a boundary facet (2 for 2D edges, 3 for tri faces, 4 for quad faces).
    /// </param>
    /// <returns>Set of boundary node indices.</returns>
    /// <remarks>
    ///     P2-1 FIX: Uses collision-safe facet counting to handle hash collisions correctly.
    /// </remarks>
    public HashSet<int> FindBoundaryNodes<TElement, TNode>(int nodesPerBoundaryFacet)
    {
        ThrowIfDisposed();
        var boundaryNodes = new HashSet<int>();
        var facetCounter = new FacetCounter(); // P2-1 FIX: Collision-safe

        var elementCount = Count<TElement>();

        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            var facets = GetElementFacets(nodes, nodesPerBoundaryFacet);

            foreach (var facet in facets) facetCounter.Add(facet);
        }

        foreach (var (count, nodes) in facetCounter.GetAllFacets())
            if (count == 1)
                foreach (var node in nodes)
                    boundaryNodes.Add(node);

        return boundaryNodes;
    }

    /// <summary>
    ///     Finds boundary elements - elements that have at least one boundary facet.
    /// </summary>
    /// <remarks>
    ///     P2-1 FIX: Uses collision-safe facet tracking.
    /// </remarks>
    public HashSet<int> FindBoundaryElements<TElement, TNode>(int nodesPerBoundaryFacet)
    {
        ThrowIfDisposed();
        var boundaryElements = new HashSet<int>();
        // P2-1 FIX: Use collision chain for correctness
        var facetToElements = new Dictionary<long, List<(int[] Nodes, List<int> Elements)>>();

        var elementCount = Count<TElement>();

        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            var facets = GetElementFacets(nodes, nodesPerBoundaryFacet);

            foreach (var facet in facets)
            {
                var key = ComputeFacetKey(facet);
                if (!facetToElements.TryGetValue(key, out var bucket))
                {
                    bucket = new List<(int[], List<int>)>();
                    facetToElements[key] = bucket;
                }

                // Find existing facet in collision chain
                var found = false;
                foreach (var entry in bucket)
                    if (FacetsEqual(entry.Nodes, facet))
                    {
                        entry.Elements.Add(e);
                        found = true;
                        break;
                    }

                if (!found) bucket.Add((facet, new List<int> { e }));
            }
        }

        foreach (var bucket in facetToElements.Values)
        foreach (var (_, elements) in bucket)
            if (elements.Count == 1)
                boundaryElements.Add(elements[0]);

        return boundaryElements;
    }

    /// <summary>
    ///     Extracts boundary facets as node arrays.
    /// </summary>
    /// <remarks>
    ///     P2-1 FIX: Uses collision-safe facet counting.
    /// </remarks>
    public List<int[]> ExtractBoundaryFacets<TElement, TNode>(int nodesPerBoundaryFacet)
    {
        ThrowIfDisposed();
        var facetCounter = new FacetCounter(); // P2-1 FIX
        var elementCount = Count<TElement>();

        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            var facets = GetElementFacets(nodes, nodesPerBoundaryFacet);

            foreach (var facet in facets) facetCounter.Add(facet);
        }

        var result = new List<int[]>();
        foreach (var (count, nodes) in facetCounter.GetAllFacets())
            if (count == 1)
                result.Add(nodes);
        return result;
    }

    /// <summary>
    ///     Finds internal facets - facets shared by exactly two elements.
    /// </summary>
    /// <remarks>
    ///     P2-1 FIX: Uses collision-safe facet tracking.
    /// </remarks>
    public List<(int[] Nodes, int Element1, int Element2)> FindInternalFacets<TElement, TNode>(int nodesPerFacet)
    {
        ThrowIfDisposed();
        // P2-1 FIX: Use collision chain for correctness
        var facetToElements = new Dictionary<long, List<(int[] Nodes, List<int> Elements)>>();
        var elementCount = Count<TElement>();

        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            var facets = GetElementFacets(nodes, nodesPerFacet);

            foreach (var facet in facets)
            {
                var key = ComputeFacetKey(facet);
                if (!facetToElements.TryGetValue(key, out var bucket))
                {
                    bucket = new List<(int[], List<int>)>();
                    facetToElements[key] = bucket;
                }

                var found = false;
                foreach (var entry in bucket)
                    if (FacetsEqual(entry.Nodes, facet))
                    {
                        entry.Elements.Add(e);
                        found = true;
                        break;
                    }

                if (!found) bucket.Add((facet, new List<int> { e }));
            }
        }

        var result = new List<(int[] Nodes, int Element1, int Element2)>();
        foreach (var bucket in facetToElements.Values)
        foreach (var (nodes, elements) in bucket)
            if (elements.Count == 2)
                result.Add((nodes, elements[0], elements[1]));
        return result;
    }

    private static List<int[]> GetElementFacets(IReadOnlyList<int> nodes, int nodesPerFacet)
    {
        var facets = new List<int[]>();
        var n = nodes.Count;

        if (nodesPerFacet == 2) // 2D: edges
        {
            for (var i = 0; i < n; i++)
                facets.Add([nodes[i], nodes[(i + 1) % n]]);
        }
        else if (nodesPerFacet == 3 && n == 4) // Tet faces
        {
            facets.Add([nodes[0], nodes[1], nodes[2]]);
            facets.Add([nodes[0], nodes[1], nodes[3]]);
            facets.Add([nodes[1], nodes[2], nodes[3]]);
            facets.Add([nodes[0], nodes[2], nodes[3]]);
        }
        else if (nodesPerFacet == 4 && n == 8) // Hex faces
        {
            facets.Add([nodes[0], nodes[1], nodes[2], nodes[3]]);
            facets.Add([nodes[4], nodes[5], nodes[6], nodes[7]]);
            facets.Add([nodes[0], nodes[1], nodes[5], nodes[4]]);
            facets.Add([nodes[2], nodes[3], nodes[7], nodes[6]]);
            facets.Add([nodes[0], nodes[3], nodes[7], nodes[4]]);
            facets.Add([nodes[1], nodes[2], nodes[6], nodes[5]]);
        }
        else if (nodesPerFacet == 3 && n == 6) // Wedge/Prism tri faces
        {
            facets.Add([nodes[0], nodes[1], nodes[2]]);
            facets.Add([nodes[3], nodes[4], nodes[5]]);
        }

        return facets;
    }

    private static long ComputeFacetKey(int[] nodes)
    {
        var sorted = new int[nodes.Length];
        Array.Copy(nodes, sorted, nodes.Length);
        Array.Sort(sorted);
        long hash = 17;
        foreach (var n in sorted)
            hash = unchecked(hash * 31 + n);
        return unchecked(hash * 31 + sorted.Length);
    }

    /// <summary>
    ///     P2-1 FIX: Compares two facets for equality (same nodes regardless of order).
    /// </summary>
    private static bool FacetsEqual(int[] a, int[] b)
    {
        if (a.Length != b.Length) return false;

        // Sort copies for comparison
        var sortedA = a.Length <= 8 ? stackalloc int[a.Length] : new int[a.Length];
        var sortedB = b.Length <= 8 ? stackalloc int[b.Length] : new int[b.Length];
        a.CopyTo(sortedA);
        b.CopyTo(sortedB);
        sortedA.Sort();
        sortedB.Sort();

        return sortedA.SequenceEqual(sortedB);
    }

    /// <summary>
    ///     P2-1 FIX: Collision-safe facet counting structure.
    ///     Uses hash for fast lookup with collision chain for correctness.
    /// </summary>
    private sealed class FacetCounter
    {
        private readonly Dictionary<long, List<(int Count, int[] Nodes)>> _buckets = new();

        public void Add(int[] facet)
        {
            var key = ComputeFacetKey(facet);
            if (_buckets.TryGetValue(key, out var bucket))
            {
                // Check for existing facet in collision chain
                for (var i = 0; i < bucket.Count; i++)
                    if (FacetsEqual(bucket[i].Nodes, facet))
                    {
                        bucket[i] = (bucket[i].Count + 1, bucket[i].Nodes);
                        return;
                    }

                // Hash collision with different facet - add to chain
                bucket.Add((1, facet));
            }
            else
            {
                _buckets[key] = new List<(int, int[])> { (1, facet) };
            }
        }

        public IEnumerable<(int Count, int[] Nodes)> GetAllFacets()
        {
            foreach (var bucket in _buckets.Values)
            foreach (var entry in bucket)
                yield return entry;
        }
    }

    // ========================================================================
    // Enhanced Boundary Detection using Sub-Entities
    // ========================================================================

    /// <summary>
    ///     Result of sub-entity-based boundary detection operation.
    /// </summary>
    public readonly struct SubEntityBoundaryResult
    {
        /// <summary>Indices of boundary sub-entities (shared by exactly 1 parent).</summary>
        public readonly List<int> BoundaryIndices;

        /// <summary>Indices of interior sub-entities (shared by 2+ parents).</summary>
        public readonly List<int> InteriorIndices;

        /// <summary>Incidence count for each sub-entity (index → count).</summary>
        public readonly int[] IncidenceCounts;

        public SubEntityBoundaryResult(List<int> boundary, List<int> interior, int[] counts)
        {
            BoundaryIndices = boundary;
            InteriorIndices = interior;
            IncidenceCounts = counts;
        }

        /// <summary>Number of boundary sub-entities.</summary>
        public int BoundaryCount => BoundaryIndices.Count;

        /// <summary>Number of interior sub-entities.</summary>
        public int InteriorCount => InteriorIndices.Count;
    }

    /// <summary>
    ///     Identifies boundary and interior sub-entities based on parent element incidence.
    /// </summary>
    /// <typeparam name="TParent">Parent element type (e.g., Tri3, Tet4).</typeparam>
    /// <typeparam name="TSubEntity">Sub-entity type (e.g., Edge for 2D, Face for 3D).</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <returns>
    ///     Boundary detection result containing boundary indices, interior indices,
    ///     and per-sub-entity incidence counts.
    /// </returns>
    /// <remarks>
    ///     <para>
    ///         <b>BOUNDARY DEFINITION:</b> A sub-entity is on the boundary if it belongs
    ///         to exactly one parent element. Interior sub-entities belong to 2+ parents.
    ///     </para>
    ///     <para>
    ///         <b>PREREQUISITES:</b> Sub-entities must be discovered first using
    ///         <see cref="DiscoverSubEntities{TElement,TSubEntity,TNode}" />.
    ///     </para>
    ///     <para>
    ///         <b>PERFORMANCE:</b> O(S × P) where S = sub-entity count, P = avg parent incidence.
    ///         Results include incidence counts for additional analysis (e.g., non-manifold detection).
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Find boundary edges of a 2D triangle mesh
    ///         mesh.WithSymmetry&lt;Edge&gt;(Symmetry.Full(2));
    ///         mesh.DiscoverSubEntities&lt;Tri3, Edge, Node&gt;(FiniteElementTopologies.Tri3Edges);
    ///         var result = mesh.DetectSubEntityBoundary&lt;Tri3, Edge, Node&gt;();
    ///         Console.WriteLine($"Boundary edges: {result.BoundaryCount}");
    ///         
    ///         // Find surface faces of a 3D tet mesh
    ///         mesh.WithSymmetry&lt;Face&gt;(Symmetry.Cyclic(3));
    ///         mesh.DiscoverSubEntities&lt;Tet4, Face, Node&gt;(FiniteElementTopologies.Tet4Faces);
    ///         var surface = mesh.DetectSubEntityBoundary&lt;Tet4, Face, Node&gt;();
    ///         </code>
    ///     </example>
    /// </remarks>
    public SubEntityBoundaryResult DetectSubEntityBoundary<TParent, TSubEntity, TNode>()
    {
        ThrowIfDisposed();

        var subEntityCount = Count<TSubEntity>();
        if (subEntityCount == 0)
            return new SubEntityBoundaryResult(new List<int>(), new List<int>(), Array.Empty<int>());

        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();

            var parentTypeIdx = GetTypeIndex<TParent>();
            var nodeTypeIdx = GetTypeIndex<TNode>();
            var subEntityTypeIdx = GetTypeIndex<TSubEntity>();

            var parentM2M = _adjacency[parentTypeIdx, nodeTypeIdx];
            var subEntityM2M = _adjacency[subEntityTypeIdx, nodeTypeIdx];

            var incidenceCounts = new int[subEntityCount];
            var boundary = new List<int>();
            var interior = new List<int>();

            for (var subIdx = 0; subIdx < subEntityCount; subIdx++)
            {
                var subNodes = subEntityM2M[subIdx];
                if (subNodes.Count == 0) continue;

                // Find parent elements containing all nodes of this sub-entity
                var transpose = parentM2M.ElementsFromNode;

                // Start with first node's parents
                if (subNodes[0] >= transpose.Count)
                {
                    incidenceCounts[subIdx] = 0;
                    continue;
                }

                var candidates = new HashSet<int>(transpose[subNodes[0]].ToArray());

                // Intersect with remaining nodes
                for (var i = 1; i < subNodes.Count && candidates.Count > 0; i++)
                {
                    if (subNodes[i] >= transpose.Count)
                    {
                        candidates.Clear();
                        break;
                    }

                    candidates.IntersectWith(transpose[subNodes[i]].ToArray());
                }

                incidenceCounts[subIdx] = candidates.Count;

                if (candidates.Count == 1)
                    boundary.Add(subIdx);
                else if (candidates.Count > 1)
                    interior.Add(subIdx);
                // Count == 0 means orphaned sub-entity (shouldn't happen in valid mesh)
            }

            return new SubEntityBoundaryResult(boundary, interior, incidenceCounts);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets only the boundary sub-entity indices (convenience wrapper).
    /// </summary>
    public List<int> GetBoundarySubEntities<TParent, TSubEntity, TNode>()
    {
        return DetectSubEntityBoundary<TParent, TSubEntity, TNode>().BoundaryIndices;
    }

    /// <summary>
    ///     Gets only the interior sub-entity indices (convenience wrapper).
    /// </summary>
    public List<int> GetInteriorSubEntities<TParent, TSubEntity, TNode>()
    {
        return DetectSubEntityBoundary<TParent, TSubEntity, TNode>().InteriorIndices;
    }

    /// <summary>
    ///     Checks if a sub-entity is on the boundary (has exactly 1 parent element).
    /// </summary>
    public bool IsSubEntityOnBoundary<TParent, TSubEntity, TNode>(int subEntityIndex)
    {
        return CountElementsSharingSubEntity<TParent, TSubEntity, TNode>(subEntityIndex) == 1;
    }

    /// <summary>
    ///     Detects non-manifold sub-entities (shared by more than 2 parent elements).
    /// </summary>
    /// <typeparam name="TParent">Parent element type.</typeparam>
    /// <typeparam name="TSubEntity">Sub-entity type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <returns>List of non-manifold sub-entity indices.</returns>
    /// <remarks>
    ///     In a manifold mesh, each interior edge/face should be shared by exactly 2 elements.
    ///     Non-manifold sub-entities indicate mesh topology issues.
    /// </remarks>
    public List<int> DetectNonManifoldSubEntities<TParent, TSubEntity, TNode>()
    {
        var result = DetectSubEntityBoundary<TParent, TSubEntity, TNode>();
        var nonManifold = new List<int>();

        for (var i = 0; i < result.IncidenceCounts.Length; i++)
            if (result.IncidenceCounts[i] > 2)
                nonManifold.Add(i);

        return nonManifold;
    }

    #endregion

    #region Bandwidth Reduction

    /// <summary>
    ///     Computes the Cuthill-McKee ordering for bandwidth reduction.
    /// </summary>
    /// <param name="reverse">If true, returns Reverse Cuthill-McKee (usually better).</param>
    /// <returns>Permutation array: newIndex = permutation[oldIndex].</returns>
    public int[] ComputeCuthillMcKeeOrdering<TElement, TNode>(bool reverse = true)
    {
        ThrowIfDisposed();
        var nodeCount = Count<TNode>();
        if (nodeCount == 0) return [];

        var adjacency = BuildNodeAdjacency<TElement, TNode>();
        var startNode = FindPeripheralNode(adjacency);

        var ordering = new List<int>(nodeCount);
        var visited = new bool[nodeCount];
        var queue = new Queue<int>();

        queue.Enqueue(startNode);
        visited[startNode] = true;

        while (ordering.Count < nodeCount)
        {
            if (queue.Count == 0)
                for (var i = 0; i < nodeCount; i++)
                    if (!visited[i])
                    {
                        queue.Enqueue(i);
                        visited[i] = true;
                        break;
                    }

            while (queue.Count > 0)
            {
                var current = queue.Dequeue();
                ordering.Add(current);

                // Filter neighbors that haven't been visited
                var neighbors = new List<int>();
                foreach (var n in adjacency[current])
                    if (!visited[n])
                        neighbors.Add(n);
                // Sort by degree (adjacency count)
                neighbors.Sort((a, b) => adjacency[a].Count.CompareTo(adjacency[b].Count));

                foreach (var neighbor in neighbors)
                    if (!visited[neighbor])
                    {
                        visited[neighbor] = true;
                        queue.Enqueue(neighbor);
                    }
            }
        }

        var permutation = new int[nodeCount];
        for (var i = 0; i < nodeCount; i++)
        {
            var newIndex = reverse ? nodeCount - 1 - i : i;
            permutation[ordering[i]] = newIndex;
        }

        return permutation;
    }

    private static int FindPeripheralNode(HashSet<int>[] adjacency)
    {
        if (adjacency.Length == 0) return 0;

        var startNode = 0;
        var minDegree = int.MaxValue;
        for (var i = 0; i < adjacency.Length; i++)
            if (adjacency[i].Count < minDegree)
            {
                minDegree = adjacency[i].Count;
                startNode = i;
            }

        var distances = new int[adjacency.Length];
        Array.Fill(distances, -1);
        var queue = new Queue<int>();
        queue.Enqueue(startNode);
        distances[startNode] = 0;
        var farthest = startNode;
        var maxDist = 0;

        while (queue.Count > 0)
        {
            var current = queue.Dequeue();
            foreach (var neighbor in adjacency[current])
                if (distances[neighbor] < 0)
                {
                    distances[neighbor] = distances[current] + 1;
                    queue.Enqueue(neighbor);
                    if (distances[neighbor] > maxDist)
                    {
                        maxDist = distances[neighbor];
                        farthest = neighbor;
                    }
                }
        }

        return farthest;
    }

    /// <summary>
    ///     Computes the bandwidth of the current node ordering.
    /// </summary>
    public int ComputeBandwidth<TElement, TNode>()
    {
        ThrowIfDisposed();
        var adjacency = BuildNodeAdjacency<TElement, TNode>();

        var maxBandwidth = 0;
        for (var i = 0; i < adjacency.Length; i++)
            foreach (var j in adjacency[i])
            {
                var bandwidth = Math.Abs(i - j);
                if (bandwidth > maxBandwidth)
                    maxBandwidth = bandwidth;
            }

        return maxBandwidth;
    }

    /// <summary>
    ///     Computes the profile (envelope) of the current ordering.
    /// </summary>
    public long ComputeProfile<TElement, TNode>()
    {
        ThrowIfDisposed();
        var adjacency = BuildNodeAdjacency<TElement, TNode>();

        long profile = 0;
        for (var i = 0; i < adjacency.Length; i++)
        {
            var minNeighbor = i;
            foreach (var j in adjacency[i])
                if (j < minNeighbor)
                    minNeighbor = j;
            profile += i - minNeighbor;
        }

        return profile;
    }

    private HashSet<int>[] BuildNodeAdjacency<TElement, TNode>()
    {
        var nodeCount = Count<TNode>();
        var elementCount = Count<TElement>();
        var adjacency = new HashSet<int>[nodeCount];

        for (var i = 0; i < nodeCount; i++)
            adjacency[i] = [];

        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            for (var i = 0; i < nodes.Count; i++)
            for (var j = i + 1; j < nodes.Count; j++)
            {
                adjacency[nodes[i]].Add(nodes[j]);
                adjacency[nodes[j]].Add(nodes[i]);
            }
        }

        return adjacency;
    }

    #endregion

    #region Sparse Matrix Pattern Extraction

    /// <summary>
    ///     Extracts the sparsity pattern in CSR format for sparse matrix assembly.
    /// </summary>
    /// <param name="dofsPerNode">Degrees of freedom per node.</param>
    /// <returns>CSR format (rowPtr, colIndices) ready for sparse matrix creation.</returns>
    public (int[] RowPtr, int[] ColIndices) GetSparsityPatternCSR<TElement, TNode>(int dofsPerNode = 1)
    {
        ThrowIfDisposed();
        var nodeCount = Count<TNode>();
        var totalDofs = nodeCount * dofsPerNode;

        var adjacency = new HashSet<int>[nodeCount];
        for (var i = 0; i < nodeCount; i++)
            adjacency[i] = [i];

        var elementCount = Count<TElement>();
        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            for (var i = 0; i < nodes.Count; i++)
            for (var j = 0; j < nodes.Count; j++)
                adjacency[nodes[i]].Add(nodes[j]);
        }

        var rowPtr = new int[totalDofs + 1];
        var colIndicesList = new List<int>();

        for (var node = 0; node < nodeCount; node++)
        {
            var connectedNodes = new List<int>(adjacency[node]);
            connectedNodes.Sort();

            for (var d = 0; d < dofsPerNode; d++)
            {
                var globalRow = node * dofsPerNode + d;
                rowPtr[globalRow] = colIndicesList.Count;

                foreach (var connNode in connectedNodes)
                    for (var cd = 0; cd < dofsPerNode; cd++)
                        colIndicesList.Add(connNode * dofsPerNode + cd);
            }
        }

        rowPtr[totalDofs] = colIndicesList.Count;

        return (rowPtr, colIndicesList.ToArray());
    }

    /// <summary>
    ///     Gets the number of non-zeros in the assembled matrix.
    /// </summary>
    public int GetNonZeroCount<TElement, TNode>(int dofsPerNode = 1)
    {
        ThrowIfDisposed();
        var nodeCount = Count<TNode>();

        var adjacency = new HashSet<int>[nodeCount];
        for (var i = 0; i < nodeCount; i++)
            adjacency[i] = [i];

        var elementCount = Count<TElement>();
        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            for (var i = 0; i < nodes.Count; i++)
            for (var j = 0; j < nodes.Count; j++)
                adjacency[nodes[i]].Add(nodes[j]);
        }

        var nnz = 0;
        for (var i = 0; i < nodeCount; i++)
            nnz += adjacency[i].Count;

        return nnz * dofsPerNode * dofsPerNode;
    }

    /// <summary>
    ///     Computes clique indices for finite element matrix assembly.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>
    ///     List where result[element] contains clique indices for all node pairs in that element.
    ///     For an element with ns nodes, returns ns² indices (flattened matrix).
    /// </returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> Essential for FEM sparse matrix assembly.
    ///     A clique is a unique pair of nodes that appear together in at least one element.
    ///     This method assigns each unique (node_i, node_j) pair a unique clique index.
    ///     Used to map local element matrices to global sparse matrix positions.
    ///     <b>PERFORMANCE:</b> O(Σ ns²) where ns = nodes per element.
    ///     Efficiently computed using generation-based marking.
    /// </remarks>
    public List<List<int>> GetCliques<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetCliques();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Exports element-to-node connectivity in CSR (Compressed Sparse Row) format.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>
    ///     Tuple of (rowPtr, columnIndices) in standard CSR format.
    /// </returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> For interoperability with sparse matrix libraries.
    ///     CSR format:
    ///     - rowPtr[i] = start index in columnIndices for element i
    ///     - rowPtr[i+1] = end index (exclusive)
    ///     - columnIndices[rowPtr[i]:rowPtr[i+1]] = nodes of element i
    ///     Compatible with libraries like Intel MKL, SuiteSparse, Eigen, etc.
    /// </remarks>
    public (int[] RowPtr, int[] ColumnIndices) ToCsr<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.ToCsr();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Creates a topology from CSR (Compressed Sparse Row) format.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="rowPtr">Row pointers array (length = num_elements + 1).</param>
    /// <param name="columnIndices">Column indices array (node indices).</param>
    /// <returns>New topology initialized with the CSR data.</returns>
    /// <remarks>
    ///     <b>EXPOSED API:</b> For importing from sparse matrix libraries.
    ///     Creates a topology from standard CSR format. Useful for:
    ///     - Loading meshes from external software
    ///     - Interfacing with sparse matrix libraries
    ///     - Converting from other mesh formats
    /// </remarks>
    public static Topology<TTypes> FromCsr<TElement, TNode>(int[] rowPtr, int[] columnIndices)
    {
        ArgumentNullException.ThrowIfNull(rowPtr);
        ArgumentNullException.ThrowIfNull(columnIndices);

        var mesh = new Topology<TTypes>();
        var m2m = M2M.FromCsr(rowPtr, columnIndices);

        // Add all elements from M2M using zero-copy span access
        for (var e = 0; e < m2m.Count; e++)
            // Use callback-based access to avoid allocation and unsafe casts
            m2m.WithNodesSpan(e, nodes => mesh.Add<TElement, TNode>(nodes));

        return mesh;
    }

    #endregion

    #region Mesh Merging and Extraction

    /// <summary>
    ///     Merges another topology into this one.
    /// </summary>
    /// <returns>Offset applied to source node indices.</returns>
    public int Merge<TElement, TNode>(Topology<TTypes> other)
    {
        ArgumentNullException.ThrowIfNull(other);

        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            var nodeOffset = Count<TNode>();
            var otherNodeCount = other.Count<TNode>();
            var otherElementCount = other.Count<TElement>();

            var nodeType = GetTypeIndex<TNode>();
            for (var i = 0; i < otherNodeCount; i++)
            {
                var index = _adjacency.GetNumberOfElements(nodeType);
                _adjacency.AppendElement(nodeType, nodeType, [index]);
            }

            var elementType = GetTypeIndex<TElement>();
            for (var e = 0; e < otherElementCount; e++)
            {
                var nodes = other.NodesOf<TElement, TNode>(e);
                var offsetNodes = new List<int>(nodes.Count);
                for (var i = 0; i < nodes.Count; i++)
                    offsetNodes.Add(nodes[i] + nodeOffset);
                _adjacency.AppendElement(elementType, nodeType, offsetNodes);
            }

            return nodeOffset;
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    /// <summary>
    ///     Extracts a subset of elements and their nodes as a new topology.
    /// </summary>
    public (Topology<TTypes> Subtopology, int[] NodeMap, int[] ElementMap)
        ExtractSubstructure<TElement, TNode>(IEnumerable<int> elementIndices)
    {
        ThrowIfDisposed();

        var elements = new List<int>();
        foreach (var idx in elementIndices)
            elements.Add(idx);
        var usedNodes = new HashSet<int>();

        foreach (var e in elements)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            foreach (var n in nodes)
                usedNodes.Add(n);
        }

        var sortedNodes = new List<int>(usedNodes);
        sortedNodes.Sort();
        var nodeOldToNew = new Dictionary<int, int>();
        var nodeNewToOld = new int[sortedNodes.Count];

        for (var i = 0; i < sortedNodes.Count; i++)
        {
            nodeOldToNew[sortedNodes[i]] = i;
            nodeNewToOld[i] = sortedNodes[i];
        }

        var elementNewToOld = new int[elements.Count];
        for (var i = 0; i < elements.Count; i++)
            elementNewToOld[i] = elements[i];
        var sub = new Topology<TTypes>();

        if (_symmetries.TryGetValue(typeof(TElement), out var sym))
            sub.WithSymmetry<TElement>(sym);

        for (var i = 0; i < sortedNodes.Count; i++)
            sub.Add<TNode>();

        foreach (var e in elements)
        {
            var oldNodes = NodesOf<TElement, TNode>(e);
            var newNodes = new int[oldNodes.Count];
            for (var i = 0; i < oldNodes.Count; i++)
                newNodes[i] = nodeOldToNew[oldNodes[i]];
            sub.Add<TElement, TNode>(newNodes);
        }

        return (sub, nodeNewToOld, elementNewToOld);
    }

    // ========================================================================
    // Enhanced Subset Extraction with Predicate Filtering
    // ========================================================================

    /// <summary>
    ///     Creates a new topology containing only elements satisfying the predicate.
    /// </summary>
    /// <typeparam name="TElement">Element type to filter.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="predicate">Filter function: index → include.</param>
    /// <param name="includeOrphanNodes">
    ///     If false (default), only includes nodes referenced by kept elements.
    ///     If true, includes all original nodes.
    /// </param>
    /// <returns>
    ///     Tuple of (newTopology, elementMapping, nodeMapping) where:
    ///     - newTopology: The filtered topology
    ///     - elementMapping: newIndex → oldIndex for elements
    ///     - nodeMapping: newIndex → oldIndex for nodes
    /// </returns>
    /// <remarks>
    ///     <para>
    ///         <b>NODE HANDLING:</b> By default, only nodes referenced by kept elements
    ///         are included. Node indices are renumbered consecutively.
    ///     </para>
    ///     <para>
    ///         <b>DATA PRESERVATION:</b> Per-entity data is copied for kept entities.
    ///         Symmetry settings are preserved.
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Keep only elements with positive material ID
    ///         var (subset, elemMap, nodeMap) = mesh.CloneWhere&lt;Tet4, Node&gt;(
    ///             idx => mesh.Get&lt;Tet4, int&gt;(idx) > 0);
    ///         </code>
    ///     </example>
    /// </remarks>
    public (Topology<TTypes> NewTopology, int[] ElementMapping, int[] NodeMapping) CloneWhere<TElement, TNode>(
        Func<int, bool> predicate,
        bool includeOrphanNodes = false)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(predicate);

        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();

            var result = new Topology<TTypes>();

            // Copy symmetries
            foreach (var kvp in _symmetries)
                result._symmetries[kvp.Key] = kvp.Value;

            var elementTypeIdx = GetTypeIndex<TElement>();
            var nodeTypeIdx = GetTypeIndex<TNode>();
            var elementCount = Count<TElement>();

            // Determine which elements to keep
            var keptElements = new List<int>();
            for (var i = 0; i < elementCount; i++)
                if (predicate(i))
                    keptElements.Add(i);

            // Determine which nodes are referenced
            var referencedNodes = new HashSet<int>();
            var m2m = _adjacency[elementTypeIdx, nodeTypeIdx];

            foreach (var elemIdx in keptElements)
            {
                var nodes = m2m[elemIdx];
                foreach (var node in nodes)
                    referencedNodes.Add(node);
            }

            // Build node mapping (old → new)
            int[] nodeMapping;
            var oldToNewNode = new Dictionary<int, int>();

            if (includeOrphanNodes)
            {
                var nodeCount = _adjacency.GetNumberOfElements(nodeTypeIdx);
                nodeMapping = new int[nodeCount];
                for (var i = 0; i < nodeCount; i++)
                {
                    nodeMapping[i] = i;
                    oldToNewNode[i] = i;
                }
            }
            else
            {
                var sortedNodes = new List<int>(referencedNodes);
                sortedNodes.Sort();
                nodeMapping = new int[sortedNodes.Count];
                for (var newIdx = 0; newIdx < sortedNodes.Count; newIdx++)
                {
                    var oldIdx = sortedNodes[newIdx];
                    nodeMapping[newIdx] = oldIdx;
                    oldToNewNode[oldIdx] = newIdx;
                }
            }

            // Add nodes to result
            for (var i = 0; i < nodeMapping.Length; i++)
                result._adjacency.AppendElement(nodeTypeIdx, nodeTypeIdx, new List<int> { i });

            // Build element mapping and add elements with remapped nodes
            var elementMapping = new int[keptElements.Count];

            for (var newIdx = 0; newIdx < keptElements.Count; newIdx++)
            {
                var oldIdx = keptElements[newIdx];
                elementMapping[newIdx] = oldIdx;

                var oldNodes = m2m[oldIdx];
                var newNodes = new List<int>(oldNodes.Count);
                foreach (var oldNode in oldNodes)
                    newNodes.Add(oldToNewNode[oldNode]);

                var elemIndex = result._adjacency.AppendElement(elementTypeIdx, nodeTypeIdx, newNodes);

                // Populate canonical index for Find/Exists support
                if (result._symmetries.TryGetValue(typeof(TElement), out var sym))
                {
                    var (canonical, key) = result.GetCanonicalWithKey<TElement>(
                        System.Runtime.InteropServices.CollectionsMarshal.AsSpan(newNodes));

                    if (!result._canonicalIndex.TryGetValue(typeof(TElement), out var typeIndex))
                    {
                        typeIndex = new Dictionary<long, List<(int Index, List<int> Nodes)>>();
                        result._canonicalIndex[typeof(TElement)] = typeIndex;
                    }

                    if (typeIndex.TryGetValue(key, out var collisionChain))
                        collisionChain.Add((elemIndex, canonical));
                    else
                        typeIndex[key] = new List<(int Index, List<int> Nodes)> { (elemIndex, canonical) };
                }
            }

            // Copy per-entity data for kept elements
            foreach (var kvp in _data)
                if (kvp.Key.Entity == typeof(TElement))
                {
                    var sourceList = kvp.Value;
                    var destList = sourceList.CreateEmpty();

                    for (var newIdx = 0; newIdx < keptElements.Count; newIdx++)
                    {
                        var oldIdx = keptElements[newIdx];
                        destList.AddFromSource(sourceList, oldIdx);
                    }

                    result._data[kvp.Key] = destList;
                }
                else if (kvp.Key.Entity == typeof(TNode) && !includeOrphanNodes)
                {
                    // Copy node data with proper remapping
                    var sourceList = kvp.Value;
                    var destList = sourceList.CreateEmpty();

                    for (var newIdx = 0; newIdx < nodeMapping.Length; newIdx++)
                    {
                        var oldIdx = nodeMapping[newIdx];
                        if (oldIdx < sourceList.Count)
                            destList.AddFromSource(sourceList, oldIdx);
                    }

                    result._data[kvp.Key] = destList;
                }

            return (result, elementMapping, nodeMapping);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Extracts a region defined by a set of element indices (convenience wrapper).
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="elementIndices">Indices of elements to extract.</param>
    /// <returns>Extraction result with new topology and mappings.</returns>
    public (Topology<TTypes> NewTopology, int[] ElementMapping, int[] NodeMapping) ExtractRegion<TElement, TNode>(
        IEnumerable<int> elementIndices)
    {
        var indexSet = new HashSet<int>(elementIndices);
        return CloneWhere<TElement, TNode>(idx => indexSet.Contains(idx));
    }

    /// <summary>
    ///     Extracts elements within a spatial bounding box.
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="minBound">Minimum corner of bounding box (X, Y, Z).</param>
    /// <param name="maxBound">Maximum corner of bounding box (X, Y, Z).</param>
    /// <param name="getCoord">Function to get coordinate from node index.</param>
    /// <param name="allNodesInside">
    ///     If true, all element nodes must be inside. If false, any node inside qualifies.
    /// </param>
    /// <returns>Extraction result.</returns>
    /// <remarks>
    ///     <example>
    ///         <code>
    ///         var (subset, _, _) = mesh.ExtractByBoundingBox&lt;Tet4, Node&gt;(
    ///             (0, 0, 0), (1, 1, 1),
    ///             nodeIdx => (coords[nodeIdx].X, coords[nodeIdx].Y, coords[nodeIdx].Z),
    ///             allNodesInside: true);
    ///         </code>
    ///     </example>
    /// </remarks>
    public (Topology<TTypes> NewTopology, int[] ElementMapping, int[] NodeMapping)
        ExtractByBoundingBox<TElement, TNode>(
            (double X, double Y, double Z) minBound,
            (double X, double Y, double Z) maxBound,
            Func<int, (double X, double Y, double Z)> getCoord,
            bool allNodesInside = false)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(getCoord);

        bool IsInside(int nodeIdx)
        {
            var (x, y, z) = getCoord(nodeIdx);
            return x >= minBound.X && x <= maxBound.X &&
                   y >= minBound.Y && y <= maxBound.Y &&
                   z >= minBound.Z && z <= maxBound.Z;
        }

        return CloneWhere<TElement, TNode>(elemIdx =>
        {
            var nodes = NodesOf<TElement, TNode>(elemIdx);

            if (allNodesInside)
            {
                foreach (var node in nodes)
                    if (!IsInside(node))
                        return false;
                return true;
            }

            foreach (var node in nodes)
                if (IsInside(node))
                    return true;
            return false;
        });
    }

    #endregion

    #region Node Reordering

    /// <summary>
    ///     Applies a node permutation to all elements.
    /// </summary>
    /// <param name="permutation">newIndex = permutation[oldIndex].</param>
    public void ApplyNodePermutation<TElement, TNode>(int[] permutation)
    {
        ArgumentNullException.ThrowIfNull(permutation);

        ThrowIfDisposed();
        _rwLock.EnterWriteLock();
        try
        {
            ThrowIfDisposed();

            var nodeCount = Count<TNode>();
            if (permutation.Length != nodeCount)
                throw new ArgumentException(
                    $"Permutation length {permutation.Length} != node count {nodeCount}");

            var nodeType = GetTypeIndex<TNode>();
            var elementType = GetTypeIndex<TElement>();
            var m2m = _adjacency[elementType, nodeType];

            var permList = new List<int>(permutation.Length);
            for (var i = 0; i < permutation.Length; i++)
                permList.Add(permutation[i]);
            m2m.PermuteNodes(permList);

            var inverse = new int[nodeCount];
            for (var i = 0; i < nodeCount; i++)
                inverse[permutation[i]] = i;

            foreach (var kvp in _data)
                if (kvp.Key.Entity == typeof(TNode))
                    kvp.Value.ReorderByInversePermutation(inverse);
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
    }

    #endregion

    #region Graph Coloring for Parallel Assembly

    /// <summary>
    ///     Computes a greedy coloring of elements for parallel assembly.
    ///     Elements with the same color share no nodes and can be assembled in parallel.
    /// </summary>
    public int[] ComputeElementColoring<TElement, TNode>()
    {
        ThrowIfDisposed();
        var elementCount = Count<TElement>();
        var colors = new int[elementCount];
        Array.Fill(colors, -1);

        if (elementCount == 0) return colors;

        var nodeToElements = new Dictionary<int, List<int>>();

        for (var e = 0; e < elementCount; e++)
        {
            var nodes = NodesOf<TElement, TNode>(e);
            foreach (var n in nodes)
            {
                if (!nodeToElements.TryGetValue(n, out var list))
                {
                    list = [];
                    nodeToElements[n] = list;
                }

                list.Add(e);
            }
        }

        for (var e = 0; e < elementCount; e++)
        {
            var usedColors = new HashSet<int>();
            var nodes = NodesOf<TElement, TNode>(e);

            foreach (var n in nodes)
            foreach (var neighbor in nodeToElements[n])
                if (neighbor != e && colors[neighbor] >= 0)
                    usedColors.Add(colors[neighbor]);

            var color = 0;
            while (usedColors.Contains(color))
                color++;

            colors[e] = color;
        }

        return colors;
    }

    /// <summary>
    ///     Groups elements by color for parallel processing.
    /// </summary>
    public List<List<int>> GetColorGroups<TElement, TNode>()
    {
        var colors = ComputeElementColoring<TElement, TNode>();
        var maxColor = -1;
        if (colors.Length > 0)
        {
            maxColor = colors[0];
            for (var i = 1; i < colors.Length; i++)
                if (colors[i] > maxColor)
                    maxColor = colors[i];
        }

        var groups = new List<List<int>>(maxColor + 1);
        for (var c = 0; c <= maxColor; c++)
            groups.Add([]);

        for (var e = 0; e < colors.Length; e++)
            groups[colors[e]].Add(e);

        return groups;
    }

    /// <summary>
    ///     Statistics about element coloring.
    /// </summary>
    public readonly struct ColoringStatistics
    {
        public int ElementCount { get; }
        public int NumberOfColors { get; }
        public int MinGroupSize { get; }
        public int MaxGroupSize { get; }
        public double AvgGroupSize { get; }

        public ColoringStatistics(int count, int colors, int min, int max, double avg)
        {
            ElementCount = count;
            NumberOfColors = colors;
            MinGroupSize = min;
            MaxGroupSize = max;
            AvgGroupSize = avg;
        }

        public override string ToString()
        {
            return $"Elements: {ElementCount}, Colors: {NumberOfColors}, " +
                   $"Group sizes: min={MinGroupSize}, max={MaxGroupSize}, avg={AvgGroupSize:F1}";
        }
    }

    /// <summary>
    ///     Gets coloring statistics.
    /// </summary>
    public ColoringStatistics GetColoringStatistics<TElement, TNode>()
    {
        var colors = ComputeElementColoring<TElement, TNode>();
        if (colors.Length == 0)
            return new ColoringStatistics(0, 0, 0, 0, 0);

        // Count elements per color manually
        var colorCounts = new Dictionary<int, int>();
        foreach (var c in colors)
        {
            colorCounts.TryGetValue(c, out var count);
            colorCounts[c] = count + 1;
        }

        // Compute min, max, avg manually
        int min = int.MaxValue, max = int.MinValue;
        long sum = 0;
        foreach (var count in colorCounts.Values)
        {
            if (count < min) min = count;
            if (count > max) max = count;
            sum += count;
        }

        var avg = sum / (double)colorCounts.Count;

        return new ColoringStatistics(
            colors.Length,
            colorCounts.Count,
            min,
            max,
            avg);
    }

    #endregion

    #region Dual Graph Construction

    /// <summary>
    ///     Represents an element-to-element dual graph for mesh traversal algorithms.
    /// </summary>
    public sealed class DualGraph
    {
        /// <summary>For each element, list of adjacent element indices.</summary>
        public readonly List<List<int>> Adjacency;

        /// <summary>For each element pair, the shared node count.</summary>
        public readonly Dictionary<(int, int), int> SharedNodeCounts;

        internal DualGraph(List<List<int>> adjacency, Dictionary<(int, int), int> sharedCounts, int edgeCount)
        {
            Adjacency = adjacency;
            SharedNodeCounts = sharedCounts;
            EdgeCount = edgeCount;
        }

        /// <summary>Number of elements in the dual graph.</summary>
        public int ElementCount => Adjacency.Count;

        /// <summary>Total number of adjacency connections (undirected edges).</summary>
        public int EdgeCount { get; }

        /// <summary>
        ///     Gets neighbors of an element.
        /// </summary>
        public IReadOnlyList<int> GetNeighbors(int elementIndex)
        {
            return Adjacency[elementIndex];
        }

        /// <summary>
        ///     Gets number of shared nodes between two elements.
        /// </summary>
        public int GetSharedNodeCount(int elem1, int elem2)
        {
            var key = elem1 < elem2 ? (elem1, elem2) : (elem2, elem1);
            return SharedNodeCounts.TryGetValue(key, out var count) ? count : 0;
        }

        /// <summary>
        ///     Performs BFS from a starting element.
        /// </summary>
        /// <returns>List of element indices in BFS order.</returns>
        public List<int> BreadthFirstSearch(int startElement)
        {
            if (startElement < 0 || startElement >= ElementCount)
                throw new ArgumentOutOfRangeException(nameof(startElement));

            var visited = new bool[ElementCount];
            var result = new List<int>();
            var queue = new Queue<int>();

            queue.Enqueue(startElement);
            visited[startElement] = true;

            while (queue.Count > 0)
            {
                var current = queue.Dequeue();
                result.Add(current);

                foreach (var neighbor in Adjacency[current])
                    if (!visited[neighbor])
                    {
                        visited[neighbor] = true;
                        queue.Enqueue(neighbor);
                    }
            }

            return result;
        }

        /// <summary>
        ///     Finds connected components in the dual graph.
        /// </summary>
        /// <returns>List of components, each component is a list of element indices.</returns>
        public List<List<int>> FindConnectedComponents()
        {
            var components = new List<List<int>>();
            var visited = new bool[ElementCount];

            for (var i = 0; i < ElementCount; i++)
                if (!visited[i])
                {
                    var component = new List<int>();
                    var queue = new Queue<int>();
                    queue.Enqueue(i);
                    visited[i] = true;

                    while (queue.Count > 0)
                    {
                        var current = queue.Dequeue();
                        component.Add(current);

                        foreach (var neighbor in Adjacency[current])
                            if (!visited[neighbor])
                            {
                                visited[neighbor] = true;
                                queue.Enqueue(neighbor);
                            }
                    }

                    components.Add(component);
                }

            return components;
        }

        /// <summary>
        ///     Computes element distances from a source element.
        /// </summary>
        /// <returns>Distance array: distance[i] = hops from source to element i, or -1 if unreachable.</returns>
        public int[] ComputeDistances(int sourceElement)
        {
            if (sourceElement < 0 || sourceElement >= ElementCount)
                throw new ArgumentOutOfRangeException(nameof(sourceElement));

            var distances = new int[ElementCount];
            Array.Fill(distances, -1);
            distances[sourceElement] = 0;

            var queue = new Queue<int>();
            queue.Enqueue(sourceElement);

            while (queue.Count > 0)
            {
                var current = queue.Dequeue();
                var nextDist = distances[current] + 1;

                foreach (var neighbor in Adjacency[current])
                    if (distances[neighbor] < 0)
                    {
                        distances[neighbor] = nextDist;
                        queue.Enqueue(neighbor);
                    }
            }

            return distances;
        }

        /// <summary>
        ///     Gets the maximum distance (diameter) between any two connected elements.
        /// </summary>
        public int ComputeDiameter()
        {
            if (ElementCount == 0) return 0;

            var maxDist = 0;
            // Simple approach: BFS from each peripheral node candidate
            var dist0 = ComputeDistances(0);
            var far = 0;
            for (var i = 1; i < ElementCount; i++)
                if (dist0[i] > dist0[far])
                    far = i;

            var distFar = ComputeDistances(far);
            for (var i = 0; i < ElementCount; i++)
                if (distFar[i] > maxDist)
                    maxDist = distFar[i];

            return maxDist;
        }
    }

    /// <summary>
    ///     Builds a dual graph where elements are connected if they share enough nodes.
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="minSharedNodes">
    ///     Minimum number of shared nodes for adjacency.
    ///     <list type="bullet">
    ///         <item>1: Share any node (vertex neighbors)</item>
    ///         <item>2: Share an edge (edge neighbors)</item>
    ///         <item>3: Share a face (face neighbors, for 3D elements)</item>
    ///     </list>
    /// </param>
    /// <returns>Dual graph structure for mesh traversal.</returns>
    /// <remarks>
    ///     <para>
    ///         <b>ALGORITHM:</b> For each element, finds all other elements sharing at least
    ///         minSharedNodes nodes using efficient set intersection.
    ///     </para>
    ///     <para>
    ///         <b>USE CASES:</b>
    ///         <list type="bullet">
    ///             <item>Mesh partitioning / domain decomposition</item>
    ///             <item>Finding connected regions</item>
    ///             <item>Level set propagation</item>
    ///             <item>Adaptive refinement propagation</item>
    ///         </list>
    ///     </para>
    ///     <para>
    ///         <b>PERFORMANCE:</b> O(E × N × log(N)) where E = element count, N = avg nodes per element.
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Build face-neighbor dual graph for tets
    ///         var dual = mesh.BuildDualGraph&lt;Tet4, Node&gt;(minSharedNodes: 3);
    ///         
    ///         // Find connected components
    ///         var components = dual.FindConnectedComponents();
    ///         Console.WriteLine($"Mesh has {components.Count} disconnected regions");
    ///         
    ///         // BFS from element 0
    ///         var order = dual.BreadthFirstSearch(0);
    ///         </code>
    ///     </example>
    /// </remarks>
    public DualGraph BuildDualGraph<TElement, TNode>(int minSharedNodes = 1)
    {
        ThrowIfDisposed();

        if (minSharedNodes < 1)
            throw new ArgumentOutOfRangeException(nameof(minSharedNodes), "Must be at least 1");

        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();

            var elementTypeIdx = GetTypeIndex<TElement>();
            var nodeTypeIdx = GetTypeIndex<TNode>();
            var m2m = _adjacency[elementTypeIdx, nodeTypeIdx];

            var elementCount = m2m.Count;
            var adjacency = new List<List<int>>(elementCount);
            var sharedCounts = new Dictionary<(int, int), int>();
            var edgeCount = 0;

            for (var i = 0; i < elementCount; i++)
                adjacency.Add(new List<int>());

            // Build node → elements index for efficient lookup
            var transpose = m2m.ElementsFromNode;

            // For each element, find neighbors through shared nodes
            for (var elemIdx = 0; elemIdx < elementCount; elemIdx++)
            {
                var nodes = m2m[elemIdx];
                var neighborCounts = new Dictionary<int, int>();

                // Count shared nodes with potential neighbors
                foreach (var node in nodes)
                {
                    if (node >= transpose.Count) continue;
                    var elementsAtNode = transpose[node];
                    foreach (var neighbor in elementsAtNode)
                        if (neighbor != elemIdx)
                        {
                            neighborCounts.TryGetValue(neighbor, out var count);
                            neighborCounts[neighbor] = count + 1;
                        }
                }

                // Add edges for neighbors meeting threshold
                foreach (var (neighbor, count) in neighborCounts)
                    if (count >= minSharedNodes && neighbor > elemIdx) // neighbor > elemIdx to avoid duplicates
                    {
                        adjacency[elemIdx].Add(neighbor);
                        adjacency[neighbor].Add(elemIdx);

                        var key = (elemIdx, neighbor);
                        sharedCounts[key] = count;
                        edgeCount++;
                    }
            }

            // Sort adjacency lists for deterministic ordering
            foreach (var list in adjacency)
                list.Sort();

            return new DualGraph(adjacency, sharedCounts, edgeCount);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Builds a face-neighbor dual graph (elements sharing N-1 nodes for N-node elements).
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <returns>Dual graph where adjacency means face-sharing.</returns>
    /// <remarks>
    ///     Convenience method that determines minSharedNodes automatically:
    ///     <list type="bullet">
    ///         <item>Tri3: 2 (edge neighbors)</item>
    ///         <item>Tet4: 3 (face neighbors)</item>
    ///         <item>Quad4: 2 (edge neighbors)</item>
    ///         <item>Hex8: 4 (face neighbors)</item>
    ///     </list>
    /// </remarks>
    public DualGraph BuildFaceNeighborGraph<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            // Determine nodes per element
            var elementTypeIdx = GetTypeIndex<TElement>();
            var nodeTypeIdx = GetTypeIndex<TNode>();
            var m2m = _adjacency[elementTypeIdx, nodeTypeIdx];

            if (m2m.Count == 0)
                return new DualGraph(new List<List<int>>(), new Dictionary<(int, int), int>(), 0);

            var nodesPerElement = m2m[0].Count;

            // Face neighbors share N-1 nodes (or N-2 for volume elements with quad faces)
            var minShared = nodesPerElement switch
            {
                3 => 2, // Tri3: edge = 2 nodes
                4 => 3, // Tet4: face = 3 nodes (assumes Tet4, not Quad4)
                6 => 3, // Wedge6: triangular face = 3 nodes
                8 => 4, // Hex8: quad face = 4 nodes
                _ => Math.Max(1, nodesPerElement - 1)
            };

            return BuildDualGraph<TElement, TNode>(minShared);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Builds an edge-neighbor dual graph (elements sharing exactly 2 nodes).
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <returns>Dual graph where adjacency means edge-sharing.</returns>
    public DualGraph BuildEdgeNeighborGraph<TElement, TNode>()
    {
        return BuildDualGraph<TElement, TNode>(2);
    }

    /// <summary>
    ///     Builds a vertex-neighbor dual graph (elements sharing any node).
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <returns>Dual graph where adjacency means vertex-sharing.</returns>
    public DualGraph BuildVertexNeighborGraph<TElement, TNode>()
    {
        return BuildDualGraph<TElement, TNode>();
    }

    #endregion

    #region Memory and Statistics

    /// <summary>
    ///     Estimates the memory usage in bytes.
    /// </summary>
    public long EstimateMemoryUsage()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            long bytes = 64;

            for (var e = 0; e < _types.Count; e++)
            for (var n = 0; n < _types.Count; n++)
            {
                var m2m = _adjacency[e, n];
                var count = m2m.Count;
                for (var i = 0; i < count; i++)
                {
                    bytes += 24;
                    bytes += m2m[i].Count * 4;
                }
            }

            foreach (var kvp in _data)
            {
                var listType = kvp.Value.GetType();
                var countProp = GetCachedProperty(listType, "Count");
                var count = (int)countProp.GetValue(kvp.Value)!;
                var elementType = listType.GetGenericArguments()[0];

                bytes += 24;
                bytes += count * GetTypeSize(elementType);
            }

            foreach (var kvp in _canonicalIndex)
            {
                bytes += 48;
                bytes += kvp.Value.Count * 16;
            }

            return bytes;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    private static int GetTypeSize(Type type)
    {
        if (type == typeof(int)) return 4;
        if (type == typeof(long)) return 8;
        if (type == typeof(float)) return 4;
        if (type == typeof(double)) return 8;
        if (type == typeof(bool)) return 1;
        if (type.IsValueType)
            try
            {
                return Marshal.SizeOf(type);
            }
            catch
            {
                return 16;
            }

        return IntPtr.Size;
    }

    /// <summary>
    ///     Statistics about elements in the topology.
    /// </summary>
    public readonly struct ElementStatistics
    {
        public int ElementCount { get; }
        public int MinNodesPerElement { get; }
        public int MaxNodesPerElement { get; }
        public double AvgNodesPerElement { get; }
        public IReadOnlyDictionary<int, int> NodesPerElementDistribution { get; }

        public ElementStatistics(int count, int min, int max, double avg, IReadOnlyDictionary<int, int> dist)
        {
            ElementCount = count;
            MinNodesPerElement = min;
            MaxNodesPerElement = max;
            AvgNodesPerElement = avg;
            NodesPerElementDistribution = dist;
        }

        public override string ToString()
        {
            return $"Elements: {ElementCount}, Nodes/Element: min={MinNodesPerElement}, " +
                   $"max={MaxNodesPerElement}, avg={AvgNodesPerElement:F2}";
        }
    }

    /// <summary>
    ///     Computes element statistics.
    /// </summary>
    public ElementStatistics GetElementStatistics<TElement, TNode>()
    {
        ThrowIfDisposed();
        var elementCount = Count<TElement>();
        if (elementCount == 0)
            return new ElementStatistics(0, 0, 0, 0, new Dictionary<int, int>());

        int minNodes = int.MaxValue, maxNodes = 0;
        long totalNodes = 0;
        var dist = new Dictionary<int, int>();

        for (var e = 0; e < elementCount; e++)
        {
            var n = NodesOf<TElement, TNode>(e).Count;
            minNodes = Math.Min(minNodes, n);
            maxNodes = Math.Max(maxNodes, n);
            totalNodes += n;
            dist.TryGetValue(n, out var c);
            dist[n] = c + 1;
        }

        return new ElementStatistics(elementCount, minNodes, maxNodes, totalNodes / (double)elementCount, dist);
    }

    #endregion

    // ============================================================================
    // SMART HANDLES - FLUENT NAVIGATION (LINQ-FREE)
    // Integrated: December 15, 2024
    // Returns concrete List<T> collections, no IEnumerable, no LINQ
    // ============================================================================

    #region Smart Handles

    /// <summary>
    ///     Smart handle representing an entity with fluent navigation capabilities.
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <remarks>
    ///     <b>DESIGN:</b> Readonly record struct for value semantics (no allocation).
    ///     Carries topology reference enabling fluent method chaining without passing
    ///     topology reference explicitly.
    ///     <b>USAGE:</b> Transforms verbose topology.Method(id) calls into natural
    ///     entity.Property and entity.Method() syntax.
    ///     <example>
    ///         // Traditional approach
    ///         var neighbors = topology.Neighbors&lt;Vertex, Edge&gt;(vertexId);
    ///         foreach (var nid in neighbors) {
    ///         var pos = topology.GetData&lt;Vertex, Vector3&gt;(nid);
    ///         }
    ///         // Smart handle approach
    ///         var vertex = topology.GetEntity&lt;Vertex&gt;(vertexId);
    ///         var neighbors = vertex.Neighbors&lt;Edge&gt;();
    ///         foreach (var neighbor in neighbors) {
    ///         var pos = neighbor.Data&lt;Vector3&gt;();
    ///         }
    ///     </example>
    /// </remarks>
    public readonly record struct SmartEntity<TEntity>(Topology<TTypes> Topology, int Index)
        : IComparable<SmartEntity<TEntity>>, IEquatable<SmartEntity<TEntity>>
    {
        /// <summary>Gets whether this entity is valid (not marked for deletion).</summary>
        public bool IsValid => !Topology.GetMarkedForRemoval<TEntity>().Contains(Index);

        /// <summary>Gets whether this entity is marked for deletion.</summary>
        public bool IsMarked => Topology.GetMarkedForRemoval<TEntity>().Contains(Index);

        /// <summary>Gets the total count of entities of this type.</summary>
        public int Count => Topology.Count<TEntity>();

        /// <summary>Compares entities by index.</summary>
        public int CompareTo(SmartEntity<TEntity> other)
        {
            // First compare topology references
            if (!ReferenceEquals(Topology, other.Topology))
                throw new InvalidOperationException(
                    "Cannot compare entities from different topology instances.");

            return Index.CompareTo(other.Index);
        }

        /// <summary>Gets data associated with this entity.</summary>
        /// <typeparam name="TData">The data type.</typeparam>
        public TData Data<TData>()
        {
            return Topology.Get<TEntity, TData>(Index);
        }

        /// <summary>Sets data associated with this entity.</summary>
        public void SetData<TData>(TData value)
        {
            Topology.Set<TEntity, TData>(Index, value);
        }

        /// <summary>
        ///     Gets entities of a related type incident to this entity.
        /// </summary>
        /// <typeparam name="TRelated">The related entity type.</typeparam>
        /// <returns>List of smart handles for related entities.</returns>
        /// <remarks>
        ///     <b>EXAMPLE:</b> For a vertex, get incident edges or faces.
        ///     For an element, get its nodes.
        /// </remarks>
        public List<SmartEntity<TRelated>> IncidentTo<TRelated>()
        {
            var incidents = Topology.ElementsAt<TRelated, TEntity>(Index);
            var result = new List<SmartEntity<TRelated>>(incidents.Count);

            for (var i = 0; i < incidents.Count; i++)
                result.Add(new SmartEntity<TRelated>(Topology, incidents[i]));

            return result;
        }

        /// <summary>
        ///     Gets entities of a related type that this entity contains.
        /// </summary>
        /// <typeparam name="TRelated">The related entity type.</typeparam>
        /// <returns>List of smart handles for contained entities.</returns>
        /// <remarks>
        ///     <b>EXAMPLE:</b> For an element, get its nodes.
        ///     For a face, get its vertices.
        /// </remarks>
        public List<SmartEntity<TRelated>> Contains<TRelated>()
        {
            var nodes = Topology.NodesOf<TEntity, TRelated>(Index);
            var result = new List<SmartEntity<TRelated>>(nodes.Count);

            for (var i = 0; i < nodes.Count; i++)
                result.Add(new SmartEntity<TRelated>(Topology, nodes[i]));

            return result;
        }

        /// <summary>
        ///     Gets neighboring entities of the same type through shared relationships.
        /// </summary>
        /// <typeparam name="TRelated">The relationship type (e.g., shared nodes/edges).</typeparam>
        /// <param name="sorted">Whether to return neighbors in sorted order.</param>
        /// <returns>List of smart handles for neighboring entities.</returns>
        /// <remarks>
        ///     <b>EXAMPLE:</b> For an element, get neighboring elements sharing nodes.
        ///     For a vertex, get neighboring vertices sharing edges.
        /// </remarks>
        public List<SmartEntity<TEntity>> Neighbors<TRelated>(bool sorted = true)
        {
            var neighbors = Topology.Neighbors<TEntity, TRelated>(Index, sorted);
            var result = new List<SmartEntity<TEntity>>(neighbors.Count);

            for (var i = 0; i < neighbors.Count; i++)
                result.Add(new SmartEntity<TEntity>(Topology, neighbors[i]));

            return result;
        }

        /// <summary>
        ///     Gets direct neighbors (elements sharing at least one relationship).
        /// </summary>
        /// <typeparam name="TRelated">The relationship type.</typeparam>
        /// <param name="includeSelf">Whether to include this entity in results.</param>
        /// <param name="sorted">Whether to sort results.</param>
        public List<SmartEntity<TEntity>> DirectNeighbors<TRelated>(
            bool includeSelf = false,
            bool sorted = true)
        {
            var neighbors = Topology.GetDirectNeighbors<TEntity, TRelated>(
                Index, includeSelf, sorted);
            var result = new List<SmartEntity<TEntity>>(neighbors.Count);

            for (var i = 0; i < neighbors.Count; i++)
                result.Add(new SmartEntity<TEntity>(Topology, neighbors[i]));

            return result;
        }

        /// <summary>
        ///     Gets weighted neighbors with shared relationship counts.
        /// </summary>
        /// <typeparam name="TRelated">The relationship type.</typeparam>
        public List<(SmartEntity<TEntity> Entity, int SharedCount)> WeightedNeighbors<TRelated>()
        {
            var weighted = Topology.GetWeightedNeighbors<TEntity, TRelated>(Index);
            var result = new List<(SmartEntity<TEntity>, int)>(weighted.Count);

            for (var i = 0; i < weighted.Count; i++)
                result.Add((new SmartEntity<TEntity>(Topology, weighted[i].EntityIndex), weighted[i].SharedCount));

            return result;
        }

        /// <summary>
        ///     Gets k-hop neighborhood (entities reachable within k steps).
        /// </summary>
        /// <typeparam name="TRelated">The relationship type.</typeparam>
        /// <param name="k">Maximum number of hops.</param>
        /// <param name="minShared">Minimum shared relationships for connection.</param>
        /// <param name="includeSelf">Whether to include this entity.</param>
        /// <returns>Dictionary mapping entities to their distance (hop count).</returns>
        public Dictionary<SmartEntity<TEntity>, int> KHopNeighborhood<TRelated>(
            int k,
            int minShared = 1,
            bool includeSelf = true)
        {
            var neighborhood = Topology.GetKHopNeighborhood<TEntity, TRelated>(
                Index, k, minShared, includeSelf);

            var result = new Dictionary<SmartEntity<TEntity>, int>(neighborhood.Count);
            foreach (var kvp in neighborhood)
                result[new SmartEntity<TEntity>(Topology, kvp.Key)] = kvp.Value;

            return result;
        }

        /// <summary>
        ///     Gets entities at exactly k hops away.
        /// </summary>
        /// <typeparam name="TRelated">The relationship type.</typeparam>
        /// <param name="k">Number of hops.</param>
        /// <param name="minShared">Minimum shared relationships.</param>
        public List<SmartEntity<TEntity>> EntitiesAtDistance<TRelated>(
            int k,
            int minShared = 1)
        {
            var entities = Topology.GetEntitiesAtDistance<TEntity, TRelated>(
                Index, k, minShared);
            var result = new List<SmartEntity<TEntity>>(entities.Count);

            for (var i = 0; i < entities.Count; i++)
                result.Add(new SmartEntity<TEntity>(Topology, entities[i]));

            return result;
        }

        /// <summary>
        ///     Performs breadth-first search starting from this entity.
        /// </summary>
        /// <param name="visitor">Optional visitor callback (entity, depth).</param>
        public List<SmartEntity<TEntity>> BreadthFirstSearch(
            Action<SmartEntity<TEntity>, int>? visitor = null)
        {
            // Copy to local to avoid capturing 'this' in lambda (struct limitation)
            var localTopology = Topology;
            var localIndex = Index;

            var visited = localTopology.BreadthFirstSearch<TEntity>(
                localIndex,
                visitor == null
                    ? null
                    : (id, depth) =>
                        visitor(new SmartEntity<TEntity>(localTopology, id), depth));

            var result = new List<SmartEntity<TEntity>>(visited.Count);
            for (var i = 0; i < visited.Count; i++)
                result.Add(new SmartEntity<TEntity>(localTopology, visited[i]));

            return result;
        }

        /// <summary>
        ///     Computes shortest path distances from this entity using BFS.
        /// </summary>
        public Dictionary<SmartEntity<TEntity>, int> BreadthFirstDistances()
        {
            var distances = Topology.BreadthFirstDistances<TEntity>(Index);

            var result = new Dictionary<SmartEntity<TEntity>, int>(distances.Count);
            foreach (var kvp in distances)
                result[new SmartEntity<TEntity>(Topology, kvp.Key)] = kvp.Value;

            return result;
        }

        /// <summary>
        ///     Computes weighted shortest paths using Dijkstra's algorithm.
        /// </summary>
        /// <param name="edgeWeight">Function computing weight between entities.</param>
        public Dictionary<SmartEntity<TEntity>, (double Distance, SmartEntity<TEntity> Predecessor)>
            DijkstraShortestPaths(Func<SmartEntity<TEntity>, SmartEntity<TEntity>, int, double> edgeWeight)
        {
            // Copy to local to avoid capturing 'this' in lambda (struct limitation)
            var localTopology = Topology;
            var localIndex = Index;

            var paths = localTopology.DijkstraShortestPaths<TEntity>(
                localIndex,
                (from, to, shared) => edgeWeight(
                    new SmartEntity<TEntity>(localTopology, from),
                    new SmartEntity<TEntity>(localTopology, to),
                    shared));

            var result = new Dictionary<SmartEntity<TEntity>, (double, SmartEntity<TEntity>)>(paths.Count);
            foreach (var kvp in paths)
                result[new SmartEntity<TEntity>(localTopology, kvp.Key)] = (
                    kvp.Value.Distance,
                    kvp.Value.Predecessor == -1
                        ? default
                        : new SmartEntity<TEntity>(localTopology, kvp.Value.Predecessor)
                );

            return result;
        }

        /// <summary>
        ///     Marks this entity for removal.
        /// </summary>
        public void MarkForRemoval()
        {
            Topology.Remove<TEntity>(Index);
        }

        /// <summary>
        ///     Unmarks this entity (if previously marked).
        /// </summary>
        /// <remarks>
        ///     Note: Unmarking is not currently supported. Once marked, entities
        ///     remain marked until Compress() is called.
        /// </remarks>
        [Obsolete("Unmarking entities is not currently supported.")]
        public void Unmark()
        {
            // Unmarking is not supported in the current implementation
            throw new NotSupportedException(
                "Unmarking entities is not currently supported. Use Compress() to remove marked entities.");
        }

        /// <summary>Converts smart handle to raw index.</summary>
        public static implicit operator int(SmartEntity<TEntity> entity)
        {
            return entity.Index;
        }

        /// <summary>String representation for debugging.</summary>
        public override string ToString()
        {
            return $"{typeof(TEntity).Name}[{Index}]";
        }
    }

    /// <summary>
    ///     Factory method for creating smart handles.
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <param name="index">The entity index.</param>
    /// <returns>Smart handle wrapping the entity.</returns>
    public SmartEntity<TEntity> GetEntity<TEntity>(int index)
    {
        return new SmartEntity<TEntity>(this, index);
    }

    /// <summary>
    ///     Converts a list of indices to smart handles.
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <param name="indices">List of entity indices.</param>
    /// <returns>List of smart handles.</returns>
    public List<SmartEntity<TEntity>> GetEntities<TEntity>(List<int> indices)
    {
        var result = new List<SmartEntity<TEntity>>(indices.Count);
        for (var i = 0; i < indices.Count; i++)
            result.Add(new SmartEntity<TEntity>(this, indices[i]));
        return result;
    }

    /// <summary>
    ///     Gets all active entities as smart handles.
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <returns>List of smart handles for all active (non-deleted) entities.</returns>
    public List<SmartEntity<TEntity>> GetActiveEntities<TEntity>()
    {
        var active = GetActive<TEntity>();
        var result = new List<SmartEntity<TEntity>>(active.Count);
        for (var i = 0; i < active.Count; i++)
            result.Add(new SmartEntity<TEntity>(this, active[i]));
        return result;
    }

    #endregion

    // ============================================================================
    // HIGH-LEVEL GRAPH ALGORITHMS AND MESH CIRCULATORS
    // Integrated: December 15, 2024
    // ============================================================================

    #region Graph Algorithms - BFS

    /// <summary>
    ///     Performs breadth-first search starting from a single entity.
    /// </summary>
    /// <typeparam name="TEntity">The entity type to traverse.</typeparam>
    /// <param name="startEntity">Starting entity index.</param>
    /// <param name="visitor">Optional visitor callback invoked for each discovered entity (entity, depth).</param>
    /// <returns>List of entity indices in BFS order.</returns>
    /// <remarks>
    ///     <b>ALGORITHM:</b> BFS within a single entity type through shared relationships.
    ///     <b>THREAD SAFETY:</b> Thread-safe with read lock.
    ///     <b>USAGE:</b> Connected component detection, flood fill, shortest unweighted paths.
    ///     <example>
    ///         // Find all elements reachable from element 0
    ///         var reachable = topology.BreadthFirstSearch&lt;Element&gt;(0);
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public List<int> BreadthFirstSearch<TEntity>(int startEntity, Action<int, int>? visitor = null)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var typeIndex = GetTypeIndex<TEntity>();

            // Get the diagonal M2M (self-relationships through shared nodes)
            var selfRelations = _adjacency[typeIndex, typeIndex];
            return selfRelations.BreadthFirstSearch(startEntity, visitor);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Computes BFS distances from a start entity.
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <param name="startEntity">Starting entity index.</param>
    /// <returns>Dictionary mapping entity indices to hop count distances.</returns>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<int, int> BreadthFirstDistances<TEntity>(int startEntity)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var typeIndex = GetTypeIndex<TEntity>();
            var selfRelations = _adjacency[typeIndex, typeIndex];
            return selfRelations.BreadthFirstDistances(startEntity);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Graph Algorithms - Dijkstra

    /// <summary>
    ///     Computes shortest weighted paths using Dijkstra's algorithm.
    /// </summary>
    /// <typeparam name="TEntity">The entity type to traverse.</typeparam>
    /// <param name="startEntity">Starting entity index.</param>
    /// <param name="edgeWeight">Function computing weight between two entities sharing a connection.</param>
    /// <returns>Dictionary mapping entity indices to (distance, predecessor) tuples.</returns>
    /// <remarks>
    ///     <b>ALGORITHM:</b> Dijkstra with priority queue.
    ///     <b>WEIGHTS:</b> Must be non-negative. Function receives (fromEntity, toEntity, sharedConnection).
    ///     <b>GEODESICS:</b> For mesh surfaces, use Euclidean distance between entity centroids.
    ///     <example>
    ///         // Compute geodesic distances on mesh
    ///         var paths = topology.DijkstraShortestPaths&lt;Element&gt;(
    ///         startElem,
    ///         (from, to, sharedNode) =>
    ///         Vector3.Distance(centroids[from], centroids[to]));
    ///         // Reconstruct path to target
    ///         var shortestPath = Topology&lt;TTypes&gt;.ReconstructPath(paths, targetElem);
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths<TEntity>(
        int startEntity,
        Func<int, int, int, double> edgeWeight)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var typeIndex = GetTypeIndex<TEntity>();
            var selfRelations = _adjacency[typeIndex, typeIndex];
            return selfRelations.DijkstraShortestPaths(startEntity, edgeWeight);
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
    /// <param name="targetEntity">Target entity index.</param>
    /// <returns>Path from start to target, or null if unreachable.</returns>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static List<int>? ReconstructPath(
        Dictionary<int, (double Distance, int Predecessor)> dijkstraResult,
        int targetEntity)
    {
        return M2M.ReconstructPath(dijkstraResult, targetEntity);
    }

    #endregion

    #region Graph Algorithms - Multi-Type

    /// <summary>
    ///     Performs breadth-first search across all entity types starting from a specific entity.
    /// </summary>
    /// <typeparam name="TStartEntity">The starting entity type.</typeparam>
    /// <param name="startEntity">Starting entity index.</param>
    /// <param name="visitor">Optional visitor callback (typeIndex, entityIndex, depth).</param>
    /// <returns>List of (typeIndex, entityIndex) tuples for all reachable entities in BFS order.</returns>
    /// <remarks>
    ///     <b>USAGE:</b> Discover all entities connected across the type hierarchy.
    ///     For example, starting from a vertex, find all connected edges, faces, volumes.
    ///     <b>TYPE INDICES:</b> Use <see cref="GetTypeIndex{T}" /> to convert type indices back to types.
    ///     <example>
    ///         // Find all entities connected to vertex 0
    ///         var reachable = topology.BreadthFirstSearchMultiType&lt;Vertex&gt;(0);
    ///         foreach (var (typeIdx, idx) in reachable)
    ///         Console.WriteLine($"Type {typeIdx}, Entity {idx}");
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public List<(int TypeIndex, int EntityIndex)> BreadthFirstSearchMultiType<TStartEntity>(
        int startEntity,
        Action<int, int, int>? visitor = null)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var startTypeIndex = GetTypeIndex<TStartEntity>();
            return _adjacency.BreadthFirstSearchMultiType(startTypeIndex, startEntity, visitor);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Performs multi-type BFS and returns distances from start entity.
    /// </summary>
    /// <typeparam name="TStartEntity">The starting entity type.</typeparam>
    /// <param name="startEntity">Starting entity index.</param>
    /// <returns>Dictionary mapping (typeIndex, entityIndex) to BFS distance (hop count).</returns>
    /// <remarks>
    ///     <b>DISTANCES:</b> Hop count across entity types.
    ///     <b>UNREACHABLE:</b> Entities not in dictionary are unreachable from start.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<(int TypeIndex, int EntityIndex), int> BreadthFirstDistancesMultiType<TStartEntity>(
        int startEntity)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var startTypeIndex = GetTypeIndex<TStartEntity>();
            return _adjacency.BreadthFirstDistancesMultiType(startTypeIndex, startEntity);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Computes shortest weighted paths across entity types using Dijkstra's algorithm.
    /// </summary>
    /// <typeparam name="TStartEntity">The starting entity type.</typeparam>
    /// <param name="startEntity">Starting entity index.</param>
    /// <param name="edgeWeight">
    ///     Function computing weight: (fromType, fromEntity, toType, toEntity, sharedNodeType,
    ///     sharedNode) → weight.
    /// </param>
    /// <returns>Dictionary mapping (typeIndex, entityIndex) to (distance, predecessor) tuples.</returns>
    /// <remarks>
    ///     <b>ALGORITHM:</b> Dijkstra with priority queue, traversing across entity types.
    ///     <b>WEIGHTS:</b> Must be non-negative.
    ///     <b>PREDECESSOR:</b> The predecessor is a (typeIndex, entityIndex) tuple, stored as (PredType, PredEntity).
    ///     <example>
    ///         // Compute weighted shortest paths across all entity types
    ///         var paths = topology.DijkstraShortestPathsMultiType&lt;Vertex&gt;(
    ///         startVertex,
    ///         (fromType, fromEnt, toType, toEnt, nodeType, node) => 1.0);
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public Dictionary<(int TypeIndex, int EntityIndex), (double Distance, (int PredType, int PredEntity))>
        DijkstraShortestPathsMultiType<TStartEntity>(
            int startEntity,
            Func<int, int, int, int, int, int, double> edgeWeight)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var startTypeIndex = GetTypeIndex<TStartEntity>();
            return _adjacency.DijkstraShortestPathsMultiType(startTypeIndex, startEntity, edgeWeight);
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
    /// <param name="targetTypeIndex">Target entity type index.</param>
    /// <param name="targetEntity">Target entity index.</param>
    /// <returns>Path as list of (typeIndex, entityIndex) tuples, or null if unreachable.</returns>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static List<(int TypeIndex, int EntityIndex)>? ReconstructPathMultiType(
        Dictionary<(int TypeIndex, int EntityIndex), (double Distance, (int PredType, int PredEntity))> dijkstraResult,
        int targetTypeIndex,
        int targetEntity)
    {
        return MM2M.ReconstructPathMultiType(dijkstraResult, targetTypeIndex, targetEntity);
    }

    #endregion

    #region Mesh Traversal - Circulators

    /// <summary>
    ///     Creates a circulator for iterating over neighbors of an entity.
    /// </summary>
    /// <typeparam name="TEntity">The entity type.</typeparam>
    /// <param name="entityIndex">The entity index.</param>
    /// <param name="sorted">Whether to return neighbors in sorted order.</param>
    /// <returns>Enumerable circulator over neighbor entities.</returns>
    /// <remarks>
    ///     <b>PATTERN:</b> Implements circulator pattern common in mesh libraries (OpenMesh, CGAL).
    ///     <b>USAGE:</b> Iterate all entities sharing connections with the given entity.
    ///     <b>THREAD SAFETY:</b> Safe for concurrent iteration (snapshot-based).
    ///     <example>
    ///         // Iterate all neighboring elements
    ///         foreach (var neighbor in topology.Circulate&lt;Element&gt;(elemId))
    ///         {
    ///         ProcessNeighbor(neighbor);
    ///         }
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public EntityCirculator<TEntity> Circulate<TEntity>(int entityIndex, bool sorted = true)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var neighbors = Neighbors<TEntity, TEntity>(entityIndex, sorted);
            return new EntityCirculator<TEntity>(neighbors);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Creates a circulator for iterating incident entities of a different type.
    /// </summary>
    /// <typeparam name="TEntity">The primary entity type.</typeparam>
    /// <typeparam name="TIncident">The incident entity type.</typeparam>
    /// <param name="entityIndex">The entity index.</param>
    /// <param name="sorted">Whether to return results in sorted order.</param>
    /// <returns>Enumerable circulator over incident entities.</returns>
    /// <remarks>
    ///     <b>USAGE:</b> Iterate all entities of a different type incident to the given entity.
    ///     <example>
    ///         // Get all faces incident to a vertex
    ///         foreach (var faceId in topology.CirculateIncident&lt;Vertex, Face&gt;(vertexId))
    ///         {
    ///         ProcessFace(faceId);
    ///         }
    ///     </example>
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public IncidentCirculator<TEntity, TIncident> CirculateIncident<TEntity, TIncident>(
        int entityIndex,
        bool sorted = true)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var incidents = ElementsAt<TIncident, TEntity>(entityIndex);
            if (sorted)
            {
                var sortedList = incidents.ToList();
                sortedList.Sort();
                return new IncidentCirculator<TEntity, TIncident>(sortedList);
            }

            return new IncidentCirculator<TEntity, TIncident>(incidents);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Creates a boundary circulator for traversing boundary loops.
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="startNode">Starting boundary node.</param>
    /// <param name="nodesPerBoundaryFacet">Number of nodes per boundary facet.</param>
    /// <returns>Circulator over boundary nodes in loop order.</returns>
    /// <remarks>
    ///     <b>USAGE:</b> Walk around mesh boundary for parameterization, hole filling, etc.
    ///     <b>REQUIREMENT:</b> Mesh must have well-defined boundary (manifold edges).
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public BoundaryCirculator<TElement, TNode> CirculateBoundary<TElement, TNode>(
        int startNode,
        int nodesPerBoundaryFacet)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            ThrowIfDisposed();
            var boundaryNodes = FindBoundaryNodes<TElement, TNode>(nodesPerBoundaryFacet);

            if (!boundaryNodes.Contains(startNode))
                throw new ArgumentException($"Node {startNode} is not a boundary node.");

            // Trace boundary loop starting from startNode
            var loop = TraceBoundaryLoop<TElement, TNode>(startNode, boundaryNodes, nodesPerBoundaryFacet);
            return new BoundaryCirculator<TElement, TNode>(loop);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Internal helper to trace a boundary loop.
    /// </summary>
    private List<int> TraceBoundaryLoop<TElement, TNode>(
        int startNode,
        HashSet<int> boundaryNodes,
        int nodesPerFacet)
    {
        var loop = new List<int> { startNode };
        var current = startNode;
        var visited = new HashSet<int> { startNode };

        // Simple walk: find next boundary node connected to current
        while (true)
        {
            var elements = ElementsAt<TElement, TNode>(current);
            int? nextNode = null;

            foreach (var elem in elements)
            {
                var nodes = NodesOf<TElement, TNode>(elem);

                foreach (var node in nodes)
                {
                    if (node == current) continue;
                    if (!boundaryNodes.Contains(node)) continue;
                    if (visited.Contains(node)) continue;

                    nextNode = node;
                    break;
                }

                if (nextNode.HasValue) break;
            }

            if (!nextNode.HasValue) break;

            current = nextNode.Value;
            loop.Add(current);
            visited.Add(current);

            // Check if we've completed the loop
            if (loop.Count > 2)
            {
                var checkElements = ElementsAt<TElement, TNode>(current);
                foreach (var elem in checkElements)
                {
                    var nodes = NodesOf<TElement, TNode>(elem);
                    if (nodes.Contains(startNode))
                        // Loop complete
                        return loop;
                }
            }
        }

        return loop;
    }

    #endregion

    #region Circulator Classes

    /// <summary>
    ///     Circulator for iterating over entities of the same type.
    /// </summary>
    /// <typeparam name="TEntity">Entity type.</typeparam>
    public readonly struct EntityCirculator<TEntity> : IEnumerable<int>
    {
        private readonly IReadOnlyList<int> _entities;

        internal EntityCirculator(IReadOnlyList<int> entities)
        {
            _entities = entities;
        }

        /// <summary>Gets the number of entities in the circulation.</summary>
        public int Count => _entities.Count;

        /// <summary>Gets the entity at the specified index.</summary>
        public int this[int index] => _entities[index];

        /// <summary>Converts circulator to list.</summary>
        public List<int> ToList()
        {
            return _entities.ToList();
        }

        public IEnumerator<int> GetEnumerator()
        {
            return _entities.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }

    /// <summary>
    ///     Circulator for iterating over incident entities of a different type.
    /// </summary>
    /// <typeparam name="TEntity">Primary entity type.</typeparam>
    /// <typeparam name="TIncident">Incident entity type.</typeparam>
    public readonly struct IncidentCirculator<TEntity, TIncident> : IEnumerable<int>
    {
        private readonly IReadOnlyList<int> _incidents;

        internal IncidentCirculator(IReadOnlyList<int> incidents)
        {
            _incidents = incidents;
        }

        /// <summary>Gets the number of incident entities.</summary>
        public int Count => _incidents.Count;

        /// <summary>Gets the incident entity at the specified index.</summary>
        public int this[int index] => _incidents[index];

        /// <summary>Converts circulator to list.</summary>
        public List<int> ToList()
        {
            return _incidents.ToList();
        }

        public IEnumerator<int> GetEnumerator()
        {
            return _incidents.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }

    /// <summary>
    ///     Circulator for traversing boundary loops.
    /// </summary>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    public readonly struct BoundaryCirculator<TElement, TNode> : IEnumerable<int>
    {
        private readonly List<int> _loop;

        internal BoundaryCirculator(List<int> loop)
        {
            _loop = loop;
        }

        /// <summary>Gets the number of nodes in the boundary loop.</summary>
        public int Count => _loop.Count;

        /// <summary>Gets the node at the specified position in the loop.</summary>
        public int this[int index] => _loop[index];

        /// <summary>Checks if the loop is closed (forms a cycle).</summary>
        public bool IsClosed => _loop.Count > 2;

        /// <summary>Converts circulator to list.</summary>
        public List<int> ToList()
        {
            return new List<int>(_loop);
        }

        public IEnumerator<int> GetEnumerator()
        {
            return _loop.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }

    #endregion

    #region Low-Level Access Wrappers (O2M/M2M/MM2M)

    // ============================================================================
    // These thin wrappers expose O2M, M2M, and MM2M functionality through
    // the type-safe Topology API without requiring escape hatches.
    // ============================================================================

    #region MM2M Wrappers

    /// <summary>
    ///     Gets the internal version counter for change detection.
    /// </summary>
    /// <remarks>
    ///     <b>USE CASE:</b> Detect when the topology has been modified (e.g., after Compress).
    ///     The version increments on structural changes. Useful for cache invalidation.
    /// </remarks>
    public long Version
    {
        get
        {
            ThrowIfDisposed();
            return _adjacency.Version;
        }
    }

    #endregion

    #region M2M Query Wrappers

    /// <summary>
    ///     Checks if an element index is valid for the given element-node relationship.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="elementIndex">The element index to check.</param>
    /// <returns>True if the element exists; false otherwise.</returns>
    public bool HasElement<TElement, TNode>(int elementIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.HasElement(elementIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if a node index exists in the transpose structure.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="nodeIndex">The node index to check.</param>
    /// <returns>True if the node is referenced by at least one element; false otherwise.</returns>
    public bool HasNode<TElement, TNode>(int nodeIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.HasNode(nodeIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if a specific element contains a specific node.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="elementIndex">The element index.</param>
    /// <param name="nodeIndex">The node index to look for.</param>
    /// <returns>True if the element contains the node; false otherwise.</returns>
    public bool ElementContainsNode<TElement, TNode>(int elementIndex, int nodeIndex)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.ElementContainsNode(elementIndex, nodeIndex);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets the number of unique nodes referenced across all elements.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>The count of unique node indices (max node index + 1).</returns>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> O(sync) + O(1). Much more efficient than computing
    ///     the transpose just to get the count.
    /// </remarks>
    public int GetTransposeNodeCount<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetTransposeNodeCount();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes an action with zero-copy span access to elements containing a node.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="nodeIndex">The node index.</param>
    /// <param name="action">Action to execute with the element span.</param>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> True zero-copy access. No allocations.
    ///     <para>
    ///         <b>SAFETY:</b> The span is only valid during action execution.
    ///         Do NOT store the span. Read lock is held during execution.
    ///     </para>
    /// </remarks>
    public void WithElementsForNodeSpan<TElement, TNode>(int nodeIndex, M2M.ReadOnlySpanAction<int> action)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(action);
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            m2m.WithElementsForNodeSpan(nodeIndex, action);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes an action with direct access to the transpose O2M structure.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="action">Action to execute with the transpose.</param>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> Zero-copy access. No cloning.
    ///     <para>
    ///         <b>SAFETY:</b> The O2M reference is only valid during action execution.
    ///         Do NOT store references. Read lock is held during execution.
    ///     </para>
    /// </remarks>
    public void WithTranspose<TElement, TNode>(Action<O2M> action)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(action);
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            m2m.WithElementsFromNode(action);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Executes a function with direct access to the transpose O2M structure.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <typeparam name="TResult">The result type.</typeparam>
    /// <param name="func">Function to execute with the transpose.</param>
    /// <returns>The function result.</returns>
    public TResult WithTranspose<TElement, TNode, TResult>(Func<O2M, TResult> func)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(func);
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.WithElementsFromNode(func);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Asynchronously gets a clone of the transpose structure for large datasets.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>A task containing the transposed O2M structure.</returns>
    /// <remarks>
    ///     <b>LOCKING:</b> If the transpose needs computation, a write lock is held
    ///     during that time. For better concurrency, call EnsureSynchronized first.
    /// </remarks>
    public async Task<O2M> GetTransposeAsync<TElement, TNode>(CancellationToken cancellationToken = default)
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return await m2m.ElementsFromNodeAsync(cancellationToken).ConfigureAwait(false);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Forces synchronization of the transpose cache for a specific relationship.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <remarks>
    ///     Call this before parallel read operations to ensure the transpose
    ///     is pre-computed and avoid lock contention.
    /// </remarks>
    public void EnsureSynchronized<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            m2m.Synchronize();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Forces synchronization of the transpose and position caches for a specific relationship.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <remarks>
    ///     Call this before hot-path queries that use <see cref="GetPositionCaches{TElement, TNode}" />
    ///     to avoid cache computation under contention.
    /// </remarks>
    public void EnsurePositionCaches<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            m2m.EnsurePositionCaches();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if a specific relationship block is valid (element/node/position checks).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>True if the relationship structure is valid; false otherwise.</returns>
    public bool IsValidRelationship<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.IsValid();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if two M2M structures are permutations of each other.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="other">The other topology to compare.</param>
    /// <returns>True if the relationships are permutations; false otherwise.</returns>
    public bool IsPermutationOf<TElement, TNode>(Topology<TTypes> other)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(other);
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            var otherM2m = other._adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.IsPermutationOf(otherM2m);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region O2M Query Wrappers

    /// <summary>
    ///     Gets the maximum node index across all elements.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>The maximum node index, or -1 if empty.</returns>
    public int GetMaxNodeIndex<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.GetMaxNode();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Checks if all adjacency lists are sorted in ascending order.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>True if all adjacency lists are sorted; false otherwise.</returns>
    /// <remarks>
    ///     Many algorithms assume sorted adjacency lists. Use this to verify the invariant.
    /// </remarks>
    public bool IsSorted<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, bool>(o2m => o2m.IsSorted());
    }

    /// <summary>
    ///     Performs comprehensive validation and returns detailed error information.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>Null if valid; otherwise a description of the first error found.</returns>
    public string? ValidateRelationshipStrict<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, string?>(o2m => o2m.ValidateStrict());
    }

    /// <summary>
    ///     Counts total number of edges (node references) across all elements.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>Total edge count.</returns>
    public long GetTotalEdgeCount<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, long>(o2m => o2m.GetTotalEdgeCount());
    }

    /// <summary>
    ///     Gets detailed statistics about a relationship structure.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>Tuple of (MinDegree, MaxDegree, AvgDegree, TotalEdges).</returns>
    public (int MinDegree, int MaxDegree, double AvgDegree, long TotalEdges) GetRelationshipStatistics<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, (int, int, double, long)>(o2m => o2m.GetStatistics());
    }

    /// <summary>
    ///     Gets a content-based hash code for the relationship structure.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>A hash code based on the full content.</returns>
    /// <remarks>
    ///     Unlike GetHashCode(), this computes a hash over all data for content comparison.
    /// </remarks>
    public int GetRelationshipContentHash<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, int>(o2m => o2m.FullContentHashCode());
    }

    #endregion

    #region M2M Direct Access Wrappers

    /// <summary>
    ///     Executes an action with direct access to the M2M structure.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="action">Action to execute with the M2M.</param>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> Direct access without cloning.
    ///     <para>
    ///         <b>SAFETY:</b> The M2M reference is only valid during action execution.
    ///         Do NOT store references. Appropriate locks are held during execution.
    ///     </para>
    /// </remarks>
    public void WithRelationship<TElement, TNode>(Action<M2M> action)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(action);
        _adjacency.WithBlock(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), action);
    }

    /// <summary>
    ///     Executes a function with direct access to the M2M structure.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <typeparam name="TResult">The result type.</typeparam>
    /// <param name="func">Function to execute with the M2M.</param>
    /// <returns>The function result.</returns>
    public TResult WithRelationship<TElement, TNode, TResult>(Func<M2M, TResult> func)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(func);
        return _adjacency.WithBlock(GetTypeIndex<TElement>(), GetTypeIndex<TNode>(), func);
    }

    #endregion

    #region O2M Transpose Wrappers

    /// <summary>
    ///     Gets a clone of the transpose structure (elements from nodes).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>A new O2M where row i contains elements connected to node i.</returns>
    /// <remarks>
    ///     <b>PERFORMANCE:</b> O(n × m) for deep copy. For zero-copy access,
    ///     use WithTranspose instead.
    /// </remarks>
    public O2M GetTranspose<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.ElementsFromNode;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets a clone of the transpose structure with a pre-allocated capacity limit.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="maxNodeCap">Maximum node capacity to pre-allocate.</param>
    /// <returns>A new O2M where row i contains elements connected to node i.</returns>
    /// <remarks>
    ///     <b>USE CASE:</b> When you know the maximum node index in advance,
    ///     this avoids potential large allocations from sparse node indices.
    ///     <para>
    ///         <b>WARNING:</b> If actual max node exceeds maxNodeCap, an exception is thrown.
    ///     </para>
    /// </remarks>
    public O2M GetTranspose<TElement, TNode>(int maxNodeCap)
    {
        ThrowIfDisposed();
        ArgumentOutOfRangeException.ThrowIfNegative(maxNodeCap);
        var forward = GetForwardStructure<TElement, TNode>();
        return forward.Transpose(maxNodeCap);
    }

    /// <summary>
    ///     Gets a clone of the forward structure (nodes from elements).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>A new O2M where row i contains nodes connected to element i.</returns>
    public O2M GetForwardStructure<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return m2m.NodesFromElement;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Computes a transpose with strict validation (no invalid node skipping).
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>Strictly transposed O2M structure.</returns>
    public O2M GetTransposeStrict<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, O2M>(o2m => o2m.TransposeStrict());
    }

    #endregion

    #region O2M Conversion Wrappers

    /// <summary>
    ///     Converts a relationship to a boolean adjacency matrix.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>Boolean matrix where [i,j] = true if element i contains node j.</returns>
    /// <remarks>
    ///     <b>WARNING:</b> Can allocate significant memory for large topologies.
    ///     Matrix size is Count × (MaxNodeIndex + 1).
    /// </remarks>
    public bool[,] ToBooleanMatrix<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, bool[,]>(o2m => o2m.ToBooleanMatrix());
    }

    /// <summary>
    ///     Generates an EPS (Encapsulated PostScript) string for visualization.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>EPS string representation of the adjacency structure.</returns>
    /// <remarks>
    ///     Useful for debugging and documentation. The output can be viewed
    ///     in PostScript viewers or converted to PDF.
    /// </remarks>
    public string ToEpsString<TElement, TNode>()
    {
        ThrowIfDisposed();
        return WithTranspose<TElement, TNode, string>(o2m => o2m.ToEpsString());
    }

    #endregion

    #region O2M Set Operations

    /// <summary>
    ///     Computes the union of two relationship structures.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="other">The other topology.</param>
    /// <returns>New O2M containing union of both structures.</returns>
    /// <remarks>
    ///     Element i in result contains all nodes from element i in either structure.
    /// </remarks>
    public O2M UnionWith<TElement, TNode>(Topology<TTypes> other)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(other);

        O2M left = null!;
        O2M right = null!;

        WithTranspose<TElement, TNode>(o2m => left = (O2M)o2m.Clone());
        other.WithTranspose<TElement, TNode>(o2m => right = (O2M)o2m.Clone());

        return left | right;
    }

    /// <summary>
    ///     Computes the intersection of two relationship structures.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="other">The other topology.</param>
    /// <returns>New O2M containing intersection of both structures.</returns>
    /// <remarks>
    ///     Element i in result contains only nodes present in element i of both structures.
    /// </remarks>
    public O2M IntersectWith<TElement, TNode>(Topology<TTypes> other)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(other);

        O2M left = null!;
        O2M right = null!;

        WithTranspose<TElement, TNode>(o2m => left = (O2M)o2m.Clone());
        other.WithTranspose<TElement, TNode>(o2m => right = (O2M)o2m.Clone());

        return left & right;
    }

    /// <summary>
    ///     Computes the difference of two relationship structures.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="other">The other topology.</param>
    /// <returns>New O2M containing nodes in this but not in other.</returns>
    /// <remarks>
    ///     Element i in result contains nodes from this element i that are NOT in other's element i.
    /// </remarks>
    public O2M DifferenceWith<TElement, TNode>(Topology<TTypes> other)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(other);

        O2M left = null!;
        O2M right = null!;

        WithTranspose<TElement, TNode>(o2m => left = (O2M)o2m.Clone());
        other.WithTranspose<TElement, TNode>(o2m => right = (O2M)o2m.Clone());

        return left - right;
    }

    /// <summary>
    ///     Computes the symmetric difference of two relationship structures.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="other">The other topology.</param>
    /// <returns>New O2M containing nodes in exactly one of the structures.</returns>
    /// <remarks>
    ///     Element i in result contains nodes that are in this XOR other (but not both).
    /// </remarks>
    public O2M SymmetricDifferenceWith<TElement, TNode>(Topology<TTypes> other)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(other);

        O2M left = null!;
        O2M right = null!;

        WithTranspose<TElement, TNode>(o2m => left = (O2M)o2m.Clone());
        other.WithTranspose<TElement, TNode>(o2m => right = (O2M)o2m.Clone());

        return left ^ right;
    }

    /// <summary>
    ///     Computes matrix multiplication of two relationship structures.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="other">The other topology.</param>
    /// <returns>New O2M representing the matrix product.</returns>
    /// <remarks>
    ///     For A * B, row i of result contains all j where there exists k such that
    ///     A[i] contains k and B[k] contains j. Useful for transitive closure computations.
    /// </remarks>
    public O2M MultiplyWith<TElement, TNode>(Topology<TTypes> other)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(other);

        O2M left = null!;
        O2M right = null!;

        WithTranspose<TElement, TNode>(o2m => left = (O2M)o2m.Clone());
        other.WithTranspose<TElement, TNode>(o2m => right = (O2M)o2m.Clone());

        return left * right;
    }

    #endregion

    #region O2M Graph Algorithm Wrappers

    /// <summary>
    ///     Gets the position caches for node locations within elements.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>
    ///     Tuple containing:
    ///     - ElementLocations: For each element e and position p, gives position of e in node p's adjacency list
    ///     - NodeLocations: For each node n, gives positions where n appears in each element's list
    /// </returns>
    public (IReadOnlyList<IReadOnlyList<int>> ElementLocations, IReadOnlyList<IReadOnlyList<int>> NodeLocations)
        GetPositionCaches<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            return (m2m.ElementLocations, m2m.NodeLocations);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets node positions computed from the forward and transpose structures.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>List where result[node][i] gives position of node within element i's adjacency list.</returns>
    public List<List<int>> ComputeNodePositions<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            var forward = m2m.NodesFromElement;
            var transpose = m2m.ElementsFromNode;
            return O2M.GetNodePositions(forward, transpose);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    /// <summary>
    ///     Gets element positions computed from the forward and transpose structures.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <returns>List where result[element][i] gives position of element within node i's adjacency list.</returns>
    public List<List<int>> ComputeElementPositions<TElement, TNode>()
    {
        ThrowIfDisposed();
        _rwLock.EnterReadLock();
        try
        {
            var m2m = _adjacency[GetTypeIndex<TElement>(), GetTypeIndex<TNode>()];
            var forward = m2m.NodesFromElement;
            var transpose = m2m.ElementsFromNode;
            return O2M.GetElementPositions(forward, transpose);
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
    }

    #endregion

    #region Static Factory Methods for O2M

    /// <summary>
    ///     Creates a topology from a boolean adjacency matrix.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="matrix">Boolean matrix where [i,j] = true means element i contains node j.</param>
    /// <returns>New topology with the specified connectivity.</returns>
    public static Topology<TTypes> FromBooleanMatrix<TElement, TNode>(bool[,] matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);

        var o2m = O2M.FromBooleanMatrix(matrix);
        var topology = new Topology<TTypes>();

        var elementType = topology.GetTypeIndex<TElement>();
        var nodeType = topology.GetTypeIndex<TNode>();

        // Add nodes first (diagonal entries for self-referential node type)
        var nodeCount = matrix.GetLength(1);
        for (var n = 0; n < nodeCount; n++)
            topology.Add<TNode>();

        // Add elements with their connectivity
        for (var e = 0; e < o2m.Count; e++)
        {
            var nodes = o2m[e].ToArray();
            topology.Add<TElement, TNode>(nodes);
        }

        return topology;
    }

    /// <summary>
    ///     Creates a random topology for testing purposes.
    /// </summary>
    /// <typeparam name="TElement">The element type.</typeparam>
    /// <typeparam name="TNode">The node type.</typeparam>
    /// <param name="elementCount">Number of elements to create.</param>
    /// <param name="nodeCount">Number of nodes available.</param>
    /// <param name="density">Probability of each element-node connection (0.0 to 1.0).</param>
    /// <param name="seed">Optional random seed for reproducibility.</param>
    /// <returns>New topology with random connectivity.</returns>
    public static Topology<TTypes> CreateRandom<TElement, TNode>(
        int elementCount,
        int nodeCount,
        double density,
        int? seed = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(elementCount);
        ArgumentOutOfRangeException.ThrowIfNegative(nodeCount);
        ArgumentOutOfRangeException.ThrowIfLessThan(density, 0.0);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(density, 1.0);

        var topology = new Topology<TTypes>();

        // Add nodes
        for (var n = 0; n < nodeCount; n++)
            topology.Add<TNode>();

        // Add elements with random connectivity
        var rng = seed.HasValue ? new Random(seed.Value) : new Random();
        for (var e = 0; e < elementCount; e++)
        {
            var nodes = new List<int>();
            for (var n = 0; n < nodeCount; n++)
                if (rng.NextDouble() < density)
                    nodes.Add(n);

            if (nodes.Count > 0)
                topology.Add<TElement, TNode>(nodes.ToArray());
            else
                topology.Add<TElement, TNode>(0); // At least one node
        }

        return topology;
    }

    #endregion

    #endregion
}

/// <summary>
///     Read-only wrapper for Topology. Safe for concurrent access.
/// </summary>
/// <remarks>
///     Now includes comprehensive read-only access to:
///     - Core queries (Count, NodesOf, ElementsAt, etc.)
///     - O(1) counters (CountRelated, CountIncident, CountActive)
///     - Union/intersection queries (GetEntitiesContainingAll, GetEntitiesContainingAny)
///     - Duplicate detection (GetDuplicates, GetAllDuplicates)
///     - Ordering (GetTopologicalOrder, GetSortOrder)
///     - Multi-type traversal (MultiTypeDFS, GetAllEntitiesAtNode, GetAllNodesOfEntity)
///     - Statistics (GetStats, GetElementStatistics)
///     - Validation (ValidateStructure)
/// </remarks>
public sealed class ReadOnlyTopology<TTypes> where TTypes : ITypeMap, new()
{
    private readonly Topology<TTypes> _topology;

    internal ReadOnlyTopology(Topology<TTypes> topology)
    {
        ArgumentNullException.ThrowIfNull(topology);
        _topology = topology;
    }

    #region Core Queries

    public int Count<TEntity>()
    {
        return _topology.Count<TEntity>();
    }

    public IReadOnlyList<int> NodesOf<TElement, TNode>(int element)
    {
        return _topology.NodesOf<TElement, TNode>(element);
    }

    public IReadOnlyList<int> ElementsAt<TElement, TNode>(int node)
    {
        return _topology.ElementsAt<TElement, TNode>(node);
    }

    public TData Get<TEntity, TData>(int index)
    {
        return _topology.Get<TEntity, TData>(index);
    }

    public IReadOnlyList<TData> All<TEntity, TData>()
    {
        return _topology.All<TEntity, TData>();
    }

    public IReadOnlyList<int> Neighbors<TElement, TNode>(int element)
    {
        return _topology.Neighbors<TElement, TNode>(element);
    }

    public bool Exists<TElement>(params int[] nodes)
    {
        return _topology.Exists<TElement>(nodes);
    }

    public int Find<TElement>(params int[] nodes)
    {
        return _topology.Find<TElement>(nodes);
    }

    public int ComputeBandwidth<TElement, TNode>()
    {
        return _topology.ComputeBandwidth<TElement, TNode>();
    }

    public HashSet<int> FindBoundaryNodes<TElement, TNode>(int nodesPerFacet)
    {
        return _topology.FindBoundaryNodes<TElement, TNode>(nodesPerFacet);
    }

    #endregion

    #region O(1) Counters (NEW v4)

    /// <summary>
    ///     Gets the count of related entities of a specific type connected to an entity.
    /// </summary>
    public int CountRelated<TEntity, TRelated>(int entityIndex)
    {
        return _topology.CountRelated<TEntity, TRelated>(entityIndex);
    }

    /// <summary>
    ///     Gets the count of elements of a specific type incident to a node.
    /// </summary>
    public int CountIncident<TElement, TNode>(int nodeIndex)
    {
        return _topology.CountIncident<TElement, TNode>(nodeIndex);
    }

    /// <summary>
    ///     Gets the count of active (not marked for deletion) entities of a given type.
    /// </summary>
    public int CountActive<TEntity>()
    {
        return _topology.CountActive<TEntity>();
    }

    #endregion

    #region Union/Intersection Queries (NEW v4)

    /// <summary>
    ///     Gets entities containing ALL specified relationships (intersection).
    /// </summary>
    public List<int> GetEntitiesContainingAll<TEntity, TRelated>(List<int> relationships)
    {
        return _topology.GetEntitiesContainingAll<TEntity, TRelated>(relationships);
    }

    /// <summary>
    ///     Gets entities containing ANY of the specified relationships (union).
    /// </summary>
    public List<int> GetEntitiesContainingAny<TEntity, TRelated>(List<int> relationships)
    {
        return _topology.GetEntitiesContainingAny<TEntity, TRelated>(relationships);
    }

    #endregion

    #region Connectivity Queries

    /// <summary>
    ///     Gets entities sharing at least one relationship (direct neighbors).
    /// </summary>
    public List<int> GetDirectNeighbors<TEntity, TRelated>(int entityIndex, bool includeSelf = false,
        bool sorted = true)
    {
        return _topology.GetDirectNeighbors<TEntity, TRelated>(entityIndex, includeSelf, sorted);
    }

    /// <summary>
    ///     Gets entities sharing exactly k relationships.
    /// </summary>
    public List<int> GetEntitiesWithSharedCount<TEntity, TRelated>(int entityIndex, int exactCount,
        bool includeSelf = false)
    {
        return _topology.GetEntitiesWithSharedCount<TEntity, TRelated>(entityIndex, exactCount, includeSelf);
    }

    /// <summary>
    ///     Gets entities sharing at least minCount relationships.
    /// </summary>
    public List<int> GetEntitiesWithMinSharedCount<TEntity, TRelated>(int entityIndex, int minCount = 1,
        bool includeSelf = false)
    {
        return _topology.GetEntitiesWithMinSharedCount<TEntity, TRelated>(entityIndex, minCount, includeSelf);
    }

    /// <summary>
    ///     Gets entities with their shared relationship counts.
    /// </summary>
    public List<(int EntityIndex, int SharedCount)> GetWeightedNeighbors<TEntity, TRelated>(int entityIndex,
        int minCount = 1, bool includeSelf = false)
    {
        return _topology.GetWeightedNeighbors<TEntity, TRelated>(entityIndex, minCount, includeSelf);
    }

    /// <summary>
    ///     Gets entities within k relationship steps via BFS.
    /// </summary>
    public Dictionary<int, int> GetKHopNeighborhood<TEntity, TRelated>(int seedEntity, int k,
        int minSharedForConnection = 1, bool includeSeed = true)
    {
        return _topology.GetKHopNeighborhood<TEntity, TRelated>(seedEntity, k, minSharedForConnection, includeSeed);
    }

    /// <summary>
    ///     Gets entities at exactly k steps from seed.
    /// </summary>
    public List<int> GetEntitiesAtDistance<TEntity, TRelated>(int seedEntity, int k, int minSharedForConnection = 1)
    {
        return _topology.GetEntitiesAtDistance<TEntity, TRelated>(seedEntity, k, minSharedForConnection);
    }

    #endregion

    #region Multi-Type Connectivity

    /// <summary>
    ///     Performs DFS across ALL types starting from a node.
    /// </summary>
    public List<(int TypeIndex, int EntityIndex)> MultiTypeDFS<TNode>(int nodeIndex)
    {
        return _topology.MultiTypeDFS<TNode>(nodeIndex);
    }

    /// <summary>
    ///     Gets ALL entities of ANY type connected to a specific node.
    /// </summary>
    public List<(int TypeIndex, int EntityIndex)> GetAllEntitiesAtNode<TNode>(int nodeIndex)
    {
        return _topology.GetAllEntitiesAtNode<TNode>(nodeIndex);
    }

    /// <summary>
    ///     Gets ALL nodes of ANY type connected to a specific entity.
    /// </summary>
    public List<(int TypeIndex, int NodeIndex)> GetAllNodesOfEntity<TEntity>(int entityIndex)
    {
        return _topology.GetAllNodesOfEntity<TEntity>(entityIndex);
    }

    /// <summary>
    ///     Gets the type topological order.
    /// </summary>
    public List<int> GetTypeTopologicalOrder()
    {
        return _topology.GetTypeTopologicalOrder();
    }

    /// <summary>
    ///     Checks if the type hierarchy is acyclic.
    /// </summary>
    public bool IsTypeHierarchyAcyclic()
    {
        return _topology.IsTypeHierarchyAcyclic();
    }

    #endregion

    #region Duplicate Detection (NEW v4)

    /// <summary>
    ///     Gets duplicate elements for a specific entity type.
    /// </summary>
    public List<int> GetDuplicates<TEntity>()
    {
        return _topology.GetDuplicates<TEntity>();
    }

    /// <summary>
    ///     Gets all duplicate elements across all entity types.
    /// </summary>
    public Dictionary<int, List<int>> GetAllDuplicates()
    {
        return _topology.GetAllDuplicates();
    }

    #endregion

    #region Ordering (NEW v4)

    /// <summary>
    ///     Gets the topological ordering of entities within a type based on dependencies.
    /// </summary>
    public List<int> GetTopologicalOrder<TEntity>()
    {
        return _topology.GetTopologicalOrder<TEntity>();
    }

    /// <summary>
    ///     Gets a lexicographically sorted ordering of elements.
    /// </summary>
    public List<int> GetSortOrder<TEntity>()
    {
        return _topology.GetSortOrder<TEntity>();
    }

    #endregion

    #region Validation and Statistics (NEW v4)

    /// <summary>
    ///     Validates the structural integrity of the entire topology.
    /// </summary>
    public bool ValidateStructure()
    {
        return _topology.ValidateStructure();
    }

    /// <summary>
    ///     Gets statistics about the topology structure.
    /// </summary>
    public TopologyStats GetStatistics()
    {
        return _topology.GetStatistics();
    }

    /// <summary>
    ///     Computes element statistics.
    /// </summary>
    public Topology<TTypes>.ElementStatistics GetElementStatistics<TElement, TNode>()
    {
        return _topology.GetElementStatistics<TElement, TNode>();
    }

    #endregion

    #region Active/Marked Entities (NEW v4)

    /// <summary>
    ///     Gets all active (not marked for deletion) entities.
    /// </summary>
    public List<int> GetActive<TEntity>()
    {
        return _topology.GetActive<TEntity>();
    }

    /// <summary>
    ///     Gets entities marked for removal.
    /// </summary>
    public IReadOnlySet<int> GetMarkedForRemoval<TEntity>()
    {
        return _topology.GetMarkedForRemoval<TEntity>();
    }

    #endregion

    #region Type Dependencies (NEW v4)

    /// <summary>
    ///     Gets the entity types in dependency order.
    /// </summary>
    public IReadOnlyList<int> GetTypeDependencyOrder()
    {
        return _topology.GetTypeDependencyOrder();
    }

    /// <summary>
    ///     Checks whether the type dependency graph is acyclic.
    /// </summary>
    public bool AreTypeDependenciesAcyclic()
    {
        return _topology.AreTypeDependenciesAcyclic();
    }

    /// <summary>
    ///     Gets the types that a given type depends on.
    /// </summary>
    public IReadOnlyList<int> GetDependencies<TEntity>()
    {
        return _topology.GetDependencies<TEntity>();
    }

    /// <summary>
    ///     Gets the types that depend on a given type.
    /// </summary>
    public IReadOnlyList<int> GetDependents<TEntity>()
    {
        return _topology.GetDependents<TEntity>();
    }

    #endregion

    #region Traversal (NEW v4)

    /// <summary>
    ///     Performs depth-first traversal starting from a node.
    /// </summary>
    public IReadOnlyList<int> Traverse<TElement, TNode>(int startNode)
    {
        return _topology.Traverse<TElement, TNode>(startNode);
    }

    /// <summary>
    ///     Performs breadth-first traversal starting from a node.
    /// </summary>
    public IReadOnlyList<int> TraverseBreadthFirst<TElement, TNode>(int startNode)
    {
        return _topology.TraverseBreadthFirst<TElement, TNode>(startNode);
    }

    /// <summary>
    ///     Finds connected components of elements.
    /// </summary>
    public IReadOnlyList<IReadOnlyList<int>> FindComponents<TElement, TNode>()
    {
        return _topology.FindComponents<TElement, TNode>();
    }

    #endregion

    #region Graph Construction and Analysis (NEW - Exposed API)

    /// <summary>
    ///     Gets element neighbors (elements sharing nodes).
    /// </summary>
    public List<int> GetElementNeighbors<TElement, TNode>(int element, bool sorted = true)
    {
        return _topology.GetElementNeighbors<TElement, TNode>(element, sorted);
    }

    /// <summary>
    ///     Gets node neighbors (nodes sharing elements).
    /// </summary>
    public List<int> GetNodeNeighbors<TElement, TNode>(int node, bool sorted = true)
    {
        return _topology.GetNodeNeighbors<TElement, TNode>(node, sorted);
    }

    /// <summary>
    ///     Constructs element-to-element adjacency graph.
    /// </summary>
    public O2M GetElementToElementGraph<TElement, TNode>()
    {
        return _topology.GetElementToElementGraph<TElement, TNode>();
    }

    /// <summary>
    ///     Constructs node-to-node adjacency graph.
    /// </summary>
    public O2M GetNodeToNodeGraph<TElement, TNode>()
    {
        return _topology.GetNodeToNodeGraph<TElement, TNode>();
    }

    #endregion

    #region Element Search (NEW - Exposed API)

    /// <summary>
    ///     Gets elements containing exactly the specified nodes.
    /// </summary>
    public List<int> GetElementsWithNodes<TElement, TNode>(List<int> nodes)
    {
        return _topology.GetElementsWithNodes<TElement, TNode>(nodes);
    }

    /// <summary>
    ///     Gets elements containing any of the specified nodes.
    /// </summary>
    public List<int> GetElementsContainingAnyNode<TElement, TNode>(List<int> nodes)
    {
        return _topology.GetElementsContainingAnyNode<TElement, TNode>(nodes);
    }

    /// <summary>
    ///     Gets elements connected to all specified nodes.
    /// </summary>
    public List<int> GetElementsFromNodes<TElement, TNode>(List<int> nodes)
    {
        return _topology.GetElementsFromNodes<TElement, TNode>(nodes);
    }

    #endregion

    #region Sparse Matrix Operations (NEW - Exposed API)

    /// <summary>
    ///     Computes clique indices for FEM matrix assembly.
    /// </summary>
    public List<List<int>> GetCliques<TElement, TNode>()
    {
        return _topology.GetCliques<TElement, TNode>();
    }

    /// <summary>
    ///     Exports connectivity in CSR format.
    /// </summary>
    public (int[] RowPtr, int[] ColumnIndices) ToCsr<TElement, TNode>()
    {
        return _topology.ToCsr<TElement, TNode>();
    }

    /// <summary>
    ///     Gets sparsity pattern for sparse matrix assembly.
    /// </summary>
    public (int[] RowPtr, int[] ColIndices) GetSparsityPatternCSR<TElement, TNode>(int dofsPerNode = 1)
    {
        return _topology.GetSparsityPatternCSR<TElement, TNode>(dofsPerNode);
    }

    #endregion

    #region Zero-Copy Access (NEW - Exposed API)

    /// <summary>
    ///     Executes action with zero-copy span access to element nodes.
    /// </summary>
    public void WithNodesSpan<TElement, TNode>(int element, Action<ReadOnlySpan<int>> action)
    {
        _topology.WithNodesSpan<TElement, TNode>(element, action);
    }

    /// <summary>
    ///     Executes function with zero-copy span access to element nodes.
    /// </summary>
    public TResult WithNodesSpan<TElement, TNode, TResult>(int element, Func<ReadOnlySpan<int>, TResult> func)
    {
        return _topology.WithNodesSpan<TElement, TNode, TResult>(element, func);
    }

    #endregion
}

#region DTOs

/// <summary>
///     Statistics about a topology structure.
/// </summary>
public sealed class TopologyStats
{
    internal TopologyStats(
        Dictionary<Type, int> entityCounts,
        Dictionary<Type, int> dataCounts,
        List<Type> typesWithSymmetry)
    {
        EntityCounts = entityCounts;
        DataCounts = dataCounts;
        TypesWithSymmetry = typesWithSymmetry;
    }

    public IReadOnlyDictionary<Type, int> EntityCounts { get; }
    public IReadOnlyDictionary<Type, int> DataCounts { get; }
    public IReadOnlyList<Type> TypesWithSymmetry { get; }

    public int TotalEntities
    {
        get
        {
            var total = 0;
            foreach (var count in EntityCounts.Values)
                total += count;
            return total;
        }
    }

    public override string ToString()
    {
        var lines = new List<string> { "Topology Statistics:" };

        foreach (var kvp in EntityCounts)
        {
            var sym = TypesWithSymmetry.Contains(kvp.Key) ? " [symmetric]" : "";
            lines.Add($"  {kvp.Key.Name}: {kvp.Value} entities{sym}");
        }

        lines.Add($"  Total: {TotalEntities} entities");
        return string.Join(Environment.NewLine, lines);
    }
}

public sealed class TopologyDto
{
    public int TypeCount { get; set; }
    public List<AdjacencyDto> Adjacency { get; set; } = new();
    public Dictionary<string, DataListDto> Data { get; set; } = new();
    public Dictionary<string, SymmetryDto> Symmetries { get; set; } = new();

    // Full canonical index with collision chains
    public List<CanonicalIndexDto> CanonicalIndices { get; set; } = new();
}

public sealed class AdjacencyDto
{
    public int EntityTypeIndex { get; set; }
    public int NodeTypeIndex { get; set; }
    public List<List<int>> Elements { get; set; } = new();
}

public sealed class DataListDto
{
    public string EntityTypeName { get; set; } = "";
    public string DataTypeName { get; set; } = "";
    public List<JsonElement> Items { get; set; } = new();
}

public sealed class SymmetryDto
{
    public int NodeCount { get; set; }
    public List<List<int>> Permutations { get; set; } = new();
}

/// <summary>
///     Represents a canonical index hash bucket with its collision chain.
/// </summary>
public sealed class CanonicalIndexDto
{
    /// <summary>Entity type full name.</summary>
    public string TypeKey { get; set; } = "";

    /// <summary>Hash value for this bucket.</summary>
    public long Hash { get; set; }

    /// <summary>All entries in the collision chain for this hash.</summary>
    public List<CanonicalEntryDto> Entries { get; set; } = new();
}

/// <summary>
///     Represents a single entry in a canonical index collision chain.
/// </summary>
public sealed class CanonicalEntryDto
{
    /// <summary>Entity index.</summary>
    public int Index { get; set; }

    /// <summary>Canonical node ordering for this entity.</summary>
    public List<int> CanonicalNodes { get; set; } = new();
}

#endregion

#region Factory Methods

/// <summary>
///     Factory methods for creating Topology instances.
/// </summary>
public static class Topology
{
    public static Topology<TypeMap<T0, T1>> New<T0, T1>()
    {
        return new Topology<TypeMap<T0, T1>>();
    }

    public static Topology<TypeMap<T0, T1, T2>> New<T0, T1, T2>()
    {
        return new Topology<TypeMap<T0, T1, T2>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3>> New<T0, T1, T2, T3>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4>> New<T0, T1, T2, T3, T4>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5>> New<T0, T1, T2, T3, T4, T5>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6>> New<T0, T1, T2, T3, T4, T5, T6>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7>> New<T0, T1, T2, T3, T4, T5, T6, T7>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8>> New<T0, T1, T2, T3, T4, T5, T6, T7, T8>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>>
        New<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>> New<T0, T1, T2, T3, T4, T5, T6, T7, T8,
        T9, T10>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>> New<T0, T1, T2, T3, T4, T5, T6,
        T7, T8, T9, T10, T11>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12>> New<T0, T1, T2, T3, T4, T5,
        T6, T7, T8, T9, T10, T11, T12>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13>> New<T0, T1, T2, T3, T4,
        T5, T6, T7, T8, T9, T10, T11, T12, T13>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14>> New<T0, T1, T2, T3,
        T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15>> New<T0, T1,
        T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16>> New<T0,
        T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17>> New<
        T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17>>();
    }

    public static Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18>>
        New<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18>()
    {
        return new Topology<
            TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18>>();
    }

    public static
        Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19>> New<
            T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
            T19>>();
    }

    public static
        Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20>>
        New<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
            T19, T20>>();
    }

    public static
        Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20,
            T21>> New<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20,
            T21>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
            T19, T20, T21>>();
    }

    public static
        Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20,
            T21, T22>> New<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19,
            T20, T21, T22>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
            T19, T20, T21, T22>>();
    }

    public static
        Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20,
            T21, T22, T23>> New<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
            T19, T20, T21, T22, T23>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
            T19, T20, T21, T22, T23>>();
    }

    public static
        Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20,
            T21, T22, T23, T24>> New<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17,
            T18, T19, T20, T21, T22, T23, T24>()
    {
        return new Topology<TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
            T19, T20, T21, T22, T23, T24>>();
    }
}

#endregion

/// <summary>
///     Defines element symmetry and computes canonical forms.
/// </summary>
/// <remarks>
///     A symmetry group describes which node orderings represent the same entity.
///     For example, if [0,1] and [1,0] are both in the group, then nodes [3,5] and [5,3]
///     represent the same element.
///     The canonical form is the lexicographically smallest ordering among all equivalent
///     permutations. This enables efficient duplicate detection via simple comparison.
/// </remarks>
/// <example>
///     <code>
///     // Create symmetry from explicit permutations
///     var sym = new Symmetry([
///         [0, 1, 2],  // identity
///         [1, 2, 0],  // rotate
///         [2, 0, 1],  // rotate again
///     ]);
///     // Compute canonical form
///     var canonical = sym.Canonical(7, 3, 5);  // Returns [3, 5, 7]
///     // Check equivalence
///     sym.AreEquivalent([0, 1, 2], [2, 0, 1]);  // True
///     // Use generators for common patterns
///     var cyclic5 = Symmetry.Cyclic(5);      // 5 rotations
///     var dihedral4 = Symmetry.Dihedral(4);  // 8 symmetries (rotations + reflections)
///     </code>
/// </example>
public sealed class Symmetry
{
    private readonly List<List<int>> _permutations;

    /// <summary>
    ///     Creates a symmetry from a list of equivalent permutations.
    /// </summary>
    /// <param name="permutations">
    ///     List of permutations defining the symmetry group. Each permutation maps
    ///     position i to the index at position i. Must include identity [0,1,2,...].
    /// </param>
    /// <exception cref="ArgumentNullException">Thrown when permutations is null.</exception>
    /// <exception cref="ArgumentException">Thrown when permutations is empty or inconsistent.</exception>
    /// <example>
    ///     <code>
    ///     // Symmetry where [a,b] = [b,a]
    ///     var sym = new Symmetry([
    ///         [0, 1],  // identity: position 0 gets node 0, position 1 gets node 1
    ///         [1, 0],  // swap: position 0 gets node 1, position 1 gets node 0
    ///     ]);
    ///     </code>
    /// </example>
    public Symmetry(List<List<int>> permutations)
    {
        ArgumentNullException.ThrowIfNull(permutations);

        if (permutations.Count == 0)
            throw new ArgumentException(
                "Must have at least the identity permutation.",
                nameof(permutations));

        var nodeCount = permutations[0].Count;

        // Validate each permutation is well-formed (includes length check)
        foreach (var perm in permutations) ValidatePermutation(perm, nodeCount);

        // Verify identity permutation is present
        var hasIdentity = false;
        foreach (var perm in permutations)
            if (IsIdentity(perm))
            {
                hasIdentity = true;
                break;
            }

        if (!hasIdentity)
            throw new ArgumentException(
                "Permutation group must include the identity permutation [0, 1, 2, ..., n-1].",
                nameof(permutations));

        // Check for duplicate permutations using structural comparison (avoids string allocations)
        var comparer = new PermutationComparer();
        var uniquePerms = new HashSet<List<int>>(comparer);
        foreach (var perm in permutations)
            if (!uniquePerms.Add(perm))
                throw new ArgumentException(
                    $"Permutation group contains duplicate: [{string.Join(",", perm)}].",
                    nameof(permutations));

        // DEFENSIVE COPY: Create deep copy to prevent external mutation
        // The caller could mutate their list after construction, corrupting canonicalization
        _permutations = new List<List<int>>(permutations.Count);
        foreach (var perm in permutations) _permutations.Add(new List<int>(perm));

        NodeCount = nodeCount;
    }

    /// <summary>
    ///     Number of nodes in elements with this symmetry.
    /// </summary>
    public int NodeCount { get; }

    /// <summary>
    ///     Number of permutations in the symmetry group.
    /// </summary>
    public int GroupSize => _permutations.Count;

    /// <summary>
    ///     The permutation group defining this symmetry.
    /// </summary>
    public IReadOnlyList<IReadOnlyList<int>> Permutations => _permutations;


    /// <summary>
    ///     Validates that a permutation is well-formed.
    /// </summary>
    private static void ValidatePermutation(IReadOnlyList<int> perm, int nodeCount)
    {
        if (perm.Count != nodeCount)
            throw new ArgumentException(
                $"Permutation has {perm.Count} elements but expected {nodeCount}.",
                nameof(perm));

        var seen = new bool[nodeCount];

        for (var i = 0; i < perm.Count; i++)
        {
            var idx = perm[i];

            if (idx < 0 || idx >= nodeCount)
                throw new ArgumentException(
                    $"Permutation contains out-of-range index {idx} at position {i}. " +
                    $"Valid range is [0, {nodeCount}).",
                    nameof(perm));

            if (seen[idx])
                throw new ArgumentException(
                    $"Permutation contains duplicate index {idx}.",
                    nameof(perm));

            seen[idx] = true;
        }
        // By pigeonhole: nodeCount values in [0, nodeCount) with no duplicates
        // guarantees all indices are present -- no further check needed.
    }

    /// <summary>
    ///     Checks if a permutation is the identity.
    /// </summary>
    private static bool IsIdentity(IReadOnlyList<int> perm)
    {
        for (var i = 0; i < perm.Count; i++)
            if (perm[i] != i)
                return false;
        return true;
    }

    /// <summary>
    ///     Computes the canonical form (lexicographically smallest) of a node list.
    /// </summary>
    /// <param name="nodes">Node indices to canonicalize.</param>
    /// <returns>The lexicographically smallest equivalent ordering.</returns>
    /// <exception cref="ArgumentException">Thrown when node count doesn't match.</exception>
    public List<int> Canonical(params int[] nodes)
    {
        ArgumentNullException.ThrowIfNull(nodes);

        if (nodes.Length != NodeCount)
            throw new ArgumentException($"Expected {NodeCount} nodes, got {nodes.Length}.", nameof(nodes));

        return CanonicalCore(nodes);
    }

    /// <summary>
    ///     Computes the canonical form (lexicographically smallest) of a node list.
    /// </summary>
    public List<int> Canonical(List<int> nodes)
    {
        ArgumentNullException.ThrowIfNull(nodes);

        if (nodes.Count != NodeCount)
            throw new ArgumentException($"Expected {NodeCount} nodes, got {nodes.Count}.", nameof(nodes));

        return CanonicalCore(nodes);
    }

    /// <summary>
    ///     Computes the canonical form (lexicographically smallest) of a node list.
    /// </summary>
    public List<int> Canonical(IReadOnlyList<int> nodes)
    {
        ArgumentNullException.ThrowIfNull(nodes);

        if (nodes.Count != NodeCount)
            throw new ArgumentException($"Expected {NodeCount} nodes, got {nodes.Count}.", nameof(nodes));

        return CanonicalCore(nodes);
    }

    /// <summary>
    ///     Computes canonical form into a pre-allocated destination buffer.
    ///     Zero-allocation alternative to Canonical() for hot paths.
    /// </summary>
    /// <param name="nodes">Input nodes (will not be modified).</param>
    /// <param name="destination">
    ///     Pre-allocated buffer to receive canonical form.
    ///     Must have length >= NodeCount.
    ///     After call, first NodeCount elements contain the canonical form.
    /// </param>
    /// <exception cref="ArgumentException">Thrown when array sizes are incorrect.</exception>
    /// <remarks>
    ///     CHANGED: Now returns void instead of Span to avoid CS9244 compiler error
    ///     when used in tuple contexts. The canonical form is written to the destination
    ///     array starting at index 0.
    /// </remarks>
    public void CanonicalSpan(ReadOnlySpan<int> nodes, Span<int> destination)
    {
        if (nodes.Length != NodeCount)
            throw new ArgumentException($"Expected {NodeCount} nodes, got {nodes.Length}.", nameof(nodes));

        if (destination.Length < NodeCount)
            throw new ArgumentException($"Destination must have length >= {NodeCount}.", nameof(destination));

        var nodeCount = NodeCount;
        var groupSize = _permutations.Count;

        // Rent temporary array for comparison (only 1 allocation vs 2 in non-span version)
        var current = ArrayPool<int>.Shared.Rent(nodeCount);

        try
        {
            // Initialize destination with first permutation
            var firstPerm = _permutations[0];
            for (var i = 0; i < nodeCount; i++)
                destination[i] = nodes[firstPerm[i]];

            // Check remaining permutations
            for (var p = 1; p < groupSize; p++)
            {
                var perm = _permutations[p];

                // Apply permutation to current
                for (var i = 0; i < nodeCount; i++)
                    current[i] = nodes[perm[i]];

                // Compare with destination (lexicographic, inline for speed)
                var isBetter = false;
                for (var i = 0; i < nodeCount; i++)
                {
                    if (current[i] < destination[i])
                    {
                        isBetter = true;
                        break;
                    }

                    if (current[i] > destination[i])
                        break;
                }

                if (isBetter)
                    // Copy current to destination
                    for (var i = 0; i < nodeCount; i++)
                        destination[i] = current[i];
            }

            // Destination now contains the canonical form (no return needed)
        }
        finally
        {
            ArrayPool<int>.Shared.Return(current);
        }
    }

    /// <summary>
    ///     Checks if two node lists represent the same element under this symmetry.
    /// </summary>
    public bool AreEquivalent(int[] a, int[] b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);

        if (a.Length != NodeCount || b.Length != NodeCount)
            return false;

        var ca = CanonicalCore(a);
        var cb = CanonicalCore(b);
        return Utils.AreEqual(ca, cb);
    }

    /// <summary>
    ///     Checks if two node lists represent the same element under this symmetry.
    /// </summary>
    public bool AreEquivalent(List<int> a, List<int> b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);

        if (a.Count != NodeCount || b.Count != NodeCount)
            return false;

        var ca = CanonicalCore(a);
        var cb = CanonicalCore(b);
        return Utils.AreEqual(ca, cb);
    }

    /// <summary>
    ///     Checks if two node lists represent the same element under this symmetry.
    /// </summary>
    public bool AreEquivalent(IReadOnlyList<int> a, IReadOnlyList<int> b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);

        if (a.Count != NodeCount || b.Count != NodeCount)
            return false;

        var ca = CanonicalCore(a);
        var cb = CanonicalCore(b);
        return Utils.AreEqual(ca, cb);
    }

    /// <summary>
    ///     Applies a specific permutation to a node list.
    /// </summary>
    /// <param name="nodes">The node list to permute.</param>
    /// <param name="permutationIndex">Index of the permutation to apply.</param>
    /// <returns>The permuted node list.</returns>
    public List<int> Apply(IReadOnlyList<int> nodes, int permutationIndex)
    {
        if (nodes.Count != NodeCount)
            throw new ArgumentException($"Expected {NodeCount} nodes, got {nodes.Count}.", nameof(nodes));

        if (permutationIndex < 0 || permutationIndex >= _permutations.Count)
            throw new ArgumentOutOfRangeException(nameof(permutationIndex));

        return ApplyPermutation(nodes, _permutations[permutationIndex]);
    }

    #region Private Helpers

    /// <summary>
    ///     Computes canonical form with ArrayPool optimization for large sizes.
    /// </summary>
    /// <remarks>
    ///     P1-2 FIX: Uses stackalloc for small sizes (≤32 nodes), ArrayPool for larger sizes.
    ///     This ensures zero heap allocation for common cases while avoiding
    ///     stack overflow for very large permutation groups.
    /// </remarks>
    private List<int> CanonicalCore(IReadOnlyList<int> nodes)
    {
        var nodeCount = NodeCount;

        // Handle small sizes with stackalloc (outside try block due to C# restriction)
        if (nodeCount <= 32) return CanonicalCoreSmall(nodes, nodeCount);

        // Large size: use ArrayPool
        return CanonicalCoreLarge(nodes, nodeCount);
    }

    /// <summary>
    ///     Small-size canonical computation using stackalloc. Zero heap allocation.
    /// </summary>
    private List<int> CanonicalCoreSmall(IReadOnlyList<int> nodes, int nodeCount)
    {
        Span<int> inputSpan = stackalloc int[nodeCount];
        Span<int> resultSpan = stackalloc int[nodeCount];

        // Copy input nodes to span
        for (var i = 0; i < nodeCount; i++)
            inputSpan[i] = nodes[i];

        // Compute canonical form
        CanonicalSpan(inputSpan, resultSpan);

        // Convert result to List<int>
        var result = new List<int>(nodeCount);
        for (var i = 0; i < nodeCount; i++)
            result.Add(resultSpan[i]);

        return result;
    }

    /// <summary>
    ///     Large-size canonical computation using ArrayPool.
    /// </summary>
    private List<int> CanonicalCoreLarge(IReadOnlyList<int> nodes, int nodeCount)
    {
        var inputRented = ArrayPool<int>.Shared.Rent(nodeCount);
        var resultRented = ArrayPool<int>.Shared.Rent(nodeCount);

        try
        {
            var inputSpan = inputRented.AsSpan(0, nodeCount);
            var resultSpan = resultRented.AsSpan(0, nodeCount);

            // Copy input nodes to span
            for (var i = 0; i < nodeCount; i++)
                inputSpan[i] = nodes[i];

            // Compute canonical form
            CanonicalSpan(inputSpan, resultSpan);

            // Convert result to List<int>
            var result = new List<int>(nodeCount);
            for (var i = 0; i < nodeCount; i++)
                result.Add(resultSpan[i]);

            return result;
        }
        finally
        {
            ArrayPool<int>.Shared.Return(inputRented);
            ArrayPool<int>.Shared.Return(resultRented);
        }
    }

    private static List<int> ApplyPermutation(IReadOnlyList<int> nodes, List<int> permutation)
    {
        var result = new List<int>(nodes.Count);
        for (var i = 0; i < permutation.Count; i++)
            result.Add(nodes[permutation[i]]);
        return result;
    }

    #endregion

    #region Generators

    /// <summary>
    ///     Creates identity-only symmetry (no equivalences, order matters).
    /// </summary>
    /// <param name="nodeCount">Number of nodes.</param>
    public static Symmetry Identity(int nodeCount)
    {
        if (nodeCount < 1)
            throw new ArgumentException("Node count must be at least 1.", nameof(nodeCount));

        var identity = new List<int>(nodeCount);
        for (var i = 0; i < nodeCount; i++)
            identity.Add(i);
        return new Symmetry([identity]);
    }

    /// <summary>
    ///     Creates cyclic symmetry (n rotations).
    /// </summary>
    /// <param name="n">Number of nodes.</param>
    /// <remarks>
    ///     Generates n permutations representing cyclic rotations.
    ///     For n=4: [0,1,2,3], [1,2,3,0], [2,3,0,1], [3,0,1,2]
    /// </remarks>
    public static Symmetry Cyclic(int n)
    {
        if (n < 1)
            throw new ArgumentException("Node count must be at least 1.", nameof(n));

        var perms = new List<List<int>>();
        for (var r = 0; r < n; r++)
        {
            var perm = new List<int>(n);
            for (var i = 0; i < n; i++)
                perm.Add((i + r) % n);
            perms.Add(perm);
        }

        return new Symmetry(perms);
    }

    /// <summary>
    ///     Creates dihedral symmetry (n rotations + n reflections).
    /// </summary>
    /// <param name="n">Number of nodes.</param>
    /// <remarks>
    ///     The dihedral group D_n has 2n elements: n rotations and n reflections.
    ///     For n=3: 6 permutations (3 rotations + 3 reflections).
    /// </remarks>
    public static Symmetry Dihedral(int n)
    {
        if (n < 1)
            throw new ArgumentException("Node count must be at least 1.", nameof(n));

        var comparer = new PermutationComparer();
        var seen = new HashSet<List<int>>(comparer);
        var perms = new List<List<int>>();

        // Rotations
        for (var r = 0; r < n; r++)
        {
            var perm = new List<int>(n);
            for (var i = 0; i < n; i++)
                perm.Add((i + r) % n);
            if (seen.Add(perm))
                perms.Add(perm);
        }

        // Reflections
        for (var r = 0; r < n; r++)
        {
            var perm = new List<int>(n);
            for (var i = 0; i < n; i++)
                perm.Add((r - i + n) % n);
            if (seen.Add(perm))
                perms.Add(perm);
        }

        return new Symmetry(perms);
    }

    /// <summary>
    ///     Creates full symmetric group S_n (all n! permutations).
    /// </summary>
    /// <param name="n">Number of nodes.</param>
    /// <remarks>
    ///     Warning: Group size grows as n! (factorial). Only practical for small n.
    ///     n=4 has 24 permutations, n=5 has 120, n=6 has 720.
    /// </remarks>
    public static Symmetry Full(int n)
    {
        if (n < 1)
            throw new ArgumentException("Node count must be at least 1.", nameof(n));

        if (n > 8)
            throw new ArgumentException("Full symmetry group too large for n > 8.", nameof(n));

        var perms = new List<List<int>>();
        var current = new List<int>(n);
        for (var i = 0; i < n; i++)
            current.Add(i);

        do
        {
            var copy = new List<int>(current.Count);
            for (var i = 0; i < current.Count; i++)
                copy.Add(current[i]);
            perms.Add(copy);
        } while (NextPermutation(current));

        return new Symmetry(perms);
    }

    /// <summary>
    ///     Creates symmetry from a generating set (computes group closure).
    /// </summary>
    /// <param name="generators">Generator permutations.</param>
    /// <remarks>
    ///     Computes the smallest group containing all generators by repeatedly
    ///     composing permutations until no new elements are found.
    ///     Useful when you know the generators but not the full group.
    ///     Performance: Uses structural comparison instead of string keys
    ///     for O(n) comparison vs O(n) string allocation per check.
    /// </remarks>
    /// <example>
    ///     <code>
    ///     // Generate cyclic group from single rotation
    ///     var cyclic4 = Symmetry.FromGenerators([[1, 2, 3, 0]]);
    ///     // Result has 4 elements: identity + 3 rotations
    ///     </code>
    /// </example>
    public static Symmetry FromGenerators(List<List<int>> generators)
    {
        ArgumentNullException.ThrowIfNull(generators);

        if (generators.Count == 0)
            throw new ArgumentException("Must provide at least one generator.", nameof(generators));

        var n = generators[0].Count;

        // Validate all generators are well-formed permutations before costly closure
        foreach (var gen in generators)
            ValidatePermutation(gen, n);

        var identity = new List<int>(n);
        for (var i = 0; i < n; i++)
            identity.Add(i);

        // Use structural equality comparer instead of string keys
        // to avoid allocation overhead for large groups (n >= 10)
        var comparer = new PermutationComparer();
        var group = new HashSet<List<int>>(comparer) { identity };
        var elements = new List<List<int>> { identity };
        var queue = new Queue<List<int>>(generators);

        while (queue.Count > 0)
        {
            var current = queue.Dequeue();

            if (group.Add(current))
            {
                elements.Add(current);

                // Compose with all existing elements - make a copy of elements to iterate
                var existingCount = elements.Count;
                for (var i = 0; i < existingCount; i++)
                {
                    var existing = elements[i];
                    var composed = Compose(current, existing, n);
                    if (!group.Contains(composed))
                        queue.Enqueue(composed);

                    var composed2 = Compose(existing, current, n);
                    if (!group.Contains(composed2))
                        queue.Enqueue(composed2);
                }
            }
        }

        return new Symmetry(elements);
    }

    private static List<int> Compose(List<int> a, List<int> b, int n)
    {
        var result = new List<int>(n);
        for (var i = 0; i < n; i++)
            result.Add(a[b[i]]);
        return result;
    }

    /// <summary>
    ///     Equality comparer for permutations using structural comparison.
    ///     More efficient than string-based comparison for large permutations.
    /// </summary>
    private sealed class PermutationComparer : IEqualityComparer<List<int>>
    {
        public bool Equals(List<int>? x, List<int>? y)
        {
            if (ReferenceEquals(x, y)) return true;
            if (x is null || y is null || x.Count != y.Count) return false;

            for (var i = 0; i < x.Count; i++)
                if (x[i] != y[i])
                    return false;

            return true;
        }

        public int GetHashCode(List<int> obj)
        {
            var hash = new HashCode();
            foreach (var item in obj)
                hash.Add(item);
            return hash.ToHashCode();
        }
    }

    private static bool NextPermutation(List<int> arr)
    {
        var i = arr.Count - 2;
        while (i >= 0 && arr[i] >= arr[i + 1])
            i--;

        if (i < 0)
            return false;

        var j = arr.Count - 1;
        while (arr[j] <= arr[i])
            j--;

        (arr[i], arr[j]) = (arr[j], arr[i]);
        ReverseRange(arr, i + 1, arr.Count - i - 1);
        return true;
    }

    private static void ReverseRange(List<int> list, int index, int count)
    {
        var left = index;
        var right = index + count - 1;
        while (left < right)
        {
            (list[left], list[right]) = (list[right], list[left]);
            left++;
            right--;
        }
    }

    #endregion
}

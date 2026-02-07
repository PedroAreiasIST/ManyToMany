using System.Diagnostics;
using System.Runtime.CompilerServices;

namespace Numerical;

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
    ///     Tracks whether element location cache has been computed.
    ///     <b>INVARIANT:</b> Must be updated atomically with <see cref="_nodeLocComputed" />.
    /// </summary>
    private bool _elemLocComputed;

    /// <summary>
    ///     Tracks whether node location cache has been computed.
    ///     <b>INVARIANT:</b> Must be updated atomically with <see cref="_elemLocComputed" />.
    /// </summary>
    private bool _nodeLocComputed;

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
        _elemLocComputed = false;
        _nodeLocComputed = false;
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
        _elemLocComputed = false;
        _nodeLocComputed = false;
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
        _elemLocComputed = false;
        _nodeLocComputed = false;
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
        _elemLocComputed = false;
        _nodeLocComputed = false;
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
        _elemLocComputed = false;
        _nodeLocComputed = false;
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
            _elemLocComputed = false;
            _nodeLocComputed = false;
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
        if (_isInSync && _elemLocComputed && _nodeLocComputed)
            return; // Read lock held, all caches valid
        _rwLock.ExitReadLock();

        // Slow path: compute under write lock, then downgrade
        _rwLock.EnterWriteLock();
        try
        {
            if (!_isInSync && _batchNesting == 0)
                SynchronizeTranspose();
            if (!_elemLocComputed || !_nodeLocComputed)
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
    ///     Ensures position caches are computed (requires transpose to be synchronized first).
    /// </summary>
    private void EnsurePositionCachesComputed()
    {
        EnsureSynchronized(); // First ensure transpose is ready

        // Fast path: check if already computed
        _rwLock.EnterReadLock();
        try
        {
            if (_elemLocComputed && _nodeLocComputed)
                return;
        }
        finally
        {
            _rwLock.ExitReadLock();
        }

        // Slow path: compute position caches
        _rwLock.EnterWriteLock();
        try
        {
            // Double-check
            if (_elemLocComputed && _nodeLocComputed)
                return;

            ComputePositionCaches();
        }
        finally
        {
            _rwLock.ExitWriteLock();
        }
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
        _elemLocComputed = false;
        _nodeLocComputed = false;
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

        _elemLocComputed = true;
        _nodeLocComputed = true;
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
        _elemLocComputed = false;
        _nodeLocComputed = false;

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
        ThrowIfDisposed();
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
            if (_elemLocComputed && _nodeLocComputed)
            {
                cloned._elemeloc = _elemeloc; // ReadOnlyCollection, safe to share
                cloned._nodeloc = _nodeloc; // ReadOnlyCollection, safe to share
                cloned._elemLocComputed = true;
                cloned._nodeLocComputed = true;
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
        // If disposal was previously incomplete, allow retry
        if (_disposalIncomplete)
        {
            _disposalIncomplete = false;
            // Fall through to retry disposal
        }
        else
        {
            // Atomic check-and-set to ensure disposal happens only once
            if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0)
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
                    // The object remains marked as disposed to prevent further use
                    _disposalIncomplete = true;
                    
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

            // Lock was acquired - safe to dispose
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
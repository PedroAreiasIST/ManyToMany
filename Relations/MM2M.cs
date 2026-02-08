using System.Diagnostics;
using System.Runtime.CompilerServices;

namespace Numerical;

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

        _listOfMarked = new Dictionary<int, HashSet<int>>(numberOfTypes);
        for (var i = 0; i < numberOfTypes; i++)
            _listOfMarked[i] = new HashSet<int>();
    }

    #endregion

    #region Fields

    private readonly Dictionary<int, HashSet<int>> _listOfMarked;
    private readonly M2M[,] _mat;
    private readonly ReaderWriterLockSlim _rwLock = new(LockRecursionPolicy.SupportsRecursion);

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
                var result = new Dictionary<int, IReadOnlySet<int>>(_listOfMarked.Count);
                foreach (var kvp in _listOfMarked)
                    result[kvp.Key] = new HashSet<int>(kvp.Value);
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
        var maxNode = 0;
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

        return maxNode;
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
            foreach (var kvp in _listOfMarked)
                if (kvp.Value.Count > 0)
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

            // PHASE 4: Atomic swap - all clones succeeded, now swap into place
            // Dispose old matrices and replace with new ones
            for (var i = 0; i < NumberOfTypes; i++)
            for (var j = 0; j < NumberOfTypes; j++)
            {
                var old = _mat[i, j];
                _mat[i, j] = clonedMat[i, j];
                old.Dispose();
            }

            // Increment version to signal stale reference invalidation
            // Any consumers with cached M2M references should check Version before use
            Interlocked.Increment(ref _version);

            // PHASE 5: Cleanup (only after successful compression)
            // Clear all marked lists
            foreach (var kvp in _listOfMarked)
                kvp.Value.Clear();

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
            return AreTypesAcyclic();
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

            return !hasDeps || typeDeps.IsAcyclic();
        }
        finally
        {
            _rwLock.ExitReadLock();
        }
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
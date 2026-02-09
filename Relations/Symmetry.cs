using System.Buffers;

namespace Numerical;

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

        // Check for duplicate permutations
        var uniquePerms = new HashSet<string>();
        foreach (var perm in permutations)
        {
            var key = string.Join(",", perm);
            if (!uniquePerms.Add(key))
                throw new ArgumentException(
                    $"Permutation group contains duplicate: [{key}].",
                    nameof(permutations));
        }

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
        if (nodes.Length != NodeCount)
            throw new ArgumentException($"Expected {NodeCount} nodes, got {nodes.Length}.", nameof(nodes));

        return CanonicalCore(nodes);
    }

    /// <summary>
    ///     Computes the canonical form (lexicographically smallest) of a node list.
    /// </summary>
    public List<int> Canonical(List<int> nodes)
    {
        if (nodes.Count != NodeCount)
            throw new ArgumentException($"Expected {NodeCount} nodes, got {nodes.Count}.", nameof(nodes));

        return CanonicalCore(nodes);
    }

    /// <summary>
    ///     Computes the canonical form (lexicographically smallest) of a node list.
    /// </summary>
    public List<int> Canonical(IReadOnlyList<int> nodes)
    {
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

        var seen = new HashSet<string>();
        var perms = new List<List<int>>();

        // Rotations
        for (var r = 0; r < n; r++)
        {
            var perm = new List<int>(n);
            for (var i = 0; i < n; i++)
                perm.Add((i + r) % n);
            var key = string.Join(",", perm);
            if (seen.Add(key))
                perms.Add(perm);
        }

        // Reflections
        for (var r = 0; r < n; r++)
        {
            var perm = new List<int>(n);
            for (var i = 0; i < n; i++)
                perm.Add((r - i + n) % n);
            var key = string.Join(",", perm);
            if (seen.Add(key))
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
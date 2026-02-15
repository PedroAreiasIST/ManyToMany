using System.Buffers;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Numerics;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Text;
using System.Text.Json.Serialization;
using System.Text.RegularExpressions;
using Microsoft.Extensions.ObjectPool;
using static System.Runtime.CompilerServices.MethodImplOptions;
using VectorT = System.Numerics.Vector;

namespace Numerical;

// Helper class for object pooling
file class HashSetPoolPolicy : IPooledObjectPolicy<HashSet<int>>
{
    public HashSet<int> Create()
    {
        return new HashSet<int>(256);
    }

    public bool Return(HashSet<int> obj)
    {
        if (obj.Count > 10000) return false;
        obj.Clear();
        return true;
    }
}

/// <summary>
///     High-performance Compressed Sparse Row (CSR) matrix with complete functionality.
///     FINAL PATCHED VERSION: All native backends (PARDISO, CUDA)
///     now use dynamic NativeLibrary resolution for maximum portability.
///     NOTE: MKL Sparse BLAS has been removed; MKL is only used for PARDISO solver.
///     REVISION 2.0 - ADDITIONAL AUDIT FIXES:
///     - Improved ArrayPool usage in MultiplyTransposedParallel for large matrices
///     - Added shared static ParallelOptions to reduce allocations
///     - Enhanced error messages with actionable guidance throughout
///     - Made threshold constants public for user tuning
///     - Added null validation in Diagonal() factory method
///     - Improved thread-safety documentation for Equals method
///     - Better validation for skipValidation matrices in SIMD operations
///     PRIOR AUDIT FIXES:
///     - Fixed lock ordering deadlock risk in Equals method
///     - Fixed integer overflow check placement in MultiplySymbolic
///     - Added bounds checking in SIMD operations
///     - Fixed GPU disposal race condition
///     - Improved modification counter semantics
///     - Enhanced error messages throughout
///     - Added duplicate detection in validation
///     - Improved BiCGSTAB breakdown detection
///     - Added GPU initialization check
///     - Enhanced parallel tuning options
/// </summary>
public sealed class CSR : IFormattable, IEquatable<CSR>, ICloneable, IDisposable
{
    // Constants for tolerances and thresholds - public for user tuning
    /// <summary>Default tolerance for near-zero comparisons.</summary>
    public const double DEFAULT_TOLERANCE = 1e-14;

    /// <summary>Minimum rows before parallel processing is beneficial.</summary>
    public const int MIN_ROWS_FOR_PARALLEL = 1000;

    /// <summary>Minimum rows before SIMD processing is beneficial.</summary>
    public const int MIN_ROWS_FOR_SIMD = 5000;

    /// <summary>Minimum rows for GPU acceleration consideration.</summary>
    public const int MIN_ROWS_FOR_GPU = 50000;

    /// <summary>Minimum non-zeros for GPU acceleration consideration.</summary>
    public const int MIN_NNZ_FOR_GPU = 1000000;

    private const int HASHSET_POOL_MAX_SIZE = 10000;
    private const int SMALL_ARRAY_MERGE_THRESHOLD = 10000;

    // Use ParallelConfig.Options for all parallel operations (centralized control)

    // Object pooling for temporary allocations
    private static readonly ObjectPool<HashSet<int>> hashSetPool =
        new DefaultObjectPool<HashSet<int>>(new HashSetPoolPolicy());

    [JsonInclude] private readonly int[] columnIndices;

    // Flag to track if matrix was constructed with skipValidation
    private readonly bool constructedWithSkipValidation;

    private readonly object disposeLock = new();
    [JsonInclude] private readonly int ncols;
    [JsonInclude] private readonly int nrows;

    [JsonInclude] private readonly int[] rowPointers;

    // Values array - the only mutable part, protected by lock
    private readonly object syncLock = new();

    // Caching using Lazy<T>
    private Lazy<CSR>? cachedTranspose;
    private volatile int cachedTransposeModCount;

    // GPU acceleration
    private CuSparseBackend? gpuAccelerator;
    private volatile bool isDisposed;
    private volatile bool isGpuInitialized; // FIXED: Issue #14 - Track GPU state

    // Modification tracking - FIXED: Issue #6 - Removed unchecked for better semantics
    private volatile int modificationCount;
    private double[]? values;

    #region Constructors

    public CSR(int[] rowPointers, int[] columnIndices, int nCols, double[]? values = null, bool skipValidation = false)
    {
        this.rowPointers = rowPointers ?? throw new ArgumentNullException(nameof(rowPointers));
        this.columnIndices = columnIndices ?? throw new ArgumentNullException(nameof(columnIndices));

        if (nCols < 0) throw new ArgumentException("Number of columns must be non-negative", nameof(nCols));

        ncols = nCols;
        nrows = this.rowPointers.Length - 1;

        if (nrows < 0)
            throw new ArgumentException("Row pointers array must have at least 1 element", nameof(rowPointers));

        constructedWithSkipValidation = skipValidation;
        if (!skipValidation) ValidateCSRStructure(this.rowPointers, this.columnIndices, nrows, ncols);

        // FIXED: Issue #20 - Better error messages
        if (values != null && values.Length != this.columnIndices.Length)
            throw new ArgumentException(
                $"Values array length {values.Length} must match column indices length {this.columnIndices.Length}. " +
                $"Expected {this.columnIndices.Length} elements.",
                nameof(values));

        this.values = values;
        modificationCount = 0;
    }

    // FIXED: Issue #16 - Added duplicate detection in validation
    // NOTE: Columns do NOT need to be sorted - we use marker-based lookup instead of binary search
    private static void ValidateCSRStructure(int[] rowPtrs, int[] colIndices, int nrows, int ncols)
    {
        for (var i = 0; i < nrows; i++)
        {
            if (rowPtrs[i] < 0)
                throw new ArgumentException($"Row pointer at position {i} is negative: {rowPtrs[i]}", nameof(rowPtrs));
            if (rowPtrs[i] > rowPtrs[i + 1])
                throw new ArgumentException(
                    $"Row pointers must be monotonically increasing. Found rowPtrs[{i}]={rowPtrs[i]} > rowPtrs[{i + 1}]={rowPtrs[i + 1]}",
                    nameof(rowPtrs));
        }

        if (rowPtrs[nrows] != colIndices.Length)
            throw new ArgumentException(
                $"Last row pointer ({rowPtrs[nrows]}) must equal column indices length ({colIndices.Length})",
                nameof(rowPtrs));

        // Validate column indices and check for duplicates using marker array
        // This allows unsorted columns while still detecting duplicates
        var marker = ArrayPool<int>.Shared.Rent(ncols);
        try
        {
            Array.Fill(marker, -1, 0, ncols);

            for (var i = 0; i < nrows; i++)
            {
                var rowStart = rowPtrs[i];
                var rowEnd = rowPtrs[i + 1];

                for (var j = rowStart; j < rowEnd; j++)
                {
                    var col = colIndices[j];
                    if (col < 0 || col >= ncols)
                        throw new ArgumentException(
                            $"Column index {col} at position {j} (row {i}) is out of range [0, {ncols})",
                            nameof(colIndices));

                    // Check for duplicates using marker
                    if (marker[col] == i)
                        throw new ArgumentException(
                            $"Duplicate column index {col} found in row {i}. " +
                            $"Each column may appear at most once per row.",
                            nameof(colIndices));

                    marker[col] = i;
                }
            }
        }
        finally
        {
            ArrayPool<int>.Shared.Return(marker);
        }
    }

    // FIXED: Issue #19 - Changed to warning (false) instead of compile error (true)
    [Obsolete(
        "Use constructor with explicit column count to avoid dimension errors. This will become an error in the next major version.",
        false)]
    public CSR(int[] rowPointers, int[] columnIndices, double[]? values = null)
        : this(rowPointers, columnIndices, InferColumnCount(columnIndices), values)
    {
    }

    /// <summary>
    ///     Creates a CSR matrix from a list of rows, where each row is a list of column indices.
    /// </summary>
    /// <param name="rows">List of rows, each containing column indices of non-zero entries.</param>
    /// <param name="sorted">Ignored. Kept for backward compatibility. Column indices do not need to be sorted.</param>
    /// <param name="enableGpu">If true, attempts to initialize GPU acceleration.</param>
    public CSR(List<List<int>> rows, bool sorted = false, bool enableGpu = false)
    {
        ArgumentNullException.ThrowIfNull(rows);

        nrows = rows.Count;
        ncols = InferColumnCount(rows);

        rowPointers = new int[nrows + 1];
        var nnz = 0;

        for (var i = 0; i < nrows; i++)
        {
            rowPointers[i] = nnz;
            nnz += rows[i].Count;
        }

        rowPointers[nrows] = nnz;
        columnIndices = new int[nnz];

        var pos = 0;
        for (var i = 0; i < nrows; i++)
        {
            var rowSpan = CollectionsMarshal.AsSpan(rows[i]);
            for (var j = 0; j < rowSpan.Length; j++)
                columnIndices[pos++] = rowSpan[j];
        }

        ValidateCSRStructure(rowPointers, columnIndices, nrows, ncols);
        constructedWithSkipValidation = false;

        if (enableGpu)
            try
            {
                InitializeGpu();
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"GPU initialization failed: {ex.Message}. Continuing with CPU only.");
                gpuAccelerator = null;
                isGpuInitialized = false;
            }
    }

    private static int InferColumnCount(int[] columnIndices)
    {
        ArgumentNullException.ThrowIfNull(columnIndices);

        if (columnIndices.Length == 0)
            return 0;

        var maxColumn = columnIndices.Max();
        if (maxColumn == int.MaxValue)
            throw new ArgumentException(
                "Column index contains Int32.MaxValue. Column count cannot be inferred safely; use constructor with explicit nCols.",
                nameof(columnIndices));

        return maxColumn + 1;
    }

    private static int InferColumnCount(List<List<int>> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);

        var maxColumn = -1;
        for (var i = 0; i < rows.Count; i++)
        {
            var row = rows[i] ?? throw new ArgumentException(
                $"Row {i} is null. Rows collection must not contain null entries.",
                nameof(rows));

            var rowSpan = CollectionsMarshal.AsSpan(row);
            for (var j = 0; j < rowSpan.Length; j++)
                if (rowSpan[j] > maxColumn)
                    maxColumn = rowSpan[j];
        }

        if (maxColumn == int.MaxValue)
            throw new ArgumentException(
                "Column index contains Int32.MaxValue. Column count cannot be inferred safely; use constructor with explicit nCols.",
                nameof(rows));

        return maxColumn + 1;
    }

    #endregion

    #region Properties and Accessors

    public int Rows
    {
        get
        {
            ThrowIfDisposed();
            return nrows;
        }
    }

    public int Columns
    {
        get
        {
            ThrowIfDisposed();
            return ncols;
        }
    }

    public int NonZeroCount
    {
        get
        {
            ThrowIfDisposed();
            return columnIndices.Length;
        }
    }

    /// <summary>Density of the matrix: fraction of non-zero entries (nnz / (rows * cols)).</summary>
    public double Sparsity
    {
        get
        {
            ThrowIfDisposed();
            return nrows > 0 && ncols > 0 ? NonZeroCount / (nrows * (double)ncols) : 0.0;
        }
    }

    public ReadOnlySpan<int> RowPointers
    {
        get
        {
            ThrowIfDisposed();
            return rowPointers;
        }
    }

    public ReadOnlySpan<int> ColumnIndices
    {
        get
        {
            ThrowIfDisposed();
            return columnIndices;
        }
    }

    internal int[] RowPointersArray
    {
        get
        {
            ThrowIfDisposed();
            return rowPointers;
        }
    }

    internal int[] ColumnIndicesArray
    {
        get
        {
            ThrowIfDisposed();
            return columnIndices;
        }
    }

    // FIXED: Issue #21 - Better validation in Values property setter
    public double[]? Values
    {
        get
        {
            ThrowIfDisposed();
            lock (syncLock)
            {
                return values;
            }
        }
        set
        {
            ThrowIfDisposed();

            // FIXED: Issue #20 - Improved error message
            if (value != null && value.Length != columnIndices.Length)
                throw new ArgumentException(
                    $"Values array length {value.Length} must match non-zero count {NonZeroCount}. " +
                    $"Expected exactly {NonZeroCount} elements.",
                    nameof(value));

            // FIXED: Issue #7 - Keep lock through GPU update to prevent race condition
            lock (syncLock)
            {
                values = value;
                cachedTranspose = null;
                // Use unchecked to wrap around on overflow - this is intentional
                // for long-running systems that may exceed int.MaxValue modifications
                unchecked
                {
                    Interlocked.Increment(ref modificationCount);
                }

                // GPU update inside lock to prevent disposal race
                if (isGpuInitialized && gpuAccelerator != null && value != null)
                    try
                    {
                        gpuAccelerator.UpdateValues(value);
                    }
                    catch (Exception ex)
                    {
                        Debug.WriteLine($"GPU update failed: {ex.Message}");
                        // Don't throw - allow CPU fallback
                    }
            }
        }
    }

    public bool HasValues
    {
        get
        {
            ThrowIfDisposed();
            lock (syncLock)
            {
                return values != null;
            }
        }
    }

    #endregion

    #region Helper Methods

    [MethodImpl(AggressiveInlining)]
    private double[] GetValuesOrThrow()
    {
        lock (syncLock)
        {
            return values ?? throw new InvalidOperationException(
                "Matrix must have values for this operation. Call the Values property setter to assign values first.");
        }
    }

    [MethodImpl(AggressiveInlining)]
    internal double[] GetValuesInternal()
    {
        lock (syncLock)
        {
            return values ?? throw new InvalidOperationException("Matrix must have values");
        }
    }

    [MethodImpl(AggressiveInlining)]
    internal void ThrowIfDisposed([CallerMemberName] string? caller = null)
    {
        if (isDisposed)
            throw new ObjectDisposedException(nameof(CSR), $"Cannot call {caller} on a disposed CSR matrix");
    }

    #endregion

    #region Matrix-Vector Operations

    public double[] Multiply(double[] vector)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better null and length checks
        if (vector == null)
            throw new ArgumentNullException(nameof(vector), "Input vector cannot be null");
        if (vector.Length != ncols)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix columns {ncols}",
                nameof(vector));

        var localValues = GetValuesOrThrow();
        var result = new double[nrows];

        for (var i = 0; i < nrows; i++)
        {
            double sum = 0;
            for (var j = rowPointers[i]; j < rowPointers[i + 1]; j++)
                sum += localValues[j] * vector[columnIndices[j]];
            result[i] = sum;
        }

        return result;
    }

    /// <summary>
    ///     Multiply into a pre-allocated result span, avoiding allocation.
    ///     result = A * vector
    /// </summary>
    public void Multiply(ReadOnlySpan<double> vector, Span<double> result)
    {
        ThrowIfDisposed();

        if (vector.Length != ncols)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix columns {ncols}",
                nameof(vector));
        if (result.Length < nrows)
            throw new ArgumentException(
                $"Result span length {result.Length} must be at least {nrows}",
                nameof(result));

        var localValues = GetValuesOrThrow();

        for (var i = 0; i < nrows; i++)
        {
            double sum = 0;
            for (var j = rowPointers[i]; j < rowPointers[i + 1]; j++)
                sum += localValues[j] * vector[columnIndices[j]];
            result[i] = sum;
        }
    }

    public double[] MultiplyParallel(double[] vector)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better error messages
        if (vector == null)
            throw new ArgumentNullException(nameof(vector), "Input vector cannot be null");
        if (vector.Length != ncols)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix columns {ncols}",
                nameof(vector));

        var localValues = GetValuesOrThrow();
        var result = new double[nrows];

        if (nrows < MIN_ROWS_FOR_PARALLEL) return Multiply(vector);

        // Use shared ParallelOptions to avoid allocation
        Parallel.For(0, nrows, ParallelConfig.Options, i =>
        {
            double sum = 0;
            for (var j = rowPointers[i]; j < rowPointers[i + 1]; j++)
                sum += localValues[j] * vector[columnIndices[j]];
            result[i] = sum;
        });

        return result;
    }

    public double[] MultiplyTransposedParallel(double[] vector)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better error messages
        if (vector == null)
            throw new ArgumentNullException(nameof(vector), "Input vector cannot be null");
        if (vector.Length != nrows)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix rows {nrows} for transposed multiply",
                nameof(vector));

        var localValues = GetValuesOrThrow();
        var result = new double[ncols];

        if (nrows < MIN_ROWS_FOR_PARALLEL)
        {
            for (var i = 0; i < nrows; i++)
            {
                var xi = vector[i];
                for (var j = rowPointers[i]; j < rowPointers[i + 1]; j++)
                    result[columnIndices[j]] += localValues[j] * xi;
            }

            return result;
        }

        // FIXED: Issue #1 - Use Parallel.For with proper thread-local state pattern
        // This avoids the bug of using Thread.CurrentThread.ManagedThreadId which is not bounded.
        // The Parallel.For overload with localInit/localFinally ensures:
        // 1. Each thread gets its own accumulation array (no contention during computation)
        // 2. Thread-local arrays are properly initialized and finalized
        // 3. Final merge uses locking only once per thread, not per element

        var lockObj = new object();

        Parallel.For(0, nrows, ParallelConfig.Options,
            // localInit: Each thread gets its own result array
            () =>
            {
                // Use ArrayPool for large result arrays to reduce GC pressure
                if (ncols > 1024)
                {
                    var rented = ArrayPool<double>.Shared.Rent(ncols);
                    Array.Clear(rented, 0, ncols);
                    return (array: rented, isRented: true);
                }

                return (array: new double[ncols], isRented: false);
            },
            // body: Process one row, accumulate into thread-local array
            (i, loopState, localResult) =>
            {
                var xi = vector[i];
                var localArray = localResult.array;
                for (var j = rowPointers[i]; j < rowPointers[i + 1]; j++)
                    localArray[columnIndices[j]] += localValues[j] * xi;
                return localResult;
            },
            // localFinally: Merge thread-local results into global result
            localResult =>
            {
                var localArray = localResult.array;
                lock (lockObj)
                {
                    for (var col = 0; col < ncols; col++)
                        result[col] += localArray[col];
                }

                // Return rented array to pool
                if (localResult.isRented)
                    ArrayPool<double>.Shared.Return(localArray, true);
            }
        );

        return result;
    }

    #endregion

    #region SIMD Operations

    public double[] MultiplySIMD(double[] vector)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better error messages
        if (vector == null)
            throw new ArgumentNullException(nameof(vector), "Input vector cannot be null");
        if (vector.Length != ncols)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix columns {ncols}",
                nameof(vector));

        var localValues = GetValuesOrThrow();
        var result = new double[nrows];

        // FIXED: Validate column indices at runtime if matrix was constructed with skipValidation
        // This prevents memory corruption from out-of-bounds access in unsafe SIMD code
        if (constructedWithSkipValidation) ValidateColumnIndicesForMultiply(vector.Length);

        Parallel.For(0, nrows, ParallelConfig.Options, i =>
        {
            var kStart = rowPointers[i];
            var kEnd = rowPointers[i + 1];
            double sum;

            if (Avx512F.IsSupported && kEnd - kStart >= 8)
                sum = ComputeRowAVX512(kStart, kEnd, vector, localValues);
            else if (Avx2.IsSupported && kEnd - kStart >= 4)
                sum = ComputeRowAVX2(kStart, kEnd, vector, localValues);
            else if (VectorT.IsHardwareAccelerated && kEnd - kStart >= Vector<double>.Count)
                sum = ComputeRowVectorized(kStart, kEnd, vector, localValues);
            else
                sum = ComputeRowScalar(kStart, kEnd, vector, localValues);

            result[i] = sum;
        });

        return result;
    }

    /// <summary>
    ///     Validates that all column indices are within bounds for the given vector length.
    ///     Called before SIMD operations on matrices constructed with skipValidation.
    /// </summary>
    private void ValidateColumnIndicesForMultiply(int vectorLength)
    {
        for (var i = 0; i < columnIndices.Length; i++)
        {
            var col = columnIndices[i];
            if (col < 0 || col >= vectorLength)
                throw new InvalidOperationException(
                    $"Column index {col} at position {i} is out of bounds for vector of length {vectorLength}. " +
                    $"Matrix was constructed with skipValidation=true and contains invalid indices. " +
                    "Reconstruct the matrix with skipValidation=false to get detailed validation errors, " +
                    "or fix the source data that generated the invalid indices.");
        }
    }

    [MethodImpl(AggressiveInlining)]
    private unsafe double ComputeRowAVX512(int kStart, int kEnd, double[] vector, double[] localValues)
    {
        var k = kStart;
        var simdEnd = kStart + (kEnd - kStart) / 8 * 8;
        var vsum = Vector512<double>.Zero;

        fixed (double* pValues = localValues, pVector = vector)
        fixed (int* pCols = columnIndices)
        {
            for (; k < simdEnd; k += 8)
            {
                // FIXED: Issue #5 - Added bounds checking in debug mode
                Debug.Assert(k + 7 < columnIndices.Length, "Column index array bounds violation");
                Debug.Assert(pCols[k] >= 0 && pCols[k] < vector.Length, $"Column index {pCols[k]} out of bounds");
                Debug.Assert(pCols[k + 7] >= 0 && pCols[k + 7] < vector.Length,
                    $"Column index {pCols[k + 7]} out of bounds");

                var vdata = Avx512F.LoadVector512(pValues + k);
                var vvec = Vector512.Create(
                    pVector[pCols[k]], pVector[pCols[k + 1]],
                    pVector[pCols[k + 2]], pVector[pCols[k + 3]],
                    pVector[pCols[k + 4]], pVector[pCols[k + 5]],
                    pVector[pCols[k + 6]], pVector[pCols[k + 7]]
                );
                vsum = Avx512F.FusedMultiplyAdd(vdata, vvec, vsum);
            }
        }

        var sum = vsum[0] + vsum[1] + vsum[2] + vsum[3] + vsum[4] + vsum[5] + vsum[6] + vsum[7];
        for (; k < kEnd; k++) sum += localValues[k] * vector[columnIndices[k]];
        return sum;
    }

    [MethodImpl(AggressiveInlining)]
    private unsafe double ComputeRowAVX2(int kStart, int kEnd, double[] vector, double[] localValues)
    {
        var k = kStart;
        var simdEnd = kStart + (kEnd - kStart) / 4 * 4;
        var vsum = Vector256<double>.Zero;

        fixed (double* pValues = localValues, pVector = vector)
        fixed (int* pCols = columnIndices)
        {
            for (; k < simdEnd; k += 4)
            {
                // FIXED: Issue #5 - Added bounds checking
                Debug.Assert(k + 3 < columnIndices.Length, "Column index array bounds violation");
                Debug.Assert(pCols[k] >= 0 && pCols[k] < vector.Length, $"Column index {pCols[k]} out of bounds");

                var vdata = Avx.LoadVector256(pValues + k);
                var vvec = Vector256.Create(pVector[pCols[k]], pVector[pCols[k + 1]], pVector[pCols[k + 2]],
                    pVector[pCols[k + 3]]);
                vsum = Fma.IsSupported
                    ? Fma.MultiplyAdd(vdata, vvec, vsum)
                    : Avx.Add(vsum, Avx.Multiply(vdata, vvec));
            }
        }

        var sum = vsum[0] + vsum[1] + vsum[2] + vsum[3];
        for (; k < kEnd; k++) sum += localValues[k] * vector[columnIndices[k]];
        return sum;
    }

    [MethodImpl(AggressiveInlining)]
    private double ComputeRowVectorized(int kStart, int kEnd, double[] vector, double[] localValues)
    {
        var vectorSize = Vector<double>.Count;
        var k = kStart;
        var simdEnd = kStart + (kEnd - kStart) / vectorSize * vectorSize;
        var vsum = Vector<double>.Zero;

        // FIXED: Issue #9 - Move stackalloc outside the loop to prevent stack overflow
        // Note: Manual gather is used here for portability. Consider hardware gather (Avx2.GatherVector256)
        // for better performance on supported platforms.
        Span<double> gatherBuffer = stackalloc double[vectorSize];
        while (k < simdEnd)
        {
            var vdata = new Vector<double>(localValues, k);
            for (var i = 0; i < vectorSize; i++) gatherBuffer[i] = vector[columnIndices[k + i]];
            var vvec = new Vector<double>(gatherBuffer);
            vsum += vdata * vvec;
            k += vectorSize;
        }

        var sum = 0.0;
        for (var i = 0; i < vectorSize; i++) sum += vsum[i];
        for (; k < kEnd; k++) sum += localValues[k] * vector[columnIndices[k]];
        return sum;
    }

    [MethodImpl(AggressiveInlining)]
    private double ComputeRowScalar(int kStart, int kEnd, double[] vector, double[] localValues)
    {
        var k = kStart;
        var sum = 0.0;
        for (; k < kEnd; k++) sum += localValues[k] * vector[columnIndices[k]];
        return sum;
    }

    public double[] MultiplyAuto(double[] vector, bool preferGPU = true)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better error messages
        if (vector == null)
            throw new ArgumentNullException(nameof(vector), "Input vector cannot be null");
        if (vector.Length != ncols)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix columns {ncols}",
                nameof(vector));

        if (preferGPU && ParallelConfig.EnableGPU && isGpuInitialized && gpuAccelerator != null &&
            SparseBackendFactory.ShouldUseGPU(nrows, ncols, NonZeroCount))
            try
            {
                return MultiplyGPU(vector);
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"GPU multiply failed: {ex.Message}. Falling back to CPU.");
            }

        if (nrows >= MIN_ROWS_FOR_SIMD)
            if (Avx512F.IsSupported || Avx2.IsSupported || VectorT.IsHardwareAccelerated)
                return MultiplySIMD(vector);

        if (nrows >= MIN_ROWS_FOR_PARALLEL) return MultiplyParallel(vector);
        return Multiply(vector);
    }

    #endregion

    #region GPU Operations

    // FIXED: Issue #14 - Added check to avoid redundant GPU initialization
    public void InitializeGpu()
    {
        ThrowIfDisposed();

        // Skip if already initialized with valid state
        if (isGpuInitialized && gpuAccelerator != null)
        {
            Debug.WriteLine("GPU already initialized. Skipping redundant initialization.");
            return;
        }

        try
        {
            if (gpuAccelerator != null)
            {
                gpuAccelerator.Dispose();
                isGpuInitialized = false;
            }

            // Validates library loading via new backend
            gpuAccelerator = new CuSparseBackend(nrows, ncols, NonZeroCount);

            lock (syncLock)
            {
                if (values != null)
                    gpuAccelerator.Initialize(rowPointers, columnIndices, values);
            }

            isGpuInitialized = true;
        }
        catch (Exception ex)
        {
            Debug.WriteLine($"GPU initialization failed: {ex.Message}");
            gpuAccelerator?.Dispose();
            gpuAccelerator = null;
            isGpuInitialized = false;
            throw new InvalidOperationException(
                "Failed to initialize GPU acceleration. Ensure CUDA drivers are installed, " +
                "a compatible NVIDIA GPU is available, and the CUDA toolkit is properly configured.", ex);
        }
    }

    public double[] MultiplyGPU(double[] vector)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better error messages
        if (vector == null)
            throw new ArgumentNullException(nameof(vector), "Input vector cannot be null");
        if (vector.Length != ncols)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix columns {ncols}",
                nameof(vector));

        if (!isGpuInitialized || gpuAccelerator == null)
            throw new InvalidOperationException("GPU not initialized. Call InitializeGpu() first.");

        return gpuAccelerator.Multiply(vector);
    }

    public double[] MultiplyTransposedGPU(double[] vector)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better error messages
        if (vector == null)
            throw new ArgumentNullException(nameof(vector), "Input vector cannot be null");
        if (vector.Length != nrows)
            throw new ArgumentException(
                $"Vector length {vector.Length} does not match matrix rows {nrows} for transposed multiply",
                nameof(vector));

        if (!isGpuInitialized || gpuAccelerator == null)
            throw new InvalidOperationException("GPU not initialized. Call InitializeGpu() first.");

        return gpuAccelerator.MultiplyTransposed(vector);
    }

    public void UpdateGPUValues(double[] newValues)
    {
        ThrowIfDisposed();

        // FIXED: Issue #20 - Better error messages
        if (newValues == null)
            throw new ArgumentNullException(nameof(newValues), "New values array cannot be null");
        if (newValues.Length != NonZeroCount)
            throw new ArgumentException(
                $"Values array length {newValues.Length} must match non-zero count {NonZeroCount}",
                nameof(newValues));

        if (!isGpuInitialized || gpuAccelerator == null)
            throw new InvalidOperationException("GPU not initialized. Call InitializeGpu() first.");

        gpuAccelerator.UpdateValues(newValues);
    }

    #endregion

    #region Matrix-Matrix Operations

    public static CSR operator *(CSR a, CSR b)
    {
        a.ThrowIfDisposed();
        b.ThrowIfDisposed();
        return Multiply3Phase(a, b);
    }

    public static (int[] rowPtr, int nnz) MultiplySymbolic(CSR A, CSR B)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();

        if (A.ncols != B.nrows)
            throw new ArgumentException(
                $"Dimension mismatch: A columns ({A.ncols}) must equal B rows ({B.nrows}) for matrix multiplication");

        var rowCounts = new int[A.nrows];
        Parallel.For(0, A.nrows, ParallelConfig.Options, () => new int[B.ncols], (i, loopState, workspace) =>
        {
            var currentRow = i + 1;
            var count = 0;
            for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
            {
                var ca = A.columnIndices[ka];
                for (var kb = B.rowPointers[ca]; kb < B.rowPointers[ca + 1]; kb++)
                {
                    var cb = B.columnIndices[kb];
                    if (workspace[cb] != currentRow)
                    {
                        count++;
                        workspace[cb] = currentRow;
                    }
                }
            }

            rowCounts[i] = count;
            return workspace;
        }, _ => { });

        var rowPtr = new int[A.nrows + 1];
        rowPtr[0] = 0;
        long cumulativeNnz = 0;

        for (var i = 0; i < A.nrows; i++)
        {
            cumulativeNnz += rowCounts[i];

            // FIXED: Issue #4 - Check overflow BEFORE cast
            if (cumulativeNnz > int.MaxValue)
                throw new InvalidOperationException(
                    $"Result matrix too large: non-zero count {cumulativeNnz} exceeds int.MaxValue ({int.MaxValue}). " +
                    "Consider using a different algorithm, reducing matrix density, or using a 64-bit sparse library.");

            rowPtr[i + 1] = (int)cumulativeNnz;
        }

        return (rowPtr, rowPtr[A.nrows]);
    }

    public static int[] MultiplyIndices(CSR A, CSR B, int[] rowPtr)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();
        var colIdx = new int[rowPtr[^1]];

        Parallel.For(0, A.nrows, ParallelConfig.Options, () => new int[B.ncols], (i, loopState, workspace) =>
        {
            var currentRow = i + 1;
            var pos = rowPtr[i];
            for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
            {
                var ca = A.columnIndices[ka];
                for (var kb = B.rowPointers[ca]; kb < B.rowPointers[ca + 1]; kb++)
                {
                    var cb = B.columnIndices[kb];
                    if (workspace[cb] != currentRow)
                    {
                        colIdx[pos++] = cb;
                        workspace[cb] = currentRow;
                    }
                }
            }

            return workspace;
        }, _ => { });
        return colIdx;
    }

    // Each thread writes directly to finalValues for its own row range.
    // Since rows are partitioned across threads, each thread's output positions
    // (rowPtr[startRow]..rowPtr[endRow]) are non-overlapping, eliminating the
    // need for thread-local copies and the merge phase.
    public static double[] MultiplyValues(CSR A, CSR B, int[] rowPtr, int[] colIdx)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();
        var aValues = A.GetValuesInternal();
        var bValues = B.GetValuesInternal();
        var finalValues = new double[colIdx.Length];
        var numThreads = ParallelConfig.MaxDegreeOfParallelism;

        Parallel.For(0, numThreads, ParallelConfig.Options, tid =>
        {
            var rowsPerThread = (A.nrows + numThreads - 1) / numThreads;
            var startRow = tid * rowsPerThread;
            var endRow = Math.Min(startRow + rowsPerThread, A.nrows);
            var avgNnzPerRow = colIdx.Length / Math.Max(1, A.nrows);
            var colToPos = new Dictionary<int, int>(avgNnzPerRow * 2);

            for (var i = startRow; i < endRow; i++)
            {
                colToPos.Clear();
                for (var k = rowPtr[i]; k < rowPtr[i + 1]; k++) colToPos[colIdx[k]] = k;

                for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                {
                    var aval = aValues[ka];
                    var ca = A.columnIndices[ka];
                    for (var kb = B.rowPointers[ca]; kb < B.rowPointers[ca + 1]; kb++)
                    {
                        var cb = B.columnIndices[kb];
                        var pos = colToPos[cb];
                        finalValues[pos] += aval * bValues[kb];
                    }
                }
            }
        });

        return finalValues;
    }

    public static CSR Multiply3Phase(CSR A, CSR B)
    {
        var (rowPtr, _) = MultiplySymbolic(A, B);
        var colIdx = MultiplyIndices(A, B, rowPtr);

        double[]? values = null;
        if (A.HasValues && B.HasValues) values = MultiplyValues(A, B, rowPtr, colIdx);

        // FIXED: Issue #17 - Optional zero-dropping could be added here
        // For now, keeping as-is to preserve exact numerical behavior
        return new CSR(rowPtr, colIdx, B.ncols, values, true);
    }

    public static CSR MultiplySymbolicOnly(CSR A, CSR B)
    {
        var (rowPtr, _) = MultiplySymbolic(A, B);
        var colIdx = MultiplyIndices(A, B, rowPtr);
        return new CSR(rowPtr, colIdx, B.ncols, skipValidation: true);
    }

    public static CSR TransposeAndMultiply(CSR A, CSR B)
    {
        if (A.nrows != B.nrows)
            throw new ArgumentException(
                $"Row count mismatch: A has {A.nrows} rows, B has {B.nrows} rows");

        A.ThrowIfDisposed();
        B.ThrowIfDisposed();

        var aValues = A.GetValuesOrThrow();
        var bValues = B.GetValuesOrThrow();
        var resultRows = A.ncols;
        var resultCols = B.ncols;

        var aTransposeStart = new int[resultRows];
        var aTransposeCount = new int[resultRows];
        for (var k = 0; k < A.columnIndices.Length; k++) aTransposeCount[A.columnIndices[k]]++;

        aTransposeStart[0] = 0;
        for (var i = 1; i < resultRows; i++) aTransposeStart[i] = aTransposeStart[i - 1] + aTransposeCount[i - 1];

        var aTotalNnz = aTransposeStart[resultRows - 1] + aTransposeCount[resultRows - 1];
        var aTransposeCol = new int[aTotalNnz];
        var aTransposeVal = new double[aTotalNnz];
        var aTransposePos = new int[resultRows];

        for (var row = 0; row < A.nrows; row++)
        for (var k = A.rowPointers[row]; k < A.rowPointers[row + 1]; k++)
        {
            var col = A.columnIndices[k];
            var idx = aTransposeStart[col] + aTransposePos[col];
            aTransposeCol[idx] = row;
            aTransposeVal[idx] = aValues[k];
            aTransposePos[col]++;
        }

        var resultRowCounts = new int[resultRows];
        Parallel.For(0, resultRows, ParallelConfig.Options, () => new int[resultCols], (i, loopState, workspace) =>
        {
            var timestamp = i + 1;
            var count = 0;
            var start = i < resultRows ? aTransposeStart[i] : aTotalNnz;
            var end = i < resultRows - 1 ? aTransposeStart[i + 1] : aTotalNnz;

            for (var ka = start; ka < end; ka++)
            {
                var k = aTransposeCol[ka];
                for (var kb = B.rowPointers[k]; kb < B.rowPointers[k + 1]; kb++)
                {
                    var j = B.columnIndices[kb];
                    if (workspace[j] != timestamp)
                    {
                        workspace[j] = timestamp;
                        count++;
                    }
                }
            }

            resultRowCounts[i] = count;
            return workspace;
        }, _ => { });

        var resultRowPtrs = new int[resultRows + 1];
        resultRowPtrs[0] = 0;
        for (var i = 0; i < resultRows; i++) resultRowPtrs[i + 1] = resultRowPtrs[i] + resultRowCounts[i];

        var resultNnz = resultRowPtrs[resultRows];
        var resultColIndices = new int[resultNnz];
        var resultValues = new double[resultNnz];

        Parallel.For(0, resultRows, ParallelConfig.Options, () => new int[resultCols], (i, loopState, workspace) =>
        {
            var timestamp = i + 1;
            var pos = resultRowPtrs[i];
            var start = i < resultRows ? aTransposeStart[i] : aTotalNnz;
            var end = i < resultRows - 1 ? aTransposeStart[i + 1] : aTotalNnz;

            for (var ka = start; ka < end; ka++)
            {
                var k = aTransposeCol[ka];
                for (var kb = B.rowPointers[k]; kb < B.rowPointers[k + 1]; kb++)
                {
                    var j = B.columnIndices[kb];
                    if (workspace[j] != timestamp)
                    {
                        resultColIndices[pos++] = j;
                        workspace[j] = timestamp;
                    }
                }
            }

            return workspace;
        }, _ => { });

        // FIXED: MEM-C1 - Eliminated massive per-thread allocation.
        // Each thread now writes directly to resultValues for its own row range,
        // since rows are partitioned and output positions are non-overlapping.
        var numThreads = ParallelConfig.MaxDegreeOfParallelism;
        Parallel.For(0, numThreads, ParallelConfig.Options, tid =>
        {
            var rowsPerThread = (resultRows + numThreads - 1) / numThreads;
            var rowStart = tid * rowsPerThread;
            var rowEnd = Math.Min(rowStart + rowsPerThread, resultRows);

            for (var i = rowStart; i < rowEnd; i++)
            {
                var posMap = new Dictionary<int, int>(resultRowCounts[i]);
                for (var p = resultRowPtrs[i]; p < resultRowPtrs[i + 1]; p++) posMap[resultColIndices[p]] = p;

                var start = i < resultRows ? aTransposeStart[i] : aTotalNnz;
                var end = i < resultRows - 1 ? aTransposeStart[i + 1] : aTotalNnz;

                for (var ka = start; ka < end; ka++)
                {
                    var k = aTransposeCol[ka];
                    var aVal = aTransposeVal[ka];
                    for (var kb = B.rowPointers[k]; kb < B.rowPointers[k + 1]; kb++)
                    {
                        var pos = posMap[B.columnIndices[kb]];
                        resultValues[pos] += aVal * bValues[kb];
                    }
                }
            }
        });

        return new CSR(resultRowPtrs, resultColIndices, resultCols, resultValues, true);
    }

    #endregion

    #region Matrix Addition Operations

    public static int[] AddSymbolic(CSR A, CSR B)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();

        if (A.nrows != B.nrows || A.ncols != B.ncols)
            throw new ArgumentException(
                $"Dimension mismatch: A is {A.nrows}x{A.ncols}, B is {B.nrows}x{B.ncols}");

        var rowPtr = new int[A.nrows + 1];

        if (A.nrows >= MIN_ROWS_FOR_PARALLEL)
        {
            // Parallel: compute row counts independently, then prefix sum
            var rowCounts = new int[A.nrows];
            Parallel.For(0, A.nrows, ParallelConfig.Options,
                () => ArrayPool<int>.Shared.Rent(A.ncols),
                (i, loopState, workspace) =>
                {
                    var nnzThisRow = 0;
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    {
                        var col = A.columnIndices[ka];
                        if (workspace[col] != i + 1)
                        {
                            workspace[col] = i + 1;
                            nnzThisRow++;
                        }
                    }

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (workspace[col] != i + 1)
                        {
                            workspace[col] = i + 1;
                            nnzThisRow++;
                        }
                    }

                    rowCounts[i] = nnzThisRow;
                    return workspace;
                },
                workspace => ArrayPool<int>.Shared.Return(workspace, true));

            rowPtr[0] = 0;
            for (var i = 0; i < A.nrows; i++)
                rowPtr[i + 1] = rowPtr[i] + rowCounts[i];
        }
        else
        {
            var workspace = ArrayPool<int>.Shared.Rent(A.ncols);
            try
            {
                Array.Clear(workspace, 0, workspace.Length);
                var totalNnz = 0;
                rowPtr[0] = 0;

                for (var i = 0; i < A.nrows; i++)
                {
                    var nnzThisRow = 0;
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    {
                        var col = A.columnIndices[ka];
                        if (workspace[col] == 0)
                        {
                            workspace[col] = 1;
                            nnzThisRow++;
                        }
                    }

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (workspace[col] == 0)
                        {
                            workspace[col] = 1;
                            nnzThisRow++;
                        }
                    }

                    // Cleanup
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++) workspace[A.columnIndices[ka]] = 0;
                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++) workspace[B.columnIndices[kb]] = 0;

                    totalNnz += nnzThisRow;
                    rowPtr[i + 1] = totalNnz;
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(workspace);
            }
        }

        return rowPtr;
    }

    public static int[] AddIndices(CSR A, CSR B, int[] rowPtr)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();
        var colIdx = new int[rowPtr[^1]];

        if (A.nrows >= MIN_ROWS_FOR_PARALLEL)
        {
            // Parallel: each row writes to its own non-overlapping segment in colIdx
            Parallel.For(0, A.nrows, ParallelConfig.Options,
                () => ArrayPool<int>.Shared.Rent(A.ncols),
                (i, loopState, workspace) =>
                {
                    var pos = rowPtr[i];
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    {
                        var col = A.columnIndices[ka];
                        if (workspace[col] != i + 1)
                        {
                            colIdx[pos++] = col;
                            workspace[col] = i + 1;
                        }
                    }

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (workspace[col] != i + 1)
                        {
                            colIdx[pos++] = col;
                            workspace[col] = i + 1;
                        }
                    }

                    return workspace;
                },
                workspace => ArrayPool<int>.Shared.Return(workspace, true));
        }
        else
        {
            var workspace = ArrayPool<int>.Shared.Rent(A.ncols);
            try
            {
                Array.Clear(workspace, 0, workspace.Length);
                var pos = 0;
                for (var i = 0; i < A.nrows; i++)
                {
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    {
                        var col = A.columnIndices[ka];
                        if (workspace[col] == 0)
                        {
                            colIdx[pos] = col;
                            workspace[col] = pos + 1;
                            pos++;
                        }
                    }

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (workspace[col] == 0)
                        {
                            colIdx[pos] = col;
                            workspace[col] = pos + 1;
                            pos++;
                        }
                    }

                    for (var k = rowPtr[i]; k < rowPtr[i + 1]; k++) workspace[colIdx[k]] = 0;
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(workspace);
            }
        }

        return colIdx;
    }

    public static double[] AddValues(CSR A, CSR B, double alpha, double beta, int[] rowPtr, int[] colIdx)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();
        var aValues = A.GetValuesInternal();
        var bValues = B.GetValuesInternal();
        var values = new double[colIdx.Length];

        if (A.nrows >= MIN_ROWS_FOR_PARALLEL)
        {
            // Parallel: each row writes to its own non-overlapping segment in values.
            // Uses workspace as a column→position map (storing position+1, with 0 meaning "unseen").
            // This gives O(1) lookup for duplicate columns, avoiding the previous O(nnz_per_row) linear search.
            Parallel.For(0, A.nrows, ParallelConfig.Options,
                () => ArrayPool<int>.Shared.Rent(A.ncols),
                (i, loopState, workspace) =>
                {
                    // First, build position map from the colIdx that AddIndices produced for this row.
                    // workspace[col] = position_in_values + 1 (0 means not yet seen).
                    // We use row-local timestamp (i+1) shifted approach: store actual output position.
                    var rowStart = rowPtr[i];
                    var rowEnd = rowPtr[i + 1];

                    // Clear workspace entries for columns in this row's output segment
                    for (var k = rowStart; k < rowEnd; k++)
                        workspace[colIdx[k]] = k + 1; // +1 so that 0 means "no mapping"

                    // Accumulate A's contributions
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    {
                        var col = A.columnIndices[ka];
                        var mappedPos = workspace[col] - 1; // O(1) lookup
                        values[mappedPos] += alpha * aValues[ka];
                    }

                    // Accumulate B's contributions
                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        var mappedPos = workspace[col] - 1; // O(1) lookup
                        values[mappedPos] += beta * bValues[kb];
                    }

                    // Reset workspace entries for this row (only touch columns we used)
                    for (var k = rowStart; k < rowEnd; k++)
                        workspace[colIdx[k]] = 0;

                    return workspace;
                },
                workspace => ArrayPool<int>.Shared.Return(workspace, true));
        }
        else
        {
            var workspace = ArrayPool<int>.Shared.Rent(A.ncols);
            try
            {
                Array.Clear(workspace, 0, workspace.Length);
                var pos = 0;
                for (var i = 0; i < A.nrows; i++)
                {
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    {
                        var col = A.columnIndices[ka];
                        if (workspace[col] == 0)
                        {
                            values[pos] = alpha * aValues[ka];
                            workspace[col] = pos + 1;
                            pos++;
                        }
                        else
                        {
                            values[workspace[col] - 1] += alpha * aValues[ka];
                        }
                    }

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (workspace[col] == 0)
                        {
                            values[pos] = beta * bValues[kb];
                            workspace[col] = pos + 1;
                            pos++;
                        }
                        else
                        {
                            values[workspace[col] - 1] += beta * bValues[kb];
                        }
                    }

                    for (var k = rowPtr[i]; k < rowPtr[i + 1]; k++) workspace[colIdx[k]] = 0;
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(workspace);
            }
        }

        return values;
    }

    public static CSR Add(CSR A, CSR B, double alpha = 1.0, double beta = 1.0)
    {
        var rowPtr = AddSymbolic(A, B);
        var colIdx = AddIndices(A, B, rowPtr);

        double[]? values = null;
        if (A.HasValues && B.HasValues)
            values = AddValues(A, B, alpha, beta, rowPtr, colIdx);

        return new CSR(rowPtr, colIdx, A.ncols, values, true);
    }

    public static CSR operator +(CSR a, CSR b)
    {
        return Add(a, b);
    }

    public static CSR operator -(CSR a, CSR b)
    {
        return Add(a, b, 1.0, -1.0);
    }

    public static CSR operator *(double alpha, CSR a)
    {
        a.ThrowIfDisposed();

        double[]? newValues = null;
        if (a.HasValues)
        {
            var aValues = a.GetValuesInternal();
            newValues = new double[aValues.Length];
            for (var i = 0; i < aValues.Length; i++) newValues[i] = alpha * aValues[i];
        }

        return new CSR((int[])a.rowPointers.Clone(), (int[])a.columnIndices.Clone(), a.ncols, newValues, true);
    }

    /// <summary>Scalar multiplication: C = A * alpha (commutative with double * CSR).</summary>
    public static CSR operator *(CSR a, double alpha)
    {
        return alpha * a;
    }

    #endregion

    #region Matrix Intersection Operations

    /// <summary>
    ///     Computes the row pointers for the intersection of two CSR matrices.
    ///     The intersection contains only entries where both matrices have non-zeros.
    /// </summary>
    /// <param name="A">First matrix.</param>
    /// <param name="B">Second matrix.</param>
    /// <returns>Row pointers array for the intersection matrix.</returns>
    public static int[] IntersectionSymbolic(CSR A, CSR B)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();

        if (A.nrows != B.nrows || A.ncols != B.ncols)
            throw new ArgumentException(
                $"Dimension mismatch: A is {A.nrows}x{A.ncols}, B is {B.nrows}x{B.ncols}");

        var rowPtr = new int[A.nrows + 1];

        if (A.nrows >= MIN_ROWS_FOR_PARALLEL)
        {
            // Parallel: compute row counts independently, then prefix sum
            var rowCounts = new int[A.nrows];
            Parallel.For(0, A.nrows, ParallelConfig.Options,
                () => ArrayPool<int>.Shared.Rent(A.ncols),
                (i, loopState, workspace) =>
                {
                    // Mark all columns present in row i of A using timestamp trick
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        workspace[A.columnIndices[ka]] = i + 1;

                    // Count columns that are also present in row i of B
                    var nnzThisRow = 0;
                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                        if (workspace[B.columnIndices[kb]] == i + 1)
                            nnzThisRow++;

                    rowCounts[i] = nnzThisRow;
                    return workspace;
                },
                workspace => ArrayPool<int>.Shared.Return(workspace, true));

            rowPtr[0] = 0;
            for (var i = 0; i < A.nrows; i++)
                rowPtr[i + 1] = rowPtr[i] + rowCounts[i];
        }
        else
        {
            var marker = ArrayPool<int>.Shared.Rent(A.ncols);
            try
            {
                Array.Fill(marker, -1, 0, A.ncols);
                var totalNnz = 0;
                rowPtr[0] = 0;

                for (var i = 0; i < A.nrows; i++)
                {
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        marker[A.columnIndices[ka]] = i;

                    var nnzThisRow = 0;
                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                        if (marker[B.columnIndices[kb]] == i)
                            nnzThisRow++;

                    totalNnz += nnzThisRow;
                    rowPtr[i + 1] = totalNnz;
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(marker);
            }
        }

        return rowPtr;
    }

    /// <summary>
    ///     Computes the column indices for the intersection of two CSR matrices.
    /// </summary>
    /// <param name="A">First matrix.</param>
    /// <param name="B">Second matrix.</param>
    /// <param name="rowPtr">Row pointers from IntersectionSymbolic.</param>
    /// <returns>Column indices array for the intersection matrix.</returns>
    public static int[] IntersectionIndices(CSR A, CSR B, int[] rowPtr)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();

        var colIdx = new int[rowPtr[^1]];

        if (A.nrows >= MIN_ROWS_FOR_PARALLEL)
        {
            // Parallel: each row writes to its own non-overlapping segment in colIdx
            Parallel.For(0, A.nrows, ParallelConfig.Options,
                () => ArrayPool<int>.Shared.Rent(A.ncols),
                (i, loopState, workspace) =>
                {
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        workspace[A.columnIndices[ka]] = i + 1;

                    var pos = rowPtr[i];
                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (workspace[col] == i + 1)
                            colIdx[pos++] = col;
                    }

                    return workspace;
                },
                workspace => ArrayPool<int>.Shared.Return(workspace, true));
        }
        else
        {
            var marker = ArrayPool<int>.Shared.Rent(A.ncols);
            try
            {
                Array.Fill(marker, -1, 0, A.ncols);
                var pos = 0;

                for (var i = 0; i < A.nrows; i++)
                {
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        marker[A.columnIndices[ka]] = i;

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (marker[col] == i)
                            colIdx[pos++] = col;
                    }
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(marker);
            }
        }

        return colIdx;
    }

    /// <summary>
    ///     Computes the values for the intersection of two CSR matrices using element-wise (Hadamard) product.
    /// </summary>
    /// <param name="A">First matrix.</param>
    /// <param name="B">Second matrix.</param>
    /// <param name="rowPtr">Row pointers from IntersectionSymbolic.</param>
    /// <param name="colIdx">Column indices from IntersectionIndices.</param>
    /// <returns>Values array for the intersection matrix.</returns>
    public static double[] IntersectionValues(CSR A, CSR B, int[] rowPtr, int[] colIdx)
    {
        A.ThrowIfDisposed();
        B.ThrowIfDisposed();

        var aValues = A.GetValuesInternal();
        var bValues = B.GetValuesInternal();
        var values = new double[colIdx.Length];

        if (A.nrows >= MIN_ROWS_FOR_PARALLEL)
        {
            // Parallel: each row writes to its own non-overlapping segment in values
            Parallel.For(0, A.nrows, ParallelConfig.Options,
                () => ArrayPool<int>.Shared.Rent(A.ncols),
                (i, loopState, workspace) =>
                {
                    // Store index into A's values for each column (using timestamp+offset encoding)
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        workspace[A.columnIndices[ka]] = ka + 1; // +1 so 0 = "not present"

                    var pos = rowPtr[i];
                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        var kaEncoded = workspace[col];
                        if (kaEncoded > 0) // Column exists in A for this row
                            values[pos++] = aValues[kaEncoded - 1] * bValues[kb];
                    }

                    // Clear markers for this row
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        workspace[A.columnIndices[ka]] = 0;

                    return workspace;
                },
                workspace => ArrayPool<int>.Shared.Return(workspace, true));
        }
        else
        {
            var markerIdx = ArrayPool<int>.Shared.Rent(A.ncols);
            try
            {
                Array.Fill(markerIdx, -1, 0, A.ncols);
                var pos = 0;

                for (var i = 0; i < A.nrows; i++)
                {
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        markerIdx[A.columnIndices[ka]] = ka;

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        var ka = markerIdx[col];
                        if (ka != -1)
                            values[pos++] = aValues[ka] * bValues[kb];
                    }

                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                        markerIdx[A.columnIndices[ka]] = -1;
                }
            }
            finally
            {
                ArrayPool<int>.Shared.Return(markerIdx);
            }
        }

        return values;
    }

    /// <summary>
    ///     Computes the intersection of two CSR matrices.
    ///     The result contains entries only where both matrices have non-zeros,
    ///     with values computed as the element-wise (Hadamard) product.
    ///     If either matrix lacks values, returns a symbolic-only result.
    /// </summary>
    /// <param name="A">First matrix.</param>
    /// <param name="B">Second matrix.</param>
    /// <returns>A new CSR matrix representing the intersection.</returns>
    public static CSR Intersection(CSR A, CSR B)
    {
        var rowPtr = IntersectionSymbolic(A, B);
        var colIdx = IntersectionIndices(A, B, rowPtr);

        double[]? values = null;
        if (A.HasValues && B.HasValues)
            values = IntersectionValues(A, B, rowPtr, colIdx);

        return new CSR(rowPtr, colIdx, A.ncols, values, true);
    }

    /// <summary>
    ///     Computes the symbolic intersection of two CSR matrices (sparsity pattern only, no values).
    /// </summary>
    /// <param name="A">First matrix.</param>
    /// <param name="B">Second matrix.</param>
    /// <returns>A new CSR matrix representing the intersection pattern without values.</returns>
    public static CSR IntersectionSymbolicOnly(CSR A, CSR B)
    {
        var rowPtr = IntersectionSymbolic(A, B);
        var colIdx = IntersectionIndices(A, B, rowPtr);
        return new CSR(rowPtr, colIdx, A.ncols, skipValidation: true);
    }

    /// <summary>
    ///     Intersection operator. Returns a matrix with entries only where both operands have non-zeros.
    ///     Values are computed as the element-wise (Hadamard) product: C[i,j] = A[i,j] * B[i,j].
    ///     If either operand lacks values, returns a symbolic-only result.
    /// </summary>
    /// <param name="a">First matrix.</param>
    /// <param name="b">Second matrix.</param>
    /// <returns>The intersection of the two matrices.</returns>
    public static CSR operator &(CSR a, CSR b)
    {
        return Intersection(a, b);
    }

    #endregion

    #region Transpose Operations

    public CSR Transpose()
    {
        // CRITICAL FIX 2: Move disposal check inside lock to prevent race condition
        // Previous code checked disposal before lock, allowing Dispose() to nullify
        // cache between the check and accessing cache.Value
        lock (syncLock)
        {
            ThrowIfDisposed();

            var currentMod = modificationCount;
            var cache = cachedTranspose;

            if (cache != null && cachedTransposeModCount == currentMod)
                // Safe: we hold syncLock, so Dispose() can't nullify cache
                try
                {
                    return cache.Value;
                }
                catch (Exception ex)
                {
                    // If Lazy initialization failed, clear cache and rethrow with context
                    cachedTranspose = null;
                    throw new InvalidOperationException(
                        "Failed to compute matrix transpose. See inner exception for details.", ex);
                }

            // Cache is stale or doesn't exist - create new one
            cachedTranspose = new Lazy<CSR>(() => TransposeParallelUncached(),
                LazyThreadSafetyMode.ExecutionAndPublication);
            cachedTransposeModCount = currentMod;

            try
            {
                return cachedTranspose.Value;
            }
            catch (Exception ex)
            {
                // Clear failed cache
                cachedTranspose = null;
                throw new InvalidOperationException(
                    "Failed to compute matrix transpose. See inner exception for details.", ex);
            }
        }
    }

    /// <summary>
    ///     Computes the transpose along with position tracking.
    ///     For each entry <c>k</c> in the transposed matrix, <c>positionTracking[k]</c> gives
    ///     the 0-based index of that entry within its original row in the source matrix.
    ///     For example, if transposed entry <c>k</c> came from the 3rd non-zero in its original row,
    ///     <c>positionTracking[k] == 2</c>.
    /// </summary>
    public (CSR transpose, int[] positionTracking) TransposeWithPositions()
    {
        ThrowIfDisposed();
        var colCounts = new int[ncols];
        Parallel.For(0, nrows, ParallelConfig.Options, i =>
        {
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
                Interlocked.Increment(ref colCounts[columnIndices[k]]);
        });

        var rowPtrT = new int[ncols + 1];
        rowPtrT[0] = 0;
        for (var c = 0; c < ncols; c++) rowPtrT[c + 1] = rowPtrT[c] + colCounts[c];

        var nnz = rowPtrT[ncols];
        var colIdxT = new int[nnz];
        double[]? localValues;
        lock (syncLock)
        {
            localValues = values;
        }

        var valuesT = localValues != null ? new double[nnz] : null;
        var positions = new int[nnz];
        var next = new int[ncols];
        for (var c = 0; c < ncols; c++) next[c] = rowPtrT[c];

        Parallel.For(0, nrows, ParallelConfig.Options, i =>
        {
            var kOrigRow = 0;
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            {
                var col = columnIndices[k];
                var posT = Interlocked.Increment(ref next[col]) - 1;
                colIdxT[posT] = i;
                if (valuesT != null && localValues != null) valuesT[posT] = localValues[k];
                positions[posT] = kOrigRow;
                kOrigRow++;
            }
        });
        return (new CSR(rowPtrT, colIdxT, nrows, valuesT, true), positions);
    }

    private CSR TransposeParallelUncached()
    {
        var colCounts = new int[ncols];
        Parallel.For(0, nrows, ParallelConfig.Options, i =>
        {
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
                Interlocked.Increment(ref colCounts[columnIndices[k]]);
        });

        var rowPtrT = new int[ncols + 1];
        rowPtrT[0] = 0;
        for (var c = 0; c < ncols; c++) rowPtrT[c + 1] = rowPtrT[c] + colCounts[c];

        var nnz = rowPtrT[ncols];
        var colIdxT = new int[nnz];
        double[]? localValues;
        lock (syncLock)
        {
            localValues = values;
        }

        var valuesT = localValues != null ? new double[nnz] : null;
        var next = new int[ncols];
        for (var c = 0; c < ncols; c++) next[c] = rowPtrT[c];

        Parallel.For(0, nrows, ParallelConfig.Options, i =>
        {
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            {
                var col = columnIndices[k];
                var pos = Interlocked.Increment(ref next[col]) - 1;
                colIdxT[pos] = i;
                if (valuesT != null && localValues != null) valuesT[pos] = localValues[k];
            }
        });
        return new CSR(rowPtrT, colIdxT, nrows, valuesT, true);
    }

    public CSR TransposeParallel()
    {
        ThrowIfDisposed();
        return TransposeParallelUncached();
    }

    #endregion

    #region PARDISO Solver Integration

    public double[] SolvePardiso(double[] rhs, int matrixType = 11, bool refresh = false)
    {
        ThrowIfDisposed();

        if (rhs == null)
            throw new ArgumentNullException(nameof(rhs), "Right-hand side vector cannot be null");
        if (rhs.Length != nrows)
            throw new ArgumentException(
                $"RHS length {rhs.Length} does not match matrix rows {nrows}",
                nameof(rhs));

        using var solver = new PardisoSolver(matrixType);
        return solver.Solve(this, rhs, refresh);
    }

    public double[] SolvePardisoMultiple(double[] rhs, int nrhs, int matrixType = 11, bool refresh = false)
    {
        ThrowIfDisposed();

        if (rhs == null)
            throw new ArgumentNullException(nameof(rhs), "Right-hand side vector cannot be null");
        // FIXED: Integer overflow risk in nrows * nrhs - use long arithmetic
        var expectedLength = (long)nrows * nrhs;
        if (expectedLength > int.MaxValue)
            throw new ArgumentException(
                $"Product nrows ({nrows}) * nrhs ({nrhs}) = {expectedLength} exceeds int.MaxValue",
                nameof(nrhs));
        if (rhs.Length != (int)expectedLength)
            throw new ArgumentException(
                $"RHS must have {expectedLength} elements for {nrhs} right-hand sides (got {rhs.Length})",
                nameof(rhs));

        using var solver = new PardisoSolver(matrixType);
        return solver.SolveMultiple(this, rhs, nrhs, refresh);
    }

    #endregion

    #region Helper Utility Methods

    public static int[] CountsToRowPointers(int[] counts)
    {
        var rowPtr = new int[counts.Length + 1];
        rowPtr[0] = 0;

        // MEDIUM FIX M11: Add overflow detection for very large sparse matrices
        long accumulated = 0;
        for (var i = 0; i < counts.Length; i++)
        {
            accumulated += counts[i];
            if (accumulated > int.MaxValue)
                throw new OverflowException(
                    $"Matrix too large: total non-zeros ({accumulated:N0}) exceeds int.MaxValue. " +
                    "Consider using 64-bit indices or splitting the matrix.");
            rowPtr[i + 1] = (int)accumulated;
        }

        return rowPtr;
    }

    public static int[] RowPointersToCounts(int[] rowPtr)
    {
        var n = rowPtr.Length - 1;
        var counts = new int[n];
        for (var i = 0; i < n; i++) counts[i] = rowPtr[i + 1] - rowPtr[i];
        return counts;
    }

    public static int GetMaximumColumn(int[] colIdx)
    {
        if (colIdx.Length == 0) return -1;
        var max = colIdx[0];
        for (var i = 1; i < colIdx.Length; i++)
            if (colIdx[i] > max)
                max = colIdx[i];
        return max;
    }

    #endregion

    #region Matrix Utility Methods

    public CSR ExtractSymmetricPart(bool upper = true)
    {
        ThrowIfDisposed();
        var localValues = GetValuesOrThrow();
        var rowPtr = new int[nrows + 1];
        var counts = new int[nrows];

        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
        {
            var j = columnIndices[k];
            if (upper ? j >= i : j <= i) counts[i]++;
        }

        rowPtr[0] = 0;
        // MEDIUM FIX M11: Add overflow detection
        long accumulated = 0;
        for (var i = 0; i < nrows; i++)
        {
            accumulated += counts[i];
            if (accumulated > int.MaxValue)
                throw new OverflowException(
                    $"Symmetric extraction failed: resulting matrix has {accumulated:N0} non-zeros, " +
                    "which exceeds int.MaxValue");
            rowPtr[i + 1] = (int)accumulated;
        }

        var nnz = rowPtr[nrows];
        var colIdx = new int[nnz];
        var vals = new double[nnz];
        Array.Clear(counts, 0, nrows);

        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
        {
            var j = columnIndices[k];
            if (upper ? j >= i : j <= i)
            {
                var pos = rowPtr[i] + counts[i];
                colIdx[pos] = j;
                vals[pos] = localValues[k];
                counts[i]++;
            }
        }

        return new CSR(rowPtr, colIdx, ncols, vals, true);
    }

    public CSR ExtractSubmatrix(int[] rows, int[] cols)
    {
        ThrowIfDisposed();

        if (rows == null) throw new ArgumentNullException(nameof(rows));
        if (cols == null) throw new ArgumentNullException(nameof(cols));

        // Validate row and column indices
        for (var i = 0; i < rows.Length; i++)
            if (rows[i] < 0 || rows[i] >= nrows)
                throw new ArgumentOutOfRangeException(nameof(rows),
                    $"Row index {rows[i]} at position {i} is out of range [0, {nrows})");

        for (var i = 0; i < cols.Length; i++)
            if (cols[i] < 0 || cols[i] >= ncols)
                throw new ArgumentOutOfRangeException(nameof(cols),
                    $"Column index {cols[i]} at position {i} is out of range [0, {ncols})");

        var colSet = hashSetPool.Get();
        try
        {
            foreach (var col in cols) colSet.Add(col);
            var colMap = new Dictionary<int, int>();
            for (var j = 0; j < cols.Length; j++) colMap[cols[j]] = j;

            double[]? localValues;
            lock (syncLock)
            {
                localValues = values;
            }

            // First pass: count non-zeros per row to build row pointers
            var resultRowPtrs = new int[rows.Length + 1];
            for (var i = 0; i < rows.Length; i++)
            {
                var count = 0;
                for (var k = rowPointers[rows[i]]; k < rowPointers[rows[i] + 1]; k++)
                    if (colSet.Contains(columnIndices[k]))
                        count++;
                resultRowPtrs[i + 1] = resultRowPtrs[i] + count;
            }

            var totalNnz = resultRowPtrs[rows.Length];
            var resultColIdx = new int[totalNnz];
            var resultVals = localValues != null ? new double[totalNnz] : null;

            // Second pass: fill column indices and values
            for (var i = 0; i < rows.Length; i++)
            {
                var pos = resultRowPtrs[i];
                for (var k = rowPointers[rows[i]]; k < rowPointers[rows[i] + 1]; k++)
                {
                    var col = columnIndices[k];
                    if (colSet.Contains(col))
                    {
                        resultColIdx[pos] = colMap[col];
                        if (resultVals != null && localValues != null)
                            resultVals[pos] = localValues[k];
                        pos++;
                    }
                }
            }

            return new CSR(resultRowPtrs, resultColIdx, cols.Length, resultVals, true);
        }
        finally
        {
            hashSetPool.Return(colSet);
        }
    }

    public CSR ScaleRows(double[] factors)
    {
        ThrowIfDisposed();

        if (factors == null)
            throw new ArgumentNullException(nameof(factors), "Scale factors array cannot be null");
        if (factors.Length != nrows)
            throw new ArgumentException(
                $"Factors length {factors.Length} must match matrix rows {nrows}",
                nameof(factors));

        var localValues = GetValuesOrThrow();
        var newValues = new double[localValues.Length];
        for (var i = 0; i < nrows; i++)
        {
            var factor = factors[i];
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++) newValues[k] = factor * localValues[k];
        }

        return new CSR((int[])rowPointers.Clone(), (int[])columnIndices.Clone(), ncols, newValues, true);
    }

    public CSR ScaleColumns(double[] factors)
    {
        ThrowIfDisposed();

        if (factors == null)
            throw new ArgumentNullException(nameof(factors), "Scale factors array cannot be null");
        if (factors.Length != ncols)
            throw new ArgumentException(
                $"Factors length {factors.Length} must match matrix columns {ncols}",
                nameof(factors));

        var localValues = GetValuesOrThrow();
        var newValues = new double[localValues.Length];
        for (var k = 0; k < localValues.Length; k++) newValues[k] = factors[columnIndices[k]] * localValues[k];
        return new CSR((int[])rowPointers.Clone(), (int[])columnIndices.Clone(), ncols, newValues, true);
    }

    public double FrobeniusNorm()
    {
        ThrowIfDisposed();
        var localValues = GetValuesOrThrow();

        if (localValues.Length >= MIN_ROWS_FOR_PARALLEL)
        {
            var lockObj = new object();
            var totalSum = 0.0;
            Parallel.ForEach(
                Partitioner.Create(0, localValues.Length),
                ParallelConfig.Options,
                () => 0.0,
                (range, loopState, partialSum) =>
                {
                    var localSum = partialSum;
                    for (var i = range.Item1; i < range.Item2; i++)
                        localSum += localValues[i] * localValues[i];
                    return localSum;
                },
                partialSum =>
                {
                    lock (lockObj) { totalSum += partialSum; }
                });
            return Math.Sqrt(totalSum);
        }

        var sum = 0.0;
        foreach (var val in localValues) sum += val * val;
        return Math.Sqrt(sum);
    }

    public double InfinityNorm()
    {
        ThrowIfDisposed();
        var localValues = GetValuesOrThrow();
        var maxNorm = 0.0;
        for (var i = 0; i < nrows; i++)
        {
            var rowSum = 0.0;
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++) rowSum += Math.Abs(localValues[k]);
            maxNorm = Math.Max(maxNorm, rowSum);
        }

        return maxNorm;
    }

    public double OneNorm()
    {
        ThrowIfDisposed();
        if (ncols == 0 || nrows == 0) return 0.0;
        var localValues = GetValuesOrThrow();
        var colSums = new double[ncols];
        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            colSums[columnIndices[k]] += Math.Abs(localValues[k]);
        return colSums.Max();
    }

    public double[,] ToDense()
    {
        ThrowIfDisposed();
        var localValues = GetValuesOrThrow();
        var dense = new double[nrows, ncols];
        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            dense[i, columnIndices[k]] = localValues[k];
        return dense;
    }

    /// <summary>True if the matrix is square (Rows == Columns).</summary>
    public bool IsSquare
    {
        get
        {
            ThrowIfDisposed();
            return nrows == ncols;
        }
    }

    /// <summary>
    ///     Density of the matrix: fraction of non-zero entries (nnz / (rows * cols)).
    ///     Note: Both <see cref="Density"/> and <see cref="Sparsity"/> return the same value
    ///     (the non-zero ratio). The <c>Sparsity</c> name is retained for backward compatibility.
    /// </summary>
    public double Density
    {
        get
        {
            ThrowIfDisposed();
            return Sparsity;
        }
    }

    /// <summary>
    ///     Returns the sum of diagonal entries (trace). Requires a square matrix with values.
    ///     Time complexity: O(nnz)
    /// </summary>
    public double Trace()
    {
        ThrowIfDisposed();
        if (nrows != ncols)
            throw new InvalidOperationException("Trace requires a square matrix");

        var localValues = GetValuesOrThrow();
        var trace = 0.0;
        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            if (columnIndices[k] == i)
            {
                trace += localValues[k];
                break;
            }

        return trace;
    }

    /// <summary>
    ///     Extracts the diagonal of the matrix as a vector.
    ///     For non-square matrices, extracts min(Rows, Columns) diagonal entries.
    ///     Missing diagonal entries are returned as 0.0.
    ///     Time complexity: O(nnz)
    /// </summary>
    public double[] DiagonalVector()
    {
        ThrowIfDisposed();
        var localValues = GetValuesOrThrow();
        var diagSize = Math.Min(nrows, ncols);
        var diag = new double[diagSize];

        for (var i = 0; i < diagSize; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            if (columnIndices[k] == i)
            {
                diag[i] = localValues[k];
                break;
            }

        return diag;
    }

    /// <summary>
    ///     Checks if the matrix is structurally symmetric (same sparsity pattern for A and A^T).
    ///     Does not check values — only the pattern.
    ///     Time complexity: O(nnz * log(nnz_per_row)) worst case.
    /// </summary>
    public bool IsStructurallySymmetric()
    {
        ThrowIfDisposed();
        if (nrows != ncols) return false;

        // For each entry (i, j), check that (j, i) also exists
        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
        {
            var j = columnIndices[k];
            // Search for column i in row j
            var found = false;
            for (var kk = rowPointers[j]; kk < rowPointers[j + 1]; kk++)
                if (columnIndices[kk] == i)
                {
                    found = true;
                    break;
                }

            if (!found) return false;
        }

        return true;
    }

    /// <summary>
    ///     Checks if the matrix is numerically symmetric (A[i,j] == A[j,i] within tolerance).
    ///     Requires values array. Time complexity: O(nnz * log(nnz_per_row)) worst case.
    /// </summary>
    public bool IsSymmetric(double tolerance = DEFAULT_TOLERANCE)
    {
        ThrowIfDisposed();
        if (nrows != ncols) return false;

        var localValues = GetValuesOrThrow();

        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
        {
            var j = columnIndices[k];
            var aij = localValues[k];

            // Search for entry (j, i) and compare values
            var found = false;
            for (var kk = rowPointers[j]; kk < rowPointers[j + 1]; kk++)
                if (columnIndices[kk] == i)
                {
                    if (Math.Abs(aij - localValues[kk]) > tolerance)
                        return false;
                    found = true;
                    break;
                }

            if (!found && Math.Abs(aij) > tolerance)
                return false;
        }

        return true;
    }

    public bool Validate(out string? errorMessage)
    {
        ThrowIfDisposed();
        try
        {
            ValidateCSRStructure(rowPointers, columnIndices, nrows, ncols);
            errorMessage = null;
            return true;
        }
        catch (ArgumentException ex)
        {
            errorMessage = ex.Message;
            return false;
        }
    }

    public MatrixStatistics GetStatistics()
    {
        ThrowIfDisposed();
        if (nrows == 0)
            return new MatrixStatistics(0, ncols, 0, 1.0, 0, 0, 0.0);

        var minNnz = int.MaxValue;
        var maxNnz = 0;
        long totalNnz = 0;
        for (var i = 0; i < nrows; i++)
        {
            var nnz = rowPointers[i + 1] - rowPointers[i];
            minNnz = Math.Min(minNnz, nnz);
            maxNnz = Math.Max(maxNnz, nnz);
            totalNnz += nnz;
        }

        return new MatrixStatistics(nrows, ncols, (int)totalNnz, Sparsity, minNnz, maxNnz, totalNnz / (double)nrows);
    }

    public void PrintStructure(int maxRows = 10, int maxCols = 10)
    {
        ThrowIfDisposed();
        Console.WriteLine($"CSR Matrix: {nrows}x{ncols}, nnz={NonZeroCount}");
        double[]? localValues;
        lock (syncLock)
        {
            localValues = values;
        }

        var rowsToShow = Math.Min(maxRows, nrows);
        for (var i = 0; i < rowsToShow; i++)
        {
            Console.Write($"Row {i}: ");
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            {
                var col = columnIndices[k];
                if (col < maxCols) Console.Write(localValues != null ? $"({col},{localValues[k]:F2}) " : $"{col} ");
            }

            Console.WriteLine();
        }

        if (nrows > maxRows) Console.WriteLine($"... ({nrows - maxRows} more rows)");
    }

    #endregion

    #region Matrix Permutation Operations

    /// <summary>
    ///     Permute rows of the matrix according to the given permutation.
    ///     Result: C[p[i],:] = A[i,:]
    ///     Time complexity: O(nnz)
    ///     Works with structural matrices (values can be null).
    /// </summary>
    /// <param name="permutation">Row permutation array where permutation[i] = new_row_index</param>
    /// <returns>Permuted matrix</returns>
    public CSR PermuteRows(int[] permutation)
    {
        ThrowIfDisposed();

        if (permutation == null)
            throw new ArgumentNullException(nameof(permutation));
        if (permutation.Length != nrows)
            throw new ArgumentException(
                $"Permutation length {permutation.Length} must equal number of rows {nrows}",
                nameof(permutation));

        // Validate permutation
        var check = new bool[nrows];
        for (var i = 0; i < nrows; i++)
        {
            if (permutation[i] < 0 || permutation[i] >= nrows)
                throw new ArgumentException(
                    $"Invalid permutation index {permutation[i]} at position {i}",
                    nameof(permutation));
            if (check[permutation[i]])
                throw new ArgumentException(
                    $"Duplicate permutation index {permutation[i]}",
                    nameof(permutation));
            check[permutation[i]] = true;
        }

        // Count non-zeros in each new row
        var newRowCounts = new int[nrows];
        for (var i = 0; i < nrows; i++)
        {
            var newRow = permutation[i];
            newRowCounts[newRow] = rowPointers[i + 1] - rowPointers[i];
        }

        // Build new row pointers
        var newRowPtr = new int[nrows + 1];
        newRowPtr[0] = 0;
        for (var i = 0; i < nrows; i++)
            newRowPtr[i + 1] = newRowPtr[i] + newRowCounts[i];

        var nnz = newRowPtr[nrows];
        var newColIdx = new int[nnz];

        // Handle values only if present
        double[]? localValues;
        lock (syncLock)
        {
            localValues = values;
        }

        var newValues = localValues != null ? new double[nnz] : null;

        // Copy data to new positions
        for (var i = 0; i < nrows; i++)
        {
            var newRow = permutation[i];
            var destPos = newRowPtr[newRow];
            var srcStart = rowPointers[i];
            var srcEnd = rowPointers[i + 1];
            var count = srcEnd - srcStart;

            // Always copy column indices
            Array.Copy(columnIndices, srcStart, newColIdx, destPos, count);

            // Copy values only if both source and destination arrays exist
            if (newValues != null && localValues != null)
                Array.Copy(localValues, srcStart, newValues, destPos, count);
        }

        return new CSR(newRowPtr, newColIdx, ncols, newValues, true);
    }

    /// <summary>
    ///     Permute columns of the matrix according to the given permutation.
    ///     Result: C[:,p[j]] = A[:,j]
    ///     Time complexity: O(nnz)
    ///     Works with structural matrices (values can be null).
    /// </summary>
    /// <param name="permutation">Column permutation array</param>
    /// <returns>Permuted matrix</returns>
    public CSR PermuteColumns(int[] permutation)
    {
        ThrowIfDisposed();

        if (permutation == null)
            throw new ArgumentNullException(nameof(permutation));
        if (permutation.Length != ncols)
            throw new ArgumentException(
                $"Permutation length {permutation.Length} must equal number of columns {ncols}",
                nameof(permutation));

        // Validate permutation
        var check = new bool[ncols];
        for (var i = 0; i < ncols; i++)
        {
            if (permutation[i] < 0 || permutation[i] >= ncols)
                throw new ArgumentException(
                    $"Invalid permutation index {permutation[i]} at position {i}",
                    nameof(permutation));
            if (check[permutation[i]])
                throw new ArgumentException(
                    $"Duplicate permutation index {permutation[i]}",
                    nameof(permutation));
            check[permutation[i]] = true;
        }

        // Column permutation: remap column indices
        var newRowPtr = (int[])rowPointers.Clone();
        var newColIdx = new int[columnIndices.Length];

        for (var i = 0; i < columnIndices.Length; i++)
            newColIdx[i] = permutation[columnIndices[i]];

        // Handle values only if present
        double[]? newValues = null;
        lock (syncLock)
        {
            if (values != null)
                newValues = (double[])values.Clone();
        }

        return new CSR(newRowPtr, newColIdx, ncols, newValues, true);
    }

    /// <summary>
    ///     Apply symmetric permutation: C = P·A·P^T where P is the permutation matrix.
    ///     This is more efficient than doing row and column permutations separately.
    ///     Time complexity: O(nnz)
    ///     Works with structural matrices (values can be null).
    /// </summary>
    /// <param name="permutation">Symmetric permutation array</param>
    /// <returns>Symmetrically permuted matrix</returns>
    public CSR PermuteSymmetric(int[] permutation)
    {
        ThrowIfDisposed();

        if (permutation == null)
            throw new ArgumentNullException(nameof(permutation));
        if (permutation.Length != nrows)
            throw new ArgumentException(
                $"Permutation length {permutation.Length} must equal matrix dimension {nrows}",
                nameof(permutation));
        if (nrows != ncols)
            throw new InvalidOperationException(
                "Symmetric permutation requires a square matrix");

        // Validate permutation
        var check = new bool[nrows];
        for (var i = 0; i < nrows; i++)
        {
            if (permutation[i] < 0 || permutation[i] >= nrows)
                throw new ArgumentException(
                    $"Invalid permutation index {permutation[i]} at position {i}",
                    nameof(permutation));
            if (check[permutation[i]])
                throw new ArgumentException(
                    $"Duplicate permutation index {permutation[i]}",
                    nameof(permutation));
            check[permutation[i]] = true;
        }

        // First permute rows (C = P·A)
        var tempMatrix = PermuteRows(permutation);

        // Then permute columns (C = (P·A)·P^T)
        var result = tempMatrix.PermuteColumns(permutation);

        return result;
    }

    #endregion

    #region Format Conversion Methods

    /// <summary>
    ///     Convert CSR format to Coordinate (COO) format.
    ///     Returns (row_indices, column_indices, values).
    ///     For structural matrices (no values), returns empty values array.
    ///     Time complexity: O(nnz)
    ///     Space complexity: O(nnz)
    /// </summary>
    /// <returns>Tuple of (rows, cols, values) arrays in COO format</returns>
    public (int[] rows, int[] cols, double[] values) ToCOO()
    {
        ThrowIfDisposed();

        var nnz = rowPointers[nrows];
        var rows = new int[nnz];
        var cols = new int[nnz];

        double[]? localValues;
        lock (syncLock)
        {
            localValues = values;
        }

        // Return empty array if no values (structural matrix)
        var vals = localValues != null ? new double[nnz] : Array.Empty<double>();

        var pos = 0;
        for (var i = 0; i < nrows; i++)
        for (var j = rowPointers[i]; j < rowPointers[i + 1]; j++)
        {
            rows[pos] = i;
            cols[pos] = columnIndices[j];
            if (localValues != null)
                vals[pos] = localValues[j];
            pos++;
        }

        return (rows, cols, vals);
    }

    /// <summary>
    ///     Build a CSR matrix from COO (coordinate) arrays.
    ///     Duplicate coordinates are summed when <paramref name="sumDuplicates" /> is true.
    /// </summary>
    /// <param name="rows">Row indices for each entry.</param>
    /// <param name="cols">Column indices for each entry.</param>
    /// <param name="nRows">Total matrix row count.</param>
    /// <param name="nCols">Total matrix column count.</param>
    /// <param name="vals">Optional values; when null, a structural matrix is created.</param>
    /// <param name="sumDuplicates">If true, duplicate coordinates are merged by summing their values.</param>
    /// <returns>New CSR matrix.</returns>
    public static CSR FromCOO(
        int[] rows,
        int[] cols,
        int nRows,
        int nCols,
        double[]? vals = null,
        bool sumDuplicates = true)
    {
        ArgumentNullException.ThrowIfNull(rows);
        ArgumentNullException.ThrowIfNull(cols);

        if (nRows < 0)
            throw new ArgumentException("Row count must be non-negative", nameof(nRows));
        if (nCols < 0)
            throw new ArgumentException("Column count must be non-negative", nameof(nCols));
        if (rows.Length != cols.Length)
            throw new ArgumentException("Rows and columns must have the same length");
        if (vals != null && vals.Length != rows.Length)
            throw new ArgumentException("Values length must match rows/cols length", nameof(vals));

        if (vals != null)
        {
            var perRow = new Dictionary<int, double>[nRows];

            for (var k = 0; k < rows.Length; k++)
            {
                var row = rows[k];
                var col = cols[k];
                if (row < 0 || row >= nRows)
                    throw new ArgumentOutOfRangeException(nameof(rows),
                        $"Row index {row} at position {k} is out of range [0, {nRows})");
                if (col < 0 || col >= nCols)
                    throw new ArgumentOutOfRangeException(nameof(cols),
                        $"Column index {col} at position {k} is out of range [0, {nCols})");

                perRow[row] ??= new Dictionary<int, double>();
                var rowMap = perRow[row]!;

                if (rowMap.TryGetValue(col, out var current))
                {
                    if (!sumDuplicates)
                        throw new ArgumentException(
                            $"Duplicate coordinate ({row}, {col}) found at position {k}. Set {nameof(sumDuplicates)}=true to merge duplicates.");
                    rowMap[col] = current + vals[k];
                }
                else
                {
                    rowMap[col] = vals[k];
                }
            }

            var rowPtr = new int[nRows + 1];
            for (var i = 0; i < nRows; i++)
                rowPtr[i + 1] = rowPtr[i] + (perRow[i]?.Count ?? 0);

            var nnz = rowPtr[nRows];
            var colIdx = new int[nnz];
            var valuesOut = new double[nnz];

            for (var i = 0; i < nRows; i++)
            {
                var rowMap = perRow[i];
                if (rowMap == null || rowMap.Count == 0)
                    continue;

                var offset = rowPtr[i];
                var count = rowMap.Count;

                // Copy keys and values into output arrays, then sort in-place
                var p = 0;
                foreach (var (col, value) in rowMap)
                {
                    colIdx[offset + p] = col;
                    valuesOut[offset + p] = value;
                    p++;
                }

                // Sort columns (and corresponding values) in-place by column index
                Array.Sort(colIdx, valuesOut, offset, count);
            }

            return new CSR(rowPtr, colIdx, nCols, valuesOut, true);
        }
        else
        {
            var perRow = new HashSet<int>[nRows];
            for (var k = 0; k < rows.Length; k++)
            {
                var row = rows[k];
                var col = cols[k];
                if (row < 0 || row >= nRows)
                    throw new ArgumentOutOfRangeException(nameof(rows),
                        $"Row index {row} at position {k} is out of range [0, {nRows})");
                if (col < 0 || col >= nCols)
                    throw new ArgumentOutOfRangeException(nameof(cols),
                        $"Column index {col} at position {k} is out of range [0, {nCols})");

                perRow[row] ??= new HashSet<int>();
                if (!perRow[row]!.Add(col) && !sumDuplicates)
                    throw new ArgumentException(
                        $"Duplicate coordinate ({row}, {col}) found at position {k}. Set {nameof(sumDuplicates)}=true to merge duplicates.");
            }

            var rowPtr = new int[nRows + 1];
            for (var i = 0; i < nRows; i++)
                rowPtr[i + 1] = rowPtr[i] + (perRow[i]?.Count ?? 0);

            var nnz = rowPtr[nRows];
            var colIdx = new int[nnz];

            for (var i = 0; i < nRows; i++)
            {
                var rowSet = perRow[i];
                if (rowSet == null || rowSet.Count == 0)
                    continue;

                var offset = rowPtr[i];
                var p = 0;
                foreach (var col in rowSet)
                    colIdx[offset + p++] = col;

                // Sort column indices in-place
                Array.Sort(colIdx, offset, rowSet.Count);
            }

            return new CSR(rowPtr, colIdx, nCols, null, true);
        }
    }

    /// <summary>
    ///     Export matrix to Matrix Market coordinate format string.
    ///     Follows the Matrix Market format specification.
    ///     For structural matrices (no values), exports pattern format.
    /// </summary>
    /// <param name="symmetric">If true, only stores upper triangular part for symmetric matrices</param>
    /// <param name="comment">Optional comment line to include in header</param>
    /// <returns>Matrix Market format string</returns>
    public string ToMatrixMarket(bool symmetric = false, string? comment = null)
    {
        ThrowIfDisposed();

        var (rows, cols, vals) = ToCOO();
        var sb = new StringBuilder();

        // Check if this is a structural matrix (no values)
        var isPattern = vals.Length == 0;

        // Header line
        if (isPattern)
            sb.AppendLine(symmetric
                ? "%%MatrixMarket matrix coordinate pattern symmetric"
                : "%%MatrixMarket matrix coordinate pattern general");
        else
            sb.AppendLine(symmetric
                ? "%%MatrixMarket matrix coordinate real symmetric"
                : "%%MatrixMarket matrix coordinate real general");

        // Optional comment
        if (!string.IsNullOrEmpty(comment))
            sb.AppendLine($"% {comment}");

        // Count entries to write
        int nnzToWrite;
        if (isPattern)
        {
            // For pattern matrices, count structural entries
            nnzToWrite = rowPointers[nrows] - rowPointers[0];
            if (symmetric)
            {
                nnzToWrite = 0;
                for (var i = 0; i < nrows; i++)
                for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
                    if (columnIndices[k] >= i)
                        nnzToWrite++;
            }
        }
        else
        {
            // For valued matrices
            nnzToWrite = vals.Length;
            if (symmetric)
            {
                // Count only upper triangular entries
                nnzToWrite = 0;
                for (var i = 0; i < vals.Length; i++)
                    if (cols[i] >= rows[i])
                        nnzToWrite++;
            }
        }

        // Dimensions line (Matrix Market uses 1-based indexing)
        sb.AppendLine($"{nrows} {ncols} {nnzToWrite}");

        // Data lines
        if (isPattern)
            // Pattern format: only row and column indices
            for (var i = 0; i < nrows; i++)
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            {
                var col = columnIndices[k];
                if (symmetric && col < i)
                    continue; // Skip lower triangle for symmetric

                // Matrix Market uses 1-based indexing
                sb.AppendLine($"{i + 1} {col + 1}");
            }
        else
            // Real format: row, column, and value
            for (var i = 0; i < vals.Length; i++)
            {
                if (symmetric && cols[i] < rows[i])
                    continue; // Skip lower triangle for symmetric

                // Matrix Market uses 1-based indexing
                sb.AppendLine($"{rows[i] + 1} {cols[i] + 1} {vals[i]:E15}");
            }

        return sb.ToString();
    }

    #endregion

    #region Triangular Solve Methods

    /// <summary>
    ///     Solve L·x = b where L is lower triangular (CSR format).
    ///     Handles unsorted column indices by scanning the entire row.
    ///     REQUIRES values array to be present (throws if structural matrix).
    ///     Time complexity: O(nnz)
    /// </summary>
    /// <param name="b">Right-hand side vector</param>
    /// <param name="unitDiagonal">If true, assumes diagonal entries are 1.0</param>
    /// <returns>Solution vector x</returns>
    public double[] SolveLowerTriangular(double[] b, bool unitDiagonal = false)
    {
        ThrowIfDisposed();

        if (b == null)
            throw new ArgumentNullException(nameof(b));
        if (b.Length != nrows)
            throw new ArgumentException(
                $"RHS length {b.Length} must equal rows {nrows}",
                nameof(b));
        if (nrows != ncols)
            throw new InvalidOperationException("Triangular solve requires square matrix");

        // Triangular solve requires values
        var localValues = GetValuesOrThrow();

        var x = new double[nrows];
        Array.Copy(b, x, nrows);

        for (var i = 0; i < nrows; i++)
        {
            var rowStart = rowPointers[i];
            var rowEnd = rowPointers[i + 1];

            // Subtract contributions from previously computed x[j] where j < i
            double diag = 0;
            var foundDiag = false;

            for (var k = rowStart; k < rowEnd; k++)
            {
                var col = columnIndices[k];
                if (col < i)
                {
                    // Lower triangular contribution
                    x[i] -= localValues[k] * x[col];
                }
                else if (col == i)
                {
                    // Found diagonal
                    diag = localValues[k];
                    foundDiag = true;
                }
                // col > i: upper triangular part, ignore for lower triangular solve
            }

            // Divide by diagonal
            if (!unitDiagonal)
            {
                if (!foundDiag)
                    throw new InvalidOperationException($"Missing diagonal entry at row {i}");
                if (Math.Abs(diag) < DEFAULT_TOLERANCE)
                    throw new InvalidOperationException(
                        $"Zero or near-zero diagonal at row {i}: {diag}");

                x[i] /= diag;
            }
        }

        return x;
    }

    /// <summary>
    ///     Solve U·x = b where U is upper triangular (CSR format).
    ///     Handles unsorted column indices by scanning the entire row.
    ///     REQUIRES values array to be present (throws if structural matrix).
    ///     Time complexity: O(nnz)
    /// </summary>
    /// <param name="b">Right-hand side vector</param>
    /// <param name="unitDiagonal">If true, assumes diagonal entries are 1.0</param>
    /// <returns>Solution vector x</returns>
    public double[] SolveUpperTriangular(double[] b, bool unitDiagonal = false)
    {
        ThrowIfDisposed();

        if (b == null)
            throw new ArgumentNullException(nameof(b));
        if (b.Length != nrows)
            throw new ArgumentException(
                $"RHS length {b.Length} must equal rows {nrows}",
                nameof(b));
        if (nrows != ncols)
            throw new InvalidOperationException("Triangular solve requires square matrix");

        // Triangular solve requires values
        var localValues = GetValuesOrThrow();

        var x = (double[])b.Clone();

        for (var i = nrows - 1; i >= 0; i--)
        {
            var rowStart = rowPointers[i];
            var rowEnd = rowPointers[i + 1];

            // Subtract contributions from previously computed x[j] where j > i
            double diag = 0;
            var foundDiag = false;

            for (var k = rowStart; k < rowEnd; k++)
            {
                var col = columnIndices[k];
                if (col == i)
                {
                    // Found diagonal
                    diag = localValues[k];
                    foundDiag = true;
                }
                else if (col > i)
                {
                    // Upper triangular contribution
                    x[i] -= localValues[k] * x[col];
                }
                // col < i: lower triangular part, ignore for upper triangular solve
            }

            // Divide by diagonal
            if (!unitDiagonal)
            {
                if (!foundDiag)
                    throw new InvalidOperationException($"Missing diagonal entry at row {i}");
                if (Math.Abs(diag) < DEFAULT_TOLERANCE)
                    throw new InvalidOperationException(
                        $"Zero or near-zero diagonal at row {i}: {diag}");

                x[i] /= diag;
            }
        }

        return x;
    }

    /// <summary>
    ///     Returns a matrix with entries |a_ij| ≤ tolerance removed.
    ///     Structural matrices are returned unchanged.
    /// </summary>
    /// <param name="tolerance">Non-negative pruning tolerance.</param>
    /// <returns>New matrix with small values removed.</returns>
    public CSR DropZeros(double tolerance = DEFAULT_TOLERANCE)
    {
        ThrowIfDisposed();

        if (tolerance < 0)
            throw new ArgumentException("Tolerance must be non-negative", nameof(tolerance));

        double[]? localValues;
        lock (syncLock)
        {
            localValues = values;
        }

        if (localValues == null)
            return (CSR)Clone();

        var rowPtr = new int[nrows + 1];
        var nnz = 0;

        for (var i = 0; i < nrows; i++)
        {
            for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
                if (Math.Abs(localValues[k]) > tolerance)
                    nnz++;

            rowPtr[i + 1] = nnz;
        }

        if (nnz == localValues.Length)
            return (CSR)Clone();

        var colIdx = new int[nnz];
        var outValues = new double[nnz];
        var pos = 0;

        for (var i = 0; i < nrows; i++)
        for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++)
            if (Math.Abs(localValues[k]) > tolerance)
            {
                colIdx[pos] = columnIndices[k];
                outValues[pos] = localValues[k];
                pos++;
            }

        return new CSR(rowPtr, colIdx, ncols, outValues, true);
    }

    #endregion

    #region Factory Methods

    public static CSR Identity(int n)
    {
        if (n < 0)
            throw new ArgumentOutOfRangeException(nameof(n), "Identity matrix size must be non-negative");

        var rowPtr = new int[n + 1];
        var colIdx = new int[n];
        var vals = new double[n];
        for (var i = 0; i <= n; i++) rowPtr[i] = i;
        for (var i = 0; i < n; i++)
        {
            colIdx[i] = i;
            vals[i] = 1.0;
        }

        return new CSR(rowPtr, colIdx, n, vals, true);
    }

    public static CSR Diagonal(double[] diagonal)
    {
        if (diagonal == null)
            throw new ArgumentNullException(nameof(diagonal), "Diagonal array cannot be null");

        var n = diagonal.Length;
        var rowPtr = new int[n + 1];
        var colIdx = new int[n];
        var vals = new double[n];
        for (var i = 0; i <= n; i++) rowPtr[i] = i;
        for (var i = 0; i < n; i++)
        {
            colIdx[i] = i;
            vals[i] = diagonal[i];
        }

        return new CSR(rowPtr, colIdx, n, vals, true);
    }

    public static CSR Random(int rows, int cols, double sparsity = 0.1, int seed = 42, bool symmetric = false)
    {
        if (rows <= 0)
            throw new ArgumentException("Rows must be positive", nameof(rows));
        if (cols <= 0)
            throw new ArgumentException("Columns must be positive", nameof(cols));
        if (sparsity < 0 || sparsity > 1)
            throw new ArgumentException("Sparsity must be between 0 and 1", nameof(sparsity));
        if (symmetric && rows != cols)
            throw new ArgumentException("Symmetric matrix requires rows == cols", nameof(symmetric));

        var random = new Random(seed);
        var nnzPerRow = Math.Max(1, (int)(cols * sparsity));

        // First pass: build row structure and count total nnz
        var rowPointers = new int[rows + 1];
        var rowColLists = new int[rows][];
        var rowValLists = new double[rows][];

        for (var i = 0; i < rows; i++)
        {
            var colsSet = new HashSet<int>();
            var rowCols = new List<int>(nnzPerRow);
            var rowVals = new List<double>(nnzPerRow);

            var maxAttempts = nnzPerRow * 10;
            var attempts = 0;
            while (colsSet.Count < nnzPerRow && attempts < maxAttempts)
            {
                var col = random.Next(cols);
                if (colsSet.Add(col))
                {
                    rowCols.Add(col);
                    rowVals.Add(random.NextDouble());
                }

                attempts++;
            }

            rowColLists[i] = rowCols.ToArray();
            rowValLists[i] = rowVals.ToArray();
            rowPointers[i + 1] = rowPointers[i] + rowColLists[i].Length;
        }

        // Build flat arrays directly
        var totalNnz = rowPointers[rows];
        var columnIndices = new int[totalNnz];
        var values = new double[totalNnz];

        for (var i = 0; i < rows; i++)
        {
            var start = rowPointers[i];
            Array.Copy(rowColLists[i], 0, columnIndices, start, rowColLists[i].Length);
            Array.Copy(rowValLists[i], 0, values, start, rowValLists[i].Length);
        }

        var ncols = rows > 0 && totalNnz > 0 ? columnIndices.Max() + 1 : 0;
        ncols = Math.Max(ncols, cols);
        var result = new CSR(rowPointers, columnIndices, ncols, values, true);

        if (symmetric && rows == cols) result = Add(result, result.Transpose(), 0.5, 0.5);
        return result;
    }

    /// <summary>
    ///     Creates a CSR sparsity pattern from mesh topology (node-to-node adjacency).
    /// </summary>
    /// <typeparam name="TTypes">TypeMap containing element and node types.</typeparam>
    /// <typeparam name="TElement">Element type (e.g., Tri3, Tet4).</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="topology">Source mesh topology.</param>
    /// <param name="dofsPerNode">DOFs per node (1 for scalar, 2-3 for vector problems).</param>
    /// <param name="includeDiagonal">Include diagonal entries (recommended).</param>
    /// <returns>CSR matrix with sparsity pattern (values array is null).</returns>
    /// <remarks>
    ///     <para>
    ///         <b>USE CASE:</b> Building sparsity patterns for FEM stiffness matrices.
    ///         Two nodes are adjacent if they belong to the same element.
    ///     </para>
    ///     <para>
    ///         <b>DOF NUMBERING:</b> For dofsPerNode > 1, DOFs are numbered consecutively:
    ///         Node 0 → DOFs [0, dofsPerNode-1], Node 1 → DOFs [dofsPerNode, 2*dofsPerNode-1], etc.
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Scalar problem (heat conduction)
    ///         var pattern = CSR.FromTopology&lt;MyMesh, Tri3, Node&gt;(mesh);
    ///         
    ///         // 3D elasticity (3 DOFs per node)
    ///         var pattern = CSR.FromTopology&lt;MyMesh, Tet4, Node&gt;(mesh, dofsPerNode: 3);
    ///         </code>
    ///     </example>
    /// </remarks>
    public static CSR FromTopology<TTypes, TElement, TNode>(
        Topology<TTypes> topology,
        int dofsPerNode = 1,
        bool includeDiagonal = true)
        where TTypes : ITypeMap, new()
        where TNode : struct
        where TElement : struct
    {
        if (topology == null)
            throw new ArgumentNullException(nameof(topology));
        if (dofsPerNode <= 0)
            throw new ArgumentOutOfRangeException(nameof(dofsPerNode), "Must be positive");

        var nodeCount = topology.Count<TNode>();
        var elementCount = topology.Count<TElement>();

        if (nodeCount == 0)
            return new CSR(new int[1], Array.Empty<int>(), 0);

        // Build node-level adjacency using arrays instead of HashSet for better performance
        // Use a marker array to track which nodes are already added to each row
        var marker = ArrayPool<int>.Shared.Rent(nodeCount);
        try
        {
            // Initialize marker to -1 (not visited)
            for (var i = 0; i < nodeCount; i++)
                marker[i] = -1;

            // First pass: count non-zeros per node row
            var nodeCounts = new int[nodeCount];

            for (var e = 0; e < elementCount; e++)
            {
                var nodes = topology.NodesOf<TElement, TNode>(e);
                var n = nodes.Count;

                for (var i = 0; i < n; i++)
                {
                    var ni = nodes[i];

                    // Add diagonal if needed
                    if (includeDiagonal && marker[ni] != ni)
                    {
                        marker[ni] = ni;
                        nodeCounts[ni]++;
                    }

                    // Add all other nodes in this element
                    for (var j = 0; j < n; j++)
                    {
                        if (i == j) continue;
                        var nj = nodes[j];
                        if (marker[nj] != ni)
                        {
                            marker[nj] = ni;
                            nodeCounts[ni]++;
                        }
                    }
                }
            }

            // Reset marker
            for (var i = 0; i < nodeCount; i++)
                marker[i] = -1;

            // Build row pointers for node-level adjacency
            var nodeRowPtrs = new int[nodeCount + 1];
            for (var i = 0; i < nodeCount; i++)
                nodeRowPtrs[i + 1] = nodeRowPtrs[i] + nodeCounts[i];

            var nodeNnz = nodeRowPtrs[nodeCount];
            var nodeColIdx = new int[nodeNnz];
            var rowPos = new int[nodeCount];
            for (var i = 0; i < nodeCount; i++)
                rowPos[i] = nodeRowPtrs[i];

            // Second pass: fill column indices
            for (var e = 0; e < elementCount; e++)
            {
                var nodes = topology.NodesOf<TElement, TNode>(e);
                var n = nodes.Count;

                for (var i = 0; i < n; i++)
                {
                    var ni = nodes[i];

                    if (includeDiagonal && marker[ni] != ni)
                    {
                        marker[ni] = ni;
                        nodeColIdx[rowPos[ni]++] = ni;
                    }

                    for (var j = 0; j < n; j++)
                    {
                        if (i == j) continue;
                        var nj = nodes[j];
                        if (marker[nj] != ni)
                        {
                            marker[nj] = ni;
                            nodeColIdx[rowPos[ni]++] = nj;
                        }
                    }
                }
            }

            // Sort column indices within each row (required for many sparse operations)
            for (var i = 0; i < nodeCount; i++)
            {
                var start = nodeRowPtrs[i];
                var end = nodeRowPtrs[i + 1];
                Array.Sort(nodeColIdx, start, end - start);
            }

            // If single DOF per node, we're done
            if (dofsPerNode == 1)
                return new CSR(nodeRowPtrs, nodeColIdx, nodeCount, null, true);

            // Expand to DOF level
            var totalDofsLong = (long)nodeCount * dofsPerNode;
            if (totalDofsLong > int.MaxValue)
                throw new OverflowException(
                    $"DOF expansion overflows: {nodeCount} nodes * {dofsPerNode} DOFs/node = {totalDofsLong:N0}, " +
                    "which exceeds int.MaxValue. Consider reducing mesh size or using a partitioned approach.");
            var totalDofs = (int)totalDofsLong;
            var dofRowPtrs = new int[totalDofs + 1];

            // Each node-node connection becomes a dofsPerNode × dofsPerNode block
            for (var node = 0; node < nodeCount; node++)
            {
                var nnzPerNodeRow = nodeRowPtrs[node + 1] - nodeRowPtrs[node];
                var nnzPerDofRow = nnzPerNodeRow * dofsPerNode;

                for (var d = 0; d < dofsPerNode; d++)
                {
                    var dofRow = node * dofsPerNode + d;
                    dofRowPtrs[dofRow + 1] = nnzPerDofRow;
                }
            }

            // Prefix sum
            for (var i = 1; i <= totalDofs; i++)
                dofRowPtrs[i] += dofRowPtrs[i - 1];

            var totalNnz = dofRowPtrs[totalDofs];
            var dofColIdx = new int[totalNnz];

            // Fill DOF-level column indices
            for (var node = 0; node < nodeCount; node++)
            {
                var nodeStart = nodeRowPtrs[node];
                var nodeEnd = nodeRowPtrs[node + 1];

                for (var d = 0; d < dofsPerNode; d++)
                {
                    var dofRow = node * dofsPerNode + d;
                    var pos = dofRowPtrs[dofRow];

                    for (var k = nodeStart; k < nodeEnd; k++)
                    {
                        var adjNode = nodeColIdx[k];
                        for (var dd = 0; dd < dofsPerNode; dd++)
                            dofColIdx[pos++] = adjNode * dofsPerNode + dd;
                    }
                }
            }

            return new CSR(dofRowPtrs, dofColIdx, totalDofs, null, true);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(marker);
        }
    }

    /// <summary>
    ///     Creates a CSR matrix representing the element-to-element adjacency (dual graph).
    /// </summary>
    /// <typeparam name="TTypes">TypeMap containing element and node types.</typeparam>
    /// <typeparam name="TElement">Element type.</typeparam>
    /// <typeparam name="TNode">Node type.</typeparam>
    /// <param name="topology">Source mesh topology.</param>
    /// <param name="minSharedNodes">Minimum shared nodes for adjacency (1=vertex, 2=edge, 3=face).</param>
    /// <param name="includeDiagonal">Include diagonal entries.</param>
    /// <returns>CSR matrix representing element adjacency.</returns>
    /// <remarks>
    ///     <para>
    ///         <b>USE CASE:</b> Domain decomposition, mesh partitioning, element coloring.
    ///     </para>
    ///     <example>
    ///         <code>
    ///         // Face neighbors only (for tet mesh)
    ///         var dual = CSR.ElementAdjacencyFromTopology&lt;MyMesh, Tet4, Node&gt;(mesh, minSharedNodes: 3);
    ///         </code>
    ///     </example>
    /// </remarks>
    public static CSR ElementAdjacencyFromTopology<TTypes, TElement, TNode>(
        Topology<TTypes> topology,
        int minSharedNodes = 1,
        bool includeDiagonal = true)
        where TTypes : ITypeMap, new()
        where TNode : struct
        where TElement : struct
    {
        if (topology == null)
            throw new ArgumentNullException(nameof(topology));
        if (minSharedNodes <= 0)
            throw new ArgumentOutOfRangeException(nameof(minSharedNodes), "Must be positive");

        var dual = topology.BuildDualGraph<TElement, TNode>(minSharedNodes);

        var n = dual.ElementCount;
        if (n == 0)
            return new CSR(new int[1], Array.Empty<int>(), 0);

        // Count entries per row
        var rowPointers = new int[n + 1];
        for (var i = 0; i < n; i++)
        {
            var count = dual.Adjacency[i].Count;
            if (includeDiagonal) count++;
            rowPointers[i + 1] = count;
        }

        // Prefix sum
        for (var i = 1; i <= n; i++)
            rowPointers[i] += rowPointers[i - 1];

        var nnz = rowPointers[n];
        var columnIndices = new int[nnz];

        // Fill column indices
        var pos = 0;
        for (var i = 0; i < n; i++)
        {
            var neighbors = dual.Adjacency[i];

            if (includeDiagonal)
            {
                // Insert diagonal in sorted position
                var diagonalInserted = false;
                for (var j = 0; j < neighbors.Count; j++)
                {
                    if (!diagonalInserted && i < neighbors[j])
                    {
                        columnIndices[pos++] = i;
                        diagonalInserted = true;
                    }

                    columnIndices[pos++] = neighbors[j];
                }

                if (!diagonalInserted)
                    columnIndices[pos++] = i;
            }
            else
            {
                for (var j = 0; j < neighbors.Count; j++)
                    columnIndices[pos++] = neighbors[j];
            }
        }

        return new CSR(rowPointers, columnIndices, n, null, true);
    }

    #endregion

    #region IEquatable/IFormattable Implementation

    /// <summary>
    ///     Determines whether this CSR matrix equals another CSR matrix.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         THREAD SAFETY NOTE: This method takes snapshots of the values arrays under individual locks
    ///         to avoid deadlock. However, this means Equals is NOT atomic - if another thread modifies
    ///         either matrix's values between the two snapshots, the comparison may be inconsistent.
    ///     </para>
    ///     <para>
    ///         For thread-safe equality comparison under concurrent modification, callers should use
    ///         external synchronization (e.g., lock both matrices before calling Equals).
    ///     </para>
    ///     <para>
    ///         Comparison order: dimensions → row pointers → column indices → values (if present).
    ///         Returns true only if all components match exactly.
    ///     </para>
    /// </remarks>
    public bool Equals(CSR? other)
    {
        ThrowIfDisposed();
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        other.ThrowIfDisposed();

        if (nrows != other.nrows || ncols != other.ncols) return false;
        if (!rowPointers.SequenceEqual(other.rowPointers)) return false;
        if (!columnIndices.SequenceEqual(other.columnIndices)) return false;

        // Take snapshots under individual locks to avoid nested locking deadlock
        // This is safe because values arrays are effectively immutable once assigned
        // (we replace the entire array rather than modifying elements)
        double[]? thisValues, otherValues;

        lock (syncLock)
        {
            thisValues = values;
        }

        lock (other.syncLock)
        {
            otherValues = other.values;
        }

        if (thisValues == null && otherValues == null) return true;
        if (thisValues == null || otherValues == null) return false;
        return thisValues.SequenceEqual(otherValues);
    }

    public override bool Equals(object? obj)
    {
        return Equals(obj as CSR);
    }

    /// <summary>
    ///     Returns a hash code based on dimensions and structural data.
    ///     Uses only structural invariants (dimensions, nnz, row pointers) for fast computation.
    ///     Note: Equals uses exact value comparison via SequenceEqual.
    /// </summary>
    public override int GetHashCode()
    {
        ThrowIfDisposed();
        var hash = new HashCode();
        hash.Add(nrows);
        hash.Add(ncols);
        hash.Add(columnIndices.Length); // NonZeroCount

        // Include a sample of row pointers for better distribution
        // without being too expensive
        var step = Math.Max(1, nrows / 8);
        for (var i = 0; i < nrows; i += step)
            hash.Add(rowPointers[i + 1] - rowPointers[i]);

        return hash.ToHashCode();
    }

    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        ThrowIfDisposed();
        return $"CSR Matrix [{nrows}x{ncols}], nnz={NonZeroCount}";
    }

    public override string ToString()
    {
        return ToString(null, null);
    }

    public object Clone()
    {
        ThrowIfDisposed();
        lock (syncLock)
        {
            var newValues = values != null ? (double[])values.Clone() : null;
            // Preserve skipValidation state - if original was validated, clone is also valid
            return new CSR((int[])rowPointers.Clone(), (int[])columnIndices.Clone(), ncols, newValues,
                true); // Safe since we're cloning validated data
        }
    }

    #endregion

    #region Dispose Pattern

    public void Dispose()
    {
        Dispose(true);
        GC.SuppressFinalize(this);
    }

    private void Dispose(bool disposing)
    {
        // Quick check without lock
        if (isDisposed) return;

        lock (disposeLock)
        {
            // Double-check inside lock
            if (isDisposed) return;

            // CRITICAL FIX: Set disposed flag FIRST to prevent new operations
            // This allows us to safely clear state without syncLock (preventing deadlock)
            isDisposed = true;

            if (disposing)
            {
                // GPU cleanup doesn't need locks
                var gpu = gpuAccelerator;
                gpuAccelerator = null;
                isGpuInitialized = false;

                try
                {
                    gpu?.Dispose();
                }
                catch (Exception ex)
                {
                    Debug.WriteLine($"GPU cleanup failed: {ex.Message}");
                }

                // CRITICAL FIX: Clear cache and values without syncLock
                // This is safe because isDisposed=true prevents concurrent access
                // All public methods call ThrowIfDisposed() which checks this flag
                cachedTranspose = null;
                values = null;
            }
        }
    }

    ~CSR()
    {
        Dispose(false);
    }

    #endregion
}

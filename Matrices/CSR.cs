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
        : this(rowPointers, columnIndices, columnIndices.Length > 0 ? columnIndices.Max() + 1 : 0, values)
    {
    }

    public CSR(List<List<int>> rows, bool sorted = false, bool enableGpu = false)
    {
        ArgumentNullException.ThrowIfNull(rows);

        for (var i = 0; i < rows.Count; i++)
            if (rows[i] == null)
                throw new ArgumentException($"Row {i} is null. Rows collection must not contain null entries.",
                    nameof(rows));

        nrows = rows.Count;
        ncols = rows.Count > 0 ? rows.Max(r => r.Count > 0 ? r.Max() + 1 : 0) : 0;

        rowPointers = new int[nrows + 1];
        var nnz = 0;

        for (var i = 0; i < nrows; i++)
        {
            rowPointers[i] = nnz;
            // NOTE: Sorting is no longer required - columns can be in any order
            // The 'sorted' parameter is kept for API compatibility but ignored
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

        // Keep this constructor's validation behavior aligned with the primary constructor.
        // This is required for safety because SIMD code paths use unsafe pointers.
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

        var numThreads = ParallelConfig.MaxDegreeOfParallelism;
        var threadLocalResults = new double[numThreads][];
        Parallel.For(0, numThreads, ParallelConfig.Options, tid =>
        {
            threadLocalResults[tid] = new double[resultNnz];
            var localResult = threadLocalResults[tid];
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
                        localResult[pos] += aVal * bValues[kb];
                    }
                }
            }
        });

        if (resultNnz < SMALL_ARRAY_MERGE_THRESHOLD)
            for (var pos = 0; pos < resultNnz; pos++)
            {
                var sum = 0.0;
                for (var tid = 0; tid < numThreads; tid++) sum += threadLocalResults[tid][pos];
                resultValues[pos] = sum;
            }
        else
            Parallel.For(0, numThreads, ParallelConfig.Options, tid =>
            {
                var chunkSize = (resultNnz + numThreads - 1) / numThreads;
                var start = tid * chunkSize;
                var end = Math.Min(start + chunkSize, resultNnz);
                for (var pos = start; pos < end; pos++)
                {
                    var sum = 0.0;
                    for (var t = 0; t < numThreads; t++) sum += threadLocalResults[t][pos];
                    resultValues[pos] = sum;
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
            // Parallel: each row writes to its own non-overlapping segment in values
            Parallel.For(0, A.nrows, ParallelConfig.Options,
                () => ArrayPool<int>.Shared.Rent(A.ncols),
                (i, loopState, workspace) =>
                {
                    // Build position map for this row's output segment
                    var pos = rowPtr[i];
                    for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    {
                        var col = A.columnIndices[ka];
                        if (workspace[col] != i + 1)
                        {
                            values[pos] = alpha * aValues[ka];
                            workspace[col] = i + 1;
                            pos++;
                        }
                        else
                        {
                            // Find position of this column in the output segment
                            for (var k = rowPtr[i]; k < pos; k++)
                                if (colIdx[k] == col)
                                {
                                    values[k] += alpha * aValues[ka];
                                    break;
                                }
                        }
                    }

                    for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    {
                        var col = B.columnIndices[kb];
                        if (workspace[col] != i + 1)
                        {
                            values[pos] = beta * bValues[kb];
                            workspace[col] = i + 1;
                            pos++;
                        }
                        else
                        {
                            // Find position of this column in the output segment
                            for (var k = rowPtr[i]; k < pos; k++)
                                if (colIdx[k] == col)
                                {
                                    values[k] += beta * bValues[kb];
                                    break;
                                }
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
        var marker = ArrayPool<int>.Shared.Rent(A.ncols);
        try
        {
            Array.Fill(marker, -1, 0, A.ncols);
            var totalNnz = 0;
            rowPtr[0] = 0;

            for (var i = 0; i < A.nrows; i++)
            {
                // Mark all columns present in row i of A
                for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    marker[A.columnIndices[ka]] = i;

                // Count columns that are also present in row i of B
                var nnzThisRow = 0;
                for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                    if (marker[B.columnIndices[kb]] == i)
                        nnzThisRow++;

                totalNnz += nnzThisRow;
                rowPtr[i + 1] = totalNnz;
            }

            return rowPtr;
        }
        finally
        {
            ArrayPool<int>.Shared.Return(marker);
        }
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
        var marker = ArrayPool<int>.Shared.Rent(A.ncols);
        try
        {
            Array.Fill(marker, -1, 0, A.ncols);
            var pos = 0;

            for (var i = 0; i < A.nrows; i++)
            {
                // Mark all columns present in row i of A
                for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    marker[A.columnIndices[ka]] = i;

                // Add columns that are also present in row i of B
                for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                {
                    var col = B.columnIndices[kb];
                    if (marker[col] == i)
                        colIdx[pos++] = col;
                }
            }

            return colIdx;
        }
        finally
        {
            ArrayPool<int>.Shared.Return(marker);
        }
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
        var markerIdx = ArrayPool<int>.Shared.Rent(A.ncols);
        try
        {
            Array.Fill(markerIdx, -1, 0, A.ncols);
            var pos = 0;

            for (var i = 0; i < A.nrows; i++)
            {
                // Store index into A's values for each column in row i
                for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    markerIdx[A.columnIndices[ka]] = ka;

                // Compute Hadamard product for intersecting entries
                for (var kb = B.rowPointers[i]; kb < B.rowPointers[i + 1]; kb++)
                {
                    var col = B.columnIndices[kb];
                    var ka = markerIdx[col];
                    if (ka != -1) // Column exists in A for this row
                        values[pos++] = aValues[ka] * bValues[kb];
                }

                // Clear markers for this row
                for (var ka = A.rowPointers[i]; ka < A.rowPointers[i + 1]; ka++)
                    markerIdx[A.columnIndices[ka]] = -1;
            }

            return values;
        }
        finally
        {
            ArrayPool<int>.Shared.Return(markerIdx);
        }
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
                kOrigRow++;
                var col = columnIndices[k];
                var posT = Interlocked.Increment(ref next[col]) - 1;
                colIdxT[posT] = i;
                if (valuesT != null && localValues != null) valuesT[posT] = localValues[k];
                positions[posT] = kOrigRow;
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
        if (rhs.Length != nrows * nrhs)
            throw new ArgumentException(
                $"RHS must have {nrows * nrhs} elements for {nrhs} right-hand sides (got {rhs.Length})",
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
        var colSet = hashSetPool.Get();
        try
        {
            foreach (var col in cols) colSet.Add(col);
            var colMap = new Dictionary<int, int>();
            for (var j = 0; j < cols.Length; j++) colMap[cols[j]] = j;

            var resultRows = new List<List<int>>(rows.Length);
            List<double>? resultVals;
            double[]? localValues;
            lock (syncLock)
            {
                localValues = values;
                resultVals = localValues != null ? new List<double>() : null;
            }

            foreach (var row in rows)
            {
                List<int> newRow = [];
                for (var k = rowPointers[row]; k < rowPointers[row + 1]; k++)
                {
                    var col = columnIndices[k];
                    if (colSet.Contains(col))
                    {
                        newRow.Add(colMap[col]);
                        if (resultVals != null && localValues != null) resultVals.Add(localValues[k]);
                    }
                }

                resultRows.Add(newRow);
            }

            var result = new CSR(resultRows, true);
            if (resultVals != null) result.Values = resultVals.ToArray();
            return result;
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
                if (Math.Abs(diag) < 1e-14)
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
                if (Math.Abs(diag) < 1e-14)
                    throw new InvalidOperationException(
                        $"Zero or near-zero diagonal at row {i}: {diag}");

                x[i] /= diag;
            }
        }

        return x;
    }

    #endregion

    #region Factory Methods

    public static CSR Identity(int n)
    {
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
            var totalDofs = nodeCount * dofsPerNode;
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

    public override int GetHashCode()
    {
        ThrowIfDisposed();
        return HashCode.Combine(nrows, ncols, NonZeroCount);
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

file class DiagonalCache
{
    public DiagonalCache(CSR matrix)
    {
        var matrixValues = matrix.GetValuesInternal();
        if (matrix.Rows != matrix.Columns)
            throw new InvalidOperationException("Matrix must be square for diagonal extraction");

        var n = matrix.Rows;
        InverseDiagonal = new double[n];
        var rowPtrs = matrix.RowPointersArray;
        var colIndices = matrix.ColumnIndicesArray;

        // FIXED: Use relative tolerance based on matrix scale
        // First pass: find the maximum absolute diagonal value
        var maxAbsDiag = 0.0;
        for (var i = 0; i < n; i++)
        for (var k = rowPtrs[i]; k < rowPtrs[i + 1]; k++)
            if (colIndices[k] == i)
            {
                maxAbsDiag = Math.Max(maxAbsDiag, Math.Abs(matrixValues[k]));
                break;
            }

        // Use relative tolerance for zero detection
        var zeroTol = Math.Max(maxAbsDiag * 1e-14, 1e-300);

        for (var i = 0; i < n; i++)
        {
            var found = false;
            for (var k = rowPtrs[i]; k < rowPtrs[i + 1]; k++)
                if (colIndices[k] == i)
                {
                    if (Math.Abs(matrixValues[k]) < zeroTol)
                        throw new InvalidOperationException(
                            $"Zero or near-zero diagonal element at position {i}: value = {matrixValues[k]}, " +
                            $"threshold = {zeroTol} (relative to max diagonal {maxAbsDiag})");
                    InverseDiagonal[i] = 1.0 / matrixValues[k];
                    found = true;
                    break;
                }

            if (!found)
                throw new InvalidOperationException($"Missing diagonal entry at row {i}");
        }
    }

    public double[] InverseDiagonal { get; }
}

public static class CSRIterativeSolvers
{
    // FIXED: Issue #18 - Relative tolerance for breakdown detection
    public static double[] DiagonallyPreconditionedBiCGSTAB(
        this CSR matrix,
        double[] b,
        double tolerance = 1e-10,
        int maxIterations = 1000,
        double[]? x0 = null)
    {
        if (matrix == null) throw new ArgumentNullException(nameof(matrix));
        matrix.ThrowIfDisposed();
        if (matrix.Rows != matrix.Columns)
            throw new InvalidOperationException("Matrix must be square for iterative solvers");

        var n = matrix.Rows;
        var diagInv = new DiagonalCache(matrix).InverseDiagonal;
        var x = new double[n];
        if (x0 != null) Array.Copy(x0, x, n);

        var r = new double[n];
        Array.Copy(b, r, n);
        if (x0 != null)
        {
            var Ax = matrix.Multiply(x);
            for (var i = 0; i < n; i++) r[i] -= Ax[i];
        }

        var r0_hat = (double[])r.Clone();
        var rho = 1.0;
        var alpha = 1.0;
        var omega = 1.0;
        var v = new double[n];
        var p = new double[n];
        var s = new double[n];
        var t = new double[n];
        var phat = new double[n];
        var shat = new double[n];

        // FIXED: Issue #18 - Compute relative breakdown tolerance based on RHS norm
        // This prevents false breakdowns for systems with large/small values
        var bNormSq = 0.0;
        for (var i = 0; i < n; i++) bNormSq += b[i] * b[i];
        var bNorm = Math.Sqrt(bNormSq);

        // Use relative tolerance scaled by squared norm (since we're checking dot products)
        // Minimum absolute tolerance to handle zero RHS case
        var breakdownTol = Math.Max(tolerance * tolerance * bNormSq, 1e-300);

        for (var iter = 0; iter < maxIterations; iter++)
        {
            var rho_prev = rho;
            rho = 0.0;
            for (var i = 0; i < n; i++) rho += r0_hat[i] * r[i];

            // FIXED: Issue #18 - Use relative tolerance
            if (Math.Abs(rho) < breakdownTol)
                throw new InvalidOperationException(
                    $"BiCGSTAB breakdown: rho = {rho} is below threshold {breakdownTol}");

            if (iter == 0)
            {
                Array.Copy(r, p, n);
            }
            else
            {
                var beta = rho / rho_prev * (alpha / omega);
                for (var i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }

            for (var i = 0; i < n; i++) phat[i] = diagInv[i] * p[i];
            Array.Copy(matrix.Multiply(phat), v, n);
            alpha = rho;
            var temp = 0.0;
            for (var i = 0; i < n; i++) temp += r0_hat[i] * v[i];

            // FIXED: Issue #18 - Use relative tolerance
            if (Math.Abs(temp) < breakdownTol)
                throw new InvalidOperationException(
                    $"BiCGSTAB breakdown: (r0_hat, v) = {temp} is below threshold {breakdownTol}");
            alpha /= temp;

            for (var i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
            var sNorm = 0.0;
            for (var i = 0; i < n; i++) sNorm += s[i] * s[i];
            if (Math.Sqrt(sNorm) < tolerance)
            {
                for (var i = 0; i < n; i++) x[i] += alpha * phat[i];
                return x;
            }

            for (var i = 0; i < n; i++) shat[i] = diagInv[i] * s[i];
            Array.Copy(matrix.Multiply(shat), t, n);
            var ts = 0.0;
            var tt = 0.0;
            for (var i = 0; i < n; i++)
            {
                ts += t[i] * s[i];
                tt += t[i] * t[i];
            }

            // FIXED: Issue #18 - Use relative tolerance
            if (Math.Abs(tt) < breakdownTol)
                throw new InvalidOperationException(
                    $"BiCGSTAB breakdown: (t, t) = {tt} is below threshold {breakdownTol}");
            omega = ts / tt;

            for (var i = 0; i < n; i++)
            {
                x[i] += alpha * phat[i] + omega * shat[i];
                r[i] = s[i] - omega * t[i];
            }

            var rNorm = 0.0;
            for (var i = 0; i < n; i++) rNorm += r[i] * r[i];
            if (Math.Sqrt(rNorm) < tolerance) return x;
        }

        throw new InvalidOperationException($"BiCGSTAB did not converge in {maxIterations} iterations");
    }

    /// <summary>
    ///     Non-throwing variant of DiagonallyPreconditionedBiCGSTAB.
    ///     Returns a SolverResult with convergence information instead of throwing on failure.
    /// </summary>
    public static SolverResult TrySolve(
        this CSR matrix,
        double[] b,
        double tolerance = 1e-10,
        int maxIterations = 1000,
        double[]? x0 = null)
    {
        if (matrix == null) throw new ArgumentNullException(nameof(matrix));
        matrix.ThrowIfDisposed();
        if (matrix.Rows != matrix.Columns)
            throw new InvalidOperationException("Matrix must be square for iterative solvers");

        var n = matrix.Rows;
        DiagonalCache diagCache;
        try
        {
            diagCache = new DiagonalCache(matrix);
        }
        catch (InvalidOperationException ex)
        {
            return new SolverResult(
                new double[n], 0, double.NaN, false,
                $"Preconditioner failed: {ex.Message}");
        }

        var diagInv = diagCache.InverseDiagonal;
        var x = new double[n];
        if (x0 != null) Array.Copy(x0, x, n);

        var r = new double[n];
        Array.Copy(b, r, n);
        if (x0 != null)
        {
            var Ax = matrix.Multiply(x);
            for (var i = 0; i < n; i++) r[i] -= Ax[i];
        }

        var r0_hat = (double[])r.Clone();
        var rho = 1.0;
        var alpha = 1.0;
        var omega = 1.0;
        var v = new double[n];
        var p = new double[n];
        var s = new double[n];
        var t = new double[n];
        var phat = new double[n];
        var shat = new double[n];

        var bNormSq = 0.0;
        for (var i = 0; i < n; i++) bNormSq += b[i] * b[i];
        var bNorm = Math.Sqrt(bNormSq);
        var breakdownTol = Math.Max(tolerance * tolerance * bNormSq, 1e-300);

        for (var iter = 0; iter < maxIterations; iter++)
        {
            var rho_prev = rho;
            rho = 0.0;
            for (var i = 0; i < n; i++) rho += r0_hat[i] * r[i];

            if (Math.Abs(rho) < breakdownTol)
                return new SolverResult(x, iter + 1, ComputeResidualNorm(r, n), false,
                    $"BiCGSTAB breakdown: rho = {rho} is below threshold {breakdownTol}");

            if (iter == 0)
            {
                Array.Copy(r, p, n);
            }
            else
            {
                var beta = rho / rho_prev * (alpha / omega);
                for (var i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }

            for (var i = 0; i < n; i++) phat[i] = diagInv[i] * p[i];
            Array.Copy(matrix.Multiply(phat), v, n);
            alpha = rho;
            var temp = 0.0;
            for (var i = 0; i < n; i++) temp += r0_hat[i] * v[i];

            if (Math.Abs(temp) < breakdownTol)
                return new SolverResult(x, iter + 1, ComputeResidualNorm(r, n), false,
                    $"BiCGSTAB breakdown: (r0_hat, v) = {temp} is below threshold {breakdownTol}");
            alpha /= temp;

            for (var i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
            var sNorm = 0.0;
            for (var i = 0; i < n; i++) sNorm += s[i] * s[i];
            if (Math.Sqrt(sNorm) < tolerance)
            {
                for (var i = 0; i < n; i++) x[i] += alpha * phat[i];
                return new SolverResult(x, iter + 1, Math.Sqrt(sNorm), true);
            }

            for (var i = 0; i < n; i++) shat[i] = diagInv[i] * s[i];
            Array.Copy(matrix.Multiply(shat), t, n);
            var ts = 0.0;
            var tt = 0.0;
            for (var i = 0; i < n; i++)
            {
                ts += t[i] * s[i];
                tt += t[i] * t[i];
            }

            if (Math.Abs(tt) < breakdownTol)
                return new SolverResult(x, iter + 1, ComputeResidualNorm(r, n), false,
                    $"BiCGSTAB breakdown: (t, t) = {tt} is below threshold {breakdownTol}");
            omega = ts / tt;

            for (var i = 0; i < n; i++)
            {
                x[i] += alpha * phat[i] + omega * shat[i];
                r[i] = s[i] - omega * t[i];
            }

            var rNorm = ComputeResidualNorm(r, n);
            if (rNorm < tolerance)
                return new SolverResult(x, iter + 1, rNorm, true);
        }

        return new SolverResult(x, maxIterations, ComputeResidualNorm(r, n), false,
            $"BiCGSTAB did not converge in {maxIterations} iterations");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double ComputeResidualNorm(double[] r, int n)
    {
        var norm = 0.0;
        for (var i = 0; i < n; i++) norm += r[i] * r[i];
        return Math.Sqrt(norm);
    }

    public static double[] SmallestEigenvalues(this CSR matrix, int m, double tolerance = 1e-8, int maxIterations = 100)
    {
        throw new NotImplementedException("Full eigenvalue computation requires specialized libraries like ARPACK.");
    }
}

public record MatrixStatistics(
    int Rows,
    int Columns,
    int NonZeros,
    double Sparsity,
    int MinNnzPerRow,
    int MaxNnzPerRow,
    double AvgNnzPerRow);

public record SolverResult(
    double[] Solution,
    int Iterations,
    double ResidualNorm,
    bool Converged,
    string? Message = null);

public class SolverException : Exception
{
    public SolverException(string message) : base(message)
    {
    }

    public SolverException(string message, Exception inner) : base(message, inner)
    {
    }
}

/// <summary>
///     GPU backend using NVIDIA cuSPARSE - DYNAMIC RESOLUTION VERSION
///     Portability: Dynamically resolves CUDA 12/11 on Windows and Linux.
/// </summary>
internal sealed class CuSparseBackend : IDisposable
{
    // Constants
    private const int CUSPARSE_INDEX_32I = 1;
    private const int CUSPARSE_INDEX_BASE_ZERO = 0;
    private const int CUDA_R_64F = 6;
    private const int CUSPARSE_OPERATION_NON_TRANSPOSE = 0;
    private const int CUSPARSE_OPERATION_TRANSPOSE = 1;
    private const int CUSPARSE_SPMV_CSR_ALG2 = 3;
    private const int cudaMemcpyHostToDevice = 1;
    private const int cudaMemcpyDeviceToHost = 2;

    // --- 2. Function Pointers (Static Loading) ---
    private static IntPtr _libCuSparse;
    private static IntPtr _libCudaRt;

    private static CudaErrorDelegate? _cusparseCreate;
    private static CudaErrorSinglePtrDelegate? _cusparseDestroy;
    private static CreateCsrDelegate? _createCsr;
    private static CudaErrorSinglePtrDelegate? _destroySpMat;
    private static CreateDnVecDelegate? _createDnVec;
    private static CudaErrorSinglePtrDelegate? _destroyDnVec;
    private static SpMVBufferSizeDelegate? _spmvBufferSize;
    private static SpMVDelegate? _spmv;

    private static MallocDelegate? _cudaMalloc;
    private static FreeDelegate? _cudaFree;
    private static MemcpyDelegate? _cudaMemcpy;
    private static SynchronizeDelegate? _cudaDeviceSynchronize;

    // Complete native library configuration:
    // - Deep system scanning + pattern-based discovery
    // - Automatic latest version selection
    // - Package manager auto-install fallback (Linux/macOS)
    // - Works with ANY CUDA/MKL version forever
    private static PardisoSolver.INativeLibraryConfig _libraryConfig = new PardisoSolver.NativeLibraryConfig();
    private readonly object _disposeLock = new(); // FIXED: Thread-safe disposal

    private readonly int nrows, ncols, nnz;
    private ulong allocatedBufferSize;
    private IntPtr cusparseHandle = IntPtr.Zero;
    private IntPtr d_rowPtr = IntPtr.Zero, d_colInd = IntPtr.Zero, d_val = IntPtr.Zero;
    private IntPtr d_x = IntPtr.Zero, d_y = IntPtr.Zero, d_buffer = IntPtr.Zero;
    private volatile bool isDisposed;
    private bool isInitialized;
    private IntPtr spMatDescr = IntPtr.Zero;

    static CuSparseBackend()
    {
        LoadLibraries();
    }

    public CuSparseBackend(int rows, int cols, int nonZeros)
    {
        if (!IsAvailable)
            throw new PlatformNotSupportedException("CUDA libraries (cudart/cusparse) could not be loaded.");
        nrows = rows;
        ncols = cols;
        nnz = nonZeros;
    }

    public static bool IsAvailable { get; private set; }

    public void Dispose()
    {
        // FIXED: Thread-safe disposal to prevent double-free from finalizer race
        lock (_disposeLock)
        {
            if (isDisposed) return;
            isDisposed = true; // Set early to prevent re-entry
            isInitialized = false; // Prevent use of freed GPU resources
        }

        // Actual cleanup outside lock (CUDA calls may block)
        try
        {
            if (_cusparseDestroy != null && spMatDescr != IntPtr.Zero) _destroySpMat!(spMatDescr);
            if (_cudaFree != null)
            {
                if (d_buffer != IntPtr.Zero) _cudaFree(d_buffer);
                if (d_x != IntPtr.Zero) _cudaFree(d_x);
                if (d_y != IntPtr.Zero) _cudaFree(d_y);
                if (d_val != IntPtr.Zero) _cudaFree(d_val);
                if (d_colInd != IntPtr.Zero) _cudaFree(d_colInd);
                if (d_rowPtr != IntPtr.Zero) _cudaFree(d_rowPtr);
            }

            if (_cusparseDestroy != null && cusparseHandle != IntPtr.Zero) _cusparseDestroy(cusparseHandle);
        }
        catch (Exception ex)
        {
            Debug.WriteLine($"CUDA cleanup exception: {ex.Message}");
        }

        GC.SuppressFinalize(this);
    }

    /// <summary>
    ///     Sets the native library configuration.
    ///     Must be called before any CuSparseBackend instances are created.
    /// </summary>
    /// <param name="config">Configuration specifying library names and search paths</param>
    public static void SetLibraryConfig(PardisoSolver.INativeLibraryConfig config)
    {
        _libraryConfig = config ?? throw new ArgumentNullException(nameof(config));
    }

    private static void LoadLibraries()
    {
        try
        {
            // TURN-KEY: Use robust library loader with comprehensive auto-discovery
            var rtNames = _libraryConfig.GetCudaRuntimeLibraries();
            var spNames = _libraryConfig.GetCuSparseLibraries();
            var searchPaths = _libraryConfig.GetSearchPaths();

            if (rtNames == null || spNames == null)
            {
                Debug.WriteLine("[CuSparse] No CUDA libraries configured for this platform");
                IsAvailable = false;
                return;
            }

            Debug.WriteLine("[CuSparse] Searching for CUDA libraries...");

            // Load CUDA Runtime
            _libCudaRt = PardisoSolver.RobustNativeLibraryLoader.TryLoadLibrary(rtNames, searchPaths, "CUDA Runtime");

            if (_libCudaRt != IntPtr.Zero)
            {
                TryGetExport(_libCudaRt, "cudaMalloc", out _cudaMalloc);
                TryGetExport(_libCudaRt, "cudaFree", out _cudaFree);
                TryGetExport(_libCudaRt, "cudaMemcpy", out _cudaMemcpy);
                TryGetExport(_libCudaRt, "cudaDeviceSynchronize", out _cudaDeviceSynchronize);
            }

            // Load cuSPARSE
            _libCuSparse = PardisoSolver.RobustNativeLibraryLoader.TryLoadLibrary(spNames, searchPaths, "cuSPARSE");

            if (_libCuSparse != IntPtr.Zero)
            {
                TryGetExport(_libCuSparse, "cusparseCreate", out _cusparseCreate);
                TryGetExport(_libCuSparse, "cusparseDestroy", out _cusparseDestroy);
                TryGetExport(_libCuSparse, "cusparseCreateCsr", out _createCsr);
                TryGetExport(_libCuSparse, "cusparseDestroySpMat", out _destroySpMat);
                TryGetExport(_libCuSparse, "cusparseCreateDnVec", out _createDnVec);
                TryGetExport(_libCuSparse, "cusparseDestroyDnVec", out _destroyDnVec);
                TryGetExport(_libCuSparse, "cusparseSpMV_bufferSize", out _spmvBufferSize);
                TryGetExport(_libCuSparse, "cusparseSpMV", out _spmv);
            }

            IsAvailable = _cudaMalloc != null && _cusparseCreate != null;

            if (IsAvailable)
                Debug.WriteLine("[CuSparse] ✓ CUDA backend available");
            else
                Debug.WriteLine("[CuSparse] × CUDA backend not available");
        }
        catch (Exception ex)
        {
            Debug.WriteLine($"[CuSparse] Error: {ex.Message}");
            IsAvailable = false;
        }
    }

    private static bool TryGetExport<T>(IntPtr handle, string name, out T? del) where T : Delegate
    {
        if (NativeLibrary.TryGetExport(handle, name, out var addr))
        {
            del = Marshal.GetDelegateForFunctionPointer<T>(addr);
            return true;
        }

        del = null;
        return false;
    }

    public void Initialize(int[] rowPointers, int[] columnIndices, double[] values)
    {
        if (isInitialized) throw new InvalidOperationException("Already initialized");
        CheckCuda(_cusparseCreate!(ref cusparseHandle), "cusparseCreate");

        var rowPtrSize = (ulong)(nrows + 1) * sizeof(int);
        var colIndSize = (ulong)nnz * sizeof(int);
        var valSize = (ulong)nnz * sizeof(double);
        var vecSize = (ulong)Math.Max(nrows, ncols) * sizeof(double);

        CheckCuda(_cudaMalloc!(ref d_rowPtr, rowPtrSize), "malloc rowPtr");
        CheckCuda(_cudaMalloc!(ref d_colInd, colIndSize), "malloc colInd");
        CheckCuda(_cudaMalloc!(ref d_val, valSize), "malloc val");
        CheckCuda(_cudaMalloc!(ref d_x, vecSize), "malloc x");
        CheckCuda(_cudaMalloc!(ref d_y, vecSize), "malloc y");

        unsafe
        {
            fixed (int* pR = rowPointers, pC = columnIndices)
            fixed (double* pV = values)
            {
                CheckCuda(_cudaMemcpy!(d_rowPtr, (IntPtr)pR, rowPtrSize, cudaMemcpyHostToDevice), "memcpy rows");
                CheckCuda(_cudaMemcpy!(d_colInd, (IntPtr)pC, colIndSize, cudaMemcpyHostToDevice), "memcpy cols");
                CheckCuda(_cudaMemcpy!(d_val, (IntPtr)pV, valSize, cudaMemcpyHostToDevice), "memcpy vals");
            }
        }

        CheckCuda(
            _createCsr!(ref spMatDescr, nrows, ncols, nnz, d_rowPtr, d_colInd, d_val, CUSPARSE_INDEX_32I,
                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F), "createCsr");
        isInitialized = true;
    }

    public double[] Multiply(double[] x)
    {
        if (!isInitialized) throw new InvalidOperationException("Not initialized");
        unsafe
        {
            fixed (double* px = x)
            {
                CheckCuda(_cudaMemcpy!(d_x, (IntPtr)px, (ulong)ncols * sizeof(double), cudaMemcpyHostToDevice),
                    "memcpy x");
            }
        }

        IntPtr vecX = IntPtr.Zero, vecY = IntPtr.Zero;
        try
        {
            CheckCuda(_createDnVec!(ref vecX, ncols, d_x, CUDA_R_64F), "createDnVec X");
            CheckCuda(_createDnVec!(ref vecY, nrows, d_y, CUDA_R_64F), "createDnVec Y");
            double alpha = 1.0, beta = 0.0;
            ulong bufferSize = 0;
            unsafe
            {
                CheckCuda(
                    _spmvBufferSize!(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, (IntPtr)(&alpha), spMatDescr,
                        vecX, (IntPtr)(&beta), vecY, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2, ref bufferSize), "bufferSize");
                if (bufferSize > allocatedBufferSize)
                {
                    if (d_buffer != IntPtr.Zero) _cudaFree!(d_buffer);
                    d_buffer = IntPtr.Zero;
                    allocatedBufferSize = 0;
                    CheckCuda(_cudaMalloc!(ref d_buffer, bufferSize), "malloc buffer");
                    allocatedBufferSize = bufferSize;
                }

                CheckCuda(
                    _spmv!(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, (IntPtr)(&alpha), spMatDescr, vecX,
                        (IntPtr)(&beta), vecY, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2, d_buffer), "spmv");
            }

            CheckCuda(_cudaDeviceSynchronize!(), "sync");
        }
        finally
        {
            if (vecX != IntPtr.Zero) _destroyDnVec!(vecX);
            if (vecY != IntPtr.Zero) _destroyDnVec!(vecY);
        }

        var result = new double[nrows];
        unsafe
        {
            fixed (double* pr = result)
            {
                CheckCuda(_cudaMemcpy((IntPtr)pr, d_y, (ulong)nrows * sizeof(double), cudaMemcpyDeviceToHost),
                    "memcpy result");
            }
        }

        return result;
    }

    public double[] MultiplyTransposed(double[] x)
    {
        if (!isInitialized) throw new InvalidOperationException("Not initialized");
        unsafe
        {
            fixed (double* px = x)
            {
                CheckCuda(_cudaMemcpy!(d_x, (IntPtr)px, (ulong)nrows * sizeof(double), cudaMemcpyHostToDevice),
                    "memcpy x");
            }
        }

        IntPtr vecX = IntPtr.Zero, vecY = IntPtr.Zero;
        try
        {
            CheckCuda(_createDnVec!(ref vecX, nrows, d_x, CUDA_R_64F), "createDnVec X");
            CheckCuda(_createDnVec!(ref vecY, ncols, d_y, CUDA_R_64F), "createDnVec Y");
            double alpha = 1.0, beta = 0.0;
            ulong bufferSize = 0;
            unsafe
            {
                CheckCuda(
                    _spmvBufferSize!(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, (IntPtr)(&alpha), spMatDescr, vecX,
                        (IntPtr)(&beta), vecY, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2, ref bufferSize), "bufferSize");
                if (bufferSize > allocatedBufferSize)
                {
                    if (d_buffer != IntPtr.Zero) _cudaFree!(d_buffer);
                    d_buffer = IntPtr.Zero;
                    allocatedBufferSize = 0;
                    CheckCuda(_cudaMalloc!(ref d_buffer, bufferSize), "malloc buffer");
                    allocatedBufferSize = bufferSize;
                }

                CheckCuda(
                    _spmv!(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, (IntPtr)(&alpha), spMatDescr, vecX,
                        (IntPtr)(&beta), vecY, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2, d_buffer), "spmv");
            }

            CheckCuda(_cudaDeviceSynchronize!(), "sync");
        }
        finally
        {
            if (vecX != IntPtr.Zero) _destroyDnVec!(vecX);
            if (vecY != IntPtr.Zero) _destroyDnVec!(vecY);
        }

        var result = new double[ncols];
        unsafe
        {
            fixed (double* pr = result)
            {
                CheckCuda(_cudaMemcpy((IntPtr)pr, d_y, (ulong)ncols * sizeof(double), cudaMemcpyDeviceToHost),
                    "memcpy result");
            }
        }

        return result;
    }

    public void UpdateValues(double[] newValues)
    {
        if (!isInitialized) throw new InvalidOperationException("Not initialized");
        unsafe
        {
            fixed (double* pV = newValues)
            {
                CheckCuda(_cudaMemcpy!(d_val, (IntPtr)pV, (ulong)nnz * sizeof(double), cudaMemcpyHostToDevice),
                    "update values");
            }
        }
    }

    private void CheckCuda(int error, string msg)
    {
        if (error != 0) throw new Exception($"CUDA Error ({msg}): {error}");
    }

    ~CuSparseBackend()
    {
        Dispose();
    }

    // --- 1. Delegate Definitions ---
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int CudaErrorDelegate(ref IntPtr handle);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int CudaErrorSinglePtrDelegate(IntPtr handle);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int CreateCsrDelegate(ref IntPtr spMatDescr, long rows, long cols, long nnz, IntPtr csrRowOffsets,
        IntPtr csrColInd, IntPtr csrValues, int csrRowOffsetsType, int csrColIndType, int idxBase, int valueType);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int CreateDnVecDelegate(ref IntPtr dnVecDescr, long size, IntPtr values, int valueType);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int SpMVBufferSizeDelegate(IntPtr handle, int opA, IntPtr alpha, IntPtr matA, IntPtr vecX,
        IntPtr beta, IntPtr vecY, int computeType, int alg, ref ulong bufferSize);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int SpMVDelegate(IntPtr handle, int opA, IntPtr alpha, IntPtr matA, IntPtr vecX, IntPtr beta,
        IntPtr vecY, int computeType, int alg, IntPtr externalBuffer);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int MallocDelegate(ref IntPtr devPtr, ulong size);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int FreeDelegate(IntPtr devPtr);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int MemcpyDelegate(IntPtr dst, IntPtr src, ulong count, int kind);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int SynchronizeDelegate();
}

public static class SparseBackendFactory
{
    private const int MIN_ROWS_FOR_GPU = 50_000;
    private const int MIN_NNZ_FOR_GPU = 1_000_000;

    public static bool ShouldUseGPU(int nrows, int ncols, int nnz)
    {
        return ParallelConfig.EnableGPU && CuSparseBackend.IsAvailable &&
               (nrows >= MIN_ROWS_FOR_GPU || nnz >= MIN_NNZ_FOR_GPU);
    }
}

public sealed class HybridScheduler : IDisposable
{
    public enum BackendType
    {
        Auto,
        CPU_Sequential,
        CPU_Parallel,
        CPU_SIMD,
        GPU_CUDA,
        Hybrid_CPUandGPU
    }

    public enum OperationType
    {
        SpMV,
        SpMV_Transpose,
        SpMM,
        SpAdd,
        Solve
    }

    private readonly bool gpuAvailable;
    private readonly object gpuLock = new();
    private readonly CSR matrix;
    private readonly PerformanceMonitor perfMonitor;
    private CuSparseBackend? gpuBackend;
    private volatile bool gpuInitialized;
    private BackendType preferredBackend;

    public HybridScheduler(CSR matrix, BackendType preferredBackend = BackendType.Auto, bool eagerGpuInit = true)
    {
        this.matrix = matrix;
        this.preferredBackend = preferredBackend;
        perfMonitor = new PerformanceMonitor();
        gpuAvailable = SparseBackendFactory.ShouldUseGPU(matrix.Rows, matrix.Columns, matrix.NonZeroCount);

        if (preferredBackend == BackendType.Auto) AutoSelectBackend();
        if (eagerGpuInit && gpuAvailable)
            if (this.preferredBackend == BackendType.GPU_CUDA || this.preferredBackend == BackendType.Hybrid_CPUandGPU)
                try
                {
                    InitializeGPU();
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"GPU Eager init failed: {ex.Message}");
                }
    }

    public void Dispose()
    {
        lock (gpuLock)
        {
            gpuBackend?.Dispose();
            gpuBackend = null;
            gpuInitialized = false;
        }
    }

    public void InitializeGPU()
    {
        if (!gpuAvailable) throw new InvalidOperationException("GPU backend not available.");
        if (!gpuInitialized)
            lock (gpuLock)
            {
                if (!gpuInitialized)
                {
                    gpuBackend = new CuSparseBackend(matrix.Rows, matrix.Columns, matrix.NonZeroCount);
                    if (matrix.Values != null)
                        gpuBackend.Initialize(matrix.RowPointersArray, matrix.ColumnIndicesArray, matrix.Values);
                    gpuInitialized = true;
                }
            }
    }

    private void AutoSelectBackend()
    {
        var nrows = matrix.Rows;
        if (gpuAvailable && nrows >= 100_000) preferredBackend = BackendType.GPU_CUDA;
        else if (nrows >= 5_000) preferredBackend = BackendType.CPU_SIMD;
        else if (nrows >= 1_000) preferredBackend = BackendType.CPU_Parallel;
        else preferredBackend = BackendType.CPU_Sequential;
    }

    public double[] Execute(OperationType opType, double[] vector)
    {
        var sw = Stopwatch.StartNew();
        double[] result;
        try
        {
            result = opType switch
            {
                OperationType.SpMV => ExecuteSpMV(vector),
                OperationType.SpMV_Transpose => ExecuteSpMV_Transpose(vector),
                _ => throw new NotImplementedException()
            };
        }
        finally
        {
            sw.Stop();
            perfMonitor.Record(preferredBackend, opType, sw.Elapsed.TotalMilliseconds);
        }

        return result;
    }

    private double[] ExecuteSpMV(double[] vector)
    {
        return preferredBackend switch
        {
            BackendType.CPU_Sequential => matrix.Multiply(vector),
            BackendType.CPU_Parallel => matrix.MultiplyParallel(vector),
            BackendType.CPU_SIMD => matrix.MultiplySIMD(vector),
            BackendType.GPU_CUDA => ExecuteOnGPU(vector),
            _ => matrix.MultiplyParallel(vector)
        };
    }

    private double[] ExecuteSpMV_Transpose(double[] vector)
    {
        return preferredBackend switch
        {
            BackendType.CPU_Parallel => matrix.MultiplyTransposedParallel(vector),
            BackendType.CPU_SIMD => matrix.MultiplyTransposedParallel(vector),
            BackendType.GPU_CUDA => ExecuteTransposeOnGPU(vector),
            _ => matrix.MultiplyTransposedParallel(vector)
        };
    }

    private double[] ExecuteOnGPU(double[] vector)
    {
        if (!gpuInitialized) InitializeGPU();
        return gpuBackend!.Multiply(vector);
    }

    private double[] ExecuteTransposeOnGPU(double[] vector)
    {
        if (!gpuInitialized) InitializeGPU();
        return gpuBackend!.MultiplyTransposed(vector);
    }

    internal class PerformanceMonitor
    {
        private readonly ConcurrentDictionary<(BackendType, OperationType), List<double>> measurements = new();

        public void Record(BackendType backend, OperationType opType, double milliseconds)
        {
            var key = (backend, opType);
            var list = measurements.GetOrAdd(key, _ => new List<double>());
            lock (list)
            {
                list.Add(milliseconds);
            }
        }
    }
}

/// <summary>
///     Intel MKL PARDISO direct solver for sparse linear systems.
///     High-performance multithreaded solver for Ax = b.
///     This is a standalone class for advanced users who need direct control
///     over PARDISO. Basic users should use FEMSystem.Solve() instead.
///     Usage:
///     <code>
///     CSR matrix = GetSparseMatrix();
///     double[] rhs = GetRHS();
/// 
///     using var solver = new PardisoSolver();
///     double[] solution = solver.Solve(matrix, rhs);
///     </code>
/// </summary>
internal class PardisoSolver : IDisposable
{
    private const int PARDISO_PT_SIZE = 64;
    private const int PARDISO_IPARM_SIZE = 64;

    private static PardisoInitDelegate? _pardisoinit;
    private static PardisoDelegate? _pardiso;
    private static bool _isAvailable;
    private readonly int[] iparm = new int[64];
    private readonly IntPtr[] pt = new IntPtr[64];
    private volatile bool isDisposed;
    private bool isInitialized;
    private int matrixSize;
    private int matrixType = 11;

    static PardisoSolver()
    {
        LoadPardisoFunctions();
    }

    // --- Constructor and Methods ---

    public PardisoSolver(int mtype = 11)
    {
        if (!_isAvailable)
            throw new SolverException("Intel MKL runtime for PARDISO not found or failed initialization.");
        matrixType = mtype;

        // FIXED: Issue #3 - Proper initialization
        try
        {
            var mt = matrixType;
            _pardisoinit!(pt, ref mt, iparm);
        }
        catch (Exception ex)
        {
            throw new SolverException($"PARDISO initialization failed during instance creation: {ex.Message}");
        }
    }

    public void Dispose()
    {
        Dispose(true);
        GC.SuppressFinalize(this);
    }

    private static void LoadPardisoFunctions()
    {
        try
        {
            // ULTIMATE+ FUTURE-PROOF: Deep system scanning + latest version + auto-install fallback
            var config = new NativeLibraryConfig();
            var mklNames = config.GetMKLLibraries();
            var searchPaths = config.GetSearchPaths();

            if (mklNames == null)
            {
                Debug.WriteLine("[PARDISO] No MKL libraries configured for this platform");
                _isAvailable = false;
                return;
            }

            Debug.WriteLine("[PARDISO] Searching for PARDISO solver...");

            // Load MKL (PARDISO is part of MKL)
            var mklHandle = RobustNativeLibraryLoader.TryLoadLibrary(mklNames, searchPaths, "Intel MKL/PARDISO");

            if (mklHandle != IntPtr.Zero)
            {
                if (TryGetExport(mklHandle, "pardisoinit", out _pardisoinit) &&
                    TryGetExport(mklHandle, "pardiso", out _pardiso))
                {
                    _isAvailable = true;
                    Debug.WriteLine("[PARDISO] ✓ PARDISO solver available");
                }
                else
                {
                    Debug.WriteLine("[PARDISO] × PARDISO functions not found");
                    _isAvailable = false;
                }
            }
            else
            {
                Debug.WriteLine("[PARDISO] × MKL library not found");
                _isAvailable = false;
            }
        }
        catch (Exception ex)
        {
            Debug.WriteLine($"[PARDISO] Error: {ex.Message}");
            _isAvailable = false;
        }
    }

    private static bool TryGetExport<T>(IntPtr handle, string name, out T? del) where T : Delegate
    {
        if (NativeLibrary.TryGetExport(handle, name, out var address))
        {
            del = Marshal.GetDelegateForFunctionPointer<T>(address);
            return true;
        }

        del = null;
        return false;
    }

    public double[] Solve(CSR matrix, double[] rhs, bool refresh = false)
    {
        var n = matrix.Rows;
        var solution = new double[n];
        var ia = matrix.RowPointersArray;
        var ja = matrix.ColumnIndicesArray;
        var vals = matrix.GetValuesInternal();

        int maxfct = 1, mnum = 1, nrhs = 1, msglvl = 0, error = 0;
        var idum = new int[1];
        int phase;

        if (refresh || !isInitialized || matrixSize != n)
        {
            // FIXED: Issue #3 - Proper cleanup before reinitializing
            if (isInitialized)
            {
                // Release old factorization
                phase = -1;
                CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs,
                    iparm, ref msglvl, rhs, solution, ref error);

                // Clear pt array properly
                for (var i = 0; i < 64; i++) pt[i] = IntPtr.Zero;
                isInitialized = false;
            }

            // Reinitialize
            Array.Clear(iparm, 0, 64);
            CallPardisoInit(pt, ref matrixType, iparm);
            iparm[34] = 1; // 0-based indexing

            // Symbolic factorization
            phase = 11;
            CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs, iparm,
                ref msglvl, rhs, solution, ref error);
            if (error != 0)
                throw new SolverException($"PARDISO symbolic factorization failed with error code: {error}");

            isInitialized = true;
            matrixSize = n;
        }

        // Numerical factorization
        phase = 22;
        CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs, iparm,
            ref msglvl, rhs, solution, ref error);
        if (error != 0) throw new SolverException($"PARDISO numerical factorization failed with error code: {error}");

        // Solve
        phase = 33;
        CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs, iparm,
            ref msglvl, rhs, solution, ref error);
        if (error != 0) throw new SolverException($"PARDISO solve failed with error code: {error}");

        return solution;
    }

    public double[] SolveMultiple(CSR matrix, double[] rhs, int nrhs, bool refresh = false)
    {
        var n = matrix.Rows;
        var solution = new double[n * nrhs];
        var ia = matrix.RowPointersArray;
        var ja = matrix.ColumnIndicesArray;
        var vals = matrix.GetValuesInternal();

        int maxfct = 1, mnum = 1, msglvl = 0, error = 0;
        var idum = new int[1];
        int phase;

        if (refresh || !isInitialized)
        {
            // FIXED: Issue #3 - Proper cleanup
            if (isInitialized)
            {
                phase = -1;
                CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs,
                    iparm, ref msglvl, rhs, solution, ref error);
                for (var i = 0; i < 64; i++) pt[i] = IntPtr.Zero;
                isInitialized = false;
            }

            Array.Clear(iparm, 0, 64);
            CallPardisoInit(pt, ref matrixType, iparm);
            iparm[34] = 1;

            phase = 11;
            CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs, iparm,
                ref msglvl, rhs, solution, ref error);
            if (error != 0)
                throw new SolverException($"PARDISO symbolic factorization failed with error code: {error}");
            isInitialized = true;
        }

        phase = 22;
        CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs, iparm,
            ref msglvl, rhs, solution, ref error);
        if (error != 0) throw new SolverException($"PARDISO numerical factorization failed with error code: {error}");

        phase = 33;
        CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, vals, ia, ja, idum, ref nrhs, iparm,
            ref msglvl, rhs, solution, ref error);
        if (error != 0) throw new SolverException($"PARDISO solve failed with error code: {error}");

        return solution;
    }

    protected virtual void Dispose(bool disposing)
    {
        if (isDisposed) return;

        // MEDIUM FIX M7: Set disposed flag FIRST to ensure idempotency
        // If cleanup throws in finalizer, we won't retry endlessly
        isDisposed = true;

        if (isInitialized && disposing)
            try
            {
                int phase = -1, maxfct = 1, mnum = 1, nrhs = 1, msglvl = 0, error = 0;
                var dummy = new double[1];
                var idum = new int[1];
                var n = matrixSize;
                CallPardiso(pt, ref maxfct, ref mnum, ref matrixType, ref phase, ref n, dummy, idum, idum, idum,
                    ref nrhs, iparm, ref msglvl, dummy, dummy, ref error);

                // Only clear if cleanup succeeded
                isInitialized = false;
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"PARDISO cleanup exception: {ex.Message}");
                // Don't rethrow from Dispose - swallow exceptions per .NET guidelines
            }
    }

    ~PardisoSolver()
    {
        Dispose(false);
    }

    // DYNAMIC WRAPPERS - Use the static delegates resolved by NativeLibrary
    private void CallPardisoInit(IntPtr[] pt, ref int mtype, int[] iparm)
    {
        _pardisoinit!(pt, ref mtype, iparm);
    }

    private void CallPardiso(IntPtr[] pt, ref int maxfct, ref int mnum, ref int mtype, ref int phase, ref int n,
        double[] a, int[] ia, int[] ja, int[] perm, ref int nrhs, int[] iparm, ref int msglvl, double[] b, double[] x,
        ref int error)
    {
        _pardiso!(pt, ref maxfct, ref mnum, ref mtype, ref phase, ref n, a, ia, ja, perm, ref nrhs, iparm, ref msglvl,
            b, x, ref error);
    }

    // --- Dynamic Loading Delegates and State ---

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate void PardisoInitDelegate(IntPtr[] pt, ref int mtype, int[] iparm);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate void PardisoDelegate(IntPtr[] pt, ref int maxfct, ref int mnum, ref int mtype, ref int phase,
        ref int n, double[] a, int[] ia, int[] ja, int[] perm, ref int nrhs, int[] iparm, ref int msglvl, double[] b,
        double[] x, ref int error);

    #region Optimized Sparse Matrix Kernels

    /// <summary>
    ///     Ultra-optimized sparse matrix-vector multiplication kernels
    /// </summary>
    internal static class OptimizedSpMV
    {
        [MethodImpl(AggressiveOptimization)]
        public static unsafe void MultiplyVector(
            int* rowPtr, int* colInd, double* values,
            double* x, double* y,
            int numRows)
        {
            // Optimized SpMV: y = A * x
            // where A is in CSR format

            if (Avx512F.IsSupported)
                MultiplyVector_AVX512(rowPtr, colInd, values, x, y, numRows);
            else if (Avx2.IsSupported)
                MultiplyVector_AVX2(rowPtr, colInd, values, x, y, numRows);
            else
                MultiplyVector_Scalar(rowPtr, colInd, values, x, y, numRows);
        }

        [MethodImpl(AggressiveOptimization)]
        private static unsafe void MultiplyVector_AVX512(
            int* rowPtr, int* colInd, double* values,
            double* x, double* y,
            int numRows)
        {
            // Process multiple rows at once for better instruction-level parallelism
            const int UNROLL = 4;
            var row = 0;

            for (; row + UNROLL - 1 < numRows; row += UNROLL)
            {
                // Process 4 rows simultaneously
                ProcessRow_AVX512(rowPtr, colInd, values, x, y, row + 0);
                ProcessRow_AVX512(rowPtr, colInd, values, x, y, row + 1);
                ProcessRow_AVX512(rowPtr, colInd, values, x, y, row + 2);
                ProcessRow_AVX512(rowPtr, colInd, values, x, y, row + 3);
            }

            // Handle remaining rows
            for (; row < numRows; row++) ProcessRow_AVX512(rowPtr, colInd, values, x, y, row);
        }

        [MethodImpl(AggressiveInlining | AggressiveOptimization)]
        private static unsafe void ProcessRow_AVX512(
            int* rowPtr, int* colInd, double* values,
            double* x, double* y, int row)
        {
            var rowStart = rowPtr[row];
            var rowEnd = rowPtr[row + 1];
            var nnz = rowEnd - rowStart;

            if (nnz == 0)
            {
                y[row] = 0.0;
                return;
            }

            // Use AVX-512 for reduction (8 doubles at a time)
            var sum = Vector512<double>.Zero;
            var j = rowStart;

            // Allocate gather buffer once outside loop to avoid stack overflow
            Span<double> xGathered = stackalloc double[8];

            // Main vectorized loop: process 8 non-zeros at a time
            for (; j + 7 < rowEnd; j += 8)
            {
                // Prefetch next iteration (heuristic: 16 elements ahead)
                if (j + 16 < rowEnd)
                {
                    Sse.Prefetch0(values + j + 16);
                    Sse.Prefetch0(colInd + j + 16);
                }

                // Load 8 values
                var vals = Avx512F.LoadVector512(values + j);

                // Manual gather of x-values (GatherVector512 not yet available in .NET)
                // This is the bottleneck in SpMV - irregular memory access
                for (var g = 0; g < 8; g++) xGathered[g] = x[colInd[j + g]];

                fixed (double* pGathered = xGathered)
                {
                    var xVals = Avx512F.LoadVector512(pGathered);
                    // Multiply and accumulate
                    sum = Avx512F.FusedMultiplyAdd(vals, xVals, sum);
                }
            }

            // Horizontal reduction of sum vector
            var total = HorizontalSum_AVX512(sum);

            // Process remaining elements (< 8)
            for (; j < rowEnd; j++) total += values[j] * x[colInd[j]];

            y[row] = total;
        }

        [MethodImpl(AggressiveInlining)]
        private static double HorizontalSum_AVX512(Vector512<double> v)
        {
            // Efficient horizontal sum of 8 doubles
            var low = v.GetLower(); // Lower 4 doubles
            var high = v.GetUpper(); // Upper 4 doubles

            var sum4 = Avx.Add(low, high); // Now 4 doubles

            // Continue reduction to scalar
            var low2 = Avx.ExtractVector128(sum4, 0);
            var high2 = Avx.ExtractVector128(sum4, 1);
            var sum2 = Sse2.Add(low2, high2);

            var high1 = Sse2.UnpackHigh(sum2, sum2);
            var sum1 = Sse2.Add(sum2, high1);

            return sum1.ToScalar();
        }

        [MethodImpl(AggressiveOptimization)]
        private static unsafe void MultiplyVector_AVX2(
            int* rowPtr, int* colInd, double* values,
            double* x, double* y,
            int numRows)
        {
            const int UNROLL = 4;
            var row = 0;

            for (; row + UNROLL - 1 < numRows; row += UNROLL)
            {
                ProcessRow_AVX2(rowPtr, colInd, values, x, y, row + 0);
                ProcessRow_AVX2(rowPtr, colInd, values, x, y, row + 1);
                ProcessRow_AVX2(rowPtr, colInd, values, x, y, row + 2);
                ProcessRow_AVX2(rowPtr, colInd, values, x, y, row + 3);
            }

            for (; row < numRows; row++) ProcessRow_AVX2(rowPtr, colInd, values, x, y, row);
        }

        [MethodImpl(AggressiveInlining | AggressiveOptimization)]
        private static unsafe void ProcessRow_AVX2(
            int* rowPtr, int* colInd, double* values,
            double* x, double* y, int row)
        {
            var rowStart = rowPtr[row];
            var rowEnd = rowPtr[row + 1];

            if (rowStart == rowEnd)
            {
                y[row] = 0.0;
                return;
            }

            // AVX2: process 4 doubles at a time
            var sum = Vector256<double>.Zero;
            var j = rowStart;

            // Main loop: 4 at a time
            for (; j + 3 < rowEnd; j += 4)
            {
                // Load 4 values
                var vals = Avx.LoadVector256(values + j);

                // Manual gather (AVX2 doesn't have double gather)
                var xVals = Vector256.Create(
                    x[colInd[j]],
                    x[colInd[j + 1]],
                    x[colInd[j + 2]],
                    x[colInd[j + 3]]
                );

                // FMA
                sum = Fma.IsSupported
                    ? Fma.MultiplyAdd(vals, xVals, sum)
                    : Avx.Add(Avx.Multiply(vals, xVals), sum);
            }

            // Horizontal sum
            var total = HorizontalSum_AVX2(sum);

            // Remaining elements
            for (; j < rowEnd; j++) total += values[j] * x[colInd[j]];

            y[row] = total;
        }

        [MethodImpl(AggressiveInlining)]
        private static double HorizontalSum_AVX2(Vector256<double> v)
        {
            var low = Avx.ExtractVector128(v, 0);
            var high = Avx.ExtractVector128(v, 1);
            var sum = Sse2.Add(low, high);

            var high1 = Sse2.UnpackHigh(sum, sum);
            var final = Sse2.Add(sum, high1);

            return final.ToScalar();
        }

        [MethodImpl(AggressiveOptimization)]
        private static unsafe void MultiplyVector_Scalar(
            int* rowPtr, int* colInd, double* values,
            double* x, double* y,
            int numRows)
        {
            // Optimized scalar fallback - still better than naive
            for (var row = 0; row < numRows; row++)
            {
                var rowStart = rowPtr[row];
                var rowEnd = rowPtr[row + 1];

                var sum = 0.0;

                // Unroll by 4 for better ILP
                var j = rowStart;
                for (; j + 3 < rowEnd; j += 4)
                {
                    sum += values[j] * x[colInd[j]];
                    sum += values[j + 1] * x[colInd[j + 1]];
                    sum += values[j + 2] * x[colInd[j + 2]];
                    sum += values[j + 3] * x[colInd[j + 3]];
                }

                for (; j < rowEnd; j++) sum += values[j] * x[colInd[j]];

                y[row] = sum;
            }
        }

        /// <summary>
        ///     Optimized SpMV with scaling: y = alpha * A * x + beta * y
        /// </summary>
        [MethodImpl(AggressiveOptimization)]
        public static unsafe void MultiplyVectorScaled(
            int* rowPtr, int* colInd, double* values,
            double* x, double* y,
            int numRows,
            double alpha, double beta)
        {
            if (alpha == 0.0)
            {
                if (beta == 0.0)
                    // y = 0
                    for (var i = 0; i < numRows; i++)
                        y[i] = 0.0;
                else if (beta != 1.0)
                    // y = beta * y
                    ScaleVector(y, numRows, beta);
                // else: y = y (no-op)
                return;
            }

            if (beta == 0.0)
            {
                // y = alpha * A * x
                MultiplyVector(rowPtr, colInd, values, x, y, numRows);
                if (alpha != 1.0)
                    ScaleVector(y, numRows, alpha);
            }
            else
            {
                // General case: y = alpha * A * x + beta * y
                // Use temporary buffer to avoid read-modify-write conflicts
                Span<double> temp = stackalloc double[Math.Min(numRows, 1024)];

                for (var blockStart = 0; blockStart < numRows; blockStart += temp.Length)
                {
                    var blockSize = Math.Min(temp.Length, numRows - blockStart);

                    fixed (double* pTemp = temp)
                    {
                        // Compute temp = A[block] * x
                        MultiplyVector(rowPtr + blockStart, colInd, values,
                            x, pTemp, blockSize);

                        // y[block] = alpha * temp + beta * y[block]
                        for (var i = 0; i < blockSize; i++)
                            y[blockStart + i] = alpha * pTemp[i] + beta * y[blockStart + i];
                    }
                }
            }
        }

        [MethodImpl(AggressiveInlining | AggressiveOptimization)]
        private static unsafe void ScaleVector(double* v, int n, double alpha)
        {
            var i = 0;

            if (Avx512F.IsSupported)
            {
                var valpha = Vector512.Create(alpha);
                for (; i + 7 < n; i += 8)
                {
                    var vec = Avx512F.LoadVector512(v + i);
                    vec = Avx512F.Multiply(vec, valpha);
                    Avx512F.Store(v + i, vec);
                }
            }
            else if (Avx.IsSupported)
            {
                var valpha = Vector256.Create(alpha);
                for (; i + 3 < n; i += 4)
                {
                    var vec = Avx.LoadVector256(v + i);
                    vec = Avx.Multiply(vec, valpha);
                    Avx.Store(v + i, vec);
                }
            }

            for (; i < n; i++) v[i] *= alpha;
        }
    }

    /// <summary>
    ///     Optimized sparse matrix-matrix multiply and other CSR operations
    /// </summary>
    internal static class OptimizedSpMM
    {
        /// <summary>
        ///     Sparse matrix-matrix multiply: C = A * B where A is sparse CSR, B is dense
        ///     Highly optimized for this common case in FEA
        ///     FIXED: Now correctly handles column-major storage with arbitrary leading dimensions.
        ///     The algorithm processes one output element at a time using dot products, which
        ///     naturally handles strided column-major access patterns.
        ///     Storage layout: Column-major (Fortran/BLAS convention)
        ///     B[row, col] is at B[row + col * ldb]
        ///     C[row, col] is at C[row + col * ldc]
        /// </summary>
        [MethodImpl(AggressiveOptimization)]
        public static unsafe void MultiplyDense(
            int* rowPtr, int* colInd, double* aValues,
            double* B, int ldb,
            double* C, int ldc,
            int m, int n, int k)
        {
            // C(m x n) = A(m x k) * B(k x n)
            // A is sparse CSR, B and C are dense column-major

            // Initialize C to zero
            for (var j = 0; j < n; j++)
            for (var i = 0; i < m; i++)
                C[i + j * ldc] = 0.0;

            // For each row of A (row i of output C)
            for (var i = 0; i < m; i++)
            {
                var rowStart = rowPtr[i];
                var rowEnd = rowPtr[i + 1];
                var nnzRow = rowEnd - rowStart;

                if (nnzRow == 0) continue;

                // For each column j of output C
                // C[i,j] = sum over non-zeros in row i of A: A[i,k] * B[k,j]
                for (var j = 0; j < n; j++)
                {
                    var sum = 0.0;

                    // Accumulate dot product of sparse row i of A with column j of B
                    // Use SIMD for the accumulation when we have enough non-zeros
                    if (Avx2.IsSupported && nnzRow >= 4)
                    {
                        var sumVec = Vector256<double>.Zero;
                        var jj = rowStart;

                        // Process 4 non-zeros at a time
                        for (; jj + 3 < rowEnd; jj += 4)
                        {
                            // Load 4 values from A
                            var aVec = Avx.LoadVector256(aValues + jj);

                            // Gather 4 values from column j of B (strided access)
                            // B[colInd[jj], j], B[colInd[jj+1], j], etc.
                            var bVec = Vector256.Create(
                                B[colInd[jj] + j * ldb],
                                B[colInd[jj + 1] + j * ldb],
                                B[colInd[jj + 2] + j * ldb],
                                B[colInd[jj + 3] + j * ldb]
                            );

                            // FMA: sumVec += aVec * bVec
                            sumVec = Fma.IsSupported
                                ? Fma.MultiplyAdd(aVec, bVec, sumVec)
                                : Avx.Add(Avx.Multiply(aVec, bVec), sumVec);
                        }

                        // Horizontal sum of vector
                        var low = Avx.ExtractVector128(sumVec, 0);
                        var high = Avx.ExtractVector128(sumVec, 1);
                        var sum2 = Sse2.Add(low, high);
                        var high1 = Sse2.UnpackHigh(sum2, sum2);
                        var final = Sse2.Add(sum2, high1);
                        sum = final.ToScalar();

                        // Handle remaining non-zeros
                        for (; jj < rowEnd; jj++) sum += aValues[jj] * B[colInd[jj] + j * ldb];
                    }
                    else
                    {
                        // Scalar path for small nnz or no AVX2
                        for (var jj = rowStart; jj < rowEnd; jj++) sum += aValues[jj] * B[colInd[jj] + j * ldb];
                    }

                    C[i + j * ldc] = sum;
                }
            }
        }

        /// <summary>
        ///     Alternative SpMM using outer product formulation: C += A[:,k] * B[k,:]
        ///     Better cache behavior when B has many columns and rows of A are dense.
        ///     Correctly handles column-major storage.
        /// </summary>
        [MethodImpl(AggressiveOptimization)]
        public static unsafe void MultiplyDenseOuterProduct(
            int* rowPtr, int* colInd, double* aValues,
            double* B, int ldb,
            double* C, int ldc,
            int m, int n, int k)
        {
            // C(m x n) = A(m x k) * B(k x n)
            // Initialize C to zero
            for (var j = 0; j < n; j++)
            for (var i = 0; i < m; i++)
                C[i + j * ldc] = 0.0;

            // Outer product formulation: for each row i of A
            for (var i = 0; i < m; i++)
            {
                var rowStart = rowPtr[i];
                var rowEnd = rowPtr[i + 1];

                // For each non-zero A[i, k_idx]
                for (var jj = rowStart; jj < rowEnd; jj++)
                {
                    var k_idx = colInd[jj];
                    var aVal = aValues[jj];

                    // C[i, :] += aVal * B[k_idx, :]
                    // This is a SAXPY on row i of C using row k_idx of B
                    // Column-major: elements of row k_idx of B are at B[k_idx + j*ldb] for j=0..n-1
                    for (var j = 0; j < n; j++) C[i + j * ldc] += aVal * B[k_idx + j * ldb];
                }
            }
        }
    }

    #endregion

    #region Native Library Loading

// ============================================================================
// CUDA FUNCTION DELEGATES (for verification)
// ============================================================================

/// <summary>
///     Delegate for cudaGetDeviceCount function.
///     Returns cudaError_t (int): 0 = success, 35 = insufficient driver, 100 = no device
/// </summary>
[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    internal delegate int CudaGetDeviceCountDelegate(ref int count);

    /// <summary>
    ///     Delegate for cudaRuntimeGetVersion function.
    ///     Returns cudaError_t (int): 0 = success
    /// </summary>
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    internal delegate int CudaGetVersionDelegate(ref int version);

// ============================================================================
// NATIVE LIBRARY CONFIGURATION INTERFACE
// ============================================================================

/// Configuration interface for native library paths.
/// Allows applications to specify custom library locations or versions.
/// </summary>
public interface INativeLibraryConfig
    {
        /// <summary>
        ///     Gets the ordered list of CUDA Runtime library names to try loading.
        ///     Returns null if CUDA should not be used.
        /// </summary>
        string[]? GetCudaRuntimeLibraries();

        /// <summary>
        ///     Gets the ordered list of cuSPARSE library names to try loading.
        ///     Returns null if cuSPARSE should not be used.
        /// </summary>
        string[]? GetCuSparseLibraries();

        /// <summary>
        ///     Gets the ordered list of Intel MKL library names to try loading.
        ///     Returns null if MKL should not be used.
        /// </summary>
        string[]? GetMKLLibraries();

        /// <summary>
        ///     Gets the ordered list of directories to search for native libraries.
        ///     Libraries are searched in the order paths are returned.
        /// </summary>
        string[]? GetSearchPaths();
    }

// ============================================================================
// NATIVE LIBRARY CONFIGURATION - Discovery and Auto-Install
// ============================================================================

/// <summary>
///     COMPLETE native library configuration - everything in one place.
///     Features:
///     - Deep recursive system scanning (Windows/Linux/macOS)
///     - Pattern-based discovery (works with ANY CUDA/MKL version)
///     - Automatic latest version selection
///     - Directory structure resilient
///     - Automatic MKL installation via package manager (Linux/macOS)
///     - Zero maintenance forever
/// </summary>
public sealed class NativeLibraryConfig : INativeLibraryConfig
    {
        private static readonly object _cacheLock = new();
        private static string[]? _cachedCudaRt;
        private static string[]? _cachedCuSparse;
        private static string[]? _cachedMKL;
        private static string[]? _cachedSearchPaths;
        private static bool _mklInstallAttempted;

        /// <summary>
        ///     If true, automatically attempts MKL installation via package manager if not found.
        ///     Default: true (can be disabled for air-gapped systems or CI/CD)
        /// </summary>
        public bool EnableAutoInstall { get; set; } = true;

        /// <summary>
        ///     If true, shows interactive prompts before installing MKL.
        ///     Default: true (set false for silent/CI environments)
        /// </summary>
        public bool InteractiveInstall { get; set; } = true;

        /// <summary>
        ///     If true, performs deep recursive scanning of system directories to find libraries.
        ///     Default: false (PERFORMANCE FIX: deep scanning can cause massive IO spikes on startup)
        ///     When disabled, only checks:
        ///     - Known/standard installation paths (CUDA_PATH, MKLROOT, etc.)
        ///     - Environment variables (PATH, LD_LIBRARY_PATH)
        ///     - Standard system library directories
        ///     Enable this only if libraries are installed in non-standard locations.
        /// </summary>
        public bool EnableDeepScanning { get; set; } = false;

        public string[]? GetCudaRuntimeLibraries()
        {
            lock (_cacheLock)
            {
                if (_cachedCudaRt == null)
                {
                    // Check for environment variable override
                    // Example: CUDA_FORCE_VERSION=11.0 or CUDA_FORCE_PATH=/path/to/libcudart.so
                    var forceVersion = Environment.GetEnvironmentVariable("CUDA_FORCE_VERSION");
                    var forcePath = Environment.GetEnvironmentVariable("CUDA_FORCE_PATH");

                    if (!string.IsNullOrEmpty(forcePath) && File.Exists(forcePath))
                    {
                        Debug.WriteLine($"[NativeLibrary] Using CUDA_FORCE_PATH: {forcePath}");
                        _cachedCudaRt = new[] { forcePath };
                        return _cachedCudaRt;
                    }

                    if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                        _cachedCudaRt = DiscoverLibrariesInSystem(LibraryPatterns.CudaRuntimeWindows, "CUDA Runtime");
                    else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                        _cachedCudaRt = DiscoverLibrariesInSystem(LibraryPatterns.CudaRuntimeLinux, "CUDA Runtime");

                    // Apply version filter if specified
                    if (!string.IsNullOrEmpty(forceVersion) && _cachedCudaRt != null)
                    {
                        var filtered = _cachedCudaRt.Where(lib => lib.Contains(forceVersion)).ToArray();
                        if (filtered.Length > 0)
                        {
                            Debug.WriteLine($"[NativeLibrary] Filtering CUDA to version {forceVersion}");
                            _cachedCudaRt = filtered;
                        }
                        else
                        {
                            Debug.WriteLine($"[NativeLibrary] Warning: CUDA_FORCE_VERSION={forceVersion} not found");
                        }
                    }
                }

                return _cachedCudaRt;
            }
        }

        public string[]? GetCuSparseLibraries()
        {
            lock (_cacheLock)
            {
                if (_cachedCuSparse == null)
                {
                    if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                        _cachedCuSparse = DiscoverLibrariesInSystem(LibraryPatterns.CuSparseWindows, "cuSPARSE");
                    else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                        _cachedCuSparse = DiscoverLibrariesInSystem(LibraryPatterns.CuSparseLinux, "cuSPARSE");
                }

                return _cachedCuSparse;
            }
        }

        public string[]? GetMKLLibraries()
        {
            lock (_cacheLock)
            {
                // Check for environment variable override
                var forceVersion = Environment.GetEnvironmentVariable("MKL_FORCE_VERSION");
                var forcePath = Environment.GetEnvironmentVariable("MKL_FORCE_PATH");

                if (!string.IsNullOrEmpty(forcePath) && File.Exists(forcePath))
                {
                    Debug.WriteLine($"[NativeLibrary] Using MKL_FORCE_PATH: {forcePath}");
                    _cachedMKL = new[] { forcePath };
                    return _cachedMKL;
                }

                // Try discovery first
                if (_cachedMKL == null)
                {
                    if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                        _cachedMKL = DiscoverLibrariesInSystem(LibraryPatterns.MKLWindows, "Intel MKL");
                    else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                        _cachedMKL = DiscoverLibrariesInSystem(LibraryPatterns.MKLLinux, "Intel MKL");
                    else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                        _cachedMKL = DiscoverLibrariesInSystem(LibraryPatterns.MKLMacOS, "Intel MKL");

                    // Apply version filter if specified
                    if (!string.IsNullOrEmpty(forceVersion) && _cachedMKL != null)
                    {
                        var filtered = _cachedMKL.Where(lib => lib.Contains(forceVersion)).ToArray();
                        if (filtered.Length > 0)
                        {
                            Debug.WriteLine($"[NativeLibrary] Filtering MKL to version {forceVersion}");
                            _cachedMKL = filtered;
                        }
                        else
                        {
                            Debug.WriteLine($"[NativeLibrary] Warning: MKL_FORCE_VERSION={forceVersion} not found");
                        }
                    }
                }

                // ENHANCED: Check if discovered libraries actually work
                // If they all fail verification, try auto-install
                var shouldTryAutoInstall = false;

                if (_cachedMKL != null && _cachedMKL.Length > 0)
                {
                    // Quick verification: try to load first library to see if any work
                    var searchPaths = GetSearchPaths();
                    if (searchPaths != null)
                    {
                        var anyWorking = false;
                        foreach (var libName in _cachedMKL.Take(3)) // Test first 3 versions
                        {
                            foreach (var searchPath in searchPaths.Take(10)) // Test first 10 paths
                                try
                                {
                                    var fullPath = Path.Combine(searchPath, libName);
                                    if (File.Exists(fullPath))
                                        if (NativeLibrary.TryLoad(fullPath, out var testHandle))
                                        {
                                            // Check if it works
                                            var works =
                                                NativeLibrary.TryGetExport(testHandle, "mkl_get_version", out _) ||
                                                NativeLibrary.TryGetExport(testHandle, "MKL_Get_Max_Threads", out _);
                                            NativeLibrary.Free(testHandle);

                                            if (works)
                                            {
                                                anyWorking = true;
                                                break;
                                            }
                                        }
                                }
                                catch
                                {
                                }

                            if (anyWorking) break;
                        }

                        if (!anyWorking)
                        {
                            Debug.WriteLine("[NativeLibrary] MKL libraries found but all failed verification!");
                            shouldTryAutoInstall = true;
                        }
                    }
                }

                // If MKL not found OR all found libraries failed verification, try auto-install
                if ((_cachedMKL == null || _cachedMKL.Length == 0 || shouldTryAutoInstall) &&
                    EnableAutoInstall &&
                    !_mklInstallAttempted)
                {
                    _mklInstallAttempted = true;

                    if (shouldTryAutoInstall)
                        Debug.WriteLine(
                            "[NativeLibrary] All MKL libraries failed verification, attempting to install working version...");
                    else
                        Debug.WriteLine(
                            "[NativeLibrary] MKL not found via scanning, attempting automatic installation...");

                    if (TryAutoInstallMKL())
                    {
                        // Re-scan after installation
                        Debug.WriteLine("[NativeLibrary] Re-scanning for MKL after installation...");

                        // Clear cache to force re-scan
                        _cachedSearchPaths = null;

                        // Try discovery again
                        if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                            _cachedMKL = DiscoverLibrariesInSystem(LibraryPatterns.MKLWindows, "Intel MKL");
                        else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                            _cachedMKL = DiscoverLibrariesInSystem(LibraryPatterns.MKLLinux, "Intel MKL");
                        else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                            _cachedMKL = DiscoverLibrariesInSystem(LibraryPatterns.MKLMacOS, "Intel MKL");

                        if (_cachedMKL != null && _cachedMKL.Length > 0)
                            Debug.WriteLine("[NativeLibrary] ✓ MKL found after installation!");
                    }
                }

                return _cachedMKL;
            }
        }

        public string[]? GetSearchPaths()
        {
            lock (_cacheLock)
            {
                if (_cachedSearchPaths == null) _cachedSearchPaths = BuildComprehensiveSearchPaths();
                return _cachedSearchPaths;
            }
        }

        /// <summary>
        ///     Builds comprehensive search paths across all common locations.
        /// </summary>
        private string[] BuildComprehensiveSearchPaths()
        {
            var paths = new List<string>();
            var checkedPaths = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

            void AddPath(string path, string description)
            {
                if (!string.IsNullOrEmpty(path) && Directory.Exists(path))
                    try
                    {
                        var normalized = Path.GetFullPath(path);
                        if (checkedPaths.Add(normalized))
                        {
                            paths.Add(normalized);
                            Debug.WriteLine($"[NativeLibrary] Search path: {description} → {normalized}");
                        }
                    }
                    catch
                    {
                    }
            }

            // PRIORITY 1: Current working directory
            try
            {
                AddPath(Directory.GetCurrentDirectory(), "pwd");
            }
            catch
            {
            }

            // PRIORITY 2: Application base directory
            try
            {
                AddPath(AppDomain.CurrentDomain.BaseDirectory, "app base");
            }
            catch
            {
            }

            // PRIORITY 3: Assembly location
            try
            {
                var assemblyLocation = Assembly.GetExecutingAssembly().Location;
                if (!string.IsNullOrEmpty(assemblyLocation))
                {
                    var assemblyDir = Path.GetDirectoryName(assemblyLocation);
                    if (!string.IsNullOrEmpty(assemblyDir))
                        AddPath(assemblyDir, "assembly dir");
                }
            }
            catch
            {
            }

            // PRIORITY 4: NuGet runtime paths
            var rid = GetRuntimeIdentifier();
            var baseDirs = new List<string>();
            try
            {
                baseDirs.Add(Directory.GetCurrentDirectory());
            }
            catch
            {
            }

            try
            {
                baseDirs.Add(AppDomain.CurrentDomain.BaseDirectory);
            }
            catch
            {
            }

            foreach (var baseDir in baseDirs)
                if (!string.IsNullOrEmpty(baseDir))
                    AddPath(Path.Combine(baseDir, "runtimes", rid, "native"), $"NuGet runtime ({rid})");

            // PRIORITY 5: Environment variables
            var cudaPath = Environment.GetEnvironmentVariable("CUDA_PATH");
            if (!string.IsNullOrEmpty(cudaPath))
            {
                var cudaBinPath = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                    ? Path.Combine(cudaPath, "bin")
                    : Path.Combine(cudaPath, "lib64");
                AddPath(cudaBinPath, "CUDA_PATH env");
            }

            var mklRoot = Environment.GetEnvironmentVariable("MKLROOT");
            if (!string.IsNullOrEmpty(mklRoot))
            {
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                {
                    AddPath(Path.Combine(mklRoot, "redist", "intel64"), "MKLROOT env");
                    AddPath(Path.Combine(mklRoot, "bin", "intel64"), "MKLROOT env/bin");
                }
                else
                {
                    AddPath(Path.Combine(mklRoot, "lib", "intel64"), "MKLROOT env");
                }
            }

            // PRIORITY 6: Comprehensive system directory scanning
            AddSystemPaths(paths, checkedPaths, AddPath);

            // PRIORITY 7: PATH/LD_LIBRARY_PATH environment variable
            var pathEnv = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                ? Environment.GetEnvironmentVariable("PATH")
                : Environment.GetEnvironmentVariable("LD_LIBRARY_PATH");

            if (!string.IsNullOrEmpty(pathEnv))
            {
                var separator = RuntimeInformation.IsOSPlatform(OSPlatform.Windows) ? ';' : ':';
                var pathDirs = pathEnv.Split(separator, StringSplitOptions.RemoveEmptyEntries);

                foreach (var dir in pathDirs.Take(100)) AddPath(dir, "From PATH/LD_LIBRARY_PATH");
            }

            Debug.WriteLine($"[NativeLibrary] Total search paths: {paths.Count}");
            return paths.ToArray();
        }

        private void AddSystemPaths(List<string> paths, HashSet<string> checkedPaths, Action<string, string> addPath)
        {
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                AddWindowsSystemPaths(addPath);
            else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                AddLinuxSystemPaths(addPath);
            else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX)) AddMacOSSystemPaths(addPath);
        }

        private void AddWindowsSystemPaths(Action<string, string> addPath)
        {
            // ============================================================================
            // ULTIMATE FUTURE-PROOF PATH DISCOVERY FOR WINDOWS
            // ============================================================================
            // Same comprehensive approach as Linux, adapted for Windows conventions
            // ============================================================================

            // 1. STANDARD WINDOWS SYSTEM PATHS
            addPath(@"C:\Windows\System32", "Windows System32");
            addPath(@"C:\Windows\SysWOW64", "Windows SysWOW64 (32-bit)");

            // 2. PROGRAM FILES - ALL VARIANTS (handles different Windows versions/locales)
            var programFilesPaths = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

            try
            {
                programFilesPaths.Add(Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles));
            }
            catch
            {
            }

            try
            {
                programFilesPaths.Add(Environment.GetFolderPath(Environment.SpecialFolder.ProgramFilesX86));
            }
            catch
            {
            }

            // Fallback hardcoded paths
            programFilesPaths.Add(@"C:\Program Files");
            programFilesPaths.Add(@"C:\Program Files (x86)");

            // Also check other drives (D:, E:, etc.) - some enterprise installs use these
            for (var drive = 'D'; drive <= 'E'; drive++)
            {
                programFilesPaths.Add($@"{drive}:\Program Files");
                programFilesPaths.Add($@"{drive}:\Program Files (x86)");
            }

            // 3. DEEP SCAN: Program Files
            foreach (var programFiles in programFilesPaths)
            {
                if (string.IsNullOrEmpty(programFiles) || !Directory.Exists(programFiles))
                    continue;

                // Scan for NVIDIA/CUDA
                ScanForCUDAWindows(programFiles, addPath);

                // Scan for Intel/MKL
                ScanForMKLWindows(programFiles, addPath);

                // Deep scan for ANY library directories
                // PERFORMANCE FIX: Only do deep scanning when explicitly enabled
                // Deep scanning can cause massive IO spikes on startup (scanning Program Files recursively)
                if (EnableDeepScanning) ScanForLibraryDirectories(programFiles, 3, addPath, "Program Files");
            }

            // 4. CUDA DISCOVERY - Check ALL possible Windows CUDA locations
            var cudaBasePaths = new[]
            {
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA",
                @"C:\Program Files (x86)\NVIDIA GPU Computing Toolkit\CUDA",
                @"C:\NVIDIA\CUDA",
                @"C:\cuda",
                @"D:\NVIDIA\CUDA",
                @"D:\cuda"
            };

            foreach (var basePath in cudaBasePaths)
                if (Directory.Exists(basePath))
                    try
                    {
                        // Find all version directories
                        var versionDirs = Directory.GetDirectories(basePath)
                            .OrderByDescending(d => ExtractVersionFromPath(d));

                        foreach (var versionDir in versionDirs)
                            // Check multiple possible subdirectories
                        foreach (var subdir in new[] { "bin", "lib\\x64", "lib", "bin\\x64" })
                        {
                            var fullPath = Path.Combine(versionDir, subdir);
                            if (Directory.Exists(fullPath))
                                addPath(fullPath, $"CUDA {Path.GetFileName(versionDir)}");
                        }
                    }
                    catch
                    {
                    }

            // 5. MKL DISCOVERY - Check ALL possible Windows MKL locations
            var mklBasePaths = new[]
            {
                @"C:\Program Files\Intel\oneAPI\mkl",
                @"C:\Program Files (x86)\Intel\oneAPI\mkl",
                @"C:\Program Files\Intel\MKL",
                @"C:\Intel\oneAPI\mkl",
                @"C:\Intel\MKL"
            };

            foreach (var basePath in mklBasePaths)
                if (Directory.Exists(basePath))
                    try
                    {
                        var versionDirs = Directory.GetDirectories(basePath)
                            .OrderByDescending(d => ExtractVersionFromPath(d));

                        foreach (var versionDir in versionDirs)
                        foreach (var subdir in new[]
                                     { @"redist\intel64", @"bin\intel64", @"lib\intel64", "bin", "lib" })
                        {
                            var fullPath = Path.Combine(versionDir, subdir);
                            if (Directory.Exists(fullPath))
                                addPath(fullPath, $"MKL {Path.GetFileName(versionDir)}");
                        }
                    }
                    catch
                    {
                    }

            // 6. ENVIRONMENT VARIABLE PATHS
            var envVars = new[]
            {
                "CUDA_PATH",
                "CUDA_HOME",
                "MKL_ROOT",
                "MKLROOT",
                "INTEL_ROOT"
            };

            foreach (var envVar in envVars)
            {
                var value = Environment.GetEnvironmentVariable(envVar);
                if (!string.IsNullOrEmpty(value) && Directory.Exists(value))
                    foreach (var subdir in new[] { "bin", @"lib\x64", "lib", @"bin\x64" })
                    {
                        var fullPath = Path.Combine(value, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, $"Env: {envVar}");
                    }
            }

            // 7. PATH ENVIRONMENT VARIABLE
            // Parse PATH and add directories that might contain our libraries
            var pathEnv = Environment.GetEnvironmentVariable("PATH");
            if (!string.IsNullOrEmpty(pathEnv))
            {
                var paths = pathEnv.Split(';', StringSplitOptions.RemoveEmptyEntries);
                foreach (var path in paths)
                    try
                    {
                        if (Directory.Exists(path))
                        {
                            var lowerPath = path.ToLowerInvariant();
                            // Only add paths that look like they contain vendor libraries
                            if (lowerPath.Contains("cuda") || lowerPath.Contains("nvidia") ||
                                lowerPath.Contains("mkl") || lowerPath.Contains("intel"))
                                addPath(path, "PATH env");
                        }
                    }
                    catch
                    {
                    }
            }

            // 8. NUGET PACKAGE PATHS
            // Check for libraries installed via NuGet (common for Intel MKL)
            // PERFORMANCE FIX: Only do deep scanning when explicitly enabled
            if (EnableDeepScanning)
            {
                var nugetPaths = new[]
                {
                    Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), @".nuget\packages"),
                    @"C:\Program Files\dotnet\shared"
                };

                foreach (var nugetBase in nugetPaths)
                    if (Directory.Exists(nugetBase))
                        try
                        {
                            ScanForLibraryDirectories(nugetBase, 4, addPath, "NuGet");
                        }
                        catch
                        {
                        }
            }
        }

        private void ScanForCUDAWindows(string programFiles, Action<string, string> addPath)
        {
            try
            {
                var nvidiaBase = Path.Combine(programFiles, "NVIDIA GPU Computing Toolkit");
                if (Directory.Exists(nvidiaBase))
                {
                    var cudaBase = Path.Combine(nvidiaBase, "CUDA");
                    if (Directory.Exists(cudaBase))
                    {
                        var versionDirs = Directory.GetDirectories(cudaBase)
                            .Where(d => Path.GetFileName(d).StartsWith("v", StringComparison.OrdinalIgnoreCase))
                            .OrderByDescending(d => ExtractVersionFromPath(d));

                        foreach (var versionDir in versionDirs)
                            addPath(Path.Combine(versionDir, "bin"), $"CUDA {Path.GetFileName(versionDir)}");
                    }
                }

                var nvidiaDirs = Directory.GetDirectories(programFiles, "NVIDIA*", SearchOption.TopDirectoryOnly);
                foreach (var nvidiaDir in nvidiaDirs) ScanDirectoryForLibraries(nvidiaDir, "cudart*.dll", 2, addPath);
            }
            catch
            {
            }
        }

        private void ScanForMKLWindows(string programFiles, Action<string, string> addPath)
        {
            try
            {
                var intelBase = Path.Combine(programFiles, "Intel");
                if (Directory.Exists(intelBase))
                {
                    var oneAPIBase = Path.Combine(intelBase, "oneAPI");
                    if (Directory.Exists(oneAPIBase))
                    {
                        var mklBase = Path.Combine(oneAPIBase, "mkl");
                        if (Directory.Exists(mklBase))
                        {
                            var versionDirs = Directory.GetDirectories(mklBase)
                                .OrderByDescending(d => ExtractVersionFromPath(d));

                            foreach (var versionDir in versionDirs)
                            {
                                addPath(Path.Combine(versionDir, "redist", "intel64"),
                                    $"MKL {Path.GetFileName(versionDir)}");
                                addPath(Path.Combine(versionDir, "bin", "intel64"),
                                    $"MKL {Path.GetFileName(versionDir)} bin");
                            }
                        }
                    }

                    ScanDirectoryForLibraries(intelBase, "mkl_rt*.dll", 4, addPath);
                }
            }
            catch
            {
            }
        }

        private void AddLinuxSystemPaths(Action<string, string> addPath)
        {
            // ============================================================================
            // ULTIMATE FUTURE-PROOF PATH DISCOVERY FOR LINUX
            // ============================================================================
            // This scans ALL possible locations where libraries could be installed,
            // ensuring compatibility for years to come regardless of distribution,
            // package manager, or installation method.
            // ============================================================================

            // 1. STANDARD SYSTEM LIBRARY PATHS (always check these first)
            var standardPaths = new[]
            {
                "/lib", // Base system libraries
                "/lib64", // 64-bit system libraries
                "/usr/lib", // User-space libraries
                "/usr/lib64", // 64-bit user libraries
                "/usr/local/lib", // Locally installed libraries
                "/usr/local/lib64" // 64-bit local libraries
            };

            foreach (var path in standardPaths)
                if (Directory.Exists(path))
                    addPath(path, "Linux system lib");

            // 2. MULTIARCH PATHS (Debian/Ubuntu convention - all possible variants)
            var multarchSuffixes = new[]
            {
                "x86_64-linux-gnu", // 64-bit x86
                "aarch64-linux-gnu", // ARM64
                "i386-linux-gnu", // 32-bit x86 (legacy)
                "powerpc64le-linux-gnu", // IBM POWER
                "s390x-linux-gnu" // IBM Z
            };

            foreach (var suffix in multarchSuffixes)
            {
                var paths = new[]
                {
                    $"/lib/{suffix}",
                    $"/usr/lib/{suffix}",
                    $"/usr/local/lib/{suffix}"
                };

                foreach (var path in paths)
                    if (Directory.Exists(path))
                        addPath(path, $"Multiarch ({suffix})");
            }

            // 3. DEEP SCAN: /opt (vendors often install here)
            // Scan up to 4 levels deep to find lib/lib64 directories
            // PERFORMANCE FIX: Only do deep scanning when explicitly enabled
            if (EnableDeepScanning) ScanForLibraryDirectories("/opt", 4, addPath, "opt");

            // 4. DEEP SCAN: /usr/local (local installations)
            // PERFORMANCE FIX: Only do deep scanning when explicitly enabled
            if (EnableDeepScanning) ScanForLibraryDirectories("/usr/local", 3, addPath, "usr/local");

            // 5. CUDA DISCOVERY - Check ALL possible CUDA locations
            var cudaBasePaths = new[]
            {
                "/usr/local/cuda", // Standard NVIDIA installer location
                "/opt/cuda", // Alternative location
                "/usr/cuda", // Another alternative
                "/opt/nvidia/cuda" // Enterprise installations
            };

            foreach (var basePath in cudaBasePaths)
                if (Directory.Exists(basePath))
                    // Check standard subdirectories
                    foreach (var subdir in new[] { "lib64", "lib", "lib/x64", "targets/x86_64-linux/lib" })
                    {
                        var fullPath = Path.Combine(basePath, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, $"CUDA {basePath}");
                    }

            // Scan for versioned CUDA directories (cuda-12.6, cuda-12.5, etc.)
            try
            {
                var searchDirs = new[] { "/usr/local", "/opt", "/usr" };
                foreach (var searchDir in searchDirs)
                {
                    if (!Directory.Exists(searchDir)) continue;

                    var cudaDirs = Directory.GetDirectories(searchDir, "cuda*", SearchOption.TopDirectoryOnly)
                        .OrderByDescending(d => ExtractVersionFromPath(d));

                    foreach (var cudaDir in cudaDirs)
                    foreach (var subdir in new[] { "lib64", "lib", "targets/x86_64-linux/lib" })
                    {
                        var fullPath = Path.Combine(cudaDir, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, $"CUDA {Path.GetFileName(cudaDir)}");
                    }
                }
            }
            catch
            {
            }

            // 6. INTEL MKL DISCOVERY - Check ALL possible MKL locations
            var mklBasePaths = new[]
            {
                "/opt/intel/oneapi/mkl", // oneAPI installation
                "/opt/intel/mkl", // Standalone MKL
                "/opt/intel/compilers_and_libraries/linux/mkl", // Legacy Parallel Studio
                "/usr/local/intel/mkl", // Custom installation
                "/opt/intel/compilers_and_libraries_*/linux/mkl" // Pattern for versions
            };

            foreach (var pattern in mklBasePaths)
                // Handle wildcards
                if (pattern.Contains("*"))
                    try
                    {
                        var baseDir = Path.GetDirectoryName(pattern.Replace("*", ""));
                        if (Directory.Exists(baseDir))
                        {
                            var searchPattern = Path.GetFileName(pattern);
                            var dirs = Directory.GetDirectories(baseDir, searchPattern.Replace("*", "*"));
                            foreach (var dir in dirs)
                            {
                                var libPath = Path.Combine(dir, "lib/intel64");
                                if (Directory.Exists(libPath))
                                    addPath(libPath, "MKL Legacy");
                            }
                        }
                    }
                    catch
                    {
                    }
                else if (Directory.Exists(pattern))
                    // Check multiple possible lib subdirectories
                    foreach (var subdir in new[] { "latest/lib/intel64", "latest/lib", "lib/intel64", "lib" })
                    {
                        var fullPath = Path.Combine(pattern, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, "MKL Linux");
                    }

            // Deep scan /opt/intel for ANY MKL libraries
            ScanDirectoryForLibraries("/opt/intel", "libmkl_rt.so*", 5, addPath);

            // 7. LDCONFIG QUERY - Ask the system what it knows about
            // This catches anything we might have missed
            QueryLdconfigForLibraries(addPath);

            // 8. ENVIRONMENT VARIABLE PATHS
            // LD_LIBRARY_PATH and other environment variables are handled by GetSearchPaths()
            // but let's also check for specific vendor variables
            var envVars = new[]
            {
                "CUDA_HOME",
                "CUDA_PATH",
                "MKL_ROOT",
                "MKLROOT",
                "INTEL_ROOT",
                "PARDISO_LIC_PATH"
            };

            foreach (var envVar in envVars)
            {
                var value = Environment.GetEnvironmentVariable(envVar);
                if (!string.IsNullOrEmpty(value) && Directory.Exists(value))
                    foreach (var subdir in new[] { "lib64", "lib", "lib/intel64", "bin" })
                    {
                        var fullPath = Path.Combine(value, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, $"Env: {envVar}");
                    }
            }
        }

        private void AddMacOSSystemPaths(Action<string, string> addPath)
        {
            // ============================================================================
            // ULTIMATE FUTURE-PROOF PATH DISCOVERY FOR MACOS
            // ============================================================================
            // Same comprehensive approach as Linux, adapted for macOS conventions
            // ============================================================================

            // 1. STANDARD MACOS SYSTEM PATHS
            var standardPaths = new[]
            {
                "/usr/lib",
                "/usr/local/lib",
                "/opt/local/lib", // MacPorts
                "/opt/homebrew/lib", // Homebrew on Apple Silicon
                "/usr/local/opt", // Homebrew on Intel
                "/Library/Frameworks" // System frameworks
            };

            foreach (var path in standardPaths)
                if (Directory.Exists(path))
                    addPath(path, "macOS system lib");

            // 2. HOMEBREW DISCOVERY (Intel vs Apple Silicon)
            var homebrewPaths = new[]
            {
                "/opt/homebrew", // Apple Silicon (M1/M2/M3)
                "/usr/local/Homebrew", // Intel Mac
                "/usr/local" // Intel Mac (legacy)
            };

            foreach (var brewBase in homebrewPaths)
                if (Directory.Exists(brewBase))
                    // Check standard Homebrew library locations
                    foreach (var subdir in new[] { "lib", "opt", "Cellar" })
                    {
                        var fullPath = Path.Combine(brewBase, subdir);
                        if (Directory.Exists(fullPath))
                        {
                            addPath(fullPath, $"Homebrew {subdir}");

                            // Also scan Cellar for versioned packages
                            // PERFORMANCE FIX: Only do deep scanning when explicitly enabled
                            if (subdir == "Cellar" && EnableDeepScanning)
                                ScanForLibraryDirectories(fullPath, 3, addPath, "Homebrew Cellar");
                        }
                    }

            // 3. DEEP SCAN: /opt (vendor installations)
            // PERFORMANCE FIX: Only do deep scanning when explicitly enabled
            if (EnableDeepScanning) ScanForLibraryDirectories("/opt", 4, addPath, "opt");

            // 4. DEEP SCAN: /usr/local
            // PERFORMANCE FIX: Only do deep scanning when explicitly enabled
            if (EnableDeepScanning) ScanForLibraryDirectories("/usr/local", 3, addPath, "usr/local");

            // 5. CUDA DISCOVERY - Check ALL possible macOS CUDA locations
            // Note: CUDA on macOS is deprecated after macOS 10.13, but check anyway for legacy systems
            var cudaBasePaths = new[]
            {
                "/Developer/NVIDIA/CUDA-*",
                "/usr/local/cuda*",
                "/opt/cuda",
                "/Library/Frameworks/CUDA.framework"
            };

            foreach (var pattern in cudaBasePaths)
                try
                {
                    if (pattern.Contains("*"))
                    {
                        var baseDir = Path.GetDirectoryName(pattern.Replace("*", "")) ?? "";
                        if (Directory.Exists(baseDir))
                        {
                            var searchPattern = Path.GetFileName(pattern);
                            var dirs = Directory.GetDirectories(baseDir, searchPattern);
                            foreach (var dir in dirs.OrderByDescending(d => ExtractVersionFromPath(d)))
                            foreach (var subdir in new[] { "lib", "lib64" })
                            {
                                var fullPath = Path.Combine(dir, subdir);
                                if (Directory.Exists(fullPath))
                                    addPath(fullPath, $"CUDA {Path.GetFileName(dir)}");
                            }
                        }
                    }
                    else if (Directory.Exists(pattern))
                    {
                        foreach (var subdir in new[] { "lib", "lib64", "Libraries" })
                        {
                            var fullPath = Path.Combine(pattern, subdir);
                            if (Directory.Exists(fullPath))
                                addPath(fullPath, "CUDA macOS");
                        }
                    }
                }
                catch
                {
                }

            // 6. MKL DISCOVERY - Check ALL possible macOS MKL locations
            var mklBasePaths = new[]
            {
                "/opt/intel/oneapi/mkl",
                "/opt/intel/mkl",
                "/opt/intel/compilers_and_libraries/mac/mkl",
                "/usr/local/intel/mkl",
                "/Library/Frameworks/Intel_MKL.framework"
            };

            foreach (var basePath in mklBasePaths)
                if (Directory.Exists(basePath))
                    try
                    {
                        // Check for versioned directories
                        var versionDirs = Directory.GetDirectories(basePath)
                            .OrderByDescending(d => ExtractVersionFromPath(d));

                        if (versionDirs.Any())
                            foreach (var versionDir in versionDirs)
                            foreach (var subdir in new[] { "lib", "lib/intel64", "Libraries" })
                            {
                                var fullPath = Path.Combine(versionDir, subdir);
                                if (Directory.Exists(fullPath))
                                    addPath(fullPath, $"MKL {Path.GetFileName(versionDir)}");
                            }
                        else
                            // No version subdirectories, check for lib directly
                            foreach (var subdir in new[] { "latest/lib", "lib", "lib/intel64", "Libraries" })
                            {
                                var fullPath = Path.Combine(basePath, subdir);
                                if (Directory.Exists(fullPath))
                                    addPath(fullPath, "MKL macOS");
                            }
                    }
                    catch
                    {
                    }

            // Deep scan /opt/intel for ANY MKL libraries
            if (Directory.Exists("/opt/intel"))
                ScanDirectoryForLibraries("/opt/intel", "libmkl_rt*.dylib", 5, addPath);

            // 7. ENVIRONMENT VARIABLE PATHS
            var envVars = new[]
            {
                "CUDA_HOME",
                "CUDA_PATH",
                "MKL_ROOT",
                "MKLROOT",
                "INTEL_ROOT"
            };

            foreach (var envVar in envVars)
            {
                var value = Environment.GetEnvironmentVariable(envVar);
                if (!string.IsNullOrEmpty(value) && Directory.Exists(value))
                    foreach (var subdir in new[] { "lib", "lib64", "lib/intel64", "Libraries" })
                    {
                        var fullPath = Path.Combine(value, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, $"Env: {envVar}");
                    }
            }

            // 8. DYLD_LIBRARY_PATH and DYLD_FALLBACK_LIBRARY_PATH
            var dyldVars = new[] { "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH" };
            foreach (var dyldVar in dyldVars)
            {
                var value = Environment.GetEnvironmentVariable(dyldVar);
                if (!string.IsNullOrEmpty(value))
                {
                    var paths = value.Split(':', StringSplitOptions.RemoveEmptyEntries);
                    foreach (var path in paths)
                        if (Directory.Exists(path))
                            addPath(path, $"Env: {dyldVar}");
                }
            }

            // 9. FRAMEWORK SEARCH
            // macOS uses Frameworks for some libraries
            var frameworkPaths = new[]
            {
                "/Library/Frameworks",
                "/System/Library/Frameworks",
                Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "Library/Frameworks")
            };

            foreach (var frameworkPath in frameworkPaths)
                if (Directory.Exists(frameworkPath))
                    try
                    {
                        var frameworks = Directory.GetDirectories(frameworkPath, "*.framework");
                        foreach (var framework in frameworks)
                        {
                            var name = Path.GetFileName(framework);
                            // Only check frameworks that might contain our libraries
                            if (name.Contains("CUDA", StringComparison.OrdinalIgnoreCase) ||
                                name.Contains("MKL", StringComparison.OrdinalIgnoreCase) ||
                                name.Contains("Intel", StringComparison.OrdinalIgnoreCase))
                            {
                                var libPath = Path.Combine(framework, "Libraries");
                                if (Directory.Exists(libPath))
                                    addPath(libPath, $"Framework: {name}");
                            }
                        }
                    }
                    catch
                    {
                    }
        }

        private void ScanDirectoryForLibraries(string baseDir, string filePattern, int maxDepth,
            Action<string, string> addPath)
        {
            if (maxDepth <= 0 || !Directory.Exists(baseDir))
                return;

            try
            {
                var matchingFiles = Directory.GetFiles(baseDir, filePattern, SearchOption.TopDirectoryOnly);
                if (matchingFiles.Length > 0)
                {
                    addPath(baseDir, $"Found {filePattern}");
                    return;
                }

                var subdirs = Directory.GetDirectories(baseDir);
                foreach (var subdir in subdirs)
                    ScanDirectoryForLibraries(subdir, filePattern, maxDepth - 1, addPath);
            }
            catch
            {
            }
        }

        /// <summary>
        ///     Recursively scans for lib, lib64, and library directories.
        ///     This finds libraries regardless of where they're installed.
        /// </summary>
        private void ScanForLibraryDirectories(string baseDir, int maxDepth, Action<string, string> addPath,
            string source)
        {
            if (maxDepth <= 0 || !Directory.Exists(baseDir))
                return;

            try
            {
                // Check if current directory contains library files
                var hasLibraries = Directory.GetFiles(baseDir, "*.so*", SearchOption.TopDirectoryOnly).Length > 0 ||
                                   Directory.GetFiles(baseDir, "*.a", SearchOption.TopDirectoryOnly).Length > 0;

                if (hasLibraries) addPath(baseDir, source);

                // Look for standard library directory names
                var libDirNames = new[] { "lib", "lib64", "lib32", "libs", "library" };
                foreach (var libDirName in libDirNames)
                {
                    var libPath = Path.Combine(baseDir, libDirName);
                    if (Directory.Exists(libPath))
                    {
                        addPath(libPath, $"{source}/{libDirName}");

                        // Also check subdirectories like lib64/intel64
                        try
                        {
                            var libSubdirs = Directory.GetDirectories(libPath);
                            foreach (var subdir in libSubdirs)
                            {
                                var subdirHasLibs = Directory.GetFiles(subdir, "*.so*", SearchOption.TopDirectoryOnly)
                                    .Length > 0;
                                if (subdirHasLibs)
                                    addPath(subdir, $"{source}/{libDirName}/{Path.GetFileName(subdir)}");
                            }
                        }
                        catch
                        {
                        }
                    }
                }

                // Recurse into subdirectories (but skip common non-library directories)
                var skipDirs = new HashSet<string>(StringComparer.OrdinalIgnoreCase)
                {
                    "bin", "sbin", "include", "share", "doc", "docs", "man",
                    "src", "source", "samples", "examples", "test", "tests"
                };

                var subdirs = Directory.GetDirectories(baseDir);
                foreach (var subdir in subdirs)
                {
                    var dirName = Path.GetFileName(subdir);
                    if (!skipDirs.Contains(dirName))
                        ScanForLibraryDirectories(subdir, maxDepth - 1, addPath, source);
                }
            }
            catch
            {
            }
        }

        /// <summary>
        ///     Queries ldconfig to find libraries already registered with the system.
        ///     This catches anything we might have missed in directory scanning.
        /// </summary>
        private void QueryLdconfigForLibraries(Action<string, string> addPath)
        {
            if (!RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                return;

            try
            {
                var startInfo = new ProcessStartInfo
                {
                    FileName = "ldconfig",
                    Arguments = "-p",
                    RedirectStandardOutput = true,
                    RedirectStandardError = true,
                    UseShellExecute = false,
                    CreateNoWindow = true
                };

                using var process = Process.Start(startInfo);
                if (process == null)
                    return;

                var output = process.StandardOutput.ReadToEnd();
                process.WaitForExit(2000); // 2 second timeout

                if (process.ExitCode != 0)
                    return;

                // Parse ldconfig output for library paths
                // Format: "libname.so.version (arch) => /path/to/lib"
                var seenPaths = new HashSet<string>();
                var lines = output.Split('\n');

                foreach (var line in lines)
                    if (line.Contains("=>"))
                    {
                        var parts = line.Split(new[] { "=>" }, StringSplitOptions.None);
                        if (parts.Length == 2)
                        {
                            var libPath = parts[1].Trim();
                            if (!string.IsNullOrEmpty(libPath) && File.Exists(libPath))
                            {
                                var directory = Path.GetDirectoryName(libPath);
                                if (!string.IsNullOrEmpty(directory) && !seenPaths.Contains(directory))
                                {
                                    seenPaths.Add(directory);

                                    // Check if this looks like a CUDA, MKL, or other vendor library
                                    var lowerPath = directory.ToLowerInvariant();
                                    var source = "ldconfig";

                                    if (lowerPath.Contains("cuda"))
                                        source = "ldconfig (CUDA)";
                                    else if (lowerPath.Contains("mkl") || lowerPath.Contains("intel"))
                                        source = "ldconfig (MKL)";
                                    else if (lowerPath.Contains("nvidia"))
                                        source = "ldconfig (NVIDIA)";

                                    addPath(directory, source);
                                }
                            }
                        }
                    }
            }
            catch
            {
            }
        }

        private string[]? DiscoverLibrariesInSystem(Regex pattern, string libraryType)
        {
            var searchPaths = GetSearchPaths();
            if (searchPaths == null || searchPaths.Length == 0)
                return null;

            var discoveredLibraries = new List<LibraryInfo>();

            foreach (var searchPath in searchPaths)
                try
                {
                    if (!Directory.Exists(searchPath))
                        continue;

                    // Get all files that match the pattern
                    var files = Directory.GetFiles(searchPath, "*.*", SearchOption.TopDirectoryOnly)
                        .Where(f => pattern.IsMatch(Path.GetFileName(f)));

                    foreach (var file in files)
                    {
                        var fileName = Path.GetFileName(file);
                        var fileVersion = ExtractVersionFromFileName(fileName);
                        var dirVersion = ExtractVersionFromPath(searchPath);

                        // CRITICAL: Resolve symlinks to find the REAL version!
                        var resolvedPath = ResolveSymlink(file);
                        var isSymlink = resolvedPath != null && resolvedPath != file;
                        var resolvedVersion = isSymlink && resolvedPath != null
                            ? ExtractVersionFromFileName(Path.GetFileName(resolvedPath))
                            : fileVersion;

                        var info = new LibraryInfo
                        {
                            FullPath = file,
                            FileName = fileName,
                            DirectoryPath = searchPath,
                            FileVersion = fileVersion,
                            DirectoryVersion = dirVersion,
                            ResolvedPath = resolvedPath,
                            ResolvedVersion = resolvedVersion,
                            IsSymlink = isSymlink
                        };

                        discoveredLibraries.Add(info);
                    }

                    // CRITICAL FIX: Also look for versioned files without symlinks
                    // e.g., if user has libcudart.so.11.0 but no libcudart.so symlink
                    // This handles incomplete installations that PhD students might have
                    if (libraryType.Contains("CUDA") || libraryType.Contains("cuda"))
                    {
                        // Look for libcudart.so.*.* (any version)
                        var cudaFiles = Directory.GetFiles(searchPath, "libcudart.so.*");
                        foreach (var file in cudaFiles)
                        {
                            var fileName = Path.GetFileName(file);
                            // Add if not already in list
                            if (!discoveredLibraries.Any(lib => lib.FullPath == file))
                            {
                                var fileVersion = ExtractVersionFromFileName(fileName);
                                var dirVersion = ExtractVersionFromPath(searchPath);
                                var resolvedPath = ResolveSymlink(file);
                                var isSymlink = resolvedPath != null && resolvedPath != file;
                                var resolvedVersion = isSymlink && resolvedPath != null
                                    ? ExtractVersionFromFileName(Path.GetFileName(resolvedPath))
                                    : fileVersion;

                                discoveredLibraries.Add(new LibraryInfo
                                {
                                    FullPath = file,
                                    FileName = fileName,
                                    DirectoryPath = searchPath,
                                    FileVersion = fileVersion,
                                    DirectoryVersion = dirVersion,
                                    ResolvedPath = resolvedPath,
                                    ResolvedVersion = resolvedVersion,
                                    IsSymlink = isSymlink
                                });
                            }
                        }

                        // Also look for cusparse
                        var cusparseFiles = Directory.GetFiles(searchPath, "libcusparse.so.*");
                        foreach (var file in cusparseFiles)
                        {
                            var fileName = Path.GetFileName(file);
                            if (!discoveredLibraries.Any(lib => lib.FullPath == file))
                            {
                                var fileVersion = ExtractVersionFromFileName(fileName);
                                var dirVersion = ExtractVersionFromPath(searchPath);
                                var resolvedPath = ResolveSymlink(file);
                                var isSymlink = resolvedPath != null && resolvedPath != file;
                                var resolvedVersion = isSymlink && resolvedPath != null
                                    ? ExtractVersionFromFileName(Path.GetFileName(resolvedPath))
                                    : fileVersion;

                                discoveredLibraries.Add(new LibraryInfo
                                {
                                    FullPath = file,
                                    FileName = fileName,
                                    DirectoryPath = searchPath,
                                    FileVersion = fileVersion,
                                    DirectoryVersion = dirVersion,
                                    ResolvedPath = resolvedPath,
                                    ResolvedVersion = resolvedVersion,
                                    IsSymlink = isSymlink
                                });
                            }
                        }
                    }

                    // Same for MKL
                    if (libraryType.Contains("MKL") || libraryType.Contains("mkl"))
                    {
                        var mklFiles = Directory.GetFiles(searchPath, "libmkl_rt.so.*");
                        foreach (var file in mklFiles)
                        {
                            var fileName = Path.GetFileName(file);
                            if (!discoveredLibraries.Any(lib => lib.FullPath == file))
                            {
                                var fileVersion = ExtractVersionFromFileName(fileName);
                                var dirVersion = ExtractVersionFromPath(searchPath);
                                var resolvedPath = ResolveSymlink(file);
                                var isSymlink = resolvedPath != null && resolvedPath != file;
                                var resolvedVersion = isSymlink && resolvedPath != null
                                    ? ExtractVersionFromFileName(Path.GetFileName(resolvedPath))
                                    : fileVersion;

                                discoveredLibraries.Add(new LibraryInfo
                                {
                                    FullPath = file,
                                    FileName = fileName,
                                    DirectoryPath = searchPath,
                                    FileVersion = fileVersion,
                                    DirectoryVersion = dirVersion,
                                    ResolvedPath = resolvedPath,
                                    ResolvedVersion = resolvedVersion,
                                    IsSymlink = isSymlink
                                });
                            }
                        }
                    }
                }
                catch
                {
                }

            if (discoveredLibraries.Count == 0)
                return null;

            // CRITICAL: Sort by RESOLVED version (after following symlinks) to pick the LATEST!
            // Example: cuda/11.8/libcudart.so -> REDIST/cuda/12.6/libcudart.so.12.6.77
            // We want to use 12.6.77, not 11.8!
            var sorted = discoveredLibraries
                .OrderByDescending(lib => lib.ResolvedVersion) // HIGHEST version first!
                .ThenByDescending(lib => lib.DirectoryVersion) // Then by directory version
                .ThenByDescending(lib => lib.FileVersion) // Then by file version
                .ThenBy(lib => lib.FileName.Count(c => c == '.')) // Prefer simpler names (symlinks)
                .ThenBy(lib => lib.FileName)
                .ToList();

            Debug.WriteLine($"[NativeLibrary] Discovered {libraryType} libraries ({sorted.Count}):");
            foreach (var lib in sorted.Take(10))
            {
                var symlinkInfo = lib.IsSymlink ? $" -> {Path.GetFileName(lib.ResolvedPath ?? "")}" : "";
                Debug.WriteLine(
                    $"[NativeLibrary]   ✓ {lib.FileName} (resolved v{lib.ResolvedVersion:F5}){symlinkInfo}");
            }

            return sorted.Select(lib => lib.FileName).Distinct().ToArray();
        }

        private double ExtractVersionFromFileName(string fileName)
        {
            // Extract version from filenames like:
            // libcudart.so.12.6.77 -> 12.677
            // libcudart.so.11.0 -> 11.0
            // libmkl_rt.so.2 -> 2.0

            var match = Regex.Match(fileName, @"\.so\.(\d+)\.(\d+)\.(\d+)");
            if (match.Success)
            {
                // Full version: 12.6.77 -> 12.00677
                var major = int.Parse(match.Groups[1].Value);
                var minor = int.Parse(match.Groups[2].Value);
                var patch = int.Parse(match.Groups[3].Value);
                return major + minor / 100.0 + patch / 100000.0;
            }

            match = Regex.Match(fileName, @"\.so\.(\d+)\.(\d+)");
            if (match.Success)
            {
                // Major.minor: 12.6 -> 12.6
                var major = int.Parse(match.Groups[1].Value);
                var minor = int.Parse(match.Groups[2].Value);
                return major + minor / 100.0;
            }

            match = Regex.Match(fileName, @"\.so\.(\d+)");
            if (match.Success)
                // Major only: 12 -> 12.0
                return double.Parse(match.Groups[1].Value);

            match = Regex.Match(fileName, @"[\._](\d+(?:\.\d+)?)(?:[\._]|$)");
            if (match.Success && double.TryParse(match.Groups[1].Value, out var version))
                return version;

            return 0;
        }

        private string? ResolveSymlink(string path)
        {
            // Resolve symlinks to find the actual file
            // This is critical for finding the REAL version!
            // Example: cuda/11.8/lib/libcudart.so -> ../../REDIST/cuda/12.6/lib/libcudart.so.12.6.77

            if (!File.Exists(path))
                return null;

            try
            {
                // WINDOWS: Handle symlinks, junction points, and hard links
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                {
                    var fileInfo = new FileInfo(path);

                    // Check if it's a reparse point (symlink or junction)
                    if (fileInfo.Attributes.HasFlag(FileAttributes.ReparsePoint))
                    {
                        // Use LinkTarget (available in .NET 6+)
                        var target = fileInfo.LinkTarget;
                        if (!string.IsNullOrEmpty(target))
                        {
                            // Resolve relative paths
                            var resolvedPath = Path.IsPathRooted(target)
                                ? target
                                : Path.GetFullPath(Path.Combine(fileInfo.DirectoryName ?? "", target));

                            if (File.Exists(resolvedPath))
                                // Recursively resolve in case of chained symlinks
                                return ResolveSymlink(resolvedPath);
                        }
                    }

                    // If not a symlink or can't resolve, return original
                    return path;
                }

                // LINUX/MACOS: Use readlink command
                else
                {
                    var fileInfo = new FileInfo(path);

                    // Check if it's a symbolic link by checking attributes
                    if (fileInfo.Attributes.HasFlag(FileAttributes.ReparsePoint))
                    {
                        // Resolve the symlink using readlink -f (follows entire chain)
                        var startInfo = new ProcessStartInfo
                        {
                            FileName = "readlink",
                            Arguments = $"-f \"{path}\"",
                            RedirectStandardOutput = true,
                            UseShellExecute = false,
                            CreateNoWindow = true
                        };

                        using var process = Process.Start(startInfo);
                        if (process != null)
                        {
                            var resolved = process.StandardOutput.ReadToEnd().Trim();
                            // FIXED: Add timeout to prevent indefinite hang
                            if (!process.WaitForExit(5000)) // 5 second timeout for readlink
                            {
                                try
                                {
                                    process.Kill();
                                }
                                catch
                                {
                                }

                                return path;
                            }

                            if (!string.IsNullOrEmpty(resolved) && File.Exists(resolved))
                                return resolved;
                        }
                    }

                    // If not a symlink or can't resolve, return original
                    return path;
                }
            }
            catch
            {
                return path;
            }
        }

        private double ExtractVersionFromPath(string path)
        {
            var parts = path.Split(new[] { '/', '\\' }, StringSplitOptions.RemoveEmptyEntries);

            for (var i = parts.Length - 1; i >= 0; i--)
            {
                var part = parts[i];

                var match1 = Regex.Match(part, @"^v(\d+(?:\.\d+)?)$", RegexOptions.IgnoreCase);
                if (match1.Success && double.TryParse(match1.Groups[1].Value, out var ver1))
                    return ver1;

                var match2 = Regex.Match(part, @"^cuda-(\d+(?:\.\d+)?)$", RegexOptions.IgnoreCase);
                if (match2.Success && double.TryParse(match2.Groups[1].Value, out var ver2))
                    return ver2;

                var match3 = Regex.Match(part, @"^(\d+(?:\.\d+)?)$");
                if (match3.Success && double.TryParse(match3.Groups[1].Value, out var ver3))
                    return ver3;
            }

            return 0;
        }

        private static string GetRuntimeIdentifier()
        {
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                return RuntimeInformation.ProcessArchitecture == Architecture.X64 ? "win-x64" : "win-arm64";
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                return RuntimeInformation.ProcessArchitecture == Architecture.X64 ? "linux-x64" : "linux-arm64";
            if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                return RuntimeInformation.ProcessArchitecture == Architecture.X64 ? "osx-x64" : "osx-arm64";
            return "unknown";
        }

        /// <summary>
        ///     SECURITY FIX: Chocolatey auto-installation has been disabled.
        ///     The previous implementation downloaded and executed a PowerShell script from the internet
        ///     using iex (Invoke-Expression), which is a Remote Code Execution (RCE) vulnerability.
        ///     An attacker controlling the network (DNS spoofing/MITM) could execute arbitrary code
        ///     with Administrator privileges on the host machine.
        ///     Instead, this method now returns false and logs instructions for manual installation.
        /// </summary>
        /// <returns>Always returns false. Users must install Chocolatey manually if desired.</returns>
        private bool TryInstallChocolatey()
        {
            // SECURITY: Do NOT auto-install from remote URLs
            // The previous implementation was vulnerable to MITM/DNS spoofing attacks

            Debug.WriteLine("[NativeLibrary] SECURITY: Auto-installation of Chocolatey has been disabled.");
            Debug.WriteLine("[NativeLibrary] To install Chocolatey manually, visit: https://chocolatey.org/install");
            Debug.WriteLine("[NativeLibrary] After installation, run: choco install intel-mkl");

            if (InteractiveInstall)
            {
                Console.WriteLine();
                Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                Console.WriteLine("║  Manual Installation Required                                  ║");
                Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                Console.WriteLine();
                Console.WriteLine("Automatic Chocolatey installation has been disabled for security reasons.");
                Console.WriteLine();
                Console.WriteLine("To install Chocolatey manually:");
                Console.WriteLine("  1. Open PowerShell as Administrator");
                Console.WriteLine("  2. Visit https://chocolatey.org/install for the official installation script");
                Console.WriteLine("  3. After installing Chocolatey, run: choco install intel-mkl");
                Console.WriteLine();
                Console.WriteLine("Alternative: Install MKL via NuGet:");
                Console.WriteLine("  dotnet add package Intel.MKL.redist.win");
                Console.WriteLine();
            }

            return false;
        }

        private void RefreshEnvironmentPath()
        {
            // Refresh PATH environment variable to include newly installed tools
            var path = Environment.GetEnvironmentVariable("PATH", EnvironmentVariableTarget.Machine) ?? "";
            var userPath = Environment.GetEnvironmentVariable("PATH", EnvironmentVariableTarget.User) ?? "";
            Environment.SetEnvironmentVariable("PATH", path + ";" + userPath, EnvironmentVariableTarget.Process);
        }

        private bool TryAutoInstallMKL()
        {
            try
            {
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                    return HandleWindowsMKLInstall();
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                    return HandleLinuxMKLInstall();
                if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                    return HandleMacOSMKLInstall();

                return false;
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[NativeLibrary] Auto-install error: {ex.Message}");
                return false;
            }
        }

        private bool HandleWindowsMKLInstall()
        {
            // Windows: Try NuGet first, then Chocolatey (auto-installing if needed)

            // Check if dotnet CLI is available for NuGet install
            var hasDotnet = CommandExists("dotnet");
            var hasChoco = CommandExists("choco");

            // If neither dotnet nor Chocolatey exists, offer to install Chocolatey
            if (!hasDotnet && !hasChoco)
            {
                if (InteractiveInstall)
                {
                    Console.WriteLine();
                    Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                    Console.WriteLine("║  Package Manager Not Found                                     ║");
                    Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                    Console.WriteLine();
                    Console.WriteLine("To auto-install MKL, we need a package manager.");
                    Console.WriteLine();
                    Console.WriteLine("Would you like to install Chocolatey (Windows package manager)?");
                    Console.WriteLine("  This will enable automatic MKL installation.");
                    Console.WriteLine("  Requires: Administrator privileges");
                    Console.WriteLine("  Time: ~30 seconds");
                    Console.WriteLine();
                    Console.Write("Install Chocolatey now? (y/n): ");

                    var response = Console.ReadLine()?.Trim().ToLower();
                    if (response == "y" || response == "yes")
                    {
                        Console.WriteLine();
                        Console.WriteLine("Installing Chocolatey...");

                        if (TryInstallChocolatey())
                        {
                            Console.WriteLine("✓ Chocolatey installed successfully!");
                            hasChoco = true;
                            // Continue to MKL installation below
                        }
                        else
                        {
                            Console.WriteLine("✗ Chocolatey installation failed.");
                            Console.WriteLine();
                            Console.WriteLine("Please install Intel MKL manually:");
                            Console.WriteLine(
                                "  Option 1: Install .NET SDK, then: dotnet add package Intel.MKL.redist.win");
                            Console.WriteLine("  Option 2: Install Chocolatey, then: choco install intel-mkl");
                            Console.WriteLine("  Option 3: Download from intel.com");
                            Console.WriteLine();
                            return false;
                        }
                    }
                    else
                    {
                        Console.WriteLine();
                        Console.WriteLine("Skipping Chocolatey installation.");
                        Console.WriteLine();
                        Console.WriteLine("Please install Intel MKL manually:");
                        Console.WriteLine(
                            "  Option 1: Install .NET SDK, then: dotnet add package Intel.MKL.redist.win");
                        Console.WriteLine("  Option 2: Download from intel.com");
                        Console.WriteLine();
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }

            if (InteractiveInstall)
            {
                Console.WriteLine();
                Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                Console.WriteLine("║  Intel MKL Not Found - Auto-Install Available                 ║");
                Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                Console.WriteLine();

                if (hasDotnet)
                {
                    Console.WriteLine("Option 1: Install via NuGet (recommended)");
                    Console.WriteLine("  Command: dotnet add package Intel.MKL.redist.win");
                }

                if (hasChoco)
                {
                    Console.WriteLine("Option 2: Install via Chocolatey");
                    Console.WriteLine("  Command: choco install intel-mkl");
                }

                Console.WriteLine();
                Console.Write("Install Intel MKL now? (y/n): ");

                var response = Console.ReadLine()?.Trim().ToLower();
                if (response != "y" && response != "yes")
                {
                    Console.WriteLine();
                    Console.WriteLine("Skipping installation.");
                    return false;
                }

                Console.WriteLine();
            }

            // Try NuGet first (cleaner, project-local install)
            if (hasDotnet)
            {
                if (InteractiveInstall)
                    Console.WriteLine("Attempting NuGet installation...");

                var csprojPath = FindProjectFile();
                if (csprojPath != null)
                {
                    var command = new[] { "dotnet", "add", csprojPath, "package", "Intel.MKL.redist.win" };
                    if (ExecutePackageManagerInstall(command))
                    {
                        if (InteractiveInstall)
                            Console.WriteLine("✓ MKL installed via NuGet!");
                        return true;
                    }
                }

                if (InteractiveInstall)
                    Console.WriteLine("NuGet install failed or .csproj not found.");
            }

            // Fallback to Chocolatey (system-wide install)
            if (hasChoco)
            {
                if (InteractiveInstall)
                    Console.WriteLine("Attempting Chocolatey installation...");

                var command = new[] { "choco", "install", "-y", "intel-mkl" };
                if (ExecutePackageManagerInstall(command))
                {
                    if (InteractiveInstall)
                        Console.WriteLine("✓ MKL installed via Chocolatey!");
                    return true;
                }
            }

            if (InteractiveInstall)
                Console.WriteLine("Auto-installation failed. Please install manually.");

            return false;
        }

        private string? FindProjectFile()
        {
            // Look for .csproj in current directory and up to 3 levels up
            var currentDir = Directory.GetCurrentDirectory();

            for (var level = 0; level < 4; level++)
            {
                var csprojFiles = Directory.GetFiles(currentDir, "*.csproj");
                if (csprojFiles.Length > 0)
                    return csprojFiles[0];

                var parent = Directory.GetParent(currentDir);
                if (parent == null)
                    break;

                currentDir = parent.FullName;
            }

            return null;
        }

        private bool HandleLinuxMKLInstall()
        {
            var hasApt = CommandExists("apt-get");
            var hasDnf = CommandExists("dnf");
            var hasYum = CommandExists("yum");

            if (!hasApt && !hasDnf && !hasYum)
                return false;

            string[] installCommand;
            if (hasApt)
                installCommand = new[] { "sudo", "apt-get", "install", "-y", "intel-mkl" };
            else if (hasDnf)
                installCommand = new[] { "sudo", "dnf", "install", "-y", "intel-mkl" };
            else
                installCommand = new[] { "sudo", "yum", "install", "-y", "intel-mkl" };

            if (InteractiveInstall)
            {
                Console.WriteLine();
                Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                Console.WriteLine("║  Intel MKL Not Found - Auto-Install Available                 ║");
                Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                Console.WriteLine();
                Console.WriteLine($"This will run: {string.Join(" ", installCommand)}");
                Console.WriteLine();
                Console.Write("Install Intel MKL now? (y/n): ");

                var response = Console.ReadLine()?.Trim().ToLower();
                if (response != "y" && response != "yes")
                    return false;

                Console.WriteLine();
                Console.WriteLine("Installing Intel MKL...");
            }

            return ExecutePackageManagerInstall(installCommand);
        }

        private bool HandleMacOSMKLInstall()
        {
            var hasBrew = CommandExists("brew");

            // If Homebrew not installed, offer to install it
            if (!hasBrew)
            {
                if (InteractiveInstall)
                {
                    Console.WriteLine();
                    Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                    Console.WriteLine("║  Homebrew Not Found                                            ║");
                    Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                    Console.WriteLine();
                    Console.WriteLine("To auto-install MKL, we need Homebrew (macOS package manager).");
                    Console.WriteLine();
                    Console.WriteLine("Would you like to install Homebrew?");
                    Console.WriteLine("  This will enable automatic MKL installation.");
                    Console.WriteLine("  Requires: No sudo (installs to user directory)");
                    Console.WriteLine("  Time: ~2-5 minutes");
                    Console.WriteLine();
                    Console.Write("Install Homebrew now? (y/n): ");

                    var response = Console.ReadLine()?.Trim().ToLower();
                    if (response == "y" || response == "yes")
                    {
                        Console.WriteLine();
                        Console.WriteLine("Installing Homebrew...");
                        Console.WriteLine("(This may take a few minutes)");

                        if (TryInstallHomebrew())
                        {
                            Console.WriteLine("✓ Homebrew installed successfully!");
                            hasBrew = true;
                            // Continue to MKL installation below
                        }
                        else
                        {
                            Console.WriteLine("✗ Homebrew installation failed.");
                            Console.WriteLine();
                            Console.WriteLine("Please install Intel MKL manually:");
                            Console.WriteLine("  Option 1: Install Homebrew, then: brew install intel-mkl");
                            Console.WriteLine("  Option 2: Download from intel.com");
                            Console.WriteLine();
                            return false;
                        }
                    }
                    else
                    {
                        Console.WriteLine();
                        Console.WriteLine("Skipping Homebrew installation.");
                        Console.WriteLine();
                        Console.WriteLine("Please install Intel MKL manually:");
                        Console.WriteLine("  Option 1: Download from intel.com");
                        Console.WriteLine();
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }

            // Now we have Homebrew, proceed with MKL installation
            if (InteractiveInstall)
            {
                Console.WriteLine();
                Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                Console.WriteLine("║  Intel MKL Not Found - Auto-Install Available                 ║");
                Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                Console.WriteLine();
                Console.WriteLine("This will run: brew install intel-mkl");
                Console.WriteLine();
                Console.Write("Install Intel MKL now? (y/n): ");

                var response = Console.ReadLine()?.Trim().ToLower();
                if (response != "y" && response != "yes")
                    return false;

                Console.WriteLine();
                Console.WriteLine("Installing Intel MKL via Homebrew...");
            }

            return ExecutePackageManagerInstall(new[] { "brew", "install", "intel-mkl" });
        }

        private bool TryInstallHomebrew()
        {
            try
            {
                // Homebrew installation script (official method)
                // No sudo required on modern macOS (installs to user directory)
                var installScript =
                    "/bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"";

                var process = new Process
                {
                    StartInfo = new ProcessStartInfo
                    {
                        FileName = "/bin/bash",
                        Arguments = $"-c \"{installScript}\"",
                        RedirectStandardOutput = false, // Let user see progress
                        RedirectStandardError = false,
                        UseShellExecute = false,
                        CreateNoWindow = false
                    }
                };

                process.Start();

                // FIXED: Add timeout to prevent indefinite hang (10 minutes for Homebrew install)
                const int timeoutMs = 10 * 60 * 1000; // 10 minutes
                if (!process.WaitForExit(timeoutMs))
                {
                    Debug.WriteLine("[NativeLibrary] Homebrew installation timed out after 10 minutes");
                    try
                    {
                        process.Kill();
                    }
                    catch
                    {
                    }

                    return false;
                }

                if (process.ExitCode == 0)
                {
                    // Refresh PATH to include Homebrew
                    // On Apple Silicon: /opt/homebrew/bin
                    // On Intel: /usr/local/bin
                    var brewPath = RuntimeInformation.ProcessArchitecture == Architecture.Arm64
                        ? "/opt/homebrew/bin"
                        : "/usr/local/bin";

                    var currentPath = Environment.GetEnvironmentVariable("PATH") ?? "";
                    if (!currentPath.Contains(brewPath))
                        Environment.SetEnvironmentVariable("PATH", brewPath + ":" + currentPath);

                    return true;
                }

                return false;
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[NativeLibrary] Homebrew install error: {ex.Message}");
                return false;
            }
        }

        private bool ExecutePackageManagerInstall(string[] command)
        {
            try
            {
                var process = new Process
                {
                    StartInfo = new ProcessStartInfo
                    {
                        FileName = command[0],
                        Arguments = string.Join(" ", command.Skip(1)),
                        RedirectStandardOutput = !InteractiveInstall,
                        RedirectStandardError = !InteractiveInstall,
                        UseShellExecute = InteractiveInstall,
                        CreateNoWindow = !InteractiveInstall
                    }
                };

                process.Start();

                // FIXED: Add timeout to prevent indefinite hang (15 minutes for package install)
                const int timeoutMs = 15 * 60 * 1000; // 15 minutes
                if (!process.WaitForExit(timeoutMs))
                {
                    Debug.WriteLine(
                        $"[NativeLibrary] Package manager command timed out after 15 minutes: {command[0]}");
                    try
                    {
                        process.Kill();
                    }
                    catch
                    {
                    }

                    return false;
                }

                if (process.ExitCode == 0)
                {
                    Debug.WriteLine("[NativeLibrary] ✓ MKL installed successfully");
                    if (InteractiveInstall)
                    {
                        Console.WriteLine();
                        Console.WriteLine("✓ Intel MKL installed successfully!");
                        Console.WriteLine();
                    }

                    return true;
                }

                return false;
            }
            catch
            {
                return false;
            }
        }

        private bool CommandExists(string command)
        {
            try
            {
                var process = new Process
                {
                    StartInfo = new ProcessStartInfo
                    {
                        FileName = RuntimeInformation.IsOSPlatform(OSPlatform.Windows) ? "where" : "which",
                        Arguments = command,
                        RedirectStandardOutput = true,
                        UseShellExecute = false,
                        CreateNoWindow = true
                    }
                };

                process.Start();

                // FIXED: Add timeout to prevent indefinite hang (10 seconds for which/where)
                if (!process.WaitForExit(10000))
                {
                    try
                    {
                        process.Kill();
                    }
                    catch
                    {
                    }

                    return false;
                }

                return process.ExitCode == 0;
            }
            catch
            {
                return false;
            }
        }

        // Library file patterns - matches ANY version
        private static class LibraryPatterns
        {
            public static readonly Regex CudaRuntimeWindows = new(@"^cudart64(_\d+)?\.dll$", RegexOptions.IgnoreCase);

            public static readonly Regex CudaRuntimeLinux =
                new(@"^libcudart\.so(\.\d+(\.\d+)?)?$", RegexOptions.IgnoreCase);

            public static readonly Regex CuSparseWindows = new(@"^cusparse64(_\d+)?\.dll$", RegexOptions.IgnoreCase);

            public static readonly Regex CuSparseLinux =
                new(@"^libcusparse\.so(\.\d+(\.\d+)?)?$", RegexOptions.IgnoreCase);

            public static readonly Regex MKLWindows = new(@"^mkl_rt(\.\d+)?\.dll$", RegexOptions.IgnoreCase);
            public static readonly Regex MKLLinux = new(@"^libmkl_rt\.so(\.\d+)?$", RegexOptions.IgnoreCase);
            public static readonly Regex MKLMacOS = new(@"^libmkl_rt(\.\d+)?\.dylib$", RegexOptions.IgnoreCase);
        }

        private class LibraryInfo
        {
            public string FullPath { get; set; } = "";
            public string FileName { get; set; } = "";
            public string DirectoryPath { get; set; } = "";
            public double FileVersion { get; set; }
            public double DirectoryVersion { get; set; }
            public string? ResolvedPath { get; set; } // Path after resolving symlinks
            public double ResolvedVersion { get; set; } // Version from resolved file
            public bool IsSymlink { get; set; }
        }
    }

// ============================================================================
// NATIVE LIBRARY STATUS - Availability Checking and Diagnostics
// ============================================================================

    /// <summary>
    ///     Provides information about native library availability and performance.
    ///     Use this to check which backends are available and get diagnostic information.
    /// </summary>
    public static class NativeLibraryStatus
    {
        private static bool _initialized;
        private static readonly object _initLock = new();

        private static readonly LibraryInfo _cudaInfo = new() { LibraryName = "CUDA Runtime" };
        private static readonly LibraryInfo _cuSparseInfo = new() { LibraryName = "cuSPARSE" };
        private static readonly LibraryInfo _mklInfo = new() { LibraryName = "Intel MKL" };
        private static readonly LibraryInfo _pardisoInfo = new() { LibraryName = "PARDISO" };

        /// <summary>
        ///     Checks if CUDA runtime is available for GPU acceleration.
        /// </summary>
        public static bool IsCudaAvailable
        {
            get
            {
                EnsureInitialized();
                return _cudaInfo.IsAvailable;
            }
        }

        /// <summary>
        ///     Checks if cuSPARSE is available for GPU sparse matrix operations.
        /// </summary>
        public static bool IsCuSparseAvailable
        {
            get
            {
                EnsureInitialized();
                return _cuSparseInfo.IsAvailable;
            }
        }

        /// <summary>
        ///     Checks if Intel MKL is available for CPU acceleration (Sparse BLAS).
        /// </summary>
        public static bool IsMklAvailable
        {
            get
            {
                EnsureInitialized();
                return _mklInfo.IsAvailable;
            }
        }

        /// <summary>
        ///     Checks if PARDISO solver is available (requires MKL).
        /// </summary>
        public static bool IsPardisoAvailable
        {
            get
            {
                EnsureInitialized();
                return _pardisoInfo.IsAvailable;
            }
        }

        /// <summary>
        ///     Gets detailed information about CUDA availability.
        /// </summary>
        public static LibraryInfo GetCudaInfo()
        {
            EnsureInitialized();
            return _cudaInfo;
        }

        /// <summary>
        ///     Gets detailed information about cuSPARSE availability.
        /// </summary>
        public static LibraryInfo GetCuSparseInfo()
        {
            EnsureInitialized();
            return _cuSparseInfo;
        }

        /// <summary>
        ///     Gets detailed information about Intel MKL availability.
        /// </summary>
        public static LibraryInfo GetMklInfo()
        {
            EnsureInitialized();
            return _mklInfo;
        }

        /// <summary>
        ///     Gets detailed information about PARDISO solver availability.
        /// </summary>
        public static LibraryInfo GetPardisoInfo()
        {
            EnsureInitialized();
            return _pardisoInfo;
        }

        /// <summary>
        ///     Gets a comprehensive status report of all native libraries.
        /// </summary>
        public static string GetStatusReport()
        {
            EnsureInitialized();

            var report = new StringBuilder();

            report.AppendLine("═══════════════════════════════════════════════════════════");
            report.AppendLine("  Native Library Status Report");
            report.AppendLine("═══════════════════════════════════════════════════════════");
            report.AppendLine();

            // GPU Acceleration
            report.AppendLine("GPU ACCELERATION:");
            AppendLibraryStatus(report, _cudaInfo);
            AppendLibraryStatus(report, _cuSparseInfo);
            report.AppendLine();

            // CPU Acceleration
            report.AppendLine("CPU ACCELERATION:");
            AppendLibraryStatus(report, _mklInfo);
            AppendLibraryStatus(report, _pardisoInfo);
            report.AppendLine();

            // Performance Summary
            report.AppendLine("PERFORMANCE SUMMARY:");
            if (_cudaInfo.IsAvailable && _cuSparseInfo.IsAvailable)
            {
                report.AppendLine("  ✓ GPU acceleration ENABLED (CUDA)");
                report.AppendLine("    Expected performance: Excellent for large sparse matrices");
            }
            else
            {
                report.AppendLine("  × GPU acceleration DISABLED");
                report.AppendLine("    Reason: CUDA not available");
            }

            if (_mklInfo.IsAvailable)
            {
                report.AppendLine("  ✓ CPU acceleration ENABLED (Intel MKL)");
                report.AppendLine("    Expected performance: Very good for medium-sized problems");
            }
            else
            {
                report.AppendLine("  × CPU acceleration DISABLED");
                report.AppendLine("    Reason: Intel MKL not available");
            }

            if (_pardisoInfo.IsAvailable)
            {
                report.AppendLine("  ✓ Direct solver ENABLED (PARDISO)");
                report.AppendLine("    Expected performance: Excellent for direct solutions");
            }
            else
            {
                report.AppendLine("  × Direct solver DISABLED");
                report.AppendLine("    Reason: PARDISO not available");
            }

            if (!_cudaInfo.IsAvailable && !_mklInfo.IsAvailable)
            {
                report.AppendLine();
                report.AppendLine("WARNING: No hardware acceleration available!");
                report.AppendLine("         Performance will be significantly reduced.");
                report.AppendLine("         Consider installing Intel MKL or CUDA.");
            }

            report.AppendLine("═══════════════════════════════════════════════════════════");

            return report.ToString();
        }

        /// <summary>
        ///     Prints a status report to the console.
        /// </summary>
        public static void PrintStatusReport()
        {
            Console.WriteLine(GetStatusReport());
        }

        /// <summary>
        ///     Gets a short summary suitable for logging.
        /// </summary>
        public static string GetShortSummary()
        {
            EnsureInitialized();

            var parts = new List<string>();

            if (_cudaInfo.IsAvailable)
                parts.Add($"CUDA {_cudaInfo.Version}");

            if (_mklInfo.IsAvailable)
                parts.Add($"MKL {_mklInfo.Version}");

            if (parts.Count == 0)
                return "No acceleration available";

            return string.Join(", ", parts);
        }

        private static void AppendLibraryStatus(StringBuilder sb, LibraryInfo info)
        {
            var status = info.IsAvailable ? "✓" : "×";
            var available = info.IsAvailable ? "AVAILABLE" : "NOT AVAILABLE";

            sb.AppendLine($"  {status} {info.LibraryName}: {available}");

            if (info.IsAvailable)
            {
                if (!string.IsNullOrEmpty(info.Version) && info.Version != "Unknown")
                    sb.AppendLine($"      Version: {info.Version}");

                if (!string.IsNullOrEmpty(info.Source))
                    sb.AppendLine($"      Source: {info.Source}");

                if (!string.IsNullOrEmpty(info.Location))
                    sb.AppendLine($"      Location: {info.Location}");
            }
            else
            {
                if (info.MissingDependencies.Length > 0)
                    sb.AppendLine($"      Missing: {string.Join(", ", info.MissingDependencies)}");
            }
        }

        private static void EnsureInitialized()
        {
            if (_initialized) return;

            lock (_initLock)
            {
                if (_initialized) return;

                // Initialize library info by checking backends
                CheckCuda();
                CheckMkl();
                CheckPardiso();

                _initialized = true;
            }
        }

        private static void CheckCuda()
        {
            try
            {
                // Check if CuSparseBackend is available (from CSR.cs)
                var type = Type.GetType("Numerical.CSR+CuSparseBackend, Matrices");
                if (type != null)
                {
                    var isAvailableProp = type.GetProperty("IsAvailable", BindingFlags.Public | BindingFlags.Static);
                    if (isAvailableProp != null)
                    {
                        var isAvailable = (bool)(isAvailableProp.GetValue(null) ?? false);

                        _cudaInfo.IsAvailable = isAvailable;
                        _cuSparseInfo.IsAvailable = isAvailable;

                        if (isAvailable)
                        {
                            // Try to detect version and location
                            var config = new NativeLibraryConfig();
                            var cudaLibs = config.GetCudaRuntimeLibraries();

                            if (cudaLibs != null && cudaLibs.Length > 0)
                            {
                                var version = ExtractVersionFromLibraryName(cudaLibs[0]);
                                _cudaInfo.Version = version;
                                _cuSparseInfo.Version = version;

                                // Try to find actual location
                                var searchPaths = config.GetSearchPaths();
                                if (searchPaths != null)
                                    foreach (var path in searchPaths)
                                    {
                                        var fullPath = Path.Combine(path, cudaLibs[0]);
                                        if (File.Exists(fullPath))
                                        {
                                            _cudaInfo.Location = fullPath;
                                            _cuSparseInfo.Location = Path.GetDirectoryName(fullPath) ?? "";

                                            // Determine source
                                            if (fullPath.Contains("runtimes") && fullPath.Contains("native"))
                                                _cudaInfo.Source = "NuGet package";
                                            else if (fullPath.Contains("NVIDIA"))
                                                _cudaInfo.Source = "NVIDIA Toolkit";
                                            else
                                                _cudaInfo.Source = "System installation";

                                            _cuSparseInfo.Source = _cudaInfo.Source;
                                            break;
                                        }
                                    }
                            }
                        }
                        else
                        {
                            _cudaInfo.MissingDependencies = new[] { "CUDA Toolkit not installed" };
                            _cuSparseInfo.MissingDependencies = new[] { "cuSPARSE (part of CUDA)" };
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[LibraryStatus] Error checking CUDA: {ex.Message}");
            }
        }

        private static void CheckMkl()
        {
            try
            {
                // MKL Sparse BLAS has been removed; MKL is only used for PARDISO solver.
                // Check MKL availability via PARDISO solver
                var config = new NativeLibraryConfig();
                var mklLibs = config.GetMKLLibraries();

                if (mklLibs != null && mklLibs.Length > 0)
                {
                    var searchPaths = config.GetSearchPaths();
                    if (searchPaths != null)
                        foreach (var libName in mklLibs)
                        foreach (var path in searchPaths)
                        {
                            var fullPath = Path.Combine(path, libName);
                            if (File.Exists(fullPath))
                            {
                                _mklInfo.IsAvailable = true;
                                _mklInfo.Version = ExtractVersionFromLibraryName(libName);
                                _mklInfo.Location = fullPath;

                                // Determine source
                                if (fullPath.Contains("runtimes") && fullPath.Contains("native"))
                                    _mklInfo.Source = "NuGet package (Intel.oneAPI.MKL.redist) - PARDISO only";
                                else if (fullPath.Contains("oneapi") || fullPath.Contains("Intel"))
                                    _mklInfo.Source = "Intel oneAPI installation - PARDISO only";
                                else
                                    _mklInfo.Source = "System installation - PARDISO only";

                                return;
                            }
                        }
                }

                _mklInfo.IsAvailable = false;
                _mklInfo.MissingDependencies = new[] { "Intel MKL not installed (required for PARDISO solver)" };
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[LibraryStatus] Error checking MKL: {ex.Message}");
            }
        }

        private static void CheckPardiso()
        {
            try
            {
                // Check if PardisoSolver is available (from CSR.cs)
                var type = Type.GetType("Numerical.PardisoSolver, Matrices");
                if (type != null)
                {
                    // Try to get IsAvailable field/property
                    var isAvailableField = type.GetField("_isAvailable", BindingFlags.NonPublic | BindingFlags.Static);
                    if (isAvailableField != null)
                    {
                        var isAvailable = (bool)(isAvailableField.GetValue(null) ?? false);

                        _pardisoInfo.IsAvailable = isAvailable;

                        if (isAvailable)
                        {
                            _pardisoInfo.Version = _mklInfo.Version; // PARDISO is part of MKL
                            _pardisoInfo.Location = _mklInfo.Location;
                            _pardisoInfo.Source = _mklInfo.Source + " (PARDISO included)";
                        }
                        else
                        {
                            _pardisoInfo.MissingDependencies = new[] { "Intel MKL (includes PARDISO)" };
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[LibraryStatus] Error checking PARDISO: {ex.Message}");
            }
        }

        private static string ExtractVersionFromLibraryName(string libraryName)
        {
            // Extract version from library filename
            // cudart64_15.dll → "15.x"
            // libcudart.so.12.0 → "12.0"
            // mkl_rt.2.dll → "2024.x"

            var match = Regex.Match(libraryName, @"[\._](\d+(?:\.\d+)?)");
            if (match.Success)
            {
                var version = match.Groups[1].Value;

                // For MKL, version 2 typically means 2024.x
                if (libraryName.Contains("mkl") && version == "2")
                    return "2024.x";

                return version;
            }

            return "Unknown";
        }

        /// <summary>
        ///     Information about a specific native library backend.
        /// </summary>
        public class LibraryInfo
        {
            public bool IsAvailable { get; set; }
            public string LibraryName { get; set; } = "";
            public string Version { get; set; } = "Unknown";
            public string Source { get; set; } = "Not found";
            public string Location { get; set; } = "";
            public string[] MissingDependencies { get; set; } = Array.Empty<string>();
        }
    }

// ============================================================================
// NUGET LIBRARY CHECKER - Automatic Package Restoration
// ============================================================================

    /// <summary>
    ///     Checks for NuGet-provided native libraries and provides guidance for automatic restoration.
    ///     This ensures MKL libraries from Intel.oneAPI.MKL.redist packages are available.
    /// </summary>
    public static class NuGetLibraryChecker
    {
        /// <summary>
        ///     Checks if NuGet-provided MKL libraries are present.
        ///     If not, provides clear guidance on how to restore them automatically.
        /// </summary>
        public static bool EnsureMklNuGetPackage()
        {
            if (IsMklNuGetPackagePresent())
            {
                Debug.WriteLine("[NuGet] Intel MKL NuGet package detected");
                return true;
            }

            Debug.WriteLine("[NuGet] Intel MKL NuGet package not found");

            // Check if we're in development environment
            if (IsInDevelopmentEnvironment())
            {
                ShowNuGetRestoreGuidance();
                return false;
            }

            // In production, packages should already be restored
            Debug.WriteLine("[NuGet] Running in production mode - NuGet packages should be pre-installed");
            return false;
        }

        /// <summary>
        ///     Checks if MKL NuGet package libraries are present in the expected location.
        /// </summary>
        public static bool IsMklNuGetPackagePresent()
        {
            try
            {
                var baseDir = AppDomain.CurrentDomain.BaseDirectory;
                var rid = GetRuntimeIdentifier();

                // Expected NuGet location: runtimes/{rid}/native/
                var nugetPath = Path.Combine(baseDir, "runtimes", rid, "native");

                if (!Directory.Exists(nugetPath))
                {
                    Debug.WriteLine($"[NuGet] NuGet native path not found: {nugetPath}");
                    return false;
                }

                // Check for MKL library files
                var mklFiles = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                    ? new[] { "mkl_rt.2.dll", "mkl_rt.1.dll", "mkl_rt.dll" }
                    : new[] { "libmkl_rt.so.2", "libmkl_rt.so.1", "libmkl_rt.so" };

                foreach (var file in mklFiles)
                {
                    var fullPath = Path.Combine(nugetPath, file);
                    if (File.Exists(fullPath))
                    {
                        Debug.WriteLine($"[NuGet] Found MKL library: {fullPath}");
                        return true;
                    }
                }

                Debug.WriteLine("[NuGet] No MKL libraries found in NuGet location");
                return false;
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[NuGet] Error checking for NuGet packages: {ex.Message}");
                return false;
            }
        }

        /// <summary>
        ///     Gets the current runtime identifier.
        /// </summary>
        private static string GetRuntimeIdentifier()
        {
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                return RuntimeInformation.ProcessArchitecture == Architecture.X64 ? "win-x64" : "win-arm64";
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                return RuntimeInformation.ProcessArchitecture == Architecture.X64 ? "linux-x64" : "linux-arm64";
            if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                return RuntimeInformation.ProcessArchitecture == Architecture.X64 ? "osx-x64" : "osx-arm64";
            return "unknown";
        }

        /// <summary>
        ///     Checks if running in a development environment (source code available).
        /// </summary>
        public static bool IsInDevelopmentEnvironment()
        {
            try
            {
                var baseDir = AppDomain.CurrentDomain.BaseDirectory;

                // Check for .csproj file in current or parent directories
                var currentDir = new DirectoryInfo(baseDir);
                while (currentDir != null)
                {
                    if (currentDir.GetFiles("*.csproj").Any())
                        return true;

                    currentDir = currentDir.Parent;
                }

                // Check for bin/Debug or bin/Release pattern (typical dev build)
                if (baseDir.Contains("bin\\Debug") || baseDir.Contains("bin/Debug") ||
                    baseDir.Contains("bin\\Release") || baseDir.Contains("bin/Release"))
                    return true;

                return false;
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        ///     Shows guidance for restoring NuGet packages automatically.
        /// </summary>
        private static void ShowNuGetRestoreGuidance()
        {
            Console.WriteLine();
            Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║  Intel MKL NuGet Package Not Found                             ║");
            Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();
            Console.WriteLine("The Intel MKL libraries are provided via NuGet packages but are not");
            Console.WriteLine("currently available. This is typically caused by missing NuGet restore.");
            Console.WriteLine();
            Console.WriteLine("AUTOMATIC FIX:");
            Console.WriteLine("  Run: dotnet restore");
            Console.WriteLine();
            Console.WriteLine("This will automatically download and install:");
            Console.WriteLine("  - Intel.oneAPI.MKL.redist.win (Windows)");
            Console.WriteLine("  - Intel.oneAPI.MKL.redist.lin (Linux)");
            Console.WriteLine();
            Console.WriteLine("Expected download size: ~250 MB");
            Console.WriteLine("Location: runtimes/{platform}/native/");
            Console.WriteLine();
            Console.WriteLine("After running 'dotnet restore', rebuild and run your application.");
            Console.WriteLine();
            Console.WriteLine("NOTE: This is a ONE-TIME setup. Packages will be cached locally.");
            Console.WriteLine();
        }

        /// <summary>
        ///     Attempts to automatically restore NuGet packages by invoking dotnet CLI.
        ///     Returns true if restore succeeded.
        /// </summary>
        public static bool TryAutoRestoreNuGetPackages()
        {
            try
            {
                if (!IsInDevelopmentEnvironment())
                {
                    Debug.WriteLine("[NuGet] Not in development environment - skipping auto-restore");
                    return false;
                }

                // Find the .csproj file
                var projectFile = FindProjectFile();
                if (projectFile == null)
                {
                    Debug.WriteLine("[NuGet] Could not find .csproj file for auto-restore");
                    return false;
                }

                Console.WriteLine();
                Console.WriteLine("Attempting automatic NuGet package restore...");
                Console.WriteLine($"Project: {Path.GetFileName(projectFile)}");
                Console.WriteLine();

                var process = new Process
                {
                    StartInfo = new ProcessStartInfo
                    {
                        FileName = "dotnet",
                        Arguments = $"restore \"{projectFile}\"",
                        WorkingDirectory = Path.GetDirectoryName(projectFile),
                        RedirectStandardOutput = true,
                        RedirectStandardError = true,
                        UseShellExecute = false,
                        CreateNoWindow = false
                    }
                };

                process.Start();

                // Show output in real-time with timeout protection
                // FIXED: Add timeout to prevent indefinite hang (5 minutes for NuGet restore)
                const int timeoutMs = 5 * 60 * 1000; // 5 minutes
                var outputTask = Task.Run(() =>
                {
                    while (!process.StandardOutput.EndOfStream)
                    {
                        var line = process.StandardOutput.ReadLine();
                        Console.WriteLine(line);
                    }
                });

                if (!process.WaitForExit(timeoutMs))
                {
                    Console.WriteLine();
                    Console.WriteLine("✗ NuGet restore timed out after 5 minutes");
                    try
                    {
                        process.Kill();
                    }
                    catch
                    {
                    }

                    return false;
                }

                // Give the output task a moment to finish
                outputTask.Wait(1000);

                if (process.ExitCode == 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("✓ NuGet packages restored successfully!");
                    Console.WriteLine("  Please restart the application to use the restored libraries.");
                    Console.WriteLine();
                    return true;
                }

                var error = process.StandardError.ReadToEnd();
                Console.WriteLine();
                Console.WriteLine($"✗ NuGet restore failed (exit code {process.ExitCode})");
                if (!string.IsNullOrEmpty(error))
                    Console.WriteLine($"  Error: {error}");
                Console.WriteLine();
                return false;
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[NuGet] Auto-restore error: {ex.Message}");
                Console.WriteLine($"✗ Auto-restore failed: {ex.Message}");
                Console.WriteLine("  Please run 'dotnet restore' manually.");
                return false;
            }
        }

        private static string? FindProjectFile()
        {
            try
            {
                var baseDir = new DirectoryInfo(AppDomain.CurrentDomain.BaseDirectory);

                // Search up to 5 levels up for .csproj file
                var currentDir = baseDir;
                for (var i = 0; i < 5 && currentDir != null; i++)
                {
                    var projectFiles = currentDir.GetFiles("*.csproj");
                    if (projectFiles.Length > 0) return projectFiles[0].FullName;

                    currentDir = currentDir.Parent;
                }

                return null;
            }
            catch
            {
                return null;
            }
        }
    }

    /// <summary>
    ///     Enhanced Ultimate+ configuration with automatic NuGet awareness.
    ///     Checks for NuGet packages first, then falls back to system installation.
    /// </summary>
    public sealed class FinalUltimateNativeLibraryConfig : INativeLibraryConfig
    {
        private static bool _nugetCheckPerformed;
        private readonly NativeLibraryConfig _baseConfig;

        public FinalUltimateNativeLibraryConfig()
        {
            _baseConfig = new NativeLibraryConfig();

            // Auto-detect development vs production
            EnableNuGetAutoRestore = NuGetLibraryChecker.IsInDevelopmentEnvironment();
        }

        /// <summary>
        ///     If true, automatically attempts to restore NuGet packages if MKL not found.
        ///     Default: true for development, false for production
        /// </summary>
        public bool EnableNuGetAutoRestore { get; set; }

        /// <summary>
        ///     If true, enables automatic MKL installation via package manager.
        ///     Default: true
        /// </summary>
        public bool EnableAutoInstall
        {
            get => _baseConfig.EnableAutoInstall;
            set => _baseConfig.EnableAutoInstall = value;
        }

        /// <summary>
        ///     If true, shows interactive prompts to user.
        ///     Default: true
        /// </summary>
        public bool InteractiveInstall
        {
            get => _baseConfig.InteractiveInstall;
            set => _baseConfig.InteractiveInstall = value;
        }

        public string[]? GetCudaRuntimeLibraries()
        {
            return _baseConfig.GetCudaRuntimeLibraries();
        }

        public string[]? GetCuSparseLibraries()
        {
            return _baseConfig.GetCuSparseLibraries();
        }

        public string[]? GetMKLLibraries()
        {
            // First-time check for NuGet packages
            if (!_nugetCheckPerformed)
            {
                _nugetCheckPerformed = true;

                if (!NuGetLibraryChecker.IsMklNuGetPackagePresent())
                {
                    Debug.WriteLine("[Final] MKL NuGet package not found");

                    if (EnableNuGetAutoRestore)
                    {
                        Console.WriteLine();
                        Console.WriteLine("Intel MKL NuGet package not found.");
                        Console.Write("Attempt automatic NuGet restore? (y/n): ");

                        var response = Console.ReadLine()?.Trim().ToLower();
                        if (response == "y" || response == "yes")
                            if (NuGetLibraryChecker.TryAutoRestoreNuGetPackages())
                                // Packages installed successfully - restart required
                                throw new InvalidOperationException(
                                    "Intel MKL NuGet packages have been installed successfully. " +
                                    "The application must be restarted for the changes to take effect.");
                    }
                    else
                    {
                        Debug.WriteLine("[Final] NuGet auto-restore disabled");
                    }
                }
            }

            // Use base config (which includes deep scanning and auto-install)
            return _baseConfig.GetMKLLibraries();
        }

        public string[]? GetSearchPaths()
        {
            return _baseConfig.GetSearchPaths();
        }
    }

// ============================================================================
// ROBUST NATIVE LIBRARY LOADER - Safe Library Loading with Fallback
// ============================================================================

    /// <summary>
    ///     Robust native library loading with multiple path attempts and detailed diagnostics.
    /// </summary>
    internal static class RobustNativeLibraryLoader
    {
        /// <summary>
        ///     Tries to load a native library from multiple names and search paths.
        /// </summary>
        /// <param name="libraryNames">Ordered list of library names to try (e.g., cudart64_15.dll, cudart64_14.dll)</param>
        /// <param name="searchPaths">Ordered list of directories to search</param>
        /// <param name="libraryType">Description for logging (e.g., "CUDA Runtime")</param>
        /// <returns>Library handle if successful, IntPtr.Zero if all attempts fail</returns>
        public static IntPtr TryLoadLibrary(string[]? libraryNames, string[]? searchPaths, string libraryType)
        {
            if (libraryNames == null || libraryNames.Length == 0)
            {
                Debug.WriteLine($"[NativeLibrary] No library names provided for {libraryType}");
                return IntPtr.Zero;
            }

            if (searchPaths == null || searchPaths.Length == 0)
            {
                Debug.WriteLine($"[NativeLibrary] No search paths provided for {libraryType}");
                return IntPtr.Zero;
            }

            Debug.WriteLine($"[NativeLibrary] Attempting to load {libraryType}...");
            Debug.WriteLine($"[NativeLibrary]   Library names: {string.Join(", ", libraryNames.Take(5))}");
            Debug.WriteLine($"[NativeLibrary]   Search paths: {searchPaths.Length} paths");

            // CRITICAL: Try each version in order (latest first) until one WORKS
            // This handles corrupt installations, incompatibilities, missing dependencies
            var attemptedPaths = new List<string>();
            var foundButAllFailed = false; // Track if we found libraries but all failed verification

            // Try each library name
            foreach (var libName in libraryNames)
            {
                // Strategy 1: Try without explicit path (use system search)
                // This works if symlinks exist (libcudart.so -> libcudart.so.11.0)
                try
                {
                    if (NativeLibrary.TryLoad(libName, out var handle))
                    {
                        foundButAllFailed = true; // Found at least one library

                        // SUCCESS! But verify it actually works
                        if (VerifyLibraryWorks(handle, libraryType))
                        {
                            Debug.WriteLine($"[NativeLibrary] ✓ Loaded {libraryType}: {libName} (system search)");
                            return handle;
                        }

                        Debug.WriteLine($"[NativeLibrary] ✗ {libName} loaded but verification failed, trying next...");
                        NativeLibrary.Free(handle);
                    }
                }
                catch (Exception ex)
                {
                    Debug.WriteLine($"[NativeLibrary] ✗ Failed to load {libName}: {ex.Message}");
                }

                // Strategy 2: Try with explicit full path
                // This works even without symlinks (direct load of libcudart.so.11.0)
                foreach (var searchPath in searchPaths)
                    try
                    {
                        var fullPath = Path.Combine(searchPath, libName);
                        if (File.Exists(fullPath) && !attemptedPaths.Contains(fullPath))
                        {
                            attemptedPaths.Add(fullPath);
                            foundButAllFailed = true; // Found at least one library

                            if (NativeLibrary.TryLoad(fullPath, out var handle))
                            {
                                // SUCCESS! But verify it actually works
                                if (VerifyLibraryWorks(handle, libraryType))
                                {
                                    Debug.WriteLine($"[NativeLibrary] ✓ Loaded {libraryType}: {fullPath}");
                                    return handle;
                                }

                                Debug.WriteLine(
                                    $"[NativeLibrary] ✗ {fullPath} loaded but verification failed, trying next...");
                                NativeLibrary.Free(handle);
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        Debug.WriteLine($"[NativeLibrary] ✗ Failed to load {libName} from {searchPath}: {ex.Message}");
                    }

                // Strategy 3: AGGRESSIVE - Look for ANY file starting with the base name
                // e.g., libName = "libcudart.so" -> try "libcudart.so.11.0", "libcudart.so.11.8", etc.
                // This handles missing symlinks for PhD students who just apt-get installed
                var baseName = libName.Split('.')[0]; // "libcudart" from "libcudart.so"
                foreach (var searchPath in searchPaths)
                    try
                    {
                        if (!Directory.Exists(searchPath))
                            continue;

                        // Find all files matching the base name
                        var pattern = baseName + ".so*";
                        var candidates = Directory.GetFiles(searchPath, pattern)
                            .Where(f => !attemptedPaths.Contains(f))
                            .OrderBy(f => f.Count(c => c == '.')) // Try simpler names first
                            .ThenByDescending(f => f); // Then newest versions

                        foreach (var candidate in candidates)
                        {
                            attemptedPaths.Add(candidate);
                            foundButAllFailed = true; // Found at least one library

                            try
                            {
                                if (NativeLibrary.TryLoad(candidate, out var handle))
                                {
                                    // SUCCESS! But verify it actually works
                                    if (VerifyLibraryWorks(handle, libraryType))
                                    {
                                        Debug.WriteLine(
                                            $"[NativeLibrary] ✓ Loaded {libraryType}: {candidate} (aggressive search)");
                                        return handle;
                                    }

                                    Debug.WriteLine(
                                        $"[NativeLibrary] ✗ {candidate} loaded but verification failed, trying next...");
                                    NativeLibrary.Free(handle);
                                }
                            }
                            catch (Exception ex)
                            {
                                Debug.WriteLine($"[NativeLibrary] ✗ Failed to load {candidate}: {ex.Message}");
                            }
                        }
                    }
                    catch
                    {
                    }
            }

            // CRITICAL: If we found libraries but ALL failed verification, try auto-install
            if (foundButAllFailed)
            {
                Debug.WriteLine($"[NativeLibrary] ✗ Failed to load {libraryType} from any location");
                Debug.WriteLine("[NativeLibrary]   Note: Libraries were found but all failed verification");
                Debug.WriteLine(
                    "[NativeLibrary]   Possible causes: corrupt installations, incompatible versions, missing dependencies");
                if (libraryType.Contains("MKL", StringComparison.OrdinalIgnoreCase))
                    Debug.WriteLine(
                        "[NativeLibrary]   Recommendation: Try 'sudo apt-get install --reinstall intel-mkl' or similar");
            }
            else
            {
                Debug.WriteLine($"[NativeLibrary] ✗ Failed to load {libraryType} from any location");
            }

            Debug.WriteLine($"[NativeLibrary]   Attempted {attemptedPaths.Count} different files");
            return IntPtr.Zero;
        }

        /// <summary>
        ///     Verifies that a loaded library actually works by calling a basic function.
        ///     This catches corrupt installations, incompatibilities, missing dependencies.
        /// </summary>
        private static bool VerifyLibraryWorks(IntPtr handle, string libraryType)
        {
            if (handle == IntPtr.Zero)
                return false;

            try
            {
                // For CUDA libraries, try to ACTUALLY CALL a function, not just check symbols
                if (libraryType.Contains("CUDA", StringComparison.OrdinalIgnoreCase))
                {
                    // Strategy 1: Try cudaGetDeviceCount (most reliable test)
                    if (NativeLibrary.TryGetExport(handle, "cudaGetDeviceCount", out var funcPtr))
                        try
                        {
                            // Actually call the function to verify it works!
                            // cudaGetDeviceCount has signature: cudaError_t cudaGetDeviceCount(int* count)
                            var cudaGetDeviceCount =
                                Marshal.GetDelegateForFunctionPointer<CudaGetDeviceCountDelegate>(funcPtr);
                            var deviceCount = 0;
                            var result = cudaGetDeviceCount(ref deviceCount);

                            // cudaSuccess = 0
                            // cudaErrorNoDevice = 100 (no GPU, but driver/library OK)
                            // cudaErrorInsufficientDriver = 35 (driver too old)
                            // cudaErrorInitializationError = 3 (major problem)
                            if (result == 0)
                            {
                                // Success! Library works and we have GPUs
                                Debug.WriteLine($"[NativeLibrary] ✓ CUDA library works! Found {deviceCount} GPU(s)");
                                return true;
                            }

                            if (result == 100)
                            {
                                // No GPU present, but library and driver are OK
                                Debug.WriteLine(
                                    $"[NativeLibrary] ✓ CUDA library works (no GPUs present, error code {result})");
                                return true;
                            }

                            if (result == 35)
                            {
                                // Driver too old for this CUDA version
                                Debug.WriteLine($"[NativeLibrary] ✗ CUDA driver insufficient (error code {result})");
                                return false;
                            }

                            // Some other error - library might be corrupt
                            Debug.WriteLine($"[NativeLibrary] ✗ CUDA library call failed (error code {result})");
                            return false;
                        }
                        catch (Exception ex)
                        {
                            Debug.WriteLine($"[NativeLibrary] ✗ CUDA function call failed: {ex.Message}");
                            return false;
                        }

                    // Strategy 2: Try cudaRuntimeGetVersion (less reliable but works without GPU)
                    if (NativeLibrary.TryGetExport(handle, "cudaRuntimeGetVersion", out funcPtr))
                        try
                        {
                            // Actually call the function
                            var cudaRuntimeGetVersion =
                                Marshal.GetDelegateForFunctionPointer<CudaGetVersionDelegate>(funcPtr);
                            var version = 0;
                            var result = cudaRuntimeGetVersion(ref version);

                            if (result == 0)
                            {
                                Debug.WriteLine($"[NativeLibrary] ✓ CUDA library works! Runtime version: {version}");
                                return true;
                            }

                            Debug.WriteLine($"[NativeLibrary] ✗ cudaRuntimeGetVersion failed (error code {result})");
                            return false;
                        }
                        catch (Exception ex)
                        {
                            Debug.WriteLine($"[NativeLibrary] ✗ cudaRuntimeGetVersion call failed: {ex.Message}");
                            return false;
                        }

                    // If we can't find or call basic CUDA functions, library is broken
                    Debug.WriteLine("[NativeLibrary] ✗ Library loaded but CUDA symbols not found or not callable");
                    return false;
                }

                // For MKL libraries, try to get a known MKL function
                if (libraryType.Contains("MKL", StringComparison.OrdinalIgnoreCase))
                {
                    // Try mkl_get_version or similar
                    if (NativeLibrary.TryGetExport(handle, "mkl_get_version", out var funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "MKL_Get_Max_Threads", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "mkl_get_max_threads", out funcPtr)) // lowercase variant
                        return true;

                    Debug.WriteLine("[NativeLibrary] ✗ Library loaded but MKL symbols not found");
                    return false;
                }

                // For PARDISO libraries, try to get the main pardiso function
                if (libraryType.Contains("PARDISO", StringComparison.OrdinalIgnoreCase))
                {
                    // PARDISO main function (different naming conventions across platforms)
                    // Linux/macOS: pardiso, pardiso_ (Fortran)
                    // Windows: PARDISO, pardiso (may vary by compiler)
                    if (NativeLibrary.TryGetExport(handle, "pardiso", out var funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "pardiso_", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "PARDISO", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "PARDISO_", out funcPtr))
                        return true;

                    // Try MKL PARDISO (part of MKL)
                    if (NativeLibrary.TryGetExport(handle, "pardisoinit", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "pardisoinit_", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "PARDISOINIT", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "pardiso_64", out funcPtr))
                        return true;

                    Debug.WriteLine("[NativeLibrary] ✗ Library loaded but PARDISO symbols not found");
                    return false;
                }

                // For LAPACK libraries, try to get a known LAPACK function
                if (libraryType.Contains("LAPACK", StringComparison.OrdinalIgnoreCase) ||
                    libraryType.Contains("BLAS", StringComparison.OrdinalIgnoreCase))
                {
                    // Try common LAPACK functions (dgesv = solve linear system)
                    // Naming conventions vary by platform and compiler:
                    // - Linux/macOS GCC/gfortran: dgesv_
                    // - Windows Intel Fortran: DGESV or dgesv
                    // - Windows MSVC: DGESV
                    if (NativeLibrary.TryGetExport(handle, "dgesv", out var funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "dgesv_", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "DGESV", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "DGESV_", out funcPtr))
                        return true;

                    // Try BLAS function (dgemm = matrix multiply)
                    if (NativeLibrary.TryGetExport(handle, "dgemm", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "dgemm_", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "DGEMM", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "DGEMM_", out funcPtr))
                        return true;

                    // Try cblas_ prefix (common on some systems)
                    if (NativeLibrary.TryGetExport(handle, "cblas_dgemm", out funcPtr) ||
                        NativeLibrary.TryGetExport(handle, "cblas_dgesv", out funcPtr))
                        return true;

                    Debug.WriteLine("[NativeLibrary] ✗ Library loaded but LAPACK/BLAS symbols not found");
                    return false;
                }

                // For other libraries, if we can load them, assume they work
                // (we don't know what functions to check)
                return true;
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"[NativeLibrary] ✗ Verification failed: {ex.Message}");
                return false;
            }
        }
    }

    #endregion
}

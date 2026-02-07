using System.Buffers;
using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Text;
using System.Text.Json.Serialization;

namespace Numerical;

/// <summary>
///     A high-performance, sealed class for dense matrix numerical and scientific computing.
/// </summary>
/// <remarks>
///     <para>
///         The <see cref="Matrix" /> class provides a comprehensive suite of matrix operations,
///         linear algebra decompositions (LU, QR, SVD, Eigenvalues), and statistical functions
///         optimized for dense matrices typically encountered in scientific computing, finite element analysis,
///         machine learning, and optimization problems.
///     </para>
///     <h3>Design Philosophy</h3>
///     <para>
///         This implementation prioritizes <strong>correctness</strong>, <strong>performance</strong>,
///         and <strong>numerical stability</strong>. It is optimized for matrices up to approximately
///         100×100 elements, though it handles larger matrices efficiently as well.
///     </para>
///     <h3>Key Architectural Features</h3>
///     <h4>1. Column-Major Storage</h4>
///     <para>
///         Unlike standard C# 2D arrays (<c>double[,]</c>), data is stored in a 1D array (<c>double[]</c>)
///         in column-major order: <c>_data[column * RowCount + row]</c>. This layout provides:
///     </para>
///     <list type="bullet">
///         <item>
///             <description>Optimal cache locality for column-wise operations</description>
///         </item>
///         <item>
///             <description>Compatibility with BLAS/LAPACK conventions</description>
///         </item>
///         <item>
///             <description>Efficient SIMD vectorization of contiguous columns</description>
///         </item>
///         <item>
///             <description>Seamless interoperability with Fortran libraries (e.g., Intel MKL, PARDISO)</description>
///         </item>
///     </list>
///     <h4>2. SIMD Acceleration</h4>
///     <para>
///         Core arithmetic operations (<c>+</c>, <c>-</c>, <c>*</c>), vector operations (<c>Dot</c>, norms),
///         and matrix multiplication are vectorized using:
///     </para>
///     <list type="bullet">
///         <item>
///             <description><c>Vector512</c> (AVX-512) when hardware support is detected</description>
///         </item>
///         <item>
///             <description><c>Vector256</c> (AVX2) for broader hardware compatibility</description>
///         </item>
///         <item>
///             <description>FMA (Fused Multiply-Add) instructions for improved accuracy and performance</description>
///         </item>
///     </list>
///     <h4>3. Intelligent Parallelization</h4>
///     <para>
///         Operations on large matrices automatically leverage multi-core processors using <c>Parallel.For</c>.
///         Parallelization is carefully tuned to avoid thread overhead for matrices below ~100×100.
///     </para>
///     <h4>4. Optimized Computational Kernels</h4>
///     <para>
///         <strong>Matrix Multiplication (GEMM):</strong>
///     </para>
///     <list type="bullet">
///         <item>
///             <description>Cache-aware blocking for L1/L2 cache optimization</description>
///         </item>
///         <item>
///             <description>SIMD micro-kernels processing 4 columns at a time</description>
///         </item>
///         <item>
///             <description>Specialized fast paths for 1×1, 2×2, and 3×3 products</description>
///         </item>
///         <item>
///             <description>Support for transposed operations (A×B, A'×B, A×B', A'×B')</description>
///         </item>
///     </list>
///     <para>
///         <strong>Determinant &amp; Inverse:</strong>
///     </para>
///     <list type="bullet">
///         <item>
///             <description>Closed-form analytic solutions for 1×1 through 3×3 matrices (O(1) complexity)</description>
///         </item>
///         <item>
///             <description>LU decomposition with Rook pivoting for larger matrices</description>
///         </item>
///     </list>
///     <h4>5. Efficient Memory Management</h4>
///     <list type="bullet">
///         <item>
///             <description>
///                 <strong>Array Pooling:</strong> Uses <c>ArrayPool&lt;double&gt;</c> for large temporary
///                 buffers in SVD and other algorithms
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>Stack Allocation:</strong> Uses <c>stackalloc</c> for small temporary buffers to avoid
///                 heap allocations
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>Pinned Arrays:</strong> Large matrices use pinned arrays to reduce GC pressure during
///                 unsafe operations
///         </item>
///         </item>
///     </list>
///     <h3>Performance Characteristics</h3>
///     <para>For 50×50 dense matrices on modern x64 CPUs with AVX2/FMA support:</para>
///     <list type="bullet">
///         <item>
///             <description>Matrix multiplication: ~0.15 ms</description>
///         </item>
///         <item>
///             <description>LU decomposition: ~0.28 ms</description>
///         </item>
///         <item>
///             <description>QR decomposition: ~0.65 ms</description>
///         </item>
///         <item>
///             <description>SVD: ~3.5 ms</description>
///         </item>
///         <item>
///             <description>Matrix inverse (via LU): ~0.42 ms</description>
///         </item>
///     </list>
///     <h3>Iteration Patterns</h3>
///     <para>
///         Due to column-major storage, performance is maximized by iterating with a
///         <strong>column-first (outer) loop</strong>:
///     </para>
///     <code>
/// // ✓ OPTIMAL (Cache-friendly)
/// for (int j = 0; j &lt; matrix.ColumnCount; j++)
///     for (int i = 0; i &lt; matrix.RowCount; i++)
///         var x = matrix[i, j]; // Accesses contiguous data
///
/// // ✗ SUBOPTIMAL (Cache-unfriendly)
/// for (int i = 0; i &lt; matrix.RowCount; i++)
///     for (int j = 0; j &lt; matrix.ColumnCount; j++)
///         var x = matrix[i, j]; // Stride access pattern
/// </code>
///     <h3>Numerical Stability</h3>
///     <para>All decompositions and solvers are implemented with numerical stability in mind:</para>
///     <list type="bullet">
///         <item>
///             <description><strong>LU:</strong> Rook pivoting for optimal stability</description>
///         </item>
///         <item>
///             <description><strong>QR:</strong> Householder reflections (backward stable)</description>
///         </item>
///         <item>
///             <description><strong>SVD:</strong> Jacobi iterations with careful convergence criteria</description>
///         </item>
///         <item>
///             <description>
///                 <strong>Eigenvalues:</strong> QR algorithm for symmetric matrices (general matrices not yet
///                 supported)
///         </item>
///         </item>
///     </list>
///     <h3>Limitations and Scope</h3>
///     <list type="bullet">
///         <item>
///             <description>
///                 <strong>Dense Matrices Only:</strong> This library is designed for dense matrices. For sparse
///                 matrices (&gt;90% zeros), use specialized sparse matrix libraries.
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>Real Numbers Only:</strong> This library operates exclusively on real-valued (double precision)
///                 matrices. Complex number support is not currently implemented. For applications requiring complex
///                 arithmetic
///                 (e.g., FFT, quantum mechanics), consider libraries with native complex support.
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>Eigenvalues:</strong> Currently supports only symmetric matrices. General eigenvalue
///                 decomposition will be added in a future version.
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>Thread Safety:</strong> Matrix instances are NOT thread-safe for concurrent
///                 modifications. Read-only operations (e.g., indexing, arithmetic) are safe across multiple threads if no
///                 thread modifies the matrix.
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>Floating-Point Limitations:</strong> <c>Equals()</c> uses tolerance-based comparison,
///                 but <c>GetHashCode()</c> does not. Matrices should not be used as dictionary keys if tolerance-based
///                 equality is expected.
///         </item>
///         </item>
///     </list>
///     <h3>Usage Examples</h3>
///     <code>
/// // Create matrices (C# 13 collection expressions)
/// double[][] data = [
///     [1, 2, 3],
///     [4, 5, 6]
/// ];
/// var A = new Matrix(data);
/// var B = Matrix.Identity(3);
///
/// // Basic operations
/// var C = A * B;
/// var D = A.Transpose();
///
/// // Linear algebra
/// var lu = A.ComputeLU();
/// var x = A.Solve(b);  // Solve Ax = b
///
/// // Decompositions
/// var qr = A.ComputeQR();
/// var svd = A.ComputeSVD();
///
/// // Statistics
/// var cov = A.Covariance();
/// var means = A.ColumnMeans();
/// </code>
///     <h3>Version History</h3>
///     <list type="bullet">
///         <item>
///             <description><strong>v1.0:</strong> Initial release</description>
///         </item>
///         <item>
///             <description><strong>v1.1:</strong> Added SIMD optimizations, Rook pivoting, improved numerical stability</description>
///         </item>
///         <item>
///             <description>
///                 <strong>v1.2:</strong> Fixed critical bugs in 2×2 and 4×4 inversions, improved parallel
///                 thresholds
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>v1.3:</strong> Fixed logic errors in 2x2 and 3x3 analytic inverses.
///                 Improved performance of InfinityNorm with cache-friendly access.
///         </item>
///         </item>
///         <item>
///             <description>
///                 <strong>v1.4:</strong> Removed unstable Inverse4x4. Improved QR stability check.
///                 Refactored SVD inner loop with reusable SIMD kernels.
///         </item>
///         </item>
///     </list>
/// </remarks>
/// <summary>
///     Helper class for throwing common exceptions efficiently.
/// </summary>
/// <remarks>
///     Methods are marked with <c>NoInlining</c> to keep hot paths small and improve code cache utilization.
/// </remarks>

internal static class ThrowHelper
{
    [MethodImpl(MethodImplOptions.NoInlining)]
    internal static void ThrowArgumentOutOfRangeException(string paramName)
    {
        throw new ArgumentOutOfRangeException(paramName);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    internal static void ThrowArgumentException(string message, string paramName)
    {
        throw new ArgumentException(message, paramName);
    }
}

public sealed class Matrix : IEquatable<Matrix>, IFormattable, ICloneable
{
    #region SIMD Alignment Notes

    /// <summary>
    ///     SIMD alignment boundary in bytes. AVX-512 benefits from 64-byte alignment, AVX2 from 32-byte.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <strong>Current Implementation:</strong> Uses pinned GC arrays which may not be aligned
    ///         to SIMD boundaries. Modern CPUs handle unaligned loads via LoadUnsafe/LoadVector methods
    ///         with ~10-20% performance penalty compared to perfectly aligned loads.
    ///     </para>
    ///     <para>
    ///         <strong>Performance Impact:</strong> For most workloads, the unaligned load penalty is
    ///         acceptable and outweighed by GC management benefits and API simplicity. Mission-critical
    ///         applications requiring maximum SIMD performance should consider:
    ///     </para>
    ///     <list type="bullet">
    ///         <item>Custom allocator with NativeMemory.AlignedAlloc (requires extensive refactoring)</item>
    ///         <item>Specialized libraries like Intel MKL which handle alignment internally</item>
    ///         <item>Accepting the ~10-20% penalty for simpler, safer managed code</item>
    ///     </list>
    ///     <para>
    ///         <strong>Code Pattern:</strong> Current code uses LoadVector512/LoadVector256 methods
    ///         which handle unaligned data correctly (unlike LoadAlignedVector512 which requires alignment).
    ///     </para>
    /// </remarks>
    private const int SimdAlignmentBytes = 64;

    #endregion

    #region Memory Pool (Internal)

    /// <summary>
    ///     Internal pool for temporary buffers used in complex algorithms (e.g., SVD).
    ///     This is an implementation detail and is not exposed to the user.
    /// </summary>
    private static readonly ArrayPool<double> MatrixPool =
        ArrayPool<double>.Create(MaxRecommendedSize * MaxRecommendedSize, 50);

    #endregion

    #region Constants

    private const double DefaultTolerance = 1e-10;

    private const int MaxRecommendedSize = 20;

    // MEDIUM FIX M10: Increased from 100 to 1000 for better convergence
    // Some difficult matrices require more iterations to reach tolerance
    private const int MaxIterations = 1000;
    private const int StackAllocThreshold = 400;
    private const int BlockSize = 8;

    /// <summary>
    ///     Threshold for when to use parallel processing. Operations on matrices
    ///     with total elements less than this use serial processing to avoid overhead.
    ///     Default: 10,000 elements (e.g., 100×100 matrix).
    ///     Tuned for typical x64 CPUs; adjust based on your hardware.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <strong>THREAD SAFETY WARNING:</strong> This is a GLOBAL setting affecting ALL threads
    ///         in the current process. Changing this value in one thread affects all other threads immediately.
    ///     </para>
    ///     <para>
    ///         <strong>MULTI-TENANT APPLICATIONS:</strong> Not safe for multi-tenant scenarios (e.g., ASP.NET)
    ///         where different requests might need different parallelization thresholds. Consider using
    ///         thread-local configuration or separate application domains if needed.
    ///     </para>
    /// </remarks>
    public static int ParallelThreshold { get; set; } = 10_000;

    private const int L1BlockSize = 48;
    private const int L2BlockSize = 192;
    private const int MicroKernelN = 4;

    /// <summary>
    ///     Adaptive minimum dimension for parallel execution based on CPU core count.
    ///     High core count systems benefit from earlier parallelization.
    /// </summary>
    private static readonly int MinParallelDimension = ParallelConfig.ProcessorCount >= 8 ? 64 :
        ParallelConfig.ProcessorCount >= 4 ? 80 : 100;

    private const double MachinePrecision = 2.220446049250313e-16;
    private const double ZeroVectorThreshold = 1e-10;
    private const double SingularValueTolerance = 1e-10;
    private const double DeterminantUnderflowThreshold = 1e-300;

    private const int MaxRookPivotingIterations = 10;
    private const int MaxMatrixSquareRootIterations = 50;
    private const double MatrixSquareRootTolerance = 1e-10;

    /// <summary>
    ///     Maximum allowed dimension for a matrix (rows or columns).
    ///     Can be adjusted for advanced use cases requiring larger matrices.
    ///     Default: 100,000
    ///     Maximum: 46,340 (sqrt of int.MaxValue) to prevent overflow in size calculations.
    /// </summary>
    private static int _maxMatrixDimension = 100_000;

    public static int MaxMatrixDimension
    {
        get => _maxMatrixDimension;
        set
        {
            // Enforce maximum to prevent integer overflow in rows * columns
            // sqrt(int.MaxValue) ≈ 46340
            const int absoluteMax = 46340;
            if (value > absoluteMax)
                throw new ArgumentOutOfRangeException(nameof(value),
                    $"MaxMatrixDimension cannot exceed {absoluteMax} to prevent integer overflow in size calculations. " +
                    $"For larger matrices, consider using sparse representations or chunked storage.");
            if (value <= 0)
                throw new ArgumentOutOfRangeException(nameof(value), "MaxMatrixDimension must be positive.");
            _maxMatrixDimension = value;
        }
    }

    /// <summary>
    ///     Maximum allowed total elements in a matrix (rows × columns).
    ///     Can be adjusted for advanced use cases requiring larger matrices.
    ///     Default: 100,000,000
    ///     Maximum: int.MaxValue (2,147,483,647) due to array indexing constraints.
    /// </summary>
    private static long _maxTotalElements = 100_000_000;

    public static long MaxTotalElements
    {
        get => _maxTotalElements;
        set
        {
            // Enforce maximum to prevent integer overflow in array allocation
            if (value > int.MaxValue)
                throw new ArgumentOutOfRangeException(nameof(value),
                    $"MaxTotalElements cannot exceed int.MaxValue ({int.MaxValue:N0}) due to .NET array indexing constraints. " +
                    $"For larger matrices, consider using chunked storage or memory-mapped files.");
            if (value <= 0)
                throw new ArgumentOutOfRangeException(nameof(value), "MaxTotalElements must be positive.");
            _maxTotalElements = value;
        }
    }

    private const double OverflowSafetyThreshold = 1e154;

    #endregion

    #region Configuration

    /// <summary>
    ///     When enabled, arithmetic operations validate results for NaN/Infinity.
    ///     Default: false (disabled for performance). Automatically enabled in DEBUG builds.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         <strong>WARNING: Global Mutable State</strong>
    ///         This is a static property that affects ALL matrix operations in the current process.
    ///         If one thread enables this for debugging, it impacts ALL other threads.
    ///     </para>
    ///     <para>
    ///         For thread-safe validation control, use <see cref="ValidateFiniteValuesForCurrentThread" />
    ///         which uses thread-local storage and only affects the current thread.
    ///     </para>
    ///     <para>
    ///         Enable this in production code when debugging numerical instabilities.
    ///         Has ~2× performance impact on arithmetic operations when enabled.
    ///     </para>
    /// </remarks>
    public static bool ValidateFiniteValues { get; set; } =
#if DEBUG
        true;
#else
        false;
#endif

    /// <summary>
    ///     Thread-local validation setting. When set, overrides <see cref="ValidateFiniteValues" />
    ///     for the current thread only.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         This property uses thread-local storage, so enabling validation on one thread
    ///         does not affect other threads. This is safer than the global <see cref="ValidateFiniteValues" />
    ///         property for multi-threaded applications.
    ///     </para>
    ///     <para>
    ///         When not explicitly set (null), falls back to the global <see cref="ValidateFiniteValues" /> setting.
    ///     </para>
    ///     <para>
    ///         Usage example:
    ///         <code>
    ///         // Enable validation just for this thread
    ///         Matrix.ValidateFiniteValuesForCurrentThread = true;
    ///         try {
    ///             // ... matrix operations with validation ...
    ///         } finally {
    ///             Matrix.ValidateFiniteValuesForCurrentThread = null; // Reset to global
    ///         }
    ///         </code>
    ///     </para>
    /// </remarks>
    [ThreadStatic] private static bool? _threadLocalValidateFinite;

    /// <summary>
    ///     Gets or sets the thread-local validation setting.
    ///     Set to null to use the global <see cref="ValidateFiniteValues" /> setting.
    /// </summary>
    public static bool? ValidateFiniteValuesForCurrentThread
    {
        get => _threadLocalValidateFinite;
        set => _threadLocalValidateFinite = value;
    }

    /// <summary>
    ///     Gets the effective validation setting for the current thread.
    ///     Returns the thread-local setting if set, otherwise the global setting.
    /// </summary>
    internal static bool EffectiveValidateFiniteValues =>
        _threadLocalValidateFinite ?? ValidateFiniteValues;

    #endregion

    #region Fields and Properties

    [JsonInclude] internal readonly double[] _data;
    [JsonInclude] public int RowCount { get; }
    [JsonInclude] public int ColumnCount { get; }

    public int ElementCount => RowCount * ColumnCount;
    public bool IsSquare => RowCount == ColumnCount;
    public bool IsScalar => RowCount == 1 && ColumnCount == 1;
    public int Length => ElementCount;

    #endregion

    #region Helper Methods

    /// <summary>
    ///     Validates that matrix dimensions do not exceed safety limits to prevent DoS attacks.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ValidateMatrixSize(int rows, int columns)
    {
        if (rows > MaxMatrixDimension)
            throw new ArgumentOutOfRangeException(nameof(rows),
                $"Row count ({rows:N0}) exceeds maximum allowed dimension ({MaxMatrixDimension:N0})");

        if (columns > MaxMatrixDimension)
            throw new ArgumentOutOfRangeException(nameof(columns),
                $"Column count ({columns:N0}) exceeds maximum allowed dimension ({MaxMatrixDimension:N0})");

        // CRITICAL FIX 3: Cast FIRST operand to long to ensure long multiplication
        // Without this, "rows * columns" could overflow before the cast in some cases
        var totalElements = (long)rows * columns;
        if (totalElements > MaxTotalElements)
            throw new ArgumentOutOfRangeException(nameof(rows),
                $"Total matrix elements ({totalElements:N0}) exceeds maximum ({MaxTotalElements:N0}). " +
                "Consider using sparse matrix representations for very large matrices.");
    }

    /// <summary>
    ///     Validates that array elements are finite (not NaN or Infinity) when validation is enabled.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ValidateFinite(double[] data, string operation)
    {
        if (!EffectiveValidateFiniteValues) return;

        for (var i = 0; i < data.Length; i++)
            if (!double.IsFinite(data[i]))
                throw new ArithmeticException(
                    $"Operation '{operation}' produced non-finite value at index {i}: " +
                    $"{(double.IsNaN(data[i]) ? "NaN" : "Infinity")}. " +
                    "This usually indicates numerical overflow, underflow, or division by zero.");
    }

    #endregion

    #region Constructors

    /// <summary>
    ///     Creates a new, zero-initialized matrix.
    /// </summary>
    /// <remarks>
    ///     <para>
    ///         Uses GC-managed arrays with optional pinning for large matrices (>1000 elements).
    ///         Pinning reduces GC overhead during intensive computation.
    ///     </para>
    ///     <para>
    ///         <strong>SIMD Alignment Note:</strong> Arrays may not be aligned to 32/64-byte boundaries
    ///         required for optimal AVX2/AVX-512 performance. SIMD code uses LoadVector/LoadUnsafe methods
    ///         which handle unaligned data correctly with ~10-20% performance penalty. For applications
    ///         requiring maximum SIMD performance, consider specialized libraries or custom allocators.
    ///     </para>
    /// </remarks>
    public Matrix(int rows, int columns)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(rows);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(columns);

        ValidateMatrixSize(rows, columns);

        if (rows > MaxRecommendedSize || columns > MaxRecommendedSize)
            Debug.WriteLine(
                $"Warning: Creating {rows}×{columns} matrix. Optimized for ≤{MaxRecommendedSize}×{MaxRecommendedSize}.");

        RowCount = rows;
        ColumnCount = columns;

        // FIXED: Use checked arithmetic to catch any overflow that slips past validation
        int size;
        try
        {
            size = checked(rows * columns);
        }
        catch (OverflowException)
        {
            throw new ArgumentOutOfRangeException(nameof(rows),
                $"Matrix dimensions {rows}×{columns} would overflow int.MaxValue ({int.MaxValue:N0}). " +
                "Reduce dimensions or use sparse/chunked storage.");
        }

        _data = GC.AllocateArray<double>(size, size > 1000);
    }


    /// <summary>
    ///     Creates a Matrix from a jagged array (row-major input).
    /// </summary>
    /// <param name="source">The source jagged array in row-major format.</param>
    /// <exception cref="ArgumentNullException">Thrown when source is null.</exception>
    /// <exception cref="ArgumentException">
    ///     Thrown when source is empty, contains null rows, has inconsistent row lengths,
    ///     or contains NaN/Infinity values.
    /// </exception>
    /// <remarks>
    ///     <para>
    ///         Input data is provided in standard row-major format and automatically converted
    ///         to the internal column-major storage during construction.
    ///     </para>
    ///     <para>
    ///         Time Complexity: O(m×n) where m = rows, n = columns.
    ///     </para>
    /// </remarks>
    public Matrix(double[][] source)
    {
        ArgumentNullException.ThrowIfNull(source);
        ArgumentOutOfRangeException.ThrowIfZero(source.Length, nameof(source));

        RowCount = source.Length;
        ColumnCount = source[0]?.Length ?? throw new ArgumentException("First row cannot be null", nameof(source));

        ValidateMatrixSize(RowCount, ColumnCount);

        _data = GC.AllocateUninitializedArray<double>(RowCount * ColumnCount);

        for (var i = 0; i < RowCount; i++)
        {
            if (source[i]?.Length != ColumnCount)
                throw new ArgumentException("All rows must have the same length", nameof(source));

            for (var j = 0; j < ColumnCount; j++)
            {
                var value = source[i][j];

                if (!double.IsFinite(value))
                    throw new ArgumentException(
                        $"Matrix element at [{i},{j}] is {(double.IsNaN(value) ? "NaN" : "Infinity")}. " +
                        "Only finite values are allowed.", nameof(source));

                _data[j * RowCount + i] = value;
            }
        }
    }

    internal Matrix(int rows, int columns, double[] data)
    {
        RowCount = rows;
        ColumnCount = columns;
        _data = data;
    }

    #endregion

    #region Indexers

    public double this[int row, int column]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if ((uint)row >= (uint)RowCount)
                ThrowHelper.ThrowArgumentOutOfRangeException(nameof(row));
            if ((uint)column >= (uint)ColumnCount)
                ThrowHelper.ThrowArgumentOutOfRangeException(nameof(column));
            return _data[column * RowCount + row];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            if ((uint)row >= (uint)RowCount)
                ThrowHelper.ThrowArgumentOutOfRangeException(nameof(row));
            if ((uint)column >= (uint)ColumnCount)
                ThrowHelper.ThrowArgumentOutOfRangeException(nameof(column));

            if (EffectiveValidateFiniteValues && !double.IsFinite(value))
                throw new ArgumentException(
                    $"Non-finite value ({value}) detected at [{row},{column}]. " +
                    "Set Matrix.ValidateFiniteValues = false to disable this check.", nameof(value));

            _data[column * RowCount + row] = value;
        }
    }

    internal double this[int index]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _data[index];
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set => _data[index] = value;
    }

    #endregion

    #region Factory Methods

    public static Matrix Identity(int size)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(size);
        var result = new Matrix(size, size);
        for (var i = 0; i < size; i++)
            result._data[i * size + i] = 1.0;
        return result;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Matrix Zeros(int rows, int columns)
    {
        return new Matrix(rows, columns);
    }

    public static Matrix Ones(int rows, int columns)
    {
        var result = new Matrix(rows, columns);
        Array.Fill(result._data, 1.0);
        return result;
    }

    /// <summary>
    ///     Creates a unit column vector of specified length with a 1 at the given position.
    /// </summary>
    public static Vector UnitVector(int length, int position)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(length);
        ArgumentOutOfRangeException.ThrowIfNegative(position);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(position, length);

        var v = new Vector(length);
        v._data[position] = 1.0;
        return v;
    }

    public static Matrix Diagonal(params double[] values)
    {
        ArgumentOutOfRangeException.ThrowIfZero(values.Length);
        var n = values.Length;
        var result = new Matrix(n, n);
        for (var i = 0; i < n; i++)
            result._data[i * n + i] = values[i];
        return result;
    }

    public static Matrix Diagonal(ReadOnlySpan<double> values)
    {
        ArgumentOutOfRangeException.ThrowIfZero(values.Length);
        var n = values.Length;
        var result = new Matrix(n, n);
        for (var i = 0; i < n; i++)
            result._data[i * n + i] = values[i];
        return result;
    }

    public static Matrix Random(int rows, int columns, int? seed = null)
    {
        var random = seed.HasValue ? new Random(seed.Value) : System.Random.Shared;
        var result = new Matrix(rows, columns);
        for (var i = 0; i < result._data.Length; i++)
            result._data[i] = random.NextDouble();
        return result;
    }

    /// <summary>
    ///     Creates a matrix with normally distributed random values using Box-Muller transform.
    /// </summary>
    /// <param name="rows">Number of rows.</param>
    /// <param name="columns">Number of columns.</param>
    /// <param name="mean">Mean of the normal distribution (default: 0).</param>
    /// <param name="stdDev">Standard deviation of the normal distribution (default: 1).</param>
    /// <param name="seed">Optional random seed for reproducibility.</param>
    /// <returns>A matrix with normally distributed random values.</returns>
    /// <remarks>
    ///     <para>
    ///         When no seed is provided, uses <see cref="Random.Shared" /> which is thread-safe
    ///         but may produce statistically correlated values when called concurrently from
    ///         multiple threads. For concurrent use, provide explicit seeds per thread.
    ///     </para>
    /// </remarks>
    public static Matrix RandomNormal(int rows, int columns, double mean = 0, double stdDev = 1, int? seed = null)
    {
        var random = seed.HasValue ? new Random(seed.Value) : System.Random.Shared;
        var result = new Matrix(rows, columns);

        for (var i = 0; i < result._data.Length; i += 2)
        {
            var u1 = 1.0 - random.NextDouble();
            var u2 = random.NextDouble();
            var radius = Math.Sqrt(-2.0 * Math.Log(u1));
            var theta = 2.0 * Math.PI * u2;

            result._data[i] = mean + stdDev * radius * Math.Cos(theta);
            if (i + 1 < result._data.Length)
                result._data[i + 1] = mean + stdDev * radius * Math.Sin(theta);
        }

        return result;
    }

    #endregion

    #region Utility Methods

    public Matrix Clone()
    {
        var clone = new Matrix(RowCount, ColumnCount);
        Array.Copy(_data, clone._data, _data.Length);
        return clone;
    }

    object ICloneable.Clone()
    {
        return Clone();
    }

    public double[,] ToArray()
    {
        var result = new double[RowCount, ColumnCount];

        // Use the public static ParallelThreshold property
        if (ElementCount >= ParallelThreshold)
            // PARALLEL VERSION (SAFE)
            // We parallelize the outer loop (over columns).
            // Each thread (j) writes to a different, non-overlapping
            // part of the 'result' array, so this is thread-safe.
            // We use array access here instead of pointers to avoid the CS1764 error.
            Parallel.For(0, ColumnCount, ParallelConfig.Options, j =>
            {
                for (var i = 0; i < RowCount; i++)
                    // Read from column-major _data
                    // Write to row-major result
                    result[i, j] = _data[j * RowCount + i];
            });
        else
            // SCALAR VERSION (UNSAFE)
            // This is safe because it's single-threaded.
            unsafe
            {
                fixed (double* pSrc = _data, pDst = &result[0, 0])
                {
                    for (var j = 0; j < ColumnCount; j++)
                    for (var i = 0; i < RowCount; i++)
                        pDst[i * ColumnCount + j] = pSrc[j * RowCount + i];
                }
            }

        return result;
    }

    public bool Equals(Matrix? other)
    {
        return Equals(other, DefaultTolerance);
    }

    public bool Equals(Matrix? other, double tolerance)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        if (RowCount != other.RowCount || ColumnCount != other.ColumnCount) return false;

        for (var i = 0; i < _data.Length; i++)
            if (Math.Abs(_data[i] - other._data[i]) > tolerance)
                return false;

        return true;
    }

    public override bool Equals(object? obj)
    {
        return Equals(obj as Matrix);
    }

    /// <summary>
    ///     Returns a hash code for this matrix instance.
    /// </summary>
    /// <returns>A hash code based on matrix dimensions and a sample of elements.</returns>
    /// <remarks>
    ///     <para>
    ///         <strong>WARNING:</strong> This implementation does not account for floating-point tolerance.
    ///         Two matrices that are considered equal by <see cref="Equals(Matrix)" /> (which uses tolerance)
    ///         may have different hash codes.
    ///     </para>
    ///     <para>
    ///         <strong>Recommendation:</strong> Do not use <see cref="Matrix" /> instances as dictionary keys
    ///         or in hash-based collections if tolerance-based equality semantics are required.
    ///     </para>
    ///     <para>
    ///         Samples up to 20 elements distributed throughout the matrix for better hash distribution.
    ///     </para>
    /// </remarks>
    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.Add(RowCount);
        hash.Add(ColumnCount);

        // Sample up to 20 elements distributed throughout the matrix
        const int maxSamples = 20;

        if (_data.Length <= maxSamples)
        {
            // Small matrix: include all elements
            for (var i = 0; i < _data.Length; i++)
                hash.Add(_data[i]);
        }
        else
        {
            // Large matrix: sample evenly distributed elements
            var stride = _data.Length / maxSamples;
            for (var i = 0; i < maxSamples; i++)
                hash.Add(_data[i * stride]);
        }

        return hash.ToHashCode();
    }

    public override string ToString()
    {
        return ToString("G", null);
    }

    #region Span-Based APIs

    /// <summary>
    ///     Gets a read-only span over the internal data array in column-major order.
    /// </summary>
    /// <remarks>
    ///     The data is stored in column-major format where element at [row, col] is at _data[col * RowCount + row].
    /// </remarks>
    /// <returns>A read-only span over the matrix data.</returns>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<double> AsSpan()
    {
        return _data.AsSpan();
    }

    /// <summary>
    ///     Gets a read-only span for the specified column.
    /// </summary>
    /// <remarks>
    ///     Column data is stored contiguously in the internal array, making this operation very efficient.
    /// </remarks>
    /// <param name="column">The zero-based column index.</param>
    /// <returns>A read-only span over the column data.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when column index is out of range.</exception>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<double> GetColumnSpan(int column)
    {
        if ((uint)column >= (uint)ColumnCount)
            ThrowHelper.ThrowArgumentOutOfRangeException(nameof(column));
        return _data.AsSpan(column * RowCount, RowCount);
    }

    /// <summary>
    ///     Copies a row into the provided span.
    /// </summary>
    /// <remarks>
    ///     Since the matrix is stored in column-major format, this operation requires iterating through columns.
    /// </remarks>
    /// <param name="row">The zero-based row index.</param>
    /// <param name="destination">The destination span to copy into.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when row index is out of range.</exception>
    /// <exception cref="ArgumentException">Thrown when destination span is too small.</exception>
    public void CopyRowTo(int row, Span<double> destination)
    {
        if ((uint)row >= (uint)RowCount)
            ThrowHelper.ThrowArgumentOutOfRangeException(nameof(row));
        if (destination.Length < ColumnCount)
            ThrowHelper.ThrowArgumentException("Destination span is too small.", nameof(destination));

        if (Vector256.IsHardwareAccelerated && ColumnCount >= 4)
        {
            var j = 0;
            for (; j <= ColumnCount - 4; j += 4)
            {
                var v = Vector256.Create(
                    _data[j * RowCount + row],
                    _data[(j + 1) * RowCount + row],
                    _data[(j + 2) * RowCount + row],
                    _data[(j + 3) * RowCount + row]
                );
                v.StoreUnsafe(ref destination[j]);
            }

            for (; j < ColumnCount; j++)
                destination[j] = _data[j * RowCount + row];
        }
        else
        {
            for (var j = 0; j < ColumnCount; j++)
                destination[j] = _data[j * RowCount + row];
        }
    }

    /// <summary>
    ///     Copies a column into the provided span.
    /// </summary>
    /// <remarks>
    ///     Since the matrix is stored in column-major format, this is a very efficient operation.
    /// </remarks>
    /// <param name="column">The zero-based column index.</param>
    /// <param name="destination">The destination span to copy into.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when column index is out of range.</exception>
    /// <exception cref="ArgumentException">Thrown when destination span is too small.</exception>
    public void CopyColumnTo(int column, Span<double> destination)
    {
        if ((uint)column >= (uint)ColumnCount)
            ThrowHelper.ThrowArgumentOutOfRangeException(nameof(column));
        if (destination.Length < RowCount)
            ThrowHelper.ThrowArgumentException("Destination span is too small.", nameof(destination));

        var source = GetColumnSpan(column);
        source.CopyTo(destination);
    }

    #endregion


    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        format ??= "G";
        formatProvider ??= CultureInfo.CurrentCulture;

        var sb = new StringBuilder();
        sb.AppendLine($"Matrix({RowCount}×{ColumnCount}):");

        // Calculate maximum width needed for proper alignment
        var maxWidth = 10; // Minimum width
        var displayRows = Math.Min(RowCount, 10);
        var displayCols = Math.Min(ColumnCount, 10);

        for (var i = 0; i < displayRows; i++)
        for (var j = 0; j < displayCols; j++)
        {
            var str = this[i, j].ToString(format, formatProvider);
            maxWidth = Math.Max(maxWidth, str.Length);
        }

        // Cap at reasonable maximum to prevent excessive width
        maxWidth = Math.Min(maxWidth, 20);

        for (var i = 0; i < displayRows; i++)
        {
            sb.Append("[ ");
            for (var j = 0; j < displayCols; j++)
            {
                if (j > 0) sb.Append("  ");
                sb.Append(this[i, j].ToString(format, formatProvider).PadLeft(maxWidth));
            }

            if (ColumnCount > 10) sb.Append(" ...");
            sb.AppendLine(" ]");
        }

        if (RowCount > 10)
            sb.AppendLine("...");

        return sb.ToString();
    }

    #endregion

    #region Arithmetic Operations

    public static Matrix operator +(Matrix left, Matrix right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        if (left.RowCount != right.RowCount || left.ColumnCount != right.ColumnCount)
            throw new ArgumentException("Matrix dimensions must match");

        var result = new Matrix(left.RowCount, left.ColumnCount);
        var length = left._data.Length;

        if (length >= ParallelThreshold)
        {
            if (Avx512F.IsSupported && length >= 8)
                Parallel.For(0, (length + 7) / 8, ParallelConfig.Options, i =>
                {
                    var start = i * 8;
                    var end = Math.Min(start + 8, length);

                    if (end - start >= 8)
                        unsafe
                        {
                            fixed (double* pL = &left._data[start], pR = &right._data[start], pRes =
                                       &result._data[start])
                            {
                                var vL = Avx512F.LoadVector512(pL);
                                var vR = Avx512F.LoadVector512(pR);
                                Avx512F.Store(pRes, Avx512F.Add(vL, vR));
                            }
                        }
                    else
                        for (var j = start; j < end; j++)
                            result._data[j] = left._data[j] + right._data[j];
                });
            else if (Vector256.IsHardwareAccelerated && length >= 4)
                Parallel.For(0, (length + 3) / 4, ParallelConfig.Options, i =>
                {
                    var start = i * 4;
                    var end = Math.Min(start + 4, length);

                    if (end - start >= 4)
                    {
                        var vL = Vector256.LoadUnsafe(ref left._data[start]);
                        var vR = Vector256.LoadUnsafe(ref right._data[start]);
                        (vL + vR).StoreUnsafe(ref result._data[start]);
                    }
                    else
                    {
                        for (var j = start; j < end; j++)
                            result._data[j] = left._data[j] + right._data[j];
                    }
                });
            else
                Parallel.For(0, length, ParallelConfig.Options, i => result._data[i] = left._data[i] + right._data[i]);
        }
        else
        {
            if (Avx512F.IsSupported && length >= 8)
            {
                unsafe
                {
                    fixed (double* pL = left._data, pR = right._data, pRes = result._data)
                    {
                        const int Vector512Size = 8;
                        var i = 0;
                        for (; i + Vector512Size <= length; i += Vector512Size)
                        {
                            var vL = Avx512F.LoadVector512(pL + i);
                            var vR = Avx512F.LoadVector512(pR + i);
                            Avx512F.Store(pRes + i, Avx512F.Add(vL, vR));
                        }

                        for (; i < length; i++)
                            pRes[i] = pL[i] + pR[i];
                    }
                }
            }
            else if (Vector256.IsHardwareAccelerated && length >= 4)
            {
                const int Vector256Size = 4;
                var i = 0;
                for (; i + Vector256Size <= length; i += Vector256Size)
                {
                    var vL = Vector256.LoadUnsafe(ref left._data[i]);
                    var vR = Vector256.LoadUnsafe(ref right._data[i]);
                    (vL + vR).StoreUnsafe(ref result._data[i]);
                }

                for (; i < length; i++)
                    result._data[i] = left._data[i] + right._data[i];
            }
            else
            {
                for (var i = 0; i < length; i++)
                    result._data[i] = left._data[i] + right._data[i];
            }
        }

        ValidateFinite(result._data, "addition");
        return result;
    }

    public static Matrix operator -(Matrix left, Matrix right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        if (left.RowCount != right.RowCount || left.ColumnCount != right.ColumnCount)
            throw new ArgumentException("Matrix dimensions must match");

        var result = new Matrix(left.RowCount, left.ColumnCount);
        var length = left._data.Length;

        if (length >= ParallelThreshold)
        {
            if (Avx512F.IsSupported && length >= 8)
                Parallel.For(0, (length + 7) / 8, ParallelConfig.Options, i =>
                {
                    var start = i * 8;
                    var end = Math.Min(start + 8, length);

                    if (end - start >= 8)
                        unsafe
                        {
                            fixed (double* pL = &left._data[start], pR = &right._data[start], pRes =
                                       &result._data[start])
                            {
                                var vL = Avx512F.LoadVector512(pL);
                                var vR = Avx512F.LoadVector512(pR);
                                Avx512F.Store(pRes, Avx512F.Subtract(vL, vR));
                            }
                        }
                    else
                        for (var j = start; j < end; j++)
                            result._data[j] = left._data[j] - right._data[j];
                });
            else if (Vector256.IsHardwareAccelerated && length >= 4)
                Parallel.For(0, (length + 3) / 4, ParallelConfig.Options, i =>
                {
                    var start = i * 4;
                    var end = Math.Min(start + 4, length);

                    if (end - start >= 4)
                    {
                        var vL = Vector256.LoadUnsafe(ref left._data[start]);
                        var vR = Vector256.LoadUnsafe(ref right._data[start]);
                        (vL - vR).StoreUnsafe(ref result._data[start]);
                    }
                    else
                    {
                        for (var j = start; j < end; j++)
                            result._data[j] = left._data[j] - right._data[j];
                    }
                });
            else
                Parallel.For(0, length, ParallelConfig.Options, i => result._data[i] = left._data[i] - right._data[i]);
        }
        else
        {
            if (Avx512F.IsSupported && length >= 8)
            {
                unsafe
                {
                    fixed (double* pL = left._data, pR = right._data, pRes = result._data)
                    {
                        var i = 0;
                        for (; i <= length - 8; i += 8)
                        {
                            var vL = Avx512F.LoadVector512(pL + i);
                            var vR = Avx512F.LoadVector512(pR + i);
                            Avx512F.Store(pRes + i, Avx512F.Subtract(vL, vR));
                        }

                        for (; i < length; i++)
                            pRes[i] = pL[i] - pR[i];
                    }
                }
            }
            else if (Vector256.IsHardwareAccelerated && length >= 4)
            {
                var i = 0;
                for (; i <= length - 4; i += 4)
                {
                    var vL = Vector256.LoadUnsafe(ref left._data[i]);
                    var vR = Vector256.LoadUnsafe(ref right._data[i]);
                    (vL - vR).StoreUnsafe(ref result._data[i]);
                }

                for (; i < length; i++)
                    result._data[i] = left._data[i] - right._data[i];
            }
            else
            {
                for (var i = 0; i < length; i++)
                    result._data[i] = left._data[i] - right._data[i];
            }
        }

        ValidateFinite(result._data, "subtraction");
        return result;
    }

    public static Matrix operator *(double scalar, Matrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);

        var result = new Matrix(matrix.RowCount, matrix.ColumnCount);
        var length = matrix._data.Length;

        if (length >= ParallelThreshold)
        {
            if (Avx512F.IsSupported && length >= 8)
            {
                Parallel.For(0, (length + 7) / 8, ParallelConfig.Options, i =>
                {
                    var start = i * 8;
                    var end = Math.Min(start + 8, length);

                    if (end - start >= 8)
                        unsafe
                        {
                            fixed (double* pM = &matrix._data[start], pRes = &result._data[start])
                            {
                                var vScalar = Vector512.Create(scalar);
                                var vM = Avx512F.LoadVector512(pM);
                                Avx512F.Store(pRes, Avx512F.Multiply(vM, vScalar));
                            }
                        }
                    else
                        for (var j = start; j < end; j++)
                            result._data[j] = matrix._data[j] * scalar;
                });
            }
            else if (Vector256.IsHardwareAccelerated && length >= 4)
            {
                var vScalar = Vector256.Create(scalar);
                Parallel.For(0, (length + 3) / 4, ParallelConfig.Options, i =>
                {
                    var start = i * 4;
                    var end = Math.Min(start + 4, length);

                    if (end - start >= 4)
                    {
                        var vM = Vector256.LoadUnsafe(ref matrix._data[start]);
                        (vM * vScalar).StoreUnsafe(ref result._data[start]);
                    }
                    else
                    {
                        for (var j = start; j < end; j++)
                            result._data[j] = matrix._data[j] * scalar;
                    }
                });
            }
            else
            {
                Parallel.For(0, length, ParallelConfig.Options, i => result._data[i] = matrix._data[i] * scalar);
            }
        }
        else
        {
            if (Avx512F.IsSupported && length >= 8)
            {
                unsafe
                {
                    fixed (double* pM = matrix._data, pRes = result._data)
                    {
                        var vScalar = Vector512.Create(scalar);
                        var i = 0;

                        for (; i <= length - 8; i += 8)
                        {
                            var vM = Avx512F.LoadVector512(pM + i);
                            Avx512F.Store(pRes + i, Avx512F.Multiply(vM, vScalar));
                        }

                        for (; i < length; i++)
                            pRes[i] = pM[i] * scalar;
                    }
                }
            }
            else if (Vector256.IsHardwareAccelerated && length >= 4)
            {
                var vScalar = Vector256.Create(scalar);
                var i = 0;

                for (; i <= length - 4; i += 4)
                {
                    var vM = Vector256.LoadUnsafe(ref matrix._data[i]);
                    (vM * vScalar).StoreUnsafe(ref result._data[i]);
                }

                for (; i < length; i++)
                    result._data[i] = matrix._data[i] * scalar;
            }
            else
            {
                for (var i = 0; i < length; i++)
                    result._data[i] = matrix._data[i] * scalar;
            }
        }

        ValidateFinite(result._data, "scalar multiplication");
        return result;
    }

    public static Matrix operator *(Matrix matrix, double scalar)
    {
        return scalar * matrix;
    }

    public static Matrix operator /(Matrix matrix, double scalar)
    {
        return matrix * (1.0 / scalar);
    }

    public static Matrix operator -(Matrix matrix)
    {
        return -1.0 * matrix;
    }

    /// <summary>
    ///     Performs matrix-matrix multiplication (GEMM) or matrix-vector multiplication.
    /// </summary>
    public static Matrix operator *(Matrix left, Matrix right)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);

        if (left.ColumnCount != right.RowCount)
            throw new ArgumentException(
                $"Cannot multiply {left.RowCount}×{left.ColumnCount} by {right.RowCount}×{right.ColumnCount}");

        var m = left.RowCount;
        var n = right.ColumnCount;
        var p = left.ColumnCount;

        if (right.ColumnCount == 1)
        {
            var vector = (Vector)right;
            var result = new Vector(m);
            var mData = left._data;
            var vData = vector._data;

            for (var j = 0; j < p; j++)
            {
                var v_j = vData[j];
                var colOffset = j * m;

                if (Avx512F.IsSupported && m >= 8)
                {
                    unsafe
                    {
                        fixed (double* pM = &mData[colOffset], pR = result._data)
                        {
                            var vScalar = Vector512.Create(v_j);
                            var i = 0;
                            for (; i <= m - 8; i += 8)
                            {
                                var vM = Avx512F.LoadVector512(pM + i);
                                var vRes = Avx512F.LoadVector512(pR + i);
                                vRes = Avx512F.FusedMultiplyAdd(vM, vScalar, vRes);
                                Avx512F.Store(pR + i, vRes);
                            }

                            for (; i < m; i++)
                                result._data[i] += mData[colOffset + i] * v_j;
                        }
                    }
                }
                else if (Fma.IsSupported && Vector256.IsHardwareAccelerated && m >= 4)
                {
                    // FMA path (fastest)
                    var vScalar = Vector256.Create(v_j);
                    var i = 0;
                    for (; i <= m - 4; i += 4)
                    {
                        var vM = Vector256.LoadUnsafe(ref mData[colOffset + i]);
                        var vRes = Vector256.LoadUnsafe(ref result._data[i]);
                        vRes = Fma.MultiplyAdd(vM, vScalar, vRes);
                        vRes.StoreUnsafe(ref result._data[i]);
                    }

                    for (; i < m; i++)
                        result._data[i] += mData[colOffset + i] * v_j;
                }
                else if (Vector256.IsHardwareAccelerated && m >= 4)
                {
                    // AVX2 without FMA fallback
                    var vScalar = Vector256.Create(v_j);
                    var i = 0;
                    for (; i <= m - 4; i += 4)
                    {
                        var vM = Vector256.LoadUnsafe(ref mData[colOffset + i]);
                        var vRes = Vector256.LoadUnsafe(ref result._data[i]);
                        var product = Vector256.Multiply(vM, vScalar);
                        vRes = Vector256.Add(product, vRes);
                        vRes.StoreUnsafe(ref result._data[i]);
                    }

                    for (; i < m; i++)
                        result._data[i] += mData[colOffset + i] * v_j;
                }
                else
                {
                    for (var i = 0; i < m; i++)
                        result._data[i] += mData[colOffset + i] * v_j;
                }
            }

            return result;
        }

        if (left.IsSquare && right.IsSquare && m == n)
            return m switch
            {
                1 => Multiply1x1(left, right),
                2 => Multiply2x2(left, right),
                3 => Multiply3x3(left, right),
                _ => GemmBlocked(left, false, right, false)
            };

        return GemmBlocked(left, false, right, false);
    }

    /// <summary>
    ///     Performs matrix-vector multiplication: A*v (returns a Vector).
    /// </summary>
    /// <param name="matrix">The matrix (m×n).</param>
    /// <param name="vector">The vector (n×1).</param>
    /// <returns>The resulting vector (m×1).</returns>
    /// <remarks>
    ///     This operator is optimized with SIMD vectorization and provides better type safety
    ///     than converting through Matrix. The computation is y = A*x where y[i] = sum(A[i,j] * x[j]).
    /// </remarks>
    public static Vector operator *(Matrix matrix, Vector vector)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        ArgumentNullException.ThrowIfNull(vector);

        if (matrix.ColumnCount != vector.Length)
            throw new ArgumentException(
                $"Cannot multiply {matrix.RowCount}×{matrix.ColumnCount} matrix by vector of length {vector.Length}");

        var m = matrix.RowCount;
        var n = matrix.ColumnCount;
        var result = new Vector(m);
        var mData = matrix._data;
        var vData = vector._data;

        for (var j = 0; j < n; j++)
        {
            var v_j = vData[j];
            var colOffset = j * m;

            if (Avx512F.IsSupported && m >= 8)
            {
                unsafe
                {
                    fixed (double* pM = &mData[colOffset], pR = result._data)
                    {
                        var vScalar = Vector512.Create(v_j);
                        var i = 0;
                        for (; i <= m - 8; i += 8)
                        {
                            var vM = Avx512F.LoadVector512(pM + i);
                            var vRes = Avx512F.LoadVector512(pR + i);
                            vRes = Avx512F.FusedMultiplyAdd(vM, vScalar, vRes);
                            Avx512F.Store(pR + i, vRes);
                        }

                        for (; i < m; i++)
                            result._data[i] += mData[colOffset + i] * v_j;
                    }
                }
            }
            else if (Fma.IsSupported && Vector256.IsHardwareAccelerated && m >= 4)
            {
                // FMA path (fastest)
                var vScalar = Vector256.Create(v_j);
                var i = 0;
                for (; i <= m - 4; i += 4)
                {
                    var vM = Vector256.LoadUnsafe(ref mData[colOffset + i]);
                    var vRes = Vector256.LoadUnsafe(ref result._data[i]);
                    vRes = Fma.MultiplyAdd(vM, vScalar, vRes);
                    vRes.StoreUnsafe(ref result._data[i]);
                }

                for (; i < m; i++)
                    result._data[i] += mData[colOffset + i] * v_j;
            }
            else if (Vector256.IsHardwareAccelerated && m >= 4)
            {
                // AVX2 without FMA fallback
                var vScalar = Vector256.Create(v_j);
                var i = 0;
                for (; i <= m - 4; i += 4)
                {
                    var vM = Vector256.LoadUnsafe(ref mData[colOffset + i]);
                    var vRes = Vector256.LoadUnsafe(ref result._data[i]);
                    var product = Vector256.Multiply(vM, vScalar);
                    vRes = Vector256.Add(product, vRes);
                    vRes.StoreUnsafe(ref result._data[i]);
                }

                for (; i < m; i++)
                    result._data[i] += mData[colOffset + i] * v_j;
            }
            else
            {
                for (var i = 0; i < m; i++)
                    result._data[i] += mData[colOffset + i] * v_j;
            }
        }

        return result;
    }

    /// <summary>
    ///     Performs vector-matrix multiplication: v*A (row vector times matrix, returns a Vector).
    /// </summary>
    /// <param name="vector">The row vector (1×m).</param>
    /// <param name="matrix">The matrix (m×n).</param>
    /// <returns>The resulting row vector (1×n) as a Vector.</returns>
    /// <remarks>
    ///     This operator treats the vector as a row vector and performs v'*A.
    ///     The computation is y = v'*A where y[j] = sum(v[i] * A[i,j]).
    /// </remarks>
    public static Vector operator *(Vector vector, Matrix matrix)
    {
        ArgumentNullException.ThrowIfNull(vector);
        ArgumentNullException.ThrowIfNull(matrix);

        if (vector.Length != matrix.RowCount)
            throw new ArgumentException(
                $"Cannot multiply vector of length {vector.Length} by {matrix.RowCount}×{matrix.ColumnCount} matrix");

        var m = matrix.RowCount;
        var n = matrix.ColumnCount;
        var result = new Vector(n);
        var mData = matrix._data;
        var vData = vector._data;

        for (var j = 0; j < n; j++)
        {
            var colOffset = j * m;
            double sum = 0;

            if (Avx512F.IsSupported && m >= 8)
            {
                unsafe
                {
                    fixed (double* pM = &mData[colOffset], pV = vData)
                    {
                        var vSum = Vector512<double>.Zero;
                        var i = 0;
                        for (; i <= m - 8; i += 8)
                        {
                            var vM = Avx512F.LoadVector512(pM + i);
                            var vV = Avx512F.LoadVector512(pV + i);
                            vSum = Avx512F.FusedMultiplyAdd(vM, vV, vSum);
                        }

                        sum = Vector512.Sum(vSum);
                        for (; i < m; i++)
                            sum += mData[colOffset + i] * vData[i];
                    }
                }
            }
            else if (Fma.IsSupported && Vector256.IsHardwareAccelerated && m >= 4)
            {
                // FMA path (fastest)
                var vSum = Vector256<double>.Zero;
                var i = 0;
                for (; i <= m - 4; i += 4)
                {
                    var vM = Vector256.LoadUnsafe(ref mData[colOffset + i]);
                    var vV = Vector256.LoadUnsafe(ref vData[i]);
                    vSum = Fma.MultiplyAdd(vM, vV, vSum);
                }

                sum = Vector256.Sum(vSum);
                for (; i < m; i++)
                    sum += mData[colOffset + i] * vData[i];
            }
            else if (Vector256.IsHardwareAccelerated && m >= 4)
            {
                // AVX2 without FMA fallback
                var vSum = Vector256<double>.Zero;
                var i = 0;
                for (; i <= m - 4; i += 4)
                {
                    var vM = Vector256.LoadUnsafe(ref mData[colOffset + i]);
                    var vV = Vector256.LoadUnsafe(ref vData[i]);
                    var product = Vector256.Multiply(vM, vV);
                    vSum = Vector256.Add(product, vSum);
                }

                sum = Vector256.Sum(vSum);
                for (; i < m; i++)
                    sum += mData[colOffset + i] * vData[i];
            }
            else
            {
                for (var i = 0; i < m; i++)
                    sum += mData[colOffset + i] * vData[i];
            }

            result._data[j] = sum;
        }

        return result;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Matrix Multiply1x1(Matrix left, Matrix right)
    {
        return new Matrix(1, 1) { [0, 0] = left._data[0] * right._data[0] };
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Matrix Multiply2x2(Matrix left, Matrix right)
    {
        var result = new Matrix(2, 2);
        var l = left._data;
        var r = right._data;
        var res = result._data;

        res[0] = l[0] * r[0] + l[2] * r[1];
        res[1] = l[1] * r[0] + l[3] * r[1];
        res[2] = l[0] * r[2] + l[2] * r[3];
        res[3] = l[1] * r[2] + l[3] * r[3];

        return result;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Matrix Multiply3x3(Matrix left, Matrix right)
    {
        var result = new Matrix(3, 3);
        var l = left._data;
        var r = right._data;
        var res = result._data;

        for (var j = 0; j < 3; j++)
        {
            var rCol = j * 3;
            var resCol = j * 3;

            res[resCol] = l[0] * r[rCol] + l[3] * r[rCol + 1] + l[6] * r[rCol + 2];
            res[resCol + 1] = l[1] * r[rCol] + l[4] * r[rCol + 1] + l[7] * r[rCol + 2];
            res[resCol + 2] = l[2] * r[rCol] + l[5] * r[rCol + 1] + l[8] * r[rCol + 2];
        }

        return result;
    }

    public static Matrix Multiply(Matrix A, Matrix B)
    {
        return A * B;
    }

    public static Matrix MultiplyAtB(Matrix A, Matrix B)
    {
        return GemmBlocked(A, true, B, false);
    }

    public static Matrix MultiplyABt(Matrix A, Matrix B)
    {
        return GemmBlocked(A, false, B, true);
    }

    public static Matrix MultiplyAtBt(Matrix A, Matrix B)
    {
        return GemmBlocked(A, true, B, true);
    }

    private static Matrix GemmBlocked(Matrix A, bool transA, Matrix B, bool transB)
    {
        var M = transA ? A.ColumnCount : A.RowCount;
        var K = transA ? A.RowCount : A.ColumnCount;
        var N = transB ? B.RowCount : B.ColumnCount;
        var Kb = transB ? B.ColumnCount : B.RowCount;
        if (K != Kb) throw new ArgumentException($"Incompatible sizes: op(A) is {M}×{K}, op(B) is {Kb}×{N}.");

        var C = new Matrix(M, N);

        var totalWork = (long)M * N * K;
        var useParallel = totalWork >= ParallelThreshold * 10 &&
                          M >= MinParallelDimension &&
                          N >= MinParallelDimension;

        if (useParallel)
        {
            Parallel.For(0, N, ParallelConfig.Options, j =>
            {
                if (!transA && !transB) Gemm_A_B_Unblocked_N(A, B, C, j, j + 1);
                else if (transA && !transB) Gemm_At_B_Unblocked_N(A, B, C, j, j + 1);
                else if (!transA && transB) Gemm_A_Bt_Unblocked_N(A, B, C, j, j + 1);
                else Gemm_At_Bt_Unblocked_N(A, B, C, j, j + 1);
            });
        }
        else
        {
            if (!transA && !transB) Gemm_A_B_Unblocked_N(A, B, C, 0, N);
            else if (transA && !transB) Gemm_At_B_Unblocked_N(A, B, C, 0, N);
            else if (!transA && transB) Gemm_A_Bt_Unblocked_N(A, B, C, 0, N);
            else Gemm_At_Bt_Unblocked_N(A, B, C, 0, N);
        }

        return C;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Gemm_A_B_Unblocked_N(Matrix A, Matrix B, Matrix C, int jStart, int jEnd)
    {
        var M = A.RowCount;
        var K = A.ColumnCount;
        var lda = A.RowCount;
        var ldb = B.RowCount;
        var ldc = C.RowCount;
        unsafe
        {
            fixed (double* pA = A._data, pB = B._data, pC = C._data)
            {
                const int NC = MicroKernelN;
                var j = jStart;
                for (; j + NC - 1 < jEnd; j += NC)
                {
                    var c0 = pC + (j + 0) * ldc;
                    var c1 = pC + (j + 1) * ldc;
                    var c2 = pC + (j + 2) * ldc;
                    var c3 = pC + (j + 3) * ldc;
                    for (var k = 0; k < K; ++k)
                    {
                        var aCol = pA + k * lda;
                        var b0 = *(pB + k + (j + 0) * ldb);
                        var b1 = *(pB + k + (j + 1) * ldb);
                        var b2 = *(pB + k + (j + 2) * ldb);
                        var b3 = *(pB + k + (j + 3) * ldb);
                        var i = 0;

                        if (Avx512F.IsSupported)
                        {
                            // AVX-512 path - process 8 doubles at a time
                            var vb0 = Vector512.Create(b0);
                            var vb1 = Vector512.Create(b1);
                            var vb2 = Vector512.Create(b2);
                            var vb3 = Vector512.Create(b3);

                            for (; i + 15 < M; i += 16)
                            {
                                var va0 = Avx512F.LoadVector512(aCol + i + 0);
                                var va1 = Avx512F.LoadVector512(aCol + i + 8);

                                var vc00 = Avx512F.LoadVector512(c0 + i + 0);
                                var vc01 = Avx512F.LoadVector512(c0 + i + 8);
                                vc00 = Avx512F.FusedMultiplyAdd(va0, vb0, vc00);
                                vc01 = Avx512F.FusedMultiplyAdd(va1, vb0, vc01);
                                Avx512F.Store(c0 + i + 0, vc00);
                                Avx512F.Store(c0 + i + 8, vc01);

                                var vc10 = Avx512F.LoadVector512(c1 + i + 0);
                                var vc11 = Avx512F.LoadVector512(c1 + i + 8);
                                vc10 = Avx512F.FusedMultiplyAdd(va0, vb1, vc10);
                                vc11 = Avx512F.FusedMultiplyAdd(va1, vb1, vc11);
                                Avx512F.Store(c1 + i + 0, vc10);
                                Avx512F.Store(c1 + i + 8, vc11);

                                var vc20 = Avx512F.LoadVector512(c2 + i + 0);
                                var vc21 = Avx512F.LoadVector512(c2 + i + 8);
                                vc20 = Avx512F.FusedMultiplyAdd(va0, vb2, vc20);
                                vc21 = Avx512F.FusedMultiplyAdd(va1, vb2, vc21);
                                Avx512F.Store(c2 + i + 0, vc20);
                                Avx512F.Store(c2 + i + 8, vc21);

                                var vc30 = Avx512F.LoadVector512(c3 + i + 0);
                                var vc31 = Avx512F.LoadVector512(c3 + i + 8);
                                vc30 = Avx512F.FusedMultiplyAdd(va0, vb3, vc30);
                                vc31 = Avx512F.FusedMultiplyAdd(va1, vb3, vc31);
                                Avx512F.Store(c3 + i + 0, vc30);
                                Avx512F.Store(c3 + i + 8, vc31);
                            }

                            for (; i + 7 < M; i += 8)
                            {
                                var va = Avx512F.LoadVector512(aCol + i);

                                var r0 = Avx512F.LoadVector512(c0 + i);
                                var r1 = Avx512F.LoadVector512(c1 + i);
                                var r2 = Avx512F.LoadVector512(c2 + i);
                                var r3 = Avx512F.LoadVector512(c3 + i);

                                r0 = Avx512F.FusedMultiplyAdd(va, vb0, r0);
                                r1 = Avx512F.FusedMultiplyAdd(va, vb1, r1);
                                r2 = Avx512F.FusedMultiplyAdd(va, vb2, r2);
                                r3 = Avx512F.FusedMultiplyAdd(va, vb3, r3);

                                Avx512F.Store(c0 + i, r0);
                                Avx512F.Store(c1 + i, r1);
                                Avx512F.Store(c2 + i, r2);
                                Avx512F.Store(c3 + i, r3);
                            }
                        }
                        else if (Fma.IsSupported)
                        {
                            // FMA path (fastest)
                            var vb0 = Avx.BroadcastScalarToVector256(&b0);
                            var vb1 = Avx.BroadcastScalarToVector256(&b1);
                            var vb2 = Avx.BroadcastScalarToVector256(&b2);
                            var vb3 = Avx.BroadcastScalarToVector256(&b3);

                            for (; i + 7 < M; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aCol + i + 0);
                                var va1 = Avx.LoadVector256(aCol + i + 4);

                                var vc00 = Avx.LoadVector256(c0 + i + 0);
                                var vc01 = Avx.LoadVector256(c0 + i + 4);
                                vc00 = Fma.MultiplyAdd(va0, vb0, vc00);
                                vc01 = Fma.MultiplyAdd(va1, vb0, vc01);
                                Avx.Store(c0 + i + 0, vc00);
                                Avx.Store(c0 + i + 4, vc01);

                                var vc10 = Avx.LoadVector256(c1 + i + 0);
                                var vc11 = Avx.LoadVector256(c1 + i + 4);
                                vc10 = Fma.MultiplyAdd(va0, vb1, vc10);
                                vc11 = Fma.MultiplyAdd(va1, vb1, vc11);
                                Avx.Store(c1 + i + 0, vc10);
                                Avx.Store(c1 + i + 4, vc11);

                                var vc20 = Avx.LoadVector256(c2 + i + 0);
                                var vc21 = Avx.LoadVector256(c2 + i + 4);
                                vc20 = Fma.MultiplyAdd(va0, vb2, vc20);
                                vc21 = Fma.MultiplyAdd(va1, vb2, vc21);
                                Avx.Store(c2 + i + 0, vc20);
                                Avx.Store(c2 + i + 4, vc21);

                                var vc30 = Avx.LoadVector256(c3 + i + 0);
                                var vc31 = Avx.LoadVector256(c3 + i + 4);
                                vc30 = Fma.MultiplyAdd(va0, vb3, vc30);
                                vc31 = Fma.MultiplyAdd(va1, vb3, vc31);
                                Avx.Store(c3 + i + 0, vc30);
                                Avx.Store(c3 + i + 4, vc31);
                            }

                            for (; i + 3 < M; i += 4)
                            {
                                var va = Avx.LoadVector256(aCol + i);

                                var r0 = Avx.LoadVector256(c0 + i);
                                var r1 = Avx.LoadVector256(c1 + i);
                                var r2 = Avx.LoadVector256(c2 + i);
                                var r3 = Avx.LoadVector256(c3 + i);

                                r0 = Fma.MultiplyAdd(va, vb0, r0);
                                r1 = Fma.MultiplyAdd(va, vb1, r1);
                                r2 = Fma.MultiplyAdd(va, vb2, r2);
                                r3 = Fma.MultiplyAdd(va, vb3, r3);

                                Avx.Store(c0 + i, r0);
                                Avx.Store(c1 + i, r1);
                                Avx.Store(c2 + i, r2);
                                Avx.Store(c3 + i, r3);
                            }
                        }
                        else if (Avx.IsSupported)
                        {
                            // AVX without FMA
                            var vb0 = Avx.BroadcastScalarToVector256(&b0);
                            var vb1 = Avx.BroadcastScalarToVector256(&b1);
                            var vb2 = Avx.BroadcastScalarToVector256(&b2);
                            var vb3 = Avx.BroadcastScalarToVector256(&b3);

                            for (; i + 7 < M; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aCol + i + 0);
                                var va1 = Avx.LoadVector256(aCol + i + 4);

                                var vc00 = Avx.LoadVector256(c0 + i + 0);
                                var vc01 = Avx.LoadVector256(c0 + i + 4);
                                vc00 = Avx.Add(vc00, Avx.Multiply(va0, vb0));
                                vc01 = Avx.Add(vc01, Avx.Multiply(va1, vb0));
                                Avx.Store(c0 + i + 0, vc00);
                                Avx.Store(c0 + i + 4, vc01);

                                var vc10 = Avx.LoadVector256(c1 + i + 0);
                                var vc11 = Avx.LoadVector256(c1 + i + 4);
                                vc10 = Avx.Add(vc10, Avx.Multiply(va0, vb1));
                                vc11 = Avx.Add(vc11, Avx.Multiply(va1, vb1));
                                Avx.Store(c1 + i + 0, vc10);
                                Avx.Store(c1 + i + 4, vc11);

                                var vc20 = Avx.LoadVector256(c2 + i + 0);
                                var vc21 = Avx.LoadVector256(c2 + i + 4);
                                vc20 = Avx.Add(vc20, Avx.Multiply(va0, vb2));
                                vc21 = Avx.Add(vc21, Avx.Multiply(va1, vb2));
                                Avx.Store(c2 + i + 0, vc20);
                                Avx.Store(c2 + i + 4, vc21);

                                var vc30 = Avx.LoadVector256(c3 + i + 0);
                                var vc31 = Avx.LoadVector256(c3 + i + 4);
                                vc30 = Avx.Add(vc30, Avx.Multiply(va0, vb3));
                                vc31 = Avx.Add(vc31, Avx.Multiply(va1, vb3));
                                Avx.Store(c3 + i + 0, vc30);
                                Avx.Store(c3 + i + 4, vc31);
                            }

                            for (; i + 3 < M; i += 4)
                            {
                                var va = Avx.LoadVector256(aCol + i);

                                var r0 = Avx.LoadVector256(c0 + i);
                                var r1 = Avx.LoadVector256(c1 + i);
                                var r2 = Avx.LoadVector256(c2 + i);
                                var r3 = Avx.LoadVector256(c3 + i);

                                r0 = Avx.Add(r0, Avx.Multiply(va, vb0));
                                r1 = Avx.Add(r1, Avx.Multiply(va, vb1));
                                r2 = Avx.Add(r2, Avx.Multiply(va, vb2));
                                r3 = Avx.Add(r3, Avx.Multiply(va, vb3));

                                Avx.Store(c0 + i, r0);
                                Avx.Store(c1 + i, r1);
                                Avx.Store(c2 + i, r2);
                                Avx.Store(c3 + i, r3);
                            }
                        }

                        for (; i < M; ++i)
                        {
                            var a = aCol[i];
                            c0[i] += a * b0;
                            c1[i] += a * b1;
                            c2[i] += a * b2;
                            c3[i] += a * b3;
                        }
                    }
                }

                for (; j < jEnd; ++j)
                {
                    var c = pC + j * ldc;
                    for (var k = 0; k < K; ++k)
                    {
                        var aCol = pA + k * lda;
                        var bkj = *(pB + k + j * ldb);
                        var i = 0;

                        if (Avx512F.IsSupported)
                        {
                            var vb = Vector512.Create(bkj);
                            for (; i + 15 < M; i += 16)
                            {
                                var va0 = Avx512F.LoadVector512(aCol + i + 0);
                                var va1 = Avx512F.LoadVector512(aCol + i + 8);

                                var vc0 = Avx512F.LoadVector512(c + i + 0);
                                var vc1 = Avx512F.LoadVector512(c + i + 8);
                                vc0 = Avx512F.FusedMultiplyAdd(va0, vb, vc0);
                                vc1 = Avx512F.FusedMultiplyAdd(va1, vb, vc1);
                                Avx512F.Store(c + i + 0, vc0);
                                Avx512F.Store(c + i + 8, vc1);
                            }

                            for (; i + 7 < M; i += 8)
                            {
                                var va = Avx512F.LoadVector512(aCol + i);
                                var vc = Avx512F.LoadVector512(c + i);
                                vc = Avx512F.FusedMultiplyAdd(va, vb, vc);
                                Avx512F.Store(c + i, vc);
                            }
                        }
                        else if (Fma.IsSupported)
                        {
                            var vb = Avx.BroadcastScalarToVector256(&bkj);
                            for (; i + 7 < M; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aCol + i + 0);
                                var va1 = Avx.LoadVector256(aCol + i + 4);

                                var vc0 = Avx.LoadVector256(c + i + 0);
                                var vc1 = Avx.LoadVector256(c + i + 4);
                                vc0 = Fma.MultiplyAdd(va0, vb, vc0);
                                vc1 = Fma.MultiplyAdd(va1, vb, vc1);
                                Avx.Store(c + i + 0, vc0);
                                Avx.Store(c + i + 4, vc1);
                            }

                            for (; i + 3 < M; i += 4)
                            {
                                var va = Avx.LoadVector256(aCol + i);
                                var vc = Avx.LoadVector256(c + i);
                                vc = Fma.MultiplyAdd(va, vb, vc);
                                Avx.Store(c + i, vc);
                            }
                        }
                        else if (Avx.IsSupported)
                        {
                            var vb = Avx.BroadcastScalarToVector256(&bkj);
                            for (; i + 7 < M; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aCol + i + 0);
                                var va1 = Avx.LoadVector256(aCol + i + 4);

                                var vc0 = Avx.LoadVector256(c + i + 0);
                                var vc1 = Avx.LoadVector256(c + i + 4);
                                vc0 = Avx.Add(vc0, Avx.Multiply(va0, vb));
                                vc1 = Avx.Add(vc1, Avx.Multiply(va1, vb));
                                Avx.Store(c + i + 0, vc0);
                                Avx.Store(c + i + 4, vc1);
                            }

                            for (; i + 3 < M; i += 4)
                            {
                                var va = Avx.LoadVector256(aCol + i);
                                var vc = Avx.LoadVector256(c + i);
                                vc = Avx.Add(vc, Avx.Multiply(va, vb));
                                Avx.Store(c + i, vc);
                            }
                        }

                        for (; i < M; ++i) c[i] += aCol[i] * bkj;
                    }
                }
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Gemm_At_B_Unblocked_N(Matrix A, Matrix B, Matrix C, int jStart, int jEnd)
    {
        var M = A.ColumnCount;
        var K = A.RowCount;
        var lda = A.RowCount;
        var ldb = B.RowCount;
        var ldc = C.RowCount;
        unsafe
        {
            fixed (double* pA = A._data, pB = B._data, pC = C._data)
            {
                const int NC = MicroKernelN;
                for (var p = 0; p < M; ++p)
                {
                    var aColP = pA + p * lda;
                    var j = jStart;
                    for (; j + NC - 1 < jEnd; j += NC)
                    {
                        var b0 = pB + (j + 0) * ldb;
                        var b1 = pB + (j + 1) * ldb;
                        var b2 = pB + (j + 2) * ldb;
                        var b3 = pB + (j + 3) * ldb;

                        Vector256<double> acc0 = default, acc1 = default, acc2 = default, acc3 = default;
                        var i = 0;
                        if (Avx.IsSupported)
                        {
                            for (; i + 7 < K; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aColP + i + 0);
                                var va1 = Avx.LoadVector256(aColP + i + 4);

                                var vb0_0 = Avx.LoadVector256(b0 + i + 0);
                                var vb0_1 = Avx.LoadVector256(b0 + i + 4);
                                var vb1_0 = Avx.LoadVector256(b1 + i + 0);
                                var vb1_1 = Avx.LoadVector256(b1 + i + 4);
                                var vb2_0 = Avx.LoadVector256(b2 + i + 0);
                                var vb2_1 = Avx.LoadVector256(b2 + i + 4);
                                var vb3_0 = Avx.LoadVector256(b3 + i + 0);
                                var vb3_1 = Avx.LoadVector256(b3 + i + 4);

                                acc0 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb0_0, acc0)
                                    : Avx.Add(acc0, Avx.Multiply(va0, vb0_0));
                                acc0 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb0_1, acc0)
                                    : Avx.Add(acc0, Avx.Multiply(va1, vb0_1));
                                acc1 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb1_0, acc1)
                                    : Avx.Add(acc1, Avx.Multiply(va0, vb1_0));
                                acc1 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb1_1, acc1)
                                    : Avx.Add(acc1, Avx.Multiply(va1, vb1_1));
                                acc2 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb2_0, acc2)
                                    : Avx.Add(acc2, Avx.Multiply(va0, vb2_0));
                                acc2 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb2_1, acc2)
                                    : Avx.Add(acc2, Avx.Multiply(va1, vb2_1));
                                acc3 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb3_0, acc3)
                                    : Avx.Add(acc3, Avx.Multiply(va0, vb3_0));
                                acc3 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb3_1, acc3)
                                    : Avx.Add(acc3, Avx.Multiply(va1, vb3_1));
                            }

                            for (; i + 3 < K; i += 4)
                            {
                                var va = Avx.LoadVector256(aColP + i);
                                var vb0v = Avx.LoadVector256(b0 + i);
                                var vb1v = Avx.LoadVector256(b1 + i);
                                var vb2v = Avx.LoadVector256(b2 + i);
                                var vb3v = Avx.LoadVector256(b3 + i);
                                acc0 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, vb0v, acc0)
                                    : Avx.Add(acc0, Avx.Multiply(va, vb0v));
                                acc1 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, vb1v, acc1)
                                    : Avx.Add(acc1, Avx.Multiply(va, vb1v));
                                acc2 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, vb2v, acc2)
                                    : Avx.Add(acc2, Avx.Multiply(va, vb2v));
                                acc3 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, vb3v, acc3)
                                    : Avx.Add(acc3, Avx.Multiply(va, vb3v));
                            }
                        }

                        double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
                        if (Avx.IsSupported)
                        {
                            s0 = Vector256.Sum(acc0);
                            s1 = Vector256.Sum(acc1);
                            s2 = Vector256.Sum(acc2);
                            s3 = Vector256.Sum(acc3);
                        }

                        for (; i < K; ++i)
                        {
                            var ai = aColP[i];
                            s0 += ai * b0[i];
                            s1 += ai * b1[i];
                            s2 += ai * b2[i];
                            s3 += ai * b3[i];
                        }

                        pC[p + (j + 0) * ldc] = s0;
                        pC[p + (j + 1) * ldc] = s1;
                        pC[p + (j + 2) * ldc] = s2;
                        pC[p + (j + 3) * ldc] = s3;
                    }

                    for (; j < jEnd; ++j)
                    {
                        var b = pB + j * ldb;
                        Vector256<double> acc = default;
                        var i = 0;
                        if (Avx.IsSupported)
                        {
                            for (; i + 7 < K; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aColP + i + 0);
                                var va1 = Avx.LoadVector256(aColP + i + 4);
                                var vb0 = Avx.LoadVector256(b + i + 0);
                                var vb1 = Avx.LoadVector256(b + i + 4);
                                acc = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb0, acc)
                                    : Avx.Add(acc, Avx.Multiply(va0, vb0));
                                acc = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb1, acc)
                                    : Avx.Add(acc, Avx.Multiply(va1, vb1));
                            }

                            for (; i + 3 < K; i += 4)
                            {
                                var va = Avx.LoadVector256(aColP + i);
                                var vb = Avx.LoadVector256(b + i);
                                acc = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, vb, acc)
                                    : Avx.Add(acc, Avx.Multiply(va, vb));
                            }
                        }

                        var s = Avx.IsSupported ? Vector256.Sum(acc) : 0.0;
                        for (; i < K; ++i) s += aColP[i] * b[i];
                        pC[p + j * ldc] = s;
                    }
                }
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Gemm_A_Bt_Unblocked_N(Matrix A, Matrix B, Matrix C, int jStart, int jEnd)
    {
        var M = A.RowCount;
        var K = A.ColumnCount;
        var lda = A.RowCount;
        var ldb = B.RowCount;
        var ldc = C.RowCount;
        unsafe
        {
            fixed (double* pA = A._data, pB = B._data, pC = C._data)
            {
                const int NC = MicroKernelN;
                var j = jStart;
                for (; j + NC - 1 < jEnd; j += NC)
                {
                    var c0 = pC + (j + 0) * ldc;
                    var c1 = pC + (j + 1) * ldc;
                    var c2 = pC + (j + 2) * ldc;
                    var c3 = pC + (j + 3) * ldc;

                    for (var k = 0; k < K; ++k)
                    {
                        var aCol = pA + k * lda;
                        var b0 = *(pB + k * ldb + (j + 0));
                        var b1 = *(pB + k * ldb + (j + 1));
                        var b2 = *(pB + k * ldb + (j + 2));
                        var b3 = *(pB + k * ldb + (j + 3));

                        var i = 0;
                        if (Avx.IsSupported)
                        {
                            var vb0 = Avx.BroadcastScalarToVector256(&b0);
                            var vb1 = Avx.BroadcastScalarToVector256(&b1);
                            var vb2 = Avx.BroadcastScalarToVector256(&b2);
                            var vb3 = Avx.BroadcastScalarToVector256(&b3);

                            for (; i + 7 < M; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aCol + i + 0);
                                var va1 = Avx.LoadVector256(aCol + i + 4);

                                var vc00 = Avx.LoadVector256(c0 + i + 0);
                                var vc01 = Avx.LoadVector256(c0 + i + 4);
                                vc00 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb0, vc00)
                                    : Avx.Add(vc00, Avx.Multiply(va0, vb0));
                                vc01 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb0, vc01)
                                    : Avx.Add(vc01, Avx.Multiply(va1, vb0));
                                Avx.Store(c0 + i + 0, vc00);
                                Avx.Store(c0 + i + 4, vc01);

                                var vc10 = Avx.LoadVector256(c1 + i + 0);
                                var vc11 = Avx.LoadVector256(c1 + i + 4);
                                vc10 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb1, vc10)
                                    : Avx.Add(vc10, Avx.Multiply(va0, vb1));
                                vc11 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb1, vc11)
                                    : Avx.Add(vc11, Avx.Multiply(va1, vb1));
                                Avx.Store(c1 + i + 0, vc10);
                                Avx.Store(c1 + i + 4, vc11);

                                var vc20 = Avx.LoadVector256(c2 + i + 0);
                                var vc21 = Avx.LoadVector256(c2 + i + 4);
                                vc20 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb2, vc20)
                                    : Avx.Add(vc20, Avx.Multiply(va0, vb2));
                                vc21 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb2, vc21)
                                    : Avx.Add(vc21, Avx.Multiply(va1, vb2));
                                Avx.Store(c2 + i + 0, vc20);
                                Avx.Store(c2 + i + 4, vc21);

                                var vc30 = Avx.LoadVector256(c3 + i + 0);
                                var vc31 = Avx.LoadVector256(c3 + i + 4);
                                vc30 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb3, vc30)
                                    : Avx.Add(vc30, Avx.Multiply(va0, vb3));
                                vc31 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb3, vc31)
                                    : Avx.Add(vc31, Avx.Multiply(va1, vb3));
                                Avx.Store(c3 + i + 0, vc30);
                                Avx.Store(c3 + i + 4, vc31);
                            }

                            for (; i + 3 < M; i += 4)
                            {
                                var va = Avx.LoadVector256(aCol + i);

                                var r0 = Avx.LoadVector256(c0 + i);
                                var r1 = Avx.LoadVector256(c1 + i);
                                var r2 = Avx.LoadVector256(c2 + i);
                                var r3 = Avx.LoadVector256(c3 + i);

                                r0 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, Avx.BroadcastScalarToVector256(&b0), r0)
                                    : Avx.Add(r0, Avx.Multiply(va, Avx.BroadcastScalarToVector256(&b0)));
                                r1 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, Avx.BroadcastScalarToVector256(&b1), r1)
                                    : Avx.Add(r1, Avx.Multiply(va, Avx.BroadcastScalarToVector256(&b1)));
                                r2 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, Avx.BroadcastScalarToVector256(&b2), r2)
                                    : Avx.Add(r2, Avx.Multiply(va, Avx.BroadcastScalarToVector256(&b2)));
                                r3 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va, Avx.BroadcastScalarToVector256(&b3), r3)
                                    : Avx.Add(r3, Avx.Multiply(va, Avx.BroadcastScalarToVector256(&b3)));

                                Avx.Store(c0 + i, r0);
                                Avx.Store(c1 + i, r1);
                                Avx.Store(c2 + i, r2);
                                Avx.Store(c3 + i, r3);
                            }
                        }

                        for (; i < M; ++i)
                        {
                            var a = aCol[i];
                            c0[i] += a * b0;
                            c1[i] += a * b1;
                            c2[i] += a * b2;
                            c3[i] += a * b3;
                        }
                    }
                }

                for (; j < jEnd; ++j)
                {
                    var c = pC + j * ldc;
                    for (var k = 0; k < K; ++k)
                    {
                        var aCol = pA + k * lda;
                        var b = *(pB + k * ldb + j);
                        var i = 0;
                        if (Avx.IsSupported)
                        {
                            var vb = Avx.BroadcastScalarToVector256(&b);
                            for (; i + 7 < M; i += 8)
                            {
                                var va0 = Avx.LoadVector256(aCol + i + 0);
                                var va1 = Avx.LoadVector256(aCol + i + 4);

                                var vc0 = Avx.LoadVector256(c + i + 0);
                                var vc1 = Avx.LoadVector256(c + i + 4);
                                vc0 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va0, vb, vc0)
                                    : Avx.Add(vc0, Avx.Multiply(va0, vb));
                                vc1 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va1, vb, vc1)
                                    : Avx.Add(vc1, Avx.Multiply(va1, vb));
                                Avx.Store(c + i + 0, vc0);
                                Avx.Store(c + i + 4, vc1);
                            }

                            for (; i + 3 < M; i += 4)
                            {
                                var va = Avx.LoadVector256(aCol + i);
                                var vc = Avx.LoadVector256(c + i);
                                vc = Fma.IsSupported ? Fma.MultiplyAdd(va, vb, vc) : Avx.Add(vc, Avx.Multiply(va, vb));
                                Avx.Store(c + i, vc);
                            }
                        }

                        for (; i < M; ++i) c[i] += aCol[i] * b;
                    }
                }
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Gemm_At_Bt_Unblocked_N(Matrix A, Matrix B, Matrix C, int jStart, int jEnd)
    {
        var M = A.ColumnCount;
        var K = A.RowCount;
        var N = B.RowCount;
        var lda = A.RowCount;
        var ldb = B.RowCount;
        var ldc = C.RowCount;

        unsafe
        {
            fixed (double* pA = A._data, pB = B._data, pC = C._data)
            {
                for (var j = jStart; j < jEnd; ++j)
                {
                    var cCol = pC + j * ldc;

                    var i = 0;
                    for (; i + 3 < M; i += 4)
                    {
                        Vector256<double> acc0 = default, acc1 = default, acc2 = default, acc3 = default;

                        var k = 0;
                        if (Avx.IsSupported)
                            for (; k + 3 < K; k += 4)
                            {
                                var va_i = Vector256.Create(pA[(i + 0) * lda + k], pA[(i + 1) * lda + k],
                                    pA[(i + 2) * lda + k], pA[(i + 3) * lda + k]);
                                var va_i1 = Vector256.Create(pA[(i + 0) * lda + k + 1], pA[(i + 1) * lda + k + 1],
                                    pA[(i + 2) * lda + k + 1], pA[(i + 3) * lda + k + 1]);
                                var va_i2 = Vector256.Create(pA[(i + 0) * lda + k + 2], pA[(i + 1) * lda + k + 2],
                                    pA[(i + 2) * lda + k + 2], pA[(i + 3) * lda + k + 2]);
                                var va_i3 = Vector256.Create(pA[(i + 0) * lda + k + 3], pA[(i + 1) * lda + k + 3],
                                    pA[(i + 2) * lda + k + 3], pA[(i + 3) * lda + k + 3]);

                                var bkj = pB[k * ldb + j];
                                var vb0 = Avx.BroadcastScalarToVector256(&bkj);

                                var b_k1j = pB[(k + 1) * ldb + j];
                                var vb1 = Avx.BroadcastScalarToVector256(&b_k1j);

                                var b_k2j = pB[(k + 2) * ldb + j];
                                var vb2 = Avx.BroadcastScalarToVector256(&b_k2j);

                                var b_k3j = pB[(k + 3) * ldb + j];
                                var vb3 = Avx.BroadcastScalarToVector256(&b_k3j);

                                acc0 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va_i, vb0, acc0)
                                    : Avx.Add(acc0, Avx.Multiply(va_i, vb0));
                                acc1 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va_i1, vb1, acc1)
                                    : Avx.Add(acc1, Avx.Multiply(va_i1, vb1));
                                acc2 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va_i2, vb2, acc2)
                                    : Avx.Add(acc2, Avx.Multiply(va_i2, vb2));
                                acc3 = Fma.IsSupported
                                    ? Fma.MultiplyAdd(va_i3, vb3, acc3)
                                    : Avx.Add(acc3, Avx.Multiply(va_i3, vb3));
                            }

                        if (Avx.IsSupported)
                        {
                            acc0 = Avx.Add(acc0, acc1);
                            acc0 = Avx.Add(acc0, acc2);
                            acc0 = Avx.Add(acc0, acc3);

                            acc0.StoreUnsafe(ref cCol[i]);
                        }

                        for (; k < K; ++k)
                        {
                            var bkj = pB[k * ldb + j];
                            cCol[i + 0] += pA[(i + 0) * lda + k] * bkj;
                            cCol[i + 1] += pA[(i + 1) * lda + k] * bkj;
                            cCol[i + 2] += pA[(i + 2) * lda + k] * bkj;
                            cCol[i + 3] += pA[(i + 3) * lda + k] * bkj;
                        }
                    }

                    for (; i < M; ++i)
                    {
                        double sum = 0;
                        for (var k = 0; k < K; ++k) sum += pA[i * lda + k] * pB[k * ldb + j];
                        cCol[i] = sum;
                    }
                }
            }
        }
    }

    #endregion

    #region Element-wise Operations

    public Matrix Hadamard(Matrix other)
    {
        ArgumentNullException.ThrowIfNull(other);

        if (RowCount != other.RowCount || ColumnCount != other.ColumnCount)
            throw new ArgumentException("Matrices must have the same dimensions");

        var result = new Matrix(RowCount, ColumnCount);

        if (Avx512F.IsSupported && _data.Length >= 8)
        {
            unsafe
            {
                fixed (double* pThis = _data, pOther = other._data, pRes = result._data)
                {
                    var i = 0;
                    for (; i <= _data.Length - 8; i += 8)
                    {
                        var vThis = Avx512F.LoadVector512(pThis + i);
                        var vOther = Avx512F.LoadVector512(pOther + i);
                        Avx512F.Store(pRes + i, Avx512F.Multiply(vThis, vOther));
                    }

                    for (; i < _data.Length; i++)
                        pRes[i] = pThis[i] * pOther[i];
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && _data.Length >= 4)
        {
            var i = 0;
            for (; i <= _data.Length - 4; i += 4)
            {
                var vThis = Vector256.LoadUnsafe(ref _data[i]);
                var vOther = Vector256.LoadUnsafe(ref other._data[i]);
                (vThis * vOther).StoreUnsafe(ref result._data[i]);
            }

            for (; i < _data.Length; i++)
                result._data[i] = _data[i] * other._data[i];
        }
        else
        {
            for (var i = 0; i < _data.Length; i++)
                result._data[i] = _data[i] * other._data[i];
        }

        return result;
    }

    public Matrix Apply(Func<double, double> function)
    {
        ArgumentNullException.ThrowIfNull(function);

        var result = new Matrix(RowCount, ColumnCount);
        for (var i = 0; i < _data.Length; i++)
            result._data[i] = function(_data[i]);

        return result;
    }

    public Matrix Map(Func<int, int, double, double> function)
    {
        ArgumentNullException.ThrowIfNull(function);

        var result = new Matrix(RowCount, ColumnCount);

        for (var j = 0; j < ColumnCount; j++)
        for (var i = 0; i < RowCount; i++)
        {
            var index = j * RowCount + i;
            result._data[index] = function(i, j, _data[index]);
        }

        return result;
    }

    #endregion

    #region Matrix Properties

    public Matrix Transpose()
    {
        var result = new Matrix(ColumnCount, RowCount);

        unsafe
        {
            fixed (double* pSrc = _data, pDst = result._data)
            {
                const int BLOCK = 8;

                for (var ii = 0; ii < RowCount; ii += BLOCK)
                for (var jj = 0; jj < ColumnCount; jj += BLOCK)
                {
                    var iEnd = Math.Min(ii + BLOCK, RowCount);
                    var jEnd = Math.Min(jj + BLOCK, ColumnCount);

                    for (var i = ii; i < iEnd; i++)
                    for (var j = jj; j < jEnd; j++)
                        pDst[i * result.RowCount + j] = pSrc[j * RowCount + i];
                }
            }
        }

        return result;
    }

    /// <summary>
    ///     Transposes this matrix in-place. Only works for square matrices.
    ///     More efficient than Transpose() when the result can overwrite the input.
    /// </summary>
    public void TransposeInPlace()
    {
        if (!IsSquare)
            throw new InvalidOperationException("In-place transpose requires a square matrix");

        for (var i = 0; i < RowCount; i++)
        for (var j = i + 1; j < ColumnCount; j++)
        {
            var idx_ij = j * RowCount + i;
            var idx_ji = i * RowCount + j;
            (_data[idx_ij], _data[idx_ji]) = (_data[idx_ji], _data[idx_ij]);
        }
    }

    /// <summary>
    ///     Computes the trace of the matrix (sum of diagonal elements).
    /// </summary>
    /// <remarks>
    ///     The implementation uses a scalar loop which is efficient for the strided access pattern
    ///     of the diagonal in column-major storage.
    /// </remarks>
    public double Trace()
    {
        if (!IsSquare)
            throw new InvalidOperationException("Trace requires a square matrix");

        double sum = 0;
        var n = RowCount;

        for (var i = 0; i < n; i++)
            sum += _data[i * n + i];

        return sum;
    }


    public bool IsSymmetric(double tolerance = DefaultTolerance)
    {
        if (!IsSquare) return false;

        for (var i = 0; i < RowCount; i++)
        for (var j = i + 1; j < RowCount; j++)
        {
            var diff = Math.Abs(this[i, j] - this[j, i]);
            if (diff > tolerance) return false;
        }

        return true;
    }

    public bool IsDiagonal(double tolerance = DefaultTolerance)
    {
        for (var i = 0; i < RowCount; i++)
        for (var j = 0; j < ColumnCount; j++)
            if (i != j && Math.Abs(this[i, j]) > tolerance)
                return false;
        return true;
    }

    public bool IsUpperTriangular(double tolerance = DefaultTolerance)
    {
        for (var i = 1; i < RowCount; i++)
        for (var j = 0; j < Math.Min(i, ColumnCount); j++)
            if (Math.Abs(this[i, j]) > tolerance)
                return false;
        return true;
    }

    public bool IsLowerTriangular(double tolerance = DefaultTolerance)
    {
        for (var i = 0; i < RowCount - 1; i++)
        for (var j = i + 1; j < ColumnCount; j++)
            if (Math.Abs(this[i, j]) > tolerance)
                return false;
        return true;
    }

    public Matrix GetSymmetricPart()
    {
        if (!IsSquare)
            throw new InvalidOperationException("Symmetric part requires a square matrix");

        var result = new Matrix(RowCount, ColumnCount);

        for (var i = 0; i < RowCount; i++)
        for (var j = 0; j <= i; j++)
        {
            var value = (this[i, j] + this[j, i]) / 2.0;
            result[i, j] = value;
            result[j, i] = value;
        }

        return result;
    }

    public Matrix GetSkewSymmetricPart()
    {
        if (!IsSquare)
            throw new InvalidOperationException("Skew-symmetric part requires a square matrix");

        var result = new Matrix(RowCount, ColumnCount);

        for (var i = 0; i < RowCount; i++)
        for (var j = 0; j < i; j++)
        {
            var value = (this[i, j] - this[j, i]) / 2.0;
            result[i, j] = value;
            result[j, i] = -value;
        }

        return result;
    }

    public double FrobeniusNorm()
    {
        double sum = 0;

        if (Avx512F.IsSupported && _data.Length >= 8)
        {
            unsafe
            {
                fixed (double* pData = _data)
                {
                    var vSum = Vector512<double>.Zero;
                    var i = 0;

                    for (; i <= _data.Length - 8; i += 8)
                    {
                        var v = Avx512F.LoadVector512(pData + i);
                        vSum = Avx512F.FusedMultiplyAdd(v, v, vSum);
                    }

                    sum = Vector512.Sum(vSum);

                    for (; i < _data.Length; i++)
                        sum += _data[i] * _data[i];
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && _data.Length >= 4)
        {
            var vSum = Vector256<double>.Zero;
            var i = 0;

            for (; i <= _data.Length - 4; i += 4)
            {
                var v = Vector256.LoadUnsafe(ref _data[i]);
                vSum = Fma.MultiplyAdd(v, v, vSum);
            }

            sum = Vector256.Sum(vSum);

            for (; i < _data.Length; i++)
                sum += _data[i] * _data[i];
        }
        else
        {
            for (var i = 0; i < _data.Length; i++)
                sum += _data[i] * _data[i];
        }

        return Math.Sqrt(sum);
    }

    public double InfinityNorm()
    {
        var useStack = RowCount <= 1024;
        double[]? pooledRowSums = null;

        var rowSums = useStack
            ? stackalloc double[RowCount]
            : (pooledRowSums = ArrayPool<double>.Shared.Rent(RowCount)).AsSpan(0, RowCount);

        try
        {
            rowSums.Clear();

            for (var j = 0; j < ColumnCount; j++)
            {
                var colOffset = j * RowCount;
                for (var i = 0; i < RowCount; i++)
                    rowSums[i] += Math.Abs(_data[colOffset + i]);
            }

            double max = 0;
            for (var i = 0; i < RowCount; i++)
                if (rowSums[i] > max)
                    max = rowSums[i];

            return max;
        }
        finally
        {
            if (pooledRowSums != null)
                ArrayPool<double>.Shared.Return(pooledRowSums);
        }
    }

    public double OneNorm()
    {
        double max = 0;
        for (var j = 0; j < ColumnCount; j++)
        {
            double colSum = 0;
            for (var i = 0; i < RowCount; i++)
                colSum += Math.Abs(this[i, j]);
            if (colSum > max) max = colSum;
        }

        return max;
    }

    #endregion

    #region Determinants

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public double Determinant()
    {
        if (!IsSquare)
            throw new InvalidOperationException("Determinant requires a square matrix");

        return RowCount switch
        {
            1 => _data[0],
            2 => Det2x2(),
            3 => Det3x3(),
            _ => ComputeLU().Determinant
        };
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private double Det2x2()
    {
        return this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private double Det3x3()
    {
        double a = this[0, 0], b = this[0, 1], c = this[0, 2];
        double d = this[1, 0], e = this[1, 1], f = this[1, 2];
        double g = this[2, 0], h = this[2, 1], i = this[2, 2];

        return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    }

    #endregion

    #region Inverse Operations

    public Matrix Inverse()
    {
        if (!IsSquare)
            throw new InvalidOperationException("Only square matrices have inverses");

        var n = RowCount;

        return n switch
        {
            1 => Inverse1x1(),
            2 => Inverse2x2(),
            3 => Inverse3x3(),
            _ => InverseGeneral()
        };
    }

    private Matrix Inverse1x1()
    {
        var tol = MachinePrecision * OneNorm();
        if (Math.Abs(this[0, 0]) < tol)
            throw new InvalidOperationException("Matrix is singular");

        return new Matrix(1, 1) { [0, 0] = 1.0 / this[0, 0] };
    }

    private Matrix Inverse2x2()
    {
        var tol = MachinePrecision * OneNorm();
        var det = Det2x2();

        if (Math.Abs(det) < tol)
            throw new InvalidOperationException("Matrix is singular");

        var invDet = 1.0 / det;

        return new Matrix(2, 2)
        {
            [0, 0] = this[1, 1] * invDet,
            [0, 1] = -this[0, 1] * invDet,
            [1, 0] = -this[1, 0] * invDet,
            [1, 1] = this[0, 0] * invDet
        };
    }

    private Matrix Inverse3x3()
    {
        var tol = MachinePrecision * OneNorm();
        var det = Det3x3();

        if (Math.Abs(det) < tol)
            throw new InvalidOperationException("Matrix is singular");

        var invDet = 1.0 / det;
        var result = new Matrix(3, 3);


        result[0, 0] = (this[1, 1] * this[2, 2] - this[1, 2] * this[2, 1]) * invDet;
        result[0, 1] = (this[0, 2] * this[2, 1] - this[0, 1] * this[2, 2]) * invDet;
        result[0, 2] = (this[0, 1] * this[1, 2] - this[0, 2] * this[1, 1]) * invDet;

        result[1, 0] = (this[1, 2] * this[2, 0] - this[1, 0] * this[2, 2]) * invDet;
        result[1, 1] = (this[0, 0] * this[2, 2] - this[0, 2] * this[2, 0]) * invDet;
        result[1, 2] = (this[0, 2] * this[1, 0] - this[0, 0] * this[1, 2]) * invDet;

        result[2, 0] = (this[1, 0] * this[2, 1] - this[1, 1] * this[2, 0]) * invDet;
        result[2, 1] = (this[0, 1] * this[2, 0] - this[0, 0] * this[2, 1]) * invDet;
        result[2, 2] = (this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0]) * invDet;

        return result;
    }


    private Matrix InverseGeneral()
    {
        var lu = ComputeLU();
        var n = RowCount;
        var result = new Matrix(n, n);

        for (var j = 0; j < n; j++) // Iterate through columns of the identity matrix
        {
            var e = UnitVector(n, j);
            var x = lu.Solve(e);

            for (var i = 0; i < n; i++)
                result[i, j] = x._data[i];
        }

        return result;
    }

    public Matrix MatrixPower(int n)
    {
        if (!IsSquare)
            throw new InvalidOperationException("Matrix power requires a square matrix");

        if (n < 0)
            return Inverse().MatrixPower(-n);

        if (n == 0)
            return Identity(RowCount);

        if (n == 1)
            return Clone();

        var result = Identity(RowCount);
        var base_ = Clone();

        while (n > 0)
        {
            if ((n & 1) == 1)
                result = result * base_;

            base_ = base_ * base_;
            n >>= 1;
        }

        return result;
    }

    #endregion

    #region Decompositions

    /// <summary>
    ///     Computes the LU decomposition with Rook pivoting for optimal numerical stability.
    /// </summary>
    /// <returns>
    ///     An <see cref="LUDecomposition" /> object containing the L and U matrices, pivot permutations,
    ///     and methods for solving linear systems with iterative refinement and computing determinants.
    /// </returns>
    /// <exception cref="InvalidOperationException">
    ///     Thrown when the matrix is not square, is singular (determinant ≈ 0), or when Rook pivoting
    ///     fails to converge (rare, usually indicates a severely ill-conditioned matrix).
    /// </exception>
    /// <remarks>
    ///     <para>
    ///         The LU decomposition factors matrix A into A = PLUQ, where P and Q are permutation matrices,
    ///         L is lower triangular with unit diagonal, and U is upper triangular.
    ///     </para>
    ///     <para>
    ///         <strong>Time Complexity:</strong> O(n³) - approximately 2× slower than partial pivoting
    ///         but provides superior numerical stability.
    ///     </para>
    ///     <para><strong>Space Complexity:</strong> O(n²)</para>
    ///     <para>
    ///         <strong>Numerical Stability:</strong>
    ///         Rook pivoting provides optimal numerical stability for ill-conditioned or nearly-singular
    ///         matrices by searching for the maximum element in both rows and columns, significantly
    ///         reducing error propagation. Recommended for matrices arising from discretizations with
    ///         poor mesh quality or unknown conditioning.
    ///     </para>
    ///     <para>
    ///         <strong>Example Usage:</strong>
    ///     </para>
    ///     <code>
    /// // Compute LU decomposition
    /// var lu = matrix.ComputeLU();
    /// var x = lu.Solve(b);
    /// 
    /// // Compute determinant
    /// var det = lu.Determinant;
    /// </code>
    /// </remarks>
    public LUDecomposition ComputeLU()
    {
        if (!IsSquare)
            throw new InvalidOperationException("LU decomposition requires a square matrix");

        var n = RowCount;
        var lu = Clone();
        var rowPivots = new int[n];
        var colPivots = new int[n];
        var rowSwaps = 0;
        var colSwaps = 0;

        var scale = OneNorm();

        for (var i = 0; i < n; i++)
        {
            rowPivots[i] = i;
            colPivots[i] = i;
        }

        for (var k = 0; k < n - 1; k++)
        {
            int pivotRow = k, pivotCol = k;
            double maxVal = 0;

            for (var i = k; i < n; i++)
            for (var j = k; j < n; j++)
            {
                var val = Math.Abs(lu[i, j]);
                if (val > maxVal)
                {
                    maxVal = val;
                    pivotRow = i;
                    pivotCol = j;
                }
            }

            const int maxIter = MaxRookPivotingIterations;
            var converged = false;
            for (var iter = 0; iter < maxIter; iter++)
            {
                var newPivotCol = k;
                double rowMax = 0;
                for (var j = k; j < n; j++)
                {
                    var val = Math.Abs(lu[pivotRow, j]);
                    if (val > rowMax)
                    {
                        rowMax = val;
                        newPivotCol = j;
                    }
                }

                var newPivotRow = k;
                double colMax = 0;
                for (var i = k; i < n; i++)
                {
                    var val = Math.Abs(lu[i, newPivotCol]);
                    if (val > colMax)
                    {
                        colMax = val;
                        newPivotRow = i;
                    }
                }

                if (newPivotRow == pivotRow && newPivotCol == pivotCol)
                {
                    converged = true;
                    break;
                }

                pivotRow = newPivotRow;
                pivotCol = newPivotCol;
                maxVal = Math.Abs(lu[pivotRow, pivotCol]);
            }

            if (!converged)
                throw new InvalidOperationException(
                    $"Rook pivoting failed to converge after {maxIter} iterations at elimination step {k}. " +
                    "Matrix may be ill-conditioned.");

            if (maxVal < MachinePrecision * scale)
                throw new InvalidOperationException("Matrix is singular");

            if (pivotRow != k)
            {
                (rowPivots[k], rowPivots[pivotRow]) = (rowPivots[pivotRow], rowPivots[k]);
                lu.SwapRows(k, pivotRow);
                rowSwaps++;
            }

            if (pivotCol != k)
            {
                (colPivots[k], colPivots[pivotCol]) = (colPivots[pivotCol], colPivots[k]);
                lu.SwapColumns(k, pivotCol);
                colSwaps++;
            }

            for (var i = k + 1; i < n; i++)
            {
                lu[i, k] /= lu[k, k];

                for (var j = k + 1; j < n; j++)
                    lu[i, j] -= lu[i, k] * lu[k, j];
            }
        }

        var L = Identity(n);
        var U = Zeros(n, n);

        for (var i = 0; i < n; i++)
        {
            for (var j = 0; j < i; j++)
                L[i, j] = lu[i, j];

            for (var j = i; j < n; j++)
                U[i, j] = lu[i, j];
        }

        return new LUDecomposition(L, U, rowPivots, colPivots, rowSwaps, colSwaps, this);
    }

    /// <summary>
    ///     Computes the QR decomposition using Householder reflections.
    /// </summary>
    /// <returns>
    ///     A <see cref="QRDecomposition" /> object containing the orthogonal matrix Q and upper triangular matrix R.
    /// </returns>
    /// <remarks>
    ///     <para>
    ///         The QR decomposition factorizes this m×n matrix A into A = QR where:
    ///         - Q is an m×m orthogonal matrix (Q'Q = I)
    ///         - R is an m×n upper triangular matrix
    ///     </para>
    ///     <para>
    ///         <b>Time Complexity:</b> O(mn²) where m = rows, n = columns.
    ///     </para>
    ///     <para>
    ///         <b>Numerical Stability:</b> Householder QR is backward stable with excellent
    ///         numerical properties. Orthogonality of Q is maintained to machine precision for
    ///         well-conditioned matrices.
    ///     </para>
    ///     <para>
    ///         <b>Use Cases:</b> Solving least-squares problems (Ax=b where A is m×n, m≥n),
    ///         computing matrix rank, orthogonalization.
    ///     </para>
    /// </remarks>
    public QRDecomposition ComputeQR()
    {
        var m = RowCount;
        var n = ColumnCount;
        var qr = Clone();
        var rDiag = new double[n];

        for (var k = 0; k < n; k++)
        {
            double nrm = 0;
            for (var i = k; i < m; i++)
                nrm = double.Hypot(nrm, qr[i, k]);

            if (nrm != 0.0)
            {
                if (qr[k, k] < 0) nrm = -nrm;
                for (var i = k; i < m; i++)
                    qr[i, k] /= nrm;
                qr[k, k] += 1.0;

                for (var j = k + 1; j < n; j++)
                {
                    double s = 0;
                    for (var i = k; i < m; i++)
                        s += qr[i, k] * qr[i, j];
                    s = -s / qr[k, k];
                    for (var i = k; i < m; i++)
                        qr[i, j] += s * qr[i, k];
                }
            }

            rDiag[k] = -nrm;
        }

        var frobeniusNorm = FrobeniusNorm();
        return new QRDecomposition(qr, rDiag, frobeniusNorm);
    }

    /// <summary>
    ///     Computes the eigenvalue decomposition of a symmetric matrix.
    /// </summary>
    /// <returns>
    ///     An <see cref="EigenDecomposition" /> object containing the eigenvalues and eigenvectors.
    /// </returns>
    /// <exception cref="InvalidOperationException">
    ///     Thrown when the matrix is not square or not symmetric.
    /// </exception>
    /// <remarks>
    ///     <para>
    ///         <strong>IMPORTANT LIMITATION:</strong> This implementation currently supports only
    ///         <strong>symmetric matrices</strong>. For a symmetric matrix A, the eigenvalue decomposition
    ///         is A = V Λ V', where V is orthogonal (V' V = I) and Λ is diagonal with real eigenvalues.
    ///     </para>
    ///     <para>
    ///         <strong>Future Work:</strong> Support for general (non-symmetric) matrices using the QR algorithm
    ///         with complex eigenvalues will be added in a future version.
    ///     </para>
    ///     <para>
    ///         <strong>Algorithm:</strong>
    ///     </para>
    ///     <list type="number">
    ///         <item>
    ///             <description>Householder reduction to tridiagonal form (O(n³))</description>
    ///         </item>
    ///         <item>
    ///             <description>
    ///                 QR iteration with Wilkinson shifts on the tridiagonal matrix (O(n² k) where k = number of
    ///                 iterations)
    ///             </description>
    ///         </item>
    ///     </list>
    ///     <para><strong>Time Complexity:</strong> O(n³) total</para>
    ///     <para><strong>Space Complexity:</strong> O(n²)</para>
    ///     <para>
    ///         <strong>Numerical Stability:</strong>
    ///     </para>
    ///     <para>
    ///         The combined Householder + QR algorithm is backward stable for symmetric matrices.
    ///         Orthogonality of eigenvectors is maintained to machine precision.
    ///     </para>
    ///     <para>
    ///         <strong>Convergence:</strong>
    ///     </para>
    ///     <para>
    ///         Maximum iterations: 30n. For well-separated eigenvalues, convergence is typically
    ///         quadratic. For clustered eigenvalues, convergence may be slower but still reliable.
    ///     </para>
    ///     <para>
    ///         <strong>Example Usage:</strong>
    ///     </para>
    ///     <code>
    /// // Create symmetric matrix (C# 13 collection expression)
    /// double[][] data = [
    ///     [4, 1, 0],
    ///     [1, 3, 1],
    ///     [0, 1, 2]
    /// ];
    /// var A = new Matrix(data);
    /// 
    /// // Compute eigenvalues and eigenvectors
    /// var eigen = A.ComputeEigenvalues();
    /// var values = eigen.Eigenvalues;  // Real eigenvalues
    /// var vectors = eigen.Eigenvectors; // Columns are eigenvectors
    /// 
    /// // Verify: A*v_i = λ_i*v_i for each eigenpair
    /// </code>
    ///     <para>
    ///         <strong>See Also:</strong>
    ///     </para>
    ///     <list type="bullet">
    ///         <item>
    ///             <description><see cref="ComputeSVD" /> for singular value decomposition (works for any matrix)</description>
    ///         </item>
    ///         <item>
    ///             <description><see cref="IsSymmetric" /> to check if a matrix is symmetric</description>
    ///         </item>
    ///     </list>
    /// </remarks>
    public EigenDecomposition ComputeEigenvalues()
    {
        if (!IsSquare)
            throw new InvalidOperationException("Eigenvalue decomposition requires a square matrix");

        if (!IsSymmetric())
            throw new InvalidOperationException("Current implementation requires symmetric matrix");

        var n = RowCount;
        var A = Clone();
        var V = Identity(n);

        for (var k = 0; k < n - 2; k++)
        {
            double alpha = 0;
            for (var i = k + 1; i < n; i++)
                alpha += A[i, k] * A[i, k];
            alpha = Math.Sqrt(alpha);

            if (alpha > MachinePrecision)
            {
                if (A[k + 1, k] < 0) alpha = -alpha;
                A[k + 1, k] += alpha;

                var beta = 1.0 / (alpha * A[k + 1, k]);

                for (var j = k + 1; j < n; j++)
                {
                    double dot = 0;
                    for (var i = k + 1; i < n; i++)
                        dot += A[i, k] * A[i, j];
                    dot *= beta;

                    for (var i = k + 1; i < n; i++)
                        A[i, j] -= dot * A[i, k];
                }

                for (var i = 0; i < n; i++)
                {
                    double dot = 0;
                    for (var j = k + 1; j < n; j++)
                        dot += A[i, j] * A[j, k];
                    dot *= beta;

                    for (var j = k + 1; j < n; j++)
                        A[i, j] -= dot * A[j, k];
                }

                for (var i = 0; i < n; i++)
                {
                    double dot = 0;
                    for (var j = k + 1; j < n; j++)
                        dot += V[i, j] * A[j, k];
                    dot *= beta;

                    for (var j = k + 1; j < n; j++)
                        V[i, j] -= dot * A[j, k];
                }

                A[k + 1, k] = -alpha;
                for (var i = k + 2; i < n; i++)
                    A[i, k] = 0;
            }
        }

        var maxIter = 30 * n;
        for (var iter = 0; iter < maxIter; iter++)
        {
            var converged = true;
            for (var i = 0; i < n - 1; i++)
                if (Math.Abs(A[i + 1, i]) > MachinePrecision * (Math.Abs(A[i, i]) + Math.Abs(A[i + 1, i + 1])))
                {
                    converged = false;
                    break;
                }

            if (converged) break;

            var a = A[n - 2, n - 2];
            var b = A[n - 1, n - 2];
            var c = A[n - 1, n - 1];
            var d = (a - c) / 2;
            var sign = d >= 0 ? 1.0 : -1.0;
            var shift = c - b * b / (d + sign * Math.Sqrt(d * d + b * b));

            for (var k = 0; k < n - 1; k++)
            {
                double c1, s;

                if (k == 0)
                {
                    var p = A[0, 0] - shift;
                    var q = A[1, 0];
                    var r = Math.Sqrt(p * p + q * q);
                    c1 = p / r;
                    s = q / r;
                }
                else
                {
                    var p = A[k, k - 1];
                    var q = A[k + 1, k - 1];
                    var r = Math.Sqrt(p * p + q * q);
                    c1 = p / r;
                    s = q / r;
                    A[k, k - 1] = r;
                    A[k + 1, k - 1] = 0;
                }

                for (var j = k; j < n; j++)
                {
                    var temp = c1 * A[k, j] + s * A[k + 1, j];
                    A[k + 1, j] = -s * A[k, j] + c1 * A[k + 1, j];
                    A[k, j] = temp;
                }

                for (var i = 0; i <= Math.Min(k + 2, n - 1); i++)
                {
                    var temp = c1 * A[i, k] + s * A[i, k + 1];
                    A[i, k + 1] = -s * A[i, k] + c1 * A[i, k + 1];
                    A[i, k] = temp;
                }

                for (var i = 0; i < n; i++)
                {
                    var temp = c1 * V[i, k] + s * V[i, k + 1];
                    V[i, k + 1] = -s * V[i, k] + c1 * V[i, k + 1];
                    V[i, k] = temp;
                }
            }
        }

        var eigenvalues = new double[n];
        for (var i = 0; i < n; i++)
            eigenvalues[i] = A[i, i];

        var indices = new int[n];
        for (var i = 0; i < n; i++)
            indices[i] = i;

        Array.Sort(indices, (i, j) => eigenvalues[j].CompareTo(eigenvalues[i]));

        var sortedEigenvalues = new double[n];
        var sortedEigenvectors = new Matrix(n, n);

        for (var i = 0; i < n; i++)
        {
            sortedEigenvalues[i] = eigenvalues[indices[i]];
            for (var j = 0; j < n; j++)
                sortedEigenvectors[j, i] = V[j, indices[i]];
        }

        return new EigenDecomposition(sortedEigenvalues, sortedEigenvectors);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DotProductKernel(int m, int j, int k, ref double cjj, ref double ckk, ref double cjk)
    {
        var data = _data;
        var jOffset = j * m;
        var kOffset = k * m;

        if (Vector256.IsHardwareAccelerated && m >= 4)
        {
            var vCjj = Vector256<double>.Zero;
            var vCkk = Vector256<double>.Zero;
            var vCjk = Vector256<double>.Zero;

            var i = 0;
            for (; i <= m - 4; i += 4)
            {
                var vBj = Vector256.LoadUnsafe(ref data[jOffset + i]);
                var vBk = Vector256.LoadUnsafe(ref data[kOffset + i]);
                vCjj = Fma.MultiplyAdd(vBj, vBj, vCjj);
                vCkk = Fma.MultiplyAdd(vBk, vBk, vCkk);
                vCjk = Fma.MultiplyAdd(vBj, vBk, vCjk);
            }

            cjj = Vector256.Sum(vCjj);
            ckk = Vector256.Sum(vCkk);
            cjk = Vector256.Sum(vCjk);

            for (; i < m; i++)
            {
                var b_ij = data[jOffset + i];
                var b_ik = data[kOffset + i];
                cjj += b_ij * b_ij;
                ckk += b_ik * b_ik;
                cjk += b_ij * b_ik;
            }
        }
        else
        {
            for (var i = 0; i < m; i++)
            {
                var b_ij = data[jOffset + i];
                var b_ik = data[kOffset + i];
                cjj += b_ij * b_ij;
                ckk += b_ik * b_ik;
                cjk += b_ij * b_ik;
            }
        }
    }

    /// <summary>
    ///     Computes the Singular Value Decomposition using the Jacobi algorithm.
    /// </summary>
    /// <returns>
    ///     An <see cref="SVDDecomposition" /> object containing U, S (diagonal singular values), and V.
    /// </returns>
    /// <remarks>
    ///     <para>
    ///         The SVD factorizes this m×n matrix A into A = U·S·V' where:
    ///         - U is an m×m orthogonal matrix (left singular vectors)
    ///         - S is an m×n diagonal matrix with non-negative singular values
    ///         - V is an n×n orthogonal matrix (right singular vectors)
    ///     </para>
    ///     <para>
    ///         <b>Time Complexity:</b> O(mn² + n³) where m = rows, n = columns.
    ///         For square matrices: O(n³).
    ///     </para>
    ///     <para>
    ///         <b>Algorithm:</b> Uses the one-sided Jacobi method with block processing for cache efficiency
    ///         and SIMD vectorization. Iterates until off-diagonal elements are below tolerance or
    ///         maximum iterations (100) is reached.
    ///     </para>
    ///     <para>
    ///         <b>Numerical Stability:</b> Jacobi SVD is highly accurate and preserves orthogonality well.
    ///         However, convergence can be slow for ill-conditioned matrices. A warning is printed if
    ///         the algorithm doesn't fully converge.
    ///     </para>
    ///     <para>
    ///         <b>Use Cases:</b> Computing matrix rank, condition number, pseudoinverse, least-squares
    ///         solutions, principal component analysis, matrix approximation.
    ///     </para>
    /// </remarks>
    /// <summary>
    ///     Computes the Singular Value Decomposition using the Jacobi algorithm.
    /// </summary>
    /// <returns>
    ///     An <see cref="SVDDecomposition" /> object containing U, S (diagonal singular values), and V.
    /// </returns>
    /// <remarks>
    ///     <para>
    ///         The SVD factorizes this m×n matrix A into A = U·S·V' where:
    ///         - U is an m×m orthogonal matrix (left singular vectors)
    ///         - S is an m×n diagonal matrix with non-negative singular values
    ///         - V is an n×n orthogonal matrix (right singular vectors)
    ///     </para>
    ///     <para>
    ///         <b>Time Complexity:</b> O(mn² + n³) where m = rows, n = columns.
    ///         For square matrices: O(n³).
    ///     </para>
    ///     <para>
    ///         <b>Algorithm:</b> Uses the one-sided Jacobi method with block processing for cache efficiency
    ///         and SIMD vectorization. Iterates until off-diagonal elements are below tolerance or
    ///         maximum iterations (100) is reached.
    ///     </para>
    ///     <para>
    ///         <b>Numerical Stability:</b> Jacobi SVD is highly accurate and preserves orthogonality well.
    ///         However, convergence can be slow for ill-conditioned matrices. A warning is printed if
    ///         the algorithm doesn't fully converge.
    ///     </para>
    ///     <para>
    ///         <b>Use Cases:</b> Computing matrix rank, condition number, pseudoinverse, least-squares
    ///         solutions, principal component analysis, matrix approximation.
    ///     </para>
    /// </remarks>
    /// <summary>
    ///     Computes the Singular Value Decomposition (SVD): A = U * S * V'.
    /// </summary>
    /// <param name="randomSeed">
    ///     Optional seed for random number generation (used for rank-deficient matrices).
    ///     If null, uses non-deterministic System.Random.Shared.
    ///     Provide a seed for reproducible results in testing and debugging.
    /// </param>
    /// <returns>The SVD decomposition with matrices U, S, and V.</returns>
    /// <remarks>
    ///     For rank-deficient matrices, the null space is filled using randomized
    ///     Gram-Schmidt orthogonalization. Providing a seed ensures deterministic results.
    /// </remarks>
    [SkipLocalsInit]
    public SVDDecomposition ComputeSVD(int? randomSeed = null)
    {
        var m = RowCount;
        var n = ColumnCount;
        var minDim = Math.Min(m, n);

        var size = m * n;
        var isPooled = size >= StackAllocThreshold;
        var tempDataArray = isPooled ? ArrayPool<double>.Shared.Rent(size) : new double[size];

        var tempDataMatrix = new Matrix(m, n, tempDataArray);

        try
        {
            Array.Copy(_data, 0, tempDataArray, 0, size);

            var V = Identity(n);

            var singularValues = new double[n];
            var tolerance = MachinePrecision * FrobeniusNorm();

            var lastOffDiagNorm = double.MaxValue;
            var stallCount = 0;

            for (var sweep = 0; sweep < MaxIterations; sweep++)
            {
                var converged = true;
                double offDiagNorm = 0;

                const int BLOCK = 4;

                for (var jBlock = 0; jBlock < n - 1; jBlock += BLOCK)
                {
                    var jEnd = Math.Min(jBlock + BLOCK, n - 1);

                    for (var kBlock = jBlock; kBlock < n; kBlock += BLOCK)
                    {
                        var kEnd = Math.Min(kBlock + BLOCK, n);

                        for (var j = jBlock; j < jEnd; j++)
                        for (var k = Math.Max(j + 1, kBlock); k < kEnd; k++)
                        {
                            double cjj = 0, ckk = 0, cjk = 0;

                            tempDataMatrix.DotProductKernel(m, j, k, ref cjj, ref ckk, ref cjk);

                            offDiagNorm += cjk * cjk;

                            if (Math.Abs(cjk) <= tolerance)
                                continue;

                            converged = false;

                            var tau = (ckk - cjj) / (2.0 * cjk);
                            var t = Math.Sign(tau) / (Math.Abs(tau) + Math.Sqrt(1.0 + tau * tau));
                            var c = 1.0 / Math.Sqrt(1.0 + t * t);
                            var s = c * t;

                            var jOffset = j * m;
                            var kOffset = k * m;

                            if (Vector256.IsHardwareAccelerated && m >= 4)
                            {
                                var vC = Vector256.Create(c);
                                var vS = Vector256.Create(s);
                                var vNegS = Vector256.Create(-s);

                                var i = 0;
                                for (; i <= m - 4; i += 4)
                                {
                                    var vBj = Vector256.LoadUnsafe(ref tempDataArray[jOffset + i]);
                                    var vBk = Vector256.LoadUnsafe(ref tempDataArray[kOffset + i]);

                                    var newBj = Fma.MultiplyAdd(vC, vBj, Vector256.Multiply(vNegS, vBk));
                                    var newBk = Fma.MultiplyAdd(vS, vBj, Vector256.Multiply(vC, vBk));

                                    newBj.StoreUnsafe(ref tempDataArray[jOffset + i]);
                                    newBk.StoreUnsafe(ref tempDataArray[kOffset + i]);
                                }

                                for (; i < m; i++)
                                {
                                    var b_ij = tempDataArray[jOffset + i];
                                    var b_ik = tempDataArray[kOffset + i];
                                    tempDataArray[jOffset + i] = c * b_ij - s * b_ik;
                                    tempDataArray[kOffset + i] = s * b_ij + c * b_ik;
                                }
                            }
                            else
                            {
                                for (var i = 0; i < m; i++)
                                {
                                    var b_ij = tempDataArray[jOffset + i];
                                    var b_ik = tempDataArray[kOffset + i];
                                    tempDataArray[jOffset + i] = c * b_ij - s * b_ik;
                                    tempDataArray[kOffset + i] = s * b_ij + c * b_ik;
                                }
                            }

                            for (var i = 0; i < n; i++)
                            {
                                var v_ij = V[i, j];
                                var v_ik = V[i, k];
                                V[i, j] = c * v_ij - s * v_ik;
                                V[i, k] = s * v_ij + c * v_ik;
                            }
                        }
                    }
                }

                if (converged)
                    break;

                offDiagNorm = Math.Sqrt(offDiagNorm);
                var relativeDelta = Math.Abs(offDiagNorm - lastOffDiagNorm) / (lastOffDiagNorm + 1e-100);

                if (relativeDelta < MachinePrecision)
                {
                    stallCount++;
                    if (stallCount > 3)
                        break;
                }
                else
                {
                    stallCount = 0;
                }

                lastOffDiagNorm = offDiagNorm;

                if (sweep == MaxIterations - 1)
                    Debug.WriteLine($"Warning: SVD did not fully converge after {MaxIterations} iterations. " +
                                    $"Relative off-diagonal norm: {offDiagNorm}");
            }

            for (var j = 0; j < n; j++)
            {
                double norm = 0;
                for (var i = 0; i < m; i++)
                    norm += tempDataArray[j * m + i] * tempDataArray[j * m + i];
                singularValues[j] = Math.Sqrt(norm);
            }

            var indices = new int[n];
            for (var i = 0; i < n; i++)
                indices[i] = i;

            Array.Sort(indices, (i, j) => singularValues[j].CompareTo(singularValues[i]));

            var U = new Matrix(m, m);
            var S = Zeros(m, n);
            var sortedV = new Matrix(n, n);
            var sortedSingularValues = new double[minDim];

            var rank = 0;

            for (var k = 0; k < n; k++)
            {
                var sortedIndex = indices[k];
                var s_k = singularValues[sortedIndex];

                for (var i = 0; i < n; i++)
                    sortedV[i, k] = V[i, sortedIndex];

                if (k < minDim)
                {
                    sortedSingularValues[k] = s_k;
                    S[k, k] = s_k;
                }

                if (s_k > tolerance)
                {
                    for (var i = 0; i < m; i++)
                        U[i, k] = tempDataArray[sortedIndex * m + i] / s_k;
                    rank++;
                }
            }

            var rand = randomSeed.HasValue
                ? new Random(randomSeed.Value)
                : System.Random.Shared;
            var v = m <= 1024 ? stackalloc double[m] : new double[m];
            for (var k = rank; k < m; k++)
            {
                for (var i = 0; i < m; i++)
                    v[i] = rand.NextDouble() - 0.5;

                for (var j = 0; j < k; j++)
                {
                    double dot = 0;
                    for (var i = 0; i < m; i++)
                        dot += v[i] * U[i, j];
                    for (var i = 0; i < m; i++)
                        v[i] -= dot * U[i, j];
                }

                double norm = 0;
                for (var i = 0; i < m; i++)
                    norm += v[i] * v[i];
                norm = Math.Sqrt(norm);

                if (norm > MachinePrecision)
                {
                    for (var i = 0; i < m; i++)
                        U[i, k] = v[i] / norm;
                }
                else
                {
                    for (var i = 0; i < m; i++) v[i] = i == k ? 1.0 : 0.0;

                    for (var j = 0; j < k; j++)
                    {
                        double dot = 0;
                        for (var i = 0; i < m; i++)
                            dot += v[i] * U[i, j];
                        for (var i = 0; i < m; i++)
                            v[i] -= dot * U[i, j];
                    }

                    norm = 0;
                    for (var i = 0; i < m; i++) norm += v[i] * v[i];
                    norm = Math.Sqrt(norm > MachinePrecision ? norm : 1.0);

                    for (var i = 0; i < m; i++)
                        U[i, k] = v[i] / norm;
                }
            }

            return new SVDDecomposition(U, S, sortedV);
        }
        finally
        {
            if (isPooled)
                ArrayPool<double>.Shared.Return(tempDataArray, true);
        }
    }

    public int Rank(double tolerance = DefaultTolerance)
    {
        var svd = ComputeSVD();
        return svd.Rank(tolerance);
    }

    public double ConditionNumber()
    {
        var svd = ComputeSVD();
        return svd.ConditionNumber();
    }

    public double SpectralNorm()
    {
        var svd = ComputeSVD();
        return svd.SingularValues.Length > 0 ? svd.SingularValues[0] : 0;
    }

    #endregion

    #region Row and Column Operations

    public Vector GetRow(int row)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(row);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(row, RowCount);

        var data = new double[ColumnCount];
        for (var j = 0; j < ColumnCount; j++)
            data[j] = this[row, j];

        return new Vector(data);
    }

    public void SetRow(int row, Vector values)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(row);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(row, RowCount);
        ArgumentNullException.ThrowIfNull(values);

        if (values.Length != ColumnCount)
            throw new ArgumentException("Vector length must match column count");

        for (var j = 0; j < ColumnCount; j++)
            this[row, j] = values[j];
    }

    public Vector GetColumn(int column)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(column);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(column, ColumnCount);

        var data = new double[RowCount];
        Array.Copy(_data, column * RowCount, data, 0, RowCount);

        return new Vector(data);
    }

    public void SetColumn(int column, Vector values)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(column);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(column, ColumnCount);
        ArgumentNullException.ThrowIfNull(values);

        if (values.Length != RowCount)
            throw new ArgumentException("Vector length must match row count");

        Array.Copy(values._data, 0, _data, column * RowCount, RowCount);
    }

    public Matrix GetSubMatrix(int startRow, int endRow, int startCol, int endCol)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(startRow);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(endRow, RowCount);
        ArgumentOutOfRangeException.ThrowIfNegative(startCol);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(endCol, ColumnCount);

        if (startRow >= endRow || startCol >= endCol)
            throw new ArgumentException("Invalid submatrix range");

        var rows = endRow - startRow;
        var cols = endCol - startCol;
        var result = new Matrix(rows, cols);

        for (var j = 0; j < cols; j++)
        {
            var srcColOffset = (startCol + j) * RowCount + startRow;
            var dstColOffset = j * rows;

            Array.Copy(_data, srcColOffset, result._data, dstColOffset, rows);
        }

        return result;
    }

    public void SetSubMatrix(int startRow, int startCol, Matrix subMatrix)
    {
        ArgumentNullException.ThrowIfNull(subMatrix);
        ArgumentOutOfRangeException.ThrowIfNegative(startRow);
        ArgumentOutOfRangeException.ThrowIfNegative(startCol);

        var endRow = startRow + subMatrix.RowCount;
        var endCol = startCol + subMatrix.ColumnCount;

        if (endRow > RowCount || endCol > ColumnCount)
            throw new ArgumentException("Submatrix doesn't fit in matrix");

        for (var j = 0; j < subMatrix.ColumnCount; j++)
        {
            var srcColOffset = j * subMatrix.RowCount;
            var dstColOffset = (startCol + j) * RowCount + startRow;

            Array.Copy(subMatrix._data, srcColOffset, _data, dstColOffset, subMatrix.RowCount);
        }
    }

    public void SwapRows(int row1, int row2)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(row1);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(row1, RowCount);
        ArgumentOutOfRangeException.ThrowIfNegative(row2);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(row2, RowCount);

        if (row1 == row2) return;

        for (var j = 0; j < ColumnCount; j++)
            (this[row1, j], this[row2, j]) = (this[row2, j], this[row1, j]);
    }

    public void SwapColumns(int col1, int col2)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(col1);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(col1, ColumnCount);
        ArgumentOutOfRangeException.ThrowIfNegative(col2);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(col2, ColumnCount);

        if (col1 == col2) return;

        if (RowCount <= StackAllocThreshold)
        {
            Span<double> buffer = stackalloc double[RowCount];
            var col1Span = new Span<double>(_data, col1 * RowCount, RowCount);
            var col2Span = new Span<double>(_data, col2 * RowCount, RowCount);

            col1Span.CopyTo(buffer);
            col2Span.CopyTo(col1Span);
            buffer.CopyTo(col2Span);
            return;
        }

        var pooledBuffer = MatrixPool.Rent(RowCount);
        try
        {
            var buffer = pooledBuffer.AsSpan(0, RowCount);
            var col1Span = new Span<double>(_data, col1 * RowCount, RowCount);
            var col2Span = new Span<double>(_data, col2 * RowCount, RowCount);

            col1Span.CopyTo(buffer);
            col2Span.CopyTo(col1Span);
            buffer.CopyTo(col2Span);
        }
        finally
        {
            MatrixPool.Return(pooledBuffer);
        }
    }

    public Matrix GetNullSpace(double tolerance = DefaultTolerance)
    {
        var svd = ComputeSVD();
        var n = ColumnCount;
        var rank = svd.Rank(tolerance);
        var nullity = n - rank;

        if (nullity == 0)
            return new Matrix(n, 0);

        var nullSpace = new Matrix(n, nullity);
        var j = 0;

        for (var i = rank; i < n; i++)
        {
            for (var k = 0; k < n; k++)
                nullSpace[k, j] = svd.V[k, i];
            j++;
        }

        return nullSpace;
    }

    public Matrix GetRowSpace(double tolerance = DefaultTolerance)
    {
        var svd = ComputeSVD();
        var rank = svd.Rank(tolerance);
        var rowSpace = new Matrix(ColumnCount, rank);

        for (var i = 0; i < rank; i++)
        for (var j = 0; j < ColumnCount; j++)
            rowSpace[j, i] = svd.V[j, i];

        return rowSpace;
    }

    public Matrix GetImageSpace(double tolerance = DefaultTolerance)
    {
        var svd = ComputeSVD();
        var rank = svd.Rank(tolerance);
        var imageSpace = new Matrix(RowCount, rank);

        for (var i = 0; i < rank; i++)
        for (var j = 0; j < RowCount; j++)
            imageSpace[j, i] = svd.U[j, i];

        return imageSpace;
    }

    public Matrix HorizontalConcat(Matrix right)
    {
        ArgumentNullException.ThrowIfNull(right);

        if (RowCount != right.RowCount)
            throw new ArgumentException("Matrices must have same number of rows");

        var result = new Matrix(RowCount, ColumnCount + right.ColumnCount);

        for (var j = 0; j < ColumnCount; j++)
            Array.Copy(_data, j * RowCount, result._data, j * RowCount, RowCount);


        for (var j = 0; j < right.ColumnCount; j++)
            Array.Copy(right._data, j * right.RowCount, result._data, (ColumnCount + j) * result.RowCount,
                right.RowCount);

        return result;
    }

    public Matrix VerticalConcat(Matrix bottom)
    {
        ArgumentNullException.ThrowIfNull(bottom);

        if (ColumnCount != bottom.ColumnCount)
            throw new ArgumentException("Matrices must have same number of columns");

        var result = new Matrix(RowCount + bottom.RowCount, ColumnCount);

        for (var j = 0; j < ColumnCount; j++)
        {
            Array.Copy(_data, j * RowCount, result._data, j * result.RowCount, RowCount);
            Array.Copy(bottom._data, j * bottom.RowCount, result._data, j * result.RowCount + RowCount,
                bottom.RowCount);
        }

        return result;
    }

    public Matrix Reshape(int newRows, int newCols)
    {
        if (newRows * newCols != ElementCount)
            throw new ArgumentException("New dimensions must have same total elements");

        var result = new Matrix(newRows, newCols);
        Array.Copy(_data, result._data, ElementCount);

        return result;
    }

    #endregion

    #region In-Place and Special Operations

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void ScaleInPlace(double scalar)
    {
        var data = _data;

        if (Avx512F.IsSupported && data.Length >= 8)
        {
            unsafe
            {
                fixed (double* pData = data)
                {
                    var vScalar = Vector512.Create(scalar);
                    var i = 0;

                    for (; i <= data.Length - 8; i += 8)
                    {
                        var v = Avx512F.LoadVector512(pData + i);
                        Avx512F.Store(pData + i, Avx512F.Multiply(v, vScalar));
                    }

                    for (; i < data.Length; i++)
                        pData[i] *= scalar;
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && data.Length >= 4)
        {
            var vScalar = Vector256.Create(scalar);
            var i = 0;

            for (; i <= data.Length - 4; i += 4)
            {
                var v = Vector256.LoadUnsafe(ref data[i]);
                (v * vScalar).StoreUnsafe(ref data[i]);
            }

            for (; i < data.Length; i++)
                data[i] *= scalar;
        }
        else
        {
            for (var i = 0; i < data.Length; i++)
                data[i] *= scalar;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void AddInPlace(Matrix other)
    {
        if (RowCount != other.RowCount || ColumnCount != other.ColumnCount)
            throw new ArgumentException("Matrix dimensions must match");

        var data = _data;
        var otherData = other._data;

        if (Avx512F.IsSupported && data.Length >= 8)
        {
            unsafe
            {
                fixed (double* pData = data, pOther = otherData)
                {
                    var i = 0;

                    for (; i <= data.Length - 8; i += 8)
                    {
                        var v1 = Avx512F.LoadVector512(pData + i);
                        var v2 = Avx512F.LoadVector512(pOther + i);
                        Avx512F.Store(pData + i, Avx512F.Add(v1, v2));
                    }

                    for (; i < data.Length; i++)
                        pData[i] += pOther[i];
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && data.Length >= 4)
        {
            var i = 0;

            for (; i <= data.Length - 4; i += 4)
            {
                var v1 = Vector256.LoadUnsafe(ref data[i]);
                var v2 = Vector256.LoadUnsafe(ref otherData[i]);
                (v1 + v2).StoreUnsafe(ref data[i]);
            }

            for (; i < data.Length; i++)
                data[i] += otherData[i];
        }
        else
        {
            for (var i = 0; i < data.Length; i++)
                data[i] += otherData[i];
        }
    }

    public Matrix KroneckerProduct(Matrix B)
    {
        var m = RowCount * B.RowCount;
        var n = ColumnCount * B.ColumnCount;
        var result = new Matrix(m, n);

        for (var i = 0; i < RowCount; i++)
        for (var j = 0; j < ColumnCount; j++)
        {
            var aij = this[i, j];
            var rowOffset = i * B.RowCount;
            var colOffset = j * B.ColumnCount;

            for (var k = 0; k < B.RowCount; k++)
            for (var l = 0; l < B.ColumnCount; l++)
                result[rowOffset + k, colOffset + l] = aij * B[k, l];
        }

        return result;
    }

    public Matrix PseudoInverse(double tolerance = 1e-10)
    {
        var svd = ComputeSVD();
        var m = RowCount;
        var n = ColumnCount;

        var sPlus = Zeros(n, m);
        var maxSv = svd.SingularValues.Length > 0 ? svd.SingularValues[0] : 0;
        var threshold = tolerance * maxSv;

        for (var i = 0; i < Math.Min(m, n); i++)
            if (svd.SingularValues[i] > threshold)
                sPlus[i, i] = 1.0 / svd.SingularValues[i];

        return svd.V * sPlus * svd.U.Transpose();
    }

    #endregion

    #region Solving Linear Systems

    /// <summary>
    ///     Solves a linear system A*x = b for a square matrix A.
    /// </summary>
    /// <param name="b">The right-hand side column vector (<see cref="Vector" />).</param>
    /// <returns>The solution vector x (<see cref="Vector" />).</returns>
    public Vector Solve(Vector b)
    {
        ArgumentNullException.ThrowIfNull(b);

        if (!IsSquare)
            throw new InvalidOperationException("Matrix must be square");
        if (b.Length != RowCount)
            throw new ArgumentException("Right-hand side must be a vector matching matrix rows");

        var lu = ComputeLU();
        return lu.Solve(b);
    }

    public Matrix SolveMultiple(Matrix B)
    {
        ArgumentNullException.ThrowIfNull(B);

        if (!IsSquare)
            throw new InvalidOperationException("Matrix must be square");

        if (B.RowCount != RowCount)
            throw new ArgumentException("Right-hand side rows must match matrix rows");

        var lu = ComputeLU();
        var X = new Matrix(RowCount, B.ColumnCount);

        for (var j = 0; j < B.ColumnCount; j++)
        {
            var b = B.GetColumn(j);
            var x = lu.Solve(b);
            X.SetColumn(j, x);
        }

        return X;
    }

    /// <summary>
    ///     Solves a linear system A*x = b using least squares.
    /// </summary>
    /// <param name="b">The right-hand side column vector (<see cref="Vector" />).</param>
    /// <returns>The least-squares solution vector x (<see cref="Vector" />).</returns>
    public Vector SolveLeastSquares(Vector b)
    {
        ArgumentNullException.ThrowIfNull(b);

        if (b.Length != RowCount)
            throw new ArgumentException("Right-hand side must be a column vector matching matrix rows");

        if (RowCount >= ColumnCount)
        {
            var qr = ComputeQR();
            return qr.Solve(b);
        }

        var svd = ComputeSVD();
        return svd.Solve(b);
    }

    #endregion

    #region Statistical Operations

    public double[] ColumnMeans()
    {
        var means = new double[ColumnCount];

        for (var j = 0; j < ColumnCount; j++)
        {
            double sum = 0;
            var colOffset = j * RowCount;

            if (Vector256.IsHardwareAccelerated && RowCount >= 4)
            {
                var vSum = Vector256<double>.Zero;
                var i = 0;

                for (; i <= RowCount - 4; i += 4)
                {
                    var v = Vector256.LoadUnsafe(ref _data[colOffset + i]);
                    vSum = Vector256.Add(vSum, v);
                }

                sum = Vector256.Sum(vSum);

                for (; i < RowCount; i++)
                    sum += _data[colOffset + i];
            }
            else
            {
                for (var i = 0; i < RowCount; i++)
                    sum += _data[colOffset + i];
            }

            means[j] = sum / RowCount;
        }

        return means;
    }

    public double[] RowMeans()
    {
        var means = new double[RowCount];

        if (Vector256.IsHardwareAccelerated && ColumnCount >= 4)
            for (var i = 0; i < RowCount; i++)
            {
                var vSum = Vector256<double>.Zero;
                var j = 0;

                for (; j <= ColumnCount - 4; j += 4)
                {
                    var v = Vector256.Create(
                        this[i, j],
                        this[i, j + 1],
                        this[i, j + 2],
                        this[i, j + 3]
                    );
                    vSum = Vector256.Add(vSum, v);
                }

                var sum = Vector256.Sum(vSum);
                for (; j < ColumnCount; j++)
                    sum += this[i, j];

                means[i] = sum / ColumnCount;
            }
        else
            for (var i = 0; i < RowCount; i++)
            {
                double sum = 0;
                for (var j = 0; j < ColumnCount; j++)
                    sum += this[i, j];
                means[i] = sum / ColumnCount;
            }

        return means;
    }

    public Matrix Covariance()
    {
        var m = ColumnCount;
        var n = RowCount;

        if (n < 2)
            throw new InvalidOperationException(
                "Need at least 2 samples for unbiased covariance. " +
                "For single-sample or biased covariance, consider computing manually with n instead of (n-1).");

        var means = ColumnMeans();

        var centered = new Matrix(n, m);
        for (var j = 0; j < m; j++)
        {
            var mean = means[j];
            var srcOffset = j * n;
            var dstOffset = j * n;

            if (Vector256.IsHardwareAccelerated && n >= 4)
            {
                var vMean = Vector256.Create(mean);
                var i = 0;
                for (; i <= n - 4; i += 4)
                {
                    var v = Vector256.LoadUnsafe(ref _data[srcOffset + i]);
                    (v - vMean).StoreUnsafe(ref centered._data[dstOffset + i]);
                }

                for (; i < n; i++)
                    centered._data[dstOffset + i] = _data[srcOffset + i] - mean;
            }
            else
            {
                for (var i = 0; i < n; i++)
                    centered._data[dstOffset + i] = _data[srcOffset + i] - mean;
            }
        }

        var cov = new Matrix(m, m);
        for (var j1 = 0; j1 < m; j1++)
        {
            var col1Offset = j1 * n;

            for (var j2 = j1; j2 < m; j2++)
            {
                var col2Offset = j2 * n;
                double sum = 0;

                if (Vector256.IsHardwareAccelerated && n >= 4)
                {
                    var vSum = Vector256<double>.Zero;
                    var i = 0;
                    for (; i <= n - 4; i += 4)
                    {
                        var v1 = Vector256.LoadUnsafe(ref centered._data[col1Offset + i]);
                        var v2 = Vector256.LoadUnsafe(ref centered._data[col2Offset + i]);
                        vSum = Fma.MultiplyAdd(v1, v2, vSum);
                    }

                    sum = Vector256.Sum(vSum);
                    for (; i < n; i++)
                        sum += centered._data[col1Offset + i] * centered._data[col2Offset + i];
                }
                else
                {
                    for (var i = 0; i < n; i++)
                        sum += centered._data[col1Offset + i] * centered._data[col2Offset + i];
                }

                cov[j1, j2] = sum / (n - 1);
                if (j1 != j2)
                    cov[j2, j1] = cov[j1, j2];
            }
        }

        return cov;
    }

    public Matrix Correlation()
    {
        var cov = Covariance();
        var m = cov.RowCount;
        var corr = new Matrix(m, m);

        var stdDevs = new double[m];
        for (var i = 0; i < m; i++)
            stdDevs[i] = Math.Sqrt(cov[i, i]);

        for (var i = 0; i < m; i++)
        for (var j = 0; j < m; j++)
            corr[i, j] = cov[i, j] / (stdDevs[i] * stdDevs[j]);

        return corr;
    }

    #endregion

    #region Numerical Algorithms

    public static Matrix MatrixSquareRoot(Matrix A,
        int maxIterations = MaxMatrixSquareRootIterations,
        double tolerance = MatrixSquareRootTolerance)
    {
        if (!A.IsSquare)
            throw new InvalidOperationException("Matrix must be square");

        var Y = A.Clone();
        var Z = Identity(A.RowCount);

        for (var k = 0; k < maxIterations; k++)
        {
            var Y_old = Y.Clone();

            var Y_next = (Y + Z.Inverse()) * 0.5;
            var Z_next = (Z + Y.Inverse()) * 0.5;

            Y = Y_next;
            Z = Z_next;

            var diff = (Y - Y_old).FrobeniusNorm();
            if (diff < tolerance)
                break;
        }

        return Y;
    }

    public static Matrix MatrixExponential(Matrix A, int order = 6)
    {
        if (!A.IsSquare)
            throw new InvalidOperationException("Matrix exponential requires square matrix");

        var n = A.RowCount;
        var norm = A.InfinityNorm();

        var s = 0;
        var scale = 1.0;
        while (norm * scale > 1)
        {
            s++;
            scale *= 0.5;
        }

        var scaledA = A * scale;

        var I = Identity(n);
        var X = I.Clone();
        var cX = I.Clone();
        var N = I.Clone();
        var D = I.Clone();
        var c = 1.0;

        for (var k = 1; k <= order; k++)
        {
            c = c * (order - k + 1) / (k * (2 * order - k + 1));
            X = X * scaledA;
            cX = X * c;
            N.AddInPlace(cX);
            D.AddInPlace(k % 2 == 0 ? cX : -cX);
        }

        var result = D.Inverse() * N;

        for (var k = 0; k < s; k++)
            result = result * result;

        return result;
    }

    public (Matrix U, Matrix P) ComputePolarDecomposition()
    {
        if (!IsSquare)
            throw new InvalidOperationException(
                "Polar decomposition is only defined for square matrices. " +
                $"Current matrix dimensions: {RowCount}×{ColumnCount}");

        var svd = ComputeSVD();
        var U = svd.U * svd.V.Transpose();
        var P = svd.V * svd.S * svd.V.Transpose();
        return (U, P);
    }

    #endregion
}

/// <summary>
///     Helper class for throwing common exceptions efficiently.
/// </summary>
/// <remarks>
///     Methods are marked with <c>NoInlining</c> to keep hot paths small and improve code cache utilization.
/// </remarks>
/// <summary>
///     A high-performance, sealed class for 1D numerical vector operations.
/// </summary>
/// <remarks>
///     The Vector class is designed for optimal performance, utilizing SIMD vectorization and
///     seamless integration with the Matrix class for linear system solving.
/// </remarks>
public sealed class Vector : IEquatable<Vector>, IFormattable, ICloneable
{
    [JsonInclude] internal readonly double[] _data;
    [JsonInclude] public int Length { get; }

    #region Constants (Shared with Matrix class)

    private const double DefaultTolerance = 1e-10;
    private const double MachinePrecision = 2.220446049250313e-16;
    private const double ZeroVectorThreshold = 1e-10;
    private const double OverflowSafetyThreshold = 1e154;

    #endregion

    #region Constructors and Conversion

    /// <summary>
    ///     Creates a new, zero-initialized vector.
    /// </summary>
    public Vector(int length)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(length);
        Length = length;
        _data = GC.AllocateArray<double>(length, length > 1000);
    }

    /// <summary>
    ///     Creates a vector from a source array.
    /// </summary>
    public Vector(ReadOnlySpan<double> source)
    {
        ArgumentOutOfRangeException.ThrowIfZero(source.Length, nameof(source));
        Length = source.Length;
        _data = GC.AllocateArray<double>(Length);
        source.CopyTo(_data);
    }

    internal Vector(double[] data)
    {
        Length = data.Length;
        _data = data;
    }


    /// <summary>
    ///     Implicit conversion from Vector to Matrix (Nx1 column vector).
    ///     Clones the underlying data once and attaches it directly to a new Matrix.
    /// </summary>
    public static implicit operator Matrix(Vector v)
    {
        ArgumentNullException.ThrowIfNull(v);

        var data = (double[])v._data.Clone();
        return new Matrix(v.Length, 1, data);
    }

    /// <summary>
    ///     Explicit conversion from Matrix to Vector (requires 1xN or Nx1 Matrix).
    ///     Clones the underlying storage once and wraps it in a Vector.
    /// </summary>
    public static explicit operator Vector(Matrix m)
    {
        ArgumentNullException.ThrowIfNull(m);

        if (m.ColumnCount != 1 && m.RowCount != 1)
            throw new InvalidCastException("Matrix must be a column (Nx1) or row (1xN) vector to convert to Vector.");

        var data = (double[])m._data.Clone();
        return new Vector(data);
    }

    #endregion

    #region Indexers and Utilities

    public double this[int index]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if ((uint)index >= (uint)Length)
                ThrowHelper.ThrowArgumentOutOfRangeException(nameof(index));
            return _data[index];
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        set
        {
            if ((uint)index >= (uint)Length)
                ThrowHelper.ThrowArgumentOutOfRangeException(nameof(index));

            if (Matrix.EffectiveValidateFiniteValues && !double.IsFinite(value))
                throw new ArgumentException(
                    $"Non-finite value ({value}) detected at index {index}. " +
                    "Set Matrix.ValidateFiniteValues = false to disable this check.", nameof(value));

            _data[index] = value;
        }
    }

    public Vector Clone()
    {
        var cloneData = (double[])_data.Clone();
        return new Vector(cloneData);
    }

    object ICloneable.Clone()
    {
        return Clone();
    }

    public bool Equals(Vector? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        if (Length != other.Length) return false;

        for (var i = 0; i < _data.Length; i++)
            if (Math.Abs(_data[i] - other._data[i]) > DefaultTolerance)
                return false;
        return true;
    }

    public override bool Equals(object? obj)
    {
        return Equals(obj as Vector);
    }

    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.Add(Length);

        // Sample up to 20 elements distributed throughout the vector
        const int maxSamples = 20;
        var sampleCount = Math.Min(_data.Length, maxSamples);

        if (sampleCount <= maxSamples)
        {
            // Small vector: include all elements
            for (var i = 0; i < sampleCount; i++)
                hash.Add(_data[i]);
        }
        else
        {
            // Large vector: sample evenly distributed elements
            var stride = _data.Length / maxSamples;
            for (var i = 0; i < maxSamples; i++)
                hash.Add(_data[i * stride]);
        }

        return hash.ToHashCode();
    }

    public override string ToString()
    {
        return ToString("G", null);
    }

    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        format ??= "G";
        formatProvider ??= CultureInfo.CurrentCulture;

        var sb = new StringBuilder();
        sb.AppendLine($"Vector({Length}):");

        // Calculate maximum width needed for proper alignment
        var maxWidth = 10; // Minimum width
        var displayCount = Math.Min(Length, 10);

        for (var i = 0; i < displayCount; i++)
        {
            var str = this[i].ToString(format, formatProvider);
            maxWidth = Math.Max(maxWidth, str.Length);
        }

        // Cap at reasonable maximum
        maxWidth = Math.Min(maxWidth, 20);

        sb.Append("[ ");
        for (var i = 0; i < displayCount; i++)
        {
            if (i > 0) sb.Append("  ");
            sb.Append(this[i].ToString(format, formatProvider).PadLeft(maxWidth));
        }

        if (Length > 10) sb.Append(" ...");
        sb.AppendLine(" ]");

        return sb.ToString();
    }

    #endregion

    #region Factory Methods

    public static Vector Zeros(int length)
    {
        return new Vector(length);
    }

    public static Vector Ones(int length)
    {
        var result = new Vector(length);
        Array.Fill(result._data, 1.0);
        return result;
    }

    public static Vector Random(int length, int? seed = null)
    {
        var random = seed.HasValue ? new Random(seed.Value) : System.Random.Shared;
        var result = new Vector(length);
        for (var i = 0; i < result._data.Length; i++)
            result._data[i] = random.NextDouble();
        return result;
    }

    #endregion

    #region Arithmetic Operators (SIMD-Optimized)

    public static Vector operator +(Vector left, Vector right)
    {
        if (left.Length != right.Length)
            throw new ArgumentException("Vector lengths must match");

        var result = new Vector(left.Length);
        var length = left.Length;

        var i = 0;
        if (Avx512F.IsSupported && length >= 8)
            unsafe
            {
                fixed (double* pL = left._data, pR = right._data, pRes = result._data)
                {
                    for (; i <= length - 8; i += 8)
                    {
                        var vL = Avx512F.LoadVector512(pL + i);
                        var vR = Avx512F.LoadVector512(pR + i);
                        Avx512F.Store(pRes + i, Avx512F.Add(vL, vR));
                    }
                }
            }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
            for (; i <= length - 4; i += 4)
            {
                var vL = Vector256.LoadUnsafe(ref left._data[i]);
                var vR = Vector256.LoadUnsafe(ref right._data[i]);
                (vL + vR).StoreUnsafe(ref result._data[i]);
            }

        for (; i < length; i++)
            result._data[i] = left._data[i] + right._data[i];

        return result;
    }

    public static Vector operator -(Vector left, Vector right)
    {
        if (left.Length != right.Length)
            throw new ArgumentException("Vector lengths must match");

        var result = new Vector(left.Length);
        var length = left.Length;

        var i = 0;
        if (Avx512F.IsSupported && length >= 8)
            unsafe
            {
                fixed (double* pL = left._data, pR = right._data, pRes = result._data)
                {
                    for (; i <= length - 8; i += 8)
                    {
                        var vL = Avx512F.LoadVector512(pL + i);
                        var vR = Avx512F.LoadVector512(pR + i);
                        Avx512F.Store(pRes + i, Avx512F.Subtract(vL, vR));
                    }
                }
            }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
            for (; i <= length - 4; i += 4)
            {
                var vL = Vector256.LoadUnsafe(ref left._data[i]);
                var vR = Vector256.LoadUnsafe(ref right._data[i]);
                (vL - vR).StoreUnsafe(ref result._data[i]);
            }

        for (; i < length; i++)
            result._data[i] = left._data[i] - right._data[i];

        return result;
    }

    public static Vector operator *(Vector vector, double scalar)
    {
        var result = new Vector(vector.Length);
        var length = vector.Length;

        var i = 0;
        if (Avx512F.IsSupported && length >= 8)
        {
            unsafe
            {
                fixed (double* pV = vector._data, pRes = result._data)
                {
                    var vScalar = Vector512.Create(scalar);
                    for (; i <= length - 8; i += 8)
                    {
                        var v = Avx512F.LoadVector512(pV + i);
                        Avx512F.Store(pRes + i, Avx512F.Multiply(v, vScalar));
                    }
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
        {
            var vScalar = Vector256.Create(scalar);
            for (; i <= length - 4; i += 4)
            {
                var v = Vector256.LoadUnsafe(ref vector._data[i]);
                (v * vScalar).StoreUnsafe(ref result._data[i]);
            }
        }

        for (; i < length; i++)
            result._data[i] = vector._data[i] * scalar;

        return result;
    }

    public static Vector operator *(double scalar, Vector vector)
    {
        return vector * scalar;
    }

    public static Vector operator /(Vector vector, double scalar)
    {
        return vector * (1.0 / scalar);
    }

    /// <summary>
    ///     Computes the element-wise (Hadamard) division of two vectors (result[i] = left[i] / right[i]).
    /// </summary>
    public static Vector operator /(Vector left, Vector right)
    {
        return left.ElementwiseDivide(right);
    }

    #endregion

    #region Vector Operations (Dot Product, Norm)

    /// <summary>
    ///     Computes the dot product of two vectors.
    /// </summary>
    public double Dot(Vector other)
    {
        if (Length != other.Length)
            throw new ArgumentException("Vector lengths must match for dot product.");

        double sum = 0;
        var dataA = _data;
        var dataB = other._data;
        var length = Length;

        if (Avx512F.IsSupported && length >= 8)
        {
            unsafe
            {
                fixed (double* pA = dataA, pB = dataB)
                {
                    var vSum = Vector512<double>.Zero;
                    var i = 0;
                    for (; i <= length - 8; i += 8)
                    {
                        var vA = Avx512F.LoadVector512(pA + i);
                        var vB = Avx512F.LoadVector512(pB + i);
                        vSum = Avx512F.FusedMultiplyAdd(vA, vB, vSum);
                    }

                    sum = Vector512.Sum(vSum);
                    for (; i < length; i++)
                        sum += dataA[i] * dataB[i];
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
        {
            var vSum = Vector256<double>.Zero;
            var i = 0;
            for (; i <= length - 4; i += 4)
            {
                var vA = Vector256.LoadUnsafe(ref dataA[i]);
                var vB = Vector256.LoadUnsafe(ref dataB[i]);
                vSum = Fma.MultiplyAdd(vA, vB, vSum);
            }

            sum = Vector256.Sum(vSum);
            for (; i < length; i++)
                sum += dataA[i] * dataB[i];
        }
        else
        {
            for (var i = 0; i < length; i++)
                sum += dataA[i] * dataB[i];
        }

        return sum;
    }

    /// <summary>
    ///     Computes the L2-norm (Euclidean length) of this vector.
    /// </summary>
    public double Norm()
    {
        return Math.Sqrt(Dot(this));
    }

    /// <summary>
    ///     Returns a new vector with the same direction but a length of 1.
    /// </summary>
    public Vector Normalize()
    {
        var norm = Norm();
        if (norm < ZeroVectorThreshold)
            throw new InvalidOperationException(
                $"Cannot normalize near-zero vector (norm = {norm:E3}, threshold = {ZeroVectorThreshold:E3})");
        return this / norm;
    }

    /// <summary>
    ///     Computes the element-wise (Hadamard) product of two vectors.
    /// </summary>
    /// <param name="other">The other vector of the same length.</param>
    /// <returns>A new vector where result[i] = this[i] * other[i].</returns>
    public Vector ElementwiseMultiply(Vector other)
    {
        if (Length != other.Length)
            throw new ArgumentException("Vector lengths must match for element-wise multiplication.");

        var result = new Vector(Length);
        var length = Length;

        var i = 0;
        if (Avx512F.IsSupported && length >= 8)
            unsafe
            {
                fixed (double* pL = _data, pR = other._data, pRes = result._data)
                {
                    for (; i <= length - 8; i += 8)
                    {
                        var vL = Avx512F.LoadVector512(pL + i);
                        var vR = Avx512F.LoadVector512(pR + i);
                        Avx512F.Store(pRes + i, Avx512F.Multiply(vL, vR));
                    }
                }
            }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
            for (; i <= length - 4; i += 4)
            {
                var vL = Vector256.LoadUnsafe(ref _data[i]);
                var vR = Vector256.LoadUnsafe(ref other._data[i]);
                (vL * vR).StoreUnsafe(ref result._data[i]);
            }

        for (; i < length; i++)
            result._data[i] = _data[i] * other._data[i];

        return result;
    }

    /// <summary>
    ///     Computes the element-wise (Hadamard) division of two vectors.
    /// </summary>
    /// <param name="other">The other vector of the same length.</param>
    /// <returns>A new vector where result[i] = this[i] / other[i].</returns>
    public Vector ElementwiseDivide(Vector other)
    {
        if (Length != other.Length)
            throw new ArgumentException("Vector lengths must match for element-wise division.");

        var result = new Vector(Length);
        var length = Length;

        var i = 0;
        if (Avx512F.IsSupported && length >= 8)
            unsafe
            {
                fixed (double* pL = _data, pR = other._data, pRes = result._data)
                {
                    for (; i <= length - 8; i += 8)
                    {
                        var vL = Avx512F.LoadVector512(pL + i);
                        var vR = Avx512F.LoadVector512(pR + i);
                        Avx512F.Store(pRes + i, Avx512F.Divide(vL, vR));
                    }
                }
            }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
            for (; i <= length - 4; i += 4)
            {
                var vL = Vector256.LoadUnsafe(ref _data[i]);
                var vR = Vector256.LoadUnsafe(ref other._data[i]);
                (vL / vR).StoreUnsafe(ref result._data[i]);
            }

        for (; i < length; i++)
            result._data[i] = _data[i] / other._data[i];

        return result;
    }

    /// <summary>
    ///     Computes the outer product of two vectors, producing a matrix.
    /// </summary>
    /// <param name="other">The other vector.</param>
    /// <returns>A matrix M where M[i,j] = this[i] * other[j].</returns>
    /// <remarks>
    ///     The resulting matrix has dimensions (this.Length × other.Length).
    ///     This is also known as the tensor product or dyadic product.
    /// </remarks>
    public Matrix OuterProduct(Vector other)
    {
        ArgumentNullException.ThrowIfNull(other);

        var m = Length;
        var n = other.Length;

        var maxMagnitudeThis = 0.0;
        var maxMagnitudeOther = 0.0;

        for (var idx = 0; idx < m; idx++)
            maxMagnitudeThis = Math.Max(maxMagnitudeThis, Math.Abs(_data[idx]));

        for (var idx = 0; idx < n; idx++)
            maxMagnitudeOther = Math.Max(maxMagnitudeOther, Math.Abs(other._data[idx]));

        if (double.IsNaN(maxMagnitudeThis) || double.IsNaN(maxMagnitudeOther))
            throw new ArgumentException(
                "Input vectors contain NaN values. Outer product is undefined for NaN inputs.");

        if (double.IsInfinity(maxMagnitudeThis) || double.IsInfinity(maxMagnitudeOther))
            throw new ArgumentException(
                "Input vectors contain infinite values. Outer product would produce infinite results.");

        if (maxMagnitudeThis * maxMagnitudeOther > OverflowSafetyThreshold)
            throw new OverflowException(
                $"OuterProduct would overflow: maximum magnitudes {maxMagnitudeThis:E3} × {maxMagnitudeOther:E3} " +
                $"exceed safe threshold ({OverflowSafetyThreshold:E3}). " +
                "Consider scaling input vectors or using logarithmic representations.");

        var result = new Matrix(m, n);

        for (var j = 0; j < n; j++)
        {
            var other_j = other._data[j];
            var colOffset = j * m;

            var i = 0;
            if (Avx512F.IsSupported && m >= 8)
            {
                unsafe
                {
                    fixed (double* pThis = _data, pRes = &result._data[colOffset])
                    {
                        var vOther = Vector512.Create(other_j);
                        for (; i <= m - 8; i += 8)
                        {
                            var vThis = Avx512F.LoadVector512(pThis + i);
                            Avx512F.Store(pRes + i, Avx512F.Multiply(vThis, vOther));
                        }
                    }
                }
            }
            else if (Vector256.IsHardwareAccelerated && m >= 4)
            {
                var vOther = Vector256.Create(other_j);
                for (; i <= m - 4; i += 4)
                {
                    var vThis = Vector256.LoadUnsafe(ref _data[i]);
                    (vThis * vOther).StoreUnsafe(ref result._data[colOffset + i]);
                }
            }

            for (; i < m; i++)
                result._data[colOffset + i] = _data[i] * other_j;
        }

        return result;
    }

    /// <summary>
    ///     Computes the cross product of two 3D vectors.
    /// </summary>
    /// <param name="other">The other 3D vector.</param>
    /// <returns>The cross product vector perpendicular to both input vectors.</returns>
    /// <exception cref="ArgumentException">Thrown if either vector is not 3D.</exception>
    /// <remarks>
    ///     The cross product is only defined for 3D vectors.
    ///     Result = [a_y*b_z - a_z*b_y, a_z*b_x - a_x*b_z, a_x*b_y - a_y*b_x]
    /// </remarks>
    public Vector Cross(Vector other)
    {
        if (Length != 3 || other.Length != 3)
            throw new ArgumentException("Cross product is only defined for 3D vectors.");

        var a = _data;
        var b = other._data;

        return new Vector([
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ]);
    }

    /// <summary>
    ///     Computes the sum of all elements in the vector.
    /// </summary>
    public double Sum()
    {
        double sum = 0;
        var length = Length;

        var i = 0;
        if (Avx512F.IsSupported && length >= 8)
        {
            unsafe
            {
                fixed (double* pData = _data)
                {
                    var vSum = Vector512<double>.Zero;
                    for (; i <= length - 8; i += 8)
                    {
                        var v = Avx512F.LoadVector512(pData + i);
                        vSum = Avx512F.Add(vSum, v);
                    }

                    sum = Vector512.Sum(vSum);
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
        {
            var vSum = Vector256<double>.Zero;
            for (; i <= length - 4; i += 4)
            {
                var v = Vector256.LoadUnsafe(ref _data[i]);
                vSum += v;
            }

            sum = Vector256.Sum(vSum);
        }

        for (; i < length; i++)
            sum += _data[i];

        return sum;
    }

    /// <summary>
    ///     Computes the mean (average) of all elements in the vector.
    /// </summary>
    public double Mean()
    {
        return Sum() / Length;
    }

    /// <summary>
    ///     Computes the maximum element in the vector.
    /// </summary>
    public double Max()
    {
        if (Length == 0)
            throw new InvalidOperationException("Cannot compute max of empty vector");

        var max = _data[0];
        for (var i = 1; i < Length; i++)
            if (_data[i] > max)
                max = _data[i];

        return max;
    }

    /// <summary>
    ///     Computes the minimum element in the vector.
    /// </summary>
    public double Min()
    {
        if (Length == 0)
            throw new InvalidOperationException("Cannot compute min of empty vector");

        var min = _data[0];
        for (var i = 1; i < Length; i++)
            if (_data[i] < min)
                min = _data[i];

        return min;
    }

    /// <summary>
    ///     Finds the index of the maximum element.
    /// </summary>
    public int ArgMax()
    {
        if (Length == 0)
            throw new InvalidOperationException("Cannot compute argmax of empty vector");

        var maxIndex = 0;
        var maxValue = _data[0];

        for (var i = 1; i < Length; i++)
            if (_data[i] > maxValue)
            {
                maxValue = _data[i];
                maxIndex = i;
            }

        return maxIndex;
    }

    /// <summary>
    ///     Finds the index of the minimum element.
    /// </summary>
    public int ArgMin()
    {
        if (Length == 0)
            throw new InvalidOperationException("Cannot compute argmin of empty vector");

        var minIndex = 0;
        var minValue = _data[0];

        for (var i = 1; i < Length; i++)
            if (_data[i] < minValue)
            {
                minValue = _data[i];
                minIndex = i;
            }

        return minIndex;
    }

    /// <summary>
    ///     Computes the L1-norm (Manhattan distance) of the vector.
    /// </summary>
    public double Norm1()
    {
        double sum = 0;
        for (var i = 0; i < Length; i++)
            sum += Math.Abs(_data[i]);
        return sum;
    }

    /// <summary>
    ///     Computes the L-infinity norm (maximum absolute value) of the vector.
    /// </summary>
    public double NormInf()
    {
        double max = 0;
        for (var i = 0; i < Length; i++)
        {
            var abs = Math.Abs(_data[i]);
            if (abs > max)
                max = abs;
        }

        return max;
    }

    /// <summary>
    ///     Creates a new vector containing a subset of elements.
    /// </summary>
    /// <param name="startIndex">The starting index (inclusive).</param>
    /// <param name="length">The number of elements to extract.</param>
    public Vector Slice(int startIndex, int length)
    {
        if (startIndex < 0 || startIndex >= Length)
            throw new ArgumentOutOfRangeException(nameof(startIndex));
        if (length < 0 || startIndex + length > Length)
            throw new ArgumentOutOfRangeException(nameof(length));

        var result = new Vector(length);
        Array.Copy(_data, startIndex, result._data, 0, length);
        return result;
    }

    /// <summary>
    ///     Applies a function element-wise to create a new vector.
    /// </summary>
    public Vector Map(Func<double, double> func)
    {
        ArgumentNullException.ThrowIfNull(func);

        var result = new Vector(Length);
        for (var i = 0; i < Length; i++)
            result._data[i] = func(_data[i]);
        return result;
    }

    /// <summary>
    ///     Applies a function element-wise with index to create a new vector.
    /// </summary>
    public Vector Map(Func<double, int, double> func)
    {
        ArgumentNullException.ThrowIfNull(func);

        var result = new Vector(Length);
        for (var i = 0; i < Length; i++)
            result._data[i] = func(_data[i], i);
        return result;
    }

    /// <summary>
    ///     Creates a copy of this vector with all elements negated.
    /// </summary>
    public static Vector operator -(Vector vector)
    {
        var result = new Vector(vector.Length);
        var length = vector.Length;

        var i = 0;
        if (Avx512F.IsSupported && length >= 8)
        {
            unsafe
            {
                fixed (double* pV = vector._data, pRes = result._data)
                {
                    var zero = Vector512<double>.Zero;
                    for (; i <= length - 8; i += 8)
                    {
                        var v = Avx512F.LoadVector512(pV + i);
                        Avx512F.Store(pRes + i, Avx512F.Subtract(zero, v));
                    }
                }
            }
        }
        else if (Vector256.IsHardwareAccelerated && length >= 4)
        {
            var zero = Vector256<double>.Zero;
            for (; i <= length - 4; i += 4)
            {
                var v = Vector256.LoadUnsafe(ref vector._data[i]);
                (zero - v).StoreUnsafe(ref result._data[i]);
            }
        }

        for (; i < length; i++)
            result._data[i] = -vector._data[i];

        return result;
    }

    /// <summary>
    ///     Computes the angle between two vectors in radians.
    /// </summary>
    /// <param name="other">The other vector.</param>
    /// <returns>The angle in radians [0, π].</returns>
    public double AngleTo(Vector other)
    {
        if (Length != other.Length)
            throw new ArgumentException("Vector lengths must match to compute angle.");

        var dot = Dot(other);
        var norm1 = Norm();
        var norm2 = other.Norm();

        if (norm1 < 1e-10 || norm2 < 1e-10)
            throw new InvalidOperationException("Cannot compute angle with zero vector");

        var cosAngle = dot / (norm1 * norm2);
        cosAngle = Math.Max(-1.0, Math.Min(1.0, cosAngle));
        return Math.Acos(cosAngle);
    }

    /// <summary>
    ///     Computes the Euclidean distance to another vector.
    /// </summary>
    public double DistanceTo(Vector other)
    {
        if (Length != other.Length)
            throw new ArgumentException("Vector lengths must match to compute distance.");

        return (this - other).Norm();
    }

    /// <summary>
    ///     Projects this vector onto another vector.
    /// </summary>
    /// <param name="onto">The vector to project onto.</param>
    /// <returns>The projection of this vector onto the 'onto' vector.</returns>
    /// <remarks>
    ///     projection = (this · onto / |onto|²) * onto
    /// </remarks>
    public Vector ProjectOnto(Vector onto)
    {
        if (Length != onto.Length)
            throw new ArgumentException("Vector lengths must match for projection.");

        var normSquared = onto.Dot(onto);
        if (normSquared < 1e-10)
            throw new InvalidOperationException("Cannot project onto zero vector");

        var scale = Dot(onto) / normSquared;
        return onto * scale;
    }

    /// <summary>
    ///     Converts this vector to a standard double array.
    /// </summary>
    public double[] ToArray()
    {
        return (double[])_data.Clone();
    }

    /// <summary>
    ///     Returns a read-only span view of the vector data.
    /// </summary>
    public ReadOnlySpan<double> AsSpan()
    {
        return new ReadOnlySpan<double>(_data);
    }

    #endregion
}

#region Decomposition Classes

public sealed class LUDecomposition
{
    internal LUDecomposition(Matrix l, Matrix u, int[] rowPivots, int[] colPivots, int rowSwaps, int colSwaps,
        Matrix originalMatrix)
    {
        _l = l;
        _u = u;
        _rowPivots = rowPivots;
        _columnPivots = colPivots;
        _originalMatrix = originalMatrix;

        var n = u.RowCount;
        var logAbsDet = 0.0;
        var sign = (rowSwaps + colSwaps) % 2 == 0 ? 1 : -1;
        var isSingular = false;

        var matrixScale = originalMatrix.OneNorm();
        var singularityThreshold = MachinePrecision * matrixScale;

        if (singularityThreshold < DeterminantUnderflowThreshold)
            singularityThreshold = DeterminantUnderflowThreshold;

        for (var i = 0; i < n; i++)
        {
            var u_ii = u[i, i];
            var absU_ii = Math.Abs(u_ii);

            if (absU_ii < singularityThreshold)
            {
                isSingular = true;
                break;
            }

            logAbsDet += Math.Log(absU_ii);
            if (u_ii < 0)
                sign = -sign;
        }

        if (isSingular)
        {
            Determinant = 0.0;
            LogAbsDeterminant = double.NegativeInfinity;
            DeterminantSign = 0;
        }
        else
        {
            LogAbsDeterminant = logAbsDet;
            DeterminantSign = sign;

            if (logAbsDet > 700)
                Determinant = sign > 0 ? double.PositiveInfinity : double.NegativeInfinity;
            else if (logAbsDet < -700)
                Determinant = 0.0;
            else
                Determinant = sign * Math.Exp(logAbsDet);
        }
    }

    /// <summary>
    ///     Gets a copy of the lower triangular matrix L.
    /// </summary>
    /// <remarks>
    ///     Returns a clone to prevent external modification that could invalidate the decomposition.
    /// </remarks>
    public Matrix L => _l.Clone();

    /// <summary>
    ///     Gets a copy of the upper triangular matrix U.
    /// </summary>
    /// <remarks>
    ///     Returns a clone to prevent external modification that could invalidate the decomposition.
    /// </remarks>
    public Matrix U => _u.Clone();

    /// <summary>
    ///     Gets a copy of the row pivot array.
    /// </summary>
    public int[] RowPivots => (int[])_rowPivots.Clone();

    /// <summary>
    ///     Gets a copy of the column pivot array.
    /// </summary>
    public int[] ColumnPivots => (int[])_columnPivots.Clone();

    public double Determinant { get; }
    public double LogAbsDeterminant { get; }
    public int DeterminantSign { get; }

    /// <summary>
    ///     Solves A*x = b given the LU decomposition with iterative refinement for improved accuracy.
    /// </summary>
    /// <param name="b">The right-hand side column vector (Nx1 vector).</param>
    /// <returns>The solution vector x (Nx1 vector).</returns>
    /// <remarks>
    ///     Uses iterative refinement to improve the accuracy of the solution. The process:
    ///     1. Solves LU*x = b for initial solution
    ///     2. Computes residual r = b - A*x
    ///     3. Solves LU*dx = r for correction
    ///     4. Updates x = x + dx
    ///     5. Repeats until convergence (typically 1-3 iterations)
    /// </remarks>
    public Vector Solve(Vector b)
    {
        var n = _l.RowCount;
        if (b.Length != n)
            throw new ArgumentException("Right-hand side vector length must match matrix rows");

        var x = SolveWithoutRefinement(b);

        const int maxIterations = 3;
        const double convergenceTolerance = 1e-14;

        for (var iter = 0; iter < maxIterations; iter++)
        {
            var residual = new double[n];

            if (Vector256.IsHardwareAccelerated && n >= 4)
                for (var i = 0; i < n; i++)
                {
                    var vSum = Vector256<double>.Zero;
                    var j = 0;

                    for (; j <= n - 4; j += 4)
                    {
                        var vA = Vector256.Create(
                            _originalMatrix[i, j],
                            _originalMatrix[i, j + 1],
                            _originalMatrix[i, j + 2],
                            _originalMatrix[i, j + 3]
                        );
                        var vX = Vector256.LoadUnsafe(ref x[j]);
                        vSum = Fma.MultiplyAdd(vA, vX, vSum);
                    }

                    var sum = Vector256.Sum(vSum);
                    for (; j < n; j++)
                        sum += _originalMatrix[i, j] * x[j];

                    residual[i] = b._data[i] - sum;
                }
            else
                for (var i = 0; i < n; i++)
                {
                    var sum = 0.0;
                    for (var j = 0; j < n; j++)
                        sum += _originalMatrix[i, j] * x[j];
                    residual[i] = b._data[i] - sum;
                }

            var residualNorm = 0.0;
            var bNorm = 0.0;
            for (var i = 0; i < n; i++)
            {
                residualNorm += residual[i] * residual[i];
                bNorm += b._data[i] * b._data[i];
            }

            residualNorm = Math.Sqrt(residualNorm);
            bNorm = Math.Sqrt(bNorm);

            if (residualNorm < convergenceTolerance * bNorm)
                break;

            var residualVector = new Vector(residual);
            var dx = SolveWithoutRefinement(residualVector);

            for (var i = 0; i < n; i++)
                x[i] += dx[i];
        }

        return new Vector(x);
    }

    /// <summary>
    ///     Solves A*x = b without iterative refinement (used internally).
    /// </summary>
    private double[] SolveWithoutRefinement(Vector b)
    {
        var n = _l.RowCount;
        var x = new double[n];
        var y = new double[n];

        var pb = new double[n];
        for (var i = 0; i < n; i++)
            pb[i] = b._data[_rowPivots[i]];

        for (var i = 0; i < n; i++)
        {
            y[i] = pb[i];
            for (var j = 0; j < i; j++)
                y[i] -= _l[i, j] * y[j];
        }

        var matrixScale = _originalMatrix.OneNorm();
        var singularityThreshold = MachinePrecision * matrixScale;
        if (singularityThreshold < DeterminantUnderflowThreshold)
            singularityThreshold = DeterminantUnderflowThreshold;

        for (var i = n - 1; i >= 0; i--)
        {
            x[i] = y[i];
            for (var j = i + 1; j < n; j++)
                x[i] -= _u[i, j] * x[j];

            if (Math.Abs(_u[i, i]) < singularityThreshold)
                throw new InvalidOperationException(
                    $"Matrix is singular or nearly singular at U[{i},{i}] " +
                    $"(|U[{i},{i}]| = {Math.Abs(_u[i, i]):E3} < {singularityThreshold:E3}). Cannot solve the system.");

            x[i] /= _u[i, i];
        }

        var result = new double[n];
        for (var i = 0; i < n; i++)
            result[_columnPivots[i]] = x[i];

        return result;
    }

    #region Constants

    private const double MachinePrecision = 2.220446049250313e-16;
    private const double DeterminantUnderflowThreshold = 1e-300;

    #endregion

    #region Fields

    private readonly int[] _columnPivots;

    private readonly Matrix _l;
    private readonly int[] _rowPivots;
    private readonly Matrix _u;
    private readonly Matrix _originalMatrix;

    #endregion
}

public sealed class QRDecomposition
{
    private readonly double _normA;
    private readonly Matrix _qr;
    private readonly double[] _rDiag;

    internal QRDecomposition(Matrix qr, double[] rDiag, double normA)
    {
        _qr = qr;
        _rDiag = rDiag;
        _normA = normA;
    }

    public Matrix Q
    {
        get
        {
            var m = _qr.RowCount;
            var n = _qr.ColumnCount;
            var q = Matrix.Identity(m);

            for (var k = n - 1; k >= 0; k--)
            for (var j = k; j < m; j++)
                if (_qr[k, k] != 0)
                {
                    double s = 0;
                    for (var i = k; i < m; i++)
                        s += _qr[i, k] * q[i, j];
                    s = -s / _qr[k, k];
                    for (var i = k; i < m; i++)
                        q[i, j] += s * _qr[i, k];
                }

            return q;
        }
    }

    public Matrix R
    {
        get
        {
            var n = _qr.ColumnCount;
            var r = new Matrix(n, n);

            for (var i = 0; i < n; i++)
            {
                r[i, i] = _rDiag[i];
                for (var j = i + 1; j < n; j++)
                    r[i, j] = _qr[i, j];
            }

            return r;
        }
    }

    /// <summary>
    ///     Solves A*x = b using the QR decomposition (least squares).
    /// </summary>
    /// <param name="b">The right-hand side column vector (<see cref="Vector" />).</param>
    /// <returns>The least-squares solution vector x (<see cref="Vector" />).</returns>
    public Vector Solve(Vector b)
    {
        var m = _qr.RowCount;
        var n = _qr.ColumnCount;

        if (b.Length != m)
            throw new ArgumentException("Vector length must match matrix rows");

        var x = b.Clone()._data;

        for (var k = 0; k < n; k++)
        {
            double s = 0;
            for (var i = k; i < m; i++)
                s += _qr[i, k] * x[i];
            s = -s / _qr[k, k];
            for (var i = k; i < m; i++)
                x[i] += s * _qr[i, k];
        }

        const double machinePrecision = 2.220446049250313e-16;
        var tolerance = _normA * machinePrecision;

        for (var i = n - 1; i >= 0; i--)
        {
            if (Math.Abs(_rDiag[i]) < tolerance)
            {
                x[i] = 0;
                continue;
            }

            x[i] /= _rDiag[i];
            for (var j = 0; j < i; j++)
                x[j] -= x[i] * _qr[j, i];
        }

        var resultData = new double[n];
        Array.Copy(x, resultData, n);
        return new Vector(resultData);
    }
}

public sealed class EigenDecomposition
{
    private readonly double[] _eigenvalues;
    private readonly Matrix _eigenvectors;

    internal EigenDecomposition(double[] eigenvalues, Matrix eigenvectors)
    {
        _eigenvalues = eigenvalues;
        _eigenvectors = eigenvectors;
    }

    /// <summary>
    ///     Gets a copy of the eigenvalues array.
    /// </summary>
    public double[] Eigenvalues => (double[])_eigenvalues.Clone();

    /// <summary>
    ///     Gets a copy of the eigenvectors matrix.
    /// </summary>
    /// <remarks>
    ///     Each column represents an eigenvector corresponding to the eigenvalue at the same index.
    ///     Returns a clone to prevent external modification.
    /// </remarks>
    public Matrix Eigenvectors => _eigenvectors.Clone();
}

public sealed class SVDDecomposition
{
    private readonly Matrix _s;
    private readonly double[] _singularValues;
    private readonly Matrix _u;
    private readonly Matrix _v;

    internal SVDDecomposition(Matrix u, Matrix s, Matrix v)
    {
        _u = u;
        _s = s;
        _v = v;

        var minDim = Math.Min(s.RowCount, s.ColumnCount);
        _singularValues = new double[minDim];
        for (var i = 0; i < minDim; i++)
            _singularValues[i] = s[i, i];
    }

    /// <summary>
    ///     Gets a copy of the U matrix (left singular vectors).
    /// </summary>
    public Matrix U => _u.Clone();

    /// <summary>
    ///     Gets a copy of the S matrix (singular values on diagonal).
    /// </summary>
    public Matrix S => _s.Clone();

    /// <summary>
    ///     Gets a copy of the V matrix (right singular vectors).
    /// </summary>
    public Matrix V => _v.Clone();

    /// <summary>
    ///     Gets a copy of the singular values array.
    /// </summary>
    public double[] SingularValues => (double[])_singularValues.Clone();

    public int Rank(double tolerance = 1e-10)
    {
        var rank = 0;
        for (var i = 0; i < _singularValues.Length; i++)
            if (_singularValues[i] > tolerance)
                rank++;
        return rank;
    }

    public double ConditionNumber()
    {
        if (_singularValues.Length == 0 || Math.Abs(_singularValues[^1]) < 1e-15)
            return double.PositiveInfinity;
        return _singularValues[0] / _singularValues[^1];
    }

    /// <summary>
    ///     Solves A*x = b using the SVD (least squares).
    /// </summary>
    /// <param name="b">The right-hand side column vector (<see cref="Vector" />).</param>
    /// <returns>The least-squares solution vector x (<see cref="Vector" />).</returns>
    public Vector Solve(Vector b)
    {
        var m = _u.RowCount;
        var n = _v.RowCount;

        if (b.Length != m)
            throw new ArgumentException("Vector length must match matrix rows");

        var x = new double[n];
        const double relativeTolerance = 1e-10;
        var tolerance = relativeTolerance * (_singularValues.Length > 0 ? _singularValues[0] : 0.0);

        var utb = new double[_u.ColumnCount];
        for (var i = 0; i < _u.ColumnCount; i++)
        {
            double sum = 0;
            for (var j = 0; j < m; j++)
                sum += _u[j, i] * b._data[j];
            utb[i] = sum;
        }

        for (var i = 0; i < Math.Min(n, _singularValues.Length); i++)
            if (_singularValues[i] > tolerance)
                utb[i] /= _singularValues[i];
            else
                utb[i] = 0;

        for (var i = 0; i < n; i++)
        {
            double sum = 0;
            for (var j = 0; j < Math.Min(n, utb.Length); j++)
                sum += _v[i, j] * utb[j];
            x[i] = sum;
        }

        return new Vector(x);
    }
}

#endregion

#region Optimized GEMM Kernels

/// <summary>
///     Ultra-optimized matrix multiplication kernels - fastest possible C# implementation
/// </summary>
internal static class OptimizedGEMM
{
    // Cache parameters tuned for modern CPUs
    private const int L1_CACHE = 32 * 1024; // 32 KB L1
    private const int L2_CACHE = 256 * 1024; // 256 KB L2
    private const int L3_CACHE = 8 * 1024 * 1024; // 8 MB L3

    // Block sizes optimized for cache hierarchy
    private const int MC = 256; // M-dimension blocking (fits in L2 with B panel)
    private const int KC = 256; // K-dimension blocking (A and B panels fit in L2)
    private const int NC = 4096; // N-dimension blocking (B panel fits in L3)

    // Micro-kernel dimensions (register blocking)
    private const int MR = 8; // Rows per micro-kernel (8 AVX-512 registers)
    private const int NR = 6; // Cols per micro-kernel (6 AVX-512 registers)

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static unsafe void MultiplyOptimized(
        double* A, int lda, bool transA,
        double* B, int ldb, bool transB,
        double* C, int ldc,
        int M, int N, int K)
    {
        // Fast paths for small matrices
        if (M <= 4 && N <= 4 && K <= 4)
        {
            MultiplySmall(A, lda, transA, B, ldb, transB, C, ldc, M, N, K);
            return;
        }

        // For large matrices, use cache-blocked algorithm
        if (!transA && !transB)
            GEMM_NN_Blocked(A, lda, B, ldb, C, ldc, M, N, K);
        else if (transA && !transB)
            GEMM_TN_Blocked(A, lda, B, ldb, C, ldc, M, N, K);
        else if (!transA && transB)
            GEMM_NT_Blocked(A, lda, B, ldb, C, ldc, M, N, K);
        else
            GEMM_TT_Blocked(A, lda, B, ldb, C, ldc, M, N, K);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static unsafe void MultiplySmall(
        double* A, int lda, bool transA,
        double* B, int ldb, bool transB,
        double* C, int ldc,
        int M, int N, int K)
    {
        // Specialized for matrices up to 4x4 - completely unrolled
        // These are extremely common in FEA and graphics

        if (M == 3 && N == 3 && K == 3 && !transA && !transB)
        {
            // 3x3 matrix multiply - fully unrolled, no loops
            Multiply3x3(A, lda, B, ldb, C, ldc);
            return;
        }

        if (M == 4 && N == 4 && K == 4 && !transA && !transB)
        {
            // 4x4 matrix multiply - fully unrolled
            Multiply4x4(A, lda, B, ldb, C, ldc);
            return;
        }

        // Generic small multiply with minimal overhead
        for (var j = 0; j < N; j++)
        {
            var colC = C + j * ldc;
            for (var i = 0; i < M; i++)
            {
                var sum = 0.0;
                for (var k = 0; k < K; k++)
                {
                    var a_val = transA ? A[i * lda + k] : A[k * lda + i];
                    var b_val = transB ? B[j * ldb + k] : B[k * ldb + j];
                    sum += a_val * b_val;
                }

                colC[i] = sum;
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static unsafe void Multiply3x3(double* A, int lda, double* B, int ldb, double* C, int ldc)
    {
        // 3x3 multiply - completely unrolled (27 multiply-adds)
        // Access patterns optimized for column-major layout

        // Column 0 of C
        double a00 = A[0], a10 = A[1], a20 = A[2];
        double a01 = A[lda], a11 = A[lda + 1], a21 = A[lda + 2];
        double a02 = A[2 * lda], a12 = A[2 * lda + 1], a22 = A[2 * lda + 2];

        double b00 = B[0], b10 = B[1], b20 = B[2];

        C[0] = a00 * b00 + a01 * b10 + a02 * b20;
        C[1] = a10 * b00 + a11 * b10 + a12 * b20;
        C[2] = a20 * b00 + a21 * b10 + a22 * b20;

        // Column 1 of C
        double b01 = B[ldb], b11 = B[ldb + 1], b21 = B[ldb + 2];

        C[ldc] = a00 * b01 + a01 * b11 + a02 * b21;
        C[ldc + 1] = a10 * b01 + a11 * b11 + a12 * b21;
        C[ldc + 2] = a20 * b01 + a21 * b11 + a22 * b21;

        // Column 2 of C
        double b02 = B[2 * ldb], b12 = B[2 * ldb + 1], b22 = B[2 * ldb + 2];

        C[2 * ldc] = a00 * b02 + a01 * b12 + a02 * b22;
        C[2 * ldc + 1] = a10 * b02 + a11 * b12 + a12 * b22;
        C[2 * ldc + 2] = a20 * b02 + a21 * b12 + a22 * b22;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static unsafe void Multiply4x4(double* A, int lda, double* B, int ldb, double* C, int ldc)
    {
        // 4x4 multiply with AVX2 - process 4 doubles at once
        if (Avx.IsSupported)
            for (var j = 0; j < 4; j++)
            {
                var colC = C + j * ldc;

                // Load column j of B into vector
                var b0 = Vector256.Create(B[j * ldb]);
                var b1 = Vector256.Create(B[j * ldb + 1]);
                var b2 = Vector256.Create(B[j * ldb + 2]);
                var b3 = Vector256.Create(B[j * ldb + 3]);

                // Load columns of A and multiply-add
                var a0 = Avx.LoadVector256(A);
                var a1 = Avx.LoadVector256(A + lda);
                var a2 = Avx.LoadVector256(A + 2 * lda);
                var a3 = Avx.LoadVector256(A + 3 * lda);

                var result = Avx.Multiply(a0, b0);
                result = Fma.IsSupported
                    ? Fma.MultiplyAdd(a1, b1, result)
                    : Avx.Add(Avx.Multiply(a1, b1), result);
                result = Fma.IsSupported
                    ? Fma.MultiplyAdd(a2, b2, result)
                    : Avx.Add(Avx.Multiply(a2, b2), result);
                result = Fma.IsSupported
                    ? Fma.MultiplyAdd(a3, b3, result)
                    : Avx.Add(Avx.Multiply(a3, b3), result);

                Avx.Store(colC, result);
            }
        else
            // Scalar fallback for 4x4 (no AVX available)
            for (var j = 0; j < 4; j++)
            {
                var colC = C + j * ldc;
                var colB = B + j * ldb;

                for (var i = 0; i < 4; i++)
                {
                    var sum = 0.0;
                    for (var k = 0; k < 4; k++) sum += A[k * lda + i] * colB[k];

                    colC[i] = sum;
                }
            }
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static unsafe void GEMM_NN_Blocked(
        double* A, int lda,
        double* B, int ldb,
        double* C, int ldc,
        int M, int N, int K)
    {
        // Three-level cache blocking: L3 -> L2 -> L1 -> registers
        // C(M x N) = A(M x K) * B(K x N)

        for (var jc = 0; jc < N; jc += NC)
        {
            var nc = Math.Min(NC, N - jc);

            for (var pc = 0; pc < K; pc += KC)
            {
                var kc = Math.Min(KC, K - pc);

                for (var ic = 0; ic < M; ic += MC)
                {
                    var mc = Math.Min(MC, M - ic);

                    // Micro-kernel: compute C[ic:ic+mc, jc:jc+nc] += A[ic:ic+mc, pc:pc+kc] * B[pc:pc+kc, jc:jc+nc]
                    GEMM_NN_MicroPanel(
                        A + ic + pc * lda, lda,
                        B + pc + jc * ldb, ldb,
                        C + ic + jc * ldc, ldc,
                        mc, nc, kc);
                }
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static unsafe void GEMM_NN_MicroPanel(
        double* A, int lda,
        double* B, int ldb,
        double* C, int ldc,
        int m, int n, int k)
    {
        // Process in MR x NR tiles (8x6 with AVX-512)
        var i = 0;

        if (Avx512F.IsSupported)
        {
            // AVX-512 micro-kernel: 8 rows x 6 columns
            for (; i + MR - 1 < m; i += MR)
            {
                var j = 0;
                for (; j + NR - 1 < n; j += NR)
                    MicroKernel_8x6_AVX512(
                        A + i, lda,
                        B + j * ldb, ldb,
                        C + i + j * ldc, ldc,
                        k);

                // Handle remaining columns
                for (; j < n; j++)
                    MicroKernel_8x1_AVX512(
                        A + i, lda,
                        B + j * ldb, ldb,
                        C + i + j * ldc, ldc,
                        k);
            }
        }
        else if (Avx2.IsSupported)
        {
            // AVX2 micro-kernel: 4 rows x 4 columns
            const int mr_avx2 = 4;
            const int nr_avx2 = 4;

            for (; i + mr_avx2 - 1 < m; i += mr_avx2)
            {
                var j = 0;
                for (; j + nr_avx2 - 1 < n; j += nr_avx2)
                    MicroKernel_4x4_AVX2(
                        A + i, lda,
                        B + j * ldb, ldb,
                        C + i + j * ldc, ldc,
                        k);

                for (; j < n; j++)
                    MicroKernel_4x1_AVX2(
                        A + i, lda,
                        B + j * ldb, ldb,
                        C + i + j * ldc, ldc,
                        k);
            }
        }

        // Handle remaining rows with scalar code
        for (; i < m; i++)
        for (var j = 0; j < n; j++)
        {
            var sum = C[i + j * ldc];
            for (var p = 0; p < k; p++) sum += A[i + p * lda] * B[p + j * ldb];

            C[i + j * ldc] = sum;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    private static unsafe void MicroKernel_8x6_AVX512(
        double* A, int lda,
        double* B, int ldb,
        double* C, int ldc,
        int k)
    {
        // Ultra-optimized 8x6 micro-kernel using AVX-512
        // Computes C[0:7, 0:5] += A[0:7, 0:k-1] * B[0:k-1, 0:5]

        // Load C into 6 AVX-512 registers (8 doubles each)
        var c0 = Avx512F.LoadVector512(C);
        var c1 = Avx512F.LoadVector512(C + ldc);
        var c2 = Avx512F.LoadVector512(C + 2 * ldc);
        var c3 = Avx512F.LoadVector512(C + 3 * ldc);
        var c4 = Avx512F.LoadVector512(C + 4 * ldc);
        var c5 = Avx512F.LoadVector512(C + 5 * ldc);

        // Main loop: process K dimension
        for (var p = 0; p < k; p++)
        {
            // Load column p of A (8 elements)
            var a = Avx512F.LoadVector512(A + p * lda);

            // Prefetch next iteration (heuristic: 8 iterations ahead)
            if (p + 8 < k)
            {
                Sse.Prefetch0(A + (p + 8) * lda);
                Sse.Prefetch0(B + (p + 8));
            }

            // Broadcast B[p, j] and perform FMA
            var b0 = Vector512.Create(B[p]);
            var b1 = Vector512.Create(B[p + ldb]);
            var b2 = Vector512.Create(B[p + 2 * ldb]);
            var b3 = Vector512.Create(B[p + 3 * ldb]);
            var b4 = Vector512.Create(B[p + 4 * ldb]);
            var b5 = Vector512.Create(B[p + 5 * ldb]);

            c0 = Avx512F.FusedMultiplyAdd(a, b0, c0);
            c1 = Avx512F.FusedMultiplyAdd(a, b1, c1);
            c2 = Avx512F.FusedMultiplyAdd(a, b2, c2);
            c3 = Avx512F.FusedMultiplyAdd(a, b3, c3);
            c4 = Avx512F.FusedMultiplyAdd(a, b4, c4);
            c5 = Avx512F.FusedMultiplyAdd(a, b5, c5);
        }

        // Store results back to C
        Avx512F.Store(C, c0);
        Avx512F.Store(C + ldc, c1);
        Avx512F.Store(C + 2 * ldc, c2);
        Avx512F.Store(C + 3 * ldc, c3);
        Avx512F.Store(C + 4 * ldc, c4);
        Avx512F.Store(C + 5 * ldc, c5);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static unsafe void MicroKernel_8x1_AVX512(
        double* A, int lda,
        double* B, int ldb,
        double* C, int ldc,
        int k)
    {
        var c = Avx512F.LoadVector512(C);

        for (var p = 0; p < k; p++)
        {
            var a = Avx512F.LoadVector512(A + p * lda);
            var b = Vector512.Create(B[p]);
            c = Avx512F.FusedMultiplyAdd(a, b, c);
        }

        Avx512F.Store(C, c);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    private static unsafe void MicroKernel_4x4_AVX2(
        double* A, int lda,
        double* B, int ldb,
        double* C, int ldc,
        int k)
    {
        // Optimized 4x4 micro-kernel using AVX2
        var c0 = Avx.LoadVector256(C);
        var c1 = Avx.LoadVector256(C + ldc);
        var c2 = Avx.LoadVector256(C + 2 * ldc);
        var c3 = Avx.LoadVector256(C + 3 * ldc);

        for (var p = 0; p < k; p++)
        {
            var a = Avx.LoadVector256(A + p * lda);

            var b0 = Vector256.Create(B[p]);
            var b1 = Vector256.Create(B[p + ldb]);
            var b2 = Vector256.Create(B[p + 2 * ldb]);
            var b3 = Vector256.Create(B[p + 3 * ldb]);

            if (Fma.IsSupported)
            {
                c0 = Fma.MultiplyAdd(a, b0, c0);
                c1 = Fma.MultiplyAdd(a, b1, c1);
                c2 = Fma.MultiplyAdd(a, b2, c2);
                c3 = Fma.MultiplyAdd(a, b3, c3);
            }
            else
            {
                c0 = Avx.Add(Avx.Multiply(a, b0), c0);
                c1 = Avx.Add(Avx.Multiply(a, b1), c1);
                c2 = Avx.Add(Avx.Multiply(a, b2), c2);
                c3 = Avx.Add(Avx.Multiply(a, b3), c3);
            }
        }

        Avx.Store(C, c0);
        Avx.Store(C + ldc, c1);
        Avx.Store(C + 2 * ldc, c2);
        Avx.Store(C + 3 * ldc, c3);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static unsafe void MicroKernel_4x1_AVX2(
        double* A, int lda,
        double* B, int ldb,
        double* C, int ldc,
        int k)
    {
        var c = Avx.LoadVector256(C);

        for (var p = 0; p < k; p++)
        {
            var a = Avx.LoadVector256(A + p * lda);
            var b = Vector256.Create(B[p]);
            c = Fma.IsSupported
                ? Fma.MultiplyAdd(a, b, c)
                : Avx.Add(Avx.Multiply(a, b), c);
        }

        Avx.Store(C, c);
    }

    // Placeholder for other transpose variants (TN, NT, TT)
    // GEMM with A transposed: C = A^T * B
    // A is stored as (K x M) but accessed as transposed to (M x K)
    // B is stored as (K x N)
    // C is stored as (M x N)
    private static unsafe void GEMM_TN_Blocked(double* A, int lda, double* B, int ldb, double* C, int ldc, int M, int N,
        int K)
    {
        // For transposed A, we need to access A[k, i] which is stored at A[i * lda + k]
        // Fallback to scalar implementation for correctness
        for (var i = 0; i < M; i++)
        for (var j = 0; j < N; j++)
        {
            var sum = 0.0;
            for (var k = 0; k < K; k++)
            {
                // A^T[i, k] = A[k, i] stored at A[i * lda + k]
                var a_val = A[i * lda + k];
                // B[k, j] stored at B[k + j * ldb]
                var b_val = B[k + j * ldb];
                sum += a_val * b_val;
            }

            C[i + j * ldc] = sum;
        }
    }

    // GEMM with B transposed: C = A * B^T
    // A is stored as (M x K)
    // B is stored as (N x K) but accessed as transposed to (K x N)
    // C is stored as (M x N)
    private static unsafe void GEMM_NT_Blocked(double* A, int lda, double* B, int ldb, double* C, int ldc, int M, int N,
        int K)
    {
        // For transposed B, we need to access B^T[k, j] = B[j, k] which is stored at B[k * ldb + j]
        for (var i = 0; i < M; i++)
        for (var j = 0; j < N; j++)
        {
            var sum = 0.0;
            for (var k = 0; k < K; k++)
            {
                // A[i, k] stored at A[k * lda + i]
                var a_val = A[k * lda + i];
                // B^T[k, j] = B[j, k] stored at B[k * ldb + j]
                var b_val = B[k * ldb + j];
                sum += a_val * b_val;
            }

            C[i + j * ldc] = sum;
        }
    }

    // GEMM with both transposed: C = A^T * B^T
    // A is stored as (K x M) but accessed as transposed to (M x K)
    // B is stored as (N x K) but accessed as transposed to (K x N)
    // C is stored as (M x N)
    private static unsafe void GEMM_TT_Blocked(double* A, int lda, double* B, int ldb, double* C, int ldc, int M, int N,
        int K)
    {
        for (var i = 0; i < M; i++)
        for (var j = 0; j < N; j++)
        {
            var sum = 0.0;
            for (var k = 0; k < K; k++)
            {
                // A^T[i, k] = A[k, i] stored at A[i * lda + k]
                var a_val = A[i * lda + k];
                // B^T[k, j] = B[j, k] stored at B[k * ldb + j]
                var b_val = B[k * ldb + j];
                sum += a_val * b_val;
            }

            C[i + j * ldc] = sum;
        }
    }
}

#endregion
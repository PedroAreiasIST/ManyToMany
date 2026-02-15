using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using static System.Runtime.CompilerServices.MethodImplOptions;

namespace Numerical;

/// <summary>
///     Ultra-optimized sparse matrix-vector multiplication kernels.
///     Provides AVX-512, AVX2, and scalar fallback paths.
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
        var final2 = Sse2.Add(sum, high1);

        return final2.ToScalar();
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
///     Optimized sparse matrix-matrix multiply and other CSR operations.
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
                    var final2 = Sse2.Add(sum2, high1);
                    sum = final2.ToScalar();

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

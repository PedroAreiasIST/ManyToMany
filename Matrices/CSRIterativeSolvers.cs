using System.Runtime.CompilerServices;

namespace Numerical;

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
    private const double DEFAULT_TOLERANCE = CSR.DEFAULT_TOLERANCE;

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
        ArgumentNullException.ThrowIfNull(matrix);
        matrix.ThrowIfDisposed();

        if (matrix.Rows != matrix.Columns)
            throw new InvalidOperationException("Eigenvalue computation requires a square matrix.");
        if (m <= 0)
            throw new ArgumentOutOfRangeException(nameof(m), "Number of requested eigenvalues must be positive.");
        if (m > matrix.Rows)
            throw new ArgumentOutOfRangeException(nameof(m),
                $"Cannot request {m} eigenvalues from a {matrix.Rows}x{matrix.Columns} matrix.");
        if (tolerance <= 0)
            throw new ArgumentOutOfRangeException(nameof(tolerance), "Tolerance must be positive.");
        if (maxIterations <= 0)
            throw new ArgumentOutOfRangeException(nameof(maxIterations), "Maximum iterations must be positive.");

        var n = matrix.Rows;
        var basis = new double[m][];
        var current = new double[m][];
        var eigenvalues = new double[m];

        // Deterministic initialization to keep behavior reproducible.
        for (var k = 0; k < m; k++)
        {
            var v = new double[n];
            for (var i = 0; i < n; i++) v[i] = Math.Sin((k + 1) * (i + 1));
            basis[k] = v;
            current[k] = new double[n];
        }

        OrthogonalizeAndNormalize(basis, m, n);

        for (var iteration = 0; iteration < maxIterations; iteration++)
        {
            var maxLambdaDelta = 0.0;

            for (var k = 0; k < m; k++)
            {
                var solve = matrix.TrySolve(basis[k], tolerance * 0.1, Math.Max(200, matrix.Rows * 2));
                if (!solve.Converged)
                    throw new SolverException(
                        $"Inverse iteration failed while computing eigenvalue #{k + 1}: {solve.Message ?? "solver did not converge"}");

                Array.Copy(solve.Solution, current[k], n);
            }

            OrthogonalizeAndNormalize(current, m, n);

            for (var k = 0; k < m; k++)
            {
                var Av = matrix.Multiply(current[k]);
                var numerator = Dot(current[k], Av, n);
                var denominator = Dot(current[k], current[k], n);

                if (Math.Abs(denominator) <= DEFAULT_TOLERANCE)
                    throw new SolverException("Encountered near-zero eigenvector norm during eigenvalue estimation.");

                var lambda = numerator / denominator;
                maxLambdaDelta = Math.Max(maxLambdaDelta, Math.Abs(lambda - eigenvalues[k]));
                eigenvalues[k] = lambda;
            }

            // Keep iteration basis updated.
            for (var k = 0; k < m; k++) Array.Copy(current[k], basis[k], n);

            if (maxLambdaDelta < tolerance) break;
        }

        Array.Sort(eigenvalues);
        return eigenvalues;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Dot(double[] a, double[] b, int n)
    {
        var sum = 0.0;
        for (var i = 0; i < n; i++) sum += a[i] * b[i];
        return sum;
    }

    private static void OrthogonalizeAndNormalize(double[][] vectors, int count, int length)
    {
        for (var i = 0; i < count; i++)
        {
            var vi = vectors[i];

            for (var j = 0; j < i; j++)
            {
                var vj = vectors[j];
                var projection = Dot(vi, vj, length);
                for (var p = 0; p < length; p++) vi[p] -= projection * vj[p];
            }

            var norm = Math.Sqrt(Dot(vi, vi, length));
            if (norm <= DEFAULT_TOLERANCE)
                throw new SolverException("Failed to build a stable orthonormal basis for eigenvalue computation.");

            var invNorm = 1.0 / norm;
            for (var p = 0; p < length; p++) vi[p] *= invNorm;
        }
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

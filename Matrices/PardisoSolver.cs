using System.Diagnostics;
using System.Runtime.InteropServices;

namespace Numerical;

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
}

using System.Diagnostics;
using System.Runtime.InteropServices;

namespace Numerical;

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
    private static INativeLibraryConfig _libraryConfig = new NativeLibraryConfig();
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
    public static void SetLibraryConfig(INativeLibraryConfig config)
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
            _libCudaRt = RobustNativeLibraryLoader.TryLoadLibrary(rtNames, searchPaths, "CUDA Runtime");

            if (_libCudaRt != IntPtr.Zero)
            {
                TryGetExport(_libCudaRt, "cudaMalloc", out _cudaMalloc);
                TryGetExport(_libCudaRt, "cudaFree", out _cudaFree);
                TryGetExport(_libCudaRt, "cudaMemcpy", out _cudaMemcpy);
                TryGetExport(_libCudaRt, "cudaDeviceSynchronize", out _cudaDeviceSynchronize);
            }

            // Load cuSPARSE
            _libCuSparse = RobustNativeLibraryLoader.TryLoadLibrary(spNames, searchPaths, "cuSPARSE");

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

    // FIXED: MEM-C2 - Finalizer must NOT call CUDA P/Invoke functions.
    // The CUDA context may be torn down during process shutdown, causing crashes.
    // Only clean up managed resources; native CUDA resources are abandoned if not disposed.
    ~CuSparseBackend()
    {
        lock (_disposeLock)
        {
            isDisposed = true;
            isInitialized = false;
        }
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

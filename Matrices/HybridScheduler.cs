using System.Collections.Concurrent;
using System.Diagnostics;

namespace Numerical;

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
                    Debug.WriteLine($"GPU Eager init failed: {ex.Message}");
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

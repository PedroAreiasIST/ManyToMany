using System.Diagnostics;
using System.Reflection;
using System.Runtime.InteropServices;

namespace Numerical;

/// <summary>
///     Global configuration for parallel operations.
///     Extended with GPU and MKL thread control.
/// </summary>
public static class ParallelConfig
{
    // ========================================================================
    // PRIVATE FIELDS
    // ========================================================================

    private static readonly object _lock = new();
    private static int _maxDegreeOfParallelism = Environment.ProcessorCount;

    // GPU caching
    private static bool? _isGPUAvailable;
    private static bool _skipGPUCheck; // Skip GPU check entirely if it causes crashes

    // MKL state
    private static int _mklNumThreads = Environment.ProcessorCount;
    private static IntPtr _mklHandle = IntPtr.Zero;
    private static MklSetNumThreadsDelegate? _mkl_set_num_threads;

    private static MklGetMaxThreadsDelegate? _mkl_get_max_threads;

    // P0-3 FIX: Split into two flags to allow retry after Cleanup()
    private static volatile bool _mklLoadAttempted;
    private static volatile bool _mklLoadedSuccessfully;

    // Debug logging

    // ========================================================================
    // STATIC CONSTRUCTOR
    // ========================================================================

    static ParallelConfig()
    {
        // Don't apply MKL threads on startup - wait until first use
        // This prevents segfaults if MKL is not available
        // ApplyMKLThreads() will be called on first MKL access
    }

    // ========================================================================
    // CPU PARALLELISM CONTROL
    // ========================================================================

    /// <summary>
    ///     Gets or sets whether debug output is enabled for diagnostics.
    ///     When true, diagnostic messages are written to System.Diagnostics.Debug.
    ///     Default: false
    /// </summary>
    public static bool EnableDebugOutput { get; set; }

    /// <summary>
    ///     Maximum degree of parallelism for all parallel operations.
    ///     Set to Environment.ProcessorCount by default.
    ///     Minimum value is 1 (no parallelism).
    /// </summary>
    public static int MaxDegreeOfParallelism
    {
        get => _maxDegreeOfParallelism;
        set
        {
            var newValue = Math.Max(1, value);
            if (newValue != _maxDegreeOfParallelism)
            {
                _maxDegreeOfParallelism = newValue;
                // P2.2 FIX: Invalidate cached ParallelOptions when config changes
                _cachedOptions = null;
            }
        }
    }

    // P2.2 FIX: Cache ParallelOptions to avoid allocation on every hot-path call
    private static volatile ParallelOptions? _cachedOptions;

    /// <summary>
    ///     Gets ParallelOptions configured with the global MaxDegreeOfParallelism.
    ///     P2.2 FIX: Cached to avoid allocation per call in hot paths.
    /// </summary>
    public static ParallelOptions Options
    {
        get
        {
            var opts = _cachedOptions;
            if (opts != null && opts.MaxDegreeOfParallelism == _maxDegreeOfParallelism)
                return opts;
            opts = new ParallelOptions { MaxDegreeOfParallelism = _maxDegreeOfParallelism };
            _cachedOptions = opts;
            return opts;
        }
    }

    // ========================================================================
    // GPU CONTROL
    // ========================================================================

    /// <summary>
    ///     Gets or sets whether to skip GPU availability checking entirely.
    ///     Set this to true if GPU detection causes crashes.
    ///     Default: false
    /// </summary>
    public static bool SkipGPUCheck
    {
        get => _skipGPUCheck;
        set
        {
            _skipGPUCheck = value;
            if (value)
                // Clear cached value so it won't be checked
                _isGPUAvailable = false;
        }
    }

    /// <summary>
    ///     Gets or sets whether GPU acceleration is enabled globally.
    ///     When false, all operations use CPU regardless of GPU availability.
    ///     Default: true
    /// </summary>
    public static bool EnableGPU { get; set; } = true;

    /// <summary>
    ///     Returns true if GPU should be used (enabled AND available).
    /// </summary>
    public static bool UseGPU => EnableGPU && IsGPUAvailable;

    /// <summary>
    ///     Returns true if CUDA/cuSPARSE is available on this system.
    ///     Result is cached after first check.
    /// </summary>
    public static bool IsGPUAvailable
    {
        get
        {
            // If user explicitly disabled GPU checking, return false immediately
            if (_skipGPUCheck)
                return false;

            if (_isGPUAvailable.HasValue)
                return _isGPUAvailable.Value;

            lock (_lock)
            {
                // Re-check after acquiring lock
                if (_isGPUAvailable.HasValue)
                    return _isGPUAvailable.Value;

                var available = false;
                try
                {
                    // CRITICAL: Reflection can trigger static constructors that may crash
                    // Wrap the entire operation in multiple layers of protection
                    try
                    {
                        // P2-2 FIX: First try simple name (works if in same assembly)
                        var type = Type.GetType("Numerical.CuSparseBackend");

                        // P2-2 FIX: If not found, search all currently loaded assemblies
                        // This handles the case where CuSparseBackend is in a separate DLL
                        if (type == null)
                            try
                            {
                                foreach (var assembly in AppDomain.CurrentDomain.GetAssemblies())
                                {
                                    // Skip dynamic assemblies (can cause issues)
                                    if (assembly.IsDynamic)
                                        continue;

                                    type = assembly.GetType("Numerical.CuSparseBackend", false);
                                    if (type != null)
                                    {
                                        LogDebug($"Found CuSparseBackend in assembly: {assembly.FullName}");
                                        break;
                                    }
                                }
                            }
                            catch (Exception ex)
                            {
                                LogDebug($"Assembly search for GPU type failed: {ex.Message}");
                            }

                        if (type != null)
                            try
                            {
                                var prop = type.GetProperty("IsAvailable",
                                    BindingFlags.Public | BindingFlags.Static);
                                if (prop != null)
                                    try
                                    {
                                        var result = prop.GetValue(null);
                                        available = result is bool b && b;
                                    }
                                    catch (TargetInvocationException ex)
                                    {
                                        LogDebug(
                                            $"GPU property getter threw: {ex.InnerException?.Message ?? ex.Message}");
                                        available = false;
                                    }
                            }
                            catch (Exception ex)
                            {
                                LogDebug($"GPU property access failed: {ex.Message}");
                                available = false;
                            }
                    }
                    catch (TypeLoadException ex)
                    {
                        LogDebug($"GPU type load failed: {ex.Message}");
                        available = false;
                    }
                    catch (FileNotFoundException ex)
                    {
                        LogDebug($"GPU assembly not found: {ex.Message}");
                        available = false;
                    }
                }
                catch (Exception ex)
                {
                    LogDebug($"GPU availability check failed: {ex.Message}");
                    available = false;
                }

                _isGPUAvailable = available;
                return available;
            }
        }
    }

    // ========================================================================
    // MKL THREAD CONTROL
    // ========================================================================

    /// <summary>
    ///     Gets or sets the number of threads MKL uses internally.
    ///     Default: Environment.ProcessorCount
    ///     Valid range: 1 to 2*ProcessorCount (allows hyperthreading scenarios)
    /// </summary>
    public static int MKLNumThreads
    {
        get => _mklNumThreads;
        set
        {
            lock (_lock)
            {
                _mklNumThreads = Math.Clamp(value, 1, Environment.ProcessorCount * 2);
                ApplyMKLThreads();
            }
        }
    }

    /// <summary>
    ///     Gets the current number of threads MKL is using.
    ///     Returns null if MKL is not available.
    /// </summary>
    public static int? MKLCurrentThreads
    {
        get
        {
            EnsureMKLLoaded();

            // Only call if we have a valid function pointer
            if (_mkl_get_max_threads != null)
                try
                {
                    return _mkl_get_max_threads.Invoke();
                }
                catch (Exception ex)
                {
                    LogDebug($"Failed to get MKL threads: {ex.Message}");
                    return null;
                }

            return null;
        }
    }

    /// <summary>
    ///     Returns true if MKL is loaded and available.
    /// </summary>
    public static bool IsMKLAvailable => _mklLoadedSuccessfully && _mkl_set_num_threads != null;

    // ========================================================================
    // CONVENIENCE METHODS
    // ========================================================================

    /// <summary>
    ///     Shorthand for Environment.ProcessorCount.
    /// </summary>
    public static int ProcessorCount => Environment.ProcessorCount;

    /// <summary>
    ///     Writes a debug message if debug output is enabled.
    /// </summary>
    private static void LogDebug(string message)
    {
        if (EnableDebugOutput)
            Debug.WriteLine($"[ParallelConfig] {message}");
    }

    /// <summary>
    ///     Applies the current MKL thread count setting to the loaded MKL library.
    ///     Does nothing if MKL is not loaded or not available.
    /// </summary>
    private static void ApplyMKLThreads()
    {
        EnsureMKLLoaded();

        // Only call if we have a valid function pointer
        if (_mkl_set_num_threads != null)
            try
            {
                _mkl_set_num_threads.Invoke(_mklNumThreads);
            }
            catch (Exception ex)
            {
                LogDebug($"Failed to set MKL threads: {ex.Message}");
            }
    }

    /// <summary>
    ///     Ensures MKL library is loaded and function pointers are initialized.
    ///     Uses lazy loading with double-checked locking for thread safety.
    ///     Searches for platform-specific MKL library names and loads the first available.
    /// </summary>
    private static void EnsureMKLLoaded()
    {
        // P0-3 FIX: Use _mklLoadAttempted to allow retry after Cleanup()
        if (_mklLoadAttempted) return;

        lock (_lock)
        {
            if (_mklLoadAttempted) return;

            try
            {
                string[] mklNames;
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                    mklNames = new[] { "mkl_rt.2.dll", "mkl_rt.dll" };
                else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                    mklNames = new[] { "libmkl_rt.2.dylib", "libmkl_rt.dylib" };
                else // Linux
                    mklNames = new[] { "libmkl_rt.so.2", "libmkl_rt.so" };

                var loadedSuccessfully = false;

                foreach (var name in mklNames)
                    try
                    {
                        if (NativeLibrary.TryLoad(name, out var handle))
                        {
                            // Free previous handle if we're retrying with a different library
                            if (_mklHandle != IntPtr.Zero && _mklHandle != handle)
                                try
                                {
                                    NativeLibrary.Free(_mklHandle);
                                    LogDebug("Freed previous MKL handle before loading new one");
                                }
                                catch (Exception ex)
                                {
                                    LogDebug($"Failed to free previous MKL handle: {ex.Message}");
                                }

                            _mklHandle = handle;

                            // Platform-specific function name priority (lowercase more common on Linux/macOS)
                            var setFuncNames = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                                ? new[] { "MKL_Set_Num_Threads", "mkl_set_num_threads" }
                                : new[] { "mkl_set_num_threads", "MKL_Set_Num_Threads" };

                            var getFuncNames = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                                ? new[] { "MKL_Get_Max_Threads", "mkl_get_max_threads" }
                                : new[] { "mkl_get_max_threads", "MKL_Get_Max_Threads" };

                            foreach (var funcName in setFuncNames)
                                if (NativeLibrary.TryGetExport(handle, funcName, out var setPtr))
                                    try
                                    {
                                        _mkl_set_num_threads =
                                            Marshal.GetDelegateForFunctionPointer<MklSetNumThreadsDelegate>(setPtr);
                                        break;
                                    }
                                    catch (Exception ex)
                                    {
                                        LogDebug($"Failed to create delegate for {funcName}: {ex.Message}");
                                    }

                            foreach (var funcName in getFuncNames)
                                if (NativeLibrary.TryGetExport(handle, funcName, out var getPtr))
                                    try
                                    {
                                        _mkl_get_max_threads =
                                            Marshal.GetDelegateForFunctionPointer<MklGetMaxThreadsDelegate>(getPtr);
                                        break;
                                    }
                                    catch (Exception ex)
                                    {
                                        LogDebug($"Failed to create delegate for {funcName}: {ex.Message}");
                                    }

                            if (_mkl_set_num_threads != null)
                            {
                                LogDebug($"MKL loaded successfully: {name}");
                                loadedSuccessfully = true;
                                break;
                            }

                            // Couldn't get function pointers, free this handle
                            try
                            {
                                NativeLibrary.Free(handle);
                                _mklHandle = IntPtr.Zero;
                            }
                            catch (Exception ex)
                            {
                                LogDebug($"Failed to free unusable MKL handle: {ex.Message}");
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        LogDebug($"Failed to load {name}: {ex.Message}");
                    }

                if (!loadedSuccessfully)
                {
                    LogDebug("MKL not found or could not be loaded");
                    _mklHandle = IntPtr.Zero;
                    _mkl_set_num_threads = null;
                    _mkl_get_max_threads = null;
                }

                // P0-3 FIX: Set both flags - attempted always true, successful only if loaded
                _mklLoadAttempted = true;
                _mklLoadedSuccessfully = loadedSuccessfully;
            }
            catch (Exception ex)
            {
                LogDebug($"MKL loading failed: {ex.GetType().Name}: {ex.Message}");
                _mklHandle = IntPtr.Zero;
                _mkl_set_num_threads = null;
                _mkl_get_max_threads = null;
                _mklLoadAttempted = true;
                _mklLoadedSuccessfully = false;
            }
        }
    }

    /// <summary>
    ///     Cleans up native MKL resources.
    ///     Call this before application exit if you need deterministic cleanup.
    ///     Note: This is typically not necessary as the OS will clean up on process exit.
    ///     After Cleanup(), MKL can be reloaded on next use (P0-3 fix).
    /// </summary>
    public static void Cleanup()
    {
        lock (_lock)
        {
            if (_mklHandle != IntPtr.Zero)
                try
                {
                    NativeLibrary.Free(_mklHandle);
                    LogDebug("MKL handle freed successfully");
                }
                catch (Exception ex)
                {
                    LogDebug($"MKL cleanup failed: {ex.GetType().Name}: {ex.Message}");
                }
                finally
                {
                    _mklHandle = IntPtr.Zero;
                    _mkl_set_num_threads = null;
                    _mkl_get_max_threads = null;
                    // P0-3 FIX: Reset load state to allow retry after cleanup
                    _mklLoadAttempted = false;
                    _mklLoadedSuccessfully = false;
                    LogDebug("MKL state reset - reload will be attempted on next use");
                }
        }
    }

    /// <summary>
    ///     Sets both CPU and MKL thread counts to the same value.
    /// </summary>
    /// <param name="numThreads">Number of threads (must be between 1 and 2*ProcessorCount)</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown if numThreads is out of valid range</exception>
    public static void SetAllThreads(int numThreads)
    {
        var maxAllowed = Environment.ProcessorCount * 2;
        if (numThreads < 1 || numThreads > maxAllowed)
            throw new ArgumentOutOfRangeException(nameof(numThreads),
                $"Thread count must be between 1 and {maxAllowed}, got {numThreads}");

        MaxDegreeOfParallelism = numThreads;
        MKLNumThreads = numThreads;
    }

    /// <summary>
    ///     Resets all settings to defaults.
    /// </summary>
    public static void Reset()
    {
        MaxDegreeOfParallelism = Environment.ProcessorCount;
        EnableGPU = true;
        MKLNumThreads = Environment.ProcessorCount;
    }

    /// <summary>
    ///     Gets a summary of current configuration.
    /// </summary>
    public static string GetSummary()
    {
        var gpuAvail = IsGPUAvailable; // Cache to avoid multiple reflection calls
        var mklThreads = MKLCurrentThreads;
        var mklStatus = IsMKLAvailable
            ? $"MKL={_mklNumThreads} (actual={mklThreads?.ToString() ?? "?"})"
            : "MKL=unavailable";

        return $"ParallelConfig: CPU={MaxDegreeOfParallelism}/{ProcessorCount}, " +
               $"GPU={EnableGPU} (Avail={gpuAvail}), " +
               $"{mklStatus}";
    }

    // ========================================================================
    // DELEGATES
    // ========================================================================

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate void MklSetNumThreadsDelegate(int numThreads);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int MklGetMaxThreadsDelegate();
}
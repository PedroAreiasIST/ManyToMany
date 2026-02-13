using System.Diagnostics;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Text.RegularExpressions;

namespace Numerical;

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

/// <summary>
///     Configuration interface for native library paths.
///     Allows applications to specify custom library locations or versions.
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
// LIBRARY STATUS - Consolidated from LibraryAvailability.cs
// ============================================================================

/// <summary>
///     Detailed status of all native libraries.
///     Provides a unified view of library availability and configuration.
/// </summary>
public class LibraryStatus
{
    public bool CudaAvailable { get; set; }
    public int CudaVersion { get; set; }
    public int CudaErrorCode { get; set; }
    public bool CuSparseAvailable { get; set; }
    public bool MKLAvailable { get; set; }
    public bool PardisoAvailable { get; set; }

    public bool HasGpuAcceleration => CudaAvailable && CuSparseAvailable;
    public bool HasCpuAcceleration => MKLAvailable;
    public bool HasAnyAcceleration => HasGpuAcceleration || HasCpuAcceleration;

    public string CudaVersionString
    {
        get
        {
            if (!CudaAvailable) return "N/A";
            var major = CudaVersion / 1000;
            var minor = CudaVersion % 1000 / 10;
            return $"{major}.{minor}";
        }
    }

    public override string ToString()
    {
        return $"CUDA: {(CudaAvailable ? $"v{CudaVersionString}" : "No")}, " +
               $"cuSPARSE: {(CuSparseAvailable ? "Yes" : "No")}, " +
               $"MKL: {(MKLAvailable ? "Yes" : "No")}, " +
               $"PARDISO: {(PardisoAvailable ? "Yes" : "No")}";
    }
}

// ============================================================================
// LIBRARY AVAILABILITY - Simple Runtime Checker (consolidated from LibraryAvailability.cs)
// ============================================================================

/// <summary>
///     Runtime checker for native library availability.
///     This provides a simple API for checking library availability without caching complexity.
///     For detailed diagnostics, use <see cref="NativeLibraryStatus" />.
/// </summary>
public static class LibraryAvailability
{
    #region Public API

    /// <summary>
    ///     Check if CUDA Runtime is available and loadable
    /// </summary>
    public static bool IsCudaRuntimeAvailable()
    {
        return CheckLibrary(GetCudaRuntimePaths(), "cudaRuntimeGetVersion");
    }

    /// <summary>
    ///     Check if cuSPARSE is available and loadable
    /// </summary>
    public static bool IsCuSparseAvailable()
    {
        return CheckLibrary(GetCuSparsePaths(), "cusparseCreate");
    }

    /// <summary>
    ///     Check if Intel MKL is available and loadable
    /// </summary>
    public static bool IsMKLAvailable()
    {
        return CheckLibrary(GetMKLPaths(), new[] { "mkl_get_max_threads", "MKL_Get_Max_Threads" });
    }

    /// <summary>
    ///     Check if PARDISO is available (usually via MKL)
    /// </summary>
    public static bool IsPardisoAvailable()
    {
        return CheckLibrary(GetMKLPaths(), new[] { "pardiso", "pardiso_", "pardiso_64", "PARDISO" });
    }

    /// <summary>
    ///     Get CUDA version if available.
    ///     Returns (available, version, errorCode)
    /// </summary>
    public static (bool available, int version, int errorCode) GetCudaVersion()
    {
        foreach (var path in GetCudaRuntimePaths())
        {
            if (!File.Exists(path)) continue;

            try
            {
                if (NativeLibrary.TryLoad(path, out var handle))
                {
                    if (NativeLibrary.TryGetExport(handle, "cudaRuntimeGetVersion", out var funcPtr))
                    {
                        var getVersion = Marshal.GetDelegateForFunctionPointer<CudaGetVersionDelegate>(funcPtr);
                        var version = 0;
                        var result = getVersion(ref version);
                        NativeLibrary.Free(handle);

                        // Only accept success (0) - error codes 100 (NoDevice), 35 (InsufficientDriver),
                        // 3 (InitializationError) mean the library loaded but CUDA cannot work
                        if (result == 0)
                            return (true, version, result);

                        return (false, 0, result);
                    }

                    NativeLibrary.Free(handle);
                }
            }
            catch
            {
            }
        }

        return (false, 0, -1);
    }

    /// <summary>
    ///     Get detailed status of all libraries
    /// </summary>
    public static LibraryStatus GetDetailedStatus()
    {
        var (cudaAvail, cudaVer, cudaErr) = GetCudaVersion();

        return new LibraryStatus
        {
            CudaAvailable = cudaAvail,
            CudaVersion = cudaVer,
            CudaErrorCode = cudaErr,
            CuSparseAvailable = IsCuSparseAvailable(),
            MKLAvailable = IsMKLAvailable(),
            PardisoAvailable = IsPardisoAvailable()
        };
    }

    /// <summary>
    ///     Print a formatted status report to console
    /// </summary>
    public static void PrintStatus()
    {
        var status = GetDetailedStatus();

        Console.WriteLine("═══════════════════════════════════════════════════════════");
        Console.WriteLine("  Runtime Library Availability Check");
        Console.WriteLine("═══════════════════════════════════════════════════════════");
        Console.WriteLine();

        // Platform info
        Console.WriteLine($"Platform:      {GetPlatformName()}");
        Console.WriteLine();

        // CUDA
        Console.Write("CUDA Runtime:  ");
        if (status.CudaAvailable)
        {
            var major = status.CudaVersion / 1000;
            var minor = status.CudaVersion % 1000 / 10;
            Console.Write($"✓ AVAILABLE (v{major}.{minor})");

            if (status.CudaErrorCode == 0)
                Console.WriteLine(" - Fully functional");
            else if (status.CudaErrorCode == 100)
                Console.WriteLine(" - No GPU detected");
            else if (status.CudaErrorCode == 35)
                Console.WriteLine(" - Driver needs update");
            else
                Console.WriteLine($" - Error code {status.CudaErrorCode}");
        }
        else
        {
            Console.WriteLine("✗ NOT AVAILABLE");
        }

        // cuSPARSE
        Console.Write("cuSPARSE:      ");
        Console.WriteLine(status.CuSparseAvailable ? "✓ AVAILABLE" : "✗ NOT AVAILABLE");

        // MKL
        Console.Write("Intel MKL:     ");
        Console.WriteLine(status.MKLAvailable ? "✓ AVAILABLE" : "✗ NOT AVAILABLE");

        // PARDISO
        Console.Write("PARDISO:       ");
        Console.WriteLine(status.PardisoAvailable ? "✓ AVAILABLE" : "✗ NOT AVAILABLE");

        Console.WriteLine();

        // Performance summary
        if (status.CudaAvailable && status.CuSparseAvailable)
            Console.WriteLine("GPU Acceleration: ENABLED");
        else
            Console.WriteLine("GPU Acceleration: DISABLED");

        if (status.MKLAvailable)
            Console.WriteLine("CPU Acceleration: ENABLED (Intel MKL)");
        else
            Console.WriteLine("CPU Acceleration: DISABLED");

        Console.WriteLine();
        Console.WriteLine("═══════════════════════════════════════════════════════════");
    }

    #endregion

    #region Helper Methods

    private static string GetPlatformName()
    {
        if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            return "Windows";
        if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
            return "Linux";
        if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
            return "macOS";
        return "Unknown";
    }

    private static bool CheckLibrary(IEnumerable<string> paths, string requiredFunction)
    {
        return CheckLibrary(paths, new[] { requiredFunction });
    }

    private static bool CheckLibrary(IEnumerable<string> paths, string[] requiredFunctions)
    {
        foreach (var path in paths)
        {
            if (!File.Exists(path)) continue;

            try
            {
                if (NativeLibrary.TryLoad(path, out var handle))
                {
                    var hasFunction = requiredFunctions.Any(func =>
                        NativeLibrary.TryGetExport(handle, func, out _));

                    NativeLibrary.Free(handle);

                    if (hasFunction) return true;
                }
            }
            catch
            {
            }
        }

        return false;
    }

    /// <summary>
    ///     Get environment variable value or null if not set
    /// </summary>
    private static string? GetEnvVar(string name)
    {
        return Environment.GetEnvironmentVariable(name);
    }

    /// <summary>
    ///     Expand path with environment variable
    /// </summary>
    private static string ExpandPath(string path)
    {
        return Environment.ExpandEnvironmentVariables(path);
    }

    #endregion

    #region Path Discovery (Simple hardcoded paths for quick checks)

    private static IEnumerable<string> GetCudaRuntimePaths()
    {
        var paths = new List<string>();

        if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
        {
            var cudaPath = GetEnvVar("CUDA_PATH");
            if (!string.IsNullOrEmpty(cudaPath))
            {
                paths.Add(Path.Combine(cudaPath, "lib64", "libcudart.so"));
                paths.Add(Path.Combine(cudaPath, "lib", "libcudart.so"));
            }

            paths.AddRange(new[]
            {
                "/usr/lib/x86_64-linux-gnu/libcudart.so",
                "/usr/lib/x86_64-linux-gnu/libcudart.so.12",
                "/usr/lib/x86_64-linux-gnu/libcudart.so.11",
                "/usr/lib64/libcudart.so",
                "/usr/lib64/libcudart.so.12",
                "/usr/lib64/libcudart.so.11",
                "/usr/local/cuda/lib64/libcudart.so",
                "/usr/local/cuda-12/lib64/libcudart.so",
                "/usr/local/cuda-11/lib64/libcudart.so",
                "/opt/cuda/lib64/libcudart.so"
            });
        }
        else if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
        {
            var cudaPath = GetEnvVar("CUDA_PATH");
            if (!string.IsNullOrEmpty(cudaPath))
            {
                paths.Add(Path.Combine(cudaPath, "bin", "cudart64_12.dll"));
                paths.Add(Path.Combine(cudaPath, "bin", "cudart64_11.dll"));
            }

            paths.AddRange(new[]
            {
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\bin\cudart64_12.dll",
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.5\bin\cudart64_12.dll",
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.0\bin\cudart64_12.dll",
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin\cudart64_11.dll",
                @"C:\Windows\System32\cudart64_12.dll",
                @"C:\Windows\System32\cudart64_11.dll"
            });
        }

        return paths.Where(p => !string.IsNullOrEmpty(p));
    }

    private static IEnumerable<string> GetCuSparsePaths()
    {
        var paths = new List<string>();

        if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
        {
            var cudaPath = GetEnvVar("CUDA_PATH");
            if (!string.IsNullOrEmpty(cudaPath))
            {
                paths.Add(Path.Combine(cudaPath, "lib64", "libcusparse.so"));
                paths.Add(Path.Combine(cudaPath, "lib", "libcusparse.so"));
            }

            paths.AddRange(new[]
            {
                "/usr/lib/x86_64-linux-gnu/libcusparse.so",
                "/usr/lib/x86_64-linux-gnu/libcusparse.so.12",
                "/usr/lib/x86_64-linux-gnu/libcusparse.so.11",
                "/usr/lib64/libcusparse.so",
                "/usr/lib64/libcusparse.so.12",
                "/usr/local/cuda/lib64/libcusparse.so",
                "/usr/local/cuda-12/lib64/libcusparse.so",
                "/opt/cuda/lib64/libcusparse.so"
            });
        }
        else if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
        {
            var cudaPath = GetEnvVar("CUDA_PATH");
            if (!string.IsNullOrEmpty(cudaPath))
            {
                paths.Add(Path.Combine(cudaPath, "bin", "cusparse64_12.dll"));
                paths.Add(Path.Combine(cudaPath, "bin", "cusparse64_11.dll"));
            }

            paths.AddRange(new[]
            {
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\bin\cusparse64_12.dll",
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.5\bin\cusparse64_12.dll",
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.0\bin\cusparse64_12.dll",
                @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin\cusparse64_11.dll",
                @"C:\Windows\System32\cusparse64_12.dll",
                @"C:\Windows\System32\cusparse64_11.dll"
            });
        }

        return paths.Where(p => !string.IsNullOrEmpty(p));
    }

    private static IEnumerable<string> GetMKLPaths()
    {
        var paths = new List<string>();

        if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
        {
            var mklRoot = GetEnvVar("MKLROOT");
            if (!string.IsNullOrEmpty(mklRoot))
            {
                paths.Add(Path.Combine(mklRoot, "lib", "intel64", "libmkl_rt.so"));
                paths.Add(Path.Combine(mklRoot, "lib", "libmkl_rt.so"));
            }

            paths.AddRange(new[]
            {
                "/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_rt.so",
                "/opt/intel/oneapi/mkl/2025.0/lib/libmkl_rt.so",
                "/opt/intel/oneapi/mkl/2024.2/lib/libmkl_rt.so",
                "/opt/intel/mkl/lib/intel64/libmkl_rt.so",
                "/usr/lib/x86_64-linux-gnu/libmkl_rt.so",
                "/usr/lib64/libmkl_rt.so",
                "/usr/local/lib/libmkl_rt.so",
                ExpandPath("$HOME/anaconda3/lib/libmkl_rt.so"),
                ExpandPath("$HOME/miniconda3/lib/libmkl_rt.so")
            });
        }
        else if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
        {
            var mklRoot = GetEnvVar("MKLROOT");
            if (!string.IsNullOrEmpty(mklRoot))
            {
                paths.Add(Path.Combine(mklRoot, "redist", "intel64", "mkl_rt.2.dll"));
                paths.Add(Path.Combine(mklRoot, "redist", "intel64", "mkl_rt.dll"));
            }

            paths.AddRange(new[]
            {
                @"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\redist\intel64\mkl_rt.2.dll",
                @"C:\Program Files (x86)\Intel\oneAPI\mkl\2025.0\redist\intel64\mkl_rt.2.dll",
                @"C:\Program Files (x86)\Intel\oneAPI\mkl\2024.2\redist\intel64\mkl_rt.2.dll",
                @"C:\Program Files\Intel\oneAPI\mkl\latest\redist\intel64\mkl_rt.2.dll",
                @"C:\Program Files (x86)\Intel\MKL\redist\intel64\mkl_rt.dll",
                ExpandPath(@"%USERPROFILE%\anaconda3\Library\bin\mkl_rt.2.dll"),
                ExpandPath(@"%USERPROFILE%\miniconda3\Library\bin\mkl_rt.2.dll")
            });
        }
        else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
        {
            var mklRoot = GetEnvVar("MKLROOT");
            if (!string.IsNullOrEmpty(mklRoot)) paths.Add(Path.Combine(mklRoot, "lib", "libmkl_rt.dylib"));

            paths.AddRange(new[]
            {
                "/opt/intel/oneapi/mkl/latest/lib/libmkl_rt.dylib",
                "/opt/intel/mkl/lib/libmkl_rt.dylib",
                "/opt/homebrew/opt/intel-oneapi-mkl/lib/libmkl_rt.dylib",
                "/opt/homebrew/lib/libmkl_rt.dylib",
                "/usr/local/opt/intel-oneapi-mkl/lib/libmkl_rt.dylib",
                "/usr/local/lib/libmkl_rt.dylib",
                ExpandPath("$HOME/anaconda3/lib/libmkl_rt.dylib"),
                ExpandPath("$HOME/miniconda3/lib/libmkl_rt.dylib")
            });
        }

        return paths.Where(p => !string.IsNullOrEmpty(p) && !p.Contains("$"));
    }

    #endregion
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
                var forceVersion = Environment.GetEnvironmentVariable("CUDA_FORCE_VERSION");
                var forcePath = Environment.GetEnvironmentVariable("CUDA_FORCE_PATH");

                if (!string.IsNullOrEmpty(forcePath) && File.Exists(forcePath))
                {
                    Debug.WriteLine($"[NativeLibrary] Using CUDA_FORCE_PATH: {forcePath}");
                    _cachedCudaRt = new[] { forcePath };
                    return _cachedCudaRt;
                }

                Debug.WriteLine("[NativeLibrary] Starting CUDA Runtime discovery...");

                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                    _cachedCudaRt = DiscoverLibrariesInSystem(LibraryPatterns.CudaRuntimeWindows, "CUDA Runtime");
                else if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                    _cachedCudaRt = DiscoverLibrariesInSystem(LibraryPatterns.CudaRuntimeLinux, "CUDA Runtime");

                if (_cachedCudaRt == null || _cachedCudaRt.Length == 0)
                    Debug.WriteLine("[NativeLibrary] ✗ CUDA Runtime discovery returned NOTHING!");
                else
                    Debug.WriteLine($"[NativeLibrary] ✓ CUDA Runtime discovery found {_cachedCudaRt.Length} libraries");

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

            // Check if discovered libraries actually work
            var shouldTryAutoInstall = false;

            if (_cachedMKL != null && _cachedMKL.Length > 0)
            {
                var searchPaths = GetSearchPaths();
                if (searchPaths != null)
                {
                    var anyWorking = false;
                    foreach (var libName in _cachedMKL.Take(3))
                    {
                        foreach (var searchPath in searchPaths.Take(10))
                            try
                            {
                                var fullPath = Path.Combine(searchPath, libName);
                                if (File.Exists(fullPath))
                                    if (NativeLibrary.TryLoad(fullPath, out var testHandle))
                                    {
                                        var works = NativeLibrary.TryGetExport(testHandle, "mkl_get_version", out _) ||
                                                    NativeLibrary.TryGetExport(testHandle, "MKL_Get_Max_Threads",
                                                        out _);
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
                    Debug.WriteLine("[NativeLibrary] MKL not found via scanning, attempting automatic installation...");

                if (TryAutoInstallMKL())
                {
                    Debug.WriteLine("[NativeLibrary] Re-scanning for MKL after installation...");
                    _cachedSearchPaths = null;

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
            if (_cachedSearchPaths == null)
            {
                var paths = BuildComprehensiveSearchPaths();

                // SAFETY: Ensure critical paths are always included on Linux
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Linux))
                {
                    var pathList = new List<string>(paths);
                    var criticalPaths = new[]
                    {
                        "/usr/lib/x86_64-linux-gnu",
                        "/usr/lib64",
                        "/usr/lib"
                    };

                    foreach (var critical in criticalPaths)
                        if (Directory.Exists(critical) && !pathList.Contains(critical))
                        {
                            Debug.WriteLine($"[NativeLibrary] SAFETY: Adding missing critical path: {critical}");
                            pathList.Insert(0, critical);
                        }

                    _cachedSearchPaths = pathList.ToArray();
                }
                else
                {
                    _cachedSearchPaths = paths;
                }
            }

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

        // PRIORITY 1: Application base directory
        try
        {
            AddPath(AppDomain.CurrentDomain.BaseDirectory, "app base");
        }
        catch
        {
        }

        // PRIORITY 2: Assembly location
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

        // PRIORITY 6: System directory scanning
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

        // PRIORITY 8: Current working directory (lowest priority to avoid DLL hijacking)
        try
        {
            AddPath(Directory.GetCurrentDirectory(), "pwd");
        }
        catch
        {
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
        // Standard Windows system paths
        addPath(@"C:\Windows\System32", "Windows System32");
        addPath(@"C:\Windows\SysWOW64", "Windows SysWOW64 (32-bit)");

        // Program Files variants
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

        programFilesPaths.Add(@"C:\Program Files");
        programFilesPaths.Add(@"C:\Program Files (x86)");

        foreach (var programFiles in programFilesPaths)
        {
            if (string.IsNullOrEmpty(programFiles) || !Directory.Exists(programFiles))
                continue;

            ScanForCUDAWindows(programFiles, addPath);
            ScanForMKLWindows(programFiles, addPath);

            if (EnableDeepScanning) ScanForLibraryDirectories(programFiles, 3, addPath, "Program Files");
        }

        // CUDA discovery
        var cudaBasePaths = new[]
        {
            @"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA",
            @"C:\Program Files (x86)\NVIDIA GPU Computing Toolkit\CUDA",
            @"C:\NVIDIA\CUDA",
            @"C:\cuda"
        };

        foreach (var basePath in cudaBasePaths)
            if (Directory.Exists(basePath))
                try
                {
                    var versionDirs = Directory.GetDirectories(basePath)
                        .OrderByDescending(d => ExtractVersionFromPath(d));

                    foreach (var versionDir in versionDirs)
                    foreach (var subdir in new[] { "bin", @"lib\x64", "lib", @"bin\x64" })
                    {
                        var fullPath = Path.Combine(versionDir, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, $"CUDA {Path.GetFileName(versionDir)}");
                    }
                }
                catch
                {
                }

        // MKL discovery
        var mklBasePaths = new[]
        {
            @"C:\Program Files\Intel\oneAPI\mkl",
            @"C:\Program Files (x86)\Intel\oneAPI\mkl",
            @"C:\Program Files\Intel\MKL",
            @"C:\Intel\oneAPI\mkl"
        };

        foreach (var basePath in mklBasePaths)
            if (Directory.Exists(basePath))
                try
                {
                    var versionDirs = Directory.GetDirectories(basePath)
                        .OrderByDescending(d => ExtractVersionFromPath(d));

                    foreach (var versionDir in versionDirs)
                    foreach (var subdir in new[] { @"redist\intel64", @"bin\intel64", @"lib\intel64", "bin", "lib" })
                    {
                        var fullPath = Path.Combine(versionDir, subdir);
                        if (Directory.Exists(fullPath))
                            addPath(fullPath, $"MKL {Path.GetFileName(versionDir)}");
                    }
                }
                catch
                {
                }

        // Environment variable paths
        var envVars = new[] { "CUDA_PATH", "CUDA_HOME", "MKL_ROOT", "MKLROOT", "INTEL_ROOT" };

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
            }
        }
        catch
        {
        }
    }

    private void AddLinuxSystemPaths(Action<string, string> addPath)
    {
        // Standard system library paths
        var standardPaths = new[]
        {
            "/lib", "/lib64", "/usr/lib", "/usr/lib64",
            "/usr/local/lib", "/usr/local/lib64"
        };

        foreach (var path in standardPaths)
            if (Directory.Exists(path))
                addPath(path, "Linux system lib");

        // Multiarch paths (Debian/Ubuntu)
        var multarchSuffixes = new[] { "x86_64-linux-gnu", "aarch64-linux-gnu" };

        foreach (var suffix in multarchSuffixes)
        {
            var multiarchPaths = new[] { $"/lib/{suffix}", $"/usr/lib/{suffix}", $"/usr/local/lib/{suffix}" };

            foreach (var path in multiarchPaths)
                if (Directory.Exists(path))
                    addPath(path, $"Multiarch ({suffix})");
        }

        // CUDA discovery
        var cudaBasePaths = new[]
        {
            "/usr/local/cuda", "/opt/cuda", "/usr/cuda", "/opt/nvidia/cuda"
        };

        foreach (var basePath in cudaBasePaths)
            if (Directory.Exists(basePath))
                foreach (var subdir in new[] { "lib64", "lib", "lib/x64" })
                {
                    var fullPath = Path.Combine(basePath, subdir);
                    if (Directory.Exists(fullPath))
                        addPath(fullPath, $"CUDA {basePath}");
                }

        // Versioned CUDA directories
        try
        {
            var searchDirs = new[] { "/usr/local", "/opt", "/usr" };
            foreach (var searchDir in searchDirs)
            {
                if (!Directory.Exists(searchDir)) continue;

                var cudaDirs = Directory.GetDirectories(searchDir, "cuda*", SearchOption.TopDirectoryOnly)
                    .OrderByDescending(d => ExtractVersionFromPath(d));

                foreach (var cudaDir in cudaDirs)
                foreach (var subdir in new[] { "lib64", "lib" })
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

        // MKL discovery
        var mklBasePaths = new[]
        {
            "/opt/intel/oneapi/mkl",
            "/opt/intel/mkl",
            "/opt/intel/compilers_and_libraries/linux/mkl",
            "/usr/local/intel/mkl"
        };

        foreach (var basePath in mklBasePaths)
            if (Directory.Exists(basePath))
                foreach (var subdir in new[] { "latest/lib/intel64", "latest/lib", "lib/intel64", "lib" })
                {
                    var fullPath = Path.Combine(basePath, subdir);
                    if (Directory.Exists(fullPath))
                        addPath(fullPath, "MKL Linux");
                }

        // Deep scanning if enabled
        if (EnableDeepScanning)
        {
            ScanForLibraryDirectories("/opt", 4, addPath, "opt");
            ScanForLibraryDirectories("/usr/local", 3, addPath, "usr/local");
        }

        // Query ldconfig
        QueryLdconfigForLibraries(addPath);

        // Environment variables
        var envVars = new[] { "CUDA_HOME", "CUDA_PATH", "MKL_ROOT", "MKLROOT", "INTEL_ROOT" };

        foreach (var envVar in envVars)
        {
            var value = Environment.GetEnvironmentVariable(envVar);
            if (!string.IsNullOrEmpty(value) && Directory.Exists(value))
                foreach (var subdir in new[] { "lib64", "lib", "lib/intel64" })
                {
                    var fullPath = Path.Combine(value, subdir);
                    if (Directory.Exists(fullPath))
                        addPath(fullPath, $"Env: {envVar}");
                }
        }
    }

    private void AddMacOSSystemPaths(Action<string, string> addPath)
    {
        // Standard macOS paths
        var standardPaths = new[]
        {
            "/usr/lib", "/usr/local/lib", "/opt/local/lib",
            "/opt/homebrew/lib", "/usr/local/opt", "/Library/Frameworks"
        };

        foreach (var path in standardPaths)
            if (Directory.Exists(path))
                addPath(path, "macOS system lib");

        // Homebrew paths
        var homebrewPaths = new[] { "/opt/homebrew", "/usr/local/Homebrew", "/usr/local" };

        foreach (var brewBase in homebrewPaths)
            if (Directory.Exists(brewBase))
                foreach (var subdir in new[] { "lib", "opt", "Cellar" })
                {
                    var fullPath = Path.Combine(brewBase, subdir);
                    if (Directory.Exists(fullPath))
                    {
                        addPath(fullPath, $"Homebrew {subdir}");

                        if (subdir == "Cellar" && EnableDeepScanning)
                            ScanForLibraryDirectories(fullPath, 3, addPath, "Homebrew Cellar");
                    }
                }

        // MKL discovery
        var mklBasePaths = new[]
        {
            "/opt/intel/oneapi/mkl", "/opt/intel/mkl",
            "/usr/local/intel/mkl", "/Library/Frameworks/Intel_MKL.framework"
        };

        foreach (var basePath in mklBasePaths)
            if (Directory.Exists(basePath))
                try
                {
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
                        foreach (var subdir in new[] { "latest/lib", "lib", "Libraries" })
                        {
                            var fullPath = Path.Combine(basePath, subdir);
                            if (Directory.Exists(fullPath))
                                addPath(fullPath, "MKL macOS");
                        }
                }
                catch
                {
                }

        // Environment variables
        var envVars = new[] { "CUDA_HOME", "CUDA_PATH", "MKL_ROOT", "MKLROOT" };

        foreach (var envVar in envVars)
        {
            var value = Environment.GetEnvironmentVariable(envVar);
            if (!string.IsNullOrEmpty(value) && Directory.Exists(value))
                foreach (var subdir in new[] { "lib", "lib64", "lib/intel64" })
                {
                    var fullPath = Path.Combine(value, subdir);
                    if (Directory.Exists(fullPath))
                        addPath(fullPath, $"Env: {envVar}");
                }
        }
    }

    private void ScanForLibraryDirectories(string baseDir, int maxDepth, Action<string, string> addPath, string source)
    {
        if (maxDepth <= 0 || !Directory.Exists(baseDir))
            return;

        try
        {
            var hasLibraries = Directory.GetFiles(baseDir, "*.so*", SearchOption.TopDirectoryOnly).Length > 0 ||
                               Directory.GetFiles(baseDir, "*.a", SearchOption.TopDirectoryOnly).Length > 0;

            if (hasLibraries) addPath(baseDir, source);

            var libDirNames = new[] { "lib", "lib64", "lib32", "libs", "library" };
            foreach (var libDirName in libDirNames)
            {
                var libPath = Path.Combine(baseDir, libDirName);
                if (Directory.Exists(libPath)) addPath(libPath, $"{source}/{libDirName}");
            }

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
            if (!process.WaitForExit(2000))
            {
                try
                {
                    process.Kill();
                }
                catch
                {
                }

                return;
            }

            if (process.ExitCode != 0)
                return;

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
                                var lowerPath = directory.ToLowerInvariant();
                                var source = "ldconfig";

                                if (lowerPath.Contains("cuda"))
                                    source = "ldconfig (CUDA)";
                                else if (lowerPath.Contains("mkl") || lowerPath.Contains("intel"))
                                    source = "ldconfig (MKL)";

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

                var files = Directory.GetFiles(searchPath, "*", SearchOption.TopDirectoryOnly)
                    .Where(f => pattern.IsMatch(Path.GetFileName(f)));

                foreach (var file in files)
                {
                    var fileName = Path.GetFileName(file);
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

                // Handle versioned files without symlinks
                if (libraryType.Contains("CUDA Runtime", StringComparison.OrdinalIgnoreCase))
                    AddVersionedFiles(searchPath, "libcudart.so.*", discoveredLibraries);
                else if (libraryType.Contains("cuSPARSE", StringComparison.OrdinalIgnoreCase))
                    AddVersionedFiles(searchPath, "libcusparse.so.*", discoveredLibraries);
                else if (libraryType.Contains("MKL") || libraryType.Contains("mkl"))
                    AddVersionedFiles(searchPath, "libmkl_rt.so.*", discoveredLibraries);
            }
            catch
            {
            }

        if (discoveredLibraries.Count == 0)
            return null;

        // Sort by resolved version (highest first)
        var sorted = discoveredLibraries
            .OrderByDescending(lib => lib.ResolvedVersion)
            .ThenByDescending(lib => lib.DirectoryVersion)
            .ThenByDescending(lib => lib.FileVersion)
            .ThenBy(lib => lib.IsSymlink ? 0 : 1)
            .ToList();

        var result = sorted.Select(lib => lib.FullPath).Distinct().ToArray();

        Debug.WriteLine($"[NativeLibrary] Discovered {libraryType}: {string.Join(", ", result.Take(5))}");

        return result;
    }

    private void AddVersionedFiles(string searchPath, string pattern, List<LibraryInfo> discoveredLibraries)
    {
        try
        {
            var files = Directory.GetFiles(searchPath, pattern);
            foreach (var file in files)
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
        catch
        {
        }
    }

    private string? ResolveSymlink(string path)
    {
        try
        {
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            {
                var fileInfo = new FileInfo(path);
                if (fileInfo.Attributes.HasFlag(FileAttributes.ReparsePoint))
                {
                    var target = fileInfo.LinkTarget;
                    if (!string.IsNullOrEmpty(target))
                    {
                        var resolvedPath = Path.IsPathRooted(target)
                            ? target
                            : Path.GetFullPath(Path.Combine(fileInfo.DirectoryName ?? "", target));

                        if (File.Exists(resolvedPath))
                            return ResolveSymlink(resolvedPath);
                    }
                }

                return path;
            }
            else
            {
                var fileInfo = new FileInfo(path);

                if (fileInfo.Attributes.HasFlag(FileAttributes.ReparsePoint))
                {
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
                        if (!process.WaitForExit(5000))
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

    private double ExtractVersionFromFileName(string fileName)
    {
        // Extract version from patterns like libcudart.so.12.0 or cudart64_12.dll
        var match = Regex.Match(fileName, @"[\._](\d+(?:\.\d+(?:\.\d+)?)?)");
        if (match.Success && double.TryParse(match.Groups[1].Value.Split('.')[0], out var version))
        {
            // Try to get more precision
            var fullVersion = match.Groups[1].Value;
            var parts = fullVersion.Split('.');
            if (parts.Length >= 2 && double.TryParse($"{parts[0]}.{parts[1]}", out var preciseVersion))
                return preciseVersion;
            return version;
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

    // Library file patterns - matches ANY version
    private static class LibraryPatterns
    {
        public static readonly Regex CudaRuntimeWindows = new(@"^cudart64(_\d+)?\.dll$", RegexOptions.IgnoreCase);

        public static readonly Regex CudaRuntimeLinux =
            new(@"^libcudart\.so(\.\d+(\.\d+)?)?$", RegexOptions.IgnoreCase);

        public static readonly Regex CuSparseWindows = new(@"^cusparse64(_\d+)?\.dll$", RegexOptions.IgnoreCase);
        public static readonly Regex CuSparseLinux = new(@"^libcusparse\.so(\.\d+(\.\d+)?)?$", RegexOptions.IgnoreCase);

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
        public string? ResolvedPath { get; set; }
        public double ResolvedVersion { get; set; }
        public bool IsSymlink { get; set; }
    }

    #region Auto-Install Methods

    /// <summary>
    ///     SECURITY FIX: Chocolatey auto-installation has been disabled.
    /// </summary>
    private bool TryInstallChocolatey()
    {
        Debug.WriteLine("[NativeLibrary] SECURITY: Auto-installation of Chocolatey has been disabled.");
        Debug.WriteLine("[NativeLibrary] To install Chocolatey manually, visit: https://chocolatey.org/install");

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
            Console.WriteLine();
        }

        return false;
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
        var hasDotnet = CommandExists("dotnet");
        var hasChoco = CommandExists("choco");

        if (!hasDotnet && !hasChoco)
        {
            if (InteractiveInstall)
            {
                Console.WriteLine();
                Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                Console.WriteLine("║  Package Manager Not Found                                     ║");
                Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                Console.WriteLine();
                Console.WriteLine("Please install Intel MKL manually:");
                Console.WriteLine("  Option 1: Install .NET SDK, then: dotnet add package Intel.MKL.redist.win");
                Console.WriteLine("  Option 2: Download from intel.com");
                Console.WriteLine();
            }

            return false;
        }

        if (InteractiveInstall)
        {
            Console.WriteLine();
            Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║  Intel MKL Not Found - Auto-Install Available                 ║");
            Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();
            Console.Write("Install Intel MKL now? (y/n): ");

            var response = Console.ReadLine()?.Trim().ToLower();
            if (response != "y" && response != "yes")
                return false;
        }

        // Try NuGet first
        if (hasDotnet)
        {
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
        }

        // Fallback to Chocolatey
        if (hasChoco)
        {
            var command = new[] { "choco", "install", "-y", "intel-mkl" };
            if (ExecutePackageManagerInstall(command))
            {
                if (InteractiveInstall)
                    Console.WriteLine("✓ MKL installed via Chocolatey!");
                return true;
            }
        }

        return false;
    }

    private string? FindProjectFile()
    {
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
            Console.Write("Install Intel MKL now? (y/n): ");

            var response = Console.ReadLine()?.Trim().ToLower();
            if (response != "y" && response != "yes")
                return false;
        }

        return ExecutePackageManagerInstall(installCommand);
    }

    private bool HandleMacOSMKLInstall()
    {
        var hasBrew = CommandExists("brew");

        if (!hasBrew)
        {
            if (InteractiveInstall)
            {
                Console.WriteLine();
                Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
                Console.WriteLine("║  Homebrew Not Found                                            ║");
                Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
                Console.WriteLine();
                Console.WriteLine("To install MKL, please install Homebrew first: https://brew.sh");
                Console.WriteLine("Then run: brew install intel-mkl");
                Console.WriteLine();
            }

            return false;
        }

        if (InteractiveInstall)
        {
            Console.WriteLine();
            Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
            Console.WriteLine("║  Intel MKL Not Found - Auto-Install Available                 ║");
            Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
            Console.WriteLine();
            Console.Write("Install Intel MKL now? (y/n): ");

            var response = Console.ReadLine()?.Trim().ToLower();
            if (response != "y" && response != "yes")
                return false;
        }

        return ExecutePackageManagerInstall(new[] { "brew", "install", "intel-mkl" });
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

            const int timeoutMs = 15 * 60 * 1000; // 15 minutes
            if (!process.WaitForExit(timeoutMs))
            {
                Debug.WriteLine($"[NativeLibrary] Package manager command timed out after 15 minutes: {command[0]}");
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

    #endregion
}

// ============================================================================
// NATIVE LIBRARY STATUS - Availability Checking and Diagnostics
// ============================================================================

/// <summary>
///     Provides detailed information about native library availability and performance.
///     Use this for comprehensive diagnostics; use <see cref="LibraryAvailability" /> for quick checks.
/// </summary>
public static class NativeLibraryStatus
{
    private static volatile bool _initialized;
    private static readonly object _initLock = new();

    private static readonly DetailedLibraryInfo _cudaInfo = new() { LibraryName = "CUDA Runtime" };
    private static readonly DetailedLibraryInfo _cuSparseInfo = new() { LibraryName = "cuSPARSE" };
    private static readonly DetailedLibraryInfo _mklInfo = new() { LibraryName = "Intel MKL" };
    private static readonly DetailedLibraryInfo _pardisoInfo = new() { LibraryName = "PARDISO" };

    /// <summary>Checks if CUDA runtime is available for GPU acceleration.</summary>
    public static bool IsCudaAvailable
    {
        get
        {
            EnsureInitialized();
            return _cudaInfo.IsAvailable;
        }
    }

    /// <summary>Checks if cuSPARSE is available for GPU sparse matrix operations.</summary>
    public static bool IsCuSparseAvailable
    {
        get
        {
            EnsureInitialized();
            return _cuSparseInfo.IsAvailable;
        }
    }

    /// <summary>Checks if Intel MKL is available for CPU acceleration.</summary>
    public static bool IsMklAvailable
    {
        get
        {
            EnsureInitialized();
            return _mklInfo.IsAvailable;
        }
    }

    /// <summary>Checks if PARDISO solver is available (requires MKL).</summary>
    public static bool IsPardisoAvailable
    {
        get
        {
            EnsureInitialized();
            return _pardisoInfo.IsAvailable;
        }
    }

    /// <summary>Gets detailed information about CUDA availability.</summary>
    public static DetailedLibraryInfo GetCudaInfo()
    {
        EnsureInitialized();
        return _cudaInfo;
    }

    /// <summary>Gets detailed information about cuSPARSE availability.</summary>
    public static DetailedLibraryInfo GetCuSparseInfo()
    {
        EnsureInitialized();
        return _cuSparseInfo;
    }

    /// <summary>Gets detailed information about Intel MKL availability.</summary>
    public static DetailedLibraryInfo GetMklInfo()
    {
        EnsureInitialized();
        return _mklInfo;
    }

    /// <summary>Gets detailed information about PARDISO solver availability.</summary>
    public static DetailedLibraryInfo GetPardisoInfo()
    {
        EnsureInitialized();
        return _pardisoInfo;
    }

    /// <summary>Gets a comprehensive status report of all native libraries.</summary>
    public static string GetStatusReport()
    {
        EnsureInitialized();

        var report = new StringBuilder();

        report.AppendLine("═══════════════════════════════════════════════════════════");
        report.AppendLine("  Native Library Status Report");
        report.AppendLine("═══════════════════════════════════════════════════════════");
        report.AppendLine();

        report.AppendLine("GPU ACCELERATION:");
        AppendLibraryStatus(report, _cudaInfo);
        AppendLibraryStatus(report, _cuSparseInfo);
        report.AppendLine();

        report.AppendLine("CPU ACCELERATION:");
        AppendLibraryStatus(report, _mklInfo);
        AppendLibraryStatus(report, _pardisoInfo);
        report.AppendLine();

        report.AppendLine("PERFORMANCE SUMMARY:");
        if (_cudaInfo.IsAvailable && _cuSparseInfo.IsAvailable)
            report.AppendLine("  ✓ GPU acceleration ENABLED (CUDA)");
        else
            report.AppendLine("  × GPU acceleration DISABLED");

        if (_mklInfo.IsAvailable)
            report.AppendLine("  ✓ CPU acceleration ENABLED (Intel MKL)");
        else
            report.AppendLine("  × CPU acceleration DISABLED");

        if (_pardisoInfo.IsAvailable)
            report.AppendLine("  ✓ Direct solver ENABLED (PARDISO)");
        else
            report.AppendLine("  × Direct solver DISABLED");

        if (!_cudaInfo.IsAvailable && !_mklInfo.IsAvailable)
        {
            report.AppendLine();
            report.AppendLine("WARNING: No hardware acceleration available!");
            report.AppendLine("         Consider installing Intel MKL or CUDA.");
        }

        report.AppendLine("═══════════════════════════════════════════════════════════");

        return report.ToString();
    }

    /// <summary>Prints a status report to the console.</summary>
    public static void PrintStatusReport()
    {
        Console.WriteLine(GetStatusReport());
    }

    /// <summary>Gets a short summary suitable for logging.</summary>
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

    private static void AppendLibraryStatus(StringBuilder sb, DetailedLibraryInfo info)
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

            // Use LibraryAvailability for quick checks
            var status = LibraryAvailability.GetDetailedStatus();

            _cudaInfo.IsAvailable = status.CudaAvailable;
            _cuSparseInfo.IsAvailable = status.CuSparseAvailable;
            _mklInfo.IsAvailable = status.MKLAvailable;
            _pardisoInfo.IsAvailable = status.PardisoAvailable;

            if (status.CudaAvailable)
            {
                _cudaInfo.Version = status.CudaVersionString;
                _cuSparseInfo.Version = status.CudaVersionString;
            }
            else
            {
                _cudaInfo.MissingDependencies = new[] { "CUDA Toolkit not installed" };
                _cuSparseInfo.MissingDependencies = new[] { "cuSPARSE (part of CUDA)" };
            }

            if (!status.MKLAvailable) _mklInfo.MissingDependencies = new[] { "Intel MKL not installed" };

            if (!status.PardisoAvailable) _pardisoInfo.MissingDependencies = new[] { "Intel MKL (includes PARDISO)" };

            _initialized = true;
        }
    }

    /// <summary>
    ///     Information about a specific native library backend.
    /// </summary>
    public class DetailedLibraryInfo
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
/// </summary>
public static class NuGetLibraryChecker
{
    /// <summary>Checks if NuGet-provided MKL libraries are present.</summary>
    public static bool EnsureMklNuGetPackage()
    {
        if (IsMklNuGetPackagePresent())
        {
            Debug.WriteLine("[NuGet] Intel MKL NuGet package detected");
            return true;
        }

        Debug.WriteLine("[NuGet] Intel MKL NuGet package not found");

        if (IsInDevelopmentEnvironment())
        {
            ShowNuGetRestoreGuidance();
            return false;
        }

        Debug.WriteLine("[NuGet] Running in production mode - NuGet packages should be pre-installed");
        return false;
    }

    /// <summary>Checks if MKL NuGet package libraries are present in the expected location.</summary>
    public static bool IsMklNuGetPackagePresent()
    {
        try
        {
            var baseDir = AppDomain.CurrentDomain.BaseDirectory;
            var rid = GetRuntimeIdentifier();
            var nugetPath = Path.Combine(baseDir, "runtimes", rid, "native");

            if (!Directory.Exists(nugetPath))
                return false;

            var mklFiles = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                ? new[] { "mkl_rt.2.dll", "mkl_rt.1.dll", "mkl_rt.dll" }
                : new[] { "libmkl_rt.so.2", "libmkl_rt.so.1", "libmkl_rt.so" };

            foreach (var file in mklFiles)
            {
                var fullPath = Path.Combine(nugetPath, file);
                if (File.Exists(fullPath))
                    return true;
            }

            return false;
        }
        catch
        {
            return false;
        }
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

    /// <summary>Checks if running in a development environment.</summary>
    public static bool IsInDevelopmentEnvironment()
    {
        try
        {
            var baseDir = AppDomain.CurrentDomain.BaseDirectory;

            var currentDir = new DirectoryInfo(baseDir);
            while (currentDir != null)
            {
                if (currentDir.GetFiles("*.csproj").Any())
                    return true;
                currentDir = currentDir.Parent;
            }

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

    private static void ShowNuGetRestoreGuidance()
    {
        Console.WriteLine();
        Console.WriteLine("╔════════════════════════════════════════════════════════════════╗");
        Console.WriteLine("║  Intel MKL NuGet Package Not Found                             ║");
        Console.WriteLine("╚════════════════════════════════════════════════════════════════╝");
        Console.WriteLine();
        Console.WriteLine("AUTOMATIC FIX:");
        Console.WriteLine("  Run: dotnet restore");
        Console.WriteLine();
    }

    /// <summary>Attempts to automatically restore NuGet packages.</summary>
    public static bool TryAutoRestoreNuGetPackages()
    {
        try
        {
            if (!IsInDevelopmentEnvironment())
                return false;

            var projectFile = FindProjectFile();
            if (projectFile == null)
                return false;

            Console.WriteLine("Attempting automatic NuGet package restore...");

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

            const int timeoutMs = 5 * 60 * 1000;
            if (!process.WaitForExit(timeoutMs))
            {
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

            if (process.ExitCode == 0)
            {
                Console.WriteLine("✓ NuGet packages restored successfully!");
                return true;
            }

            return false;
        }
        catch (Exception ex)
        {
            Debug.WriteLine($"[NuGet] Auto-restore error: {ex.Message}");
            return false;
        }
    }

    private static string? FindProjectFile()
    {
        try
        {
            var baseDir = new DirectoryInfo(AppDomain.CurrentDomain.BaseDirectory);
            var currentDir = baseDir;

            for (var i = 0; i < 5 && currentDir != null; i++)
            {
                var projectFiles = currentDir.GetFiles("*.csproj");
                if (projectFiles.Length > 0)
                    return projectFiles[0].FullName;

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

// ============================================================================
// ROBUST NATIVE LIBRARY LOADER - Safe Library Loading with Fallback
// ============================================================================

/// <summary>
///     Robust native library loading with multiple path attempts and detailed diagnostics.
/// </summary>
public static class RobustNativeLibraryLoader
{
    /// <summary>
    ///     Tries to load a native library from multiple names and search paths.
    /// </summary>
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

        var attemptedPaths = new List<string>();
        var foundButAllFailed = false;

        foreach (var libName in libraryNames)
        {
            // Strategy 1: System search
            try
            {
                if (NativeLibrary.TryLoad(libName, out var handle))
                {
                    foundButAllFailed = true;

                    if (VerifyLibraryWorks(handle, libraryType))
                    {
                        Debug.WriteLine($"[NativeLibrary] ✓ Loaded {libraryType}: {libName} (system search)");
                        return handle;
                    }

                    NativeLibrary.Free(handle);
                }
            }
            catch
            {
            }

            // Strategy 2: Explicit full path
            foreach (var searchPath in searchPaths)
                try
                {
                    var fullPath = Path.Combine(searchPath, libName);
                    if (File.Exists(fullPath) && !attemptedPaths.Contains(fullPath))
                    {
                        attemptedPaths.Add(fullPath);
                        foundButAllFailed = true;

                        if (NativeLibrary.TryLoad(fullPath, out var handle))
                        {
                            if (VerifyLibraryWorks(handle, libraryType))
                            {
                                Debug.WriteLine($"[NativeLibrary] ✓ Loaded {libraryType}: {fullPath}");
                                return handle;
                            }

                            NativeLibrary.Free(handle);
                        }
                    }
                }
                catch
                {
                }

            // Strategy 3: Aggressive - look for versioned files
            var baseName = libName.Split('.')[0];
            foreach (var searchPath in searchPaths)
                try
                {
                    if (!Directory.Exists(searchPath))
                        continue;

                    var pattern = baseName + ".so*";
                    var candidates = Directory.GetFiles(searchPath, pattern)
                        .Where(f => !attemptedPaths.Contains(f))
                        .OrderBy(f => f.Count(c => c == '.'))
                        .ThenByDescending(f => f);

                    foreach (var candidate in candidates)
                    {
                        attemptedPaths.Add(candidate);
                        foundButAllFailed = true;

                        try
                        {
                            if (NativeLibrary.TryLoad(candidate, out var handle))
                            {
                                if (VerifyLibraryWorks(handle, libraryType))
                                {
                                    Debug.WriteLine(
                                        $"[NativeLibrary] ✓ Loaded {libraryType}: {candidate} (aggressive search)");
                                    return handle;
                                }

                                NativeLibrary.Free(handle);
                            }
                        }
                        catch
                        {
                        }
                    }
                }
                catch
                {
                }
        }

        if (foundButAllFailed)
        {
            Debug.WriteLine($"[NativeLibrary] ✗ Failed to load {libraryType} from any location");
            Debug.WriteLine("[NativeLibrary]   Note: Libraries were found but all failed verification");
        }
        else
        {
            Debug.WriteLine($"[NativeLibrary] ✗ Failed to load {libraryType} from any location");
        }

        return IntPtr.Zero;
    }

    /// <summary>
    ///     Verifies that a loaded library actually works by checking symbols exist.
    /// </summary>
    private static bool VerifyLibraryWorks(IntPtr handle, string libraryType)
    {
        if (handle == IntPtr.Zero)
            return false;

        try
        {
            if (libraryType.Contains("CUDA", StringComparison.OrdinalIgnoreCase))
            {
                if (NativeLibrary.TryGetExport(handle, "cudaRuntimeGetVersion", out var funcPtr))
                    try
                    {
                        var cudaRuntimeGetVersion =
                            Marshal.GetDelegateForFunctionPointer<CudaGetVersionDelegate>(funcPtr);
                        var version = 0;
                        var result = cudaRuntimeGetVersion(ref version);

                        // Only accept success (0) - error codes mean CUDA cannot work
                        if (result == 0)
                            return true;
                    }
                    catch
                    {
                    }

                // Fallback: check for basic CUDA symbols
                if (NativeLibrary.TryGetExport(handle, "cudaGetDeviceCount", out _) ||
                    NativeLibrary.TryGetExport(handle, "cudaMalloc", out _))
                    return true;

                return false;
            }

            if (libraryType.Contains("cuSPARSE", StringComparison.OrdinalIgnoreCase))
            {
                if (NativeLibrary.TryGetExport(handle, "cusparseCreate", out _) ||
                    NativeLibrary.TryGetExport(handle, "cusparseCreateCsr", out _) ||
                    NativeLibrary.TryGetExport(handle, "cusparseSpMV", out _))
                    return true;
                return false;
            }

            if (libraryType.Contains("MKL", StringComparison.OrdinalIgnoreCase))
            {
                if (NativeLibrary.TryGetExport(handle, "mkl_get_version", out _) ||
                    NativeLibrary.TryGetExport(handle, "MKL_Get_Max_Threads", out _) ||
                    NativeLibrary.TryGetExport(handle, "mkl_get_max_threads", out _))
                    return true;
                return false;
            }

            if (libraryType.Contains("PARDISO", StringComparison.OrdinalIgnoreCase))
            {
                if (NativeLibrary.TryGetExport(handle, "pardiso", out _) ||
                    NativeLibrary.TryGetExport(handle, "pardiso_", out _) ||
                    NativeLibrary.TryGetExport(handle, "PARDISO", out _) ||
                    NativeLibrary.TryGetExport(handle, "pardisoinit", out _) ||
                    NativeLibrary.TryGetExport(handle, "pardiso_64", out _))
                    return true;
                return false;
            }

            // Unknown library type - if we can load it, assume it works
            return true;
        }
        catch
        {
            return false;
        }
    }
}
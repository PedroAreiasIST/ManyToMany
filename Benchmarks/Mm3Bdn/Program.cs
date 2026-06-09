// Mm3Bdn/Program.cs
// A standalone BenchmarkDotNet harness for MM3's connectivity-query path. It
// measures exactly what the cross-library Python harness times for MM3 --
// build the topology, then traverse every cell's face-neighbours and every
// vertex's incident cells once -- but with BenchmarkDotNet's rigorous warm-up,
// multi-iteration statistics and allocation tracking, so the numbers can be
// independently audited. Meshes are generated in-process and match the Python
// generators, so nothing external (Python, OpenMesh, VTK, ...) is needed.
using System;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Running;
using Numerical;   // O2M, Topology<>, TypeMap<>, Node, Tri3, Tet4

namespace Mm3Bdn;

/// <summary>Synthetic mesh generators, mirroring mesh_connectivity_benchmark.py.</summary>
internal static class MeshGen
{
    /// <summary>Triangulated n×n grid hitting ~targetFaces triangles (matches _grid_triangles).</summary>
    public static (int nVerts, int[][] cells) GridTriangles(int targetFaces)
    {
        int n = Math.Max(2, (int)Math.Round(Math.Sqrt(targetFaces / 2.0)) + 1);
        int nv = n * n;
        var cells = new int[2 * (n - 1) * (n - 1)][];
        int t = 0;
        for (int r = 0; r < n - 1; r++)
        for (int c = 0; c < n - 1; c++)
        {
            int v00 = r * n + c, v10 = r * n + c + 1, v01 = (r + 1) * n + c, v11 = (r + 1) * n + c + 1;
            cells[t++] = new[] { v00, v10, v11 };
            cells[t++] = new[] { v00, v11, v01 };
        }
        return (nv, cells);
    }

    /// <summary>Structured box, 6 tets per hex sharing the main diagonal (matches _box_tets).</summary>
    public static (int nVerts, int[][] cells) BoxTets(int targetTets)
    {
        int cellsCount = Math.Max(1, targetTets / 6);
        int n = Math.Max(1, (int)Math.Round(Math.Pow(cellsCount, 1.0 / 3.0)));
        int nx = n, ny = n, nz = Math.Max(1, cellsCount / (nx * ny));
        int NX = nx + 1, NY = ny + 1, NZ = nz + 1;
        int Vid(int i, int j, int k) => i + NX * j + NX * NY * k;   // Fortran order, matches numpy order='F'
        int[][] hexsplit =
        {
            new[] { 0, 1, 3, 7 }, new[] { 0, 3, 2, 7 }, new[] { 0, 2, 6, 7 },
            new[] { 0, 6, 4, 7 }, new[] { 0, 4, 5, 7 }, new[] { 0, 5, 1, 7 },
        };
        var cells = new int[nx * ny * nz * 6][];
        int t = 0;
        for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
        {
            int[] cor =
            {
                Vid(i, j, k), Vid(i + 1, j, k), Vid(i, j + 1, k), Vid(i + 1, j + 1, k),
                Vid(i, j, k + 1), Vid(i + 1, j, k + 1), Vid(i, j + 1, k + 1), Vid(i + 1, j + 1, k + 1),
            };
            foreach (var h in hexsplit)
                cells[t++] = new[] { cor[h[0]], cor[h[1]], cor[h[2]], cor[h[3]] };
        }
        return (NX * NY * NZ, cells);
    }
}

[MemoryDiagnoser]
[SimpleJob(warmupCount: 3, iterationCount: 7)]
public class TriangleConnectivity
{
    // Face counts matching the README columns (grid-2K / 100K / 500K).
    [Params(2_000, 100_000, 500_000)] public int Faces;

    private int[][] _cells = null!;
    private int _nVerts;
    private Topology<TypeMap<Node, Tri3>> _prebuilt = null!;

    [GlobalSetup]
    public void Setup()
    {
        (_nVerts, _cells) = MeshGen.GridTriangles(Faces);
        _prebuilt = Build();   // reused by QueryOnly
    }

    private Topology<TypeMap<Node, Tri3>> Build()
    {
        var topo = new Topology<TypeMap<Node, Tri3>>();
        topo.AddRange<Node, int>(new ReadOnlySpan<int>(new int[_nVerts]));   // nodes 0..n-1, one lock
        topo.AddRangeParallel<Tri3, Node>(_cells, 4096);                     // bulk parallel element add
        return topo;
    }

    private static long Traverse(Topology<TypeMap<Node, Tri3>> topo)
    {
        long s = 0;
        O2M ee = topo.GetElementToElementGraph<Tri3, Node>();   // element -> edge/node-adjacent elements
        O2M ve = topo.GetTranspose<Tri3, Node>();               // node -> incident elements
        for (int e = 0; e < ee.Count; e++) s += ee[e].Length;
        for (int v = 0; v < ve.Count; v++) s += ve[v].Length;
        return s;
    }

    [Benchmark(Description = "build + query")]
    public long BuildAndQuery() => Traverse(Build());

    [Benchmark(Description = "query only (adjacency on prebuilt topology)")]
    public long QueryOnly() => Traverse(_prebuilt);
}

[MemoryDiagnoser]
[SimpleJob(warmupCount: 3, iterationCount: 7)]
public class TetConnectivity
{
    // Tet counts matching the README column (box-6K) plus larger volume meshes.
    [Params(6_000, 200_000, 1_200_000)] public int Tets;

    private int[][] _cells = null!;
    private int _nVerts;
    private Topology<TypeMap<Node, Tet4>> _prebuilt = null!;

    [GlobalSetup]
    public void Setup()
    {
        (_nVerts, _cells) = MeshGen.BoxTets(Tets);
        _prebuilt = Build();
    }

    private Topology<TypeMap<Node, Tet4>> Build()
    {
        var topo = new Topology<TypeMap<Node, Tet4>>();
        topo.AddRange<Node, int>(new ReadOnlySpan<int>(new int[_nVerts]));
        topo.AddRangeParallel<Tet4, Node>(_cells, 4096);
        return topo;
    }

    private static long Traverse(Topology<TypeMap<Node, Tet4>> topo)
    {
        long s = 0;
        O2M ee = topo.GetElementToElementGraph<Tet4, Node>();
        O2M ve = topo.GetTranspose<Tet4, Node>();
        for (int e = 0; e < ee.Count; e++) s += ee[e].Length;
        for (int v = 0; v < ve.Count; v++) s += ve[v].Length;
        return s;
    }

    [Benchmark(Description = "build + query")]
    public long BuildAndQuery() => Traverse(Build());

    [Benchmark(Description = "query only (adjacency on prebuilt topology)")]
    public long QueryOnly() => Traverse(_prebuilt);
}

internal static class Program
{
    private static void Main(string[] args)
    {
        // Default to running every benchmark non-interactively.
        if (args is null || args.Length == 0)
            args = new[] { "--filter", "*" };
        BenchmarkSwitcher.FromTypes(new[] { typeof(TriangleConnectivity), typeof(TetConnectivity) }).Run(args);
    }
}

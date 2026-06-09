// Mm3Bench/Program.cs
// Benchmark entry point used by ../mesh_connectivity_benchmark.py to time the
// ManyToMany (MM3) library on the same mesh as the other libraries.
//
// Contract:  dotnet run -c Release --project Benchmarks/Mm3Bench -- <meshfile> <tri|tet> [repeats]
//            prints, as the final stdout line:  BUILD_MS QUERY_MS
//
// Methodology (kept fair against the native/array libraries the harness also times):
//   * The mesh file is read ONCE; file I/O is never inside the timed region.
//   * MM3 is the only JIT-compiled entrant. To compare it against the always
//     AOT-compiled native libraries (OpenMesh/CGAL/VTK) and the compiled libigl
//     wheel on equal footing, it is run through one DISCARDED warm-up iteration so
//     the recorded samples execute fully-optimized code, then `repeats` timed
//     iterations; the reported BUILD_MS / QUERY_MS are the MEDIAN of those. This is
//     the same steady-state regime in which the in-process Python libraries are
//     measured (their median of `repeats` warm calls).
//   * Topology construction uses the bulk AddRangeParallel API -- the realistic way
//     to build a large mesh -- instead of one locked Add call per element.
//
// Mesh formats (written by the Python side):
//   .obj      lines "v x y z" and "f i j k" (1-based)             [tri]
//   .tetmesh  "<nV> <nT>" then nV "x y z" then nT "i j k l" (0-based) [tet]
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Numerical;   // O2M, Topology<>, TypeMap<>, Node, Tri3, Tet4 all live here

internal static class Program
{
    private static int Main(string[] args)
    {
        if (args.Length < 2) { Console.Error.WriteLine("usage: <meshfile> <tri|tet> [repeats]"); return 2; }
        string path = args[0];
        string kind = args[1].ToLowerInvariant();
        int repeats = args.Length > 2 ? int.Parse(args[2], CultureInfo.InvariantCulture) : 5;
        if (repeats < 1) repeats = 1;

        // Read the mesh ONCE. Connectivity is all the topology needs; vertex count
        // is enough to register the nodes (coordinates do not affect the queries).
        int nVerts;
        int[][] cells;
        if (kind == "tri") ReadObj(path, out nVerts, out cells);
        else ReadTet(path, out nVerts, out cells);

        var builds = new double[repeats];
        var queries = new double[repeats];
        // r == -1 is a discarded warm-up so every recorded sample runs JIT-warmed.
        for (int r = -1; r < repeats; r++)
        {
            double b, q;
            long check = kind == "tri"
                ? BuildAndQueryTri(cells, nVerts, out b, out q)
                : BuildAndQueryTet(cells, nVerts, out b, out q);
            if (r >= 0) { builds[r] = b; queries[r] = q; }
            if (check < 0) Console.Error.WriteLine("unreachable"); // keep result live
        }

        Array.Sort(builds); Array.Sort(queries);
        Console.WriteLine($"{builds[repeats / 2].ToString("F3", CultureInfo.InvariantCulture)} " +
                          $"{queries[repeats / 2].ToString("F3", CultureInfo.InvariantCulture)}");
        return 0;
    }

    private static long BuildAndQueryTri(int[][] cells, int nVerts, out double buildMs, out double queryMs)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var topo = new Topology<TypeMap<Node, Tri3>>();
        topo.AddRange<Node, int>(new ReadOnlySpan<int>(new int[nVerts]));  // nodes 0..nVerts-1, one lock
        topo.AddRangeParallel<Tri3, Node>(cells, 4096);                    // bulk parallel element add
        sw.Stop(); buildMs = sw.Elapsed.TotalMilliseconds;

        sw.Restart();
        long s = 0;
        O2M ee = topo.GetElementToElementGraph<Tri3, Node>();   // element -> edge/node-adjacent elements
        O2M ve = topo.GetTranspose<Tri3, Node>();               // node -> incident elements
        for (int e = 0; e < ee.Count; e++) s += ee[e].Length;   // O2M row indexer => ReadOnlySpan<int>
        for (int v = 0; v < ve.Count; v++) s += ve[v].Length;
        sw.Stop(); queryMs = sw.Elapsed.TotalMilliseconds;

        // Self-check: the bulk path must build exactly the requested mesh.
        if (ee.Count != cells.Length || ve.Count != nVerts)
            throw new InvalidOperationException(
                $"tri build mismatch: elements {ee.Count}/{cells.Length}, nodes {ve.Count}/{nVerts}");
        return s;
    }

    private static long BuildAndQueryTet(int[][] cells, int nVerts, out double buildMs, out double queryMs)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var topo = new Topology<TypeMap<Node, Tet4>>();
        topo.AddRange<Node, int>(new ReadOnlySpan<int>(new int[nVerts]));
        topo.AddRangeParallel<Tet4, Node>(cells, 4096);
        sw.Stop(); buildMs = sw.Elapsed.TotalMilliseconds;

        sw.Restart();
        long s = 0;
        O2M ee = topo.GetElementToElementGraph<Tet4, Node>();
        O2M ve = topo.GetTranspose<Tet4, Node>();
        for (int e = 0; e < ee.Count; e++) s += ee[e].Length;
        for (int v = 0; v < ve.Count; v++) s += ve[v].Length;
        sw.Stop(); queryMs = sw.Elapsed.TotalMilliseconds;

        if (ee.Count != cells.Length || ve.Count != nVerts)
            throw new InvalidOperationException(
                $"tet build mismatch: elements {ee.Count}/{cells.Length}, nodes {ve.Count}/{nVerts}");
        return s;
    }

    private static void ReadObj(string path, out int nVerts, out int[][] cells)
    {
        var verts = 0;
        var cellList = new List<int[]>();
        foreach (var line in File.ReadLines(path))
        {
            if (line.Length < 2) continue;
            if (line[0] == 'v' && line[1] == ' ')
            {
                verts++;
            }
            else if (line[0] == 'f' && line[1] == ' ')
            {
                var p = line.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
                int A(int i) => int.Parse(p[i].Split('/')[0], CultureInfo.InvariantCulture) - 1;
                cellList.Add(new[] { A(1), A(2), A(3) });
            }
        }
        nVerts = verts;
        cells = cellList.ToArray();
    }

    private static void ReadTet(string path, out int nVerts, out int[][] cells)
    {
        using var rd = new StreamReader(path);
        var first = rd.ReadLine()!.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
        int nV = int.Parse(first[0], CultureInfo.InvariantCulture);
        int nT = int.Parse(first[1], CultureInfo.InvariantCulture);
        for (int i = 0; i < nV; i++) rd.ReadLine();   // skip coordinates (not needed for connectivity)
        cells = new int[nT][];
        for (int i = 0; i < nT; i++)
        {
            var p = rd.ReadLine()!.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
            cells[i] = new[] { int.Parse(p[0], CultureInfo.InvariantCulture),
                               int.Parse(p[1], CultureInfo.InvariantCulture),
                               int.Parse(p[2], CultureInfo.InvariantCulture),
                               int.Parse(p[3], CultureInfo.InvariantCulture) };
        }
        nVerts = nV;
    }
}

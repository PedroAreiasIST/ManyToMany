// Mm3Bench/Program.cs
// Benchmark entry point used by ../mesh_connectivity_benchmark.py to time the
// ManyToMany (MM3) library on the same mesh as the other libraries.
//
// Contract:  dotnet run -c Release --project Benchmarks/Mm3Bench -- <meshfile> <tri|tet> <repeats>
//            prints, as the final stdout line:  BUILD_MS QUERY_MS
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
        int repeats = args.Length > 2 ? int.Parse(args[2], CultureInfo.InvariantCulture) : 1;

        var builds = new double[repeats];
        var queries = new double[repeats];
        for (int r = 0; r < repeats; r++)
        {
            if (kind == "tri") RunTri(path, out builds[r], out queries[r]);
            else RunTet(path, out builds[r], out queries[r]);
        }
        Array.Sort(builds); Array.Sort(queries);
        Console.WriteLine($"{builds[builds.Length / 2].ToString("F3", CultureInfo.InvariantCulture)} " +
                          $"{queries[queries.Length / 2].ToString("F3", CultureInfo.InvariantCulture)}");
        return 0;
    }

    private static void RunTri(string path, out double buildMs, out double queryMs)
    {
        ReadObj(path, out var verts, out var cells);
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var topo = new Topology<TypeMap<Node, Tri3>>();
        var vid = new int[verts.Count];
        for (int i = 0; i < verts.Count; i++) vid[i] = topo.Add<Node>();
        foreach (var c in cells) topo.Add<Tri3, Node>(vid[c[0]], vid[c[1]], vid[c[2]]);
        sw.Stop(); buildMs = sw.Elapsed.TotalMilliseconds;

        sw.Restart();
        long s = 0;
        O2M ee = topo.GetElementToElementGraph<Tri3, Node>();   // element -> element (shared node/edge)
        O2M ve = topo.GetTranspose<Tri3, Node>();               // node -> incident elements
        for (int e = 0; e < ee.Count; e++) s += ee[e].Length;   // O2M row indexer => ReadOnlySpan<int>
        for (int v = 0; v < ve.Count; v++) s += ve[v].Length;
        sw.Stop(); queryMs = sw.Elapsed.TotalMilliseconds;
        if (s == long.MinValue) Console.Error.WriteLine("unreachable"); // keep s live
    }

    private static void RunTet(string path, out double buildMs, out double queryMs)
    {
        ReadTet(path, out var verts, out var cells);
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var topo = new Topology<TypeMap<Node, Tet4>>();
        var vid = new int[verts.Count];
        for (int i = 0; i < verts.Count; i++) vid[i] = topo.Add<Node>();
        foreach (var c in cells) topo.Add<Tet4, Node>(vid[c[0]], vid[c[1]], vid[c[2]], vid[c[3]]);
        sw.Stop(); buildMs = sw.Elapsed.TotalMilliseconds;

        sw.Restart();
        long s = 0;
        O2M ee = topo.GetElementToElementGraph<Tet4, Node>();
        O2M ve = topo.GetTranspose<Tet4, Node>();
        for (int e = 0; e < ee.Count; e++) s += ee[e].Length;
        for (int v = 0; v < ve.Count; v++) s += ve[v].Length;
        sw.Stop(); queryMs = sw.Elapsed.TotalMilliseconds;
        if (s == long.MinValue) Console.Error.WriteLine("unreachable");
    }

    private static void ReadObj(string path, out List<double[]> verts, out List<int[]> cells)
    {
        verts = new List<double[]>(); cells = new List<int[]>();
        foreach (var line in File.ReadLines(path))
        {
            if (line.Length < 2) continue;
            if (line[0] == 'v' && line[1] == ' ')
            {
                var p = line.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
                verts.Add(new[] { double.Parse(p[1], CultureInfo.InvariantCulture),
                                  double.Parse(p[2], CultureInfo.InvariantCulture),
                                  double.Parse(p[3], CultureInfo.InvariantCulture) });
            }
            else if (line[0] == 'f' && line[1] == ' ')
            {
                var p = line.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
                int A(int i) => int.Parse(p[i].Split('/')[0], CultureInfo.InvariantCulture) - 1;
                cells.Add(new[] { A(1), A(2), A(3) });
            }
        }
    }

    private static void ReadTet(string path, out List<double[]> verts, out List<int[]> cells)
    {
        using var rd = new StreamReader(path);
        var first = rd.ReadLine()!.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
        int nV = int.Parse(first[0], CultureInfo.InvariantCulture);
        int nT = int.Parse(first[1], CultureInfo.InvariantCulture);
        verts = new List<double[]>(nV); cells = new List<int[]>(nT);
        for (int i = 0; i < nV; i++)
        {
            var p = rd.ReadLine()!.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
            verts.Add(new[] { double.Parse(p[0], CultureInfo.InvariantCulture),
                              double.Parse(p[1], CultureInfo.InvariantCulture),
                              double.Parse(p[2], CultureInfo.InvariantCulture) });
        }
        for (int i = 0; i < nT; i++)
        {
            var p = rd.ReadLine()!.Split((char[])null!, StringSplitOptions.RemoveEmptyEntries);
            cells.Add(new[] { int.Parse(p[0], CultureInfo.InvariantCulture),
                              int.Parse(p[1], CultureInfo.InvariantCulture),
                              int.Parse(p[2], CultureInfo.InvariantCulture),
                              int.Parse(p[3], CultureInfo.InvariantCulture) });
        }
    }
}

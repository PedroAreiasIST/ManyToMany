// FlatBuildVerify.cs — the threaded tuned rebuild against the sequential one.
//
// Copyright (C) 2026 Pedro Miguel de Almeida Areias
//
// This file is part of ManyToMany.
//
// ManyToMany is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ManyToMany is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ManyToMany. If not, see <https://www.gnu.org/licenses/>.
//
//
//   verify-flatbuild n
//
// On a Kuhn box, the same box with its elements in random order and a mixed mesh
// (boundary triangles, bars, point elements), the threaded rebuild must give the
// sequential transpose, positions, pattern, destinations, facets and boundary list
// exactly, a proper colouring, and the same colouring on every run. The colouring
// depends on the element count only: the printed hash must not change with the CPUs
// the process may use (run it under taskset with different CPU sets).

using Step2;

namespace Step2.Bench;

public static class FlatBuildVerify
{
    public static void Run(int n)
    {
        var (tets, _, nN) = Kuhn.Box(n);
        var shuffled = tets.ToArray();
        new Random(n).Shuffle(shuffled);
        var (mixed, _, nM) = CutVerify.MixedBox(n);
        Check($"Kuhn n={n}", ElementSpec.Tetrahedron, nN, tets);
        Check($"Kuhn n={n}, shuffled", ElementSpec.Tetrahedron, nN, shuffled);
        Check($"mixed n={n}", ElementSpec.MixedSimplex, nM, mixed.ToArray());
    }

    private static void Check(string what, ElementSpec spec, int nN, int[][] rows)
    {
        var seq = FlatBuild.Build(spec, nN, rows, parallel: false);
        Proper(what, seq);
        int[]? colour = null;
        for (var run = 0; run < 3; run++)
        {
            var par = FlatBuild.Build(spec, nN, rows, parallel: true);
            Same(what, "TPtr", seq.TPtr, par.TPtr);
            Same(what, "TPairs", seq.TPairs, par.TPairs);
            Same(what, "EP", seq.EP, par.EP);
            Same(what, "PPtr", seq.PPtr, par.PPtr);
            Same(what, "PCol", seq.PCol, par.PCol);
            Same(what, "PCnt", seq.PCnt, par.PCnt);
            Same(what, "Dest", seq.Dest, par.Dest);
            Same(what, "FA", seq.FA, par.FA);
            Same(what, "Boundary", seq.Boundary, par.Boundary);
            Proper(what, par);
            if (colour == null) colour = par.Color;
            else Same(what, $"Color (run {run})", colour, par.Color);
        }

        ulong h = 1469598103934665603UL;
        foreach (var c in colour!)
            unchecked
            {
                h = (h ^ (uint)c) * 1099511628211UL;
            }

        Console.WriteLine($"verify-flatbuild {what}: {rows.Length} elements, {Environment.ProcessorCount} threads: the " +
                          $"threaded rebuild equals the sequential one; colouring proper, identical in three runs, " +
                          $"{colour.Max() + 1} colours (sequential greedy {seq.ColorCount}), hash {h:x16}");
    }

    private static void Same(string what, string field, int[] a, int[] b)
    {
        if (a.Length != b.Length) throw new Exception($"{what}: {field} lengths {a.Length} and {b.Length}");
        for (var i = 0; i < a.Length; i++)
            if (a[i] != b[i])
                throw new Exception($"{what}: {field}[{i}] = {a[i]} (sequential) and {b[i]} (threaded)");
    }

    private static void Proper(string what, FlatStructures s)
    {
        for (var e = 0; e < s.NE; e++)
        {
            if (s.A[e].Length == 0) continue;
            if (s.Color[e] < 0 || s.Color[e] >= s.ColorCount) throw new Exception($"{what}: element {e} colour {s.Color[e]}");
            foreach (var v in s.A[e])
                for (var q = s.TPtr[v]; q < s.TPtr[v + 1]; q++)
                {
                    var g = s.TPairs[2 * q];
                    if (g != e && s.Color[g] == s.Color[e])
                        throw new Exception($"{what}: elements {e} and {g} share node {v} and colour {s.Color[e]}");
                }
        }
    }
}

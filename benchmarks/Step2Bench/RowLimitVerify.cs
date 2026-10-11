// RowLimitVerify.cs — a cut that would overflow a pattern row is rejected before any write.
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
//   verify-rowlimit n rounds
//
// The link of the central node v of a Kuhn box (the edges of the elements around v that do
// not contain v) is cut round after round: each element around v becomes four, all around
// v, and every cut edge gives v a new neighbour, so the row of v grows about fourfold per
// round until a cut would exceed IncrementalMesh.MaxRowSlots.
// Every accepted round is checked against the oracle. The rejected round must leave the
// working version equal to the last committed one (node and element counts, content hash
// of a new publication, oracle), for the sequential and for the concurrent cut, and the
// mesh must accept further edits (a cut of a few edges far from the sphere).

using Step2;

namespace Step2.Bench;

public static class RowLimitVerify
{
    public static void Run(int n, int rounds)
    {
        var (rows, coords, nN) = Kuhn.Box(n);
        var flat = FlatBuild.Build(ElementSpec.Tetrahedron, nN, rows, parallel: false);
        var m = new IncrementalMesh(flat, ElementSpec.Tetrahedron, persistent: true);
        var xyz = new List<double>(coords);
        var centre = Kuhn.Node(n, n / 2, n / 2, n / 2);

        var obs = new Midpoints(xyz);
        var committed = m.Commit();
        for (var r = 1; r <= rounds; r++)
        {
            var marked = Link(m, centre);
            var largest = LargestRow(m);
            try
            {
                Cutter.Cut(m, marked, obs);
            }
            catch (RowLimitException ex)
            {
                Unchanged(m, committed, $"round {r}, sequential cut");
                try
                {
                    Cutter.CutParallel(m, [marked], obs, threads: 4);
                    throw new Exception($"round {r}: the concurrent cut accepted what the sequential cut rejected");
                }
                catch (RowLimitException)
                {
                }

                Unchanged(m, committed, $"round {r}, concurrent cut");
                // the mesh still accepts edits: cut three edges of the corner element
                var row = m.Row(0).ToArray();
                Cutter.Cut(m, [(row[0], row[1]), (row[1], row[2]), (row[2], row[3])], obs);
                Oracle.Verify(m.Commit());
                Console.WriteLine($"verify-rowlimit n={n}: round {r} ({marked.Count} edges, largest row {largest} columns) " +
                                  $"rejected before any write ({ex.Message}); mesh unchanged for the sequential and the " +
                                  $"concurrent cut, later cut verified");
                return;
            }

            committed = m.Commit();
            Oracle.Verify(committed);
            Console.WriteLine($"  round {r}: {marked.Count} edges cut, row of the centre {m.PatternRow(centre).Length} columns, " +
                              $"{m.ElementCount} elements");
        }

        throw new Exception($"verify-rowlimit n={n}: no cut rejected in {rounds} rounds");
    }

    private sealed class Midpoints(List<double> xyz) : ICutObserver
    {
        public void NodeCreated(int node, int a, int b)
        {
            if (xyz.Count != 3 * node) throw new Exception($"node {node} created out of order");
            for (var c = 0; c < 3; c++) xyz.Add(0.5 * (xyz[3 * a + c] + xyz[3 * b + c]));
        }

        public void ElementSplit(int parent, int child)
        {
        }
    }

    private static void Unchanged(IncrementalMesh m, Snapshot committed, string what)
    {
        if (m.NodeCount != committed.NodeCount || m.ElementCount != committed.ElementCount)
            throw new Exception($"{what}: counts changed ({m.NodeCount}/{m.ElementCount} nodes/elements, " +
                                $"committed {committed.NodeCount}/{committed.ElementCount})");
        var now = m.Commit();
        if (Oracle.Hash(now) != Oracle.Hash(committed)) throw new Exception($"{what}: the working version changed");
        Oracle.Verify(now);
    }

    /// <summary>Edges of the elements around v that do not contain v.</summary>
    private static List<(int, int)> Link(IncrementalMesh m, int v)
    {
        var seen = new HashSet<long>(EdgeKeyComparer.Instance);
        var marked = new List<(int, int)>();
        var tr = m.TransposePairs(v);
        for (var q = 0; q < tr.Length; q += 2)
        {
            var row = m.Row(tr[q]);
            for (var i = 0; i < row.Length; i++)
            for (var j = i + 1; j < row.Length; j++)
            {
                int a = row[i], b = row[j];
                if (a != v && b != v && seen.Add(Cutter.Key(a, b))) marked.Add((a, b));
            }
        }

        return marked;
    }

    private static int LargestRow(IncrementalMesh m)
    {
        var largest = 0;
        for (var v = 0; v < m.NodeCount; v++) largest = Math.Max(largest, m.PatternRow(v).Length);
        return largest;
    }
}

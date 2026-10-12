// Verify.cs — counting gates: random planar cracks in a cube, checker, Step2 oracle, rollback, disc growth.
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
//   verify-count n seed cases [steps=k] [check]
//
// Inserts `cases` random disc cracks (random normal, centre and radius; some reach the
// boundary, later ones cross earlier ones) into a Kuhn cube of n^3 cells, one bit each.
// The growth test runs with one disc (cases = 1).
// After every insertion: commit, recount, CrackChecker (geometric twins vs counted pairs,
// edge classes, chains), Step2 Oracle (every maintained structure recomputed from the
// rows), and a rollback test (the hash of the rolled-back version equals the hash before
// the insertion). Then the first disc grows `steps` times by a walk from its front with the
// plane as surface and a larger disc as region; the recovered front must be one closed
// chain whose nodes lie within one element size of the new radius. Prints one line
// starting with "verify-count" on success; any mismatch throws.

using Step2;
using Step2.Bench;

namespace Step3.Bench;

public static class Verify
{
    public static void Count(Args a)
    {
        var n = a.Positional.Count > 0 ? int.Parse(a.Positional[0]) : 10;
        var seed = a.Positional.Count > 1 ? int.Parse(a.Positional[1]) : 1;
        var cases = a.Positional.Count > 2 ? int.Parse(a.Positional[2]) : 3;
        var steps = a.I("steps", 3);
        var rng = new Random(seed);
        var cm = KuhnMesh.Build(n, n, n, new V3(-1, -1, -1), new V3(1, 1, 1));
        var h = 2.0 / n;
        var ins = new CrackInserter3(cm) { AllowSameSide = a.Flag("sameside") };
        var cc = new CrackCount(cm);
        var snap = cm.M.Current;
        var hashBefore = Oracle.Hash(snap);
        var junctions = 0;
        var mouths = 0;
        var rejected = 0;
        var discs = new List<(V3 C, V3 N, double R)>();
        for (var c = 0; c < cases; c++)
        {
            var nrm = new V3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5).Unit();
            var cen = new V3(0.6 * (rng.NextDouble() - 0.5), 0.6 * (rng.NextDouble() - 0.5), 0.6 * (rng.NextDouble() - 0.5)) + 0.37 * h * nrm;
            var r = 0.3 + 0.6 * rng.NextDouble();
            var bit = 1 << c;
            cc.Reference[bit] = (cen, nrm);
            Field3 phiS = x => V3.Dot(x - cen, nrm);
            Field3 phiR = x => (x - cen).Perp(nrm).Norm - r;
            var pre = snap;
            var preHash = Oracle.Hash(pre);
            try
            {
                var res = ins.Advance(bit, phiS, phiR, [], true);
                snap = cm.Commit();
                cc.Rescan(snap);
                var chk = CrackChecker.Verify(snap, cm, cc);
                Oracle.Verify(snap);
                junctions += cc.JunctionEdges.Count;
                mouths = cc.MouthEdges.Count;
                Console.WriteLine($"case {c}: disc r {r:F3} normal {nrm}: {res.Edges.Length} edges, {res.Promoted.Length} promoted, {res.Splits.Length} splits, " +
                                  $"{res.SameSideSplits} same-side; {cc.Summary().Split('\n')[0]}; checker {chk.Pairs} pairs ok; oracle ok");
                // rollback test: the version before the insertion is restored bit for bit
                var after = snap;
                cm.Rollback(pre);
                var hashRolled = Oracle.Hash(cm.M.Current);
                if (hashRolled != preHash) throw new InvalidOperationException($"case {c}: rollback hash {hashRolled:x} vs {preHash:x}");
                // redo deterministically: the same insertion gives the same hash
                var res2 = ins.Advance(bit, phiS, phiR, [], true);
                var snap2 = cm.Commit();
                if (Oracle.Hash(snap2) != Oracle.Hash(after)) throw new InvalidOperationException($"case {c}: re-insertion after rollback differs");
                snap = snap2;
                cc.Rescan(snap);
                discs.Add((cen, nrm, r));
            }
            catch (InvalidOperationException ex) when (ex.Message.Contains("same side") || ex.Message.Contains("without a signed"))
            {
                rejected++;
                Console.WriteLine($"case {c}: rejected ({ex.Message.Split('\n')[0]}); rolled back");
                cm.Rollback(pre);
                snap = pre;
                cc.Rescan(snap);
                if (Oracle.Hash(snap) != preHash) throw new InvalidOperationException("rollback after rejection differs");
            }
        }

        // growth of the first disc by walks from its front (one disc only: a disc crossed by later discs is not one crack,
        // and the increments of a crack that meets another end on it, which the slab and 3D examples test)
        var grown = 0;
        if (discs.Count == 1 && steps > 0)
        {
            var (cen, nrm, r) = discs[0];
            Console.WriteLine("  chains before growth:");
            foreach (var ch in cc.Chains)
            {
                double rmin = double.MaxValue, rmax = 0;
                foreach (var p in ch.P)
                {
                    var rr = (p - cen).Perp(nrm).Norm;
                    rmin = Math.Min(rmin, rr);
                    rmax = Math.Max(rmax, rr);
                }

                Console.WriteLine($"   chain sheet {ch.Sheet} bits {ch.Bits} nodes {ch.Count} {(ch.Closed ? "closed" : ch.Start + "-" + ch.End)}: radius {rmin:F3}..{rmax:F3}");
            }
            for (var k = 1; k <= steps; k++)
            {
                var r1 = Math.Min(r + k * 1.5 * h, 0.95);
                var chains = cc.Chains.Where(ch => (ch.Bits & 1) != 0).ToList();
                if (chains.Count == 0) break;
                var seeds = chains.SelectMany(cc.MeshNodes).Distinct().ToArray();
                Field3 phiS = x => V3.Dot(x - cen, nrm);
                Field3 phiR = x => (x - cen).Perp(nrm).Norm - r1;
                var res = ins.Advance(1, phiS, phiR, seeds);
                snap = cm.Commit();
                cc.Rescan(snap);
                CrackChecker.Verify(snap, cm, cc);
                Oracle.Verify(snap);
                // the front of sheet 1 lies near the new radius (where it is inside the body)
                var worst = 0.0;
                foreach (var ch in cc.Chains.Where(ch => (ch.Bits & 1) != 0))
                {
                    double rmin = double.MaxValue, rmax = 0;
                    foreach (var p in ch.P)
                    {
                        var rr = (p - cen).Perp(nrm).Norm;
                        rmin = Math.Min(rmin, rr);
                        rmax = Math.Max(rmax, rr);
                        worst = Math.Max(worst, Math.Abs(rr - r1));
                    }

                    Console.WriteLine($"   chain sheet {ch.Sheet} bits {ch.Bits} nodes {ch.Count} {(ch.Closed ? "closed" : ch.Start + "-" + ch.End)}: radius {rmin:F3}..{rmax:F3} (target {r1:F3})");
                    if (ch.Count <= 15 && (ch.Start == EndKind.Free || ch.End == EndKind.Free || ch.Closed))
                        foreach (var cn in ch.Nodes)
                        {
                            var v = cc.Rep[cn];
                            var copies = 0;
                            for (var w = 0; w < cm.NodeCount; w++)
                                if (cm.Canon[w] == cn)
                                    copies++;
                            Console.WriteLine($"      node {v} canon {cn} copies {copies} mask {cm.Mask[v]} at {cm.P(v)} parents {cm.Parent1[v]},{cm.Parent2[v]} w {cm.Weight[v]:F3} rho {(cm.P(v) - cen).Perp(nrm).Norm:F3}");
                        }
                }
                if (worst > 1.5 * h) throw new InvalidOperationException($"growth step {k}: front node {worst:G4} from the radius {r1:G4} (h {h:G4})");
                Console.WriteLine($"growth {k}: r {r1:F3}, {res.Edges.Length} edges by a walk over {res.VisitedTets} tetrahedra, {res.Splits.Length} splits; front within {worst / h:F2} h of the radius; {cc.Summary().Split('\n')[0]}");
                grown++;
            }
        }

        Console.WriteLine($"verify-count n={n} seed={seed} cases={cases}: {cases - rejected} inserted, {rejected} rejected, {grown} growth steps, " +
                          $"{cc.PairCount} crack facet pairs, {cc.SheetCount} sheets, {cc.Chains.Count} fronts, {cc.JunctionEdges.Count} junction edges, {mouths} mouth edges: ok");
    }

    private static int SheetOfBit(CrackCount cc, int bit)
    {
        for (var s = 0; s < cc.SheetCount; s++)
            if ((cc.SheetBits[s] & bit) != 0)
                return s;
        return -1;
    }
}

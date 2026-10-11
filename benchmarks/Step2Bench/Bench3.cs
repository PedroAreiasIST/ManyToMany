// Bench3.cs — end-to-end crack step: maintenance vs assembly vs solve.
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
// Each crack step is followed by the assembly of a P1 operator (stiffness plus mass,
// so that the cracked body needs no Dirichlet condition) and a fixed number of
// Jacobi-preconditioned conjugate-gradient iterations. The incremental path
// refreshes the flat copy of the pattern (FlatPattern) and assembles through its
// destination table into the gapped pattern; the MM3 path rebuilds transpose,
// positions, cliques and the node-node pattern, then assembles into a compact CSR
// addressed by the cliques. Assembly and solver are sequential and identical on both
// paths. Also measured on the incremental path: the iterations reading the pattern
// through the snapshot, and two threaded assemblies, by colour classes of the
// maintained colouring and by a gather over the transpose rows (bit-identical to the
// sequential element loop).

using System.Collections.Concurrent;
using Numerical;

namespace Step2.Bench;

internal static partial class Bench
{
    private static void ElementMatrix(MeshState st, ReadOnlySpan<int> row, Span<double> ke)
    {
        Span<double> j = stackalloc double[9];
        for (var c = 0; c < 3; c++)
        {
            j[c * 3 + 0] = st.X(row[1], c) - st.X(row[0], c);
            j[c * 3 + 1] = st.X(row[2], c) - st.X(row[0], c);
            j[c * 3 + 2] = st.X(row[3], c) - st.X(row[0], c);
        }

        var det = j[0] * (j[4] * j[8] - j[5] * j[7]) - j[1] * (j[3] * j[8] - j[5] * j[6]) +
                  j[2] * (j[3] * j[7] - j[4] * j[6]);
        var vol = Math.Abs(det) / 6.0;
        // inverse transpose rows give the gradients of the barycentric coordinates 1..3
        Span<double> inv = stackalloc double[9];
        inv[0] = (j[4] * j[8] - j[5] * j[7]) / det;
        inv[1] = (j[2] * j[7] - j[1] * j[8]) / det;
        inv[2] = (j[1] * j[5] - j[2] * j[4]) / det;
        inv[3] = (j[5] * j[6] - j[3] * j[8]) / det;
        inv[4] = (j[0] * j[8] - j[2] * j[6]) / det;
        inv[5] = (j[2] * j[3] - j[0] * j[5]) / det;
        inv[6] = (j[3] * j[7] - j[4] * j[6]) / det;
        inv[7] = (j[1] * j[6] - j[0] * j[7]) / det;
        inv[8] = (j[0] * j[4] - j[1] * j[3]) / det;
        Span<double> g = stackalloc double[12];
        for (var c = 0; c < 3; c++)
        {
            g[3 + c] = inv[0 * 3 + c];
            g[6 + c] = inv[1 * 3 + c];
            g[9 + c] = inv[2 * 3 + c];
            g[c] = -(g[3 + c] + g[6 + c] + g[9 + c]);
        }

        for (var a = 0; a < 4; a++)
        for (var b = 0; b < 4; b++)
        {
            var k = vol * (g[3 * a] * g[3 * b] + g[3 * a + 1] * g[3 * b + 1] + g[3 * a + 2] * g[3 * b + 2]);
            ke[4 * a + b] = k + vol / 20.0 * (a == b ? 2.0 : 1.0);
        }
    }

    private static void AssembleGapped(Snapshot s, FlatPattern fp, MeshState x, double[] vals)
    {
        Array.Clear(vals, 0, s.TotalSlots);
        Span<double> ke = stackalloc double[16];
        var rowStart = fp.RowStart;
        for (var e = 0; e < s.ElementCount; e++)
        {
            var row = s.Row(e);
            if (row.Length == 0) continue;
            ElementMatrix(x, row, ke);
            var dest = s.Destinations(e);
            for (var i = 0; i < 4; i++)
            {
                var rs = rowStart[row[i]];
                for (var j = 0; j < 4; j++) vals[rs + dest[4 * i + j]] += ke[4 * i + j];
            }
        }
    }

    /// <summary>
    ///     Threaded assembly by colour classes: the elements of one colour share no node, so they
    ///     write disjoint pattern rows. The classes are sorted out of the maintained colours first
    ///     (counting sort, included in the time).
    /// </summary>
    private static void AssembleColour(Snapshot s, FlatPattern fp, MeshState x, double[] vals)
    {
        Array.Clear(vals, 0, s.TotalSlots);
        int nE = s.ElementCount, nc = s.ColorCount;
        var ptr = new int[nc + 1];
        for (var e = 0; e < nE; e++)
            if (s.Color(e) >= 0)
                ptr[s.Color(e) + 1]++;
        for (var c = 0; c < nc; c++) ptr[c + 1] += ptr[c];
        var byColour = new int[ptr[nc]];
        var next = (int[])ptr.Clone();
        for (var e = 0; e < nE; e++)
            if (s.Color(e) >= 0)
                byColour[next[s.Color(e)]++] = e;
        var rowStart = fp.RowStart;
        for (var c = 0; c < nc; c++)
            Parallel.ForEach(Partitioner.Create(ptr[c], ptr[c + 1], 2048), r =>
            {
                Span<double> ke = stackalloc double[16];
                for (var k = r.Item1; k < r.Item2; k++)
                {
                    var e = byColour[k];
                    var row = s.Row(e);
                    ElementMatrix(x, row, ke);
                    var dest = s.Destinations(e);
                    for (var i = 0; i < 4; i++)
                    {
                        var rs = rowStart[row[i]];
                        for (var j = 0; j < 4; j++) vals[rs + dest[4 * i + j]] += ke[4 * i + j];
                    }
                }
            });
    }

    /// <summary>
    ///     Threaded assembly by rows: row v sums the element rows of the pairs (e, p) of T(v) in
    ///     element order, the order of the sequential element loop, so every value equals that
    ///     loop's bit for bit. Each element matrix is evaluated once per node.
    /// </summary>
    private static void AssembleGather(Snapshot s, FlatPattern fp, MeshState x, double[] vals)
    {
        Array.Clear(vals, 0, s.TotalSlots);
        var rowStart = fp.RowStart;
        Parallel.ForEach(Partitioner.Create(0, s.NodeCount, 1024), r =>
        {
            Span<double> ke = stackalloc double[16];
            for (var v = r.Item1; v < r.Item2; v++)
            {
                var tr = s.TransposePairs(v);
                var rs = rowStart[v];
                for (var q = 0; q < tr.Length; q += 2)
                {
                    var e = tr[q];
                    var p = tr[q + 1];
                    ElementMatrix(x, s.Row(e), ke);
                    var dest = s.Destinations(e);
                    for (var j = 0; j < 4; j++) vals[rs + dest[4 * p + j]] += ke[4 * p + j];
                }
            }
        });
    }

    /// <summary>Product with the pattern read from its flat copy.</summary>
    private static void SpmvFlat(FlatPattern fp, double[] vals, double[] xv, double[] y)
    {
        int[] rowStart = fp.RowStart, rowLen = fp.RowLen, col = fp.Col;
        for (var v = 0; v < fp.NodeCount; v++)
        {
            var rs = rowStart[v];
            double acc = 0;
            if (rs >= 0)
                for (int t = rs, end = rs + rowLen[v]; t < end; t++)
                    acc += vals[t] * xv[col[t]];
            y[v] = acc;
        }
    }

    /// <summary>Product with the pattern read through the snapshot (two directory levels).</summary>
    private static void SpmvGapped(Snapshot s, double[] vals, double[] xv, double[] y)
    {
        for (var v = 0; v < s.NodeCount; v++)
        {
            var cols = s.PatternColumns(v);
            if (cols.IsEmpty)
            {
                y[v] = 0;
                continue;
            }

            var rs = s.PatternRowStart(v);
            double acc = 0;
            for (var t = 0; t < cols.Length; t++) acc += vals[rs + t] * xv[cols[t]];
            y[v] = acc;
        }
    }

    private static double[] DiagGapped(Snapshot s, double[] vals)
    {
        var d = new double[s.NodeCount];
        for (var v = 0; v < s.NodeCount; v++)
        {
            var cols = s.PatternColumns(v);
            var rs = s.PatternRowStart(v);
            for (var t = 0; t < cols.Length; t++)
                if (cols[t] == v)
                    d[v] += vals[rs + t];
            if (d[v] == 0) d[v] = 1;
        }

        return d;
    }

    private static void AssembleCsr(int[][] rows, List<List<int>> cliques, int[] rowPtr, MeshState x, double[] vals)
    {
        Array.Clear(vals);
        Span<double> ke = stackalloc double[16];
        for (var e = 0; e < rows.Length; e++)
        {
            var row = rows[e];
            if (row.Length == 0) continue;
            ElementMatrix(x, row, ke);
            var cl = cliques[e];
            for (var i = 0; i < 4; i++)
            {
                var rs = rowPtr[row[i]];
                for (var j = 0; j < 4; j++) vals[rs + cl[j + 4 * i]] += ke[4 * i + j];
            }
        }
    }

    private static void SpmvCsr(int[] rowPtr, int[] col, double[] vals, double[] xv, double[] y)
    {
        for (var v = 0; v < rowPtr.Length - 1; v++)
        {
            double acc = 0;
            for (var t = rowPtr[v]; t < rowPtr[v + 1]; t++) acc += vals[t] * xv[col[t]];
            y[v] = acc;
        }
    }

    /// <summary>Jacobi-preconditioned CG, fixed iteration count; returns the final residual norm.</summary>
    private static double Cg(int nN, Action<double[], double[]> apply, double[] diag, int iters)
    {
        var b = new double[nN];
        for (var i = 0; i < nN; i++) b[i] = 1.0;
        var xv = new double[nN];
        var r = (double[])b.Clone();
        var z = new double[nN];
        for (var i = 0; i < nN; i++) z[i] = r[i] / diag[i];
        var p = (double[])z.Clone();
        var q = new double[nN];
        double rz = 0;
        for (var i = 0; i < nN; i++) rz += r[i] * z[i];
        for (var it = 0; it < iters; it++)
        {
            apply(p, q);
            double pq = 0;
            for (var i = 0; i < nN; i++) pq += p[i] * q[i];
            var alpha = rz / pq;
            for (var i = 0; i < nN; i++)
            {
                xv[i] += alpha * p[i];
                r[i] -= alpha * q[i];
                z[i] = r[i] / diag[i];
            }

            double rz2 = 0;
            for (var i = 0; i < nN; i++) rz2 += r[i] * z[i];
            var beta = rz2 / rz;
            rz = rz2;
            for (var i = 0; i < nN; i++) p[i] = z[i] + beta * p[i];
        }

        double rr = 0;
        for (var i = 0; i < nN; i++) rr += r[i] * r[i];
        return Math.Sqrt(rr);
    }

    private static void EndToEnd(int n, int iters)
    {
        var st = Build(n, persistent: true, mirror: true);
        var mesh = st.Mesh;
        var crack = new CrackPlane(n);
        var every = Math.Max(1, crack.Steps / 8);
        Console.WriteLine($"endtoend n={n}: {mesh.ElementCount:N0} elements, {crack.Steps} crack steps, CG {iters} it.");
        using var w = Csv($"endtoend_{n}.csv",
            "step,t_maint_ms,t_copy_ms,chunks_copied,chunks_total,t_asm_ms,t_asm_gather_ms,t_asm_colour_ms,t_solve_ms," +
            "t_solve_snapshot_ms,res_inc,t_mm3_rebuild_ms,t_mm3_asm_ms,t_mm3_solve_ms,res_mm3");
        var fp = new FlatPattern();
        var tm = new List<double>();
        var tc = new List<double>();
        var ta = new List<double>();
        var tag = new List<double>();
        var tac = new List<double>();
        var ts = new List<double>();
        var tss = new List<double>();
        var copied = new List<double>();
        var tr = new List<double>();
        var tra = new List<double>();
        var trs = new List<double>();
        for (var s = 0; s < crack.Steps; s++)
        {
            FullGC();
            var tMaint = TimeMs(() =>
            {
                crack.Step(st);
                mesh.Commit();
            });
            st.FlushMirror();
            var snap = mesh.Current;
            var tCopy = TimeMs(() => fp.Refresh(snap));
            var vals = new double[snap.TotalSlots];
            FullGC();
            var tAsm = TimeMs(() => AssembleGapped(snap, fp, st, vals));
            var check = new double[snap.TotalSlots];
            FullGC();
            var tAsmGather = TimeMs(() => AssembleGather(snap, fp, st, check));
            if (!check.AsSpan().SequenceEqual(vals)) throw new Exception($"step {s}: gathered assembly differs from the element loop");
            FullGC();
            var tAsmColour = TimeMs(() => AssembleColour(snap, fp, st, check));
            double colourDiff = 0, scale = 0;
            for (var t = 0; t < vals.Length; t++)
            {
                colourDiff = Math.Max(colourDiff, Math.Abs(check[t] - vals[t]));
                scale = Math.Max(scale, Math.Abs(vals[t]));
            }

            if (colourDiff > 1e-12 * scale) throw new Exception($"step {s}: colour-class assembly differs by {colourDiff:E2}");
            check = null;
            var diag = DiagGapped(snap, vals);
            double res = 0, resSnap = 0;
            FullGC();
            var tSolve = TimeMs(() => res = Cg(snap.NodeCount, (xin, yout) => SpmvFlat(fp, vals, xin, yout), diag, iters));
            FullGC();
            var tSolveSnap = TimeMs(() =>
                resSnap = Cg(snap.NodeCount, (xin, yout) => SpmvGapped(snap, vals, xin, yout), diag, iters));
            if (res != resSnap) throw new Exception($"step {s}: flat and snapshot products differ");
            tm.Add(tMaint);
            tc.Add(tCopy);
            ta.Add(tAsm);
            tag.Add(tAsmGather);
            tac.Add(tAsmColour);
            ts.Add(tSolve);
            tss.Add(tSolveSnap);
            copied.Add(fp.ChunksCopied);

            double tR = double.NaN, tRa = double.NaN, tRs = double.NaN, resM = double.NaN;
            if (s % every == 0 || s == crack.Steps - 1)
            {
                var mirror = st.Mirror!;
                object? built = null;
                FullGC();
                tR = TimeMs(() => built = mirror.Rebuild());
                var (_, _, cl, nn) = ((O2M, List<List<int>>, List<List<int>>, O2M))built!;
                var nN = nn.Count;
                var rowPtr = new int[nN + 1];
                for (var v = 0; v < nN; v++) rowPtr[v + 1] = rowPtr[v] + nn[v].Length;
                var col = new int[rowPtr[nN]];
                for (var v = 0; v < nN; v++) nn[v].CopyTo(col.AsSpan(rowPtr[v]));
                var rows = CurrentRows(mesh);
                var cvals = new double[rowPtr[nN]];
                FullGC();
                tRa = TimeMs(() => AssembleCsr(rows, cl, rowPtr, st, cvals));
                var cdiag = new double[nN];
                for (var v = 0; v < nN; v++)
                {
                    for (var t = rowPtr[v]; t < rowPtr[v + 1]; t++)
                        if (col[t] == v) cdiag[v] += cvals[t];
                    if (cdiag[v] == 0) cdiag[v] = 1;
                }

                FullGC();
                tRs = TimeMs(() => resM = Cg(nN, (xin, yout) => SpmvCsr(rowPtr, col, cvals, xin, yout), cdiag, iters));
                tr.Add(tR);
                tra.Add(tRa);
                trs.Add(tRs);
            }

            w.WriteLine($"{s},{tMaint:F4},{tCopy:F4},{fp.ChunksCopied},{fp.ChunksTotal},{tAsm:F3},{tAsmGather:F3}," +
                        $"{tAsmColour:F3},{tSolve:F3},{tSolveSnap:F3},{res:E6},{tR:F3},{tRa:F3},{tRs:F3},{resM:E6}");
        }

        var inc = Median(tm) + Median(tc) + Median(ta) + Median(ts);
        var mm3 = Median(tr) + Median(tra) + Median(trs);
        var line = $"endtoend n={n}: incremental step = maint {Median(tm):F2} + copy {Median(tc):F3} + asm {Median(ta):F1} + " +
                   $"solve {Median(ts):F1} = {inc:F1} ms (maint + copy {100 * (Median(tm) + Median(tc)) / inc:F2} %); " +
                   $"MM3 step = rebuild {Median(tr):F1} + asm {Median(tra):F1} + solve {Median(trs):F1} = {mm3:F1} ms " +
                   $"(rebuild {100 * Median(tr) / mm3:F1} %); end-to-end ratio {mm3 / inc:F2}x; solve flat/CSR " +
                   $"{Median(ts) / Median(trs):F2}, snapshot/CSR {Median(tss) / Median(trs):F2} (snapshot {Median(tss):F1}); " +
                   $"threaded assembly: gather {Median(tag):F1}, colour classes {Median(tac):F1} ({Environment.ProcessorCount} threads); " +
                   $"chunks copied per refresh {Median(copied):F0} of {fp.ChunksTotal}";
        Console.WriteLine(line);
        File.AppendAllText(Path.Combine(DataDir, "endtoend_summary.txt"), line + Environment.NewLine);
    }
}

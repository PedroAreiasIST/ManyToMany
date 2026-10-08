// Bench3.cs — end-to-end crack step: maintenance vs assembly vs solve.
//
// Each crack step is followed by the assembly of a P1 operator (stiffness plus mass,
// so that the cracked body needs no Dirichlet condition) and a fixed number of
// Jacobi-preconditioned conjugate-gradient iterations. The incremental path
// assembles through its destination table into the gapped pattern; the MM3 path
// rebuilds transpose, positions, cliques and the node-node pattern, then assembles
// into a compact CSR addressed by the cliques. The solver is identical.

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

    private static void AssembleGapped(Snapshot s, MeshState x, double[] vals)
    {
        Array.Clear(vals, 0, s.TotalSlots);
        Span<double> ke = stackalloc double[16];
        for (var e = 0; e < s.ElementCount; e++)
        {
            var row = s.Row(e);
            if (row.Length == 0) continue;
            ElementMatrix(x, row, ke);
            var dest = s.Destinations(e);
            for (var i = 0; i < 4; i++)
            {
                var rs = s.PatternRowStart(row[i]);
                for (var j = 0; j < 4; j++) vals[rs + dest[4 * i + j]] += ke[4 * i + j];
            }
        }
    }

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
            "step,t_maint_ms,t_asm_ms,t_solve_ms,res_inc,t_mm3_rebuild_ms,t_mm3_asm_ms,t_mm3_solve_ms,res_mm3");
        var tm = new List<double>();
        var ta = new List<double>();
        var ts = new List<double>();
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
            var vals = new double[snap.TotalSlots];
            FullGC();
            var tAsm = TimeMs(() => AssembleGapped(snap, st, vals));
            var diag = DiagGapped(snap, vals);
            double res = 0;
            FullGC();
            var tSolve = TimeMs(() =>
                res = Cg(snap.NodeCount, (xin, yout) => SpmvGapped(snap, vals, xin, yout), diag, iters));
            tm.Add(tMaint);
            ta.Add(tAsm);
            ts.Add(tSolve);

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

            w.WriteLine($"{s},{tMaint:F4},{tAsm:F3},{tSolve:F3},{res:E6},{tR:F3},{tRa:F3},{tRs:F3},{resM:E6}");
        }

        var inc = Median(tm) + Median(ta) + Median(ts);
        var mm3 = Median(tr) + Median(tra) + Median(trs);
        var line = $"endtoend n={n}: incremental step = maint {Median(tm):F2} + asm {Median(ta):F1} + solve {Median(ts):F1} " +
                   $"= {inc:F1} ms (maint {100 * Median(tm) / inc:F2} %); MM3 step = rebuild {Median(tr):F1} + asm " +
                   $"{Median(tra):F1} + solve {Median(trs):F1} = {mm3:F1} ms (rebuild {100 * Median(tr) / mm3:F1} %); " +
                   $"end-to-end ratio {mm3 / inc:F2}x";
        Console.WriteLine(line);
        File.AppendAllText(Path.Combine(DataDir, "endtoend_summary.txt"), line + Environment.NewLine);
    }
}

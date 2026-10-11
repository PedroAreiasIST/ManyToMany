// CornerRun.cs — driver of the corner-crack example (see CornerCrack.cs).
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
//   corner [n=20] [W=1] [H=1] [beta=30] [a0=0.1] [b0=0.1] [hmin=0.0125] [gr=0.6] [zlo=-0.45] [zhi=0.05] [rz=1.0]
//          [steps=12] [da=0.07] [m=3] [rq=] [M=64] [tol=1e-9] [out=<dir>] [ens=<dir>] [ensevery=1]
//          [patch] [penny] [nofs] [oracle] [est]
//
// penny: an embedded circular crack of radius a0 at the centre of the block (W = 2, H = 1 for the faces at 5 a0), the
// front closed by repeating its first knot; K at 16 interior knots against Kassir and Sih, as in Fracture.Run.

using System.Diagnostics;
using Step2;

namespace Step2.Bench;

public static partial class Fracture
{
    public static void Corner(string[] a)
    {
        double D(string key, double dflt) =>
            a.Where(x => x.StartsWith(key + "=")).Select(x => double.Parse(x[(key.Length + 1)..], Ci)).DefaultIfEmpty(dflt).First();
        string S(string key, string dflt) =>
            a.Where(x => x.StartsWith(key + "=")).Select(x => x[(key.Length + 1)..]).DefaultIfEmpty(dflt).First();
        var penny = a.Contains("penny");
        var patch = a.Contains("patch");
        var freeSurface = !a.Contains("nofs") && !penny;
        var n = (int)D("n", 20);
        var W = D("W", penny ? 2.0 : 1.0);
        var H = D("H", 1.0);
        var beta = D("beta", 30) * Math.PI / 180;
        var a0 = D("a0", penny ? 0.2 : 0.10);
        var b0 = D("b0", a0);
        var hmin = D("hmin", 0.0125);
        var gr = D("gr", penny ? 0.0 : 0.6);
        double zlo = D("zlo", -0.45), zhi = D("zhi", 0.05), rz = D("rz", 1.0);
        var maxSteps = (int)D("steps", 12);
        var daRel = D("dar", 0.25);   // advance of the largest knot as a fraction of the mean front radius
        var daFix = D("da", 0);
        var mExp = D("m", 3.0);
        var M = (int)D("M", 64);
        var tolCg = D("tol", 1e-9);
        var rqRel = D("rqr", 0.15);
        var kappa = D("kappa", 0);   // mesh size kappa (rho + rho0) with the distance rho from the corner (0: hmin (1 + gr rho))
        var rho0 = D("rho0", 0.05);
        var kinkMax = D("kinkmax", 45) * Math.PI / 180;
        var smoothPasses = (int)D("smooth", 3);
        var curv = D("curv", 12);
        var fsmooth = (int)D("fsmooth", 3);
        var borderMargin = D("border", 0.08);
        var dasmooth = (int)D("dasmooth", 6);
        var outDir = S("out", DataDirOf());
        var ensDir = S("ens", "");
        var ensEvery = (int)D("ensevery", 1);
        var nx = (int)Math.Round(W * n);
        var ny = nx;
        var nz = (int)Math.Round(2 * H * n);
        var delta0 = 0.37 * (kappa > 0 ? kappa * rho0 : hmin);
        var cx = penny ? W / 2 : 0.0;
        var cy = penny ? W / 2 : 0.0;
        var pl = new Plane(beta, delta0, cx, cy);
        var sw = Stopwatch.StartNew();
        Directory.CreateDirectory(outDir);

        // ---- mesh size: fine in the zone the crack can reach, graded outside
        const double grading = 0.35;
        double Size(double x, double y, double z)
        {
            if (penny)
            {
                var (p, q, h) = pl.Local(x, y, z);
                var rho = Math.Sqrt(p * p + q * q);
                var d = Math.Sqrt(Math.Max(0, rho - a0) * Math.Max(0, rho - a0) + h * h);
                return hmin + grading * Math.Max(0, d - 0.06);
            }

            var r = Math.Sqrt(x * x + y * y);
            var dz = Math.Max(0, Math.Max(zlo - z, z - zhi));
            var dr = Math.Max(0, r - rz);
            var dd = Math.Sqrt(dz * dz + dr * dr);
            var hz = kappa > 0 ? kappa * (Math.Min(r, rz) + rho0) : hmin * (1 + gr * Math.Min(r, rz));
            return hz + grading * Math.Max(0, dd - 0.03);
        }

        if (a.Contains("est"))
        {
            // volume integral of 6 / size^3 on a grid: an estimate of the tetrahedra before bisection
            double tot = 0;
            const int g = 160;
            for (var i = 0; i < g; i++)
            for (var j = 0; j < g; j++)
            for (var k = 0; k < g; k++)
            {
                var x = W * (i + 0.5) / g;
                var y = W * (j + 0.5) / g;
                var z = -H + 2 * H * (k + 0.5) / g;
                var sz = Size(x, y, z);
                tot += 6 / (sz * sz * sz) * (W / g) * (W / g) * (2 * H / g);
            }

            Console.WriteLine($"  estimate: {tot:F0} tetrahedra (for 1.5 h as the bisection limit divide by about 2..4)");
            return;
        }

        var (m, st, passes, bisections) = RefinedCorner(nx, ny, nz, W, H, Size);
        var quality = Quality(m, st);
        Console.WriteLine($"  mesh: block {W} x {W} x {2 * H}, h0 = {1.0 / n:F4}, hmin = {hmin:F5}, {passes} bisection passes " +
                          $"({bisections} bisections), {Live(m)} tetrahedra, {m.NodeCount} nodes; radius ratio min " +
                          $"{quality.Min:F3}, mean {quality.Mean:F3} ({sw.Elapsed.TotalSeconds:F1} s)");
        if (a.Contains("meshonly")) return;
        var xyz = new double[3 * m.NodeCount];
        for (var v = 0; v < m.NodeCount; v++)
        for (var c = 0; c < 3; c++)
            xyz[3 * v + c] = st.X(v, c);
        var ins = new CrackInserter(m, xyz);
        m.CompactPattern();

        // ---- initial front and crack surface
        var up = new[] { pl.N[0], pl.N[1], pl.N[2] };
        var (pc, qc, _) = pl.Local(0, 0, delta0 / Math.Cos(beta)); // the corner in local coordinates
        if (penny) (pc, qc) = (0, 0);
        double[] P0;
        Func<int, (double, double, double)> n0;
        int nk;
        if (penny)
        {
            nk = M + 1;
            P0 = new double[3 * nk];
            for (var j = 0; j < nk; j++)
            {
                var psi = 2 * Math.PI * j / M;
                var (gx, gy, gz) = pl.Global(a0 * Math.Cos(psi), a0 * Math.Sin(psi), 0);
                P0[3 * j] = gx;
                P0[3 * j + 1] = gy;
                P0[3 * j + 2] = gz;
            }

            n0 = j =>
            {
                var psi = 2 * Math.PI * j / M;
                var gx = Math.Cos(psi) * pl.E1[0] + Math.Sin(psi) * pl.E2[0];
                var gy = Math.Cos(psi) * pl.E1[1] + Math.Sin(psi) * pl.E2[1];
                var gz = Math.Cos(psi) * pl.E1[2] + Math.Sin(psi) * pl.E2[2];
                return (gx, gy, gz);
            };
        }
        else
        {
            nk = M;
            P0 = new double[3 * nk];
            for (var j = 0; j < nk; j++)
            {
                var psi = 0.5 * Math.PI * j / (nk - 1);
                var (gx, gy, gz) = pl.Global(pc + a0 * Math.Cos(psi), qc + b0 * Math.Sin(psi), 0);
                P0[3 * j] = gx;
                P0[3 * j + 1] = gy;
                P0[3 * j + 2] = gz;
            }

            n0 = j =>
            {
                var psi = 0.5 * Math.PI * j / (nk - 1);
                double ex = Math.Cos(psi) / a0, ey = Math.Sin(psi) / b0;
                var l = Math.Sqrt(ex * ex + ey * ey);
                ex /= l;
                ey /= l;
                return (ex * pl.E1[0] + ey * pl.E2[0], ex * pl.E1[1] + ey * pl.E2[1], ex * pl.E1[2] + ey * pl.E2[2]);
            };
        }

        var front = MakeFront(P0, null, up, n0);
        var surface = new CSurface(-0.9, -0.3, 1.8, 1.6, Math.Max((kappa > 0 ? kappa * rho0 : hmin) / 2, 0.0015));
        double[] LocalFront(CFront f)
        {
            var o = new double[3 * f.M];
            for (var j = 0; j < f.M; j++)
            {
                var (p, q, h) = pl.Local(f.P[3 * j], f.P[3 * j + 1], f.P[3 * j + 2]);
                o[3 * j] = p;
                o[3 * j + 1] = q;
                o[3 * j + 2] = h;
            }

            return o;
        }

        double[] Polygon(double[] fl, int cnt)
        {
            var poly = new List<double>();
            if (penny)
            {
                for (var j = 0; j < cnt; j++)
                {
                    poly.Add(fl[3 * j]);
                    poly.Add(fl[3 * j + 1]);
                }
            }
            else
            {
                const double L = 0.8;
                // y = 0 end first (index 0), x = 0 end last; extended beyond the faces so that the body lies inside
                poly.Add(fl[0]); poly.Add(fl[1] - L);
                for (var j = 0; j < cnt; j++)
                {
                    poly.Add(fl[3 * j]);
                    poly.Add(fl[3 * j + 1]);
                }

                poly.Add(fl[3 * (cnt - 1)] - L); poly.Add(fl[3 * (cnt - 1) + 1]);
                poly.Add(Math.Min(pc, fl[0]) - L); poly.Add(Math.Min(qc, fl[1]) - L);
            }

            return poly.ToArray();
        }

        void AddFan(double[] fl, int cnt)
        {
            Span<double> c0 = stackalloc double[3] { pc, qc, 0 };
            for (var j = 0; j + 1 < cnt; j++)
                surface.AddTriangle(c0, fl.AsSpan(3 * j, 3), fl.AsSpan(3 * j + 3, 3));
        }

        var fl0 = LocalFront(front);
        AddFan(fl0, front.M);
        surface.Finish(fl0, Polygon(fl0, front.M), 14);

        // The level set of step k is the chain  W_k = W_{k-1} behind the front of step k-1, a smooth blend of W_{k-1} and the
        // surface S_k over three grid cells ahead of it: behind every earlier front it reproduces exactly the field in
        // which the nodes there were cut, so no node already classified changes side.
        var levels = new List<(CSurface S, CSurface? Fprev)> { (surface, null) };
        double Wk(int k, double p, double q, double margin)
        {
            if (k == 0) return levels[0].S.W(p, q);
            // frozen where the region function of the previous step was negative (distance to its front below the margin)
            var d = levels[k].Fprev!.SignedDistance(p, q) - margin;
            if (d <= 0) return Wk(k - 1, p, q, margin);
            var lam = Math.Min(1.0, d / (3 * levels[k].S.Cell));
            lam = lam * lam * (3 - 2 * lam);
            return (1 - lam) * Wk(k - 1, p, q, margin) + lam * levels[k].S.W(p, q);
        }

        CrackInserter.Field phiS = (x, y, z) =>
        {
            var (p, q, h) = pl.Local(x, y, z);
            return h - Wk(levels.Count - 1, p, q, 0.25 * Size(x, y, z));
        };
        CrackInserter.Field phiR = (x, y, z) =>
        {
            var (p, q, _) = pl.Local(x, y, z);
            return surface.SignedDistance(p, q) - 0.25 * Size(x, y, z);
        };

        var tag = penny ? "penny" : "corner";
        using var wk = new StreamWriter(Path.Combine(outDir, $"{tag}_K.csv"));
        using var wf = new StreamWriter(Path.Combine(outDir, $"{tag}_fronts.csv"));
        using var ws = new StreamWriter(Path.Combine(outDir, $"{tag}_steps.csv"));
        wk.WriteLine("step,knot,s,x,y,z,KI,KII,KIII,Keq,theta_deg,da");
        wf.WriteLine("step,knot,x,y,z");
        ws.WriteLine("step,da,t_insert_s,t_publish_ms,t_assemble_s,cg_iters,t_solve_s,t_K_s,elements,nodes,unknowns,crack_edges," +
                     "split_nodes,front_nodes,visited,max_theta_deg,front_length,rq");
        var fp = new FlatPattern();
        double[]? uPrev = null;
        var ensSteps = new List<int>();
        var ensTimes = new List<double>();
        var accepted = 0;
        double traveled = 0;
        Func<double, double, double, (CFront F, CSurface S, double Adv)?>? retryBuild = null;
        Snapshot? lastSnap = null;
        int[]? lastFrontNodes = null;
        double lastAdv = 0;
        var retries = 0;
        while (true)
        {
            var tIns = Stopwatch.StartNew();
            CrackStep? cstep = null;
            if (!patch)
            {
                var tries = 0;
                while (true)
                {
                    try
                    {
                        cstep = ins.Advance(phiS, phiR, global: accepted == 0);
                        break;
                    }
                    catch (InvalidOperationException ex) when (retryBuild != null && tries < 8)
                    {
                        tries++;
                        var mt = System.Text.RegularExpressions.Regex.Match(ex.Message, @"node (\d+)");
                        if (mt.Success)
                        {
                            var vbad = int.Parse(mt.Groups[1].Value);
                            var (pp, qq, hh) = pl.Local(ins.X[vbad], ins.Y[vbad], ins.Z[vbad]);
                            Console.WriteLine($"    node {vbad} at ({ins.X[vbad]:F4}, {ins.Y[vbad]:F4}, {ins.Z[vbad]:F4}), local (p {pp:F4}, q {qq:F4}, h {hh:F4}), " +
                                              $"phi_s {phiS(ins.X[vbad], ins.Y[vbad], ins.Z[vbad]):E2}, phi_r {phiR(ins.X[vbad], ins.Y[vbad], ins.Z[vbad]):E2}");
                        }

                        m.Rollback(lastSnap!);
                        ins.RestoreFront(lastFrontNodes!);
                        double[] kSch = [1, 0.8, 1, 0.8, 0.6, 0.6, 0.4, 0.4], aSch = [1, 0.9, 0.9, 0.8, 0.8, 0.7, 0.7, 0.6],
                            shSch = [0.06, -0.06, 0.12, -0.12, 0.2, -0.2, 0.3, -0.3], tauSch = [0.05, 0.08, 0.08, 0.12, 0.12, 0.16, 0.16, 0.2];
                        var ks = kSch[tries - 1];
                        var asc = aSch[tries - 1];
                        var shiftRel = shSch[tries - 1];
                        ins.Tau = tauSch[tries - 1];
                        Console.WriteLine($"    insertion failed ({ex.Message.Split('\n')[0]}): rolled back to the committed version, " +
                                          $"retry {tries} with kink x{ks:F2}, advance x{asc:F2}, level set moved by {shiftRel:F2} h, promotion tolerance {ins.Tau:F2}");
                        traveled -= lastAdv;
                        var b = retryBuild(ks, asc, shiftRel * Size(front.P[3 * (front.M / 2)], front.P[3 * (front.M / 2) + 1], front.P[3 * (front.M / 2) + 2]));
                        if (b == null) throw;
                        front = b.Value.F;
                        surface = b.Value.S;
                        levels[^1] = (b.Value.S, levels[^1].Fprev);
                        lastAdv = b.Value.Adv;
                        traveled += lastAdv;
                        retries++;
                    }
                }
            }
            tIns.Stop();
            var tPub = Stopwatch.StartNew();
            var s = m.Commit();
            tPub.Stop();
            if (a.Contains("oracle"))
            {
                var tO = Stopwatch.StartNew();
                Oracle.Verify(s);
                Console.WriteLine($"  oracle: all structures of the committed version match ({tO.Elapsed.TotalSeconds:F1} s)");
            }

            var (X, f, fixedDof) = CornerProblem(s, ins, nx, ny, H);
            fp.Refresh(s);
            var tAsm = Stopwatch.StartNew();
            var vals = Assemble(s, fp, X);
            tAsm.Stop();
            double[]? u0 = null;
            if (uPrev != null)
            {
                u0 = new double[3 * s.NodeCount];
                var nPrev = uPrev.Length / 3;
                Array.Copy(uPrev, u0, uPrev.Length);
                var p1 = s.View(ins.Parent1);
                var p2 = s.View(ins.Parent2);
                var wt = s.View(ins.Weight);
                for (var v = nPrev; v < s.NodeCount; v++)
                for (var k = 0; k < 3; k++)
                    u0[3 * v + k] = (1 - wt[v]) * u0[3 * p1[v] + k] + wt[v] * u0[3 * p2[v] + k];
            }

            var tSol = Stopwatch.StartNew();
            var (u, iters, rel) = Pcg(fp, vals, f, fixedDof, tolCg, 60000, u0);
            tSol.Stop();
            if (!(rel <= 1e-7)) throw new Exception($"step {accepted}: CG did not converge ({iters} iterations, residual {rel:E2})");
            vals = null;
            Console.WriteLine($"  step {accepted}: {Live(m)} tetrahedra, {s.NodeCount} nodes, solve {3 * s.NodeCount} unknowns, " +
                              $"CG {iters} it {tSol.Elapsed.TotalSeconds:F1} s, residual {rel:E1}");

            if (patch)
            {
                double err = 0, scale = 0;
                for (var v = 0; v < s.NodeCount; v++)
                {
                    if (s.TransposePairs(v).Length == 0) continue;
                    double ex = -Nu * Sigma / E * X[3 * v], ey = -Nu * Sigma / E * X[3 * v + 1], ez = Sigma / E * (X[3 * v + 2] + H);
                    err = Math.Max(err, Math.Max(Math.Abs(u[3 * v] - ex), Math.Max(Math.Abs(u[3 * v + 1] - ey), Math.Abs(u[3 * v + 2] - ez))));
                    scale = Math.Max(scale, Math.Abs(ez));
                }

                Console.WriteLine($"corner patch test: max nodal error {err / scale:E2} relative to the exact linear field");
                return;
            }

            // ---- K along the front
            var lengthEff = penny ? a0 : front.Length / (0.5 * Math.PI);
            var mid = front.M / 2;
            var hLoc = Size(front.P[3 * mid], front.P[3 * mid + 1], front.P[3 * mid + 2]);
            var rq = a.Any(x => x.StartsWith("rq=")) ? D("rq", 0.04) : Math.Clamp(Math.Max(rqRel * lengthEff, 6 * hLoc), 0.03, 0.30);
            var daBase = daFix > 0 ? daFix : daRel * lengthEff;
            var delta = Math.Max(rq, 2.5 * front.Length / (front.M - 1));
            var tK = Stopwatch.StartNew();
            var tube = BuildTube(s, X, front, rq);
            var vn = new double[3 * s.NodeCount];
            var anyV = new bool[s.NodeCount];
            var K = new (double I, double II, double III, double J)[front.M];
            var knots = penny
                ? Enumerable.Range(0, front.M).Where(j => front.S[j] >= delta && front.S[j] <= front.Length - delta).ToArray()
                : Enumerable.Range(0, front.M).ToArray();
            foreach (var j in knots)
                K[j] = CInteraction(s, X, u, front, tube, front.S[j], delta, freeSurface, W, vn, anyV);
            tK.Stop();

            if (penny)
            {
                var kRef = 2 * Sigma * Math.Sqrt(a0 / Math.PI);
                double sn = Sigma * Math.Cos(beta) * Math.Cos(beta), tau = Sigma * Math.Sin(beta) * Math.Cos(beta);
                double errI = 0, errII = 0, errIII = 0;
                const int points = 16;
                for (var kq = 0; kq < points; kq++)
                {
                    var psiT = 2 * Math.PI * (kq + 0.5) / points;
                    var j = (int)Math.Round(psiT / (2 * Math.PI) * M);
                    if (!knots.Contains(j)) continue;
                    var th = 2 * Math.PI * j / M;
                    var refI = 2 * sn * Math.Sqrt(a0 / Math.PI);
                    var refII = 4 / (2 - Nu) * tau * Math.Sqrt(a0 / Math.PI) * Math.Cos(th);
                    var refIII = 4 * (1 - Nu) / (2 - Nu) * tau * Math.Sqrt(a0 / Math.PI) * Math.Sin(th);
                    Console.WriteLine($"  psi {th * 180 / Math.PI,6:F1}: K_I/K0 {K[j].I / kRef,8:F4} (ref {refI / kRef,7:F4})  " +
                                      $"K_II/K0 {K[j].II / kRef,8:F4} (ref {refII / kRef,7:F4})  K_III/K0 {K[j].III / kRef,8:F4} (ref {refIII / kRef,7:F4})");
                    errI = Math.Max(errI, Math.Abs(K[j].I - refI) / kRef);
                    errII = Math.Max(errII, Math.Abs(Math.Abs(K[j].II) - Math.Abs(refII)) / kRef);
                    errIII = Math.Max(errIII, Math.Abs(Math.Abs(K[j].III) - Math.Abs(refIII)) / kRef);
                }

                Console.WriteLine($"penny (corner machinery) beta={beta * 180 / Math.PI:F0}: max |K - K_ref|/K0 = {errI:P2} (I), {errII:P2} (|II|), " +
                                  $"{errIII:P2} (|III|); rq={rq:F3}, delta={delta:F3}, {Live(m)} tetrahedra, {3 * s.NodeCount} unknowns " +
                                  $"({sw.Elapsed.TotalSeconds:F0} s)");
                return;
            }

            // ---- Erdogan-Sih directions and Paris advances
            // K_I and K_II used for the direction and the advance: within one hat width of the free faces they are
            // held at the value of the knot one hat width inside (the corner singularity and the clipped hat make
            // the values at the ends unreliable); K_III does not enter (two-dimensional criterion in the plane N, B)
            var kuI = K.Select(t => t.I).ToArray();
            var kuII = K.Select(t => t.II).ToArray();
            int jLo = Enumerable.Range(0, front.M).First(j => front.S[j] >= delta);
            int jHi = Enumerable.Range(0, front.M).Last(j => front.S[j] <= front.Length - delta);
            for (var j = 0; j < jLo; j++) { kuI[j] = K[jLo].I; kuII[j] = K[jLo].II; }
            for (var j = jHi + 1; j < front.M; j++) { kuI[j] = K[jHi].I; kuII[j] = K[jHi].II; }
            var th0 = new double[front.M];
            var keq = new double[front.M];
            for (var j = 0; j < front.M; j++)
            {
                th0[j] = ErdoganSih(kuI[j], kuII[j]);
                keq[j] = Math.Sqrt(Math.Max(0, kuI[j] * kuI[j] + kuII[j] * kuII[j]));
            }

            // smoothing along the front (the ends are on the free faces)
            for (var pass = 0; pass < Math.Max(smoothPasses, dasmooth); pass++)
            {
                var t2 = (double[])th0.Clone();
                var k2 = (double[])keq.Clone();
                for (var j = 1; j + 1 < front.M; j++)
                {
                    if (pass < smoothPasses) t2[j] = 0.25 * th0[j - 1] + 0.5 * th0[j] + 0.25 * th0[j + 1];
                    k2[j] = 0.25 * keq[j - 1] + 0.5 * keq[j] + 0.25 * keq[j + 1];
                }

                th0 = t2;
                keq = k2;
            }

            var kinkLimit = Math.Min(kinkMax, daBase / (curv * hLoc)); // radius of curvature of the path at least curv element sizes
            for (var j = 0; j < front.M; j++) th0[j] = Math.Clamp(th0[j], -kinkLimit, kinkLimit);
            var maxTh = th0.Max(t => Math.Abs(t)) * 180 / Math.PI;
            var kMax = keq.Max();
            var daJ = keq.Select(k => daBase * Math.Pow(k / kMax, mExp)).ToArray();
            for (var j = 0; j < front.M; j++)
            {
                var (kI, kII, kIII, _) = K[j];
                wk.WriteLine(string.Format(Ci, "{0},{1},{2:F5},{3:F5},{4:F5},{5:F5},{6:E5},{7:E5},{8:E5},{9:E5},{10:F3},{11:F5}", accepted, j,
                    front.S[j], front.P[3 * j], front.P[3 * j + 1], front.P[3 * j + 2], kI, kII, kIII, keq[j], th0[j] * 180 / Math.PI, daJ[j]));
                wf.WriteLine(string.Format(Ci, "{0},{1},{2:R},{3:R},{4:R}", accepted, j, front.P[3 * j], front.P[3 * j + 1], front.P[3 * j + 2]));
            }

            wk.Flush();
            wf.Flush();
            ws.WriteLine(string.Format(Ci, "{0},{1:F5},{2:F2},{3:F3},{4:F2},{5},{6:F1},{7:F1},{8},{9},{10},{11},{12},{13},{14},{15:F2},{16:F4},{17:F4}",
                accepted, daBase, tIns.Elapsed.TotalSeconds, tPub.Elapsed.TotalMilliseconds, tAsm.Elapsed.TotalSeconds, iters,
                tSol.Elapsed.TotalSeconds, tK.Elapsed.TotalSeconds, Live(m), s.NodeCount, 3 * s.NodeCount, cstep!.Edges.Length,
                cstep.Splits.Length, cstep.Front.Length, cstep.VisitedTets, maxTh, front.Length, rq));
            ws.Flush();
            Console.WriteLine($"    K_I {K.Min(t => t.I):F4}..{K.Max(t => t.I):F4}, |K_II| max {K.Max(t => Math.Abs(t.II)):F4}, |K_III| max " +
                              $"{K.Max(t => Math.Abs(t.III)):F4}, kink up to {maxTh:F1} deg, front length {front.Length:F3}, rq {rq:F3}, " +
                              $"insert {tIns.Elapsed.TotalSeconds:F1} s ({cstep.Edges.Length} crack edges, {cstep.Splits.Length} splits), K {tK.Elapsed.TotalSeconds:F1} s");

            // ---- output of this step (EnSight Gold, ASCII)
            if (ensDir.Length > 0 && accepted % ensEvery == 0)
            {
                var tE = Stopwatch.StartNew();
                var nNodes = s.NodeCount;
                var p1 = s.View(ins.Parent1);
                var p2 = s.View(ins.Parent2);
                var opening = new double[nNodes];
                for (var v = 0; v < nNodes; v++)
                {
                    if (p1[v] == v || p1[v] != p2[v]) continue;
                    var (p, q, _) = pl.Local(X[3 * v], X[3 * v + 1], X[3 * v + 2]);
                    const double dq = 1e-3;
                    var dwp = (surface.W(p + dq, q) - surface.W(p - dq, q)) / (2 * dq);
                    var dwq = (surface.W(p, q + dq) - surface.W(p, q - dq)) / (2 * dq);
                    double nx_ = pl.N[0] - dwp * pl.E1[0] - dwq * pl.E2[0], ny_ = pl.N[1] - dwp * pl.E1[1] - dwq * pl.E2[1],
                        nz_ = pl.N[2] - dwp * pl.E1[2] - dwq * pl.E2[2];
                    var nl = Math.Sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_);
                    var o = p1[v];
                    opening[v] = ((u[3 * v] - u[3 * o]) * nx_ + (u[3 * v + 1] - u[3 * o + 1]) * ny_ + (u[3 * v + 2] - u[3 * o + 2]) * nz_) / nl;
                }

                var vm = new double[s.ElementCount];
                Parallel.For(0, s.ElementCount, e =>
                {
                    var row = s.Row(e);
                    if (row.Length != 4) return;
                    Span<double> g = stackalloc double[12];
                    Gradients(X, row, g);
                    Span<double> G = stackalloc double[9];
                    G.Clear();
                    for (var aa = 0; aa < 4; aa++)
                    for (var i = 0; i < 3; i++)
                    for (var j = 0; j < 3; j++)
                        G[3 * i + j] += u[3 * row[aa] + i] * g[3 * aa + j];
                    var tr = G[0] + G[4] + G[8];
                    double s11 = Mu * 2 * G[0] + Lambda * tr, s22 = Mu * 2 * G[4] + Lambda * tr, s33 = Mu * 2 * G[8] + Lambda * tr;
                    double s12 = Mu * (G[1] + G[3]), s13 = Mu * (G[2] + G[6]), s23 = Mu * (G[5] + G[7]);
                    vm[e] = Math.Sqrt(0.5 * ((s11 - s22) * (s11 - s22) + (s22 - s33) * (s22 - s33) + (s33 - s11) * (s33 - s11)) +
                                      3 * (s12 * s12 + s13 * s13 + s23 * s23));
                });
                var nComp = WriteEnsightStep(ensDir, "corner", true, accepted, s, ins, X, u, front, W, H, opening, vm);
                Console.WriteLine($"    upper crack face: {nComp} connected component(s)" + (nComp > 1 ? "  <-- WARNING: more than one" : ""));
                WriteEnsightStep(Path.Combine(ensDir, "crack"), "crack", false, accepted, s, ins, X, u, front, W, H, opening, vm);
                ensSteps.Add(accepted);
                ensTimes.Add(traveled);
                WriteEnsightCase(ensDir, "corner", true, ensSteps.ToArray(), ensTimes.ToArray());
                WriteEnsightCase(Path.Combine(ensDir, "crack"), "crack", false, ensSteps.ToArray(), ensTimes.ToArray());
                Console.WriteLine($"    EnSight Gold step {accepted} written to {ensDir} ({tE.Elapsed.TotalSeconds:F1} s)");
            }

            uPrev = u;
            if (accepted >= maxSteps) break;

            // ---- new front: P'_j = P_j + da_j (cos theta N_j + sin theta B_j), clipped to the block. Kept as a function of the
            // kink and advance scales so that an insertion that fails can be retried with a gentler step after a rollback.
            var frontOld = front;
            var surfOld = surface;
            var th0Old = th0;
            var daOld = daJ;
            lastSnap = s;
            lastFrontNodes = (int[])ins.Front.Clone();
            (CFront F, CSurface S, double Adv)? Build(double kScale, double aScale, double hShift = 0)
            {
                var fo = frontOld;
                var lim = Math.Max(3 * Size(fo.P[3 * (fo.M / 2)], fo.P[3 * (fo.M / 2) + 1], fo.P[3 * (fo.M / 2) + 2]), borderMargin);
                var thk = th0Old.Select(t => t * kScale).ToArray();
                double[] Trial(double fct)
                {
                    var P2 = new double[3 * fo.M];
                    for (var j = 0; j < fo.M; j++)
                    {
                        var d = fct * daOld[j];
                        var ct = Math.Cos(thk[j]);
                        var stt = Math.Sin(thk[j]);
                        for (var c = 0; c < 3; c++)
                            P2[3 * j + c] = fo.P[3 * j + c] + d * (ct * fo.N[3 * j + c] + stt * fo.B[3 * j + c]);
                    }

                    for (var j = 0; j < fo.M; j++)
                    {
                        P2[3 * j] = Math.Max(0, P2[3 * j]);         // the free faces x = 0 and y = 0: the front slides on them
                        P2[3 * j + 1] = Math.Max(0, P2[3 * j + 1]);
                    }

                    // light smoothing of the new polyline (a concave shoulder would fold under unequal advances)
                    for (var pass = 0; pass < fsmooth; pass++)
                    {
                        var Q = (double[])P2.Clone();
                        for (var j = 1; j + 1 < fo.M; j++)
                        for (var c = 0; c < 3; c++)
                            Q[3 * j + c] = 0.25 * P2[3 * j - 3 + c] + 0.5 * P2[3 * j + c] + 0.25 * P2[3 * j + 3 + c];
                        P2 = Q;
                    }

                    for (var j = 0; j < fo.M; j++)
                    {
                        P2[3 * j] = Math.Max(0, P2[3 * j]);
                        P2[3 * j + 1] = Math.Max(0, P2[3 * j + 1]);
                    }

                    P2[1] = 0;                              // the y = 0 end stays on its face
                    P2[3 * (fo.M - 1)] = 0;                 // the x = 0 end stays on its face
                    return P2;
                }

                bool Inside(double[] P2) => Enumerable.Range(0, fo.M).All(j =>
                    P2[3 * j] <= W - lim && P2[3 * j + 1] <= W - lim && P2[3 * j + 2] >= -H + lim && P2[3 * j + 2] <= H - lim);

                var factor = aScale;
                var trial = Trial(factor);
                while (!Inside(trial) && factor > 0.05 * aScale)
                {
                    factor *= 0.8;
                    trial = Trial(factor);
                }

                if (!Inside(trial))
                {
                    Console.WriteLine($"  the front is within {lim:F3} of the border: stop");
                    return null;
                }

                var ns = surfOld.Clone();
                // Hermite strips: from P_j with the previous advance direction N_j to P'_j with the new direction d_j, so
                // that the crack surface has no crease at the old front (the chord length scales both tangents)
                var dvec = new double[3 * fo.M];
                for (var j = 0; j < fo.M; j++)
                {
                    var ct = Math.Cos(thk[j]);
                    var stt = Math.Sin(thk[j]);
                    for (var c = 0; c < 3; c++) dvec[3 * j + c] = ct * fo.N[3 * j + c] + stt * fo.B[3 * j + c];
                }

                const int nsub = 6;
                var sub = new double[nsub + 1][];
                for (var k = 0; k <= nsub; k++)
                {
                    var t = (double)k / nsub;
                    double h00 = 2 * t * t * t - 3 * t * t + 1, h10 = t * t * t - 2 * t * t + t, h01 = -2 * t * t * t + 3 * t * t,
                        h11 = t * t * t - t * t;
                    var gf = new double[3 * fo.M];
                    for (var j = 0; j < fo.M; j++)
                    {
                        double chord = Math.Sqrt(Math.Pow(trial[3 * j] - fo.P[3 * j], 2) + Math.Pow(trial[3 * j + 1] - fo.P[3 * j + 1], 2) +
                                                 Math.Pow(trial[3 * j + 2] - fo.P[3 * j + 2], 2));
                        for (var c = 0; c < 3; c++)
                            gf[3 * j + c] = h00 * fo.P[3 * j + c] + h10 * chord * fo.N[3 * j + c] + h01 * trial[3 * j + c] +
                                            h11 * chord * dvec[3 * j + c];
                    }

                    var fl = new double[3 * fo.M];
                    for (var j = 0; j < fo.M; j++)
                    {
                        var (p, q, h) = pl.Local(gf[3 * j], gf[3 * j + 1], gf[3 * j + 2]);
                        fl[3 * j] = p;
                        fl[3 * j + 1] = q;
                        fl[3 * j + 2] = h + hShift * t; // the retry perturbation, ramped from 0 at the old front
                    }

                    sub[k] = fl;
                }

                for (var k = 0; k < nsub; k++)
                for (var j = 0; j + 1 < fo.M; j++)
                {
                    var f0 = sub[k];
                    var f1 = sub[k + 1];
                    ns.AddTriangle(f0.AsSpan(3 * j, 3), f0.AsSpan(3 * j + 3, 3), f1.AsSpan(3 * j, 3));
                    ns.AddTriangle(f0.AsSpan(3 * j + 3, 3), f1.AsSpan(3 * j + 3, 3), f1.AsSpan(3 * j, 3));
                }

                var dres = new double[3 * M];
                var res = Resample(trial, M, dvec, dres);
                var newFront = MakeFront(res, fo, up, null, dres);
                var flRes = LocalFront(newFront);
                for (var j = 0; j < newFront.M; j++) flRes[3 * j + 2] += hShift;
                ns.Finish(flRes, Polygon(flRes, newFront.M), 14);
                return (newFront, ns, factor * daBase);
            }

            ins.Tau = 0.05;
            var built = Build(1, 1);
            if (built == null) break;
            retryBuild = Build;
            levels.Add((built.Value.S, surfOld));
            front = built.Value.F;
            surface = built.Value.S;
            lastAdv = built.Value.Adv;
            traveled += lastAdv;
            accepted++;
        }

        Console.WriteLine($"corner: {accepted} growth steps, {retries} insertion retries after rollback, {Live(m)} tetrahedra ({sw.Elapsed.TotalSeconds:F0} s)");
    }
}

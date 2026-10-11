// SurfaceCrack.cs — semi-elliptical surface crack in a finite plate: Newman-Raju stress intensity factors and
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
// Paris-law fatigue growth on the maintained structures.
//
//   surface n div a c steps damax m tolK rq [M] [t=<thickness>] [h=<half-length>] [patch] [nofs] [rebuild] [oracle] [dump=<dir>] [u=<file>]
//
// Plate (bar) x in [-1,1] (b = 1), y in [-t,0], z in [-h,h] (default t = 1, h = 2), E = 1, nu = 0.3, remote tension
// sigma = 1 along z on the faces z = +-h, rigid motions removed by three corner constraints (the uniform field is
// exact). The crack lies in the plane z = 0.37 hmin, centred on the face y = -t: the front is the star-shaped
// curve rho = R(psi) of FrontCurve about the centre (0, -t), of which the arc psi in [0, pi] lies in the body, and
// meets the free surface at psi = 0, pi. Mesh, fitting, crack insertion, assembly, solution and interaction
// integral are those of Fracture.cs (the front curve is closed, the half outside the body has no elements; the
// hat of the interaction integral is normalized by its length inside the body).
// Reference: Newman and Raju (1981), K_I = sigma sqrt(pi a / Q) F_s(a/c, a/t, c/b, phi) at the parametric angle phi
// of the front point, a = R(pi/2), c = R(0). Growth: da_j = da (K_j / K_max)^m along the normal, K mirrored about
// the plane psi = 0 so that the front stays symmetric.

using Step2;

namespace Step2.Bench;

public static partial class Fracture
{
    private const double BarB = 1.0;

    /// <summary>Newman-Raju (1981) K_I (sigma = 1) of a semi-elliptical surface crack, a/c in (0, 2].</summary>
    public static double NewmanRaju(double a, double c, double phi, double t, double b = BarB)
    {
        var r = a / c;
        var q = 1 + 1.464 * Math.Pow(r <= 1 ? r : 1 / r, 1.65);
        double m1, m2, m3, g, fphi;
        if (r <= 1)
        {
            m1 = 1.13 - 0.09 * r;
            m2 = -0.54 + 0.89 / (0.2 + r);
            m3 = 0.5 - 1 / (0.65 + r) + 14 * Math.Pow(1 - r, 24);
            g = 1 + (0.1 + 0.35 * (a / t) * (a / t)) * Math.Pow(1 - Math.Sin(phi), 2);
            fphi = Math.Pow(r * r * Math.Cos(phi) * Math.Cos(phi) + Math.Sin(phi) * Math.Sin(phi), 0.25);
        }
        else
        {
            m1 = Math.Sqrt(c / a) * (1 + 0.04 * c / a);
            m2 = 0.2 * Math.Pow(c / a, 4);
            m3 = -0.11 * Math.Pow(c / a, 4);
            g = 1 + (0.1 + 0.35 * (c / a) * (a / t) * (a / t)) * Math.Pow(1 - Math.Sin(phi), 2);
            fphi = Math.Pow(Math.Pow(c / a, 2) * Math.Sin(phi) * Math.Sin(phi) + Math.Cos(phi) * Math.Cos(phi), 0.25);
        }

        var fw = Math.Sqrt(1 / Math.Cos(Math.PI * c / (2 * b) * Math.Sqrt(a / t)));
        var f = (m1 + m2 * Math.Pow(a / t, 2) + m3 * Math.Pow(a / t, 4)) * g * fphi * fw;
        return Sigma * Math.Sqrt(Math.PI * a / q) * f;
    }

    /// <summary>Kuhn mesh of nx x ny x nz cubes of the bar [-1,1] x [-t,0] x [-h,h] (rows in Maubach order).</summary>
    private static (int[][] Rows, List<double> Coords, int NodeCount) BarBox(int nx, int ny, int nz, double thick, double half)
    {
        int[][] perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
        int Node(int i, int j, int k) => i + (nx + 1) * (j + (ny + 1) * k);
        var nN = (nx + 1) * (ny + 1) * (nz + 1);
        var coords = new List<double>(3 * nN);
        for (var k = 0; k <= nz; k++)
        for (var j = 0; j <= ny; j++)
        for (var i = 0; i <= nx; i++)
        {
            coords.Add(-1 + 2.0 * i / nx);
            coords.Add(-thick + thick * j / ny);
            coords.Add(-half + 2 * half * k / nz);
        }

        var rows = new int[6 * nx * ny * nz][];
        var t = 0;
        var c = new int[3];
        for (var k = 0; k < nz; k++)
        for (var j = 0; j < ny; j++)
        for (var i = 0; i < nx; i++)
            foreach (var perm in perms)
            {
                c[0] = i;
                c[1] = j;
                c[2] = k;
                var row = new int[4];
                row[0] = Node(c[0], c[1], c[2]);
                for (var s = 0; s < 3; s++)
                {
                    c[perm[s]]++;
                    row[s + 1] = Node(c[0], c[1], c[2]);
                }

                rows[t++] = row;
            }

        return (rows, coords, nN);
    }

    /// <summary>The bar refined by recursive bisection until no element is longer than 1.5 h(x) (as RefinedBox).</summary>
    private static (IncrementalMesh M, MeshState St, int Passes, long Bisections) RefinedBar(int nx, double thick, double half,
        double hmin, double rq, Func<double, double, double, double> distance)
    {
        var (rows, coords, nN) = BarBox(nx, (int)Math.Round(thick * nx / 2), (int)Math.Round(half * nx), thick, half);
        var flat = FlatBuild.Initial(ElementSpec.Tetrahedron, nN, rows);
        var m = new IncrementalMesh(flat, ElementSpec.Tetrahedron, persistent: true);
        var st = new MeshState(m, coords, null);
        var bis = new Bisection();
        const double grading = 0.35;
        var passes = 0;
        while (true)
        {
            var marked = new List<int>();
            for (var e = 0; e < m.ElementCount; e++)
            {
                var row = m.Row(e);
                if (row.Length != 4) continue;
                double longest = 0;
                for (var i = 0; i < 4; i++)
                for (var j = i + 1; j < 4; j++)
                {
                    double dx = st.X(row[i], 0) - st.X(row[j], 0), dy = st.X(row[i], 1) - st.X(row[j], 1),
                        dz = st.X(row[i], 2) - st.X(row[j], 2);
                    longest = Math.Max(longest, Math.Sqrt(dx * dx + dy * dy + dz * dz));
                }

                var (cx, cy, cz) = st.Centroid(row);
                if (longest > 1.5 * (hmin + grading * Math.Max(0, distance(cx, cy, cz) - rq))) marked.Add(e);
            }

            if (marked.Count == 0) break;
            bis.Refine(st, marked, maxLevel: 1000);
            passes++;
        }

        m.Commit();
        return (m, st, passes, bis.Bisections);
    }

    /// <summary>Nodal coordinates, tractions sigma on z = +-h and three corner constraints of the bar.</summary>
    private static (double[] X, double[] F, bool[] Fixed) BarProblem(Snapshot s, CrackInserter ins, int nx, double thick, double hh)
    {
        var nNodes = s.NodeCount;
        var X = new double[3 * nNodes];
        var xv = s.View(ins.X);
        var yv = s.View(ins.Y);
        var zv = s.View(ins.Z);
        for (var v = 0; v < nNodes; v++)
        {
            X[3 * v] = xv[v];
            X[3 * v + 1] = yv[v];
            X[3 * v + 2] = zv[v];
        }

        var f = new double[3 * nNodes];
        var t = new int[3];
        for (var i = 0; i < s.BoundaryCount; i++)
        {
            var id = s.BoundaryFacet(i);
            int e = id / 4, fl = id % 4;
            var row = s.Row(e);
            var k = 0;
            for (var p = 0; p < 4; p++)
                if (p != fl)
                    t[k++] = row[p];
            var zTop = t.All(v => Math.Abs(X[3 * v + 2] - hh) < 1e-12);
            var zBot = t.All(v => Math.Abs(X[3 * v + 2] + hh) < 1e-12);
            if (!zTop && !zBot) continue;
            var area = TriangleArea(X, t[0], t[1], t[2]);
            foreach (var v in t) f[3 * v + 2] += (zTop ? Sigma : -Sigma) * area / 3;
        }

        var ny = (int)Math.Round(thick * nx / 2);
        int Node(int i, int j, int k) => i + (nx + 1) * (j + (ny + 1) * k);
        var fixedDof = new bool[3 * nNodes];
        int cA = Node(0, 0, 0), cB = Node(nx, 0, 0), cC = Node(0, ny, 0);
        fixedDof[3 * cA] = fixedDof[3 * cA + 1] = fixedDof[3 * cA + 2] = true;
        fixedDof[3 * cB + 1] = fixedDof[3 * cB + 2] = true;
        fixedDof[3 * cC + 2] = true;
        return (X, f, fixedDof);
    }

    /// <summary>Parametric angle of the front point at polar angle psi of a front with semi-axes c (along p) and a.</summary>
    private static double Parametric(FrontCurve curve, double psi, double a, double c)
    {
        var (r, _, _) = curve.Eval(psi);
        return Math.Atan2(r * Math.Sin(psi) / a, r * Math.Cos(psi) / c);
    }

    public static void Surface(string[] arg)
    {
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        double D(int i, double d) => arg.Length > i && double.TryParse(arg[i], System.Globalization.NumberStyles.Float, ci, out var dv) ? dv : d;
        int n = arg.Length > 0 ? int.Parse(arg[0]) : 20, div = arg.Length > 1 ? int.Parse(arg[1]) : 8;
        double A0 = D(2, 0.15), C0 = D(3, 0.30);
        var maxSteps = arg.Length > 4 ? int.Parse(arg[4]) : 12;
        double daMax = D(5, 0.02), mExp = D(6, 3.0), tolK = D(7, 0.04);
        var rq = D(8, 0.05);
        var M = arg.Length > 9 && int.TryParse(arg[9], out var mk) ? mk : 48;
        var patch = arg.Contains("patch");
        double Opt(string key, double dflt) => arg.Where(x => x.StartsWith(key + "=")).Select(x => double.Parse(x[(key.Length + 1)..], ci)).DefaultIfEmpty(dflt).First();
        double tPl = Opt("t", 1.0), hPl = Opt("h", 2.0);
        var freeSurface = !arg.Contains("nofs"); // free-surface term of the interaction integral ("nofs": domain form only)
        var h0 = 2.0 / n;
        var hmin = h0 / div;
        var pl = new Plane(0, 0.37 * hmin, 0, -tPl);
        var daMin = hmin;
        var sw = System.Diagnostics.Stopwatch.StartNew();
        const double rOut = 0.52;
        var rIn = A0 - 0.01;

        double Region(double x, double y, double z)
        {
            var rho = Math.Sqrt(x * x + (y + tPl) * (y + tPl));
            var dr = Math.Max(0, Math.Max(rIn - rho, rho - rOut));
            var h = z - pl.Offset;
            return Math.Sqrt(dr * dr + h * h);
        }

        var (m, st, passes, bisections) = RefinedBar(n, tPl, hPl, hmin, rq, Region);
        var quality = Quality(m, st);
        Console.WriteLine($"  mesh: n={n}, hmin={hmin:F5}, {passes} bisection passes ({bisections} bisections), " +
                          $"{Live(m)} tetrahedra, {m.NodeCount} nodes; radius ratio min {quality.Min:F3}, mean {quality.Mean:F3}");
        var xyz = new double[3 * m.NodeCount];
        for (var v = 0; v < m.NodeCount; v++)
        for (var c = 0; c < 3; c++)
            xyz[3 * v + c] = st.X(v, c);
        var ins = new CrackInserter(m, xyz);
        m.CompactPattern(); // slots, hence destinations, unchanged

        var dir = DataDirOf();
        var tag = $"{n}_{div}" + (tPl != 1.0 || hPl != 2.0 ? $"_t{tPl}_h{hPl}" : "");
        using var wf = new StreamWriter(Path.Combine(dir, $"surface_fronts_{tag}.csv"));
        using var wk = new StreamWriter(Path.Combine(dir, $"surface_K_{tag}.csv"));
        using var ws = new StreamWriter(Path.Combine(dir, $"surface_steps_{tag}.csv"));
        wf.WriteLine("step,psi_deg,R,x,y,a,c");
        wk.WriteLine("step,psi_deg,phi_deg,KI,KII,KIII,KI_NR,a,c,J_over_JK");
        ws.WriteLine("step,attempt,da,accepted,err,t_fit_ms,t_insert_ms,t_publish_ms,t_rollback_ms,t_rebuild_ms," +
                     "t_assemble_ms,t_solve_ms,cg_iters,t_K_ms,elements,nodes,unknowns,moved,fit_cut,crack_edges," +
                     "split_nodes,front_nodes,visited,t_copy_ms,chunks_copied,chunks_total,f_dot_u");
        var fp = new FlatPattern();

        var curve = FrontCurve.Ellipse(C0, A0, M);
        FrontCurve? oldCurve = null;
        double[]? uPrev = null;
        double[]? kPrev = null;
        var da = daMax;
        Snapshot? prev = null;
        int[] prevFront = [];
        var accepted = 0;
        var attempt = 0;
        while (accepted <= maxSteps)
        {
            attempt++;
            var tFit = System.Diagnostics.Stopwatch.StartNew();
            var fit = Fit(m, ins, pl, curve, oldCurve, hmin, oldCurve == null ? null : ins.Front);
            tFit.Stop();
            var tIns = System.Diagnostics.Stopwatch.StartNew();
            CrackStep? step = null;
            if (!patch) step = Insert(ins, pl, curve, hmin, global: oldCurve == null);
            tIns.Stop();
            var tPub = System.Diagnostics.Stopwatch.StartNew();
            var s = m.Commit();
            tPub.Stop();
            if (arg.Contains("oracle"))
            {
                // every committed version equals the naive recomputation of all structures (also with the front on the free surface)
                var tO = System.Diagnostics.Stopwatch.StartNew();
                Oracle.Verify(s);
                Console.WriteLine($"  oracle: all structures of the committed version match ({tO.Elapsed.TotalSeconds:F1} s)");
            }

            if (accepted == 0 && attempt == 1)
                Console.WriteLine($"  front fitting: {fit.Moved} nodes moved, {fit.Cut} edges cut, {fit.Split} split elements" +
                                  (step == null ? "" : $"; crack: {step.Edges.Length} crack edges, {step.Splits.Length} split nodes, front {step.Front.Length} nodes"));

            // rebuild of every structure from the rows (the baseline the maintenance replaces), only on request
            var tReb = System.Diagnostics.Stopwatch.StartNew();
            if (arg.Contains("rebuild"))
            {
                var rows = new List<int[]>(s.ElementCount);
                for (var e = 0; e < s.ElementCount; e++) rows.Add(s.Row(e).ToArray());
                tReb.Restart();
                GC.KeepAlive(FlatBuild.Build(ElementSpec.Tetrahedron, s.NodeCount, rows, parallel: true));
            }

            tReb.Stop();

            var (X, f, fixedDof) = BarProblem(s, ins, n, tPl, hPl);
            var tCopy = System.Diagnostics.Stopwatch.StartNew();
            fp.Refresh(s);
            tCopy.Stop();
            var tAsm = System.Diagnostics.Stopwatch.StartNew();
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

            var tSol = System.Diagnostics.Stopwatch.StartNew();
            var (u, iters, rel) = Pcg(fp, vals, f, fixedDof, 1e-10, 60000, u0);
            tSol.Stop();
            if (!(rel <= 1e-8)) throw new Exception($"step {accepted}: CG did not converge ({iters} iterations, residual {rel:E2})");
            vals = null;
            double fu = 0;
            for (var i = 0; i < f.Length; i++) fu += f[i] * u[i];
            Console.WriteLine($"  solve: {3 * s.NodeCount} unknowns, assembly {tAsm.Elapsed.TotalSeconds:F2} s, CG {iters} iterations, " +
                              $"relative residual {rel:E1}, {tSol.Elapsed.TotalSeconds:F1} s");

            if (patch)
            {
                double err = 0, scale = 0;
                for (var v = 0; v < s.NodeCount; v++)
                {
                    if (s.TransposePairs(v).Length == 0) continue;
                    double ex = -Nu * Sigma / E * (X[3 * v] + 1), ey = -Nu * Sigma / E * (X[3 * v + 1] + tPl),
                        ez = Sigma / E * (X[3 * v + 2] + hPl);
                    err = Math.Max(err, Math.Max(Math.Abs(u[3 * v] - ex), Math.Max(Math.Abs(u[3 * v + 1] - ey),
                        Math.Abs(u[3 * v + 2] - ez))));
                    scale = Math.Max(scale, Math.Abs(ez));
                }

                Console.WriteLine($"surface patch test n={n}: max nodal error {err / scale:E2} relative to the exact linear field " +
                                  $"({sw.Elapsed.TotalSeconds:F1} s)");
                if (err / scale > 1e-7) throw new Exception("patch test failed");
                return;
            }

            var ext = arg.FirstOrDefault(x => x.StartsWith("u="));
            if (ext != null && accepted == 0)
            {
                var lines = File.ReadAllLines(ext[2..]);
                if (lines.Length != s.NodeCount) throw new Exception($"{ext}: {lines.Length} rows for {s.NodeCount} nodes");
                for (var v = 0; v < s.NodeCount; v++)
                {
                    var c = lines[v].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    for (var k = 0; k < 3; k++) u[3 * v + k] = double.Parse(c[k], ci);
                }

                Console.WriteLine($"  displacements read from {ext[2..]}");
            }

            var dump = arg.FirstOrDefault(x => x.StartsWith("dump="));
            if (dump != null && accepted == 0) DumpMesh(dump[5..], s, X, u, fixedDof, f);

            // K at the knots of the half of the front inside the body, mirrored about the plane psi = 0
            var tK = System.Diagnostics.Stopwatch.StartNew();
            var meanR = curve.Radii.Average();
            var delta = Math.Max(rq / meanR, 2 * Math.PI / 64);
            var K = new (double I, double II, double III)[M];
            var Jr = new double[M];
            for (var j = 0; j <= M / 2; j++)
            {
                var r = Interaction(s, X, u, pl, curve, 2 * Math.PI * j / M, rq, delta, 0, Math.PI, freeSurface);
                K[j] = (r.KI, r.KII, r.KIII);
                Jr[j] = r.J / ((1 - Nu * Nu) / E * (r.KI * r.KI + r.KII * r.KII) + r.KIII * r.KIII / (2 * Mu));
                if (j > 0 && j < M / 2) K[M - j] = (r.KI, -r.KII, -r.KIII);
            }

            tK.Stop();
            var kI = K.Select(t => t.I).ToArray();
            var kMean = kI.Average();
            double err2 = 0;
            if (kPrev != null)
            {
                var pMean = kPrev.Average();
                for (var j = 0; j < M; j++) err2 = Math.Max(err2, Math.Abs(kI[j] / kMean - kPrev[j] / pMean));
            }

            var accept = kPrev == null || err2 <= tolK || da <= daMin;
            double tRb = 0;
            if (!accept)
            {
                var tR = System.Diagnostics.Stopwatch.StartNew();
                m.Rollback(prev!);
                ins.RestoreFront(prevFront);
                tR.Stop();
                tRb = tR.Elapsed.TotalMilliseconds;
                if (arg.Contains("oracle"))
                {
                    Oracle.Verify(m.Current);
                    if (Oracle.Hash(m.Current) != Oracle.Hash(prev!)) throw new Exception("rollback did not restore the committed version");
                    Console.WriteLine("  oracle after rollback: structures match and equal the committed version");
                }
            }

            ws.WriteLine(string.Format(ci,
                "{0},{1},{2:R},{3},{4:F5},{5:F3},{6:F3},{7:F4},{8:F4},{9:F2},{10:F2},{11:F1},{12},{13:F1},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23:F3},{24},{25},{26:R}",
                accepted, attempt, da, accept ? 1 : 0, err2, tFit.Elapsed.TotalMilliseconds, tIns.Elapsed.TotalMilliseconds,
                tPub.Elapsed.TotalMilliseconds, tRb, tReb.Elapsed.TotalMilliseconds, tAsm.Elapsed.TotalMilliseconds,
                tSol.Elapsed.TotalMilliseconds, iters, tK.Elapsed.TotalMilliseconds, Live(m), s.NodeCount, 3 * s.NodeCount,
                fit.Moved, fit.Cut, step!.Edges.Length, step.Splits.Length, step.Front.Length, fit.Visited + step.VisitedTets,
                tCopy.Elapsed.TotalMilliseconds, fp.ChunksCopied, fp.ChunksTotal, fu));
            ws.Flush();
            var aD = curve.R(Math.PI / 2);
            var cD = curve.R(0);
            Console.WriteLine($"  step {accepted,2} attempt {attempt,2}: a {aD:F4} c {cD:F4} a/c {aD / cD:F3} a/t {aD / tPl:F3}, da {da:F4}, " +
                              $"err {err2:F4} -> {(accept ? "accepted" : "rolled back")}; fit {tFit.Elapsed.TotalMilliseconds:F1} ms, " +
                              $"insert {tIns.Elapsed.TotalMilliseconds:F1} ms, commit {tPub.Elapsed.TotalMilliseconds:F2} ms, " +
                              $"rollback {tRb:F3} ms, rebuild {tReb.Elapsed.TotalMilliseconds:F0} ms, assemble {tAsm.Elapsed.TotalMilliseconds:F0} ms, " +
                              $"CG {iters} it {tSol.Elapsed.TotalSeconds:F1} s, K {tK.Elapsed.TotalSeconds:F1} s; {Live(m)} tetrahedra, " +
                              $"{step.Edges.Length} crack edges, {step.Splits.Length} split nodes, front {step.Front.Length}");
            if (!accept)
            {
                da = Math.Max(daMin, da / 2);
                curve = Advance(oldCurve!, kPrev!, da, mExp);
                continue;
            }

            // accepted front and K against Newman and Raju; the end points (free surface) are reported apart
            double eAll = 0, eInt = 0, e2 = 0, e3 = 0, eDeep = 0, eSurf = 0;
            for (var j = 0; j < M; j++)
            {
                var psi = 2 * Math.PI * j / M;
                var (x, y, _, _, _, _) = curve.Point(psi);
                wf.WriteLine(string.Format(ci, "{0},{1:F4},{2:R},{3:R},{4:R},{5:R},{6:R}", accepted, psi * 180 / Math.PI, curve.Radii[j], x, y, aD, cD));
                var phi = Parametric(curve, psi, aD, cD);
                var nr = NewmanRaju(aD, cD, Math.Abs(phi), tPl);
                wk.WriteLine(string.Format(ci, "{0},{1:F4},{2:F4},{3:E6},{4:E6},{5:E6},{6:E6},{7:R},{8:R},{9:F4}", accepted, psi * 180 / Math.PI,
                    phi * 180 / Math.PI, K[j].I, K[j].II, K[j].III, nr, aD, cD, Jr[j <= M / 2 ? j : M - j]));
                if (j > M / 2) continue;
                var rel2 = Math.Abs(K[j].I - nr) / nr;
                eAll = Math.Max(eAll, rel2);
                if (j > 0 && j < M / 2) eInt = Math.Max(eInt, rel2);
                if (j == 0) eSurf = rel2;
                if (j == M / 4) eDeep = rel2;
                e2 = Math.Max(e2, Math.Abs(K[j].II) / nr);
                e3 = Math.Max(e3, Math.Abs(K[j].III) / nr);
            }

            Console.WriteLine($"    K_I against Newman-Raju: max error {eAll:P2} over the arc ({eInt:P2} without the surface points), " +
                              $"deepest point {eDeep:P2}, surface point {eSurf:P2}; max |K_II|/K {e2:P2}, |K_III|/K {e3:P2}");
            wf.Flush();
            wk.Flush();
            accepted++;
            prev = s;
            prevFront = (int[])ins.Front.Clone();
            uPrev = u;
            kPrev = kI;
            oldCurve = curve;
            da = Math.Min(daMax, 1.5 * da);
            if (curve.Radii.Max() + da > rOut - 0.02) break;
            curve = Advance(curve, kI, da, mExp);
        }

        Console.WriteLine($"surface n={n} div={div}: {accepted} accepted steps, {attempt} attempts ({sw.Elapsed.TotalSeconds:F0} s)");
    }
}

// Fracture.cs — stress intensity factors on the maintained structures: penny cracks and Paris-law growth.
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
//   penny n hmin_div beta_deg [a] [rq] [patch] [dump=<dir>] [u=<file>]
//   growth n hmin_div a b steps damax m tolK rq
//
// Body [-1,1]^3 (Kuhn box, spacing h0 = 2/n), E = 1, nu = 0.3, uniaxial tension sigma = 1 along z
// on the faces z = +-1, rigid motions removed by six point constraints at three corners. A planar
// crack lies in the plane n . x = c, n = (sin beta, 0, cos beta), with the star-shaped front
// rho = R(psi) of FrontCurve (circle for penny, ellipse and its Paris-law advances for growth).
//   mesh       recursive (Maubach) bisection until no element is longer than 1.5 h(x),
//              h = hmin + 0.35 max(0, d - rq), d the distance to the region of the fronts (hmin = h0 / div);
//   fitting    in a band of 3 hmin: nodes within 0.2 hmin of rho = R(psi) move onto it radially (no
//              neighbour moved, no element of the star inverted or below 10 % of its volume), the other
//              crossing edges are cut at the crossing (Cutter); after the first increment the nodes are
//              found by a walk from the previous front;
//   crack      CrackInserter, region rho < R(psi) + hmin / 4;
//   FE         P1 elasticity assembled through the destination table (3 x 3 blocks per pattern slot),
//              gathered row by row on all threads (bit-identical to the sequential element loop);
//              block-Jacobi CG to 1e-10 ||f|| on the flat copy of the pattern (FlatPattern, refreshed
//              per publication), warm start on stable node indices in growth;
//   K          domain interaction integral, global form with V = q N (N the in-plane front normal at the
//              foot point), Gosz-Moran auxiliary fields (2D gradient in the frame of the foot point,
//              traction-free asymptotic stress) with the non-equilibrium and incompatibility terms.
// References: Sneddon (beta = 0) K_I = 2 sigma sqrt(a/pi); Kassir and Sih (inclined), sigma_n = sigma
// cos^2 beta, tau = sigma sin beta cos beta, K_I = 2 sigma_n sqrt(a/pi), K_II = 4/(2 - nu) tau sqrt(a/pi)
// cos theta, K_III = 4 (1 - nu)/(2 - nu) tau sqrt(a/pi) sin theta; Irwin for the ellipse. With "patch" no
// crack is inserted and the solution must be the exact linear field; "dump" writes the mesh and the
// solution (SimPlas cross-check, figures); "u" replaces the solution by another solver's.

using System.Collections.Concurrent;
using Step2;

namespace Step2.Bench;

public static partial class Fracture
{
    private const double E = 1.0, Nu = 0.3, Sigma = 1.0;

    /// <summary>Elements per block of the interaction-integral sum (STEP2_KBLOCK); the contributing elements lie
    ///     near one front point, so small blocks spread them over the threads.</summary>
    private static readonly int KBlock = int.TryParse(Environment.GetEnvironmentVariable("STEP2_KBLOCK"), out var kb) ? kb : 64;
    /// <summary>Fraction of the tube radius with q = 1 (STEP2_QPLATEAU, default 0: q falls linearly from the front).
    ///     A plateau keeps the elements next to the front, where the finite element stress is least accurate, out of
    ///     the gradient terms of the domain integral.</summary>
    private static readonly double QPlateau = double.TryParse(Environment.GetEnvironmentVariable("STEP2_QPLATEAU"),
        System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out var qp0) ? qp0 : 0.0;

    private static readonly double Lambda = E * Nu / ((1 + Nu) * (1 - 2 * Nu)), Mu = E / (2 * (1 + Nu));

    /// <summary>
    ///     Crack plane n . x = offset: normal N and in-plane axes E1, E2 with E1 x E2 = N; Local gives the in-plane
    ///     coordinates (p, q) and the height h above the plane.
    /// </summary>
    private sealed class Plane(double beta, double offset = 0, double cx = 0, double cy = 0)
    {
        public readonly double[] N = [Math.Sin(beta), 0, Math.Cos(beta)];
        public readonly double[] E1 = [Math.Cos(beta), 0, -Math.Sin(beta)];
        public readonly double[] E2 = [0, 1, 0];
        public readonly double Offset = offset;

        // (cx, cy): global x and y of the in-plane origin (the centre of the front; on the face for a surface crack)
        public (double P, double Q, double H) Local(double x, double y, double z)
        {
            x -= cx;
            y -= cy;
            return (E1[0] * x + E1[1] * y + E1[2] * z, E2[0] * x + E2[1] * y + E2[2] * z,
                N[0] * x + N[1] * y + N[2] * z - Offset);
        }

        /// <summary>Global point of the local coordinates.</summary>
        public (double X, double Y, double Z) Global(double p, double q, double h) =>
            (p * E1[0] + q * E2[0] + (h + Offset) * N[0] + cx, p * E1[1] + q * E2[1] + (h + Offset) * N[1] + cy,
                p * E1[2] + q * E2[2] + (h + Offset) * N[2]);
    }

    /// <summary>Places the node of a cut edge at the stored root (parameter from the smaller label).</summary>
    private sealed class FitObserver(CrackInserter ins, Dictionary<long, double> roots) : ICutObserver
    {
        public void NodeCreated(int node, int a, int b)
        {
            int lo = Math.Min(a, b), hi = Math.Max(a, b);
            var t = roots[Cutter.Key(lo, hi)];
            ins.X.Set(node, ins.X[lo] + t * (ins.X[hi] - ins.X[lo]));
            ins.Y.Set(node, ins.Y[lo] + t * (ins.Y[hi] - ins.Y[lo]));
            ins.Z.Set(node, ins.Z[lo] + t * (ins.Z[hi] - ins.Z[lo]));
            ins.Parent1.Set(node, lo);
            ins.Parent2.Set(node, hi);
            ins.Weight.Set(node, t);
        }

        public void ElementSplit(int parent, int child) => ins.Origin.Set(child, ins.Origin[parent]);
    }

    /// <summary>
    ///     Kuhn box [-1,1]^3 refined by recursive (Maubach) bisection until no element is longer than
    ///     1.5 h(x), h = hmin + g max(0, d(x) - rq), d the distance to the region of the fronts.
    /// </summary>
    private static (IncrementalMesh M, MeshState St, int Passes, long Bisections) RefinedBox(int n, double hmin,
        double rq, Func<double, double, double, double> distance)
    {
        var (rows, coords01, nN) = Kuhn.Box(n);
        var coords = coords01.Select(c => 2 * c - 1).ToList();
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

    /// <summary>
    ///     Fits the mesh to the front curve in a band of 3 hmin around the crack plane: nodes within 0.2 hmin of
    ///     the surface rho = R(psi) move onto it radially (a node only if no neighbour moved), and the edges that
    ///     cross it farther than 0.1 of their length from both ends are cut at the crossing. Only nodes off the
    ///     crack surface take part; after the first increment they are found by a walk from the previous front
    ///     nodes through the band between 2 hmin behind it and 2 hmin beyond the new curve, so the cost follows
    ///     the advance.
    /// </summary>
    private static (int Moved, int Cut, int Split, int Visited) Fit(IncrementalMesh m, CrackInserter ins, Plane pl,
        FrontCurve curve, FrontCurve? old, double hmin, int[]? seeds)
    {
        var band = 3 * hmin;
        (double P, double Q, double H) L(int v) => pl.Local(ins.X[v], ins.Y[v], ins.Z[v]);
        double PhiR(int v)
        {
            var (p, q, _) = L(v);
            return Math.Sqrt(p * p + q * q) - curve.R(Math.Atan2(q, p));
        }

        bool Candidate(int v)
        {
            if (ins.OnSurface[v] != 0) return false;
            var (_, _, h) = L(v);
            return Math.Abs(h) < band;
        }

        // the walk stays between 2 hmin behind the previous front and 2 hmin beyond the new one
        bool InWalk(int v)
        {
            if (PhiR(v) >= 2 * hmin) return false;
            if (old == null) return true;
            var (p, q, _) = L(v);
            return Math.Sqrt(p * p + q * q) > old.R(Math.Atan2(q, p)) - 2 * hmin;
        }

        var nodes = new List<int>();
        if (seeds == null)
        {
            for (var v = 0; v < m.NodeCount; v++)
                if (m.TransposePairs(v).Length > 0 && Candidate(v))
                    nodes.Add(v);
        }
        else
        {
            var seen = new HashSet<int>();
            var stack = new Stack<int>();
            foreach (var s0 in seeds)
            foreach (var w in m.PatternColumns(s0))
                if (seen.Add(w) && Candidate(w) && InWalk(w))
                    stack.Push(w);
            while (stack.Count > 0)
            {
                var v = stack.Pop();
                nodes.Add(v);
                foreach (var w in m.PatternColumns(v))
                    if (seen.Add(w) && Candidate(w) && InWalk(w))
                        stack.Push(w);
            }

            nodes.Sort();
        }

        // signed volume of element e with node v at (vx, vy, vz)
        double Volume(int e, int v, double vx, double vy, double vz)
        {
            var r = m.Row(e);
            Span<double> c = stackalloc double[12];
            for (var i = 0; i < 4; i++)
            {
                var w = r[i];
                c[3 * i] = w == v ? vx : ins.X[w];
                c[3 * i + 1] = w == v ? vy : ins.Y[w];
                c[3 * i + 2] = w == v ? vz : ins.Z[w];
            }

            double ax = c[3] - c[0], ay = c[4] - c[1], az = c[5] - c[2];
            double bx = c[6] - c[0], by = c[7] - c[1], bz = c[8] - c[2];
            double cx = c[9] - c[0], cy = c[10] - c[1], cz = c[11] - c[2];
            return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
        }

        var moved = new HashSet<int>();
        foreach (var v in nodes)
        {
            if (Math.Abs(PhiR(v)) >= 0.2 * hmin) continue;
            var nb = false;
            foreach (var w in m.PatternColumns(v)) nb |= moved.Contains(w);
            if (nb) continue;
            var (p, q, h) = L(v);
            var f = curve.R(Math.Atan2(q, p)) / Math.Sqrt(p * p + q * q);
            var (nx, ny, nz) = pl.Global(f * p, f * q, h);
            // no element of the star may change orientation or lose more than 90 % of its volume
            var ok = true;
            var tr = m.TransposePairs(v);
            for (var t = 0; t < tr.Length && ok; t += 2)
            {
                var e = tr[t];
                var v0 = Volume(e, v, ins.X[v], ins.Y[v], ins.Z[v]);
                var v1 = Volume(e, v, nx, ny, nz);
                ok = v0 * v1 > 0 && Math.Abs(v1) > 0.1 * Math.Abs(v0);
            }

            if (!ok) continue;
            moved.Add(v);
            ins.X.Set(v, nx);
            ins.Y.Set(v, ny);
            ins.Z.Set(v, nz);
        }

        CrackInserter.Field fr = (x, y, z) =>
        {
            var (p, q, _) = pl.Local(x, y, z);
            return Math.Sqrt(p * p + q * q) - curve.R(Math.Atan2(q, p));
        };
        var inSet = new HashSet<int>(nodes);
        var roots = new Dictionary<long, double>(EdgeKeyComparer.Instance);
        var edges = new List<(int, int)>();
        foreach (var v in nodes)
        {
            var cols = m.PatternColumns(v);
            var cnts = m.PatternCounts(v);
            for (var t = 0; t < cols.Length; t++)
            {
                var w = cols[t];
                if (cnts[t] <= 0 || w <= v || !inSet.Contains(w)) continue;
                double r0 = PhiR(v), r1 = PhiR(w);
                if (!((r0 < -1e-12 && r1 > 1e-12) || (r0 > 1e-12 && r1 < -1e-12))) continue;
                var tl = r0 / (r0 - r1);
                if (tl < 0.1 || tl > 0.9) continue;
                roots[Cutter.Key(v, w)] = CrackInserter.RootOf(ins.X[v], ins.Y[v], ins.Z[v], ins.X[w], ins.Y[w], ins.Z[w], fr);
                edges.Add((v, w));
            }
        }

        edges.Sort();
        var rec = Cutter.Cut(m, edges, new FitObserver(ins, roots));
        return (moved.Count, rec.Edges.Length, rec.ShellElements, nodes.Count);
    }

    /// <summary>One crack increment up to the curve (nodes within 0.25 hmin outside it count as inside).</summary>
    private static CrackStep Insert(CrackInserter ins, Plane pl, FrontCurve curve, double hmin, bool global)
    {
        CrackInserter.Field phiS = (x, y, z) => pl.Local(x, y, z).H;
        CrackInserter.Field phiR = (x, y, z) =>
        {
            var (p, q, _) = pl.Local(x, y, z);
            return Math.Sqrt(p * p + q * q) - curve.R(Math.Atan2(q, p)) - 0.25 * hmin;
        };
        return ins.Advance(phiS, phiR, global);
    }

    /// <summary>Nodal coordinates, tractions sigma on z = +-1 and the six point constraints.</summary>
    private static (double[] X, double[] F, bool[] Fixed) Problem(Snapshot s, CrackInserter ins, int n)
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
            var zTop = t.All(v => Math.Abs(X[3 * v + 2] - 1) < 1e-12);
            var zBot = t.All(v => Math.Abs(X[3 * v + 2] + 1) < 1e-12);
            if (!zTop && !zBot) continue;
            var area = TriangleArea(X, t[0], t[1], t[2]);
            foreach (var v in t) f[3 * v + 2] += (zTop ? Sigma : -Sigma) * area / 3;
        }

        var fixedDof = new bool[3 * nNodes];
        int cA = Kuhn.Node(n, 0, 0, 0), cB = Kuhn.Node(n, n, 0, 0), cC = Kuhn.Node(n, 0, n, 0);
        fixedDof[3 * cA] = fixedDof[3 * cA + 1] = fixedDof[3 * cA + 2] = true;
        fixedDof[3 * cB + 1] = fixedDof[3 * cB + 2] = true;
        fixedDof[3 * cC + 2] = true;
        return (X, f, fixedDof);
    }

    public static void Run(string[] a)
    {
        int n = int.Parse(a[0]), div = a.Length > 1 ? int.Parse(a[1]) : 8;
        var beta = (a.Length > 2 ? double.Parse(a[2]) : 0.0) * Math.PI / 180;
        var R = a.Length > 3 ? double.Parse(a[3]) : 0.2;
        var h0 = 2.0 / n;
        var hmin = h0 / div;
        var rq = a.Length > 4 ? double.Parse(a[4]) : Math.Max(0.2 * R, 6 * hmin);
        var patch = a.Contains("patch");
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var pl = new Plane(beta);
        var curve = FrontCurve.Circle(R, 64);

        double FrontDistance(double x, double y, double z)
        {
            var (p, q, h) = pl.Local(x, y, z);
            var rho = Math.Sqrt(p * p + q * q);
            return Math.Sqrt((rho - R) * (rho - R) + h * h);
        }

        var (m, st, passes, bisections) = RefinedBox(n, hmin, rq, FrontDistance);
        var quality = Quality(m, st);
        Console.WriteLine($"  mesh: n={n} (h0={h0:F4}), hmin={hmin:F5}, {passes} bisection passes ({bisections} " +
                          $"bisections), {Live(m)} tetrahedra, {m.NodeCount} nodes; radius ratio min {quality.Min:F3}, " +
                          $"mean {quality.Mean:F3}");

        var xyz = new double[3 * m.NodeCount];
        for (var v = 0; v < m.NodeCount; v++)
        for (var c = 0; c < 3; c++)
            xyz[3 * v + c] = st.X(v, c);
        var ins = new CrackInserter(m, xyz);
        var fit = Fit(m, ins, pl, curve, null, hmin, null);
        Console.WriteLine($"  front fitting: {fit.Moved} nodes moved onto the cylinder, {fit.Cut} edges cut, " +
                          $"{fit.Split} split elements");

        CrackStep? step = null;
        if (!patch)
        {
            step = Insert(ins, pl, curve, hmin, global: true);
            Console.WriteLine($"  crack: a={R}, beta={beta * 180 / Math.PI:F1} deg, {step.Edges.Length} crack edges, " +
                              $"{step.Splits.Length} split nodes, front {step.Front.Length} nodes");
        }

        // one solve: the pattern is compacted first (slots, and hence destinations, unchanged), so that the
        // 3 x 3 blocks of the abandoned segments of the refinement do not count against the array size
        m.CompactPattern();
        var s = m.Commit();
        var nNodes = s.NodeCount;
        var (X, f, fixedDof) = Problem(s, ins, n);
        var fp = new FlatPattern();
        var tCopy = System.Diagnostics.Stopwatch.StartNew();
        fp.Refresh(s);
        tCopy.Stop();
        var tAsm = System.Diagnostics.Stopwatch.StartNew();
        var vals = Assemble(s, fp, X);
        tAsm.Stop();
        var tSol = System.Diagnostics.Stopwatch.StartNew();
        var (u, iters, rel) = Pcg(fp, vals, f, fixedDof, 1e-10, 20000, null);
        tSol.Stop();
        Console.WriteLine($"  solve: {3 * nNodes} unknowns, pattern copy {tCopy.Elapsed.TotalMilliseconds:F1} ms, assembly " +
                          $"{tAsm.Elapsed.TotalSeconds:F2} s, CG {iters} iterations, relative residual {rel:E1}, " +
                          $"{tSol.Elapsed.TotalSeconds:F1} s");

        if (patch)
        {
            double err = 0, scale = 0;
            for (var v = 0; v < nNodes; v++)
            {
                if (s.TransposePairs(v).Length == 0) continue;
                double ex = -Nu * Sigma / E * (X[3 * v] + 1), ey = -Nu * Sigma / E * (X[3 * v + 1] + 1),
                    ez = Sigma / E * (X[3 * v + 2] + 1);
                err = Math.Max(err, Math.Max(Math.Abs(u[3 * v] - ex), Math.Max(Math.Abs(u[3 * v + 1] - ey),
                    Math.Abs(u[3 * v + 2] - ez))));
                scale = Math.Max(scale, Math.Abs(ez));
            }

            Console.WriteLine($"penny patch test n={n}: max nodal error {err / scale:E2} relative to the exact linear " +
                              $"field ({sw.Elapsed.TotalSeconds:F1} s)");
            if (err / scale > 1e-7) throw new Exception("patch test failed");
            return;
        }

        var front = step!.Front.Select(v =>
        {
            var (p, q, _) = pl.Local(X[3 * v], X[3 * v + 1], X[3 * v + 2]);
            return Math.Sqrt(p * p + q * q);
        }).ToArray();
        Console.WriteLine($"  front: {front.Length} nodes, radius min {front.Min():F5}, mean {front.Average():F5}, " +
                          $"max {front.Max():F5}");

        // optional: the displacements of another solver on the same mesh (SimPlas), as "u=<file>"
        var ext = a.FirstOrDefault(x => x.StartsWith("u="));
        if (ext != null)
        {
            var lines = File.ReadAllLines(ext[2..]);
            if (lines.Length != nNodes) throw new Exception($"{ext}: {lines.Length} rows for {nNodes} nodes");
            for (var v = 0; v < nNodes; v++)
            {
                var c = lines[v].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                for (var k = 0; k < 3; k++) u[3 * v + k] = double.Parse(c[k], System.Globalization.CultureInfo.InvariantCulture);
            }

            Console.WriteLine($"  displacements read from {ext[2..]}");
        }

        var dump = a.FirstOrDefault(x => x.StartsWith("dump="));
        if (dump != null) DumpMesh(dump[5..], s, X, u, fixedDof, f);

        const int points = 16;
        var delta = Math.Max(rq / R, 2 * Math.PI / 64);
        var kRef = 2 * Sigma * Math.Sqrt(R / Math.PI);
        double sn = Sigma * Math.Cos(beta) * Math.Cos(beta), tau = Sigma * Math.Sin(beta) * Math.Cos(beta);
        var tag = ext != null ? "_ext" : "";
        using var w = new StreamWriter(Path.Combine(DataDirOf(),
            string.Format(System.Globalization.CultureInfo.InvariantCulture, "penny_{0}_{1}_{2}_rq{3:F3}{4}.csv", n, div,
                (int)Math.Round(beta * 180 / Math.PI), rq, tag)));
        w.WriteLine("theta_deg,KI,KII,KIII,KI_ref,KII_ref,KIII_ref,J,J_from_K");
        double errI = 0, errII = 0, errIII = 0;
        var tK = System.Diagnostics.Stopwatch.StartNew();
        for (var k = 0; k < points; k++)
        {
            var th = 2 * Math.PI * k / points;
            var (kI, kII, kIII, J) = Interaction(s, X, u, pl, curve, th, rq, delta);
            var refI = 2 * sn * Math.Sqrt(R / Math.PI);
            var refII = 4 / (2 - Nu) * tau * Math.Sqrt(R / Math.PI) * Math.Cos(th);
            var refIII = 4 * (1 - Nu) / (2 - Nu) * tau * Math.Sqrt(R / Math.PI) * Math.Sin(th);
            var jK = (1 - Nu * Nu) / E * (kI * kI + kII * kII) + kIII * kIII / (2 * Mu);
            w.WriteLine($"{th * 180 / Math.PI:F2},{kI:E6},{kII:E6},{kIII:E6},{refI:E6},{refII:E6},{refIII:E6},{J:E6},{jK:E6}");
            Console.WriteLine($"  theta {th * 180 / Math.PI,6:F1}: K_I/K0 {kI / kRef,8:F4} (ref {refI / kRef,7:F4})  " +
                              $"K_II/K0 {kII / kRef,8:F4} (ref {refII / kRef,7:F4})  K_III/K0 {kIII / kRef,8:F4} " +
                              $"(ref {refIII / kRef,7:F4})  J/J(K) {J / jK:F4}");
            errI = Math.Max(errI, Math.Abs(kI - refI) / kRef);
            errII = Math.Max(errII, Math.Abs(Math.Abs(kII) - Math.Abs(refII)) / kRef);
            errIII = Math.Max(errIII, Math.Abs(Math.Abs(kIII) - Math.Abs(refIII)) / kRef);
        }

        Console.WriteLine($"  K at {points} points: {tK.Elapsed.TotalSeconds:F2} s");
        Console.WriteLine($"penny n={n} div={div} beta={beta * 180 / Math.PI:F0}{tag}: max |K - K_ref| / K0 = {errI:P2} (I), " +
                          $"{errII:P2} (|II|), {errIII:P2} (|III|); rq={rq:F4}, delta={delta * 180 / Math.PI:F1} deg, " +
                          $"{Live(m)} tetrahedra, {3 * nNodes} unknowns ({sw.Elapsed.TotalSeconds:F1} s)");
    }

    /// <summary>
    ///     Writes the mesh and the solution for an independent solver and the figures: nodes (x y z), tetrahedra,
    ///     boundary triangles with their element, fixed DOFs, nodal loads and displacements.
    /// </summary>
    private static void DumpMesh(string dir, Snapshot s, double[] X, double[] u, bool[] fixedDof, double[] f)
    {
        Directory.CreateDirectory(dir);
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        using (var w = new StreamWriter(Path.Combine(dir, "nodes.txt")))
            for (var v = 0; v < s.NodeCount; v++)
                w.WriteLine(string.Format(ci, "{0:R} {1:R} {2:R}", X[3 * v], X[3 * v + 1], X[3 * v + 2]));
        using (var w = new StreamWriter(Path.Combine(dir, "tets.txt")))
            for (var e = 0; e < s.ElementCount; e++)
            {
                var r = s.Row(e);
                if (r.Length == 4) w.WriteLine($"{r[0]} {r[1]} {r[2]} {r[3]}");
            }

        using (var w = new StreamWriter(Path.Combine(dir, "boundary.txt")))
            for (var i = 0; i < s.BoundaryCount; i++)
            {
                var id = s.BoundaryFacet(i);
                int e = id / 4, fl = id % 4;
                var r = s.Row(e).ToArray();
                var t = Enumerable.Range(0, 4).Where(p => p != fl).Select(p => r[p]).ToArray();
                w.WriteLine($"{t[0]} {t[1]} {t[2]} {e}");
            }

        using (var w = new StreamWriter(Path.Combine(dir, "fixed.txt")))
            for (var d = 0; d < fixedDof.Length; d++)
                if (fixedDof[d])
                    w.WriteLine($"{d / 3} {d % 3}");
        using (var w = new StreamWriter(Path.Combine(dir, "loads.txt")))
            for (var d = 0; d < f.Length; d++)
                if (f[d] != 0)
                    w.WriteLine(string.Format(ci, "{0} {1} {2:R}", d / 3, d % 3, f[d]));
        using (var w = new StreamWriter(Path.Combine(dir, "u.txt")))
            for (var v = 0; v < s.NodeCount; v++)
                w.WriteLine(string.Format(ci, "{0:R} {1:R} {2:R}", u[3 * v], u[3 * v + 1], u[3 * v + 2]));
    }

    /// <summary>
    ///     Fatigue growth of an elliptical crack (semi-axes a along x, b along y, plane z = 0): Irwin's K_I at the
    ///     first front, then Paris-law increments da_j = da (K_j / K_max)^m along the front normal, each a local
    ///     fit, crack increment, publication, warm-started solve and K evaluation. A step whose normalized K
    ///     changes by more than tolK is rolled back to the committed version before it and retried with da / 2.
    ///
    ///   growth n div a b steps damax m tolK rq
    /// </summary>
    public static void Growth(string[] a)
    {
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        double D(int i, double d) => a.Length > i ? double.Parse(a[i], ci) : d;
        int n = a.Length > 0 ? int.Parse(a[0]) : 20, div = a.Length > 1 ? int.Parse(a[1]) : 8;
        double A0 = D(2, 0.25), B0 = D(3, 0.15);
        var maxSteps = a.Length > 4 ? int.Parse(a[4]) : 12;
        double daMax = D(5, 0.03), mExp = D(6, 4.0), tolK = D(7, 0.02);
        var h0 = 2.0 / n;
        var hmin = h0 / div;
        var rq = D(8, 0.06);
        const int M = 32;
        double rIn = B0 - 0.01, rOut = 0.40;
        // the plane z = 0.37 hmin is not a plane of the lattice, so the increments cut crack edges
        var pl = new Plane(0, 0.37 * hmin);
        var daMin = hmin;
        var sw = System.Diagnostics.Stopwatch.StartNew();
        double Region(double x, double y, double z)
        {
            var rho = Math.Sqrt(x * x + y * y);
            var dr = Math.Max(0, Math.Max(rIn - rho, rho - rOut));
            var h = z - pl.Offset;
            return Math.Sqrt(dr * dr + h * h);
        }

        var (m, st, passes, bisections) = RefinedBox(n, hmin, rq, Region);
        var quality = Quality(m, st);
        Console.WriteLine($"  mesh: n={n}, hmin={hmin:F5}, {passes} bisection passes ({bisections} bisections), " +
                          $"{Live(m)} tetrahedra, {m.NodeCount} nodes; radius ratio min {quality.Min:F3}, mean {quality.Mean:F3}");
        var xyz = new double[3 * m.NodeCount];
        for (var v = 0; v < m.NodeCount; v++)
        for (var c = 0; c < 3; c++)
            xyz[3 * v + c] = st.X(v, c);
        var ins = new CrackInserter(m, xyz);

        var dir = DataDirOf();
        using var wf = new StreamWriter(Path.Combine(dir, "growth_fronts.csv"));
        using var wk = new StreamWriter(Path.Combine(dir, "growth_K.csv"));
        using var ws = new StreamWriter(Path.Combine(dir, "growth_steps.csv"));
        wf.WriteLine("step,theta_deg,R,x,y");
        wk.WriteLine("step,theta_deg,KI,KII,KIII,KI_irwin");
        ws.WriteLine("step,attempt,da,accepted,err,t_fit_ms,t_insert_ms,t_publish_ms,t_rollback_ms,t_rebuild_ms," +
                     "t_assemble_ms,t_solve_ms,cg_iters,t_K_ms,elements,nodes,unknowns,moved,fit_cut,crack_edges," +
                     "split_nodes,front_nodes,visited,t_copy_ms,chunks_copied,chunks_total");
        var fp = new FlatPattern();

        var curve = FrontCurve.Ellipse(A0, B0, M);
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
            var step = Insert(ins, pl, curve, hmin, global: oldCurve == null);
            tIns.Stop();
            var tPub = System.Diagnostics.Stopwatch.StartNew();
            var s = m.Commit();
            tPub.Stop();
            if (Environment.GetEnvironmentVariable("STEP2_FRONTDBG") == "1")
            {
                // where the front nodes lie relative to the curve: d = rho - R(psi) in units of hmin
                var bins = new int[6];
                foreach (var v in step.Front)
                {
                    var (p, q, h) = pl.Local(ins.X[v], ins.Y[v], ins.Z[v]);
                    var d = (Math.Sqrt(p * p + q * q) - curve.R(Math.Atan2(q, p))) / hmin;
                    bins[d < -2 ? 0 : d < -1 ? 1 : d < -0.5 ? 2 : d < 0.5 ? 3 : d < 1 ? 4 : 5]++;
                }

                Console.WriteLine($"    [front] {step.Front.Length} nodes by (rho - R)/hmin: <-2 {bins[0]}, -2..-1 {bins[1]}, " +
                                  $"-1..-0.5 {bins[2]}, -0.5..0.5 {bins[3]}, 0.5..1 {bins[4]}, >1 {bins[5]}");
            }

            // rebuild of every structure from the rows (the baseline the maintenance replaces)
            var rows = new List<int[]>(s.ElementCount);
            for (var e = 0; e < s.ElementCount; e++) rows.Add(s.Row(e).ToArray());
            var tReb = System.Diagnostics.Stopwatch.StartNew();
            GC.KeepAlive(FlatBuild.Build(ElementSpec.Tetrahedron, s.NodeCount, rows, parallel: true));
            tReb.Stop();
            rows = null;

            var (X, f, fixedDof) = Problem(s, ins, n);
            var tCopy = System.Diagnostics.Stopwatch.StartNew();
            fp.Refresh(s); // copies the chunks the step changed
            tCopy.Stop();
            var tAsm = System.Diagnostics.Stopwatch.StartNew();
            var vals = Assemble(s, fp, X);
            tAsm.Stop();
            // warm start on the stable node indices; new nodes from their parents
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
            var (u, iters, rel) = Pcg(fp, vals, f, fixedDof, 1e-10, 40000, u0);
            tSol.Stop();
            if (!(rel <= 1e-8)) throw new Exception($"step {accepted}: CG did not converge ({iters} iterations, residual {rel:E2})");
            vals = null;

            var tK = System.Diagnostics.Stopwatch.StartNew();
            var meanR = curve.Radii.Average();
            var delta = Math.Max(rq / meanR, 2 * Math.PI / 64);
            var K = new (double I, double II, double III)[M];
            for (var j = 0; j < M; j++)
            {
                var r = Interaction(s, X, u, pl, curve, 2 * Math.PI * j / M, rq, delta);
                K[j] = (r.KI, r.KII, r.KIII);
            }

            tK.Stop();
            var kI = K.Select(t => t.I).ToArray();
            var kMean = kI.Average();
            double err = 0;
            if (kPrev != null)
            {
                var pMean = kPrev.Average();
                for (var j = 0; j < M; j++) err = Math.Max(err, Math.Abs(kI[j] / kMean - kPrev[j] / pMean));
            }

            var accept = kPrev == null || err <= tolK || da <= daMin;
            double tRb = 0;
            if (!accept)
            {
                var tR = System.Diagnostics.Stopwatch.StartNew();
                m.Rollback(prev!);
                ins.RestoreFront(prevFront);
                tR.Stop();
                tRb = tR.Elapsed.TotalMilliseconds;
            }

            ws.WriteLine(string.Format(ci,
                "{0},{1},{2:R},{3},{4:F5},{5:F3},{6:F3},{7:F4},{8:F4},{9:F2},{10:F2},{11:F1},{12},{13:F1},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23:F4},{24},{25}",
                accepted, attempt, da, accept ? 1 : 0, err, tFit.Elapsed.TotalMilliseconds, tIns.Elapsed.TotalMilliseconds,
                tPub.Elapsed.TotalMilliseconds, tRb, tReb.Elapsed.TotalMilliseconds, tAsm.Elapsed.TotalMilliseconds,
                tSol.Elapsed.TotalMilliseconds, iters, tK.Elapsed.TotalMilliseconds, Live(m), s.NodeCount, 3 * s.NodeCount,
                fit.Moved, fit.Cut, step.Edges.Length, step.Splits.Length, step.Front.Length, fit.Visited + step.VisitedTets,
                tCopy.Elapsed.TotalMilliseconds, fp.ChunksCopied, fp.ChunksTotal));
            ws.Flush();
            Console.WriteLine($"  step {accepted,2} attempt {attempt,2}: da {da:F4}, err {err:F4} -> {(accept ? "accepted" : "rolled back")}; " +
                              $"fit {tFit.Elapsed.TotalMilliseconds:F1} ms, insert {tIns.Elapsed.TotalMilliseconds:F1} ms, " +
                              $"commit {tPub.Elapsed.TotalMilliseconds:F2} ms, rollback {tRb:F3} ms, rebuild {tReb.Elapsed.TotalMilliseconds:F0} ms, " +
                              $"pattern copy {tCopy.Elapsed.TotalMilliseconds:F2} ms ({fp.ChunksCopied}/{fp.ChunksTotal} chunks), " +
                              $"assemble {tAsm.Elapsed.TotalMilliseconds:F0} ms, CG {iters} it {tSol.Elapsed.TotalSeconds:F1} s, K {tK.Elapsed.TotalSeconds:F1} s; " +
                              $"{Live(m)} tetrahedra, {step.Edges.Length} crack edges, {step.Splits.Length} split nodes, front {step.Front.Length}");
            if (!accept)
            {
                da = Math.Max(daMin, da / 2);
                curve = Advance(oldCurve!, kPrev!, da, mExp);
                continue;
            }

            // record the accepted front and K
            for (var j = 0; j < M; j++)
            {
                var th = 2 * Math.PI * j / M;
                var (x, y, _, _, _, _) = curve.Point(th);
                wf.WriteLine(string.Format(ci, "{0},{1:F4},{2:R},{3:R},{4:R}", accepted, th * 180 / Math.PI, curve.Radii[j], x, y));
                var irwin = accepted == 0 ? Irwin(A0, B0, th) : double.NaN;
                wk.WriteLine(string.Format(ci, "{0},{1:F4},{2:E6},{3:E6},{4:E6},{5:E6}", accepted, th * 180 / Math.PI,
                    K[j].I, K[j].II, K[j].III, irwin));
            }

            wf.Flush();
            wk.Flush();
            if (accepted == 0)
            {
                double e0 = 0;
                for (var j = 0; j < M; j++)
                {
                    var r = Irwin(A0, B0, 2 * Math.PI * j / M);
                    e0 = Math.Max(e0, Math.Abs(K[j].I - r) / r);
                }

                Console.WriteLine($"  initial front: max |K_I - K_Irwin| / K_Irwin = {e0:P2}");
            }

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

        Console.WriteLine($"growth n={n} div={div}: {accepted} accepted steps, {attempt} attempts ({sw.Elapsed.TotalSeconds:F0} s)");
    }

    /// <summary>Paris-law advance of the knots: da_j = da (K_j / K_max)^m along the normal, as a radial increment.</summary>
    private static FrontCurve Advance(FrontCurve c, double[] k, double da, double mExp)
    {
        var kMax = k.Max();
        var r = new double[c.M];
        for (var j = 0; j < c.M; j++)
        {
            var th = 2 * Math.PI * j / c.M;
            var (nx, ny, _) = c.Normal(th);
            var cosA = nx * Math.Cos(th) + ny * Math.Sin(th);
            r[j] = c.Radii[j] + da * Math.Pow(k[j] / kMax, mExp) / cosA;
        }

        return new FrontCurve(r);
    }

    /// <summary>Irwin's K_I of an elliptical crack (semi-axes a > b) under remote tension Sigma, at polar angle psi.</summary>
    private static double Irwin(double a, double b, double psi)
    {
        var phi = Math.Atan2(a * Math.Sin(psi), b * Math.Cos(psi)); // parametric angle: (a cos phi, b sin phi)
        var k2 = 1 - b * b / (a * a);
        double ek = 0;
        const int nq = 4000;
        for (var i = 0; i <= nq; i++)
        {
            var t = 0.5 * Math.PI * i / nq;
            ek += (i == 0 || i == nq ? 0.5 : 1.0) * Math.Sqrt(1 - k2 * Math.Sin(t) * Math.Sin(t));
        }

        ek *= 0.5 * Math.PI / nq;
        return Sigma * Math.Sqrt(Math.PI * b) / ek *
               Math.Pow(Math.Sin(phi) * Math.Sin(phi) + b * b / (a * a) * Math.Cos(phi) * Math.Cos(phi), 0.25);
    }

    private static string DataDirOf() =>
        Environment.GetEnvironmentVariable("STEP2_DATA") ??
        Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "../../../../../benchmarks/data"));

    private static int Live(IncrementalMesh m)
    {
        var c = 0;
        for (var e = 0; e < m.ElementCount; e++)
            if (m.Row(e).Length == 4)
                c++;
        return c;
    }

    private static double TriangleArea(double[] X, int a, int b, int c)
    {
        double ux = X[3 * b] - X[3 * a], uy = X[3 * b + 1] - X[3 * a + 1], uz = X[3 * b + 2] - X[3 * a + 2];
        double vx = X[3 * c] - X[3 * a], vy = X[3 * c + 1] - X[3 * a + 1], vz = X[3 * c + 2] - X[3 * a + 2];
        double cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx;
        return 0.5 * Math.Sqrt(cx * cx + cy * cy + cz * cz);
    }

    /// <summary>Gradients of the barycentric coordinates (g[3a + c]) and the volume of a tetrahedron.</summary>
    private static double Gradients(double[] X, ReadOnlySpan<int> row, Span<double> g)
    {
        Span<double> j = stackalloc double[9];
        for (var c = 0; c < 3; c++)
        {
            j[c * 3 + 0] = X[3 * row[1] + c] - X[3 * row[0] + c];
            j[c * 3 + 1] = X[3 * row[2] + c] - X[3 * row[0] + c];
            j[c * 3 + 2] = X[3 * row[3] + c] - X[3 * row[0] + c];
        }

        var det = j[0] * (j[4] * j[8] - j[5] * j[7]) - j[1] * (j[3] * j[8] - j[5] * j[6]) +
                  j[2] * (j[3] * j[7] - j[4] * j[6]);
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
        for (var c = 0; c < 3; c++)
        {
            g[3 + c] = inv[0 * 3 + c];
            g[6 + c] = inv[1 * 3 + c];
            g[9 + c] = inv[2 * 3 + c];
            g[c] = -(g[3 + c] + g[6 + c] + g[9 + c]);
        }

        return Math.Abs(det) / 6.0;
    }

    /// <summary>
    ///     Stiffness into 3 x 3 blocks per pattern slot, addressed by the destination table and gathered row by
    ///     row on all threads: row v sums the blocks of the pairs (e, a) of T(v) in element order, the order of a
    ///     sequential loop over the elements, so every value equals that loop's bit for bit.
    /// </summary>
    private static double[] Assemble(Snapshot s, FlatPattern fp, double[] X)
    {
        if (9L * s.TotalSlots > Array.MaxLength)
            throw new InvalidOperationException($"{s.TotalSlots} pattern slots: 3 x 3 blocks exceed the array size");
        var vals = new double[9 * s.TotalSlots];
        var rowStart = fp.RowStart;
        Parallel.ForEach(Partitioner.Create(0, s.NodeCount, 1024), r =>
        {
            Span<double> g = stackalloc double[12];
            for (var v = r.Item1; v < r.Item2; v++)
            {
                var tr = s.TransposePairs(v);
                var rs = rowStart[v];
                for (var q = 0; q < tr.Length; q += 2)
                {
                    var e = tr[q];
                    var row = s.Row(e);
                    if (row.Length != 4) continue;
                    var a = tr[q + 1];
                    var vol = Gradients(X, row, g);
                    var dest = s.Destinations(e);
                    for (var b = 0; b < 4; b++)
                    {
                        var at = 9 * (rs + dest[4 * a + b]);
                        var gg = g[3 * a] * g[3 * b] + g[3 * a + 1] * g[3 * b + 1] + g[3 * a + 2] * g[3 * b + 2];
                        for (var i = 0; i < 3; i++)
                        for (var jj = 0; jj < 3; jj++)
                            vals[at + 3 * i + jj] += vol * (Lambda * g[3 * a + i] * g[3 * b + jj] +
                                                            Mu * g[3 * a + jj] * g[3 * b + i] + (i == jj ? Mu * gg : 0));
                    }
                }
            }
        });

        return vals;
    }

    private static void Apply(FlatPattern fp, double[] vals, double[] x, double[] y, bool[] fixedDof)
    {
        int[] rowStart = fp.RowStart, rowLen = fp.RowLen, col = fp.Col;
        Parallel.ForEach(Partitioner.Create(0, fp.NodeCount, 1024), r =>
        {
            for (var v = r.Item1; v < r.Item2; v++)
            {
                double y0 = 0, y1 = 0, y2 = 0;
                var rs = rowStart[v];
                if (rs >= 0)
                    for (int t = rs, end = rs + rowLen[v]; t < end; t++)
                    {
                        var at = 9 * t;
                        var c = 3 * col[t];
                        double x0 = x[c], x1 = x[c + 1], x2 = x[c + 2];
                        y0 += vals[at] * x0 + vals[at + 1] * x1 + vals[at + 2] * x2;
                        y1 += vals[at + 3] * x0 + vals[at + 4] * x1 + vals[at + 5] * x2;
                        y2 += vals[at + 6] * x0 + vals[at + 7] * x1 + vals[at + 8] * x2;
                    }

                y[3 * v] = fixedDof[3 * v] ? 0 : y0;
                y[3 * v + 1] = fixedDof[3 * v + 1] ? 0 : y1;
                y[3 * v + 2] = fixedDof[3 * v + 2] ? 0 : y2;
            }
        });
    }

    /// <summary>Block-Jacobi preconditioned CG on the free unknowns (prescribed values zero).</summary>
    private static (double[] U, int Iters, double Rel) Pcg(FlatPattern fp, double[] vals, double[] f, bool[] fixedDof,
        double tol, int maxIt, double[]? x0)
    {
        var nN = fp.NodeCount;
        var nd = 3 * nN;
        // inverse diagonal blocks, fixed unknowns decoupled
        var minv = new double[9 * nN];
        for (var v = 0; v < nN; v++)
        {
            var d = new double[9];
            var rs = fp.RowStart[v];
            if (rs >= 0)
                for (int t = rs, end = rs + fp.RowLen[v]; t < end; t++)
                    if (fp.Col[t] == v)
                        for (var q = 0; q < 9; q++)
                            d[q] = vals[9 * t + q];
            for (var i = 0; i < 3; i++)
            {
                if (!fixedDof[3 * v + i] && d[4 * i] != 0) continue;
                fixedDof[3 * v + i] = true; // also nodes without elements
                for (var q = 0; q < 3; q++) d[3 * i + q] = d[3 * q + i] = 0;
                d[4 * i] = 1;
            }

            Invert3(d, minv.AsSpan(9 * v, 9));
        }

        var u = new double[nd];
        var r = (double[])f.Clone();
        var r0 = 0.0;
        for (var i = 0; i < nd; i++)
        {
            if (fixedDof[i]) r[i] = 0;
            r0 += r[i] * r[i];
        }

        r0 = Math.Sqrt(r0);
        if (x0 != null)
        {
            // warm start: r = f - A x0 on the free unknowns
            for (var i = 0; i < nd; i++) u[i] = fixedDof[i] ? 0 : x0[i];
            var ax = new double[nd];
            Apply(fp, vals, u, ax, fixedDof);
            for (var i = 0; i < nd; i++) r[i] -= ax[i];
        }
        var z = new double[nd];
        var p = new double[nd];
        var q2 = new double[nd];
        void Precond(double[] src, double[] dst)
        {
            Parallel.For(0, nN, v =>
            {
                var m = minv.AsSpan(9 * v, 9);
                double a0 = src[3 * v], a1 = src[3 * v + 1], a2 = src[3 * v + 2];
                dst[3 * v] = fixedDof[3 * v] ? 0 : m[0] * a0 + m[1] * a1 + m[2] * a2;
                dst[3 * v + 1] = fixedDof[3 * v + 1] ? 0 : m[3] * a0 + m[4] * a1 + m[5] * a2;
                dst[3 * v + 2] = fixedDof[3 * v + 2] ? 0 : m[6] * a0 + m[7] * a1 + m[8] * a2;
            });
        }

        double DotV(double[] a, double[] b)
        {
            double acc = 0;
            for (var i = 0; i < nd; i++) acc += a[i] * b[i];
            return acc;
        }

        Precond(r, z);
        Array.Copy(z, p, nd);
        var rz = DotV(r, z);
        var it = 0;
        var rel = 1.0;
        for (; it < maxIt; it++)
        {
            Apply(fp, vals, p, q2, fixedDof);
            var alpha = rz / DotV(p, q2);
            for (var i = 0; i < nd; i++)
            {
                u[i] += alpha * p[i];
                r[i] -= alpha * q2[i];
            }

            rel = Math.Sqrt(DotV(r, r)) / r0;
            if (rel < tol) break;
            Precond(r, z);
            var rz2 = DotV(r, z);
            var beta = rz2 / rz;
            rz = rz2;
            for (var i = 0; i < nd; i++) p[i] = z[i] + beta * p[i];
        }

        return (u, it + 1, rel);
    }

    private static void Invert3(double[] a, Span<double> o)
    {
        var det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6]) +
                  a[2] * (a[3] * a[7] - a[4] * a[6]);
        o[0] = (a[4] * a[8] - a[5] * a[7]) / det;
        o[1] = (a[2] * a[7] - a[1] * a[8]) / det;
        o[2] = (a[1] * a[5] - a[2] * a[4]) / det;
        o[3] = (a[5] * a[6] - a[3] * a[8]) / det;
        o[4] = (a[0] * a[8] - a[2] * a[6]) / det;
        o[5] = (a[2] * a[3] - a[0] * a[5]) / det;
        o[6] = (a[3] * a[7] - a[4] * a[6]) / det;
        o[7] = (a[1] * a[6] - a[0] * a[7]) / det;
        o[8] = (a[0] * a[4] - a[1] * a[3]) / det;
    }

    /// <summary>Radius ratio 3 r_in / r_circ of the live tetrahedra (1 for the regular one).</summary>
    private static (double Min, double Mean) Quality(IncrementalMesh m, MeshState st)
    {
        double min = 1, sum = 0;
        var cnt = 0;
        var p = new double[12];
        for (var e = 0; e < m.ElementCount; e++)
        {
            var row = m.Row(e);
            if (row.Length != 4) continue;
            for (var i = 0; i < 4; i++)
            {
                p[3 * i] = st.X(row[i], 0);
                p[3 * i + 1] = st.X(row[i], 1);
                p[3 * i + 2] = st.X(row[i], 2);
            }

            var q = RadiusRatio(p);
            min = Math.Min(min, q);
            sum += q;
            cnt++;
        }

        return (min, sum / cnt);
    }

    private static double RadiusRatio(double[] p)
    {
        // r_in = 3V / sum of face areas; r_circ from the circumcentre solve
        double ax = p[3] - p[0], ay = p[4] - p[1], az = p[5] - p[2];
        double bx = p[6] - p[0], by = p[7] - p[1], bz = p[8] - p[2];
        double cx = p[9] - p[0], cy = p[10] - p[1], cz = p[11] - p[2];
        var det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
        var vol = Math.Abs(det) / 6;
        double Area(int i, int j, int k)
        {
            double ux = p[3 * j] - p[3 * i], uy = p[3 * j + 1] - p[3 * i + 1], uz = p[3 * j + 2] - p[3 * i + 2];
            double vx = p[3 * k] - p[3 * i], vy = p[3 * k + 1] - p[3 * i + 1], vz = p[3 * k + 2] - p[3 * i + 2];
            double qx = uy * vz - uz * vy, qy = uz * vx - ux * vz, qz = ux * vy - uy * vx;
            return 0.5 * Math.Sqrt(qx * qx + qy * qy + qz * qz);
        }

        var rin = 3 * vol / (Area(0, 1, 2) + Area(0, 1, 3) + Area(0, 2, 3) + Area(1, 2, 3));
        double a2 = ax * ax + ay * ay + az * az, b2 = bx * bx + by * by + bz * bz, c2 = cx * cx + cy * cy + cz * cz;
        // circumcentre o (relative to p0): 2 [a b c]^T o = [a2 b2 c2]
        var ox = (a2 * (by * cz - bz * cy) - ay * (b2 * cz - bz * c2) + az * (b2 * cy - by * c2)) / (2 * det);
        var oy = (ax * (b2 * cz - bz * c2) - a2 * (bx * cz - bz * cx) + az * (bx * c2 - b2 * cx)) / (2 * det);
        var oz = (ax * (by * c2 - b2 * cy) - ay * (bx * c2 - b2 * cx) + a2 * (bx * cy - by * cx)) / (2 * det);
        var rc = Math.Sqrt(ox * ox + oy * oy + oz * oz);
        return 3 * rin / rc;
    }

    // ------------------------------------------------------------------ interaction integral

    /// <summary>
    ///     Plane-strain (modes I, II) and antiplane (mode III) near-tip displacements, unit K, local
    ///     frame. <paramref name="side"/> is the side of the crack the element lies on: behind the
    ///     front the angle is continued across the cut on that side, since nodes joined to the
    ///     surface keep their coordinates and an element may reach slightly across the plane.
    /// </summary>
    private static void Aux(int mode, double y1, double y2, int side, Span<double> uOut)
    {
        var r = Math.Sqrt(y1 * y1 + y2 * y2);
        var phi = Math.Atan2(y2, y1);
        if (side > 0 && phi < -Math.PI / 2) phi += 2 * Math.PI;
        else if (side < 0 && phi > Math.PI / 2) phi -= 2 * Math.PI;
        var kappa = 3 - 4 * Nu;
        var sr = Math.Sqrt(r / (2 * Math.PI));
        double c = Math.Cos(phi / 2), s = Math.Sin(phi / 2);
        uOut.Clear();
        switch (mode)
        {
            case 0:
                uOut[0] = sr / (2 * Mu) * c * (kappa - 1 + 2 * s * s);
                uOut[1] = sr / (2 * Mu) * s * (kappa + 1 - 2 * c * c);
                break;
            case 1:
                uOut[0] = sr / (2 * Mu) * s * (kappa + 1 + 2 * c * c);
                uOut[1] = -sr / (2 * Mu) * c * (kappa - 1 - 2 * s * s);
                break;
            default:
                uOut[2] = 2 / Mu * sr * s;
                break;
        }
    }

    /// <summary>
    ///     Auxiliary gradient H[i, k] (global) of the curvilinear field in the form of Gosz and Moran (2002): the
    ///     plane-strain / antiplane gradient in the frame (N, n, T) of the foot point of x on the front curve (signed
    ///     normal distance r1, height r2 above the plane), with no derivative along the front, rotated to global axes.
    ///     The auxiliary stress C : sym(H) is the two-dimensional asymptotic stress, traction-free on the crack faces;
    ///     H is not a gradient (the curvature terms are missing), which the incompatibility term accounts for.
    /// </summary>
    private static void AuxH(int mode, double x, double y, double z, Plane pl, FrontCurve c, int side, Span<double> Hg)
    {
        var (p, q, h) = pl.Local(x, y, z);
        var (psi, r1) = c.Project(p, q);
        var (nx, ny, _) = c.Normal(psi);
        var r2 = h;
        var d = 1e-6 * Math.Max(Math.Sqrt(r1 * r1 + r2 * r2), 1e-12);
        Span<double> up = stackalloc double[3], um = stackalloc double[3];
        Span<double> Hl = stackalloc double[9];
        Hl.Clear();
        for (var j = 0; j < 2; j++)
        {
            Aux(mode, r1 + (j == 0 ? d : 0), r2 + (j == 1 ? d : 0), side, up);
            Aux(mode, r1 - (j == 0 ? d : 0), r2 - (j == 1 ? d : 0), side, um);
            for (var i = 0; i < 3; i++) Hl[3 * i + j] = (up[i] - um[i]) / (2 * d);
        }

        // basis B = [N, n, T] (columns), T = N x n; Hg = B Hl B^T
        Span<double> B = stackalloc double[9];
        for (var k = 0; k < 3; k++)
        {
            B[3 * k] = nx * pl.E1[k] + ny * pl.E2[k];
            B[3 * k + 1] = pl.N[k];
            B[3 * k + 2] = ny * pl.E1[k] - nx * pl.E2[k];
        }

        for (var i = 0; i < 3; i++)
        for (var k = 0; k < 3; k++)
        {
            double acc = 0;
            for (var aa = 0; aa < 3; aa++)
            for (var bb = 0; bb < 2; bb++)
                acc += B[3 * i + aa] * Hl[3 * aa + bb] * B[3 * k + bb];
            Hg[3 * i + k] = acc;
        }
    }

    /// <summary>
    ///     Derivatives of the auxiliary gradient, dH[9 j + 3 i + k] = d H_ik / d x_j, and the divergence of the
    ///     auxiliary stress, by central differences of <see cref="AuxH"/>.
    /// </summary>
    private static void AuxDerivatives(int mode, double x, double y, double z, Plane pl, FrontCurve c, int side,
        Span<double> dH, Span<double> div)
    {
        var (p, q, h) = pl.Local(x, y, z);
        var (_, r1) = c.Project(p, q);
        var d = 1e-3 * Math.Max(Math.Sqrt(r1 * r1 + h * h), 1e-9);
        Span<double> Hp = stackalloc double[9], Hm = stackalloc double[9];
        div.Clear();
        for (var j = 0; j < 3; j++)
        {
            AuxH(mode, x + (j == 0 ? d : 0), y + (j == 1 ? d : 0), z + (j == 2 ? d : 0), pl, c, side, Hp);
            AuxH(mode, x - (j == 0 ? d : 0), y - (j == 1 ? d : 0), z - (j == 2 ? d : 0), pl, c, side, Hm);
            for (var t = 0; t < 9; t++) dH[9 * j + t] = (Hp[t] - Hm[t]) / (2 * d);
            double trp = Hp[0] + Hp[4] + Hp[8], trm = Hm[0] + Hm[4] + Hm[8];
            for (var i = 0; i < 3; i++)
            {
                var sp = Mu * (Hp[3 * i + j] + Hp[3 * j + i]) + (i == j ? Lambda * trp : 0);
                var sm = Mu * (Hm[3 * i + j] + Hm[3 * j + i]) + (i == j ? Lambda * trm : 0);
                div[i] += (sp - sm) / (2 * d);
            }
        }
    }

    private static (double KI, double KII, double KIII, double J) Interaction(Snapshot s, double[] X, double[] u,
        Plane pl, FrontCurve curve, double th, double rq, double delta, double psiLo = double.NegativeInfinity,
        double psiHi = double.PositiveInfinity, bool freeSurface = false)
    {
        double[] nrm = pl.N;

        // virtual extension field V = q N at the nodes: N the in-plane outward normal of the front at the node's
        // foot point, q a tube of radius rq around the front times a hat of half-width delta in the curve
        // parameter (global form of the domain integral, Shih et al. 1986)
        var nN = s.NodeCount;
        var vn = new double[3 * nN];
        var any = new bool[nN];
        for (var v = 0; v < nN; v++)
        {
            var (p, q, h) = pl.Local(X[3 * v], X[3 * v + 1], X[3 * v + 2]);
            if (Math.Abs(h) >= rq) continue;
            var rho = Math.Sqrt(p * p + q * q);
            if (Math.Abs(rho - curve.R(Math.Atan2(q, p))) >= 2 * rq) continue;
            var (psi, dn) = curve.Project(p, q);
            var d = Math.Sqrt(dn * dn + h * h);
            var dpsi = Math.IEEERemainder(psi - th, 2 * Math.PI);
            var tube = QPlateau <= 0 ? Math.Max(0, 1 - d / rq)
                : d <= QPlateau * rq ? 1.0 : Math.Max(0, 1 - (d - QPlateau * rq) / (rq - QPlateau * rq));
            var qv = tube * Math.Max(0, 1 - Math.Abs(dpsi) / delta);
            if (qv <= 0) continue;
            any[v] = true;
            var (nx, ny, _) = curve.Normal(psi);
            for (var k = 0; k < 3; k++) vn[3 * v + k] = qv * (nx * pl.E1[k] + ny * pl.E2[k]);
        }

        // 4-point rule (degree 2); fixed blocks of elements in parallel, their sums of (I_I, I_II, I_III, J) added
        // in block order, so the result does not depend on the thread schedule
        const double qa = 0.5854101966249685, qb = 0.1381966011250105;
        var block = KBlock;
        var partial = new double[(s.ElementCount + block - 1) / block][];
        Parallel.For(0, partial.Length, b =>
        {
            var part = new double[4];
            for (int e = b * block, hi = Math.Min(s.ElementCount, e + block); e < hi; e++) Term(e, part);
            partial[b] = part;
        });
        var sum = new double[4];
        foreach (var part in partial)
            for (var k = 0; k < 4; k++)
                sum[k] += part[k];

        void Term(int e, double[] part)
        {
            Span<double> g = stackalloc double[12];
            Span<double> G = stackalloc double[9];   // FE displacement gradient G[i, k] = du_i / dx_k
            Span<double> S = stackalloc double[9];   // FE stress
            Span<double> GV = stackalloc double[9];  // GV[k, j] = dV_k / dx_j
            Span<double> H = stackalloc double[9];   // auxiliary gradient
            Span<double> dv = stackalloc double[3];
            Span<double> dH = stackalloc double[27];
            Span<double> Vq = stackalloc double[3];
            var row = s.Row(e);
            if (row.Length != 4) return;
            var touched = false;
            foreach (var v in row) touched |= any[v];
            if (!touched) return;
            var vol = Gradients(X, row, g);
            double cxs = 0, cys = 0, czs = 0;
            foreach (var v in row)
            {
                cxs += X[3 * v];
                cys += X[3 * v + 1];
                czs += X[3 * v + 2];
            }

            var side = pl.Local(cxs / 4, cys / 4, czs / 4).H >= 0 ? 1 : -1;
            G.Clear();
            GV.Clear();
            for (var a = 0; a < 4; a++)
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
            {
                G[3 * i + j] += u[3 * row[a] + i] * g[3 * a + j];
                GV[3 * i + j] += vn[3 * row[a] + i] * g[3 * a + j];
            }

            var tr = G[0] + G[4] + G[8];
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
                S[3 * i + j] = Mu * (G[3 * i + j] + G[3 * j + i]) + (i == j ? Lambda * tr : 0);
            double W = 0, divV = GV[0] + GV[4] + GV[8];
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
                W += 0.5 * S[3 * i + j] * 0.5 * (G[3 * i + j] + G[3 * j + i]);
            // J: (sigma_ij u_i,k - W delta_jk) V_k,j
            double jv = -W * divV;
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
            for (var k = 0; k < 3; k++)
                jv += S[3 * i + j] * G[3 * i + k] * GV[3 * k + j];
            part[3] += jv * vol;

            for (var qp = 0; qp < 4; qp++)
            {
                double xx = 0, yy = 0, zz = 0;
                Vq.Clear();
                for (var a = 0; a < 4; a++)
                {
                    var wa = a == qp ? qa : qb;
                    xx += wa * X[3 * row[a]];
                    yy += wa * X[3 * row[a] + 1];
                    zz += wa * X[3 * row[a] + 2];
                    for (var k = 0; k < 3; k++) Vq[k] += wa * vn[3 * row[a] + k];
                }

                for (var mode = 0; mode < 3; mode++)
                {
                    AuxH(mode, xx, yy, zz, pl, curve, side, H);
                    var htr = H[0] + H[4] + H[8];
                    double acc = 0, wint = 0;
                    for (var i = 0; i < 3; i++)
                    for (var j = 0; j < 3; j++)
                    {
                        var sa = Mu * (H[3 * i + j] + H[3 * j + i]) + (i == j ? Lambda * htr : 0);
                        wint += S[3 * i + j] * 0.5 * (H[3 * i + j] + H[3 * j + i]);
                        for (var k = 0; k < 3; k++)
                            acc += (S[3 * i + j] * H[3 * i + k] + sa * G[3 * i + k]) * GV[3 * k + j];
                    }

                    acc -= wint * divV;
                    // div(F V) - F : grad V for an equilibrated FE field (constant in the element):
                    // non-equilibrium div(sigma^aux)_i u_i,k V_k and incompatibility sigma_ij (H_ik,j - H_ij,k) V_k
                    AuxDerivatives(mode, xx, yy, zz, pl, curve, side, dH, dv);
                    for (var k = 0; k < 3; k++)
                    {
                        double ck = 0;
                        for (var i = 0; i < 3; i++)
                        {
                            ck += dv[i] * G[3 * i + k];
                            for (var j = 0; j < 3; j++)
                                ck += S[3 * i + j] * (dH[9 * j + 3 * i + k] - dH[9 * k + 3 * i + j]);
                        }

                        acc += ck * Vq[k];
                    }
                    part[mode] += acc * vol / 4;
                }
            }

            return;
        }

        // Free surface inside the tube: the auxiliary field is not traction-free there, and the divergence theorem
        // leaves the boundary term  B = int_{Sf} (sigma^aux_ij m_j u_i,k - W^int m_k) V_k dS  (sigma_ij m_j = 0 for the
        // solution); the interaction integral is D - B. Faces of the crack (normal along the crack normal) and the
        // loaded faces are not free surfaces of the front.
        if (freeSurface) FreeSurfaceTerm(sum);

        void FreeSurfaceTerm(double[] acc)
        {
            Span<double> g = stackalloc double[12];
            Span<double> G = stackalloc double[9];
            Span<double> S = stackalloc double[9];
            Span<double> H = stackalloc double[9];
            Span<double> V = stackalloc double[3];
            Span<int> t = stackalloc int[3];
            var bsum = new double[3];
            for (var i = 0; i < s.BoundaryCount; i++)
            {
                var id = s.BoundaryFacet(i);
                int e = id / 4, fl = id % 4;
                var row = s.Row(e);
                if (row.Length != 4) continue;
                var k0 = 0;
                for (var p = 0; p < 4; p++)
                    if (p != fl)
                        t[k0++] = row[p];
                if (!(any[t[0]] || any[t[1]] || any[t[2]])) continue;
                // outward normal and area from the facet nodes, oriented away from the opposite node
                double ux = X[3 * t[1]] - X[3 * t[0]], uy = X[3 * t[1] + 1] - X[3 * t[0] + 1], uz = X[3 * t[1] + 2] - X[3 * t[0] + 2];
                double wx = X[3 * t[2]] - X[3 * t[0]], wy = X[3 * t[2] + 1] - X[3 * t[0] + 1], wz = X[3 * t[2] + 2] - X[3 * t[0] + 2];
                double mx = uy * wz - uz * wy, my = uz * wx - ux * wz, mz = ux * wy - uy * wx;
                var area = 0.5 * Math.Sqrt(mx * mx + my * my + mz * mz);
                mx /= 2 * area;
                my /= 2 * area;
                mz /= 2 * area;
                var opp = row[fl];
                if (mx * (X[3 * opp] - X[3 * t[0]]) + my * (X[3 * opp + 1] - X[3 * t[0] + 1]) + mz * (X[3 * opp + 2] - X[3 * t[0] + 2]) > 0)
                {
                    mx = -mx;
                    my = -my;
                    mz = -mz;
                }

                var cn = Math.Abs(mx * pl.N[0] + my * pl.N[1] + mz * pl.N[2]);
                if (cn > 0.5) continue; // a crack face or a loaded face
                Gradients(X, row, g);
                G.Clear();
                for (var a = 0; a < 4; a++)
                for (var ii = 0; ii < 3; ii++)
                for (var j = 0; j < 3; j++)
                    G[3 * ii + j] += u[3 * row[a] + ii] * g[3 * a + j];
                var tr = G[0] + G[4] + G[8];
                for (var ii = 0; ii < 3; ii++)
                for (var j = 0; j < 3; j++)
                    S[3 * ii + j] = Mu * (G[3 * ii + j] + G[3 * j + ii]) + (ii == j ? Lambda * tr : 0);
                double cxs = 0, cys = 0, czs = 0;
                foreach (var v in row)
                {
                    cxs += X[3 * v];
                    cys += X[3 * v + 1];
                    czs += X[3 * v + 2];
                }

                var side = pl.Local(cxs / 4, cys / 4, czs / 4).H >= 0 ? 1 : -1;
                for (var qp = 0; qp < 3; qp++)
                {
                    double xx = 0, yy = 0, zz = 0;
                    V.Clear();
                    for (var a = 0; a < 3; a++)
                    {
                        var wa = a == qp ? 2.0 / 3 : 1.0 / 6;
                        xx += wa * X[3 * t[a]];
                        yy += wa * X[3 * t[a] + 1];
                        zz += wa * X[3 * t[a] + 2];
                        for (var k = 0; k < 3; k++) V[k] += wa * vn[3 * t[a] + k];
                    }

                    for (var mode = 0; mode < 3; mode++)
                    {
                        AuxH(mode, xx, yy, zz, pl, curve, side, H);
                        var htr = H[0] + H[4] + H[8];
                        double wint = 0, mv = mx * V[0] + my * V[1] + mz * V[2], b = 0;
                        for (var ii = 0; ii < 3; ii++)
                        {
                            double sam = 0;
                            for (var j = 0; j < 3; j++)
                            {
                                var sa = Mu * (H[3 * ii + j] + H[3 * j + ii]) + (ii == j ? Lambda * htr : 0);
                                sam += sa * (j == 0 ? mx : j == 1 ? my : mz);
                                wint += S[3 * ii + j] * 0.5 * (H[3 * ii + j] + H[3 * j + ii]);
                            }

                            for (var k = 0; k < 3; k++) b += sam * G[3 * ii + k] * V[k];
                        }

                        bsum[mode] += (b - wint * mv) * area / 3;
                    }
                }
            }

            for (var mode = 0; mode < 3; mode++) acc[mode] -= bsum[mode];
        }

        // J(s) times the integral of q along the front (V . N = q on the front) = volume integral
        var L = curve.HatLength(th, delta, 400, psiLo, psiHi);
        var ep = E / (1 - Nu * Nu);
        return (ep / 2 * sum[0] / L, ep / 2 * sum[1] / L, Mu * sum[2] / L, sum[3] / L);
    }
}

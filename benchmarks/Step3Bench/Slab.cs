// Slab.cs — plane-strain slab examples: centre, edge, collinear, inclined and T-junction cracks.
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
// A plate [-W/2, W/2] x [-H/2, H/2] of one Kuhn cell in thickness (t = h0) with u_z = 0 on
// both faces reproduces plane strain, so the two-dimensional solutions apply: Tada's centre
// and edge cracks, Erdogan's two collinear cracks, the inclined crack of Erdogan and Sih. The
// fronts are straight lines through the thickness, evaluated in through mode (one factor per
// front). Tension sigma on y = +-H/2; in-plane rigid motions removed by point constraints on
// the free faces x = +-W/2 at mid-height. The mesh is refined to h_min in the bands given by
// the example and graded outside with slope 0.35.
//
//   slab kind=centre|edge|collinear|inclined|tjunction [a= b= beta= W= H= n= div= rq= plateau=]
//        [steps= da= damin= m= tol= crit=mts|schoellmann|richard] [out=dir] [ens=dir] [check] [nofs]

using System.Diagnostics;
using System.Globalization;
using Step2;

namespace Step3.Bench;

public static class Slab
{
    private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

    public sealed class LineCrack
    {
        public V3 P0, P1; // endpoints in the mid-plane (z free)
        public int Bit;
        public string Name = "";
        public V3[]? Polyline; // kinked crack: the vertices from P0 to P1 (P0 and P1 included)
        public V3 Dir => (P1 - P0).Unit();
        public V3 Normal => V3.Cross(Dir, new V3(0, 0, 1)).Unit();
        public double Length => (P1 - P0).Norm;
    }

    public static void Run(Args a)
    {
        var kind = a.S("kind", "centre");
        double W = a.D("W", 8), H = a.D("H", 8);
        var n = a.I("n", 80);
        var div = a.I("div", 8);
        var h0 = W / n;
        var hmin = h0 / div;
        var ny = (int)Math.Round(H / h0);
        var sigma = a.D("sigma", 1.0);
        var mat = new Material(a.D("E", 1.0), a.D("nu", 0.3));
        var steps = a.I("steps", 0);
        var outDir = a.S("out", "");
        if (outDir.Length > 0) Directory.CreateDirectory(outDir);
        var yOff = 0.37 * hmin; // crack lines off the lattice
        var cracks = new List<LineCrack>();
        var bands = new List<(V3 Lo, V3 Hi)>(); // boxes at h_min
        var aLen = a.D("a", 0.2);
        var bLen = a.D("b", 0.4);
        var beta = a.D("beta", 45) * Math.PI / 180;
        var grow = steps > 0;
        var pad = grow ? a.D("pad", 1.0) * aLen : 0; // room for the path
        switch (kind)
        {
            case "centre":
                cracks.Add(new LineCrack { P0 = new V3(-aLen, yOff, 0), P1 = new V3(aLen, yOff, 0), Bit = 1, Name = "centre" });
                bands.Add((new V3(-aLen - pad, yOff - 0.4 * aLen - pad, -1), new V3(aLen + pad, yOff + 0.4 * aLen + pad, 1)));
                break;
            case "edge":
                // edge crack of length a from the face x = -W/2
                cracks.Add(new LineCrack { P0 = new V3(-W / 2 - 2 * h0, yOff, 0), P1 = new V3(-W / 2 + aLen, yOff, 0), Bit = 1, Name = "edge" });
                bands.Add((new V3(-W / 2, yOff - 0.4 * aLen - pad, -1), new V3(-W / 2 + aLen + 0.4 * aLen + pad, yOff + 0.4 * aLen + pad, 1)));
                break;
            case "collinear":
                // cracks on a <= |x| <= b
                cracks.Add(new LineCrack { P0 = new V3(-bLen, yOff, 0), P1 = new V3(-aLen, yOff, 0), Bit = 1, Name = "left" });
                cracks.Add(new LineCrack { P0 = new V3(aLen, yOff, 0), P1 = new V3(bLen, yOff, 0), Bit = 2, Name = "right" });
                var lig = bLen - aLen;
                bands.Add((new V3(-bLen - pad, yOff - 0.5 * lig - 0.5 * pad, -1), new V3(bLen + pad, yOff + 0.5 * lig + 0.5 * pad, 1)));
                break;
            case "inclined":
            {
                // crack of half-length a through the origin, angle beta between the crack and the load axis (y)
                var e = new V3(Math.Sin(beta), Math.Cos(beta), 0);
                var c = new V3(0.37 * hmin, yOff, 0);
                cracks.Add(new LineCrack { P0 = c - aLen * e, P1 = c + aLen * e, Bit = 1, Name = "inclined" });
                var r = aLen + pad;
                bands.Add((new V3(-r, -0.75 * r, -1), new V3(r, 0.75 * r, 1)));
                break;
            }
            case "kinked":
            {
                // main crack of half-length a at angle beta, then a kink of length ell = kl * a at the angle kink (deg, from the
                // crack direction toward its normal) at the tip P1; static check against Cotterell and Rice (small kink)
                var e = new V3(Math.Sin(beta), Math.Cos(beta), 0);
                var nrm0 = V3.Cross(e, new V3(0, 0, 1)).Unit();
                var c = new V3(0.37 * hmin, yOff, 0);
                var kink = a.D("kink", -53.13) * Math.PI / 180;
                var ell = a.D("kl", 0.2) * aLen;
                var d = (Math.Cos(kink) * e + Math.Sin(kink) * nrm0).Unit();
                V3 p0 = c - aLen * e, p1 = c + aLen * e, p2 = p1 + ell * d;
                cracks.Add(new LineCrack { P0 = p0, P1 = p2, Bit = 1, Name = "kinked", Polyline = [p0, p1, p2] });
                var r = aLen + ell + pad;
                bands.Add((new V3(-r, -0.75 * r, -1), new V3(r, 0.75 * r, 1)));
                break;
            }
            case "tjunction":
            {
                // horizontal crack of half-length b at height y1 (bit 1); inclined crack of half-length a below it (bit 2) growing up
                var y1 = a.D("y1", 0.5);
                cracks.Add(new LineCrack { P0 = new V3(-bLen, y1 + yOff, 0), P1 = new V3(bLen, y1 + yOff, 0), Bit = 1, Name = "barrier" });
                var e = new V3(Math.Sin(beta), Math.Cos(beta), 0);
                var c = new V3(0.37 * hmin, yOff, 0);
                cracks.Add(new LineCrack { P0 = c - aLen * e, P1 = c + aLen * e, Bit = 2, Name = "inclined" });
                bands.Add((new V3(-bLen - pad, -aLen - pad, -1), new V3(bLen + pad, y1 + 0.3 * aLen + pad, 1)));
                break;
            }
            default:
                throw new ArgumentException($"unknown kind {kind}");
        }

        var sw = Stopwatch.StartNew();
        var cm = KuhnMesh.Build(n, ny, 1, new V3(-W / 2, -H / 2, 0), new V3(W / 2, H / 2, h0));
        var refiner = new MaubachRefiner(cm);
        double Size(V3 x)
        {
            var d = double.MaxValue;
            foreach (var (lo, hi) in bands)
            {
                var dx = Math.Max(0, Math.Max(lo.X - x.X, x.X - hi.X));
                var dy = Math.Max(0, Math.Max(lo.Y - x.Y, x.Y - hi.Y));
                d = Math.Min(d, Math.Sqrt(dx * dx + dy * dy));
            }

            return hmin + 0.35 * d;
        }

        // size: h_min within r_q of the crack lines (and in the growth bands), graded outside with slope 0.35
        var rq = a.D("rq", aLen / 3);
        double SizeRq(V3 x)
        {
            var d = double.MaxValue;
            foreach (var cr in cracks)
            {
                Geo.SegmentSegment(cr.P0, cr.P1, x, x, out _, out _);
                var dd = Geo.SegmentSegment(new V3(cr.P0.X, cr.P0.Y, x.Z), new V3(cr.P1.X, cr.P1.Y, x.Z), x, x, out _, out _);
                d = Math.Min(d, Math.Max(0, dd - rq));
            }

            if (grow) d = Math.Min(d, Size(x) - hmin);
            return hmin + 0.35 * d;
        }

        refiner.RefineToSize(SizeRq);
        var snap0 = cm.Commit();
        Console.WriteLine($"slab {kind}: W {W} H {H} t {h0:G4} h0 {h0:G4} hmin {hmin:G5}; mesh {snap0.ElementCount} tets {snap0.NodeCount} nodes after {refiner.Bisections} bisections, {sw.ElapsedMilliseconds} ms");

        var opt = new GrowthOptions
        {
            Through = !a.Flag("nothrough"), DaMax = a.D("da", 2 * hmin), DaMin = a.D("damin", hmin), TolRelative = a.D("tol", 0.06),
            Rq = rq, RqElements = a.D("rqel", 0), KnotSpacingElements = 2, MinKnots = 6, Fit = !a.Flag("nofit"),
            Paris = new ParisLaw { M = a.D("m", 3), C = a.D("C", 0) }, Dct = !a.Flag("nodct"), SolverTol = a.D("stol", 1e-10)
        };
        opt.Sif.Plateau = a.D("plateau", 0.5);
        opt.Sif.SurfaceFollowing = !a.Flag("planar");
        opt.Sif.FreeSurfaceTerm = !a.Flag("nofs");
        opt.Sif.CrackFaceTerm = !a.Flag("nocf");
        opt.Advance.Direction = a.S("crit", "mts") switch
        {
            "schoellmann" => DirectionCriterion.Schoellmann, "richard" => DirectionCriterion.Richard, "planar" => DirectionCriterion.Planar,
            _ => DirectionCriterion.MaximumTangentialStress
        };
        opt.Advance.KinkMax = a.D("kinkmax", 70) * Math.PI / 180;
        opt.Advance.CurvatureElements = a.D("curv", 6);
        opt.Advance.SnapElements = a.D("snap", 1.0);
        opt.Advance.OvershootElements = a.D("over", 0.6);

        BoundaryConditions Bc(Snapshot s, double[] X)
        {
            var f = new double[3 * s.NodeCount];
            var fixedDof = new bool[3 * s.NodeCount];
            var kinds = new Dictionary<int, FacetKind>();
            const double tol = 1e-9;
            BoundaryConditions.Traction(s, X, (c, nrm) => Math.Abs(c.Y - H / 2) < tol, new V3(0, sigma, 0), f, kinds);
            BoundaryConditions.Traction(s, X, (c, nrm) => Math.Abs(c.Y + H / 2) < tol, new V3(0, -sigma, 0), f, kinds);
            BoundaryConditions.Mark(s, X, (c, nrm) => Math.Abs(c.Z) < tol || Math.Abs(c.Z - h0) < tol, FacetKind.Symmetry, kinds);
            int best1 = -1, best2 = -1;
            double d1 = double.MaxValue, d2 = double.MaxValue;
            for (var v = 0; v < s.NodeCount; v++)
            {
                if (s.TransposePairs(v).Length == 0) continue;
                double z = X[3 * v + 2];
                if (Math.Abs(z) < tol || Math.Abs(z - h0) < tol) fixedDof[3 * v + 2] = true;
                var q1 = (X[3 * v] + W / 2) * (X[3 * v] + W / 2) + X[3 * v + 1] * X[3 * v + 1] + z * z;
                var q2 = (X[3 * v] - W / 2) * (X[3 * v] - W / 2) + X[3 * v + 1] * X[3 * v + 1] + z * z;
                if (q1 < d1) { d1 = q1; best1 = v; }
                if (q2 < d2) { d2 = q2; best2 = v; }
            }

            fixedDof[3 * best1] = fixedDof[3 * best1 + 1] = true;
            fixedDof[3 * best2 + 1] = true;
            return new BoundaryConditions { Force = f, FixedDof = fixedDof, FacetKinds = kinds };
        }

        var drv = new GrowthDriver(cm, mat, Bc, opt);
        foreach (var cr in cracks)
        {
            var e = cr.Dir;
            var nrm = cr.Normal;
            var mid = 0.5 * (cr.P0 + cr.P1);
            var half = 0.5 * cr.Length;
            var margin = hmin / 4;
            HashSet<int>? onWall = null;
            // polyline cracks: the extruded polyline as an oriented soup with the arc length as vertex value
            TriangleSoup? poly = null;
            var polyLength = 0.0;
            if (cr.Polyline != null)
            {
                poly = new TriangleSoup();
                var sAcc = 0.0;
                for (var k = 0; k + 1 < cr.Polyline.Length; k++)
                {
                    V3 q0 = cr.Polyline[k], q1 = cr.Polyline[k + 1];
                    var seg = (q1 - q0).Norm;
                    var up = new V3(0, 0, 1);
                    var nk = V3.Cross((q1 - q0).Unit(), up);
                    V3 a0 = q0 - h0 * up, a1 = q0 + 2 * h0 * up, b0 = q1 - h0 * up, b1 = q1 + 2 * h0 * up;
                    poly.AddOriented(a0, b0, b1, nk, k, sAcc, sAcc + seg, sAcc + seg);
                    poly.AddOriented(a0, b1, a1, nk, k, sAcc, sAcc + seg, sAcc);
                    sAcc += seg;
                }

                poly.Build();
                polyLength = sAcc;
            }

            if (opt.Fit)
            {
                onWall = new HashSet<int>();
                // walls at the tips inside the plate: planes normal to the crack line through the thickness
                var tips = cr.Polyline != null
                    ? new[] { (cr.Polyline[0], (cr.Polyline[1] - cr.Polyline[0]).Unit()), (cr.Polyline[^1], (cr.Polyline[^1] - cr.Polyline[^2]).Unit()) }
                    : new[] { (cr.P0, e), (cr.P1, e) };
                foreach (var (tip, dirTip) in tips)
                {
                    if (Math.Abs(tip.X) > W / 2 || Math.Abs(tip.Y) > H / 2) continue;
                    var wallHalf = 2 * hmin;
                    var nTip = V3.Cross(dirTip, new V3(0, 0, 1)).Unit();
                    V3 p0 = new(tip.X, tip.Y, -h0), p1 = new(tip.X, tip.Y, 2 * h0);
                    var wall = FrontFitter.Wall([p0, p1], [nTip, nTip], wallHalf, false);
                    var fit = FrontFitter.Fit(cm, wall, hmin);
                    onWall.UnionWith(fit.WallNodes);
                    Console.WriteLine($"fit tip ({tip.X:F3},{tip.Y:F3}): {fit.Moved}/{fit.Candidates} moved, {fit.Rejected} rejected, {fit.EdgesCut} edges cut");
                }

                margin = opt.FitMarginElements * hmin;
            }

            Field3 phiS, phiR;
            V3 refNormal = nrm, refPoint = mid;
            if (poly != null)
            {
                var pl = poly;
                var len = polyLength;
                phiS = x => pl.SignedDistance(x).Phi;
                phiR = x =>
                {
                    var (tri, cp, _, _) = pl.Nearest(x);
                    var sv = pl.ValueAt(tri, cp);
                    return Math.Max(sv - len, -sv) - margin;
                };
                refNormal = poly.Normal(0);
                refPoint = 0.5 * (cr.Polyline![0] + cr.Polyline[1]);
            }
            else
            {
                phiS = x => V3.Dot(x - cr.P0, nrm);
                phiR = x => Math.Abs(V3.Dot(x - mid, e)) - half - margin;
            }

            drv.InsertInitial(cr.Bit, phiS, phiR, refPoint, refNormal, onWall);
        }

        if (a.Flag("check")) Check(drv);
        drv.Solve();
        drv.Current = drv.EvaluateFronts();
        Report(drv, 0, kind, a, sigma, W, cracks, outDir);
        var ensDir = a.S("ens", "");
        var ensSteps = new List<int>();
        var ensTimes = new List<double>();
        void Ens(int step)
        {
            if (ensDir.Length == 0) return;
            EnsightGold.WriteStep(ensDir, "slab", step, drv.Snap, cm, drv.Count, drv.X, drv.U, mat, true, (new V3(-W / 2, -H / 2, 0), new V3(W / 2, H / 2, h0)), $"slab {kind}");
            ensSteps.Add(step);
            ensTimes.Add(step);
            EnsightGold.WriteCase(ensDir, "slab", true, ensSteps, ensTimes);
        }

        Ens(0);
        for (var st = 1; st <= steps; st++)
        {
            var rec = drv.Step();
            if (a.Flag("check")) Check(drv);
            Report(drv, st, kind, a, sigma, W, cracks, outDir, rec);
            Ens(st);
            if (drv.Count.Chains.Count == 0)
            {
                Console.WriteLine("no front left");
                break;
            }
        }

        Console.WriteLine($"total {sw.Elapsed.TotalSeconds:F1} s");
    }

    private static void Check(GrowthDriver drv)
    {
        var r = CrackChecker.Verify(drv.Snap, drv.Mesh, drv.Count);
        Console.WriteLine($"check: {r.Pairs} crack facet pairs, {r.Front} front, {r.Mouth} mouth, {r.Junction} junction edges agree with the count");
    }

    /// <summary>Prints the factors of every front with the reference of the example, and appends them to out/steps.txt.</summary>
    private static void Report(GrowthDriver drv, int step, string kind, Args a, double sigma, double W, List<LineCrack> cracks, string outDir, StepRecord? rec = null)
    {
        var sb = new System.Text.StringBuilder();
        var cc = drv.Count;
        sb.Append($"step {step}: {cc.SheetCount} sheets, {cc.Chains.Count} fronts");
        if (rec != null) sb.Append($", da {rec.DaMax:G4}, dPi {rec.EnergyAfter - rec.EnergyBefore:G5}, -int G dA {-rec.SweptWork:G5}, events [{string.Join("; ", rec.Events)}]");
        Console.WriteLine(sb.ToString());
        // tips: order the fronts by x of their mid knot
        var fronts = drv.Current.OrderBy(f => f.Curve.P[f.M / 2].X).ToList();
        var refs = ReferenceFor(kind, a, sigma, W, cracks, fronts);
        for (var i = 0; i < fronts.Count; i++)
        {
            var f = fronts[i];
            var p = f.Curve.P[f.M / 2];
            var k = f.K[f.M / 2];
            var dct = f.DctKI == null ? "" : $" dct {f.DctKI.Where(x => !double.IsNaN(x)).DefaultIfEmpty(double.NaN).Average():G5} {f.DctKII!.Where(x => !double.IsNaN(x)).DefaultIfEmpty(double.NaN).Average():G5}";
            var g = f.G[f.M / 2];
            var line = $"  front {i} sheet {f.Chain.Sheet} at ({p.X:F4},{p.Y:F4}): KI {k.KI:G6} KII {k.KII:G6} KIII {k.KIII:G6} J/G {(g > 0 ? k.J / g : double.NaN):F4} kink {f.Kink[f.M / 2] * 180 / Math.PI:F2} deg Kv {f.Kv[f.M / 2]:G5}{dct}";
            if (refs != null && i < refs.Count && refs[i].HasValue)
            {
                var (rI, rII, label) = refs[i]!.Value;
                line += $" | ref {label}: KI {rI:G6} ({100 * (k.KI - rI) / rI:+0.00;-0.00} %)";
                if (rII != 0) line += $" KII {rII:G6} ({100 * (k.KII - rII) / rII:+0.00;-0.00} %) kink_mts {Reference.MtsKink(rI, rII) * 180 / Math.PI:F2} deg";
            }

            Console.WriteLine(line);
        }

        if (outDir.Length == 0) return;
        using (var w = new StreamWriter(Path.Combine(outDir, "steps.txt"), true))
        {
            for (var i = 0; i < fronts.Count; i++)
            {
                var f = fronts[i];
                var p = f.Curve.P[f.M / 2];
                var k = f.K[f.M / 2];
                var rI = refs != null && i < refs.Count && refs[i].HasValue ? refs[i]!.Value.KI : double.NaN;
                var rII = refs != null && i < refs.Count && refs[i].HasValue ? refs[i]!.Value.KII : double.NaN;
                w.WriteLine(string.Join(" ", new[]
                {
                    step.ToString(Ci), i.ToString(Ci), f.Chain.Sheet.ToString(Ci), p.X.ToString("G8", Ci), p.Y.ToString("G8", Ci), k.KI.ToString("G8", Ci), k.KII.ToString("G8", Ci),
                    k.KIII.ToString("G8", Ci), k.J.ToString("G8", Ci), f.G[f.M / 2].ToString("G8", Ci), (f.Kink[f.M / 2] * 180 / Math.PI).ToString("G6", Ci),
                    rI.ToString("G8", Ci), rII.ToString("G8", Ci), (f.DctKI?[f.M / 2] ?? double.NaN).ToString("G8", Ci), (f.DctKII?[f.M / 2] ?? double.NaN).ToString("G8", Ci),
                    (rec?.DaMax ?? 0).ToString("G6", Ci), drv.Energy.ToString("G10", Ci), (rec?.SweptWork ?? 0).ToString("G8", Ci), cc.SheetCount.ToString(Ci), cc.Chains.Count.ToString(Ci),
                    drv.Snap.ElementCount.ToString(Ci), (rec?.Iterations ?? 0).ToString(Ci), (rec?.MsSolve ?? 0).ToString("F0", Ci)
                }));
            }
        }

        drv.WriteFronts(Path.Combine(outDir, $"fronts_{step:D3}.txt"));
        drv.WriteFactors(Path.Combine(outDir, $"factors_{step:D3}.txt"));
        if (rec != null)
            File.AppendAllText(Path.Combine(outDir, "events.txt"), $"{step} {string.Join("; ", rec.Events)}\n");
    }

    /// <summary>Reference factors (KI, KII, label) per front ordered by x, or null entries when no reference applies.</summary>
    private static List<(double KI, double KII, string Label)?>? ReferenceFor(string kind, Args a, double sigma, double W, List<LineCrack> cracks, List<FrontFactors> fronts)
    {
        var list = new List<(double, double, string)?>();
        switch (kind)
        {
            case "centre":
            {
                // half-length from the tips (mid knots of the two fronts)
                if (fronts.Count != 2) return null;
                var half = 0.5 * (fronts[1].Curve.P[fronts[1].M / 2].X - fronts[0].Curve.P[fronts[0].M / 2].X);
                var k = sigma * Math.Sqrt(Math.PI * half) * Reference.CentreCrackStrip(half, W);
                list.Add((k, 0, $"Tada centre 2a/W={2 * half / W:F3}"));
                list.Add((k, 0, $"Tada centre 2a/W={2 * half / W:F3}"));
                return list;
            }
            case "edge":
            {
                if (fronts.Count != 1) return null;
                var len = fronts[0].Curve.P[fronts[0].M / 2].X + W / 2;
                list.Add((sigma * Math.Sqrt(Math.PI * len) * Reference.EdgeCrackStrip(len, W), 0, $"Tada edge a/W={len / W:F3}"));
                return list;
            }
            case "collinear":
            {
                if (fronts.Count == 4)
                {
                    double xo1 = fronts[0].Curve.P[fronts[0].M / 2].X, xi1 = fronts[1].Curve.P[fronts[1].M / 2].X,
                        xi2 = fronts[2].Curve.P[fronts[2].M / 2].X, xo2 = fronts[3].Curve.P[fronts[3].M / 2].X;
                    var ai = 0.5 * (xi2 - xi1);
                    var bo = 0.5 * (xo2 - xo1);
                    var (inner, outer) = Reference.CollinearCracks(ai, bo, sigma);
                    var lbl = $"Erdogan a/b={ai / bo:F3}";
                    list.Add((outer, 0, lbl + " outer"));
                    list.Add((inner, 0, lbl + " inner"));
                    list.Add((inner, 0, lbl + " inner"));
                    list.Add((outer, 0, lbl + " outer"));
                    return list;
                }

                if (fronts.Count == 2)
                {
                    var half = 0.5 * (fronts[1].Curve.P[fronts[1].M / 2].X - fronts[0].Curve.P[fronts[0].M / 2].X);
                    var k = sigma * Math.Sqrt(Math.PI * half) * Reference.CentreCrackStrip(half, W);
                    list.Add((k, 0, $"Tada centre (coalesced) 2a/W={2 * half / W:F3}"));
                    list.Add((k, 0, $"Tada centre (coalesced) 2a/W={2 * half / W:F3}"));
                    return list;
                }

                return null;
            }
            case "kinked":
            {
                // Cotterell and Rice (1980) first-order kinked-crack factors at the kinked tip (the front at larger x along the kink)
                if (fronts.Count != 2) return null;
                var beta = a.D("beta", 45) * Math.PI / 180;
                var kink = a.D("kink", -53.13) * Math.PI / 180;
                var aMain = a.D("a", 0.2);
                var (kI, kII) = Reference.InclinedCrack(aMain, sigma, beta);
                // sign of K_II at the kinked tip P1 in the frame (N = e, B = normal): the inclined formula gives the magnitude;
                // the kink angle is measured from N toward B, and the MTS kink has the sign opposite to K_II
                var kIIs = -Math.Sign(kink) * kII;
                double c1 = Math.Cos(kink / 2), c3 = Math.Cos(1.5 * kink), s1 = Math.Sin(kink / 2), s3 = Math.Sin(1.5 * kink);
                double c11 = (3 * c1 + c3) / 4, c12 = -3 * (s1 + s3) / 4, c21 = (s1 + s3) / 4, c22 = (c1 + 3 * c3) / 4;
                var kIk = c11 * kI + c12 * kIIs;
                var kIIk = c21 * kI + c22 * kIIs;
                // the straight tip P0 (first front by x for beta < 90) keeps the inclined values
                var e = new V3(Math.Sin(beta), Math.Cos(beta), 0);
                var kinkedFirst = fronts[0].Curve.P[fronts[0].M / 2].X > fronts[1].Curve.P[fronts[1].M / 2].X == (e.X > 0);
                var s0 = Math.Sign(fronts[kinkedFirst ? 1 : 0].K[fronts[0].M / 2].KII);
                var straight = (kI, (s0 == 0 ? 1 : s0) * kII, $"inclined tip beta={a.D("beta", 45)}");
                var kinked = (kIk, kIIk, $"Cotterell-Rice kink={a.D("kink", -53.13)} l/a={a.D("kl", 0.2)}");
                list.Add(kinkedFirst ? kinked : straight);
                list.Add(kinkedFirst ? straight : kinked);
                return list;
            }
            case "inclined":
            {
                if (fronts.Count != 2) return null;
                var beta = a.D("beta", 45) * Math.PI / 180;
                var half = 0.5 * (fronts[1].Curve.P[fronts[1].M / 2] - fronts[0].Curve.P[fronts[0].M / 2]).Norm;
                var (kI, kII) = Reference.InclinedCrack(half, sigma, beta);
                // sign of K_II in the frame (N, B) of each tip: N outward along the crack, B the crack normal (reference seeds (cos beta, -sin beta))
                // the computed sign is reported as is; the reference carries the magnitude with the sign of the first evaluation
                var s0 = Math.Sign(fronts[0].K[fronts[0].M / 2].KII);
                var s1 = Math.Sign(fronts[1].K[fronts[1].M / 2].KII);
                list.Add((kI, (s0 == 0 ? 1 : s0) * kII, $"inclined beta={a.D("beta", 45)}"));
                list.Add((kI, (s1 == 0 ? 1 : s1) * kII, $"inclined beta={a.D("beta", 45)}"));
                return list;
            }
            default:
                return null;
        }
    }
}

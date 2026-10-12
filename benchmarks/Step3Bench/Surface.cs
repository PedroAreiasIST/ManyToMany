// Surface.cs — semi-elliptical surface crack of Newman and Raju growing through the plate (breakthrough).
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
//   surface [a=0.2 c=0.2 t=0.5 L=1 n=20 div=4 rq= steps= da= m= out= ens= check nofit]
//
// Plate [-1, 1] x [-t, 0] x [-L, L] (half-width b = 1), tension 1 along z, crack of depth a and
// half-length c at the face y = -t in the plane z = 0.37 h_min (Step2, Section 8.4). Before
// breakthrough the factors are compared with the equation of Newman and Raju (1981) at the
// parametric angle; the free-surface term acts at the mouth on y = -t. Paris growth drives the
// deepest point to the back face y = 0: the increment snaps through it, the counting reports
// mouth edges on the back face and the single open front becomes two through-thickness fronts,
// which are then compared with Tada's centre crack of half-length c in the strip of width 2.

using Step2;

namespace Step3.Bench;

public static class Surface
{
    public static void Run(Args a)
    {
        var n = a.I("n", 20);
        var div = a.I("div", 4);
        var h0 = 2.0 / n;
        var hmin = h0 / div;
        double aD = a.D("a", 0.2), cH = a.D("c", 0.2), t = a.D("t", 0.5), L = a.D("L", 1);
        var rq = a.D("rq", aD / 3); // a/3 as for the penny: a tube of 4 h at a/8 reaches the vertex field and puts the surface point 13 % low
        var sigma = 1.0;
        var mat = new Material(1, 0.3);
        var steps = a.I("steps", 0);
        var outDir = a.S("out", "");
        if (outDir.Length > 0) Directory.CreateDirectory(outDir);
        var zOff = 0.37 * hmin;
        V3 lo = new(-1, -t, -L), hi = new(1, 0, L);
        var pad = steps > 0 ? a.D("pad", 1.0) * aD : 0;
        // distance to the front band: the semi-ellipse (and its growth ring) in the plane z = zOff
        double Dist(V3 x)
        {
            var rho = Math.Sqrt(x.X * x.X / (cH * cH) + (x.Y + t) * (x.Y + t) / (aD * aD)); // 1 on the front
            var scale = Math.Min(aD, cH);
            var dr = Math.Max(0, Math.Max(rho - 1 - pad / scale, 1 - pad / scale - rho)) * scale;
            return Math.Sqrt(dr * dr + (x.Z - zOff) * (x.Z - zOff));
        }

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var cm = Box3D.Mesh(lo, hi, h0, hmin, rq, Dist, out var refiner);
        Console.WriteLine($"surface: a {aD} c {cH} t {t} L {L} h0 {h0:G4} hmin {hmin:G5} rq {rq:G4}; mesh {cm.M.ElementCount} tets after {refiner.Bisections} bisections, {sw.ElapsedMilliseconds} ms");
        var opt = Box3D.Options(a, hmin, rq);
        var drv = new GrowthDriver(cm, mat, Box3D.TensionZ(lo, hi, sigma), opt);
        drv.Inserter.RejectSameSide = !a.Flag("samesideok"); // samesideok: a star whose components all lie on one side is not split instead of rejecting
        var nrm = new V3(0, 0, 1);
        var margin = hmin / 4;
        HashSet<int>? onWall = null;
        HashSet<int>? nearWall = null;
        if (opt.Fit)
        {
            const int m = 96;
            var p = new V3[m];
            var b = new V3[m];
            for (var k = 0; k < m; k++)
            {
                var ph = -0.15 + (Math.PI + 0.3) * k / (m - 1); // slightly beyond the face at both ends
                p[k] = new V3(cH * Math.Cos(ph), -t + aD * Math.Sin(ph), zOff);
                b[k] = nrm;
            }

            var wall = FrontFitter.Wall(p, b, 2 * hmin, false);
            var fit = FrontFitter.Fit(cm, wall, hmin);
            onWall = fit.WallNodes;
            nearWall = fit.Unmoved;
            margin = opt.FitMarginElements * hmin;
            Console.WriteLine($"fit: {fit.Moved}/{fit.Candidates} moved, {fit.Rejected} rejected, {fit.EdgesCut} edges cut");
        }

        Field3 phiS = x => x.Z - zOff;
        Field3 phiR = x =>
        {
            var rho = Math.Sqrt(x.X * x.X / (cH * cH) + (x.Y + t) * (x.Y + t) / (aD * aD));
            return (rho - 1) * Math.Min(aD, cH) - margin;
        };
        {
            // intended front: the semi-ellipse inside the plate, N the outward in-plane normal of the ellipse
            const int mc = 192;
            var pc = new V3[mc];
            var nc = new V3[mc];
            var bc = new V3[mc];
            for (var k = 0; k < mc; k++)
            {
                var ph = Math.PI * k / (mc - 1);
                pc[k] = new V3(cH * Math.Cos(ph), -t + aD * Math.Sin(ph), zOff);
                nc[k] = new V3(aD * Math.Cos(ph), cH * Math.Sin(ph), 0).Unit();
                bc[k] = nrm;
            }

            drv.InsertInitial(1, phiS, phiR, new V3(0, -t + aD / 2, zOff), nrm, onWall, hmin, nearWall, [new FrontCurve3(pc, nc, bc, false, EndKind.Mouth, EndKind.Mouth)]);
        }
        if (a.Flag("check")) Check(drv);
        drv.Solve();
        drv.Current = drv.EvaluateFronts();
        Report(drv, 0, aD, cH, t, sigma, outDir);
        var ensDir = a.S("ens", "");
        var ensSteps = new List<int>();
        void Ens(int step)
        {
            if (ensDir.Length == 0) return;
            EnsightGold.WriteStep(ensDir, "surface", step, drv.Snap, cm, drv.Count, drv.X, drv.U, mat, a.Flag("body"), (lo, hi), "surface crack");
            ensSteps.Add(step);
            EnsightGold.WriteCase(ensDir, "surface", a.Flag("body"), ensSteps, ensSteps.Select(s => (double)s).ToList());
        }

        Ens(0);
        for (var st = 1; st <= steps; st++)
        {
            var rec = drv.Step();
            if (a.Flag("check")) Check(drv);
            Report(drv, st, aD, cH, t, sigma, outDir, rec);
            Ens(st);
        }

        Console.WriteLine($"total {sw.Elapsed.TotalSeconds:F1} s");
    }

    private static void Check(GrowthDriver drv)
    {
        var r = CrackChecker.Verify(drv.Snap, drv.Mesh, drv.Count);
        Console.WriteLine($"check: {r.Pairs} crack facet pairs, {r.Front} front, {r.Mouth} mouth, {r.Junction} junction edges agree with the count");
    }

    private static void Report(GrowthDriver drv, int step, double a0, double c0, double t, double sigma, string outDir, StepRecord? rec = null)
    {
        var cc = drv.Count;
        var head = $"step {step}: {cc.SheetCount} sheets, {cc.Chains.Count} fronts, {cc.PairCount} crack facets";
        if (rec != null) head += $", da {rec.DaMax:G4}, dPi {rec.EnergyAfter - rec.EnergyBefore:G5}, -int G dA {-rec.SweptWork:G5}, events [{string.Join("; ", rec.Events)}]";
        Console.WriteLine(head);
        // current a (depth) and c (half-length) from the front nodes
        double aCur = 0, cCur = 0;
        foreach (var ch in cc.Chains)
        foreach (var p in ch.P)
        {
            aCur = Math.Max(aCur, p.Y + t);
            cCur = Math.Max(cCur, Math.Abs(p.X));
        }

        var through = cc.Chains.Count >= 2 && cc.Chains.All(ch => ch.Start == EndKind.Mouth && ch.End == EndKind.Mouth) && aCur > t - 1e-6;
        (double, double, double)? Ref(FrontFactors f, int j)
        {
            if (through)
            {
                var k = sigma * Math.Sqrt(Math.PI * cCur) * Reference.CentreCrackStrip(cCur, 2.0);
                return (k, 0, 0);
            }

            var p = f.Curve.P[j];
            var phi = Math.Atan2((p.Y + t) / aCur, p.X / cCur);
            if (aCur / cCur > 1.001) return null; // the equation is fitted for a/c <= 1
            return (Reference.NewmanRaju(aCur, cCur, t, 1.0, sigma, phi), 0, 0);
        }

        Console.WriteLine($"  a {aCur:G4} c {cCur:G4} a/c {aCur / cCur:F3} a/t {aCur / t:F3}{(through ? " THROUGH" : "")}");
        for (var c = 0; c < drv.Current.Count; c++)
        {
            var f = drv.Current[c];
            var kI = f.K.Select(k => k.KI).ToArray();
            var line = $"  front {c} {(f.Curve.Closed ? "closed" : f.Chain.Start + "-" + f.Chain.End)} knots {f.M} length {f.Curve.Length:G4}: KI {kI.Min():G5}..{kI.Max():G5} |KII| max {f.K.Max(k => Math.Abs(k.KII)):G4} |KIII| max {f.K.Max(k => Math.Abs(k.KIII)):G4} J/G mean {f.K.Select((k, j) => k.J / f.G[j]).Average():F4}";
            if (Ref(f, 0) != null)
            {
                var (mx, mean) = Box3D.Compare(f, j => Ref(f, j)!.Value.Item1);
                var jd = f.M / 2;
                var jS = 0;
                line += through
                    ? $" | Tada through crack 2c/W={cCur:F3}: largest diff {100 * mx:F2} %, mean {100 * mean:F2} %"
                    : $" | Newman-Raju: largest diff {100 * mx:F2} %, mean {100 * mean:F2} %; deepest {100 * (f.K[jd].KI - Ref(f, jd)!.Value.Item1) / Ref(f, jd)!.Value.Item1:+0.00;-0.00} %, surface {100 * (f.K[jS].KI - Ref(f, jS)!.Value.Item1) / Ref(f, jS)!.Value.Item1:+0.00;-0.00} %";
            }

            Console.WriteLine(line);
        }

        if (outDir.Length == 0) return;
        Box3D.WriteKnots(Path.Combine(outDir, "knots.txt"), step, drv, Ref);
        drv.WriteFronts(Path.Combine(outDir, $"fronts_{step:D3}.txt"));
        drv.WriteCrackFacets(Path.Combine(outDir, $"facets_{step:D3}.txt"));
        if (rec != null)
            File.AppendAllText(Path.Combine(outDir, "steps.txt"),
                $"{step} {rec.DaMax:G6} {aCur:G6} {cCur:G6} {rec.EnergyBefore:G10} {rec.EnergyAfter:G10} {rec.SweptWork:G8} {cc.SheetCount} {cc.Chains.Count} {cc.MouthEdges.Count} {drv.Snap.ElementCount} {rec.Iterations} {rec.MsSolve:F0} {string.Join(";", rec.Events)}\n");
    }
}

// Penny.cs — penny cracks in a cube: one inclined (Sneddon, Kassir and Sih) or two coplanar (coalescence).
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
//   penny kind=one|two [a=0.2 beta=0 gap=0.2 n=20 div=4 rq= steps= da= m= plateau= out= ens= check nofit]
//
// Cube [-1, 1]^3, E = 1, nu = 0.3, tension 1 along z, six point constraints (Step2, Section 8).
// "one": a penny crack of radius a through the origin with the normal (sin beta, 0, cos beta):
// K_I against Sneddon (beta = 0) or K_I, K_II, K_III against Kassir and Sih; a closed front
// recovered by counting, hats along the resampled front, no free-surface term involved.
// "two": two coplanar penny cracks of radius a with centres (-(a + gap/2), 0) and (a + gap/2, 0)
// in the plane z = 0.37 h_min; the inner points see a raised K_I; Paris growth brings the inner
// fronts together, the increments snap across the ligament, and the counting reports the
// coalescence: two sheets and two closed fronts become one sheet with one closed front.

using Step2;

namespace Step3.Bench;

public static class Penny
{
    private static readonly System.Globalization.CultureInfo Ci = System.Globalization.CultureInfo.InvariantCulture;

    public static void Run(Args a)
    {
        var kind = a.S("kind", "one");
        var n = a.I("n", 20);
        var div = a.I("div", 4);
        var h0 = 2.0 / n;
        var hmin = h0 / div;
        var aR = a.D("a", 0.2);
        var beta = a.D("beta", 0) * Math.PI / 180;
        var gap = a.D("gap", 0.2);
        var rq = a.D("rq", aR / 3);
        var sigma = 1.0;
        var mat = new Material(1, 0.3);
        var steps = a.I("steps", 0);
        var outDir = a.S("out", "");
        if (outDir.Length > 0) Directory.CreateDirectory(outDir);
        var zOff = 0.37 * hmin;
        V3 lo = new(-1, -1, -1), hi = new(1, 1, 1);
        var pad = steps > 0 ? a.D("pad", 1.0) * aR : 0;

        // crack discs: centre, normal, radius; the region for the refinement is the disc plus the growth annulus
        var discs = new List<(V3 C, V3 N, double R, int Bit)>();
        if (kind == "one")
        {
            // the centre (option centre=x,y,z in units of h_min) is 0.37 h above the lattice plane z = 0; an offset in x and y
            // as well made the first increment fail at every retry (2026-10-12), so x = y = 0 is the default
            var cs = a.S("centre", "0,0,0.37").Split(',');
            var c0 = new V3(double.Parse(cs[0], Ci) * hmin, double.Parse(cs[1], Ci) * hmin, double.Parse(cs[2], Ci) * hmin);
            discs.Add((c0, new V3(Math.Sin(beta), 0, Math.Cos(beta)).Unit(), aR, 1));
        }
        else if (kind == "list")
        {
            // any number of discs, inserted in order: discs=cx,cy,cz,nx,ny,nz,r:cx,... (":" or ";" between discs; absolute coordinates; a centre on a
            // face gives an edge crack, a later disc crossing an earlier one a junction). 2026-10-12, crack interaction example.
            var bit = 1;
            foreach (var d in a.S("discs", "").Split(new[] { ';', ':' }, StringSplitOptions.RemoveEmptyEntries))
            {
                var v = d.Split(',').Select(t => double.Parse(t.Trim(), Ci)).ToArray();
                if (v.Length != 7) throw new ArgumentException("discs: each entry is cx,cy,cz,nx,ny,nz,r");
                discs.Add((new V3(v[0], v[1], v[2]), new V3(v[3], v[4], v[5]).Unit(), v[6], bit++));
            }

            if (discs.Count == 0) throw new ArgumentException("kind=list needs discs=");
        }
        else
        {
            // second disc: radius a2 (default a) and centre offset dy along y (default 0): unequal, offset coalescence
            var a2 = a.D("a2", aR);
            var dy = a.D("dy", 0);
            discs.Add((new V3(-(aR + gap / 2), -dy / 2, zOff), new V3(0, 0, 1), aR, 1));
            discs.Add((new V3(a2 + gap / 2, dy / 2, zOff), new V3(0, 0, 1), a2, 2));
        }

        double Dist(V3 x)
        {
            var d = double.MaxValue;
            foreach (var (c, nrm, r, _) in discs)
            {
                var e = x - c;
                var hgt = V3.Dot(e, nrm);
                var rad = e.Perp(nrm).Norm;
                // distance to the annulus r - pad .. r + pad of the disc plane (the band the fronts cross)
                var dr = Math.Max(0, Math.Max(rad - (r + pad), (r - pad) - rad));
                d = Math.Min(d, Math.Sqrt(dr * dr + hgt * hgt));
            }

            return d;
        }

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var cm = Box3D.Mesh(lo, hi, h0, hmin, rq, Dist, out var refiner);
        Console.WriteLine($"penny {kind}: a {aR} beta {a.D("beta", 0)} gap {gap} h0 {h0:G4} hmin {hmin:G5} rq {rq:G4}; mesh {cm.M.ElementCount} tets after {refiner.Bisections} bisections, {sw.ElapsedMilliseconds} ms");
        var opt = Box3D.Options(a, hmin, rq);
        var drv = new GrowthDriver(cm, mat, Box3D.TensionZ(lo, hi, sigma), opt);
        drv.Inserter.RejectSameSide = !a.Flag("samesideok"); // samesideok: a star whose components all lie on one side is not split instead of rejecting
        foreach (var (c, nrm, r, bit) in discs)
        {
            var margin = hmin / 4;
            HashSet<int>? onWall = null;
            HashSet<int>? nearWall = null;
            if (opt.Fit)
            {
                // wall: the cylinder rho = r around the disc axis, as a polygon of 128 segments
                var e1 = Math.Abs(nrm.X) < 0.9 ? V3.Cross(nrm, new V3(1, 0, 0)).Unit() : V3.Cross(nrm, new V3(0, 1, 0)).Unit();
                var e2 = V3.Cross(nrm, e1);
                const int m = 128;
                var p = new V3[m];
                var b = new V3[m];
                for (var k = 0; k < m; k++)
                {
                    var th = 2 * Math.PI * k / m;
                    p[k] = c + r * (Math.Cos(th) * e1 + Math.Sin(th) * e2);
                    b[k] = nrm;
                }

                var wall = FrontFitter.Wall(p, b, 2 * hmin, true);
                var fit = FrontFitter.Fit(cm, wall, hmin);
                drv.DebugWall = wall;
                onWall = fit.WallNodes;
                nearWall = fit.Unmoved;
                margin = opt.FitMarginElements * hmin;
                Console.WriteLine($"fit disc bit {bit}: {fit.Moved}/{fit.Candidates} moved, {fit.Rejected} rejected, {fit.EdgesCut} edges cut");
            }

            Field3 phiS = x => V3.Dot(x - c, nrm);
            Field3 phiR = x => (x - c).Perp(nrm).Norm - r - margin;
            // intended front: the circle with outward N and the normal B
            {
                var e1 = Math.Abs(nrm.X) < 0.9 ? V3.Cross(nrm, new V3(1, 0, 0)).Unit() : V3.Cross(nrm, new V3(0, 1, 0)).Unit();
                var e2 = V3.Cross(nrm, e1);
                const int mc = 256;
                var pc = new V3[mc];
                var nc = new V3[mc];
                var bc = new V3[mc];
                for (var k = 0; k < mc; k++)
                {
                    var th = 2 * Math.PI * k / mc;
                    nc[k] = Math.Cos(th) * e1 + Math.Sin(th) * e2;
                    pc[k] = c + r * nc[k];
                    bc[k] = nrm;
                }

                drv.InsertInitial(bit, phiS, phiR, c, nrm, onWall, hmin, nearWall, [new FrontCurve3(pc, nc, bc, true)]);
            }
        }

        if (a.Flag("check")) Check(drv);
        drv.Solve();
        drv.Current = drv.EvaluateFronts();
        Report(drv, 0, kind, aR, beta, sigma, mat.Nu, discs, outDir);
        var ensDir = a.S("ens", "");
        var ensSteps = new List<int>();
        void Ens(int step)
        {
            if (ensDir.Length == 0) return;
            EnsightGold.WriteStep(ensDir, "penny", step, drv.Snap, cm, drv.Count, drv.X, drv.U, mat, a.Flag("body"), (lo, hi), $"penny {kind}");
            ensSteps.Add(step);
            EnsightGold.WriteCase(ensDir, "penny", a.Flag("body"), ensSteps, ensSteps.Select(s => (double)s).ToList());
        }

        Ens(0);
        for (var st = 1; st <= steps; st++)
        {
            var rec = drv.Step();
            if (a.Flag("check")) Check(drv);
            Report(drv, st, kind, aR, beta, sigma, mat.Nu, discs, outDir, rec);
            Ens(st);
        }

        Console.WriteLine($"total {sw.Elapsed.TotalSeconds:F1} s");
    }

    private static void Check(GrowthDriver drv)
    {
        var r = CrackChecker.Verify(drv.Snap, drv.Mesh, drv.Count);
        Console.WriteLine($"check: {r.Pairs} crack facet pairs, {r.Front} front, {r.Mouth} mouth, {r.Junction} junction edges agree with the count");
    }

    private static void Report(GrowthDriver drv, int step, string kind, double aR, double beta, double sigma, double nu,
        List<(V3 C, V3 N, double R, int Bit)> discs, string outDir, StepRecord? rec = null)
    {
        var cc = drv.Count;
        var head = $"step {step}: {cc.SheetCount} sheets, {cc.Chains.Count} fronts, {cc.PairCount} crack facets";
        if (rec != null) head += $", da {rec.DaMax:G4}, dPi {rec.EnergyAfter - rec.EnergyBefore:G5}, -int G dA {-rec.SweptWork:G5}, events [{string.Join("; ", rec.Events)}]";
        Console.WriteLine(head);
        var k0 = Reference.Sneddon(aR, sigma);
        (double, double, double)? Ref(FrontFactors f, int j)
        {
            if (kind != "one" || step > 0) return null;
            var (c, nrm, r, _) = discs[0];
            // theta from the shear direction: the shear on the crack plane is along the projection of z onto the plane
            var shear = new V3(0, 0, 1).Perp(nrm);
            var p = f.Curve.P[j] - c;
            var e1 = shear.Norm > 1e-12 ? shear.Unit() : (Math.Abs(nrm.Z) > 0.9 ? new V3(1, 0, 0) : new V3(0, 0, 1).Perp(nrm).Unit());
            var e2 = V3.Cross(nrm, e1);
            var th = Math.Atan2(V3.Dot(p, e2), V3.Dot(p, e1));
            var (kI, kII, kIII) = Reference.KassirSih(r, sigma, beta, th, nu);
            return (kI, kII, kIII);
        }

        for (var c = 0; c < drv.Current.Count; c++)
        {
            var f = drv.Current[c];
            var kI = f.K.Select(k => k.KI).ToArray();
            var line = $"  front {c} sheet {f.Chain.Sheet} {(f.Curve.Closed ? "closed" : "open")} knots {f.M} length {f.Curve.Length:G4}: KI {kI.Min():G5}..{kI.Max():G5} (mean {kI.Average():G5}, {100 * kI.Average() / k0:F1} % of Sneddon) |KII| max {f.K.Max(k => Math.Abs(k.KII)):G4} |KIII| max {f.K.Max(k => Math.Abs(k.KIII)):G4} J/G mean {f.K.Select((k, j) => k.J / f.G[j]).Average():F4}";
            if (kind == "one" && step == 0)
            {
                double eI = 0, eII = 0, eIII = 0;
                for (var j = 0; j < f.M; j++)
                {
                    var r = Ref(f, j)!.Value;
                    eI = Math.Max(eI, Math.Abs(f.K[j].KI - r.Item1) / k0);
                    eII = Math.Max(eII, Math.Abs(Math.Abs(f.K[j].KII) - Math.Abs(r.Item2)) / k0);
                    eIII = Math.Max(eIII, Math.Abs(Math.Abs(f.K[j].KIII) - Math.Abs(r.Item3)) / k0);
                }

                line += $" | largest error vs {(beta == 0 ? "Sneddon" : "Kassir-Sih")} in % of K0: KI {100 * eI:F2} KII {100 * eII:F2} KIII {100 * eIII:F2}";
            }

            if (f.DctKI != null) line += $" | dct KI mean {f.DctKI.Where(x => !double.IsNaN(x)).DefaultIfEmpty(double.NaN).Average():G5}";
            Console.WriteLine(line);
        }

        if (outDir.Length == 0) return;
        Box3D.WriteKnots(Path.Combine(outDir, "knots.txt"), step, drv, Ref);
        drv.WriteFronts(Path.Combine(outDir, $"fronts_{step:D3}.txt"));
        drv.WriteCrackFacets(Path.Combine(outDir, $"facets_{step:D3}.txt"));
        if (rec != null)
            File.AppendAllText(Path.Combine(outDir, "steps.txt"),
                $"{step} {rec.DaMax:G6} {rec.EnergyBefore:G10} {rec.EnergyAfter:G10} {rec.SweptWork:G8} {cc.SheetCount} {cc.Chains.Count} {cc.PairCount} {drv.Snap.ElementCount} {rec.Iterations} {rec.MsSolve:F0} {rec.MsFactors:F0} {rec.MsInsert:F0} {string.Join(";", rec.Events)}\n");
    }
}

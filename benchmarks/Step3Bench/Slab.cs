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
    /// <summary>Maximum tangential stress equivalent factor of a plane state (the hoop factor at the MTS angle).</summary>
    private static double MtsEquivalent(double ki, double kii)
    {
        var th = Math.Abs(kii) < 1e-12 * Math.Max(1, Math.Abs(ki)) ? 0 : 2 * Math.Atan((ki - Math.Sign(kii) * Math.Sqrt(ki * ki + 8 * kii * kii)) / (4 * kii));
        var c = Math.Cos(th / 2);
        return c * (ki * c * c - 1.5 * kii * Math.Sin(th));
    }

    private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;
    private const double HoleX = -4.0, HoleR = 0.25;
    private static readonly double[] HoleY = [2.75, 0.75, -1.25];

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
                // crack of half-length a through the origin, angle beta between the crack and the load axis (y); the centre
                // is offset so that the crack line avoids the lattice lines and diagonals
                var e = new V3(Math.Sin(beta), Math.Cos(beta), 0);
                var c = new V3(0.37 * hmin, 0.11 * hmin, 0);
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
                var c = new V3(0.37 * hmin, 0.11 * hmin, 0);
                var kink = a.D("kink", 53.13) * Math.PI / 180;
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
                var c = new V3(0.37 * hmin, 0.11 * hmin, 0);
                cracks.Add(new LineCrack { P0 = c - aLen * e, P1 = c + aLen * e, Bit = 2, Name = "inclined" });
                bands.Add((new V3(-bLen - pad, -aLen - pad, -1), new V3(bLen + pad, y1 + 0.3 * aLen + pad, 1)));
                break;
            }
            case "beam":
            {
                // three-point-bending beam with three holes (Bittencourt et al. 1996; Ingraffea and Grigoriu): W = 20, H = 8 (in),
                // supports at x = -+9 on the bottom face, load P at x = 0 on the top face, holes of radius 0.25 at x = -4,
                // y = 2.75, 0.75, -1.25; notch of length a from the bottom face at x = -d (cases d = 5, a = 1.5 and d = 6, a = 1)
                var d = a.D("d", 5.0);
                var xN = -d + 0.37 * hmin;
                cracks.Add(new LineCrack { P0 = new V3(xN, -H / 2 - 2 * h0, 0), P1 = new V3(xN, -H / 2 + aLen, 0), Bit = 1, Name = "notch" });
                bands.Add((new V3(xN - 0.4 * aLen - pad, -H / 2, -1), new V3(-2.5, 3.6, 1)));
                foreach (var hy in HoleY)
                    bands.Add((new V3(HoleX - HoleR - 3 * hmin, hy - HoleR - 3 * hmin, -1), new V3(HoleX + HoleR + 3 * hmin, hy + HoleR + 3 * hmin, 1)));
                break;
            }
            case "nooru":
            {
                // Nooru-Mohamed double-edge-notched panel (Areias et al. 2009, Sect. 4.6): W = H = 200, notches of depth a
                // at mid-height from both faces; bit 1 left, bit 2 right; the load programme is in Bc (bc=nooru)
                cracks.Add(new LineCrack { P0 = new V3(-W / 2 - 2 * h0, yOff, 0), P1 = new V3(-W / 2 + aLen, yOff, 0), Bit = 1, Name = "left" });
                cracks.Add(new LineCrack { P0 = new V3(W / 2 - aLen, yOff, 0), P1 = new V3(W / 2 + 2 * h0, yOff, 0), Bit = 2, Name = "right" });
                bands.Add((new V3(-W / 2, -0.25 * H, -1), new V3(W / 2, 0.25 * H, 1)));
                break;
            }
            case "echelon":
            {
                // two parallel cracks of half-length a under tension (Baud and Reuschle 1997, Fig. 7): normal spacing c between
                // the crack lines, longitudinal distance d between the inner tips; bit 1 lower left, bit 2 upper right
                var c = a.D("c", 0.5 * aLen);
                var d = a.D("d", aLen);
                cracks.Add(new LineCrack { P0 = new V3(-2 * aLen - d / 2, -c / 2 + yOff, 0), P1 = new V3(-d / 2, -c / 2 + yOff, 0), Bit = 1, Name = "lower" });
                cracks.Add(new LineCrack { P0 = new V3(d / 2, c / 2 + yOff, 0), P1 = new V3(2 * aLen + d / 2, c / 2 + yOff, 0), Bit = 2, Name = "upper" });
                bands.Add((new V3(-2 * aLen - d / 2 - pad, -c / 2 - 0.6 * aLen - pad, -1), new V3(2 * aLen + d / 2 + pad, c / 2 + 0.6 * aLen + pad, 1)));
                break;
            }
            default:
                throw new ArgumentException($"unknown kind {kind}");
        }

        var bc = a.S("bc", kind == "beam" ? "bend" : kind == "nooru" ? "nooru" : "tension");
        var loadCase = 0; // nooru: 0 = unit frame displacement in x, 1 = unit frame displacement in y
        var gap = a.D("gap", 5.0);
        var useBands = grow || kind == "beam";
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

            if (useBands) d = Math.Min(d, Size(x) - hmin);
            return hmin + 0.35 * d;
        }

        refiner.RefineToSize(SizeRq);
        var snap0 = cm.Commit();
        if (kind == "beam")
        {
            // the holes: the nodes within half an element of a circle are moved radially onto it (star kept valid), then the
            // tetrahedra whose centroid lies inside the circle are removed, which leaves a polygonal hole through the moved
            // nodes (smooth=0 keeps the stair-step boundary of the plain removal)
            var removed = 0;
            var snapped = 0;
            if (a.D("smooth", 1) > 0)
            {
                var m = cm.M;
                for (var v = 0; v < snap0.NodeCount; v++)
                {
                    if (m.TransposePairs(v).Length == 0) continue;
                    var p = cm.P(v);
                    foreach (var hy in HoleY)
                    {
                        var r = Math.Sqrt((p.X - HoleX) * (p.X - HoleX) + (p.Y - hy) * (p.Y - hy));
                        if (r < 1e-12 || Math.Abs(r - HoleR) > a.D("smooth", 1) * 0.5 * hmin) continue;
                        var target = new V3(HoleX + (p.X - HoleX) * HoleR / r, hy + (p.Y - hy) * HoleR / r, p.Z);
                        var ok = true;
                        var tr = m.TransposePairs(v);
                        for (var q = 0; q < tr.Length && ok; q += 2)
                        {
                            var row = m.Row(tr[q]);
                            if (row.Length != 4) continue;
                            V3[] pts = [cm.P(row[0]), cm.P(row[1]), cm.P(row[2]), cm.P(row[3])];
                            var before = Geo.Volume6(pts[0], pts[1], pts[2], pts[3]);
                            pts[row.IndexOf(v)] = target;
                            var after = Geo.Volume6(pts[0], pts[1], pts[2], pts[3]);
                            if (Math.Sign(after) != Math.Sign(before) || Math.Abs(after) < 0.15 * Math.Abs(before)) ok = false;
                        }

                        if (!ok) continue;
                        cm.X.Set(v, target.X);
                        cm.Y.Set(v, target.Y);
                        snapped++;
                        break;
                    }
                }

                snap0 = cm.Commit();
            }

            for (var e = 0; e < snap0.ElementCount; e++)
            {
                var row = snap0.Row(e);
                if (row.Length != 4) continue;
                var c = V3.Zero;
                foreach (var v in row) c += cm.P(v) / 4;
                foreach (var hy in HoleY)
                    if ((c.X - HoleX) * (c.X - HoleX) + (c.Y - hy) * (c.Y - hy) < HoleR * HoleR * (1 - 1e-9))
                    {
                        cm.M.RemoveElement(e);
                        removed++;
                        break;
                    }
            }

            snap0 = cm.Commit();
            Console.WriteLine($"holes: {snapped} nodes moved onto the circles, {removed} tetrahedra removed");
        }
        Console.WriteLine($"slab {kind}: W {W} H {H} t {h0:G4} h0 {h0:G4} hmin {hmin:G5}; mesh {snap0.ElementCount} tets {snap0.NodeCount} nodes after {refiner.Bisections} bisections, {sw.ElapsedMilliseconds} ms");

        var opt = new GrowthOptions
        {
            Through = !a.Flag("nothrough"), DaMax = a.D("da", 2 * hmin), DaMin = a.D("damin", hmin), TolRelative = a.D("tol", 0.06),
            Rq = rq, RqElements = a.D("rqel", 0), KnotSpacingElements = 2, MinKnots = 6, Fit = !a.Flag("nofit"),
            Paris = new ParisLaw { M = a.D("m", 3), C = a.D("C", 0) }, Dct = !a.Flag("nodct"), SolverTol = a.D("stol", 1e-10)
        };
        opt.Sif.Plateau = a.D("plateau", 0.5);
        opt.Sif.SurfaceFollowing = a.Flag("surf");
        opt.RejectDefects = !a.Flag("nodefect"); // nodefect: the topological acceptance test of a step is off (as in Box3D)
        opt.Advance.HermiteStrip = a.Flag("hermite");
        opt.Sif.FreeSurfaceTerm = !a.Flag("nofs");
        opt.Sif.CrackFaceTerm = !a.Flag("nocf");
        opt.Advance.Direction = a.S("crit", "mts") switch
        {
            "schoellmann" => DirectionCriterion.Schoellmann, "richard" => DirectionCriterion.Richard, "planar" => DirectionCriterion.Planar,
            _ => DirectionCriterion.MaximumTangentialStress
        };
        opt.Advance.KinkMax = a.D("kinkmax", 70) * Math.PI / 180;
        opt.Advance.CurvatureElements = a.D("curv", 0);
        opt.Advance.SnapElements = a.D("snap", 1.0);
        opt.Advance.OvershootElements = a.D("over", 0.6);

        BoundaryConditions Bc(Snapshot s, double[] X)
        {
            var f = new double[3 * s.NodeCount];
            var fixedDof = new bool[3 * s.NodeCount];
            var kinds = new Dictionary<int, FacetKind>();
            const double tol = 1e-9;
            BoundaryConditions.Mark(s, X, (c, nrm) => Math.Abs(c.Z) < tol || Math.Abs(c.Z - h0) < tol, FacetKind.Symmetry, kinds);
            switch (bc)
            {
                case "shear":
                    // bottom face clamped, shear traction tau = sigma along x on the top face (Yau et al.; Moes et al. 1999)
                    BoundaryConditions.Traction(s, X, (c, nrm) => Math.Abs(c.Y - H / 2) < tol, new V3(sigma, 0, 0), f, kinds);
                    for (var v = 0; v < s.NodeCount; v++)
                    {
                        if (s.TransposePairs(v).Length == 0) continue;
                        double z = X[3 * v + 2];
                        if (Math.Abs(z) < tol || Math.Abs(z - h0) < tol) fixedDof[3 * v + 2] = true;
                        if (Math.Abs(X[3 * v + 1] + H / 2) < tol) fixedDof[3 * v] = fixedDof[3 * v + 1] = true;
                    }

                    break;
                case "nooru":
                {
                    // bottom face and the lower part of the right face fixed; the top face and the upper part of the left face
                    // are glued to a rigid frame (no rotation) whose displacement is prescribed: (1, 0) for load case 0,
                    // (0, 1) for case 1; the shear force P and the vertical displacement v of the programme are combined
                    // from the two cases in CombineNooru
                    var pres = new double[3 * s.NodeCount];
                    for (var v = 0; v < s.NodeCount; v++)
                    {
                        if (s.TransposePairs(v).Length == 0) continue;
                        double x = X[3 * v], y = X[3 * v + 1], z = X[3 * v + 2];
                        if (Math.Abs(z) < tol || Math.Abs(z - h0) < tol) fixedDof[3 * v + 2] = true;
                        if (Math.Abs(y + H / 2) < tol || (Math.Abs(x - W / 2) < tol && y < -gap))
                            fixedDof[3 * v] = fixedDof[3 * v + 1] = true;
                        else if (Math.Abs(y - H / 2) < tol || (Math.Abs(x + W / 2) < tol && y > gap))
                        {
                            fixedDof[3 * v] = fixedDof[3 * v + 1] = true;
                            pres[3 * v] = loadCase == 0 ? 1 : 0;
                            pres[3 * v + 1] = loadCase == 0 ? 0 : 1;
                        }
                    }

                    return new BoundaryConditions { Force = f, FixedDof = fixedDof, FacetKinds = kinds, Prescribed = pres };
                }
                case "bend":
                {
                    // three-point bending: load P = sigma spread over |x| <= 0.25 on the top face, supports (u_y = 0) over
                    // |x -+ xs| <= h0/2 on the bottom face, u_x = 0 at the left support
                    var xs = a.D("xs", 9.0);
                    var pw = a.D("pw", 0.25);
                    BoundaryConditions.Traction(s, X, (c, nrm) => Math.Abs(c.Y - H / 2) < tol && Math.Abs(c.X) <= pw, new V3(0, -sigma / (2 * pw * h0), 0), f, kinds);
                    for (var v = 0; v < s.NodeCount; v++)
                    {
                        if (s.TransposePairs(v).Length == 0) continue;
                        double x = X[3 * v], y = X[3 * v + 1], z = X[3 * v + 2];
                        if (Math.Abs(z) < tol || Math.Abs(z - h0) < tol) fixedDof[3 * v + 2] = true;
                        if (Math.Abs(y + H / 2) < tol && Math.Abs(Math.Abs(x) - xs) <= h0 / 2 + tol)
                        {
                            fixedDof[3 * v + 1] = true;
                            if (x < 0) fixedDof[3 * v] = true;
                        }
                    }

                    break;
                }
                default:
                {
                    BoundaryConditions.Traction(s, X, (c, nrm) => Math.Abs(c.Y - H / 2) < tol, new V3(0, sigma, 0), f, kinds);
                    BoundaryConditions.Traction(s, X, (c, nrm) => Math.Abs(c.Y + H / 2) < tol, new V3(0, -sigma, 0), f, kinds);
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
                    break;
                }
            }

            return new BoundaryConditions { Force = f, FixedDof = fixedDof, FacetKinds = kinds };
        }

        var drv = new GrowthDriver(cm, mat, Bc, opt);
        if (kind == "nooru")
        {
            drv.Bc2Func = (s, X) =>
            {
                loadCase = 1;
                try { return Bc(s, X); }
                finally { loadCase = 0; }
            };
            drv.Combine = CombineNooru;
        }

        // Load programme of the Nooru-Mohamed test: the shear force P on the frame is kept and the vertical displacement v
        // grows until the maximum tangential stress factor of a tip reaches K_Ic. With u = u_x A + v B (A, B the unit
        // cases), P = u_x R_A + v R_B with R the x-reactions of the frame, so u_x(v) = (P - v R_B) / R_A, and K at a tip
        // is u_x K_A + v K_B; v is found by bisection, u_x follows.
        (double cA, double cB) CombineNooru(GrowthDriver d)
        {
            var X = d.X;
            var n = X.Length / 3;
            double RA = 0, RB = 0;
            for (var v = 0; v < n; v++)
            {
                double x = X[3 * v], y = X[3 * v + 1];
                if (Math.Abs(y - H / 2) < 1e-9 || (Math.Abs(x + W / 2) < 1e-9 && y > gap))
                {
                    RA += d.Reaction[3 * v];
                    RB += d.Reaction2[3 * v];
                }
            }

            var P = a.D("P", 10000.0) / a.D("thick", 50.0) * h0; // force on the slab of thickness h0
            var kic = a.D("kic", Math.Sqrt(mat.E * a.D("gf", 0.11) / (1 - mat.Nu * mat.Nu)));
            double Keq(double v, out double ux)
            {
                ux = (P - v * RB) / RA;
                var best = 0.0;
                for (var i = 0; i < d.CaseFactors[0].Count; i++)
                {
                    var ka = d.CaseFactors[0][i].K;
                    var kb = d.CaseFactors[1][i].K;
                    double ki = 0, kii = 0;
                    for (var j = 0; j < ka.Length; j++)
                    {
                        ki += (ux * ka[j].KI + v * kb[j].KI) / ka.Length;
                        kii += (ux * ka[j].KII + v * kb[j].KII) / ka.Length;
                    }

                    best = Math.Max(best, MtsEquivalent(ki, kii));
                }

                return best;
            }

            var k0 = Keq(0, out _);
            double vv;
            if (k0 >= kic) vv = 0;
            else
            {
                double lo = 0, hi = 1e-3;
                while (Keq(hi, out _) < kic && hi < 1e6) hi *= 2;
                for (var it = 0; it < 80; it++)
                {
                    var mid = 0.5 * (lo + hi);
                    if (Keq(mid, out _) < kic) lo = mid;
                    else hi = mid;
                }

                vv = hi;
            }

            var keq = Keq(vv, out var uxF);
            Console.WriteLine($"  nooru: R_xA {RA:G5} R_xB {RB:G5} P {P:G5}: u_x {uxF:G5} v {vv:G5} Keq {keq:G5} (KIc {kic:G5}, shear alone {k0:G5})");
            return (uxF, vv);
        }
        drv.Inserter.RejectSameSide = !a.Flag("samesideok"); // samesideok: a star whose components all lie on one side is not split instead of rejecting
        foreach (var cr in cracks)
        {
            var e = cr.Dir;
            var nrm = cr.Normal;
            var mid = 0.5 * (cr.P0 + cr.P1);
            var half = 0.5 * cr.Length;
            var margin = hmin / 4;
            HashSet<int>? onWall = null;
            HashSet<int>? nearWall = null;
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
                    drv.DebugWall = wall;
                    onWall.UnionWith(fit.WallNodes);
                    (nearWall ??= new HashSet<int>()).UnionWith(fit.Unmoved);
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
                var pts = cr.Polyline!;
                phiS = x => pl.SignedDistance(x).Phi;
                // arc coordinate along the polyline (in the plane), extrapolated beyond the ends
                phiR = x =>
                {
                    var best = double.MaxValue;
                    var sBest = 0.0;
                    var acc = 0.0;
                    for (var k = 0; k + 1 < pts.Length; k++)
                    {
                        V3 q0 = new(pts[k].X, pts[k].Y, 0), q1 = new(pts[k + 1].X, pts[k + 1].Y, 0);
                        var seg = q1 - q0;
                        var l2 = seg.Norm2;
                        var xp = new V3(x.X, x.Y, 0);
                        var t = V3.Dot(xp - q0, seg) / l2;
                        var tc = k == 0 ? Math.Min(t, 1) : k + 2 == pts.Length ? Math.Max(t, 0) : Math.Clamp(t, 0, 1);
                        var dd = (xp - q0 - tc * seg).Norm;
                        if (dd < best)
                        {
                            best = dd;
                            sBest = acc + tc * Math.Sqrt(l2);
                        }

                        acc += Math.Sqrt(l2);
                    }

                    return Math.Max(sBest - len, -sBest) - margin;
                };
                refNormal = poly.Normal(0);
                refPoint = 0.5 * (cr.Polyline![0] + cr.Polyline[1]);
            }
            else
            {
                phiS = x => V3.Dot(x - cr.P0, nrm);
                phiR = x => Math.Abs(V3.Dot(x - mid, e)) - half - margin;
            }

            // intended fronts: straight lines through the thickness at the tips inside the plate
            var intended = new List<FrontCurve3>();
            {
                var tipsI = cr.Polyline != null
                    ? new[] { (cr.Polyline[0], (cr.Polyline[0] - cr.Polyline[1]).Unit()), (cr.Polyline[^1], (cr.Polyline[^1] - cr.Polyline[^2]).Unit()) }
                    : new[] { (cr.P0, -e), (cr.P1, e) };
                foreach (var (tip, outward) in tipsI)
                {
                    if (Math.Abs(tip.X) > W / 2 || Math.Abs(tip.Y) > H / 2) continue;
                    var nTip = V3.Cross(outward, new V3(0, 0, 1)).Unit();
                    var pts2 = new V3[9];
                    var ns = new V3[9];
                    var bs = new V3[9];
                    for (var k = 0; k < 9; k++)
                    {
                        pts2[k] = new V3(tip.X, tip.Y, h0 * k / 8.0);
                        ns[k] = outward;
                        bs[k] = nTip;
                    }

                    intended.Add(new FrontCurve3(pts2, ns, bs, false, EndKind.Mouth, EndKind.Mouth));
                }
            }

            drv.InsertInitial(cr.Bit, phiS, phiR, refPoint, refNormal, onWall, hmin, nearWall, intended);
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
        drv.WriteCrackFacets(Path.Combine(outDir, $"facets_{step:D3}.txt"));
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
                if (a.S("bc", "tension") == "shear")
                {
                    // edge crack a = 3.5 in the 7 x 16 plate, bottom clamped, unit shear on the top: K_I = 34.0, K_II = 4.55
                    // (Yau, Wang and Corten 1980, as used by Moes, Dolbow and Belytschko 1999); the reference holds at the first front only
                    list.Add((34.0 * sigma, 4.55 * sigma, "Yau et al. shear plate (first front)"));
                    return list;
                }

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
                var kink = a.D("kink", 53.13) * Math.PI / 180;
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
                var kinked = (kIk, kIIk, $"Cotterell-Rice kink={a.D("kink", 53.13)} l/a={a.D("kl", 0.2)}");
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

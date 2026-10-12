// Box3D.cs — shared setup of the three-dimensional examples: box mesh, tension along z, corner constraints.
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

using System.Globalization;
using Step2;

namespace Step3.Bench;

public static class Box3D
{
    public static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

    /// <summary>Kuhn box refined to h = hmin + 0.35 max(0, d - rq) with d the distance given by <paramref name="dist"/>.</summary>
    public static CrackMesh Mesh(V3 lo, V3 hi, double h0, double hmin, double rq, Func<V3, double> dist, out MaubachRefiner refiner)
    {
        int nx = (int)Math.Round((hi.X - lo.X) / h0), ny = (int)Math.Round((hi.Y - lo.Y) / h0), nz = (int)Math.Round((hi.Z - lo.Z) / h0);
        var cm = KuhnMesh.Build(nx, ny, nz, lo, hi);
        refiner = new MaubachRefiner(cm);
        refiner.RefineToSize(x => hmin + 0.35 * Math.Max(0, dist(x) - rq));
        cm.Commit();
        return cm;
    }

    /// <summary>Tension sigma along z on the faces z = lo.Z and z = hi.Z; six point constraints at three corners of z = lo.Z.</summary>
    public static Func<Snapshot, double[], BoundaryConditions> TensionZ(V3 lo, V3 hi, double sigma)
    {
        return (s, X) =>
        {
            var f = new double[3 * s.NodeCount];
            var fixedDof = new bool[3 * s.NodeCount];
            var kinds = new Dictionary<int, FacetKind>();
            const double tol = 1e-9;
            BoundaryConditions.Traction(s, X, (c, _) => Math.Abs(c.Z - hi.Z) < tol, new V3(0, 0, sigma), f, kinds);
            BoundaryConditions.Traction(s, X, (c, _) => Math.Abs(c.Z - lo.Z) < tol, new V3(0, 0, -sigma), f, kinds);
            int Nearest(V3 p)
            {
                var best = -1;
                var bd = double.MaxValue;
                for (var v = 0; v < s.NodeCount; v++)
                {
                    if (s.TransposePairs(v).Length == 0) continue;
                    var d = (V3.Of(X, v) - p).Norm2;
                    if (d < bd)
                    {
                        bd = d;
                        best = v;
                    }
                }

                return best;
            }

            int a = Nearest(new V3(lo.X, lo.Y, lo.Z)), b = Nearest(new V3(hi.X, lo.Y, lo.Z)), c = Nearest(new V3(lo.X, hi.Y, lo.Z));
            fixedDof[3 * a] = fixedDof[3 * a + 1] = fixedDof[3 * a + 2] = true;
            fixedDof[3 * b + 1] = fixedDof[3 * b + 2] = true;
            fixedDof[3 * c + 2] = true;
            return new BoundaryConditions { Force = f, FixedDof = fixedDof, FacetKinds = kinds };
        };
    }

    public static GrowthOptions Options(Args a, double hmin, double rq)
    {
        var opt = new GrowthOptions
        {
            DaMax = a.D("da", 2 * hmin), DaMin = a.D("damin", hmin), TolRelative = a.D("tol", 0.04), Rq = rq, RqElements = a.D("rqel", 0),
            KnotSpacingElements = a.D("knot", 2), MinKnots = a.I("knots", 16), Paris = new ParisLaw { M = a.D("m", 3), C = a.D("C", 0) },
            Dct = !a.Flag("nodct"), SolverTol = a.D("stol", 1e-10), Fit = !a.Flag("nofit"), Through = a.Flag("through")
        };
        opt.Sif.SurfaceFollowing = a.Flag("surf");
        opt.RejectDefects = !a.Flag("nodefect");
        opt.GrowBits = a.I("growbits", 0); // growbits=1: only crack 1 grows
        opt.Timing = a.Flag("timing"); // time the from-scratch recomputation after every committed increment
        opt.Advance.HermiteStrip = a.Flag("hermite");
        opt.Sif.Plateau = a.D("plateau", 0.5);
        opt.Sif.FreeSurfaceTerm = !a.Flag("nofs");
        opt.Sif.CrackFaceTerm = !a.Flag("nocf");
        opt.Advance.Direction = a.S("crit", "mts") switch
        {
            "schoellmann" => DirectionCriterion.Schoellmann, "richard" => DirectionCriterion.Richard, "planar" => DirectionCriterion.Planar,
            _ => DirectionCriterion.MaximumTangentialStress
        };
        opt.Advance.KinkMax = a.D("kinkmax", 45) * Math.PI / 180;
        opt.Advance.CurvatureElements = a.D("curv", 0);
        opt.Advance.SnapElements = a.D("snap", 1.0);
        opt.Advance.OvershootElements = a.D("over", 0.6);
        opt.Advance.SmoothKink = a.I("smooth", 3);
        opt.Advance.SmoothAdvance = a.I("dasmooth", 6);
        opt.Advance.SmoothFront = a.I("fsmooth", 3);
        return opt;
    }

    /// <summary>Appends one line per knot of every front: step, chain, knot, x, y, z, KI, KII, KIII, J, G, Kv, kink, ref.</summary>
    public static void WriteKnots(string path, int step, GrowthDriver drv, Func<FrontFactors, int, (double KI, double KII, double KIII)?> reference)
    {
        using var w = new StreamWriter(path, true);
        for (var c = 0; c < drv.Current.Count; c++)
        {
            var f = drv.Current[c];
            for (var j = 0; j < f.M; j++)
            {
                var k = f.K[j];
                var p = f.Curve.P[j];
                var r = reference(f, j);
                w.WriteLine(string.Join(" ", new[]
                {
                    step.ToString(Ci), c.ToString(Ci), j.ToString(Ci), f.Curve.S[j].ToString("G8", Ci), p.X.ToString("G8", Ci), p.Y.ToString("G8", Ci), p.Z.ToString("G8", Ci),
                    k.KI.ToString("G8", Ci), k.KII.ToString("G8", Ci), k.KIII.ToString("G8", Ci), k.J.ToString("G8", Ci), f.G[j].ToString("G8", Ci), f.Kv[j].ToString("G8", Ci),
                    (f.Kink[j] * 180 / Math.PI).ToString("G6", Ci), (r?.KI ?? double.NaN).ToString("G8", Ci), (r?.KII ?? double.NaN).ToString("G8", Ci), (r?.KIII ?? double.NaN).ToString("G8", Ci),
                    (f.DctKI?[j] ?? double.NaN).ToString("G8", Ci), (f.DctKII?[j] ?? double.NaN).ToString("G8", Ci), (f.DctKIII?[j] ?? double.NaN).ToString("G8", Ci)
                }));
            }
        }
    }

    /// <summary>Largest relative difference of KI from the reference over the knots of a front, and the mean.</summary>
    public static (double Max, double Mean) Compare(FrontFactors f, Func<int, double> refKI)
    {
        double max = 0, sum = 0;
        for (var j = 0; j < f.M; j++)
        {
            var r = refKI(j);
            var e = Math.Abs(f.K[j].KI - r) / Math.Abs(r);
            max = Math.Max(max, e);
            sum += e;
        }

        return (max, sum / f.M);
    }
}

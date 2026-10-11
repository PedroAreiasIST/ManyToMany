// GrowthDriver.cs — one growth step: factors, advance, increments, commit, recount, solve, acceptance, energy balance.
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
// A step evaluates the factors of every front chain (interaction integral, with the
// displacement-correlation estimate and the J-integral as checks), builds the advance of
// each chain by the direction criterion and the Paris law relative to the largest
// equivalent factor of the step, inserts one increment per chain (walk from the chain's
// nodes), commits, recounts the crack facets and fronts, reports the topological events,
// assembles and solves from the previous displacements, and evaluates the new factors. The
// step is rolled back and repeated with half the advance when the relative factors of a
// chain change by more than a tolerance, and repeated with a modified increment (promotion
// tolerance, level-set shift) when an insertion is rejected. The energy balance compares
// the decrease of the potential energy with the swept-area integral of the energy release
// rate averaged between the two fronts.

using System.Diagnostics;
using System.Globalization;
using Step2;

namespace Step3;

public sealed class GrowthOptions
{
    public double DaMax = 0.02, DaMin = 0.005, DaGrowth = 1.5;
    public double TolRelative = 0.04;
    public int MaxRetries = 6;
    public double Rq = 0.05;              // tube radius (max(Rq, RqElements * h))
    public double RqElements = 6;
    public double KnotSpacingElements = 2; // knots of the resampled front per local element size
    public int MinKnots = 8;
    public SifOptions Sif = new();
    public AdvanceOptions Advance = new();
    public ParisLaw Paris = new();
    public bool Dct = true;
    public double[] DctRadiiElements = [1.5, 2.0, 2.5, 3.0, 4.0];
    public double SolverTol = 1e-10;
    public int SolverMaxIt = 20000;
    public bool Verbose = true;
    public TextWriter Log = Console.Out;
    public bool Through;                  // plane-strain slab: one factor per chain
    public bool Fit = true;               // fit the mesh to the new front before each increment
    public double FitWidthElements = 2;   // half-width of the front wall in element sizes
    public double FitMoveElements = 0.2;
    public double FitMarginElements = 0.01; // region margin beyond the fitted wall (roots on the wall are included by the wall rule)
    public double FaceTolElements = 0;    // optional geometric crack-face confirmation (0 = off; the topological rule applies)
    public double CoincidenceElements = 0.1; // other-crack nodes within this distance of the increment surface join it
}

public sealed class ChainReport
{
    public required FrontFactors F;
    public double KvMax, DaMax, KinkMean, Advance;
    public double DctKI = double.NaN, DctKII = double.NaN, DctKIII = double.NaN;
    public double JOverG = double.NaN;
}

public sealed class StepRecord
{
    public int Step;
    public double DaMax, Cycles, TotalCycles;
    public int Rollbacks, Retries;
    public List<string> Events = new();
    public List<ChainReport> Before = new();
    public double EnergyBefore, EnergyAfter, SweptWork;
    public int Elements, Nodes, Unknowns, Iterations, CrackFacets, Sheets, Chains;
    public double MsInsert, MsCommit, MsScan, MsAssemble, MsSolve, MsFactors;
    public string Summary = "";
}

public sealed class GrowthDriver
{
    public readonly CrackMesh Mesh;
    public readonly Material Mat;
    public readonly CrackInserter3 Inserter;
    public readonly CrackCount Count;
    public readonly GrowthOptions Opt;
    private readonly Func<Snapshot, double[], BoundaryConditions> _bc;
    public Snapshot Snap;
    public double[] X = [], U = [];
    public BoundaryConditions? Bc;
    public readonly FlatPattern Pattern = new();
    public double Energy = double.NaN;
    public double TotalCycles;
    public double Da;
    public int StepNo;
    public List<FrontFactors> Current = new();  // factors of the current version, evaluated by the last step
    public readonly List<StepRecord> History = new();
    private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

    public GrowthDriver(CrackMesh cm, Material mat, Func<Snapshot, double[], BoundaryConditions> bc, GrowthOptions opt)
    {
        Mesh = cm;
        Mat = mat;
        _bc = bc;
        Opt = opt;
        Inserter = new CrackInserter3(cm);
        Count = new CrackCount(cm);
        Snap = cm.M.Current;
        Da = opt.DaMax;
        if (opt.Through) opt.Advance.Through = true;
    }

    private void Log(string s)
    {
        if (Opt.Verbose) Opt.Log.WriteLine(s);
    }

    /// <summary>Inserts a crack by a global scan (first increment of a crack), commits and recounts.</summary>
    public IncrementResult InsertInitial(int bit, Field3 phiS, Field3 phiR, V3 refPoint, V3 refNormal, HashSet<int>? onWall = null)
    {
        Count.Reference[bit] = (refPoint, refNormal);
        var sw = Stopwatch.StartNew();
        Inserter.OnWall = onWall;
        var r = Inserter.Advance(bit, phiS, phiR, [], true);
        Inserter.OnWall = null;
        Snap = Mesh.Commit();
        Count.Rescan(Snap);
        Log($"insert bit {bit}: {r.Edges.Length} edges cut, {r.Promoted.Length} promoted, {r.Splits.Length} splits, {r.FrontNodes.Length} front nodes, {sw.ElapsedMilliseconds} ms");
        Log(Count.Summary());
        return r;
    }

    /// <summary>Assembles and solves the current version from the previous displacements.</summary>
    public StepRecord Solve(StepRecord? rec = null)
    {
        rec ??= new StepRecord { Step = StepNo };
        var sw = Stopwatch.StartNew();
        Pattern.Refresh(Snap);
        X = CrackMesh.Coordinates(Snap, Mesh);
        Bc = _bc(Snap, X);
        var vals = P1Elasticity.Assemble(Snap, Pattern, X, Mat);
        rec.MsAssemble = sw.Elapsed.TotalMilliseconds;
        sw.Restart();
        double[]? x0 = null;
        if (U.Length > 0)
        {
            // warm start: existing nodes keep their values, new ones interpolate from their parents
            x0 = new double[3 * Snap.NodeCount];
            Array.Copy(U, x0, Math.Min(U.Length, x0.Length));
            var p1 = Snap.View(Mesh.Parent1);
            var p2 = Snap.View(Mesh.Parent2);
            var wt = Snap.View(Mesh.Weight);
            for (var v = U.Length / 3; v < Snap.NodeCount; v++)
            {
                int a = p1[v], b = p2[v];
                var t = wt[v];
                if (a < 0 || b < 0) continue;
                for (var c = 0; c < 3; c++)
                    x0[3 * v + c] = a == b ? x0[3 * a + c] : (1 - t) * x0[3 * a + c] + t * x0[3 * b + c];
            }
        }

        var res = P1Elasticity.Pcg(Pattern, vals, Bc.Force, Bc.FixedDof, Opt.SolverTol, Opt.SolverMaxIt, x0);
        U = res.U;
        rec.MsSolve = sw.Elapsed.TotalMilliseconds;
        rec.Iterations = res.Iterations;
        Energy = BoundaryConditions.PotentialEnergy(Bc.Force, U);
        rec.EnergyAfter = Energy;
        rec.Elements = Snap.ElementCount;
        rec.Nodes = Snap.NodeCount;
        rec.Unknowns = Bc.FixedDof.Count(f => !f);
        if (!res.Converged) Log($"WARNING solver not converged: {res.Iterations} iterations, relative {res.Relative:E2}");
        Log($"solve: {Snap.ElementCount} elements, {Snap.NodeCount} nodes, {res.Iterations} iterations, rel {res.Relative:E2}, energy {Energy:G8}, {rec.MsAssemble:F0}+{rec.MsSolve:F0} ms");
        return rec;
    }

    private FacetKind FacetKindOf(int id)
    {
        if (Count.PairOfFacet.ContainsKey(id)) return FacetKind.Crack;
        return Bc != null && Bc.FacetKinds.TryGetValue(id, out var k) ? k : FacetKind.Free;
    }

    /// <summary>
    ///     Local element size at a chain: the mean longest edge of the uncut (Maubach-tagged) tetrahedra whose centroid
    ///     lies within <paramref name="radius"/> of the chain; the front edges themselves are edges of template children
    ///     and are shorter than the elements of the band. Falls back to the children when no uncut element is near.
    /// </summary>
    public double LocalSize(FrontChain c, double radius)
    {
        var bis = Snap.View(Mesh.Bisectable);
        double lo0 = double.MaxValue, lo1 = lo0, lo2 = lo0, hi0 = -lo0, hi1 = hi0, hi2 = hi0;
        foreach (var p in c.P)
        {
            lo0 = Math.Min(lo0, p.X); hi0 = Math.Max(hi0, p.X);
            lo1 = Math.Min(lo1, p.Y); hi1 = Math.Max(hi1, p.Y);
            lo2 = Math.Min(lo2, p.Z); hi2 = Math.Max(hi2, p.Z);
        }

        var curve = new FrontCurve3(c.P, c.N, c.B, c.Closed, c.Start, c.End);
        double sum = 0, sumAll = 0;
        int cnt = 0, cntAll = 0;
        for (var e = 0; e < Snap.ElementCount; e++)
        {
            var row = Snap.Row(e);
            if (row.Length != 4) continue;
            var cen = Geo.Centroid(X, row);
            if (cen.X < lo0 - radius || cen.X > hi0 + radius || cen.Y < lo1 - radius || cen.Y > hi1 + radius ||
                cen.Z < lo2 - radius || cen.Z > hi2 + radius) continue;
            if (curve.Nearest(cen).D2 > radius * radius) continue;
            var l2 = 0.0;
            for (var i = 0; i < 4; i++)
            for (var j = i + 1; j < 4; j++)
                l2 = Math.Max(l2, (V3.Of(X, row[i]) - V3.Of(X, row[j])).Norm2);
            var l = Math.Sqrt(l2);
            sumAll += l;
            cntAll++;
            if (bis[e] == 0) continue;
            sum += l;
            cnt++;
        }

        if (cnt > 0) return sum / cnt;
        return cntAll > 0 ? sumAll / cntAll : c.Length / Math.Max(1, c.Count - 1);
    }

    /// <summary>Factors of every front chain of the current version.</summary>
    public List<FrontFactors> EvaluateFronts(StepRecord? rec = null)
    {
        var sw = Stopwatch.StartNew();
        var list = new List<FrontFactors>();
        foreach (var chain in Count.Chains)
        {
            var h = LocalSize(chain, Math.Max(Opt.Rq, 2 * (chain.Length / Math.Max(1, chain.Count - 1))));
            var rq = Math.Max(Opt.Rq, Opt.RqElements * h);
            var curve = FrontCurve3.Resample(chain, Opt.KnotSpacingElements * h, Opt.MinKnots);
            var sheetSoup = Count.Soup(chain.Sheet);
            var ii = new InteractionIntegral(Snap, X, U, Mat, curve, rq, Opt.Sif, FacetKindOf, sheetSoup);
            if (Opt.Verbose) Log($"  front sheet {chain.Sheet}: h {h:G4} rq {rq:G4} knots {curve.M} tube {ii.TubeNodes} nodes {ii.TubeElements} elements, hat {ii.Delta:G4}");
            SifResult[] k;
            if (Opt.Through)
            {
                var one = ii.EvaluateThrough();
                k = new SifResult[curve.M];
                for (var j = 0; j < curve.M; j++)
                    k[j] = new SifResult { KI = one.KI, KII = one.KII, KIII = one.KIII, J = one.J, HatLength = one.HatLength, S = curve.S[j], FreeTerm = one.FreeTerm, CrackTerm = one.CrackTerm, Elements = one.Elements };
            }
            else k = ii.EvaluateAll();

            var bit = chain.Bits & -chain.Bits; // lowest bit
            var f = new FrontFactors
            {
                Chain = chain, Curve = curve, K = k, Kink = new double[curve.M], Twist = new double[curve.M], Kv = new double[curve.M],
                G = new double[curve.M], Rq = rq, H = h, Bit = bit
            };
            FrontAdvance.Evaluate(f, Opt.Advance, Mat.Nu, ii.Delta);
            for (var j = 0; j < curve.M; j++) f.G[j] = k[j].G(Mat);
            if (Opt.Dct)
            {
                var soup = sheetSoup;
                var radii = Opt.DctRadiiElements.Select(r => r * h).ToArray();
                f.DctKI = new double[curve.M];
                f.DctKII = new double[curve.M];
                f.DctKIII = new double[curve.M];
                for (var j = 0; j < curve.M; j++)
                {
                    var (a, b, c, _) = DisplacementCorrelation.At(curve, j, Count, soup, Snap, X, U, Mat, radii);
                    f.DctKI[j] = a;
                    f.DctKII[j] = b;
                    f.DctKIII[j] = c;
                }
            }

            list.Add(f);
        }

        if (rec != null) rec.MsFactors = sw.Elapsed.TotalMilliseconds;
        return list;
    }

    private static double Mean(IEnumerable<double> v)
    {
        var a = v.Where(x => !double.IsNaN(x)).ToArray();
        return a.Length == 0 ? double.NaN : a.Average();
    }

    private ChainReport Report(FrontFactors f)
    {
        var r = new ChainReport { F = f, KvMax = f.KvMax, KinkMean = f.Kink.Average() };
        if (f.DctKI != null)
        {
            r.DctKI = Mean(f.DctKI);
            r.DctKII = Mean(f.DctKII!);
            r.DctKIII = Mean(f.DctKIII!);
        }

        var g = f.G.Where(x => x > 0).ToArray();
        var j = f.K.Select(k => k.J).ToArray();
        if (g.Length > 0) r.JOverG = j.Average() / g.Average();
        return r;
    }

    /// <summary>Relative-factor change of a chain between two evaluations (same chain topology assumed).</summary>
    private static double RelativeChange(FrontFactors a, FrontFactors b)
    {
        if (a.M != b.M) return double.NaN;
        double ma = a.Kv.Average(), mb = b.Kv.Average();
        if (ma <= 0 || mb <= 0) return double.NaN;
        var d = 0.0;
        for (var j = 0; j < a.M; j++) d = Math.Max(d, Math.Abs(a.Kv[j] / ma - b.Kv[j] / mb));
        return d;
    }

    /// <summary>Checks that the current version is solved and its fronts evaluated (the first step needs it).</summary>
    public void Prepare()
    {
        if (U.Length != 3 * Snap.NodeCount) Solve();
        if (Current.Count != Count.Chains.Count) Current = EvaluateFronts();
    }

    /// <summary>One growth step with acceptance control; returns the record (also appended to History).</summary>
    public StepRecord Step()
    {
        Prepare();
        StepNo++;
        var rec = new StepRecord { Step = StepNo, EnergyBefore = Energy };
        var before = Current;
        foreach (var f in before) rec.Before.Add(Report(f));
        var kvMax = before.Count == 0 ? 0 : before.Max(f => f.KvMax);
        if (kvMax <= 0) throw new InvalidOperationException("no front to grow");
        var preSnap = Snap;
        var preU = U;
        var preEnergy = Energy;
        var prevSummary = new CountSummary(Count);
        var attempts = 0;
        while (true)
        {
            attempts++;
            var ok = TryStep(rec, before, kvMax, out var newFactors, out var failure);
            if (ok)
            {
                // acceptance: relative factors of chains that kept their topology
                var worst = 0.0;
                if (newFactors.Count == before.Count)
                    for (var i = 0; i < before.Count; i++)
                    {
                        var d = RelativeChange(before[i], newFactors[i]);
                        if (!double.IsNaN(d)) worst = Math.Max(worst, d);
                    }

                if (worst > Opt.TolRelative && Da > Opt.DaMin * 1.0001)
                {
                    Log($"step {StepNo}: relative factors changed by {worst:F3} > {Opt.TolRelative}: rollback, da {Da:G4} -> {Math.Max(Opt.DaMin, Da / 2):G4}");
                    Rollback(preSnap, preU, preEnergy, prevSummary);
                    Da = Math.Max(Opt.DaMin, Da / 2);
                    rec.Rollbacks++;
                    continue;
                }

                Current = newFactors;
                rec.DaMax = Da;
                var cycles = Opt.Paris.Cycles(Da, kvMax);
                TotalCycles += cycles;
                rec.Cycles = cycles;
                rec.TotalCycles = TotalCycles;
                rec.CrackFacets = Count.PairCount;
                rec.Sheets = Count.SheetCount;
                rec.Chains = Count.Chains.Count;
                Da = Math.Min(Opt.DaMax, Da * Opt.DaGrowth);
                rec.Summary = $"step {StepNo}: da {rec.DaMax:G4} rel-change {worst:F3} energy {rec.EnergyBefore:G7} -> {rec.EnergyAfter:G7} (dPi {rec.EnergyAfter - rec.EnergyBefore:G5}, -int G dA {-rec.SweptWork:G5}) events [{string.Join("; ", rec.Events)}]";
                Log(rec.Summary);
                History.Add(rec);
                return rec;
            }

            rec.Retries++;
            Log($"step {StepNo}: insertion rejected ({failure}); retry {attempts}");
            Rollback(preSnap, preU, preEnergy, prevSummary);
            if (attempts > Opt.MaxRetries)
            {
                if (Da > Opt.DaMin * 1.0001)
                {
                    Da = Math.Max(Opt.DaMin, Da / 2);
                    attempts = 0;
                    Log($"step {StepNo}: advance halved to {Da:G4}");
                    continue;
                }

                throw new InvalidOperationException($"step {StepNo}: insertion failed after {Opt.MaxRetries} retries at the minimum advance: {failure}");
            }
        }
    }

    private sealed class CountSummary(CrackCount c)
    {
        public readonly int Pairs = c.PairCount, Sheets = c.SheetCount, Junction = c.JunctionEdges.Count, Mouth = c.MouthEdges.Count;
        public readonly int Chains = c.Chains.Count, Closed = c.Chains.Count(ch => ch.Closed);
        public readonly List<(int Sheet, HashSet<(int, int, int)> Tris)> SheetTris = Enumerable.Range(0, c.SheetCount)
            .Select(s => (s, new HashSet<(int, int, int)>(Enumerable.Range(0, c.PairCount).Where(i => c.SheetOfPair[i] == s)
                .Select(i => Sorted(c.Tri[i]))))).ToList();

        public static (int, int, int) Sorted((int, int, int) t)
        {
            var (a, b, cc) = t;
            if (a > b) (a, b) = (b, a);
            if (b > cc) (b, cc) = (cc, b);
            if (a > b) (a, b) = (b, a);
            return (a, b, cc);
        }
    }

    private void Rollback(Snapshot pre, double[] preU, double preEnergy, CountSummary prev)
    {
        Mesh.Rollback(pre);
        Snap = pre;
        Count.Rescan(Snap);
        U = preU;
        Energy = preEnergy;
    }

    private readonly double[] _shiftSchedule = [0, 0.06, -0.06, 0.15, -0.15, 0.3, -0.3];
    private readonly double[] _tauSchedule = [0.05, 0.1, 0.05, 0.1, 0.2, 0.2, 0.2];

    private bool TryStep(StepRecord rec, List<FrontFactors> before, double kvMax, out List<FrontFactors> newFactors, out string failure)
    {
        newFactors = new List<FrontFactors>();
        failure = "";
        var sw = Stopwatch.StartNew();
        var attempt = rec.Retries;
        var shiftIdx = Math.Min(attempt, _shiftSchedule.Length - 1);
        Inserter.Tau = _tauSchedule[Math.Min(attempt, _tauSchedule.Length - 1)];
        var scale = attempt >= _shiftSchedule.Length ? 0.8 : 1.0;
        var outer = new AdvanceTarget { Soup = Count.OuterSoup(Snap), Name = "outer" };
        var sheetSoups = new Dictionary<int, TriangleSoup>();
        for (var s = 0; s < Count.SheetCount; s++) sheetSoups[s] = Count.Soup(s);
        var prevSummary = new CountSummary(Count);
        _lastAdvance.Clear();
        try
        {
            foreach (var f in before)
            {
                var chain = f.Chain;
                var h = f.H;
                var targets = new List<AdvanceTarget> { outer };
                foreach (var kv in sheetSoups)
                    if (kv.Key != chain.Sheet)
                        targets.Add(new AdvanceTarget { Soup = kv.Value, Name = $"sheet{kv.Key}", ExcludeFarSide = true });
                var others = Count.Chains.Where(c => !ReferenceEquals(c, chain)).ToList();
                var near = Count.Soup(chain.Sheet, p => chain.P.Any(q => V3.Dist(p, q) < 2 * f.Rq + Da));
                var geo = FrontAdvance.Build(f, kv => scale * Opt.Paris.Advance(Da, kv, kvMax), Opt.Advance, h, targets, others, near, scale * Da);
                var shift = _shiftSchedule[shiftIdx] * h;
                foreach (var tg in targets)
                    if (tg.ExcludeFarSide)
                    {
                        // side of the chain relative to the target, from the mid knot; a chain within the margin of the
                        // other sheet's plane (coplanar cracks about to coalesce) has no far side to exclude
                        var (t, cq, _, _) = tg.Soup.Nearest(chain.P[chain.Count / 2]);
                        if (t < 0) continue;
                        var d = V3.Dot(chain.P[chain.Count / 2] - cq, tg.Soup.Normal(t));
                        if (Math.Abs(d) > 2 * geo.Behind) geo.Excluded.Add((tg.Soup, d >= 0 ? 1 : -1));
                    }

                if (Opt.Fit)
                {
                    // wall through the new front along the rotated normal
                    var nb = new V3[chain.Count];
                    for (var i = 0; i < chain.Count; i++)
                        nb[i] = (Math.Cos(geo.DenseKink[i]) * chain.B[i] - Math.Sin(geo.DenseKink[i]) * chain.N[i]).Unit();
                    var wall = FrontFitter.Wall(geo.NewDense, nb, Opt.FitWidthElements * h, chain.Closed);
                    var fit = FrontFitter.Fit(Mesh, wall, h, Opt.FitMoveElements);
                    Inserter.OnWall = fit.WallNodes;
                    geo.Margin = Opt.FitMarginElements * h;
                    Log($"  fit: {fit.Moved} of {fit.Candidates} nodes moved ({fit.Rejected} rejected), {fit.EdgesCut} edges cut");
                }

                Field3 phiS = x => geo.PhiS(x) + shift;
                Field3 phiR = geo.PhiR;
                if (Environment.GetEnvironmentVariable("STEP3_DEBUG") == "1")
                {
                    var i0 = chain.Count / 2;
                    var dir0 = (geo.NewDense[i0] - chain.P[i0]).Unit();
                    var da0 = geo.DenseAdvance[i0];
                    Log($"  debug chain mid knot {chain.P[i0]} N {chain.N[i0]} B {chain.B[i0]} new {geo.NewDense[i0]} da {da0:G4} soup {geo.Surface.Count} tris");
                    foreach (var sfrac in new[] { -0.3, 0.25, 0.5, 0.75, 0.99, 1.2 })
                    foreach (var bfrac in new[] { -0.5, 0.0, 0.5 })
                    {
                        var x = chain.P[i0] + sfrac * da0 * dir0 + bfrac * h * chain.B[i0];
                        var (tri, cp, d2, _) = geo.Surface.Nearest(x);
                        Log($"    s {sfrac,5} b {bfrac,4}: phiS {geo.PhiS(x),9:G4} phiR {geo.PhiR(x),9:G4} tri {tri} tag {(tri >= 0 ? geo.Surface.Tag(tri) : -9)} d {Math.Sqrt(d2):G4}");
                    }
                }

                var seeds = Count.MeshNodes(chain);
                Inserter.FaceTolerance = Opt.FaceTolElements > 0 ? Opt.FaceTolElements * h : double.PositiveInfinity;
                Inserter.CoincidenceTolerance = Opt.CoincidenceElements * h;
                var r = Inserter.Advance(f.Bit, phiS, phiR, seeds);
                Inserter.OnWall = null;
                _lastAdvance.Add(geo.DenseAdvance);
                if (Environment.GetEnvironmentVariable("STEP3_DEBUG") == "2") LagDiagnostic(chain, geo, f.Bit, phiS, phiR, h);
                Log($"  chain sheet {chain.Sheet} bit {f.Bit}: Kv max {f.KvMax:G5}, kink mean {f.Kink.Average() * 180 / Math.PI:F1} deg, advance {geo.DenseAdvance.Max():G4}, snapped {geo.Snapped}; cut {r.Edges.Length} edges, {r.Promoted.Length} promoted, {r.Splits.Length} splits, visited {r.VisitedTets}");
                // swept work: sum over knots of G (before) times advance times arc weight (half of the trapezoid; the other half after)
                var swept = 0.0;
                for (var j = 0; j < f.M; j++)
                {
                    var w = f.Curve.Closed ? f.Curve.Length / f.M : (j == 0 || j == f.M - 1 ? 0.5 : 1) * f.Curve.Length / (f.M - 1);
                    // advance at the knot from the dense advance by arc length
                    var si = Math.Clamp((int)Math.Round(f.Curve.S[j] / chain.Length * (chain.Count - 1)), 0, chain.Count - 1);
                    swept += f.G[j] * geo.DenseAdvance[si] * w;
                }

                rec.SweptWork += 0.5 * swept;
            }

            rec.MsInsert = sw.Elapsed.TotalMilliseconds;
            sw.Restart();
            Snap = Mesh.Commit();
            rec.MsCommit = sw.Elapsed.TotalMilliseconds;
            sw.Restart();
            Count.Rescan(Snap);
            rec.MsScan = Count.ScanMs;
            rec.Events = EventsFrom(prevSummary);
            foreach (var ev in rec.Events) Log("  event: " + ev);
            if (Environment.GetEnvironmentVariable("STEP3_DEBUG") == "3")
                foreach (var ch in Count.Chains.Where(c => c.Count <= 12 && (c.Closed || c.Start == EndKind.Free || c.End == EndKind.Free)))
                {
                    Log($"  small chain sheet {ch.Sheet} {(ch.Closed ? "closed" : ch.Start + "-" + ch.End)} {ch.Count} nodes:");
                    foreach (var c in ch.Nodes)
                    {
                        var v = Count.Rep[c];
                        var copies = 0;
                        for (var w = 0; w < Mesh.NodeCount; w++)
                            if (Mesh.Canon[w] == c)
                                copies++;
                        Log($"    node {v} canon {c} copies {copies} mask {Mesh.Mask[v]} at {Mesh.P(v)} parents {Mesh.Parent1[v]},{Mesh.Parent2[v]} w {Mesh.Weight[v]:F3}");
                    }
                }
            MergeCoalescedBits();
            Log("  " + Count.Summary().Replace("\n", "\n  "));
            Solve(rec);
            newFactors = EvaluateFronts(rec);
            // second half of the swept work with the new factors on the chains that persisted
            if (newFactors.Count == before.Count)
            {
                var swept = 0.0;
                for (var i = 0; i < before.Count; i++)
                {
                    var f = newFactors[i];
                    var fb = before[i];
                    if (f.M != fb.M) continue;
                    for (var j = 0; j < f.M; j++)
                    {
                        var w = f.Curve.Closed ? f.Curve.Length / f.M : (j == 0 || j == f.M - 1 ? 0.5 : 1) * f.Curve.Length / (f.M - 1);
                        var si = Math.Clamp((int)Math.Round(fb.Curve.S[j] / fb.Chain.Length * (fb.Chain.Count - 1)), 0, fb.Chain.Count - 1);
                        swept += f.G[j] * _lastAdvance[i][si] * w;
                    }
                }

                rec.SweptWork += 0.5 * swept;
            }
            else rec.SweptWork *= 2; // topology changed: use the before-estimate only

            return true;
        }
        catch (InvalidOperationException ex)
        {
            failure = ex.Message.Split('\n')[0];
            return false;
        }
    }

    private readonly List<double[]> _lastAdvance = new();

    /// <summary>Debug: front nodes of the increment that lag behind the wall, with the crack-edge classification of their stars.</summary>
    private void LagDiagnostic(FrontChain chain, IncrementGeometry geo, int bit, Field3 phiS, Field3 phiR, double h)
    {
        var tmp = Mesh.Commit();
        var cc = new CrackCount(Mesh);
        foreach (var kv in Count.Reference) cc.Reference[kv.Key] = kv.Value;
        cc.Rescan(tmp);
        var wall = FrontFitter.Wall(geo.NewDense, chain.B, 1.0, chain.Closed);
        var lag = 0;
        foreach (var ch in cc.Chains)
        {
            if (ch.Sheet != chain.Sheet) continue;
            for (var j = 0; j < ch.Count; j++)
            {
                var d = wall.Distance(ch.P[j]);
                if (d < 0.05 * h || d > geo.DenseAdvance.Max() + h) continue;
                lag++;
                if (lag > 6) continue;
                var v = cc.Rep[ch.Nodes[j]];
                Log($"    lag node {v} canon {Mesh.Canon[v]} mask {Mesh.Mask[v]} at {ch.P[j]} dist to wall {d:G4} ({d / h:F2} h) kind {(ch.Closed ? "closed" : ch.Start + "-" + ch.End)}");
                var tr = Mesh.M.TransposePairs(v);
                for (var q = 0; q < tr.Length && q < 24; q += 2)
                {
                    var row = Mesh.M.Row(tr[q]).ToArray();
                    if (row.Length != 4) continue;
                    var sb = new System.Text.StringBuilder($"      tet {tr[q]}:");
                    foreach (var w in row)
                    {
                        var sg = (Mesh.Mask[w] & bit) != 0 ? 0 : phiS(Mesh.P(w)) >= 0 ? 1 : -1;
                        sb.Append($" {w}[m{Mesh.Mask[w]} s{sg} wall{wall.Distance(Mesh.P(w)) / h:F2}]");
                    }

                    foreach (var en in Numerical.SubdivisionTemplates.EdgeNodes)
                    {
                        int a = row[en[0]], b = row[en[1]];
                        int sa = (Mesh.Mask[a] & bit) != 0 ? 0 : phiS(Mesh.P(a)) >= 0 ? 1 : -1;
                        int sbb = (Mesh.Mask[b] & bit) != 0 ? 0 : phiS(Mesh.P(b)) >= 0 ? 1 : -1;
                        if (sa * sbb != -1) continue;
                        var t = Inserter.Root(Math.Min(a, b), Math.Max(a, b), phiS);
                        var x = V3.Lerp(Mesh.P(Math.Min(a, b)), Mesh.P(Math.Max(a, b)), t);
                        sb.Append($" edge({a},{b}) t {t:F3} phiR {phiR(x):G3}");
                    }

                    Log(sb.ToString());
                }
            }
        }

        Log($"    lag diagnostic: {lag} front nodes farther than 0.05 h from the wall");
    }

    private List<string> EventsFrom(CountSummary prev)
    {
        var ev = new List<string>();
        // sheet mapping through canonical triangles
        var curOf = new Dictionary<(int, int, int), int>();
        for (var i = 0; i < Count.PairCount; i++) curOf[CountSummary.Sorted(Count.Tri[i])] = Count.SheetOfPair[i];
        var prevToCur = new Dictionary<int, HashSet<int>>();
        var curToPrev = new Dictionary<int, HashSet<int>>();
        foreach (var (sheet, tris) in prev.SheetTris)
        foreach (var t in tris)
        {
            if (!curOf.TryGetValue(t, out var cs)) continue;
            if (!prevToCur.TryGetValue(sheet, out var a)) prevToCur[sheet] = a = new HashSet<int>();
            a.Add(cs);
            if (!curToPrev.TryGetValue(cs, out var b)) curToPrev[cs] = b = new HashSet<int>();
            b.Add(sheet);
        }

        foreach (var kv in curToPrev)
            if (kv.Value.Count > 1)
                ev.Add($"coalescence: sheets {string.Join(",", kv.Value.OrderBy(x => x))} -> sheet {kv.Key}");
        foreach (var kv in prevToCur)
            if (kv.Value.Count > 1)
                ev.Add($"split: sheet {kv.Key} -> sheets {string.Join(",", kv.Value.OrderBy(x => x))}");
        if (prev.Junction == 0 && Count.JunctionEdges.Count > 0) ev.Add($"junction: {Count.JunctionEdges.Count} edges");
        else if (Count.JunctionEdges.Count > prev.Junction) ev.Add($"junction +{Count.JunctionEdges.Count - prev.Junction} edges");
        if (prev.Mouth == 0 && Count.MouthEdges.Count > 0) ev.Add($"breakthrough: {Count.MouthEdges.Count} mouth edges");
        else if (Count.MouthEdges.Count > prev.Mouth) ev.Add($"mouth +{Count.MouthEdges.Count - prev.Mouth} edges");
        var closed = Count.Chains.Count(c => c.Closed);
        if (closed != prev.Closed) ev.Add($"closed fronts {prev.Closed} -> {closed}");
        if (Count.Chains.Count != prev.Chains) ev.Add($"fronts {prev.Chains} -> {Count.Chains.Count}");
        if (Count.SheetCount != prev.Sheets && curToPrev.All(kv => kv.Value.Count == 1) && prevToCur.All(kv => kv.Value.Count == 1))
            ev.Add($"sheets {prev.Sheets} -> {Count.SheetCount}");
        return ev;
    }

    /// <summary>After a coalescence the merged sheet carries two crack bits: every node of the sheet receives both.</summary>
    private void MergeCoalescedBits()
    {
        var changed = false;
        for (var sheet = 0; sheet < Count.SheetCount; sheet++)
        {
            var bits = Count.SheetBits[sheet];
            if ((bits & (bits - 1)) == 0) continue; // one bit
            foreach (var v in Count.SheetNodes(sheet))
                if ((Mesh.Mask[v] & bits) != 0 && (Mesh.Mask[v] & bits) != bits)
                {
                    Mesh.Mask.Set(v, Mesh.Mask[v] | bits);
                    changed = true;
                }
        }

        if (changed)
        {
            Snap = Mesh.Commit();
            Count.Rescan(Snap);
        }
    }

    /// <summary>Writes the factors of the current fronts: one line per knot (chain, knot, s, x, y, z, KI, KII, KIII, J, G, Kv, kink, DCT KI, KII, KIII).</summary>
    public void WriteFactors(string path)
    {
        using var w = new StreamWriter(path);
        w.WriteLine("# chain knot s x y z KI KII KIII J G Kv kink_deg dctKI dctKII dctKIII");
        for (var c = 0; c < Current.Count; c++)
        {
            var f = Current[c];
            for (var j = 0; j < f.M; j++)
            {
                var k = f.K[j];
                var p = f.Curve.P[j];
                w.WriteLine(string.Join(" ", new[]
                {
                    c.ToString(Ci), j.ToString(Ci), f.Curve.S[j].ToString("G8", Ci), p.X.ToString("G8", Ci), p.Y.ToString("G8", Ci), p.Z.ToString("G8", Ci),
                    k.KI.ToString("G8", Ci), k.KII.ToString("G8", Ci), k.KIII.ToString("G8", Ci), k.J.ToString("G8", Ci), f.G[j].ToString("G8", Ci),
                    f.Kv[j].ToString("G8", Ci), (f.Kink[j] * 180 / Math.PI).ToString("G6", Ci),
                    (f.DctKI?[j] ?? double.NaN).ToString("G8", Ci), (f.DctKII?[j] ?? double.NaN).ToString("G8", Ci), (f.DctKIII?[j] ?? double.NaN).ToString("G8", Ci)
                }));
            }
        }
    }

    /// <summary>Writes the fronts of the current version (chain, node, x, y, z).</summary>
    public void WriteFronts(string path)
    {
        using var w = new StreamWriter(path);
        w.WriteLine("# chain node x y z closed");
        for (var c = 0; c < Count.Chains.Count; c++)
        {
            var ch = Count.Chains[c];
            for (var j = 0; j < ch.Count; j++)
                w.WriteLine($"{c} {j} {ch.P[j].X.ToString("G8", Ci)} {ch.P[j].Y.ToString("G8", Ci)} {ch.P[j].Z.ToString("G8", Ci)} {(ch.Closed ? 1 : 0)}");
            w.WriteLine();
        }
    }
}

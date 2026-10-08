// Program.cs — verification and measurement driver for the step2 paper.
//
//   verify-fuzz  n batches batch seed   randomized edits, oracle after every batch,
//                                       snapshot immutability and rollback checks
//   verify-work  n                      crack / bisection / erosion with the oracle after
//                                       every step and conformity checks for bisection
//   (measurement modes live in Bench.cs)

using System.Diagnostics;
using Step2;
using Step2.Bench;

var mode = args.Length > 0 ? args[0] : "verify-work";
var rest = args.Skip(1).ToArray();
int Arg(int i, int dflt) => rest.Length > i ? int.Parse(rest[i]) : dflt;

switch (mode)
{
    case "verify-fuzz":
        Verify.Fuzz(Arg(0, 3), Arg(1, 2000), Arg(2, 5), Arg(3, 12345));
        break;
    case "verify-work":
        Verify.Workloads(Arg(0, 6));
        break;
    case "verify-uncut":
        UncutVerify.Run(Arg(0, 4), Arg(1, 10), Arg(2, 1));
        break;
    case "verify-crack":
        CrackVerify.Run(Arg(0, 8), Arg(1, 6), Arg(2, 0));
        break;
    case "verify-cut":
        CutVerify.Run(Arg(0, 4), Arg(1, 9), Arg(2, 1));
        break;
    case "penny":
        Fracture.Run(rest);
        break;
    case "growth":
        Fracture.Growth(rest);
        break;
    case "interfaces":
        Interfaces.Run(rest);
        break;
    case "verify-cutpar":
        ParallelCutVerify.Run(Arg(0, 6), Arg(1, 8), Arg(2, 1), Arg(3, 8));
        break;
    default:
        Bench.Run(mode, rest);
        break;
}

internal static class Verify
{
    public static MeshState NewState(int n, bool persistent, bool mirror = false)
    {
        var (rows, coords, nN) = Kuhn.Box(n);
        var flat = FlatBuild.Build(ElementSpec.Tetrahedron, nN, rows, parallel: false);
        var mesh = new IncrementalMesh(flat, ElementSpec.Tetrahedron, persistent);
        return new MeshState(mesh, coords, mirror ? new Mm3Mirror(rows) : null);
    }

    public static void Fuzz(int n, int batches, int batchSize, int seed)
    {
        var sw = Stopwatch.StartNew();
        var st = NewState(n, persistent: true);
        var rnd = new RandomEdits(seed);
        var pick = new Random(seed + 1);
        Oracle.Verify(st.Mesh.Current);
        var history = new List<(Snapshot S, ulong H)>();
        long edits = 0, rollbacks = 0, immutabilityChecks = 0;
        for (var b = 0; b < batches; b++)
        {
            for (var i = 0; i < batchSize; i++)
            {
                rnd.One(st, allowNonManifold: true);
                edits++;
            }

            var snap = st.Mesh.Commit();
            Oracle.Verify(snap);
            history.Add((snap, Oracle.Hash(snap) ^ AttrHash(st, snap)));
            if (history.Count > 64) history.RemoveAt(pick.Next(history.Count - 1));

            // published versions never change
            var (old, oldHash) = history[pick.Next(history.Count)];
            if ((Oracle.Hash(old) ^ AttrHash(st, old)) != oldHash)
                throw new Exception($"batch {b}: snapshot modified after publication");
            immutabilityChecks++;

            // occasionally roll back to an older version and continue from it
            if (b % 37 == 36)
            {
                // any version, including one on a branch abandoned by an earlier rollback:
                // topology and attributes (coordinates, tags) are restored together
                var (target, targetHash) = history[pick.Next(history.Count)];
                st.Mesh.Rollback(target);
                if ((Oracle.Hash(st.Mesh.Current) ^ AttrHash(st, st.Mesh.Current)) != targetHash)
                    throw new Exception("rollback changed the version");
                // the working state after rollback reads the restored attributes
                for (var v = 0; v < st.Mesh.NodeCount; v++)
                    if (st.X(v, 0) != st.Mesh.Current.View(st.Xa)[v])
                        throw new Exception("rollback did not restore coordinates");
                Oracle.Verify(st.Mesh.Current);
                rollbacks++;
            }
        }

        // the oracle also re-checks the working state after compaction of the pattern
        st.Mesh.CompactPattern();
        Oracle.Verify(st.Mesh.Commit());
        Console.WriteLine($"verify-fuzz n={n}: {batches} batches, {edits} edits, {rollbacks} rollbacks, " +
                          $"{immutabilityChecks} immutability checks, {st.Mesh.Relocations} relocations, " +
                          $"{st.Mesh.SlotReuses} slot reuses: all structures match the oracle " +
                          $"({sw.Elapsed.TotalSeconds:F1} s)");
    }

    /// <summary>Hash of the coordinate and tag attributes as published in a version.</summary>
    public static ulong AttrHash(MeshState st, Snapshot s)
    {
        ulong h = 0xcbf29ce484222325UL;
        var x = s.View(st.Xa);
        var tag = s.View(st.Tag);
        unchecked
        {
            for (var v = 0; v < s.NodeCount; v++) h = (h ^ (ulong)BitConverter.DoubleToInt64Bits(x[v])) * 1099511628211UL;
            for (var e = 0; e < s.ElementCount; e++) h = (h ^ (uint)tag[e]) * 1099511628211UL;
        }

        return h;
    }

    public static void Workloads(int n)
    {
        var sw = Stopwatch.StartNew();

        // crack insertion by node duplication
        {
            var st = NewState(n, persistent: true);
            var crack = new CrackPlane(n);
            for (var s = 0; s < crack.Steps; s++)
            {
                crack.Step(st);
                Oracle.Verify(st.Mesh.Commit());
            }

            Console.WriteLine($"crack: {crack.Steps} steps verified; boundary facets {st.Mesh.BoundaryCount}");
        }

        // recursive bisection in a moving ball, conformity checked through the boundary facets
        {
            var st = NewState(n, persistent: true);
            var bis = new Bisection();
            for (var s = 0; s < 8; s++)
            {
                var c = 0.2 + 0.6 * s / 7.0;
                var marked = Bisection.InBall(st, c, 0.5, 0.5, 0.25);
                bis.Refine(st, marked, maxLevel: 6);
                var snap = st.Mesh.Commit();
                Oracle.Verify(snap);
                CheckBoundaryOnBox(st, snap);
            }

            Console.WriteLine($"bisection: 8 steps verified and conforming; {bis.Bisections} bisections, " +
                              $"{st.Mesh.ElementCount} elements");
        }

        // erosion by a moving ball
        {
            var st = NewState(n, persistent: true);
            for (var s = 0; s < 8; s++)
            {
                var c = 0.1 + 0.8 * s / 7.0;
                Erosion.Step(st, Bisection.InBall(st, c, c, 0.5, 0.2));
                Oracle.Verify(st.Mesh.Commit());
            }

            Console.WriteLine($"erosion: 8 steps verified; boundary facets {st.Mesh.BoundaryCount}");
        }

        Console.WriteLine($"verify-work n={n}: all workloads match the oracle ({sw.Elapsed.TotalSeconds:F1} s)");
    }

    /// <summary>A conforming refinement of the unit box has boundary facets on the box faces only.</summary>
    public static void CheckBoundaryOnBox(MeshState st, Snapshot s)
    {
        for (var i = 0; i < s.BoundaryCount; i++)
        {
            var id = s.BoundaryFacet(i);
            int e = id / 4, f = id % 4;
            var row = s.Row(e);
            var onFace = false;
            for (var c = 0; c < 3 && !onFace; c++)
            foreach (var w in new[] { 0.0, 1.0 })
            {
                var all = true;
                for (var p = 0; p < 4; p++)
                    if (p != f && Math.Abs(st.X(row[p], c) - w) > 1e-12)
                    {
                        all = false;
                        break;
                    }

                if (all)
                {
                    onFace = true;
                    break;
                }
            }

            if (!onFace) throw new Exception($"non-conforming: interior boundary facet ({e},{f})");
        }
    }
}

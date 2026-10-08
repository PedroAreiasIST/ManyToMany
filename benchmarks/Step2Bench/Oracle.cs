// Oracle.cs — naive from-scratch evaluation and comparison against a snapshot.
//
// The oracle shares no code with FlatBuild or IncrementalMesh: it recomputes every
// derived structure from the element rows with dictionaries and brute-force facet
// matching, then compares field by field. Destinations and the pattern are checked
// through their defining invariants (slot order is implementation-defined): every
// destination addresses a live slot holding the right column, the live slots of a
// row are exactly the node's neighbours with the right pair counts, and no column
// occurs twice in a row. An assembled operator with integer element matrices is
// compared exactly against direct summation.

namespace Step2.Bench;

public static class Oracle
{
    public sealed class Mismatch(string message) : Exception(message);

    public static void Verify(Snapshot s, bool checkAssembly = true)
    {
        int K = s.NodesPerElement, Nf = s.FacetsPerElement;
        int nE = s.ElementCount, nN = s.NodeCount;

        // ---- transpose (element, position) lists, ascending by construction
        var tref = new List<(int E, int P)>[nN];
        for (var v = 0; v < nN; v++) tref[v] = new List<(int, int)>();
        for (var e = 0; e < nE; e++)
        {
            var row = s.Row(e);
            for (var p = 0; p < row.Length; p++) tref[row[p]].Add((e, p));
        }

        for (var v = 0; v < nN; v++)
        {
            var pairs = s.TransposePairs(v);
            if (pairs.Length != 2 * tref[v].Count)
                throw new Mismatch($"transpose row {v}: length {pairs.Length / 2} vs {tref[v].Count}");
            for (var q = 0; q < tref[v].Count; q++)
                if (pairs[2 * q] != tref[v][q].E || pairs[2 * q + 1] != tref[v][q].P)
                    throw new Mismatch($"transpose row {v} entry {q}: ({pairs[2 * q]},{pairs[2 * q + 1]}) vs {tref[v][q]}");
        }

        // ---- element positions
        for (var e = 0; e < nE; e++)
        {
            var row = s.Row(e);
            for (var p = 0; p < K; p++)
            {
                var expected = p < row.Length ? tref[row[p]].IndexOf((e, p)) : -1;
                if (s.ElementPosition(e, p) != expected)
                    throw new Mismatch($"EP({e},{p}) = {s.ElementPosition(e, p)} vs {expected}");
            }
        }

        // ---- node-node pair counts
        var pairCount = new Dictionary<int, int>[nN];
        for (var v = 0; v < nN; v++) pairCount[v] = new Dictionary<int, int>();
        for (var e = 0; e < nE; e++)
        {
            var row = s.Row(e);
            foreach (var a in row)
            foreach (var b in row)
                pairCount[a][b] = pairCount[a].GetValueOrDefault(b) + 1;
        }

        for (var v = 0; v < nN; v++)
        {
            var cols = s.PatternColumns(v);
            var cnts = s.PatternCounts(v);
            var seen = new HashSet<int>();
            var live = 0;
            for (var t = 0; t < cols.Length; t++)
            {
                if (!seen.Add(cols[t])) throw new Mismatch($"pattern row {v}: column {cols[t]} twice");
                if (cnts[t] < 0) throw new Mismatch($"pattern row {v}: negative count");
                if (cnts[t] == 0) continue;
                live++;
                if (pairCount[v].GetValueOrDefault(cols[t]) != cnts[t])
                    throw new Mismatch($"pattern ({v},{cols[t]}): count {cnts[t]} vs {pairCount[v].GetValueOrDefault(cols[t])}");
            }

            if (live != pairCount[v].Count)
                throw new Mismatch($"pattern row {v}: {live} live slots vs {pairCount[v].Count} neighbours");
        }

        // ---- destinations
        for (var e = 0; e < nE; e++)
        {
            var row = s.Row(e);
            if (row.Length == 0)
            {
                for (var t = 0; t < K * K; t++)
                    if (s.Destinations(e)[t] != -1)
                        throw new Mismatch($"destination of deleted element {e} not cleared");
                continue;
            }

            for (var i = 0; i < K; i++)
            for (var j = 0; j < K; j++)
                if ((i >= row.Length || j >= row.Length) && s.Destination(e, i, j) != -1)
                    throw new Mismatch($"destination ({e},{i},{j}) beyond the row not cleared");

            for (var i = 0; i < row.Length; i++)
            {
                var cols = s.PatternColumns(row[i]);
                var cnts = s.PatternCounts(row[i]);
                for (var j = 0; j < row.Length; j++)
                {
                    var slot = s.Destination(e, i, j);
                    if (slot < 0 || slot >= cols.Length || cols[slot] != row[j] || cnts[slot] <= 0)
                        throw new Mismatch($"destination ({e},{i},{j}) = {slot} invalid");
                }
            }
        }

        // ---- facet adjacency and boundary set
        var boundaryRef = new HashSet<int>();
        for (var e = 0; e < nE; e++)
        {
            var row = s.Row(e);
            for (var f = 0; f < Nf; f++)
            {
                var expected = -1;
                // a row of length k has k facets when k >= 2 (facet f opposite local node f)
                if (row.Length > 0 && f < (row.Length < 2 ? 0 : row.Length))
                {
                    // brute force: candidates from the first facet node, filtered by containment;
                    // facets are matched between rows of the same length (same dimension)
                    var facet = new List<int>();
                    for (var p = 0; p < row.Length; p++)
                        if (p != f) facet.Add(row[p]);
                    foreach (var (g, _) in tref[facet[0]])
                    {
                        if (g == e) continue;
                        var gr = s.Row(g);
                        if (gr.Length != row.Length) continue;
                        var all = true;
                        foreach (var x in facet)
                            if (gr.IndexOf(x) < 0)
                            {
                                all = false;
                                break;
                            }

                        if (all && (expected < 0 || g < expected)) expected = g;
                    }

                    if (expected < 0) boundaryRef.Add(e * Nf + f);
                }

                if (s.FacetNeighbour(e, f) != expected)
                    throw new Mismatch($"FA({e},{f}) = {s.FacetNeighbour(e, f)} vs {expected}");
            }
        }

        if (s.BoundaryCount != boundaryRef.Count)
            throw new Mismatch($"boundary count {s.BoundaryCount} vs {boundaryRef.Count}");
        var bset = new HashSet<int>();
        for (var i = 0; i < s.BoundaryCount; i++)
            if (!bset.Add(s.BoundaryFacet(i)))
                throw new Mismatch("boundary facet listed twice");
        if (!bset.SetEquals(boundaryRef)) throw new Mismatch("boundary facet sets differ");

        // ---- colouring: proper on the vertex-adjacency graph
        for (var e = 0; e < nE; e++)
        {
            var row = s.Row(e);
            if (row.Length == 0)
            {
                if (s.Color(e) != -1) throw new Mismatch($"deleted element {e} keeps a colour");
                continue;
            }

            if (s.Color(e) < 0) throw new Mismatch($"element {e} uncoloured");
            foreach (var v in row)
            foreach (var (g, _) in tref[v])
                if (g != e && s.Color(g) == s.Color(e))
                    throw new Mismatch($"elements {e} and {g} share node {v} and colour {s.Color(e)}");
        }

        if (checkAssembly) VerifyAssembly(s, tref);
    }

    /// <summary>
    ///     Assembles integer-valued element matrices through the destinations and compares
    ///     with direct summation; integer values make the comparison exact.
    /// </summary>
    private static void VerifyAssembly(Snapshot s, List<(int E, int P)>[] tref)
    {
        int K = s.NodesPerElement, nE = s.ElementCount, nN = s.NodeCount;
        var values = new double[s.TotalSlots];
        var direct = new Dictionary<long, double>(EdgeKeyComparer.Instance);
        for (var e = 0; e < nE; e++)
        {
            var row = s.Row(e);
            if (row.Length == 0) continue;
            for (var i = 0; i < row.Length; i++)
            {
                var rs = s.PatternRowStart(row[i]);
                for (var j = 0; j < row.Length; j++)
                {
                    double ke = 1 + (e * 7 + i * 3 + j) % 11;
                    values[rs + s.Destination(e, i, j)] += ke;
                    var key = (long)row[i] * nN + row[j];
                    direct[key] = direct.GetValueOrDefault(key) + ke;
                }
            }
        }

        var checkedEntries = 0;
        for (var v = 0; v < nN; v++)
        {
            var cols = s.PatternColumns(v);
            var cnts = s.PatternCounts(v);
            var rs = s.PatternRowStart(v);
            for (var t = 0; t < cols.Length; t++)
            {
                var expected = cnts[t] > 0 ? direct.GetValueOrDefault((long)v * nN + cols[t]) : 0.0;
                if (values[rs + t] != expected)
                    throw new Mismatch($"assembled ({v},{cols[t]}) = {values[rs + t]} vs {expected}");
                if (cnts[t] > 0) checkedEntries++;
            }
        }

        if (checkedEntries != direct.Count)
            throw new Mismatch($"assembled entries {checkedEntries} vs {direct.Count}");
    }

    /// <summary>Order-independent content hash of a snapshot (immutability checks).</summary>
    public static ulong Hash(Snapshot s)
    {
        ulong h = 1469598103934665603UL;
        void Mix(int x)
        {
            unchecked
            {
                h ^= (uint)x;
                h *= 1099511628211UL;
            }
        }

        Mix(s.ElementCount);
        Mix(s.NodeCount);
        for (var e = 0; e < s.ElementCount; e++)
        {
            foreach (var x in s.Row(e)) Mix(x);
            for (var p = 0; p < s.NodesPerElement; p++) Mix(s.ElementPosition(e, p));
            foreach (var x in s.Destinations(e)) Mix(x);
            for (var f = 0; f < s.FacetsPerElement; f++) Mix(s.FacetNeighbour(e, f));
            Mix(s.Color(e));
        }

        for (var v = 0; v < s.NodeCount; v++)
        {
            foreach (var x in s.TransposePairs(v)) Mix(x);
            foreach (var x in s.PatternColumns(v)) Mix(x);
            foreach (var x in s.PatternCounts(v)) Mix(x);
        }

        return h;
    }
}

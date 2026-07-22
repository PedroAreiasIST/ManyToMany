// Unit tests for the M2M layer -- the paper's BO2M (bidirectional one-to-many):
// a primary one-to-many map plus its lazily synchronized transpose and position
// caches, under a reader-writer lock.
using Numerical;

namespace UnitTests;

public class M2MTests
{
    private static M2M FromRows(params int[][] rows)
    {
        var o = new O2M(rows.Length);
        foreach (var r in rows) o.AppendElementCopy(r);
        return new M2M(o);
    }

    [Fact]
    public void ForwardLookupReturnsSuppliedOrder()
    {
        using var m = FromRows(new[] { 2, 4, 0 }, new[] { 0, 3 });
        Assert.Equal(new[] { 2, 4, 0 }, m.GetNodesForElement(0).ToArray());
    }

    [Fact]
    public void ReverseLookupUsesTheCachedTranspose()
    {
        using var m = FromRows(new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 0 });
        Assert.Equal(new[] { 0, 2 }, m.GetElementsForNode(0).ToArray());
        Assert.Equal(new[] { 0, 1 }, m.GetElementsForNode(1).ToArray());
        Assert.Equal(new[] { 1, 2 }, m.GetElementsForNode(2).ToArray());
    }

    [Fact]
    public void ReverseLookupRowsAreAscending()
    {
        var rng = new System.Random(5);
        var o = new O2M(2000);
        var seen = new HashSet<int>();
        var buf = new List<int>();
        for (var i = 0; i < 2000; i++)
        {
            seen.Clear(); buf.Clear();
            while (buf.Count < 4) { var v = rng.Next(300); if (seen.Add(v)) buf.Add(v); }
            o.AppendElementCopy(System.Runtime.InteropServices.CollectionsMarshal.AsSpan(buf));
        }
        using var m = new M2M(o);
        for (var n = 0; n <= o.GetMaxNode(); n++)
        {
            var row = m.GetElementsForNode(n);
            for (var i = 1; i < row.Count; i++)
                Assert.True(row[i - 1] < row[i], $"reverse row {n} not ascending");
        }
    }

    [Fact]
    public void ModificationInvalidatesAndLazilyRebuildsTheCache()
    {
        // The paper's lazy-synchronization contract: a structural modification
        // makes the derived caches stale; the next reverse query rebuilds them.
        using var m = FromRows(new[] { 0, 1 });
        Assert.Equal(new[] { 0 }, m.GetElementsForNode(1).ToArray());   // builds caches

        m.AppendElement(new List<int> { 1, 2 });                         // invalidates

        Assert.Equal(new[] { 0, 1 }, m.GetElementsForNode(1).ToArray()); // rebuilt
        Assert.Equal(new[] { 1 }, m.GetElementsForNode(2).ToArray());
    }

    [Fact]
    public void RepeatedQueriesBetweenModificationsAreStable()
    {
        using var m = FromRows(new[] { 0, 1 }, new[] { 1, 2 });
        var first = m.GetElementsForNode(1).ToArray();
        for (var i = 0; i < 5; i++)
            Assert.Equal(first, m.GetElementsForNode(1).ToArray());
    }

    [Fact]
    public void BatchUpdateDefersSynchronizationUntilScopeCloses()
    {
        using var m = FromRows(new[] { 0, 1 });
        using (m.BeginBatchUpdate())
        {
            m.AppendElement(new List<int> { 1, 2 });
            m.AppendElement(new List<int> { 2, 3 });
        }
        // After the scope closes, the next query resynchronizes and sees everything.
        Assert.Equal(new[] { 0, 1 }, m.GetElementsForNode(1).ToArray());
        Assert.Equal(new[] { 1, 2 }, m.GetElementsForNode(2).ToArray());
        Assert.Equal(new[] { 2 }, m.GetElementsForNode(3).ToArray());
    }

    [Fact]
    public void MembershipTestAgreesWithTheStoredRow()
    {
        using var m = FromRows(new[] { 4, 7 });
        Assert.True(m.ElementContainsNode(0, 7));
        Assert.False(m.ElementContainsNode(0, 5));
    }

    [Fact]
    public void ElementNeighboursAreThoseSharingAtLeastOneNode()
    {
        //  e0=[0,1]  e1=[1,2]  e2=[3,4]   -> e0~e1, e2 isolated
        using var m = FromRows(new[] { 0, 1 }, new[] { 1, 2 }, new[] { 3, 4 });
        var n0 = m.GetElementNeighbors(0);
        Assert.Contains(1, n0);
        Assert.DoesNotContain(2, n0);
        Assert.Empty(m.GetElementNeighbors(2).Where(x => x != 2));
    }

    [Fact]
    public void ElementToElementGraphIsReflexiveAndSymmetric()
    {
        using var m = FromRows(new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 });
        var ee = m.GetElementsToElements();
        for (var e = 0; e < ee.Count; e++)
        {
            Assert.Contains(e, ee[e].ToArray());
            foreach (var f in ee[e].ToArray())
                Assert.Contains(e, ee[f].ToArray());
        }
    }

    [Fact]
    public void NodeToNodeGraphIsTheAssemblySparsityPattern()
    {
        // Two nodes are coupled iff they appear together in some element.
        using var m = FromRows(new[] { 0, 1 }, new[] { 2, 3 });
        var nn = m.GetNodesToNodes();
        Assert.Contains(1, nn[0].ToArray());
        Assert.DoesNotContain(2, nn[0].ToArray());
        Assert.Contains(3, nn[2].ToArray());
    }

    [Fact]
    public void ConcurrentReadersObserveConsistentResults()
    {
        using var m = FromRows(new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 3, 0 });
        var expected = m.GetElementsForNode(1).ToArray();
        Parallel.For(0, 200, _ =>
        {
            Assert.Equal(expected, m.GetElementsForNode(1).ToArray());
        });
    }
}

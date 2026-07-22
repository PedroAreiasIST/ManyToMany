// Unit tests for the O2M layer (the paper's one-to-many sparse structure).
// Each test names the contract it pins down, so a change in behaviour that
// would invalidate a statement in the accompanying paper fails here first.
using System.Runtime.InteropServices;
using Numerical;

namespace UnitTests;

public class O2MTests
{
    private static O2M FromRows(params int[][] rows)
    {
        var o = new O2M(rows.Length);
        foreach (var r in rows) o.AppendElementCopy(r);
        return o;
    }

    private static int[][] ToJagged(O2M o)
    {
        var res = new int[o.Count][];
        for (var i = 0; i < o.Count; i++) res[i] = o[i].ToArray();
        return res;
    }

    /// <summary>
    ///     Random structure with rows that are duplicate-free but NOT sorted — the
    ///     shape of a real primary connectivity, where each element lists distinct
    ///     nodes in an orientation-carrying order.
    /// </summary>
    private static O2M Random(int rows, int perRow, int maxNode, int seed)
    {
        var rng = new System.Random(seed);
        var o = new O2M(rows);
        var seen = new HashSet<int>();
        var buf = new List<int>(perRow);
        for (var i = 0; i < rows; i++)
        {
            seen.Clear();
            buf.Clear();
            while (buf.Count < perRow)
            {
                var v = rng.Next(maxNode);
                if (seen.Add(v)) buf.Add(v);
            }
            o.AppendElementCopy(CollectionsMarshal.AsSpan(buf));
        }
        return o;
    }

    // ---------------------------------------------------------------- storage

    [Fact]
    public void PrimaryRowsPreserveSuppliedOrder()
    {
        // The paper relies on this: for elements the node order encodes orientation.
        var a = FromRows(new[] { 2, 4, 0, 6 }, new[] { 0, 3 });
        Assert.Equal(new[] { 2, 4, 0, 6 }, a[0].ToArray());
        Assert.Equal(new[] { 0, 3 }, a[1].ToArray());
    }

    [Fact]
    public void AppendReturnsIndexAndGrowsCount()
    {
        var a = new O2M();
        Assert.Equal(0, a.AppendElementCopy(new[] { 1, 2 }));
        Assert.Equal(1, a.AppendElementCopy(new[] { 3 }));
        Assert.Equal(2, a.Count);
    }

    // -------------------------------------------------------------- transpose

    [Fact]
    public void TransposeInvertsTheRelation()
    {
        // Worked example of the paper's motivating table.
        var a = FromRows(new[] { 2, 4, 0, 6 }, new[] { 0, 3 }, new[] { 1, 3 }, new[] { 3, 6, 1, 5, 4, 2 });
        var t = a.Transpose();

        for (var e = 0; e < a.Count; e++)
            foreach (var n in a[e].ToArray())
                Assert.Contains(e, t[n].ToArray());

        for (var n = 0; n < t.Count; n++)
            foreach (var e in t[n].ToArray())
                Assert.Contains(n, a[e].ToArray());
    }

    [Fact]
    public void TransposeRowsAreSortedByConstruction()
    {
        // Theorem "Transpose Correctness and Sorting" in the paper.
        foreach (var (rows, perRow, maxNode) in new[] { (50, 4, 30), (5000, 4, 512), (20000, 6, 4096) })
        {
            var t = Random(rows, perRow, maxNode, seed: 7).Transpose();
            for (var n = 0; n < t.Count; n++)
            {
                var row = t[n].ToArray();
                for (var i = 1; i < row.Length; i++)
                    Assert.True(row[i - 1] < row[i], $"transpose row {n} not strictly ascending");
            }
        }
    }

    [Fact]
    public void TransposeIsInvolutiveOnSortedInput()
    {
        var a = FromRows(new[] { 0, 1 }, new[] { 1, 2 }, new[] { 0, 2 });
        var back = a.Transpose().Transpose();
        Assert.Equal(ToJagged(a), ToJagged(back));
    }

    // --------------------------------------------------- set-operation orders

    [Fact]
    public void UnionPreservesPreEstablishedOrder()
    {
        // The reviewer's example: left-operand precedence, NOT sorted output.
        var u = FromRows(new[] { 2 }) | FromRows(new[] { 1 });
        Assert.Equal(new[] { 2, 1 }, u[0].ToArray());

        var v = FromRows(new[] { 5, 3, 9 }) | FromRows(new[] { 9, 1, 3 });
        Assert.Equal(new[] { 5, 3, 9, 1 }, v[0].ToArray());
    }

    [Fact]
    public void UnionRemovesDuplicatesWithinAndAcrossOperands()
    {
        var u = FromRows(new[] { 4, 4, 2 }) | FromRows(new[] { 2, 7, 7 });
        Assert.Equal(new[] { 4, 2, 7 }, u[0].ToArray());
    }

    [Fact]
    public void UnionPlusOperatorMatchesPipeOperator()
    {
        var a = Random(500, 4, 200, 11);
        var b = Random(500, 4, 200, 12);
        Assert.Equal(ToJagged(a | b), ToJagged(a + b));
    }

    [Theory]
    [InlineData(4)]      // short rows: hash path
    [InlineData(80)]     // long rows: sort-merge path
    public void IntersectionReturnsAscendingRows(int perRow)
    {
        var a = Random(400, perRow, 300, 21);
        var b = Random(400, perRow, 300, 22);
        var c = a & b;
        for (var i = 0; i < c.Count; i++)
        {
            var row = c[i].ToArray();
            for (var k = 1; k < row.Length; k++)
                Assert.True(row[k - 1] < row[k], $"intersection row {i} not ascending (perRow={perRow})");
        }
    }

    [Fact]
    public void IntersectionContainsExactlyTheCommonEntries()
    {
        var a = Random(300, 5, 60, 31);
        var b = Random(300, 5, 60, 32);
        var c = a & b;
        for (var i = 0; i < c.Count; i++)
        {
            var expected = a[i].ToArray().Intersect(b[i].ToArray()).OrderBy(x => x).ToArray();
            Assert.Equal(expected, c[i].ToArray());
        }
    }

    [Theory]
    [InlineData(4)]
    [InlineData(80)]
    public void DifferenceReturnsAscendingRows(int perRow)
    {
        var a = Random(400, perRow, 300, 41);
        var b = Random(400, perRow, 300, 42);
        var c = a - b;
        for (var i = 0; i < c.Count; i++)
        {
            var row = c[i].ToArray();
            for (var k = 1; k < row.Length; k++)
                Assert.True(row[k - 1] < row[k], $"difference row {i} not ascending (perRow={perRow})");
        }
    }

    [Fact]
    public void DifferenceOfRowWithoutCounterpartIsStillAscending()
    {
        // left has more rows than right: the tail rows are copied through.
        var a = FromRows(new[] { 5, 1, 3 }, new[] { 9, 2 });
        var b = FromRows(new[] { 1 });
        var c = a - b;
        Assert.Equal(new[] { 3, 5 }, c[0].ToArray());
        Assert.Equal(new[] { 2, 9 }, c[1].ToArray());
    }

    [Fact]
    public void SymmetricDifferenceKeepsEntriesOfExactlyOneOperand()
    {
        var x = FromRows(new[] { 1, 2, 3 }) ^ FromRows(new[] { 3, 4 });
        Assert.Equal(new[] { 1, 2, 4 }, x[0].ToArray().OrderBy(v => v).ToArray());
    }

    // ---------------------------------------------------------- multiplication

    [Fact]
    public void MultiplicationComposesTheTwoRelations()
    {
        var a = FromRows(new[] { 0, 1 }, new[] { 2 });          // element -> node
        var b = FromRows(new[] { 10 }, new[] { 11 }, new[] { 12 }); // node -> face
        var c = a * b;
        Assert.Equal(new[] { 10, 11 }, c[0].ToArray().OrderBy(v => v).ToArray());
        Assert.Equal(new[] { 12 }, c[1].ToArray());
    }

    [Fact]
    public void ProductRowsAreDuplicateFree()
    {
        var a = Random(600, 4, 120, 51);
        var p = a * a.Transpose();
        for (var i = 0; i < p.Count; i++)
        {
            var row = p[i].ToArray();
            Assert.Equal(row.Length, row.Distinct().Count());
        }
    }

    [Fact]
    public void ElementToElementProductIsReflexiveAndSymmetric()
    {
        var a = Random(400, 4, 100, 61);
        var ee = a * a.Transpose();
        for (var e = 0; e < ee.Count; e++)
        {
            Assert.Contains(e, ee[e].ToArray());                 // shares nodes with itself
            foreach (var f in ee[e].ToArray())
                Assert.Contains(e, ee[f].ToArray());             // symmetry
        }
    }

    [Fact]
    public void MultiplicationRejectsIncompatibleDimensions()
    {
        var a = FromRows(new[] { 0, 7 });
        var b = FromRows(new[] { 1 });     // only one row, but 'a' references node 7
        Assert.Throws<InvalidOperationException>(() => _ = a * b);
    }

    // ------------------------------------------- pooled-buffer reuse (regression)

    [Fact]
    public void RepeatedOperationsReuseBuffersWithoutCorruption()
    {
        // The set-operation kernels rent marker/touched buffers from a shared pool;
        // calls after the first receive dirty buffers. Results must not drift.
        foreach (var threshold in new[] { int.MaxValue, 1 })   // sequential and parallel kernels
        {
            var a = Random(4000, 4, 512, 71);
            var b = Random(4000, 4, 512, 72);
            a.ParallelizationThreshold = threshold;
            b.ParallelizationThreshold = threshold;

            var firstUnion = ToJagged(a | b);
            var firstXor = ToJagged(a ^ b);
            var firstProd = ToJagged(a * a.Transpose());

            for (var rep = 0; rep < 3; rep++)
            {
                Assert.Equal(firstUnion, ToJagged(a | b));
                Assert.Equal(firstXor, ToJagged(a ^ b));
                Assert.Equal(firstProd, ToJagged(a * a.Transpose()));
            }
        }
    }

    // ------------------------------------------------------------ maintenance

    [Fact]
    public void CompressElementsKeepsSelectedRowsInOrder()
    {
        var a = FromRows(new[] { 0, 1 }, new[] { 2, 3 }, new[] { 4, 5 });
        a.CompressElements(new List<int> { 2, 0 });
        Assert.Equal(2, a.Count);
        Assert.Equal(new[] { 4, 5 }, a[0].ToArray());
        Assert.Equal(new[] { 0, 1 }, a[1].ToArray());
    }

    [Fact]
    public void PermuteNodesRelabelsEveryReference()
    {
        var a = FromRows(new[] { 0, 1 }, new[] { 1, 2 });
        a.PermuteNodes(new List<int> { 10, 11, 12 });
        Assert.Equal(new[] { 10, 11 }, a[0].ToArray());
        Assert.Equal(new[] { 11, 12 }, a[1].ToArray());
    }

    [Fact]
    public void MaxNodeTracksContent()
    {
        var a = FromRows(new[] { 3, 9 }, new[] { 1 });
        Assert.Equal(9, a.GetMaxNode());
    }

    [Fact]
    public void ValidityDetectsDuplicateNodesWithinAnElement()
    {
        Assert.True(FromRows(new[] { 0, 1, 2 }).IsValid());
        Assert.False(FromRows(new[] { 0, 1, 1 }).IsValid());
    }

    // ------------------------------------------------------- topological sort

    [Fact]
    public void TopologicalOrderRespectsDependencies()
    {
        // 0 -> 1 -> 2 (row i lists the successors of vertex i)
        var dag = FromRows(new[] { 1 }, new[] { 2 }, System.Array.Empty<int>());
        var order = dag.GetTopOrder();
        var pos = new Dictionary<int, int>();
        for (var i = 0; i < order.Count; i++) pos[order[i]] = i;
        Assert.True(pos[0] < pos[1]);
        Assert.True(pos[1] < pos[2]);
    }

    [Fact]
    public void TopologicalOrderRaisesOnCycle()
    {
        // Documented behaviour: cycles raise rather than returning an empty list.
        var cyclic = FromRows(new[] { 1 }, new[] { 0 });
        Assert.Throws<InvalidOperationException>(() => cyclic.GetTopOrder());
    }

    // ----------------------------------------------------------------- CSR

    [Fact]
    public void CsrRoundTripPreservesContent()
    {
        var a = Random(200, 4, 90, 81);
        var (rowPtr, colIdx) = a.ToCsr();
        Assert.Equal(a.Count + 1, rowPtr.Length);
        Assert.Equal(rowPtr[^1], colIdx.Length);
        var back = O2M.FromCsr(rowPtr, colIdx);
        Assert.Equal(ToJagged(a), ToJagged(back));
    }
}

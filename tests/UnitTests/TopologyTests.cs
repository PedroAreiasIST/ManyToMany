// Unit tests for the MM2M layer (the paper's MBO2M multi-type matrix) and for
// the type-safe Topology API, including the finite-element support operations
// described in the paper.
using Numerical;

namespace UnitTests;

// Entity marker types for the type-safe API.
public readonly struct Vertex;
public readonly struct Cell;

public class MM2MTests
{
    [Fact]
    public void RelationsAreKeptPerTypePair()
    {
        using var mm = new MM2M(2);
        mm.AppendElement(1, 0, new List<int> { 0, 1 });   // type1 -> type0
        mm.AppendElement(1, 0, new List<int> { 1, 2 });
        Assert.Equal(2, mm[1, 0].Count);
        Assert.Equal(0, mm[0, 1].Count);                  // untouched pair stays empty
    }

    [Fact]
    public void ReverseQueriesWorkThroughTheSelectedPair()
    {
        using var mm = new MM2M(2);
        mm.AppendElement(1, 0, new List<int> { 0, 1 });
        mm.AppendElement(1, 0, new List<int> { 1, 2 });
        Assert.Equal(new[] { 0, 1 }, mm[1, 0].GetElementsForNode(1).ToArray());
    }

    [Fact]
    public void MarkAndCompressRemovesMarkedEntities()
    {
        using var mm = new MM2M(2);
        mm.AppendElement(1, 0, new List<int> { 0, 1 });
        mm.AppendElement(1, 0, new List<int> { 2, 3 });
        mm.AppendElement(1, 0, new List<int> { 4, 5 });

        mm.MarkToErase(1, 1);            // mark the middle element of type 1
        mm.Compress();

        Assert.Equal(2, mm[1, 0].Count);
        Assert.Equal(new[] { 0, 1 }, mm[1, 0].GetNodesForElement(0).ToArray());
        Assert.Equal(new[] { 4, 5 }, mm[1, 0].GetNodesForElement(1).ToArray());
    }

    [Fact]
    public void VersionAdvancesOnCompression()
    {
        using var mm = new MM2M(2);
        mm.AppendElement(1, 0, new List<int> { 0, 1 });
        mm.AppendElement(1, 0, new List<int> { 2, 3 });
        var before = mm.Version;
        mm.MarkToErase(1, 0);
        mm.Compress();
        Assert.True(mm.Version > before, "Compress must publish a new version");
    }
}

public class TopologyTests
{
    /// <summary>Structured n-by-n quadrilateral grid; returns the topology.</summary>
    private static Topology<TypeMap<Vertex, Cell>> Grid(int n)
    {
        var topo = new Topology<TypeMap<Vertex, Cell>>();
        for (var i = 0; i < (n + 1) * (n + 1); i++) topo.Add<Vertex>();
        for (var r = 0; r < n; r++)
        for (var c = 0; c < n; c++)
        {
            int v00 = r * (n + 1) + c, v10 = v00 + 1;
            int v01 = (r + 1) * (n + 1) + c, v11 = v01 + 1;
            topo.Add<Cell, Vertex>(v00, v10, v11, v01);
        }
        return topo;
    }

    [Fact]
    public void TypeSafeQueriesRoundTrip()
    {
        var topo = Grid(2);
        Assert.Equal(9, topo.Count<Vertex>());
        Assert.Equal(4, topo.Count<Cell>());

        foreach (var v in topo.NodesOf<Cell, Vertex>(0))
            Assert.Contains(0, topo.ElementsAt<Cell, Vertex>(v));
    }

    [Fact]
    public void CellNodeOrderIsPreserved()
    {
        var topo = new Topology<TypeMap<Vertex, Cell>>();
        for (var i = 0; i < 4; i++) topo.Add<Vertex>();
        topo.Add<Cell, Vertex>(2, 0, 3, 1);
        Assert.Equal(new[] { 2, 0, 3, 1 }, topo.NodesOf<Cell, Vertex>(0).ToArray());
    }

    [Fact]
    public void TransposeExposesIncidentCells()
    {
        var topo = Grid(2);
        var t = topo.GetTranspose<Cell, Vertex>();
        // The centre vertex of a 2x2 grid belongs to all four cells.
        Assert.Equal(4, t[4].Length);
    }

    [Fact]
    public void ElementGraphConnectsCellsSharingVertices()
    {
        var topo = Grid(2);
        var ee = topo.GetElementToElementGraph<Cell, Vertex>();
        for (var e = 0; e < ee.Count; e++)
            Assert.Equal(4, ee[e].Length);   // all four cells meet at the centre
    }

    [Fact]
    public void BoundaryNodesExcludeTheInterior()
    {
        var topo = Grid(2);
        var boundary = topo.FindBoundaryNodes<Cell, Vertex>(2);
        Assert.DoesNotContain(4, boundary);          // centre vertex is interior
        Assert.Contains(0, boundary);                // corner is on the boundary
        Assert.Equal(8, boundary.Count);             // 9 vertices minus the centre
    }

    [Fact]
    public void CuthillMcKeeReturnsAPermutation()
    {
        var topo = Grid(4);
        var perm = topo.ComputeCuthillMcKeeOrdering<Cell, Vertex>();
        Assert.Equal(topo.Count<Vertex>(), perm.Length);
        Assert.Equal(Enumerable.Range(0, perm.Length), perm.OrderBy(x => x));
    }

    [Fact]
    public void ColoringSeparatesCellsThatShareAVertex()
    {
        // The parallel-assembly guarantee stated in the paper.
        var topo = Grid(4);
        var colors = topo.ComputeElementColoring<Cell, Vertex>();
        var ee = topo.GetElementToElementGraph<Cell, Vertex>();
        for (var e = 0; e < ee.Count; e++)
            foreach (var f in ee[e].ToArray())
                if (f != e)
                    Assert.True(colors[e] != colors[f], $"cells {e},{f} share a vertex but share a colour");
    }

    [Fact]
    public void ColorGroupsPartitionTheCells()
    {
        var topo = Grid(3);
        var groups = topo.GetColorGroups<Cell, Vertex>();
        var all = groups.SelectMany(g => g).OrderBy(x => x).ToArray();
        Assert.Equal(Enumerable.Range(0, topo.Count<Cell>()), all);
    }

    [Fact]
    public void SparsityPatternRowsAreSortedAndSymmetric()
    {
        var topo = Grid(3);
        var (rowPtr, colIdx) = topo.GetSparsityPatternCSR<Cell, Vertex>();
        var n = topo.Count<Vertex>();
        Assert.Equal(n + 1, rowPtr.Length);
        Assert.Equal(rowPtr[^1], colIdx.Length);

        var pairs = new HashSet<(int, int)>();
        for (var r = 0; r < n; r++)
        {
            for (var k = rowPtr[r] + 1; k < rowPtr[r + 1]; k++)
                Assert.True(colIdx[k - 1] < colIdx[k], $"CSR row {r} not ascending");
            for (var k = rowPtr[r]; k < rowPtr[r + 1]; k++)
                pairs.Add((r, colIdx[k]));
        }
        foreach (var (r, c) in pairs)
            Assert.Contains((c, r), pairs);       // structurally symmetric
    }

    [Fact]
    public void SparsityPatternExpandsWithDegreesOfFreedom()
    {
        var topo = Grid(2);
        var (rp1, ci1) = topo.GetSparsityPatternCSR<Cell, Vertex>();
        var (rp2, ci2) = topo.GetSparsityPatternCSR<Cell, Vertex>(2);
        Assert.Equal(2 * (rp1.Length - 1) + 1, rp2.Length);
        Assert.Equal(4 * ci1.Length, ci2.Length);   // d x d block per node pair
    }

    [Fact]
    public void SharedCountIdentifiesEdgeNeighbours()
    {
        var topo = Grid(2);
        // In a quad grid, two cells sharing an edge have exactly two common vertices.
        var edgeNeighbours = topo.GetEntitiesWithSharedCount<Cell, Vertex>(0, 2);
        Assert.Equal(2, edgeNeighbours.Count);
    }

    [Fact]
    public void KHopNeighbourhoodGrowsWithDepth()
    {
        var topo = Grid(4);
        var hop1 = topo.GetKHopNeighborhood<Cell, Vertex>(0, 1);
        var hop2 = topo.GetKHopNeighborhood<Cell, Vertex>(0, 2);
        Assert.True(hop2.Count > hop1.Count);
        foreach (var kv in hop1) Assert.True(hop2.ContainsKey(kv.Key));
    }

    [Fact]
    public void IntegrityValidationAcceptsAWellFormedMesh()
    {
        var topo = Grid(3);
        Assert.True(topo.ValidateIntegrity<Cell, Vertex>().IsValid);
    }

    [Fact]
    public void CanonicalFormIdentifiesPermutedElements()
    {
        // Symmetry-based deduplication: the same triangle in a different node
        // order canonicalizes to the same key.
        var sym = Symmetry.Full(3);
        Span<int> a = stackalloc int[3];
        Span<int> b = stackalloc int[3];
        sym.CanonicalSpan(new[] { 7, 3, 5 }, a);
        sym.CanonicalSpan(new[] { 5, 7, 3 }, b);
        Assert.True(a.SequenceEqual(b));
        Assert.Equal(new[] { 3, 5, 7 }, a.ToArray());     // lexicographic minimum
    }

    [Fact]
    public void CyclicSymmetryDistinguishesReflections()
    {
        var cyclic = Symmetry.Cyclic(3);
        Span<int> a = stackalloc int[3];
        Span<int> b = stackalloc int[3];
        cyclic.CanonicalSpan(new[] { 1, 2, 3 }, a);       // rotations only
        cyclic.CanonicalSpan(new[] { 3, 2, 1 }, b);       // a reflection
        Assert.False(a.SequenceEqual(b), "orientation must survive cyclic canonicalization");
    }

    [Fact]
    public void SymmetryGroupSizesMatchTheStandardGroups()
    {
        Assert.Equal(1, Symmetry.Identity(4).Permutations.Count);
        Assert.Equal(4, Symmetry.Cyclic(4).Permutations.Count);
        Assert.Equal(8, Symmetry.Dihedral(4).Permutations.Count);
        Assert.Equal(24, Symmetry.Full(4).Permutations.Count);
    }
}

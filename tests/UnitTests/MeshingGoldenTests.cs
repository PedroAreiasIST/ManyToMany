// Golden-master tests for the Meshing layer.
//
// Each test runs a deterministic pipeline (edge refinement, smoothing, crack insertion)
// and reduces its complete output - element rows in order, node lineage, parent-element
// records and coordinates bit for bit - to a 64-bit FNV-1a fingerprint. The expected
// values were recorded on the code before the deduplication of the subdivision templates
// and of the SimplexRemesher boundary/neighbour helpers, so a refactor that changes any
// node label, element order or coordinate bit fails here.
using Numerical;

namespace UnitTests;

public class MeshingGoldenTests
{
    // ------------------------------------------------------------------ fingerprinting

    private sealed class Fnv
    {
        private ulong _h = 0xcbf29ce484222325UL;

        public void Add(long v)
        {
            unchecked
            {
                for (var i = 0; i < 8; i++)
                {
                    _h ^= (byte)(v >> (8 * i));
                    _h *= 1099511628211UL;
                }
            }
        }

        public void Add(double d) => Add(BitConverter.DoubleToInt64Bits(d));

        public string Hex => _h.ToString("x16");
    }

    private static void AddRows<TElem>(Fnv h, SimplexMesh m, bool withOrigin)
    {
        var n = m.Count<TElem>();
        h.Add(n);
        for (var e = 0; e < n; e++)
        {
            var row = m.NodesOf<TElem, Node>(e);
            h.Add(row.Count);
            for (var k = 0; k < row.Count; k++) h.Add(row[k]);
            if (withOrigin) h.Add(m.Get<TElem, OriginalElement>(e).Index);
        }
    }

    private static void AddMesh(Fnv h, SimplexMesh m, bool withLineage)
    {
        var nn = m.Count<Node>();
        h.Add(nn);
        if (withLineage)
            for (var v = 0; v < nn; v++)
            {
                var p = m.Get<Node, ParentNodes>(v);
                h.Add(p.Parent1);
                h.Add(p.Parent2);
            }

        AddRows<Point>(h, m, withLineage);
        AddRows<Bar2>(h, m, withLineage);
        AddRows<Tri3>(h, m, withLineage);
        AddRows<Tet4>(h, m, withLineage);
    }

    private static void AddCoords(Fnv h, double[,] x)
    {
        h.Add(x.GetLength(0));
        h.Add(x.GetLength(1));
        foreach (var v in x) h.Add(v);
    }

    private static void AddSet(Fnv h, HashSet<int> s)
    {
        h.Add(s.Count);
        foreach (var v in s.OrderBy(v => v)) h.Add(v);
    }

    // ------------------------------------------------------------------ inputs

    private static List<(int, int)> UniqueEdges<TElem>(SimplexMesh m)
    {
        var seen = new HashSet<(int, int)>();
        var list = new List<(int, int)>();
        for (var e = 0; e < m.Count<TElem>(); e++)
        {
            var r = m.NodesOf<TElem, Node>(e);
            for (var i = 0; i < r.Count; i++)
            for (var j = i + 1; j < r.Count; j++)
            {
                var key = r[i] < r[j] ? (r[i], r[j]) : (r[j], r[i]);
                if (seen.Add(key)) list.Add((r[i], r[j]));
            }
        }

        return list;
    }

    private static List<(int, int)> Pick(List<(int, int)> edges, double p, int seed)
    {
        var rnd = new Random(seed);
        var marked = edges.Where(_ => rnd.NextDouble() < p).ToList();
        // shuffle so that midpoint labels are not ordered like the edges
        for (var i = marked.Count - 1; i > 0; i--)
        {
            var j = rnd.Next(i + 1);
            (marked[i], marked[j]) = (marked[j], marked[i]);
        }

        return marked;
    }

    private static void Perturb(double[,] x, int dims, double amp, int seed)
    {
        var rnd = new Random(seed);
        for (var v = 0; v < x.GetLength(0); v++)
        for (var c = 0; c < dims; c++)
        {
            var w = x[v, c];
            if (w > 1e-12 && w < 1 - 1e-12) x[v, c] = w + amp * (2 * rnd.NextDouble() - 1);
        }
    }

    private static readonly int[][] Perms4 = Permutations([0, 1, 2, 3]).ToArray();

    private static IEnumerable<int[]> Permutations(int[] a)
    {
        if (a.Length <= 1)
        {
            yield return a;
            yield break;
        }

        for (var i = 0; i < a.Length; i++)
        {
            var rest = a.Where((_, k) => k != i).ToArray();
            foreach (var p in Permutations(rest)) yield return [a[i], .. p];
        }
    }

    // ------------------------------------------------------------------ refinement

    /// <summary>
    ///     Every marked-edge pattern of a single tetrahedron, under every relabelling of
    ///     its vertices and with shuffled midpoint labels: all branches of the subdivision
    ///     templates, including every global-label diagonal choice.
    /// </summary>
    [Fact]
    public void RefineSingleTetAllPatternsMatchesGolden()
    {
        var h = new Fnv();
        int[][] local = [[0, 1], [1, 2], [0, 2], [2, 3], [3, 0], [1, 3]];
        for (var mask = 0; mask < 64; mask++)
        for (var ip = 0; ip < Perms4.Length; ip++)
        {
            var p = Perms4[ip];
            var mesh = new SimplexMesh();
            for (var v = 0; v < 4; v++) mesh.AddNode(v);
            mesh.AddTetrahedron(p[0], p[1], p[2], p[3]);
            var marked = new List<(int, int)>();
            for (var e = 0; e < 6; e++)
                if ((mask & (1 << e)) != 0)
                    marked.Add((p[local[e][0]], p[local[e][1]]));
            var rnd = new Random(31 * mask + ip);
            for (var i = marked.Count - 1; i > 0; i--)
            {
                var j = rnd.Next(i + 1);
                (marked[i], marked[j]) = (marked[j], marked[i]);
            }

            var (refined, remap) = MeshRefinement.Refine(mesh, marked);
            AddMesh(h, refined, withLineage: true);
            foreach (var r in remap) h.Add(r);
        }

        Assert.Equal("e2cb564b71cfde86", h.Hex);
    }

    /// <summary>The sizing count of the templates equals the number of children they emit.</summary>
    [Fact]
    public void TetChildrenCountsMatchTheTemplates()
    {
        var rows = new List<int[]>();
        Span<int> mid = stackalloc int[6];
        for (var mask = 0; mask < 64; mask++)
        {
            var nmarked = 0;
            for (var e = 0; e < 6; e++)
            {
                mid[e] = (mask & (1 << e)) != 0 ? 10 + e : 0;
                if (mid[e] > 0) nmarked++;
            }

            rows.Clear();
            SubdivisionTemplates.SplitTet([0, 1, 2, 3], mid, rows);
            Assert.Equal(SubdivisionTemplates.TetChildren(nmarked, mid), rows.Count);
        }
    }

    /// <summary>Every marked-edge pattern of a triangle under every vertex relabelling.</summary>
    [Fact]
    public void RefineSingleTriangleAllPatternsMatchesGolden()
    {
        var h = new Fnv();
        int[][] perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
        for (var mask = 0; mask < 8; mask++)
        for (var ip = 0; ip < perms.Length; ip++)
        for (var order = 0; order < 2; order++)
        {
            var p = perms[ip];
            var mesh = new SimplexMesh();
            for (var v = 0; v < 3; v++) mesh.AddNode(v);
            mesh.AddTriangle(p[0], p[1], p[2]);
            mesh.AddBar(p[0], p[1]);
            mesh.AddPoint(p[1]);
            var marked = new List<(int, int)>();
            if ((mask & 1) != 0) marked.Add((p[0], p[1]));
            if ((mask & 2) != 0) marked.Add((p[1], p[2]));
            if ((mask & 4) != 0) marked.Add((p[2], p[0]));
            if (order == 1) marked.Reverse();
            var (refined, remap) = MeshRefinement.Refine(mesh, marked);
            AddMesh(h, refined, withLineage: true);
            foreach (var r in remap) h.Add(r);
        }

        Assert.Equal("6afe93a9319edf65", h.Hex);
    }

    /// <summary>Random markings of a structured tetrahedral box, with and without closure.</summary>
    [Fact]
    public void RefineBoxMeshMatchesGolden()
    {
        var h = new Fnv();
        foreach (var (p, seed, closure) in new[] { (0.15, 1, false), (0.4, 2, false), (0.75, 3, false), (1.0, 4, false), (0.1, 5, true) })
        {
            var (mesh, x) = SimplexRemesher.CreateBoxMesh(3, 3, 3);
            var marked = Pick(UniqueEdges<Tet4>(mesh), p, seed);
            var (refined, remap) = MeshRefinement.Refine(mesh, marked, enforceClosureForTets: closure);
            AddMesh(h, refined, withLineage: true);
            foreach (var r in remap) h.Add(r);
            AddCoords(h, MeshRefinement.InterpolateCoordinates(refined, x));
        }

        Assert.Equal("9b75039591a0a72c", h.Hex);
    }

    /// <summary>Random markings of a structured triangulation carrying bars and points.</summary>
    [Fact]
    public void RefineRectangleWithBarsAndPointsMatchesGolden()
    {
        var h = new Fnv();
        foreach (var (p, seed) in new[] { (0.2, 11), (0.5, 12), (1.0, 13) })
        {
            var (mesh, x) = SimplexRemesher.CreateRectangularMesh(5, 4, 0, 1, 0, 1);
            for (var i = 0; i < 5; i++) mesh.AddBar(i, i + 1);
            mesh.AddPoint(0);
            mesh.AddPoint(5);
            var marked = Pick(UniqueEdges<Tri3>(mesh), p, seed);
            var (refined, remap) = MeshRefinement.Refine(mesh, marked);
            AddMesh(h, refined, withLineage: true);
            foreach (var r in remap) h.Add(r);
            AddCoords(h, MeshRefinement.InterpolateCoordinates(refined, x));
        }

        Assert.Equal("a8d65776cff60746", h.Hex);
    }

    // ------------------------------------------------------------------ remesher helpers

    [Fact]
    public void BoundaryNodesMatchGolden()
    {
        var h = new Fnv();
        var (tri, _) = SimplexRemesher.CreateRectangularMesh(6, 5, 0, 1, 0, 1);
        AddSet(h, SimplexRemesher.FindBoundaryNodes(tri));
        AddSet(h, SimplexRemesher.IdentifyBoundaryNodes(tri));
        var (tet, _) = SimplexRemesher.CreateBoxMesh(3, 4, 2);
        AddSet(h, SimplexRemesher.FindBoundaryNodes3D(tet));
        AddSet(h, SimplexRemesher.IdentifyBoundaryNodes(tet));
        Assert.Equal("72d25d43333a7665", h.Hex);
    }

    [Fact]
    public void SmoothingMatchesGolden()
    {
        var h = new Fnv();
        var (tri, x2) = SimplexRemesher.CreateRectangularMesh(7, 6, 0, 1, 0, 1);
        Perturb(x2, 2, 0.04, 21);
        AddCoords(h, SimplexRemesher.LaplacianSmoothing(tri, x2, iterations: 3, relaxation: 0.6));
        AddCoords(h, SimplexRemesher.LaplacianSmoothing(tri, x2, iterations: 2, fixedNodes: [8, 9],
            keepBoundaryFixed: false));
        AddCoords(h, SimplexRemesher.CVTSmoothing(tri, x2, iterations: 3, relaxation: 0.7));
        var (tet, x3) = SimplexRemesher.CreateBoxMesh(3, 3, 3);
        Perturb(x3, 3, 0.03, 22);
        AddCoords(h, SimplexRemesher.LaplacianSmoothing(tet, x3, iterations: 3, relaxation: 0.5));
        Assert.Equal("c70f65b654649ce4", h.Hex);
    }

    /// <summary>
    ///     Two successive 2D crack insertions with smoothing: the second reads the edges the
    ///     first left discovered on its result, so the state of the returned mesh is pinned too.
    /// </summary>
    [Fact]
    public void CrackInsertion2DMatchesGolden()
    {
        var h = new Fnv();
        var (mesh, x) = SimplexRemesher.CreateRectangularMesh(12, 12, 0, 1, 0, 1);
        Perturb(x, 2, 0.012, 31);
        var (c1, x1) = SimplexRemesher.CreateCrackFromSignedField(mesh, x,
            (px, py, _) => py - 0.513, (px, _, _) => px < 0.55 ? -1.0 : px - 0.55, enableSmoothing: true);
        Assert.True(c1.Count<Tri3>() > mesh.Count<Tri3>() && c1.Count<Node>() > mesh.Count<Node>());
        AddMesh(h, c1, withLineage: false);
        AddCoords(h, x1);
        h.Add(c1.Count<Edge>());
        var (c2, x2) = SimplexRemesher.CreateCrackFromSignedField(c1, x1,
            (_, py, _) => py - 0.257, (px, _, _) => px < 0.3 ? 0.3 - px : px > 0.7 ? px - 0.7 : -1.0,
            enableSmoothing: true);
        Assert.True(c2.Count<Node>() > c1.Count<Node>());
        AddMesh(h, c2, withLineage: false);
        AddCoords(h, x2);
        Assert.Equal("bf7b61985a23c655", h.Hex);
    }

    [Fact]
    public void CrackInsertion3DMatchesGolden()
    {
        var h = new Fnv();
        var (mesh, x) = SimplexRemesher.CreateBoxMesh(5, 5, 5);
        Perturb(x, 3, 0.01, 41);
        var (c, xc) = SimplexRemesher.CreateCrackFromSignedField3D(mesh, x,
            (_, _, pz) => pz - 0.511, (px, _, _) => px < 0.6 ? -1.0 : px - 0.6, enableSmoothing: true);
        Assert.True(c.Count<Tet4>() > mesh.Count<Tet4>() && c.Count<Node>() > mesh.Count<Node>());
        AddMesh(h, c, withLineage: false);
        AddCoords(h, xc);
        Assert.Equal("351cc9f42a3c5950", h.Hex);
    }
}

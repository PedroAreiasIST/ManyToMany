// IncrementalMesh.cs — delta maintenance of the MM3 derived connectivity structures.
//
// State (all persistent, see PArray):
//   A      element -> ordered node row (empty row = deleted element)
//   T      node -> (element, local position) pairs sorted by element  [transpose + node positions]
//   EP     (e, p) -> index of e in T[A[e][p]]                         [element positions]
//   P      gapped CSR over nodes: per row start/len/cap, per slot column and pair count
//   Dest   (e, i, j) -> slot of A[e][j] in pattern row A[e][i]        [assembly destinations]
//   FA     (e, f) -> smallest other element containing facet f, or -1 [facet adjacency]
//   B      sparse set of boundary facets (FA = -1)
//   Color  greedy colouring of the vertex-adjacency graph
//
// An edit replaces the row of one element (creation: empty old row; deletion: empty
// new row). Apply(e, O, N) updates every structure from O and N alone:
//   - transpose: e leaves the rows of O \ N and enters the rows of N \ O at its sorted
//     place; nodes in O n N only update the stored local position;
//   - element positions: entries behind an insertion or removal shift by one, and
//     are found through the node positions stored in the same transpose row;
//   - pattern: every pair of O is decremented through Dest (one access per pair),
//     every pair of N is incremented at its slot; a pair keeps its slot while its
//     count is positive, so destinations of untouched elements never change;
//   - facets: entries of elements sharing an old or a new facet of e are recomputed
//     by intersecting the transpose rows of the facet nodes;
//   - colour: e takes the smallest colour unused by its current vertex neighbours.

using System.Runtime.CompilerServices;

namespace Step2;

public sealed class IncrementalMesh
{
    private const int RowShift = 8;   // 256 row references per chunk
    private const int IntShift = 10;  // 1024 ints per chunk
    private const int IntChunk = 1 << IntShift;

    public readonly ElementSpec Spec;
    private readonly int K, KK, Nf;
    private readonly int[][] _facetLocal;
    private readonly int[][][] _facetByLen;
    private readonly bool _mixed;
    private readonly Epoch _epoch = new();
    private readonly bool _persistent;
    private readonly int _slackMin;
    private readonly double _slackFactor;

    private readonly PArray<int[]> _A;
    private readonly PArray<int[]> _T;
    private readonly PArray<int> _EP;
    private readonly PArray<int> _Dest;
    private readonly PArray<int> _FA;
    private readonly PArray<int> _Color;
    private readonly PArray<int> _RowStart;
    private readonly PArray<int> _RowLen;
    private readonly PArray<int> _RowCap;
    private readonly PArray<int> _Col;
    private readonly PArray<int> _Cnt;
    private readonly PArray<int> _BList;
    private readonly PArray<int> _BPos;

    private int _nE, _nN, _total, _bCount, _colorCount;
    private long _garbage;
    private Snapshot? _current;

    private readonly List<IAttribute> _attributes = new();
    // per-thread scratch (several writer threads during a concurrent cut)
    private sealed class Scratch
    {
        public readonly List<(int G, int F)> Pending = new();
        public readonly List<int> Sharers = new();
        public readonly List<int> Mixed = new();
        public int[] ColourUsed = new int[128];
        public int ColourStamp;
    }

    [ThreadStatic] private static Scratch? t_scratch;
    private static Scratch S => t_scratch ??= new Scratch();
    private readonly object _segmentLock = new();

    // concurrent mode: counters and boundary-set changes are kept per writer thread and
    // merged by EndConcurrent (a facet is touched by one writer only, so replaying each
    // thread's log in order gives the set of a sequential application)
    private const int CEdits = 0, CShifts = 1, CReuses = 2, CGarbage = 3, CRelocations = 4, CFacets = 5;
    private readonly ThreadLocal<long[]> _localCounts = new(() => new long[6], trackAllValues: true);
    private readonly ThreadLocal<List<int>> _boundaryLog = new(() => new List<int>(), trackAllValues: true);

    private bool Concurrent => _epoch.Concurrent;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void Count(int counter, long d = 1)
    {
        if (Concurrent)
        {
            _localCounts.Value![counter] += d;
            return;
        }

        switch (counter)
        {
            case CEdits: Edits += d; break;
            case CShifts: PositionShifts += d; break;
            case CReuses: SlotReuses += d; break;
            case CGarbage: _garbage += d; break;
            case CRelocations: Relocations += d; break;
            case CFacets: FacetRecomputes += d; break;
        }
    }

    // instrumentation
    public long Edits, Relocations, SlotReuses, Compactions, FacetRecomputes, PositionShifts;

    /// <summary>Phase timers (Stopwatch ticks), active when STEP2_PROFILE=1.</summary>
    public static readonly bool Profile = Environment.GetEnvironmentVariable("STEP2_PROFILE") == "1";
    public static readonly long[] PhaseTicks = new long[8];
    public static readonly string[] PhaseNames =
        ["facets (old)", "pattern (old pairs)", "transpose + positions", "pattern (new pairs) + destinations",
         "facets (new)", "colour", "facets (neighbours)", "row, checks"];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static long Tick() => Profile ? System.Diagnostics.Stopwatch.GetTimestamp() : 0;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Lap(int phase, ref long t)
    {
        if (!Profile) return;
        var now = System.Diagnostics.Stopwatch.GetTimestamp();
        PhaseTicks[phase] += now - t;
        t = now;
    }

    public IncrementalMesh(FlatStructures s, ElementSpec spec, bool persistent,
        int slackMin = 4, double slackFactor = 0.25)
    {
        Spec = spec;
        K = spec.K;
        KK = K * K;
        Nf = spec.Nf;
        _facetLocal = spec.FacetLocal;
        _facetByLen = spec.FacetLocalByLength;
        _mixed = spec.Mixed;
        _persistent = persistent;
        _slackMin = slackMin;
        _slackFactor = slackFactor;
        _nE = s.NE;
        _nN = s.NN;

        _A = new PArray<int[]>(_epoch, RowShift, Array.Empty<int>(), _nE);
        for (var e = 0; e < _nE; e++) _A.Set(e, (int[])s.A[e].Clone());

        _T = new PArray<int[]>(_epoch, RowShift, Array.Empty<int>(), _nN);
        for (var v = 0; v < _nN; v++)
        {
            var b = 2 * s.TPtr[v];
            var n = 2 * (s.TPtr[v + 1] - s.TPtr[v]);
            if (n > 0) _T.Set(v, s.TPairs.AsSpan(b, n).ToArray());
        }

        _EP = PArray<int>.FromFlat(_epoch, IntShift, -1, s.EP, _nE * K);
        _Dest = PArray<int>.FromFlat(_epoch, IntShift, -1, s.Dest, _nE * KK);
        _FA = PArray<int>.FromFlat(_epoch, IntShift, -1, s.FA, _nE * Nf);
        _Color = PArray<int>.FromFlat(_epoch, IntShift, -1, s.Color, _nE);
        _colorCount = s.ColorCount;

        // pattern rows with slack, no row crossing a chunk boundary
        var rowStart = new int[_nN];
        var rowLen = new int[_nN];
        var rowCap = new int[_nN];
        var total = 0;
        long garbage = 0;
        for (var v = 0; v < _nN; v++)
        {
            var len = s.PPtr[v + 1] - s.PPtr[v];
            rowLen[v] = len;
            if (len == 0)
            {
                rowStart[v] = -1;
                continue;
            }

            var cap = CapacityFor(len);
            var within = total & (IntChunk - 1);
            if (within + cap > IntChunk)
            {
                garbage += IntChunk - within;
                total += IntChunk - within;
            }

            rowStart[v] = total;
            rowCap[v] = cap;
            total += cap;
        }

        var col = new int[total];
        var cnt = new int[total];
        for (var v = 0; v < _nN; v++)
        {
            if (rowLen[v] == 0) continue;
            Array.Copy(s.PCol, s.PPtr[v], col, rowStart[v], rowLen[v]);
            Array.Copy(s.PCnt, s.PPtr[v], cnt, rowStart[v], rowLen[v]);
        }

        _total = total;
        _garbage = garbage;
        _RowStart = PArray<int>.FromFlat(_epoch, IntShift, -1, rowStart, _nN);
        _RowLen = PArray<int>.FromFlat(_epoch, IntShift, 0, rowLen, _nN);
        _RowCap = PArray<int>.FromFlat(_epoch, IntShift, 0, rowCap, _nN);
        _Col = PArray<int>.FromFlat(_epoch, IntShift, 0, col, total);
        _Cnt = PArray<int>.FromFlat(_epoch, IntShift, 0, cnt, total);

        // boundary sparse set
        var bpos = new int[_nE * Nf];
        Array.Fill(bpos, -1);
        for (var i = 0; i < s.Boundary.Length; i++) bpos[s.Boundary[i]] = i;
        _bCount = s.Boundary.Length;
        _BList = PArray<int>.FromFlat(_epoch, IntShift, -1, s.Boundary, Math.Max(_bCount, 16));
        _BPos = PArray<int>.FromFlat(_epoch, IntShift, -1, bpos, _nE * Nf);

        Commit();
    }

    private int CapacityFor(int len) =>
        Math.Min(IntChunk, len + Math.Max(_slackMin, (int)Math.Ceiling(len * _slackFactor)));

    // ================================================================ reads (working version)

    public int ElementCount => _nE;
    public int NodeCount => _nN;
    public bool Persistent => _persistent;
    public int CurrentEpoch => _epoch.Value;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<int> Row(int e) => _A[e];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<int> TransposePairs(int v) => _T[v];

    public int ElementPosition(int e, int p) => _EP[e * K + p];
    public int Destination(int e, int i, int j) => _Dest[e * KK + i * K + j];
    public int FacetNeighbour(int e, int f) => _FA[e * Nf + f];
    public int Color(int e) => _Color[e];
    public int ColorCount => _colorCount;
    public int BoundaryCount => _bCount;
    public int TotalSlots => _total;
    public long GarbageSlots => _garbage;

    public (int Start, int Length) PatternRow(int v) => (_RowStart[v], _RowLen[v]);
    public ReadOnlySpan<int> PatternColumns(int v) =>
        _RowStart[v] < 0 ? ReadOnlySpan<int>.Empty : _Col.Segment(_RowStart[v], _RowLen[v]);
    public ReadOnlySpan<int> PatternCounts(int v) =>
        _RowStart[v] < 0 ? ReadOnlySpan<int>.Empty : _Cnt.Segment(_RowStart[v], _RowLen[v]);

    /// <summary>Last published version (lock-free read).</summary>
    public Snapshot Current => Volatile.Read(ref _current)!;

    /// <summary>
    ///     Registers a per-element or per-node array versioned with the mesh. Register before
    ///     publishing the versions it must belong to: a rollback to a version published before
    ///     the registration leaves the attribute unchanged.
    /// </summary>
    public Attribute<T> AddAttribute<T>(AttributeKind kind, T fill)
    {
        var length = kind == AttributeKind.Element ? _nE : _nN;
        var shift = typeof(T).IsValueType ? IntShift : RowShift;
        var a = new Attribute<T>(new PArray<T>(_epoch, shift, fill, length), kind, _attributes.Count);
        _attributes.Add(a);
        return a;
    }

    // ================================================================ edits

    public int AddNode()
    {
        var v = _nN++;
        _T.Resize(_nN);
        _RowStart.Resize(_nN);
        _RowLen.Resize(_nN);
        _RowCap.Resize(_nN);
        foreach (var a in _attributes)
            if (a.Kind == AttributeKind.Node) a.Resize(_nN);
        return v;
    }

    public int AddElement(ReadOnlySpan<int> nodes)
    {
        var e = _nE++;
        _A.Resize(_nE);
        _EP.Resize(_nE * K);
        _Dest.Resize(_nE * KK);
        _FA.Resize(_nE * Nf);
        _Color.Resize(_nE);
        _BPos.Resize(_nE * Nf);
        foreach (var a in _attributes)
            if (a.Kind == AttributeKind.Element) a.Resize(_nE);
        Apply(e, Array.Empty<int>(), nodes.ToArray());
        return e;
    }

    public void ReplaceElement(int e, ReadOnlySpan<int> nodes) => Apply(e, _A[e], nodes.ToArray());

    public void RemoveElement(int e)
    {
        if (_A[e].Length == 0) return;
        Apply(e, _A[e], Array.Empty<int>());
    }

    /// <summary>Appends <paramref name="count"/> deleted (empty) elements and returns the first index.</summary>
    public int ReserveElements(int count)
    {
        var first = _nE;
        _nE += count;
        _A.Resize(_nE);
        _EP.Resize(_nE * K);
        _Dest.Resize(_nE * KK);
        _FA.Resize(_nE * Nf);
        _Color.Resize(_nE);
        _BPos.Resize(_nE * Nf);
        foreach (var a in _attributes)
            if (a.Kind == AttributeKind.Element) a.Resize(_nE);
        return first;
    }

    /// <summary>Gives a row to an empty element (reserved or deleted); same as ReplaceElement.</summary>
    public void CreateElement(int e, ReadOnlySpan<int> nodes)
    {
        if (_A[e].Length != 0) throw new InvalidOperationException($"element {e} is not empty");
        Apply(e, Array.Empty<int>(), nodes.ToArray());
    }

    /// <summary>
    ///     Enters concurrent mode: makes every directory private to the working version and
    ///     reserves pattern slots, since no array may grow while several threads write. Edits
    ///     of elements whose nodes no other thread touches may then run in parallel (see
    ///     Cutter.CutParallel).
    /// </summary>
    public void BeginConcurrent(int extraSlots)
    {
        if (_total + extraSlots > _Col.Length)
        {
            _Col.Resize(_total + extraSlots);
            _Cnt.Resize(_total + extraSlots);
        }

        _A.PrepareWrite(); _T.PrepareWrite(); _EP.PrepareWrite(); _Dest.PrepareWrite(); _FA.PrepareWrite();
        _Color.PrepareWrite(); _RowStart.PrepareWrite(); _RowLen.PrepareWrite(); _RowCap.PrepareWrite();
        _Col.PrepareWrite(); _Cnt.PrepareWrite(); _BList.PrepareWrite(); _BPos.PrepareWrite();
        foreach (var a in _attributes) a.PrepareWrite();
        _epoch.Concurrent = true;
    }

    /// <summary>Leaves concurrent mode: merges the per-thread counters and boundary changes.</summary>
    public void EndConcurrent()
    {
        _epoch.Concurrent = false;
        foreach (var log in _boundaryLog.Values)
        {
            foreach (var x in log)
                if (x >= 0) BoundaryAddCore(x);
                else BoundaryRemoveCore(~x);
            log.Clear();
        }

        foreach (var c in _localCounts.Values)
        {
            Edits += c[CEdits];
            PositionShifts += c[CShifts];
            SlotReuses += c[CReuses];
            _garbage += c[CGarbage];
            Relocations += c[CRelocations];
            FacetRecomputes += c[CFacets];
            Array.Clear(c);
        }
    }

    /// <summary>
    ///     Upper bound of the pattern slots (row segments plus chunk padding) that growing the
    ///     row of node v to <paramref name="finalLength"/> entries can allocate.
    /// </summary>
    public int SlotsToGrow(int v, int finalLength)
    {
        int cap, total = 0;
        if (_RowStart[v] < 0)
        {
            if (finalLength == 0) return 0;
            cap = CapacityFor(8);
            total += 2 * cap;
        }
        else cap = _RowCap[v];

        while (cap < finalLength && cap < IntChunk)
        {
            cap = Math.Min(IntChunk, Math.Max(8, 2 * cap));
            total += 2 * cap;
        }

        return total;
    }

    private void Validate(int[] row)
    {
        if (row.Length == 0) return;
        if (_mixed ? row.Length > K : row.Length != K)
            throw new ArgumentException(_mixed ? $"row must have 1 to {K} nodes" : $"row must have {K} nodes");
        for (var p = 0; p < row.Length; p++)
        {
            if ((uint)row[p] >= (uint)_nN) throw new ArgumentOutOfRangeException(nameof(row), $"node {row[p]}");
            for (var q = 0; q < p; q++)
                if (row[q] == row[p])
                    throw new ArgumentException($"node {row[p]} repeated in row");
        }
    }

    /// <summary>
    ///     True when facet f has the same node set before and after the edit; its adjacency
    ///     entry, its boundary status and the entries of its sharers are then unchanged
    ///     (e stays in the transpose rows of all its nodes), so the facet is skipped.
    /// </summary>
    private bool FacetUnchanged(int[] O, int[] N, int f)
    {
        if (O.Length == 0 || O.Length != N.Length) return false;
        var loc = _facetByLen[N.Length][f];
        foreach (var l in loc)
            if (Array.IndexOf(O, N[l], 0, O.Length) < 0 || !OnFacet(O, loc, N[l]))
                return false;
        return true;
    }

    private static bool OnFacet(int[] row, int[] loc, int v)
    {
        foreach (var l in loc)
            if (row[l] == v)
                return true;
        return false;
    }

    private void Apply(int e, int[] O, int[] N)
    {
        var tk = Tick();
        Validate(N);
        Count(CEdits);
        var pending = S.Pending;
        pending.Clear();
        int lo = O.Length, ln = N.Length;
        int nfO = Spec.FacetCount(lo), nfN = Spec.FacetCount(ln);
        var sameLen = lo == ln;
        Span<bool> same = stackalloc bool[Nf];
        for (var f = 0; f < nfN; f++) same[f] = sameLen && FacetUnchanged(O, N, f);
        // keep[p]: node at local position p is the same before and after the edit
        Span<bool> keep = stackalloc bool[K];
        for (var p = 0; p < K; p++) keep[p] = sameLen && p < lo && O[p] == N[p];
        Lap(7, ref tk);

        // ------------------------------------------------ removal of the old incidences
        if (O.Length > 0)
        {
            for (var f = 0; f < nfO; f++)
            {
                if (same[f]) continue;
                var id = e * Nf + f;
                if (_FA[id] < 0) BoundaryRemove(id);
                CollectSharers(O, f, e, pending);
            }

            Lap(0, ref tk);

            for (var i = 0; i < lo; i++)
            {
                var rs = _RowStart[O[i]];
                var db = e * KK + i * K;
                for (var j = 0; j < lo; j++)
                {
                    if (keep[i] && keep[j]) continue; // pair, count and slot unchanged
                    var at = rs + _Dest[db + j];
                    _Cnt.Set(at, _Cnt[at] - 1);
                }
            }

            Lap(1, ref tk);
            for (var p = 0; p < O.Length; p++)
                if (Array.IndexOf(N, O[p]) < 0)
                    TransposeRemove(O[p], e);
            Lap(2, ref tk);
        }

        _A.Set(e, N);
        Lap(7, ref tk);

        // ------------------------------------------------ insertion of the new incidences
        if (N.Length > 0)
        {
            for (var p = 0; p < ln; p++)
            {
                if (keep[p]) continue; // e stays at the same index of T[v] with the same position
                var v = N[p];
                if (Array.IndexOf(O, v) >= 0) TransposeRelocate(v, e, p);
                else TransposeInsert(v, e, p);
            }

            Lap(2, ref tk);

            for (var i = 0; i < ln; i++)
            {
                var db = e * KK + i * K;
                for (var j = 0; j < ln; j++)
                {
                    if (keep[i] && keep[j]) continue;
                    var slot = AcquireSlot(N[i], N[j]);
                    var at = _RowStart[N[i]] + slot;
                    _Cnt.Set(at, _Cnt[at] + 1);
                    _Dest.Set(db + j, slot);
                }
            }

            Lap(3, ref tk);

            if (lo > ln) ClearTail(e, ln, nfN);

            for (var f = 0; f < nfN; f++)
            {
                if (same[f]) continue;
                var id = e * Nf + f;
                var nb = MinOtherSharer(N, f, e);
                _FA.Set(id, nb);
                if (nb < 0) BoundaryAdd(id);
                CollectSharers(N, f, e, pending);
            }

            Lap(4, ref tk);

            // Neighbours through nodes kept from O were neighbours before and already differ;
            // only neighbours through new nodes can conflict with the current colour.
            var cur = O.Length > 0 ? _Color[e] : -1;
            if (cur < 0 || ConflictsThroughNewNodes(e, O, N, cur))
            {
                var c = SmallestFreeColour(e, N);
                _Color.Set(e, c);
                if (c + 1 > _colorCount)
                {
                    if (!Concurrent) _colorCount = c + 1;
                    else
                    {
                        int cur2;
                        while ((cur2 = Volatile.Read(ref _colorCount)) < c + 1 &&
                               Interlocked.CompareExchange(ref _colorCount, c + 1, cur2) != cur2)
                        {
                        }
                    }
                }
            }

            Lap(5, ref tk);
        }
        else
        {
            for (var p = 0; p < K; p++) _EP.Set(e * K + p, -1);
            for (var t = 0; t < KK; t++) _Dest.Set(e * KK + t, -1);
            for (var f = 0; f < Nf; f++) _FA.Set(e * Nf + f, -1);
            _Color.Set(e, -1);
        }

        Lap(7, ref tk);
        foreach (var (g, fg) in pending) RecomputeFacet(g, fg);
        Lap(6, ref tk);
    }

    /// <summary>Clears the table entries beyond a shorter new row (mixed meshes, type change).</summary>
    private void ClearTail(int e, int ln, int nfN)
    {
        for (var p = ln; p < K; p++) _EP.Set(e * K + p, -1);
        for (var i = 0; i < K; i++)
        for (var j = 0; j < K; j++)
            if (i >= ln || j >= ln)
                _Dest.Set(e * KK + i * K + j, -1);
        for (var f = nfN; f < Nf; f++) _FA.Set(e * Nf + f, -1);
    }

    // ------------------------------------------------ transpose and positions

    private void TransposeInsert(int v, int e, int p)
    {
        var row = _T[v];
        var idx = PairRow.IndexOf(row, e);
        if (idx >= 0) throw new InvalidOperationException($"element {e} already in row {v}");
        var ins = ~idx;
        var n = row.Length >> 1;
        var nr = new int[row.Length + 2];
        Array.Copy(row, 0, nr, 0, 2 * ins);
        nr[2 * ins] = e;
        nr[2 * ins + 1] = p;
        Array.Copy(row, 2 * ins, nr, 2 * ins + 2, row.Length - 2 * ins);
        _T.Set(v, nr);
        _EP.Set(e * K + p, ins);
        for (var q = ins + 1; q <= n; q++) _EP.Set(nr[2 * q] * K + nr[2 * q + 1], q);
        Count(CShifts, n - ins);
    }

    private void TransposeRemove(int v, int e)
    {
        var row = _T[v];
        var idx = PairRow.IndexOf(row, e);
        if (idx < 0) throw new InvalidOperationException($"element {e} missing from row {v}");
        var n = row.Length >> 1;
        var nr = new int[row.Length - 2];
        Array.Copy(row, 0, nr, 0, 2 * idx);
        Array.Copy(row, 2 * idx + 2, nr, 2 * idx, row.Length - 2 * idx - 2);
        _T.Set(v, nr);
        for (var q = idx; q < n - 1; q++) _EP.Set(nr[2 * q] * K + nr[2 * q + 1], q);
        Count(CShifts, n - 1 - idx);
    }

    private void TransposeRelocate(int v, int e, int p)
    {
        var row = _T[v];
        var idx = PairRow.IndexOf(row, e);
        if (idx < 0) throw new InvalidOperationException($"element {e} missing from row {v}");
        if (row[2 * idx + 1] != p)
        {
            var nr = (int[])row.Clone();
            nr[2 * idx + 1] = p;
            _T.Set(v, nr);
        }

        _EP.Set(e * K + p, idx);
    }

    // ------------------------------------------------ pattern

    private int AcquireSlot(int v, int m)
    {
        var rs = _RowStart[v];
        if (rs < 0)
        {
            AllocateRow(v, CapacityFor(8));
            rs = _RowStart[v];
        }

        var len = _RowLen[v];
        var idx = _Col.Segment(rs, len).IndexOf(m);
        if (idx >= 0) return idx;

        var t = _Cnt.Segment(rs, len).IndexOf(0);
        if (t >= 0)
        {
            _Col.Set(rs + t, m);
            Count(CReuses);
            return t;
        }

        if (len == _RowCap[v])
        {
            Relocate(v);
            rs = _RowStart[v];
        }

        _Col.Set(rs + len, m);
        _Cnt.Set(rs + len, 0);
        _RowLen.Set(v, len + 1);
        return len;
    }

    private int AllocateSegment(int cap)
    {
        if (!Concurrent) return AllocateSegmentCore(cap);
        lock (_segmentLock)
        {
            var start = AllocateSegmentCore(cap);
            if (_total > _Col.Length) throw new InvalidOperationException("pattern reserve exhausted");
            return start;
        }
    }

    private int AllocateSegmentCore(int cap)
    {
        var start = _total;
        var within = start & (IntChunk - 1);
        if (within + cap > IntChunk)
        {
            Count(CGarbage, IntChunk - within);
            start += IntChunk - within;
        }

        _total = start + cap;
        if (_total > _Col.Length && !Concurrent)
        {
            var grow = Math.Max(_total, _Col.Length + (_Col.Length >> 3));
            _Col.Resize(grow);
            _Cnt.Resize(grow);
        }

        return start;
    }

    private void AllocateRow(int v, int cap)
    {
        var start = AllocateSegment(cap);
        _RowStart.Set(v, start);
        _RowLen.Set(v, 0);
        _RowCap.Set(v, cap);
    }

    private void Relocate(int v)
    {
        var rs = _RowStart[v];
        var len = _RowLen[v];
        var cap = _RowCap[v];
        var ncap = Math.Min(IntChunk, Math.Max(8, 2 * cap));
        if (ncap <= len) throw new InvalidOperationException($"pattern row {v} exceeds {IntChunk} slots");
        var ns = AllocateSegment(ncap);
        var src = _Col.Segment(rs, len).ToArray();
        var srcCnt = _Cnt.Segment(rs, len).ToArray();
        src.CopyTo(_Col.WritableSegment(ns, len));
        srcCnt.CopyTo(_Cnt.WritableSegment(ns, len));
        Count(CGarbage, cap);
        _RowStart.Set(v, ns);
        _RowCap.Set(v, ncap);
        Count(CRelocations);
    }

    /// <summary>
    ///     Rewrites the pattern rows contiguously (slots unchanged, so destinations stay valid).
    ///     Global; intended to run rarely, when abandoned segments exceed half the slots.
    /// </summary>
    public void CompactPattern()
    {
        var rowStart = new int[_nN];
        var rowCap = new int[_nN];
        var total = 0;
        long garbage = 0;
        for (var v = 0; v < _nN; v++)
        {
            var len = _RowLen[v];
            if (_RowStart[v] < 0)
            {
                rowStart[v] = -1;
                continue;
            }

            var cap = Math.Max(_RowCap[v] > 2 * len ? CapacityFor(len) : _RowCap[v], len);
            var within = total & (IntChunk - 1);
            if (within + cap > IntChunk)
            {
                garbage += IntChunk - within;
                total += IntChunk - within;
            }

            rowStart[v] = total;
            rowCap[v] = cap;
            total += cap;
        }

        var col = new int[total];
        var cnt = new int[total];
        for (var v = 0; v < _nN; v++)
        {
            if (rowStart[v] < 0) continue;
            var len = _RowLen[v];
            _Col.Segment(_RowStart[v], len).CopyTo(col.AsSpan(rowStart[v], len));
            _Cnt.Segment(_RowStart[v], len).CopyTo(cnt.AsSpan(rowStart[v], len));
        }

        var newCol = PArray<int>.FromFlat(_epoch, IntShift, 0, col, total);
        var newCnt = PArray<int>.FromFlat(_epoch, IntShift, 0, cnt, total);
        _Col.Restore(newCol.Freeze());
        _Cnt.Restore(newCnt.Freeze());
        for (var v = 0; v < _nN; v++)
        {
            _RowStart.Set(v, rowStart[v]);
            _RowCap.Set(v, rowCap[v]);
        }

        _total = total;
        _garbage = garbage;
        Compactions++;
    }

    // ------------------------------------------------ facets

    /// <summary>Appends (g, facet of g) for every element g != self that contains facet f of row.</summary>

    private void CollectSharers(int[] row, int f, int self, List<(int, int)> pending)
    {
        var loc = _facetByLen[row.Length][f];
        if (_mixed)
        {
            var sharers = S.Sharers;
            sharers.Clear();
            CommonOthersAny(row, loc, self, sharers);
            foreach (var g in sharers)
            {
                var gr = _A[g];
                if (gr.Length == row.Length) pending.Add((g, OppositeLocal(gr, row, loc)));
            }

            return;
        }

        if (loc.Length == 3)
        {
            var sharers = S.Sharers;
            sharers.Clear();
            PairRow.CommonOthers(_T[row[loc[0]]], _T[row[loc[1]]], _T[row[loc[2]]], self, sharers);
            foreach (var g in sharers) pending.Add((g, OppositeLocal(_A[g], row, loc)));
            return;
        }

        var three = false;
        ReadOnlySpan<int> r0 = _T[row[loc[0]]];
        ReadOnlySpan<int> r1 = _T[row[loc[1]]];
        ReadOnlySpan<int> r2 = three ? _T[row[loc[2]]] : ReadOnlySpan<int>.Empty;
        if (r1.Length < r0.Length) { var t = r0; r0 = r1; r1 = t; }
        if (three && r2.Length < r0.Length) { var t = r0; r0 = r2; r2 = t; }
        for (var q = 0; q < r0.Length; q += 2)
        {
            var g = r0[q];
            if (g == self) continue;
            if (PairRow.IndexOf(r1, g) < 0) continue;
            if (three && PairRow.IndexOf(r2, g) < 0) continue;
            pending.Add((g, OppositeLocal(_A[g], row, loc)));
        }
    }

    /// <summary>Local index in <paramref name="g"/> of the node not on the facet (simplices).</summary>
    private static int OppositeLocal(int[] g, int[] row, int[] loc)
    {
        for (var p = 0; p < g.Length; p++)
        {
            var x = g[p];
            var onFacet = false;
            foreach (var l in loc)
                if (row[l] == x)
                {
                    onFacet = true;
                    break;
                }

            if (!onFacet) return p;
        }

        throw new InvalidOperationException("degenerate element: all nodes on the facet");
    }

    /// <summary>All elements other than self containing the facet nodes, ascending (any facet size).</summary>
    private void CommonOthersAny(int[] row, int[] loc, int self, List<int> result)
    {
        switch (loc.Length)
        {
            case 3:
                PairRow.CommonOthers(_T[row[loc[0]]], _T[row[loc[1]]], _T[row[loc[2]]], self, result);
                return;
            case 2:
            {
                ReadOnlySpan<int> r0 = _T[row[loc[0]]], r1 = _T[row[loc[1]]];
                if (r1.Length < r0.Length) { var t = r0; r0 = r1; r1 = t; }
                for (var q = 0; q < r0.Length; q += 2)
                    if (r0[q] != self && PairRow.IndexOf(r1, r0[q]) >= 0)
                        result.Add(r0[q]);
                return;
            }
            case 1:
            {
                ReadOnlySpan<int> r0 = _T[row[loc[0]]];
                for (var q = 0; q < r0.Length; q += 2)
                    if (r0[q] != self)
                        result.Add(r0[q]);
                return;
            }
        }
    }


    private int MinOtherSharer(int[] row, int f, int self)
    {
        if (_mixed)
        {
            var l = _facetByLen[row.Length][f];
            var mixed = S.Mixed;
            mixed.Clear();
            CommonOthersAny(row, l, self, mixed);
            foreach (var g in mixed)
                if (_A[g].Length == row.Length)
                    return g;
            return -1;
        }

        var loc = _facetLocal[f];
        return loc.Length == 3
            ? PairRow.MinOtherCommon(_T[row[loc[0]]], _T[row[loc[1]]], _T[row[loc[2]]], self)
            : PairRow.MinOtherCommon(_T[row[loc[0]]], _T[row[loc[1]]], self);
    }

    private void RecomputeFacet(int g, int fg)
    {
        var row = _A[g];
        if (row.Length == 0) return;
        Count(CFacets);
        var nb = MinOtherSharer(row, fg, g);
        var id = g * Nf + fg;
        var old = _FA[id];
        if (old == nb) return;
        _FA.Set(id, nb);
        if (old < 0) BoundaryRemove(id);
        if (nb < 0) BoundaryAdd(id);
    }

    private void BoundaryAdd(int id)
    {
        if (Concurrent) _boundaryLog.Value!.Add(id);
        else BoundaryAddCore(id);
    }

    private void BoundaryRemove(int id)
    {
        if (Concurrent) _boundaryLog.Value!.Add(~id);
        else BoundaryRemoveCore(id);
    }

    private void BoundaryAddCore(int id)
    {
        if (_BPos[id] >= 0) return;
        if (_bCount == _BList.Length) _BList.Resize(Math.Max(16, 2 * _BList.Length));
        _BList.Set(_bCount, id);
        _BPos.Set(id, _bCount);
        _bCount++;
    }

    private void BoundaryRemoveCore(int id)
    {
        var pos = _BPos[id];
        if (pos < 0) return;
        var last = --_bCount;
        var lastId = _BList[last];
        _BList.Set(pos, lastId);
        _BPos.Set(lastId, pos);
        _BPos.Set(id, -1);
    }

    // ------------------------------------------------ colour

    private bool ConflictsThroughNewNodes(int e, int[] O, int[] N, int colour)
    {
        foreach (var v in N)
        {
            if (Array.IndexOf(O, v) >= 0) continue;
            var tr = _T[v];
            for (var q = 0; q < tr.Length; q += 2)
                if (tr[q] != e && _Color[tr[q]] == colour)
                    return true;
        }

        return false;
    }

    private int SmallestFreeColour(int e, int[] row)
    {
        var sc = S;
        var mark = ++sc.ColourStamp;
        var used = sc.ColourUsed;
        foreach (var v in row)
        {
            var tr = _T[v];
            for (var q = 0; q < tr.Length; q += 2)
            {
                var g = tr[q];
                if (g == e) continue;
                var c = _Color[g];
                if (c < 0) continue;
                if (c >= used.Length)
                {
                    Array.Resize(ref sc.ColourUsed, Math.Max(c + 1, 2 * used.Length));
                    used = sc.ColourUsed;
                }

                used[c] = mark;
            }
        }

        var k = 0;
        while (k < used.Length && used[k] == mark) k++;
        if (k == used.Length) Array.Resize(ref sc.ColourUsed, 2 * used.Length);
        return k;
    }

    // ================================================================ versions

    /// <summary>
    ///     Publishes the working state as an immutable version (persistent mode) and returns it.
    ///     In transient mode the returned snapshot aliases the working arrays and is valid
    ///     only until the next edit.
    /// </summary>
    public Snapshot Commit()
    {
        var snap = new Snapshot(K, Nf, _nE, _nN, _total, _bCount, _colorCount, _garbage, _epoch.Value,
            _A.Freeze(), _T.Freeze(), _EP.Freeze(), _Dest.Freeze(), _FA.Freeze(), _Color.Freeze(),
            _RowStart.Freeze(), _RowLen.Freeze(), _RowCap.Freeze(), _Col.Freeze(), _Cnt.Freeze(),
            _BList.Freeze(), _BPos.Freeze(), _attributes.Select(a => a.Freeze()).ToArray());
        if (_persistent) _epoch.Value++;
        Volatile.Write(ref _current, snap);
        return snap;
    }

    /// <summary>Makes a published version the working state again; O(1) in the mesh size.</summary>
    public void Rollback(Snapshot s)
    {
        if (!_persistent) throw new InvalidOperationException("rollback requires persistent mode");
        _A.Restore(s.A);
        _T.Restore(s.T);
        _EP.Restore(s.EP);
        _Dest.Restore(s.Dest);
        _FA.Restore(s.FA);
        _Color.Restore(s.ColorView);
        _RowStart.Restore(s.RowStart);
        _RowLen.Restore(s.RowLen);
        _RowCap.Restore(s.RowCap);
        _Col.Restore(s.Col);
        _Cnt.Restore(s.Cnt);
        _BList.Restore(s.BList);
        _BPos.Restore(s.BPos);
        _nE = s.ElementCount;
        _nN = s.NodeCount;
        _total = s.TotalSlots;
        _bCount = s.BoundaryCount;
        _colorCount = s.ColorCount;
        _garbage = s.GarbageSlots;
        for (var i = 0; i < Math.Min(_attributes.Count, s.AttributeViews.Length); i++)
            _attributes[i].Restore(s.AttributeViews[i]);
        _epoch.Value++;
        Volatile.Write(ref _current, s);
    }

    /// <summary>Copy-on-write counters summed over all arrays.</summary>
    public (long Chunks, long Directories) CopyCounters()
    {
        long c = 0, d = 0;
        void Add<T>(PArray<T> a)
        {
            c += a.ChunkCopies;
            d += a.DirectoryCopies;
        }

        Add(_A); Add(_T); Add(_EP); Add(_Dest); Add(_FA); Add(_Color); Add(_RowStart); Add(_RowLen);
        Add(_RowCap); Add(_Col); Add(_Cnt); Add(_BList); Add(_BPos);
        return (c, d);
    }
}

// Snapshot.cs — an immutable published version of an IncrementalMesh.
//
// A snapshot holds the directory references of every array at publication time.
// In persistent mode the writer never writes to those directories or chunks again,
// so any number of threads read a snapshot without locks while the writer prepares
// the next version; the snapshot stays valid for as long as a reader holds it.

using System.Runtime.CompilerServices;

namespace Step2;

public sealed class Snapshot
{
    private readonly int K, KK, Nf;
    public readonly int ElementCount;
    public readonly int NodeCount;
    public readonly int TotalSlots;
    public readonly int BoundaryCount;
    public readonly int ColorCount;
    public readonly long GarbageSlots;
    public readonly int Epoch;

    internal readonly PView<int[]> A, T;
    internal readonly PView<int> EP, Dest, FA, ColorView, RowStart, RowLen, RowCap, Col, Cnt, BList, BPos;
    internal readonly object[] AttributeViews;

    internal Snapshot(int k, int nf, int nE, int nN, int total, int bCount, int colorCount, long garbage, int epoch,
        PView<int[]> a, PView<int[]> t, PView<int> ep, PView<int> dest, PView<int> fa, PView<int> color,
        PView<int> rowStart, PView<int> rowLen, PView<int> rowCap, PView<int> col, PView<int> cnt,
        PView<int> bList, PView<int> bPos, object[] attributeViews)
    {
        AttributeViews = attributeViews;
        K = k;
        KK = k * k;
        Nf = nf;
        ElementCount = nE;
        NodeCount = nN;
        TotalSlots = total;
        BoundaryCount = bCount;
        ColorCount = colorCount;
        GarbageSlots = garbage;
        Epoch = epoch;
        A = a;
        T = t;
        EP = ep;
        Dest = dest;
        FA = fa;
        ColorView = color;
        RowStart = rowStart;
        RowLen = rowLen;
        RowCap = rowCap;
        Col = col;
        Cnt = cnt;
        BList = bList;
        BPos = bPos;
    }

    public int NodesPerElement => K;
    public int FacetsPerElement => Nf;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<int> Row(int e) => A[e];

    /// <summary>(element, local position) pairs of node v, sorted by element.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<int> TransposePairs(int v) => T[v];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int ElementPosition(int e, int p) => EP[e * K + p];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Destination(int e, int i, int j) => Dest[e * KK + i * K + j];

    /// <summary>The K x K destination block of element e (row-relative slots).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<int> Destinations(int e) => Dest.Segment(e * KK, KK);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int FacetNeighbour(int e, int f) => FA[e * Nf + f];

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Color(int e) => ColorView[e];

    public int PatternRowStart(int v) => RowStart[v];
    public int PatternRowLength(int v) => RowLen[v];

    public ReadOnlySpan<int> PatternColumns(int v) =>
        RowStart[v] < 0 ? ReadOnlySpan<int>.Empty : Col.Segment(RowStart[v], RowLen[v]);

    public ReadOnlySpan<int> PatternCounts(int v) =>
        RowStart[v] < 0 ? ReadOnlySpan<int>.Empty : Cnt.Segment(RowStart[v], RowLen[v]);

    public int BoundaryFacet(int i) => BList[i];

    /// <summary>The attribute as published in this version.</summary>
    public PView<T> View<T>(Attribute<T> a) => (PView<T>)AttributeViews[a.Index];
}

// PArray.cs — persistent chunked arrays with epoch-tagged copy-on-write.
//
// A PArray<T> is a two-level structure: a directory of chunk references and
// fixed-size chunks of 2^shift entries. Every chunk, and the directory itself,
// carries an owner tag. The writer works in the current epoch: an object whose
// tag equals the epoch belongs to the working version and is written in place;
// any other tag means the object is shared with a published version and is
// copied before the first write (path copying, depth two). Publishing a version
// captures the directory references (a View) and advances the epoch, so the
// published arrays are never written again. Reclamation of versions nobody holds
// is left to the garbage collector.

using System.Runtime.CompilerServices;

namespace Step2;

/// <summary>
///     The writer's epoch, shared by all arrays of one mesh. Objects tagged with
///     <see cref="Value"/> are private to the working version.
/// </summary>
public sealed class Epoch
{
    public int Value = 1;

    /// <summary>
    ///     Set while several writer threads edit disjoint parts of the working version: chunk
    ///     copies are then claimed atomically, and directories must have been made writable
    ///     (<see cref="PArray{T}.PrepareWrite"/>) and sized beforehand.
    /// </summary>
    public volatile bool Concurrent;
}

/// <summary>Immutable view of a <see cref="PArray{T}"/> as published in a version.</summary>
public readonly struct PView<T>
{
    internal readonly T[][] Chunks;
    internal readonly int[] Owners;
    internal readonly int Shift;
    internal readonly int Mask;
    public readonly int Length;

    internal PView(T[][] chunks, int[] owners, int shift, int length)
    {
        Chunks = chunks;
        Owners = owners;
        Shift = shift;
        Mask = (1 << shift) - 1;
        Length = length;
    }

    public T this[int i]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => Chunks[i >> Shift][i & Mask];
    }

    /// <summary>Contiguous segment; the caller guarantees it does not cross a chunk boundary.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<T> Segment(int start, int length) =>
        new(Chunks[start >> Shift], start & Mask, length);
}

/// <summary>Persistent chunked array (writer side).</summary>
public sealed class PArray<T>
{
    private readonly Epoch _epoch;
    private readonly int _shift;
    private readonly int _mask;
    private readonly int _chunkSize;
    private readonly T _fill;

    private T[][] _chunks;
    private int[] _owners;
    private int _dirOwner;
    private int _chunkCount;
    private int _length;

    /// <summary>Chunks copied because they were shared with a published version.</summary>
    public long ChunkCopies;

    /// <summary>Directory copies (first write after a publication).</summary>
    public long DirectoryCopies;

    public PArray(Epoch epoch, int shift, T fill, int length = 0)
    {
        _epoch = epoch;
        _shift = shift;
        _chunkSize = 1 << shift;
        _mask = _chunkSize - 1;
        _fill = fill;
        _chunks = new T[Math.Max(4, ChunksFor(length))][];
        _owners = new int[_chunks.Length];
        _dirOwner = epoch.Value;
        Resize(length);
    }

    public int Length => _length;
    public int ChunkSize => _chunkSize;
    public int Shift => _shift;

    private int ChunksFor(int length) => (length + _chunkSize - 1) >> _shift;

    public T this[int i]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => _chunks[i >> _shift][i & _mask];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Set(int i, T value)
    {
        var c = i >> _shift;
        var e = _epoch.Value;
        // acquire read: a chunk copied by another writer thread (concurrent mode) is published
        // before its tag, so the chunk reference read below is the copy
        if (_dirOwner == e && Volatile.Read(ref _owners[c]) == e)
        {
            _chunks[c][i & _mask] = value;
            return;
        }

        WritableChunk(c)[i & _mask] = value;
    }

    /// <summary>Writable segment inside one chunk (copy-on-write applied).</summary>
    public Span<T> WritableSegment(int start, int length)
    {
        var c = start >> _shift;
        return new Span<T>(WritableChunk(c), start & _mask, length);
    }

    /// <summary>Read-only segment inside one chunk.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ReadOnlySpan<T> Segment(int start, int length) =>
        new(_chunks[start >> _shift], start & _mask, length);

    [MethodImpl(MethodImplOptions.NoInlining)]
    private T[] WritableChunk(int c)
    {
        var e = _epoch.Value;
        if (_epoch.Concurrent) return WritableChunkConcurrent(c, e);
        if (_dirOwner != e)
        {
            var chunks = new T[_chunks.Length][];
            Array.Copy(_chunks, chunks, _chunkCount);
            var owners = new int[_chunks.Length];
            Array.Copy(_owners, owners, _chunkCount);
            _chunks = chunks;
            _owners = owners;
            _dirOwner = e;
            DirectoryCopies++;
        }

        if (_owners[c] != e)
        {
            _chunks[c] = (T[])_chunks[c].Clone();
            _owners[c] = e;
            ChunkCopies++;
        }

        return _chunks[c];
    }

    private const int Claimed = int.MinValue;

    /// <summary>
    ///     Copy-on-write under concurrent writers: the first thread to claim the owner tag of a
    ///     shared chunk copies it and publishes the copy before the tag; the others wait for the
    ///     tag. Writes only go to chunks owned by the epoch, so no write is lost to a copy.
    /// </summary>
    private T[] WritableChunkConcurrent(int c, int e)
    {
        if (_dirOwner != e) throw new InvalidOperationException("directory not prepared for concurrent writes");
        var spin = new SpinWait();
        while (true)
        {
            var o = Volatile.Read(ref _owners[c]);
            if (o == e) return Volatile.Read(ref _chunks[c]);
            if (o != Claimed && Interlocked.CompareExchange(ref _owners[c], Claimed, o) == o)
            {
                var copy = (T[])_chunks[c].Clone();
                Volatile.Write(ref _chunks[c], copy);
                Interlocked.Increment(ref ChunkCopies);
                Volatile.Write(ref _owners[c], e);
                return copy;
            }

            spin.SpinOnce();
        }
    }

    /// <summary>Makes the directory private to the working version (before concurrent writes).</summary>
    public void PrepareWrite()
    {
        var e = _epoch.Value;
        if (_dirOwner == e) return;
        var chunks = new T[_chunks.Length][];
        Array.Copy(_chunks, chunks, _chunkCount);
        var owners = new int[_chunks.Length];
        Array.Copy(_owners, owners, _chunkCount);
        _chunks = chunks;
        _owners = owners;
        _dirOwner = e;
        DirectoryCopies++;
    }

    /// <summary>Grows the logical length; new entries take the fill value.</summary>
    public void Resize(int length)
    {
        if (_epoch.Concurrent && length > _length)
            throw new InvalidOperationException("arrays cannot grow during concurrent writes");
        if (length <= _length)
        {
            _length = Math.Max(length, _length);
            return;
        }

        var need = ChunksFor(length);
        var e = _epoch.Value;
        if (need > _chunks.Length || _dirOwner != e)
        {
            var capacity = Math.Max(need, _chunks.Length);
            if (need > _chunks.Length) capacity = Math.Max(need, _chunks.Length * 2);
            var chunks = new T[capacity][];
            var owners = new int[capacity];
            Array.Copy(_chunks, chunks, _chunkCount);
            Array.Copy(_owners, owners, _chunkCount);
            _chunks = chunks;
            _owners = owners;
            if (_dirOwner != e) DirectoryCopies++;
            _dirOwner = e;
        }

        for (var c = _chunkCount; c < need; c++)
        {
            var chunk = new T[_chunkSize];
            if (!EqualityComparer<T>.Default.Equals(_fill, default!)) Array.Fill(chunk, _fill);
            _chunks[c] = chunk;
            _owners[c] = e;
        }

        _chunkCount = Math.Max(_chunkCount, need);
        _length = length;
    }

    /// <summary>Captures the current contents; the caller must advance the epoch before writing again.</summary>
    public PView<T> Freeze() => new(_chunks, _owners, _shift, _length);

    /// <summary>Makes a published view the working state again (rollback).</summary>
    public void Restore(PView<T> view)
    {
        _chunks = view.Chunks;
        _owners = view.Owners;
        _length = view.Length;
        _chunkCount = ChunksFor(view.Length);
        _dirOwner = 0; // never equal to a live epoch: the next write copies the directory
    }

    /// <summary>Bulk load from a flat array into fresh chunks owned by the current epoch.</summary>
    public static PArray<T> FromFlat(Epoch epoch, int shift, T fill, ReadOnlySpan<T> data, int length)
    {
        var a = new PArray<T>(epoch, shift, fill, length);
        var n = Math.Min(data.Length, length);
        for (var start = 0; start < n; start += a._chunkSize)
        {
            var len = Math.Min(a._chunkSize, n - start);
            data.Slice(start, len).CopyTo(a._chunks[start >> shift]);
        }

        return a;
    }

    /// <summary>Approximate bytes held by the working version (directory + chunks).</summary>
    public long ApproximateBytes(int elementBytes) =>
        (long)_chunks.Length * 12 + (long)_chunkCount * (_chunkSize * (long)elementBytes + 24);
}

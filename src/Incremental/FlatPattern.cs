// FlatPattern.cs — flat copy of the pattern rows of a committed version, for solvers.
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
// A snapshot reads a pattern entry through the two levels of its chunked arrays. A solver
// that sweeps the pattern hundreds of times per version reads it faster from flat arrays.
// The copy keeps the slot layout of the gapped pattern (row start, row length, column per
// slot), so the values array and the destinations of the maintained structures address it
// unchanged. It follows the committed versions: in persistent mode a chunk shared by two
// versions is the same object and is never written again, so a refresh copies only the
// chunks whose reference differs from the version the copy holds (this covers edits,
// rollbacks and compactions alike); a transient snapshot aliases the working arrays, so
// every chunk is copied.

namespace Step2;

public sealed class FlatPattern
{
    /// <summary>Row start per node (-1: no row), row length per node, column per slot.</summary>
    public int[] RowStart = [], RowLen = [], Col = [];

    public int NodeCount { get; private set; }
    public int TotalSlots { get; private set; }

    /// <summary>Chunks copied by the last refresh, and chunks the three arrays hold.</summary>
    public int ChunksCopied { get; private set; }
    public int ChunksTotal { get; private set; }

    private PView<int> _rs, _rl, _col;
    private bool _valid;

    /// <summary>Makes the copy equal to the pattern of <paramref name="s"/>.</summary>
    public void Refresh(Snapshot s)
    {
        var incremental = _valid && s.Persistent;
        int copied = 0, total = 0;
        Sync(s.RowStart, _rs, incremental, ref RowStart, ref copied, ref total);
        Sync(s.RowLen, _rl, incremental, ref RowLen, ref copied, ref total);
        Sync(s.Col, _col, incremental, ref Col, ref copied, ref total);
        _rs = s.RowStart;
        _rl = s.RowLen;
        _col = s.Col;
        _valid = s.Persistent;
        NodeCount = s.NodeCount;
        TotalSlots = s.TotalSlots;
        ChunksCopied = copied;
        ChunksTotal = total;
    }

    /// <summary>
    ///     Copies the chunks of <paramref name="cur"/> that differ from <paramref name="prev"/>. A
    ///     chunk is skipped only when it is the same object in both views and lies within both
    ///     lengths, so the tail of a partly used chunk is always copied.
    /// </summary>
    private static void Sync(PView<int> cur, PView<int> prev, bool incremental, ref int[] flat, ref int copied,
        ref int total)
    {
        var len = cur.Length;
        if (flat.Length < len) Array.Resize(ref flat, Math.Max(len, flat.Length + (flat.Length >> 3)));
        var size = 1 << cur.Shift;
        var nc = (len + size - 1) >> cur.Shift;
        var covered = incremental ? Math.Min(len, prev.Length) : 0;
        var dst = flat;
        var n = 0;
        if (nc >= 256)
        {
            var counts = new int[(nc + 255) >> 8];
            Parallel.For(0, counts.Length, b =>
            {
                for (int c = b << 8, hi = Math.Min(nc, c + 256); c < hi; c++)
                    if (CopyChunk(cur, prev, c, size, len, covered, dst))
                        counts[b]++;
            });
            n = counts.Sum();
        }
        else
            for (var c = 0; c < nc; c++)
                if (CopyChunk(cur, prev, c, size, len, covered, dst))
                    n++;

        copied += n;
        total += nc;
    }

    private static bool CopyChunk(PView<int> cur, PView<int> prev, int c, int size, int len, int covered, int[] dst)
    {
        var start = c * size;
        if (start + size <= covered && ReferenceEquals(cur.Chunks[c], prev.Chunks[c])) return false;
        Array.Copy(cur.Chunks[c], 0, dst, start, Math.Min(size, len - start));
        return true;
    }
}

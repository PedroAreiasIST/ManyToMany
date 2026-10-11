// CornerCrack.cs — corner crack of a block with Erdogan-Sih growth in 3D, on the maintained structures.
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
//   corner [key=value ...]
//
// Body [0,W] x [0,W] x [-H,H] (E = 1, nu = 0.3), remote tension sigma = 1 along z on z = +-H, rigid motions removed by
// three corner constraints (the uniform field is exact). The faces x = 0 and y = 0 are free. A quarter-elliptical
// crack lies at the corner (0, 0) in the plane n . x = delta, n = (sin beta, 0, cos beta): beta > 0 loads it in
// mixed mode. Unlike the planar examples of Fracture.cs the front is a general 3D polyline with a local frame
// (N, B, T) at every knot (N the direction of the last advance, B the normal of the crack surface, T = N x B), and the
// crack surface is a graph h = w(p, q) over the initial plane, grown by strips between successive fronts.
//   Growth     at every knot K_I, K_II, K_III (domain interaction integral in the local frame, with the surface term of
//              Section 9.4 on the free faces) give the Erdogan-Sih kink angle
//                theta0 = 2 atan[(K_I - sqrt(K_I^2 + 8 K_II^2)) / (4 K_II)]      (in the plane (N, B), from N toward B)
//              and the advance da_j = da (K_eq,j / K_eq,max)^m with K_eq^2 = K_I^2 + K_II^2 + K_III^2 / (1 - nu); K_III
//              does not enter the direction. The knots move to P_j + da_j (cos theta N + sin theta B), the ends staying on
//              their faces, and are resampled to equal arc length.
//   Insertion  CrackInserter with phi_s = h - w(p, q) and phi_r = (2D signed distance to the front in (p, q)) - h_loc / 4.
//   Output     steps CSV, K along the front, fronts, and EnSight Gold (ASCII) cases of selected steps: the body (tetra4),
//              the upper crack face (tri3, opening) and the front (bar2); variables displacement, |u|, von Mises, opening.
//   "penny"    embedded circular crack of radius a0 in the same machinery (the front closed by repeating the first knot),
//              for the check against Kassir and Sih.

using System.Globalization;
using System.Text;
using Step2;

namespace Step2.Bench;

public static partial class Fracture
{
    private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

    // ------------------------------------------------------------------ the front

    /// <summary>Polyline front with the frame (N, B) at every knot; T = N x B.</summary>
    private sealed class CFront
    {
        public readonly int M;
        public readonly double[] P, N, B, S;

        public CFront(double[] p, double[] n, double[] b)
        {
            M = p.Length / 3;
            P = p;
            N = n;
            B = b;
            S = new double[M];
            for (var j = 1; j < M; j++)
            {
                double dx = P[3 * j] - P[3 * j - 3], dy = P[3 * j + 1] - P[3 * j - 2], dz = P[3 * j + 2] - P[3 * j - 1];
                S[j] = S[j - 1] + Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }
        }

        public double Length => S[M - 1];

        private (int K, double T, double D2) Search(double x, double y, double z, int lo, int hi)
        {
            var bk = lo;
            double bt = 0, bd = double.MaxValue;
            for (var k = lo; k <= hi; k++)
            {
                double ax = P[3 * k], ay = P[3 * k + 1], az = P[3 * k + 2];
                double dx = P[3 * k + 3] - ax, dy = P[3 * k + 4] - ay, dz = P[3 * k + 5] - az;
                var dd = dx * dx + dy * dy + dz * dz;
                var t = dd > 0 ? ((x - ax) * dx + (y - ay) * dy + (z - az) * dz) / dd : 0;
                t = t < 0 ? 0 : t > 1 ? 1 : t;
                double ex = x - ax - t * dx, ey = y - ay - t * dy, ez = z - az - t * dz;
                var d2 = ex * ex + ey * ey + ez * ez;
                if (d2 < bd)
                {
                    bd = d2;
                    bk = k;
                    bt = t;
                }
            }

            return (bk, bt, bd);
        }

        public (int K, double T, double D2) Nearest(double x, double y, double z, int hint = -1, int win = 4)
        {
            if (hint < 0) return Search(x, y, z, 0, M - 2);
            int lo = Math.Max(0, hint - win), hi = Math.Min(M - 2, hint + win);
            var r = Search(x, y, z, lo, hi);
            if ((r.K == lo && lo > 0) || (r.K == hi && hi < M - 2)) r = Search(x, y, z, 0, M - 2);
            return r;
        }

        public struct FootInfo
        {
            public int K;
            public double T, R1, R2, S, D;
            public double Nx, Ny, Nz, Bx, By, Bz, Tx, Ty, Tz;
        }

        /// <summary>Foot point of x on the polyline, its orthonormal frame (N, B, T) and the coordinates r1 = (x - F) . N, r2 = (x - F) . B.</summary>
        public FootInfo Foot(double x, double y, double z, int hint = -1)
        {
            var (k, t, d2) = Nearest(x, y, z, hint);
            return FootAt(x, y, z, k, t, d2);
        }

        public FootInfo FootAt(double x, double y, double z, int k, double t, double d2)
        {
            var f = new FootInfo { K = k, T = t, D = Math.Sqrt(d2) };
            double sx = P[3 * k + 3] - P[3 * k], sy = P[3 * k + 4] - P[3 * k + 1], sz = P[3 * k + 5] - P[3 * k + 2];
            var sl = Math.Sqrt(sx * sx + sy * sy + sz * sz);
            f.S = S[k] + t * sl;
            double fx = P[3 * k] + t * sx, fy = P[3 * k + 1] + t * sy, fz = P[3 * k + 2] + t * sz;
            sx /= sl;
            sy /= sl;
            sz /= sl;
            double nx = (1 - t) * N[3 * k] + t * N[3 * k + 3], ny = (1 - t) * N[3 * k + 1] + t * N[3 * k + 4],
                nz = (1 - t) * N[3 * k + 2] + t * N[3 * k + 5];
            double bx = (1 - t) * B[3 * k] + t * B[3 * k + 3], by = (1 - t) * B[3 * k + 1] + t * B[3 * k + 4],
                bz = (1 - t) * B[3 * k + 2] + t * B[3 * k + 5];
            var nt = nx * sx + ny * sy + nz * sz;
            nx -= nt * sx;
            ny -= nt * sy;
            nz -= nt * sz;
            var nl = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            nx /= nl;
            ny /= nl;
            nz /= nl;
            // B' = Ts x N' oriented as the interpolated B
            double cx = sy * nz - sz * ny, cy = sz * nx - sx * nz, cz = sx * ny - sy * nx;
            if (cx * bx + cy * by + cz * bz < 0)
            {
                cx = -cx;
                cy = -cy;
                cz = -cz;
            }

            f.Nx = nx;
            f.Ny = ny;
            f.Nz = nz;
            f.Bx = cx;
            f.By = cy;
            f.Bz = cz;
            // T = N x B
            f.Tx = ny * cz - nz * cy;
            f.Ty = nz * cx - nx * cz;
            f.Tz = nx * cy - ny * cx;
            double ex = x - fx, ey = y - fy, ez = z - fz;
            f.R1 = ex * nx + ey * ny + ez * nz;
            f.R2 = ex * cx + ey * cy + ez * cz;
            return f;
        }

        /// <summary>Integral of a hat of half-width delta (arc length) centered at s0, clipped to the front.</summary>
        public double HatLength(double s0, double delta)
        {
            double xa = Math.Max(-delta, -s0), xb = Math.Min(delta, Length - s0);
            if (xb <= xa) return 0;
            double F(double x) => x - x * Math.Abs(x) / (2 * delta);
            return F(xb) - F(xa);
        }
    }

    // ------------------------------------------------------------------ the crack surface

    /// <summary>
    ///     The crack surface h = w(p, q) over the initial plane on a regular grid, filled by strips (triangles of
    ///     successive fronts, first come first kept) and extended by dilation beyond the front; bilinear between the
    ///     grid nodes. The signed distance in (p, q) to the front, negative on the side of the crack, gives the region.
    /// </summary>
    private sealed class CSurface
    {
        private readonly double _p0, _q0, _dp;
        private readonly int _np, _nq;
        private readonly double[] _w, _ext;
        private readonly bool[] _set, _hasExt;
        private double[] _fp = [], _fq = []; // front in (p, q)
        private double[] _poly = [];          // polygon (p, q pairs) for the sign
        private double _hMean;

        public CSurface(double p0, double q0, double p1, double q1, double dp)
        {
            _p0 = p0;
            _q0 = q0;
            _dp = dp;
            _np = (int)Math.Ceiling((p1 - p0) / dp) + 1;
            _nq = (int)Math.Ceiling((q1 - q0) / dp) + 1;
            _w = new double[_np * _nq];
            _ext = new double[_np * _nq];
            _set = new bool[_np * _nq];
            _hasExt = new bool[_np * _nq];
        }

        public double Cell => _dp;

        /// <summary>A frozen copy (grid, extension, front and polygon) of the surface as it is.</summary>
        public CSurface Clone()
        {
            var c = new CSurface(_p0, _q0, _p0 + (_np - 1) * _dp, _q0 + (_nq - 1) * _dp, _dp);
            Array.Copy(_w, c._w, _w.Length);
            Array.Copy(_ext, c._ext, _ext.Length);
            Array.Copy(_set, c._set, _set.Length);
            Array.Copy(_hasExt, c._hasExt, _hasExt.Length);
            c._fp = (double[])_fp.Clone();
            c._fq = (double[])_fq.Clone();
            c._poly = (double[])_poly.Clone();
            c._hMean = _hMean;
            return c;
        }

        public void AddTriangle(ReadOnlySpan<double> a, ReadOnlySpan<double> b, ReadOnlySpan<double> c)
        {
            double pmin = Math.Min(a[0], Math.Min(b[0], c[0])), pmax = Math.Max(a[0], Math.Max(b[0], c[0]));
            double qmin = Math.Min(a[1], Math.Min(b[1], c[1])), qmax = Math.Max(a[1], Math.Max(b[1], c[1]));
            int i0 = Math.Max(0, (int)Math.Floor((pmin - _p0) / _dp)), i1 = Math.Min(_np - 1, (int)Math.Ceiling((pmax - _p0) / _dp));
            int j0 = Math.Max(0, (int)Math.Floor((qmin - _q0) / _dp)), j1 = Math.Min(_nq - 1, (int)Math.Ceiling((qmax - _q0) / _dp));
            var det = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);
            if (Math.Abs(det) < 1e-300) return;
            const double eps = 1e-9;
            for (var j = j0; j <= j1; j++)
            for (var i = i0; i <= i1; i++)
            {
                if (_set[i + _np * j]) continue;
                double p = _p0 + i * _dp, q = _q0 + j * _dp;
                var l1 = ((b[1] - c[1]) * (p - c[0]) + (c[0] - b[0]) * (q - c[1])) / det;
                var l2 = ((c[1] - a[1]) * (p - c[0]) + (a[0] - c[0]) * (q - c[1])) / det;
                var l3 = 1 - l1 - l2;
                if (l1 < -eps || l2 < -eps || l3 < -eps) continue;
                _w[i + _np * j] = l1 * a[2] + l2 * b[2] + l3 * c[2];
                _set[i + _np * j] = true;
            }
        }

        /// <summary>Dilates the filled nodes by a few cells (the extension beyond the front) and stores the front.</summary>
        public void Finish(double[] frontPQH, double[] polygonPQ, int cells)
        {
            var m = frontPQH.Length / 3;
            _fp = new double[m];
            _fq = new double[m];
            double hs = 0;
            for (var k = 0; k < m; k++)
            {
                _fp[k] = frontPQH[3 * k];
                _fq[k] = frontPQH[3 * k + 1];
                hs += frontPQH[3 * k + 2];
            }

            _hMean = hs / m;
            _poly = polygonPQ;
            Array.Clear(_hasExt);
            var cur = new bool[_set.Length];
            Array.Copy(_set, cur, _set.Length);
            for (var it = 0; it < cells; it++)
            {
                var next = (bool[])cur.Clone();
                var add = new List<(int, double)>();
                for (var j = 0; j < _nq; j++)
                for (var i = 0; i < _np; i++)
                {
                    var id = i + _np * j;
                    if (cur[id]) continue;
                    double s = 0;
                    var cnt = 0;
                    for (var dj = -1; dj <= 1; dj++)
                    for (var di = -1; di <= 1; di++)
                    {
                        int ii = i + di, jj = j + dj;
                        if (ii < 0 || jj < 0 || ii >= _np || jj >= _nq) continue;
                        var nid = ii + _np * jj;
                        if (!cur[nid]) continue;
                        s += _set[nid] ? _w[nid] : _ext[nid];
                        cnt++;
                    }

                    if (cnt > 0) add.Add((id, s / cnt));
                }

                foreach (var (id, v) in add)
                {
                    _ext[id] = v;
                    _hasExt[id] = true;
                    next[id] = true;
                }

                cur = next;
            }

            // every remaining cell takes the value of the nearest filled one (two chamfer sweeps), so that the extension
            // of the surface never has a zero of h - w inside the body away from the crack (a constant fallback did)
            var np = _np;
            var dist = new float[_set.Length];
            for (var id = 0; id < dist.Length; id++) dist[id] = cur[id] ? 0f : float.MaxValue;
            double V(int id) => _set[id] ? _w[id] : _ext[id];
            void Relax(int id, int nid, float cost)
            {
                if (dist[nid] == float.MaxValue || dist[nid] + cost >= dist[id]) return;
                dist[id] = dist[nid] + cost;
                _ext[id] = V(nid);
                _hasExt[id] = true;
            }

            for (var j = 0; j < _nq; j++)
            for (var i = 0; i < _np; i++)
            {
                var id = i + np * j;
                if (i > 0) Relax(id, id - 1, 1);
                if (j > 0) Relax(id, id - np, 1);
                if (i > 0 && j > 0) Relax(id, id - np - 1, 1.4142f);
                if (i + 1 < _np && j > 0) Relax(id, id - np + 1, 1.4142f);
            }

            for (var j = _nq - 1; j >= 0; j--)
            for (var i = _np - 1; i >= 0; i--)
            {
                var id = i + np * j;
                if (i + 1 < _np) Relax(id, id + 1, 1);
                if (j + 1 < _nq) Relax(id, id + np, 1);
                if (i + 1 < _np && j + 1 < _nq) Relax(id, id + np + 1, 1.4142f);
                if (i > 0 && j + 1 < _nq) Relax(id, id + np - 1, 1.4142f);
            }
        }

        private double Node(int i, int j, double fallback)
        {
            i = Math.Clamp(i, 0, _np - 1);
            j = Math.Clamp(j, 0, _nq - 1);
            var id = i + _np * j;
            return _set[id] ? _w[id] : _hasExt[id] ? _ext[id] : fallback;
        }

        public double W(double p, double q)
        {
            var fx = (p - _p0) / _dp;
            var fy = (q - _q0) / _dp;
            int i = (int)Math.Floor(fx), j = (int)Math.Floor(fy);
            double tx = fx - i, ty = fy - j;
            var fb = _hMean;
            return (1 - ty) * ((1 - tx) * Node(i, j, fb) + tx * Node(i + 1, j, fb)) +
                   ty * ((1 - tx) * Node(i, j + 1, fb) + tx * Node(i + 1, j + 1, fb));
        }

        /// <summary>Signed distance in (p, q) to the front: negative on the crack side (inside the polygon).</summary>
        public double SignedDistance(double p, double q)
        {
            double best = double.MaxValue;
            for (var k = 0; k + 1 < _fp.Length; k++)
            {
                double ax = _fp[k], ay = _fq[k], dx = _fp[k + 1] - ax, dy = _fq[k + 1] - ay;
                var dd = dx * dx + dy * dy;
                var t = dd > 0 ? ((p - ax) * dx + (q - ay) * dy) / dd : 0;
                t = t < 0 ? 0 : t > 1 ? 1 : t;
                double ex = p - ax - t * dx, ey = q - ay - t * dy;
                best = Math.Min(best, ex * ex + ey * ey);
            }

            var inside = false;
            var n = _poly.Length / 2;
            for (int i = 0, j = n - 1; i < n; j = i++)
            {
                double xi = _poly[2 * i], yi = _poly[2 * i + 1], xj = _poly[2 * j], yj = _poly[2 * j + 1];
                if ((yi > q) != (yj > q) && p < (xj - xi) * (q - yi) / (yj - yi) + xi) inside = !inside;
            }

            return inside ? -Math.Sqrt(best) : Math.Sqrt(best);
        }
    }

    // ------------------------------------------------------------------ mesh and problem

    private static (int[][] Rows, List<double> Coords, int NodeCount) CornerBox(int nx, int ny, int nz, double w, double hh)
    {
        int[][] perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
        int Node(int i, int j, int k) => i + (nx + 1) * (j + (ny + 1) * k);
        var nN = (nx + 1) * (ny + 1) * (nz + 1);
        var coords = new List<double>(3 * nN);
        for (var k = 0; k <= nz; k++)
        for (var j = 0; j <= ny; j++)
        for (var i = 0; i <= nx; i++)
        {
            coords.Add(w * i / nx);
            coords.Add(w * j / ny);
            coords.Add(-hh + 2 * hh * k / nz);
        }

        var rows = new int[6 * nx * ny * nz][];
        var t = 0;
        var c = new int[3];
        for (var k = 0; k < nz; k++)
        for (var j = 0; j < ny; j++)
        for (var i = 0; i < nx; i++)
            foreach (var perm in perms)
            {
                c[0] = i;
                c[1] = j;
                c[2] = k;
                var row = new int[4];
                row[0] = Node(c[0], c[1], c[2]);
                for (var s = 0; s < 3; s++)
                {
                    c[perm[s]]++;
                    row[s + 1] = Node(c[0], c[1], c[2]);
                }

                rows[t++] = row;
            }

        return (rows, coords, nN);
    }

    /// <summary>The block refined by recursive bisection until no element is longer than 1.5 size(x).</summary>
    private static (IncrementalMesh M, MeshState St, int Passes, long Bisections) RefinedCorner(int nx, int ny, int nz,
        double w, double hh, Func<double, double, double, double> size)
    {
        var (rows, coords, nN) = CornerBox(nx, ny, nz, w, hh);
        var flat = FlatBuild.Initial(ElementSpec.Tetrahedron, nN, rows);
        var m = new IncrementalMesh(flat, ElementSpec.Tetrahedron, persistent: true);
        var st = new MeshState(m, coords, null);
        var bis = new Bisection();
        var passes = 0;
        while (true)
        {
            var marked = new List<int>();
            for (var e = 0; e < m.ElementCount; e++)
            {
                var row = m.Row(e);
                if (row.Length != 4) continue;
                double longest = 0;
                for (var i = 0; i < 4; i++)
                for (var j = i + 1; j < 4; j++)
                {
                    double dx = st.X(row[i], 0) - st.X(row[j], 0), dy = st.X(row[i], 1) - st.X(row[j], 1),
                        dz = st.X(row[i], 2) - st.X(row[j], 2);
                    longest = Math.Max(longest, Math.Sqrt(dx * dx + dy * dy + dz * dz));
                }

                var (cx, cy, cz) = st.Centroid(row);
                if (longest > 1.5 * size(cx, cy, cz)) marked.Add(e);
            }

            if (marked.Count == 0) break;
            bis.Refine(st, marked, maxLevel: 1000);
            passes++;
        }

        m.Commit();
        return (m, st, passes, bis.Bisections);
    }

    /// <summary>Nodal coordinates, tractions sigma on z = +-H and the three corner constraints of the block.</summary>
    private static (double[] X, double[] F, bool[] Fixed) CornerProblem(Snapshot s, CrackInserter ins, int nx, int ny, double hh)
    {
        var nNodes = s.NodeCount;
        var X = new double[3 * nNodes];
        var xv = s.View(ins.X);
        var yv = s.View(ins.Y);
        var zv = s.View(ins.Z);
        for (var v = 0; v < nNodes; v++)
        {
            X[3 * v] = xv[v];
            X[3 * v + 1] = yv[v];
            X[3 * v + 2] = zv[v];
        }

        var f = new double[3 * nNodes];
        var t = new int[3];
        for (var i = 0; i < s.BoundaryCount; i++)
        {
            var id = s.BoundaryFacet(i);
            int e = id / 4, fl = id % 4;
            var row = s.Row(e);
            var k = 0;
            for (var p = 0; p < 4; p++)
                if (p != fl)
                    t[k++] = row[p];
            var zTop = t.All(v => Math.Abs(X[3 * v + 2] - hh) < 1e-12);
            var zBot = t.All(v => Math.Abs(X[3 * v + 2] + hh) < 1e-12);
            if (!zTop && !zBot) continue;
            var area = TriangleArea(X, t[0], t[1], t[2]);
            foreach (var v in t) f[3 * v + 2] += (zTop ? Sigma : -Sigma) * area / 3;
        }

        int Node(int i, int j, int k) => i + (nx + 1) * (j + (ny + 1) * k);
        var fixedDof = new bool[3 * nNodes];
        int cA = Node(0, 0, 0), cB = Node(nx, 0, 0), cC = Node(0, ny, 0);
        fixedDof[3 * cA] = fixedDof[3 * cA + 1] = fixedDof[3 * cA + 2] = true;
        fixedDof[3 * cB + 1] = fixedDof[3 * cB + 2] = true;
        fixedDof[3 * cC + 2] = true;
        return (X, f, fixedDof);
    }

    // ------------------------------------------------------------------ auxiliary fields in the local frame

    private static void CAuxH(int mode, double x, double y, double z, CFront fr, int hint, int side, Span<double> Hg)
    {
        var ft = fr.Foot(x, y, z, hint);
        double r1 = ft.R1, r2 = ft.R2;
        var d = 1e-6 * Math.Max(Math.Sqrt(r1 * r1 + r2 * r2), 1e-12);
        Span<double> up = stackalloc double[3], um = stackalloc double[3];
        Span<double> Hl = stackalloc double[9];
        Hl.Clear();
        for (var j = 0; j < 2; j++)
        {
            Aux(mode, r1 + (j == 0 ? d : 0), r2 + (j == 1 ? d : 0), side, up);
            Aux(mode, r1 - (j == 0 ? d : 0), r2 - (j == 1 ? d : 0), side, um);
            for (var i = 0; i < 3; i++) Hl[3 * i + j] = (up[i] - um[i]) / (2 * d);
        }

        Span<double> Bm = stackalloc double[9];
        Bm[0] = ft.Nx; Bm[3] = ft.Ny; Bm[6] = ft.Nz;
        Bm[1] = ft.Bx; Bm[4] = ft.By; Bm[7] = ft.Bz;
        Bm[2] = ft.Tx; Bm[5] = ft.Ty; Bm[8] = ft.Tz;
        for (var i = 0; i < 3; i++)
        for (var k = 0; k < 3; k++)
        {
            double acc = 0;
            for (var aa = 0; aa < 3; aa++)
            for (var bb = 0; bb < 2; bb++)
                acc += Bm[3 * i + aa] * Hl[3 * aa + bb] * Bm[3 * k + bb];
            Hg[3 * i + k] = acc;
        }
    }

    private static void CAuxDerivatives(int mode, double x, double y, double z, CFront fr, int hint, int side,
        Span<double> dH, Span<double> div)
    {
        var ft = fr.Foot(x, y, z, hint);
        var d = 1e-3 * Math.Max(Math.Sqrt(ft.R1 * ft.R1 + ft.R2 * ft.R2), 1e-9);
        Span<double> Hp = stackalloc double[9], Hm = stackalloc double[9];
        div.Clear();
        for (var j = 0; j < 3; j++)
        {
            CAuxH(mode, x + (j == 0 ? d : 0), y + (j == 1 ? d : 0), z + (j == 2 ? d : 0), fr, ft.K, side, Hp);
            CAuxH(mode, x - (j == 0 ? d : 0), y - (j == 1 ? d : 0), z - (j == 2 ? d : 0), fr, ft.K, side, Hm);
            for (var t = 0; t < 9; t++) dH[9 * j + t] = (Hp[t] - Hm[t]) / (2 * d);
            double trp = Hp[0] + Hp[4] + Hp[8], trm = Hm[0] + Hm[4] + Hm[8];
            for (var i = 0; i < 3; i++)
            {
                var sp = Mu * (Hp[3 * i + j] + Hp[3 * j + i]) + (i == j ? Lambda * trp : 0);
                var sm = Mu * (Hm[3 * i + j] + Hm[3 * j + i]) + (i == j ? Lambda * trm : 0);
                div[i] += (sp - sm) / (2 * d);
            }
        }
    }

    // ------------------------------------------------------------------ interaction integral along a 3D front

    private sealed class Tube
    {
        public int[] Nodes = [];
        public int[] Elements = [];
        public int[] NodeSeg = [];
        public double[] NodeS = [], NodeQ = [], NodeDir = [];
        public int[] ElemHint = [];
        public double[] ElemSMin = [], ElemSMax = [];
    }

    private static Tube BuildTube(Snapshot s, double[] X, CFront fr, double rq)
    {
        var nN = s.NodeCount;
        var tb = new Tube { NodeSeg = new int[nN], NodeS = new double[nN], NodeQ = new double[nN], NodeDir = new double[3 * nN] };
        Array.Fill(tb.NodeSeg, -1);
        // bounding box of the front inflated by rq
        double x0 = double.MaxValue, y0 = x0, z0 = x0, x1 = -x0, y1 = x1, z1 = x1;
        for (var k = 0; k < fr.M; k++)
        {
            x0 = Math.Min(x0, fr.P[3 * k]); x1 = Math.Max(x1, fr.P[3 * k]);
            y0 = Math.Min(y0, fr.P[3 * k + 1]); y1 = Math.Max(y1, fr.P[3 * k + 1]);
            z0 = Math.Min(z0, fr.P[3 * k + 2]); z1 = Math.Max(z1, fr.P[3 * k + 2]);
        }

        x0 -= rq; y0 -= rq; z0 -= rq; x1 += rq; y1 += rq; z1 += rq;
        var inTube = new bool[nN];
        Parallel.For(0, nN, v =>
        {
            if (s.TransposePairs(v).Length == 0) return;
            double x = X[3 * v], y = X[3 * v + 1], z = X[3 * v + 2];
            if (x < x0 || x > x1 || y < y0 || y > y1 || z < z0 || z > z1) return;
            var (k, t, d2) = fr.Nearest(x, y, z);
            if (d2 >= rq * rq) return;
            var ft = fr.FootAt(x, y, z, k, t, d2);
            var d = ft.D;
            var tube = QPlateau <= 0 ? Math.Max(0, 1 - d / rq)
                : d <= QPlateau * rq ? 1.0 : Math.Max(0, 1 - (d - QPlateau * rq) / (rq - QPlateau * rq));
            if (tube <= 0) return;
            tb.NodeSeg[v] = k;
            tb.NodeS[v] = ft.S;
            tb.NodeQ[v] = tube;
            tb.NodeDir[3 * v] = ft.Nx;
            tb.NodeDir[3 * v + 1] = ft.Ny;
            tb.NodeDir[3 * v + 2] = ft.Nz;
            inTube[v] = true;
        });
        var nodes = new List<int>();
        for (var v = 0; v < nN; v++)
            if (inTube[v])
                nodes.Add(v);
        tb.Nodes = nodes.ToArray();
        var els = new List<int>();
        var hint = new List<int>();
        var smin = new List<double>();
        var smax = new List<double>();
        for (var e = 0; e < s.ElementCount; e++)
        {
            var row = s.Row(e);
            if (row.Length != 4) continue;
            var any = false;
            int h = -1;
            double lo = double.MaxValue, hi = -1;
            foreach (var v in row)
                if (inTube[v])
                {
                    any = true;
                    if (h < 0) h = tb.NodeSeg[v];
                    lo = Math.Min(lo, tb.NodeS[v]);
                    hi = Math.Max(hi, tb.NodeS[v]);
                }

            if (!any) continue;
            els.Add(e);
            hint.Add(h);
            smin.Add(lo);
            smax.Add(hi);
        }

        tb.Elements = els.ToArray();
        tb.ElemHint = hint.ToArray();
        tb.ElemSMin = smin.ToArray();
        tb.ElemSMax = smax.ToArray();
        return tb;
    }

    /// <summary>Interaction integral at the knot with arc length s0: (K_I, K_II, K_III, J).</summary>
    private static (double KI, double KII, double KIII, double J) CInteraction(Snapshot s, double[] X, double[] u, CFront fr,
        Tube tb, double s0, double delta, bool freeSurface, double wBox, double[] vn, bool[] any)
    {
        // virtual extension V = q N at the tube nodes: tube factor times a hat of half-width delta in arc length
        var lo = s0 - delta;
        var hi = s0 + delta;
        var sNodes = tb.Nodes;
        foreach (var v in sNodes)
        {
            var ds = Math.Abs(tb.NodeS[v] - s0);
            var q = ds >= delta ? 0 : tb.NodeQ[v] * (1 - ds / delta);
            if (q <= 0)
            {
                any[v] = false;
                continue;
            }

            any[v] = true;
            vn[3 * v] = q * tb.NodeDir[3 * v];
            vn[3 * v + 1] = q * tb.NodeDir[3 * v + 1];
            vn[3 * v + 2] = q * tb.NodeDir[3 * v + 2];
        }

        const double qa = 0.5854101966249685, qb = 0.1381966011250105;
        var sel = new List<int>();
        for (var i = 0; i < tb.Elements.Length; i++)
            if (tb.ElemSMax[i] >= lo && tb.ElemSMin[i] <= hi)
                sel.Add(i);
        var block = KBlock;
        var nb = (sel.Count + block - 1) / block;
        var partial = new double[nb][];
        Parallel.For(0, nb, b =>
        {
            var part = new double[4];
            for (int k = b * block, end = Math.Min(sel.Count, k + block); k < end; k++) Term(sel[k], part);
            partial[b] = part;
        });
        var sum = new double[4];
        foreach (var part in partial)
            for (var k = 0; k < 4; k++)
                sum[k] += part[k];

        void Term(int idx, double[] part)
        {
            int e = tb.Elements[idx], hint = tb.ElemHint[idx];
            var row = s.Row(e);
            var touched = false;
            foreach (var v in row) touched |= any[v];
            if (!touched) return;
            Span<double> g = stackalloc double[12];
            Span<double> G = stackalloc double[9];
            Span<double> S = stackalloc double[9];
            Span<double> GV = stackalloc double[9];
            Span<double> H = stackalloc double[9];
            Span<double> dv = stackalloc double[3];
            Span<double> dH = stackalloc double[27];
            Span<double> Vq = stackalloc double[3];
            var vol = Gradients(X, row, g);
            double cxs = 0, cys = 0, czs = 0;
            foreach (var v in row)
            {
                cxs += X[3 * v];
                cys += X[3 * v + 1];
                czs += X[3 * v + 2];
            }

            cxs /= 4;
            cys /= 4;
            czs /= 4;
            var side = fr.Foot(cxs, cys, czs, hint).R2 >= 0 ? 1 : -1;
            G.Clear();
            GV.Clear();
            for (var a = 0; a < 4; a++)
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
            {
                G[3 * i + j] += u[3 * row[a] + i] * g[3 * a + j];
                GV[3 * i + j] += (any[row[a]] ? vn[3 * row[a] + i] : 0) * g[3 * a + j];
            }

            var tr = G[0] + G[4] + G[8];
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
                S[3 * i + j] = Mu * (G[3 * i + j] + G[3 * j + i]) + (i == j ? Lambda * tr : 0);
            double W = 0, divV = GV[0] + GV[4] + GV[8];
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
                W += 0.5 * S[3 * i + j] * 0.5 * (G[3 * i + j] + G[3 * j + i]);
            var jv = -W * divV;
            for (var i = 0; i < 3; i++)
            for (var j = 0; j < 3; j++)
            for (var k = 0; k < 3; k++)
                jv += S[3 * i + j] * G[3 * i + k] * GV[3 * k + j];
            part[3] += jv * vol;

            for (var qp = 0; qp < 4; qp++)
            {
                double xx = 0, yy = 0, zz = 0;
                Vq.Clear();
                for (var a = 0; a < 4; a++)
                {
                    var wa = a == qp ? qa : qb;
                    xx += wa * X[3 * row[a]];
                    yy += wa * X[3 * row[a] + 1];
                    zz += wa * X[3 * row[a] + 2];
                    if (any[row[a]])
                        for (var k = 0; k < 3; k++)
                            Vq[k] += wa * vn[3 * row[a] + k];
                }

                for (var mode = 0; mode < 3; mode++)
                {
                    CAuxH(mode, xx, yy, zz, fr, hint, side, H);
                    var htr = H[0] + H[4] + H[8];
                    double acc = 0, wint = 0;
                    for (var i = 0; i < 3; i++)
                    for (var j = 0; j < 3; j++)
                    {
                        var sa = Mu * (H[3 * i + j] + H[3 * j + i]) + (i == j ? Lambda * htr : 0);
                        wint += S[3 * i + j] * 0.5 * (H[3 * i + j] + H[3 * j + i]);
                        for (var k = 0; k < 3; k++)
                            acc += (S[3 * i + j] * H[3 * i + k] + sa * G[3 * i + k]) * GV[3 * k + j];
                    }

                    acc -= wint * divV;
                    CAuxDerivatives(mode, xx, yy, zz, fr, hint, side, dH, dv);
                    for (var k = 0; k < 3; k++)
                    {
                        double ck = 0;
                        for (var i = 0; i < 3; i++)
                        {
                            ck += dv[i] * G[3 * i + k];
                            for (var j = 0; j < 3; j++)
                                ck += S[3 * i + j] * (dH[9 * j + 3 * i + k] - dH[9 * k + 3 * i + j]);
                        }

                        acc += ck * Vq[k];
                    }

                    part[mode] += acc * vol / 4;
                }
            }
        }

        if (freeSurface) FreeFaces(sum);

        // free faces x = 0, y = 0, x = W, y = W inside the tube: the surface term of Section 9.4
        void FreeFaces(double[] acc)
        {
            Span<double> g = stackalloc double[12];
            Span<double> G = stackalloc double[9];
            Span<double> S = stackalloc double[9];
            Span<double> H = stackalloc double[9];
            Span<double> V = stackalloc double[3];
            var t = new int[3];
            var bsum = new double[3];
            for (var i = 0; i < s.BoundaryCount; i++)
            {
                var id = s.BoundaryFacet(i);
                int e = id / 4, fl = id % 4;
                var row = s.Row(e);
                if (row.Length != 4) continue;
                var k0 = 0;
                for (var p = 0; p < 4; p++)
                    if (p != fl)
                        t[k0++] = row[p];
                if (!(any[t[0]] || any[t[1]] || any[t[2]])) continue;
                const double tol = 1e-9;
                var onX0 = t.All(v => Math.Abs(X[3 * v]) < tol);
                var onY0 = t.All(v => Math.Abs(X[3 * v + 1]) < tol);
                var onX1 = t.All(v => Math.Abs(X[3 * v] - wBox) < tol);
                var onY1 = t.All(v => Math.Abs(X[3 * v + 1] - wBox) < tol);
                if (!(onX0 || onY0 || onX1 || onY1)) continue;
                double mx = onX0 ? -1 : onX1 ? 1 : 0, my = onY0 ? -1 : onY1 ? 1 : 0;
                double ux = X[3 * t[1]] - X[3 * t[0]], uy = X[3 * t[1] + 1] - X[3 * t[0] + 1], uz = X[3 * t[1] + 2] - X[3 * t[0] + 2];
                double wx = X[3 * t[2]] - X[3 * t[0]], wy = X[3 * t[2] + 1] - X[3 * t[0] + 1], wz = X[3 * t[2] + 2] - X[3 * t[0] + 2];
                double cx = uy * wz - uz * wy, cy = uz * wx - ux * wz, cz = ux * wy - uy * wx;
                var area = 0.5 * Math.Sqrt(cx * cx + cy * cy + cz * cz);
                if (area <= 0) continue;
                Gradients(X, row, g);
                G.Clear();
                for (var a = 0; a < 4; a++)
                for (var ii = 0; ii < 3; ii++)
                for (var j = 0; j < 3; j++)
                    G[3 * ii + j] += u[3 * row[a] + ii] * g[3 * a + j];
                var tr = G[0] + G[4] + G[8];
                for (var ii = 0; ii < 3; ii++)
                for (var j = 0; j < 3; j++)
                    S[3 * ii + j] = Mu * (G[3 * ii + j] + G[3 * j + ii]) + (ii == j ? Lambda * tr : 0);
                double cxs = 0, cys = 0, czs = 0;
                foreach (var v in row)
                {
                    cxs += X[3 * v];
                    cys += X[3 * v + 1];
                    czs += X[3 * v + 2];
                }

                var hint = tb.NodeSeg[t[0]] >= 0 ? tb.NodeSeg[t[0]] : tb.NodeSeg[t[1]] >= 0 ? tb.NodeSeg[t[1]] : tb.NodeSeg[t[2]];
                var side = fr.Foot(cxs / 4, cys / 4, czs / 4, hint).R2 >= 0 ? 1 : -1;
                for (var qp = 0; qp < 3; qp++)
                {
                    double xx = 0, yy = 0, zz = 0;
                    V.Clear();
                    for (var a = 0; a < 3; a++)
                    {
                        var wa = a == qp ? 2.0 / 3 : 1.0 / 6;
                        xx += wa * X[3 * t[a]];
                        yy += wa * X[3 * t[a] + 1];
                        zz += wa * X[3 * t[a] + 2];
                        if (any[t[a]])
                            for (var k = 0; k < 3; k++)
                                V[k] += wa * vn[3 * t[a] + k];
                    }

                    for (var mode = 0; mode < 3; mode++)
                    {
                        CAuxH(mode, xx, yy, zz, fr, hint, side, H);
                        var htr = H[0] + H[4] + H[8];
                        double wint = 0, mv = mx * V[0] + my * V[1], b = 0;
                        for (var ii = 0; ii < 3; ii++)
                        {
                            double sam = 0;
                            for (var j = 0; j < 3; j++)
                            {
                                var sa = Mu * (H[3 * ii + j] + H[3 * j + ii]) + (ii == j ? Lambda * htr : 0);
                                sam += sa * (j == 0 ? mx : j == 1 ? my : 0);
                                wint += S[3 * ii + j] * 0.5 * (H[3 * ii + j] + H[3 * j + ii]);
                            }

                            for (var k = 0; k < 3; k++) b += sam * G[3 * ii + k] * V[k];
                        }

                        bsum[mode] += (b - wint * mv) * area / 3;
                    }
                }
            }

            for (var mode = 0; mode < 3; mode++) acc[mode] -= bsum[mode];
        }

        foreach (var v in sNodes) any[v] = false;
        var l = fr.HatLength(s0, delta);
        var ep = E / (1 - Nu * Nu);
        return (ep / 2 * sum[0] / l, ep / 2 * sum[1] / l, Mu * sum[2] / l, sum[3] / l);
    }

    // ------------------------------------------------------------------ Erdogan-Sih growth

    /// <summary>Kink angle of the maximum circumferential stress criterion (Erdogan and Sih 1963), from N toward B.</summary>
    private static double ErdoganSih(double kI, double kII)
    {
        if (Math.Abs(kII) < 1e-14 * Math.Max(1e-300, Math.Abs(kI))) return 0;
        return 2 * Math.Atan((kI - Math.Sqrt(kI * kI + 8 * kII * kII)) / (4 * kII));
    }

    /// <summary>Equal-arc-length resampling of the polyline P (3 M doubles) into M knots.</summary>
    private static double[] Resample(double[] P, int m, double[]? companion = null, double[]? companionOut = null)
    {
        var n = P.Length / 3;
        var s = new double[n];
        for (var j = 1; j < n; j++)
        {
            double dx = P[3 * j] - P[3 * j - 3], dy = P[3 * j + 1] - P[3 * j - 2], dz = P[3 * j + 2] - P[3 * j - 1];
            s[j] = s[j - 1] + Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        var o = new double[3 * m];
        var seg = 0;
        for (var j = 0; j < m; j++)
        {
            var target = s[n - 1] * j / (m - 1);
            while (seg < n - 2 && s[seg + 1] < target) seg++;
            var t = s[seg + 1] > s[seg] ? (target - s[seg]) / (s[seg + 1] - s[seg]) : 0;
            for (var c = 0; c < 3; c++) o[3 * j + c] = (1 - t) * P[3 * seg + c] + t * P[3 * seg + 3 + c];
            if (companion != null && companionOut != null)
                for (var c = 0; c < 3; c++) companionOut[3 * j + c] = (1 - t) * companion[3 * seg + c] + t * companion[3 * seg + 3 + c];
        }

        return o;
    }

    /// <summary>
    ///     Frames of a new front from the old one: T from the polyline, N from the vector to the nearest point of the old
    ///     front (the direction of the advance) made orthogonal to T, B = T x N oriented toward <paramref name="up"/>.
    /// </summary>
    private static CFront MakeFront(double[] P, CFront? old, double[] up, Func<int, (double, double, double)>? initialN = null,
        double[]? nOverride = null)
    {
        var m = P.Length / 3;
        var N = new double[3 * m];
        var B = new double[3 * m];
        for (var j = 0; j < m; j++)
        {
            int a = Math.Max(0, j - 1), b = Math.Min(m - 1, j + 1);
            double tx = P[3 * b] - P[3 * a], ty = P[3 * b + 1] - P[3 * a + 1], tz = P[3 * b + 2] - P[3 * a + 2];
            var tl = Math.Sqrt(tx * tx + ty * ty + tz * tz);
            tx /= tl;
            ty /= tl;
            tz /= tl;
            double nx, ny, nz;
            if (nOverride != null)
            {
                nx = nOverride[3 * j];
                ny = nOverride[3 * j + 1];
                nz = nOverride[3 * j + 2];
            }
            else if (old == null)
            {
                (nx, ny, nz) = initialN!(j);
            }
            else
            {
                var (k, t, _) = old.Nearest(P[3 * j], P[3 * j + 1], P[3 * j + 2]);
                nx = P[3 * j] - ((1 - t) * old.P[3 * k] + t * old.P[3 * k + 3]);
                ny = P[3 * j + 1] - ((1 - t) * old.P[3 * k + 1] + t * old.P[3 * k + 4]);
                nz = P[3 * j + 2] - ((1 - t) * old.P[3 * k + 2] + t * old.P[3 * k + 5]);
            }

            var nt = nx * tx + ny * ty + nz * tz;
            nx -= nt * tx;
            ny -= nt * ty;
            nz -= nt * tz;
            var nl = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            nx /= nl;
            ny /= nl;
            nz /= nl;
            double bx = ty * nz - tz * ny, by = tz * nx - tx * nz, bz = tx * ny - ty * nx;
            if (bx * up[0] + by * up[1] + bz * up[2] < 0)
            {
                bx = -bx;
                by = -by;
                bz = -bz;
            }

            N[3 * j] = nx; N[3 * j + 1] = ny; N[3 * j + 2] = nz;
            B[3 * j] = bx; B[3 * j + 1] = by; B[3 * j + 2] = bz;
        }

        return new CFront(P, N, B);
    }

    // ------------------------------------------------------------------ EnSight Gold (ASCII)

    private static string F12(double v) => v.ToString("0.00000e+00", Ci).PadLeft(12);

    private sealed class EPart
    {
        public int Id;
        public string Name = "";
        public int[] Nodes = [];          // global node (or front knot) behind each local coordinate
        public double[] Xyz = [];         // 3 * count coordinates
        public string Type = "";
        public int[] Conn = [];           // 1-based local connectivity
        public int NPerElem;
        public bool IsBody, IsCrack;
    }

    /// <summary>
    ///     EnSight Gold ASCII files of one step in the conventions of ist_ensight_gold.py: prefix.geoNNNN with the parts
    ///     body (tetra4, optional), crack_face (tria3, the upper face), front (bar2) and outline (bar2, the 12 edges of the
    ///     block), and the node variables displacement, umag, mises and opening in prefix.dispNNNN, .umagNNNN, .misesNNNN,
    ///     .openNNNN. The case file is written by <see cref="WriteEnsightCase"/>.
    /// </summary>
    private static int WriteEnsightStep(string dir, string prefix, bool withBody, int stepNo, Snapshot s, CrackInserter ins,
        double[] X, double[] u, CFront fr, double wBox, double hBox, double[] opening, double[] vmElem)
    {
        Directory.CreateDirectory(dir);
        var tag = stepNo.ToString("D4", Ci);
        var nNodes = s.NodeCount;
        var parts = new List<EPart>();

        // body: live tetrahedra and the nodes they use
        if (withBody)
        {
            var live = new List<int>();
            for (var e = 0; e < s.ElementCount; e++)
                if (s.Row(e).Length == 4)
                    live.Add(e);
            var used = new bool[nNodes];
            foreach (var e in live)
                foreach (var v in s.Row(e))
                    used[v] = true;
            var newId = new int[nNodes];
            var order = new List<int>();
            for (var v = 0; v < nNodes; v++)
                if (used[v])
                {
                    order.Add(v);
                    newId[v] = order.Count;
                }

            var conn = new int[4 * live.Count];
            for (var k = 0; k < live.Count; k++)
            {
                var r = s.Row(live[k]);
                for (var c = 0; c < 4; c++) conn[4 * k + c] = newId[r[c]];
            }

            parts.Add(new EPart { Name = "body", Nodes = order.ToArray(), Type = "tetra4", Conn = conn, NPerElem = 4, IsBody = true });
        }

        // crack face: boundary triangles whose three nodes are surface nodes on the + side of their twin (the copies)
        var isUp = new bool[nNodes];
        for (var v = 0; v < nNodes; v++)
        {
            var p1 = ins.Parent1[v];
            if (p1 != v && p1 == ins.Parent2[v]) isUp[v] = true;
        }

        var tris = new List<int[]>();
        for (var i = 0; i < s.BoundaryCount; i++)
        {
            var id = s.BoundaryFacet(i);
            int e = id / 4, fl = id % 4;
            var row = s.Row(e);
            if (row.Length != 4) continue;
            var t = new int[3];
            var k = 0;
            for (var p = 0; p < 4; p++)
                if (p != fl)
                    t[k++] = row[p];
            if (!t.All(v => isUp[v] && ins.OnSurface[v] != 0)) continue;
            // a triangle of a free face (all nodes on x = 0, y = 0, x = W or y = W) is not part of the crack surface
            const double ft = 1e-9;
            if (t.All(v => Math.Abs(X[3 * v]) < ft) || t.All(v => Math.Abs(X[3 * v + 1]) < ft) ||
                t.All(v => Math.Abs(X[3 * v] - wBox) < ft) || t.All(v => Math.Abs(X[3 * v + 1] - wBox) < ft)) continue;
            tris.Add(t);
        }

        var cmap = new Dictionary<int, int>();
        var corder = new List<int>();
        var cconn = new int[3 * tris.Count];
        for (var k = 0; k < tris.Count; k++)
            for (var c = 0; c < 3; c++)
            {
                var v = tris[k][c];
                if (!cmap.TryGetValue(v, out var lid))
                {
                    corder.Add(v);
                    lid = corder.Count;
                    cmap[v] = lid;
                }

                cconn[3 * k + c] = lid;
            }

        parts.Add(new EPart { Name = "crack_face", Nodes = corder.ToArray(), Type = "tria3", Conn = cconn, NPerElem = 3, IsCrack = true });

        // front and block outline (bars, not attached to the mesh: zero fields)
        var fconn = new int[2 * (fr.M - 1)];
        for (var j = 0; j + 1 < fr.M; j++)
        {
            fconn[2 * j] = j + 1;
            fconn[2 * j + 1] = j + 2;
        }

        parts.Add(new EPart { Name = "front", Nodes = new int[fr.M], Xyz = (double[])fr.P.Clone(), Type = "bar2", Conn = fconn, NPerElem = 2 });
        var corners = new double[24];
        for (var c = 0; c < 8; c++)
        {
            corners[3 * c] = (c & 1) != 0 ? wBox : 0;
            corners[3 * c + 1] = (c & 2) != 0 ? wBox : 0;
            corners[3 * c + 2] = (c & 4) != 0 ? hBox : -hBox;
        }

        int[] edges = [1, 2, 1, 3, 2, 4, 3, 4, 5, 6, 5, 7, 6, 8, 7, 8, 1, 5, 2, 6, 3, 7, 4, 8];
        if (withBody) parts.Add(new EPart { Name = "outline", Nodes = new int[8], Xyz = corners, Type = "bar2", Conn = edges, NPerElem = 2 });
        for (var k = 0; k < parts.Count; k++)
        {
            parts[k].Id = k + 1;
            if (parts[k].Xyz.Length == 0)
            {
                var xs = new double[3 * parts[k].Nodes.Length];
                for (var q = 0; q < parts[k].Nodes.Length; q++)
                for (var c = 0; c < 3; c++)
                    xs[3 * q + c] = X[3 * parts[k].Nodes[q] + c];
                parts[k].Xyz = xs;
            }
        }

        var enc = new UTF8Encoding(false);
        using (var w = new StreamWriter(Path.Combine(dir, $"{prefix}.geo{tag}"), false, enc, 1 << 20))
        {
            w.WriteLine($"Corner crack, step {stepNo}");
            w.WriteLine("Erdogan-Sih growth on the maintained structures (Step2)");
            w.WriteLine("node id assign");
            w.WriteLine("element id assign");
            foreach (var part in parts)
            {
                var n = part.Xyz.Length / 3;
                w.WriteLine("part");
                w.WriteLine(part.Id.ToString(Ci).PadLeft(10));
                w.WriteLine(part.Name);
                w.WriteLine("coordinates");
                w.WriteLine(n.ToString(Ci).PadLeft(10));
                for (var c = 0; c < 3; c++)
                    for (var q = 0; q < n; q++)
                        w.WriteLine(F12(part.Xyz[3 * q + c]));
                w.WriteLine(part.Type);
                w.WriteLine((part.Conn.Length / part.NPerElem).ToString(Ci).PadLeft(10));
                var sb = new StringBuilder();
                for (var k = 0; k < part.Conn.Length / part.NPerElem; k++)
                {
                    sb.Clear();
                    for (var c = 0; c < part.NPerElem; c++) sb.Append(part.Conn[part.NPerElem * k + c].ToString(Ci).PadLeft(10));
                    w.WriteLine(sb.ToString());
                }
            }
        }

        // element-averaged von Mises at the nodes of the body
        double[] vmNode = new double[nNodes];
        if (withBody)
        {
            var cnt = new int[nNodes];
            for (var e = 0; e < s.ElementCount; e++)
            {
                var row = s.Row(e);
                if (row.Length != 4) continue;
                foreach (var v in row)
                {
                    vmNode[v] += vmElem[e];
                    cnt[v]++;
                }
            }

            for (var v = 0; v < nNodes; v++)
                if (cnt[v] > 0)
                    vmNode[v] /= cnt[v];
        }

        void Variable(string file, string label, int ncomp, Func<EPart, int, int, double> val)
        {
            using var w = new StreamWriter(Path.Combine(dir, $"{prefix}.{file}{tag}"), false, enc, 1 << 20);
            w.WriteLine(label);
            foreach (var part in parts)
            {
                w.WriteLine("part");
                w.WriteLine(part.Id.ToString(Ci).PadLeft(10));
                w.WriteLine("coordinates");
                for (var c = 0; c < ncomp; c++)
                    for (var q = 0; q < part.Nodes.Length; q++)
                        w.WriteLine(F12(val(part, q, c)));
            }
        }

        Variable("umag", "umag", 1, (part, q, _) =>
            part.IsBody || part.IsCrack
                ? Math.Sqrt(u[3 * part.Nodes[q]] * u[3 * part.Nodes[q]] + u[3 * part.Nodes[q] + 1] * u[3 * part.Nodes[q] + 1] +
                            u[3 * part.Nodes[q] + 2] * u[3 * part.Nodes[q] + 2])
                : 0);
        Variable("open", "opening", 1, (part, q, _) => part.IsCrack ? opening[part.Nodes[q]] : 0);
        Variable("disp", "displacement", 3, (part, q, c) => part.IsBody || part.IsCrack ? u[3 * part.Nodes[q] + c] : 0);
        if (withBody) Variable("mises", "mises", 1, (part, q, _) => part.IsBody ? vmNode[part.Nodes[q]] : 0);

        // connected components of the upper crack face (one for a crack without spurious cuts)
        var par = new int[corder.Count + 1];
        for (var k = 0; k < par.Length; k++) par[k] = k;
        int Find(int a)
        {
            while (par[a] != a) a = par[a] = par[par[a]];
            return a;
        }

        for (var k = 0; k < tris.Count; k++)
            for (var c = 1; c < 3; c++)
            {
                int ra = Find(cconn[3 * k]), rb = Find(cconn[3 * k + c]);
                if (ra != rb) par[ra] = rb;
            }

        var roots = new HashSet<int>();
        for (var k = 0; k < tris.Count; k++) roots.Add(Find(cconn[3 * k]));
        return roots.Count;
    }

    private static void WriteEnsightCase(string dir, string prefix, bool withBody, int[] stepNos, double[] times)
    {
        using var w = new StreamWriter(Path.Combine(dir, $"{prefix}.case"));
        w.WriteLine("FORMAT");
        w.WriteLine("type: ensight gold");
        w.WriteLine();
        w.WriteLine("GEOMETRY");
        w.WriteLine($"model: 1 {prefix}.geo****");
        w.WriteLine();
        w.WriteLine("VARIABLE");
        w.WriteLine($"vector per node: 1 displacement {prefix}.disp****");
        w.WriteLine($"scalar per node: 1 umag {prefix}.umag****");
        w.WriteLine($"scalar per node: 1 opening {prefix}.open****");
        if (withBody) w.WriteLine($"scalar per node: 1 mises {prefix}.mises****");
        w.WriteLine();
        w.WriteLine("TIME");
        w.WriteLine("time set: 1");
        w.WriteLine($"number of steps: {stepNos.Length}");
        var inc = stepNos.Length > 1 ? stepNos[1] - stepNos[0] : 1;
        var uniform = inc > 0 && Enumerable.Range(0, stepNos.Length - 1).All(i => stepNos[i + 1] - stepNos[i] == inc);
        if (uniform)
        {
            w.WriteLine($"filename start number: {stepNos[0]}");
            w.WriteLine($"filename increment: {inc}");
        }
        else
        {
            w.WriteLine("filename numbers:");
            foreach (var n in stepNos) w.WriteLine(n);
        }

        w.WriteLine("time values:");
        foreach (var t in times) w.WriteLine(t.ToString("0.00000000e+00", Ci));
    }
}

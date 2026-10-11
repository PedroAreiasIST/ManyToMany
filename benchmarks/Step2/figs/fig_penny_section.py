#!/usr/bin/env python3
"""Planar section of the deformed cracked mesh of the inclined penny example (vector, Xfig).

    python3 fig_penny_section.py <dump dir> [stem] [beta_deg] [scale]

<dump dir>: nodes.txt, tets.txt, boundary.txt, u.txt (Step2.Bench penny ... dump=<dir>). Every tetrahedron is cut
by the plane y = Y0 = 0.37 hmin (off every lattice plane, so that no section polygon degenerates); the cut points move with the
displacement interpolated along the cut edges, scaled by <scale>. Each segment of the section is drawn once: the
traces of the elements thin and grey, the segments on boundary faces (the crack faces and the box) black. The crack
faces carry their own nodes after the split, so the opening and sliding show in the deformed section. Panel (a): the
whole section; panel (b): the neighbourhood of the crack tip at theta = 0 (the front point in the section, x > 0).
"""
import math
import os
import sys

import numpy as np

SAND = "/home/pedro/Desktop/NovosPapers/engoopt2026/sand_paper/benchmarks"
sys.path.insert(0, SAND)
import figxfig as FX  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..")
Y0 = 0.37 * 0.0125      # 0.37 hmin (h0/8): off every lattice plane of the bisected Kuhn mesh
EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def section(X, T4, U, B, scale):
    """Segments of the section: (point a, point b, on a boundary face); points keyed by tet edge (sorted nodes)."""
    d = X[:, 1] - Y0
    s = np.sign(d[T4])
    cut = T4[(s.min(axis=1) < 0) & (s.max(axis=1) > 0)]
    bfaces = {tuple(sorted(f)) for f in B[:, :3].tolist()}
    pts = {}

    def point(i, j):
        k = (min(i, j), max(i, j))
        if k not in pts:
            t = d[k[0]] / (d[k[0]] - d[k[1]])
            p = X[k[0]] + t * (X[k[1]] - X[k[0]])
            u = U[k[0]] + t * (U[k[1]] - U[k[0]])
            q = p + scale * u
            pts[k] = (q[0], q[2])
        return k

    segs = {}
    for tet in cut:
        ks = [point(tet[a], tet[b]) for a, b in EDGES if (d[tet[a]] < 0) != (d[tet[b]] < 0)]
        n = len(ks)
        if n == 3:
            ring = ks
        else:                                                  # quad: order around the centroid
            c = np.mean([pts[k] for k in ks], axis=0)
            ring = sorted(ks, key=lambda k: math.atan2(pts[k][1] - c[1], pts[k][0] - c[0]))
        for i in range(n):
            ka, kb = ring[i], ring[(i + 1) % n]
            face = tuple(sorted(set(ka) | set(kb)))
            if len(face) != 3:                                 # opposite edges of a quad: not a face segment
                continue
            key = (min(ka, kb), max(ka, kb))
            segs[key] = segs.get(key, False) or face in bfaces
    return pts, segs


def draw(fig, pts, segs, x0, y0, w, window, title, heavy_pt=2, grey=60, frame=True, mark=None):
    (xa, xb), (za, zb) = window
    s = w / (xb - xa)
    h = s * (zb - za)

    def P(p):
        return (x0 + s * (p[0] - xa), y0 + s * (zb - p[1]))

    def inside(p):
        return xa <= p[0] <= xb and za <= p[1] <= zb

    f = fig.f
    for (ka, kb), heavy in segs.items():
        pa, pb = pts[ka], pts[kb]
        if not (inside(pa) or inside(pb)):
            continue
        # clip to the window (Liang-Barsky)
        t0, t1 = 0.0, 1.0
        dx, dz = pb[0] - pa[0], pb[1] - pa[1]
        ok = True
        for p_, q_ in ((-dx, pa[0] - xa), (dx, xb - pa[0]), (-dz, pa[1] - za), (dz, zb - pa[1])):
            if p_ == 0:
                ok = ok and q_ >= 0
                continue
            r = q_ / p_
            if p_ < 0:
                t0 = max(t0, r)
            else:
                t1 = min(t1, r)
        if not ok or t0 > t1:
            continue
        a = (pa[0] + t0 * dx, pa[1] + t0 * dz)
        b = (pa[0] + t1 * dx, pa[1] + t1 * dz)
        f.polyline(np.rint([P(a), P(b)]), thickness=heavy_pt if heavy else 1, pen=0 if heavy else grey,
                   depth=40 if heavy else 60)
    if frame:
        f.polyline(np.rint([P((xa, za)), P((xb, za)), P((xb, zb)), P((xa, zb)), P((xa, za))]), thickness=4,
                   pen=0, depth=30)
    if mark is not None:                                       # the window of another panel
        (ma, mb), (na, nb) = mark
        f.polyline(np.rint([P((ma, na)), P((mb, na)), P((mb, nb)), P((ma, nb)), P((ma, na))]), thickness=4,
                   pen=0, depth=20)
    f.text(x0 + 0.5 * w, y0 - 0.12 * FX.CM, title, sub=1, size=FX.TITLE_PT, depth=10)
    fig.ymax = max(fig.ymax, y0 + h)


def main():
    d = sys.argv[1]
    stem = sys.argv[2] if len(sys.argv) > 2 else "penny_section"
    beta = math.radians(float(sys.argv[3])) if len(sys.argv) > 3 else math.radians(30)
    scale = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
    X = np.loadtxt(os.path.join(d, "nodes.txt"))
    T4 = np.loadtxt(os.path.join(d, "tets.txt"), dtype=np.int64)
    B = np.loadtxt(os.path.join(d, "boundary.txt"), dtype=np.int64)
    U = np.loadtxt(os.path.join(d, "u.txt"))
    pts, segs = section(X, T4, U, B, scale)
    print("  %d section points, %d segments (%d on boundary faces)" % (len(pts), len(segs), sum(segs.values())))
    fig = FX.Figure()
    fig.f.color(60, "#959595")
    gap = 0.5 * FX.CM
    w = (FX.TEXTW * FX.CM - gap) / 2
    line = 1.25 * FX.TITLE_PT / 72.0 * 1200.0
    y0 = 0.12 * FX.CM + line
    allp = np.array(list(pts.values()))
    m = 0.02
    win = ((allp[:, 0].min() - m, allp[:, 0].max() + m), (allp[:, 1].min() - m, allp[:, 1].max() + m))
    side = max(win[0][1] - win[0][0], win[1][1] - win[1][0])
    cx, cz = 0.5 * (win[0][0] + win[0][1]), 0.5 * (win[1][0] + win[1][1])
    # (b) the tip at theta = 0: undeformed (a cos beta, -a sin beta) moved by the displacement of the nearest node
    a = 0.2
    tip = np.array([a * math.cos(beta), 0.0, -a * math.sin(beta)])
    k = np.argmin(np.linalg.norm(X - tip, axis=1))
    tx, tz = tip[0] + scale * U[k, 0], tip[2] + scale * U[k, 2]
    half = 0.06
    zoom = ((tx - half, tx + half), (tz - half, tz + half))
    draw(fig, pts, segs, 0.0, y0, w, ((-0.32, 0.32), (-0.32, 0.32)),
         r"(a) section $y=0.37h_{\min}$ near the crack, displacements $\times %g$" % scale, heavy_pt=2, frame=True, mark=zoom)
    draw(fig, pts, segs, w + gap, y0, w, zoom, r"(b) the crack tip at $\theta=0$ (square in (a))", heavy_pt=3)
    # fig2dev sets 7.5 slw for thickness 1 and 15 (t - 1) above: traces 0.15 pt, frames (t = 4) 0.45 pt,
    # crack faces and box (t = 2) 0.9 pt
    widths = {"7.500": "2.500", "45.000": "7.500"}

    def ps_filter(ps):
        out = []
        for line in ps.split("\n"):
            if line.endswith(" slw"):
                v = line.split()[0]
                line = widths.get(v, v) + " slw"
            out.append(line)
        return "\n".join(out)

    fig.export(stem, xfig_dir=os.path.join(FIGS, "xfig"), pdf_dir=os.path.join(FIGS, "pdf"), ps_filter=ps_filter,
               compact=True)


if __name__ == "__main__":
    main()

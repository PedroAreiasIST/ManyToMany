#!/usr/bin/env python3
"""Vector figure of the inclusion interfaces of the concurrent-cutting example (Step2.Bench interfaces mode).

    python3 fig_interfaces.py <interfaces.surf> [stem]

Reuses the 3D renderer of the SAND paper (engoopt2026/sand_paper/benchmarks/make_figures_yago.py: flat Blinn-Phong
shading, z-buffer culling, same-grey patches, compact vector PDF). Panel (a): every interface in the box [-1,1]^3,
isometric, patches merged; panel (b): one interface with the edges of its triangulation (the template triangles).
Writes figs/xfig/<stem>.fig and figs/pdf/<stem>.pdf.
"""
import functools
import os
import sys

import numpy as np

SAND = "/home/pedro/Desktop/NovosPapers/engoopt2026/sand_paper/benchmarks"
sys.path.insert(0, SAND)
import make_figures_yago as YG  # noqa: E402
import figxfig as FX  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..")


def read(path):
    b = open(path, "rb").read()
    nn, nt = (int(x) for x in np.frombuffer(b, "<i4", 2, 0))
    X = np.frombuffer(b, "<f8", 3 * nn, 8).reshape(nn, 3).copy()
    T = np.frombuffer(b, "<i4", 3 * nt, 8 + 24 * nn).reshape(nt, 3).copy()
    S = np.frombuffer(b, "<i4", nt, 8 + 24 * nn + 12 * nt).copy()
    return X, T, S


def box_geo(lo, hi, view):
    box = np.array([(a, b, c) for a in (lo[0], hi[0]) for b in (lo[1], hi[1]) for c in (lo[2], hi[2])])
    return dict(view=np.asarray(view, float), box=box, edges=YG.EDGES, planes=(), right_handed=True)


def main():
    path = sys.argv[1]
    stem = sys.argv[2] if len(sys.argv) > 2 else "interfaces"
    X, T, S = read(path)
    view = np.array([1.0, -0.62, 0.70])                    # from the +x, -y side and above, off the diagonal
    fig = YG.new_figure()
    gap = 0.5 * FX.CM
    w = (FX.TEXTW * FX.CM - gap) / 2
    line = 1.25 * FX.TITLE_PT / 72.0 * 1200.0
    y0 = 0.12 * FX.CM + line
    # (a) all interfaces in the box, merged patches
    YG.MERGE = True
    h1, n1 = YG.panel(fig, X, T, 0.0, y0, w, r"(a) $%d$ interfaces in $[-1,1]^3$" % len(set(S.tolist())), row=0,
                      geo=box_geo((-1, -1, -1), (1, 1, 1), view))
    # (b) the interface of one sphere, with the edges of its triangles
    k = int(np.bincount(S).argmax())
    Tk = T[S == k]
    used = np.unique(Tk)
    c = X[used].mean(axis=0)
    r = np.abs(X[used] - c).max() * 1.05
    YG.MERGE = False
    h2, n2 = YG.panel(fig, X, Tk, w + gap, y0, w, r"(b) one interface: %d triangles" % len(Tk), row=1,
                      geo=box_geo(c - r, c + r, view) | dict(edges=[]))
    fig.export(stem, xfig_dir=os.path.join(FIGS, "xfig"), pdf_dir=os.path.join(FIGS, "pdf"),
               ps_filter=functools.partial(YG.thin_mesh, widths=[YG.SEAM_PT, 0.12]), compact=True)
    print("  %s: %d + %d triangles drawn" % (stem, n1, n2))


if __name__ == "__main__":
    main()

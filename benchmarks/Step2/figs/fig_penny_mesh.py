#!/usr/bin/env python3
"""Vector figure of a cracked mesh of the penny example (Step2.Bench penny ... dump=<dir>).

    python3 fig_penny_mesh.py <dump dir> [stem] [beta_deg]

<dump dir> holds nodes.txt (x y z), tets.txt (4 nodes) and boundary.txt (3 nodes + element). Panel (a): the
tetrahedra with centroid y < 0 (a crinkle cut through the body), their exposed faces with mesh edges, showing the
refinement towards the front; panel (b): the crack surface (boundary faces off the box, drawn once, the faces of the
lower side) with its triangulation. Renderer: the SAND paper's make_figures_yago.panel (vector, no raster).
"""
import functools
import math
import os
import sys

import numpy as np

SAND = "/home/pedro/Desktop/NovosPapers/engoopt2026/sand_paper/benchmarks"
sys.path.insert(0, SAND)
import make_figures_yago as YG  # noqa: E402
import figxfig as FX  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..")


def box_geo(lo, hi, view, edges=True, planes=()):
    box = np.array([(a, b, c) for a in (lo[0], hi[0]) for b in (lo[1], hi[1]) for c in (lo[2], hi[2])])
    return dict(view=np.asarray(view, float), box=box, edges=YG.EDGES if edges else [], planes=planes,
                right_handed=True)


def main():
    d = sys.argv[1]
    stem = sys.argv[2] if len(sys.argv) > 2 else "penny_mesh"
    beta = math.radians(float(sys.argv[3])) if len(sys.argv) > 3 else 0.0
    X = np.loadtxt(os.path.join(d, "nodes.txt"))
    T4 = np.loadtxt(os.path.join(d, "tets.txt"), dtype=np.int64)
    B = np.loadtxt(os.path.join(d, "boundary.txt"), dtype=np.int64)
    # (a) crinkle cut: tetrahedra with centroid y < 0, exposed faces oriented outward
    cen = X[T4].mean(axis=1)
    keep = T4[cen[:, 1] < 0]
    faces = np.concatenate([keep[:, [1, 2, 3]], keep[:, [0, 3, 2]], keep[:, [0, 1, 3]], keep[:, [0, 2, 1]]])
    opp = np.concatenate([keep[:, 0], keep[:, 1], keep[:, 2], keep[:, 3]])
    key = np.sort(faces, axis=1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    exposed = cnt[inv.ravel()] == 1
    F = faces[exposed]
    O = opp[exposed]
    a, b, c = X[F[:, 0]], X[F[:, 1]], X[F[:, 2]]
    n = np.cross(b - a, c - a)
    flip = np.einsum("ij,ij->i", n, X[O] - a) > 0           # normal towards the opposite node: reverse
    F[flip] = F[flip][:, [0, 2, 1]]
    # (b) crack faces: boundary triangles off the box; one side (the element below the plane)
    nrm = np.array([math.sin(beta), 0.0, math.cos(beta)])
    onbox = np.zeros(len(B), bool)
    for k in range(3):
        for v in (-1.0, 1.0):
            onbox |= np.all(np.abs(X[B[:, :3], k] - v) < 1e-9, axis=1)
    C = B[~onbox]
    below = (X[T4[C[:, 3]]].mean(axis=1) @ nrm) < 0
    Cup = C[~below][:, :3].copy()
    Cf = C[below][:, :3].copy()
    a, b, c = X[Cf[:, 0]], X[Cf[:, 1]], X[Cf[:, 2]]
    up = np.einsum("ij,j->i", np.cross(b - a, c - a), nrm) < 0
    Cf[up] = Cf[up][:, [0, 2, 1]]                            # normals towards +n (seen from above)
    # (c) both crack faces displaced by a multiple of the computed displacements
    U = np.loadtxt(os.path.join(d, "u.txt"))
    cn = np.unique(np.concatenate([Cup.ravel(), Cf.ravel()]))
    um = U[cn].mean(axis=0)                                  # a translation does not show
    amp = np.linalg.norm(U[cn] - um, axis=1).max()
    scale = 0.25 * 0.2 / amp                                 # largest relative displacement: a quarter of the radius
    Xd = X + scale * (U - um)
    a, b, c = Xd[Cup[:, 0]], Xd[Cup[:, 1]], Xd[Cup[:, 2]]
    dn = np.einsum("ij,j->i", np.cross(b - a, c - a), nrm) < 0
    Cup[dn] = Cup[dn][:, [0, 2, 1]]                          # the upper face's normals towards +n
    a, b, c = Xd[Cf[:, 0]], Xd[Cf[:, 1]], Xd[Cf[:, 2]]
    Clo = Cf.copy()
    up2 = np.einsum("ij,j->i", np.cross(b - a, c - a), nrm) > 0
    Clo[up2] = Clo[up2][:, [0, 2, 1]]                        # the lower face's normals towards -n
    fig = YG.new_figure()
    gap = 0.5 * FX.CM
    w = (FX.TEXTW * FX.CM - gap) / 2
    line = 1.25 * FX.TITLE_PT / 72.0 * 1200.0
    y0 = 0.12 * FX.CM + line
    YG.MERGE = False
    view = np.array([0.55, 1.0, 0.62])                     # from the +y side: the cut faces y ~ 0 face the viewer
    planes = ((0, -1.0), (0, 1.0), (1, -1.0), (2, -1.0), (2, 1.0))     # the box faces keep their own normal
    YG.panel(fig, X, F, 0.0, y0, w, r"(a) tetrahedra with centroid $y<0$", row=0,
             geo=box_geo((-1, -1, -1), (1, 1, 1), view, planes=planes))
    used = np.unique(Cf)
    r = np.abs(X[used]).max() * 1.08
    YG.panel(fig, X, Cf, w + gap, y0, w, r"(b) crack surface: %d faces" % len(Cf), row=1,
             geo=box_geo((-r, -r, -r), (r, r, r), nrm, edges=False))
    fig.export(stem, xfig_dir=os.path.join(FIGS, "xfig"), pdf_dir=os.path.join(FIGS, "pdf"),
               ps_filter=functools.partial(YG.thin_mesh, widths=[0.08, 0.15]), compact=True)


if __name__ == "__main__":
    main()

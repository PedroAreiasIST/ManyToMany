#!/usr/bin/env python3
"""Colour figure of the semi-elliptical surface crack through the headless IST Vector Postprocessor.

    python3 fig_surface_color.py <dump dir> [stem] [scale] [t] [a] [h]

<dump dir>: nodes.txt, tets.txt, boundary.txt, u.txt (Step2.Bench surface ... dump=<dir>). The crack faces are the
boundary faces off the faces of the plate [-1,1] x [-t,0] x [-h,h]; the upper face (tetrahedra above the crack plane z = 0) is written as EnSight Gold
(figs/ensight/<stem>.case, ist_ensight_gold) with the deformed coordinates (displacements scaled by <scale>) and the
node field "opening", the normal jump [u].n against the partner node across the crack (the other node at the same
undeformed position; front nodes are not split and have none), divided by 8 (1 - nu^2) sigma a / (pi E) as in the
plot of the jump. ISTVectorPostprocessor --headless frames it from its best view and writes the vector PDF
figs/pdf/<stem>.pdf (cropped). Planar sections are drawn in xfig instead (fig_penny_section.py): the headless
export chooses its own view, which is not normal to a plane.
"""
import os
import subprocess
import sys

import numpy as np
from scipy.spatial import cKDTree

IVP = os.path.expanduser("~/RiderProjects/Numerical/ISTVectorPostprocessor")
sys.path.insert(0, IVP)
from ist_ensight_gold import write_ensight_gold  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..")
EXE = os.path.join(IVP, "csharp/ISTVectorPostprocessor.Desktop/bin/Release/net10.0/linux-x64/ISTVectorPostprocessor")
SIGMA, E, NU = 1.0, 1.0, 0.3


def crack_mask(X, B, t, h):
    """Boundary faces off the faces of the plate."""
    F = B[:, :3]
    onbox = np.zeros(len(F), bool)
    for ax, vals in ((0, (-1.0, 1.0)), (1, (-t, 0.0)), (2, (-h, h))):
        for v in vals:
            onbox |= np.all(np.abs(X[F][:, :, ax] - v) < 1e-12, axis=1)
    return ~onbox


def main():
    d = sys.argv[1]
    stem = sys.argv[2] if len(sys.argv) > 2 else "surface_opening"
    scale = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    t = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    A = float(sys.argv[5]) if len(sys.argv) > 5 else 0.15
    h = float(sys.argv[6]) if len(sys.argv) > 6 else 2.0
    X = np.loadtxt(os.path.join(d, "nodes.txt"))
    T4 = np.loadtxt(os.path.join(d, "tets.txt"), dtype=np.int64)
    B = np.loadtxt(os.path.join(d, "boundary.txt"), dtype=np.int64)
    U = np.loadtxt(os.path.join(d, "u.txt"))
    n = np.array([0.0, 0.0, 1.0])

    crack = B[crack_mask(X, B, t, h)]
    allnodes = np.unique(crack[:, :3])
    upper = X[T4[crack[:, 3]]].mean(axis=1) @ n > 0.37 * 0.1 / 6
    F = crack[upper, :3]
    nodes = np.unique(F)
    tree = cKDTree(X[allnodes])
    pairs = tree.query_ball_point(X[nodes], r=1e-12)
    partner = np.array([next((allnodes[j] for j in p if allnodes[j] != g), g) for g, p in zip(nodes, pairs)])
    ref = 8 * (1 - NU ** 2) * SIGMA * A / (np.pi * E)
    opening = np.abs((U[nodes] - U[partner]) @ n) / ref
    print("  upper crack face: %d triangles, %d nodes (%d front nodes); opening / ref: max %.4f (closed form %.4f)"
          % (len(F), len(nodes), int(np.sum(partner == nodes)), opening.max(), 0.0))

    local = np.full(len(X), -1, dtype=np.int64)
    local[nodes] = np.arange(len(nodes))
    out = os.path.join(FIGS, "ensight")
    os.makedirs(out, exist_ok=True)
    case = os.path.join(out, stem + ".case")
    write_ensight_gold(case, nodes=X[nodes] + scale * U[nodes], elements={"tria3": local[F]},
                       node_data={"opening": opening}, part_name="crack", binary=True)
    pdfdir = os.path.join(FIGS, "pdf")
    r = subprocess.run([EXE, "--headless", "-o", pdfdir, "-n", stem, "--views", "best", "--size", "520x300", "--field", "opening", case],
                       capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        sys.exit("headless export failed (%d): %s" % (r.returncode, r.stderr.strip()[-2000:]))
    raw = os.path.join(pdfdir, stem + "_best.pdf")
    subprocess.run(["pdfcrop", "--margins", "2", raw, os.path.join(pdfdir, stem + ".pdf")], check=True,
                   stdout=subprocess.DEVNULL)
    os.remove(raw)
    print("  wrote", os.path.join(pdfdir, stem + ".pdf"))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Write the plotgnu .data specs for the MM3 paper figures 8-13.

Stems keep the names the manuscript already \includegraphics, so the regenerated
plots drop straight in:
    assemblenew  assembleores  orderingsnew
    operationsnew  operationsnewcsr  operationsnewcores  cliquesnew
"""
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUT, exist_ok=True)


def write(stem, title, xlabel, ylabel, curves, *, legend="bottom right",
          xrange=("*", "*"), yrange=("*", "*"), ilog=0, xshift=0.0):
    """curves: list of (source, label, type, xcol, ycol, xmult, ymult)."""
    lines = [title, xlabel, ylabel, legend,
             f"{xrange[0]} {xrange[1]}", f"{yrange[0]} {yrange[1]}",
             f"{ilog} {xshift}"]
    for src, label, typ, xc, yc, xm, ym in curves:
        lines += [src, label, f"{typ} {xc} {yc} {xm} {ym}"]
    path = os.path.join(OUT, stem + ".data")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", path)


NEL = r"$n_{\mathrm{el}}$ [$\times 10^6$]"
NTH = r"$n_{\mathrm{threads}}$"
TMS = r"wall clock time [ms]"
LP = 2      # linespoints

# ---- Fig 8: assembling vs element count -------------------------------------
write("assemblenew",
      r"Assembling: flexible vs.\ CSR representation",
      NEL, TMS,
      [("fig08_assembly.txt", r"\texttt{List<List<int>>}", LP, 1, 2, 1.0, 1.0),
       ("fig08_assembly.txt", r"CSR", LP, 1, 3, 1.0, 1.0)],
      legend="top left")

# ---- Fig 9: assembling vs thread count --------------------------------------
write("assembleores",
      r"Bulk assembling ($4\times10^{6}$ elements): effect of threads",
      NTH, TMS,
      [("fig09_assembly_threads.txt", r"With symmetry group ($D_4$)", LP, 1, 2, 1.0, 1.0),
       ("fig09_assembly_threads.txt", r"Without symmetry group", LP, 1, 3, 1.0, 1.0)],
      legend="top right")

# ---- Fig 10: orderings ------------------------------------------------------
write("orderingsnew",
      r"Lexicographic ordering and compression/renumbering",
      NEL, TMS,
      [("fig10_orderings.txt", r"Lexicographic ordering", LP, 1, 2, 1.0, 1.0),
       ("fig10_orderings.txt", r"Compression/renumbering", LP, 1, 3, 1.0, 1.0)],
      legend="top left")

# ---- Fig 11a/b: set operations ---------------------------------------------
for stem, src, what in (("operationsnew", "fig11a_operations_o2m.txt", "O2M"),
                        ("operationsnewcsr", "fig11b_operations_csr.txt", "CSR")):
    write(stem,
          rf"Set operations on the {what} representation",
          NEL, TMS,
          [(src, r"Union ($+$)", LP, 1, 2, 1.0, 1.0),
           (src, r"Difference ($-$)", LP, 1, 3, 1.0, 1.0),
           (src, r"Intersection ($\&$)", LP, 1, 4, 1.0, 1.0)],
          legend="top left")

# ---- Fig 12: set operations vs thread count ---------------------------------
write("operationsnewcores",
      r"Set operations: effect of the number of threads",
      NTH, TMS,
      [("fig12_operations_threads.txt", r"Union, $4\times10^{6}$", LP, 1, 2, 1.0, 1.0),
       ("fig12_operations_threads.txt", r"Difference, $4\times10^{6}$", LP, 1, 3, 1.0, 1.0),
       ("fig12_operations_threads.txt", r"Intersection, $4\times10^{6}$", LP, 1, 4, 1.0, 1.0),
       ("fig12_operations_threads.txt", r"Union, $8\times10^{6}$", LP, 1, 5, 1.0, 1.0),
       ("fig12_operations_threads.txt", r"Difference, $8\times10^{6}$", LP, 1, 6, 1.0, 1.0),
       ("fig12_operations_threads.txt", r"Intersection, $8\times10^{6}$", LP, 1, 7, 1.0, 1.0)],
      legend="top right")

# ---- Fig 13: local positions and cliques ------------------------------------
write("cliquesnew",
      r"Local position searching and clique formation",
      NEL, TMS,
      [("fig13_cliques.txt", r"Element local positions", LP, 1, 2, 1.0, 1.0),
       ("fig13_cliques.txt", r"Node local positions", LP, 1, 3, 1.0, 1.0),
       ("fig13_cliques.txt", r"Clique formation", LP, 1, 4, 1.0, 1.0)],
      legend="top left")

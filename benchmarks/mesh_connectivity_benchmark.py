#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mesh_connectivity_benchmark.py
==============================

Reproducible benchmark of mesh *connectivity-query* performance across mesh
data-structure libraries, emitted in the same column layout as the comparison
table in Section 9 of the MM3 paper:

    Method | <mesh 1> | <mesh 2> | <mesh 3> | Reduction vs OpenMesh

--------------------------------------------------------------------------------
HONEST SCOPE -- read this before quoting any number it prints:
  * This measures performance ON THE MACHINE THAT RUNS IT. Results depend on the
    CPU, memory bandwidth, core/thread count, library versions, and whether the
    libraries use multiple threads internally. The numbers WILL NOT match numbers
    measured on different hardware, and they are NOT meant to reproduce a specific
    published table. This is a *reproducible harness*, not a way to regenerate
    someone else's measurements.
  * Different libraries use different connectivity models. OpenMesh / CGAL / VTK
    answer neighbour queries by per-element traversal of an explicit topology;
    libigl and MM3 build full adjacency arrays once (amortised). The script
    therefore times "build the topology + answer every neighbour query once" as a
    single comparable number, and also reports build vs query separately in the
    CSV. Treat cross-model comparisons as order-of-magnitude, not exact.
  * Surface half-edge libraries (OpenMesh, CGAL) do not represent tetrahedra, so
    they are reported as N/A on the volumetric (tet) mesh. That is expected.
--------------------------------------------------------------------------------

Libraries benchmarked if importable (each optional; a missing one => "N/A"):
    OpenMesh   ->  pip install openmesh          (half-edge, triangle surfaces)
    libigl     ->  pip install libigl            (adjacency arrays; `import igl`)
    VTK        ->  pip install vtk               (vtkPolyData / vtkUnstructuredGrid)
    CGAL       ->  pip install cgal              (HalfedgeDS; bindings are fragile)
    MM3 (C#)   ->  --mm3-project <dir>           (shells out to `dotnet run`; see
                                                  the contract in bench_mm3() below)

Meshes (default = synthetic, generated in-process so it runs OFFLINE anywhere):
    triangle rows : a perturbed, fully-triangulated grid at an EXACT target face
                    count (a valid 2-manifold-with-boundary surface mesh).
    tet row       : a structured box split into 6 tetrahedra per hex cell
                    (the same decomposition MM3's CreateBoxMesh uses).
    --download    : use the real Stanford Bunny/Dragon instead (needs network and
                    `pip install trimesh`); row labels then reflect the real sizes.

Usage
-----
    python mesh_connectivity_benchmark.py                 # synthetic, paper-ish sizes
    python mesh_connectivity_benchmark.py --quick         # tiny sizes (smoke test)
    python mesh_connectivity_benchmark.py --repeats 5     # median of 5 timed runs
    python mesh_connectivity_benchmark.py --download      # real Bunny/Dragon
    python mesh_connectivity_benchmark.py --mm3-project ./Mm3Bench
    python mesh_connectivity_benchmark.py --out results   # writes results.csv + results.md + results.tex

No third-party package is required to RUN the script; only NumPy is needed for the
synthetic mesh generators, and even that is optional (a pure-Python fallback is
used if NumPy is absent). Each benchmarked library is independently optional.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import math
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from typing import Callable, Optional

# --------------------------------------------------------------------------- #
# NumPy is used where available; otherwise a minimal fallback keeps us running.
# --------------------------------------------------------------------------- #
try:
    import numpy as np
    HAVE_NUMPY = True
except Exception:  # pragma: no cover
    np = None
    HAVE_NUMPY = False


# =========================================================================== #
#  Mesh container + synthetic generators
# =========================================================================== #
@dataclass
class Mesh:
    label: str
    kind: str                 # "tri" or "tet"
    V: "list|np.ndarray"      # (n, 3) float vertex coordinates
    C: "list|np.ndarray"      # (m, 3) tri indices or (m, 4) tet indices
    n_vertices: int
    n_cells: int


def _grid_triangles(target_faces: int):
    """A perturbed triangulated grid with ~target_faces triangles.

    An (N x N) grid of quads split along a diagonal gives 2*(N-1)^2 triangles,
    so we pick N to hit the target. A small pseudo-random z perturbation keeps
    the surface non-degenerate without depending on any RNG library.
    """
    n = max(2, int(round(math.sqrt(target_faces / 2.0))) + 1)
    if HAVE_NUMPY:
        xs = np.linspace(0.0, 1.0, n)
        gx, gy = np.meshgrid(xs, xs)
        gx = gx.ravel(); gy = gy.ravel()
        # deterministic, cheap perturbation (no RNG dependency)
        gz = 0.05 * np.sin(12.9898 * gx) * np.cos(78.233 * gy)
        V = np.column_stack([gx, gy, gz]).astype(np.float64)
        idx = np.arange(n * n).reshape(n, n)
        v00 = idx[:-1, :-1].ravel(); v10 = idx[:-1, 1:].ravel()
        v01 = idx[1:, :-1].ravel();  v11 = idx[1:, 1:].ravel()
        tri_a = np.column_stack([v00, v10, v11])
        tri_b = np.column_stack([v00, v11, v01])
        F = np.vstack([tri_a, tri_b]).astype(np.int64)
    else:  # pragma: no cover - pure-Python fallback
        V = []
        for j in range(n):
            for i in range(n):
                x = i / (n - 1); y = j / (n - 1)
                z = 0.05 * math.sin(12.9898 * x) * math.cos(78.233 * y)
                V.append((x, y, z))
        F = []
        def vid(i, j): return j * n + i
        for j in range(n - 1):
            for i in range(n - 1):
                a, b, c, d = vid(i, j), vid(i + 1, j), vid(i, j + 1), vid(i + 1, j + 1)
                F.append((a, b, d)); F.append((a, d, c))
    return Mesh(f"grid-{len(F)//1000}K-tri", "tri", V, F, len(V), len(F))


def _box_tets(target_tets: int):
    """Structured box meshed as 6 tetrahedra per hex cell (~target_tets tets)."""
    cells = max(1, target_tets // 6)
    n = max(1, int(round(cells ** (1.0 / 3.0))))
    nx = ny = n
    nz = max(1, cells // (nx * ny))
    NX, NY, NZ = nx + 1, ny + 1, nz + 1

    def vid(i, j, k): return (k * NY + j) * NX + i

    if HAVE_NUMPY:
        ii, jj, kk = np.meshgrid(np.arange(NX), np.arange(NY), np.arange(NZ), indexing="ij")
        V = np.column_stack([ii.ravel(order="F") / nx,
                             jj.ravel(order="F") / ny,
                             kk.ravel(order="F") / nz]).astype(np.float64)
        # rebuild vid consistent with the ravel order used above
        lin = np.arange(NX * NY * NZ).reshape(NX, NY, NZ, order="F")
        def vid(i, j, k): return lin[i, j, k]
    else:  # pragma: no cover
        V = [(i / nx, j / ny, k / nz)
             for k in range(NZ) for j in range(NY) for i in range(NX)]

    # 6-tet decomposition of each hex sharing the main diagonal n000-n111
    hexsplit = ((0, 1, 3, 7), (0, 3, 2, 7), (0, 2, 6, 7),
                (0, 6, 4, 7), (0, 4, 5, 7), (0, 5, 1, 7))
    corners = lambda i, j, k: [vid(i, j, k), vid(i + 1, j, k), vid(i, j + 1, k),
                               vid(i + 1, j + 1, k), vid(i, j, k + 1), vid(i + 1, j, k + 1),
                               vid(i, j + 1, k + 1), vid(i + 1, j + 1, k + 1)]
    T = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = corners(i, j, k)
                for a, b, d, e in hexsplit:
                    T.append((int(c[a]), int(c[b]), int(c[d]), int(c[e])))
    if HAVE_NUMPY:
        T = np.asarray(T, dtype=np.int64)
    return Mesh(f"box-{len(T)//1000}K-tet", "tet", V, T, len(V), len(T))


def load_real_meshes(quick: bool):
    """Optional: real Stanford Bunny / Dragon via trimesh (needs network)."""
    import trimesh  # raises if missing
    out = []
    sources = {
        "bunny": "https://github.com/alecjacobson/common-3d-test-models/raw/master/data/stanford-bunny.obj",
        "dragon": "https://github.com/alecjacobson/common-3d-test-models/raw/master/data/dragon.obj",
    }
    for name, url in sources.items():
        m = trimesh.load_remote(url, force="mesh")
        V = np.asarray(m.vertices, dtype=np.float64)
        F = np.asarray(m.faces, dtype=np.int64)
        out.append(Mesh(f"{name}-{len(F)//1000}K-tri", "tri", V, F, len(V), len(F)))
    return out


# =========================================================================== #
#  Mesh exchange writers (for the MM3 / C# bridge)
# =========================================================================== #
def write_obj(mesh: Mesh, path: str):
    with open(path, "w") as fh:
        for v in mesh.V:
            fh.write(f"v {float(v[0])} {float(v[1])} {float(v[2])}\n")
        for c in mesh.C:
            fh.write("f " + " ".join(str(int(x) + 1) for x in c) + "\n")


def write_tetmesh(mesh: Mesh, path: str):
    with open(path, "w") as fh:
        fh.write(f"{mesh.n_vertices} {mesh.n_cells}\n")
        for v in mesh.V:
            fh.write(f"{float(v[0])} {float(v[1])} {float(v[2])}\n")
        for c in mesh.C:
            fh.write(" ".join(str(int(x)) for x in c) + "\n")


# =========================================================================== #
#  Per-library benchmarks.  Each returns (build_ms, query_ms) or None if N/A.
#  The "score" used in the table is build_ms + query_ms (see run()).
# =========================================================================== #
def _importable(modname: str) -> bool:
    try:
        importlib.import_module(modname)
        return True
    except Exception:
        return False


def bench_openmesh(mesh: Mesh):
    if mesh.kind != "tri":
        return None  # half-edge surface library: no tets
    import openmesh as om
    V, F = mesh.V, mesh.C
    t0 = time.perf_counter()
    m = om.TriMesh()
    vh = [m.add_vertex([float(p[0]), float(p[1]), float(p[2])]) for p in V]
    for f in F:
        m.add_face(vh[int(f[0])], vh[int(f[1])], vh[int(f[2])])
    build = time.perf_counter() - t0
    t0 = time.perf_counter()
    s = 0
    for fh in m.faces():                 # face -> edge-adjacent faces
        for ff in m.ff(fh):
            s += ff.idx()
    for v in m.vertices():               # vertex -> incident faces
        for vf in m.vf(v):
            s += vf.idx()
    query = time.perf_counter() - t0
    return build * 1e3, query * 1e3


def bench_cgal(mesh: Mesh):
    if mesh.kind != "tri":
        return None
    # CGAL's Python (SWIG) bindings vary a lot between installs; build a
    # Polyhedron_3 via an OFF stream when possible, else report N/A.
    try:
        from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
    except Exception:
        return None
    V, F = mesh.V, mesh.C
    off = [f"OFF\n{len(V)} {len(F)} 0"]
    for v in V:
        off.append(f"{float(v[0])} {float(v[1])} {float(v[2])}")
    for f in F:
        off.append("3 " + " ".join(str(int(x)) for x in f))
    off_str = "\n".join(off) + "\n"
    with tempfile.NamedTemporaryFile("w", suffix=".off", delete=False) as tf:
        tf.write(off_str); off_path = tf.name
    try:
        t0 = time.perf_counter()
        P = Polyhedron_3(off_path)
        build = time.perf_counter() - t0
        t0 = time.perf_counter()
        s = 0
        for facet in P.facets():                 # facet -> adjacent facets
            h = facet.halfedge(); start = h
            while True:
                s += 1 if not h.opposite().is_border() else 0
                h = h.next()
                if h == start:
                    break
        query = time.perf_counter() - t0
        return build * 1e3, query * 1e3
    finally:
        try:
            os.unlink(off_path)
        except OSError:
            pass


def bench_vtk(mesh: Mesh):
    import vtk
    from vtk.util import numpy_support as nps
    V = mesh.V if HAVE_NUMPY else np.array(mesh.V)  # vtk path needs numpy
    C = mesh.C if HAVE_NUMPY else np.array(mesh.C)
    k = C.shape[1]
    t0 = time.perf_counter()
    pts = vtk.vtkPoints()
    pts.SetData(nps.numpy_to_vtk(np.ascontiguousarray(V), deep=1))
    # fast cell-array construction (VTK 9): offsets + connectivity
    n = C.shape[0]
    offsets = np.arange(0, (n + 1) * k, k, dtype=np.int64)
    conn = np.ascontiguousarray(C.ravel().astype(np.int64))
    ca = vtk.vtkCellArray()
    ca.SetData(nps.numpy_to_vtkIdTypeArray(offsets, deep=1),
               nps.numpy_to_vtkIdTypeArray(conn, deep=1))
    if mesh.kind == "tri":
        grid = vtk.vtkPolyData(); grid.SetPoints(pts); grid.SetPolys(ca)
    else:
        grid = vtk.vtkUnstructuredGrid(); grid.SetPoints(pts)
        grid.SetCells(vtk.VTK_TETRA, ca)
    grid.BuildLinks()
    build = time.perf_counter() - t0

    t0 = time.perf_counter()
    s = 0
    idl = vtk.vtkIdList(); nb = vtk.vtkIdList(); edge = vtk.vtkIdList()
    npts = grid.GetNumberOfPoints()
    for i in range(npts):                          # vertex -> incident cells
        grid.GetPointCells(i, idl); s += idl.GetNumberOfIds()
    ncells = grid.GetNumberOfCells()
    if mesh.kind == "tri":
        faces_of = [(0, 1), (1, 2), (2, 0)]        # edges of a triangle
        for cid in range(ncells):
            cell = grid.GetCell(cid)
            for a, b in faces_of:
                edge.Reset(); edge.InsertNextId(cell.GetPointId(a)); edge.InsertNextId(cell.GetPointId(b))
                grid.GetCellNeighbors(cid, edge, nb); s += nb.GetNumberOfIds()
    else:
        faces_of = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]  # faces of a tet
        for cid in range(ncells):
            cell = grid.GetCell(cid)
            for a, b, c in faces_of:
                edge.Reset()
                edge.InsertNextId(cell.GetPointId(a)); edge.InsertNextId(cell.GetPointId(b)); edge.InsertNextId(cell.GetPointId(c))
                grid.GetCellNeighbors(cid, edge, nb); s += nb.GetNumberOfIds()
    query = time.perf_counter() - t0
    return build * 1e3, query * 1e3


def bench_igl(mesh: Mesh):
    import igl
    F = mesh.C if HAVE_NUMPY else np.array(mesh.C)
    nV = mesh.n_vertices
    if mesh.kind == "tri":
        t0 = time.perf_counter()
        _ = igl.adjacency_list(F)                       # vertex-vertex
        build = time.perf_counter() - t0
        t0 = time.perf_counter()
        _tt, _tti = igl.triangle_triangle_adjacency(F)  # face-face
        _vf, _ni = igl.vertex_triangle_adjacency(F, nV) # vertex-face
        query = time.perf_counter() - t0
        return build * 1e3, query * 1e3
    else:
        # tet adjacency is not exposed in every libigl build
        fn = getattr(igl, "tet_tet_adjacency", None)
        if fn is None:
            return None
        t0 = time.perf_counter()
        out = fn(F)
        query = time.perf_counter() - t0
        return 0.0, query * 1e3


def bench_mm3(mesh: Mesh, project: Optional[str], repeats: int = 5):
    """Shell out to a .NET benchmark. Contract:

        dotnet run -c Release --project <project> -- <meshfile> <tri|tet> <repeats>

    must print a final line:   BUILD_MS QUERY_MS   (two floats, milliseconds).

    The bridge runs `repeats` timed iterations (after a discarded warm-up) inside a
    single warm process and prints the MEDIAN, so MM3 is measured in the same
    steady-state regime as the in-process Python libraries and pays .NET process
    start / JIT only once for the whole cell. Without --mm3-project, MM3 is N/A.
    """
    if not project:
        return None
    suffix = ".obj" if mesh.kind == "tri" else ".tetmesh"
    with tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False) as tf:
        path = tf.name
    try:
        (write_obj if mesh.kind == "tri" else write_tetmesh)(mesh, path)
        cmd = ["dotnet", "run", "-c", "Release", "--project", project, "--",
               path, mesh.kind, str(repeats)]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if res.returncode != 0:
            sys.stderr.write(f"[MM3] dotnet failed:\n{res.stderr[-2000:]}\n")
            return None
        last = [ln for ln in res.stdout.strip().splitlines() if ln.strip()][-1]
        b, q = last.split()
        return float(b), float(q)
    except Exception as exc:  # pragma: no cover
        sys.stderr.write(f"[MM3] {exc}\n")
        return None
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


METHODS = [
    ("OpenMesh (Half-Edge)", "openmesh", bench_openmesh),
    ("CGAL (HalfedgeDS)",    "CGAL",     bench_cgal),
    ("VTK",                  "vtk",      bench_vtk),
    ("libigl",               "igl",      bench_igl),
    # MM3 handled separately (external process)
]


# =========================================================================== #
#  Runner + reporting
# =========================================================================== #
def host_info() -> dict:
    info = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "processor": platform.processor() or "unknown",
        "logical_cpus": os.cpu_count() or 0,
    }
    try:
        import psutil
        info["ram_gb"] = round(psutil.virtual_memory().total / 1e9, 1)
        info["physical_cpus"] = psutil.cpu_count(logical=False) or info["logical_cpus"]
    except Exception:
        info["ram_gb"] = None
        info["physical_cpus"] = None
    for mod in ("numpy", "openmesh", "vtk", "igl", "CGAL", "trimesh"):
        try:
            m = importlib.import_module(mod)
            info[f"ver_{mod}"] = getattr(m, "__version__", "present")
        except Exception:
            info[f"ver_{mod}"] = "absent"
    return info


def median_run(fn: Callable[[Mesh], Optional[tuple]], mesh: Mesh, repeats: int):
    builds, queries = [], []
    for _ in range(repeats):
        r = fn(mesh)
        if r is None:
            return None
        builds.append(r[0]); queries.append(r[1])
    return statistics.median(builds), statistics.median(queries)


def run(meshes, repeats: int, mm3_project: Optional[str]):
    rows = []  # (method_name, {mesh_label: (build,query) or None})
    method_list = list(METHODS) + [("MM3", "__mm3__", None)]
    for name, mod, fn in method_list:
        is_mm3 = mod == "__mm3__"
        available = (is_mm3 and mm3_project is not None) or \
                    (not is_mm3 and _importable(mod))
        per_mesh = {}
        for mesh in meshes:
            if not available:
                per_mesh[mesh.label] = None
                continue
            sys.stderr.write(f"  {name:22s} on {mesh.label:18s} ... ")
            sys.stderr.flush()
            try:
                if is_mm3:
                    # The .NET bridge does its own warm-up + `repeats` timed runs in
                    # one process and returns the median, so MM3 is measured in the
                    # same steady-state regime as the in-process libraries above and
                    # pays process start / JIT only once per cell.
                    per_mesh[mesh.label] = bench_mm3(mesh, mm3_project, repeats)
                else:
                    per_mesh[mesh.label] = median_run(fn, mesh, repeats)
                sys.stderr.write("ok\n" if per_mesh[mesh.label] else "N/A\n")
            except Exception as exc:
                sys.stderr.write(f"error ({type(exc).__name__})\n")
                per_mesh[mesh.label] = None
        rows.append((name, per_mesh))
    return rows


def _score(cell):  # total time used for the table & reduction
    return None if cell is None else (cell[0] + cell[1])


def format_table(rows, meshes) -> str:
    labels = [m.label for m in meshes]
    ref = dict((m.label, None) for m in meshes)
    for name, per in rows:
        if name.startswith("OpenMesh"):
            ref = {lab: _score(per[lab]) for lab in labels}
    headers = ["Method"] + labels + ["Reduction vs OpenMesh"]
    table = [headers]
    for name, per in rows:
        line = [name]
        reds = []
        for lab in labels:
            sc = _score(per[lab])
            line.append("N/A" if sc is None else f"{sc:.1f}")
            r = ref.get(lab)
            if sc is not None and r not in (None, 0) and not name.startswith("OpenMesh"):
                reds.append(1.0 - sc / r)
        if name.startswith("OpenMesh"):
            line.append("--")
        elif reds:
            line.append(f"{100.0 * statistics.mean(reds):.1f}%")
        else:
            line.append("N/A")
        table.append(line)
    widths = [max(len(r[i]) for r in table) for i in range(len(headers))]
    out = []
    for ri, row in enumerate(table):
        out.append("  ".join(c.ljust(widths[i]) for i, c in enumerate(row)))
        if ri == 0:
            out.append("  ".join("-" * widths[i] for i in range(len(headers))))
    return "\n".join(out)


def _tex_escape(s: str) -> str:
    """Escape the LaTeX special characters that can appear in a table cell.

    The important one here is '%' (reduction percentages such as "33.8%"),
    which would otherwise comment out the rest of the row.
    """
    repl = {
        "\\": r"\textbackslash{}",
        "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
        "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(repl.get(ch, ch) for ch in s)


def write_outputs(rows, meshes, info, stem: str):
    labels = [m.label for m in meshes]
    with open(stem + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["method"] + [f"{l}__build_ms" for l in labels]
                   + [f"{l}__query_ms" for l in labels]
                   + [f"{l}__total_ms" for l in labels])
        for name, per in rows:
            b = [("" if per[l] is None else f"{per[l][0]:.3f}") for l in labels]
            q = [("" if per[l] is None else f"{per[l][1]:.3f}") for l in labels]
            t = [("" if per[l] is None else f"{per[l][0] + per[l][1]:.3f}") for l in labels]
            w.writerow([name] + b + q + t)
    with open(stem + ".md", "w") as fh:
        fh.write("# Mesh connectivity benchmark\n\n")
        fh.write("Host: " + ", ".join(f"{k}={v}" for k, v in info.items()
                                       if k in ("platform", "processor",
                                                "physical_cpus", "logical_cpus",
                                                "ram_gb")) + "\n\n")
        fh.write("Times are **total (build + query) milliseconds**, "
                 "median over repeats, measured on the host above. "
                 "They are host-specific and not comparable across machines.\n\n")
        body = format_table(rows, meshes).splitlines()
        # turn the fixed-width table into a markdown table
        hdr = [c.strip() for c in body[0].split("  ") if c.strip()]
        fh.write("| " + " | ".join(hdr) + " |\n")
        fh.write("| " + " | ".join("---" for _ in hdr) + " |\n")
        for line in body[2:]:
            cells = [c.strip() for c in line.split("  ") if c.strip()]
            fh.write("| " + " | ".join(cells) + " |\n")
    with open(stem + ".tex", "w") as fh:
        host = ", ".join(f"{k}={v}" for k, v in info.items()
                         if k in ("platform", "processor", "physical_cpus",
                                  "logical_cpus", "ram_gb"))
        body = format_table(rows, meshes).splitlines()
        hdr = [c.strip() for c in body[0].split("  ") if c.strip()]
        colspec = "l" + "r" * (len(hdr) - 1)   # method left, numbers right
        # A bare tabular meant to be \input{} into a table float in the paper.
        fh.write("% Mesh connectivity benchmark: total (build + query) ms, "
                 "median over repeats.\n")
        fh.write(f"% Host: {host}\n")
        fh.write("% Host-specific; not comparable across machines.\n")
        fh.write("\\begin{tabular}{" + colspec + "}\n\\hline\n")
        fh.write(" & ".join(_tex_escape(h) for h in hdr) + " \\\\\n\\hline\n")
        for line in body[2:]:
            cells = [c.strip() for c in line.split("  ") if c.strip()]
            fh.write(" & ".join(_tex_escape(c) for c in cells) + " \\\\\n")
        fh.write("\\hline\n\\end{tabular}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true",
                    help="tiny meshes for a fast smoke test")
    ap.add_argument("--repeats", type=int, default=3, help="timed runs per cell (median)")
    ap.add_argument("--download", action="store_true",
                    help="use real Stanford Bunny/Dragon (needs network + trimesh)")
    ap.add_argument("--mm3-project", default=os.environ.get("MM3_PROJECT"),
                    help="path to a .NET project implementing the MM3 bench contract")
    ap.add_argument("--tri-sizes", type=int, nargs="+",
                    help="triangle face counts (overrides defaults)")
    ap.add_argument("--tet-size", type=int, help="tet count (overrides default)")
    ap.add_argument("--out", default=None,
                    help="write <out>.csv, <out>.md and <out>.tex")
    args = ap.parse_args(argv)

    if not HAVE_NUMPY:
        sys.stderr.write("WARNING: NumPy not found; using slow pure-Python mesh "
                         "generation and the VTK/libigl paths will be unavailable.\n")

    if args.quick:
        tri_sizes = args.tri_sizes or [2_000, 8_000]
        tet_size = args.tet_size or 6_000
    else:
        tri_sizes = args.tri_sizes or [100_000, 500_000]
        tet_size = args.tet_size or 1_200_000

    info = host_info()
    print("=" * 78)
    print("Mesh connectivity benchmark  (numbers are HOST-SPECIFIC; see header)")
    print("=" * 78)
    for k in ("platform", "processor", "physical_cpus", "logical_cpus", "ram_gb",
              "python", "ver_numpy", "ver_openmesh", "ver_vtk", "ver_igl",
              "ver_CGAL"):
        print(f"  {k:16s}: {info.get(k)}")
    print("-" * 78)

    if args.download:
        try:
            meshes = load_real_meshes(args.quick)
        except Exception as exc:
            sys.stderr.write(f"--download failed ({exc}); falling back to synthetic.\n")
            meshes = [_grid_triangles(s) for s in tri_sizes] + [_box_tets(tet_size)]
    else:
        meshes = [_grid_triangles(s) for s in tri_sizes] + [_box_tets(tet_size)]

    print("Meshes:")
    for m in meshes:
        print(f"  {m.label:18s}  V={m.n_vertices:>9,}  cells={m.n_cells:>9,}  ({m.kind})")
    print("-" * 78)
    print(f"Workload: build topology + traverse every cell's face-neighbours and "
          f"every vertex's incident cells once.\nTimed cells = median of "
          f"{args.repeats} run(s). MM3 {'ENABLED' if args.mm3_project else 'N/A (pass --mm3-project)'}.")
    print("-" * 78)

    rows = run(meshes, args.repeats, args.mm3_project)

    print("\nTotal time (build + query), milliseconds:\n")
    print(format_table(rows, meshes))
    print("\nNote: cross-library numbers mix traversal-based (OpenMesh/CGAL/VTK) and "
          "array-build (libigl/MM3) models; read as order-of-magnitude.")

    if args.out:
        write_outputs(rows, meshes, info, args.out)
        print(f"\nWrote {args.out}.csv, {args.out}.md and {args.out}.tex")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

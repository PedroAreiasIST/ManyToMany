#!/usr/bin/env python3
"""Plots of the three examples (house style: gnuplot fig terminal, LaTeX labels, figtex2pdf.sh -> figs/*.pdf).

    python3 make_examples.py [penny] [conv] [cutpar] [growth]      # default: all whose data exist

Data: benchmarks/data_final (amd2), benchmarks/data_deu (Deucalion), benchmarks/logs_final.
"""
import csv
import glob
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, "..", "..", "benchmarks")
NU = 0.3
COLORS = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7"]
HEAD = """set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "{stem}.fig"
""" + "".join('set linetype %d lw 2 pt %d ps 1.1 lc rgb "%s"\n' % (i + 1, pt, c)
              for i, (c, pt) in enumerate(zip(COLORS, [7, 5, 9, 13, 15, 11, 3])))


def render(stem, body):
    gp = os.path.join(HERE, stem + ".gp")
    with open(gp, "w") as f:
        f.write(HEAD.format(stem=stem) + body + "unset output\n")
    subprocess.run(["gnuplot", stem + ".gp"], cwd=HERE, check=True)
    subprocess.run(["figtex2pdf.sh", stem + ".fig"], cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    pdf = os.path.join(HERE, stem + ".pdf")
    out = os.path.join(HERE, "..", stem + ".pdf")
    os.replace(pdf, out)
    # the LyX manuscript includes the EPS (tex2lyx resolves figs/<stem> to it): keep both in step
    subprocess.run(["pdftops", "-eps", out, os.path.join(HERE, "..", stem + ".eps")], check=True)
    print("  figs/%s.pdf, .eps" % stem)


def find(pattern):
    for root in ("data_deu", "data_final"):
        hits = sorted(glob.glob(os.path.join(B, root, "**", pattern), recursive=True))
        if hits:
            return hits[-1]
    return None


def penny_profile():
    """K along the front of the inclined crack (beta = 30 deg, h0/16) against Kassir and Sih."""
    path = find("penny_20_16_30_rq0.075.csv")
    if not path:
        print("  penny profile: no data yet")
        return
    a, beta, sigma = 0.2, math.radians(30), 1.0
    k0 = 2 * sigma * math.sqrt(a / math.pi)
    rows = list(csv.DictReader(open(path)))
    # the computation measures theta from e1; the shear on the upper face acts along -e1, and Kassir and Sih
    # measure theta from the shear direction: theta_KS = theta + 180 deg. The computed values are plotted as they
    # are, against theta_KS.
    with open(os.path.join(HERE, "penny_K.txt"), "w") as f:
        for r in sorted(rows, key=lambda r: (float(r["theta_deg"]) + 180) % 360):
            f.write("%.4f %.6f %.6f %.6f\n" % ((float(r["theta_deg"]) + 180) % 360, float(r["KI"]) / k0,
                                              float(r["KII"]) / k0, float(r["KIII"]) / k0))
    sn, tau = sigma * math.cos(beta) ** 2, sigma * math.sin(beta) * math.cos(beta)
    with open(os.path.join(HERE, "penny_Kref.txt"), "w") as f:
        for i in range(361):
            th = math.radians(i)
            f.write("%d %.6f %.6f %.6f\n" % (i, 2 * sn * math.sqrt(a / math.pi) / k0,
                                             4 / (2 - NU) * tau * math.sqrt(a / math.pi) * math.cos(th) / k0,
                                             4 * (1 - NU) / (2 - NU) * tau * math.sqrt(a / math.pi) * math.sin(th) / k0))
    render("penny_K", """set xlabel "Angle from the shear direction $\\\\theta$ [deg]"
set ylabel "$K/K_0$"
set xrange [0:360]
set xtics 0,90,360
set yrange [-0.75:1.0]
set key outside right top Left reverse nobox spacing 1.4
plot "penny_Kref.txt" u 1:2 w l lt 1 dt 1 t "$K_{\\\\mathrm{I}}$, Kassir--Sih", \\
     "penny_K.txt" u 1:2 w p lt 1 t "$K_{\\\\mathrm{I}}$, computed", \\
     "penny_Kref.txt" u 1:3 w l lt 5 t "$K_{\\\\mathrm{II}}$, Kassir--Sih", \\
     "penny_K.txt" u 1:3 w p lt 5 t "$K_{\\\\mathrm{II}}$, computed", \\
     "penny_Kref.txt" u 1:4 w l lt 3 t "$K_{\\\\mathrm{III}}$, Kassir--Sih", \\
     "penny_K.txt" u 1:4 w p lt 3 t "$K_{\\\\mathrm{III}}$, computed"
""")


def penny_convergence():
    """Largest error of K along the front against hmin / a (rq = 0.075), from the run summaries."""
    pat = re.compile(r"penny n=20 div=(\d+) beta=(\d+): max \|K - K_ref\| / K0 = ([\d.]+) % \(I\), ([\d.]+) % "
                     r"\(\|II\|\), ([\d.]+) % \(\|III\|\); rq=0\.0750")
    res = {}
    for log in glob.glob(os.path.join(B, "logs_final", "penny_*.log")) + \
            glob.glob(os.path.join(B, "logs_deu", "*.out")):
        for line in open(log, errors="replace"):
            m = pat.search(line)
            if m:
                res[(int(m.group(1)), int(m.group(2)))] = tuple(float(m.group(k)) for k in (3, 4, 5))
    divs = sorted({d for d, _ in res})
    if len(divs) < 2:
        print("  convergence: not enough data")
        return
    with open(os.path.join(HERE, "penny_conv.txt"), "w") as f:
        for d in divs:
            h = 0.1 / d / 0.2
            e0 = res.get((d, 0), (math.nan,) * 3)
            e30 = res.get((d, 30), (math.nan,) * 3)
            f.write("%.6f %s %s %s %s\n" % (h, e0[0], e30[0], e30[1], e30[2]))
    render("penny_conv", """set logscale xy
set xlabel "$h_{\\\\min}/a$"
set ylabel "$\\\\max_\\\\theta |K-K_{\\\\mathrm{ref}}|/K_0$ [\\\\%]"
set xrange [0.025:0.15]
set yrange [0.1:5]
set key outside right top Left reverse nobox spacing 1.4
plot "penny_conv.txt" u 1:2 w lp lt 1 t "$\\\\beta=0$: $K_{\\\\mathrm{I}}$", \\
     "" u 1:3 w lp lt 2 t "$\\\\beta=30^\\\\circ$: $K_{\\\\mathrm{I}}$", \\
     "" u 1:4 w lp lt 5 t "$\\\\beta=30^\\\\circ$: $K_{\\\\mathrm{II}}$", \\
     "" u 1:5 w lp lt 3 t "$\\\\beta=30^\\\\circ$: $K_{\\\\mathrm{III}}$", \\
     12*x w l lt 7 dt 2 t "slope 1", 80*x**2 w l lt 7 dt 3 t "slope 2"
""")


def cutpar():
    """Speed-up of the concurrent cut against threads (sequential cut / concurrent cut, minima of the times without collector pauses)."""
    sets = [("data_fin4", "cutpar_119_4.csv", "9950X, $10^7$ el., 64 groups"),
            ("data_fin4", "cutpar_119_8.csv", "9950X, $10^7$ el., 512 groups"),
            ("data_deu7", "cutpar_119_4.csv", "EPYC 7742, $10^7$ el., 64 groups"),
            ("data_deu7", "cutpar_119_8.csv", "EPYC 7742, $10^7$ el., 512 groups"),
            ("data_deu7", "cutpar_160_8.csv", "EPYC 7742, $2.5\\\\times10^7$ el., 512 groups")]
    plots = []
    for i, (d, name, label) in enumerate(sets):
        path = os.path.join(B, d, name)
        if not os.path.exists(path):
            continue
        rows = list(csv.DictReader(open(path)))
        seq = float(next(r for r in rows if r["impl"] == "sequential")["nogc_ms_min"])
        out = "cutpar_%d.txt" % i
        with open(os.path.join(HERE, out), "w") as f:
            for r in rows:
                if r["impl"] == "parallel":
                    f.write("%s %.4f\n" % (r["threads"], seq / float(r["nogc_ms_min"])))
        plots.append('"%s" u 1:2 w lp lt %d t "%s"' % (out, i + 1, label))
    if not plots:
        print("  cutpar: no final-code data yet")
        return
    render("cutpar", """set logscale x 2
set xlabel "Threads"
set ylabel "Speed-up over the sequential cut"
set xrange [0.8:160]
set yrange [0:14]
set key inside left top Left reverse nobox spacing 1.3
plot """ + ", \\\n     ".join(plots) + "\n")


def growth():
    """Fronts of the growing elliptical crack, K_I along the front, and the cost of a step."""
    gdir = os.environ.get("GROWTH_DIR", os.path.join(B, "data_enh", "cut2", "growth"))
    path = os.path.join(gdir, "growth_fronts.csv")
    if not os.path.exists(path):
        print("  growth: no data yet")
        return
    d = os.path.dirname(path)
    rows = list(csv.DictReader(open(path)))
    steps = sorted({int(r["step"]) for r in rows})
    with open(os.path.join(HERE, "growth_fronts.txt"), "w") as f:
        for st in steps:
            pts = [(float(r["x"]), float(r["y"])) for r in rows if int(r["step"]) == st]
            for x, y in pts + pts[:1]:
                f.write("%.6f %.6f\n" % (x, y))
            f.write("\n\n")
    with open(os.path.join(HERE, "growth_aspect.txt"), "w") as f:
        for st in steps:
            rad = {round(float(r["theta_deg"])): float(r["R"]) for r in rows if int(r["step"]) == st}
            f.write("%d %.6f %.6f %.6f\n" % (st, rad[0], rad[90], rad[90] / rad[0]))
    render("growth_fronts", """set size ratio -1
set xlabel "$x$"
set ylabel "$y$"
set xtics 0.1
set ytics 0.1
set key off
plot for [i=1:%d] "growth_fronts.txt" index i u 1:2 w l lt 3 lw 1.5, "growth_fronts.txt" index 0 u 1:2 w l lt 1 lw 2.5
""" % (len(steps) - 1))
    # K_I along the front: first, middle and last step, Irwin at the first
    krows = list(csv.DictReader(open(os.path.join(d, "growth_K.csv"))))
    sig, b0 = 1.0, 0.15
    k0 = sig * math.sqrt(math.pi * b0)
    pick = sorted({steps[0], steps[len(steps) // 2], steps[-1]})
    plots = []
    for i, st in enumerate(pick):
        out = "growth_K_%d.txt" % st
        with open(os.path.join(HERE, out), "w") as f:
            for r in krows:
                if int(r["step"]) == st:
                    f.write("%s %.6f %s\n" % (r["theta_deg"], float(r["KI"]) / k0,
                                             "%.6f" % (float(r["KI_irwin"]) / k0) if r["KI_irwin"] not in ("NaN", "") else "nan"))
        plots.append('"%s" u 1:2 w lp lt %d t "step %d"' % (out, i + 2, st))
    plots.append('"growth_K_%d.txt" u 1:3 w l lt 1 t "Irwin, step %d"' % (pick[0], pick[0]))
    render("growth_K", """set xlabel "Angle along the front $\\\\psi$ [deg]"
set ylabel "$K_{\\\\mathrm{I}}/(\\\\sigma\\\\sqrt{\\\\pi b_0})$"
set xrange [0:360]
set xtics 0,90,360
set yrange [0.6:1.2]
set key inside top center horizontal maxcols 2 Left reverse nobox spacing 1.2 width -1
plot """ + ", \\\n     ".join(plots) + "\n")
    # cost of the accepted steps
    srows = [r for r in csv.DictReader(open(os.path.join(d, "growth_steps.csv")))]
    with open(os.path.join(HERE, "growth_cost.txt"), "w") as f:
        for r in srows:
            f.write("%s %s %.4f %s %s %s %s %s %s %s %s\n" % (r["attempt"], r["accepted"],
                    float(r["t_fit_ms"]) + float(r["t_insert_ms"]), r["t_publish_ms"],
                    r["t_rollback_ms"] if float(r["t_rollback_ms"]) > 0 else "nan", r["t_rebuild_ms"],
                    r["t_assemble_ms"], r["t_solve_ms"], r["t_K_ms"], r["elements"], r.get("t_copy_ms") or "nan"))
    has_rb = any(float(r["t_rollback_ms"]) > 0 for r in srows)
    has_copy = any(r.get("t_copy_ms") for r in srows)
    cp = '     "" u 1:11 w lp lt 5 pt 4 t "flat copy refresh", \\\n' if has_copy else ""
    rb = '     "" u 1:5 w p lt 5 pt 2 ps 1.6 t "rollback", \\\n' if has_rb else ""
    render("growth_cost", """set logscale y
set format y "$10^{%T}$"
set xlabel "Attempt"
set ylabel "Time [ms]"
set yrange [1e-3:1e7]
set key inside top center horizontal maxcols 2 Left reverse nobox spacing 1.2 width -2
plot "growth_cost.txt" u 1:3 w lp lt 1 t "fit and crack increment", \\
     "" u 1:4 w lp lt 2 t "commit", \\
""" + rb + cp + """     "" u 1:6 w lp lt 3 t "recomputation of every structure", \\
     "" u 1:7 w lp lt 4 t "assembly", \\
     "" u 1:8 w lp lt 6 t "solution", \\
     "" u 1:9 w lp lt 7 t "$K$ at 32 points"
""")


def penny_cod():
    """Crack-face displacement jump of the inclined crack against the closed forms (cod.py on a dump)."""
    dump = None
    for cand in (os.path.join(B, "data_final", "dump", "penny_20_8_30"), os.path.join(B, "simplas_check", "mesh")):
        if os.path.exists(os.path.join(cand, "u.txt")):
            dump = cand
            break
    if dump is None:
        print("  cod: no dump")
        return
    subprocess.run([sys.executable, os.path.join(HERE, "cod.py"), dump, "30", os.path.join(HERE, "penny_cod.txt")],
                   check=True)
    b = math.radians(30)
    render("penny_cod", """set xlabel "$r/a$"
set ylabel "Displacement jump$/(8(1-\\\\nu^2)\\\\sigma a/(\\\\pi E))$"
set xrange [0:1]
set yrange [0:0.8]
set key outside right top Left reverse nobox spacing 1.4
sn = %.10f
st = %.10f
plot "penny_cod.txt" u 1:2 w p lt 1 pt 7 ps 0.4 t "opening, computed", \\
     sn*sqrt(1-x**2) w l lt 1 t "opening, Kassir--Sih", \\
     "penny_cod.txt" u 1:3 w p lt 5 pt 7 ps 0.4 t "sliding, computed", \\
     st*sqrt(1-x**2) w l lt 5 t "sliding, Kassir--Sih"
""" % (math.cos(b) ** 2, 2 * math.sin(b) * math.cos(b) / (2 - NU)))


if __name__ == "__main__":
    want = set(sys.argv[1:]) or {"penny", "conv", "cutpar", "growth", "cod"}
    if "penny" in want:
        penny_profile()
    if "conv" in want:
        penny_convergence()
    if "cutpar" in want:
        cutpar()
    if "growth" in want:
        growth()
    if "cod" in want:
        penny_cod()

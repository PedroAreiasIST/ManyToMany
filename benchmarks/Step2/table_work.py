#!/usr/bin/env python3
"""LaTeX rows of tab:workloads from a data directory (hot sequence per step; medians of the sampled rebuilds)."""
import csv, os, re, statistics as st, sys
D = sys.argv[1]
def med(xs): return st.median([x for x in xs if x == x])
def hot(kind, n, tag):
    p = os.path.join(D, f'work_{kind}_{n}_{tag}_hot.txt')
    return float(re.search(r'per_step_ms=([0-9.]+)', open(p).read().strip().splitlines()[-1]).group(1)) if os.path.exists(p) else None
def g(x):
    if x >= 100: return f'{x:.0f}'
    if x >= 10: return f'{x:.1f}'
    if x >= 1: return f'{x:.2f}'
    return f'{x:.3f}'
for kind in ('crack', 'bisect', 'erode'):
    for n in (12, 26, 55, 119):
        rows = list(csv.DictReader(open(os.path.join(D, f'work_{kind}_{n}_p.csv'))))
        el = float(rows[0]['elements']); e = int(med([float(r['edits']) for r in rows]))
        p, t = hot(kind, n, 'p'), hot(kind, n, 't')
        mm3 = med([float(r['t_mm3_ms']) for r in rows]); core = med([float(r['t_flat_core_par_ms']) for r in rows])
        allp = med([float(r['t_flat_all_par_ms']) for r in rows]); alls = med([float(r['t_flat_all_seq_ms']) for r in rows])
        ex = int(f'{el:.2e}'.split('e')[1]); mant = el / 10 ** ex
        print(f'{kind} & ${mant:.2f}\\times10^{{{ex}}}$ & {e} & {g(p)} & {g(t) if t else "--"} & {g(mm3)} & {g(core)} & {g(allp)} & {mm3/p:.0f}\\\\'
              f'   % seq all {g(alls)}; MM3set/pers {core/p:.0f}, all/pers {allp/p:.0f}')

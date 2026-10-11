# make_data.py — reduce the benchmark CSVs to the plot tables (all values from benchmarks/data).
import csv, glob, os, re, statistics as st
D = os.environ.get('STEP2_PLOT_DATA') or os.path.join(os.path.dirname(__file__), '..', '..', 'benchmarks', 'data')
def rows(name):
    p = os.path.join(D, name)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else None
def med(xs):
    xs = [x for x in xs if x == x]
    return st.median(xs) if xs else float('nan')
def f(r, k):
    try: return float(r[k])
    except: return float('nan')
def hot(kind, n, tag):
    p = os.path.join(D, f'work_{kind}_{n}_{tag}_hot.txt')
    if not os.path.exists(p): return float('nan')
    m = re.search(r'per_step_ms=([0-9.]+)', open(p).read().strip().splitlines()[-1])
    return float(m.group(1))
for kind in ('crack', 'bisect', 'erode'):
    out = []
    for n in (12, 26, 55, 119):
        p = rows(f'work_{kind}_{n}_p.csv')
        if not p: continue
        t = rows(f'work_{kind}_{n}_t.csv')
        el = float(p[0]['elements'])
        line = [el, med([f(r, 't_inc_min_ms') for r in p]), med([f(r, 't_inc_ms') for r in t]) if t else float('nan'),
                hot(kind, n, 'p'), hot(kind, n, 't'),
                med([f(r, 't_mm3_ms') for r in p]), med([f(r, 't_flat_core_par_ms') for r in p]),
                med([f(r, 't_flat_all_par_ms') for r in p]), med([f(r, 't_flat_all_seq_ms') for r in p]),
                med([f(r, 'edits') for r in p])]
        out.append(line)
    with open(f'scaling_{kind}.txt', 'w') as w:
        w.write('# elements inc_p_cold_min inc_t_single hot_p hot_t mm3 flat_core_par flat_all_par flat_all_seq edits_per_step\n')
        for l in out: w.write(' '.join(f'{x:.6g}' for x in l) + '\n')
c = rows('crossover_55.csv')
if c:
    with open('crossover.txt', 'w') as w:
        w.write('# batch fraction t_inc_ms t_mm3 t_flat_core_par t_flat_all_seq t_flat_all_par\n')
        for r in c: w.write(f"{r['batch']} {r['fraction']} {r['t_inc_ms']} {r['t_mm3_ms']} {r['t_flat_core_par_ms']} {r['t_flat_all_seq_ms']} {r['t_flat_all_par_ms']}\n")
i = rows('interleaved.csv')
if i:
    with open('interleaved.txt', 'w') as w:
        w.write('# elements bisections us_per_bisection_inc us_per_bisection_mm3\n')
        # time without collector pauses where the CSV records them
        for r in i: w.write(f"{r['elements']} {r['bisections']} {r.get('nogc_us_per_bisection_inc', r['us_per_bisection_inc'])} {r.get('nogc_us_per_bisection_mm3', r['us_per_bisection_mm3'])}\n")
q = rows('concurrent_55.csv')
if q:
    for impl in ('step2', 'mm3'):
        with open(f'concurrent_{impl}.txt', 'w') as w:
            w.write('# readers reads_per_s writer_step_median_ms\n')
            for r in q:
                if r['impl'] == impl: w.write(f"{r['readers']} {r['reads_per_s']} {r['writer_step_median_ms']}\n")
print('ok')

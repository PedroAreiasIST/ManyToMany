#!/usr/bin/env python3
"""Old (benchmarks/logs, 2026-10-01) against new (benchmarks/logs_enh) values of every Section 7 quantity."""
import os, re, sys
B = '/home/pedro/Desktop/NovosPapers/remesh/step2/benchmarks'
OLD, NEW = os.path.join(B, 'logs'), os.path.join(B, sys.argv[1] if len(sys.argv) > 1 else 'logs_enh')
def body(d, f):
    p = os.path.join(d, f)
    if not os.path.exists(p): return ''
    return '\n'.join(l for l in open(p).read().splitlines() if not l.startswith(('#', 'pinned')))
def work(d, kind, n, tag):
    t = body(d, f'work_{kind}_{n}_{tag}.log')
    m = re.search(r'hot sequence ([\d.]+) ms/step.*median edits/step (\d+).*MM3 ([\d.]+) ms, flat core par ([\d.]+) ms, flat all seq ([\d.]+) ms, flat all par ([\d.]+) ms', t)
    return m.groups() if m else None
print('Table workloads: kind n | pers trans MM3 MM3set all allseq edits MM3/pers  (old -> new)')
for kind in ('crack', 'bisect', 'erode'):
    for n in (12, 26, 55, 119):
        r = []
        for d in (OLD, NEW):
            p = work(d, kind, n, 'p'); t = work(d, kind, n, 't') if kind == 'crack' else None
            if not p: r.append('--'); continue
            r.append(f'{float(p[0]):.3g} {float(t[0]):.3g} ' if t else f'{float(p[0]):.3g} -- ')
            r[-1] += f'{float(p[2]):.3g} {float(p[3]):.3g} {float(p[5]):.3g} {float(p[4]):.4g} e{p[1]} r{float(p[2])/float(p[0]):.0f}'
        print(f'{kind:6s} {n:3d} | {r[0]}  ->  {r[1]}')
for f in ('mm3ref.log', 'persist_55.log', 'persist_119.log', 'memory_55.log', 'memory_119.log', 'reads_55.log',
          'reads_119.log', 'concurrent_55.log', 'interleaved.log', 'crossover_55.log', 'endtoend_55.log',
          'profile_crack_55_t.log', 'profile_crack_55_p.log', 'profile_bisect_55_p.log'):
    print(f'\n===== {f}\n--- old\n{body(OLD, f)}\n--- new\n{body(NEW, f)}')

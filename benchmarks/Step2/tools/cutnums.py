"""Numbers of Section 8.8 from the fin4 (Ryzen) and cutpar7 (EPYC) data: no-pause minima."""
import csv, os
B = '/home/pedro/Desktop/NovosPapers/remesh/step2/benchmarks'
def rd(d, f):
    return list(csv.DictReader(open(os.path.join(B, d, f))))
def sg(v, n=3):
    """n significant digits, integer part never rounded; thin-space grouping from 1000 on"""
    if v >= 1000: s = f"{round(v):,}".replace(",", "\\,")
    elif v >= 100: s = f"{round(v):d}"
    else:
        digs = n - (len(str(int(v))) if v >= 1 else 0)
        s = f"{v:.{max(digs,0)}f}" if v >= 1 else f"{v:.{n}g}"
        if v < 1 and "e" in s: s = f"{v:.2f}"
        if v < 1: s = f"{v:.2f}"
    return s
def cutpar(d, n, g):
    rows = rd(d, f'cutpar_{n}_{g}.csv')
    seq = next(r for r in rows if r['impl'] == 'sequential')
    par = {int(r['threads']): r for r in rows if r['impl'] == 'parallel'}
    return seq, par
def batch(d, n, g):
    rows = rd(d, f'batchcut_{n}_{g}.csv')
    return {(r['impl'], r['marks']): r for r in rows}
if __name__ == '__main__':
    print('--- Table cutpar (no pauses) ---')
    for n, g in ((55, 4), (55, 8), (119, 4), (119, 8)):
        seq, par = cutpar('data_fin4', n, g)
        s = float(seq['nogc_ms_min']); edges = int(seq['edges'])
        row = [sg(s)] + [sg(float(par[t]['nogc_ms_min'])) for t in (1, 2, 4, 8, 16)]
        su = {t: s / float(par[t]['nogc_ms_min']) for t in par}
        el = {t: float(seq['ms_min']) / float(par[t]['ms_min']) for t in par}
        print(n, g, edges, row, 'speedup16 %.2f (elapsed %.2f)' % (su[16], el[16]), '1thr ratio %.3f' % (float(par[1]['nogc_ms_min']) / s),
              'pauses seq %s/%s ms, 16thr %s ms' % (seq['gc_count_at_min'], seq['gc_ms_at_min'], par[16]['gc_ms_at_min']))
    print('--- Table batchcut ---')
    for n, g in ((55, 4), (55, 8), (119, 4), (119, 8)):
        b = batch('data_fin4', n, g)
        for m in ('all', 'one'):
            l, bt, rb = (float(b[(k, m)]['nogc_ms_min']) for k in ('local', 'batch', 'rebuild'))
            print(n, g, m, b[('local', m)]['edges'], sg(l), sg(bt), sg(rb), sg((bt + rb) / l), '(batch+rebuild %.0f)' % (bt + rb))
    print('--- EPYC (cutpar7, no pauses) ---')
    for d, n, g, th in (('data_deu7', 119, 4, (1, 2, 4, 8, 16, 32, 64, 128)), ('data_deu7', 119, 8, (1, 2, 4, 8, 16, 32, 64, 128)),
                        ('data_deu7', 160, 8, (1, 16, 32, 64, 128)), ('data_deu7/gc64', 119, 8, (1, 16, 32, 64, 128)), ('data_deu7/gc64b', 160, 8, (1, 16, 32, 64, 128))):
        seq, par = cutpar(d, n, g); s = float(seq['nogc_ms_min'])
        print(d, n, g, 'seq nogc %.0f elapsed %.0f' % (s, float(seq['ms_min'])), ' '.join('%d:%.2f(%.2f)' % (t, s / float(par[t]['nogc_ms_min']), float(seq['ms_min']) / float(par[t]['ms_min'])) for t in th if t in par),
              '| par128 %.0f' % float(par[128]['nogc_ms_min']))

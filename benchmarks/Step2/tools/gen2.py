import sys, re
sys.path.insert(0, '/home/pedro/.claude/skills/finalpass/scripts')
import lyxflat as lf
src, spec, out = sys.argv[1:4]
lines = lf.load(src)
paras = lf.parse(lines)
flats = [lf.flat(lines, p)[0].strip() for p in paras]
blocks, cur = [], None
for raw in open(spec, encoding='utf-8').read().split('\n'):
    if raw.startswith('@@'):
        cur = {'hint': raw[2:], 'P': None, 'OLD': None, 'NEW': None}
        blocks.append(cur)
    elif raw.startswith('P: '):
        cur['P'] = raw[3:].strip()
    elif raw.startswith('OLD: '):
        cur['OLD'] = raw[5:].strip()
    elif raw.startswith('NEW: '):
        cur['NEW'] = raw[5:].strip()
bad = 0
with open(out, 'w', encoding='utf-8') as fo:
    for b in blocks:
        if b['P']:
            hits = [f for f in flats if f.startswith(b['P'])]
            if len(hits) != 1:
                print('PREFIX', len(hits), 'hits:', b['P'][:70]); bad += 1; continue
            old = hits[0]
        else:
            old = b['OLD']
        fo.write('@@%s\nOLD: %s\nNEW: %s\n\n' % (b['hint'], old, b['NEW']))
print(len(blocks), 'blocks,', bad, 'bad')

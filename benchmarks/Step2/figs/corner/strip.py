import sys
src, dst = sys.argv[1:3]
L = open(src).read().split('\n')
out = []; lw = None; col = None; removed = 0
for l in L:
    t = l.split()
    if t:
        op = t[-1]
        if op == 'slw': lw = t[0]
        elif op == 'Cr': col = tuple(t[:3])
        elif op == 'L' and lw == '100' and col == ('0.000', '0.000', '0.000'):
            removed += 1; continue
    out.append(l)
open(dst, 'w').write('\n'.join(out))
print(src, 'removed', removed)

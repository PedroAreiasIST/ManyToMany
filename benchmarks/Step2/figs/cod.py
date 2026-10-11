#!/usr/bin/env python3
"""Crack-face displacement jump of a dumped penny run against the closed forms.

    python3 cod.py <dump dir> <beta_deg> <out.txt>

Crack-face node pairs are matched by coordinates (the copy made by the split lies on the + side); the jump
[u] = u(+) - u(-) gives the normal opening [u].n and the sliding [u].s along the shear direction
s = -(e1) (the in-plane projection of the load axis), against
  opening = 8 (1 - nu^2) sigma_n sqrt(a^2 - r^2) / (pi E),  sliding = 16 (1 - nu^2) tau sqrt(a^2 - r^2) / (pi (2 - nu) E).
Writes r/a, opening, sliding (normalized by 8 (1-nu^2) sigma sqrt... at r = 0 for beta = 0, i.e. by 4(1-nu^2) sigma a/E * 2/pi)."""
import math
import os
import sys

import numpy as np

E, NU, SIG, A = 1.0, 0.3, 1.0, 0.2
d, beta, out = sys.argv[1], math.radians(float(sys.argv[2])), sys.argv[3]
X = np.loadtxt(os.path.join(d, "nodes.txt"))
U = np.loadtxt(os.path.join(d, "u.txt"))
T4 = np.loadtxt(os.path.join(d, "tets.txt"), dtype=np.int64)
B = np.loadtxt(os.path.join(d, "boundary.txt"), dtype=np.int64)
n = np.array([math.sin(beta), 0, math.cos(beta)])
e1 = np.array([math.cos(beta), 0, -math.sin(beta)])
e2 = np.array([0.0, 1.0, 0.0])
onbox = np.zeros(len(B), bool)
for k in range(3):
    for v in (-1.0, 1.0):
        onbox |= np.all(np.abs(X[B[:, :3], k] - v) < 1e-9, axis=1)
C = B[~onbox]
side = np.sign((X[T4[C[:, 3]]].mean(axis=1)) @ n)
plus = np.unique(C[side > 0][:, :3])
minus = np.unique(C[side < 0][:, :3])
key = lambda v: tuple(np.round(X[v], 12))
mp = {key(v): v for v in minus}
pairs = [(v, mp[key(v)]) for v in plus if key(v) in mp and mp[key(v)] != v]
sn, tau = SIG * math.cos(beta) ** 2, SIG * math.sin(beta) * math.cos(beta)
ref = 8 * (1 - NU ** 2) * SIG * A / (math.pi * E)            # opening at the center for sigma_n = sigma
s_dir = -e1                                                    # shear direction on the + face
with open(out, "w") as f:
    for p, m in pairs:
        x = X[p]
        r = math.hypot(x @ e1, x @ e2)
        if r > A:
            continue
        j = U[p] - U[m]
        # sliding measured along the shear direction projected onto the radial line is not needed: report [u].s
        f.write("%.6f %.6f %.6f\n" % (r / A, (j @ n) / ref, (j @ s_dir) / ref))
print("%d pairs; reference opening at the center %.5f (normal), %.5f (sliding)" % (
    len(pairs), 8 * (1 - NU ** 2) * sn * A / (math.pi * E) / ref, 16 * (1 - NU ** 2) * tau * A / (math.pi * (2 - NU) * E) / ref))

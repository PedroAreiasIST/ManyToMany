import io,re
s=io.open('before9.tex',encoding='utf-8').read()
def rep(old,new,count=1):
    global s
    pat=r'\s+'.join(re.escape(t) for t in old.split())
    m=list(re.finditer(pat,s))
    assert len(m)==count,(old[:70],len(m))
    for mm in reversed(m): s=s[:mm.start()]+new+s[mm.end():]
# ---- abstract
a=s.index('\\begin{abstract}')+len('\\begin{abstract}'); b=s.index('\\end{abstract}')
abstract=r'''
Finite element codes with evolving topology rebuild connectivity at
a cost proportional to the mesh. We maintain the transpose, node--node
pattern, assembly destinations, facet adjacency, boundary facets and
coloring under unplanned local edits, at a cost proportional to the
edit plus a directory copy after each commit. Counting relations give
exact update rules; stable pattern slots keep the destinations of
unedited elements. Copy-on-write arrays give constant-time commit
and rollback with lock-free readers. A Steiner-point-free cutter
reproduces the batch cutter, reverts locally, grows cracks from two
level sets and cuts independent groups concurrently. An oracle verifies
every structure after each batch. On $10^{7}$ tetrahedra, a crack
advance takes 3.0~ms, against 849~ms for the rebuild of our earlier
MM3 framework and 300~ms for a tuned 16-thread rebuild of its structures;
16 threads cut independent groups there 7.3 to 7.4 times faster than
the sequential cut, excluding collector pauses. Stress intensity
factors from the interaction integral differ from the penny-crack
solutions of Sneddon and of Kassir and Sih by at most 0.29 and 0.51\,\%
of $K_{0}$ at $h_{\min}/a=1/64$ and, with a free-surface term, from
the surface-crack values of Newman and Raju by 3.1 and 1.5\,\%. Elliptical-crack
Paris-law steps roll back when relative factors change by more than
a tolerance; fitting, crack increment and commit take 0.52 to 0.66\,\%
of a step. A corner crack without a reference solution grows by the
Erdogan--Sih criterion to within 0.115 of the boundary on 6.5 million
tetrahedra. Code: ManyToMany repository.
'''
s=s[:a]+abstract+s[b:]
# ---- contribution 7
rep('''a free-surface term for the interaction integral; Paris-law growth''','''the interaction integral with the surface term that a front ending on a free surface requires; Paris-law growth''')
# ---- theory
rep('''An edit
costs $\\BigO\\bigl(k\\,d_{\\max}+k^{2}\\ell_{\\max}+(k+n_{\\Sigma})\\,k\\,d_{\\max}\\bigr)$,
where $d_{\\max}$ is the largest transpose row, $\\ell_{\\max}$ the
largest pattern row and $n_{\\Sigma}$ the number of facet sharers
recomputed, plus, for each relocated row, the copy of at most $\\ell_{\\max}$
entries and, amortized, a constant per added element, node, pattern
slot or boundary facet for the growth of the arrays. Apart from the
first writes of an epoch (Proposition~\\ref{prop:versions}) and the
growth of the arrays, the cost does not depend on $n_{E}$ or $n_{N}$.''','''An edit
costs $\\BigO\\bigl(k\\,d_{\\max}+k^{2}\\ell_{\\max}+(k+n_{\\Sigma})\\,k\\,d_{\\max}\\bigr)$,
where $d_{\\max}$ is the largest transpose row, $\\ell_{\\max}$ the
largest pattern row and $n_{\\Sigma}$ the number of facet sharers
recomputed. In addition, each relocated row costs the copy of at most
$\\ell_{\\max}$ entries, and the growth of the arrays costs an amortized
constant per added element, node, pattern slot or boundary facet.
Apart from the first writes of an epoch (Proposition~\\ref{prop:versions})
and the growth of the arrays, the cost does not depend on $n_{E}$
or $n_{N}$.''')
rep('''a relocation copies at most $\\ell_{\\max}$ entries, and an edit
relocates at most one row per node of the new row.''','''a relocation copies at most $\\ell_{\\max}$ entries, an edit
relocates at most one row per node of the new row, and the arrays
grow geometrically, which costs an amortized constant per added entry.''')
rep('the concurrent planner computes the exact bound of each\ngroup.','the concurrent planner computes this bound for every row that a\ngroup touches.')
rep('Reads in both modes take two dependent loads,','In both modes, reaching an entry takes two dependent loads,')
rep('It is the fastest rebuild measured; the incremental mesh','It is the fastest rebuild measured for each set of structures; the incremental mesh')
# ---- performance
rep('''are 0.50, 4.7 and 56~ms and 1.1, 12 and 173~ms. All comparisons
below use measurements on the same machine.''','''are 0.50, 4.7 and 56~ms and 1.1, 12 and 173~ms. Without the pauses
the same protocol measures 32 and 92~ms at $10^{6}$ quadrilaterals.
All comparisons below use measurements on the same machine.''')
rep('''after full garbage collections, and the tables report the median of
these minima over the sampled steps.''','''after full garbage collections, and the tables report the median of
these minima over the sampled steps. Rebuild times are not corrected
for collector pauses; at most 10\\,\\% of the timed regions of a run
contain a collection.''')
rep('''Above these fractions, rebuilding and initializing an incremental
mesh is cheaper.''','''Above these fractions, a rebuild is cheaper; initializing an incremental
mesh from it adds to its cost.''')
rep('''Commit and rollback cost a fraction of a microsecond at both sizes.''','''Commit without edits and rollback cost a fraction of a microsecond
at both sizes, and a commit after a crack increment costs 4 to 14~\\textmu s.''')
rep('''STA-FEM reports overall speed-ups of 1.37 to 1.61 \\citep{Acharya2026stafem}.''','''The ratio depends on the iteration count and tends to 1 as the solution
dominates. STA-FEM reports overall speed-ups of 1.37 to 1.61 on other
problems \\citep{Acharya2026stafem}.''')
rep('''takes the
same 2.1~s on 64 threads, while the sequential cut takes 26.9~s.''','''takes
2.1~s on 64 threads (2.0~s with 16 heaps), while the sequential cut
takes 26.9~s.''')
rep('''Readers of committed
versions reach $8.3\\times10^{8}$ reads per second while a writer
commits, $6\\times10^{3}$ to $3\\times10^{4}$ times the rate of the
MM3 reader--writer lock.''','''Readers of committed
versions reach $8.3\\times10^{8}$ reads per second while a writer
commits, whereas MM3 readers, which wait on the write lock during
each transpose rebuild, obtain $4.2\\times10^{3}$ to $2.7\\times10^{4}$.''')
rep('the readers\nobtain $4.2\\times10^{3}$ to $2.7\\times10^{4}$ reads per second','the readers\nobtain 4.2 to 27 thousand reads per second')
# ---- examples
rep('''$J$-integral
of the same domain, which needs no auxiliary field, differs from the
energy release rate computed from the factors by at most 4.0, 2.4,
1.0 and 0.55\\,\\% for $\\beta=0$ on the penny-crack meshes''','''$J$-integral
of the same domain, which needs no auxiliary field, differs from the
energy release rate computed from the factors by at most 4.0, 2.4,
1.0 and 0.55\\,\\% for normal loading on the penny-crack meshes''')
rep('so the finite body\ndoes not explain the remaining error.','so the finite body\nis not the dominant source of the remaining error.')
rep('''SimPlas \\citep{simplas} (4-node tetrahedra, linear analysis) solves''','''The finite element code SimPlas \\citep{simplas} (4-node tetrahedra, linear analysis) solves''')
rep('up to 4.0 \\%.','up to 4.0\\,\\%.')
rep('''on
the three finest meshes, after a rise from 4.11 to 4.74\\,\\% over
the three coarsest;''','''on
the three finest meshes, after a rise from 4.11 to 4.74\\,\\% over
the three coarsest;''')
rep('3.7, 3.6 and 3.1\\,\\% on','3.66, 3.56 and 3.14\\,\\% on')
rep('graph $h=w(p,q)$','graph $\\zeta=w(s_{1},s_{2})$')
rep('There is no reference value for the corner crack.','''There is no reference value for the corner crack: the example tests
the cutter and the versioned structures on a curved mixed-mode path,
and its factors are verified only through the penny crack above.''')
rep('''The tube radius grows from 0.036 to 0.23; at fronts 11 and 12 it
(0.213 and 0.228) exceeds the distances from the front to the faces
$x=1$ (0.19 and 0.115) and $y=1$ (0.21 and 0.13), which therefore
cut the tube.''','''The tube radius grows from 0.036 to 0.228; at fronts 11 and 12
(0.213 and 0.228) it exceeds the distances from the front to the
faces $x=1$ (0.19 and 0.115) and $y=1$ (0.21 and 0.13), which therefore
cut the tube.''')
rep('into a Kuhn mesh by cutting every edge','into a Kuhn mesh of the cube $[-1,1]^{3}$ by cutting every edge')
rep('The concurrent cut gives the rows of the sequential cut, and every','The concurrent cut gives the rows of the sequential cut (up to the permutation of the appended children), and every')
io.open('after9.tex','w',encoding='utf-8').write(s)
print('ok')

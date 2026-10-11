import io,re,sys
s=io.open('before8.tex',encoding='utf-8').read()
def rep(old,new,count=1):
    global s
    pat=r'\s+'.join(re.escape(t) for t in old.split())
    m=list(re.finditer(pat,s))
    assert len(m)==count,(old[:70],len(m))
    for mm in reversed(m): s=s[:mm.start()]+new+s[mm.end():]
# ---------- abstract
a=s.index('\\begin{abstract}')+len('\\begin{abstract}'); b=s.index('\\end{abstract}')
abstract='''
Finite element codes with evolving topology rebuild connectivity at
a cost proportional to the mesh. We maintain the transpose with positions,
the node--node pattern, assembly destinations, facet adjacency, boundary
facets and coloring under unplanned local edits, at a cost proportional
to the edit plus a directory copy after each commit. Counting relations
give exact update rules, and stable pattern slots preserve the destinations
of unedited elements. Copy-on-write arrays give constant-time commit
and rollback with lock-free readers. A Steiner-point-free cutter built
on these edits reproduces the batch cutter, reverts locally, grows
cracks from two level sets and cuts independent groups concurrently.
An oracle verifies every structure after each batch. On $10^{7}$
tetrahedra, a crack advance takes 3.0~ms, against 849~ms for the
rebuild of our earlier MM3 framework and 300~ms for a tuned 16-thread
rebuild. Sixteen threads cut independent groups 7.3 to 7.4 times
faster than the sequential cut, excluding collector pauses. Stress
intensity factors from the interaction integral differ from the penny-crack
solutions of Sneddon and of Kassir and Sih by at most 0.29 and 0.51\\,\\%
of the reference factor at $h_{\\min}/a=1/64$. With a free-surface
term, they differ from the surface-crack values of Newman and Raju
by 3.1 and 1.5\\,\\%. In Paris-law growth of an elliptical crack,
steps roll back when the relative factors change by more than a tolerance,
and fitting, crack increment and commit take 0.52 to 0.66\\,\\% of
a step. A corner crack grows by the Erdogan--Sih criterion to within
0.115 of the block boundary on 6.5 million tetrahedra. The code is
in the ManyToMany repository.
'''
s=s[:a]+abstract+s[b:]
# ---------- introduction
rep('''Mesh databases update adjacency locally: half-edge-type
structures, the modifiable arrays of Ibanez and Shephard \\citep{Ibanez2017modifiable,Ibanez2016pumi},
the patch-local GPU updates of Mahmoud et al.\\ \\citep{Mahmoud2025dynamic},
the thread-parallel adaptation of Gorman et al.\\ \\citep{Rokos2015thread},
which keeps node--node and node--element lists current on independent
sets, and TopS and ParTopS, which insert cohesive elements locally
and in parallel \\citep{Paulino2008general,Espinha2009partops}. They
maintain adjacencies, but no slotted sparsity pattern, assembly destinations
or coloring.''','''Mesh databases update adjacency locally: half-edge-type
structures, the modifiable arrays of Ibanez and Shephard \\citep{Ibanez2017modifiable,Ibanez2016pumi},
and TopS and ParTopS, which insert cohesive elements locally and
in parallel \\citep{Paulino2008general,Espinha2009partops}. Adaptation
codes do the same: Mahmoud et al.\\ \\citep{Mahmoud2025dynamic} update
patches on the GPU, and Gorman et al.\\ \\citep{Rokos2015thread} keep
node--node and node--element lists current on independent sets. None
of these maintains a slotted sparsity pattern, assembly destinations
or a coloring for assembly.''')
rep('but neither derives topology-based patterns or preserves assembly','but neither derives topology-based patterns nor preserves assembly')
rep('''Persistent representations of triangulations exist \\citep{Blelloch2001persistent}.
To our knowledge, neither approach has been applied to the inverse
incidence with positions, the sparsity pattern and the assembly destinations
of a finite element mesh.''','''Blelloch et al.\\ represent triangulations persistently \\citep{Blelloch2001persistent},
but not the solver structures derived from them. To our knowledge,
neither incremental view maintenance nor persistence has been applied
to the inverse incidence with positions, to the sparsity pattern under
unplanned edits or to the assembly destinations of a finite element
mesh.''')
rep('''Three-dimensional crack growth with stress intensity factors
from domain integrals is established in remeshing codes \\citep{Carter2000automated}
and in the extended finite element method with two level sets \\citep{Moes2002non}.''','''Three-dimensional crack growth is established with remeshing
\\citep{Carter2000automated} and with the extended finite element
method and two level sets \\citep{Moes2002non}.''')
rep('''to the edit, plus a copy of the array directories after the first
edit that follows a commit''','''to the edit, plus, at the first edit after a commit, a copy of the
array directories''')
# ---------- theory
rep('so a facet enters or leaves it in constant','so a facet enters or leaves it in amortized constant')
rep('''recomputed, plus, for each relocated row, the copy of at most twice
the longest pattern row. Apart from the first writes of an epoch (Proposition~\\ref{prop:versions}),
the cost does not depend on $n_{E}$ or $n_{N}$.''','''recomputed, plus, for each relocated row, the copy of at most $\\ell_{\\max}$
entries and, amortized, a constant per added element, node, pattern
slot or boundary facet for the growth of the arrays. Apart from the
first writes of an epoch (Proposition~\\ref{prop:versions}) and the
growth of the arrays, the cost does not depend on $n_{E}$ or $n_{N}$.''')
rep('a\nrelocation copies at most twice the longest pattern row','a relocation copies at most $\\ell_{\\max}$ entries')
rep('If no child of the cut has been edited since the cut','If no element has been edited since the cut')
rep('''The
result equals the sequential cut in rows, slots, destinations and
boundary set;''','''Up to
the permutation of the appended children, the result equals the sequential
cut in rows, slots, destinations and boundary set;''')
rep('''While these bounds are far from the limit, the check
costs $\\BigO(1)$; otherwise it uses the current length, plus the
distinct new columns, of the rows whose bounds exceed the limit.''','''In the sequential cut, the check costs $\\BigO(1)$ while
these bounds are far from the limit and otherwise uses the current
length, plus the distinct new columns, of the rows whose bounds exceed
the limit; the concurrent planner computes the exact bound of each
group.''')
rep('and its star then has more than two components','and its star may then have more than two components')
rep('each element waiting until its conflicting neighbors of smaller index are colored','each element waiting until its uncolored neighbors of smaller index are colored')
rep('''It is the fastest rebuild measured, and it initializes
the incremental mesh with the greedy coloring in element order; edits
then repair the coloring by''','''It is the fastest rebuild measured; the incremental mesh
is initialized from it with the sequential greedy coloring in element
order, and edits then repair the coloring by''')
rep('''for \\code{ReplaceElement}, \\code{AddElement}, \\code{RemoveElement}
and \\code{AddNode}.''','''for \\code{ReplaceElement}, \\code{AddElement} and \\code{RemoveElement};
\\code{AddNode} appends a node.''')
rep('compaction replaces the pattern arrays','a compaction gives later versions new pattern arrays')
rep('The tuned parallel and sequential rebuilds agree exactly on all structures','The parallel and sequential forms of the tuned rebuild agree exactly on all structures')
rep('The tuned rebuild uses 206:','The tuned rebuild of all structures uses 206:')
# ---------- performance
rep('''The pauses depend on the
state of the heap and not on the operation. The sequential cut of
$10^{7}$ elements with 64 groups takes 6.0~s, with six collections
and 1.2~s of pauses, in a process that holds only the mesh, and 4.8
to 5.0~s without collections when a second mesh is alive; without
pauses it takes 4.8 to 4.9~s in all four states (Table~\\ref{tab:gc}).''','''The pauses depend on the
state of the heap, not on the order of the marks. In a process that
holds only the mesh, the sequential cut of $10^{7}$ elements with
64 groups takes 6.0~s, including six collections and 1.2~s of pauses.
When a second mesh is alive, it takes 4.8 to 5.0~s without collections.
Without pauses it takes 4.8 to 4.9~s in all four states (Table~\\ref{tab:gc}).''')
rep(r'''in transient mode (Section~\ref{sec:versioning}), accounting for
most of the persistent-mode overhead.''',r'''in transient mode (Section~\ref{sec:versioning}), comparable to the
persistent-mode overhead of 1.6~ms per step.''')
rep(r'it is 6.4 and 6.9 times faster',r'it is 6.3 and 6.9 times faster')
rep(r'''(13.1 times for 512 groups
with 64 collector heaps; Figure~\ref{fig:cutpar}).''',r'''(Figure~\ref{fig:cutpar}).
With 64 collector heaps the concurrent cut of 512 groups takes the
same 2.1~s on 64 threads, while the sequential cut takes 26.9~s.''')
rep(r'differ by at most 1.4\,\%',r'differ by at most 1.5\,\%')
rep(r'''The readers reach $2.7\times10^{7}$
reads per second with one thread and $8.3\times10^{8}$ with fifteen,
that is $2.7\times10^{7}$ and $5.5\times10^{7}$ per thread. Without
a writer, each reader reads 76.6 and 75.5~M/s at one and fifteen
threads, and a writer that commits every millisecond without editing
leaves the rates at 75.7 and 75.1~M/s. The loss therefore comes from
the edits, not from commits; pinning the writer and the readers to
separate cores gives the same rates, and the decrease of the loss
with the number of readers is not explained. Repeat runs differ by
at most 13~\%.''',r'''The readers reach 26.7~M/s with one thread and 831~M/s with fifteen,
that is 26.7 and 55.4~M/s per thread. Without a writer, each reader
reads 76.6 and 75.5~M/s at one and fifteen threads, and a writer
that commits every millisecond without editing leaves the rates at
75.7 and 75.1~M/s. The loss therefore comes from the edits, not from
the commits. It falls from 65\,\% with one reader to 27\,\% with fifteen,
and pinning the writer and each reader to its own core gives the
same snapshot rates. Repeat runs differ by at most 15\,\%.''')
rep(r'$4.3\times10^{3}$ to $2.7\times10^{4}$ reads per second',r'$4.2\times10^{3}$ to $2.7\times10^{4}$ reads per second')
# ---------- examples
rep(r'''The errors of the second and third factors are the
same for the two smaller radii; the refined band, and with it the
mesh, grows with $r_{q}$.''',r'''The refined band grows with $r_{q}$ (3.6, 6.8 and 13.2 million tetrahedra
for $r_{q}=0.05$, 0.075 and 0.11 at $\beta=30^{\circ}$), yet the errors
of the second and third factors are equal for the two smaller radii.''')
rep(r'''For a crack of
half the radius in the same cube ($a=0.1$, faces at $10a$, $r_{q}=0.0375$,
linear factor), the errors at $h_{\min}/a=1/8$ and $1/16$ are 3.77
and 1.49\,\% ($\beta=0$) and 2.67, 1.78 and 1.80\,\% and 1.39,
1.32 and 1.12\,\% ($\beta=30^{\circ}$). They do not fall when the
faces move from $5a$ to $10a$.''',r'''For a crack of
half the radius in the same cube ($a=0.1$, faces at $10a$, $r_{q}=0.0375$,
linear factor), the errors at $h_{\min}/a=1/8$ are 3.77\,\% ($\beta=0$)
and 2.67, 1.78 and 1.80\,\% ($\beta=30^{\circ}$). At $h_{\min}/a=1/16$
they are 1.49\,\% and 1.39, 1.32 and 1.12\,\%. These exceed the errors
with the faces at $5a$ (Table~\ref{tab:penny}), so the finite body
does not explain the remaining error.''')
rep(r'''on
the three finest meshes and 4.74\,\% at $h_{\min}/a=0.083$, without
a monotone trend;''',r'''on
the three finest meshes, after a rise from 4.11 to 4.74\,\% over
the three coarsest;''')
rep(r'to within 1.5\,\% for $h_{\min}=h_{0}/8$ and 0.5\,\% for',r'to within 1.6\,\% for $h_{\min}=h_{0}/8$ and 0.6\,\% for')
rep(r'half-length within 0.53\,\%',r'half-length within 0.54\,\%')
rep(r'''The crack faces are curved, so the integral
over them, which vanishes for planar faces, is omitted; the penny-crack
verification is planar.''',r'''The integral over the crack faces vanishes only for planar
faces; it is omitted here although the faces are curved, and the
penny-crack verification is planar.''')
rep(r'''the crack surface follows the element faces within
a fraction of an element.''',r'''the crack surface, made of element faces,
lies within a fraction of an element of the level set.''')
rep(r'bring the crack within 0.114 of the boundary',r'bring the crack within 0.115 of the boundary')
rep(r'''The tube radius grows from 0.036 to 0.23; at fronts 11 and 12 it exceeds
the distances 0.19 and 0.115 from the front to the faces $x=1$ and
$y=1$, which then cut the tube.''',r'''The tube radius grows from 0.036 to 0.23; at fronts 11 and 12 it
(0.213 and 0.228) exceeds the distances from the front to the faces
$x=1$ (0.19 and 0.115) and $y=1$ (0.21 and 0.13), which therefore
cut the tube.''')
# ---------- conclusions
rep(r'''through exact update rules: pair counts update the pattern
by the difference of two cliques, positive counts preserve the destinations
of unedited elements, and copy-on-write chunks preserve committed
versions.''',r'''through exact update rules: pair counts update the pattern
by the difference of two cliques, and positive counts preserve the
destinations of unedited elements. Copy-on-write chunks preserve
committed versions.''')
rep(r'''and an element edit costs 1.5 to
2.9~\textmu s.''',r'''and an element edit costs 1.5 to
2.9~\textmu s from $10^{4}$ to $10^{7}$ elements.''')
rep(r'''Independent groups
are cut 6.4 to 7.4 times faster than sequentially on 16 cores, without
collector pauses.''',r'''Independent groups
are cut 6.3 to 7.4 times faster than sequentially on 16 cores at
$10^{6}$ and $10^{7}$ elements, without collector pauses.''')
rep(r'''Interaction-integral factors converge to the solutions of Sneddon
and of Kassir and Sih, with largest errors of 0.29 and 0.51\,\%
of $K_{0}$ at $h_{\min}=a/64$.''',r'''Interaction-integral factors differ from the solutions of Sneddon
and of Kassir and Sih by at most 0.29 and 0.51\,\% of $K_{0}$ at
$h_{\min}=a/64$, and the errors fall as $h_{\min}$ halves.''')
rep(r'to within 0.114 of the block boundary',r'to within 0.115 of the block boundary')
rep(r'come from one workstation and one node of each supercomputer.',r'come from one workstation and one A64FX and one EPYC 7742 node of Deucalion.')
io.open('after8.tex','w',encoding='utf-8').write(s)
print('patched')

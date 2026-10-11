import sys; sys.path.insert(0,'/tmp/claude-1000/-home-pedro-Desktop-NovosPapers-remesh-step2/538c0fdd-ec6a-4ebb-b66b-9b73d27a2640/scratchpad/tools')
import lyxtab as T, cutnums as C
def make():
    rows=[]
    for n,lab in ((55,'$10^{6}$'),(119,'$10^{7}$')):
        L=open(f'/home/pedro/Desktop/NovosPapers/remesh/step2/benchmarks/data_fin4/cutgc_{n}_4.csv').read().strip().split('\n')
        hdr=L[0].split(',')[1:]
        for i,l in enumerate(L[1:]):
            p=l.split(','); v=p[0]+','+p[1]; f=dict(zip(hdr,p[2:]))
            gcms=float(f['gc_ms_at_min'])
            rows.append([lab if i==0 else '', ('1' if 'one mesh' in v else '2') if i%2==0 else '', 'groups' if v.startswith('group') else 'edges',
                         C.sg(float(f['ms_min'])), f['gc_count_at_min'], C.sg(gcms) if gcms>0 else '0', C.sg(float(f['nogc_ms_min'])), format(int(float(f['heap_mb'])),',').replace(',','\\,')])
    cap=("Sequential cut in four heap states (Ryzen~9 9950X, 64 groups, minimum of five samples, ms): time with collector pauses, collections and pauses in the minimum run, "
         "time without pauses, and managed heap (MB). Meshes is the number of meshes alive, the second being an MM3 mesh of the same tetrahedra; the marks are ordered by group or by edge.")
    hdr=['elements','meshes','marks','time','collections','pauses','no pauses','heap']
    return rows, T.float_table(cap,'tab:gc',hdr,rows,['left','left','left','right','right','right','right','right'])
def replace(s):
    rows,t=make()
    k=s.index('name "tab:gc"'); a=s.rindex('\\begin_layout Standard\n\\begin_inset Float table',0,k)
    tail='\\end_inset\n\n\n\\end_layout\n\n\\end_inset\n\n\n\\end_layout\n\n'
    e=s.index(tail,s.index('</lyxtabular>',k))+len(tail)
    return s[:a]+'\\begin_layout Standard\n'+t+'\\end_layout\n\n'+s[e:]
if __name__=='__main__':
    s=open(sys.argv[1],encoding='utf-8',newline='').read()
    open(sys.argv[2],'w',encoding='utf-8',newline='').write(replace(s))

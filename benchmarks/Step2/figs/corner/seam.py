import re,sys
src,dst=sys.argv[1:3]
s=open(src).read()
n=0
def sub(m):
    global n; n+=1
    return m.group(1)+' closepath gsave fill grestore stroke }'
s2=re.sub(r'(/[TQP] \{ newpath (?:moveto|lineto|\s)+?)\s*closepath\s+fill \}',sub,s)
open(dst,'w').write(s2); print(src,'patched',n)

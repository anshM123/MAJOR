# log10 sizes of the realised direct coefficients (K=3299 table), non-rigorous diagnostics
import sys,math
fn=sys.argv[1]
f=open(fn);N,D=map(int,f.readline().split())
def l10(s):
    s=s.lstrip('-')
    return len(s)-1+math.log10(int(s[:15]))-14 if len(s)>15 else (math.log10(int(s)) if int(s)>0 else -1e9)
LA=[];LG=[];sg=[]
for line in f:
    a,r,g,t=line.split();LA.append(l10(a)-D);LG.append(l10(g)-D);sg.append(a[0]=='-')
import json
json.dump(dict(LA=LA,LG=LG,neg=sg),open('sizes.json','w'))
for k in list(range(0,20))+list(range(20,len(LA),100)):
    print(k,'%.1f %.1f'%(LA[k],LG[k]),'-' if sg[k] else '+')
for c in (0.27,0.3,0.4,0.5,0.6,0.7,0.77):
    v=[LA[k]+k*math.log10(c) for k in range(len(LA))]
    km=max(range(len(v)),key=lambda k:v[k])
    print('c',c,'max log10|A_k|c^k=%.1f at k=%d; last=%.1f'%(v[km],km,v[-1]))

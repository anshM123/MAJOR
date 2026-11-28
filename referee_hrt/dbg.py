import sys; sys.path.insert(0,'.')
from btz_cone_test import *
rng=np.random.default_rng(11)
N=3; fs=facets(N)
for t in range(300):
    n=int(rng.integers(2,8)); s=BTZ(rng,n)
    lab=[int(rng.integers(0,N+1)) for _ in range(n)]; lab2=int(rng.integers(0,N+1))
    regs=[frozenset(c) for r in range(1,N+1) for c in itertools.combinations(range(1,N+1),r)]
    v={}
    for A in regs:
        inA=[1 if lab[k] in A else 0 for k in range(n)]; b2=1 if lab2 in A else 0
        v[A]=F(s,inA,b2,1e-6)[0]
    sc=max(abs(x) for x in v.values())
    for name,d in fs:
        q=sum(c*v[A] for A,c in d.items())/sc
        if q<-1e-9:
            print(t,n,lab,lab2,name,q,{tuple(sorted(A)):round(x,3) for A,x in v.items()}); break
    else: continue
    if t>40: break

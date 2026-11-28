# Does the finite part G change when the horoball normalisation at one marked point is rescaled?
# SM Sec. S5 "Cutoff scheme" remark says G shifts by a multiple of b_p; Prop. S4.3 implies G is invariant.
import sys; sys.path.insert(0,'.')
from btz_cone_test import *
rng=np.random.default_rng(5); N=3; maxdG=0; maxdF=0
for t in range(200):
    n=int(rng.integers(2,7)); s=BTZ(rng,n)
    lab=[int(rng.integers(0,N+1)) for _ in range(n)]; lab2=int(rng.integers(0,N+1))
    regs=[frozenset(c) for r in range(1,N+1) for c in itertools.combinations(range(1,N+1),r)]
    def GF():
        out={}
        for A in regs:
            inA=[1 if lab[k] in A else 0 for k in range(n)]; b2=1 if lab2 in A else 0
            f,_=F(s,inA,b2,1e-6); E=[p for p in range(n) if inA[(p-1)%n]!=inA[p]]
            out[A]=(f, f-sum(math.log(1e6)-0.5*math.log(s.dp(p,'L')*s.dp(p,'R')) for p in E))
        return out
    g0=GF(); k=int(rng.integers(0,n)); s.cL[k]*=rng.uniform(0.1,10); s.cR[k]*=rng.uniform(0.1,10); g1=GF()
    for A in regs:
        maxdG=max(maxdG,abs(g1[A][1]-g0[A][1])); maxdF=max(maxdF,abs(g1[A][0]-g0[A][0]))
print(f"max change of F^eps under rescaling one horoball: {maxdF:.3f};  max change of G: {maxdG:.1e}")

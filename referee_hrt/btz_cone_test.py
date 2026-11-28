# Independent referee check: HRT entropies in ROTATING BTZ (pi_1 = Z), regions = unions of intervals on
# boundary 1 plus (optionally) the whole boundary 2.  HRT = min over common cut classes of (l_L+l_R)/2
# (Theorem H of the submission; for pi_1=Z this is also GGHH's working formula).
# Checks: (1) all C_3 facets / SSA / MMI / 5-party cyclic on F^eps at small eps;
#         (2) G = F - sum_p t_p b_p is eps-independent and satisfies the UNBALANCED facets too (Thm 2(i));
#         (3) control: half endpoint subtraction; (4) how often common-class minimiser != separate minimisers.
import itertools, math, random, sys
import numpy as np

def matchings(pts):
    if not pts: yield []; return
    a=pts[0]
    for k in range(1,len(pts),2):        # non-crossing: partner splits list into two even blocks
        b=pts[k]
        for m1 in matchings(pts[1:k]):
            for m2 in matchings(pts[k+1:]):
                yield [(a,b)]+m1+m2

class BTZ:
    def __init__(s, rng, n):
        s.a = rng.uniform(1.0,6.0); s.b = rng.uniform(1.0,6.0)            # chiral horizon lengths l_L != l_R: rotating
        s.n = n
        s.th = np.sort(rng.uniform(0,s.a,n))
        amp = 0.8*s.b/(2*math.pi); ph=rng.uniform(0,2*math.pi)          # wiggly boundary slice: phi increasing, equivariant
        s.phi = lambda t: (s.b/s.a)*t + amp*math.sin(2*math.pi*t/s.a+ph)
        s.cL = rng.uniform(0.3,3,n); s.cR = rng.uniform(0.3,3,n)        # random horocyclic normalisations
    def arc(s,i,j,eps):   # arc homotopic to ccw boundary interval from point i to point j; returns (lL, lR)
        ti=s.th[i]; tj=s.th[j] if s.th[j]>ti else s.th[j]+s.a
        x1,x2=math.exp(ti),math.exp(tj)
        lL=math.log((x2-x1)**2/((eps*x1*s.cL[i])*(eps*x2*s.cL[j])))
        y1,y2=math.exp(s.phi(ti)),math.exp(s.phi(tj))
        lR=math.log((y2-y1)**2/((eps*y1*s.cR[i])*(eps*y2*s.cR[j])))
        return lL,lR
    def dp(s,p,side):     # d_p of Sec. S4 in the chart with p at infinity
        n=s.n; t=s.th[p]; tm=s.th[(p-1)%n]-(s.a if p==0 else 0); tp=s.th[(p+1)%n]+(s.a if p==n-1 else 0)
        if n==1: tm=t-s.a; tp=t+s.a
        f=(lambda x:x) if side=='L' else s.phi
        c=s.cL[p] if side=='L' else s.cR[p]
        X=math.exp(f(t)); Xm=math.exp(f(tm)); Xp=math.exp(f(tp))
        return X*c*(1/(X-Xm)+1/(Xp-X))

def pieces_of(s):  # piece k between point k and k+1 (ccw); interval i->j covers pieces i..j-1
    return list(range(s.n))
def interval_pieces(n,i,j):
    out=[]; k=i
    while k!=j: out.append(k); k=(k+1)%n
    return frozenset(out)

def F(s, inA, bd2A, eps, mode='common'):
    n=s.n
    E=[p for p in range(n) if inA[(p-1)%n]!=inA[p]]
    best=math.inf; bestL=math.inf; bestR=math.inf; argc=None; argL=None; argR=None
    for M in matchings(E):
        for orient in itertools.product((0,1),repeat=len(M)):
            ivs=[]; L=0; R=0
            for (p,q),o in zip(M,orient):
                i,j=(p,q) if o==0 else (q,p)
                ivs.append(interval_pieces(n,i,j)); lL,lR=s.arc(i,j,eps); L+=lL; R+=lR
            ok=True
            for u,v in itertools.combinations(ivs,2):
                if not (u<=v or v<=u or not (u&v)): ok=False;break
            if not ok: continue
            for c0 in (0,1):
                col=[(c0+sum(1 for I in ivs if k in I))%2 for k in range(n)]
                if any(col[k]!=inA[k] for k in range(n)): continue
                h=(c0+bd2A)%2
                tl=L+h*s.a; tr=R+h*s.b
                key=(tuple(sorted(M)),orient,c0)
                val=(tl+tr)/2 if mode=='common' else max(tl,tr)
                if val<best: best=val; argc=key
                if tl<bestL: bestL=tl; argL=key
                if tr<bestR: bestR=tr; argR=key
    return best, (argc!=argL or argc!=argR)

def facets(N):
    # returns list of (name, dict region->coef, LHS-RHS >=0), regions as frozensets of parties 1..N (purifier = complement)
    P=list(range(1,N+1)); out=[]
    def S(*xs): return frozenset(xs)
    if N==3:
        A,B,C=1,2,3
        out.append(('MMI',{S(A,B):1,S(B,C):1,S(A,C):1,S(A):-1,S(B):-1,S(C):-1,S(A,B,C):-1}))
        # SA instances among A,B,C,O
        out.append(('SA AB',{S(A):1,S(B):1,S(A,B):-1})); out.append(('SA AC',{S(A):1,S(C):1,S(A,C):-1}))
        out.append(('SA BC',{S(B):1,S(C):1,S(B,C):-1}))
        out.append(('SA AO',{S(A):1,S(A,B,C):1,S(B,C):-1})); out.append(('SA BO',{S(B):1,S(A,B,C):1,S(A,C):-1}))
        out.append(('SA CO',{S(C):1,S(A,B,C):1,S(A,B):-1}))
        out.append(('SSA',{S(A,B):1,S(B,C):1,S(B):-1,S(A,B,C):-1}))
    if N==5:
        for perm in itertools.permutations(P):
            a=perm; d={}
            for i in range(5):
                k=S(a[i],a[(i+1)%5],a[(i+2)%5]); d[k]=d.get(k,0)+1
                k=S(a[i],a[(i+1)%5]); d[k]=d.get(k,0)-1
            d[S(*P)]=d.get(S(*P),0)-1
            out.append(('cyc5',d))
        for t in itertools.permutations(P,3):
            A,B,C=t
            out.append(('MMI',{S(A,B):1,S(B,C):1,S(A,C):1,S(A):-1,S(B):-1,S(C):-1,S(A,B,C):-1}))
    return out

def run(N, trials, seed):
    rng=np.random.default_rng(seed); random.seed(seed)
    fs=facets(N); worst={'F':math.inf,'G':math.inf,'ctl':math.inf,'Gdrift':0}; nctl=0; nctl2=0; ndiff=0; nreg=0
    for t in range(trials):
        n=int(rng.integers(2,8)); s=BTZ(rng,n)
        lab=[int(rng.integers(0,N+1)) for _ in range(n)]   # 0 = purifier
        lab2=int(rng.integers(0,N+1))
        tp={}
        regs=[frozenset(c) for r in range(1,N+1) for c in itertools.combinations(range(1,N+1),r)]
        vals={}
        for eps in (1e-3,1e-6):
            for A in regs:
                inA=[1 if lab[k] in A else 0 for k in range(n)]; b2=1 if lab2 in A else 0
                f,diff=F(s,inA,b2,eps)
                E=[p for p in range(n) if inA[(p-1)%n]!=inA[p]]
                t_sum=sum(math.log(1/eps)-0.5*math.log(s.dp(p,'L')*s.dp(p,'R')) for p in E)
                fm=F(s,inA,b2,eps,mode='minmax')[0]
                vals[(eps,A)]=(f,f-t_sum,f-1.5*t_sum,fm)   # F, G, control1 (over-subtraction), control2 (min_c max(lL,lR))
                if eps==1e-6: nreg+=1; ndiff+=diff
        scale=max(abs(vals[(1e-6,A)][0]) for A in regs)+1e-12
        if scale<1e-9: continue
        for A in regs: worst['Gdrift']=max(worst['Gdrift'],abs(vals[(1e-3,A)][1]-vals[(1e-6,A)][1]))
        bad=False; bad2=False
        for name,d in fs:
            qF=sum(c*vals[(1e-6,A)][0] for A,c in d.items())/scale
            qG=sum(c*vals[(1e-6,A)][1] for A,c in d.items())/scale
            qC=sum(c*vals[(1e-6,A)][2] for A,c in d.items())/scale
            worst['F']=min(worst['F'],qF); worst['G']=min(worst['G'],qG); worst['ctl']=min(worst['ctl'],qC)
            qM=sum(c*vals[(1e-6,A)][3] for A,c in d.items())/scale
            if qC<-1e-9: bad=True
            if qM<-1e-9: bad2=True
        nctl+=bad; nctl2+=bad2
    print(f"N={N} trials={trials}: min facet F^1e-6 = {worst['F']:.2e}, min facet G = {worst['G']:.2e}, "
          f"G eps-drift = {worst['Gdrift']:.1e}, controls violating: over-subtraction {nctl}/{trials}, min-max {nctl2}/{trials}, "
          f"regions with common-class minimiser != a separate minimiser: {ndiff}/{nreg}")

if __name__=='__main__':
    run(3, int(sys.argv[1]) if len(sys.argv)>1 else 300, 11)
    run(5, int(sys.argv[2]) if len(sys.argv)>2 else 60, 12)

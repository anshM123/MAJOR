# Cost model for Taylor-triangle sign certificates (as in sliver_direct_cont/direct_sign_triangle.cpp)
# under three parametrisations of a cell. Majorant v = absolute Horner majorant of the prefix on the
# complex polydisk (radii 2^ac, 2^ah in rescaled coordinates); Q = least weighted degree with
# 2 v (Q+2) 2^-Q <= 1e-10; terms = #{(i,j): ac i + ah j < Q}; est. time = 1.3e-7 s * K * terms
# (calibrated on the K=1397..1567 runs in code/sliver_direct_ext/runlogs and the author's next_* runs).
#  A: h-rectangle  c in [cl,cu], h in [hl,hu]  (G1 cap; author's driver with h0!=dh allowed)
#  C: curvilinear  r=sqrt(c) in [rl,ru], y in [yl,yu], h = y (1-r)^2/998^2
#  B: x-box        c in [cl,cu], x=h/c in [xl,xu]
import sys,math
from mpmath import mp,mpf,log10
mp.dps=20
from common2 import parts,needK,arb
def amag(z):return mpf(z.mid().str(radius=False)) if hasattr(z,'mid') else mpf(z)
_cache={}
def mags(b,K):
    if (b,K) not in _cache:
        A,Bs=parts(b,K);q=b*(b+1)
        _cache[(b,K)]=([amag(a) for a in A],[amag(q*a+t) for a,t in zip(A,Bs)])
    return _cache[(b,K)]
def polymul(P,Q):
    R={}
    for (a,b),x in P.items():
        for (c,d),y in Q.items():R[(a+c,b+d)]=R.get((a+c,b+d),0)+x*y
    return R
def pmaj(P,ru,rv,skip0=False):
    return sum(abs(x)*ru**a*rv**b for (a,b),x in P.items() if not(skip0 and a==0 and b==0))
def majorant(param,cell,b,K,ac,ah):
    A,T=mags(b,K);s=2*b+1;ru=mpf(2)**ac;rv=mpf(2)**ah
    lo1,hi1,lo2,hi2=[mpf(x) for x in cell]
    m1=(lo1+hi1)/2;d1=(hi1-lo1)/2;m2=(lo2+hi2)/2;d2=(hi2-lo2)/2
    if param=='A':   # c=m1+d1 u, h=m2+d2 v
        c={(0,0):m1,(1,0):d1};h={(0,0):m2,(0,1):d2}
    elif param=='C': # r=m1+d1 u, y=m2+d2 v
        r={(0,0):m1,(1,0):d1};omr={(0,0):1-m1,(1,0):-d1};y={(0,0):m2,(0,1):d2}
        c=polymul(r,r);h={k:v/998**2 for k,v in polymul(y,polymul(omr,omr)).items()}
    elif param=='B': # c=m1+d1 u, x=m2+d2 v, h=c x
        c={(0,0):m1,(1,0):d1};h=polymul(c,{(0,0):m2,(0,1):d2})
    h00=h.get((0,0),0);hr=pmaj(h,ru,rv,True);cm=pmaj(c,ru,rv)
    v=mpf(0)
    for k in range(K,-1,-1):
        Fk=abs(A[k]+T[k]*h00)+abs(T[k])*hr
        if k<K:
            p=k*(k+s)
            num={kk:-p*x for kk,x in h.items()}
            for kk,x in c.items():num[kk]=num.get(kk,0)+x
            dd=1-p*h00-p*hr
            if dd<=0:return None
            v=Fk+v*pmaj(num,ru,rv)/dd
        else:v=Fk
    return v
def best(param,cell,b,K,acs=range(1,9),ahs=range(1,8)):
    bestv=None
    for ac in acs:
        for ah in ahs:
            v=majorant(param,cell,b,K,ac,ah)
            if v is None:continue
            Q=max(ac,ah)+1;tgt=mpf(10)**-10
            Q=max(ac,ah)+1
            lv=float(log10(v))
            while lv+math.log10(2*(Q+2))-Q*math.log10(2)>-10:Q+=1
            terms=sum((Q-1-ac*i)//ah+1 for i in range((Q-1)//ac+1))
            if bestv is None or terms<bestv[0]:bestv=(terms,ac,ah,Q,lv)
    if bestv is None:return None
    terms,ac,ah,Q,lv=bestv
    return dict(param=param,cell=[float(x) for x in cell],b=b,K=K,ac=ac,ah=ah,Q=Q,log10v=round(lv,1),terms=terms,est_s=round(1.3e-7*K*terms,1))
if __name__=='__main__':
    param=sys.argv[1];cell=[mpf(x) for x in sys.argv[2:6]];b=int(sys.argv[6]);K=int(sys.argv[7])
    print(best(param,cell,b,K))

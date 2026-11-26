"""c-direction approximation on bands of the *ideal* y norming nodes."""
import sys,os
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'sliver_close'))
from coupled_plan import *
from reference_driver import cos_nodes, NODE_DEN
from functools import lru_cache

def group_bounds(l,u,Hl,Hu,K,ylo,yhi):
    out=[];ps=p[:K];c0=(l+u)/2;dc=(u-l)/2;H0=(Hl+Hu)/2;dH=(Hu-Hl)/2
    for rho in Rhos:
        R=(rho+1/rho)/2
        pairs=[(c0+z*dc,(H0+z*dH)*y) for z in [-R,R] for y in (ylo,yhi)]
        HH=max(abs(h) for c,h in pairs);norm=max(abs(1-c) for c,h in pairs)**3
        if max(h for c,h in pairs)*(K-1)*(K+8)>=1:continue
        f=np.zeros(K)
        for pp in [ps,k[:K]*(k[:K]+1)]:
            for c,h in pairs:f=np.maximum(f,np.abs(float(c)-pp*float(h))/(1-pp*float(h)))
        M=np.r_[0.,np.cumsum(np.log(np.maximum(f,1e-300)))]+np.logaddexp(A[:K+1],T[:K+1]+math.log(max(1e-300,float(HH))))+math.log(float(norm))+math.log(2/float(rho-1))
        out.append((rho,M))
    rho=np.array([math.log(float(r)) for r,v in out]);M=np.array([v for r,v in out]);
    def err(m):return logsumexp(np.min(M-m*rho[:,None],axis=0))
    lo,hi=1,500
    while err(hi)>math.log(1e-10):hi*=2
    while hi-lo>1:
        m=(hi+lo)//2
        if err(m)>math.log(1e-10):lo=m
        else:hi=m
    m=hi+2;choice=np.argmin(M-m*rho[:,None],axis=0)
    used=sorted(set(choice));rs=[out[i][0] for i in used];ii={v:i for i,v in enumerate(used)}
    return m,rs,[ii[i] for i in choice]

def group_check(l,u,Hl,Hu,K,ylo,yhi,m,rs,choices):
    # positive upper-rounded Decimal products; ratio maxima exact rational,
    # computed from integer corner numerator/denominator without Fraction gcds.
    with localcontext() as ctx:
        ctx.prec=45;ctx.rounding=ROUND_CEILING
        def dec(f):f=F(f);return D(f.numerator)/D(f.denominator)
        bounds=[None]*(K+1)
        c0=(l+u)/2;dc=(u-l)/2;H0=(Hl+Hu)/2;dH=(Hu-Hl)/2
        for idx,rho in enumerate(rs):
            R=(rho+1/rho)/2;pairs=[(c0+z*dc,(H0+z*dH)*y) for z in [-R,R] for y in (ylo,yhi)]
            HH=max(abs(h) for c,h in pairs);norm=max(abs(1-c) for c,h in pairs)**3
            assert max(h for c,h in pairs)*(K-1)*(K+8)<1
            coef=dec(2*rho**(-m)/(rho-1)*norm);prod=D(1);hd=dec(HH)
            # common denominators c/h
            coeffs=[(c.numerator*h.denominator,h.numerator*c.denominator,c.denominator*h.denominator) for c,h in pairs]
            last=max(j for j,t in enumerate(choices) if t==idx)
            for j in range(last+1):
                if choices[j]==idx:bounds[j]=prod*(D(AI[j])/D(S)+hd*D(TI[j])/D(S))*coef
                if j<last:
                    n0,d0=0,1
                    for p0 in (j*(j+1),j*(j+9)):
                        for cn,hn,dd in coeffs:
                            n1=abs(cn-p0*hn);d1=dd-p0*hn;assert d1>0
                            if n1*d0>n0*d1:n0,d0=n1,d1
                    prod*=D(n0)/D(d0)
        assert all(x is not None for x in bounds)
        total=sum(bounds,D(0));assert total<D('1e-10'),total
        return str(total)

def make_group(path):
    path=Path(path);data=json.load(open(path));l,u,Hl,Hu=map(F,(data['l'],data['u'],data['Hl'],data['Hu']));K=data['K'];my=data['y']['deg'];yy=cos_nodes(my);groups=[]
    # use [0,.45], [.45,.7], [.7,.9], [.9,1]. Group counts checked.
    edges=[F(0),F('.45'),F('.7'),F('.9'),F(1)]
    for ylo,yhi in zip(edges,edges[1:]):
        vals=[j for j in range(my+1) if ylo<=F(yy[j],NODE_DEN)<=yhi and not any(j in g['nodes'] for g in groups)]
        if not vals:continue
        # Ideal cos nodes within 2e-24 of stored rational nodes; clip at true [0,1].
        lo=max(F(0),min(F(yy[j],NODE_DEN) for j in vals)-F(2,10**24));hi=min(F(1),max(F(yy[j],NODE_DEN) for j in vals)+F(2,10**24))
        m,rs,ch=group_bounds(l,u,Hl,Hu,K,lo,hi);e=group_check(l,u,Hl,Hu,K,lo,hi,m,rs,ch)
        groups.append(dict(nodes=vals,lo=str(lo),hi=str(hi),deg=m,rhos=list(map(str,rs)),choices=ch,error=e))
        print(path.stem,'group',len(vals),'degree',m,'err',e,flush=True)
    assert sorted(j for g in groups for j in g['nodes'])==list(range(my+1))
    cnt=2*sum(len(g['nodes'])*(1+g['deg']) for g in groups);print('nodes',cnt,'vs',2*(data['c']['deg']+1)*(my+1),flush=True)
    dest=path.with_name(path.stem+'_group.json');json.dump(dict(cell=path.name,groups=groups,nodes=cnt),open(dest,'w'))
if __name__=='__main__':
    for path in sys.argv[1:]:make_group(path)

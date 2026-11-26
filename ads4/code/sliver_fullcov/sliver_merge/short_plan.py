"""A rigorously shortened scalar prefix per twist group."""
from group_plan import *
from reference_driver import DEN

def tail_vector(l,u,Hl,Hu,K,lo,hi):
    epsy=F(2,10**24)+F(1,DEN)/Hu
    yl,yh=max(F(0),lo-epsy),min(F(1),hi+epsy)
    corners=[(c,H*y) for c,H in [(l,Hl),(u,Hu)] for y in (yl,yh)]
    coeffs=[(c.numerator*h.denominator,h.numerator*c.denominator,c.denominator*h.denominator) for c,h in corners]
    with localcontext() as ctx:
        ctx.prec=50;ctx.rounding=ROUND_CEILING
        H=Hl*yh;hd=D(H.numerator)/D(H.denominator);prod=D(1);terms=[]
        for j in range(K+1):
            terms.append(prod*(D(AI[j])/D(S)+hd*D(TI[j])/D(S)))
            if j<K:
                n0,d0=0,1
                for pp in (j*(j+1),j*(j+9)):
                    for cn,hn,dd in coeffs:
                        n1=abs(cn-pp*hn);d1=dd-pp*hn;assert d1>0 and n1<=d1
                        if n1*d0>n0*d1:n0,d0=n1,d1
                prod*=D(n0)/D(d0)
        tail=[D(0)]*(K+1)
        for j in range(K-1,-1,-1):tail[j]=tail[j+1]+terms[j+1]
        return tail

def make_short(path):
    path=Path(path);z=json.load(open(path));groups=json.load(open(path.with_name(path.stem+'_group.json')))['groups'];l,u,Hl,Hu=map(F,(z['l'],z['u'],z['Hl'],z['Hu']));K=z['K'];new=[]
    for g in groups:
        lo,hi=F(g['lo']),F(g['hi']);tail=tail_vector(l,u,Hl,Hu,K,lo,hi);Kg=next(j for j in range(1,K+1) if tail[j]<D('1e-10'));tailg=F(tail[Kg]);m,rs,ch=group_bounds(l,u,Hl,Hu,Kg,lo,hi);ec=F(group_check(l,u,Hl,Hu,Kg,lo,hi,m,rs,ch));e=tailg+ec
        new.append(dict(nodes=g['nodes'],lo=g['lo'],hi=g['hi'],K=Kg,short_tail=str(tailg),deg=m,rhos=list(map(str,rs)),choices=ch,poly_error=str(ec),error=str(e)))
        assert e<F('2e-10');print(path.stem,'K',Kg,'degree',m,'group',len(g['nodes']),'err',float(e),flush=True)
    cnt=2*sum(len(g['nodes'])*(1+g['deg']) for g in new);print('nodes',cnt,'vs',2*(z['c']['deg']+1)*(z['y']['deg']+1),flush=True)
    json.dump(dict(cell=path.name,groups=new,nodes=cnt),open(path.with_name(path.stem+'_short_group.json'),'w'))
if __name__=='__main__':
    for path in sys.argv[1:]:make_short(path)

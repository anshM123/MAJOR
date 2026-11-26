"""Degree-matched exact trigonometric norming with a rational reference."""
from fractions import Fraction as F
from math import isqrt,ceil,comb
from pathlib import Path
import json,sys,subprocess,concurrent.futures,time,os
from norm_driver_scaled import PI0,NODE_DEN,DEN
HERE=Path(__file__).resolve().parent
COEFF=Path(os.environ.get('SLIVER_COEFF',str(HERE/'coeff_compact60.txt')))
UNIT=10**20

def cos_nodes(m):
    # Restriction of the 2m+1 equispaced trigonometric grid to the even half.
    v=[]
    for j in range(m+1):
        x=PI0*2*j/(2*m+1);z=-x*x;term=F(1);s=term
        for k in range(1,27):term*=z/((2*k-1)*(2*k));s+=term
        rem=abs(term*z/(53*54));assert rem<F(1,10**40)
        y=(1+s)/2;u=min(NODE_DEN,max(0,round(y*NODE_DEN)))
        assert abs(F(u,NODE_DEN)-y)<=F(1,NODE_DEN)
        v.append(u)
    assert v[0]==NODE_DEN and all(v[i]>v[i+1] for i in range(m))
    return v

def norm_constant(m,eps):
    # Dirichlet interpolation with 2m+1 angular samples. For the nearest
    # point use 1, for each of the other m pairs use 1/(2k-1).
    lam=1+2*sum((F(1,2*k-1) for k in range(1,m+1)),F(0))
    # Markov perturbation on x in [0,1]: ||p'|| <= 2m^2 ||p||.
    assert 1-2*m*m*eps*lam>0
    return lam/(1-2*m*m*eps*lam)

def reference(data):
    l,u,Hl,Hu=map(F,(data['l'],data['u'],data['Hl'],data['Hu']));K=data['K'];
    pts=[];keys=[]
    for b in (0,4):
        for i in range(3):
            for j in range(3):
                c0=l+(u-l)*F(i,2);C=max(ceil(l*DEN),min(int(u*DEN),int(c0*DEN)));c=F(C,DEN)
                H=Hl+(Hu-Hl)*(c-l)/(u-l);h=int(H*F(j,2)*DEN)
                pts.append(f'{b} {K} {C} {h} {DEN}\n');keys.append((b,i,j,c))
    # Only a convenient choice of reference. Its coefficients are rational and
    # no accuracy assumption about this fit is made in the certificate.
    fast=HERE/'fast_scalar'
    proc=subprocess.run([str(fast),str(COEFF)],input=''.join(pts),capture_output=True,text=True,check=True)
    vals={}
    for (b,i,j,c),v in zip(keys,proc.stdout.splitlines()):vals[b,i,j]=F(int(v),10**60)*(1-c)**3
    basis=((1,-3,2),(0,4,-4),(0,-1,2)) # Lagrange polynomials at 0,1/2,1
    ref=[];low=[]
    for b in (0,4):
        P=[[round(UNIT*sum(vals[b,i,j]*basis[i][a]*basis[j][d] for i in range(3) for j in range(3))) for d in range(3)] for a in range(3)]
        # Exact Bernstein coefficients (degree 2 separately in x,y).
        B=[[sum(F(P[a][d]*comb(i,a)*comb(j,d),comb(2,a)*comb(2,d)*UNIT) for a in range(i+1) for d in range(j+1)) for j in range(3)] for i in range(3)]
        lb=min(min(row) for row in B)
        assert lb>0
        assert sum(abs(x) for row in P for x in row)<=UNIT//1000
        ref.append(P);low.append(lb)
    return ref,low

def make_cell(path,execute=True):
    p=Path(path).resolve();data=json.load(open(p));l,u,Hl,Hu=map(F,(data['l'],data['u'],data['Hl'],data['Hu']));K=data['K'];mc=data['c']['deg'];my=data['y']['deg'];assert l<u and 0<Hu<=Hl and min(mc,my)>=2
    for c,h in [(l,Hl),(u,Hu)]:
        z=F(isqrt(c.numerator*10**44//c.denominator),10**22);assert z*z<=c and h>=(1-z)**2/998**2
    for s in (1,9):
        for j in range(K):
            pp=j*(j+s)
            for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)]:assert 1-pp*h>0 and abs(c-pp*h)<=1-pp*h
    sc=2**100;upper=sc;shifts=[0];sh=0
    for j in range(K):
        mx=max(abs(c-h*j*(j+s))/(1-h*j*(j+s)) for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)] for s in (1,9))
        assert mx<=1;upper=ceil(upper*mx)
        t=0
        while upper<=sc//2 and t<2:upper*=2;sh+=1;t+=1
        assert sc//2<upper<=sc;shifts.append(sh)
    err=(K+1)*(F(4)+Hl*((K+1)**2+K*K+8*K+67))/10**60+F(10*(K+1),2**180);assert err<F(1,10**45)
    # Actual point in [0,1] is a rational rounding of the angular node.
    Cc=norm_constant(mc,F(2,10**24)+F(1,DEN)/(u-l));Cy=norm_constant(my,F(2,10**24)+F(1,DEN)/Hu)
    assert Cc<12 and Cy<12
    ref,ref_lower=reference(data)
    # Degree-2 tensor Horner, with <=1 unit truncation per operation and
    # x,y rounded down within 2^-50. Bound via sum |coeff|<=10^-3 and
    # |partial_x R|,|partial_y R|<=2*10^-3 each.
    referr=F(4,1000*2**50)+F(8,UNIT)
    assert referr<F(1,10**15)
    yc=cos_nodes(mc);yy=cos_nodes(my);rows=[]
    for a in yc:
        c0=l+(u-l)*F(a,NODE_DEN);C=max(ceil(l*DEN),min(int(u*DEN),int(c0*DEN)));c=F(C,DEN);assert l<=c<=u and abs(c-c0)<=F(1,DEN)
        H=Hl+(Hu-Hl)*(c-l)/(u-l);rows.append((C,H.numerator,H.denominator))
    assert len(set(C for C,_,_ in rows))==len(rows)
    wd=p.parent;stem=p.stem;parts=[]
    refs=' '.join(map(str,[l.numerator,l.denominator,u.numerator,u.denominator]+[x for P in ref for row in P for x in row]))+'\n'
    for worker in range(4):
        r=rows[worker::4];inp=wd/f'{stem}_r{worker}.in';out=wd/f'{stem}_r{worker}.out';log=wd/f'{stem}_r{worker}.log';parts.append((inp,out,log,len(r)))
        inp.write_text(f'{K} {DEN} {len(yy)} {len(r)} 180\n'+' '.join(map(str,shifts))+'\n'+f'{NODE_DEN}\n'+ ' '.join(map(str,yy))+'\n'+refs+''.join(f'{C} {Hn} {Hd}\n' for C,Hn,Hd in r))
    print(stem,'grid',mc+1,my+1,'nodes',2*(mc+1)*(my+1),'norm constants',float(Cc),float(Cy),'ref lower',list(map(float,ref_lower)),flush=True)
    if not execute:return
    def worker(part):
        inp,out,log,num=part
        if out.exists() and 'COUNT ' in out.read_text():return out
        with open(inp)as f,open(out,'w')as g,open(log,'w')as h:subprocess.run([str(HERE/'scalar_ref'),str(COEFF)],stdin=f,stdout=g,stderr=h,check=True)
        return out
    st=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as exe:
        for res in exe.map(worker,parts):print('done',res.name,flush=True)
    maxabs={0:0,4:0};count=0
    for inp,out,log,num in parts:
        seen=set();ct=False
        for line in out.read_text().splitlines():
            a,*z=line.split()
            if a=='COUNT':assert not ct and int(z[0])==num*len(yy)*2;ct=True;count+=int(z[0]);continue
            b=int(a);assert b in (0,4) and b not in seen;seen.add(b);v,w=map(int,z);assert v<=w;maxabs[b]=max(maxabs[b],abs(v),abs(w))
        assert seen=={0,4} and ct
    assert count==2*(mc+1)*(my+1)
    result=[]
    for i,b in enumerate((0,4)):
        delta=F(maxabs[b],UNIT)+err+referr;E_c=F(data['c']['error']);E_y=F(data['y']['error'])
        residual=Cc*(Cy*delta+(Cy+1)*E_y)+(Cc+1)*E_c
        lower=ref_lower[i]-residual
        Hphys=lower/(1+b*(b+1)*Hl)**4-F(data['physical_tail'])-F(1,2**140)-F(1,10**14)
        result.append(dict(b=b,node_residual=str(delta),ref_lower=str(ref_lower[i]),uniform_residual=str(residual),prefix_lower=str(lower),physical_lower=str(Hphys),lower_float=float(Hphys),status='PASS_CONDITIONAL' if Hphys>0 else 'INCONCLUSIVE'))
    print(stem,[(z['b'],z['lower_float'],z['status']) for z in result],'seconds',time.monotonic()-st,flush=True)
    json.dump(dict(cell=stem,count=count,Nc=mc,Ny=my,Cc=str(Cc),Cy=str(Cy),scalar_error=str(err),reference_error=str(referr),reference=ref,results=result),open(wd/f'{stem}_ref_ledger.json','w'),indent=2)
    assert all(z['status']=='PASS_CONDITIONAL' for z in result)
if __name__=='__main__':
    for p in sys.argv[1:]:make_cell(p)

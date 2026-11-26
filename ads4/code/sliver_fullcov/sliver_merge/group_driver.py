"""Rigorous grouped c-norming, with global self-consistent y-rounding bound."""
from reference_driver import *
import hashlib

def make_cell_group(path,execute=True):
    p=Path(path).resolve();data=json.load(open(p));groups=json.load(open(p.with_name(p.stem+'_group.json')))['groups'];l,u,Hl,Hu=map(F,(data['l'],data['u'],data['Hl'],data['Hu']));K=data['K'];my=data['y']['deg'];assert l<u and 0<Hu<=Hl and my>=2
    for c,h in [(l,Hl),(u,Hu)]:
        z=F(isqrt(c.numerator*10**44//c.denominator),10**22);assert z*z<=c and h>=(1-z)**2/998**2
    for s in (1,9):
        for j in range(K):
            pp=j*(j+s)
            for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)]:assert 1-pp*h>0 and abs(c-pp*h)<=1-pp*h
    P=80
    err=(K+1)*(F(4)+Hl*((K+1)**2+K*K+8*K+67))/10**60+F(10*(K+1),2**P);assert err<F(1,10**15)
    epsc=F(2,10**24)+F(1,DEN)/(u-l);epsy=F(2,10**24)+F(1,DEN)/Hu
    Cy=norm_constant(my,F(0));gamma=2*my**2*epsy
    assert Cy<12
    ref,ref_lower=reference(data);referr=F(4,1000*2**50)+F(8,UNIT)
    yy=cos_nodes(my);parts=[];wd=p.parent;stem=p.stem;Cc=[]
    refs=' '.join(map(str,[l.numerator,l.denominator,u.numerator,u.denominator]+[x for P in ref for row in P for x in row]))+'\n'
    assert sorted(j for g in groups for j in g['nodes'])==list(range(my+1))
    for ig,g in enumerate(groups):
        lo,hi=F(g['lo']),F(g['hi']);m=g['deg'];assert m>=2;Ccg=norm_constant(m,epsc);assert Ccg<12;Cc.append(Ccg)
        yg=[yy[j] for j in g['nodes']];assert yg
        for v in yg:assert lo<=max(F(0),F(v,NODE_DEN)-F(2,10**24)) and hi>=min(F(1),F(v,NODE_DEN)+F(2,10**24))
        # A separate exact prefix scaling schedule for this y group.
        ylo=max(F(0),lo-epsy);yhi=min(F(1),hi+epsy)
        corners=[(c,H*y) for c,H in [(l,Hl),(u,Hu)] for y in (ylo,yhi)]
        sc=2**100;upper=sc;shifts=[0];sh=0
        for j in range(K):
            mx=max(abs(c-h*j*(j+s))/(1-h*j*(j+s)) for c,h in corners for s in (1,9));assert mx<=1;upper=ceil(upper*mx);t=0
            while upper<=sc//2 and t<2:upper*=2;sh+=1;t+=1
            assert sc//2<upper<=sc;shifts.append(sh)
        rows=[]
        for a in cos_nodes(m):
            c0=l+(u-l)*F(a,NODE_DEN);C=max(ceil(l*DEN),min(int(u*DEN),int(c0*DEN)));c=F(C,DEN);assert l<=c<=u and abs(c-c0)<=F(1,DEN)
            H=Hl+(Hu-Hl)*(c-l)/(u-l);rows.append((C,H.numerator,H.denominator))
        assert len(set(C for C,_,_ in rows))==len(rows)
        for w in range(4):
            r=rows[w::4];inp=wd/f'{stem}_g{ig}_{w}.in';out=wd/f'{stem}_g{ig}_{w}.out';log=wd/f'{stem}_g{ig}_{w}.log';parts.append((ig,inp,out,log,len(r)*len(yg)*2))
            payload=f'{K} {DEN} {len(yg)} {len(r)} {P}\n'+' '.join(map(str,shifts))+'\n'+f'{NODE_DEN}\n'+ ' '.join(map(str,yg))+'\n'+refs+''.join(f'{C} {Hn} {Hd}\n' for C,Hn,Hd in r)
            if inp.exists() and inp.read_text()!=payload:
                assert not out.exists() or 'COUNT ' not in out.read_text(),f'changed completed input {inp}'
            inp.write_text(payload)
    count_expected=2*sum(len(g['nodes'])*(1+g['deg']) for g in groups);assert count_expected==sum(z[-1] for z in parts)
    print(stem,'group grid',[(g['deg']+1,len(g['nodes'])) for g in groups],'nodes',count_expected,'ref lower',list(map(float,ref_lower)),flush=True)
    if not execute:return
    def worker(part):
        ig,inp,out,log,num=part
        if out.exists() and 'COUNT ' in out.read_text():return out
        with open(inp)as f,open(out,'w')as g,open(log,'w')as h:subprocess.run([str(HERE/'scalar_ref80'),str(COEFF)],stdin=f,stdout=g,stderr=h,check=True)
        return out
    st=time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as exe:
        # Put the largest tasks first to balance worker time.
        for res in exe.map(worker,sorted(parts,key=lambda z:-z[-1])):print('done',res.name,flush=True)
    maxabs=[{0:0,4:0} for g in groups];count=0;inputs={}
    for ig,inp,out,log,num in parts:
        seen=set();ct=False
        for line in out.read_text().splitlines():
            a,*z=line.split()
            if a=='COUNT':assert not ct and int(z[0])==num;ct=True;count+=int(z[0]);continue
            b=int(a);assert b in (0,4) and b not in seen;seen.add(b);v,w=map(int,z);assert v<=w;maxabs[ig][b]=max(maxabs[ig][b],abs(v),abs(w))
        assert seen=={0,4} and ct;inputs[inp.name]=hashlib.sha256(inp.read_bytes()).hexdigest()
    assert count==count_expected
    result=[];Ey=F(data['y']['error']);lam=Cy*max(Cc)*gamma;assert lam<F(1,10)
    for i,b in enumerate((0,4)):
        delta=[F(z[b],UNIT)+err+referr for z in maxabs];E=[F(g['error']) for g in groups]
        # M = sup |F-R|. At fixed ideal y, the rounded scalar y' differs
        # by <=epsy. Markov on the y approximant bounds the sample discrepancy
        # by 2 Ey + gamma(M+Ey). Then apply c norming in each group, followed
        # by exact-node y norming; move lambda*M to the left.
        T=max(C*(d+2*Ey+gamma*Ey+e)+e for C,d,e in zip(Cc,delta,E))
        residual=(Cy*(T+Ey)+Ey)/(1-lam)
        lower=ref_lower[i]-residual
        Hphys=lower/(1+b*(b+1)*Hl)**4-F(data['physical_tail'])-F(1,2**140)-F(1,10**14)
        result.append(dict(b=b,node_residuals=list(map(str,delta)),ref_lower=str(ref_lower[i]),uniform_residual=str(residual),prefix_lower=str(lower),physical_lower=str(Hphys),lower_float=float(Hphys),status='PASS_CONDITIONAL' if Hphys>0 else 'INCONCLUSIVE'))
    print(stem,[(z['b'],z['lower_float'],z['status']) for z in result],'seconds',time.monotonic()-st,flush=True)
    json.dump(dict(cell=stem,count=count,precision=P,Ny=my,Cc=list(map(str,Cc)),Cy=str(Cy),epsy=str(epsy),gamma=str(gamma),lam=str(lam),scalar_error=str(err),reference_error=str(referr),reference=ref,inputs=inputs,results=result),open(wd/f'{stem}_group_ledger.json','w'),indent=2)
    assert all(z['status']=='PASS_CONDITIONAL' for z in result)
if __name__=='__main__':
    for p in sys.argv[1:]:make_cell_group(p)

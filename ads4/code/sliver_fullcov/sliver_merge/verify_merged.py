"""Check executed grids, arithmetic, analytic bounds and interval union.
Does not redo the archived scalar runs. Use the drivers to rerun them.
"""
import sys,hashlib,json,math
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,ceil,comb
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'sliver_close'))
from coupled_plan import check_tail, maxrad
from full_preflight import residue_degree
from reference_driver import cos_nodes,norm_constant,NODE_DEN,DEN,UNIT
from group_plan import group_check
from short_plan import tail_vector
sys.set_int_max_str_digits(0)
assert hashlib.sha256((HERE/'coeff_compact60.txt').read_bytes()).hexdigest()=='811e15cdd95319309edf64f5cf862bb674f43fbc73f551a002577c4b56b1ee06'
assert maxrad<=2

def exact_schedule(K,corners):
    upper=2**100;shifts=[0];sh=0
    for j in range(K):
        mx=F(0)
        for c,h in corners:
            for s in (1,9):
                d=1-h*j*(j+s);assert d>0;mx=max(mx,abs(c-h*j*(j+s))/d)
        assert mx<=1;upper=ceil(upper*mx);t=0
        while upper<=2**99 and t<2:upper*=2;sh+=1;t+=1
        assert 2**99<upper<=2**100;shifts.append(sh)
    return shifts

def check_analytic(z,group):
    l,u,Hl,Hu=map(F,(z['l'],z['u'],z['Hl'],z['Hu']));K=z['K'];assert 0<l<u<=F('.9') and Hl>=Hu>0
    for c,h in [(l,Hl),(u,Hu)]:
        r=F(isqrt(c.numerator*10**44//c.denominator),10**22);assert r*r<=c and h>=(1-r)**2/998**2
    for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)]:
        for s in (1,9):
            for j in (0,K-1):assert 1-h*j*(j+s)>0 and abs(c-h*j*(j+s))<=1-h*j*(j+s)
    tail=max(F(check_tail(l,u,Hl,K)),F(z['physical_tail']))
    t=z['y'];assert t['method']=='residue_r4_group16';ey=max(F(residue_degree(l,u,Hl,K,t['deg'])[1]),F(t['error']));ecnew=[];shortnew=[]
    for g in group:
        Kg=g.get('K',K);m=g['deg'];rs=list(map(F,g['rhos']));ch=g['choices'];assert 1<=Kg<=K and len(ch)==Kg+1 and all(r>1 for r in rs) and all(0<=j<len(rs) for j in ch)
        ec=F(group_check(l,u,Hl,Hu,Kg,F(g['lo']),F(g['hi']),m,rs,ch))
        if 'short_tail' in g:
            st=max(F(g['short_tail']),F(tail_vector(l,u,Hl,Hu,K,F(g['lo']),F(g['hi']))[Kg]));assert st<F('1e-10') and F(g['error'])==F(g['short_tail'])+F(g['poly_error']);ec=max(ec,F(g['poly_error']))+st
        else:st=F(0)
        shortnew.append(st);ecnew.append(max(ec,F(g['error'])))
    assert tail<F('1e-11') and ey<F('1e-10') and all(e<F('2e-10') for e in ecnew)
    return l,u,Hl,Hu,K,tail,ey,ecnew,shortnew

def ref_bounds(ref):
    assert len(ref)==2
    L=[]
    for P in ref:
        assert len(P)==3 and all(len(row)==3 for row in P) and all(isinstance(x,int) for row in P for x in row)
        assert sum(abs(x) for row in P for x in row)<=UNIT//1000
        B=[[sum(F(P[a][d]*comb(i,a)*comb(j,d),comb(2,a)*comb(2,d)*UNIT) for a in range(i+1) for d in range(j+1)) for j in range(3)] for i in range(3)]
        L.append(min(min(row) for row in B));assert L[-1]>0
    return L

def verify_new(ledpath):
    ledger=json.load(open(ledpath));stem=ledger['cell'];z=json.load(open(HERE/(stem+'.json')));grouped=ledpath.name.endswith('_group_ledger.json');short=ledpath.name.endswith('_short_group_ledger.json');my=z['y']['deg']
    if grouped:groups=json.load(open(HERE/(stem+('_short_group.json' if short else '_group.json'))))['groups']
    else:groups=[dict(z['c'],lo='0',hi='1',nodes=list(range(my+1)))]
    l,u,Hl,Hu,K,tail_checked,ey_checked,ec_checked,short_checked=check_analytic(z,groups)
    assert sorted(j for g in groups for j in g['nodes'])==list(range(my+1));assert ledger['Ny']==my
    P=ledger.get('precision',180);assert 64<=P<=180
    err=(K+1)*(F(4)+Hl*((K+1)**2+K*K+8*K+67))/10**60+F(10*(K+1),2**P);assert err==F(ledger['scalar_error'])
    epsc=F(2,10**24)+F(1,DEN)/(u-l);epsy=F(2,10**24)+F(1,DEN)/Hu
    Cy=norm_constant(my,F(0) if grouped else epsy);assert Cy==F(ledger['Cy']);ref=ledger['reference'];ref_lower=ref_bounds(ref);referr=F(4,1000*2**50)+F(8,UNIT);assert referr==F(ledger['reference_error'])
    Cc=[];maxabs=[];count=0;yy=cos_nodes(my);Ey=F(z['y']['error']);Ec=[]
    for ig,g in enumerate(groups):
        m=g['deg'];Kg=g.get('K',K);lo,hi=F(g['lo']),F(g['hi']);nodes=g['nodes'];yg=[yy[j] for j in nodes];assert yg and 0<=lo<=hi<=1
        for y in yg:assert lo<=max(F(0),F(y,NODE_DEN)-F(2,10**24)) and hi>=min(F(1),F(y,NODE_DEN)+F(2,10**24))
        Cc.append(norm_constant(m,epsc));Ec.append(F(g['error']))
        if grouped:ylo,yhi=max(F(0),lo-epsy),min(F(1),hi+epsy)
        else:ylo,yhi=F(0),F(1)
        shifts=exact_schedule(Kg,[(c,H*y) for c,H in [(l,Hl),(u,Hu)] for y in (ylo,yhi)])
        expected=[]
        for a in cos_nodes(m):
            c0=l+(u-l)*F(a,NODE_DEN);C=max(ceil(l*DEN),min(int(u*DEN),int(c0*DEN)));c=F(C,DEN);assert l<=c<=u and abs(c-c0)<=F(1,DEN)
            H=Hl+(Hu-Hl)*(c-l)/(u-l);expected.append((C,H.numerator,H.denominator))
        assert len(set(C for C,_,_ in expected))==len(expected)
        mx={0:0,4:0}
        for w in range(4):
            suffix=f'_s{ig}_{w}' if short else (f'_g{ig}_{w}' if grouped else f'_r{w}');inp=HERE/(stem+suffix+'.in');out=HERE/(stem+suffix+'.out');raw=inp.read_bytes()
            if grouped:assert hashlib.sha256(raw).hexdigest()==ledger['inputs'][inp.name]
            v=list(map(int,raw.split()));rows=expected[w::4];exp=[Kg,DEN,len(yg),len(rows),P]+shifts+[NODE_DEN]+yg+[l.numerator,l.denominator,u.numerator,u.denominator]+[x for P0 in ref for row in P0 for x in row]+[a for row in rows for a in row];assert exp==v,inp.name
            lines=out.read_text().splitlines();assert len(lines)==3;seen=set();n=0
            for line in lines:
                k,*a=line.split()
                if k=='COUNT':assert not n and len(a)==1;n=int(a[0]);assert n==2*len(rows)*len(yg);count+=n
                else:
                    b=int(k);assert b in (0,4) and b not in seen;seen.add(b);a,bv=map(int,a);assert a<=bv;mx[b]=max(mx[b],abs(a),abs(bv))
            assert seen=={0,4} and n
        maxabs.append(mx)
    assert count==ledger['count']==2*sum((g['deg']+1)*len(g['nodes']) for g in groups)
    if grouped:
        assert Cc==list(map(F,ledger['Cc']));gamma=2*my**2*epsy;lam=Cy*max(Cc)*gamma;assert gamma==F(ledger['gamma']) and lam==F(ledger['lam']) and lam<F(1,10)
    else:assert len(Cc)==1 and Cc[0]==F(ledger['Cc']) and ledger['Nc']==groups[0]['deg']
    rr=[]
    for i,b in enumerate((0,4)):
        delta=[F(z[b],UNIT)+err+referr+F(g.get('short_tail','0')) for z,g in zip(maxabs,groups)]
        if grouped:
            T=max(C*(d+2*Ey+gamma*Ey+e)+e for C,d,e in zip(Cc,delta,Ec));residual=(Cy*(T+Ey)+Ey)/(1-lam)
        else:residual=Cc[0]*(Cy*delta[0]+(Cy+1)*Ey)+(Cc[0]+1)*Ec[0]
        lower=ref_lower[i]-residual;Hphys=lower/(1+b*(b+1)*Hl)**4-F(z['physical_tail'])-F(1,2**140)-F(1,10**14);r=ledger['results'][i]
        assert r['b']==b and F(r['ref_lower'])==ref_lower[i] and F(r['uniform_residual'])==residual and F(r['prefix_lower'])==lower and F(r['physical_lower'])==Hphys and Hphys>0
        assert (list(map(F,r['node_residuals']))==delta if grouped else F(r['node_residual'])==delta[0])
        if grouped:
            dnew=[d+st-F(g.get('short_tail','0')) for d,st,g in zip(delta,short_checked,groups)];Tnew=max(C*(d+2*ey_checked+gamma*ey_checked+e)+e for C,d,e in zip(Cc,dnew,ec_checked));resnew=(Cy*(Tnew+ey_checked)+ey_checked)/(1-lam)
        else:resnew=Cc[0]*(Cy*delta[0]+(Cy+1)*ey_checked)+(Cc[0]+1)*ec_checked[0]
        Hcert=(ref_lower[i]-resnew)/(1+b*(b+1)*Hl)**4-tail_checked-F(1,2**140)-F(1,10**14)
        assert 0<Hcert<=Hphys
        rr.append(dict(b=b,c_lo=str(l),c_hi=str(u),lower_exact=str(Hcert),lower=float(Hcert),source=ledpath.name,tail_checked=str(tail_checked),ey_checked=str(ey_checked),ec_checked=list(map(str,ec_checked)),node_tail_checked=list(map(str,short_checked))))
    print(stem,'verified',count,'nodes',[(r['b'],r['lower']) for r in rr],flush=True)
    return rr

def union(rows):
    ret={}
    for b in (0,4):
        it=sorted([(F(r['c_lo']),F(r['c_hi'])) for r in rows if r['b']==b and F(r['lower_exact'])>0]);merged=[]
        for l,u in it:
            assert l<u
            if merged and l<=merged[-1][1]:merged[-1]=(merged[-1][0],max(u,merged[-1][1]))
            else:merged.append((l,u))
        ret[str(b)]=[[str(l),str(u)] for l,u in merged]
    return ret

if __name__=='__main__':
    rows=[]
    files=[Path(s) for s in sys.argv[1:]] if len(sys.argv)>1 else sorted(HERE.glob('*_ref_ledger.json'))+sorted(HERE.glob('*_group_ledger.json'))
    for p in files:rows+=verify_new(p)
    print('new coverage:',union(rows),flush=True)
    json.dump(dict(status='VERIFIED_NEW_LEDGERS',new_rows=rows,union=union(rows)),open(HERE/'verified_new.json','w'),indent=2)

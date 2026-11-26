from fractions import Fraction as F
from pathlib import Path
from math import ceil
import json,hashlib
from norm_driver_scaled import ngrid,DEN,cos_nodes,NODE_DEN
root=Path(__file__).resolve().parent
for stem in ('cell_3500_4462','cell_4462_5305'):
    x=json.load(open(root/(stem+'.json')));d=json.load(open(root/(stem+'_ledger.json')))
    l,u,Hl,Hu=map(F,(x['l'],x['u'],x['Hl'],x['Hu']));K=x['K']
    Nc=ngrid(x['c']['deg'],F(4,10**24)+F(2,DEN)/(u-l));Ny=ngrid(x['y']['deg'],F(4,10**24)+F(2,DEN)/Hu)
    assert (Nc,Ny)==(d['Nc'],d['Ny'])
    yc,yy=cos_nodes(Nc),cos_nodes(Ny)
    sc=2**100;upper=sc;shifts=[0];sh=0
    for j in range(K):
        ps=[(c,h,j*(j+s)) for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)] for s in (1,9)]
        assert all(1-h*p>0 for c,h,p in ps)
        mx=max(abs(c-h*p)/(1-h*p) for c,h,p in ps);assert mx<=1
        upper=ceil(upper*mx);t=0
        while upper<=sc//2 and t<2:upper*=2;sh+=1;t+=1
        assert sc//2<upper<=sc;shifts.append(sh)
    rows=[]
    for a in yc:
        c0=l+(u-l)*F(a,NODE_DEN);C=max(ceil(l*DEN),min((u*DEN).numerator//(u*DEN).denominator,(c0*DEN).numerator//(c0*DEN).denominator));c=F(C,DEN);H=Hl+(Hu-Hl)*(c-l)/(u-l);rows.append((C,H.numerator,H.denominator))
    lo={0:10**100,4:10**100};hi={0:-10**100,4:-10**100};total=0
    for i in range(4):
        r=rows[i::4]
        inp=(root/f'{stem}_{i}.in').read_text().splitlines()
        assert list(map(int,inp[0].split()))==[K,DEN,len(yy),len(r),180]
        assert list(map(int,inp[1].split()))==shifts
        assert int(inp[2])==NODE_DEN and list(map(int,inp[3].split()))==yy
        assert [tuple(map(int,a.split())) for a in inp[4:]]==r
        lines=(root/f'{stem}_{i}.out').read_text().splitlines();assert len(lines)==3
        assert lines[-1].split()==['COUNT',str(2*len(r)*len(yy))];total+=int(lines[-1].split()[1])
        for line in lines[:2]:
            b,a,z=map(int,line.split());lo[b]=min(lo[b],a);hi[b]=max(hi[b],z)
    assert total==d['count']==2*(Nc+1)*(Ny+1)
    err=(K+1)*(F(4)+Hl*((K+1)**2+K*K+8*K+67))/10**60+F(10*(K+1),2**180)
    assert err==F(d['scalar_error'])<F(1,10**45)
    for r in d['results']:
        b=r['b'];a,z=F(lo[b],10**20)-err,F(hi[b],10**20)+err
        assert a==F(r['min_node']) and z==F(r['max_node'])
        R,delta=(a+z)/2,(z-a)/2;C=F(5,2)
        assert R==F(r['reference']) and delta==F(r['delta'])
        lower=R-C*C*delta-C*(C+1)*F(x['y']['error'])-(C+1)*F(x['c']['error'])
        phys=lower/(1+b*(b+1)*Hl)**4-F(x['physical_tail'])-F(1,2**140)-F(1,10**14)
        assert lower==F(r['prefix_lower']) and phys==F(r['physical_lower'])>0
        assert r['status']=='PASS_CONDITIONAL'
    print(stem,'PASS input, node coverage, shift schedule, output count and rational sign inequalities')
print('PARTIAL: two executed cells only. No full-sliver closure.')

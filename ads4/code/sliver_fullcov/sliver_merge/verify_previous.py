"""Exact reassembly of previous local norming certificates and interval union.
The previous direct triangle ledgers are rechecked by their original scripts.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt,ceil
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'sliver_close'))
from coupled_plan import check_tail,exact_check
from full_preflight import residue_degree
# Explicitly load the older node helper, which used a different denominator.
import importlib.util
spec=importlib.util.spec_from_file_location('oldnorm',HERE.parent/'sliver_close/norm_driver.py');oldnorm=importlib.util.module_from_spec(spec);spec.loader.exec_module(oldnorm)
sys.path.insert(0,str(HERE))
from norm_driver_scaled import ngrid,cos_nodes,NODE_DEN,DEN
from verify_merged import exact_schedule,union

def oldcell(root,stem,scaled):
    root=Path(root);data=json.load(open(root/(stem+'.json')));led=json.load(open(root/(stem+'_ledger.json')));K=data['K'];l,u,Hl,Hu=map(F,(data['l'],data['u'],data['Hl'],data['Hu']))
    for c,h in [(l,Hl),(u,Hu)]:
        t=F(isqrt(c.numerator*10**44//c.denominator),10**22);assert t*t<=c and h>=(1-t)**2/998**2
    for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)]:
        for j in (0,K-1):
            for s in (1,9):assert 1-h*j*(j+s)>0 and abs(c-h*j*(j+s))<=1-h*j*(j+s)
    tail=max(F(check_tail(l,u,Hl,K)),F(data['physical_tail']));eps=[]
    for axis in ['c','y']:
        z=data[axis];m=z['deg']
        if z.get('method')=='residue_r4_group16':e=F(residue_degree(l,u,Hl,K,m)[1])
        else:e=F(exact_check(l,u,Hl,Hu,K,axis,m,list(map(F,z['rhos'])),z['choices']))
        eps.append(max(e,F(z['error'])))
    assert tail<F('1e-11') and all(e<F('1e-10') for e in eps)
    den=DEN if scaled else oldnorm.DEN
    Nc=ngrid(data['c']['deg'],F(4,10**24)+F(2,den)/(u-l));Ny=ngrid(data['y']['deg'],F(4,10**24)+F(2,den)/Hu)
    assert (Nc,Ny)==(led['Nc'],led['Ny']);yc,yy=cos_nodes(Nc),cos_nodes(Ny)
    rows=[]
    for a in yc:
        c0=l+(u-l)*F(a,NODE_DEN);C=max(ceil(l*den),min(int(u*den),int(c0*den)));c=F(C,den);H=Hl+(Hu-Hl)*(c-l)/(u-l);rows.append((C,H.numerator,H.denominator))
    shifts=exact_schedule(K,[(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)]) if scaled else []
    lo={0:10**100,4:10**100};hi={0:-10**100,4:-10**100};total=0
    for w in range(4):
        row=rows[w::4];raw=list(map(int,(root/f'{stem}_{w}.in').read_text().split()));expected=[K,den,len(yy),len(row)]+([180]+shifts if scaled else [])+[NODE_DEN]+yy+[x for r in row for x in r];assert raw==expected
        lines=(root/f'{stem}_{w}.out').read_text().splitlines();assert len(lines)==3 and lines[-1].split()==['COUNT',str(2*len(row)*len(yy))];total+=int(lines[-1].split()[1]);seen=set()
        for ln in lines[:2]:
            b,a,z=map(int,ln.split());assert b in (0,4) and b not in seen and a<=z;seen.add(b);lo[b]=min(lo[b],a);hi[b]=max(hi[b],z)
    assert total==led['count']==2*(Nc+1)*(Ny+1)
    err=(K+1)*(F(4)+Hl*((K+1)**2+K*K+8*K+67))/10**60+(F(10*(K+1),2**180) if scaled else F(0));assert err==F(led['scalar_error'])
    out=[];C=F(5,2)
    for r in led['results']:
        b=r['b'];a,z=F(lo[b],10**20)-err,F(hi[b],10**20)+err;R,delta=(a+z)/2,(z-a)/2;assert a==F(r['min_node']) and z==F(r['max_node']) and R==F(r['reference']) and delta==F(r['delta'])
        L=R-C*C*delta-C*(C+1)*F(data['y']['error'])-(C+1)*F(data['c']['error']);H=L/(1+b*(b+1)*Hl)**4-F(data['physical_tail'])-F(1,2**140)-F(1,10**14);assert L==F(r['prefix_lower']) and H==F(r['physical_lower'])>0
        L=R-C*C*delta-C*(C+1)*eps[1]-(C+1)*eps[0];Hnew=L/(1+b*(b+1)*Hl)**4-tail-F(1,2**140)-F(1,10**14);assert 0<Hnew<=H
        out.append(dict(b=b,c_lo=str(l),c_hi=str(u),lower_exact=str(Hnew),lower=float(Hnew),source=stem+'_ledger.json',tail_checked=str(tail),ec_checked=str(eps[0]),ey_checked=str(eps[1])))
    assert sorted(r['b'] for r in out)==[0,4]
    print(stem,'verified',total,'nodes',[(r['b'],r['lower']) for r in out],flush=True);return out

if __name__=='__main__':
    # These two ledgers are produced/rechecked by the original triangle scripts.
    rows=[]
    for directory in ['sliver_next_work','physical_sliver_run']:
        d=json.load(open(HERE.parent/directory/'new_sign_ledger.json'))
        for r in d['rows']:
            # Charge the more conservative common 1e-14 coherent allowance as well.
            H=F(r['lower_exact'])-F(1,10**14);assert H>0
            rows.append(dict(b=r['b'],c_lo=r['c_lo'],c_hi=r['c_hi'],lower_exact=str(H),lower=float(H),source=directory+'/new_sign_ledger.json'))
    rows+=oldcell(HERE.parent/'sliver_close','cell_30_35',False)
    rows+=oldcell(HERE.parent/'sliver_fast','cell_3500_4462',True)
    rows+=oldcell(HERE.parent/'sliver_fast','cell_4462_5305',True)
    assert union(rows)=={'0':[['129/500','1061/2000']],'4':[['129/500','1061/2000']]}
    print('previous exact union',union(rows),flush=True)
    json.dump(dict(rows=rows,union=union(rows)),open(HERE/'verified_previous.json','w'),indent=2)

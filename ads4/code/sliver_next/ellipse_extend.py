#!/usr/bin/env python3
"""Extend a certified real-frequency bound by Hadamard on the Joukowski annulus.
The outer bound is on the *complete* J=1500 transfer, not only sampled K.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys,time
from contour_direct import *
from functools import lru_cache

def eta_for(tlo,thi,yhi,center,half):
    # ellipse containing target rectangle; rational exact ellipse membership
    X=max(abs(tlo-center),abs(thi-center))/half;Z=yhi/half
    if Z==0 and X<=1:return F(1)
    r=F('1.001')
    while True:
        A=(r+1/r)/2;B=(r-1/r)/2
        if X*X/(A*A)+Z*Z/(B*B)<=1:return r
        r+=F('.005')
        if r>10:raise ValueError('ellipse extension too far')

def choose_row(rows,p,b,tlo,thi):
    for r in rows:
        if r['p']!=p or r['b']!=b or r['status'] not in ('PASS','ENCLOSED'):continue
        c=r['config'];d,dd,h,dh,v,dv,im=parameters(c)
        if v-dv-F("1e-9")<=tlo<=thi<=v+dv+F("1e-9") and im==0 and c.get('nuAxis','real')=='real':return r
    return None

class Extend:
    def __init__(self,rows,dl,du):self.rows=rows;self.dl=dl;self.du=du;self.cache={}
    def outer(self,row,rho,p):
        key=(row['name'],row['b'],p,rho)
        if key not in self.cache:
            c=dict(row['config']);c['K']=1500
            from major_rect import major as rect
            vals=[]
            try:vals.append(major(c,p,row['b'],[F(1),F(1),(rho+1/rho)/2])[0])
            except ValueError:pass
            try:vals.append(rect(c,p,row['b'],2,rho)[0])
            except ValueError:pass
            self.cache[key]=F(min(vals),NS) if vals else None
        return self.cache[key]
    def get(self,p,b,tlo,thi,ylo,yhi):
        # Real transfer certificates available only p3; p4 use direct absolute majorant.
        row=choose_row(self.rows,p,b,tlo,thi)
        if row is None:return abs_transfer(self.dl,self.du,tlo,thi,ylo,yhi,p,b)
        c=row['config'];d,dd,h,dh,v,dv,im=parameters(c)
        assert c['fullTwist'] and d-dd<=self.dl<=self.du<=d+dd and im==0
        cap_check(c)
        uP,uQ=F(row['P_bound']),F(row['Q_bound'])
        eta=eta_for(tlo,thi,yhi,v,dv)
        if eta==1:return uP,uQ
        candidates=[]
        for rho in map(F,[2,3,4,6,8,12,16,24,32]):
            if rho<=eta:continue
            M=self.outer(row,rho,p)
            if M is None:continue
            w=min(F(1),hi(logu(eta)/logl(rho)))
            def ext(U,V):
                if V<=U:return V
                return hi(exp_u((1-w)*logu(U)+w*logu(V)))
            candidates.append((ext(uP,M),ext(uQ,M/(d-dd))))
        if not candidates:return abs_transfer(self.dl,self.du,tlo,thi,ylo,yhi,p,b)
        return min(z[0] for z in candidates),min(z[1] for z in candidates)

def right(dl,du,p,b,Y,E):
    y=F(0);total=F(0);worst=(F(0),None)
    while y<Y:
        u=min(Y,y+F('.25'));P,Q=E.get(p,b,RIGHT,RIGHT,y,u);v=hi((u-y)*osc_integrand(dl,du,p,b,RIGHT,RIGHT,y,u,P,Q));total+=v
        if v>worst[0]:worst=(v,(y,u))
        y=u
    return hi(total),worst

def top(dl,du,p,b,Y,E):
    # Cell boundaries at existing real certificate boundaries and small-frequency breaks.
    t=DELTA;total=F(0);worst=(F(0),None)
    edges=[F('1'),F('10'),F('159'),F('400'),F('700'),F('850'),F('950'),RIGHT]
    while t<RIGHT:
        step=F('.25') if t<1 else F(1) if t<10 else F(25)
        u=min(t+step,min(x for x in edges if x>t),RIGHT);P,Q=E.get(p,b,t,u,Y,Y);val=hi((u-t)*osc_integrand(dl,du,p,b,t,u,Y,Y,P,Q));total+=val
        if val>worst[0]:worst=(val,(t,u))
        t=u
    return hi(total),worst

if __name__=='__main__':
    path=sys.argv[1];rows=json.load(open(path));c=rows[0]['config'];d,dd,*_=parameters(c);dl,du=d-dd,d+dd
    Y=F(sys.argv[2] if len(sys.argv)>2 else '10.25');E=Extend(rows,dl,du);out=[]
    for p in ((int(sys.argv[4]),) if len(sys.argv)>5 else (3,4)):
      for b in ((int(sys.argv[5]),) if len(sys.argv)>5 else (0,4)):
        t=time.perf_counter();y=Y+p-3
        v,w=near_contours(dl,du,p,b,y)
        u,ut=top(dl,du,p,b,y,E);r,rt=right(dl,du,p,b,y,E);x=extra(dl,du,p,b)
        row=dict(p=p,b=b,Y=str(y),dlo=str(dl),dhi=str(du),bottom=str(v),left=str(w),top=str(u),right=str(r),extra=str(x),total=str(v+w+u+r+x),seconds=time.perf_counter()-t)
        print(json.dumps({**row,**{k+'_approx':float(F(row[k]))for k in ('bottom','left','top','right','extra','total')}}),flush=True);out.append(row)
    if len(sys.argv)>3:Path(sys.argv[3]).write_text(json.dumps(out,indent=2)+'\n')

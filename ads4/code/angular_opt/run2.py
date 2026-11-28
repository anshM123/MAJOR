import numpy as np
from lp_dual import *
def states2(g, nt=40):
    Js = sorted(set(list(range(0, 8*g)) + [int(x) for x in np.geomspace(8*g, 2000*g, 200)]))
    ts = np.concatenate([[0.0], np.geomspace(0.003*g, 300*g, nt)])
    out=[]
    for J in Js:
        for b in (0,4):
            if b==0 and J%2: continue
            if b==4 and J<4: continue
            for t in ts: out.append((b,J,J+g+t))
    return out
for g in [10,20]:
  for n in [6,8,10,12,16,20,24]:
    S = states2(g)
    R = build(S,n); sc=np.abs(R).max(axis=0); sc[sc==0]=1
    c=np.zeros(n+1); c[2]=1/sc[2]
    res=linprog(c,A_ub=-R/sc,b_ub=np.ones(len(S)),bounds=[(None,None)]*(n+1),method='highs')
    if res.status!=0: print(g,n,'status',res.status); continue
    lam=res.x/sc; L=-lam[2]
    print(g,n,'L=%.5g'%L,'tau*g^3=%.1f'%(4*np.sqrt(3/L)*g**3), 'tau*g=%.3f'%(4*np.sqrt(3/L)*g), 'check', check(g,lam,3000*g//10,nt=100)[0],flush=True)

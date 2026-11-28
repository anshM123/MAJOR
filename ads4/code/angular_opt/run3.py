import numpy as np
from lp_dual import *
from run2 import states2
g=20; n=12
S = states2(g)
R = build(S,n); sc=np.abs(R).max(axis=0); sc[sc==0]=1
c=np.zeros(n+1); c[2]=1/sc[2]
res=linprog(c,A_ub=-R/sc,b_ub=np.ones(len(S)),bounds=[(None,None)]*(n+1),method='highs')
lam=res.x/sc
print('lam*D^2/g^(2k+2):', [lam[k]*((k+1)*(k+2)*(k+3)/6)**2/g**(2*k+2) for k in range(n+1)])
for J in list(range(0,40))+[50,60,80,100,150,200,400]:
    for b in (0,4):
        if b==0 and J%2: continue
        if b==4 and J<4: continue
        vals=[]
        for t in [0,1,3,10,30,100]:
            D=J+g+t; r=ladder(b,J,D,n)
            vals.append(np.dot(lam,r)/W(b,D))
        print(b,J,' '.join('%+.3f'%v for v in vals))

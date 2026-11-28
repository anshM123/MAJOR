import numpy as np, sys
from lp_light2 import solve
from math import comb
g=1000; n=int(sys.argv[1]); kappa=float(sys.argv[2])
L, lam = solve(g, n, kappa)
D2 = lambda k: ((k+1)*(k+2)*(k+3)/6)**2
ell = np.array([lam[k]*D2(k) for k in range(n+1)])
print('L', L, 'ell', ell)
print('sum|ell|', np.abs(ell).sum(), ' sum|ell|/kappa', np.abs(ell).sum()/kappa)
M=g-2
z=np.arange(M+1)
vals=np.zeros(M+1)
for k in range(n+1):
    r=np.ones(M+1)
    for i in range(k): r*= (z-i)/(M-i)
    vals+=ell[k]*r
WB=((g+3)*(g+2)-12)/3+1.5; W0=(g+1)*g/3+1.5; W0g=(g-1)*(g-2)/3+1.5
print('heavy grid: min ell(z)', vals.min(), ' (need >= -WB=', -WB, ')  max ell(z)', vals.max(), '(need <= W0(g+2)=',W0,')', 'ell0',ell[0],'<=',W0g)
np.save('cache/robust_lam_n%d_k%.0e.npy'%(n,kappa), lam)

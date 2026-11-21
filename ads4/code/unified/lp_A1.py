"""Variant of ads4/code/kernels/lp_balance.py with Result D's A1: every operator (all J, incl. J=0,2) has twist Delta-J >= gap.
Usage: python3 lp_A1.py gap Jmax"""
import sys
sys.path.insert(0, '/home/user/MAJOR/ads4/code/kernels')
import numpy as np, mpmath as mp
from scipy.optimize import linprog
import closed_forms as cf
mp.mp.dps = 30
def run(gap, Jmax, ntheta=721):
    offs = [0, 0.02, 0.1, 0.3, 0.6, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 256, 512, 2048]
    kmax = Jmax; pp=[]; pm=[]
    def add(lst, W, J):
        base = gap + J
        for o in offs:
            D = mp.mpf(base + o)
            v = np.array([float(W(k, D, J)) for k in range(kmax+1)])
            if not np.all(np.isfinite(v)) or v.max() <= 0: continue
            sc = v.max(); lst.append((J, float(D), v/sc, sc))
    for J in range(0, Jmax+1, 2): add(pp, cf.W0, J)
    for J in range(4, Jmax+5): add(pm, cf.W4, J)
    npp, npm = len(pp), len(pm); nv = 1+npp+npm
    Aeq = np.zeros((kmax+1, nv)); beq = np.zeros(kmax+1)
    for k in range(kmax+1):
        if k == 2: Aeq[k,0] = 1
        for i,(J,D,v,sc) in enumerate(pp): Aeq[k,1+i] = v[k]
        for j,(J,D,v,sc) in enumerate(pm): Aeq[k,1+npp+j] = -v[k]
    th = np.linspace(0, np.pi, ntheta); Aub = np.zeros((ntheta, nv)); bub = np.ones(ntheta)
    for j,(J,D,v,sc) in enumerate(pm): Aub[:,1+npp+j] = -2*sum(v[k]*np.cos((2+k)*th) for k in range(kmax+1))
    cost = np.zeros(nv); cost[0] = -1
    res = linprog(cost, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(0,None)]*nv, method='highs')
    c2 = -res.fun if res.status == 0 else None
    print(f"A1 gap={gap} Jmax={Jmax}: status={res.status} max|c|^2={c2}  |c|^2 gap^4={c2*gap**4 if c2 else None}  sqrt(t4^2+aT^2)=4|c|={4*np.sqrt(c2) if c2 else None}", flush=True)
run(float(sys.argv[1]), int(sys.argv[2]))

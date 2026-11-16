"""
Exact four-operator solutions of all sum rules (S_k) with c != 0 and an arbitrary gap (KERNELS_D3.md section 7.3):
   '+-' operators:  O_4 (J=4, Delta_4 >= gap), weight q4;   O_6 (J=6, Delta_6 = D >= gap), weight q6
   '++' operators:  scalar (any Delta), weight p0;          O_2 (J=2, Delta' = Dp > 3), weight p2
   k=0: p0 + p2 = q4 + q6;  k=1: p2 W0_1(Dp,2) = q6 W4_1(D,6);  k=2: |c|^2 + p2 W0_2(Dp,2) = q6 W4_2(D,6);  k>=3: 0 = 0.
   EEC(theta) = 1 + 2 a2 cos 2theta + 2 a3 cos 3theta + 2 a4 cos 4theta >= 0, a2 = q4+q6, a3 = q6 W4_1, a4 = q6 W4_2.
Maximises |c|^2 over (q4, q6, Dp) numerically for given gap (D = Delta_4 = gap).
"""
import sys
import numpy as np
from scipy.optimize import minimize

W41 = lambda D: 352.0/((D - 6)*(D + 3))
W42 = lambda D: 26400.0/((D - 7)*(D - 6)*(D + 3)*(D + 4))
W01 = lambda Dp: 96.0/((Dp - 2)*(Dp - 1))
W02 = lambda Dp: 2400.0/((Dp - 3)*(Dp - 2)*(Dp - 1)*Dp)
TH = np.linspace(0, np.pi, 2001)

def data(gap, q4, q6, Dp):
    a2, a3, a4 = q4 + q6, q6*W41(gap), q6*W42(gap)
    p2 = a3/W01(Dp); p0 = a2 - p2
    c2 = a4 - p2*W02(Dp)
    eec = 1 + 2*a2*np.cos(2*TH) + 2*a3*np.cos(3*TH) + 2*a4*np.cos(4*TH)
    return dict(c2=c2, p0=p0, p2=p2, a=(a2, a3, a4), eecmin=eec.min())

def best(gap):
    def obj(v):
        q4, q6, Dp = np.abs(v[0]), np.abs(v[1]), 3 + np.abs(v[2])
        d = data(gap, q4, q6, Dp)
        pen = 1e3*(max(0, -d['p0'])**2 + max(0, -d['eecmin'])**2)
        return -d['c2'] + pen
    bestv = None
    for q4 in (0.1, 0.3, 0.5):
        for q6 in (0.01, 0.05, 0.2):
            for Dp in (5.0, 10.0, 0.5*gap, gap):
                r = minimize(obj, [q4, q6, Dp - 3], method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-14, maxiter=20000))
                if bestv is None or r.fun < bestv.fun:
                    bestv = r
    q4, q6, Dp = abs(bestv.x[0]), abs(bestv.x[1]), 3 + abs(bestv.x[2])
    # shrink slightly to be strictly feasible
    d = data(gap, q4, q6, Dp)
    return q4, q6, Dp, d

if __name__ == '__main__':
    print(" gap     q4       q6      Delta'    p0       p2      (a2,a3,a4)                |c|^2     |c|^2 gap^4   min EEC")
    for gap in [float(a) for a in sys.argv[1:]] or [8, 10, 15, 20, 30, 50, 100, 300, 1000]:
        q4, q6, Dp, d = best(gap)
        print(f"{gap:6.0f} {q4:8.4f} {q6:8.5f} {Dp:8.3f} {d['p0']:8.4f} {d['p2']:8.4f}  ({d['a'][0]:.4f},{d['a'][1]:.4f},{d['a'][2]:.5f})"
              f"  {d['c2']:.6g}  {d['c2']*gap**4:10.1f}   {d['eecmin']:.1e}")
    print("explicit strictly feasible example (gap = 20): q4 = 0.29, q6 = 0.16, Delta_4 = Delta_6 = 20, Delta' = 17")
    d = data(20.0, 0.29, 0.16, 17.0)
    print("   ", {k: (tuple(round(x, 6) for x in v) if isinstance(v, tuple) else round(float(v), 6)) for k, v in d.items()})
    print("asymptotics: |c|^2 -> 5400 a2/gap^4 with Delta' = (2/3)^(1/2) gap, p0 = 0, a2 -> 1/2")

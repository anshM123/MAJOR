"""
Checks of the lightcone + asymptotic-tail representation of single d = 3 scalar blocks (lcblock.LCBlock) in the corner rho_b -> 1:
 (1) leading lightcone function H_0 against the exact collinear block: H_0(x) = 4^{-(D+J)/2} x^{-(D-J)/2} k_{D+J}(z(x)),
     k_b(z) = z^{b/2} 2F1(b/2 + c1, b/2 + c2; b; z), z = 4x/(1+x)^2, evaluated with mpmath at 1 - x down to 1e-8 (independent);
 (2) derivative sums S[k,l] against the radial series with N = 24000 levels at points where the latter has converged
     (X = -N log rho_b > 90): this tests the extrapolation q = Q ... 24000 of the tail;
 (3) Casimir residual (casimir_check.py, relative to the largest term) in the deep corner, and the difference between the two
     LC settings A, B of corner.py.
Output: out_c1_block_checks.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
from lcblock import LCBlock
from blocks import CBlock
from casimir_check import casimir_residual
from corner import PAR_A, PAR_B
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def corner_points(Delta, J, ts, us):
    al, alb = Delta - J - 2, Delta + J - 1
    pts = []
    for t in ts:
        for u in us:
            s = t/(2*(al + alb*u)); s2 = s*u
            if s >= 0.985: continue
            rs = (1-s)*(1-s2)/((1+s)*(1+s2)); chi = np.log((1-s2)*(1+s)/((1+s2)*(1-s)))
            pts.append((t, u, rs, chi))
    return pts

if __name__ == "__main__":
    pr("# (1) H_0 (leading lightcone function) from the LC/tail representation (setting A) vs mpmath 2F1: relative deviation")
    for (Delta, J, c1, c2) in ((14.166224, 2, 0.625, 0.625), (105.363036, 20, 0.625, -0.375), (115.107554, 52, 1.625, 1.625), (109.825411, 82, 0.375, 0.375)):
        lb = LCBlock(Delta, J, c1, c2, **PAR_A)
        row = []
        for eps in ('1e-2', '1e-4', '1e-6', '1e-8'):
            with mp.workdps(40):
                x = 1 - mp.mpf(eps)
                z = 4*x/(1 + x)**2
                bt = mp.mpf(Delta) + J
                k = z**(bt/2)*mp.hyp2f1(bt/2 + c1, bt/2 + c2, bt, z)
                H0 = mp.mpf(4)**(-(mp.mpf(Delta) + J)/2)*x**(-(mp.mpf(Delta) - J)/2)*k
            Wv = lb.W(np.array([float(x)]))[0, 0, 0]
            row.append(f"1-x={eps}: {float(Wv/H0 - 1):+.1e}")
        pr(f"  (Delta,J)=({Delta},{J}) (c1,c2)=({c1},{c2}):  " + "  ".join(row))
    pr("# (2) S[k,l] (k+l <= 4) from LC/tail (A) vs radial series with N = 24000 at its converged points: max relative deviation")
    for (Delta, J, c1, c2) in ((14.166224, 2, 0.625, -0.375), (24.221511, 10, 0.625, 0.625), (24.221511, 10, 1.625, -0.375)):
        big = CBlock(Delta, J, c1, c2, 24000)
        pts = corner_points(Delta, J, (2.0, 5.0, 10.0, 20.0), (0.002, 0.005, 0.01, 0.03, 0.1))
        r = np.array([p[2] for p in pts]); chi = np.array([p[3] for p in pts])
        X = -24000*np.log(r*np.exp(chi))
        sel = X > 90
        Sref = big.series(r[sel], chi[sel], kmax=4, lmax=4, cut=60.0)
        del big
        lb = LCBlock(Delta, J, c1, c2, **PAR_A)
        S = lb.series(r[sel], chi[sel])
        dev = max(abs(S[p, k, l]/Sref[p, k, l] - 1) for p in range(sel.sum()) for k in range(5) for l in range(5 - k))
        pr(f"  (Delta,J)=({Delta},{J}) (c1,c2)=({c1},{c2}): {int(sel.sum())} points with 1-rho_b in [{(1-r*np.exp(chi))[sel].min():.1e}, {(1-r*np.exp(chi))[sel].max():.1e}]: max dev {dev:.1e}")
    pr("# (3) deep corner: Casimir residual (relative to the largest term) for settings A / B, and max rel. difference of S[k,l] between A and B")
    for (Delta, J, c1, c2) in ((115.107554, 52, 0.625, -0.375), (105.363036, 20, 0.375, 0.375), (14.166224, 2, 0.625, -0.375)):
        pts = corner_points(Delta, J, (1.0, 3.0, 8.0, 20.0), (1e-4, 1e-3, 1e-2, 5e-2))
        r = np.array([p[2] for p in pts]); chi = np.array([p[3] for p in pts])
        xs, xb = r*np.exp(-chi), r*np.exp(chi)
        res = {}
        for name, par in (("A", PAR_A), ("B", PAR_B)):
            lb = LCBlock(Delta, J, c1, c2, **par)
            E = lb.series(r, chi, return_E=True); S = lb.series(r, chi)
            res[name] = (S, casimir_residual(E, xs, xb, Delta, J, c1, c2))
        pr(f"  (Delta,J)=({Delta},{J}) (c1,c2)=({c1},{c2})")
        for p, (t, u, _, _) in enumerate(pts):
            d = max(abs(res['A'][0][p, k, l]/res['B'][0][p, k, l] - 1) for k in range(5) for l in range(5 - k))
            pr(f"     t={t:5.1f} u={u:7.4f} 1-rho_b={1-xb[p]:.1e}:  Casimir {res['A'][1][p]:.1e} / {res['B'][1][p]:.1e}   |S_A/S_B - 1| {d:.1e}")
    open(os.path.join(HERE, 'out_c1_block_checks.txt'), 'w').write("\n".join(out) + "\n")

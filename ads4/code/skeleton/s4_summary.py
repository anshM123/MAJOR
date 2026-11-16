"""
s4: summary of the skeleton test (SKELETON_TEST.md).  Reads the two-pair moments (cache/s1_*_std.pkl), computes the exact scalar
skeletons (skeleton.Skeleton) for Df in DFS, and reports, per case (beta, m), k, kernel (even/odd) and helicity class (same/opp):

  S0 = S(x = 0) = Q_h R_spin / (Nrat R_skel)  at nu = 0                         -> m^2 (S0 - 1): nu-independent spin factor at 1/m^2
  d(x) = (Q_h R_spin - S0 Nrat R_skel H_c) Kpred / (A E_k(x)/2)                  -> m^2 d(x): nu/m-NON-uniformity of the spin factor,
        relative to the bulk-point envelope E_k(x) = (2 - x^2)/(1 - x^2)^{1 + k/2} and the local amplitude A of the leading Regge value
        (finite at zeros of P); H_c: Kpred -> Kpred + n00 [P(eta - c/m^2) - P(eta)], the curved flat-space helicity term P(eta - c/m^2)/P(eta) applied to the preserved part
        (c = 2 J_z^2; J_z = 0 for 'same', 2 for 'opp'; SPINNING_E2E.md 6.3, inferred from flat space), c = 0 = no helicity factor.
Points whose tail estimate (dropped small-u region) exceeds TAILMAX are not used.  Fits are diagnostics only.
"""
import os, sys, glob, pickle
from common import CACHE, info
import numpy as np, mpmath as mp
import spin
from skeleton import Skeleton, Nrat
from s3_compare import envelope

DFS = (1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5)
TAILMAX = 1e-4
XS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

def Ek(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2)

_H = {}
def dPc(c, nu, Delta, J):
    """P_{-1/2+i nu}(eta - c/m^2) - P_{-1/2+i nu}(eta)  (so that the helicity factor never divides by P)"""
    if c == 0:
        return 0.0
    key = (c, round(nu, 6), Delta, J)
    if key not in _H:
        m, beta, eta = info(Delta, J)
        P = lambda e: mp.re(mp.hyp2f1(0.5 - 1j*nu, 0.5 + 1j*nu, 1, (1 - e)/2))
        _H[key] = float(P(eta - mp.mpf(c)/(m*m)) - P(eta))
    return _H[key]

def cases():
    out = []
    for f in glob.glob(os.path.join(CACHE, 's1_D*_J*_std.pkl')):
        b = os.path.basename(f)[4:-8]
        D, J = b.split('_J'); out.append((float(D), int(J)))
    out = [(D, J, float(info(D, J)[0]), float(info(D, J)[1])) for D, J in out]
    return sorted(out, key=lambda t: (round(t[3], 1), t[2]))

def analyse(Delta, J, ks=(2, 3, 4, 6), dfs=DFS, cvals=(0, 8, 'flat', 'flat2')):
    """returns dict (k, kern, h, Df, c) -> dict(S0=..., pts=[(x, m2d, tail)])"""
    m, beta, eta = info(Delta, J); m = float(m)
    sp = pickle.load(open(os.path.join(CACHE, f"s1_D{Delta:.6f}_J{J}_std.pkl"), 'rb'))
    sk = Skeleton(Delta, J, 'std')
    NR = {Df: Nrat(Delta, J, Df) for Df in dfs}
    Qh = {h: spin.Q(Delta, J, h) for h in ('same', 'opp')}
    res = {}
    for k in ks:
        nus = sorted(nu for (kk, nu) in sp if kk == k)
        if not nus or nus[0] != 0.0:
            continue
        for kern in ('even', 'odd'):
            for h in ('same', 'opp'):
                for Df in dfs:
                    R0, tail0, Kp0 = sp[(k, 0.0)][(kern, h)]
                    _, _, n000, _ = envelope(Delta, J, k, 0.0, kern, h)
                    for c in cvals:
                        cc = 0 if c in ('flat', 'flat2') else c
                        S0 = Qh[h]*R0*Kp0/(NR[Df]*sk.R(Df, k, 0.0)*(Kp0 + n000*dPc(cc, 0.0, Delta, J)))
                        pts = []
                        for nu in nus:
                            R, tail, Kp = sp[(k, nu)][(kern, h)]
                            x = nu/m
                            _, A, n00, P = envelope(Delta, J, k, nu, kern, h)
                            d = (Qh[h]*R*Kp - S0*NR[Df]*sk.R(Df, k, nu)*(Kp + n00*dPc(cc, nu, Delta, J)))/(A*Ek(k, x)/2)
                            if c in ('flat', 'flat2'):
                                d -= S0*n00*Dflat(h, J, x, c)/A
                            pts.append((x, m*m*d, tail))
                        res[(k, kern, h, Df, c)] = dict(S0=S0, pts=pts, tail0=tail0)
    sk.save()
    return res, m, float(beta)

def Dflat(h, J, x, variant='flat'):
    """flat-space helicity density minus the scalar one, at cos(theta) = 1 - 2 x^2 [I: photon analogue of FLAT_GRAVITON_BOUND.md (A.1)-(A.3)]:
       same (J_z = 0): d~^J_00 = P_J  -> 0;   opp (J_z = 2): d~^J_22 = P^{(0,4)}_{J-2}  -> P^{(0,4)}_{J-2} - P_J."""
    from scipy.special import eval_jacobi, eval_legendre
    if h == 'same':
        return 0.0
    z = 1 - 2*x*x
    pre = 1.0 if variant == 'flat' else (1 - x*x)**2           # 'flat2': d^J_22 = cos^4(theta/2) d~^J_22
    return float(pre*eval_jacobi(J - 2, 0, 4, z) - eval_legendre(J, z))

def at_x(pts, x0, tol=0.026):
    c = [p for p in pts if abs(p[0] - x0) < tol and p[2] < TAILMAX]
    if not c:
        return None
    return min(c, key=lambda p: abs(p[0] - x0))[1]

def fmt(v, w=7):
    return ' ' * (w - 1) + '-' if v is None else f"{v:{w}.1f}"

if __name__ == "__main__":
    allres = {}
    for (D, J, m, beta) in cases():
        allres[(D, J)] = analyse(D, J)
        print(f"# loaded case beta={beta:.2f} m={m:.2f} (Delta={D}, J={J})", flush=True)
    pickle.dump(allres, open(os.path.join(CACHE, 's4_allres.pkl'), 'wb'))
    for k in (6, 4, 3, 2):
        for kern in ('odd', 'even'):
            for h in ('same', 'opp'):
                for c in (0, 8, 'flat', 'flat2'):
                    if h == 'same' and c != 0:
                        continue
                    print(f"\n## k={k} kernel={kern} helicity={h}  helicity factor c={c}")
                    print("#  beta     m    Df    m^2(S0-1) | m^2 d(x) at x = " + "  ".join(f"{x:5.2f}" for x in XS) + " | max|m^2 d| (x<=0.65, tail ok)")
                    for (D, J), (res, m, beta) in allres.items():
                        for Df in DFS:
                            key = (k, kern, h, Df, c)
                            if key not in res:
                                continue
                            r = res[key]
                            ok = [p[1] for p in r['pts'] if p[2] < TAILMAX]
                            print(f"   {beta:.2f} {m:7.2f}  {Df:4.2f}  {m*m*(r['S0']-1):9.2f} | " + " ".join(fmt(at_x(r['pts'], x)) for x in XS)
                                  + f" | {max(abs(v) for v in ok) if ok else float('nan'):7.1f}")

    # ---- compact tables used in SKELETON_TEST.md
    print("\n\n######## TABLE 1: C_s,eff(X) = max_{x <= X} |m^2 (S0 - 1) + m^2 d(x)| / cosh(beta)  (tail-ok points; c = 0; odd kernel | even kernel)")
    for h, c in (('same', 0), ('opp', 0), ('opp', 8)):
        for k in (6, 4, 3, 2):
            print(f"## helicity={h} c={c} k={k}:  beta   m   | Df=2: X=0.2 0.4 0.6 | Df=2.5: X=0.2 0.4 0.6 | Df=3: X=0.2 0.4 0.6   [odd ; even]")
            for (D, J), (res, m, beta) in allres.items():
                cells = []
                for kern in ('odd', 'even'):
                    row = []
                    for Df in (2.0, 2.5, 3.0):
                        key = (k, kern, h, Df, c)
                        if key not in res:
                            row.append('   -    -    -'); continue
                        r = res[key]; base = m*m*(r['S0'] - 1)
                        vals = []
                        for X in (0.2, 0.4, 0.6):
                            pts = [p for p in r['pts'] if p[0] <= X + 1e-4 and p[2] < TAILMAX]
                            cover = max([p[0] for p in pts], default=-1)
                            vals.append(f"{max(abs(base + p[1]) for p in pts)/np.cosh(beta):5.0f}" if pts and cover >= X - 0.026 else '    -')
                        row.append(''.join(vals))
                    cells.append(' |'.join(row))
                print(f"   {beta:.2f} {m:6.1f} | " + '  ;  '.join(cells))
    print("\n######## TABLE 2: same helicity, odd kernel: Df* where S0 = 1 (linear interpolation in Df) and Df minimising max|m^2 d| (x <= 0.65, tail ok)")
    for k in (6, 4, 3, 2):
        for (D, J), (res, m, beta) in allres.items():
            s0 = [(Df, m*m*(res[(k, 'odd', 'same', Df, 0)]['S0'] - 1)) for Df in DFS if (k, 'odd', 'same', Df, 0) in res]
            if not s0:
                continue
            root = None
            for (a, fa), (b, fb) in zip(s0, s0[1:]):
                if fa == 0 or fa*fb < 0:
                    root = a + (b - a)*fa/(fa - fb); break
            mx = [(Df, max([abs(p[1]) for p in res[(k, 'odd', 'same', Df, 0)]['pts'] if p[2] < TAILMAX], default=np.nan)) for Df, _ in s0]
            best = min(mx, key=lambda t: t[1])
            xmax = max([p[0] for p in res[(k, 'odd', 'same', 2.5, 0)]['pts'] if p[2] < TAILMAX], default=0)
            print(f"   k={k} beta={beta:.2f} m={m:6.1f}:  Df*(S0=1) = {root if root is None else round(root, 3)}   argmin_Df max|m^2 d| = {best[0]} ({best[1]:.2f}; "
                  + ", ".join(f"Df={Df}: {v:.2f}" for Df, v in mx if Df in (2.0, 2.5, 3.0)) + f")   x_max(tail ok) = {xmax:.2f}")

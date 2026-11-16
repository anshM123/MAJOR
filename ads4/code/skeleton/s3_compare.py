"""
s3: exact two-pair spinning moments versus exact scalar skeletons.
For each case (Delta, J), k, nu, kernel (even / odd) and helicity h (same / opp):
   S(Df) = Q_h R_spin / (Nrat(Df) R_skel(Df))     (MFT-normalised exact spinning moment / MFT-normalised exact skeleton, each relative
                                                  to its own leading Regge formula; = 1 at leading order for every Df)
   delta(Df) = (Q_h R_spin - Nrat R_skel) * Kpred / A,   A = |n00| A_P + |Kpred - n00 P|      (envelope-normalised difference:
               finite at zeros of P; |delta| = |S - 1| away from them; n00 P = preserved part of the leading prediction,
               A_P = amplitude of P_{-1/2+i nu}(cosh beta) (Harish-Chandra form, = |P| for nu beta < 1.5))
Prints m^2 delta and S.  Usage: python3 s3_compare.py case_spec ... where case_spec = Delta:J[:grid]
"""
import sys, os, pickle
from common import CACHE, info
import numpy as np, mpmath as mp
import spin
from skeleton import Skeleton, Nrat
from predict4 import prediction
from t5_analysis_l1 import hel_vectors

DFS = (1.5, 2.0, 2.5, 3.0)

def A_P(nu, beta):
    nu = mp.mpf(nu); ch = mp.cosh(beta)
    P = mp.re(mp.hyp2f1(0.5 - 1j*nu, 0.5 + 1j*nu, 1, (1 - ch)/2))
    if nu*beta < 1.5:
        return abs(float(P)), float(P)
    Z = 2*mp.gamma(1j*nu)/(mp.sqrt(mp.pi)*mp.gamma(0.5 + 1j*nu))*(2*ch)**(-0.5 + 1j*nu)*mp.hyp2f1(0.25 - 0.5j*nu, 0.75 - 0.5j*nu, 1 - 1j*nu, 1/ch**2)
    return max(abs(float(P)), float(abs(Z))), float(P)

_env = {}
def envelope(Delta, J, k, nu, kern, h):
    key = (Delta, J, k, round(nu, 6), kern, h)
    if key not in _env:
        m, beta, eta = info(Delta, J)
        o, s = hel_vectors(Delta, J); c = s if h == 'same' else o
        w = spin.KERNELS[kern]
        tot = {}
        for qL, wl in w.items():
            for qR, wr in w.items():
                p = prediction(1, Delta, J, k, nu, (qL[0], qL[1], qR[0], qR[1]), c, c, parts=True)
                for kk, v in p.items():
                    tot[kk] = tot.get(kk, 0) + wl*wr*complex(v)
        Kp = sum(tot.values()).real
        AP, P = A_P(nu, beta)
        n00 = tot.get((0, 0), 0).real/P
        _env[key] = (Kp, abs(n00)*AP + abs(Kp - n00*P), n00, P)
    return _env[key]

def load_spin(Delta, J, grid):
    fn = os.path.join(CACHE, f"s1_D{Delta:.6f}_J{J}_{grid}.pkl")
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}

def compare(Delta, J, grid='std', ks=(4, 6), quiet=False):
    """returns list of rows: (k, nu, x, kern, h, R_spin, tail, {Df: (S, delta)})"""
    m, beta, eta = info(Delta, J); m = float(m); beta = float(beta)
    sp = load_spin(Delta, J, grid)
    sk = Skeleton(Delta, J, 'std')
    NR = {Df: Nrat(Delta, J, Df) for Df in DFS}
    Qh = {h: spin.Q(Delta, J, h) for h in ('same', 'opp')}
    rows = []
    for (k, nu), r in sorted(sp.items()):
        if k not in ks:
            continue
        for kern in ('even', 'odd'):
            for h in ('same', 'opp'):
                R, tail, Kp = r[(kern, h)]
                Kp2, A = envelope(Delta, J, k, nu, kern, h)[:2]
                out = {}
                for Df in DFS:
                    Rs = sk.R(Df, k, nu)
                    a, b = Qh[h]*R, NR[Df]*Rs
                    out[Df] = (a/b, (a - b)*Kp/A)
                rows.append((k, nu, nu/m, kern, h, R, tail, out, Qh[h], {Df: NR[Df] for Df in DFS}))
    sk.save()
    if not quiet:
        print(f"# case Delta={Delta} J={J} m={m:.3f} beta={beta:.3f} grid={grid};  Q_same={Qh['same']:.6f} Q_opp={Qh['opp']:.6f};  "
              + " ".join(f"Nrat({Df})={NR[Df]:.6f}" for Df in DFS))
        print("#  k   nu     x    kern h     R_spin      tail    | m^2 delta for Df = " + "  ".join(f"{Df:>7}" for Df in DFS) + "  |  S(Df=2)  S(Df=3)")
        for (k, nu, x, kern, h, R, tail, out, _, _) in rows:
            print(f"   {k} {nu:5.2f}  {x:.3f}  {kern:4s} {h:4s} {R:11.6f}  {tail:.1e} | " + "  ".join(f"{m*m*out[Df][1]:7.2f}" for Df in DFS)
                  + f"  |  {out[2.0][0]:8.5f} {out[3.0][0]:8.5f}")
    return rows, m, beta

if __name__ == "__main__":
    for spec in sys.argv[1:]:
        p = spec.split(':')
        compare(float(p[0]), int(p[1]), p[2] if len(p) > 2 else 'std')

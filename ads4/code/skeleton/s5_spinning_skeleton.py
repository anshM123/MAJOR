"""
s5 (follow-up): "spinning skeleton" candidates for the opposite-helicity (and, as a control, same-helicity) component.
The beta-dependence P_{-1/2+i nu}(cosh beta) of the preserved part of the leading value is replaced by a harmonic-analysis structure X_nu(beta)
of SPINNING_KERNELS.md (Theorem K) / REGGE_FROM_BLOCKS.md, normalised at nu = 0 (x -> 0) to the established leading formula (no fit):
    Kp_X(nu) = Kp(nu) - n00 P(nu) + n00 P(0) X(nu)/X(0)
    d_X(x)   = (Q_h R_spin Kp - S0 Nrat R_skel Kp_X) / (A E_k(x)/2)          (S0 at nu = 0, unchanged; A, E_k as in s4)
Candidates X:  P (= s4, control);  Phi^nu_{1,1}, Phi^nu_{2,2} (Chebyshev-Mehler, closed 2F1 form);  K^nu_{2,2}, K^nu_{2,-2} (two-pair kernel of
fourpt/moments4 with one flip at each pair, mu = +-2);  (nu^2+1/4) P^{-1};  (nu^2+1/4)(nu^2+9/4) P^{-2}  (weight-shifted associated Legendre).
Usage: python3 s5_spinning_skeleton.py   (reads cache/s1_*_std.pkl; writes cache/out_s5.txt)
"""
import os, sys, pickle
from common import CACHE, info
import numpy as np, mpmath as mp
import spin
from skeleton import Skeleton, Nrat
from s3_compare import envelope
from s4_summary import cases, Ek, TAILMAX
from moments4 import kernel_mp
mp.mp.dps = 20

def Xfun(name, nu, beta):
    nu = mp.mpf(nu); b = mp.mpf(beta); ch = mp.cosh(b)
    if name == 'P':
        return mp.re(mp.hyp2f1(0.5 - 1j*nu, 0.5 + 1j*nu, 1, (1 - ch)/2))
    if name in ('Phi11', 'Phi22'):
        mm = 1 if name == 'Phi11' else 2
        return mp.re(mp.cosh(b/2)**(2*mm)*mp.hyp2f1(0.5 + mm + 1j*nu, 0.5 + mm - 1j*nu, 1, -mp.sinh(b/2)**2))
    if name == 'K22':
        return mp.re(kernel_mp(float(nu), 2, 2, b))
    if name == 'K2m2':
        return mp.re(kernel_mp(float(nu), 2, -2, b))
    if name == 'P-1':
        return mp.re((nu**2 + 0.25)*mp.legenp(-0.5 + 1j*nu, -1, ch, type=3))
    if name == 'P-2':
        return mp.re((nu**2 + 0.25)*(nu**2 + 2.25)*mp.legenp(-0.5 + 1j*nu, -2, ch, type=3))
    raise ValueError(name)

NAMES = ('P', 'Phi11', 'Phi22', 'K22', 'P-1', 'P-2')   # K^nu_{2,-2} of moments4 equals Phi22 (inversion-referred frame), so K22 is the off-diagonal Phi_{2,-2}
XS = (0.1, 0.2, 0.3, 0.4, 0.6)

if __name__ == "__main__":
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr("# s5: spinning-skeleton candidates X replacing P in the preserved leading part (normalised at nu = 0). m^2 d_X at x = "
       + ", ".join(str(x) for x in XS) + " | max|m^2 d_X| over tail-ok x <= 0.65")
    Xc = {}
    for (D, J, m, beta) in cases():
        sp = pickle.load(open(os.path.join(CACHE, f"s1_D{D:.6f}_J{J}_std.pkl"), 'rb'))
        sk = Skeleton(D, J, 'std')
        Qh = {h: spin.Q(D, J, h) for h in ('same', 'opp')}
        for k in (6, 4):
            nus = sorted(nu for (kk, nu) in sp if kk == k)
            for h in ('opp', 'same'):
                for kern in ('odd', 'even'):
                    for Df in (2.0, 2.5):
                        NR = Nrat(D, J, Df)
                        R0, _, Kp0 = sp[(k, 0.0)][(kern, h)]
                        S0 = Qh[h]*R0/(NR*sk.R(Df, k, 0.0))
                        for name in NAMES:
                            vals = []
                            for nu in nus:
                                R, tail, Kp = sp[(k, nu)][(kern, h)]
                                if tail >= TAILMAX:
                                    continue
                                x = nu/m
                                _, A, n00, P = envelope(D, J, k, nu, kern, h)
                                for key in ((name, nu), (name, 0.0)):
                                    kk_ = key + (round(beta, 6),)
                                    if kk_ not in Xc:
                                        Xc[kk_] = float(Xfun(key[0], key[1], beta))
                                X, X0 = Xc[(name, nu, round(beta, 6))], Xc[(name, 0.0, round(beta, 6))]
                                P0 = Xc[('P', 0.0, round(beta, 6))] if ('P', 0.0, round(beta, 6)) in Xc else float(Xfun('P', 0.0, beta))
                                KpX = Kp - n00*P + n00*P0*X/X0
                                d = (Qh[h]*R*Kp - S0*NR*sk.R(Df, k, nu)*KpX)/(A*Ek(k, x)/2)
                                vals.append((x, m*m*d))
                            row = []
                            for x0 in XS:
                                c = [v for v in vals if abs(v[0] - x0) < 0.026]
                                row.append(f"{min(c, key=lambda v: abs(v[0]-x0))[1]:8.1f}" if c else "       -")
                            mx = max(abs(v[1]) for v in vals) if vals else float('nan')
                            pr(f"  b={beta:.2f} m={m:6.1f} k={k} {h:4s} {kern:4s} Df={Df} X={name:5s} |" + "".join(row) + f" | {mx:9.1f}")
        sk.save()
    open(os.path.join(CACHE, 'out_s5.txt'), 'w').write("\n".join(out) + "\n")

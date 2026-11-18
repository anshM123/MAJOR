"""nf1: flat limit of Omega^{NF,SC} (NF02, no forward terms).
(a) factored psi2, psi3 == exact expanded ones (exact rational evaluation at test points); ExactCosT vs rcommon.CosT (float) at small tau.
(b) per-k B' parts sum to the stripped NF02 action (pointwise identity of u1, generic psi; numerical check here with factored weights);
    flat min F/N_s on a grid (NOFORWARD_FUNCTIONAL.md: 0.2000 at m = M, J = 4 (++), 5 (+-)).
(c) amplification N_U/N_s; odd-k prefactor sensitivity (odd B' density -> k-continued scalar E_k).
Output: cache/nf1.log
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from fractions import Fraction as Fr
from nf_common import NF, ExactCosT, flat_parts_nf, KS, EPS, SIGMA, Dh_all
from common import Ek_scalar
from rcommon import CosT
out = []
def pr(s):
    print(s, flush=True); out.append(s)
F = NF()
pr(f"NF02: k(0) = {F.k0} = {float(F.k0):.8f}; normalised x0 = psi2'(0) = {F.x0:.6f}; psi3(0) = {float(F.s3[0])}")
for pv in (Fr(1, 3), Fr(9, 10), Fr(99, 100)):
    e2 = sum(c*pv**i for i, c in enumerate(F.s2)); f2 = F.psi2(float(pv))
    e3 = sum(c*pv**i for i, c in enumerate(F.s3)); f3 = F.psi3(float(pv))
    pr(f"(a) p={pv}: psi2 exact {float(e2):.12g} factored {float(f2):.12g}; psi3 exact {float(e3):.12g} factored {float(f3):.12g};"
       f" expanded-float psi2 = {np.polynomial.polynomial.polyval(float(pv), [float(c) for c in F.s2]):.12g}")
for k in (2, 5):
    c1 = ExactCosT(F.PHx[k]); c2 = CosT(np.array([float(c) for c in F.PHx[k]]))
    t = np.array([0.5, 3.0, 7.9, 8.1, 20.0, 60.0])
    pr(f"(a) k={k}: ExactCosT {np.array2string(c1(t), precision=8)}  float CosT {np.array2string(c2(t), precision=8)}")
from numpy.polynomial.legendre import leggauss
xg, wg = leggauss(1500); pn = (xg + 1)/2; pw = wg/2
worst = {}
for h in SIGMA:
    for mfl in (1.0, 1.0005, 1.01, 1.05, 1.2, 1.4, 1.5, 2.0, 3.0, 6.0, 20.0):
        u = 1/mfl**2
        Js, B2, B3, parts = flat_parts_nf(F, h, u, 60, pn, pw)
        sel = (Js % 2 == 0) if h == '++' else np.ones_like(Js, bool)
        Fs = B2 + B3; Ns = np.abs(B2) + np.abs(B3); NU = np.abs(parts).sum(axis=1)
        dev = np.max(np.abs(parts.sum(axis=1) - Fs)[sel]/Ns[sel])
        y = u*pn*pn; _, Dm = Dh_all(h, 60, y)
        dodd = sum(u**(k - 2)*(Dm @ (pw*F.phi(k, pn)*SIGMA[h]*(EPS[k]*Ek_scalar(k, y) - 1/(1 - y)**((k + 1)/2)))) for k in KS if k % 2)
        r = Fs[sel]/Ns[sel]; i = np.argmin(r)
        pr(f"(b,c) {h} m/Lam={mfl:7.4f}: min F/N_s = {r[i]:+.5f} (J={Js[sel][i]}); |sum_k B' - stripped|/N_s <= {dev:.1e}; max N_U/N_s = {np.max(NU[sel]/Ns[sel]):.3g};"
           f" odd-prefactor shift max/N_s = {np.max(np.abs(dodd[sel])/Ns[sel]):.3g}")
open(os.path.join(HERE, 'cache', 'nf1.log'), 'w').write("\n".join(out) + "\n")

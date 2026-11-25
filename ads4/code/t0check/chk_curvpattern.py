"""t0check/chk_curvpattern.py (INDEP_CHECK_T0 §3-4) [DV]:
  1. Lemma 5 (flat-Regge remainder): |int psi J0(bp) dp + f1/b^3| <= e(b)/b^3 at sample b (mpmath quadrature).
  2. The CONCLUSIONS of Parts K/T/L on actual values (chk_eval.curv, not bounds): for every sign pattern g, sigma = +-1, u in {0, 1},
     curv_g := (1 - c s2) 2 C2 + (1 - c s3) u sigma C3 is >= 0 for bhat >= 16 and >= -T_abs for bhat <= 16, at Lam(300), Lam(1000),
     Lam(1e4), on a bhat grid reaching beta ~ 12 (far beyond 2 Lam).  Also reports where curv_g first becomes positive."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, mpmath as mp
import chk_eval as E, t0_common as T

print('1. Lemma 5:')
for k in (2, 3):
    f = T.flat_E_consts(k)
    f1 = float(f['f1_signed']); ps = [mp.mpf(c.numerator)/c.denominator for c in T.PSI[k]]
    for b in (3.5, 10.0, 50.0, 300.0, 2000.0):
        val = mp.quad(lambda p: sum(c*p**i for i, c in enumerate(ps))*mp.besselj(0, b*p), mp.linspace(0, 1, 2 + int(b/2)))
        e = float(f['f1'])*(1 + 2/b) + float(f['g1']) + (float(f['dg1']) + float(f['Iq'].upper()))/b
        rem = float(val) + f1/b**3
        print(f'   k={k} b={b:7.1f}: Omega_flat = {float(val): .6e}, -f1/b^3 = {-f1/b**3: .6e}, |E_fl| b^3 = {abs(rem)*b**3:.4g} <= e(b) = {e:.4g}: {abs(rem)*b**3 <= e}')
print('2. actual pattern curvature:')
c = 0.199
for Dg, Tabs in ((300, 3.3036e-4), (1000, 2.9075e-5), (10000, None)):
    L = E.lam(Dg)
    bs = np.concatenate([np.linspace(0.05, 16, 160), np.linspace(16, 300, 285)[1:], np.geomspace(300, 2*L*np.sinh(6.0), 80)[1:]])
    C2 = np.array([E.curv(2, b, L) for b in bs]); C3 = np.array([E.curv(3, b, L) for b in bs])
    mn_lo = 1e9; mn_hi = 1e9; firstpos = 0
    for s2 in (1, -1):
        for s3 in (1, -1):
            for sg in (1, -1):
                for u in (0, 1):
                    cv = (1 - c*s2)*2*C2 + (1 - c*s3)*u*sg*C3
                    mn_lo = min(mn_lo, cv[bs <= 16].min()); mn_hi = min(mn_hi, (cv[bs >= 16]*L*L).min())
                    neg = bs[cv < 0]
                    if len(neg): firstpos = max(firstpos, neg.max())
    print(f'   Dg={Dg}: min over bhat<=16 of curv_g = {mn_lo:.4e} (T_abs = {Tabs}); min over bhat>=16 (to beta=12) of Lam^2 curv_g = {mn_hi:.4f}; '
          f'largest bhat with some curv_g < 0: {firstpos:.2f}')

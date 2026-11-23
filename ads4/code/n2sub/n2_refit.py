"""
n2_refit.py (N2_SUBREGGE.md §2.4): HWIN_THIMBLE §4.2's D_h-class fit of the (1100, 100) window data, redone in the convention in which
the data are defined (raw, forward-normalised; the exact weight is X).

Data: spinfunc/cache/pw_gamma.pkl (read-only): per node nu, Gamma (G), kappa image (img), y = G - img, jackknife err.
  F0  : y  = c0 D + s x D'/(2J+1) + r                                (fw_gamma.fit 'cs', unweighted; reproduces the stored fit)
  F1  : y' = c0 A[D] + s A[x D'/(2J+1)] + r A[1]                     with y' = y - (A - 1)[kappa]      (all A-columns: A[env f]/env, Re)
        -- the correct raw form: a stripped strip-analytic correction f enters the rows as A[env f]/env, and so does kappa.
  F2  : y  with A-columns but WITHOUT the (A - 1) kappa correction.
  F3  : joint C2 / C3H fit of the D_h class with ONE c0 (C3H has no de-shift; its class is sigma_h c0 D): tests the a_h-normalisation
        reading (a_h is common to k = 2, 3 by construction, l9_analysis docstring).
Also printed: |a(nu)| = |R_- + R_+|/(|R_-| + |R_+|) at the nodes (conversion Gamma_c -> a Gamma_c for the |Y|-normalised tolerance).
Usage: python3 n2_refit.py -> cache/n2_refit.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N

d = pickle.load(open(os.path.join(N.ADS4, 'code', 'spinfunc', 'cache', 'pw_gamma.pkl'), 'rb'))
lines = []


def pr(s=''):
    print(s, flush=True); lines.append(s)


def Acol(fx, nu, m, h):
    """Re A[env f]/env at real nu (f given as function of complex x)."""
    e = lambda x: N.env(h, x, 2)
    return (N.A_on_x(lambda x: e(x)*fx(x), nu, m)/e(nu/m)).real


def dprime(h, J, x, hh=1e-6):
    return (N.Dsh(h, J, x + hh) - N.Dsh(h, J, x - hh))/(2*hh)


def lsq(cols, y):
    A = np.array(cols).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ c
    dof = max(len(y) - A.shape[1], 1)
    cov = np.linalg.inv(A.T @ A)*(res @ res)/dof
    return c, np.sqrt(np.diag(cov)), np.sqrt(np.mean(res**2)), res


KEYS = [(100, 1100.0, 12.0, 34.0, 'cv', ('_w',)), (100, 1100.0, 0.0, 34.0, 'cv', ('_w',)), (100, 1100.0, 0.0, 34.0, 'dir', ('_w',)), (100, 1100.0, 0.0, 34.0, 'cv', ('_wP',)),
        (100, 1100.0, 0.0, 34.0, 'cv', ('_none',)), (64, 1064.0, 12.0, 34.0, 'cv', ('_w',)), (64, 364.0, 0.0, 34.0, 'cv', ('_w',))]
Rm, Rp = N.Rk()
def main():
  for key in KEYS:
      if key not in d:
          continue
      r = d[key]; J = key[0]; m = r['m']; nus = np.asarray(r['nus'])
      sel = (nus >= 12 - 1e-9) & (nus <= 34 + 1e-9)
      nu = nus[sel]; x = nu/m
      aabs = np.abs(Rm(nu + 0j) + Rp(nu + 0j))/(np.abs(Rm(nu + 0j)) + np.abs(Rp(nu + 0j)))
      pr(f"#### J={J} Delta={key[1]:g} m={m:.2f} est={key[4]} tags={key[5]}  nodes nu={np.round(nu, 2).tolist()}")
      pr(f"     |a(nu)| at nodes = {np.round(aabs, 3).tolist()}")
      c0s = {}
      for h in ('++', '+-'):
          k = ('C2', h)
          if k not in r:
              continue
          y = np.asarray(r[k]['y'])[sel]; err = np.asarray(r[k]['err'])[sel]
          D = N.Dsh(h, J, x).real; Dp = (dprime(h, J, x)*x/(2*J + 1)).real
          c, e, rms, res = lsq([D, Dp, np.ones_like(x)], y)
          AD = Acol(lambda z: N.Dsh(h, J, z), nu, m, h)
          ADp = Acol(lambda z: dprime(h, J, z)*z/(2*J + 1), nu, m, h)
          A1 = Acol(lambda z: 0*z + 1.0, nu, m, h)
          kap = N.kappa('C2', h, J, m, x).real
          Akap = Acol(lambda z: N.kappa('C2', h, J, m, z), nu, m, h)
          dk = Akap - kap
          y1 = y - dk
          c1, e1, rms1, res1 = lsq([AD, ADp, A1], y1)
          c2, e2, rms2, res2 = lsq([AD, ADp, A1], y)
          pr(f"  C2 {h}: y = {np.round(y).tolist()} (+- {np.round(err).tolist()})")
          pr(f"        (A-1)D/D = {np.round((AD - D)/np.where(np.abs(D) > 1e-3, D, np.nan), 3).tolist()};  (A-1)[kappa] = {np.round(dk).tolist()}")
          pr(f"     F0 (stored convention): c0 = {c[0]:+.0f}, s = {c[1]:+.0f}, r = {c[2]:+.0f} +- {e[2]:.0f}, rms resid {rms:.0f}; "
             f"max|r|+rms = {abs(c[2]) + rms:.0f}")
          pr(f"     F1 (raw form, kappa de-shifted): c0 = {c1[0]:+.0f}, s = {c1[1]:+.0f}, r = {c1[2]:+.0f} +- {e1[2]:.0f}, rms resid {rms1:.0f}; "
             f"|r|+rms = {abs(c1[2]) + rms1:.0f};  max|y'| = {np.max(np.abs(y1)):.0f};  max|a(y' - class)| = {np.max(np.abs(aabs*res1)) + 0:.0f} (+|r a|)")
          pr(f"     F2 (A-columns, kappa NOT de-shifted): c0 = {c2[0]:+.0f}, s = {c2[1]:+.0f}, r = {c2[2]:+.0f} +- {e2[2]:.0f}, rms {rms2:.0f}")
          c0s[h] = (c[0], c1[0])
          # C3H joint test (no de-shift for k = 3): y3 = sigma c0 D + s3 x D' + r3
          k3 = ('C3H', h)
          if k3 in r:
              y3 = np.asarray(r[k3]['y'])[sel]; e3 = np.asarray(r[k3]['err'])[sel]
              sg = 1.0 if h == '++' else -1.0
              if np.all(np.isfinite(y3)) and np.median(e3) < 5e3:
                  c3, ee3, rms3, _ = lsq([D, Dp, np.ones_like(x)], y3)
                  w = 1/np.maximum(e3, 1)**2
                  # weighted single-parameter test: C3H y3 - (its own s, r) ~ sigma * c0(C2)
                  pr(f"     C3H {h}: own fit c0_3 = {c3[0]:+.0f} +- {ee3[0]:.0f} (sigma_h*c0_C2 predicted {sg*c[0]:+.0f} [F0] / {sg*c1[0]:+.0f} [F1]); "
                     f"r3 = {c3[2]:+.0f}; median err {np.median(e3):.0f}")
              else:
                  pr(f"     C3H {h}: unusable (median err {np.median(e3):.0f})")
  pr('')
  # same-node comparison at (1100, 100): cv restricted to the dir node set (13.3, 16.3, 20.8, 25), F1 convention
  pr('#### (1100,100): CV vs DIR on the SAME 4 nodes (13.3, 16.3, 20.8, 25), F1 convention (1 dof)')
  for est in ('cv', 'dir'):
      r = d[(100, 1100.0, 0.0, 34.0, est, ('_w',))]; m = r['m']; nus = np.asarray(r['nus'])
      sel = np.isin(np.round(nus, 2), [13.3, 16.3, 20.8, 25.0]); nu = nus[sel]; x = nu/m
      for h in ('++', '+-'):
          y = np.asarray(r[('C2', h)]['y'])[sel]
          AD = Acol(lambda z: N.Dsh(h, 100, z), nu, m, h); ADp = Acol(lambda z: dprime(h, 100, z)*z/201, nu, m, h)
          A1 = Acol(lambda z: 0*z + 1.0, nu, m, h)
          dk = Acol(lambda z: N.kappa('C2', h, 100, m, z), nu, m, h) - N.kappa('C2', h, 100, m, x).real
          c1, e1, rms1, _ = lsq([AD, ADp, A1], y - dk)
          pr(f"   {est} {h}: c0 = {c1[0]:+.0f}, r = {c1[2]:+.0f} +- {e1[2]:.0f}, rms {rms1:.0f}")
  open(os.path.join(N.CACHE, 'n2_refit.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()

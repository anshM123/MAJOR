"""
s_cost.py (N2_SUBREGGE.md §3.3): the action of the first subleading Regge term on the H_LJ(1000) sample states, for the scalar
laboratory (CHMRS eq. (regge subleading), Dphi = 3, k = 2, validated in s_chmrs.py) -- i.e. R_win's content at O(1/m^2) where the
formula is known -- including J >= 150 where no exact data exist.

  dF_sub = 2 int_0^1 psi2(p) [delta(Lam p; beta) - 8.1 P] dp / m^2        (the -8.1 P is the measured normalisation residual of s_chmrs)
Decomposition (ring class, Lemma I as distributions at t = beta):
  delta = (1/6)[ (K4 + K2 nu^2) cosh(b) P + (1/sinh b + sinh b (nu^2 + K1)) dP/db ],   K1 = -11.75, K2 = 11, K4 = 86.75 (Dphi = 3, k = 2)
Computed two ways: (i) direct quadrature on the 0.04 nu-grid up to Lam; (ii) Lemma I: nu^2 P(cosh t) = -(d_t^2 + coth t d_t + 1/4) P, so
the action is the same operator applied to Omega_LR,2(t) at t = beta (finite differences of ls_core.OmLR).  Floors: 0.9 lam0 N_s (HLJ) and
N_flat/5 (LARGESPIN, struct_1000.json floor_rel).
Usage: python3 s_cost.py -> cache/s_cost.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N
sys.path.insert(0, os.path.join(N.ADS4, 'code', 'largespin'))
import ls_core as LC        # read-only (OmLR, hyb_fast)

Dg = 1000; L = N.lam(Dg)
XD = N.xdata(Dg); nu, ps, dn = XD['nu'], XD['psi2'], XD['dnu']
sel = nu < L; nu, ps = nu[sel], ps[sel]
K1, K2, K4 = -11.75, 11.0, 86.75
RES = -8.1
ST = json.load(open(os.path.join(N.ADS4, 'code', 'largespin', 'cache', 'struct_1000.json')))
lines = []


def pr(s=''):
    print(s, flush=True); lines.append(s)


def Om(t):
    return LC.OmLR(2, [t], Dg)[0]


pr('J  m/L  beta | dF_sub/N_s: direct (CHMRS / residual -8.1P) | Lemma I check | cost/floor: 0.9lam0 , N_flat/5 | equiv. ring-amplitude a = coeff*m (LS A-tol)')
seen = set()
for q in ST:
    if q['reg'] != 'LJ' or q['J'] in seen:
        continue
    seen.add(q['J'])
    J, m, b, Ns = q['J'], q['m'], q['beta'], q['Ns']
    P = N.Pnu(b, nu).real
    h = 1e-5
    dP = (N.Pnu(b + h, nu) - N.Pnu(b - h, nu)).real/(2*h)
    delta = ((K4 + K2*nu**2)*np.cosh(b)*P + (1/np.sinh(b) + np.sinh(b)*(nu**2 + K1))*dP)/6
    d1 = 2*(ps*delta).sum()*dn/L/m**2/Ns
    d2 = 2*(ps*RES*P).sum()*dn/L/m**2/Ns
    # Lemma I: Omega(t) and derivatives at beta
    hh = 2e-4*max(b, 0.05)
    O0, Op, Om_ = Om(b), Om(b + hh), Om(b - hh)
    O1 = (Op - Om_)/(2*hh); O2 = (Op - 2*O0 + Om_)/hh**2
    Op2, Om2 = Om(b + 2*hh), Om(b - 2*hh)
    O3 = (Op2 - 2*Op + 2*Om_ - Om2)/(2*hh**3)
    lap = lambda f0, f1, f2: -(f2 + f1/np.tanh(b) + f0/4)          # nu^2 acting on a radial function of t, evaluated at beta
    # nu^2 dP/dt  <->  d/dt of (nu^2 P) = -(d_t^3 + coth d_t^2 - d_t/sinh^2 + d_t/4)
    nu2O = lap(O0, O1, O2)
    nu2dO = -(O3 + O2/np.tanh(b) - O1/np.sinh(b)**2 + O1/4)
    li = 2*((K4*np.cosh(b)*O0 + K2*np.cosh(b)*nu2O + (1/np.sinh(b) + np.sinh(b)*K1)*O1 + np.sinh(b)*nu2dO)/6)/m**2/Ns
    st = LC.hyb_fast(Dg, J, q['mfl'], '+-')
    lam0 = st['base']/st['Ns']
    tot = d1 + d2
    pr(f"{J:3d} {q['mfl']:.3f} {b:.4f} | {d1:+.3e} / {d2:+.3e} | {li:+.3e} | {abs(tot)/(0.9*lam0):.2e} , {abs(tot)/q['floor_rel']:.2e} |"
       f" ring coeff (K4 cosh b + RES)/6..: {(K4*np.cosh(b)/6 + RES):.1f}/m^2")
open(os.path.join(N.CACHE, 's_cost.txt'), 'w').write('\n'.join(lines) + '\n')

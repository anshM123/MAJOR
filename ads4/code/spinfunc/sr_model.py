"""
sr_model (SUBLEADING_REGGE.md §4): per-state positivity of Omega^EE near the gap at large J with the Regge-resummed uniform action
   T_h(nu) = env_h(x) A(b) D_h(theta_eff) + 2[P_{-1/2+i nu}(cosh b) - A J0(nu b)]          (k = 2; k = 3 analogous with u sigma_h env3)
   theta = 2 arcsin x,  theta_eff = theta (1 - b^2 g(x)/24),  A = sqrt(b/sinh b)(1 + b^2 a(x))
 g = 1, a = 0: 'kappa' model (theta_eff = theta (b/2)/sinh(b/2) exactly in b; fixed by the exact leading-Regge limit); 'fit': g, a
 from the row fit (SUBLEADING_REGGE §2).  The correction relative to the hybrid base of ee1000_common.hyb is
   dF_h = sum p2 [T_h - env D_h(theta) - 2(P - J0((2J+1)x))] + (k=3 analogue),
 and lam = (base + dF)/N_s (diagonal; leading-Regge off-diagonal = 0 for EE).  Also: tolerance on a non-oscillating correction
 (c2 units, uniform in x) = 0.9 lam/[sens r0] and on nu < 50 alone.
Usage: python3 sr_model.py Dg J1,J2,... [npn]   -> stdout (cache/sr_model_<Dg>.txt via tee)
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import ee1000_common as E

FIT = dict(g=(1.0, 1.52, -7.288), a=(1.479, 2.153))    # ++ row fit (sr §2); used as the next-order error proxy


def corr(st, h, J, mode, k3=True):
    x, m, u, y = st['x'], st['m'], st['u'], st['y']
    l = J + 0.5; b = 2*np.arcsinh(l/m); nu = x*m
    th = 2*np.arcsin(np.minimum(x, 1 - 1e-12))
    if mode == 'kappa':
        g = 1.0; a = 0.0
    else:
        g = FIT['g'][0] + FIT['g'][1]*y + FIT['g'][2]*y*y; a = FIT['a'][0]*y + FIT['a'][1]*y*y
    te = th*(1 - b*b*g/24); A = np.sqrt(b/np.sinh(b))*(1 + b*b*a)
    xe = np.sin(te/2)
    D = E.Dc_fixed(h, J, x); De = E.Dc_fixed(h, J, xe)
    osc = A*De - D
    reg = A*j0(nu*b) - j0((2*J + 1)*x)
    s = E.Q.SIGMA[h]
    d2 = (st['p2']*(E.Q.A.env(h, x, 2)*osc - 2*reg)).sum()
    d3 = (st['p3']*u*s*(E.Q.A.env(h, x, 3)*osc - reg)).sum() if k3 else 0.0
    return d2 + d3


if __name__ == '__main__':
    Dg = int(sys.argv[1]); Js = [int(v) for v in sys.argv[2].split(',')]
    npn = int(sys.argv[3]) if len(sys.argv) > 3 else 512
    for J in Js:
        mg = E.mgap(Dg, J)
        rows = []
        for mfl in [mg] + [v for v in (1.05, 1.1, 1.15, 1.2, 1.3, 1.5, 1.75, 2.0) if v > mg + 1e-6]:
            res = []
            for h in ('++', '+-'):
                st = E.hyb(Dg, J, mfl, h, npn)
                dk = corr(st, h, J, 'kappa'); df = corr(st, h, J, 'fit')
                x = st['x']; env = E.Q.A.env(h, x, 2)
                sens = (st['p2']*env/st['m']**2).sum()/st['Ns']                      # per 1 c2 unit, uniform in x
                sl = (st['p2']*env*(x*st['m'] < 50)/st['m']**2).sum()/st['Ns']        # nu < 50 only
                l0 = st['base']/st['Ns']; lk = (st['base'] + dk)/st['Ns']; lf = (st['base'] + df)/st['Ns']
                res.append((h, l0, lk, lf, 0.9*min(lk, lf)/sens, 0.9*min(lk, lf)/sl if sl > 0 else np.inf))
            lam = min(r[2] for r in res)
            print(f"Dg={Dg} J={J:3d} m/L={mfl:.4f} bhat={(2*J+1)/mfl:6.1f} beta={2*np.arcsinh((J+.5)/(mfl*E.lam_of(Dg))):.3f}  " + " | ".join(
                f"{h}: lam0 {l0:+.3f} kappa {lk:+.3f} fit {lf:+.3f} |dk-dfit| {abs(lk-lf):.3f} tolC {tc:.0f} tolC(nu<50) {tl:.0f}"
                for h, l0, lk, lf, tc, tl in res), flush=True)

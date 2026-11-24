"""
gsm_ooid (GAP_SEAM.md §3): the identity Gamma_OO = Gamma_++ (PARITY_ODD_CLOSE §1 (D)) on the Delta_g = 2000 own rows Delta = 2000 + J,
J = 4, 6, 8, 10, from K2 points of the 3x3 p3_odd_action action at even J (gsm_rows.py oo: kind EE = w0 E (x) E integrand, C3H).
Per node: T0 = real helicity basis of the row's own even 2x2 block (l8_analysis.basis, read-only, order 0, nodes x >= 0.15);
kappa_k(x) = M_k[OO] / (T0^-1 M_k[ee] T0^-T)_++ for the Laplace orders k = 0, 2, 4.  If kappa_k(x) = kappa_0 for every k and x, the
order-by-order chain (pred is linear in the state coefficients) gives Gamma_OO = Gamma_++ exactly.  Reported:
  * max_x |kappa_0(x)/kappa_0(x_ref) - 1|  (order-0 constancy);
  * dGamma(x) = m^2 |rho(x)/rho(x_ref) - 1|,  rho = sum_k M_k[OO] / sum_k (T0^-1 M_k T0^-T)_++  (K2-summed state series; an upper
    proxy for |Gamma_OO - Gamma_++| / |D_++| in c2 units, the reference side being constant to 1e-13 at 2e4, PARITY_ODD_CLOSE).
Usage: python3 gsm_ooid.py 4,6,8,10   -> cache/gsm_ooid.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gsm_guard  # noqa: F401
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
os.chdir(SF)
import l8_analysis as L

ROWS = os.path.join(HERE, 'cache', 'rows')


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def run(J, D, lines):
    out = {}
    for kind in ('EE', 'C3H'):
        fn = os.path.join(ROWS, 'oo_p3_%s_J%d_gsm.pkl' % (kind, J))
        if not os.path.exists(fn):
            continue
        R = {(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items() if k[2] == 1.0 and k[0] == D}
        m = mphys(D, J)
        R2 = {k: dict(v, M0=np.asarray(v['M0'])[:2, :2]) for k, v in R.items()}
        B = L.basis(J, D, R2)
        if B is None:
            lines.append(f"J={J} D={D:g} {kind}: no basis ({len(R)} nodes)"); continue
        Ti = B['Ti']
        xs, k0, rho, kk = [], [], [], []
        for (DD, nu), r in sorted(R.items(), key=lambda t: t[0][1]):
            pp = [(Ti @ np.asarray(r[o])[:2, :2] @ Ti.T)[0, 0] for o in ('M0', 'M2', 'M4')]
            oo = [np.asarray(r[o])[2, 2] for o in ('M0', 'M2', 'M4')]
            xs.append(nu/m); k0.append(oo[0]/pp[0]); rho.append(sum(oo)/sum(pp)); kk.append([o/p for o, p in zip(oo, pp)])
        xs = np.array(xs); k0 = np.array(k0); rho = np.array(rho); kk = np.array(kk)
        i0 = int(np.argmin(np.abs(xs - 0.1)))
        dk0 = np.abs(k0/k0[i0] - 1); dG = m*m*np.abs(rho/rho[i0] - 1)
        dkk = np.abs(kk/kk[:, :1] - 1)                     # per-order kappa_k / kappa_0 - 1
        lines.append(f"J={J:2d} D={D:g} m={m:.1f} {kind:3s}: basis sv3 {B['sv3']:.1e} imag {B['imag']:.1e}; nodes x = "
                     + " ".join(f"{v:.3f}" for v in xs))
        lines.append(f"      kappa_0 = {k0[i0].real:+.6f}{k0[i0].imag:+.1e}i; max|kappa_0(x)/kappa_0 - 1| = {dk0.max():.2e}; "
                     f"max_k,x |kappa_k/kappa_0 - 1| = {dkk[:, 1:].max():.2e}; dGamma proxy m^2|rho/rho_ref - 1| per node: "
                     + " ".join(f"{v:.2e}" for v in dG) + f"  (max {dG.max():.3f} c2)")
        out[kind] = dict(xs=xs.tolist(), dk0=float(dk0.max()), dkk=float(dkk[:, 1:].max()), dG=dG.tolist(), dGmax=float(dG.max()))
    return out


if __name__ == '__main__':
    lines = []
    res = {}
    for J in [int(v) for v in sys.argv[1].split(',')]:
        res[J] = run(J, 2000.0 + J, lines)
    print("\n".join(lines))
    open(os.path.join(HERE, 'cache', 'gsm_ooid.txt'), 'w').write("\n".join(lines) + "\n")
    pickle.dump(res, open(os.path.join(HERE, 'cache', 'gsm_ooid.pkl'), 'wb'))

"""
End-to-end block test of the MFT normalisation (referee issue 1): leading-twist (tau = tau_0 = Dphi + 1) lightcone limit of the
MFT correlator <O_1 O_3><phi_2 phi_4> on the real section of SPINNING_E2E.md, decomposed with the exact e2e blocks
(spinblock.SpinBlock, imported, not modified) and the Fock-space coefficients p[L](N = J - l, J) of mft_exact.py.

On tau_0 the conserved structure is L = sum_n L_n S_n = 2^J sum_w c_w W_w (finite, while L_0 = 0), and by SPINNING_E2E.md sect. 1.4 an
exchange with <O phi X> = lam L contributes lam^2 (-1/2)^J W_diff.  Prediction (no parameter):
      F_MFT(x_i; z_1, z_3)  =  sum_{J >= l} p[L](J - l, J) (-1/2)^J W_diff^{(Delta = tau_0 + J, J)}  x (1 + O(u^{1/2}))
up to the continuation phases between the Euclidean and the real section (|x_ij^2|^p everywhere on the real section); these are
(a) an overall sign and (b) a possible (-1)^J (the sign between P_J(eta) and P_J(-eta), SPINNING_E2E.md sect. 1.4).  The script
prints the ratio for both choices of (b); a correct normalisation must give a ratio that is (+-1) x (1 + O(eps^{1/2})), independent of
the polarisations and of the point, for ONE of the two choices, with all J contributing.
"""
import os, sys, time
import numpy as np, mpmath as mp, sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..', 'e2e'))
from spinblock import SpinBlock, invariants, block_coords_from_uv, mdot
from mft_exact import closed_even

DPHI = 1.75

def wfull(sb, x1, x2, x3, x4, z1, z3):
    X = [np.array(v, dtype=float)[None] for v in (x1, x2, x3, x4)]
    z1 = np.array(z1, dtype=float); z3 = np.array(z3, dtype=float)
    Z1 = z1[None, :, None, :] + 0*z3[None, None, :, :]; Z3 = z3[None, None, :, :] + 0*z1[None, :, None, :]
    s, z = invariants(X[0], X[1], X[2], X[3], Z1, Z3)
    r, chi, A0, E0 = block_coords_from_uv(s)
    ls, M = sb.eval(s, z, r, chi, A0, E0)
    return (np.exp(ls - (sb.l+1)*np.log(s['s13']) - sb.Dphi*np.log(s['s24']))[:, None, None]*M)[0], s

def F_mft(l, x1, x2, x3, x4, z1, z3):
    """<O1 O3><phi2 phi4> with |x_ij^2| (real section, phases stripped): (z1.z3 x13^2 - 2 z1.x13 z3.x13)^l/|x13^2|^{2l+1}/|x24^2|^Dphi"""
    d = lambda a, b: -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    x13 = np.array(x1) - np.array(x3); x24 = np.array(x2) - np.array(x4)
    out = np.zeros((len(z1), len(z3)))
    for i, a in enumerate(z1):
        for j, b in enumerate(z3):
            out[i, j] = (d(a, b)*d(x13, x13) - 2*d(a, x13)*d(b, x13))**l/abs(d(x13, x13))**(2*l + 1)/abs(d(x24, x24))**DPHI
    return out

def config(eps, base):
    """x_12 = T (1, (1-eps) n) nearly null (s-channel lightcone x_12^2 -> 0 at fixed zbar), all other x_ij timelike, v > 1 (zbar < 0)"""
    x2, x3, x4, n, T = base
    x1 = x2 + T*np.array([1.0, *((1 - eps)*n)])
    return x1, x2, x3, x4

def main():
    bases = {'A': (np.zeros(3), np.array([5.0, 1.0, 0.2]), np.array([-4.0, 0.3, -0.5]), np.array([1.0, 0.0]), 1.0),
             'B': (np.array([0.2, -0.1, 0.3]), np.array([4.0, 0.5, 1.0]), np.array([-3.5, -0.6, 0.2]), np.array([0.6, 0.8]), 0.8)}
    print("Component z_1 = null direction of x_12 (where the tau_0 family dominates; see SPINNING_MFT.md sect. 5.1);")
    print("z_3 = that direction, (1, cos 2.6, sin 2.6), (1, -1, 0).  Ratios F_MFT / sum_J p (-1/2)^J W  [and control with (+1/2)^J].")
    for l in (1, 2):
        Jmax = 34 if l == 1 else 30
        print("== l = %d, Dphi = %s, tau_0 family J = %d..%d, p from the Fock-space closed form (exact for N <= 8/6)" % (l, DPHI, l, Jmax), flush=True)
        for bname in ('A', 'B'):
            x2, x3, x4, nn, T = bases[bname]
            n12 = [1.0, *nn]
            z1 = [n12]; z3 = [n12, [1, np.cos(2.6), np.sin(2.6)], [1, -1, 0]]
            for eps in ((1e-3, 1e-4, 1e-5, 1e-6) if l == 1 else (1e-3, 1e-4, 1e-5)):
                x1, x2_, x3_, x4_ = config(eps, bases[bname])
                F = F_mft(l, x1, x2_, x3_, x4_, z1, z3)
                S = np.zeros_like(F); Sc = np.zeros_like(F)
                for J in range(l, Jmax + 1):
                    Delta = DPHI + 1 + J
                    p = float(closed_even(l, sp.nsimplify(Delta), J, sp.Rational(7, 4)))
                    W, s = wfull(SpinBlock(l, Delta, J, DPHI, 400), x1, x2_, x3_, x4_, z1, z3)
                    S += p*(-0.5)**J*W; Sc += p*0.5**J*W
                    last = np.max(np.abs(p*W*0.5**J/S))
                u = s['s12'][0]*s['s34'][0]/(s['s13'][0]*s['s24'][0]); v = s['s14'][0]*s['s23'][0]/(s['s13'][0]*s['s24'][0])
                zb = (1 + u - v - np.sqrt((1 + u - v)**2 - 4*u))/2
                r = (F/S)[0]; rc = (F/Sc)[0]
                print("  point %s eps=%.0e u=%.2e zbar=%.3f (last J term %.0e):  ratio = %s   ratio-1 = %s   control = %s" % (
                    bname, eps, u, zb, last, " ".join("%.9f" % x for x in r), " ".join("%.1e" % (x - 1) for x in r),
                    " ".join("%.3f" % x for x in rc)), flush=True)

if __name__ == "__main__":
    import warnings; warnings.filterwarnings('ignore')
    main()

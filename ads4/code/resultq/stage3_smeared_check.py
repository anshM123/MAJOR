"""Result Q, Stage 3 (RESULT_Q_STAGE3.md §2.3): accuracy of the HYBRID smeared parts (flat + leading-Regge curvature) on the binding
near-gap states at Delta_g = 100, against the exact scalar-skeleton remainder (as cert_actionA / Result D's exact-skeleton model):

  B2 += int psi2 fac_h [ex2(nu) - E2 P_J(1-2y) - 2 Pc] dp,   B3 += sigma u int psi3 w3 [(ex3 - 2 Pc)/E3 - (P_J - ...)]  (k = 3 omitted: see below)
  y = u p^2, E2 = (2-y)/(1-y)^2, Pc = P_{-1/2+i nu}(cosh beta) - J0(bh p), fac_{++} = (1-y)^2, fac_{+-} = (1-y)^-2   (cert_actionA conventions)

Only the k = 2 remainder is computed (the k = 3 remainder needs the odd-k continuation [M]; it is reported as not computed).
Usage: OMP_NUM_THREADS=1 python3 stage3_smeared_check.py   -> cache/stage3_smeared_check.txt
"""
import os, sys, io, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
from scipy.special import j0, eval_legendre
import exactmom as E
from stage12_forward_regge import mehler_nodes

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def main():
    Dg = 100; Lam = np.sqrt(98.0 * 99.0)
    FS = json.load(open(os.path.join(CACHE, 'functionals_q.json')))
    names = ['Q0-full-robust']
    states = [(5, 1.2), (6, 1.5), (5, 1.1), (4, np.sqrt(107 / 99.) * 1.0001), (6, 1.2), (8, 1.1), (0, 1.0001)]
    x, w = np.polynomial.legendre.leggauss(48); pn = 0.5 * (x + 1); pw = 0.5 * w
    for J, mfl in states:
        m = mfl * Lam; al = 0.5 * (-(2 * J + 1) + np.sqrt((2 * J + 1) ** 2 + 4 * m * m)); Delta = al + J + 2
        u = 1 / mfl ** 2; bh = (2 * J + 1) * np.sqrt(u); beta = np.arccosh(1 + (2 * J + 1) ** 2 / (2 * m * m))
        hc = E.HeavyC(Delta, J, 3.0, N0=200)
        nus = Lam * pn
        ex2, canc, _ = hc.value(2, nus, diag=True)
        y = u * pn ** 2
        t, W = mehler_nodes(beta, 1600)
        Pl = np.cos(np.outer(nus, t)) @ W
        Pc = Pl - j0(bh * pn)
        E2 = (2 - y) / (1 - y) ** 2
        rem = ex2 - E2 * eval_legendre(J, 1 - 2 * y) - 2 * Pc
        for nm in names:
            F = FS[nm]; c2 = np.array(F['psi2'])
            ps2 = np.polynomial.polynomial.polyval(pn, c2)
            for h, fac in (('++', (1 - y) ** 2), ('+-', (1 - y) ** -2)):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                dB2 = np.sum(pw * ps2 * fac * rem)
                pr(f"J={J:2d} m/Lam={mfl:.4f} {h}: k=2 exact-skeleton remainder int psi2 fac rem = {dB2:+.3e}  (u^2 = {u * u:.3e}; "
                   f"max|rem| = {np.max(np.abs(rem)):.2e} at x = {pn[np.argmax(np.abs(rem))] * np.sqrt(u):.3f}; max cancel {canc.max():.1e})")
    open(os.path.join(CACHE, 'stage3_smeared_check.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()

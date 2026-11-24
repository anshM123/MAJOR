"""
edger/er_sval.py (EDGE_R.md §2.3): validation of the n3 contour beyond nu/Delta = 0.94 on the EXACT scalar.

The scalar crystal integral Y(Delta, J, nu) (kind 'S', pref form with the Q_J prefactor; n3 hybrid contour with the QSC rescaling of
hc_points) equals the exact scalar moment up to factors that depend on nu only, or on (Delta, J) only (n1_scalar_validate: log|Y/(qX)| =
a_nu + b_J).  Hence at fixed J the double ratio
    T(nu) = [Y(D1, nu)/Y(D2, nu)] / [ex2(D1, nu)/ex2(D2, nu)]       (ex2 = the Arb k = 2 moment, agproof/ag_exact, read-only)
must be CONSTANT in nu if Y is computed correctly.  D1 = 1000 + J (twist edge: nu/D1 up to 0.973), D2 = 1100 (nu/D2 <= 0.91, where every
method works).  T is measured at x = nu/m(D1) in [0.88, 0.972]; a jump or drift of T beyond nu/D1 = 0.94 would falsify n3 there.
(The k = 2 moment is used as the exact reference; if Y's nu-dependence differed from ex2's by more than a nu-only factor, T would
drift already below 0.94 -- that is checked too.)
Usage: python3 er_sval.py J x1,x2,...   -> cache/sval_J<J>.txt / .json
"""
import er_guard  # noqa: F401
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
sys.path.insert(0, os.path.join(HERE, '..', 'agproof'))
os.environ.setdefault('L1_DEG', '6')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import n1_pw as PW, z2_points as Z2, fw_frame as FW
import n3_hybrid as NH, n3_flow as F
import ag_exact as AE
from er_diag import mphys

OUT = []
pr = lambda s='': (print(s, flush=True), OUT.append(s))
MQ = int(os.environ.get('ER_SM', '10'))


def ystate(J, D, nu):
    c, wc, L, _ = F.saddle(max(nu, F.NUF), D, J)
    lq = Z2.set_qsc(D, J, c, wc, L)
    r = NH.integrate(nu, D, J, 'S', rho1=3.5, m=MQ, R=8, adaptive=True)
    lY = r['logN'] + lq + np.log(complex(r['Es']))
    reps = np.array([r['logN'] + lq + np.log(complex(v)) for v in r['repsS']])
    return lY, reps, abs(r['errs']/r['Es'])


def exact(J, D, nus):
    al = D - J - 2.0; ab = al + 2*J + 1
    st = dict(J=J, Delta=float(D), m=np.sqrt(al*ab), m2=al*ab, eta=(al*al + ab*ab)/(2*al*ab))
    L3 = AE.build(st, numax=max(nus)*1.02, prec=400, conf=(320, 960, 0.5))
    return np.array([float(L3.value_fast(2, float(v)).mid()) for v in nus])


if __name__ == '__main__':
    J = int(sys.argv[1]); xs = [float(v) for v in sys.argv[2].split(',')]
    FW.patch_qrec_pw(PW, SF); Z2.patch()
    D1, D2 = 1000.0 + J, 1100.0
    m1 = mphys(D1, J); nus = np.array(xs)*m1
    t0 = time.time()
    e1 = exact(J, D1, nus); e2 = exact(J, D2, nus)
    pr(__doc__.split('Usage')[0])
    pr(f"J = {J}: D1 = {D1:g}, D2 = {D2:g}; exact ex2 built [{time.time() - t0:.0f}s]; n3 scalar M = {MQ}")
    rows = []
    for x, nu, a1, a2 in zip(xs, nus, e1, e2):
        t = time.time()
        y1, r1, q1 = ystate(J, D1, nu); y2, r2, q2 = ystate(J, D2, nu)
        lT = (y1 - y2) - np.log(complex(a1/a2))
        lTr = (r1 - r2) - np.log(complex(a1/a2))
        sd = float(np.std(np.exp(lTr - lT).real, ddof=1)/np.sqrt(len(lTr)))
        rows.append(dict(x=x, nud=nu/D1, T=complex(np.exp(lT)), sd=sd, rel1=q1, rel2=q2))
        pr(f"  x = {x:.4f} (nu/D1 = {nu/D1:.4f}): T = {np.exp(lT).real:+.5f}{np.exp(lT).imag:+.5f}i  (rel. replica sd {sd:.1e}; "
           f"rel err Y1 {q1:.1e}, Y2 {q2:.1e}) [{time.time() - t:.0f}s]")
    Ts = np.array([r['T'] for r in rows])
    base = Ts[np.array([r['nud'] < 0.93 for r in rows])]
    if len(base):
        ref = base.mean()
        pr("  T / mean(T over nu/D1 < 0.93): " + ", ".join(f"{r['x']:.3f}: {abs(r['T']/ref):.4f}" for r in rows))
    json.dump(dict(J=J, rows=rows), open(os.path.join(HERE, 'cache', 'sval_J%d.json' % J), 'w'), default=str, indent=1)
    open(os.path.join(HERE, 'cache', 'sval_J%d.txt' % J), 'w').write('\n'.join(OUT) + '\n')

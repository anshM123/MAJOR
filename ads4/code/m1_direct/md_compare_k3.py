"""md_compare_k3.py -- k = 3 piece (C^_3 = Psi_3[o_H q], j' = 1): direct D' images (md_dquadJ with MD_KIND=C3) against the ladder's
k = 3 per-nu polynomial  S3_b(s) = sigma_b sum_k D_k^2 U_{b,k} q_k(s)  (md_k3).  Normalisation N3(nu) from the J = 0 reference state, for
which the ladder image is exactly nu-independent (S3 = q_0 = 1).  Every entry of I3[a, c](nu) must be N3(nu) [p S3_0 + q S3_4].
Control 'LR': U_{b,k} -> r^k (leading-Regge truncation of the same state).
Usage: python3 md_compare_k3.py Dref numax -> out/compare_k3.txt"""
import os, sys, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import md_ladder as ML
import md_k3 as K3

QS = K3.qseq(40)


def S3f(b, Delta, J, lr=False):
    if not lr:
        return sp.lambdify(ML.s, K3.S3(b, Delta, J, QS), 'numpy')
    D_ = sp.nsimplify(Delta); r = sp.Rational(J*(J + 1))/((D_ - 1)*(D_ - 2)); sig = 1 if b == 0 else -1
    return sp.lambdify(ML.s, sp.expand(sig*sum(sp.binomial(k + 3, 3)**2*r**k*QS[k] for k in range(J - b + 1))), 'numpy')


def load(Delta, J):
    out = {}
    for f in glob.glob(os.path.join(HERE, 'out', 'dq3_D%g_J%d_*.npz' % (Delta, J))):
        d = np.load(f); out['-'.join(map(str, d['ns']))] = (d['nus'], d['I']*np.exp(float(d['shift'])))
    return out


def main(Dref=300.0, numax=6.0):
    L = []
    pr = lambda t='': (print(t, flush=True), L.append(t))
    pr(f'# k = 3 test: direct C3H image vs ladder S3_b; reference J = 0 Delta = {Dref:g}; nu <= {numax:g}; residual max|res|/rms')
    R = load(Dref, 0)
    keys = sorted(set((float(os.path.basename(f).split('_')[1][1:]), int(os.path.basename(f).split('_')[2][1:]))
                      for f in glob.glob(os.path.join(HERE, 'out', 'dq3_D*_J*_*.npz'))), key=lambda t: (t[1], t[0]))
    for D, J in keys:
        if (D, J) == (Dref, 0):
            continue
        for g, (nu, I) in sorted(load(D, J).items()):
            if g not in R:
                continue
            nuR, IR = R[g]; IR = IR[0, 0]
            m = np.isin(np.round(nu, 6), np.round(nuR, 6)) & (nu <= numax)
            nu1 = nu[m]; s = nu1*nu1 + 0.25
            I2 = np.array([IR[np.argmin(np.abs(nuR - v))] for v in nu1])
            for lr in (False, True):
                bs = [0] if J < 4 else [0, 4]
                Bm = np.stack([I2*np.broadcast_to(S3f(b, D, J, lr)(s), s.shape) for b in bs], axis=1)
                res = []
                for a in range(I.shape[0]):
                    for c in range(I.shape[1]):
                        A = I[a, c, m]
                        x, *_ = np.linalg.lstsq(Bm, A, rcond=None)
                        res.append(np.max(np.abs(A - Bm @ x))/np.sqrt(np.mean(np.abs(A)**2)))
                pr(f'  Delta={D:g} J={J} grid {g} {"[LR control]" if lr else "[ladder]    "}: ' + ' '.join(f'{v:.1e}' for v in res))
    open(os.path.join(HERE, 'out', 'compare_k3.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 300.0, float(sys.argv[2]) if len(sys.argv) > 2 else 6.0)

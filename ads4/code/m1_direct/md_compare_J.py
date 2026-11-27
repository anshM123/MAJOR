"""md_compare_J.py -- spin J: the direct D' image matrix I[a, c](nu) (md_dquadJ; a, c = conserved parity-even structures of <TTO_Delta>,
<TTO~_{3-Delta}>) against the ladder (fp-H): every entry must be  N(nu) [p_ac S_0(nu) + q_ac S_4(nu)]  (S_4 absent for J < 4), with
S_b(nu; Delta, J) = sum_k D_k^2 U_{b,k} [ z_{k+2}(1/(k+3) - r/(k+1)) + (2 y_{k+2} + V_{b,k} z_{k+2})/T ]  (md_ladder.h2) and N(nu) the
nu-only normalisation, taken from the J = 0 reference state Delta_ref (md_dquad E2, same grid) as N = I_ref/P_ref.
Cross-multiplied least squares (no division by possibly vanishing images): minimise sum_nu |I_ac P_ref - I_ref (p S_0 + q S_4)|^2.
Controls: (LR) S_b -> leading-Regge truncation (T -> inf at fixed r of the same state is not a polynomial, so the control used is the
ladder with all 1/T terms dropped, U_{b,k} -> r^k); (swap) the J-ladder of a different spin.
Usage: python3 md_compare_J.py Dref numax  -> out/compare_J.txt
"""
import os, sys, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import md_ladder as ML
import md_compare_J0 as MC

Z, Y = ML.zy(40)


def S_lad(b, Delta, J, lr=False):
    if lr:
        Delta_s = sp.nsimplify(Delta); T = (Delta_s - 1)*(Delta_s - 2); C = J*(J + 1); r = sp.Rational(1)*C/T
        n = J - b; S = 0
        for k in range(n + 1):
            S += sp.binomial(k + 3, 3)**2*r**k*Z[k + 2]*(sp.Rational(1, k + 3) - r/(k + 1))
        return sp.lambdify(ML.s, sp.expand(S), 'numpy')
    return sp.lambdify(ML.s, ML.h2(b, Delta, J, Z, Y), 'numpy')


def load_J(Delta, J):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, 'out', 'dqJ_D%g_J%d_*.npz' % (Delta, J)))):
        d = np.load(f)
        out.append(('-'.join(map(str, d['ns'])), d['nus'], d['I']*np.exp(float(d['shift']))))
    return out


def analyse(Delta, J, Dref, numax, pr, lr=False):
    RR = {g: (n, I) for g, n, I in MC.best(Dref)}
    bs = [0] if J < 4 else [0, 4]
    fs = [S_lad(b, Delta, J, lr) for b in bs]
    res_all = []
    for g, nu, I in load_J(Delta, J):
        if g not in RR:
            continue
        nuR, IR = RR[g]
        m = np.isin(np.round(nu, 6), np.round(nuR, 6)) & (nu <= numax)
        nu1 = nu[m]; s = nu1*nu1 + 0.25
        I2 = np.array([IR[np.argmin(np.abs(nuR - v))] for v in nu1])
        Pr = MC.P(nu1, Dref)
        Bm = np.stack([I2*f(s) for f in fs], axis=1)
        line = []
        for a in range(I.shape[0]):
            for c in range(I.shape[1]):
                A = I[a, c, m]*Pr
                x, *_ = np.linalg.lstsq(Bm, A, rcond=None)
                r = (A - Bm @ x)/np.sqrt(np.mean(np.abs(A)**2))
                line.append(np.max(np.abs(r)))
                res_all.append((g, a, c, x, r))
        pr(f'  Delta={Delta:g} J={J} grid {g} {"[LR control]" if lr else ""}: max|res|/rms per entry = ' + ' '.join(f'{v:.1e}' for v in line))
    return res_all


def main(Dref=300.0, numax=6.0):
    L = []
    pr = lambda t='': (print(t, flush=True), L.append(t))
    pr(f'# spin-J test of (fp-H) at light Delta: direct D\' image vs ladder span; reference J=0 Delta={Dref:g}; nu <= {numax:g}')
    keys = sorted(set((float(os.path.basename(f).split('_')[1][1:]), int(os.path.basename(f).split('_')[2][1:]))
                      for f in glob.glob(os.path.join(HERE, 'out', 'dqJ_D*_J*_*.npz'))), key=lambda t: (t[1], t[0]))
    for D, J in keys:
        analyse(D, J, Dref, numax, pr)
        analyse(D, J, Dref, numax, pr, lr=True)
    open(os.path.join(HERE, 'out', 'compare_J.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 300.0, float(sys.argv[2]) if len(sys.argv) > 2 else 6.0)

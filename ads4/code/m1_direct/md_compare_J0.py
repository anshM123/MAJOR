"""md_compare_J0.py -- J = 0: direct D' per-nu image of the literal CLS-D k = 2 EE functional (md_dquad, kind E2) against the
ladder formula (fp-H) continued to the state, which at J = 0 is exactly
        P(nu; Delta) = z_2(s)/3 + 2 y_2(s)/T_Delta ,   s = nu^2 + 1/4,   T = (Delta-1)(Delta-2)          [b = 0, n = 0, r = 0, U = 1]
up to a nu-only factor (common to all states) and a Delta-only factor (B |a_{O,0,0}|^2, sin^2, shadow normalisation).
Test (normalisation-free): for states D1, D2 the cross product  I(nu;D1) P(nu;D2) - c I(nu;D2) P(nu;D1)  must vanish for ONE complex c.
c is fitted by weighted least squares on the trusted nu-window; the residual is reported relative to |I(nu;D1) P(nu;D2)|.
Also fitted: the hypothesis 'no Delta-dependence beyond leading Regge' (P -> z_2/3) as a control, which must FAIL at low Delta.
Usage: python3 md_compare_J0.py  -> out/compare_J0.txt
"""
import os, sys, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import md_ladder as ML

z, y = ML.zy(4)
z2 = sp.lambdify(ML.s, z[2], 'numpy'); y2 = sp.lambdify(ML.s, y[2], 'numpy')


def P(nu, Delta, order=1):
    s = nu*nu + 0.25; T = (Delta - 1)*(Delta - 2)
    return z2(s)/3 + (2*y2(s)/T if order >= 1 else 0.0)


def best(Delta, kind='E2'):
    fs = sorted(glob.glob(os.path.join(HERE, 'out', 'dq_D%g_%s_*.npz' % (Delta, kind))),
                key=lambda f: np.prod([int(v) for v in f.split('_')[-1][:-4].split('-')]))
    out = []
    for f in fs:
        d = np.load(f)
        out.append(('-'.join(map(str, d['ns'])), d['nus'], d['I']*np.exp(float(d['shift']))))
    return out


def fit(nu, I1, I2, P1, P2, w):
    # minimise sum w |I1 P2 - c I2 P1|^2 over complex c
    a = I1*P2; b = I2*P1
    c = np.sum(w*np.conj(b)*a)/np.sum(w*np.abs(b)**2)
    res = (a - c*b)/(np.abs(a) + np.abs(c*b))
    return c, res


def main(Dref=100.0, numax=6.0):
    L = []
    pr = lambda t='': (print(t), L.append(t))
    R = best(Dref)
    if not R:
        print('no reference'); return
    gR, nuR, IR = R[-1]
    pr(f'# J = 0 test of (fp-H) continued: direct D\' (literal CLS-D k=2 EE functional) vs ladder P = z2/3 + 2 y2/T; reference Delta = {Dref:g} (grid {gR})')
    pr(f'# residual rho(nu) = (I1 P2 - c I2 P1)/(|I1 P2| + |c I2 P1|), c fitted on nu <= {numax}; "LR" = control with P -> z2/3 (no 1/T term)')
    for f in sorted(set(glob.glob(os.path.join(HERE, 'out', 'dq_D*_E2_*.npz')))):
        pass
    Ds = sorted(set(float(os.path.basename(f).split('_')[1][1:]) for f in glob.glob(os.path.join(HERE, 'out', 'dq_D*_E2_*.npz'))))
    for D in Ds:
        if D == Dref:
            continue
        for (g, nu, I) in best(D)[-2:]:
            m = np.isin(np.round(nu, 6), np.round(nuR, 6))
            nu1 = nu[m]; I1 = I[m]; I2 = np.array([IR[np.argmin(np.abs(nuR - v))] for v in nu1])
            w = (nu1 <= numax).astype(float)
            pr(f'  Delta = {D:g} (T = {(D-1)*(D-2):g}), grid {g}:')
            for lab, order in (('ladder', 1), ('LR', 0)):
                c, res = fit(nu1, I1, I2, P(nu1, D, order), P(nu1, Dref, order), w)
                mx = np.max(np.abs(res[nu1 <= numax]))
                pr(f'     {lab:6s}: c = {c.real:+.6e}{c.imag:+.2e}i   max|rho| (nu<={numax:g}) = {mx:.2e}   rho at nu = ' +
                   ' '.join(f'{v:g}:{abs(r):.1e}' for v, r in zip(nu1, res) if v in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)))
    open(os.path.join(HERE, 'out', 'compare_J0.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 100.0, float(sys.argv[2]) if len(sys.argv) > 2 else 6.0)

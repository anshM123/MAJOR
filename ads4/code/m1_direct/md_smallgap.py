"""md_smallgap.py -- 'small-gap' cross-check: a toy k = 2 functional of reach nu <= 6 (weight w(nu) = nu^2 (36 - nu^2)^2 on [0, 6], i.e. a
functional whose Lambda is ~6, Delta_g ~ 8) applied to light J = 0 states: direct (D') vs ladder, normalised to the Delta = 300 state.
Since both images share the nu-only factor N(nu), the comparison is  int w I_D/N  vs  c_D int w P_D, with N = I_ref/P_ref and c_D the
state normalisation fitted on the per-nu data (md_compare_J0).  Reports the relative difference of the smeared actions."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import md_compare_J0 as MC
g = '12-12-8-12-12-20-12'
RR = {gg: (n, I) for gg, n, I in MC.best(300.0)}
nuR, IR = RR[g]
m = nuR <= 6.0
nu = nuR[m]; w = nu**2*(36 - nu**2)**2
L = ['# toy small-gap functional (reach nu <= 6): relative difference of smeared actions, direct vs ladder (state normalisation from the per-nu fit)']
for D in (10.0, 20.0, 30.0, 100.0):
    for gg, n2, I in MC.best(D):
        if gg != g:
            continue
        I1 = np.array([I[np.argmin(abs(n2 - v))] for v in nu])
        N = IR[m]/MC.P(nu, 300.0)
        c, _ = MC.fit(nu, I1, IR[m], MC.P(nu, D), MC.P(nu, 300.0), np.ones_like(nu))
        dl = np.trapezoid(w*(I1/N).real, nu); la = np.trapezoid(w*(c*MC.P(nu, D)).real, nu)
        lr = np.trapezoid(w*(c*MC.P(nu, D, 0)).real, nu)
        L.append(f'  Delta={D:g}: (direct - ladder)/ladder = {(dl - la)/la:+.2e};  leading-Regge control (no 1/T): {(dl - lr)/lr:+.2e}')
print('\n'.join(L)); open(os.path.join(HERE, 'out', 'smallgap.txt'), 'w').write('\n'.join(L) + '\n')

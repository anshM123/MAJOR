"""
t4: what the available coefficients say about the radius of convergence of the remainder series in y = x^2 = (nu/m)^2
(after REFEREE_UNIFORM_REMAINDER.md issue 2).  Only 4 terms (1/m^2 remainder) / 3 terms (1/m^4) exist (CHMRS order 8), so
this is a ratio-trend (Domb-Sykes) ESTIMATE, not a determination.
Domb-Sykes: c_n ~ A y_c^{-n} n^{g-1}  =>  r_n = c_{n+1}/c_n = y_c^{-1} (1 + (g-1)/n) + O(1/n^2).  With ratios r_1, r_2, r_3, a straight line
through (1/2, r_2), (1/3, r_3) gives y_c^{-1} (intercept) and g; the spread with the line through (1, r_1), (1/2, r_2) is a crude error.
Calibration: the same 4-term procedure applied to a function with known singularity, the bulk-point envelope E_k (y_c = 1).
Output: cache/out_t4.txt
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import factorial
lines = []
def pr(s):
    print(s, flush=True); lines.append(s)

def ds2(c):
    """two ratios r_1 = c_2/c_1, r_2 = c_3/c_2 (indices of y^n, n = 1..3): line through (1, r_1), (1/2, r_2)."""
    c = np.asarray(c, float)
    r1, r2 = c[2]/c[1], c[3]/c[2]
    slope = (r2 - r1)/(0.5 - 1.0)
    inter = r1 - slope
    return r1, r2, (1/inter if inter != 0 else np.inf), (1 + slope/inter if inter != 0 else np.nan)

def ds(c):
    c = np.asarray(c, float)
    r = c[1:]/c[:-1]
    out = {'ratios': r}
    if len(r) >= 3 and np.all(np.isfinite(r)):
        for lab, (i, j) in (('23', (1, 2)), ('12', (0, 1))):
            n1, n2 = i + 1, j + 1
            s = (r[j] - r[i])/(1/n2 - 1/n1)
            inter = r[i] - s/n1
            out[lab] = (1/inter if inter != 0 else np.inf, 1 + s/inter if inter != 0 else np.nan)
    return out

pr("# t4: Domb-Sykes ratio-trend estimates of y_c = x_c^2 (singularity of the remainder series in x^2); 4 coefficients only")
pr("# calibration on E_k(x) = (2-y)(1-y)^{-1-k/2} (true y_c = 1, exponent g = 1+k/2), first 4 nonconstant coefficients:")
for k in (2, 3, 4, 5):
    C = lambda j: 0.0 if j < 0 else np.prod([1 + k/2 + i for i in range(j)])/factorial(j)
    e = [2*C(j) - C(j - 1) for j in range(1, 5)]
    o = ds(e)
    a = ds2([0] + e[:3])
    pr(f"  k={k}: ratios {np.round(o['ratios'], 3)};  two-ratio estimate from y^1..y^3: y_c = {a[2]:.3f} (true 1), g = {a[3]:.2f} (true {1+k/2})")
for DF in (3, 2):
    RC = pickle.load(open(os.path.join(HERE, 'cache', 'residual_coeffs_Df%g.pkl' % DF), 'rb'))
    pr(f"# Df = {DF}: curvature-subtracted 1/m^2 remainder r1~(x, J) = sum_n c_n y^n, n = 0..3 (Domb-Sykes on n = 0..3)")
    for k in (2, 3, 4, 5):
        for J in (0, 2, 4, 8):
            JJ = J*(J + 1.0)
            c = [np.polynomial.polynomial.polyval(JJ, a) for a in RC[k]['r1']]
            o = ds(c)
            a = ds2(c)
            s23 = f"two-ratio (y^1..y^3): r = {a[0]:+.3f}, {a[1]:+.3f} -> y_c = {a[2]:+.3f}, g = {a[3]:+.2f}" 
            pr(f"  k={k} J={J}: c = {np.round(c, 2)}  ratios {np.round(o['ratios'], 3)}  {s23}")
pr("# reading [estimate]: the two-ratio procedure underestimates y_c on the calibration functions (0.72-0.86 for true 1).  For J = 0 the")
pr("#  remainder ratios are sign-regular and give 0.69-0.88 (Df = 3), 0.59-0.66 (Df = 2), i.e. after the calibration bias y_c ~ 1 (Df = 3),")
pr("#  ~ 0.8 (Df = 2): consistent with the bulk-point singularity at nu = m, not determined.  For J >= 2 the ratios change sign (JJ x^2 cross")
pr("#  terms) and nothing can be inferred from 3 ratios.  The control criterion of t3 (last ratio <= 0.8, fallback only for x <= 0.6) is a heuristic cut.")
open(os.path.join(HERE, 'cache', 'out_t4.txt'), 'w').write("\n".join(lines) + "\n")

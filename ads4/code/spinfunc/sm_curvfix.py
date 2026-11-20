"""
sm_curvfix (SMOOTH_BOUND.md §4): in-process fix of the curvature double count of the q3/w2/ee1000/z4 chains (SUBLEADING_REGGE §5).

The rows' Gamma_{k,h} (K2 and direct, C2 and C3H) contain the leading-Regge AdS curvature (sm_curvcheck.py), while every state's base
already contains it (cv2 + cv3; exact-skeleton states subtract only the scalar rho, which excludes it).  Fix ('row', default): subtract
from every row's diagonal Gamma its own curvature image before any model (spline / shape fit / 1/m interpolation / NL term) is built,
   C2 : Gamma_h -= 2 m^2 [P_{-1/2+i nu}(cosh b_row) - J0((2J+1)x)] / env_h2(x)
   C3H: Gamma_h -= m^2 sigma_h [P - J0] / env_h3(x)
(the images of cv2 = sum p2 2(P - J0), cv3 = sum p3 u sigma (P - J0) in dF = sum p2 env2 Gamma/m^2 + p3 u env3 Gamma/m^2), so that each
state carries the curvature once, from its own base, with its own beta.  Off-diagonal keys are untouched (EE leading-Regge off-diagonal 0).
Patch points (SM_CURV_STAGE): 'post' (default) wraps l9_analysis.grid as installed by w2_chain (i.e. AFTER w2's Hampel despike and
direct-row replacement, so the chain keeps exactly its original node set; the z4 non-local rows z4_chain.NLROWS are corrected too);
'pre' wraps w2_chain._A_grid (below the despike, which then acts on Gamma - curv and can drop different nodes -- nonlinear, not used).
Odd rows: q3_analysis.P4.grid (key 'C2'/'C3H' = +-).  Env SM_CURV = row (default) | none | only;  SM_CURVIMG = kappa (default) |
bessel: which image of the leading-Regge curvature is subtracted (see curv_image; 'bessel' = the hybrid's own 2(P - J0) term, which is
out of phase with the rows at x >~ 0.3 and is badly represented by the shape model -- sm_alias.py).
SM_CURV = only (diagnostic): every row's diagonal Gamma is REPLACED by its curvature image (off-diagonals zeroed), so the chain's
dF of a hybrid state is the chain model's integral of the known curvature, to be compared with the state's exact cv2 + cv3
(representation / aliasing error of the chain on a known function with the same oscillation).
Import this module before ee1000_chain / z4_chain.  Nothing on disk is modified.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
import w2_chain as W2

Q = W2.Q
MODE = os.environ.get('SM_CURV', 'row')
LOG = []
_PCACHE = {}


def _pdiff(J, m, xs):
    import remk
    key = (J, round(float(m), 9), tuple(np.round(xs, 12)))
    if key not in _PCACHE:
        b = 2*np.arcsinh((J + 0.5)/m)                 # = arccosh(1 + (2J+1)^2/(2 m^2)) exactly
        P = remk.Pders_vec(xs*m, np.cosh(b), nmax=0)[0]
        _PCACHE[key] = P - j0((2*J + 1)*xs)
    return _PCACHE[key]


IMG = os.environ.get('SM_CURVIMG', 'kappa')


def _jac(n, a, b, x):
    p0 = np.ones_like(x); p1 = 0.5*(a - b + (a + b + 2)*x)
    if n == 0:
        return p0
    for k in range(2, n + 1):
        c1 = 2*k*(k + a + b)*(2*k + a + b - 2)
        c2 = (2*k + a + b - 1)*((2*k + a + b)*(2*k + a + b - 2)*x + a*a - b*b)
        c3 = 2*(k + a - 1)*(k + b - 1)*(2*k + a + b)
        p0, p1 = p1, (c2*p1 - c3*p0)/c1
    return p1


def Dsh(h, J, x):
    """LO shapes D_++ = P_J(1-2x^2), D_+- = ((1+c)/2)^4 P^(0,8)_{J-4}(c) (stable recurrence)."""
    c = 1 - 2*np.asarray(x, float)**2
    return _jac(J, 0, 0, c) if h == '++' else ((1 + c)/2)**4*_jac(J - 4, 0, 8, c)


def curv_image(kind, h, J, m, xs, env, SIG):
    """leading-Regge curvature content of a row's Gamma (c2 units).  SM_CURVIMG = 'bessel': the hybrid's own term,
    2 m^2 (P - J0((2J+1)x))/env2 (k = 2), m^2 s (P - J0)/env3 (k = 3); its J0 carries the flat-Bessel phase 2 l x, NOT the partial-wave
    phase l theta of the rows (SUBLEADING_REGGE §0.2), so at x >~ 0.3 it is out of phase with Gamma.  'kappa' (default): the exact-in-beta
    leading-Regge uniform action of SUBLEADING_REGGE (sr_model 'kappa'), T - env D = env [A D(x_e) - D(x)] + 2 [P - A J0(nu b)]:
       C2: m^2 {A D(x_e) - D(x) + 2 [P - A J0(nu b)]/env2},   C3H: m^2 s {A D(x_e) - D(x) + [P - A J0(nu b)]/env3},
    theta_e = 2 arcsin(x) (b/2)/sinh(b/2), x_e = sin(theta_e/2), A = sqrt(b/sinh b); equal to the bessel image as x -> 0 and, integrated
    against the chain weights, to cv2 + cv3 within < 1e-3 N_s (SUBLEADING_REGGE §0.3)."""
    xs = np.asarray(xs, float)
    d = _pdiff(J, m, xs)
    if IMG == 'bessel':
        return 2*m*m*d/env(h, xs, 2) if kind == 'C2' else m*m*SIG[h]*d/env(h, xs, 3)
    b = 2*np.arcsinh((J + 0.5)/m); A = np.sqrt(b/np.sinh(b))
    te = 2*np.arcsin(np.minimum(xs, 1 - 1e-12))*(b/2)/np.sinh(b/2); xe = np.sin(te/2)
    osc = A*Dsh(h, J, xe) - Dsh(h, J, xs)
    reg = d + j0((2*J + 1)*xs) - A*j0(xs*m*b)            # P - A J0(nu b)
    if kind == 'C2':
        return m*m*(osc + 2*reg/env(h, xs, 2))
    return m*m*SIG[h]*(osc + reg/env(h, xs, 3))


def _fix_even(g, J):
    for D, r in g.items():
        xs, m = r['xs'], r['m']
        for kd in ('C2', 'C3H'):
            for h in ('++', '+-'):
                if (kd, h) in r:
                    c = curv_image(kd, h, J, m, xs, Q.A.env, Q.A.SIG)
                    r[(kd, h)] = (np.asarray(r[(kd, h)], float) - c) if MODE == 'row' else c
                    LOG.append(('even', J, float(D), kd, h, float(np.sqrt(np.mean(c**2)))))
            if MODE == 'only' and (kd, 'off') in r:
                r[(kd, 'off')] = 0*np.asarray(r[(kd, 'off')], float)
    return g


def _fix_odd(g, J):
    for D, r in g.items():
        xs, m = r['xs'], r['m']
        for kd in ('C2', 'C3H'):
            if kd in r:
                c = curv_image(kd, '+-', J, m, xs, Q.P4.env, Q.P4.SIG)
                r[kd] = (np.asarray(r[kd], float) - c) if MODE == 'row' else c
                LOG.append(('odd', J, float(D), kd, '+-', float(np.sqrt(np.mean(c**2)))))
    return g


STAGE = os.environ.get('SM_CURV_STAGE', 'post')


def _fix_nl(g, J):
    """z4's non-local k = 2 rows (built from the un-subtracted Gamma below w2's despike): subtract the image from Re Gamma_c too."""
    Z4 = sys.modules.get('z4_chain')
    if Z4 is None:
        return
    for D in g:
        rec = Z4.NLROWS.get((J, float(D)))
        if rec is None or rec.get('_curvfixed'):
            continue
        m = rec['m']; x = rec['nus']/m
        hm, hp = Z4.hhat(rec['nus'])
        for h in ('++', '+-'):
            if h not in rec:
                continue
            c2 = curv_image('C2', h, J, m, x, Q.A.env, Q.A.SIG)
            c3 = curv_image('C3H', h, J, m, x, Q.A.env, Q.A.SIG)
            rec[(h, 'c')] = (rec[(h, 'c')] - c2) if MODE == 'row' else (c2 + 0j)
            rec[h] = rec[(h, 'c')]*(hm + hp)
            rec[(h, 'C3H')] = (rec[(h, 'C3H')] - c3) if MODE == 'row' else c3
        rec['_curvfixed'] = True
        LOG.append(('nl', J, float(D), 'C2', 'NL', 0.0))


if MODE in ('row', 'only'):
    if STAGE == 'pre':          # below w2's despike (despike then acts on Gamma - curv: a different node set; nonlinear)
        _even_orig = W2._A_grid

        def _A_grid_fixed(J, K, kinds=('C2', 'C3H'), *a, **k):
            return _fix_even(_even_orig(J, K, kinds, *a, **k), J)
        W2._A_grid = _A_grid_fixed
    else:                       # 'post' (default): after w2's despike / direct-row replacement, same node set as the original chain
        _grid_orig = Q.A.grid

        def _grid_fixed(J, K, kinds=('C2', 'C3H'), *a, **k):
            g = _grid_orig(J, K, kinds, *a, **k)
            _fix_even(g, J)
            _fix_nl(g, J)
            return g
        Q.A.grid = _grid_fixed

    _odd_orig = Q.P4.grid

    def _odd_grid_fixed(J, *a, **k):
        return _fix_odd(_odd_orig(J, *a, **k), J)
    Q.P4.grid = _odd_grid_fixed


def summary():
    if not LOG:
        return "sm_curvfix: no rows patched (mode %s)" % MODE
    seen = {}
    for t in LOG:
        seen[t[:5]] = t[5]
    return "sm_curvfix (%s, image %s, stage %s): curvature image subtracted, rms per row/key: " % (MODE, IMG, STAGE) + "; ".join(
        f"{s}J{J}/{D:g}/{kd}{h}:{v:.0f}" for (s, J, D, kd, h), v in sorted(seen.items()))

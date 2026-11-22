"""
ein_run (RESULT_Q_EINSTEIN_CHECK.md): Lemma A (swapped form, D_LR, CLS-D T-shift) on spinning correlators that carry the Einstein-specific
features the free-field checks lack: Witten-diagram (D-function) transcendental structure with log/Li2 at fixed external dimension,
a single-trace SPIN-2 exchange pole (Regge intercept J0 = 2) in all three channels, and contacts of Mellin degree 0, 2, 3, 4.

Test operator  T'(x, z) = sigma(x) T_free(x, z)  in the tensor product of
   * an AdS4 theory of a scalar sigma (Delta = 2) whose O(lambda) 4-pt function G_sigma is a finite sum of D-functions (ein_dbar),
   * the free conformally coupled scalar (og_free; T_free primary, Delta = 3, J = 2).
T' is a primary with Delta = 5, J = 2.  The O(lambda) part of <T'T'T'T'> is G_sigma(X) x <TTTT>_free (every Wick class); it is microcausal
(product of Wightman functions of decoupled theories), conformally covariant and crossing symmetric.  Light-transform weight t^(Delta+J-2)
= t^5, f_k of og_free, gauge of og_free.

    OMP_NUM_THREADS=1 python3 ein_run.py [section ...]   sections: val main ee ctrl end  (default all) -> cache/ein_<sections>.txt
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import og_free as F
import ein_dbar as E

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
L = []
T0 = time.time()
CLS = ('1234', '1243', '1324', 's', 't', 'u')
PI_FREE = {'1234': '1234', '1243': '1324', '1324': '1243', 's': 't', 't': 's', 'u': 'u'}
PI_SIG = {'c0': 'c0', 'c2': 'c2', 'c3': 'c3', 'c4': 'c4', 'xs': 'xt', 'xt': 'xs', 'xu': 'xu'}
PIECES = E.pieces()
JLAB = {'c0': 'contact M=1 (J0=0)', 'c2': 'contact sigma2 (J0=2)', 'c3': 'contact sigma3 (J0=2)', 'c4': 'contact sigma2^2 (J0=4)',
        'xs': 'spin-2 exch. s', 'xt': 'spin-2 exch. t', 'xu': 'spin-2 exch. u'}
SIGMA_LEMMA_A = ('c0', 'c2', 'c3', 'xs', 'xt', 'xu')         # J0 < 4


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def rnull(rng):
    v = rng.normal(size=2); v /= np.linalg.norm(v); s = rng.uniform(0.5, 2.0)
    return np.array([s, s*v[0], s*v[1]])


_SIGCACHE = {}
ROUTE = {'r': 'ellipse', 'kappa': 1.0}     # 'contour': physical route (real spacelike axis to t = 1, then along C), order ROUTE['o']; 'spiral': grid_D
ROUTE['o'] = 'ab'


def sigma_grid(y, n, kind='sv'):
    """G_sigma pieces on the n x n contour grid (depends on dphi only through y)."""
    key = (round(y, 15), n, kind, ROUTE['r'], ROUTE['o'], ROUTE['kappa'])
    if key not in _SIGCACHE:
        t, w = F.contour(n)
        if ROUTE['r'] == 'ellipse':
            D, ms = E.grid_D_ellipse(n, y, 6, kappa=ROUTE['kappa'], kind=kind, corder=ROUTE['o'])
            t, w, _ = E.contour_ellipse(n, ROUTE['kappa'])
        elif ROUTE['r'] == 'contour':
            D, ms = E.grid_D_contour(t, y, 6, kind=kind, corder=ROUTE['o'])
        else:
            D, ms = E.grid_D(t, y, 6, kind=kind)
        T1b, T3b = np.broadcast_arrays(t[:, None], t[None, :])
        Xs = E.Xgauge(T1b, T3b, y)
        _SIGCACHE[key] = ({nm: E.evaluate(fn, Xs, D) for nm, fn in PIECES.items()}, ms)
    return _SIGCACHE[key]


def nodes(n):
    if ROUTE['r'] == 'ellipse':
        t, w, _ = E.contour_ellipse(n, ROUTE['kappa'])
        return t, w
    return F.contour(n)


def free_grid(phi1, phi3, p2, p4, n):
    z1 = F.zvec(phi1); z3 = F.zvec(phi3); y13 = 2*np.sin((phi1 - phi3)/2)**2
    t, w = nodes(n)
    Lk = F.links(t[:, None], t[None, :], z1, z3, p2, p4, y13=y13)
    return F.correlator(Lk, CLS)


def Etable(k, phi1, phi3, p2, p4, n=60, q=0, kind='sv', sig=None):
    """E[piece x class](phi1, phi3) = oint oint t1^(5+q) t3^(5+q) f_k G_sigma[piece] G_free[class]."""
    y = 2*np.sin((phi1 - phi3)/2)**2
    S, ms = sigma_grid(y, n, kind) if sig is None else sig
    G = free_grid(phi1, phi3, p2, p4, n)
    t, w = nodes(n)
    T1b, T3b = np.broadcast_arrays(t[:, None], t[None, :])
    z1 = F.zvec(phi1); z3 = F.zvec(phi3)
    h = T1b**(5 + q)*T3b**(5 + q)*F.fk(k, T1b, T3b, z1, z3, y13=y)*(w[:, None]*w[None, :])
    return {(p, c): complex(np.sum(h*S[p]*G[c])) for p in S for c in CLS}


def lemmaA(A, B, k, pieces=SIGMA_LEMMA_A, wrong=False):
    """A = E(a,b), B = E(b,a).  Returns (residual of the total, max pair residual, scale S = max single |value|)."""
    sg = 1 if k % 2 == 0 else -1
    if wrong:
        sg = -sg
    keys = [(p, c) for p in pieces for c in CLS]
    S = max(max(abs(A[x]) for x in keys), max(abs(B[x]) for x in keys))
    tot = abs(sum(A[x] for x in keys) + sg*sum(B[x] for x in keys))/S
    pair = max(abs(A[(p, c)] + sg*B[(PI_SIG[p], PI_FREE[c])]) for (p, c) in keys)/S
    return tot, pair, S


# ----------------------------------------------------------------------------------------------------------------------- (1)
def sec_val():
    from math import gamma
    pr('(1) VALIDATION of the D-function engine (ein_dbar)')

    def feyn(Dl, X, ng=160):
        x, w = np.polynomial.legendre.leggauss(ng); s = (x + 1)/2; ws = w/2
        u = np.sin(np.pi*s/2)**6; du = 3*np.sin(np.pi*s/2)**4*np.pi*np.sin(np.pi*s/2)*np.cos(np.pi*s/2)*ws
        tot = 0
        for i0 in range(0, ng, 20):
            u1 = u[i0:i0 + 20, None, None]; u2 = u[None, :, None]; u3 = u[None, None, :]
            W = du[i0:i0 + 20, None, None]*du[None, :, None]*du[None, None, :]
            al = [u1, (1 - u1)*u2, (1 - u1)*(1 - u2)*u3, (1 - u1)*(1 - u2)*(1 - u3)]
            J = (1 - u1)**2*(1 - u2); Sig = sum(Dl)/2
            Q = sum(X[(i, j)]*al[i - 1]*al[j - 1] for (i, j) in E.PAIRS)
            fv = al[0]**(Dl[0] - 1)*al[1]**(Dl[1] - 1)*al[2]**(Dl[2] - 1)*al[3]**(Dl[3] - 1)/Q**Sig*J
            tot += np.sum(W*fv)
        return gamma(sum(Dl)/2)/2*tot
    rng = np.random.default_rng(0)
    st = E.principal_state(np.array([1 + 0j]), np.array([1 + 0j]))
    pr(f'  (1a) Phi(z, zb) at u = v = 1: {E.phi_of(st)[0].real:.12f}   (Dbar_1111(1,1) = 2.343907238689..., the known value)')
    pr('  (1b) I_Delta (operator algebra from I_1111 = Phi/(2 X13 X24), raising = -d/dX_ij) vs direct Feynman-parameter integral, Euclidean points:')
    for tr in range(2):
        P = rng.normal(size=(4, 3))
        X = {(i, j): float(np.sum((P[i - 1] - P[j - 1])**2)) for (i, j) in E.PAIRS}
        U = X[(1, 2)]*X[(3, 4)]/(X[(1, 3)]*X[(2, 4)]); V = X[(1, 4)]*X[(2, 3)]/(X[(1, 3)]*X[(2, 4)])
        stt = E.principal_state(np.array([U + 0j]), np.array([V + 0j]))
        D = E.phi_derivs(stt, np.array([U + 0j]), np.array([V + 0j]), 6)
        Xs = {p: np.array([X[p]]) for p in E.PAIRS}
        r1 = E.raise_pair
        tests = {'1111': (E.I1111(), [1, 1, 1, 1]), '1122': (r1(E.I1111(), (3, 4)), [1, 1, 2, 2]), '2222': (E.base2222(), [2, 2, 2, 2]),
                 '3142': (r1(r1(r1(E.I1111(), (1, 3)), (1, 3)), (3, 4)), [3, 1, 4, 2])}
        row = []
        for nm, (f, Dl) in tests.items():
            v = E.evaluate(f, Xs, D)[0]; ref = feyn(Dl, X)
            row.append(f'{nm}: {abs(v - ref)/abs(ref):.1e}')
        B = E.base2222()
        c12 = abs(E.evaluate(E.theta(B, (1, 2)), Xs, D)[0] - E.evaluate(E.theta(B, (3, 4)), Xs, D)[0])/abs(E.evaluate(E.theta(B, (1, 2)), Xs, D)[0])
        # crossing: s-exchange at permuted points (1<->3) = t-exchange
        Xp = {(i, j): X[tuple(sorted(({1: 3, 3: 1}.get(i, i), {1: 3, 3: 1}.get(j, j))))] for (i, j) in E.PAIRS}
        Up = Xp[(1, 2)]*Xp[(3, 4)]/(Xp[(1, 3)]*Xp[(2, 4)]); Vp = Xp[(1, 4)]*Xp[(2, 3)]/(Xp[(1, 3)]*Xp[(2, 4)])
        stp = E.principal_state(np.array([Up + 0j]), np.array([Vp + 0j]))
        Dp = E.phi_derivs(stp, np.array([Up + 0j]), np.array([Vp + 0j]), 6)
        Xps = {p: np.array([Xp[p]]) for p in E.PAIRS}
        cr = []
        for a, b in (('xs', 'xt'), ('xu', 'xu'), ('c3', 'c3'), ('c4', 'c4')):
            va = E.evaluate(PIECES[a], Xps, Dp)[0]; vb = E.evaluate(PIECES[b], Xs, D)[0]
            cr.append(f'{a}(x3,x2,x1,x4)/{b}(x) - 1 = {abs(va/vb - 1):.1e}')
        pr(f'      config {tr}: rel. dev. ' + ', '.join(row) + f';  theta12 = theta34 on I_2222: {c12:.1e};  ' + '; '.join(cr))
    pr('  (1c) analytic continuation: path independence of Phi at contour points (dphi = 1), three paths from the Euclidean base (a, b) = (eps, eps)')
    y = 2*np.sin(0.5)**2
    for (ta, tb) in [(0.5 - 0.5j, 0.2 + 0.4j), (0.01 - 0.1j, 0.9 + 0.3j), (1e-6 + 1e-3j, 0.3 - 0.45j), (0.02 + 0.14j, 0.02 - 0.14j)]:
        a, b = 1/ta, 1/tb
        s1, _ = E.continue_ab(np.array([a]), np.array([b]), y)
        s2 = E.continue_path([(1.5 + 0j, 1e-3 + 0j), (1.5 + 0j, 1.5 + 0j), (a, 1.5 + 0j), (a, b)], y)
        s3 = E.continue_path([(1e-3 + 0j, 1.9 + 0j), (1.2 + 0.3j, 1.9 + 0j), (a, 1.9 + 0j), (a, b)], y)
        v1, v2, v3 = E.phi_of(s1)[0], E.phi_of(s2)[0], E.phi_of(s3)[0]
        pr(f'      t1 = {ta}, t3 = {tb}: Phi = {v1.real:+.12e}{v1.imag:+.12e}i;  path2/path1 - 1 = {abs(v2/v1 - 1):.1e}, path3: {abs(v3/v1 - 1):.1e}')
    pr('  (1d) scalar MFT with the scalar weight t^(Delta-2) (Delta = 2), k = 4: s + t must vanish (Lemma A for MFT [P]); q = 1 control')
    for n in (40, 60):
        t, w = F.contour(n); T1b, T3b = np.broadcast_arrays(t[:, None], t[None, :])
        Xs = E.Xgauge(T1b, T3b, y); U, V = E.UV_of(1/T1b, 1/T3b, y); W = w[:, None]*w[None, :]
        f4 = ((V - U)**3 - (V - U))/(4*U*V)
        for q in (0, 1):
            s_ = np.sum(T1b**q*T3b**q*f4*W*Xs[(1, 2)]**-2*Xs[(3, 4)]**-2); t_ = np.sum(T1b**q*T3b**q*f4*W*Xs[(1, 4)]**-2*Xs[(2, 3)]**-2)
            pr(f'      n = {n}, q = {q}: s = {s_.real:+.12e}, t = {t_.real:+.12e}, |s+t|/|s| = {abs(s_ + t_)/abs(s_):.1e}')
    pr('      (the scalar u-term (13 identity) and the scalar D-function pieces themselves are NOT in the convergence domain at Delta = 2:')
    pr('       |t^0 f4 G| ~ t^-1 x logs at a single endpoint; see RESULT_Q_EINSTEIN_CHECK.md §3.1.  Hence the spinning product below.)')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (2)
def sec_main():
    pr('(2) MAIN: <T\'T\'T\'T\'> = G_sigma[piece] x <TTTT>_free[class], random null (2,4) polarisations, k = 4 (Psi4) and 5 (Psi5)')
    pr('    Lemma A: E(a,b) + (-1)^k E(b,a) = 0.  Pair identity: E[p x c](a,b) + (-1)^k E[pi(p) x pi(c)](b,a) = 0 (pi = 1<->3 image).')
    pr('    S = largest single |E[p x c]| over the J0 < 4 pieces; all residuals relative to S.')
    rng = np.random.default_rng(11)
    for trial in range(3):
        p2, p4 = rnull(rng), rnull(rng); a, b = rng.uniform(-np.pi, np.pi, 2)
        pr(f'  pair {trial}: (phi1, phi3) = ({a:.4f}, {b:.4f}), dphi = {abs(a - b):.4f}')
        for n in ((60, 80) if trial == 0 else (60,)):
            for k in (4, 5):
                A = Etable(k, a, b, p2, p4, n=n); B = Etable(k, b, a, p2, p4, n=n)
                tot, pair, S = lemmaA(A, B, k)
                wtot, _, _ = lemmaA(A, B, k, wrong=True)
                per = {p: lemmaA(A, B, k, pieces=(p, PI_SIG[p]))[0] for p in ('c0', 'c2', 'c3', 'xu')}
                per['xs+xt'] = lemmaA(A, B, k, pieces=('xs', 'xt'))[0]
                c4tot, c4pair, c4S = lemmaA(A, B, k, pieces=('c4',))
                pr(f'    n = {n}, k = {k}: S = {S:.3e}; total residual {tot:.1e}; max pair residual {pair:.1e}; wrong symmetrisation {wtot:.2f}; '
                   + 'per sigma-piece: ' + ', '.join(f'{p} {v:.1e}' for p, v in per.items())
                   + f';  J0=4 contact c4: residual {c4tot:.1e} (scale {c4S:.2e})')
                if n == 60 and k == 4 and trial == 0:
                    big = sorted(((abs(A[x]), x) for x in A if x[0] in SIGMA_LEMMA_A), reverse=True)[:6]
                    pr('      largest single entries E(a,b): ' + ', '.join(f'{x[0]}x{x[1]} = {A[x].real:+.6e}{A[x].imag:+.2e}i' for _, x in big))
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (3)
def sec_ee():
    pr('(3) EE structure (T\'4 . T\'2 contracted), k = 4: E(a,b) = E(b,a), so Lemma A predicts E^EE = 0 at every angle.')
    pr('    per sigma-piece: value summed over free classes; "pair" = xs + xt; ratio = |sum over J0<4 pieces| / max single |piece x class|.')
    for d in (0.05, 0.3, 1.0, 2.0, 3.0):
        r = F.EE_combine(lambda p2, p4: Etable(4, d/2, -d/2, p2, p4, n=60))
        S = max(abs(r[x]) for x in r if x[0] in SIGMA_LEMMA_A)
        per = {p: sum(r[(p, c)] for c in CLS) for p in PIECES}
        tot = sum(per[p] for p in SIGMA_LEMMA_A)
        pr(f'  dphi = {d}: S = {S:.3e};  ' + ', '.join(f'{p} {per[p].real:+.3e}' for p in ('c0', 'c2', 'c3', 'xu'))
           + f', xs {per["xs"].real:+.6e}, xt {per["xt"].real:+.6e};  |total J0<4|/S = {abs(tot)/S:.1e};  c4 (J0=4): {per["c4"].real:+.4e} '
           f'(|c4|/max|c4 x class| = {abs(per["c4"])/max(abs(r[("c4", c)]) for c in CLS):.1e})')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (4)
def sec_ctrl():
    pr('(4) CONTROLS (one generic angle pair, random polarisations, n = 60; residual of the J0<4 total relative to S)')
    rng = np.random.default_rng(5)
    p2, p4 = rnull(rng), rnull(rng); a, b = 0.9, -1.4
    pr(f'    pair (phi1, phi3) = ({a}, {b})')
    for k in (4, 5):
        A = Etable(k, a, b, p2, p4); B = Etable(k, b, a, p2, p4)
        pr(f'    reference k = {k}: residual {lemmaA(A, B, k)[0]:.1e};  wrong symmetrisation: {lemmaA(A, B, k, wrong=True)[0]:.2f}')
    pr('  (4a) light-transform weight t^(5+q) (conformal weight Delta + J - 2 = 5 is q = 0):')
    for q in (-1, 1, 2):
        for k in (4, 5):
            A = Etable(k, a, b, p2, p4, q=q); B = Etable(k, b, a, p2, p4, q=q)
            pr(f'      q = {q:+d}, k = {k}: residual {lemmaA(A, B, k)[0]:.2e}')
    pr('  (4b) non-primary free factor (improvement xi != 1/8):')
    keep = list(F.VERTEX)
    for xi in (0.0, 0.25):
        F.VERTEX[:] = [(1 - 2*xi, (1, 1)), (-2*xi, (0, 2))]
        for k in (4, 5):
            A = Etable(k, a, b, p2, p4); B = Etable(k, b, a, p2, p4)
            pr(f'      xi = {xi}, k = {k}: residual {lemmaA(A, B, k)[0]:.2e}')
    F.VERTEX[:] = keep
    pr('  (4c) NOT single-valued sigma part (Phi without its log term: conformally covariant, same weights and exponents, but its')
    pr('       Euclidean continuation has monodromy, i.e. it is not a microcausal Wightman function):')
    y = 2*np.sin((a - b)/2)**2
    sig = sigma_grid(y, 60, kind='nsv')
    for k in (4, 5):
        A = Etable(k, a, b, p2, p4, sig=sig); B = Etable(k, b, a, p2, p4, sig=sig)
        per = {p: lemmaA(A, B, k, pieces=(p, PI_SIG[p]))[0] for p in ('c0', 'c2', 'xu')}
        pr(f'      k = {k}: total residual {lemmaA(A, B, k)[0]:.2e};  ' + ', '.join(f'{p} {v:.2e}' for p, v in per.items()))
    pr('  (4d) J0 >= k: subtraction order below the Regge intercept of the exchange (k = 2, 3: CLS-D convergence needs k > J0 = 2):')
    for k in (2, 3):
        for n in (40, 60, 80):
            A = Etable(k, a, b, p2, p4, n=n); B = Etable(k, b, a, p2, p4, n=n)
            ex = lemmaA(A, B, k, pieces=('xs', 'xt', 'xu'))
            c0 = lemmaA(A, B, k, pieces=('c0',))
            pr(f'      k = {k}, n = {n}: spin-2 exchange residual {ex[0]:.2e} (scale {ex[2]:.3e});  J0=0 contact residual {c0[0]:.1e}')
    pr(f'  [{time.time() - T0:.0f}s]')


# ----------------------------------------------------------------------------------------------------------------------- (5)
def sec_end():
    pr('(5) endpoint exponents of |t1^5 t3^5 f_k G_sigma[piece] G_free[1243]| as t -> 0 along t = -i tau (dphi = 1)')
    rng = np.random.default_rng(5)
    p2, p4 = rnull(rng), rnull(rng); a, b = 0.5, -0.5; y = 2*np.sin(0.5)**2
    z1 = F.zvec(a); z3 = F.zvec(b)
    taus = np.array([1e-2, 1e-3, 1e-4, 1e-5])
    for lab, (g1, g3) in (('single (t3 = 0.3-0.4i)', (lambda s: -1j*s, lambda s: 0.3 - 0.4j + 0*s)), ('joint t1 = t3', (lambda s: -1j*s, lambda s: -1j*s))):
        t1 = g1(taus); t3 = g3(taus)
        st, _ = E.continue_ab(1/t1, 1/t3, y, nsteps=4000)
        U, V = E.UV_of(1/t1, 1/t3, y)
        D = E.phi_derivs(st, U, V, 6); Xs = E.Xgauge(t1, t3, y)
        Gf = F.correlator(F.links(t1, t3, z1, z3, p2, p4, y13=y), ('1243', '1234'))
        for k in (2, 4):
            f = F.fk(k, t1, t3, z1, z3, y13=y)
            row = []
            for p in PIECES:
                v = np.abs(t1**5*t3**5*f*E.evaluate(PIECES[p], Xs, D)*Gf['1243'])
                row.append(f'{p} {np.polyfit(np.log(taus), np.log(v), 1)[0]:+.2f}')
            pr(f'    {lab}, k = {k}: ' + '  '.join(row))
    pr('    (a single endpoint is integrable iff the exponent > -1; the joint corner iff > -2)')
    pr(f'  [{time.time() - T0:.0f}s]')


SECTIONS = {'val': sec_val, 'main': sec_main, 'ee': sec_ee, 'ctrl': sec_ctrl, 'end': sec_end}

if __name__ == '__main__':
    secs = sys.argv[1:] or list(SECTIONS)
    pr(f'ein_run: sections {secs}')
    for s in secs:
        SECTIONS[s]()
    pr(f'total {time.time() - T0:.0f}s')
    with open(os.path.join(CACHE, 'ein_' + '_'.join(secs) + '.txt'), 'w') as fh:
        fh.write('\n'.join(L) + '\n')

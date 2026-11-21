"""
v4c_mft (RESULT_Q_V4PRIME.md §2): Lemma A in swapped (s-channel, double-commutator) form on MFT <TTTT>, with the SPINNING kernel, in the
complete gauge of v4c_core.  The swapped-form light side of Psi+_{k,nu}[f G_MFT] consists of edge terms only:
    s: the s-channel identity  <T4T3><T1T2>  (higher-order poles of the 2-pt functions x f's poles at x12^2 = x34^2 = 0),
    t: <T4T1><T3T2> = the J-RESUMMED leading family [TT]_{0,J} (all J; f's simple poles),
    u: <T4T2><T3T1> = 0 identically for light-ray polarisations off phi1 = phi3 (collinear contact term: v4c_wu.py);
   [TT]_{n>=1}: double zeros, exactly 0.
Checks:
 (1) exact (sympy, rational celestial points z = (1, (1-m^2)/(1+m^2), 2m/(1+m^2))): k = 2, 3: Res_s = 0;  k = 4: Res_s(a,b) = -Res_t(b,a);
     k = 5: Res_s(a,b) = +Res_t(b,a); for EE and for random (2,4) tensors;  EE closed forms in c = cos(phi1 - phi3) (symbolic m).
 (2) numeric (Cauchy-FFT residues, v4c_core): 300 random (phi1, phi3), all 3^4 tensor components of the (2,4) structure.
 (3) EE light side integrated against the j'=0 celestial kernel A13^{(2+delta)/2} K_nu(dphi) with a corner cutoff |dphi| > dc:
     I_s(dc), I_t(dc), I_s + I_t (the residual), nu = 0, 1, 3.
 (4) corner exponents: |A13^{(2+delta)/2} K R| ~ |dphi|^p as dphi -> 0 (dphi ~ 4 r/rho at the edge rho = 1).
Output: cache/v4c_mft.txt.   Runtime ~ 3 min (sympy part ~ 1 min), 1 core.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v4c_core as V
from v4c_wu import Kfun

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
lines = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)


# ---------------------------------------------------------------- (1) exact
t1s, t3s = sp.symbols('t1 t3')
E = sp.Matrix([1, 0, 0])


def zs(m):
    return sp.Matrix([1, (1 - m**2)/(1 + m**2), 2*m/(1 + m**2)])


def sdot(a, b):
    return -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def sI(X, v):
    return v - 2*X*sdot(X, v)/sdot(X, X)


def sres2(g):
    g = sp.factor(sp.together(g))
    r = sp.residue(g, t1s, sp.Rational(1, 2))
    return sp.factor(sp.residue(sp.factor(sp.together(r)), t3s, sp.Rational(1, 2)))


def sym_edge(M1, M3, fk, S4=None, S2=None):
    z1, z3 = zs(M1), zs(M3)
    x1 = -t1s*z1; x3 = t3s*z3; x2 = -E; x4 = E
    sq = lambda a: sp.factor(sp.expand(sdot(a, a)))
    s12, s34, s13, s14, s23, s24 = [sq(a) for a in (x1 - x2, x3 - x4, x1 - x3, x1 - x4, x2 - x3, x2 - x4)]
    up = s12*s34/(s13*s24); vp = s14*s23/(s13*s24)
    f = {'f2': vp - up, 'f3': sp.Integer(-1), 'f4': ((vp - up)**3 - (vp - up))/(4*up*vp),
         'f5': -((up - vp)**2 - 1)/(4*up*vp)}[fk]
    v = sI(x3 - x4, z3); w = sI(x1 - x2, z1); a = sI(x1 - x4, z1); b = sI(x2 - x3, z3)
    if S4 is None:
        cs, ct = sdot(v, w)**2, sdot(a, b)**2
    else:
        cs = (v.T*S4*v)[0]*(w.T*S2*w)[0]; ct = (a.T*S4*a)[0]*(b.T*S2*b)[0]
    return (sres2(t1s**3*t3s**3*f*cs/(s34**3*s12**3)), sres2(t1s**3*t3s**3*f*ct/(s14**3*s23**3)))


def part1():
    pr('(1) exact edge residues (sympy), R_X := Res_{t1=1/2} Res_{t3=1/2}[t1^3 t3^3 f W_X]  (times (-2 pi i)^2 in Psi+)')
    import random
    random.seed(11)

    def randST():
        A = sp.Matrix(3, 3, lambda i, j: sp.Rational(random.randint(-5, 5), random.randint(1, 4))); A = (A + A.T)/2
        A[1, 1] -= (-A[0, 0] + A[1, 1] + A[2, 2]); return A
    S4, S2 = randST(), randST()
    P = (sp.Rational(1, 3), sp.Rational(2, 5)); Q = (P[1], P[0])
    for fk, sig in (('f2', None), ('f3', None), ('f4', -1), ('f5', +1)):
        for lab, (A4, A2) in (('EE', (None, None)), ('random (2,4) tensor', (S4, S2))):
            Rs1, Rt1 = sym_edge(*P, fk, A4, A2); Rs2, Rt2 = sym_edge(*Q, fk, A4, A2)
            if sig is None:
                pr(f'   {fk} {lab}: R_s(P) = {Rs1}, R_s(Q) = {Rs2}  [identity edge term; [TT]_0 has no f-pole]')
            else:
                pr(f'   {fk} {lab}: R_s(P) = {Rs1}, R_t(Q) = {Rt2};  R_s(P) - ({sig:+d}) R_t(Q) = {sp.simplify(Rs1 - sig*Rt2)};  '
                   f'R_s(Q) - ({sig:+d}) R_t(P) = {sp.simplify(Rs2 - sig*Rt1)};  '
                   f'{"1<->3-SYMMETRIC" if sig == -1 else "1<->3-ANTISYMMETRIC"} part of R_s + R_t = '
                   f'{sp.simplify((Rs1 + Rt1 + Rs2 + Rt2)/2) if sig == -1 else sp.simplify((Rs1 + Rt1 - Rs2 - Rt2)/2)}')
    # EE closed forms, symbolic angle
    m = sp.symbols('m', positive=True); c = sp.symbols('c')
    for fk in ('f4', 'f5'):
        Rs, Rt = sym_edge(m, sp.Integer(0), fk)
        sub = lambda X: sp.factor(sp.simplify(X.subs(m, sp.sqrt((1 - c)/(1 + c)))))
        pr(f'   EE closed form ({fk}), all angles: R_s = {sub(Rs)},  R_t = {sub(Rt)};  R_s + R_t = {sp.simplify(Rs + Rt)},  '
           f'R_s - R_t = {sp.factor(sp.simplify(Rs - Rt))}')


# ---------------------------------------------------------------- (2) numeric, all tensor components
def part2(n=300, seed=5):
    rng = np.random.default_rng(seed)
    a = rng.uniform(-np.pi, np.pi, n); b = rng.uniform(-np.pi, np.pi, n)
    pr(f'(2) numeric (Cauchy-FFT, r = 0.2, M = 32), {n} random celestial pairs, all 81 (2,4) tensor components:')
    for k in (2, 4):
        Sab = V.edge_residue_full(a, b, k, 's'); Tba = V.edge_residue_full(b, a, k, 't')
        Tab = V.edge_residue_full(a, b, k, 't')
        scale = np.max(np.abs(Tab).reshape(n, -1), axis=1) + np.max(np.abs(Sab).reshape(n, -1), axis=1) + 1e-300
        if k == 2:
            pr(f'   k = 2: max over pairs of max|R_s|/max|R_s(k=4)| ~ {np.max(np.abs(Sab))/np.max(np.abs(V.edge_residue_full(a, b, 4, "s"))):.1e}  (identity edge term = 0)')
        else:
            d = np.max(np.abs(Sab + Tba).reshape(n, -1), axis=1)/scale
            pr(f'   k = 4: max_pairs max_comp |R_s(a,b) + R_t(b,a)| / scale = {d.max():.1e} (median {np.median(d):.1e});  '
               f'without the swap |R_s(a,b) + R_t(a,b)|/scale median = {np.median(np.max(np.abs(Sab + Tab).reshape(n, -1), axis=1)/scale):.2f}')
    # k = 5 via the EE-contracted scalar and the tensor (f5 implemented here)
    return a, b


# ---------------------------------------------------------------- (3), (4) integrated EE light side vs corner cutoff
def RtEE(c, k):
    return -4*np.pi**2*((c - 3)*(c - 1)*(c + 1)/524288.0 if k == 4 else (c - 3)*(c - 1)**2*(c + 1)/1048576.0)


def part3():
    pr('(3) EE light side x j\'=0 celestial kernel, corner cutoff |dphi| > dc (k = 4); R_s from the Cauchy-FFT residues (independent of the')
    pr('    closed form), R_t from the closed form; integrand A13^{(2+delta)/2} K_nu(dphi) R(dphi), K_nu = int dphibar [A10 A30]^{-delta/2}:')
    for nu in (0.0, 1.0, 3.0):
        d = 0.5 + 1j*nu
        xs = np.concatenate([np.logspace(-4, np.log10(np.pi*0.9999), 500)])
        Kv = np.array([Kfun(x, nu) for x in xs[::10]])
        Ki = np.interp(np.log(xs), np.log(xs[::10]), Kv.real) + 1j*np.interp(np.log(xs), np.log(xs[::10]), Kv.imag)
        c = np.cos(xs)
        Rs = np.concatenate([V.edge_residue(xs[i:i + 100]/2, -xs[i:i + 100]/2, 4, 's') for i in range(0, len(xs), 100)])
        Rt = RtEE(c, 4)
        Rt_num = V.edge_residue(xs[::50]/2, -xs[::50]/2, 4, 't')
        ker = (2*(1 - c))**((2 + d)/2)*Ki
        Is_d = 2*ker*Rs; It_d = 2*ker*Rt                       # both signs of dphi (even integrands)
        rows = []
        for dc in (1.0, 0.3, 0.1, 0.03, 0.01, 1e-3):
            msk = xs > dc
            Is = np.trapezoid(Is_d[msk], xs[msk]); It = np.trapezoid(It_d[msk], xs[msk])
            rows.append((dc, Is, It))
        pr(f'   nu = {nu}: closed-form R_t vs numeric max rel. dev. {np.max(np.abs(Rt_num - RtEE(np.cos(xs[::50]), 4))/np.abs(RtEE(np.cos(xs[::50]), 4))):.1e}')
        for dc, Is, It in rows:
            pr(f'      dc = {dc:7.0e}:  I_s = {Is.real:+.6e}{Is.imag:+.6e}i   I_t = {It.real:+.6e}{It.imag:+.6e}i   '
               f'|I_s + I_t|/|I_t| = {abs(Is + It)/abs(It):.1e}')
        # (4) corner exponent of each term
        sm = xs < 3e-2
        p_t = np.polyfit(np.log(xs[sm]), np.log(np.abs(It_d[sm])), 1)[0]
        p_s = np.polyfit(np.log(xs[sm & (xs > 3e-3)]), np.log(np.abs(Is_d[sm & (xs > 3e-3)])), 1)[0]
        pr(f'      corner exponents (dphi -> 0): |integrand_t| ~ dphi^{p_t:.3f}, |integrand_s| ~ dphi^{p_s:.3f}  '
           f'(analytic: 5/2 + 2 = 4.5, times log for nu = 0); tail of each term ~ dc^5.5')


if __name__ == '__main__':
    t0 = time.time()
    pr('# v4c_mft: swapped-form Lemma A on MFT <TTTT>, spinning CLS-D kernel, complete gauge (x0 = 0, x2 = -e, x4 = e, z0 fixed)')
    part1()
    pr(f'   [{time.time() - t0:.0f} s]')
    part2()
    part3()
    pr(f'# total {time.time() - t0:.0f} s')
    open(os.path.join(CACHE, 'v4c_mft.txt'), 'w').write('\n'.join(lines) + '\n')

"""
fn_scalar (F_NORM.md, sec. 3.1-3.3): the scalar (F) chain, CHMRS 2106.10274 conventions, d = 3.

 A1  light value of a Penedones-normalised contact under the heavy-normalised C_{k,nu} (exact algebra, own implementation):
       Bhat_{k,T}[M] = (-1)^{k/2} oint_{|S'|=oo} dS'/(2 pi i) M(S', T + 2Df - S') / [(2Df - S') (Df - S'/2)_{k/2} (Df - T'/2)_{k/2}]
                     = (-1)^{k/2} * (sum of the finite residues)          (sympy, exact rational function of T)
       C_{k,nu}[M]   = 2^{d-2-2k} pi^{3/2} Gamma((d-1)/2) a_Df / gamma_a(nu)^2  int dT/(4 pi i) Gamma(T/2)^2 gamma_{a-T}(nu) Bhat_{k,T}
                       (a = 2Df + k - 1; T-integral by Barnes' first lemma on the Newton expansion of Bhat in z = T/2; cross-checked
                        by direct mpmath quadrature on Re T = c0)
     contacts M_{a,b} = c_{2a+3b} (S^2+T^2+U^2)^a (STU)^b,  c_n = 2^{n-3} Gamma(2Df - d/2 + n) / (pi^{d/2} Gamma(Df)^2 Gamma(Df-(d-2)/2)^2)
     (c_n is CHMRS's Penedones normalisation: M_{a,b} <-> M_flat = (s^2+t^2+u^2)^a (stu)^b).  Expected (CHMRS (C2C4Light),
     (BknuContactsLeading)): -C_{4,nu}[M_{2,0}] = 4, -C_{2,nu}[M_{1,0}] = 2, C_{4,nu}[M_{1,0}] = C_{4,nu}[M_{0,1}] = 0, ...
 A3  Penedones top-degree map vs the heavy-density dictionary, on exact AdS data only (no functional):
       a J = 0 heavy exchange with Mellin residues rho_n = -1/2 K_n Q^n (exactmom 'density' convention, for which
       m^8 (-1)^J C_{4,nu}[block]/q -> 2 = the flat B_4 image of one state: CHMRS (bulkpointbknu), exactmom.HeavyC) has the
       degree-4 top coefficient  a_top = coefficient of (S^2+T^2+U^2)^2 in its large-m (EFT) expansion = -1/2 sum_n rho_n / s_n^5.
       The flat amplitude of the same exchange, with the coupling fixed by the heavy density (lambda^2/m^2 = q (-1)^J), is
       lambda^2 sum_ch 1/(m^2 - x), whose (s^2+t^2+u^2)^2 coefficient is g4_flat = q (-1)^J / (2 m^8).
       Penedones => a_top = c_4 g4_flat (1 + O(1/m^2)).  The script prints R3 = a_top / (c_4 g4_flat).
 A2' calibration at the functional level: C_{4,nu}[block] + C_{4,nu}[a_top (S^2+T^2+U^2)^2]  ->  0  (relative), nu = 0, 2, 5.

Usage: (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 fn_scalar.py)  -> cache/fn_scalar.txt   (~2 min)
Read-only imports: ../exactskel/exactmom.py, ../uniform/skel.py (q_{Delta,J}).
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
import mpmath as mp
import sympy as sp

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
OUT = []
T0 = time.time()
mp.mp.dps = 30


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


S_, T_, Sp_ = sp.symbols('S T Sp')


def c_n(n, Df, d=3):
    return mp.mpf(2)**(n - 3)*mp.gamma(2*Df - mp.mpf(d)/2 + n)/(mp.pi**(mp.mpf(d)/2)*mp.gamma(Df)**2*mp.gamma(Df - mp.mpf(d - 2)/2)**2)


def a_Df(Df, d=3):
    return (mp.pi**(mp.mpf(d - 3)/2)*mp.gamma(Df)**2*mp.gamma(Df - mp.mpf(d - 2)/2)**2/(mp.mpf(2)**(d - 5)*mp.gamma(mp.mpf(d - 1)/2)))


def gam(a, nu, d=3):
    a = mp.mpf(a)
    return mp.gamma((1 + a - mp.mpf(d)/2 - 1j*nu)/2)*mp.gamma((1 + a - mp.mpf(d)/2 + 1j*nu)/2)


def Bhat_poly(Mexpr, k, Df):
    """exact Bhat_{k,T}[M] for a polynomial M(S, T) (U = 4Df - S - T): (-1)^{k/2} x sum of finite residues in S'."""
    Dfr = sp.nsimplify(Df)
    Tp = T_ + 2*Dfr - Sp_
    Mv = Mexpr.subs({S_: Sp_, T_: Tp}, simultaneous=True)
    poch = lambda x, n: sp.Mul(*[(x + j) for j in range(n)])
    integrand = Mv/((2*Dfr - Sp_)*poch(Dfr - Sp_/2, k//2)*poch(Dfr - Tp/2, k//2))
    integrand = sp.together(integrand)
    poles = [2*Dfr + 2*j for j in range(k//2)] + [T_ - 2*j for j in range(k//2)]
    tot = 0
    for p0 in poles:
        tot += sp.residue(integrand, Sp_, p0)
    return sp.simplify(sp.Integer(-1)**(k//2)*tot)


def C_from_Bhat_poly(Bpoly, k, Df, nu, d=3):
    """C_{k,nu} for a Bhat that is a polynomial in T: Barnes first lemma termwise in the rising-factorial basis of z = T/2."""
    z = sp.Symbol('z')
    P = sp.Poly(sp.expand(Bpoly.subs(T_, 2*z)), z)
    deg = P.degree()
    # Newton expansion in rising factorials (z)_j: P(z) = sum_j b_j (z)_j  (solve triangular system at z = 0, -1, ..., -deg? use
    # falling evaluation: (z)_j at z = -i vanishes for j > i)
    coeffs = [mp.mpf(0)]*(deg + 1)
    vals = [mp.mpf(str(sp.N(P.eval(-i), 40))) for i in range(deg + 1)]
    for i in range(deg + 1):
        s = vals[i]
        for j in range(i):
            s -= coeffs[j]*mp.rf(-i, j)
        coeffs[i] = s/mp.rf(-i, i)
    a = 2*Df + k - 1
    c = (1 + a - mp.mpf(d)/2 - 1j*nu)/2; dd = (1 + a - mp.mpf(d)/2 + 1j*nu)/2
    I = 0
    for j, b in enumerate(coeffs):
        # int dz/(2 pi i) Gamma(z+j) Gamma(z) Gamma(c-z) Gamma(dd-z) = Gamma(c+j)Gamma(dd+j)Gamma(c)Gamma(dd)/Gamma(c+dd+j)
        I += b*mp.gamma(c + j)*mp.gamma(dd + j)*mp.gamma(c)*mp.gamma(dd)/mp.gamma(c + dd + j)
    pref = mp.mpf(2)**(d - 2 - 2*k)*mp.pi**1.5*mp.gamma(mp.mpf(d - 1)/2)*a_Df(Df, d)/gam(a, nu, d)**2
    return pref*I


def C_from_Bhat_quad(Bfun, k, Df, nu, c0, d=3):
    """direct quadrature on Re T = c0 (check of the Barnes evaluation)."""
    a = 2*Df + k - 1

    def f(y):
        T = c0 + 1j*y
        return (mp.gamma(T/2)**2*mp.gamma((1 + a - T - mp.mpf(d)/2 - 1j*nu)/2)*mp.gamma((1 + a - T - mp.mpf(d)/2 + 1j*nu)/2)*Bfun(T))
    I = mp.quad(f, [-mp.inf, -20, 0, 20, mp.inf])/(4*mp.pi)      # dT/(4 pi i) = i dy/(4 pi i)
    pref = mp.mpf(2)**(d - 2 - 2*k)*mp.pi**1.5*mp.gamma(mp.mpf(d - 1)/2)*a_Df(Df, d)/gam(a, nu, d)**2
    return pref*I


def contact(a_, b_, Df):
    U = 4*sp.nsimplify(Df) - S_ - T_
    return (S_**2 + T_**2 + U**2)**a_*(S_*T_*U)**b_


def secA1():
    pr('A1. Light values of Penedones-normalised contacts M_{a,b} = c_{2a+3b}(S^2+T^2+U^2)^a(STU)^b under CHMRS C_{k,nu} (exact)')
    for Df in (3.0, 2.3):
        for (a_, b_, k) in ((1, 0, 2), (2, 0, 4), (1, 0, 4), (0, 1, 4), (1, 1, 4), (3, 0, 4), (0, 2, 4), (2, 0, 2), (3, 0, 6)):
            n = 2*a_ + 3*b_
            B = Bhat_poly(contact(a_, b_, Df), k, Df)
            Bp = sp.Poly(sp.expand(B), T_) if B != 0 else None
            row = []
            for nu in (0.0, 1.0, 3.0, 7.5):
                v = -c_n(n, Df)*C_from_Bhat_poly(B, k, Df, nu) if B != 0 else mp.mpf(0)
                row.append(f'nu={nu}: {mp.nstr(mp.re(v), 15)}' + (f' (Im {mp.nstr(mp.im(v), 2)})' if abs(mp.im(v)) > 1e-20 else ''))
            # CHMRS closed forms
            nus = sp.Symbol('nu')
            dd = 3
            ref = ''
            if (a_, b_, k) in ((1, 0, 2), (2, 0, 4), (3, 0, 6)):
                ref = f'  [CHMRS: 2^a = {2**a_}]'
            if (a_, b_, k) == (1, 1, 4):
                ref = f'  [CHMRS: 2^a[nu^2 - (2Df-(d-8)/2)(2Df+(d-8)/2)] = 2(nu^2 - {float((2*Df - (dd - 8)/2)*(2*Df + (dd - 8)/2)):.4f})]'
            if (a_, b_, k) == (2, 0, 2):
                ref = (f'  [CHMRS (C2C4Light): 8 nu^4 + 4(d^2-4dDf-8d+8Df+16) nu^2 - (2Df-(d-4)/2)[16(3d-16)Df^2 - 4d(d-4)Df + (d-4)^3 - 192Df^3]'
                       f' at nu=0: {float(-(2*Df - (dd - 4)/2)*(16*(3*dd - 16)*Df**2 - 4*dd*(dd - 4)*Df + (dd - 4)**3 - 192*Df**3)):.6f}]')
            pr(f'  Df = {Df}, M_{{{a_},{b_}}}, k = {k}: Bhat_{{k,T}} poly degree {Bp.degree() if Bp else "-"};  -C = ' + ';  '.join(row) + ref)
        # quadrature cross-check of one value
        B = Bhat_poly(contact(2, 0, Df), 4, Df)
        Bf = sp.lambdify(T_, B, 'mpmath')
        vq = -c_n(4, Df)*C_from_Bhat_quad(Bf, 4, Df, 1.0, c0=1.0)
        pr(f'  Df = {Df}: quadrature check of -C_(4,nu=1)[M_(2,0)] on Re T = 1: {mp.nstr(mp.re(vq), 15)}')
    pr(f'  [{time.time() - T0:.0f}s]')


def secA3():
    import exactmom as E
    import skel
    pr('A3. Penedones top-degree map vs the heavy-density dictionary (exact AdS data, J = 0 exchange, Df = 3, d = 3)')
    pr('    a_top = -1/2 sum_n rho_n/s_n^5 (rho_n = -1/2 K_n Q^n, s_n = Delta + 2n);  g4_flat = q/(2 m^8);  R3 = a_top/(c_4 g4_flat)')
    Df = 3.0
    c4 = float(c_n(4, Df))
    rows = []
    for D in (12.5, 20.5, 40.5, 80.5, 160.5, 320.5, 640.5, 1280.5):
        B = E.BlockMellin(D, 0, Df, N0=600, nq=400)
        lK = B.logK0 + 2*B.g1
        Q = B.G[:, 0]
        sn = B.tau + 2*B.ns
        L = lK.max()
        rho = -0.5*np.exp(lK - L)*Q                        # rho_n e^{-L}
        atop = -0.5*np.sum(B.wn*rho/sn**5)                 # e^{-L} units
        lq = E.log_q(D, 0, Df)
        m2 = (D - 2)*(D - 1)
        g4f = np.exp(lq - L)/(2*m2**4)
        R3 = atop/(c4*g4f)
        # also with m^2 -> Delta^2 (leading-order identification) for comparison
        rows.append((D, R3, atop/(c4*np.exp(lq - L)/(2*D**8))))
        pr(f'    Delta = {D:6.1f}: R3 = {R3:.10f}   (with m^2 -> Delta^2: {rows[-1][2]:.6f})')
    Ds = np.array([r[0] for r in rows]); R = np.array([r[1] for r in rows])
    m2s = (Ds - 2)*(Ds - 1)
    A = np.vstack([np.ones_like(Ds), 1/m2s, 1/m2s**2]).T
    cf = np.linalg.lstsq(A[3:], R[3:], rcond=None)[0]
    pr(f'    fit R3 = r0 + r1/m^2 + r2/m^4 (Delta >= 80): r0 = {cf[0]:.8f}, r1 = {cf[1]:.4f}, r2 = {cf[2]:.2f}')
    pr(f'  [{time.time() - T0:.0f}s]')
    return c4


def secA2():
    """calibration at the functional level: C_{4,nu}[block] (exactmom.HeavyC, density) vs minus the term-wise light value of the exchange's
    EFT (1/s_n) expansion truncated at degree R:  M^(R) = -sum_{j<=R} mu_j (S^j + T^j + U^j),  mu_j = sum_n rho_n / s_n^{j+1};
    TP_R = -sum_j mu_j C_{4,nu}[S^j + T^j + U^j]  (exact Bhat + Barnes).  R = 4 is the top (borderline) contact alone."""
    import exactmom as E
    pr('A2. Calibration identity  C_{4,nu}[J=0 block] + TP_R -> 0   (TP_R: term-wise light value of the degree-<=R EFT expansion; R = 4: top')
    pr('    contact only, whose value is -4 a_top/c_4 by A1).  Entries: (C[block] + TP_R)/C[block], Df = 3')
    Df = 3.0
    c4 = float(c_n(4, Df))
    nus = [0.0001, 2.0, 5.0]
    Rmax = 9
    Cj = {}
    for nu in nus:
        for j in range(Rmax + 1):
            B = Bhat_poly(S_**j + T_**j + (4*sp.nsimplify(Df) - S_ - T_)**j, 4, Df)
            Cj[(nu, j)] = float(mp.re(C_from_Bhat_poly(B, 4, Df, nu))) if B != 0 else 0.0
    pr('    C_{4,nu}[S^j+T^j+U^j], j = 4..' + str(Rmax) + ' at nu = 0: ' + ', '.join(f'{Cj[(nus[0], j)]:.6g}' for j in range(4, Rmax + 1)) +
       '  (j <= 3: ' + ', '.join(f'{Cj[(nus[0], j)]:.1g}' for j in range(4)) + ')')
    for D in (20.5, 40.5, 80.5, 160.5):
        hc = E.HeavyC(D, 0, Df)
        v = hc.value(4, np.array(nus))                      # m^8 (-1)^J C[block_density]/q
        B = E.BlockMellin(D, 0, Df, N0=600, nq=400)
        lK = B.logK0 + 2*B.g1; L = lK.max()
        sn = B.tau + 2*B.ns
        rho = -0.5*np.exp(lK - L)*B.G[:, 0]
        mu = [np.sum(B.wn*rho/sn**(j + 1)) for j in range(Rmax + 1)]
        lq = E.log_q(D, 0, Df); m2 = hc.m2
        Cblock = v*np.exp(lq - L)/m2**4
        line = []
        for i, nu in enumerate(nus):
            parts = []
            for R in (4, 5, 6, 7, 8, 9):
                TP = -sum(mu[j]*Cj[(nu, j)] for j in range(R + 1))
                parts.append(f'R={R}: {(Cblock[i] + TP)/Cblock[i]:+.2e}')
            line.append(f'nu={nu:g} [' + ', '.join(parts) + ']')
        pr(f'    Delta = {D:5.1f} (m^8 C_dens = {", ".join(f"{x:.6f}" for x in v)}): ' + '  '.join(line))
    pr(f'  [{time.time() - T0:.0f}s]')


if __name__ == '__main__':
    pr('# fn_scalar: scalar (F) chain, CHMRS conventions, d = 3; see docstring')
    secA1()
    secA3()
    secA2()
    pr(f'# total {time.time() - T0:.0f} s')
    open(os.path.join(CACHE, 'fn_scalar.txt'), 'w').write('\n'.join(OUT) + '\n')

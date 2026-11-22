"""repro/bounds.py -- independent recomputation of the light-side constants and bound values (Results D and Q), the closed-form
error terms of DG_SCAN, and the ABJM placement.

Own implementation.  Inputs read (data only): code/robust/functionals_noforward.json (NF02), code/resultq/functionals_q3.json
(Q3-full-delta0.2).  Definitions taken from the statements (RESULT_D_STATEMENT §2-3, RESULT_Q_STATEMENT §0-1, DG_SCAN §1.3); the
definition of kappa (k/8, k = -2 mu2_6 + 2 mu3_4 + 4a) and of the H17 heuristic were read from resultq/common_q.py, resultq/lp_g4.py
and spinfunc/h3_H17_model.py (definitions only; nothing imported).

All integrals of polynomials against 1/(p^2 + a^2) are done exactly term by term in mpmath (60 digits), with exact rational
coefficients for NF02 (whose monomial coefficients cancel to ~1e-10 relative).
Usage: python3 bounds.py  -> stdout and cache/bounds.txt
"""
import json, os
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


def mom(c, j):
    """int_0^1 p^j sum c_n p^n dp (exact if c are Fractions)."""
    return sum(cn/(n + j + 1) for n, cn in enumerate(c))


def Lpole(c, Lam):
    """int_0^1 sum c_n p^n /(p^2 + a^2) dp, a = 3/(2 Lam), exact term by term (recursion I_n = 1/(n-1) - a^2 I_{n-2})."""
    a = mp.mpf(3)/(2*mp.mpf(Lam)); a2 = a*a
    I = [mp.atan(1/a)/a, mp.log(1 + 1/a2)/2]
    for n in range(2, len(c)):
        I.append(mp.mpf(1)/(n - 1) - a2*I[n - 2])
    return sum(mp.mpf(cn.numerator)/cn.denominator*I[n] if isinstance(cn, Fr) else mp.mpf(cn)*I[n] for n, cn in enumerate(c))


def lam(Dg):
    return mp.sqrt((mp.mpf(Dg) - 2)*(mp.mpf(Dg) - 1))


def expansion(c, pref):
    """pref * L = A ln(2Lam/3) + Binf + b1/Lam + (b2L ln(2Lam/3) + b2)/Lam^2 + ...  (DG_SCAN §1.1 expansion, re-derived)."""
    c1, c2, c3 = [mp.mpf(Fr(x).numerator)/Fr(x).denominator if isinstance(x, Fr) else mp.mpf(x) for x in c[1:4]]
    tail = lambda shift: sum((mp.mpf(cn.numerator)/cn.denominator if isinstance(cn, Fr) else mp.mpf(cn))/(n - shift)
                             for n, cn in enumerate(c) if n >= 4)
    Rp2 = tail(1)        # int R/p^2 = sum c_n/(n-1)
    Rp4 = tail(3)        # int R/p^4 = sum c_n/(n-3)
    A = pref*c1
    Binf = pref*(c2 + c3/2 + Rp2)
    b1 = pref*(-(mp.pi/2)*c2*mp.mpf(3)/2)
    b2L = pref*mp.mpf(9)/4*(-c3)
    b2 = pref*mp.mpf(9)/4*(c1/2 + c2 - Rp4)
    return A, Binf, b1, b2L, b2


# ------------------------------------------------------------------------------------------------------------------ NF02 (Result D)
nf = json.load(open(os.path.join(ROOT, 'robust', 'functionals_noforward.json')))['NF02']
s2 = [Fr(x) for x in nf['psi2_exact']]; s3 = [Fr(x) for x in nf['psi3_exact']]
# check psi_exact = (1-p)^5 p X(p) from X, Y
def omp_pX(Xc, r):
    poly = [Fr(0)] + [Fr(x) for x in Xc]             # p X(p)
    for _ in range(r):
        poly = [a - b for a, b in zip(poly + [Fr(0)], [Fr(0)] + poly)]
    return poly
assert omp_pX(nf['X'], nf['r'])[:len(s2)] == s2 and omp_pX(nf['Y'], nf['r'])[:len(s3)] == s3
k0 = -2*mom(s2, 6) + 2*mom(s3, 4)
psi2 = [x/k0 for x in s2]; psi3 = [x/k0 for x in s3]
pr('# repro/bounds.py: independent light-side constants (exact rational NF02, mpmath 60 digits)')
pr(f'NF02: k0 = -2 int psi2 p^6 + 2 int psi3 p^4 (unnormalised) = {float(k0):.8f}   [claimed 27.48509817]')
pr(f'      k(0) after normalisation = {float(-2*mom(psi2, 6) + 2*mom(psi3, 4)):.15f}; int psi2 = {float(mom(psi2, 0)):.4f} [claimed 10.386]')
A, Binf, b1, b2L, b2 = expansion(psi2, 16)
pr(f'      A = 16 psi2\'(0) = {float(A):.4f} [3022.34]; B_inf = {float(Binf):.4f} [-5327.21]; B_AdS = B_inf - A ln(3/2) = '
   f'{float(Binf - A*mp.log(1.5)):.3f} [-6552.66]')
pr(f'      b1 = {float(b1):.1f} [19293.4]; b2L = {float(b2L):.0f} [-13549]; b2 = {float(b2):.0f} [-14705]; '
   f'Lam0^D = 1.5 exp(-B_inf/A) = {float(1.5*mp.exp(-Binf/A)):.4f} [8.7415]')
F2m = lambda f: mp.mpf(f.numerator)/f.denominator
mu2_4 = F2m(mom(psi2, 4)); mu3_2 = F2m(mom(psi3, 2)); int2 = F2m(mom(psi2, 0))

# ------------------------------------------------------------------------------------------------------------- Q3-full-delta0.2
Qall = json.load(open(os.path.join(ROOT, 'resultq', 'functionals_q3.json')))
Q = Qall['Q3-full-delta0.2']
q2 = [mp.mpf(x) for x in Q['psi2']]; q3 = [mp.mpf(x) for x in Q['psi3']]
qmom = lambda c, j: sum(cn/(n + j + 1) for n, cn in enumerate(c))
AQ, BQ, b1Q, b2LQ, b2Q = expansion(q2, 2)
Lfin = sum(cn/(n - 1) for n, cn in enumerate(q2) if n >= 2)
kQ = -2*qmom(q2, 6) + 2*qmom(q3, 4) + 4*mp.mpf(Q['a'])
pr()
pr('Q3-full-delta0.2 (stored: A = %.6f, B = %.6f, kappa = %.6f, a = %.6f, c = %.2f)' % (Q['A'], Q['B'], Q['kappa'], Q['a'], Q['c']))
pr(f'      A = 2 psi2\'(0) = {float(AQ):.6f}; B = 2 int (psi2 - c1 p)/p^2 = {float(2*Lfin):.6f} (= B_inf {float(BQ):.6f})')
pr(f'      kappa = (-2 mu2_6 + 2 mu3_4 + 4a)/8 = {float(kQ/8):.6f} [0.1099]; int psi2 = {float(qmom(q2, 0)):.4f} [1.016]')
pr(f'      b1 = {float(b1Q):.4f} [10.314 = -(3 pi/2) c2]; b2L = {float(b2LQ):.2f} [860.3]; b2 = {float(b2Q):.1f} [-2508.6]; '
   f'Lam0^Q = {float(1.5*mp.exp(-BQ/AQ)):.4f} [4.5889]')

# ------------------------------------------------------------------------------------------------------------- tables
pr()
pr('Dg | Lam | 16 L_A | cubic bound 240 sqrt(16L_A)/Lam^4 | |g3|^2 Lam^8 <= | 2L_AdS | B_eff | c1 (D) per unit | eps_c1 (Q) | '
   'H17 band | dk(0) 1st-order | ABJM 2z3/2L')
res = {}
for Dg in (100, 300, 1000, 2000, 5000):
    L = lam(Dg)
    L16 = 16*Lpole(psi2, L)
    bound = 240*mp.sqrt(L16)/L**4
    twoL = 2*Lpole(q2, L)
    Beff = twoL - AQ*mp.log(2*L/3)
    c1rel = 1920*int2/(L**2*mp.sqrt(L16))          # relative shift of the D bound per unit c1 (linearised)
    epsQ = 2*qmom(q2, 0)*bound*L**2                          # |eps_c1| per unit c1 with t4 = D's bound
    dkH17 = 36*(6*abs(mu2_4) + 4*abs(mu3_2))/L**2            # nf5/h3 heuristic
    dk0 = -168*mu2_4/L**2                                    # first-order symbol (DG_SCAN §1.3)
    abjm = 2*mp.zeta(3)/twoL
    res[Dg] = dict(L=L, L16=L16, bound=bound, twoL=twoL, Beff=Beff)
    pr(f'{Dg} | {float(L):.3f} | {float(L16):.2f} | {float(bound):.5e} | {float(L16):.1f} | {float(twoL):.4f} | {float(Beff):.4f} | '
       f'{100*float(c1rel):.4f} % | {float(epsQ):.4f}|c1| | +-{100*float(dkH17)/2:.4f} % | {float(dk0):.4e} | {100*float(abjm):.2f} %')
pr()
pr(f'H17 heuristic: Lam^2 |dk| = 36(6|int psi2 p^4| + 4|int psi3 p^2|) = {float(36*(6*abs(mu2_4) + 4*abs(mu3_2))):.2f} [424.97]; '
   f'bound band = half: {float(18*(6*abs(mu2_4) + 4*abs(mu3_2))):.2f}/Lam^2 [212.5]')
pr(f'EE dk(0) first-order: Lam^2 dk(0) = -168 int psi2 p^4 = {float(-168*mu2_4):.3f} [-68.08]; quoted exact dk(0) Lam^2: '
   f'100: {-7.245e-3*float(lam(100))**2:.2f}, 300: {-7.666e-4*float(lam(300))**2:.2f}, 1000: {-6.830e-5*float(lam(1000))**2:.2f} '
   f'[DG_SCAN: -69.9 / -68.3 / -68.1]')
pr(f'c1 (D) coefficient 1920 int psi2 = {float(1920*int2):.1f} [19941]; eps_c1 (Q) coefficient 2 int psi2^Q = {float(2*qmom(q2, 0)):.4f} [2.032]')

# ------------------------------------------------------------------------------------------------------------- asymptotic forms
pr()
pr('DG_SCAN asymptotic forms vs exact (cubic: relative error; quartic: absolute error):')
for Dg in (300, 1000):
    L = res[Dg]['L']; ex = res[Dg]['bound']; tL = res[Dg]['twoL']
    f1 = 240*mp.sqrt(A*mp.log(L/8.7415))/L**4
    f2 = 240*mp.sqrt(A*mp.log(L/8.7415) + 19293/L)/L**4
    f3 = mp.mpf('1.3194e4')*mp.sqrt(mp.log(L/8.7415) + 6.384/L)/L**4
    f4 = 240*mp.sqrt(A*mp.log(mp.mpf(Dg)/8.7415))*(1 + mp.mpf(6)/Dg)/mp.mpf(Dg)**4
    q1 = 17.2407*mp.log(L/4.5889); q2_ = q1 + 10.314/L
    q3_ = 17.2407*mp.log(2*L/3) - 19.2782 + 10.314/L + (860.3*mp.log(2*L/3) - 2508.6)/L**2
    pr(f'  {Dg}: cubic A ln(L/8.7415): {float(f1/ex - 1):+.2e} [-0.30 % / -0.067 %]; +19293/L: {float(f2/ex - 1):+.2e} [+4.5e-5 / +3.6e-6]; '
       f'1.3194e4 sqrt(ln + 6.384/L): {float(f3/ex - 1):+.2e}; Dg form: {float(f4/ex - 1):+.2e} [-0.25 % / -0.054 %]')
    pr(f'        quartic 17.2407 ln(L/4.5889): {float(q1 - tL):+.4f} [-0.058 / -0.013]; +10.314/L: {float(q2_ - tL):+.4f} [-0.024 / -0.0031]; '
       f'+(860.3 ln - 2508.6)/L^2: {float(q3_ - tL):+.5f} [-0.0009 / -0.00002]; theorem form 17.241 ln(2L/3) + B: '
       f'{float(17.241*mp.log(2*L/3) + (-19.220 if Dg == 300 else -19.265)):.3f}')
    pr(f'        cubic with Lam^4 -> Dg^4: {float((L/Dg)**4 - 1):+.4f} [-2.2 % at 300]')

# ------------------------------------------------------------------------------------------------------------- ABJM
pr()
z3 = mp.zeta(3)
# Virasoro-Shapiro: Gamma(-s)Gamma(-t)Gamma(-u)/(Gamma(1+s)Gamma(1+t)Gamma(1+u)), s+t+u = 0, s in units of alpha'/4 = 1/M^2:
# check numerically that it equals -1/(stu) - 2 zeta(3) + O(s^2)
s, t = mp.mpf('1e-3'), mp.mpf('-0.37e-3'); u = -s - t
VS = mp.gamma(-s)*mp.gamma(-t)*mp.gamma(-u)/(mp.gamma(1 + s)*mp.gamma(1 + t)*mp.gamma(1 + u))
pr(f'ABJM: Virasoro-Shapiro check: VS + 1/(stu) = {float(VS + 1/(s*t*u)):.6f} vs -2 zeta(3) = {float(-2*z3):.6f} '
   f'(so g4/(8piG) = 2 zeta(3)/M^6 with M^2 = 4/alpha\'; zeta(3) alpha\'^3/32 = 2 zeta(3) (alpha\'/4)^3: '
   f'{float(z3/32*64):.6f} = {float(2*z3):.6f} in alpha\' = 4 units)')
pr(f'      ML = 2 sqrt(L^2/alpha\') = 2 (pi sqrt(2 lambda))^(1/2) = {float(2*mp.sqrt(mp.pi*mp.sqrt(2))):.5f} lambda^(1/4); '
   f'2^(5/4) pi^(1/2) = {float(2**mp.mpf(1.25)*mp.sqrt(mp.pi)):.5f} [4.22]')
for Dg in (300, 1000):
    tL = res[Dg]['twoL']
    pr(f'      Dg = {Dg}: 2 zeta(3) = {float(2*z3):.5f}; / 2L_AdS = {float(2*z3/tL):.5f} [{"0.0334" if Dg == 300 else "0.0259"}]; '
       f'needed lambda for Lam = ML: {float((res[Dg]["L"]/(2**mp.mpf(1.25)*mp.sqrt(mp.pi)))**4):.3e}')
pr(f'      8 pi G = 256/c_T from c_T = 32 L^2/(pi G): 8 pi * 32/pi = {8*32}')
open(os.path.join(HERE, 'cache', 'bounds.txt'), 'w').write('\n'.join(out) + '\n')
# c_T normalisation chain (standard formulas, recalled from the literature, not from the project):
#   Einstein gravity: C_T = Gamma(d+2)/((d-1) Gamma(d/2) pi^(d/2)) * L^(d-1)/(8 pi G);  free real scalar C_T = d/((d-1) S_d^2), S_d = 2 pi^(d/2)/Gamma(d/2)
d = 3
CT_grav = mp.gamma(d + 2)/((d - 1)*mp.gamma(mp.mpf(d)/2)*mp.pi**(mp.mpf(d)/2))/(8*mp.pi)      # coefficient of L^2/G
Sd = 2*mp.pi**(mp.mpf(d)/2)/mp.gamma(mp.mpf(d)/2)
CT_free = mp.mpf(d)/((d - 1)*Sd**2)
line = (f'c_T check: C_T(Einstein)/C_T(free real scalar) = {float(CT_grav/CT_free*mp.pi):.6f}/pi x L^2/G  [32/pi: c_T = 32 L^2/(pi G)]')
print(line); open(os.path.join(HERE, 'cache', 'bounds.txt'), 'a').write(line + '\n')

"""
Exact large-frequency expansion of the Ghosh-Mazumdar (arXiv:1912.06589, eq. strf1a) N = infinity thermal helicity
correlator of the regular Chern-Simons-fermion theory, at zero spatial momentum, all lambda.

    G(q3) = < J_{--}(-q) J_{++}(q) >_beta ,   beta = 1,  q3 = Matsubara frequency > 0,
    Re G = 2 g_E (shear channel),  Im G = 2 h_E (Hall channel)          (as in ../sumrule_freefield/cs_fermion_probe.py)

../sumrule_freefield/cs_fermion_probe.py FITTED the large-q3 coefficients and found the Hall constant to be
sin(pi lam)/8 * eps to ~1e-5.  Here the expansion is done ANALYTICALLY (sympy), using only
  * the large-q3 expansion of F_F(q3, mu) (exact: a convergent series in 1/q3 for q3 > 2 mu),
  * closed forms of the three 1d integrals that appear (polylogarithms of z = -exp(-mu + i pi lam)),
  * the gap equation for mu_F,
and the result is checked numerically against direct evaluation of G from the probe (imported, not modified).

Result (proved symbolically below):
    G(q3) = A3 q3^3 + A1 q3 + A0 + O(1/q3),   A2 = 0,
    A3 = sin(pi lam/2) e^{i pi lam/2} / (64 pi lam)              (vacuum; Re A3 = (2/256) sin(pi lam)/(pi lam) ~ c_T)
    A1 = mu_F^2 sin(pi lam/2) e^{i pi lam/2} / (8 pi lam)        (psibar-psi, Delta = 2: Omega^{3-Delta} = Omega)
    A0 = (e^{i pi lam} - 1) eps / 4                              (T exchange)
so, with t4 = -4 cos(pi lam), alpha_T = 4 sin(pi lam) (CDP 1707.03007):
    A0/(2 eps) = -1/8 - (t4 - i alpha_T)/32      EXACTLY, i.e.  h-constant/eps = +alpha_T/32,  g-constant/eps = -1/8 - t4/32.
"""
import os
import sys
sys.dont_write_bytecode = True  # do not leave __pycache__ in the imported (read-only) directories

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'sumrule_freefield'))
import cs_fermion_probe as P  # noqa: E402  (read-only use: G, mu_F, eps_CS)

# ---------------------------------------------------------------- symbolic part
lam, mu, u = sp.symbols('lambda mu u', positive=True)     # u = 1/q3
L2, L3, L4 = sp.symbols('L2 L3 L4', real=True)             # Im Li_n(z), z = -exp(-mu + i pi lam)
I, pi = sp.I, sp.pi
b = pi*lam/2                                               # lambda > 0 branch, sgn(lambda) = +1

# 1d integrals (psi(E) = Im log(1 + e^{-E + i pi lam}) = Im log cosh(E/2 - i b) + b):
Psi0 = -L2                                                 # Int_mu^oo psi dE
Psi2 = -(mu**2*L2 + 2*mu*L3 + 2*L4)                        # Int_mu^oo psi E^2 dE
li = 2*I*(mu*L3 + L4)                                      # Int_mu^oo x (Li2(-e^{-x+i pi lam}) - c.c.) dx
gap = {L2: -(pi/2)*(1 - lam)*mu}                           # gap equation (probe's mu_F), solved for Im Li2(z)

# q3 F_F(q3) = 2b - (4b/pi) atan(2 mu/q3) - (8/(pi q3)) Int psi/(1 + 4E^2/q3^2) dE
qFF = 2*b - (4*b/pi)*sp.atan(2*mu*u) - (8*u/pi)*(Psi0 - 4*u**2*Psi2)   # + O(u^5)
q = 1/u
e = sp.exp(I*qFF)
t1 = I/(128*pi*lam*q)*(1 - e)*(q**2 + 4*mu**2)*(q**2 - 4*mu**2 + 4*I*q*mu)
t2 = I/(96*pi*lam)*(2*mu**2*(3*q - 2*I*mu) - 3*I*q**2*mu)
t3 = mu**3/(12*pi)
t4 = -I/(4*pi**2*lam)*li
Gu3 = sp.series(sp.expand((t1 + t2 + t3 + t4)*u**3), u, 0, 4).removeO()
A = {3 - k: sp.simplify(Gu3.coeff(u, k).subs(gap)) for k in range(4)}
eps = -(1/pi)*(mu**3*(1 - 1/lam)/3 + 2*(mu*L3 + L4)/(pi*lam))   # 1211.4843, as eps_CS in the probe

pred = {3: sp.sin(b)*sp.exp(I*b)/(64*pi*lam),
        2: 0,
        1: mu**2*sp.sin(b)*sp.exp(I*b)/(8*pi*lam),
        0: (sp.exp(2*I*b) - 1)*eps/4}
print('Large-q3 expansion of G (sympy, gap equation used):')
for k in (3, 2, 1, 0):
    d = sp.simplify(sp.expand(sp.expand_complex(A[k] - pred[k])))
    print('  q3^%d coefficient - prediction = %s' % (k, d))
    assert d == 0
t4s, aTs = -4*sp.cos(pi*lam), 4*sp.sin(pi*lam)
chk = sp.simplify(sp.expand_complex(A[0]/(2*eps) - (-sp.Rational(1, 8) - (t4s - I*aTs)/32)))
print('  A0/(2 eps) - [-1/8 - (t4 - i alpha_T)/32] =', chk)
assert chk == 0
print('  => Hall constant  S_H = Im A0/(2 eps) = alpha_T/32 = sin(pi lam)/8   (exact)')
print('  => shear constant Re A0/(2 eps) = -1/8 - t4/32                       (exact; -1/8 is a lam-independent offset)')
print('  => psibar-psi term: arg A1 = pi lam/2 exactly, i.e. Hall/shear ratio b1/a1 = tan(pi lam/2)')

# ---------------------------------------------------------------- numerical checks
mp.mp.dps = 30


def num_checks(lamv):
    lamv = mp.mpf(lamv)
    m = P.mu_F(lamv)
    z = -mp.e**(-m + 1j*mp.pi*lamv)
    L2v, L3v, L4v = (mp.im(mp.polylog(n, z)) for n in (2, 3, 4))
    psi = lambda E: mp.im(mp.log1p(mp.e**(-E + 1j*mp.pi*lamv)))
    r0 = mp.quad(psi, [m, m + 5, m + 40, m + 160]) - (-L2v)
    r2 = mp.quad(lambda E: psi(E)*E**2, [m, m + 5, m + 40, m + 160]) - (-(m**2*L2v + 2*m*L3v + 2*L4v))
    rg = L2v + (mp.pi/2)*(1 - lamv)*m
    epsv = P.eps_CS(lamv)
    eps_sym = -(1/mp.pi)*(m**3*(1 - 1/lamv)/3 + 2*(m*L3v + L4v)/(mp.pi*lamv))
    bb = mp.pi*lamv/2
    A3v = mp.sin(bb)*mp.e**(1j*bb)/(64*mp.pi*lamv)
    A1v = m**2*mp.sin(bb)*mp.e**(1j*bb)/(8*mp.pi*lamv)
    A0v = (mp.e**(2j*bb) - 1)*epsv/4
    print('lambda = %s: mu_F = %s' % (mp.nstr(lamv, 3), mp.nstr(m, 15)))
    print('   closed forms vs quadrature: Psi0 %.1e, Psi2 %.1e, gap residual %.1e, eps %.1e'
          % (abs(r0), abs(r2), abs(rg), abs(epsv - eps_sym)))
    print('   n      q3          |G - A3 q^3 - A1 q - A0| * q3     (should tend to a constant)')
    for n in (25, 50, 100, 200, 400):
        qv = 2*mp.pi*n
        Gv, _ = P.G(qv, lamv)
        rem = Gv - A3v*qv**3 - A1v*qv - A0v
        print('   %-5d  %-10s  %s' % (n, mp.nstr(qv, 7), mp.nstr(abs(rem)*qv, 10)))
    print('   S_H = Im A0/(2 eps) = %s ; sin(pi lam)/8 = %s' % (mp.nstr(mp.im(A0v)/(2*epsv), 20),
                                                              mp.nstr(mp.sin(mp.pi*lamv)/8, 20)))


if __name__ == '__main__':
    for lv in ('0.05', '0.3', '0.7'):
        num_checks(lv)

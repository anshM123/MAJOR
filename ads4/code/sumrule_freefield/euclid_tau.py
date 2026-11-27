"""
Independent free-field check of the CDP d=3 shear sum rule
    S_E = (1/(pi eps)) Int_0^inf (d rho_+ + d rho_-)/w dw  =?  3/8 + t4/32,
done in IMAGINARY TIME (mixed tau-k representation, Euclidean path integral, Wick contraction with explicit
gamma-matrix traces), not with Hamiltonian matrix elements as in ../scope_sumrule/free_fields_sumrule.py.

For a bosonic Hermitian operator O at zero spatial momentum, 0 < tau < beta,
    G(tau) = Int d^2x <O(tau,x) O(0)>_beta = (1/pi) Int_0^inf dw rho(w) cosh(w(tau-beta/2))/sinh(w beta/2)
                                             + (Drude) constant,
and a tau-independent piece C of G(tau) is a zero-frequency spectral weight: rho(w)/w = pi beta C delta(w) on the
full line (the free-theory limit of the hydrodynamic peak eta/(1+w^2 tau_R^2)).  Every term of the one-loop G(tau)
at p = 0 has the form a(k) e^{-2k tau} + b(k) e^{2k tau} + c(k); from it
    pair continuum   rho(w = 2k) = pi a(k) / (2 (1 + n_B(2k)))   [KMS: b = a e^{-2 beta k} checked]
    Drude            C = Int dk c(k)
and the sum rule integral  (2/pi) Int_0^inf (rho - rho_0)/w dw + beta*C  is the zero Matsubara mode of the
vacuum-subtracted Euclidean correlator (checked separately below).  T = beta = 1 throughout.

Normalisations are cross-checked against the literature: vacuum G_0(tau) vs Osborn-Petkou C_T (free real scalar
3/(32 pi^2), free 2-component Dirac 3/(16 pi^2)), and eps from the free energy (no stress tensor used).
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 30
k, th, tau = sp.symbols('k theta tau', positive=True)
n = sp.Symbol('n')            # occupation number (n_B or n_F at energy k), substituted later
Ep, Em = sp.symbols('Ep Em')  # e^{+k tau}, e^{-k tau}
I = sp.I
kx, ky = k*sp.cos(th), k*sp.sin(th)


def ang(expr):
    """angular integral over theta, then measure k dk/(2pi)^2"""
    e = sp.expand(sp.expand_trig(expr))
    return sp.simplify(sp.integrate(e, (th, 0, 2*sp.pi))*k/(4*sp.pi**2))


def decompose(F):
    """F polynomial in Ep, Em with Ep*Em = 1 -> coefficients of e^{-2k tau}, e^{2k tau}, const"""
    P = sp.expand(F)
    P = sp.expand(P.subs(Em, 1/Ep))
    a = P.coeff(Ep, -2); b = P.coeff(Ep, 2); c = P.coeff(Ep, 0)
    rest = sp.simplify(P - a*Ep**-2 - b*Ep**2 - c)
    assert rest == 0, rest
    return sp.simplify(a), sp.simplify(b), sp.simplify(c)


# ----------------------------------------------------------------------------------------------------------
# real scalar.  mixed propagator <phi_k(tau) phi_-k(0)> = D(tau), 0<tau<beta
D = ((1 + n)*Em + n*Ep)/(2*k)
# T_ij = d_i phi d_j phi - (1/2) delta_ij (d phi)^2 - xi (d_i d_j - delta_ij d^2) phi^2, xi = 1/8.
# At p = 0 the xi-term and the delta_ij-term drop from T2 = T_xy and T1 = (T_xx - T_yy)/2 (all-spatial
# derivatives / trace); vertex of phi_k phi_-k: (i k_i)(-i k_j) = k_i k_j.  Finite-p check in finite_p.py.
VS = {'T2': kx*ky, 'T1': (kx**2 - ky**2)/2}

# ----------------------------------------------------------------------------------------------------------
# 2-component Dirac fermion, S_E = Int psibar gamma_mu d_mu psi, gamma = (sigma3, sigma1, sigma2) for (tau, x, y)
g0 = sp.Matrix([[1, 0], [0, -1]]); gx = sp.Matrix([[0, 1], [1, 0]]); gy = sp.Matrix([[0, -I], [I, 0]])
assert gx*gy + gy*gx == sp.zeros(2) and g0*gx + gx*g0 == sp.zeros(2)
DF = ((1 - n)*Em - n*Ep)/(2*k)          # T Sum e^{-i w_n tau}/(w_n^2+k^2), fermionic, 0<tau<beta
dDF = sp.diff(DF.subs({Em: sp.exp(-k*tau), Ep: sp.exp(k*tau)}), tau).subs(
    {sp.exp(-k*tau): Em, sp.exp(k*tau): Ep})
# S(tau,k) = <psi_k(tau) psibar_k(0)> = -gamma_0 d_tau DF - i gamma.k DF   (checked: (g.d)S = delta)
S_tau = -g0*dDF - I*(gx*kx + gy*ky)*DF
# S(-tau) = -S(beta - tau): replace tau -> 1 - tau, i.e. Ep -> e^k Em, Em -> e^{-k} Ep, with e^k = (1-n)/n
ek = (1 - n)/n
S_mtau = -S_tau.subs({Ep: sp.Symbol('X1'), Em: sp.Symbol('X2')}).subs({sp.Symbol('X1'): ek*Em, sp.Symbol('X2'): Ep/ek})
# check of the propagator: (gamma_0 d_tau + i gamma.k) S = 0 for tau != 0, and jump S(0+) - S(0-) = gamma_0
S_fn = S_tau.subs({Em: sp.exp(-k*tau), Ep: sp.exp(k*tau)})
chk = sp.simplify(g0*sp.diff(S_fn, tau) + I*(gx*kx + gy*ky)*S_fn)
nF_ = 1/(sp.exp(k) + 1)
jump = sp.simplify((S_tau - S_mtau).subs({Ep: 1, Em: 1}).subs(n, nF_))
# T_ij = (1/4)[psibar g_i d_j psi - d_j psibar g_i psi + (i<->j)] -> vertex (i/2)(g_i k_j + g_j k_i)
VF = {'T2': (I/2)*(gx*ky + gy*kx), 'T1': (I/2)*(gx*kx - gy*ky)}


def G_scalar(A, B):
    return 2*VS[A]*VS[B]*D**2


def G_fermion(A, B):
    # <psibar V_A psi (tau) psibar V_B psi (0)> = - tr[V_A S(tau) V_B S(-tau)]; overall phase of the Hermitian
    # operator is (+-1 or +-i): fixed below by reflection positivity G_AA(tau) > 0 (and checked by C_T).
    return -(VF[A]*S_tau*VF[B]*S_mtau).trace()


nB = lambda x: 1/(mp.e**x - 1)
nF = lambda x: 1/(mp.e**x + 1)


def analyse(name, Gfun, nsub, stat):
    out = {}
    for A, B in [('T2', 'T2'), ('T1', 'T1'), ('T1', 'T2')]:
        F = ang(Gfun(A, B))
        a, b, c = decompose(F)
        out[(A, B)] = (a, b, c)
    a, b, c = out[('T2', 'T2')]
    # sign fixing by positivity: vacuum coefficient a0 must be > 0
    a0 = sp.simplify(a.subs(n, 0))
    sgn = 1 if a0.subs(k, 1) > 0 else -1
    for key in out:
        out[key] = tuple(sp.simplify(sgn*x) for x in out[key])
    print('=== %s  (overall phase sign fixed by reflection positivity: %+d)' % (name, sgn))
    for key, (a, b, c) in out.items():
        print('  G_%s%s(tau) per k:  a(k) = %s ;  b(k) = %s ;  c(k) = %s' % (key[0], key[1], a, b, c))
    a, b, c = out[('T2', 'T2')]
    a0 = sp.simplify(a.subs(n, 0))
    # KMS check b = a e^{-2k}
    nk = nsub(k)
    kms = sp.simplify((b - a*sp.exp(-2*k)).subs(n, nk))
    print('  KMS b - a e^{-2k} =', kms)
    # T1 = T2 isotropy and cross term
    print('  G_T1T1 - G_T2T2 :', [sp.simplify(x - y) for x, y in zip(out[('T1', 'T1')], out[('T2', 'T2')])],
          '  cross G_T1T2 :', out[('T1', 'T2')])
    # vacuum spectral function rho_0(w) = pi a0(w/2)/2
    w = sp.Symbol('w', positive=True)
    rho0 = sp.simplify(sp.pi*a0.subs(k, w/2)/2)
    # thermal pair spectral function, ratio to vacuum
    Nw = 1/(sp.exp(2*k) - 1)
    rho_ratio = sp.simplify((a.subs(n, nk)/(1 + Nw))/a0)
    print('  rho_0(w) =', rho0, ' ;  rho_T/rho_0 at w = 2k:', sp.simplify(rho_ratio.rewrite(sp.exp)))
    # C_T check: G_0(tau) = Int dk a0 e^{-2k tau}
    G0 = sp.integrate(a0*sp.exp(-2*k*tau), (k, 0, sp.oo))
    # Osborn-Petkou: <T_xy(X) T_xy(0)> = C_T I_{xy,xy}/X^6, I_{xy,xy} = (I_xx I_yy + I_xy^2)/2, I = 1 - 2XX/X^2
    x, y = sp.symbols('x y', real=True)
    X2 = tau**2 + x**2 + y**2
    Ixx = 1 - 2*x**2/X2; Iyy = 1 - 2*y**2/X2; Ixy = -2*x*y/X2
    r, ph = sp.symbols('r phi', positive=True)
    integrand = ((Ixx*Iyy + Ixy**2)/2/X2**3).subs({x: r*sp.cos(ph), y: r*sp.sin(ph)})
    OP = sp.integrate(sp.integrate(sp.expand(sp.simplify(integrand*r)), (ph, 0, 2*sp.pi)), (r, 0, sp.oo))
    CT = sp.simplify(G0/OP)
    print('  vacuum G_0(tau) =', sp.simplify(G0), ' => C_T =', CT, ' = %s x 3/(32 pi^2)' % sp.nsimplify(CT/(sp.Rational(3, 32)/sp.pi**2)))
    # energy density from the free energy (independent of the stress tensor normalisation)
    if stat == 'B':
        P = -mp.quad(lambda q: q*mp.log(1 - mp.e**(-q)), [0, mp.inf])/(2*mp.pi)
    else:
        P = 2*mp.quad(lambda q: q*mp.log(1 + mp.e**(-q)), [0, mp.inf])/(2*mp.pi)   # 2 dof (particle+antiparticle)
    eps = 2*P
    # spectral pair integral (2/pi) Int (rho - rho0)/w, w = 2k
    af = sp.lambdify((k, n), a, 'mpmath'); a0f = sp.lambdify(k, a0, 'mpmath'); cf = sp.lambdify((k, n), c, 'mpmath')
    nn = nB if stat == 'B' else nF
    NB2 = lambda q: 1/(mp.e**(2*q) - 1)
    rho = lambda q: mp.pi*af(q, nn(q))/(2*(1 + NB2(q)))     # at w = 2q
    rho0f = lambda q: mp.pi*a0f(q)/2
    pair = (2/mp.pi)*mp.quad(lambda q: (rho(q) - rho0f(q))/(2*q)*2, [0, 1, mp.inf])   # dw = 2 dq
    drude = mp.quad(lambda q: cf(q, nn(q)), [0, 1, mp.inf])
    # Euclidean zero mode, per k: Int_0^1 G(tau) dtau - Int_R G_0(tau) dtau
    # (combined analytically per k first, to avoid e^{2k} cancellations)
    zmk = sp.simplify((a*(1 - sp.exp(-2*k))/(2*k) + b*(sp.exp(2*k) - 1)/(2*k) + c - a0/k).subs(n, nk))
    print('  zero-mode integrand per k (thermal - vacuum):', sp.simplify(zmk.rewrite(sp.exp)))
    zmf = sp.lambdify(k, zmk, 'mpmath')
    zm = mp.quad(lambda q: zmf(q) if q > mp.mpf('1e-12') else 0, [0, 1, 10, mp.inf])
    print('  eps = %s   (closed form %s)' % (mp.nstr(eps, 20), 'zeta(3)/pi' if stat == 'B' else '3 zeta(3)/(2 pi)'),
          mp.nstr(eps/(mp.zeta(3)/mp.pi), 15))
    print('  pair  (2/pi)Int d rho/w / eps = %s' % mp.nstr(pair/eps, 20))
    print('  Drude beta*C / eps           = %s' % mp.nstr(drude/eps, 20))
    print('  S_E = (pair + Drude)/eps     = %s' % mp.nstr((pair + drude)/eps, 20))
    print('  Euclidean zero mode / eps    = %s  (must equal S_E)' % mp.nstr(zm/eps, 20))
    # cross channel (parity-odd h) totals
    ac, bc, cc = out[('T1', 'T2')]
    print('  T1-T2 cross channel identically zero:', ac == 0 and bc == 0 and cc == 0)
    return (pair + drude)/eps, pair/eps, drude/eps


if __name__ == '__main__':
    print('fermion propagator check (g.d)S at tau != 0:', chk, ';  S(0+) - S(0-) =', jump)
    SB, pB, dB = analyse('free real scalar (improved T; xi-term drops at p = 0)', G_scalar, lambda q: 1/(sp.exp(q) - 1), 'B')
    SF, pF, dF = analyse('free 2-component Dirac fermion', G_fermion, lambda q: 1/(sp.exp(q) + 1), 'F')
    print()
    print('t4 convention (DICTIONARY_T4.md sec.1; Buchel et al. 0911.4257; CDP 1707.03007): scalar +4, Majorana/Dirac -4')
    for nm, S, t4 in [('scalar', SB, 4), ('Dirac', SF, -4)]:
        pred = mp.mpf(3)/8 + mp.mpf(t4)/32
        print('  %-7s S_E = %s   3/8 + t4/32 = %s   ratio = %s' % (nm, mp.nstr(S, 18), mp.nstr(pred, 18), mp.nstr(S/pred, 18)))
    print('  without the Drude (zero-frequency) weight: scalar %s, Dirac %s  (would NOT match)' % (mp.nstr(pB, 12), mp.nstr(pF, 12)))
    print('  with the Drude weight but no pair term     : scalar %s, Dirac %s' % (mp.nstr(dB, 12), mp.nstr(dF, 12)))

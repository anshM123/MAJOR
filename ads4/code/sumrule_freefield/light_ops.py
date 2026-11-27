"""
Light-operator caveats for the free-field check of the CDP shear sum rule (T = 1).

(1) Large-Matsubara behaviour of the vacuum-subtracted Wick-contracted correlator at p = 0:
        dG_E(w_n) = (2/pi) Int_0^inf drho(w) w/(w^2 + w_n^2) dw      (n != 0; the Drude weight sits only at n = 0)
    An operator O of dimension Delta with <O>_beta != 0 in T x T would give dG_E ~ <O> w_n^{3-Delta}.  phi^2
    (Delta = 1) would give a w_n^2 term, T itself (Delta = 3) a constant.  We show w_n^2 dG_E -> const, i.e. the
    leading thermal term is Delta = 5 (no w^2, no w^0 term) for both theories.
(2) <phi^2>_beta of the massless free scalar in 2+1 D is IR divergent (zero Matsubara mode is a 2D massless field):
        <phi^2>_beta = Int d^2k/(2pi)^2 n_B(k)/k  ~  (T/2pi) log(T/m_IR).
    The p = 0 correlators of T1, T2 contain only k^3-weighted integrals (euclid_tau.py) and are IR finite, so they
    cannot depend on <phi^2>.
(3) Free Dirac fermion: <psibar psi>_beta = 0 and <T_A(tau) psibar psi(0)>_beta = 0 at p = 0 (A = 1, 2).
(4) Improvement term at small but nonzero spatial momentum p (scalar, T_xy, p along the diagonal):
        V(q, p-q) = -(1/2)[q_x(p-q)_y + q_y(p-q)_x] + xi p_x p_y,   xi = 1/8,
        G(tau, p) = 2 Int d^2q/(2pi)^2 V^2 D(tau, E_q) D(tau, E_{p-q}),  E = sqrt(q^2 + m^2),
    split as G_cc + 2 xi G_ci + xi^2 G_ii.  G_ii is <phi^2 phi^2> times p^4: IR divergent as m -> 0 at fixed p,
    but -> 0 as p -> 0 at fixed m; G_cc -> G(tau, 0) of euclid_tau.py.
"""
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import quad

mp.mp.dps = 25
nB = lambda x: 1/(mp.e**x - 1)
nF = lambda x: 1/(mp.e**x + 1)

# ---------------------------------------------------------------- (1)
# drho from euclid_tau.py:  scalar  (w^3/512) 2 n_B(w/2);  Dirac  -(w^3/256) 2 n_F(w/2)
drho = {'scalar': lambda w: w**3/512*2*nB(w/2), 'Dirac': lambda w: -w**3/256*2*nF(w/2)}
print('(1) large-Matsubara behaviour of dG_E(w_n) (Wick contracted, p = 0, Drude excluded for n != 0)')
for nm, f in drho.items():
    mom = (2/mp.pi)*mp.quad(lambda w: w*f(w), [0, 10, mp.inf])
    row = []
    for nn in [1, 4, 16, 64, 256]:
        wn = 2*mp.pi*nn
        G = (2/mp.pi)*mp.quad(lambda w: f(w)*w/(w**2 + wn**2), [0, 10, mp.inf])
        row.append('n=%d: dG=%s, wn^2 dG=%s' % (nn, mp.nstr(G, 8), mp.nstr(G*wn**2, 10)))
    print('  %-6s' % nm, '; '.join(row))
    print('         (2/pi) Int w drho = %s  [= w_n^2 dG limit; Delta=5 data].  w_n^0 and w_n^2 coefficients: 0' % mp.nstr(mom, 12))

# ---------------------------------------------------------------- (2)
print('(2) <phi^2>_beta with IR cutoff m (thermal part, massive propagator):')
for m in [1e-1, 1e-2, 1e-3, 1e-4]:
    val = mp.quad(lambda q: q*nB(mp.sqrt(q*q + m*m))/mp.sqrt(q*q + m*m), [0, m, 1, mp.inf])/(2*mp.pi)
    print('   m = %.0e : <phi^2> = %s ;  (1/2pi) log(1/m) = %s' % (m, mp.nstr(val, 8), mp.nstr(mp.log(1/m)/(2*mp.pi), 8)))

# ---------------------------------------------------------------- (3)
k, th, tau, n = sp.symbols('k theta tau n', positive=True)
Ep, Em = sp.symbols('Ep Em')
I = sp.I
g0 = sp.Matrix([[1, 0], [0, -1]]); gx = sp.Matrix([[0, 1], [1, 0]]); gy = sp.Matrix([[0, -I], [I, 0]])
kx, ky = k*sp.cos(th), k*sp.sin(th)
DF = ((1 - n)*Em - n*Ep)/(2*k)
dDF = (-(1 - n)*Em - n*Ep)/2           # d/dtau of DF
S_tau = -g0*dDF - I*(gx*kx + gy*ky)*DF
ek = (1 - n)/n
X1, X2 = sp.symbols('X1 X2')
S_mtau = -S_tau.subs({Ep: X1, Em: X2}).subs({X1: ek*Em, X2: Ep/ek})
VF = {'T2': (I/2)*(gx*ky + gy*kx), 'T1': (I/2)*(gx*kx - gy*ky)}
one = sp.eye(2)
vev = sp.integrate(sp.expand(-(S_tau.subs({Ep: 1, Em: 1})).trace()), (th, 0, 2*sp.pi))
print('(3) Dirac: <psibar psi> integrand (angle-integrated, per k):', sp.simplify(vev))
for A in ['T1', 'T2']:
    c = sp.integrate(sp.expand(-(VF[A]*S_tau*one*S_mtau).trace()), (th, 0, 2*sp.pi))
    print('    <%s(tau) psibar psi(0)> integrand (angle-integrated, per k):' % A, sp.simplify(c))

# ---------------------------------------------------------------- (4)
xi = 1/8.


def Dm(t, E):
    nb = 1/np.expm1(E)
    return ((1 + nb)*np.exp(-E*t) + nb*np.exp(E*t))/(2*E)


def pieces(p, m, t=0.5):
    """G_cc, G_ci, G_ii at tau = t, p along the diagonal.  Polar coordinates (u, phi) centred at the midpoint
    q = p/2 + u (the integrand is symmetric under u -> -u, so phi in [0, pi) and x2); the two IR points q = 0,
    q = p sit at u = p/2, phi = pi/4 (+pi) and are given to the adaptive quadrature as break points."""
    px = py = p/np.sqrt(2)
    out = np.zeros(3)
    for j, which in enumerate(['cc', 'ci', 'ii']):
        def f(phi, u):
            qx, qy = px/2 + u*np.cos(phi), py/2 + u*np.sin(phi)
            Vc = -0.5*(qx*(py - qy) + qy*(px - qx)); Vi = px*py
            V2 = {'cc': Vc*Vc, 'ci': Vc*Vi, 'ii': Vi*Vi}[which]
            E1 = np.sqrt(qx*qx + qy*qy + m*m); E2 = np.sqrt((px - qx)**2 + (py - qy)**2 + m*m)
            return 2*V2*Dm(t, E1)*Dm(t, E2)*u/(4*np.pi**2)
        inner = lambda u: quad(f, 0, np.pi, args=(u,), points=[np.pi/4], limit=400, epsabs=1e-13, epsrel=1e-10)[0]
        tot = 0.0
        for a, b in [(0, p/2 - 20*m), (p/2 - 20*m, p/2), (p/2, p/2 + 20*m), (p/2 + 20*m, 3*p), (3*p, 60.0)]:
            tot += quad(inner, a, b, limit=400, epsabs=1e-13, epsrel=1e-9)[0]
        out[j] = 2*tot
    return out


if __name__ == '__main__':
    G0 = float(sum(mp.quad(lambda q: f(q), [0, 1, mp.inf]) for f in [
        lambda q: q**3*(1 + nB(q))**2/(32*mp.pi)*mp.e**(-q), lambda q: q**3*nB(q)**2/(32*mp.pi)*mp.e**(q),
        lambda q: q**3*nB(q)*(1 + nB(q))/(16*mp.pi)]))
    print('(4) scalar T_xy at tau = 1/2: p = 0 value (euclid_tau.py) G(1/2, 0) = %.8f' % G0)
    print('     p      m        G_cc         2 xi G_ci      xi^2 G_ii')
    for p in [0.4, 0.2, 0.1, 0.05]:
        for m in [1e-2, 1e-3]:
            cc, ci, ii = pieces(p, m)
            print('   %.2f  %.0e  %.8f  %+.3e  %+.3e' % (p, m, cc, 2*xi*ci, xi*xi*ii))

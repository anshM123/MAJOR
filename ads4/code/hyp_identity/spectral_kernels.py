"""
Free-field Kubo shear transport eta_Y(T) on R x H^2 (R=1) from first principles.

Method (see README.md, Sec. 2):
  eta_Y(T) = (pi/(2 N_Y)) * Int dE (-dn/dE) * S(E,E) / (norm)       [particle-hole, omega->0]
  S(E,E') = sum over single-particle shells |<k|O_Y|l>|^2 is evaluated by Schur's lemma on
  the discrete-series irrep of holomorphic quadratic differentials:
      sum |L_f|^2 = lambda(E,E') ||f||^2 ,  lambda = Int d^2w (1-|w|^2)^2 K(0,w)
  where K(z,w) is built from the H^2 spectral projectors (scalar: Legendre P_{-1/2+ip};
  Dirac: cosh^2(s/2) 2F1(1+ip,1-ip;1;-sinh^2(s/2)) in the conformal chi = Omega^{1/2} psi frame).

Results (this script): the shell kernels are, to all printed digits,
  scalar  J(p)       = p^3 tanh(pi p)(1+4p^2)/(48 pi^3)
  Dirac   lambda_F(p) = p coth(pi p)(1+4p^2)/(6 pi^2)
and the direct quadrature gives eta_Y(T0), c_H(T0), dlog eta/dlog T at T0 = 1/2pi.
"""
import mpmath as mp, numpy as np, time
mp.mp.dps = 20
CUTS = [0, 0.5, 1, 2, 3, 4.5, 6, 9, 12, 18, 25, 35, 50]   # geodesic distance s, u = cosh s

# ---------------- conformally coupled real scalar ----------------
def D_scalar(p, s):
    """D_E(u) = Pi_E' + (u-1) Pi_E'' , Pi_E = mu(p) P_{-1/2+ip}(u), mu = p tanh(pi p)/(2pi)."""
    x = -mp.sinh(s/2)**2
    a = mp.mpf(1)/2 - 1j*p; b = mp.mpf(1)/2 + 1j*p
    F1 = mp.hyp2f1(a, b, 1, x); F2 = mp.hyp2f1(a+1, b+1, 2, x)
    mu = p*mp.tanh(mp.pi*p)/(2*mp.pi)
    return mp.re(-mu*(mp.mpf(1)/4 + p*p)*(F1 + x*F2)/(2*(1 - x)))

def J_scalar(p):
    """J(E) = Int_1^oo du D_E(u)^2 ; lambda = 2 pi J ; S(E,E) = 2 pi N_Y J."""
    return mp.quad(lambda s: D_scalar(p, s)**2*mp.sinh(s), CUTS)

# ---------------- Dirac fermion (2-component, 3d) ----------------
def F_fermion(p, u):
    """Radial kernel F_E(u) of the component projector P^(1)(z,w) = F(sigma)/(1 - conj(w) z),
    normalised F(0) = p coth(pi p)/pi; returns F, F_u, F_uu."""
    x = (1 - u)/2; a = 1 + 1j*p; b = 1 - 1j*p
    h0 = mp.hyp2f1(a, b, 1, x)
    h1 = -(a*b/2)*mp.hyp2f1(a+1, b+1, 2, x)
    h2 = (a*b*(a+1)*(b+1)/8)*mp.hyp2f1(a+2, b+2, 3, x)
    c = (u + 1)/2
    rho = p/mp.tanh(mp.pi*p)/mp.pi
    return [rho*mp.re(v) for v in (c*h0, h0/2 + c*h1, h1 + c*h2)]

def K_fermion(p, u):
    """K_F(0,w) for the diagonal shell (A = B = F/2, the positive-energy projector)."""
    F, F1, F2 = F_fermion(p, u)
    return -(u+1)**2/4*((u-1)*(F2*F - F1**2) + F*F1) + F*F/2

def lam_fermion(p):
    """lambda_F = Int d^2w (1-|w|^2)^2 K_F(0,w) = 8 pi Int du K_F/(u+1)^4."""
    return 8*mp.pi*mp.quad(lambda s: K_fermion(p, mp.cosh(s))/(mp.cosh(s)+1)**4*mp.sinh(s), CUTS)

J_closed   = lambda p: p**3*mp.tanh(mp.pi*p)*(1+4*p**2)/(48*mp.pi**3)
lam_closed = lambda p: p/mp.tanh(mp.pi*p)*(1+4*p**2)/(6*mp.pi**2)

if __name__ == "__main__":
    print("== 1. shell kernels vs closed forms ==")
    for p in ['0.25', '0.5', '1', '1.7', '3']:
        p = mp.mpf(p)
        print(f" p={float(p):4.2f}  J/J_closed-1 = {mp.nstr(J_scalar(p)/J_closed(p)-1,3):>9}"
              f"   lamF/lamF_closed-1 = {mp.nstr(lam_fermion(p)/lam_closed(p)-1,3):>9}", flush=True)

    print("== 2. direct quadrature at T0 = 1/2pi (Gauss-Legendre in E on [0,8], N=30) ==")
    T0 = 1/(2*mp.pi)
    gB = lambda x: mp.exp(x)/mp.expm1(x)**2          # T * (-dn_B/dE)
    gF = lambda x: mp.exp(x)/(mp.exp(x)+1)**2        # T * (-dn_F/dE)
    xs, ws = np.polynomial.legendre.leggauss(30)
    E = [mp.mpf(4)*(x+1) for x in xs]; W = [mp.mpf(4)*w for w in ws]
    CT = {'scalar': mp.mpf(3)/(32*mp.pi**2), 'dirac': mp.mpf(3)/(16*mp.pi**2)}
    for name in ['scalar', 'dirac']:
        t = time.time()
        if name == 'scalar':
            ker = [mp.pi**2*J_scalar(e)/e**2 for e in E]         # eta = pi^2 Int (-n') J/E^2
            g = gB; rho = lambda e: e*mp.tanh(mp.pi*e)/(2*mp.pi)
        else:
            ker = [mp.pi/4*lam_fermion(e) for e in E]           # eta = (pi/4) Int (-f') lambda_F
            g = gF; rho = lambda e: 2*e/mp.tanh(mp.pi*e)/(2*mp.pi)  # particles + antiparticles
        dg = lambda x: mp.diff(g, x)
        eta = sum(w*k*g(e/T0)/T0 for e, w, k in zip(E, W, ker))
        Tdeta = sum(w*k*(-(g(e/T0) + (e/T0)*dg(e/T0))/T0) for e, w, k in zip(E, W, ker))
        cH = sum(w*rho(e)*e*(e/T0**2)*g(e/T0) for e, w in zip(E, W))
        print(f" {name}: eta_Y(T0)/(pi^2 C_T/48) = {mp.nstr(eta/(mp.pi**2*CT[name]/48),12)}"
              f"  c_H(T0)/(pi^3 C_T/12) = {mp.nstr(cH/(mp.pi**3*CT[name]/12),10)}"
              f"  c_H/(4 pi eta) = {mp.nstr(cH/(4*mp.pi*eta),10)}")
        print(f"   dlog eta/dlog T = {mp.nstr(Tdeta/eta,12)}  ->  V(T0) = {mp.nstr(Tdeta/eta-1,12)}"
              f"   ({time.time()-t:.0f}s)", flush=True)

"""v2 (H7_I14_LIGHTSIDE.md §12): D3 -- heavy (Laplace / crystal) versus light (Fourier / shock) plane-wave bookkeeping, phases included.

Frame (3.23) of CHMRS, d = 3, mostly-plus metric on R^{2,1}; future cone V+; H2 = unit future hyperboloid.  A leading-Regge kernel of a
dispersive functional is homogeneous, K(x, y) = (|x||y|)^{a-3} Kh(xh, yh; z's), a = 2 Delta_e + k - 1 (= 7 for <TTTT>, k = 2).
Heavy (Lemma V, REGGE_MOMENT_TTTT §3-4): the block is a plane-wave pair e^{2p.x + 2pb.y} x (constant transverse helicity tensor).
Light (CHP arXiv:1707.07689 eq. AtoB): the Regge term is int dp dpb e^{-2i p.x - 2i pb.y} B(p, pb) x (transverse tensor).

Lemma D3 (proved in §12; checked here):
 (i)  for every kernel homogeneous of degree a-3 on V+ and every constant polarisation tensor, the transform at lambda*p, Re lambda >= 0,
      is lambda^{-a} times the Laplace transform at p (radial integral Gamma(a)(lambda c)^{-a}, same angular integral).  Fourier = lambda = -i:
      phase (-i)^{-a} = e^{i pi a/2} per pair, e^{i pi a} for the pair of pairs (= -1 at a = 7).  No nu-dependence.
 (ii) Parseval on the cone: the angular (r^{-1}-coefficient) pairing of the kernel with the light Regge term equals
      e^{i pi a} int_{H2} dpb <Bh(p0, pb), L(p0, pb)>, L = the heavy per-plane-wave action.  So the light side carries exactly the same
      nu-only factor (gamma_a^2 M^2) as the heavy normalisation.

Checks:
 (A) [mpmath, 30 digits] B.2 angular identity Gamma(a) h_a(nu) = 2^{a-1} sqrt(pi) gamma_a(nu), h_a = 2 pi int_1^oo x^{-a} P_l(x) dx.
 (B) [numpy, ~1e-12] spin-2 (stress-tensor) rows of Lemma L (D1^2, D1D2, D2^2; preserving and flip), not integrated numerically before,
     at a = 7, several nu and xi:  Gamma(a) A_s(xi) vs C_a gamma_a M_s P_l(cosh xi).
 (C) [numpy] direct 3-d transform with complex lambda = e^{-i theta} (radial integral numerical): I(lambda)/I(1) = lambda^{-a} for the tensor
     structures, a = 7 and a = 6.5; theta up to 0.49 pi (Fourier = boundary value theta -> pi/2).
 (D) [numpy + mpmath] two-pair Parseval with the Einstein shock profile: position-space Regge coefficient G(t) built from
     B = Q_1(cosh L) (H2 propagator, Delta_0 = 2) by two explicit H2 convolutions, then projected on P_l: equals h_a(nu)^2 * 2 pi/(nu^2+9/4);
     the heavy action of the same kernel is h_a(nu)^2 P_l(cosh beta).  Ratio light/heavy = 2 pi/(nu^2+9/4) * (nu-independent phase).
 (E) [mpmath] the CHMRS magic-block coefficient is the heavy image of eta-multiplication: the Laplace image of the kernel eta P_l(eta) is
     (1/2)[g P_{l-1} + g~ P_{l+1}] h_7(nu)^2, g = ((2-d)/2 + i nu)/(i nu) (gamma_7(nu+i)/gamma_7(nu))^2 (l2_scalar_validate.gfac).
Output: cache/v2_d3_planewave.txt
"""
import os, sys, io, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
import numpy as np
import mpmath as mp
from numpy.polynomial.legendre import leggauss

buf = io.StringIO()
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")

# ------------------------------------------------------------------ conical functions P_{-1/2+i nu}(cosh d), complex nu allowed
_S, _WS = leggauss(400)
_S = 0.5*(_S + 1); _WS = 0.5*_WS
def Pcon(nu, d):
    """Mehler-Dirichlet P_{-1/2+i nu}(cosh d) = (sqrt2/pi) int_0^d cos(nu t) dt / sqrt(cosh d - cosh t), t = d(1-s^2); vectorised in d.
    Valid for |Im nu| < 1/2 (absolute convergence); for Im nu = +-1 use Pcon_mp."""
    d = np.asarray(d, float)[..., None]
    s = _S; t = d*(1 - s*s)
    den = 2*np.sinh((d + t)/2)*np.sinh(d*s*s/2)            # cosh d - cosh t, stable
    f = np.cos(nu*t)*2*d*s/np.sqrt(den)
    out = np.sqrt(2)/np.pi*(f*_WS).sum(-1)
    return np.where(d[..., 0] < 1e-12, 1.0, out)

def Pcon_mp(nu, x):
    return mp.legenp(-mp.mpf(1)/2 + 1j*nu, 0, x, type=3)

def gam(a, nu):        # CHMRS gamma_a(nu), d = 3
    return mp.gamma((a - mp.mpf(1)/2 - 1j*nu)/2)*mp.gamma((a - mp.mpf(1)/2 + 1j*nu)/2)

# ------------------------------------------------------------------ (A)
mp.mp.dps = 30
pr("(A) B.2 angular identity  Gamma(a) h_a(nu) / (2^{a-1} sqrt(pi) gamma_a(nu)) - 1   [mpmath, dps 30]")
def h_mp(a, nu):
    f = lambda x: x**(-a)*Pcon_mp(nu, x)
    return 2*mp.pi*mp.quad(f, [1, 1.5, 3, 10, 100, mp.inf])
HA = {}
for a in (mp.mpf(7), mp.mpf(13)/2, mp.mpf(9)):
    row = []
    for nu in (0, 0.7, 1.5, 3.0, 6.0):
        h = h_mp(a, nu); HA[(float(a), nu)] = h
        r = mp.gamma(a)*h/(2**(a - 1)*mp.sqrt(mp.pi)*gam(a, nu)) - 1
        row.append(f"nu={nu}: {mp.nstr(abs(r), 3)}")
    pr(f"   a = {mp.nstr(a, 3)}: " + "; ".join(row))

# ------------------------------------------------------------------ geometry for (B), (C)
def frame(xi):
    ph = np.array([np.cosh(xi), np.sinh(xi), 0.0])           # p-hat;  kernel centre x-hat = (1,0,0)
    npar = np.array([np.sinh(xi), np.cosh(xi), 0.0])         # tangent at p-hat along the geodesic from x-hat (parallel transported)
    eperp = np.array([0.0, 0.0, 1.0])
    ep = (npar + 1j*eperp)/np.sqrt(2); em = (npar - 1j*eperp)/np.sqrt(2)
    return ph, ep, em

def mdot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]

def boost_x(xi):
    c, s = np.cosh(xi), np.sinh(xi)
    return np.array([[c, s, 0], [s, c, 0], [0, 0, 1.0]])

def h2_grid(xi, nchi=260, nphi=96, chimax=24.0):
    """nodes on H2 in polar coordinates centred at p-hat: y = B(xi) (cosh chi, sinh chi cos phi, sinh chi sin phi)."""
    x, w = leggauss(nchi); chi = 0.5*chimax*(x + 1); wc = 0.5*chimax*w
    phi = 2*np.pi*np.arange(nphi)/nphi; wp = np.full(nphi, 2*np.pi/nphi)
    C, PH = np.meshgrid(chi, phi, indexing='ij')
    y0 = np.stack([np.cosh(C), np.sinh(C)*np.cos(PH), np.sinh(C)*np.sin(PH)], -1)
    y = y0 @ boost_x(xi).T
    W = (wc[:, None]*np.sinh(C))*wp[None, :]
    return y, W, np.cosh(C)               # -y.p_hat = cosh chi

STRUCT = {
    'D1^2': lambda D1, D2: D1*D1, 'D1D2': lambda D1, D2: D1*D2, 'D2^2': lambda D1, D2: D2*D2,
}

def dd(a, nu, b):
    return nu*nu + (a - b)**2

def Mclosed(a, nu, s, flip):
    d5, d9 = dd(a, nu, 2.5), dd(a, nu, 4.5)
    N1 = 64*a**4 - 768*a**3 + 64*a*a*nu*nu + 3280*a*a - 320*a*nu*nu - 5840*a + 16*nu**4 + 488*nu*nu + 3705
    N2 = 64*a**3 - 32*a*a*nu*nu - 584*a*a + 224*a*nu*nu + 1720*a - 16*nu**4 - 424*nu*nu - 1641
    N3 = 128*a*a - 128*a*nu*nu - 800*a + 16*nu**4 + 360*nu*nu + 1113
    if not flip:
        return {'D1^2': N1/(64*d5*d9), 'D1D2': -N2/(64*d5*d9), 'D2^2': N3/(64*d5*d9)}[s]
    return {'D1^2': 1/(4*d5*d9), 'D1D2': -1/(4*d5*d9), 'D2^2': 1/(4*d5*d9)}[s]

def P4(nu, xi):
    """P^{(4)} = sinh^4 xi d^4 P/d(cosh xi)^4 (mpmath)."""
    f = lambda x: Pcon_mp(nu, x)
    return float(mp.re(mp.sinh(xi)**4*mp.diff(f, mp.cosh(xi), 4)))

def angular(a, nu, xi, flip, grid=None):
    y, W, ch = grid if grid is not None else h2_grid(xi)
    ph, ep, em = frame(xi)
    z1 = ep; z3 = ep if flip else em
    eta = y[..., 0]                                    # -y.x_hat, x_hat = (1,0,0)
    dist = np.arccosh(np.maximum(eta, 1.0))
    P = Pcon(nu, dist)
    z1y = mdot(y, z1); z3y = mdot(y, z3)
    D1 = np.dot(np.conj(np.zeros(3)) + z1, np.array([-1, 1, 1])*z3) + z1y*z3y
    D2 = -z1y*z3y
    out = {}
    for s, f in STRUCT.items():
        out[s] = (W*ch**(-a)*P*f(D1, D2)).sum()
    return out

pr("\n(B) stress-tensor rows of Lemma L (a = 7): Gamma(a) A_s / (C_a gamma_a M_s P) - 1 ; preserving (eps+, eps-) and flip (eps+, eps+, vs P^(4))")
a = 7.0
for nu in (0.0, 0.7, 1.5, 3.0):
    for xi in (0.4, 1.1):
        g = h2_grid(xi)
        Ca = 2**(a - 1)*np.sqrt(np.pi); ga = float(mp.re(gam(a, nu)))
        Pl = float(mp.re(Pcon_mp(nu, mp.cosh(xi)))); P4v = P4(nu, xi)
        Ap = angular(a, nu, xi, False, g); Af = angular(a, nu, xi, True, g)
        rows = []
        for s in STRUCT:
            rp = float(mp.gamma(a))*Ap[s]/(Ca*ga*Mclosed(a, nu, s, False)*Pl) - 1
            rf = float(mp.gamma(a))*Af[s]/(Ca*ga*Mclosed(a, nu, s, True)*P4v)
            rows.append(f"{s}: pres {abs(rp):.1e}, flip ratio {rf.real:+.10f}{rf.imag:+.1e}i")
        pr(f"   nu={nu}, xi={xi}: " + "; ".join(rows))

# ------------------------------------------------------------------ (C) complex lambda, radial integral done numerically
pr("\n(C) direct 3-d transform at lambda = e^{-i theta}, radial integral numerical (Gauss-Legendre in the damped variable u, 3000 nodes on [0, 90]):"
   " max over structures/helicities of |I(lambda)/(lambda^{-a} I(1)) - 1|")
_x, _w = leggauss(3000); UL = 45.0*(_x + 1); WL = 45.0*_w*np.exp(-UL)      # int_0^90 du e^{-u} (...)
def full3d(a, nu, xi, lam, flip, grid):
    y, W, ch = grid
    ph, ep, em = frame(xi)
    z1 = ep; z3 = ep if flip else em
    eta = y[..., 0]; P = Pcon(nu, np.arccosh(np.maximum(eta, 1.0)))
    z1y = mdot(y, z1); z3y = mdot(y, z3)
    D1 = np.dot(z1, np.array([-1, 1, 1])*z3) + z1y*z3y; D2 = -z1y*z3y
    # radial: int_0^oo R^{a-1} e^{-lam R ch} dR; substitute R = u/(Re(lam) ch): e^{-u} e^{-i Im(lam)/Re(lam) u}
    rl = lam.real; tq = lam.imag/rl
    rad = (UL**(a - 1)*np.exp(-1j*tq*UL)*WL).sum()/(rl*ch)**a      # independent of angle except via ch
    out = {}
    for s, f in STRUCT.items():
        out[s] = (W*P*f(D1, D2)*rad).sum()
    return out
for a in (7.0, 6.5):
    for nu in (0.7, 3.0):
        xi = 0.8; g = h2_grid(xi, nchi=200, nphi=64)
        base = {fl: full3d(a, nu, xi, 1.0 + 0j, fl, g) for fl in (False, True)}
        for th in (np.pi/6, np.pi/3, 0.45*np.pi, 0.47*np.pi):
            lam = np.exp(-1j*th)
            worst = 0.0
            for fl in (False, True):
                I = full3d(a, nu, xi, lam, fl, g)
                for s in STRUCT:
                    worst = max(worst, abs(I[s]/(lam**(-a)*base[fl][s]) - 1))
            pr(f"   a={a}, nu={nu}, theta={th/np.pi:.3f} pi: {worst:.2e}")
pr("   (theta -> pi/2 is the boundary value of an analytic function of lambda, Lemma D3(i); the quadrature is uniform up to "
   "tan(theta) ~ 10.)")

# ------------------------------------------------------------------ (D) two-pair Parseval with the shock profile
pr("\n(D) two-pair Parseval with the Einstein shock (B = Q_1(cosh L)): P_light / [h_a(nu)^2 * 2 pi/(nu^2+9/4)] - 1, and the nu-only factor")
def phi_int(L, s, a):
    """int_0^{2pi} dphi (cosh L cosh s - sinh L sinh s cos phi)^{-a} = 2 pi (cosh^2 L + sinh^2 s)^{-a/2} P_{a-1}(z),
    z = cosh L cosh s / sqrt(cosh^2 L + sinh^2 s)  (exact; Legendre P of degree a-1; checked against direct phi quadrature below)."""
    L, s = np.broadcast_arrays(np.asarray(L, float), np.asarray(s, float))
    q = np.cosh(L)**2 + np.sinh(s)**2
    z = np.cosh(L)*np.cosh(s)/np.sqrt(q)
    if float(a).is_integer():
        from numpy.polynomial.legendre import legval
        Pz = legval(z, np.eye(int(a))[int(a) - 1])
    else:
        Pz = np.vectorize(lambda v: float(mp.re(mp.legenp(a - 1, 0, v, type=3))))(z)
    return 2*np.pi*q**(-a/2)*Pz

def _phi_check():
    for (L, s, a) in ((0.3, 0.7, 7.0), (2.0, 2.5, 7.0), (1.0, 0.2, 6.5)):
        f = lambda ph: (mp.cosh(L)*mp.cosh(s) - mp.sinh(L)*mp.sinh(s)*mp.cos(ph))**(-a)
        ex = mp.quad(f, [0, mp.pi, 2*mp.pi])
        pr(f"   phi-integral closed form check (L={L}, s={s}, a={a}): rel err {abs(phi_int(L, s, a)/float(ex) - 1):.1e}")

def Q1L(L):
    """Legendre Q_1(cosh L) = cosh L log coth(L/2) - 1 (L <= 1.5); series sum_k cosh(L)^{-2k}/(2k+1) for L > 1.5."""
    L = np.asarray(L, float); out = np.empty_like(L); m = L <= 1.5
    out[m] = np.cosh(L[m])*np.log(1/np.tanh(L[m]/2)) - 1
    z = 1/np.cosh(L[~m])**2; ser = np.zeros_like(z); zk = z.copy()
    for k in range(1, 60):
        ser += zk/(2*k + 1); zk = zk*z
    out[~m] = ser
    return out

def radial_nodes(n, Lmax, k=3):
    """Gauss-Legendre in u, L = Lmax u^k (clusters nodes at L = 0 for the log singularity)."""
    x, w = leggauss(n); u = 0.5*(x + 1); wu = 0.5*w
    return Lmax*u**k, Lmax*k*u**(k - 1)*wu

def two_pair(a, nus, n=200, Lmax=24.0):
    L, wL = radial_nodes(n, Lmax)
    s, ws = radial_nodes(n, Lmax, k=1)
    t, wt = radial_nodes(n, Lmax, k=1)
    # F(s) = int_{H2} dp Q1(cosh d(p, pb)) (cosh d(p, x0))^{-a},  d(pb, x0) = s
    F = (wL*np.sinh(L)*Q1L(L))[None, :]*phi_int(L[None, :], s[:, None], a)
    F = F.sum(1)
    # G(t) = int_{H2} dpb F(d(pb, x0)) (cosh d(pb, y))^{-a},  d(y, x0) = t
    G = ((ws*np.sinh(s)*F)[None, :]*phi_int(s[None, :], t[:, None], a)).sum(1)
    out = {}
    for nu in nus:
        P = Pcon(nu, t)
        out[nu] = 2*np.pi*(wt*np.sinh(t)*P*G).sum()
    return out, (t, G)
_phi_check()
for a in (7.0, 6.5):
    nus = (0.0, 0.7, 1.5, 3.0)
    t0 = time.time(); Pl, (tg, Gt) = two_pair(a, nus)
    rows = []
    for nu in nus:
        h = HA.get((a, nu))
        if h is None:
            h = h_mp(mp.mpf(a), nu)
        h = float(mp.re(h))
        pred = h*h*2*np.pi/(nu*nu + 2.25)
        rows.append(f"nu={nu}: {abs(Pl[nu]/pred - 1):.1e} (h_a^2 = {h*h:.6e})")
    pr(f"   a = {a} ({time.time() - t0:.0f}s): " + "; ".join(rows))
pr("   G(t) (the light r^{-1} coefficient, angular part) is nu-independent by construction; all nu-dependence of P_light is the"
   " projection P_l, which produces h_a(nu)^2 x (harmonic transform of the shock) = h_a^2 * 2 pi/(nu^2+9/4).")
pr("   Heavy action of the same kernel on a plane-wave pair at distance beta: Gamma(a)^2 2^{-2a} h_a(nu)^2 P_l(cosh beta) (check A + (B.2)).")
pr("   Light/heavy ratio = e^{i pi a} (2 pi/(nu^2+9/4)) x [nu-independent]: at a = 7, e^{7 i pi} = -1, no branch ambiguity.")

# ------------------------------------------------------------------ (E) magic-block coefficient = heavy image of eta-multiplication
pr("\n(E) heavy image of the kernel eta P_l(eta): one pair, int_{H2} dy (cosh d(y,pb))^{-a} eta(x,y) P_l(eta(x,y)) vs the index-(l-1), "
   "(l+1) reproduction formula; with the second pair this gives (1/2)[g P_{l-1} + g~ P_{l+1}] h_7(nu)^2")
mp.mp.dps = 20
def g_l2(nu):
    nu = mp.mpf(nu)
    lg = lambda v: mp.loggamma((1 + 7 - 1.5 - 1j*v)/2) + mp.loggamma((1 + 7 - 1.5 + 1j*v)/2)
    gp = ((2 - 3)/2 + 1j*nu)/(1j*nu)*mp.exp(2*(lg(nu + 1j) - lg(nu)))
    gm = ((2 - 3)/2 - 1j*nu)/(-1j*nu)*mp.exp(2*(lg(-nu + 1j) - lg(-nu)))
    return gp, gm
for nu in (0.7, 2.0, 4.0):
    a = 7.0
    l = -0.5 + 1j*nu
    hA = lambda v: 2**(a - 1)*mp.sqrt(mp.pi)*gam(a, v)/mp.gamma(a)          # closed form of (A), continued to complex nu
    for xi in (0.5, 1.5):
        y, W, ch = h2_grid(xi, nchi=240, nphi=96)          # y-hat in polar coordinates centred at pb = (cosh xi, sinh xi, 0)
        eta = y[..., 0]                                      # eta = -x.y with x-hat = (1,0,0)
        lhs = (W*ch**(-a)*eta*Pcon(nu, np.arccosh(np.maximum(eta, 1.0)))).sum()
        x = mp.cosh(xi)
        rhs = (l*hA(nu + 1j)*mp.legenp(l - 1, 0, x, type=3) + (l + 1)*hA(nu - 1j)*mp.legenp(l + 1, 0, x, type=3))/(2*l + 1)
        pr(f"   one pair, nu={nu}, d(x,pb)={xi}: direct {lhs:+.12e} vs [l h(nu+i) P_(l-1) + (l+1) h(nu-i) P_(l+1)]/(2l+1) = "
           f"{complex(rhs).real:+.12e}{complex(rhs).imag:+.1e}i; rel diff {abs(lhs - complex(rhs))/abs(complex(rhs)):.1e}")
    gp, gm = g_l2(nu)
    gpred = 2*l/(2*l + 1)*(gam(7, nu + 1j)/gam(7, nu))**2
    pr(f"   two pairs => normalised image (1/2)[g P_(l-1) + g~ P_(l+1)] with g = 2l/(2l+1) (gamma_7(nu+i)/gamma_7(nu))^2: "
       f"|g_l2/g_pred - 1| = {mp.nstr(abs(gp/gpred - 1), 3)}, |g~ - conj g| = {mp.nstr(abs(gm - mp.conj(gp)), 3)}")
pr("   (identity: eta P_l = [l P_{l-1} + (l+1) P_{l+1}]/(2l+1), and P_{l-1/l+1} = P at nu +- i reproduce with h_7(nu +- i);"
   " h_7(nu+i)/h_7(nu) = gamma_7(nu+i)/gamma_7(nu) by (A); 2l/(2l+1) = ((2-d)/2 + i nu)/(i nu) at d = 3.)")

open(os.path.join(HERE, 'cache', 'v2_d3_planewave.txt'), 'w').write(buf.getvalue())

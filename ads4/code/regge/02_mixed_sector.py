r"""Check 2: <J J O O> and <T T O O> in the Regge limit, d = 3 (transverse space H^2).
Input (Costa-Hansen-Penedones arXiv:1707.07689):
  (2.47)  B(p, pbar, z1, z2) ~ \int dnu S^{j(nu)-1} sum_k beta_k(nu) Dhat_k Omega_{i nu}(L),
  (2.54)  currents:  Dhat_1 = delta^mu_nu,
                     Dhat_4 = -nabla^mu nabla_nu - (nu^2 + ((d-2)/2)^2)/(d-1) delta^mu_nu,
  (2.55)  stress tensors (curly brackets = symmetrise and remove traces in (mu nu) and in (rho sigma)):
          Dhat_1 = delta^{mu}_{rho} delta^{nu}_{sigma},
          Dhat_2 = delta^{mu}_{rho} nabla^{nu} nabla_{sigma} + (nu^2 + ((d-2)/2)^2)/(d-1) delta delta,
          Dhat_3 = (1/2)(nabla^mu nabla_rho nabla^nu nabla_sigma + nabla_rho nabla^mu nabla_sigma nabla^nu)
                   + (4 nu^2 + d^2 - 2d - 2)/(2(d+1)) delta^{mu}_{rho} nabla^{nu} nabla_{sigma}.
Here: d = 3, H^2 in geodesic polar coordinates (L, phi) centred at the impact point of the scalar pair,
ds^2 = dL^2 + sinh^2 L dphi^2; F(L) any solution of F'' + coth L F' + lam F = 0, lam = nu^2 + 1/4.
Polarisations: null vectors eps_pm = (e_L pm i e_phi)/sqrt 2 (helicity pm 1 about the transverse point).
Contractions with z1 for the indices of operator 1 and z2 for those of operator 2;
z1 = eps_+, z2 = eps_- = eps_+^* is the helicity-preserving ("weight 0") configuration of CHP (2.56);
z1 = z2 = eps_pm is the helicity-flip configuration (SO(2) weight pm 2 for currents, pm 4 for stress tensors).
All statements are exact identities modulo the radial equation (computed symbolically)."""
import sympy as sp, itertools, mpmath as mp
L, lam = sp.symbols('L lam', positive=True)
f0, f1 = sp.symbols('f0 f1')            # F, F'
ch, sh = sp.cosh(L), sp.sinh(L)
def dL(expr):                           # d/dL using the radial equation
    return sp.diff(expr, L) + sp.diff(expr, f0)*f1 + sp.diff(expr, f1)*(-ch/sh*f1 - lam*f0)
# Christoffel symbols Gamma[e][c][a] (coords 0 = L, 1 = phi)
Gam = [[[0, 0], [0, -sh*ch]], [[0, ch/sh], [ch/sh, 0]]]
def cov(T, n):
    """T: dict index-tuple (a1..an) -> component of a covariant tensor (phi-independent). Returns nabla_c T_{a1..an} as dict (c,a1..an)."""
    out = {}
    for idx in itertools.product((0, 1), repeat=n+1):
        c, rest = idx[0], idx[1:]
        val = dL(T[rest]) if c == 0 else 0
        for i in range(n):
            for e in (0, 1):
                g = Gam[e][c][rest[i]]
                if g != 0:
                    val -= g*T[rest[:i] + (e,) + rest[i+1:]]
        out[idx] = sp.simplify(val)
    return out
T0 = {(): f0}
T1 = cov(T0, 0); T2 = cov(T1, 1); T3 = cov(T2, 2); T4 = cov(T3, 3)
r2 = sp.sqrt(2)
eps = {+1: (1/r2, sp.I/(r2*sh)), -1: (1/r2, -sp.I/(r2*sh))}     # contravariant components
def dot(a, b):
    return sp.simplify(a[0]*b[0] + sh**2*a[1]*b[1])
def contract(T, vecs):
    n = len(vecs); tot = 0
    for idx in itertools.product((0, 1), repeat=n):
        term = T[idx]
        for i in range(n):
            term *= vecs[i][idx[i]]
        tot += term
    return sp.simplify(sp.expand(tot))
# associated functions P^m = sinh^m L (d/dcosh L)^m F, as combinations of F, F'
def Pm(m):
    g = f0
    for _ in range(m):
        g = dL(g)/sh
    return sp.simplify(sh**m*g)
print("eps_+.eps_- =", dot(eps[1], eps[-1]), "  eps_+.eps_+ =", dot(eps[1], eps[1]))
print("\n--- currents, d = 3 ---")
for (s1, s2, name) in ((1, -1, "weight 0 (z1=eps_+, z2=eps_-)"), (1, 1, "weight +2 (z1=z2=eps_+)"), (-1, -1, "weight -2 (z1=z2=eps_-)")):
    z1, z2 = eps[s1], eps[s2]
    D1 = dot(z1, z2)*f0
    D4 = sp.simplify(-contract(T2, [z1, z2]) - lam/2*dot(z1, z2)*f0)
    print(name, ":  Dhat_1 F =", D1, " ;  Dhat_4 F =", sp.simplify(D4), " ;  Dhat_4 F + P^2/2 =", sp.simplify(D4 + Pm(2)/2))
print("\n--- stress tensors, d = 3 ---")
for (s1, s2, name) in ((1, -1, "weight 0 (z1=eps_+, z2=eps_-)"), (1, 1, "weight +4 (z1=z2=eps_+)"), (-1, -1, "weight -4 (z1=z2=eps_-)")):
    z1, z2 = eps[s1], eps[s2]
    D1 = dot(z1, z2)**2*f0
    D2 = sp.simplify(dot(z1, z2)*contract(T2, [z1, z2]) + lam/2*dot(z1, z2)**2*f0)
    # index order of T4 = order of the derivative operators (outermost first)
    D3 = sp.Rational(1, 2)*(contract(T4, [z1, z2, z1, z2]) + contract(T4, [z2, z1, z2, z1])) + lam/2*dot(z1, z2)*contract(T2, [z1, z2])
    D3 = sp.simplify(D3)
    print(name, ":  Dhat_1 F =", D1, " ;  Dhat_2 F =", D2, " ;  Dhat_3 F - P^4/4 =", sp.simplify(D3 - Pm(4)/4), " ;  Dhat_3 F =", sp.simplify(D3))
print("\n--- general statement: (eps_+ . nabla)^m F = 2^{-m/2} P^m_{-1/2+i nu}(cosh L), m = 1..4 ---")
for m, T in ((1, T1), (2, T2), (3, T3), (4, T4)):
    print(" m =", m, ":", sp.simplify(contract(T, [eps[1]]*m) - Pm(m)/r2**m), sp.simplify(contract(T, [eps[-1]]*m) - Pm(m)/r2**m))
# numerical confirmation that the symbolic P^m is mpmath's associated conical function (type 3), and of
# P^m = prod_{j<m} (nu^2 + (j+1/2)^2) * P^{-m}  (m even),  P^{m} = (1/2+i nu)_m (1/2 - i nu)_m P^{-m}
mp.mp.dps = 25
nu_, L_ = mp.mpf('1.3'), mp.mpf('0.9')
F_ = lambda x: mp.re(mp.legenp(-0.5+1j*nu_, 0, mp.cosh(x), type=3))
vals = {f0: F_(L_), f1: mp.diff(F_, L_), lam: nu_**2 + mp.mpf(1)/4, L: L_}
for m in (2, 4):
    sym = complex(sp.N(Pm(m).subs({k: sp.Float(str(v), 25) for k, v in vals.items()}), 20))
    ref = mp.legenp(-0.5+1j*nu_, m, mp.cosh(L_), type=3)
    neg = mp.legenp(-0.5+1j*nu_, -m, mp.cosh(L_), type=3)
    poly = mp.mpf(1)
    for j in range(m):
        poly *= nu_**2 + (j + mp.mpf(1)/2)**2
    print(f" m={m}: symbolic P^m = {sym.real:.15f}   mpmath P^m = {mp.nstr(mp.re(ref),16)}   prod(nu^2+(j+1/2)^2) P^-m = {mp.nstr(mp.re(poly*neg),16)}")

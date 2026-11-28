#!/usr/bin/env python3
"""
Planar (N = infinity) test of  c_S / C_T <= 16 pi zeta(3)  on large-N Chern-Simons-vector models.

Conventions (all checked in checks() below):
  p = f T^3 on R^2,  c_S = 3 f,  C_T(free real scalar) = C_s = 3/(32 pi^2).
  R := (c_S/C_T) / (16 pi zeta(3)),  R = 1 for free scalars.

Inputs (leading order in N, U(N)_k, lambda = N/k in (0,1), YM-regularised level):
  * C_T = 2N C_s sin(pi lam)/(pi lam) per fundamental boson OR fermion species
    (Aharony-Gur-Ari-Yacoby 1207.4593 eq. (sumtwo): Ntilde = 2N sin(pi lam)/(pi lam);
     Gur-Ari-Yacoby 1211.1866 for fermions; independent of lambda_4, lambda_6 at planar order
     [1211.4843 sec. 7.3]).
  * Free energies on R^2 with the table-top holonomy rho(a)=1/(2 pi lam), |a|<=pi lam
    (1211.4843 sec. 2 / 1301.6169, valid for V T^2 >> N, i.e. the thermodynamic limit):
      F = (N V T^3/lam) [ calF(mu_B, lam) - calF(mu_F, lam - 1) ]           (1211.4843 eq. bos-fer-F)
      calF(mu, lam) = (1/(2 pi^2 i)) [ mu^2/3 Li2(e^{-mu-i pi lam}) + int_mu^inf y Li2(e^{-y-i pi lam}) dy - c.c. ]
    Gap equations (1211.4843 summary sec. 7, 1305.7235 eqs. (girfreesquare),(bosfree)):
      S(c) = c/2 - Im Li2( e^{-c+i pi lam})/(pi lam)      (bosonic holonomy-averaged tadpole)
      C(c) = c/2 - Im Li2(-e^{-c+i pi lam})/(pi lam)      (fermionic)
      regular fermion   : mu_F = 2 lam C(mu_F)
      critical boson    : S(mu_B) = 0
      regular boson     : mu_B = -lamhat S(mu_B),  lamhat^2 = lam^2 + lambda_6/(8 pi^2)
      critical fermion  : (1 + ghat) mu_F = 2 lam C(mu_F),   ghat = 1/sqrt(1 - 2 pi lam lambda_6^f)
      N=1 'w' family    : c_F = 2 lam (C - x4 S), c_B^2 = lam^2(1+3 x6) S^2 + 4 x4 lam^2 (C^2 - 2 x4 S C),
                          x4=(1+w)/2, x6=w^2 (1305.7235 eq. after (generalaction)); N=2 is w=1.
One finds (derived in CSM.md, re-verified numerically here) for each species
      F_species = - N V T^3 sin(pi lam)/(pi^2 lam) * G(mu, phi),
      G(mu,phi) = (1/12) int_mu^inf t (t^2-mu^2)/(cosh t - cos phi) dt
                = Im[ mu^2/3 Li2(z) + mu Li3(z) + Li4(z) ] / sin(phi),  z = e^{-mu+i phi},
      phi = pi lam (boson), phi = pi (1-lam) (fermion).
so that R = <G>/zeta(3), averaged over species with equal weights (each species carries the same C_T).
"""
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

mp.mp.dps = 30
Z3 = float(mp.zeta(3))
PI = np.pi


# ---------------------------------------------------------------- basic functions
def ImLi(n, mu, phi):
    return float(mp.im(mp.polylog(n, mp.exp(-mu + 1j * phi))))


def G_poly(mu, phi):
    """G via polylogs (rewriting of 1211.4843 calF)."""
    return (mu ** 2 / 3 * ImLi(2, mu, phi) + mu * ImLi(3, mu, phi) + ImLi(4, mu, phi)) / np.sin(phi)


def G_int(mu, phi):
    """G via the manifestly positive integral (independent evaluation)."""
    f = lambda t: t * (t * t - mu * mu) / (np.cosh(t) - np.cos(phi)) if t < 700 else 0.0
    val, _ = quad(f, mu, np.inf, limit=400, epsabs=1e-13, epsrel=1e-12)
    return val / 12


def calF(mu, lam):
    """calF(mu,lam) of 1211.4843 directly from its definition (numerical y-integral)."""
    zc = lambda y: mp.exp(-y - 1j * PI * lam)
    t1 = mu ** 2 / 3 * mp.polylog(2, zc(mu))
    t2 = mp.quad(lambda y: y * mp.polylog(2, zc(y)), [mu, mu + 10, mp.inf])
    val = (t1 + t2) - mp.conj(t1 + t2)
    return float(mp.re(val / (2 * PI ** 2 * 1j)))


def S_tad(c, lam):
    return c / 2 - ImLi(2, c, PI * lam) / (PI * lam)


def C_tad(c, lam):
    return c / 2 - float(mp.im(mp.polylog(2, -mp.exp(-c + 1j * PI * lam)))) / (PI * lam)


def S_tad_direct(c, lam):
    """Holonomy average done by brute-force quadrature (check of S_tad)."""
    a = PI * lam
    g = lambda al: 0.5 * 2 * np.real(np.log(2 * np.sinh((c + 1j * al) / 2)))
    v, _ = quad(g, -a, a, limit=200)
    return v / (2 * a)


def C_tad_direct(c, lam):
    a = PI * lam
    g = lambda al: 0.5 * 2 * np.real(np.log(2 * np.cosh((c + 1j * al) / 2)))
    v, _ = quad(g, -a, a, limit=200)
    return v / (2 * a)


# ---------------------------------------------------------------- gap equations
def mu_RF(lam):          # regular fermion: mu = 2 lam C(mu)
    if lam == 0:
        return 0.0
    return brentq(lambda m: m - 2 * lam * C_tad(m, lam), 1e-14, 60)


def mu_CB(lam):          # critical boson: S(mu) = 0
    if lam == 0:
        return 2 * np.log((1 + np.sqrt(5)) / 2)
    return brentq(lambda m: S_tad(m, lam), 1e-12, 60)


def mu_RB(lam, lamhat):  # regular boson: mu = -lamhat S(mu)
    if lamhat == 0:
        return 0.0
    if lam == 0:         # ungauged limit: (1 + 2/lamhat) mu = -2 log(1-e^-mu)
        return brentq(lambda m: (1 + 2 / lamhat) * m + 2 * np.log(-np.expm1(-m)), 1e-14, 60)
    return brentq(lambda m: m + lamhat * S_tad(m, lam), 1e-14, 60)


def mu_CF(lam, ghat):    # critical fermion: (1+ghat) mu = 2 lam C(mu)
    if lam == 0:
        return 0.0
    return brentq(lambda m: (1 + ghat) * m - 2 * lam * C_tad(m, lam), 1e-14, 60)


def mu_N2(lam):          # N=2, one chiral: mu = 2 lam (C - S)
    if lam == 0:
        return 0.0
    return brentq(lambda m: m - 2 * lam * (C_tad(m, lam) - S_tad(m, lam)), 1e-14, 60)


def masses_w(lam, w):
    """N=1 family of 1305.7235 at the conformal point; returns (c_B, c_F) or None."""
    x4, x6 = (1 + w) / 2, w * w

    def cB_of(cF, cB_guess):
        # solve c_B^2 = lam^2(1+3x6)S^2 + 4x4 lam^2(C^2 - 2 x4 S C) for c_B at fixed c_F
        C = C_tad(cF, lam)
        h = lambda cB: cB ** 2 - (lam ** 2 * (1 + 3 * x6) * S_tad(cB, lam) ** 2
                                  + 4 * x4 * lam ** 2 * (C ** 2 - 2 * x4 * S_tad(cB, lam) * C))
        grid = np.linspace(1e-6, 30, 600)
        vals = [h(g) for g in grid]
        for i in range(len(grid) - 1):
            if vals[i] * vals[i + 1] < 0:
                return brentq(h, grid[i], grid[i + 1])
        return None

    # fixed-point iteration on (c_B, c_F)
    cB, cF = 0.5, 0.5
    for _ in range(200):
        cB_new = cB_of(cF, cB)
        if cB_new is None:
            return None
        cF_new = abs(2 * lam * (C_tad(cF, lam) - x4 * S_tad(cB_new, lam)))
        # c_F appears on both sides; solve properly
        g = lambda c: c - abs(2 * lam * (C_tad(c, lam) - x4 * S_tad(cB_new, lam)))
        try:
            cF_new = brentq(g, 1e-12, 40)
        except ValueError:
            return None
        if abs(cB_new - cB) + abs(cF_new - cF) < 1e-12:
            cB, cF = cB_new, cF_new
            break
        cB, cF = cB_new, cF_new
    return cB, cF


# ---------------------------------------------------------------- R for each family
def R_boson(mu, lam):
    return G_poly(mu, PI * lam) / Z3 if lam > 0 else G_int(mu, 1e-300) / Z3


def R_fermion(mu, lam):
    return G_poly(mu, PI * (1 - lam)) / Z3


def R_RF(lam):
    return R_fermion(mu_RF(lam), lam)


def R_CB(lam):
    return R_boson(mu_CB(lam), lam) if lam > 0 else G_int(mu_CB(0), 0.0) / Z3


def R_RB(lam, lamhat):
    m = mu_RB(lam, lamhat)
    return R_boson(m, lam) if lam > 0 else G_int(m, 0.0) / Z3


def R_CF(lam, ghat):
    return R_fermion(mu_CF(lam, ghat), lam)


def R_N2(lam):
    m = mu_N2(lam)
    if lam == 0:
        return (Z3 + 0.75 * Z3) / (2 * Z3)
    return 0.5 * (G_poly(m, PI * lam) + G_poly(m, PI * (1 - lam))) / Z3


def R_w(lam, w):
    r = masses_w(lam, w)
    if r is None:
        return None
    cB, cF = r
    return 0.5 * (G_poly(cB, PI * lam) + G_poly(cF, PI * (1 - lam))) / Z3, cB, cF


def RB_domain(lam):
    """Stable lamhat range of the regular boson (Dey et al. 1808.04415, phi1<=x6<=phi2,
    x6_Dey = lambda_6/(24 pi^2 lam^2) - 1):  lamhat in [2 lam/(2-lam), 2]."""
    return 2 * lam / (2 - lam), 2.0


# ---------------------------------------------------------------- naive (zero holonomy) control
def R_RF_naive(lam):
    """Regular fermion with rho(a)=delta(a) (GMPTWY 1110.4386, known to be wrong at finite lam).
    mu = 2 lam log(2 cosh(mu/2));  f = -(1/6pi)[mu^3(1-1/lam) + 6(mu Li2(-e^-mu) + Li3(-e^-mu))]"""
    m = brentq(lambda x: x - 2 * lam * np.log(2 * np.cosh(x / 2)), 1e-14, 200)
    li2 = float(mp.polylog(2, -mp.exp(-m)))
    li3 = float(mp.polylog(3, -mp.exp(-m)))
    f = -(m ** 3 * (1 - 1 / lam) + 6 * (m * li2 + li3)) / (6 * PI)   # per N
    CT = 2 * (3 / (32 * PI ** 2)) * np.sin(PI * lam) / (PI * lam)   # per N
    return 3 * f / CT / (16 * PI * Z3)


# ---------------------------------------------------------------- checks
def checks():
    out = []
    Cs = 3 / (32 * PI ** 2)
    # free limits computed from scratch (no CS input)
    f_real_scalar = quad(lambda x: -x * np.log(-np.expm1(-x)), 0, np.inf)[0] / (2 * PI)
    f_dirac = 2 * quad(lambda x: x * np.log1p(np.exp(-x)), 0, np.inf)[0] / (2 * PI)
    out.append(("f(real scalar)/(zeta3/2pi)", f_real_scalar / (Z3 / (2 * PI))))
    out.append(("f(Dirac)/(zeta3/2pi)  [expect 3/2]", f_dirac / (Z3 / (2 * PI))))
    out.append(("R(free Dirac) [expect 0.75]", (3 * f_dirac / (2 * Cs)) / (16 * PI * Z3)))
    # G representations agree, and with calF
    for (m, l) in [(0.3, 0.2), (1.1, 0.5), (0.05, 0.9)]:
        Gp, Gi = G_poly(m, PI * l), G_int(m, PI * l)
        cf = calF(m, l)
        out.append((f"G_poly-G_int (mu={m},lam={l})", Gp - Gi))
        out.append((f"calF + sin(pi lam)G/pi^2 (mu={m},lam={l})", cf + np.sin(PI * l) * Gp / PI ** 2))
    # tadpoles vs brute-force holonomy average
    for (c, l) in [(0.4, 0.3), (1.5, 0.8)]:
        out.append((f"S_tad-direct (c={c},lam={l})", S_tad(c, l) - S_tad_direct(c, l)))
        out.append((f"C_tad-direct (c={c},lam={l})", C_tad(c, l) - C_tad_direct(c, l)))
    # free energy of N complex scalars at lam->0: f/N = zeta3/pi -> G(0,0)=zeta3
    out.append(("G(0,0)/zeta3", G_int(0, 0) / Z3))
    out.append(("G(0,pi)/zeta3 [3/4]", G_int(0, PI) / Z3))
    # critical O(2N) at lam=0: Sachdev 4/5
    out.append(("R_CB(0) [4/5]", R_CB(0)))
    # 1211.4843 perturbative check: F_F/(N V T^3) = -[3z3/4pi - (pi^2 ln2 + 4 ln^3 2)/(6 pi) lam^2]
    l = 1e-3
    fF = np.sin(PI * l) / (PI ** 2 * l) * G_poly(mu_RF(l), PI * (1 - l))
    pert = 3 * Z3 / (4 * PI) - (PI ** 2 * np.log(2) + 4 * np.log(2) ** 3) / (6 * PI) * l ** 2
    out.append(("RF: (f_exact - f_pert)/lam^2 at lam=1e-3 [->0]", (fF - pert) / l ** 2))
    out.append(("RF: mu_F/(2 ln2 lam) at lam=1e-3 [->1]", mu_RF(l) / (2 * np.log(2) * l)))
    # bosonization  RF(lam_f) = CB(1-lam_f)
    for lf in [0.1, 0.37, 0.8, 0.99]:
        out.append((f"R_RF({lf}) - R_CB({1-lf:.2f})", R_RF(lf) - R_CB(1 - lf)))
        out.append((f"mu_RF({lf}) - mu_CB({1-lf:.2f})", mu_RF(lf) - mu_CB(1 - lf)))
    # bosonization  CF(lam_f, ghat) = RB(1-lam_f, lamhat = 2 lam_b/ghat)
    for lf, gh in [(0.2, 1.0), (0.5, 0.7), (0.7, 1.6)]:
        lb = 1 - lf
        out.append((f"R_CF({lf},g={gh}) - R_RB({lb:.2f},lh={2*lb/gh:.3f})",
                    R_CF(lf, gh) - R_RB(lb, 2 * lb / gh)))
    # primary-paper benchmarks at lam=1/2
    for l6 in [6 * PI ** 2, 12 * PI ** 2]:
        lh = np.sqrt(0.25 + l6 / (8 * PI ** 2))
        out.append((f"RB lam=1/2 lambda6={l6/PI**2:.0f}pi^2: mu, c_S/C_T", (mu_RB(0.5, lh), 16 * PI * Z3 * R_RB(0.5, lh))))
    # ungauged BMB endpoint: lamhat=2 at lam=0 -> mu=log2, ratio 7/8
    out.append(("RB(lam=0,lamhat=2): mu [log 2], R [7/8]", (mu_RB(0, 2.0), R_RB(0, 2.0))))
    # N=2 self-duality and w=1 reproduces N=2
    out.append(("R_N2(0.3)-R_N2(0.7)", R_N2(0.3) - R_N2(0.7)))
    rw = R_w(0.3, 1.0)
    out.append(("R_w(0.3,w=1) - R_N2(0.3)", rw[0] - R_N2(0.3)))
    out.append(("c_B - c_F at w=1", rw[1] - rw[2]))
    return out


if __name__ == "__main__":
    np.set_printoptions(precision=6)
    print("=== checks ===")
    for k, v in checks():
        print(f"{k:60s} {v}")

#!/usr/bin/env python3
"""REGGE_TREE.md support script (referee report 3, M2).

Checks, in a stringy Regge-pole model of the O(1/c_T) heavy double commutator on the
transverse space H_2 (d = 3), the statements of REGGE_TREE.md sec. 3:

 A. heat-kernel identity  int_0^oo e^{-2t} K_t(b) dt = Q_1(cosh b)/(2 pi) = G_2(b)
    (G_2 = H_2 propagator of the Delta = 2 transverse graviton, ~ 1/(6 pi eta^2)).
 B. model  D_eps(r, b) = (eps/r) e^{-2 eps l} K_{eps l}(b),  l = log(1/r),
    i.e. the conformal-Regge form int dnu alpha(nu) r^{1-j(nu)} Omega_nu with
    j(nu) = 2 - eps (nu^2 + 9/4) and an inelastic residue prop. to 2 - j(nu)
    (normalised so that its k = 2 Regge moment reproduces the graviton pole 1/(nu^2+9/4)).
    D_eps >= 0; fixed-eta intercept 2 - 9 eps/4; k = 2 Regge moment (CHMRS (2.39))
    Pi_2(eta) = int_{eps l0}^oo e^{-2t} K_t dt <= G_2(eta) ~ eta^{-2}   (CHMRS (4.45)).
 C. corner integrals int^H d eta eta^c * (r-moment of the majorant) for each Regge spin k,
    with three majorants: uniform r^{1-J0} at J0 = 1 (non-perturbative), uniform at
    J0 = 2 (planar, referee's count), and the tree moment bound of B.
 D. the fixed-eta tree margin 2 - j(0) = 9 eps/4 with eps = 2/Lambda^2.

Usage:  ulimit -v 3000000; OMP_NUM_THREADS=1 python3 regge_tree.py   (~1 min)
Output: cache/regge_tree.txt
"""
import os
import sys
import mpmath as mp
import numpy as np
from scipy.integrate import quad

mp.mp.dps = 15
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache", "regge_tree.txt")
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def logK(t, b):
    """log of the heat kernel of the Laplacian on H_2 (spectrum -(nu^2+1/4)), McKean's formula,
    with s = b + u^2 to remove the endpoint singularity; double precision (scipy)."""
    t = float(t); b = float(b)
    lpref = 0.5 * np.log(2) - t / 4 - 1.5 * np.log(4 * np.pi * t) - b * b / (4 * t)
    def lsinh(x):
        return x - np.log(2) + np.log1p(-np.exp(-2 * x)) if x > 1e-3 else np.log(np.sinh(x))

    def f(u):
        if u <= 0:
            return 0.0
        s = b + u * u
        # cosh s - cosh b = 2 sinh((s+b)/2) sinh(u^2/2), in logs (no cancellation / overflow)
        lden = 0.5 * (np.log(2) + lsinh((s + b) / 2) + lsinh(u * u / 2))
        return np.exp(np.log(2 * u * s) - (s * s - b * b) / (4 * t) - lden)
    w = np.sqrt(max(1.0, 40 * t / max(b, 0.05)))  # Gaussian width in u
    val = 0.0
    edges = [0, w / 20, w / 4, w, 3 * w, 10 * w, np.inf]
    for x0, x1 in zip(edges[:-1], edges[1:]):
        val += quad(f, x0, x1, limit=400, epsabs=0, epsrel=1e-12)[0]
    return lpref + np.log(val)


def G2(b):
    z = np.cosh(b)
    return (z / 2 * np.log((z + 1) / (z - 1)) - 1) / (2 * np.pi)


def tau_int(b, lo):
    """int_lo^oo e^{-2t} K_t(b) dt."""
    f = lambda t: np.exp(-2 * t + logK(t, b))
    pts = sorted(set([lo] + [x for x in (b * b / 40, b * b / 8, b * b / 2, b * b, 2 * b * b, 6 + 2 * b) if x > lo]))
    val = 0.0
    for x0, x1 in zip(pts[:-1], pts[1:]):
        val += quad(f, x0, x1, limit=400, epsrel=1e-11)[0]
    val += quad(f, pts[-1], np.inf, limit=400, epsrel=1e-11)[0]
    return val


say("REGGE_TREE support output")
say("=" * 70)

# ---------------- A. heat-kernel identity -----------------
say("A. int_0^oo e^{-2t} K_t(b) dt  vs  G_2(b) = Q_1(cosh b)/(2 pi)")
for b in (0.5, 1.0, 2.0, 4.0, 8.0):
    lhs = tau_int(b, 1e-12)
    rhs = G2(b)
    say(f"  b={b:4.1f}  int={lhs:.12e}  G2={rhs:.12e}  rel.diff={abs(lhs / rhs - 1):.1e}"
        f"   6*pi*eta^2*G2={6 * np.pi * np.cosh(b) ** 2 * rhs:.6f}")

# ---------------- B. stringy model -----------------
say("")
say("B. model D_eps(r,b) = (eps/r) e^{-2 eps l} K_{eps l}(b)  (l = log 1/r), positive (heat kernel)")
say("   (i) local fixed-eta intercept J = 1 + dlogD/dl = 2 - eps*(2 - dlogK/dtau), tau = eps*l;")
say("       (2-J)/eps -> 9/4 + 3/(2 tau) (eps-independent):")
for b in (0.5, 1.0, 3.0):
    row = []
    for tau in (5.0, 20.0, 80.0, 320.0):
        h = 1e-3 * tau
        d = (logK(tau + h, b) - logK(tau - h, b)) / (2 * h)
        row.append(f"tau={tau:5.0f}: (2-J)/eps={2 - d:.5f} [9/4+3/(2tau)={2.25 + 1.5 / tau:.5f}]")
    say(f"  b={b}: " + "; ".join(row))
for dg in (300, 1000):
    L2 = (dg - 2) * (dg - 1)
    eps = 2.0 / L2
    say(f"  Delta_g={dg}: eps=2/Lambda^2={eps:.3e}  => J(l->oo) = 2 - {9 * eps / 4:.3e}; reached only for l >~ 1/eps = {1 / eps:.2e}")
say("   (ii) k=2 Regge moment Pi_2(eta) = int_{eps l0}^oo e^{-2t}K_t dt, l0 = -log r_max(eta), r_max = eta - sqrt(eta^2-1)")
for eps in (0.05, 2.0 / 997002):
    for b in (0.5, 2.0, 6.0):
        eta = np.cosh(b)
        rmax = eta - np.sqrt(eta * eta - 1)
        l0 = -np.log(rmax)
        Pi = tau_int(b, eps * l0)
        say(f"  eps={eps:9.3e} b={b:4.1f} eta={eta:10.4f}  Pi_2={Pi:.10e}  G2={G2(b):.10e}  Pi_2/G2={Pi / G2(b):.8f}")
say("   => 0 <= Pi_2(eta) <= G_2(eta) ~ eta^{-2}/(6 pi): the CHMRS (4.45) tail (transverse graviton propagator)")

# ---------------- C. corner integrals -----------------
say("")
say("C. corner integrals I(H) = int_2^H d eta  eta^c  m(eta),  m = int_0^{1/(2 eta)} dr r^a M(r,eta)")
say("   majorants: N1 = r^0 (J0=1, finite c_T), P2 = r^{-1} (J0=2 uniform, planar), T = tree moment <= r_max^a G_2")
pieces = [  # (label, a, c, log power)
    ("k=2  C'_EE (de-shifted, cut)", 0, mp.mpf(-0.5), 1),
    ("k=3  C_3", 0, mp.mpf(-0.5), 1),
    ("k=4  C4ds, dC4ds (Q weights)", 2, mp.mpf(0.5), 0),
    ("k=4  generic weights", 2, mp.mpf(1), 0),
    ("k=5  C5ds (raw)", 2, mp.mpf(1.5), 0),
]


def moment(kind, a, eta):
    rm = 1 / (2 * eta)
    if kind == "N1":
        return rm ** (a + 1) / (a + 1)
    if kind == "P2":
        if a == 0:
            return mp.inf  # int_0 dr r^{-1}: fixed-eta divergence (J0 = 2 = k)
        return rm ** a / a
    if kind == "T":
        z = eta
        return rm ** a * (z / 2 * mp.log((z + 1) / (z - 1)) - 1) / (2 * mp.pi)
    raise ValueError


def corner(kind, a, c, lg, H):
    f = lambda e: e ** c * (mp.log(e) ** lg) * moment(kind, a, e)
    pts = [2] + [10 ** j for j in range(1, int(mp.log10(H)) + 1) if 10 ** j < H] + [H]
    return mp.quad(f, pts)


for lab, a, c, lg in pieces:
    say(f"  {lab}:  a={a}, c={mp.nstr(c, 2)}{' (x log)' if lg else ''}")
    for kind, thr in (("N1", a - c), ("P2", a - 1 - c), ("T", a + 1 - c)):
        # thr = margin in the eta exponent: convergence iff thr > 0
        if kind == "P2" and a == 0:
            say(f"     {kind}: diverges at r->0 for every eta (fixed-eta Regge: J0 = 2 = k)")
            continue
        vals = [corner(kind, a, c, lg, mp.mpf(10) ** p) for p in (2, 4, 6, 8)]
        say(f"     {kind}: I(1e2,1e4,1e6,1e8) = " + ", ".join(mp.nstr(v, 6) for v in vals)
            + f"   eta-exponent margin = {mp.nstr(thr, 3)}  -> {'converges' if thr > 0 else 'DIVERGES'}")

# ---------------- D. numbers -----------------
say("")
say("D. fixed-eta tree margin 2 - j(0) = 9 eps/4, eps = 2/(M L)^2, M L ~ Lambda (flat closed string, alpha' t/2)")
for dg in (300, 1000):
    L2 = mp.mpf((dg - 2) * (dg - 1))
    say(f"  Delta_g={dg}: Lambda^2={int(L2)}  2-j(0)={mp.nstr(mp.mpf(9) / (2 * L2), 4)}")
say("  cf. N=4 SYM (CGP 1209.4355 (95),(98)): j(0) = 2 - 2/sqrt(lambda) = 2 - 8/Delta_gap^2 with Delta_gap = 2 lambda^{1/4}")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("wrote", OUT, file=sys.stderr)

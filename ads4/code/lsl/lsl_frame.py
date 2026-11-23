"""
lsl/lsl_frame.py (LARGESPIN_L.md §2): checks of the exact frame identities behind the spin-Gamma representation.

  (F1) In the gauge x5 = 0, x1 = -e, x2 = -L e (L -> inf), x3 = (e+y)/(e+y)^2, x4 = x/(-x^2)  (CHMRS 2106.10274 eq. firstconformalframe),
       the crystal's heavy invariants (l1_crystal.build: B = m35 m45 m12/(m15 m25 m34), t = Z.W/sqrt(Z^2 W^2)) are
            B = 1/(-w^2),   t = cosh(chi_w) = w^0/|w|,    w := e + x + y.
  (F2) B and t are conformal invariants (random boosts, dilatations, translations, small special conformal maps).
  (F3) Laplace's second integral, d = 3:  B^{Delta/2} P_J(t) = (1/2pi) int_0^{2pi} dphi (-w^2)^{-(tau-1)/2} (-n_phi.w)^{-(J+1)},
       n_phi = (1, cos phi, sin phi), tau = Delta - J.
  (F4) Euler:  (-n.w)^{-(J+1)} = int_0^inf ds s^J e^{s n.w} / J!       (n.w < 0 for w future timelike).
  (F5) Riesz/Bernstein (d = 3, lam > 1/2):  (-w^2)^{-lam} = (2lam-1)/(2 pi Gamma(2lam)) int_{V+} d^3p (-p^2)^{lam-3/2} e^{p.w}.
Metric mostly plus, (t, x1, x2).  Writes cache/frame_checks.txt.  Pure numerics, no project imports.
"""
import os, sys
sys.dont_write_bytecode = True
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'cache', 'frame_checks.txt')
rng = np.random.default_rng(20261123)
E = np.array([1.0, 0.0, 0.0])


def dot(a, b):
    return -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def msq(p, q):
    r = p - q
    return -dot(r, r)


def rand_future(scale=1.0):
    v = rng.normal(size=2)*scale
    t = np.sqrt(v @ v) + abs(rng.normal())*scale + 0.05*scale
    return np.array([t, v[0], v[1]])


def invariants(x1, x2, x3, x4, x5):
    m35, m45, m12, m15, m25, m34 = msq(x3, x5), msq(x4, x5), msq(x1, x2), msq(x1, x5), msq(x2, x5), msq(x3, x4)
    B = m35*m45*m12/(m15*m25*m34)
    x15, x25, x53, x54 = x1 - x5, x2 - x5, x5 - x3, x5 - x4
    Z = x15*m25 - x25*m15
    W = x53*m45 - x54*m35
    t = dot(Z, W)/np.sqrt(dot(Z, Z)*dot(W, W))
    return B, t


def first_frame(x, y, L=1e7):
    x1 = -E; x2 = -L*E; x5 = np.zeros(3)
    ey = E + y
    x3 = ey/dot(ey, ey)
    x4 = x/(-dot(x, x))
    return x1, x2, x3, x4, x5


def boost(eta, ang):
    c, s = np.cosh(eta), np.sinh(eta)
    n = np.array([np.cos(ang), np.sin(ang)])
    M = np.eye(3)
    M[0, 0] = c; M[0, 1:] = s*n; M[1:, 0] = s*n
    M[1:, 1:] = np.eye(2) + (c - 1)*np.outer(n, n)
    return M


def sct(xv, b):
    """special conformal map x -> (x - b x^2)/(1 - 2 b.x + b^2 x^2)."""
    x2 = dot(xv, xv)
    return (xv - b*x2)/(1 - 2*dot(b, xv) + dot(b, b)*x2)


def main():
    lines = []
    # (F1)
    errB, errt, tmin = 0.0, 0.0, np.inf
    for _ in range(2000):
        x = rand_future(rng.uniform(0.05, 3)); y = rand_future(rng.uniform(0.05, 3))
        B, t = invariants(*first_frame(x, y))
        w = E + x + y
        Bw = 1/(-dot(w, w)); tw = w[0]/np.sqrt(-dot(w, w))
        errB = max(errB, abs(B/Bw - 1)); errt = max(errt, abs(abs(t)/tw - 1)); tmin = min(tmin, abs(t))
    lines.append(f"(F1) first frame, 2000 random x,y in V+ (L = 1e7): max|B (-w^2) - 1| = {errB:.2e}; max| |t| w^0/|w| ... |t|/(w0/|w|) - 1| = {errt:.2e}; min |t| = {tmin:.4f} (>= 1: t = cosh chi_w)")
    # sign of t
    x = rand_future(); y = rand_future(); B, t = invariants(*first_frame(x, y))
    lines.append(f"     sign of t on the causal region: {np.sign(t):+.0f} (the overall (+-1)^J of l1_crystal is a convention; |t| enters P_J up to (-1)^J)")
    # (F2) conformal invariance
    errB2, errt2 = 0.0, 0.0; nflip = 0
    for _ in range(500):
        x = rand_future(rng.uniform(0.05, 2)); y = rand_future(rng.uniform(0.05, 2))
        pts = list(first_frame(x, y, L=50.0))
        B0, t0 = invariants(*pts)
        M = boost(rng.normal()*0.7, rng.uniform(0, 2*np.pi)); lam_ = np.exp(rng.normal()*0.5); a = rng.normal(size=3)*0.3
        b = rng.normal(size=3)*0.003
        pts2 = [sct(M @ p*lam_ + a, b) for p in pts]
        B1, t1 = invariants(*pts2)
        errB2 = max(errB2, abs(B1/B0 - 1)); errt2 = max(errt2, abs(abs(t1)/abs(t0) - 1)); nflip += (np.sign(t1) != np.sign(t0))
    lines.append(f"(F2) 500 random conformal maps (boost, dilatation, translation, SCT |b| ~ 3e-3): max rel change B {errB2:.2e}, |t| {errt2:.2e}; sign of t flipped in {nflip} maps (a point crossed infinity: (+-1)^J patch convention)")
    # (F3) Laplace's second integral
    mp.mp.dps = 20
    print('\n'.join(lines), flush=True)
    worst = 0
    for (tau, J) in [(1000, 0), (1000, 5), (1000, 50), (1000, 344), (300, 100), (12.5, 7)]:
        for _ in range(2):
            w = E + rand_future(rng.uniform(0.01, 1.0))
            wm = mp.sqrt(-(-w[0]**2 + w[1]**2 + w[2]**2)); ch = mp.mpf(w[0])/wm
            lhs = wm**(-(tau + J)) * mp.legenp(J, 0, ch, type=3)
            f = lambda ph: wm**(-(tau - 1)) * (mp.mpf(w[0]) - mp.cos(ph)*w[1] - mp.sin(ph)*w[2])**(-(J + 1))
            rhs = mp.quad(f, mp.linspace(0, 2*mp.pi, 129))/(2*mp.pi)
            worst = max(worst, abs(rhs/lhs - 1))
    lines.append(f"(F3) (-w^2)^(-Delta/2) P_J(w0/|w|) = (1/2pi) int dphi (-w^2)^(-(tau-1)/2) (-n.w)^(-(J+1)): max rel err {float(worst):.1e} (tau,J up to (1000,344))")
    # (F4) Euler
    worst = 0
    for J in [0, 5, 50, 344, 3000]:
        a = mp.mpf('1.37')
        pk = J/a; wd = 12*mp.sqrt(J + 1)/a
        val = mp.quad(lambda s: mp.exp(J*mp.log(s) - s*a - mp.loggamma(J + 1)) if s > 0 else mp.mpf(0), [0] + [max(mp.mpf(0), pk - wd + i*wd/8) for i in range(17)] + [pk + 3*wd, mp.inf])
        worst = max(worst, abs(val/a**(-(J + 1)) - 1))
    lines.append(f"(F4) Euler int s^J e^(-a s)/J! = a^-(J+1): max rel err {float(worst):.1e} (J <= 3000)")
    # (F5) Riesz
    worst = 0
    for lam_ in [0.75, 2.0, 7.5]:
        for _ in range(2):
            w = rand_future(rng.uniform(0.3, 2.0)); wm = np.sqrt(-dot(w, w)); chi = np.arccosh(w[0]/wm)
            # rest frame of w: integrand e^{-|p| wm cosh(rel)} -> reduce to the known 1D form, then do the angular integral numerically
            def I(chi_p):
                # integral over |p| and the azimuth of p, at p-rapidity chi_p, relative to w with rapidity chi
                def az(ph):
                    c = mp.cosh(chi_p)*mp.cosh(chi) - mp.sinh(chi_p)*mp.sinh(chi)*mp.cos(ph)
                    return mp.gamma(2*lam_)/(wm*c)**(2*lam_)
                return mp.quad(az, [0, mp.pi, 2*mp.pi])*mp.sinh(chi_p)
            tot = mp.quad(I, [0, 0.5, 1.5, 4, 10, 30])
            pred = 2*mp.pi*mp.gamma(2*lam_)/(2*lam_ - 1)*wm**(-2*lam_)
            worst = max(worst, abs(tot/pred - 1))
    lines.append(f"(F5) Bernstein/Riesz representation of (-w^2)^-lam, lam in {{0.75, 2, 7.5}}, boosted w: max rel err {float(worst):.1e}")
    txt = "\n".join(lines)
    print(txt)
    with open(OUT, 'w') as fh:
        fh.write(txt + "\n")


if __name__ == '__main__':
    main()

"""
bh_spinJ (BELOW_FOLD3.md §3): closed-form spin-J contraction of the O index (the l5_spin / bg_ee null-circle Fourier x kappa pairing).
  [DV] <(A.z)^J, (B.z)^J>_kappa := sum_n f_n g_{-n}/kappa_n  =  N_J (A^2 B^2)^{J/2} P_J(A.B/sqrt(A^2 B^2)),   N_J = J!/(2J-1)!!
  [D]  polarisation:  for f = (V.z)^{J-p} prod_{i<=p} (A_i.z),  g = (W.z)^{J-q} prod_{j<=q} (B_j.z):
       <f, g> = N_J (J-p)!(J-q)!/(J!)^2  prod_i (A_i.d_a) prod_j (B_j.d_b) H(a, b) |_{a=V, b=W},   H = (a^2 b^2)^{J/2} P_J(a.b/|a||b|),
  i.e. a finite sum of P_J^{(n)}(t), n <= p + q, t = V.W/sqrt(V^2 W^2), times polynomials in the dot products -- no 2J+1 samples and no
  V^J powers.  Here: exact multilinear (nilpotent) evaluation and a pointwise check against the sampled pairing (any J).
  The structures of <T T O> have p, q <= 4 (H15^b15 H25^b25, b <= 2), so n <= 8.
Usage: python3 bh_spinJ.py  -> cache/bh_spinJ.txt
"""
import os, sys, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
import mpmath as mp
from scipy.special import comb

dot = lambda p, q: -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]


# ---------- multilinear nilpotent numbers: dict {frozenset(generators): coeff}, eps_i^2 = 0
def nmul(x, y):
    out = {}
    for s, a in x.items():
        for t, b in y.items():
            if s & t:
                continue
            u = s | t; out[u] = out.get(u, 0) + a*b
    return out


def nadd(x, y, cy=1.0):
    out = dict(x)
    for t, b in y.items():
        out[t] = out.get(t, 0) + cy*b
    return out


def nfun(x, derivs):
    """f(x0 + n) = sum_k f^(k)(x0) n^k / k!  (n nilpotent part), derivs[k] = f^(k)(x0)."""
    x0 = x.get(frozenset(), 0); n = {s: v for s, v in x.items() if s}
    out = {frozenset(): derivs[0]}; pw = {frozenset(): 1.0}; fact = 1.0
    for k in range(1, len(derivs)):
        pw = nmul(pw, n); fact *= k
        if not pw:
            break
        out = nadd(out, pw, derivs[k]/fact)
    return out


def ndot(p, q):
    return nadd(nadd(nmul(p[0], q[1]), nmul(p[1], q[0])), {}, 1.0) if False else \
        nadd(nadd({s: -0.5*v for s, v in nmul(p[0], q[1]).items()}, {s: -0.5*v for s, v in nmul(p[1], q[0]).items()}), nmul(p[2], q[2]))


def legendre_derivs(J, t, kmax):
    """P_J^(k)(t), k = 0..kmax (complex t), via P_J^(k) = (2k-1)!! C^{(k+1/2)}_{J-k}."""
    out = []
    for k in range(kmax + 1):
        if k > J:
            out.append(0.0); continue
        out.append(complex(mp.fac2(2*k - 1)*mp.gegenbauer(J - k, k + mp.mpf(1)/2, mp.mpc(t))))
    return out


def pair_closed(J, V, As, W, Bs):
    p, q = len(As), len(Bs)
    ga = [frozenset([('a', i)]) for i in range(p)]; gb = [frozenset([('b', j)]) for j in range(q)]
    a = [{frozenset(): V[c], **{ga[i]: As[i][c] for i in range(p)}} for c in range(3)]
    b = [{frozenset(): W[c], **{gb[j]: Bs[j][c] for j in range(q)}} for c in range(3)]
    al, A2, B2 = ndot(a, b), ndot(a, a), ndot(b, b)
    nmax = p + q
    S = nmul(A2, B2)
    s0 = S[frozenset()]
    rS = nfun(S, [s0**(J/2 - k)*float(mp.rf(J/2 - k + 1, k)) for k in range(nmax + 1)])          # S^{J/2}
    iS = nfun(S, [s0**(-0.5 - k)*float(mp.rf(-0.5 - k + 1, k)) for k in range(nmax + 1)])        # S^{-1/2}
    t = nmul(al, iS)
    Pt = nfun(t, legendre_derivs(J, t[frozenset()], nmax))
    H = nmul(rS, Pt)
    full = frozenset(ga + gb) if (p + q) else frozenset()
    key = frozenset().union(*ga, *gb) if (p + q) else frozenset()
    NJ = float(mp.fac(J)/mp.fac2(2*J - 1))
    return NJ*float(mp.fac(J - p)*mp.fac(J - q)/mp.fac(J)**2)*H.get(key, 0)


def pair_sampled_mp(J, V, As, W, Bs, dps=60):
    """the same sampled pairing in mpmath (dps digits): the double-precision DFT loses digits at J >~ 16 (cancellations)."""
    with mp.workdps(dps):
        n = 2*J + 1
        cv = lambda v: [mp.mpc(complex(x)) for x in v]
        V, W = cv(V), cv(W); As = [cv(A) for A in As]; Bs = [cv(B) for B in Bs]
        d = lambda p, q: -(p[0]*q[1] + p[1]*q[0])/2 + p[2]*q[2]
        Z = [[1 - mp.cos(2*mp.pi*k/n), 1 + mp.cos(2*mp.pi*k/n), mp.sin(2*mp.pi*k/n)] for k in range(n)]
        f = [d(V, z)**(J - len(As))*mp.fprod([d(A, z) for A in As]) for z in Z]
        g = [d(W, z)**(J - len(Bs))*mp.fprod([d(B, z) for B in Bs]) for z in Z]
        F = lambda h, m: mp.fsum(h[k]*mp.expjpi(-2*m*k/mp.mpf(n)) for k in range(n))/n
        tot = 0
        for m in range(-J, J + 1):
            kap = (-1)**(J + m)*mp.mpf(2)**(-J)*mp.binomial(2*J, J + m)
            tot += F(f, m)*F(g, -m)/kap
        return complex(tot)


def pair_sampled(J, V, As, W, Bs):
    n = 2*J + 1; ph = 2*np.pi*np.arange(n)/n
    Z = np.stack([1 - np.cos(ph), 1 + np.cos(ph), np.sin(ph)], axis=1)
    f = np.array([dot(V, z)**(J - len(As))*np.prod([dot(A, z) for A in As]) for z in Z])
    g = np.array([dot(W, z)**(J - len(Bs))*np.prod([dot(B, z) for B in Bs]) for z in Z])
    F = np.fft.fft(f)/n; G = np.fft.fft(g)/n
    idx = np.arange(-J, J + 1) % n
    F, G = F[idx], G[idx]
    k = np.arange(-J, J + 1); kap = (-1.0)**(J + k)*2.0**(-J)*comb(2*J, J + k)
    return np.sum(F*G[::-1]/kap)


if __name__ == '__main__':
    mp.mp.dps = 40
    rng = np.random.default_rng(3)
    lines = ['# bh_spinJ: closed-form vs sampled O-index pairing, random complex vectors; columns J, p, q, |closed/sampled - 1|']
    worst = 0
    for J in (4, 6, 8, 12, 16, 24, 28, 32):
        for p, q in ((0, 0), (1, 0), (2, 1), (2, 2), (4, 2), (4, 4)):
            if p > J or q > J:
                continue
            R = lambda: rng.normal(size=3) + 1j*rng.normal(size=3)
            V, W = R(), R(); As = [R() for _ in range(p)]; Bs = [R() for _ in range(q)]
            # near-forward kinematics as in the engine: W close to V (t ~ 1) on half of the tests
            if (p + q) % 2 == 0:
                W = V + 0.05*R()
            c, s = pair_closed(J, V, As, W, Bs), pair_sampled(J, V, As, W, Bs)
            sm = pair_sampled_mp(J, V, As, W, Bs)
            e = abs(c/sm - 1); worst = max(worst, e)
            lines.append(f'  J={J:3d} p={p} q={q}: closed vs sampled(mp, 60 digits) {e:.2e};  double-precision sampled vs mp {abs(s/sm - 1):.1e}'
                         f'   (|pair| = {abs(sm):.3e})')
    lines.append(f'  worst closed vs sampled(mp): {worst:.2e}')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'bh_spinJ.txt'), 'w').write('\n'.join(lines) + '\n')

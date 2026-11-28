"""v3: well-conditioned version of the v2 dual LP, for larger g.

The dual functional f_b(J,D) = sum_{k<=N} mu_k r_k^(b)(J,D) is rewritten in a state-independent basis
    B_n = sum_k M_nk rhat_k ,  rhat_k = r_k / binom(k+3,3)^2 = prod_{j<k} (C - y_{b+j})/(T - y_{b+j}),
    M_nk = (-n)_k (n+1)_k / (k! (2)_k)       (any invertible triangular M gives the SAME LP; this one keeps
                                               the B_n of order one, so double-precision HiGHS is well conditioned)
B_n is evaluated in 256-bit ball arithmetic (python-flint).  The LP is
    min lam  s.t.  sigma_i f(i)/W_i <= lam  for all generated states,   mu_2 = 1,
with constraint generation over the full candidate set, then an exact-rational re-check of the
certificate on the most binding states.  |c|^2_max = 3 lam.
usage: python3 v3_lp_hp.py g N [light] [odd_b0]
  light : add light columns (b=0, J=0 with 1/2<=D<g ; b=0, J=2 with 3<D<g+2): light single traces / any light state
          of spin <= 2, which enter ONLY the b=0 side.  They are given their endpoint cost W (no A2' suppression).
"""
import sys, os, time
import numpy as np
from fractions import Fraction as Fr
from math import comb
from scipy.optimize import linprog
from flint import arb, arb_mat, ctx
from common import ladder_exact, W, allowed

ctx.prec = 256
MAXIT = int(os.environ.get("MAXIT", 60)); BATCH = int(os.environ.get("BATCH", 1500))


def Mmat(N):
    M = [[Fr(0)] * (N + 1) for _ in range(N + 1)]
    for n in range(N + 1):
        t = Fr(1)
        for k in range(N + 1):
            M[n][k] = t
            t = t * (-n + k) * (n + 1 + k) / ((k + 1) * (2 + k))
    return M


def candidates(g, light=False, odd_b0=False):
    Js = sorted(set(list(range(0, 12 * g + 40)) + [int(x) for x in np.geomspace(12 * g + 40, 20000, 60)]))
    ts = [Fr(0)] + [Fr(x).limit_denominator(64) for x in np.geomspace(1e-3, 1e5, 70) * g]
    ts = sorted(set(ts))
    out = []
    for J in Js:
        for b in (0, 4):
            if not allowed(b, J, odd_b0):
                continue
            for t in ts:
                out.append((b, J, Fr(J + g) + t))
    if light:
        for D in [Fr(x).limit_denominator(64) for x in np.linspace(0.5, g - 1e-3, 200)]:
            out.append((0, 0, D))
        for D in [Fr(x).limit_denominator(64) for x in np.linspace(3.001, g + 2 - 1e-3, 200)]:
            out.append((0, 2, D))
    return out


def rhat_arb(b, J, D, N):
    C = arb(J * (J + 1)); T = (arb(D.numerator) / D.denominator - 1) * (arb(D.numerator) / D.denominator - 2)
    r = [arb(1)]
    for k in range(N):
        s = b + k
        if k + 1 > J - b:
            r.append(arb(0)); continue
        y = s * (s + 1)
        r.append(r[-1] * (C - y) / (T - y))
    return r


def basis_rows(S, N, MT):
    rows = [rhat_arb(b, J, D, N) for (b, J, D) in S]
    R = arb_mat(rows) * MT                       # (m x N+1) * (N+1 x N+1)
    out = np.empty((len(S), N + 1))
    for i, (b, J, D) in enumerate(S):
        sg = 1.0 if b == 4 else -1.0
        w = float(W(b, D))
        for n in range(N + 1):
            out[i, n] = sg * float(R[i, n].mid()) / w
    return out


def run(g, N, light=False, odd_b0=False):
    t0 = time.time()
    M = Mmat(N)
    MT = arb_mat([[arb(M[n][k].numerator) / M[n][k].denominator for n in range(N + 1)] for k in range(N + 1)])
    S = candidates(g, light, odd_b0)
    B = basis_rows(S, N, MT)
    # mu_k binom^2 = c_k = sum_n beta_n M_nk ; constraint c_2 = 100  (mu_2 = 1)
    a_eq = np.array([float(M[n][2]) for n in range(N + 1)])
    idx = [i for i, s in enumerate(S) if s[1] <= 3 * g + 10 and float(s[2]) - s[1] <= 30 * g]
    for it in range(MAXIT):
        A = np.hstack([B[idx], -np.ones((len(idx), 1))])
        c = np.zeros(N + 2); c[-1] = 1.0
        res = linprog(c, A_ub=A, b_ub=np.zeros(len(idx)), A_eq=np.append(a_eq, 0.0)[None, :], b_eq=[100.0],
                      bounds=[(None, None)] * (N + 1) + [(0, None)], method='highs')
        if res.status != 0:
            print('LP status', res.status, res.message); return None
        beta, lam = res.x[:-1], res.x[-1]
        viol = B @ beta - lam
        bad = np.where(viol > 1e-7 * lam)[0]
        if len(bad) == 0:
            break
        idx = sorted(set(idx) | set(bad[np.argsort(-viol[bad])][:BATCH]))
    viol = B @ beta - lam
    order = np.argsort(-viol)
    # exact re-check on the 6 most binding states, with mu from beta rounded to rationals
    cks = [sum(Fr(float(beta[n])) * M[n][k] for n in range(N + 1)) for k in range(N + 1)]
    mus = [cks[k] / comb(k + 3, 3) ** 2 for k in range(N + 1)]
    worst = -1e9
    for i in order[:6]:
        b, J, D = S[i]
        r = ladder_exact(b, J, D, N)
        f = sum(m * x for m, x in zip(mus, r))
        worst = max(worst, float((1 if b == 4 else -1) * f / W(b, D)) / lam - 1)
    return dict(c2=3 * lam, it=it + 1, worst=worst, binding=[S[i] for i in order[:8]], t=time.time() - t0,
                nS=len(S), mu2=float(mus[2]))


if __name__ == '__main__':
    g = int(sys.argv[1]); N = int(sys.argv[2]); fl = sys.argv[3:]
    from v2_lp import three_state
    c3, _ = three_state(g)
    o = run(g, N, light='light' in fl, odd_b0='odd_b0' in fl)
    if o:
        print(f"g={g} N={N} {fl}: LP |c|^2={o['c2']:.8e}  3-state={float(c3):.8e}  ratio={o['c2']/float(c3):.6f}  "
              f"tau g^3={4*o['c2']**0.5*g**3:.2f} (3-state {4*float(c3)**0.5*g**3:.2f})  iters={o['it']}  "
              f"exact rel. viol={o['worst']:.1e}  states={o['nS']}  {o['t']:.0f}s")
        print('   binding:', [(b, J, float(D - J - g)) for (b, J, D) in o['binding']])

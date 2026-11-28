"""v2: independent LP for  max |c|^2  subject to
     endpoint rule   sum_i W_i rho_i <= 3            (D = 0, light states dropped: see v4 for light columns)
     R_k, k=0..N     sum_i sigma_i rho_i r_k(i) = delta_{k2} |c|^2   (sigma=+1 for b=4, -1 for b=0)
     rho_i >= 0 over heavy states i=(b,J,Delta), Delta - J >= g.
Solved through its dual (a minimax / Chebyshev problem):
     |c|^2_max = 3 * min_{mu, mu_2=1} max_i sigma_i sum_k mu_k r_k(i) / W_i
with constraint generation over a very large state set (J up to 20000, twist up to 1e5 g), and an exact
rational re-check of the dual certificate on every generated state.  Also prints the exact 3-state primal.

usage: python3 v2_lp.py g N [extra flags]
   flags: odd_b0  -> also allow b=0 ladder for odd J (sensitivity test)
          dimgap  -> impose Delta >= g (+J only via unitarity) instead of twist >= g (sensitivity test)
"""
import sys, time
import numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from common import ladder_float, ladder_exact, W, allowed


def three_state(g):
    g = Fr(g)
    A, B, C = (0, 2, g + 2), (4, 4, g + 4), (4, 6, g + 6)
    rA, rC = ladder_exact(*A, 3), ladder_exact(*C, 3)
    rhoC = Fr(1); rhoA = rhoC * rC[1] / rA[1]; rhoB = rhoA - rhoC
    bud = rhoA * W(0, A[2]) + rhoB * W(4, B[2]) + rhoC * W(4, C[2])
    s = Fr(3) / bud
    rhoA, rhoB, rhoC = rhoA * s, rhoB * s, rhoC * s
    # check every rule exactly
    rB = ladder_exact(*B, 3)
    for k in range(4):
        val = rhoB * rB[k] + rhoC * rC[k] - rhoA * rA[k]
        assert (val == 0) == (k != 2), (k, val)
    c2 = rhoC * rC[2] + rhoB * rB[2] - rhoA * rA[2]
    return c2, (rhoA, rhoB, rhoC)


def candidate_states(g, odd_b0=False, dimgap=False, dense=False):
    Js = sorted(set(list(range(0, 12 * g + 40)) + [int(x) for x in np.geomspace(12 * g + 40, 20000, 150 if dense else 60)]))
    ts = np.concatenate([[0.0], np.geomspace(1e-3, 1e5, 160 if dense else 70) * g])
    out = []
    for J in Js:
        for b in (0, 4):
            if not allowed(b, J, odd_b0):
                continue
            for t in ts:
                D = (max(g, J + 1.0) + t) if dimgap else (J + g + t)
                out.append((b, J, float(D)))
    return out


def rows(S, N):
    R = np.empty((len(S), N + 1))
    for i, (b, J, D) in enumerate(S):
        R[i] = (1 if b == 4 else -1) * ladder_float(b, J, D, N) / W(b, D)
    return R


def solve_dual(R, N):
    # variables x = (mu_0..mu_N except mu_2 fixed=1 ; lam).  minimise lam
    # constraint per state: R_i . mu - lam <= 0  ->  sum_{k!=2} R_ik mu_k - lam <= -R_i2
    ks = [k for k in range(N + 1) if k != 2]
    sc = np.abs(R[:, ks]).max(axis=0); sc[sc == 0] = 1.0
    A = np.hstack([R[:, ks] / sc, -np.ones((R.shape[0], 1))])
    lam_sc = np.abs(R[:, 2]).max()
    A[:, -1] *= lam_sc
    c = np.zeros(len(ks) + 1); c[-1] = 1.0
    res = linprog(c, A_ub=A, b_ub=-R[:, 2], bounds=[(None, None)] * len(ks) + [(0, None)], method='highs')
    if res.status != 0:
        return None, None, res
    mu = np.zeros(N + 1); mu[2] = 1.0
    mu[ks] = res.x[:-1] / sc
    return mu, res.x[-1] * lam_sc, res


def exact_margin(mu, lam, st, N):
    b, J, D = st
    Dq = Fr(D).limit_denominator(10 ** 9) if not float(D).is_integer() else Fr(int(D))
    r = ladder_exact(b, J, Dq, N)
    f = sum(Fr(float(m)) * x for m, x in zip(mu, r))
    return float((1 if b == 4 else -1) * f / W(b, Dq)), float(Fr(float(lam)))


def run(g, N, odd_b0=False, dimgap=False):
    t0 = time.time()
    allS = candidate_states(g, odd_b0, dimgap)
    Rall = rows(allS, N)
    # start: low-J states only
    idx = [i for i, s in enumerate(allS) if s[1] <= 3 * g + 10 and s[2] - s[1] <= 30 * g]
    for it in range(30):
        mu, lam, res = solve_dual(Rall[idx], N)
        if mu is None:
            print('  LP failed', res.message); return None
        viol = Rall @ mu - lam
        bad = np.where(viol > 1e-9 * max(lam, 1e-300))[0]
        if len(bad) == 0:
            break
        new = list(bad[np.argsort(-viol[bad])][:400])
        idx = sorted(set(idx) | set(new))
    # binding states and exact check of the worst constraints
    viol = Rall @ mu - lam
    order = np.argsort(-viol)[:8]
    c2 = 3 * lam
    tau = 4 * np.sqrt(c2)
    ex = [exact_margin(mu, lam, allS[i], N) for i in order[:5]]
    worst_exact = max(e[0] - e[1] for e in ex)
    return dict(g=g, N=N, c2=c2, tau=tau, taug3=tau * g ** 3, c2g6=c2 * g ** 6, iters=it + 1,
                nstates=len(allS), binding=[allS[i] for i in order[:6]], worst_exact_rel=worst_exact / lam,
                time=time.time() - t0, mu=mu)


if __name__ == '__main__':
    g = int(sys.argv[1]); N = int(sys.argv[2])
    flags = sys.argv[3:]
    c2, rho = three_state(g)
    print(f"g={g}: exact 3-state primal |c|^2 = {float(c2):.8e}  |c|^2 g^6 = {float(c2)*g**6:.3f}  "
          f"tau g^3 = {4*float(c2)**0.5*g**3:.3f}")
    out = run(g, N, odd_b0='odd_b0' in flags, dimgap='dimgap' in flags)
    if out:
        print(f"  dual LP (N={N}, flags={flags}): |c|^2 = {out['c2']:.8e}  |c|^2 g^6 = {out['c2g6']:.3f}  "
              f"tau g^3 = {out['taug3']:.3f}  ratio to 3-state = {out['c2']/float(c2):.6f}  "
              f"iters={out['iters']}  exact-recheck worst rel. violation={out['worst_exact_rel']:.1e}  ({out['time']:.0f}s)")
        print('  most binding states:', [(b, J, round(D - J - g, 3)) for (b, J, D) in out['binding']])

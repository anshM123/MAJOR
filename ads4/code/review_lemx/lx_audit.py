"""lx_audit.py: step-by-step audit of the Lemma X proof (ANALYTIC_SLIVER §5.5, rt_pp.cross) along exact trajectories
(double precision; values stay < e^40 on all audited states, and the final values are cross-checked against Arb in lx_exact.py).

Indexing: step k maps (S_k, S_{k-1}) -> (S_{k+1}, S_k) with w_k = eps_k alpha/2 = u_k + i y_k.  delta = c tau^{-2/3}, yA/B = 1 -/+ delta,
dy = 2 tau/(a+1)^2 (proved bound on y_{k+1} - y_k), eta = 2(delta + u_max).
Phase I indices k <= K (y_k <= yA); c_k = P(w_k)^{-1}(S_k, S_{k-1}).  Window: steps k = K .. while y_k < yB, T_j = i^{-j} S_j,
N_j = |T_j| + |T_j - T_{j-1}|/sqrt(eta).  Phase III: indices with y_k >= yB, eigen-coordinates again.
Checks (slack = rhs - lhs in log or absolute units; negative = violated):
  P1  log[|lam+(k)|(1 + (|dlam+| + |dlam-|)/(2|r_{k+1}|)) |c_k|] - log|c_{k+1}|                 (phases I and III)
  P2  u_k/sqrt(1 - y_k^2) - log|lam+(k)|                                                         (phase I)
  P3  [log c0 + sigma log(2/(1+sqrt(1-Y^2))) + q(asin Y + q/2 log(1/(1-Y^2)))] - log(d_k^sigma |c_k|),  Y = y_k + dy (phase I)
  E1  log[(2 + e1 + e1/sqrt(eta)) |c_K|] - log N_K  and  e1 - (|lam+ - i| + |lam- - i|)  at index K
  W1  log(step factor) - log(N_{k+1}/N_k)  (first step from y_K <= yA uses eta' = eta + 2 dy)
  W2  tau(1/yA - 1/yB) + 2 - #window steps
  E2  log[(sqrt(eta)(1+e2) + e2)/(2 rB) N] - log|c|  at the phase-III entry
  P4  acosh y + u/sqrt(y^2 - 1 - u^2) - log|lam+|                                                (phase III)
  FIN log cross(sigma, tau, tau, c) - log|S^_n|  (rt_pp.cross, Arb, uniform a >= 998)
"""
import math, cmath, sys
from lx_common import d_of, n_at, LAMf
sys.path.insert(0, '/home/user/MAJOR/ads4/code/analytic')
import rt_pp


def lam(w):
    r = cmath.sqrt(1 + w*w)
    lp, lm = w + r, w - r
    if abs(lp) < abs(lm):
        lp, lm, r = lm, lp, -r
    return lp, lm, r


def traj(n, a, al):
    N = n + a
    S = [1+0j, 1+0j]          # S[j+1] = S_j, j = -1..n
    for k in range(n):
        eps = 2*(2*k + 1)/((N - k)*(N + k + 1))
        S.append(S[-2] + eps*al*S[-1])
    return S


def audit(n, a, sig, tau, c, S=None):
    N = n + a; al = complex(sig, tau); q = abs(al)/tau
    if S is None:
        S = traj(n, a, al)
    Sx = lambda j: S[j + 1]
    T = lambda j: (1j**((-j) % 4))*Sx(j)
    dl = c*tau**(-2/3); yA, yB = 1 - dl, 1 + dl
    dy = 2*tau/(a + 1)**2
    umax = sig/(a + 1)
    eta = 2*(dl + umax); se = math.sqrt(eta)
    sl = {}
    upd = lambda key, v: sl.__setitem__(key, min(sl.get(key, 1e9), v))
    eps = [2*(2*k + 1)/((N - k)*(N + k + 1)) for k in range(n)]
    ys = [e*tau/2 for e in eps]
    us = [e*sig/2 for e in eps]
    ws = [e*al/2 for e in eps]
    L = [lam(w) for w in ws]
    def coords(k):
        lp, lm, r = L[k]
        c1 = (Sx(k) - lm*Sx(k - 1))/(2*r); c2 = (-Sx(k) + lp*Sx(k - 1))/(2*r)
        return max(abs(c1), abs(c2))
    logd = lambda k: math.log(1 - k*(k + 1)/(N*(N + 1)))
    c0 = 1 + 2*abs(al)/a**2
    # phase I
    K = 0
    while K + 1 < n and ys[K + 1] <= yA:
        K += 1
    cn = [None]*(n)
    for k in range(0, K + 1):
        cn[k] = coords(k)
        lp, lm, r = L[k]
        upd('P2', us[k]/math.sqrt(1 - ys[k]**2) - math.log(abs(lp)))
        Y = min(ys[k] + dy, 1 - 1e-12)
        B = math.log(c0) + sig*math.log(2/(1 + math.sqrt(1 - Y*Y))) + q*(math.asin(Y) + q/2*math.log(1/(1 - Y*Y)))
        upd('P3', B - (sig*logd(k) + math.log(cn[k])))
        if k >= 1:
            lp0, lm0, r0 = L[k - 1]
            fac = abs(lp0)*(1 + (abs(lp - lp0) + abs(lm - lm0))/(2*abs(r)))
            upd('P1', math.log(fac*cn[k - 1]) - math.log(cn[k]))
    nwin = 0
    if K + 1 < n or ys[K] > yA:
        # window entry at index K
        e1 = 2*(umax + 1 - yA + dy + math.sqrt(1 - (yA - dy)**2 + umax**2 + 2*umax))
        lp, lm, r = L[K]
        upd('E1', e1 - (abs(lp - 1j) + abs(lm - 1j)))
        Nj = lambda j: abs(T(j)) + abs(T(j) - T(j - 1))/se
        upd('E1', math.log((2 + e1 + e1/se)*cn[K]) - math.log(Nj(K)))
        k = K
        etap = eta + 2*dy
        while k < n and ys[k] < yB:
            fac = max(1 + etap + etap/se, 1 + se) if k == K else 1 + se + eta
            upd('W1', math.log(fac) - math.log(Nj(k + 1)/Nj(k)))
            nwin += 1; k += 1
        upd('W2', tau*(1/yA - 1/yB) + 2 - nwin)
        if k < n:
            # phase III entry at index k (y_k >= yB)
            rB = math.sqrt(yB*yB - 1 - umax**2)
            e2 = umax + yB + dy - 1 + math.sqrt((yB + dy)**2 - 1 + umax**2 + 2*umax*(yB + dy))
            c3 = coords(k)
            upd('E2', math.log((se*(1 + e2) + e2)/(2*rB)*Nj(k)) - math.log(c3))
            prev = c3
            for j in range(k + 1, n):
                cj = coords(j)
                lp0, lm0, r0 = L[j - 1]; lp, lm, r = L[j]
                fac = abs(lp0)*(1 + (abs(lp - lp0) + abs(lm - lm0))/(2*abs(r)))
                upd('P1', math.log(fac*prev) - math.log(cj))
                upd('P4', math.acosh(ys[j]) + us[j]/math.sqrt(ys[j]**2 - 1 - us[j]**2) - math.log(abs(lp)))
                prev = cj
    fin = math.log(abs(Sx(n))) + sig*math.log(d_of(n, a))
    return fin, sl, nwin


if __name__ == '__main__':
    Lu = LAMf
    cases = []
    for a in (998.0, 998.49, 1000.0, 1010.0):
        for Lc in (3.65, 3.75):
            n = n_at(a, math.exp(-Lc))
            for sig in (0.0, 2.5, 5.5):
                for tau in (Lu - 1, Lu):
                    cases.append((n, a, sig, tau))
    # general lemma (phase III): tau in (a, 1150]
    for (n, a, sig, tau) in [(80000, 998.0, 2.5, 1031.0), (80000, 998.0, 0.25, 1150.0), (300000, 998.0, 5.5, 1150.0),
                             (20000, 998.0, 0.25, 1080.0), (3000, 998.0, 5.5, 1000.0), (80000, 1050.0, 2.5, 1100.0)]:
        cases.append((n, a, sig, tau))
    worst = {}
    for (n, a, sig, tau) in cases:
        S = traj(n, a, complex(sig, tau))
        for c in (0.7, 1.0, 1.5):
            fin, sl, nw = audit(n, a, sig, tau, c, S)
            bnd = math.log(rt_pp.cross(sig, tau, tau, c))
            sl['FIN'] = bnd - fin
            for kk, v in sl.items():
                worst[kk] = min(worst.get(kk, 1e9), v)
            print(f'n={n} a={a} sig={sig} tau={tau:.3f} c={c}: log|S^|={fin:.3f} cross={bnd:.3f} nwin={nw} '
                  + ' '.join(f'{k}:{v:.3g}' for k, v in sl.items()), flush=True)
    print('WORST SLACKS (negative = violated):', {k: round(v, 5) for k, v in worst.items()})

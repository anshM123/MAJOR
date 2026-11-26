"""
loc5.py: Lemma Loc''' -- the j-truncation remainder by a 2D (u = 1 - B, theta) cell decomposition of the Beta x Cauchy-circle measure,
with Chernoff (optimal-moment) B-weights and phase-aware bounds of the branch factors.  ANALYTIC_SLIVER3 §2.  [P + certified arithmetic]

Representation (TAIL_CERT §4, re-derived in REVIEW_SLIVER_L5 §2):  for G analytic in |z| < 1,
    L(G) = E_B[ (2 pi i)^{-1} oint_{|z| = rho} (1 - z)^{-a-1} G(B z) z^{-n-1} dz ] / C(N, a),   B ~ Beta(n + s, a),  rho = n/N,
hence |L(G)| <= K0 E_B E^W_theta |G(B rho e^{i theta})|,  E^W = normalised weight W(theta) = |1 - rho e^{i theta}|^{-a-1},
K0 = rho^{-n} I_W / C(N, a),  I_W = (2 pi)^{-1} int W.
  (K0)  K0 <= (pi/2) sqrt(1 + 2/a) e^{1/(6a)} <= 1.5727 for every state (a >= 998): sin(th/2) >= th/pi gives
        I_W <= t^{-a} (2 sqrt rho)^{-1} sqrt(pi) Gamma(a/2)/(2 Gamma((a+1)/2)), Robbins gives C(N,a) rho^n >= t^{-a} e^{-1/(6a)}/sqrt(2 pi a rho),
        Gautschi gives Gamma(a/2)/Gamma((a+1)/2) < 2 sqrt(a/2 + 1)/a.   (t := 1 - rho = a/N.)
  (TH)  P_W(|theta| >= th0) <= q(th0)^{-a-1}/I_W <= sqrt(2 pi (a+1)) t^{-1} (t/q(th0))^{a+1}/erf(X sqrt((a+1)/2)),
        q^2 = t^2 + 4 rho sin^2(th/2) (increasing on [0, pi]),  X = pi sqrt(rho)/t   (I_W >= (1/pi) t^{-a} rho^{-1/2} int_0^X e^{-(a+1)x^2/2} dx).
  (CH)  P(1 - B >= u) <= (a)_m / ((N + s)_m u^m) for every m (Markov, 1 - B ~ Beta(a, n + s)); with (a + i)/(N + s + i) <= (1 + x)/(K + x),
        x = i/a, K = N/a >= 1/t_hi, the summand is increasing in x, so  log P <= a int_0^{m/a} f,  f(x) = log((1 + x)/((K + x) u)),
        and with m = floor(a X0):  log P <= 998 max(F(X0 - 1/998), F(X0)) whenever both are <= 0, F(X) = int_0^X f (closed form below).
  (DR)  |V_b G(xi)| = |(1 + tau D_b) G(xi)| <= (1 + 2.6e-5) sup_{|zeta - xi| <= t/2} |G(zeta)|  (Cauchy on the circle of radius t/2;
        D_b = z(1-z)d^2/2 + [(7/2 - 2b)z + 1/2]d + 9/2 - b(b-1), tau = 1/T <= 1/N^2: tau/r^2 <= 4/a^2, tau/r <= 2/a^2, tau <= 1/a^2).
  The probability of a cell {u in [u0, u1), |theta| in [th0, th1)} is P_B * P_W <= P(1-B >= u0) * P_W(|theta| >= th0) (B and theta independent).
Bounds of |R_J(zeta)| on the delta-enlarged cell (zeta = (1 - u) rho e^{i theta} + disc(t_hi/2)), lo <= |1 - zeta| <= hi, |arg(1 - zeta)| <= A,
|Im eta| <= E (eta = 2 atanh sqrt(zeta), cosh eta = (1 + zeta)/(1 - zeta)), all Arb enclosures over the cell:
  (T)   hi < 1:  |R| <= 2 Bmax e^{nu_hi A} |c_J| hi^{J - e0} / (1 - R_t hi)          (tail series, |(1-zeta)^{-+ i nu}| = e^{+- nu arg})
  (FT)  |R| <= (pi/2) e^{nu_hi E} lo^{-p} + 2 Bmax e^{nu_hi A} sum_{j<J} |c_j| max(lo, hi)^{j - e0}   (FULLPOS3 Lemma B; |Im eta| <= pi/2)
  R_t >= sup_{j >= J} |c_{j+1}/c_j|: min(1 + nu_hi^2/(J+1)^2, 1 [if the Loc'' tail-ratio lemma applies]).
Small nu in [0, 0.4]: the B-pole handled as in Loc' (rt_loc3): sum_{+-} B Phi = (2 i nu)^{-1} int (E Phi)', bounded by (1/r') sup_{Wbox} |E Phi|.
"""
import math, os, sys, json
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.abspath(os.path.join(HERE, '..', 'analytic'))
AN2 = os.path.abspath(os.path.join(HERE, '..', 'analytic2'))
sys.path.insert(0, AN); sys.path.insert(0, AN2)
I = acb(0, 1)
A_MIN = 998
EPS_LOG = 1e-6          # safety margin on float log-accumulation


def up(x):
    return float(arb(x).upper())


def lo_(x):
    return float(arb(x).lower())


def arb_iv(a, b):
    a = arb(a); b = arb(b)
    return arb((a + b)/2, ((b - a)/2).upper()*(1 + 1e-12) + 1e-300)


def t_range(La, Lb, a_lo=A_MIN):
    """[t_lo, t_hi] containing t = a/N for every state (s in {1, 9}, a >= a_lo) with -log d in [La, Lb]:
    t^2 - (2 + e) t + d (1 + e) = 0 (d = a(2N - a + s)/(N(N + s)), e = s/N), t increasing in d and in e; e <= 9 t/a_lo <= 9 * 0.7/a_lo."""
    d_lo = arb(-Lb).exp(); d_hi = arb(-La).exp()
    def t_of(d, e):
        return ((2 + e) - ((2 + e)**2 - 4*d*(1 + e)).sqrt())/2
    tl = t_of(d_lo, arb(0)); th = t_of(d_hi, arb(9*0.7)/a_lo)
    assert th < 0.7
    return lo_(tl)*(1 - 1e-12), up(th)*(1 + 1e-12)


def K0_bound(a=A_MIN):
    aa = arb(a)
    return up(arb.pi()/2*(1 + 2/aa).sqrt()*(1/(6*aa)).exp())


def lchern(u0, t_hi, a=A_MIN):
    """log upper bound of P(1 - B >= u0), uniform over states with t <= t_hi and a >= a (Arb)."""
    if u0 <= t_hi*1.0001:
        return 0.0
    u = arb(u0); K = 1/arb(t_hi)
    Xs = float(((K*u - 1)/(1 - u)).mid())        # zero of f
    def F(X):
        X = arb(X)
        return (1 + X)*(1 + X).log() - (K + X)*(K + X).log() + K*K.log() - X*u.log()
    X0 = Xs
    X1 = max(0.0, X0 - 1.0/a)
    v = max(up(F(X1)), up(F(X0)))
    if v > 0:
        return 0.0
    return v*a


def lfrac(th0, t_lo, t_hi, a=A_MIN):
    """log upper bound of P_W(|theta| >= th0), uniform over t in [t_lo, t_hi], a >= 998 (Arb)."""
    if th0 <= 0:
        return 0.0
    th = arb(th0); tt = arb(t_hi)
    q = (tt*tt + 4*(1 - tt)*(th/2).sin()**2).sqrt()
    r = tt/q
    if up(r) > 0.999:
        return 0.0
    aa = arb(a)
    X = arb.pi()*(1 - tt).sqrt()/tt
    er = (X*((aa + 1)/2).sqrt()).erf()
    v = (2*arb.pi()*(aa + 1)).sqrt()/arb(t_lo)*r**(aa + 1)/er
    return min(0.0, float(v.log().upper()))


def cell_geom(u0, u1, th0, th1, t_lo, t_hi, dress):
    """Arb enclosures over the enlarged cell: (lo, hi, A, E) with lo <= |1-zeta| <= hi, |arg(1-zeta)| <= A, |Im eta| <= E."""
    r = arb_iv((1 - u1)*(1 - t_hi), (1 - u0)*(1 - t_lo))
    th = arb_iv(th0, th1)
    xi = r*(I*th).exp()
    dl = t_hi/2 if dress else 0.0
    zeta = xi + acb(arb(0, dl), arb(0, dl)) if dress else xi
    w = 1 - zeta
    aw = abs(w)
    lo = lo_(aw) - (dl*1e-9 if dress else 0.0)
    hi = up(aw)
    if lo <= 0:
        lo = 1e-300
    A = up(abs(w.arg()))
    A = min(A, math.pi/2)
    # Im eta: eta = 2 atanh(sqrt zeta); valid bound |Im eta| <= pi/2 always (|zeta| < 1)
    try:
        s = zeta.sqrt()
        eta = 2*s.atanh()
        E = up(abs(eta.imag))
        if not math.isfinite(E):
            E = math.pi/2
    except Exception:
        E = math.pi/2
    E = min(E, math.pi/2)
    return lo, hi, A, E


def geometry(La, Lb, dress, ratio=1.25, fine_center=True):
    """list of (lw, lo, hi, A, E) over a cell partition of [0,1] x [0,pi] in (u, |theta|)."""
    t_lo, t_hi = t_range(La, Lb)
    def breaks(top, base):
        b = [0.0]
        x = base
        while x < top:
            b.append(x); x *= ratio
        b.append(top)
        return b
    base = t_hi/4 if fine_center else t_hi
    ub = breaks(1.0, base); tb = breaks(math.pi, base)
    out = []
    for i in range(len(ub) - 1):
        wB = lchern(ub[i], t_hi)
        for k in range(len(tb) - 1):
            wT = lfrac(tb[k], t_lo, t_hi)
            lw = wB + wT
            if lw < -20000:
                lo, hi, A, E = 1e-300, 2.0 + t_hi, math.pi/2, math.pi/2    # weight negligible: crude geometry
                lw_eff = lw
            else:
                lo, hi, A, E = cell_geom(ub[i], ub[i + 1], tb[k], tb[k + 1], t_lo, t_hi, dress)
                lw_eff = lw
            out.append((lw_eff, max(lo, (0.49*0 + 1e-300)), hi, A, E, ub[i], ub[i + 1], tb[k], tb[k + 1]))
    return dict(t_lo=t_lo, t_hi=t_hi, cells=out)


# ---------------------------------------------------------------- coefficients
def lcs_real(nu_lo, nu_hi, J):
    """log upper bounds of |c_j(+-nu)|, real nu in [nu_lo, nu_hi], j = 0..J (rt_loc3.lcs, Arb-rounded)."""
    import rt_loc3 as R3
    return np.array(R3.lcs(nu_lo, nu_hi, 0.0, J + 1))


def lcs_cplx(x_hi, y, J):
    import rt_loc3 as R3
    return np.array(R3.lcs(0.0, x_hi, y, J + 1))


def logsumexp(v):
    m = np.max(v)
    return m + math.log(np.sum(np.exp(v - m)))


def lB_real(nu_lo):
    import rt_cert as RC
    return math.log(2*RC.B_abs_real(nu_lo))


def ratio_tail(nu_lo, nu_hi, J):
    """R_t >= sup_{j >= J} |c_{j+1}/c_j| (real nu):  ((j+1/2)^2 + nu^2)/((j+1) sqrt((j+1)^2 + 4 nu^2)) <= 1 + nu^2/(j+1)^2;
    and <= 1 under the Loc'' condition q = 4 nu_lo^2 - 2 nu_hi^2 - 3/2 > 0, J^2 q >= (1/4 + nu_hi^2)^2."""
    r = 1 + nu_hi**2/(J + 1)**2
    q = 4*nu_lo*nu_lo - 2*nu_hi*nu_hi - 1.5
    if q > 0 and J*J*q >= (0.25 + nu_hi*nu_hi)**2:
        r = min(r, 1.0)
    return r*(1 + 1e-12)


def bound_from_geom(G, p, J, lc, lpref, nu_hi, Rt, mode='real', yfac=None):
    """log bound of |V_b L(R_J)| (or |L(R_J)| if undressed; the dressing factor is added by the caller) for one nu-cell.
    lc: log|c_j| (j = 0..J), lpref: log of the prefactor (2 Bmax, or sup|E|/r' in small mode), nu_hi: exponent rate for e^{nu A},
    yfac(lo, hi): extra log factor for complex w (small mode)."""
    e0 = p - 0.5
    cells = G['cells']
    lw = np.array([c[0] for c in cells]); lo = np.array([c[1] for c in cells]); hi = np.array([c[2] for c in cells])
    A = np.array([c[3] for c in cells]); E = np.array([c[4] for c in cells])
    ext = np.zeros(len(cells)) if yfac is None else np.array([yfac(l, h) for l, h in zip(lo, hi)])
    # tail route
    tail = np.full(len(cells), np.inf)
    ok = hi*Rt < 0.999
    tail[ok] = lpref + nu_hi*A[ok] + ext[ok] + lc[J] + (J - e0)*np.log(hi[ok]) - np.log(1 - Rt*hi[ok])
    # F + T route
    lF = math.log(math.pi/2) + nu_hi*E - p*np.log(lo)
    js = np.arange(J)
    lcJ = lc[:J]
    # sum_{j<J} |c_j| max(lo,hi)^{j-e0}: j - e0 < 0 -> lo, else hi
    neg = js - e0 < 0
    lT = np.empty(len(cells))
    for k in range(len(cells)):
        v = np.where(neg, lcJ + (js - e0)*math.log(lo[k]), lcJ + (js - e0)*math.log(hi[k]))
        lT[k] = logsumexp(v)
    lT = lpref + nu_hi*A + ext + lT
    ft = np.logaddexp(lF, lT)
    best = np.minimum(tail, ft)
    tot = lw + best
    return logsumexp(tot) + EPS_LOG, int(np.argmax(tot))


def cell_bound(G, p, J, nu_lo, nu_hi, dress, K0=None):
    """log bound of |V_b L(R_J)(nu)| for all nu in [nu_lo, nu_hi] (real, nu_lo >= 0.4), or nu in [0, 0.4] (small mode)."""
    K0 = K0 or K0_bound()
    ldr = math.log(1 + 2.6e-5) if dress else 0.0
    if nu_hi <= 0.4 + 1e-12:
        import rt_loc3 as R3
        rp = R3.RPRIME
        lc = lcs_cplx(0.4 + rp, rp, J)
        lpref = math.log(R3.Esup_small()/rp)
        def yfac(l, h):
            return rp*max(-math.log(l), math.log(h), 0.0)
        # F part: |F| <= (pi/2) e^{0.4 E} lo^{-p} (no division); the T/tail parts carry the pole bound (1/r') sup|E Phi|
        v, k = bound_from_geom(G, p, J, lc, lpref, 0.4 + rp, 1 + (0.45**2 + 0.1)/(J + 1)**2*4, 'small', yfac)
        return math.log(K0) + ldr + v, k
    lc = lcs_real(nu_lo, nu_hi, J)
    lpref = lB_real(nu_lo)
    Rt = ratio_tail(nu_lo, nu_hi, J)
    v, k = bound_from_geom(G, p, J, lc, lpref, nu_hi, Rt)
    return math.log(K0) + ldr + v, k


if __name__ == '__main__':
    import time
    La, Lb, J = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    t0 = time.time()
    G = geometry(La, Lb, True)
    print('geometry', len(G['cells']), 'cells', G['t_lo'], G['t_hi'], f'{time.time()-t0:.1f}s', flush=True)
    for nu in (0.0, 1.0, 300.0, 600.0, 900.0, 990.0, 1029.0, 1100.0, 1149.0):
        if nu == 0.0:
            b = cell_bound(G, 3, J, 0.0, 0.4, True)
        else:
            b = cell_bound(G, 3, J, nu, nu + 1, True)
        c = G['cells'][b[1]]
        print(f'nu={nu}: log|V R_J| <= {b[0]:.1f}   worst cell u [{c[5]:.4f},{c[6]:.4f}] th [{c[7]:.4f},{c[8]:.4f}] lw {c[0]:.1f} lo {c[1]:.3g} hi {c[2]:.3g} A {c[3]:.3f} E {c[4]:.3f}  [{time.time()-t0:.1f}s]', flush=True)

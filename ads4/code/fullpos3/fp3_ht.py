"""
fullpos3/fp3_ht.py (FULLPOS3.md §4): certified high-twist bound.  For a channel b in {0,4}, a beta-cell [ba, bb] (c = tanh^2(beta/2))
and a twist threshold a0, decide (Arb) whether
     (1-c)^3 * Hhat > 0   for every heavy state with c in the cell and a = N - n >= a0,
where Hhat = H/w^3 = L(G_A + tau D_b G_A + tau sigma_b G_D)  (FULLPOS3 §4.1).  Inequality certified:
     F_cert(cell) > (1-c)^3 [ |G^r(c)| + tau|D_b G_A(c)| + sum_{m=1}^{M} h^m ||P_m||_R S(R) + Far ]          (FULLPOS3 §4.3)
Usage: fp3_ht.py b ba bb a0 [M] [m1]
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, acb, arb_poly, fmpq, fmpq_poly, ctx
import fp3_htdata as D
from fp3_htdata import A, LAM, PI, up, lo

ctx.prec = 256
PSI2 = D.Psi(2); PSI3 = D.Psi(3)
NARC = 96

# ---------------------------------------------------------------- J tables (decreasing in beta): grid
BGRID = [Fr(1, 200)*Fr(103, 100)**i for i in range(0, 258)]         # 0.005 .. ~9.9
_JT = {}
def jt(k):
    if k not in _JT:
        _JT[k] = D.jtable(PSI2 if k == 2 else PSI3, BGRID, ng=400)
    return _JT[k]

def jval(k, beta_lo):
    """J_s, J_C, J_e at the largest grid point <= beta_lo (upper bounds, decreasing in beta)."""
    ok = [i for i, g in enumerate(BGRID) if A(g) <= beta_lo]
    if not ok:
        inf = arb('inf'); return (inf, inf, inf)
    return jt(k)[max(ok)]

_RI = {}
def rint(eta_hi):
    """Rint at eta rounded up to a grid of 1e-4."""
    e = math.ceil(float(eta_hi.upper())*1e4 + 1e-9)/1e4
    if e not in _RI:
        _RI[e] = D.Rint(arb(e))[0]
    return _RI[e]

_CC = {}
def ccor(k, eta_hi):
    e = math.ceil(float(eta_hi.upper())*1e3 + 1e-9)/1e3
    if (k, e) not in _CC:
        _CC[(k, e)] = (PSI2 if k == 2 else PSI3).Ccor(arb(e))
    return _CC[(k, e)]

# ---------------------------------------------------------------- sup bounds on circles z = c e^t, |t| = R
def _piece_geom(cl, Dl, R, kind, pa, pb):
    """ball enclosures on a piece of the stadium boundary (upper half), z = cl e^t:
       kind 0: left quarter circle t = R e^{i phi}, phi in [pa, pb] subset [pi/2, pi];
       kind 1: segment t = x + iR, x in [pa, pb] subset [0, Dl];
       kind 2: right quarter circle t = Dl + R e^{i phi}, phi in [pa, pb] subset [0, pi/2]."""
    pb_ = arb((pa + pb)/2, (pb - pa)/2)
    if kind == 1:
        t = acb(pb_, R)
    elif kind == 0:
        t = acb(R*pb_.cos(), R*pb_.sin())
    else:
        t = acb(Dl + R*pb_.cos(), R*pb_.sin())
    z = acb(cl)*t.exp()
    X = (1 + z)/(1 - z)
    xi = X.acosh()
    return z, xi


_RC = {}
def rcomp(eta_hi):
    e = math.ceil(float(eta_hi.upper())*1e4 + 1e-9)/1e4
    if e not in _RC:
        _RC[e] = [v for (_k, v) in D.rint_components(arb(e))]
    return _RC[e]


def circle_sups(cl, ch, R, b, tol=2e-4):
    """component-wise upper bounds of sup|.| over the union of the circles |log(z/c)| = R, c in [cl, ch] (boundary of the stadium,
    maximum principle; upper half by conjugation symmetry).  Components (each an analytic function, FULLPOS3 §4.3):
      'As': psi2 smooth (p=0 + coth) part of G_A;  'Ae': psi2 p=1 endpoint part;  'r<i>': residual cells of G^r;
      'Ds', 'De': the same split of G_D;  'zinv': max |1/z|."""
    out = {}
    def mx(key, v):
        if key not in out or v > out[key]:
            out[key] = v
    R = arb(R); cl = arb(cl); Dl = up((arb(ch)/cl).log())
    stack = [(0, PI/2, PI), (1, arb(0), Dl), (2, arb(0), PI/2)]
    while stack:
        kind, pa, pb = stack.pop()
        z, xi = _piece_geom(cl, Dl, R, kind, pa, pb)
        ok = xi.real.is_finite() and xi.imag.is_finite()
        if ok:
            beta = lo(xi.real); eta = up(abs(xi.imag))
            ok = beta > 0 and eta < PI/2 and float(xi.imag.rad()) < tol and float(xi.real.rad()) < 0.01*float(beta.lower())
        if not ok:
            if float((pb - pa).upper()) < 1e-12:
                raise RuntimeError('stadium geometry failed at R=%s' % R)
            m = (pa + pb)/2
            stack.append((kind, pa, m)); stack.append((kind, m, pb))
            continue
        axi = up(abs(xi)); om = lo(abs(1 - z))
        Js, Jc, Je = jval(2, beta)
        k2 = 2*om**(-3)*(2**0.5/PI)
        mx('As', up(k2*(Js + ccor(2, eta)*Jc)))
        mx('Ae', up(k2*(LAM*eta).exp()*Je))
        kr = up((2/LAM)*om**(-3)*(axi/beta)*D.P_mhalf(beta))
        for i, v in enumerate(rcomp(eta)):
            mx('r%d' % i, up(kr*v))
        Js3, Jc3, Je3 = jval(3, beta)
        k3 = LAM**2*om**(-4)*(2**0.5/PI)
        mx('Ds', up(k3*(Js3 + ccor(3, eta)*Jc3)))
        mx('De', up(k3*(LAM*eta).exp()*Je3))
        mx('zinv', up(1/abs(z)))
    return out


# ---------------------------------------------------------------- exact P_m (m <= m1) as arb_poly in k
def faulhaber(q):
    """F_q(k) = sum_{j<k} j^q as fmpq_poly."""
    co = [fmpq(0)]*(q + 2)
    for r in range(q + 1):
        co[q + 1 - r] += fmpq(math.comb(q + 1, r))*fmpq.bernoulli(r)/(q + 1)
    return fmpq_poly(co)

def S_poly(l, s):
    out = fmpq_poly([0])
    for i in range(l + 1):
        out += math.comb(l, i)*fmpq(s)**(l - i)*faulhaber(l + i)
    return out

def Pm_exact(cb, s, m1):
    a = []
    for l in range(1, m1 + 1):
        Sl = S_poly(l, s)
        co = [arb(Sl[j]) for j in range(Sl.degree() + 1)]
        a.append(arb_poly(co)*((cb**(-l) - 1)/l))
    P = [arb_poly([1])]
    for m in range(1, m1 + 1):
        acc = arb_poly([0])
        for i in range(1, m + 1):
            acc += a[i - 1]*P[m - i]*i
        P.append(acc*(-arb(1)/m))
    return P

def normR(P, R):
    """||P||_R = sum_j |p_j| j! R^{-j}."""
    R = arb(R); tot = arb(0); f = arb(1)
    for j in range(P.degree() + 1):
        if j > 0:
            f = f*j/R
        tot += arb(abs(P[j]).upper())*f
    return up(tot)

# ---------------------------------------------------------------- majorant of ||P_m||_R for large m
def Bhat(y):
    y = arb(y)
    assert y < 2*PI
    return up(2 + y/2 - (y/2)/((y/2).tan()))

import numpy as np

def _logA_vec(m, rho, s, cl, logeta):
    l = np.arange(1, m + 1, dtype=float)
    yb = (2*l + 1)/rho
    Bh = 2 + yb/2 - (yb/2)/np.tan(yb/2)
    base = -l*np.log(cl) + np.log1p(-cl**l) - np.log(l) + np.log(Bh) + (l + 1)*np.log(rho) + l*np.log(rho + s) - np.log(l + 1)
    lt = base[None, :] + np.outer(logeta, l)
    mx = lt.max(axis=1)
    return mx + np.log(np.exp(lt - mx[:, None]).sum(axis=1))

def norm_major_float(m, R, cl, s):
    """float estimate (log) of the majorant and the optimising eta."""
    rho = 3*m/(math.e*R)
    if (2*m + 1)/rho >= 2*math.pi*0.99:
        return math.inf, None
    hi = cl/(rho*(rho + s))
    logeta = math.log(hi) + np.linspace(-12, 2, 141)
    LA = _logA_vec(m, rho, s, cl, logeta)
    v = -m*logeta + np.exp(np.minimum(LA, 700))
    i = int(np.argmin(v))
    lf = math.lgamma(3*m + 1) - 3*m*math.log(R*rho)
    return v[i] + max(0.0, lf), math.exp(logeta[i])

_NMC = {}
def norm_major(m, R, cb, s, eta):
    key = (m, R, s, float(cb.lower()), float(cb.upper()))
    if key not in _NMC:
        _NMC[key] = _norm_major(m, R, cb, s, eta)
    return _NMC[key]


def _norm_major(m, R, cb, s, eta):
    """Arb upper bound of ||P_m||_R at the given eta (FULLPOS3 §4.3, (M1)-(M3))."""
    cl = arb(cb.lower())
    rhoA = arb(3*m)/(arb(1).exp()*arb(R)); etaA = arb(eta); tot = arb(0)
    for l in range(1, m + 1):
        Bh = Bhat(arb(2*l + 1)/rhoA)
        tot += etaA**l*(cl**(-l) - 1)/l*Bh*rhoA**(l + 1)*(rhoA + s)**l/(l + 1)
    Pm_rho = etaA**(-m)*tot.exp()
    f3 = arb.fac_ui(3*m)/(arb(R)*rhoA)**(3*m)
    fmax = f3 if f3 > 1 else arb(1)
    return up(Pm_rho*fmax)


# ---------------------------------------------------------------- far remainder
def far_remainder(cb, hmax, nmin, M, b, tau, d_lo):
    """upper bound of |R_M| = |sum_k g_k (U_k - T_M(k))| (FULLPOS3 §4.4, author's §61 lemma, fixed for c-cells).
    U_k(h; c) = prod_{j<k} (c - h u_j)/(1 - h u_j) is analytic in |h| < 1/u_{k-1} (also past termination: exact zero at j = n), and
    T_M(k) is its Taylor polynomial in h.  On |h| = rho (x_j = rho u_j < 1), uniformly over c in [c-, c+]:
      M_k(rho) = prod_{j<k} max{(c+ + x_j)/(1 + x_j), |c- - x_j|/(1 - x_j), |c+ - x_j|/(1 - x_j)}   (Moebius extremes on the real axis)
      |U_k - T_M(k)| <= M_k(rho) (h/rho)^{M+1}/(1 - h/rho)        (h <= hmax < rho)
      |U_k - T_M(k)| <= c+^k + M_k(rho) sum_{m<=M} (h/rho)^m           (any admissible rho)
    Blocks [k0, k1): each factor is quasi-convex in x, so its max over a block is at an end point; B_k increasing; then a geometric tail
    with rho = theta/u_k, theta = 3/4."""
    s = 2*b + 1
    kap = A(abs(Fr(9, 2) - b*(b - 1)))
    KcA = (PI/2)*(2*PSI2.absexp(PI/2, 2000) + (2/LAM)*D.Rint(PI/2)[0])
    KcD = LAM**2*(PI/2)*PSI3.absexp(PI/2, 2000)
    e3 = arb(3).exp(); e4 = arb(4).exp()
    def BA(k):
        return KcA*e3*arb(k + 3)**3/27
    def Bk(k):
        return BA(k) + tau*(arb(k + 1)**2/2*BA(k + 1) + (arb(k)**2/2 + 4*k + kap)*BA(k)) + tau*KcD*e4*arb(k + 4)**4/256
    cm = arb(cb.lower()); cp = arb(cb.upper()); h = arb(hmax)
    def logfac(x):
        v = (cp + x)/(1 + x)
        v2 = abs(cm - x)/(1 - x); v3 = abs(cp - x)/(1 - x)
        m = v
        if v2 > m: m = v2
        if v3 > m: m = v3
        return up(m.log()) if m > 0 else arb(-10**6)
    rhos = [arb(10)**(arb(i)/4) for i in range(-80, 1)]          # 1e-20 .. 1
    LM = [arb(0) for _ in rhos]; valid = [True]*len(rhos)
    tot = arb(0); k0 = 0; logcp = cp.log()
    th = arb(3)/4
    KT = None
    while True:
        blk = max(1, k0//400)
        k1 = k0 + blk
        ua = arb(k0*(k0 + s)); ub = arb((k1 - 1)*(k1 - 1 + s))
        best = None
        for i, r in enumerate(rhos):
            if not valid[i]:
                continue
            if not (r*ub < 1):
                valid[i] = False; continue
            lf = logfac(r*ua); lf2 = logfac(r*ub)
            if lf2 > lf: lf = lf2
            LMblk = LM[i] + (blk*lf if lf > 0 else arb(0))        # bound of log M_k for k in [k0, k1]
            Mk = LMblk.exp()
            q = h/r
            fb = cp**k0 + Mk*((M + 1) if q <= 1 else (M + 1)*q**M)
            cand = fb
            if q < 1:
                tr = Mk*q**(M + 1)/(1 - q)
                if tr < cand:
                    cand = tr
            if best is None or cand < best:
                best = arb(cand.upper())
            LM[i] = LM[i] + blk*lf                                  # log M_{k1} <= log M_{k0} + sum over the block
        if best is None:
            raise RuntimeError('no admissible rho at k=%d' % k1)
        tot += up(arb(blk)*Bk(k1)*best)
        k0 = k1
        # switch to the analytic tail once its ratio bound is < 1 and terms are negligible
        if k0 > 1000:
            K = arb(k0)
            q_th = (cp + th)/(1 + th)
            for v in (cp, abs(cm - th)/(1 - th), abs(cp - th)/(1 - th)):
                if v > q_th: q_th = v
            hu = h*K*(K + s)/th
            rr = q_th*(1 + 1/K)**(2*M + 8)
            rr2 = cp*(1 + 1/K)**8
            if rr < 1 and rr2 < 1:
                fK = Bk(k0)*(cp**k0 + (M + 1)*q_th**k0*(hu**M if hu > 1 else arb(1)))
                if fK < arb(10)**(-400):
                    rmax = rr if rr > rr2 else rr2
                    tot += up(2*fK/(1 - rmax))
                    KT = k0
                    break
        if k0 > 5*10**8:
            raise RuntimeError('far remainder did not converge')
    return up(tot), KT


def state_params(cb_lo, a0, s):
    """n_min (lower) and h_max (upper) over the cell for a >= a0: n(n+s) = c (n+a)(n+a+s)."""
    c = arb(cb_lo); a = arb(a0)
    Aq = 1 - c; Bq = s - c*(2*a + s); Cq = -c*a*(a + s)
    n = (-Bq + (Bq*Bq - 4*Aq*Cq).sqrt())/(2*Aq)
    N = n + a
    return lo(n), up(1/(N*(N + s)))


def geometry(ba, bb, verbose=False):
    t0 = time.time()
    ba = Fr(ba); bb = Fr(bb)
    c_lo = (A(ba)/2).tanh()**2; c_hi = (A(bb)/2).tanh()**2
    Rmax = 0.92*float(lo(-(arb(c_hi.upper()).log())).lower())
    Rs = [Rmax*(1e-5)**(1 - i/40) for i in range(41)]
    sups = []
    for R in Rs:
        S = circle_sups(lo(c_lo), up(c_hi), R, 0)
        for kk in list(S):
            if not S[kk].is_finite():
                S[kk] = arb('inf')
        sups.append((R, S))
    if verbose:
        print(f'   geometry [{time.time()-t0:.0f}s]', flush=True)
    return sups, Rmax


def run(b, ba, bb, a0, M=None, m1=30, verbose=True, geo=None):
    t0 = time.time()
    s = 2*b + 1; sig = 1 if b == 0 else -1; kap = Fr(9, 2) - b*(b - 1)
    ba = Fr(ba); bb = Fr(bb)
    c_lo = (A(ba)/2).tanh()**2; c_hi = (A(bb)/2).tanh()**2
    cb = lo(c_lo).union(up(c_hi))
    d_lo = 1 - up(c_hi); w3 = (1 - lo(c_lo))**3
    # floor
    fc = [F for (x, y, F) in D.fcert_cells() if not (y <= ba or x >= bb)]
    Fcert = min(fc)
    nmin, hmax = state_params(c_lo, a0, s)
    tau = hmax                                            # tau = 1/T <= h
    if geo is None:
        geo = geometry(ba, bb, verbose)
    sups, Rmax = geo
    keysA = [k for k in sups[0][1] if k == 'As' or k == 'Ae' or k[0] == 'r']
    keysD = ['Ds', 'De']
    # per component: S2(R) = S(R) + tau Dfac(R) S(R2) for G_A components; tau S(R) for D components
    def at_least(Rq):
        for (R, S) in sups:
            if R >= Rq:
                return R, S
        return sups[-1]
    comp = {k: [] for k in keysA + keysD}
    for i, (R, S) in enumerate(sups):
        for k in keysA:
            best = None
            for (R2, S2_) in sups[i + 1:]:
                if not S2_[k].is_finite():
                    break
                dl = arb(R2 - R); zinv = S2_['zinv']
                Dfac = (1 + zinv)/2*(2/dl**2 + 1/dl) + (abs(A(Fr(7, 2) - 2*b)) + zinv/2)/dl + abs(A(kap))
                v = S[k] + tau*Dfac*S2_[k]
                if best is None or v < best:
                    best = v
            comp[k].append((R, up(best) if (best is not None and best.is_finite()) else arb('inf')))
        for k in keysD:
            comp[k].append((R, up(tau*S[k]) if S[k].is_finite() else arb('inf')))
    # m = 0 extras: G^r(c) and tau D_b G_A(c)
    beta_lo = A(ba)
    Gr = (2/LAM)*(1 - cb.upper())**(-3)*D.P_mhalf(beta_lo)*rint(arb(0))
    best_tD = arb(0)
    for k in keysA:
        bk = None
        for (R, S) in sups:
            Rr = arb(R)
            v = tau*((1 - lo(c_lo))/(2*lo(c_lo))*(2/Rr**2 + 1/Rr) + (abs(A(Fr(7, 2) - 2*b)) + 1/(2*lo(c_lo)))/Rr + abs(A(kap)))*S[k]
            if bk is None or v < bk:
                bk = up(v)
        best_tD += bk
    # local orders, exact P_m for m <= m1: sum over components of min over R
    Pex = Pm_exact(cb, s, m1)
    loc = []
    for m in range(1, m1 + 1):
        tot = arb(0)
        for k, lst in comp.items():
            bm = None
            for (R, St) in lst:
                if not St.is_finite():
                    continue
                v = hmax**m*normR(Pex[m], R)*St
                if bm is None or v < bm:
                    bm = up(v)
            tot += bm
        loc.append(up(tot))
        if verbose:
            print(f'      exact m={m}: {float(tot.upper()):.3e}', flush=True)
    # choose M: needs the far remainder; try list
    if M:
        Ms = [M]
    else:
        Ms = []
        for MM in [100, 200, 300, 450, 650, 900, 1200, 1600]:
            try:
                fr, _k = far_remainder(cb, hmax, int(float(nmin.lower())), MM, b, tau, d_lo)
            except (RuntimeError, AssertionError):
                continue
            if w3*fr < arb(Fcert)*1e-3:
                Ms = [MM]; break
        if not Ms:
            out = dict(b=b, ba=str(ba), bb=str(bb), a0=a0, Fcert=Fcert, ok=False, err=float('inf'), why='far remainder',
                       nmin=float(nmin.lower()), hmax=float(hmax.upper()), secs=time.time() - t0)
            return out
    best_total = None
    for MM in Ms:
        tail_loc = arb(0); mm_detail = []
        cl = float(cb.lower())
        for m in range(m1 + 1, MM + 1):
            # float log-majorant per radius (shared by components)
            nm = {}
            for (R, _S) in sups:
                if R > 3:
                    continue
                lv, eta = norm_major_float(m, R, cl, s)
                if eta is not None:
                    nm[R] = (lv, eta)
            bm_tot = arb(0); percomp = {}
            for k, lst in comp.items():
                cand = [(nm[R][0] + float(St.log().upper()), R) for (R, St) in lst if R in nm and St.is_finite() and St > 0]
                if not cand:
                    continue
                lvb, R = min(cand)
                St = dict(lst)[R]
                vv = hmax**m*norm_major(m, R, cb, s, nm[R][1])*St
                bm_tot += vv; percomp[k] = (float(vv.upper()), R)
            bm = up(bm_tot)
            if verbose and (m <= m1 + 2 or m % 100 == 0):
                worst = max(percomp.items(), key=lambda kv: kv[1][0])
                print(f'      m={m}: {float(bm.upper()):.3e}  worst component {worst[0]} {worst[1][0]:.2e} at R={worst[1][1]:.2e}', flush=True)
            tail_loc += bm
            if m <= m1 + 3 or m % 50 == 0:
                mm_detail.append((m, float(bm.upper())))
        try:
            far, kend = far_remainder(cb, hmax, int(float(nmin.lower())), MM, b, tau, d_lo)
        except RuntimeError as e:
            far, kend = arb('inf'), None
        err = w3*(Gr + best_tD + sum(loc, arb(0)) + tail_loc + far)
        res = dict(M=MM, err=float(err.upper()), far=float((w3*far).upper()), tail_loc=float((w3*tail_loc).upper()), kend=kend,
                   mdet=mm_detail)
        if best_total is None or err < best_total[0]:
            best_total = (err, res)
        if verbose:
            print(f'   M={MM}: err<={float(err.upper()):.3e}  far<={res["far"]:.3e}  higher local<={res["tail_loc"]:.3e}  kend={kend}', flush=True)
        if err < arb(Fcert)/2:
            break
    err, res = best_total
    ok = err < arb(Fcert)
    out = dict(b=b, ba=str(ba), bb=str(bb), a0=a0, Fcert=Fcert, ok=bool(ok),
               Gr=float((w3*Gr).upper()), tD=float((w3*best_tD).upper()), loc=[float((w3*x).upper()) for x in loc],
               nmin=float(nmin.lower()), hmax=float(hmax.upper()), secs=time.time() - t0, **res)
    return out


A0LIST = [998, 1100, 1250, 1400, 1600, 1800, 2000, 2250, 2500, 2800, 3200, 3600, 4000, 4500, 5000, 6000, 7000, 8500, 10000, 12500,
          15000, 20000, 25000, 30000, 40000, 50000, 70000, 100000, 200000, 500000, 2200000]


def scan_cell(ba, bb, out):
    """smallest a0 in A0LIST certified for each channel on the beta-cell [ba, bb] (monotone in a0: every bound increases with h_max
    and decreases with n_min)."""
    t0 = time.time()
    geo = geometry(ba, bb)
    res = {}
    for b in (0, 4):
        lo_i, hi_i = -1, len(A0LIST) - 1
        rec = {}
        r = run(b, ba, bb, A0LIST[hi_i], verbose=False, geo=geo); rec[hi_i] = r
        if not r['ok']:
            res[b] = None; continue
        while hi_i - lo_i > 1:
            mid = (lo_i + hi_i)//2
            r = run(b, ba, bb, A0LIST[mid], verbose=False, geo=geo); rec[mid] = r
            if r['ok']:
                hi_i = mid
            else:
                lo_i = mid
        res[b] = rec[hi_i]
    line = dict(ba=str(ba), bb=str(bb), a0_b0=res[0]['a0'] if res[0] else None, a0_b4=res[4]['a0'] if res[4] else None,
                rec0=res[0], rec4=res[4], secs=time.time() - t0)
    print(f"cell [{ba},{bb}]: a0(b=0) = {line['a0_b0']}, a0(b=4) = {line['a0_b4']}  "
          f"(err/Fcert = {res[0]['err']/res[0]['Fcert'] if res[0] else None}, {res[4]['err']/res[4]['Fcert'] if res[4] else None})  [{line['secs']:.0f}s]", flush=True)
    open(out, 'a').write(json.dumps(line) + '\n')
    return line


if __name__ == '__main__':
    if sys.argv[1] == 'scan':
        cells = [l.split(',') for l in sys.argv[2].split(';')]
        out = os.path.join(HERE, 'cache', sys.argv[3])
        for ba, bb in cells:
            scan_cell(ba, bb, out)
    else:
        b = int(sys.argv[1]); ba = sys.argv[2]; bb = sys.argv[3]; a0 = int(sys.argv[4])
        M = int(sys.argv[5]) if len(sys.argv) > 5 else None
        r = run(b, ba, bb, a0, M)
        print(json.dumps(r))

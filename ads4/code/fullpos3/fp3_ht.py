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
BGRID = [Fr(1, 4)*Fr(103, 100)**i for i in range(0, 125)]          # 0.25 .. ~9.7
_JT = {}
def jt(k):
    if k not in _JT:
        _JT[k] = D.jtable(PSI2 if k == 2 else PSI3, BGRID, ng=400)
    return _JT[k]

def jval(k, beta_lo):
    """J_s, J_C, J_e at the largest grid point <= beta_lo (upper bounds, decreasing in beta)."""
    i = max(i for i, g in enumerate(BGRID) if A(g) <= beta_lo)
    return jt(k)[i]

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
def _arc_geom(c0, R, pa, pb):
    """ball enclosures on the arc t = R e^{i phi}, phi in [pa, pb], z = c0 e^t."""
    phb = arb((pa + pb)/2, (pb - pa)/2)
    t = acb(R*phb.cos(), R*phb.sin())
    z = acb(c0)*t.exp()
    X = (1 + z)/(1 - z)
    xi = X.acosh()
    return z, xi


def circle_sups(c0, R, b, tol=2e-4):
    """(S_A, S_D, max|1/z|): S_A >= sup|G_A|, S_D >= sup|G_D| on the circle z = c0 e^t, |t| = R (adaptive arcs, tolerance tol on
    the enclosure width of Im xi so that e^{Lam eta} is not inflated by more than e^{Lam tol})."""
    SA = arb(0); SD = arb(0); zinv = arb(0)
    R = arb(R); c0 = arb(c0)
    stack = [(arb(0), PI)]                    # symmetric: conj(z) gives conj(xi); |G| symmetric -> phi in [0, pi]
    nseg = 0
    while stack:
        pa, pb = stack.pop()
        z, xi = _arc_geom(c0, R, pa, pb)
        ok = xi.real.is_finite() and xi.imag.is_finite()
        if ok:
            beta = lo(xi.real); eta = up(abs(xi.imag))
            ok = beta > 0 and eta < PI/2 and float(xi.imag.rad()) < tol and float(xi.real.rad()) < 0.01*float(beta.lower())
        if not ok:
            if float((pb - pa).upper()) < 1e-9:
                raise RuntimeError('arc geometry failed at R=%s' % R)
            m = (pa + pb)/2
            stack.append((pa, m)); stack.append((m, pb))
            continue
        nseg += 1
        axi = up(abs(xi)); om = lo(abs(1 - z))
        Js, Jc, Je = jval(2, beta)
        psi_part = 2*om**(-3)*(2**0.5/PI)*(Js + ccor(2, eta)*Jc + (LAM*eta).exp()*Je)
        r_part = (2/LAM)*om**(-3)*(axi/beta)*D.P_mhalf(beta)*rint(eta)
        v = up(psi_part + r_part)
        if v > SA:
            SA = v
        Js3, Jc3, Je3 = jval(3, beta)
        dpart = up(LAM**2*om**(-4)*(2**0.5/PI)*(Js3 + ccor(3, eta)*Jc3 + (LAM*eta).exp()*Je3))
        if dpart > SD:
            SD = dpart
        zi = up(1/abs(z))
        if zi > zinv:
            zinv = zi
    return SA, SD, zinv


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
    base = np.log(cl**(-l) - 1) - np.log(l) + np.log(Bh) + (l + 1)*np.log(rho) + l*np.log(rho + s) - np.log(l + 1)
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

def norm_major(m, R, cb, s, eta):
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
    """upper bound of |R_M| = |sum_{k<=n} c^k g_k rho_M(k) - sum_{k>n} c^k g_k T_M(k)|  (FULLPOS3 §4.4).
    Blocks [k0, k1) up to K_T = max(nmin+1, 6(2M+10)/d), then the analytic geometric tail."""
    s = 2*b + 1
    kap = abs(Fr(9, 2) - b*(b - 1))
    KcA = (PI/2)*(2*PSI2.absexp(PI/2, 2000) + (2/LAM)*D.Rint(PI/2)[0])
    KcD = LAM**2*(PI/2)*PSI3.absexp(PI/2, 2000)
    e3 = arb(3).exp(); e4 = arb(4).exp()
    def BA(k):
        return KcA*e3*arb(k + 3)**3/27
    def Bk(k):
        return BA(k) + tau*(arb(k + 1)**2/2*BA(k + 1) + (arb(k)**2/2 + 4*k + kap)*BA(k)) + tau*KcD*e4*arb(k + 4)**4/256
    c = arb(cb.upper()); h = arb(hmax); dd = arb(d_lo.lower())
    ts = [arb(10)**(arb(i)/4) for i in range(-16, 40)]
    rhos = [h*t for t in ts]
    logInvE = [arb(0) for _ in rhos]
    valid = [True]*len(rhos)
    tot = arb(0); k0 = 0; logc = c.log()
    KT = max(nmin + 1, int(6*(2*M + 10)/float(dd.lower())) + 1)
    while k0 < KT:
        blk = max(1, k0//500)
        k1 = min(k0 + blk, KT); blk = k1 - k0
        uj = arb((k1 - 1)*(k1 - 1 + s))
        for i, r in enumerate(rhos):
            if not valid[i]:
                continue
            if not (r*uj < c*0.999):
                valid[i] = False; continue
            phi = (1 - r*uj).log() - (1 - r*uj/c).log()
            logInvE[i] = logInvE[i] + blk*up(phi)
        best = None
        for i, r in enumerate(rhos):
            if not valid[i]:
                continue
            q = h/r; IE = logInvE[i].exp()
            TM = (M + 1)*IE if q < 1 else (M + 1)*q**M*IE
            cand = 1 + TM
            if k1 - 1 <= nmin and r >= h:
                cand2 = q**(M + 1)*IE
                if cand2 < cand:
                    cand = cand2
            if best is None or cand < best:
                best = arb(cand.upper())
        if best is None:
            raise RuntimeError('no valid rho at k=%d' % k1)
        tot += up(arb(blk)*(logc*k0).exp()*Bk(k1)*best)
        k0 = k1
    # analytic tail k >= KT: rho_k = theta c/u_k, theta = d/2: 1/E_k(rho_k) <= exp(theta (k/3 + b + 1)/(1 - theta)),
    #   f(k) = c^k B_k (1 + (M+1) max(1, h u_k/(theta c))^M exp(theta (k/3+b+1)/(1-theta)));
    #   f(k+1)/f(k) <= c e^{theta/(3(1-theta))} (1 + 1/k)^{2M+5} * (B_{k+1}/B_k (1+1/k)^{-5} <= 1)  =: r < 1  ->  tail <= f(KT)/(1-r)
    th = dd/2; K = arb(KT)
    uK = K*(K + s)
    q = h*uK/(th*c); qM = q**M if q > 1 else arb(1)
    fK = (logc*K).exp()*Bk(KT)*(1 + (M + 1)*qM*(th*(K/3 + b + 1)/(1 - th)).exp())
    rr = c*(th/(3*(1 - th))).exp()*(1 + 1/K)**(2*M + 5)
    assert rr < 1, ('tail ratio', rr)
    tail = fK/(1 - rr)
    return up(tot + tail), KT


def state_params(cb_lo, a0, s):
    """n_min (lower) and h_max (upper) over the cell for a >= a0: n(n+s) = c (n+a)(n+a+s)."""
    c = arb(cb_lo); a = arb(a0)
    Aq = 1 - c; Bq = s - c*(2*a + s); Cq = -c*a*(a + s)
    n = (-Bq + (Bq*Bq - 4*Aq*Cq).sqrt())/(2*Aq)
    N = n + a
    return lo(n), up(1/(N*(N + s)))


def run(b, ba, bb, a0, M=None, m1=8, verbose=True):
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
    c0 = ((A(ba) + A(bb))/4).tanh()**2                   # centre: c at the mid beta
    c0 = arb(c0.mid())
    dlc = up(max(abs((c0/lo(c_lo)).log()), abs((up(c_hi)/c0).log())))   # delta_c >= |log(c/c0)| on the cell
    Rlim = 0.92*float(lo(-(arb(c_hi.upper()).log())).lower())
    Rmax = Rlim - float(dlc.upper())
    Rs = [Rmax*(1e-4)**(1 - i/40) for i in range(41)]
    sups = []
    for R in Rs:
        SA, SD, zinv = circle_sups(c0, R + float(dlc.upper()), b)
        sups.append((R, SA, SD, zinv))
    # S_2(R) = S_A(R) + tau D(R) S_A(R+delta), S_3 = tau S_D ; delta: next grid radius >= 2R (or Rmax)
    def SA_at_least(Rq):
        for (R, SA, SD, zi) in sups:
            if R >= Rq:
                return R, SA
        return sups[-1][0], sups[-1][1]
    Stot = []
    cmaxA = arb(c_hi.upper())
    for (R, SA, SD, zinv) in sups:
        R2, SA2 = SA_at_least(min(2*R, Rmax))
        dl = arb(R2 - R)
        if R2 - R <= 0:
            Dfac = arb('inf')
        else:
            Dfac = (1 + zinv)/2*(2/dl**2 + 1/dl) + (abs(arb(Fr(7, 2) - 2*b)) + zinv/2)/dl + abs(A(kap))
        S2 = SA + tau*Dfac*SA2
        S3 = tau*SD
        Stot.append((R, up(S2 + S3), SA, Dfac, SA2))
    # m = 0 extras: G^r(c) and tau D_b G_A(c)
    beta_lo = A(ba)
    Gr = (2/LAM)*(1 - cb.upper())**(-3)*D.P_mhalf(beta_lo)*rint(arb(0))
    best_tD = None
    for (R, SA, SD, zinv) in sups:
        Rr = arb(R)
        v = tau*((1 - lo(c_lo))/(2*lo(c_lo))*(2/Rr**2 + 1/Rr) + (abs(arb(Fr(7, 2) - 2*b)) + 1/(2*lo(c_lo)))/Rr + abs(A(kap)))*SA
        if best_tD is None or v < best_tD:
            best_tD = up(v)
    # local orders
    Pex = Pm_exact(cb, s, m1)
    loc = []
    for m in range(1, m1 + 1):
        bm = None
        for (R, St, *_r) in Stot:
            v = hmax**m*normR(Pex[m], R)*St
            if bm is None or v < bm:
                bm = up(v)
        loc.append(bm)
    # choose M: needs the far remainder; try list
    Ms = [M] if M else [60, 120, 250, 500, 900, 1400]
    best_total = None
    for MM in Ms:
        tail_loc = arb(0); mm_detail = []
        cl = float(cb.lower())
        for m in range(m1 + 1, MM + 1):
            cand = []
            for (R, St, *_r) in Stot:
                if R > 3:
                    continue
                lv, eta = norm_major_float(m, R, cl, s)
                if eta is None:
                    continue
                cand.append((lv + m*math.log(float(hmax.upper())) + math.log(float(St.upper())), R, St, eta))
            lvb, R, St, eta = min(cand)
            if lvb < -800:                                  # negligible: still bounded rigorously below via Arb
                pass
            bm = up(hmax**m*norm_major(m, R, cb, s, eta)*St)
            tail_loc += bm
            if m <= m1 + 3 or m % 100 == 0:
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
        if err < A(Fr(Fcert))*Fr(1, 2):
            break
    err, res = best_total
    ok = err < arb(Fcert)
    out = dict(b=b, ba=str(ba), bb=str(bb), a0=a0, Fcert=Fcert, ok=bool(ok), err=float(err.upper()),
               Gr=float((w3*Gr).upper()), tD=float((w3*best_tD).upper()), loc=[float((w3*x).upper()) for x in loc],
               nmin=float(nmin.lower()), hmax=float(hmax.upper()), secs=time.time() - t0, **res)
    return out


if __name__ == '__main__':
    b = int(sys.argv[1]); ba = sys.argv[2]; bb = sys.argv[3]; a0 = int(sys.argv[4])
    M = int(sys.argv[5]) if len(sys.argv) > 5 else None
    r = run(b, ba, bb, a0, M)
    print(json.dumps(r))

"""
lemU_cert/lu_cell.py: Lemma U on a whole impact-parameter CELL u = beta_c in [u_lo, u_hi], both channels, every twist a >= 998.

Exact structure used [P, checked numerically in LEMU_CERT.md §3]:
  * R_{m,j}(c, x) = x^m d^{2m-j} R'_{m,j}(v),  v = d/c,  with R'_0 = e_0 and
        R'_{m,j+1} = j(j+s) Bv[j] + j(2j+s-1) Bv[j-1] + j(j-1) Bv[j-2],   Bv[j] = R'_{m-1,j+1} - v R'_{m-1,j};
    so R'_{m,j} is a polynomial of degree <= m in v, and R'_{m,j}(v) = Q_{m,j}(-v) with Q_{m,j} having NONNEGATIVE coefficients
    (induction: Bv = Q[j+1](-v) + (-v) Q[j](-v)).
  * Hence, with h = x d^2 and Y'_j := d^{-j} Y_j = d^3 c^j G^{(j)}(c)/j!:
        P(c, h) = sum_m h^m [A_m + q_b A_{m-1} + B_{m-1}],   A_m = sum_j R'_{m,j}(v) Y'^A_j,  B_m = sum_j R'_{m,j}(v) Y'^B_j.
  * Twist: for c in the cell, every state with a >= 998 has 0 < h <= h_*(c) <= hbar := x_*(c_lo) d_lo^2 (h_* decreasing in c);
    y := h/hbar in [0, 1].
Cell variable w = v - v0 (v0 = d0/c0 at the anchor), |w| <= delta (real).  Rigorous Taylor model of degree I0 + D in w:
  * Y'_j(w) = (d/d0)^3 (c/(c0 d0))^j sum_{i>=j} C(i,j) Y_i(c0) tau^{i-j},  tau = (c - c0)/(c0 d0)  (Taylor shift of the anchor jets):
    exact degree-D Taylor polynomial + remainder <= crude_j(rho) (delta/rho)^{D+1}/(1 - delta/rho) (Cauchy on |w| = rho, complex),
    crude_j from |Y_i(c0)| (balls) for i <= Kj and the integrated Lemma-J envelope beyond.
  * R'_{m,j}(v0 + w): exact Taylor coefficients up to I1 (recursion on truncated polynomials is exact); tail beyond I1
    <= Q_{m,j}(v0 + delta) P(Bin(m, delta/(v0+delta)) > I1) (Chernoff); orders I0 < i <= I1 enclosed by Horner on the ball [-delta, delta].
  * All remainders are carried as ball radii of the constant term (function enclosures on the real cell); arb_poly products then
    bound every cross term.  Final test:  lower(Phat_{0,0}) - sum_{l>=1}|Phat_{0,l}| - sum_{k>=1} delta^k S_k - Dmax^4 e^{-84} > 0,
    S_k = sum_l |Phat_{k,l}| (Chebyshev in y) for k <= KC, else 2 sum_m |p_{k,m}|.
Usage: lu_cell.py u_lo u_hi prec I1 I0 D [out.jsonl] [rho_fac] [precC]
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from flint import arb, acb, arb_poly, ctx
import lu_core as C
import lu_cert1 as C1

M = C1.M


def chernoff_log(m, p, k):
    """log upper bound of P(Bin(m, p) >= k) (0 if k > m)."""
    if k > m:
        return -math.inf
    a = k/m
    if a <= p:
        return 0.0
    v = a*math.log(a/p) + ((1 - a)*math.log((1 - a)/(1 - p)) if a < 1 else 0.0)
    return -m*v + 1e-9


def ballpoly_eval_abs(p, delta):
    """upper bound of sup_{|w| <= delta} |p(w)| (Horner on the ball [-delta, delta])."""
    return abs(p(arb(0, delta))).upper()


def addrad(p, r):
    if p.length() == 0:
        return arb_poly([arb(0, r)])
    c0 = p[0] + arb(0, r)
    q = arb_poly([c0]) + (p - arb_poly([p[0]]))
    return q


def series_inv1(c0, D):
    """exact Taylor coefficients of (1 + c0 w)^{-1}, degree D."""
    return arb_poly([(-c0)**k for k in range(D + 1)])


def _ser(v):
    n = int(ctx.prec*0.302) + 10
    return (v.mid().str(n, radius=False), n, str(arb(v.rad()).upper()))


def _des(t):
    """midpoint was printed with n significant digits: add |m| 10^{-(n-2)} for the decimal rounding."""
    m, n, rs = t
    a = arb(m)
    return a + arb(0, (abs(a)*arb(10)**(-(n - 2))).upper()) + arb(0, arb(rs.strip('[]')).upper() if '+/-' not in rs else arb(rs).upper())


def anchor(c0, prec, pr, Kj):
    """anchor jets with on-disk cache (cache/anchor_*.pkl); the cache stores midpoints to full precision and radii."""
    import pickle, lu_bounds as B
    key = f"anchor_{c0.mid().str(18, radius=False)}_{Kj}_{prec}_{C1.M}.pkl".replace(' ', '')
    fn = os.path.join(HERE, 'cache', key)
    if os.path.exists(fn):
        d = pickle.load(open(fn, 'rb'))
        ctx.prec = d['prec']
        YA = [_des(t) for t in d['YA']]; YD = [_des(t) for t in d['YD']]
        meta = d['meta']; meta['env'] = C1.env_floats(B.Envelope(c0))
        pr(f'anchor jets loaded from {key}')
        return YA, YD, meta
    YA, YD, meta = C1.build(c0, prec, pr, Kj=Kj)
    ctx.prec = prec
    m2 = {k: v for k, v in meta.items() if k != 'env'}
    pickle.dump(dict(prec=prec, YA=[_ser(v) for v in YA], YD=[_ser(v) for v in YD], meta=m2), open(fn + '.tmp', 'wb'))
    os.replace(fn + '.tmp', fn)
    return YA, YD, meta


def cell(u_lo, u_hi, prec, I1, I0, D, out=None, rho_fac=6.0, precC=None, KC=40, Jextra=400):
    T0 = time.time()
    precC = precC or prec
    def pr(*a):
        print(' '.join(str(x) for x in a), flush=True)
    ctx.prec = prec
    u0 = (u_lo + u_hi)/2
    c0 = arb((arb(u0)/2).tanh()**2).mid(); c0 = arb(c0)
    d0 = 1 - c0; v0 = d0/c0
    vv = lambda u: 1/((arb(u)/2).sinh()**2)
    v_hi = vv(u_lo); v_lo = vv(u_hi)
    delta = max(float((v_hi - v0).upper()), float((v0 - v_lo).upper()))*(1 + 1e-9)
    delta_a = arb(delta)
    rho = rho_fac*delta
    c_lo = (arb(u_lo)/2).tanh()**2
    pr(f'cell u in [{u_lo}, {u_hi}]  c0 = {c0.str(15)}  v0 = {v0.str(12)}  delta_w = {delta:.4e} (delta/v0 = {delta/float(v0.mid()):.3e}) '
       f' rho = {rho:.3e}  I1 = {I1} I0 = {I0} D = {D}  prec = {prec}/{precC}')
    assert float((c0*rho).upper()) < 0.5
    # ---------------- anchor jets
    Kj = 3*M + 3 + D + Jextra
    YA, YD, meta = anchor(c0, prec, pr, Kj)
    env = meta['env']; F = meta['F']
    tj = time.time() - T0
    ctx.prec = prec
    # ---------------- Y' polynomials
    t = time.time()
    inv1 = series_inv1(c0, D)
    cw = inv1*c0                                     # c(w) = c0/(1 + c0 w)
    dw = 1 - cw
    tau = (cw - c0)*(1/(c0*d0))
    taus = [arb_poly([1])]
    for k in range(1, D + 1):
        taus.append((taus[-1]*tau).truncate(D + 1))
    P3 = ((dw*(1/d0))**3).truncate(D + 1)
    c0r = float((c0*rho).upper())
    Pd3 = ((1 + float((c0*rho/d0).upper()))/(1 - c0r))**3
    Tm = float((c0*rho/(d0*(1 - c0r))).upper())
    # crude envelope of the anchor jets beyond Kj (integrated Lemma J)
    sys.path.insert(0, os.path.join(C.ADS4, 'code', 'fullpos3'))
    _p = ctx.prec
    import fp3_htdata as HT
    ctx.prec = 256
    lam_a = C.lam()
    ctx.prec = _p
    def amp_log(kind, Th, Lw, ldinv):
        old = ctx.prec; ctx.prec = 256
        try:
            th = arb(Th)
            if kind == 'A':
                v = (2/lam_a)*(lam_a*PSI2.absexp(Th, 400) + HT.Rint(th)[0])
                return float(v.log().upper()) + 3*ldinv + 0.5*Lw
            v = lam_a**2/d0*PSI3.absexp(Th, 400)
            return float(v.log().upper()) + 4*ldinv + 0.5*Lw
        finally:
            ctx.prec = old
    global PSI2, PSI3
    PSI2 = HT.Psi(2); PSI3 = HT.Psi(3)
    ampA = [amp_log('A', Th, Lw, ld) for (lr, ld, Lw, Th) in env.rows_f]
    ampD = [amp_log('D', Th, Lw, ld) for (lr, ld, Lw, Th) in env.rows_f]
    def tail_log(j, amp):
        best = math.inf
        n = Kj + 1
        for (lr, ld, Lw, Th), a in zip(env.rows_f, amp):
            r = math.exp(lr)
            q = n*Tm/((n - j)*r)
            if q >= 1:
                continue
            v = a - n*lr + math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1) + (n - j)*math.log(Tm) - math.log(1 - q)
            best = min(best, v)
        return best + 1e-6
    logYA = [float(abs(v).upper().log().upper()) if v != 0 else -math.inf for v in YA]
    logYD = [float(abs(v).upper().log().upper()) if v != 0 else -math.inf for v in YD]
    lTm = math.log(Tm); ld0 = float(d0.log().lower()); l1 = math.log(1 - c0r)
    def crude_log(j, logY, amp):
        # log sup_{|w|<=rho} |Y'_j(w)|
        terms = [logY[i] + math.lgamma(i + 1) - math.lgamma(j + 1) - math.lgamma(i - j + 1) + (i - j)*lTm for i in range(j, Kj + 1)]
        terms.append(tail_log(j, amp))
        m_ = max(terms)
        s = m_ + math.log(sum(math.exp(x - m_) for x in terms)) + 1e-6
        return 3*math.log(Pd3**(1/3)) - j*l1 - j*ld0 + s
    rem_fac = (D + 1)*math.log(delta/rho) - math.log(1 - delta/rho)
    ctx.prec = precC
    def yprime(Y, logY, amp, jmax):
        out = []
        q = arb_poly([1])
        dpow = arb(1)
        for j in range(jmax + 1):
            S = arb_poly([0])
            for k in range(D + 1):
                if j + k > Kj:
                    break
                S += taus[k]*(arb.bin_uiui(j + k, k)*Y[j + k])
            pj = ((P3*q).truncate(D + 1)*S).truncate(D + 1)*dpow
            rl = crude_log(j, logY, amp) + rem_fac
            out.append(addrad(pj, arb(rl).exp().upper()))
            q = (q*inv1).truncate(D + 1)
            dpow = dpow/d0
        return out
    JB = 3*M
    YpA = yprime(YA, logYA, ampA, JB + 2)
    YpD = yprime(YD, logYD, ampD, JB)
    # exact-function enclosures of c(w), d(w), 1/(2 c d) on the real cell
    crem = float((c0*(c0*delta_a)**(D + 1)/(1 - c0*delta_a)).upper())
    cwb = addrad(cw, crem); dwb = addrad(dw, crem)
    cd = (cwb*dwb)
    cdinv_T = arb_poly([1/(2*cd[0])])   # build T_D of 1/(2cd) by Newton on the exact series
    exact_cd = (cw*dw).truncate(D + 1)
    inv = arb_poly([1/(2*exact_cd[0])])
    for _ in range(int(math.log2(D + 1)) + 2):
        inv = (inv*(2 - (2*exact_cd*inv).truncate(D + 1))).truncate(D + 1)
    # remainder of 1/(2cd) on |w| <= delta: sup on |w| <= rho of |1/(2cd)| times (delta/rho)^{D+1}/(1-delta/rho)
    cmin = float((c0/(1 + c0*rho)).lower()); dmin = float(((d0 - c0*rho)/(1 + c0*rho)).lower())
    assert dmin > 0
    irem = 1/(2*cmin*dmin)*math.exp(rem_fac)
    invb = addrad(inv, irem)
    pr(f"Y' polynomials: {time.time() - t:.0f}s;  max log rad: A {max(float(p[0].rad().log()) for p in YpA if p[0].rad() > 0):.1f}")
    ty = time.time() - t
    res = {}
    chans = [int(x) for x in os.environ.get('LU_CHANNELS', '0,4').split(',')]
    for b in chans:
        t = time.time()
        s, q_b, sig, kb, e = C.chan(b)
        # Y'^B
        YpB = []
        for j in range(JB + 1):
            a2 = (j + 2)*(j + 1)*(dwb*dwb)
            a1 = ((dwb - cwb)*j + cwb*e + 1)*(j + 1)*dwb
            a0 = (-(cwb*dwb)*(j*(j - 1)) + (cwb*dwb)*(e*j))
            num = a2*YpA[j + 2] + a1*YpA[j + 1] + a0*YpA[j]
            v = num*invb + YpA[j]*kb + YpD[j]*sig
            low = v.truncate(D + 1); hi = v.right_shift(D + 1)
            r = ballpoly_eval_abs(hi, delta)*arb(delta)**(D + 1) if hi.length() > 0 else arb(0)
            YpB.append(addrad(low, r.upper()))
        # hbar, Dmax
        xs_lo = C.xstar(c_lo, s); hbar = arb((xs_lo*(1 - c_lo)**2).upper())
        Dmax = 1 + q_b*hbar
        # R' recursion (truncated at I1, exact), majorant Q at v0 + delta
        ctx.prec = prec
        vhat = v0 + delta_a
        pbin = float((delta_a/vhat).upper())
        prev = [arb_poly([1])]; prevQ = [arb(1)]
        Am = [YpA[0]]; Bm = [YpB[0]]
        vpol = arb_poly([v0, 1])
        lqmax = 0.0
        for m in range(1, M + 1):
            top = len(prev) - 1
            Bv = [(prev[jj + 1] - (vpol*prev[jj]).truncate(I1 + 1)) for jj in range(top)] + [-(vpol*prev[top]).truncate(I1 + 1)]
            BQ = [prevQ[jj + 1] + vhat*prevQ[jj] for jj in range(top)] + [vhat*prevQ[top]]
            cur = [arb_poly([0])]*(3*m + 1); curQ = [arb(0)]*(3*m + 1)
            for jj in range(3*m):
                acc = Bv[jj]*(jj*(jj + s)) if jj <= top else arb_poly([0])
                accQ = BQ[jj]*(jj*(jj + s)) if jj <= top else arb(0)
                if 1 <= jj <= top + 1:
                    acc = acc + Bv[jj - 1]*(jj*(2*jj + s - 1)); accQ += BQ[jj - 1]*(jj*(2*jj + s - 1))
                if 2 <= jj <= top + 2:
                    acc = acc + Bv[jj - 2]*(jj*(jj - 1)); accQ += BQ[jj - 2]*(jj*(jj - 1))
                cur[jj + 1] = acc; curQ[jj + 1] = accQ
            lch = chernoff_log(m, pbin, I1 + 1)
            ctx.prec = precC
            A_ = arb_poly([0]); B_ = arb_poly([0])
            for jj in range(1, 3*m + 1):
                rp = cur[jj]
                low = rp.truncate(I0 + 1)
                hi = rp.right_shift(I0 + 1)
                r = arb(0)
                if hi.length() > 0:
                    r = ballpoly_eval_abs(hi, delta)*arb(delta)**(I0 + 1)
                if lch > -math.inf:
                    r = r + curQ[jj]*arb(lch).exp()
                lowb = addrad(low, r.upper()) if r != 0 else low
                A_ += lowb*YpA[jj]; B_ += lowb*YpB[jj]
            Am.append(A_); Bm.append(B_)
            ctx.prec = prec
            prev = cur; prevQ = curQ
            if m % 100 == 0:
                lq = max(float(v.log().upper()) for v in curQ if v != 0)
                lr_ = max(float(abs(p_[0]).upper().log()) for p_ in cur if p_.length() > 0 and p_[0] != 0)
                pr(f'  b={b} m={m}: {time.time() - t:.0f}s  log max Q = {lq:.0f}, log max |R\'(v0)| = {lr_:.0f}')
        ctx.prec = precC
        # P(w, y) = sum_m (hbar y)^m [A_m + q_b A_{m-1} + B_{m-1}]
        Pm = [Am[0]] + [Am[m] + Am[m - 1]*q_b + Bm[m - 1] for m in range(1, M + 1)] + [Am[M]*q_b + Bm[M]]
        hp = arb(1); pk = {}
        KT = max(p_.length() for p_ in Pm)
        cols = [[None]*(M + 2) for _ in range(KT)]
        for m, p_ in enumerate(Pm):
            for k in range(KT):
                cols[k][m] = (p_[k] if k < p_.length() else arb(0))*hp
            hp = hp*hbar
        lowers = {}
        Ph0 = C.cheb(cols[0])
        P0 = Ph0[0]; rest0 = sum((abs(v) for v in Ph0[1:]), arb(0))
        tot = arb(0); Sk = []
        for k in range(1, KT):
            if k <= KC:
                Phk = C.cheb(cols[k]); S = sum((abs(v) for v in Phk), arb(0))
            else:
                S = 2*sum((abs(v) for v in cols[k]), arb(0))
            term = S*arb(delta)**k
            Sk.append(float(term.upper().log()) if term != 0 else -1e9)
            tot += term
        Eref = arb(-84).exp()*Dmax**4
        lower = P0 - rest0 - tot - Eref
        ok = bool(lower > 0)
        r = dict(u_lo=u_lo, u_hi=u_hi, b=b, c0=c0.str(20), delta_w=delta, rho=rho, I1=I1, I0=I0, D=D, prec=prec, precC=precC,
                 F_cert=F, Phat00=P0.str(10), rest0=rest0.str(5), taylor_sum=tot.str(5), lower=lower.str(10), ok=ok,
                 margin_rel=float((lower/P0).lower()) if P0 > 0 else None,
                 log_terms_k=[round(x, 1) for x in Sk[:12]] + [round(max(Sk[12:]), 1) if len(Sk) > 12 else None],
                 t_jets=tj, t_yprime=ty, t_chan=time.time() - t, secs=time.time() - T0)
        pr(json.dumps(r))
        res[b] = r
        if out:
            open(out, 'a').write(json.dumps(r) + '\n')
        ctx.prec = prec
    return res


if __name__ == '__main__':
    a = sys.argv
    u_lo, u_hi = float(a[1]), float(a[2]); prec = int(a[3]); I1, I0, D = int(a[4]), int(a[5]), int(a[6])
    out = a[7] if len(a) > 7 else os.path.join(HERE, 'cache', 'cell.jsonl')
    rf = float(a[8]) if len(a) > 8 else 6.0
    pc = int(a[9]) if len(a) > 9 else None
    cell(u_lo, u_hi, prec, I1, I0, D, out, rf, pc)

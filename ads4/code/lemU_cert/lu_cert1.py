"""
lemU_cert/lu_cert1.py [certified up to the retained inputs listed in LEMU_CERT.md]: Lemma U at ONE impact parameter c,
both channels b in {0, 4}, realised (EPSW-repaired) weight in G_A, ideal psi3 in G_D, M = 1000.

Decomposition of the integrated jets (all as Arb balls; every truncation enters as a ball radius):
  Y^A_j = (2/Lam)[ int_0^{nu_cut} psi F + int_0^{nu_cut} r_theta F + int_{nu_cut}^{NU_E} E_theta F + int_{NU_E}^oo E_theta F ]_j
    (1) psi part: Gauss-Legendre in p on [0, nu_cut/Lam], error <= (len/2)(64/15) rho^{-2N}/(rho^2-1) sup_ellipse|psi y_j|
        (Trefethen ATAP Thm 19.3) with the jet envelope of lu_bounds (Lemma J);
    (2) residual r_theta = f_theta - psi on [0, nu_cut]: |.| <= (2/Lam) sum_cells mass_cell * Ybar_j(cell) (EPSW certified panel sums,
        fp3_htdata.eps_cells, plus the theta correction);
    (3) E_theta on [nu_cut, NU_E]: Gauss-Legendre panels in mu, same error formula with sup_disc|E_theta| (lu_weight);
    (4) tail mu > NU_E: Gaussian tail of the RBF sum times the jet envelope.
  Y^D_j = (Lam^2/d) int_0^1 psi3 y^{(4)}_j dp (GL in p, error as in (1)).
Then Y^B = D_b Y^A + sigma_b Y^D, A_m, B_m by the R recursion, P(y), Chebyshev, all in ball arithmetic.
Test (ANGULAR_ROUTE §80): lower(Phat_0) - upper(sum_{l>=1}|Phat_l|) > D_max^4 E_ref, E_ref = e^{-84} (author §97; RETAINED).
Usage: lu_cert1.py L prec [out.jsonl]      or  lu_cert1.py c=<float> prec [out]
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import numpy as np
from flint import arb, acb, ctx
import lu_core as C
import lu_bounds as B

M = 1000
K = 3*M + 3
NU_E = 1059
SHIFT = 0.0       # c0 delta_w/d0 of a cell (lu_cell): jets beyond 3M+3 enter only through the Taylor shift
SLACK = 1e-6          # additive slack (nats) on every float log-bound (float64 rounding of O(1e3) quantities is < 1e-10)


def log_ybar_np(env, J, p, nur, nui):
    """float64 upper bounds (with SLACK) of log max|y_j| over |nu_r| <= nur, |nu_i| <= nui (Lemma J)."""
    rows = env.rows_f
    j = np.arange(J + 1, dtype=float)
    best = np.full(J + 1, np.inf)
    for (lr, ldinv, Lw, Th) in rows:
        b0 = p*ldinv + (0.5 + nui)*Lw + nur*Th
        best = np.minimum(best, b0 - j*lr)
    return best + SLACK + 1e-12*np.abs(best)


def env_floats(env):
    out = []
    for (r, lr, Lw, Th, dinv) in env.rows:
        # lr: log r (used with -j: need a LOWER bound of log r); others upper bounds
        out.append((float(lr.lower()) - 1e-15, float(dinv.log().upper()) + 1e-15, float(Lw.upper()) + 1e-15, float(Th.upper()) + 1e-15))
    env.rows_f = out
    return env


def lse(v):
    m = np.max(v)
    return m + math.log(np.sum(np.exp(v - m)))


def fcert_at(beta):
    d = json.load(open(os.path.join(C.ADS4, 'code', 'lrtight', 'cache', 'lrt3_rcert.json')))
    v = [r['Fcert'] for r in d['rows'] if Fr(r['ba']) <= Fr(beta) <= Fr(r['bb'])]
    return min(v) if v else None


def colsums(c, b, prec_note=''):
    """rho_j >= sum_{m=0}^{M} |R_{m,j}(c, x_*)| (arb upper bounds), ball arithmetic."""
    s = 2*b + 1; xs = C.xstar(c, s)
    d = 1 - c; d2 = d*d; ic = 1/c
    col = [arb(0)]*(3*M + 1)
    col[0] = arb(1)
    prev = [arb(1)]
    for m in range(1, M + 1):
        top = len(prev) - 1
        Bv = [prev[j + 1] - prev[j]*ic for j in range(top)] + [-prev[top]*ic]
        cur = [arb(0)]*(3*m + 1)
        for j in range(0, 3*m):
            acc = d2*(j*(j + s))*Bv[j] if j <= top else arb(0)
            if 1 <= j <= top + 1:
                acc += d*(j*(2*j + s - 1))*Bv[j - 1]
            if 2 <= j <= top + 2:
                acc += (j*(j - 1))*Bv[j - 2]
            v = xs*acc
            cur[j + 1] = v
            col[j + 1] += abs(v)
        prev = cur
    out = [float(v.upper().log().upper()) if v != 0 else -np.inf for v in col]   # log upper bounds
    return out + [-np.inf]*(K + 1 - len(out))


def omegas(c, logrho, b):
    """log upper bounds of omega^A_j (coefficient abs-sums of Y^A_j in p_m) and omega^D_j, channel b."""
    s, q_b, sig, kb, e = C.chan(b)
    cf = float(c.mid()); df = 1 - cf
    xs = float(C.xstar(c, s).mid())*(1 + 1e-12)
    dx = df*df*xs
    J = len(logrho) - 1
    lr = np.array(logrho)
    j = np.arange(J + 1, dtype=float)
    al = np.log((j + 2)*(j + 1)/(2*cf*df))
    be = np.log(np.abs(((df - cf)*j + e*cf + 1)*(j + 1))/(2*cf*df) + 1e-300)
    ga = np.log(np.abs((-cf*df*j*(j - 1) + e*cf*df*j)/(2*cf*df) + float(kb.mid())) + 1e-300)
    terms = [lr + math.log(1 + dx*q_b), math.log(dx) + ga + lr]
    t1 = np.full(J + 1, -np.inf); t1[2:] = math.log(dx) + al[:-2] + lr[:-2]
    t2 = np.full(J + 1, -np.inf); t2[1:] = math.log(dx) + be[:-1] + lr[:-1]
    oA = np.logaddexp.reduce(np.vstack(terms + [t1, t2]), axis=0) + 1e-6
    oD = lr + math.log(dx) + 1e-6
    return oA, oD


def gl_nodes(N, a, b_):
    """rigorous GL nodes/weights on [a, b_] (arb balls)."""
    out = []
    half = (b_ - a)/2; mid = (b_ + a)/2
    for k in range(N):
        x, w = arb.legendre_p_root(N, k, weight=True)
        out.append((mid + half*x, half*w))
    return out


def choose_gl(logw_vec, env, p, nu_scale, len_half, sup_w_fn, target, Ns=None):
    """choose (N, rho) for GL on an interval mapped from [-1,1]: nu = nu_c0 + nu_scale*x (nu_scale real > 0);
    returns minimal N and the per-j log error vector.  sup_w_fn(rho) = log upper bound of sup|weight| on the ellipse."""
    best = None
    for rho in [1.02, 1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2, 2.5, 3, 4, 6, 9, 14, 22]:
        a_ = (rho + 1/rho)/2; b_ = (rho - 1/rho)/2
        nur = nu_scale[0] + nu_scale[1]*a_; nui = nu_scale[1]*b_
        ly = log_ybar_np(env, K, p, nur, nui)
        base = math.log(len_half*64/15/(rho*rho - 1)) + sup_w_fn(rho)
        tot = lse(logw_vec + ly) + base
        N = max(8, math.ceil((tot - math.log(target))/(2*math.log(rho))))
        if best is None or N < best[0]:
            best = (N, rho, base + ly)
    N, rho, lv = best
    return N, rho, lv - 2*N*math.log(rho)


def build(c, prec, pr, rbudget=0.05, tgt=1e-15, Kj=None, shift=0.0, nu_e=None):
    """integrated jets Y^A, Y^D (j <= Kj) at the exact real point c, every truncation as a ball radius; returns (YA, YD, meta)."""
    global K, SHIFT, NU_E
    K0_, S0_, N0_ = K, SHIFT, NU_E
    if Kj is not None:
        K = Kj
    SHIFT = shift
    if nu_e is not None:
        NU_E = nu_e
    try:
        return _build(c, prec, pr, rbudget, tgt)
    finally:
        K, SHIFT, NU_E = K0_, S0_, N0_


def _build(c, prec, pr, rbudget, tgt):
    ctx.prec = prec
    d = 1 - c
    beta = float((2*c.sqrt().atanh()).mid())
    F = fcert_at(beta)
    LAM = C.lam()
    env = env_floats(B.Envelope(c))
    # ---- pass 1: column sums of |R| (both channels)
    t = time.time()
    lrho = {b: colsums(c, b) for b in (0, 4)}
    om = {b: omegas(c, lrho[b], b) for b in (0, 4)}
    oA = np.maximum(om[0][0], om[4][0]); oD = np.maximum(om[0][1], om[4][1])
    if K > 3*M + 3:
        # jets with j > 3M enter the cell model only through the shift sum C(j, k) tau^k Y_j, |tau| ~ SHIFT per unit w-order:
        # weight omega_pad(j) = max_{j' <= 3M} omega_{j'} C(j, j - j') SHIFT^{j - j'}  (a choice of budget, not a rigour input)
        ls = math.log(max(SHIFT, 1e-300))
        for arr in (oA, oD):
            base = arr[:3*M + 1].copy()
            jp = np.arange(3*M + 1)
            for j in range(3*M + 1, K + 1):
                k = j - jp
                lc = np.array([math.lgamma(j + 1) - math.lgamma(kk + 1) - math.lgamma(j - kk + 1) for kk in k])
                arr[j] = np.max(base + lc + k*ls) + 1.0
    pr(f'pass 1 (|R| column sums): {time.time() - t:.0f}s; max log10 rho_j = {max(max(lrho[0]), max(lrho[4]))/math.log(10):.1f}')
    # ---- residual cells and nu_cut (R1 bound)
    sys.path.insert(0, os.path.join(C.ADS4, 'code', 'fullpos3'))
    _p = ctx.prec
    import fp3_htdata as HT
    ctx.prec = _p
    cells, epsW = HT.eps_cells()
    lam_f = float(LAM.mid())
    cum = -np.inf; nu_cut = None; used = []
    budget = math.log(rbudget*F)
    for (a_, b_, w) in cells:
        if b_ > 950:
            break
        ly = log_ybar_np(env, K, 3, b_, 0.0)
        contrib = math.log(2/lam_f) + math.log(w) + lse(oA + ly)
        if np.logaddexp(cum, contrib) > budget:
            break
        cum = np.logaddexp(cum, contrib); used.append((a_, b_, w, contrib))
    # nu_cut: upper checkpoint of the last fully used cell (cells are widened by 0.1 on both sides; the next cell, which
    # overlaps [nu_cut - 0.2, nu_cut], is included as well, with its own upper edge in the envelope)
    if len(used) < 2:
        raise RuntimeError('residual budget exhausted in the first cells')
    nu_cut = Fr(str(round(used[-2][1] - 0.1, 6)))      # checkpoint value
    pr(f'R1 residual: {len(used)} EPSW cells, nu_cut = {float(nu_cut)}, bound (A-coefficient units) log10 = {cum/math.log(10):.2f} '
       f'({math.exp(cum)/F:.3e} F_cert)')
    # theta correction on [0, nu_cut]: |f_theta - f| <= 2|theta(mu - i) - 1||g|
    import lu_weight as LW
    RW = LW.RealWeight()
    old = ctx.prec; ctx.prec = 128
    nuc = LAM + 100
    th1 = (arb(1)/100).exp()*(((nuc - arb(float(nu_cut)))/10).erfc())          # >= |theta(mu - i) - 1| on [0, nu_cut]
    gs = arb(0)
    m_ = 2.5
    while m_ - 2.5 < float(nu_cut):
        Mf, _ = LW.EW.f_disc_sup(RW.W, max(m_, 2.4), 2.3, nb=48)
        gs = arb(max(gs, Mf).upper()); m_ += 4.5
    M1, M2 = arb(4.4198e5 + 2), arb(6.8782 + 1e-4)                               # EPSW log: disc sups of F1, F2 on |mu| = 0.8
    gint = M1*arb(0.2) + M2*arb(0.53) + gs*arb(float(nu_cut))
    thmass = arb((2*th1*gint).upper())
    ctx.prec = old
    # all residual cells with lower edge < nu_cut (incl. the overlapping next one) + theta mass (envelope at nu_cut + 0.1)
    rlog = -np.inf * np.ones(K + 1)
    for (a_, b_, w) in cells:
        if a_ < float(nu_cut):
            rlog = np.logaddexp(rlog, math.log(w) + log_ybar_np(env, K, 3, b_, 0.0))
    rlog = np.logaddexp(rlog, math.log(float(thmass.upper()) + 1e-300) + log_ybar_np(env, K, 3, float(nu_cut) + 0.1, 0.0))
    rlog += math.log(2/lam_f) + 1e-6
    pr(f'  theta-correction mass on [0, nu_cut] <= {float(thmass.upper()):.3e}; residual bound in P (log10) = {lse(oA + rlog)/math.log(10):.2f}')
    # ---- psi-part GL in p (A on [0, p_cut], D on [0, 1])
    pcut = C.A(nu_cut)/LAM
    pc_f = float(pcut.mid())
    abs2 = [abs(float(v)) for v in C.PSI2]; abs3 = [abs(float(v)) for v in C.PSI3]
    def supA(rho):
        pm = pc_f/2*(1 + (rho + 1/rho)/2)*(1 + 1e-12)
        return math.log(2*sum(a*pm**n for n, a in enumerate(abs2)) + 1e-300) + 1e-9
    def supD(rho):
        pm = 0.5*(1 + (rho + 1/rho)/2)*(1 + 1e-12)
        return math.log(lam_f**2/float(d.mid())*sum(a*pm**n for n, a in enumerate(abs3)) + 1e-300) + 1e-9
    NA, rhoA, errA = choose_gl(oA, env, 3, (lam_f*pc_f/2, lam_f*pc_f/2), pc_f/2, supA, tgt)
    ND, rhoD, errD = choose_gl(oD, env, 4, (lam_f/2, lam_f/2), 0.5, supD, tgt)
    pr(f'psi GL: N_A = {NA} (rho {rhoA}) on [0, {pc_f:.4f}], N_D = {ND} (rho {rhoD}); error in P (log10) A {lse(oA + errA)/math.log(10):.1f}, D {lse(oD + errD)/math.log(10):.1f}')
    # ---- E_theta panels on [nu_cut, NU_E]
    panels = []
    edges = [nu_cut]
    while edges[-1] < NU_E:
        e_ = edges[-1]
        nxt = Fr(4*(math.floor(float(e_))//4) + 4)
        if float(e_) < 12:
            nxt = min(nxt, Fr(math.ceil(float(e_)*1.4*64), 64))
        if nxt - e_ < Fr(1, 64):
            nxt += Fr(1, 64)
        edges.append(min(nxt, Fr(NU_E)))
    errE = -np.inf*np.ones(K + 1)
    t = time.time()
    for a_, b_ in zip(edges[:-1], edges[1:]):
        h = float(b_ - a_)/2; m = float(a_ + b_)/2
        best = None
        for rho in [1.5, 2, 3, 4, 6, 8]:
            R = h*(rho + 1/rho)/2
            if m - R < 0.15:
                continue
            se = float(RW.E_disc_sup(m, R*(1 + 1e-9)).upper().log().upper())
            ly = log_ybar_np(env, K, 3, m + R, h*(rho - 1/rho)/2)
            base = math.log(2/lam_f) + math.log(h*64/15/(rho*rho - 1)) + se + 1e-6
            tot = lse(oA + ly) + base
            n = max(6, math.ceil((tot - math.log(tgt/400))/(2*math.log(rho))))
            if best is None or n < best[0]:
                best = (n, rho, base + ly)
        n, rho, lv = best
        errE = np.logaddexp(errE, lv - 2*n*math.log(rho))
        panels.append((a_, b_, n))
    NE = sum(p_[2] for p_ in panels)
    pr(f'E panels: {len(panels)} panels, {NE} nodes on [{float(nu_cut)}, {NU_E}] (setup {time.time() - t:.0f}s); error in P (log10) {lse(oA + errE)/math.log(10):.1f}')
    # ---- tail mu > NU_E: |E_theta| <= 33 |B(mu - i)| <= 33 * 2S e^{1/(2 smin^2)} e^{-(mu - cmax)^2/(2 smax^2)}  (|QR_-| <= 16.3,
    #      |theta| <= 1.01 for |mu| >= 1000; Y-terms < e^{-5e5});  int_{NU_E}^oo e^{eta mu} e^{-(mu-cmax)^2/(2s^2)} <= sqrt(2pi) s
    #      e^{eta cmax + eta^2 s^2/2} e^{-x0^2/(2 s^2)}, x0 = NU_E - cmax - eta s^2 > 0.
    S_, cmax, smin, smax = HT.gauss_tail_params()
    S_, cmax, smin, smax = float(S_)*(1 + 1e-12), float(cmax), float(smin), float(smax)
    tl = np.full(K + 1, np.inf)
    j = np.arange(K + 1, dtype=float)
    for (lr, ldinv, Lw, Th) in env.rows_f:
        x0 = NU_E - cmax - Th*smax**2
        if x0 <= 0:
            continue
        lt = (math.log(33*2*S_) + 1/(2*smin**2) + 0.5*math.log(2*math.pi) + math.log(smax) + Th*cmax + Th**2*smax**2/2
              - x0**2/(2*smax**2))
        tl = np.minimum(tl, 3*ldinv + 0.5*Lw + lt - j*lr)
    tl += math.log(2/lam_f) + SLACK + 1e-12*np.abs(tl)
    pr(f'tail mu > {NU_E}: error in P (log10) {lse(oA + tl)/math.log(10):.1f}')
    radA = np.logaddexp.reduce(np.vstack([errA, rlog, errE, tl]), axis=0)
    radD = errD
    # ---- jets (ball recurrence wraps like the majorant: extra working precision, lu_core.wrap_bits)
    pj = prec + C.wrap_bits(c, K)     # a-priori majorant estimate; checked against measured losses (LEMU_CERT §2.4)
    ctx.prec = pj
    pr(f'jet precision {pj} bits')
    t = time.time()
    YA = [arb(0)]*(K + 1); YD = [arb(0)]*(K + 1)
    NOJ = bool(os.environ.get('LU_NOJETS'))   # debug only: radii without jets
    for (p, w) in ([] if NOJ else gl_nodes(NA, arb(0), pcut)):
        f = 2*C.polyval(C.PSI2, p)*w
        y = C.jets(3, LAM*p, c, K)
        for jj in range(K + 1):
            YA[jj] += f*y[jj]
    tA = time.time() - t
    t = time.time()
    for (p, w) in ([] if NOJ else gl_nodes(ND, arb(0), arb(1))):
        f = LAM*LAM*C.polyval(C.PSI3, p)*w/d
        y = C.jets(4, LAM*p, c, K)
        for jj in range(K + 1):
            YD[jj] += f*y[jj]
    tD = time.time() - t
    t = time.time()
    cache = {}
    for (a_, b_, n) in ([] if NOJ else panels):
        if n not in cache:
            cache[n] = [arb.legendre_p_root(n, k, weight=True) for k in range(n)]
        h = (C.A(b_) - C.A(a_))/2; mm = (C.A(b_) + C.A(a_))/2
        for (x, w) in cache[n]:
            mu = mm + h*x
            f = (2/LAM)*h*w*RW.E_real(mu)
            y = C.jets(3, mu, c, K)
            for jj in range(K + 1):
                YA[jj] += f*y[jj]
    tE = time.time() - t
    pr(f'jets: A {NA} nodes {tA:.0f}s, D {ND} nodes {tD:.0f}s, E {NE} nodes {tE:.0f}s')
    ctx.prec = prec
    YA = [+v for v in YA]; YD = [+v for v in YD]
    if not NOJ:
        pr(f'  max rel radius of jets: A {max(float((v.rad()/abs(v).upper()).log()) for v in YA if v.rad() > 0)/math.log(10):.1f} (log10)')
    # add the error radii
    for jj in range(K + 1):
        YA[jj] += arb(0, arb(float(radA[jj])).exp().upper())
        YD[jj] += arb(0, arb(float(radD[jj])).exp().upper())
    meta = dict(F=F, nu_cut=float(nu_cut), NA=NA, ND=ND, NE=NE, env=env)
    return YA, YD, meta


def run(L=None, cval=None, prec=1500, out=None, rbudget=0.05, tgt=1e-15):
    ctx.prec = prec
    T0 = time.time()
    log = []
    def pr(*a):
        s = ' '.join(str(x) for x in a); print(s, flush=True); log.append(s)
    if cval is None:
        c = arb((1 - arb(-L).exp()).mid())
    else:
        c = arb(cval)
    d = 1 - c
    beta = float((2*c.sqrt().atanh()).mid()); Lf = float((-(d.log())).mid())
    F = fcert_at(beta)
    LAM = C.lam()
    pr(f'c = {c.mid().str(20)}  L = {Lf:.6f}  beta_c = {beta:.6f}  F_cert = {F}  prec = {prec}')
    YA, YD, meta = build(c, prec, pr, rbudget, tgt)
    nu_cut, NA, ND, NE = meta["nu_cut"], meta["NA"], meta["ND"], meta["NE"]
    # ---- pass 2: contraction, both channels
    res = {}
    for b in (0, 4):
        t = time.time()
        s, q_b, sig, kb, e = C.chan(b)
        YB = C.dress(YA, YD, c, b)
        xs = C.xstar(c, s)
        (Am, _), (Bm, _) = C.contract(c, xs, s, [YA, YB], M)
        pm = C.poly_y(Am, Bm, c, xs, b)
        Ph = C.cheb(pm)
        P0 = Ph[0]; rest = sum((abs(v) for v in Ph[1:]), arb(0))
        Dmax = 1 + q_b*d*d*xs
        Eref = arb(-84).exp()*Dmax**4
        lower = P0 - rest - Eref
        ok = bool(lower > 0)
        r = dict(L=Lf, beta=beta, c=c.mid().str(25), b=b, M=M, prec=prec, F_cert=F, nu_cut=nu_cut, NA=NA, ND=ND, NE=NE,
                 Phat0=P0.str(12), sum_abs_rest=rest.str(12), lower=lower.str(12), ok=ok,
                 margin_rel=float((lower/P0).lower()), lower_over_F=float(lower.lower())/F if F else None,
                 rad_P0=float(P0.rad()), rad_A0=float(Am[0].rad()), A0=Am[0].str(12), H_y1=(sum(pm, arb(0))/Dmax**4).str(10),
                 t_contract=time.time() - t, secs=time.time() - T0)
        pr(json.dumps(r))
        res[b] = r
        if out:
            open(out, 'a').write(json.dumps(r) + '\n')
    return res


if __name__ == '__main__':
    a1 = sys.argv[1]; prec = int(sys.argv[2]); out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, 'cache', 'cert1.jsonl')
    if a1.startswith('c='):
        run(cval=float(a1[2:]), prec=prec, out=out)
    else:
        run(L=float(a1), prec=prec, out=out)

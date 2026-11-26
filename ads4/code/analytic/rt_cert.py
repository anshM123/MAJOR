"""
rt_cert.py: Route T' certificate (ANALYTIC_SLIVER §7): H_ang,b > 0 for every state with L = -log d in [L1, L2], a >= 998, every spin,
b in {0, 4}, via the Regge contour of the polynomial weights (psi2 with the dressing, psi3), the realised residual r = f - psi2 on the
real axis, and the j-truncation remainder by TAIL_CERT's Lemma Loc.  All numbers are Arb upper bounds (floats rounded up).

  d^3 L(G) >= M - (E2 + E3 + Er + Erho),   M >= F_cert(beta_c) (FULLPOS3 §4.1(c), LRINF_TIGHT cells)
  E2   = d^3 (2/Lam) 2 sum_{pieces} int |psi2(nu/Lam)| |B(nu)| sum_j |c_j(nu)| |d^{j-alpha}| Mdiff_j(nu) |dnu|      (alpha = 5/2 + i nu)
  E3   = d^3 tau_T Lam 2 sum_{pieces} int |psi3(nu/Lam)| |B(nu)| sum_j |c_j(nu)| |d^{j-alpha4}| Mdiff4_j(nu) |dnu|  (alpha4 = 7/2 + i nu)
  Er   = d^3 (2/Lam) 2 sum_bands int_band |r| sup_band |K_J(nu)|,     K_J = B sum_{j<J0} c_j V_b(E_{alpha-j})
  Erho = d^3 (2/Lam) int |f| |V_b(R_J0)| + d^3 tau_T Lam int |psi3| |L(R'_J0)|     (Lemma Loc, TAIL_CERT §4)
Pieces of the contour 0 -> delta -> delta + iY -> Lam + iY -> Lam:  bottom (real), left (Re nu = delta), top (Im nu = Y), right (Re nu = Lam).
Mdiff_j >= |d^{beta} V_b(E_beta) - 1| (beta = alpha - j): near-real pieces |S^ - 1| + dressing (Lemma NB / Lemma V), top/right |S^| + 1 + dressing.
Usage: python3 rt_cert.py L1 L2 [dL] [k0] [J_loc3]   (J_loc3 > 0: E_rho by Lemma Loc', rt_loc3.py)
"""
import os, sys, math, json, time
from fractions import Fraction
from flint import arb, acb, ctx, fmpq_poly, fmpq
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import rt_pp as PPm
from an_common import nf02, LAM2
I = acb(0, 1)
LAM = arb(LAM2).sqrt(); LAMf = float(LAM.mid()); LAMup = float(LAM.upper())
A_MIN = 998
DELTA = 0.25
TAU_COEF = 3.9*A_MIN**2          # tau_T <= 1/D <= d^2/(3.9 a^2)


def up(x):
    return float(arb(x).upper())


# ---------------------------------------------------------------- weights
P2, P3 = nf02()
def _poly(coefs):
    return fmpq_poly([fmpq(c.numerator, c.denominator) for c in coefs])
Q2 = _poly(P2); Q3 = _poly(P3)
Q2p = _poly(P2[1:]); Q3p = _poly(P3[1:])           # psi/p


def psi_abs(Q, z0, rad):
    """sup |psi(p)| on the disc |p - z0| <= rad (z0 acb): exact Taylor shift at a rational approximation of z0, then sum |coef| rad^k."""
    zr = fmpq(int(round(float(z0.real.mid())*2**40)), 2**40); zi = fmpq(int(round(float(z0.imag.mid())*2**40)), 2**40)
    # shift in two steps: real part exactly, imaginary part via acb
    sh = Q(fmpq_poly([zr, 1]))
    coeffs = [acb(arb(sh[i].p)/arb(sh[i].q)) for i in range(sh.degree() + 1)]
    # now expand around zr + i zi: p = zr + (i zi + u)
    w = acb(0, arb(zi.p)/arb(zi.q))
    dz = abs(z0 - acb(arb(zr.p)/arb(zr.q), arb(zi.p)/arb(zi.q))).upper()
    R = arb(rad) + arb(dz)
    # value bound: sum_k |c_k| (|w| + R)^k
    tot = arb(0)
    for k, c in enumerate(coeffs):
        tot += abs(c)*(abs(w) + R)**k
    return up(tot)


# ---------------------------------------------------------------- B and c_j
def _digamma_box(z):
    """digamma on a box, shifted to Re z > 8 by psi(z) = psi(z + K) - sum_{k<K} 1/(z + k) (Arb's reflection fails on wide boxes)."""
    K = max(0, int(math.ceil(-float(z.real.mid()) + float(z.real.rad()))) + 8)
    zk = z + K                      # Re zk >= 8: psi(zk) = log zk - 1/(2 zk) + E, |E| <= sec^4(arg/2)/(12|zk|^2) <= 1/(3*64) (|arg zk| < pi/2)
    e = arb(0, 1/150.0)
    v = zk.log() - 1/(2*zk) + acb(e, e)
    for k in range(K):
        v -= 1/(z + k)
    return v


def B_abs_cell(mid, r):
    """sup |B(nu)| on the disc |nu - mid| <= r (mid complex, disc away from the poles nu = i m):
    |B(mid)| (point evaluation) times exp(r sup_disc |d/dnu log B|),  d/dnu log B = 2i[psi(2i nu) - psi(1/2 + i nu)]."""
    m = acb(mid.real, mid.imag)
    Bm = abs((2*I*m).gamma()*((acb(0.5) + I*m).rgamma())**2)
    disc = acb(arb(mid.real, r), arb(mid.imag, r))
    D = 2*abs(_digamma_box(2*I*disc) - _digamma_box(acb(0.5) + I*disc))
    v = Bm*(arb(r)*D).exp()
    return up(v)


def B_abs_real(nu_lo):
    """sup_{nu >= nu_lo} |B(nu)| for real nu (|B|^2 = coth(pi nu)/(4 pi nu), decreasing)."""
    n = arb(nu_lo)
    return up(((arb.pi()*n).coth()/(4*arb.pi()*n)).sqrt())


def nuB_abs(nu_hi):
    """sup_{0 < nu <= nu_hi} |nu B(nu)| = sqrt(nu coth(pi nu)/(4 pi)) (increasing in nu)."""
    n = arb(nu_hi)
    return up((n*(arb.pi()*n).coth()/(4*arb.pi())).sqrt())


def cj_list(nu_box, J):
    out = [1.0]
    c = acb(1); a0 = acb(0.5) - I*nu_box; b0 = acb(1) - 2*I*nu_box
    for j in range(J):
        c = c*(a0 + j)*(a0 + j)/((b0 + j)*(j + 1))
        out.append(up(abs(c)))
    return out


def ratio_bound_box(nu_box, jmax):
    """sup over j >= jmax and nu in the box of |c_{j+1}/c_j| = |a0+j|^2/(|b0+j|(j+1)), a0 = 1/2 - i nu, b0 = 1 - 2i nu:
    |a0 + j| <= |nu| + j + 1/2,  |b0 + j| >= max(j + 1 + 2 Im_lo, 2 Re_lo) (Im_lo = min(0, inf Im nu)), and
    (|nu|+j+1)^2/((j+1) max(j+1+2Im_lo, 2Re_lo)) over j >= jmax is <= the max of its values at j = jmax and at j+1 = 2Re_lo
    (convex in j on the first branch, decreasing on the second)."""
    nmax = float(abs(nu_box).upper())
    re_lo = max(0.0, float(abs(nu_box.real).lower()))
    im_lo = min(0.0, float(nu_box.imag.lower()))
    def f(j):
        return (nmax + j + 1)**2/((j + 1)*max(j + 1 + 2*im_lo, 2*re_lo, 0.1))
    vals = [f(jmax)]
    jc = 2*re_lo - 1
    if jc > jmax:
        vals.append((1 + nmax/(2*re_lo))**2*(2*re_lo)/(2*re_lo + 2*im_lo) if 2*re_lo + 2*im_lo > 0 else f(int(jc)))
        vals.append(f(int(math.floor(jc))) ); vals.append(f(int(math.ceil(jc))))
    vals.append((1 + nmax/max(jmax + 1 + 2*im_lo, 2*re_lo, 0.1))**2*(1 + 1e-9))
    return max(vals)


def Abounds(bmod):
    A1 = bmod*(bmod + 1)/2 + 4*bmod
    A0 = bmod*(bmod + 1)/2 + 4.5*bmod + 7.5
    return A1, A0


# ---------------------------------------------------------------- j-sums
def jsum_near(nu_box, sig0, tau, J0, d, p, dress, s_list=(1, 9)):
    """sum_j |c_j| d^{j - Re alpha} Mdiff_j for near-real exponents beta_j = alpha - j, Re alpha in [lo, hi] = sig0, |Im beta| <= tau.
    Delta(beta) := |S^(beta) - 1| by Lemma NB where Re beta > 0.05, and by Lemma V (downward recurrence) below; the start is two
    consecutive shifted ranges [lo + m, hi + m] with lo + m > 0.05.  Includes the dressing and the (geometric) tail."""
    slo, shi = sig0
    nu_abs = float(abs(nu_box).upper())
    cs = cj_list(nu_box, 60)
    dt = d*d/TAU_COEF if dress else 0.0
    def NBc(lo, hi):
        return max(PPm.NB(lo, hi, tau, s, d) for s in s_list)
    jmax = 59
    m0 = max(0, int(math.ceil(0.05 - slo + 1e-9)))
    # Dl[i] = bound for the range shifted by m = m0 + 1 - i  (i = 0: m0+1, i = 1: m0, ..., i = m0 + 1 + j: shift -j)
    Dl = [NBc(slo + m0 + 1, shi + m0 + 1), NBc(slo + m0, shi + m0)]
    for i in range(2, m0 + 2 + jmax + 1):
        m = m0 + 1 - i
        if slo + m > 0.05:
            Dl.append(NBc(slo + m, shi + m)); continue
        bmod = max(abs(slo + m + 1), abs(shi + m + 1)) + tau + 1
        A, B, den, F = PPm.rec_step(bmod, 9, d)
        Dl.append((A*Dl[-1] + B*Dl[-2] + F)/den)
    def Delta(j):          # range shifted by -j
        return Dl[m0 + 1 + j]
    tot = 0.0; mds = []
    for j in range(0, jmax):
        bmod = max(abs(complex(shi - j, tau)), abs(complex(slo - j, tau)))
        A1, A0 = Abounds(bmod)
        md = Delta(j) + (dt*(A1*(1 + Delta(j - 1))/d + A0*(1 + Delta(j))) if dress else 0.0)
        mds.append(md)
        tot += cs[j]*d**(j - shi)*md
    # tail j >= jmax: |c_j| <= |c_jmax| R^{j-jmax}; M'_j = max(Mdiff_j, Mdiff_{j-1}) grows by <= 2 per step (Lemma V / REC, |beta| <= 1500)
    rb = ratio_bound_box(nu_box, jmax)*d*2.0
    assert rb < 0.5, rb
    Mp = 2*max(mds[-1], mds[-2]) + 1.0
    tot += cs[jmax]*d**(jmax - shi)*Mp/(1 - rb)
    return tot*(1 + 1e-10)


def jsum_osc(nu_box, sig_hi, J0, d, dress, base, chain_bmod):
    """sum_j |c_j| d^{j - sig_hi} Mdiff_j with Mdiff_j <= |S^(beta_j)| + 1 + dressing, |S^| from the chain base = list of bounds for the
    exponents beta_0, beta_1, ... (already computed, length >= jmax+2)."""
    nu_abs = float(abs(nu_box).upper())
    cs = cj_list(nu_box, 60)
    dt = d*d/TAU_COEF if dress else 0.0
    tot = 0.0; jmax = 59; mds = []
    for j in range(jmax):
        A1, A0 = Abounds(chain_bmod + j)
        Sj = base[j + 1]; Sj1 = base[j]     # base[0] = |S^(beta_0 + 1)|, base[k+1] = |S^(beta_k)|
        md = Sj + 1 + (dt*(A1*Sj1/d + A0*Sj) if dress else 0.0)
        mds.append(md)
        tot += cs[j]*d**(j - sig_hi)*md
    rb = ratio_bound_box(nu_box, jmax)*d*2.0
    assert rb < 0.5, rb
    Mp = 2*max(mds[-1], mds[-2]) + 1.0
    tot += cs[jmax]*d**(jmax - sig_hi)*Mp/(1 - rb)
    return tot*(1 + 1e-10)


def chain(base0, base1, sig_top, tau_abs, steps, d):
    """|S^| bounds going down: returns [M(sig_top+1), M(sig_top), M(sig_top-1), ...] (length steps+2) from M(sig_top+1) <= base1, M(sig_top) <= base0."""
    out = [base1, base0]
    for k in range(steps):
        sig = sig_top - k
        bmod = abs(complex(sig, tau_abs)) + 1
        A, B, den, F = PPm.rec_step(bmod, 9, d)
        out.append((A*out[-1] + B*out[-2])/den)
    return out


# ---------------------------------------------------------------- contour pieces
def piece_bottom(L, p, Q, Qp, dress, ncell=20):
    """nu in [0, delta] real: |psi(nu/Lam) B(nu)| <= |psi(p)/p|/Lam |nu B(nu)|."""
    d = math.exp(-L)
    tot = 0.0
    e0 = p - 0.5
    for i in range(ncell):
        a_, b_ = DELTA*i/ncell, DELTA*(i + 1)/ncell
        box = acb(arb((a_ + b_)/2, (b_ - a_)/2))
        w = psi_abs(Qp, acb((a_ + b_)/2/LAMf), (b_ - a_)/2/LAMf + 1e-12)/LAMf*nuB_abs(b_)
        js = jsum_near(box, (e0, e0), b_, 600, d, p, dress)
        tot += w*js*(b_ - a_)
    return tot


def piece_left(L, p, Q, dress, Y, dy=0.05):
    d = math.exp(-L)
    e0 = p - 0.5
    tot = 0.0
    n = int(round(Y/dy))
    for i in range(n):
        y0, y1 = i*dy, (i + 1)*dy
        box = acb(arb(DELTA), arb((y0 + y1)/2, (y1 - y0)/2))
        w = psi_abs(Q, acb(DELTA, (y0 + y1)/2)/LAM, (y1 - y0)/2/LAMf + 1e-12)*B_abs_cell(complex(DELTA, (y0 + y1)/2), (y1 - y0)/2*1.0001)
        # alpha = e0 + i nu = e0 - y + i delta: Re in [e0 - y1, e0 - y0], |Im| = delta;  |d^{j-alpha}| = d^{j - e0 + y}
        js = jsum_near(box, (e0 - y1, e0 - y0), DELTA, 600, d, p, dress)
        tot += w*js*(y1 - y0)
    return tot


def piece_top(L, p, Q, dress, k0, dt_=1.0):
    """nu = t + iY, Y = (p - 1/2) - 1/4 + k0 ... chosen so that Re alpha = 1/4 - k0 (an integer below 1/4): Y = e0 - 1/4 + k0."""
    d = math.exp(-L); d0 = d*1.0001
    e0 = p - 0.5
    Y = e0 - 0.25 + k0
    tot = 0.0
    t = DELTA
    while t < LAMf:
        step = 0.05 if t < 3 else (0.25 if t < 30 else dt_)
        t1 = min(t + step, LAMup)
        box = acb(arb((t + t1)/2, (t1 - t)/2), arb(Y))
        w = psi_abs(Q, acb((t + t1)/2, Y)/LAM, (t1 - t)/2/LAMf + 1e-12)*B_abs_cell(complex((t + t1)/2, Y), (t1 - t)/2*1.0001)
        b0 = PPm.PP(0.25, t, t1, d0); b1 = PPm.PP(1.25, t, t1, d0)
        ch = chain(b0, b1, 0.25, t1, k0 + 62, d)      # ch[0] = M(5/4), ch[1] = M(1/4), ch[2] = M(-3/4), ...
        # beta_j = alpha - j = 1/4 - k0 - j + i t: |S^(beta_j)| = ch[k0 + j + 1], |S^(beta_j + 1)| = ch[k0 + j]
        base = ch[k0:]
        js = jsum_osc(box, 0.25 - k0, 600, d, dress, base, abs(complex(0.25 - k0, t1)))
        # |d^{j - alpha}| = d^{j - 1/4 + k0}: jsum_osc uses d^{j - sig_hi} with sig_hi = 1/4 - k0
        tot += w*js*(t1 - t)
        t = t1
    return tot, Y


def piece_right(L, p, Q, dress, Y, dy=0.25):
    d = math.exp(-L); d0 = d*1.0001
    e0 = p - 0.5
    tot = 0.0
    n = int(math.ceil(Y/dy))
    for i in range(n):
        y0, y1 = i*dy, min((i + 1)*dy, Y)
        box = acb(LAM, arb((y0 + y1)/2, (y1 - y0)/2))
        w = psi_abs(Q, acb(LAM, arb((y0 + y1)/2))/LAM, (y1 - y0)/2/LAMf + 1e-9)*B_abs_cell(complex(LAMf, (y0 + y1)/2), (y1 - y0)/2*1.0001 + 1e-9)
        # beta_0 = e0 - y + i Lam, Re in [e0 - y1, e0 - y0]; base at sig_top = ceil-ish
        stop = e0 - y0
        k = max(0, int(math.ceil(-(e0 - y1))))     # steps needed below 0
        b0 = PPm.PP(max(stop + k, 0.0) + 0.0, LAMf - 1, LAMup, d0)
        b1 = PPm.PP(max(stop + k, 0.0) + 1.0, LAMf - 1, LAMup, d0)
        ch = chain(b0, b1, stop + k, LAMup, k + 62, d)
        base = ch[k:]
        js = jsum_osc(box, stop, 600, d, dress, base, abs(complex(stop, LAMup)) + 1)
        tot += w*js*(y1 - y0)
    return tot


# ---------------------------------------------------------------- residual
EPSW_BANDS = [   # (nu_lo, nu_hi, int|r| upper bound): EPSW_CERT certified ranges / running sums of certified panel bounds
    (0.0, 0.4, 1.1063e-6),
    (0.4, 40.0, 1.1063e-6 + 6.0850e-7),
    (40.0, 898.0, 1.7690e-7),
]
# [898, oo): unit cells from EPSW's own panel bounds regrouped (rt_rmass2.py) + EPSW's Gaussian tail beyond mu_end = 1059


def K_sup_real(nu_lo, nu_hi, L, ncell):
    """sup over real nu in [nu_lo, nu_hi] (nu_lo >= 0.4) of |K_J(nu)| <= |B| sum_j |c_j| d^{j-5/2} (|S^| + 1 + dressing)."""
    d = math.exp(-L); d0 = d*1.0001
    best = 0.0
    h = (nu_hi - nu_lo)/ncell
    for i in range(ncell):
        a_, b_ = nu_lo + i*h, nu_lo + (i + 1)*h
        box = acb(arb((a_ + b_)/2, (b_ - a_)/2))
        Bv = B_abs_real(a_)
        b0 = PPm.PP(2.5, a_, b_, d0); b1 = PPm.PP(3.5, a_, b_, d0)
        ch = chain(b0, b1, 2.5, b_, 62, d)
        js = jsum_osc(box, 2.5, 600, d, True, ch, abs(complex(2.5, b_)))
        best = max(best, Bv*js)
    return best


def residual(L):
    d = math.exp(-L)
    tot = 0.0
    out = []
    # [0, 0.05]: max principle on |nu| = 0.05 for g(nu) = K_J(nu) + K_J(-nu) (analytic in |nu| < 1/2); |g| <= 2 max_circle |K_J|;
    # mass int_0^0.2 |r| <= 7.8938e-7 (EPSW).  [0.05, 0.4]: real-axis bound with |B| decreasing, mass int_0^1 |r| <= 1.1063e-6 (EPSW).
    R0 = 0.05
    nb = 64; best = 0.0
    for k in range(nb):
        th = 2*math.pi*(k + 0.5)/nb
        mid = complex(R0*math.cos(th), R0*math.sin(th)); r = R0*math.pi/nb*1.5
        Bv = B_abs_cell(mid, r)
        box = acb(arb(mid.real, r), arb(mid.imag, r))
        re_hi = 2.5 + abs(mid.imag) + r; re_lo = 2.5 - abs(mid.imag) - r
        js = jsum_near(box, (re_lo, re_hi), abs(mid.real) + r, 600, d, 3, True)
        cs = cj_list(box, 60)
        js1 = sum(cs[j]*d**(j - re_hi) for j in range(60))*1.01
        best = max(best, Bv*(js + js1))
    tot += 2*best*7.8938e-7
    out.append(('[0,0.05]', 2*best*7.8938e-7))
    K = K_sup_real(0.05, 0.4, L, 14)
    tot += 2*K*1.1063e-6
    out.append(('[0.05,0.4]', 2*K*1.1063e-6))
    for (lo, hi, m) in EPSW_BANDS[1:]:
        ncell = max(4, int((hi - lo)/2))
        K = K_sup_real(max(lo, 0.4), hi, L, ncell)
        tot += 2*K*m
        out.append((f'[{lo},{hi}]', 2*K*m))
    cells = json.load(open(os.path.join(HERE, 'cache', 'rt_rmass2.json')))   # EPSW panel bounds regrouped per unit cell (rt_rmass2.py)
    Ts = sorted(int(c[0]) for c in cells if c[1] < 1e8)
    assert Ts[0] <= 898 and all(Ts[i + 1] == Ts[i] + 1 for i in range(len(Ts) - 1)), 'rt_rmass2 incomplete'
    sub = 0.0; worst = (0, 0)
    for (T0, T1, m) in cells:
        if T1 > 1e8:
            # beyond mu_end (1059): EPSW Gaussian tail mass, kernel bounded at nu <= 1150 by PP(1150) and beyond by the far sum below
            K = K_sup_real(1059.0, 1150.0, L, 8)
            sub += 2*K*m
            continue
        K = K_sup_real(float(T0), float(T1), L, 2)
        sub += 2*K*m
        if 2*K*m > worst[1]: worst = (T0, 2*K*m)
    tot += sub
    out.append((f'[898,oo) worst cell {worst[0]}', sub))
    # nu >= 1150: |K| <= |B| sum_j |c_j| d^{j-5/2} Gamma(sigma)/|Gamma(sigma + i nu)| (reflection, S^(sigma) <= 1 real), crude 1e3 e^{pi nu/2} nu^3 d^{-5/2};
    # int_T^{T+1} |f| <= tail(T) (CONTOUR); summed to 3000, Gaussian beyond.
    sys.path.insert(0, os.path.join(ADS4, 'code', 'contour'))
    import cb_common as KC, cb_lines as CL
    Wc = KC.load()
    far = 0.0
    for T in range(1150, 3000):
        intf, _ = CL.tail(Wc, 0.0, T)
        lt = float(intf.log().upper()) if intf > 0 else -1e9
        lv = lt + math.pi*(T + 1)/2 + 3*math.log(T + 1) + math.log(1e3) - 2.5*math.log(d)
        if lv > -700:
            far += math.exp(lv)
    tot += 2*far
    out.append(('[1150,oo)', 2*far))
    return d**3*(2/LAMf)*tot, out


# ---------------------------------------------------------------- remainder (Lemma Loc), F_cert
def remainder_log(L, J0):
    sys.path.insert(0, os.path.join(ADS4, 'code', 'tailcert'))
    import tc_remainder as RM
    return RM


def fcert_min(b_lo, b_hi):
    rows = json.load(open(os.path.join(ADS4, 'code', 'lrtight', 'cache', 'lrt3_rcert.json')))['rows']
    vals = [r['Fcert'] for r in rows if float(Fraction(r['ba'])) <= b_hi and float(Fraction(r['bb'])) >= b_lo]
    assert vals, (b_lo, b_hi)
    return min(vals)


def beta_of_L(L):
    d = math.exp(-L)
    return 2*math.atanh(math.sqrt(1 - d))


def remainder_terms(L):
    """log E_rho2, log E_rho3 (min over the cached J0), see rt_rem.py."""
    sys.path.insert(0, os.path.join(ADS4, 'code', 'tailcert'))
    import tc_remainder as RM
    R = json.load(open(os.path.join(HERE, 'cache', 'rt_rem.json')))
    d = math.exp(-L)
    s3 = sum(abs(float(c)) for c in P3)
    best2 = best3 = (1e300, None)
    for J0s, v in R.items():
        J0 = int(J0s); m = J0 - 4; m3 = J0 - 5
        if m > A_MIN - 3: continue
        small2 = v['lK0_2'] + L + math.log(0.4*21.33815)
        lE2 = 3*math.log(d) + math.log(2/LAMf) + math.log(2) + RM.log_loc_factor(m, L) + max(v['lI2'], small2) + math.log(2)
        small3 = v['lK0_3'] + L + math.log(0.4/LAMf)
        lE3 = 3*math.log(d) + math.log(d*d/TAU_COEF) + 2*math.log(LAMup) + math.log(s3) + math.log(2) + RM.log_loc_factor(m3, L) \
            + max(v['lI3'], small3) + math.log(2)
        if lE2 < best2[0]: best2 = (lE2, J0)
        if lE3 < best3[0]: best3 = (lE3, J0)
    return best2, best3


def remainder_terms_loc3(La, Lb, J):
    """log E_rho2, log E_rho3 by Lemma Loc' (rt_loc3.py), uniform over the L-cell [La, Lb]; same prefactors as remainder_terms.
    Needs the j-sum growth (Lemma V / REC, M'_j grows by <= 2) up to |beta| = |alpha - J| <= hypot(J + 3, 1151) + 2: asserted."""
    import rt_loc3 as R3
    d = math.exp(-La)
    A, B, den, Fi = PPm.rec_step(math.hypot(J + 3, 1151) + 2, 9, d)
    assert (A + B)/den <= 2.0, (A, B, den)
    fn = os.path.join(HERE, 'cache', 'rt_loc3.json')
    cache = json.load(open(fn)) if os.path.exists(fn) else {}
    key = f'{La:.4f}_{Lb:.4f}_{J}'
    if key not in cache:
        lI2, lI3 = R3.integrals(La, J, False, Lb)
        cache = json.load(open(fn)) if os.path.exists(fn) else {}
        cache[key] = [lI2, lI3]
        json.dump(cache, open(fn, 'w'))
    lI2, lI3 = cache[key]
    s3 = sum(abs(float(c)) for c in P3)
    lE2 = 3*math.log(d) + math.log(2/LAMf) + math.log(2) + lI2 + math.log(2)
    lE3 = 3*math.log(d) + math.log(d*d/TAU_COEF) + 2*math.log(LAMup) + math.log(s3) + math.log(2) + lI3 + math.log(2)
    return (lE2, J), (lE3, J)


if __name__ == '__main__':
    L1 = float(sys.argv[1]); L2 = float(sys.argv[2])
    dL = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    k0 = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    JLOC3 = int(sys.argv[5]) if len(sys.argv) > 5 else 0      # > 0: remainder by Lemma Loc' (rt_loc3) with J0 = JLOC3
    t00 = time.time()
    L = L1
    lines = []; rows = []
    def pr(x):
        print(x, flush=True); lines.append(x)
    pr(f'Route T\' certificate, L in [{L1}, {L2}], cells {dL}, top line at Re alpha = 1/4 - {k0}')
    while L < L2 - 1e-9:
        La, Lb = L, min(L + dL, L2)
        d = math.exp(-La)
        res = {}
        for (name, p, Q, Qp, dress) in (('psi2', 3, Q2, Q2p, True), ('psi3', 4, Q3, Q3p, False)):
            bot = piece_bottom(La, p, Q, Qp, dress)
            top, Y = piece_top(La, p, Q, dress, k0)
            left = piece_left(La, p, Q, dress, Y)
            right = piece_right(La, p, Q, dress, Y)
            fac = d**3*(2/LAMf)*2 if name == 'psi2' else d**3*(d*d/TAU_COEF)*LAMup*2
            res[name] = dict(bottom=fac*bot, left=fac*left, top=fac*top, right=fac*right, Y=Y)
        Er, Er_parts = residual(La)
        if JLOC3:
            (lE2, J2), (lE3, J3) = remainder_terms_loc3(La, Lb, JLOC3)
        else:
            (lE2, J2), (lE3, J3) = remainder_terms(La)
        F = fcert_min(beta_of_L(La), beta_of_L(Lb))
        E2 = sum(v for k, v in res['psi2'].items() if k != 'Y'); E3 = sum(v for k, v in res['psi3'].items() if k != 'Y')
        Erho = math.exp(lE2) + math.exp(lE3)
        rho = (E2 + E3 + Er + Erho)/F*(1 + 1e-9)       # float accumulation of positive upper bounds
        pr(f'L in [{La:.2f},{Lb:.2f}] (beta_c in [{beta_of_L(La):.3f},{beta_of_L(Lb):.3f}]): F_cert >= {F:.4e}')
        for nm in ('psi2', 'psi3'):
            r_ = res[nm]
            pr(f'   {nm}: bottom {r_["bottom"]/F:.2e}  left {r_["left"]/F:.2e}  top(Y={r_["Y"]}) {r_["top"]/F:.2e}  right {r_["right"]/F:.2e}')
        pr(f'   residual: {Er/F:.2e}  (' + ', '.join(f'{a}: {b*d**3*(2/LAMf)/F:.1e}' for a, b in Er_parts) + ')')
        pr(f'   Lemma-Loc{chr(39) if JLOC3 else ""} remainders: k=2 {math.exp(lE2)/F:.2e} (J0={J2}), k=3 {math.exp(lE3)/F:.2e} (J0={J3})')
        pr(f'   rho <= {rho:.4e}   {"PASS" if rho < 1 else "FAIL"}   [{time.time()-t00:.0f}s]')
        rows.append(dict(La=La, Lb=Lb, F=F, rho=rho, psi2=res['psi2'], psi3=res['psi3'], Er=Er, lE2=lE2, lE3=lE3))
        L = Lb
    tag = f'{L1}_{L2}' + (f'_loc3J{JLOC3}' if JLOC3 else '')
    json.dump(rows, open(os.path.join(HERE, 'cache', f'rt_cert_{tag}.json'), 'w'))
    open(os.path.join(HERE, 'cache', f'rt_cert_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')

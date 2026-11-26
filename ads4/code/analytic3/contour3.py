"""
contour3.py La Lb [k0]: the four Regge-contour pieces of Route T''' (ANALYTIC_SLIVER3 §3) for the L-cell [La, Lb], as ratios to F_cert.

* bottom, left (near-real exponents): rt_cert.piece_bottom / piece_left (Lemma NB, Lemma V; unchanged).  Their geometric j-tails use the
  REC growth 2 with inhomogeneity F/den <= 1 up to |beta| <= J + 3 + Y + 1 (J = 1500): ASSERTED here at d_hi (both increase with d, |beta|).
* top line nu = t + iY (Re alpha = 1/4 - k0): base |S^(1/4 + it)|, |S^(5/4 + it)| from rt_pp.PP (uniform over a >= 998, b in {0,4});
  the absolute REC chain down to Re = 1/4 - k0 - j (rt_cert.chain); j-sum with JMAX explicit terms (adaptive) and a geometric tail with
  the GENERAL growth constant G = max(2, (A+B)/den at |beta|max) and Mdiff_j <= (1 + dt(A0 + A1/d)) max(|S^_j|, |S^_{j-1}|) + 1
  (the "+1" is |kappa_inf|'s multiplier); requires ratio * d * G < 0.9 (asserted).
* right edge nu = Lam + iy: per y-cell and per a-range (998 <= a <= 4500, kernel3.a_ranges): base powers by rt2_kernel.PP_box at the
  top of the sigma-range (E', Lemma X, Lemma F with the sigma-box), Lemma B multi-frame chain (bchain) down to Re beta_j, the same j-sum.
"""
import sys, os, math, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kernel3 as K3
import rt_cert as RC
import rt_pp as PPm
import rt2_kernel as RK
import bchain as BC
import tpf3 as TF
from flint import arb, acb


def cs_and_ratio(nu_box, jmax):
    return RC.cj_list(nu_box, jmax), RC.ratio_bound_box(nu_box, jmax)


def jsum_general(nu_box, xs, d, dt, sigA_hi, jmax, Jcap, bmod0, plus_one=True):
    """sum_{j<Jcap} |c_j| d^{j - sigA_hi} Mdiff_j,  Mdiff_j = xs[j+1](1 + dt A0) + dt A1 xs[j]/d + (1 if plus_one),
    xs[j+1] >= |S^(beta_j)| (j = -1..jmax-1), explicit j < jmax; tail by REC growth G (computed at |beta| <= bmod0 + Jcap + 3)."""
    cs, ratio = cs_and_ratio(nu_box, jmax)
    tot = 0.0
    for j in range(jmax):
        bm = bmod0 + j + 1
        A1, A0 = RC.Abounds(bm)
        md = xs[j + 1]*(1 + dt*A0) + dt*A1*xs[j]/d + (1.0 if plus_one else 0.0)
        tot += cs[j]*d**(j - sigA_hi)*md
    bmx = bmod0 + Jcap + 3
    A_, B_, den_, F_ = PPm.rec_step(bmx, 9, d)
    G = max(2.0, (A_ + B_)/den_)
    rb = ratio*d*G
    assert rb < 0.9, rb
    A1m, A0m = RC.Abounds(bmx)
    Cd = 1 + dt*(A0m + A1m/d)
    Mx = G*max(xs[jmax], xs[jmax - 1])
    tot += cs[jmax]*d**(jmax - sigA_hi)*(Cd*Mx + (1.0 if plus_one else 0.0))/(1 - min(rb, 0.9))
    return tot*(1 + 1e-10)


def chain_abs(base0, base1, sig_top, tau_abs, steps, d):
    return RC.chain(base0, base1, sig_top, tau_abs, steps, d)


def piece_top(La, p, Q, dress, k0, jmax=160):
    d = math.exp(-La); d0 = d*1.0001
    e0 = p - 0.5
    Y = e0 - 0.25 + k0
    dt = d*d/RC.TAU_COEF if dress else 0.0
    tot = 0.0; t = RC.DELTA
    while t < RC.LAMf:
        step = 0.05 if t < 3 else (0.25 if t < 30 else 1.0)
        t1 = min(t + step, RC.LAMup)
        box = acb(arb((t + t1)/2, (t1 - t)/2), arb(Y))
        w = RC.psi_abs(Q, acb((t + t1)/2, Y)/RC.LAM, (t1 - t)/2/RC.LAMf + 1e-12)*RC.B_abs_cell(complex((t + t1)/2, Y), (t1 - t)/2*1.0001)
        b0 = PPm.PP(0.25, t, t1, d0); b1 = PPm.PP(1.25, t, t1, d0)
        ch = chain_abs(b0, b1, 0.25, t1, k0 + jmax + 2, d)
        xs = ch[k0:]
        js = jsum_general(box, xs, d, dt, 0.25 - k0, jmax, 1500, abs(complex(0.25 - k0, t1)))
        tot += w*js*(t1 - t)
        t = t1
    return tot, Y


def piece_right(La, Lb, p, Q, dress, Y, dy=0.25, jmax=120, ratio=1.03):
    d = math.exp(-La); d_lo = math.exp(-Lb)
    e0 = p - 0.5
    tot = 0.0
    n = int(math.ceil(Y/dy - 1e-9))
    T0 = RC.LAMf - 1e-6; T1 = RC.LAMup + 1e-6
    for i in range(n):
        y0, y1 = i*dy, min((i + 1)*dy, Y)
        box = acb(RC.LAM, arb((y0 + y1)/2, (y1 - y0)/2))
        w = RC.psi_abs(Q, acb(RC.LAM, arb((y0 + y1)/2))/RC.LAM, (y1 - y0)/2/RC.LAMf + 1e-9)*RC.B_abs_cell(complex(RC.LAMf, (y0 + y1)/2), (y1 - y0)/2*1.0001 + 1e-9)
        # beta_0 = alpha = e0 - y + i Lam, Re in [e0 - y1, e0 - y0]; chain from the top Re = s_hi = e0 - y0 + k (k shifts so that s_lo >= 0)
        k = max(0, int(math.ceil(-(e0 - y1) - 1e-12)))
        s_lo, s_hi = e0 - y1 + k, e0 - y0 + k
        best = 0.0
        for (a_lo, a_hi) in K3.a_ranges(T1, ratio):
            M0 = PP_box_sig(s_lo, s_hi, T0, T1, a_lo, a_hi, d_lo, d, dress)
            M1 = PP_box_sig(s_lo + 1, s_hi + 1, T0, T1, a_lo, a_hi, d_lo, d, dress)
            sig0 = arb((s_lo + s_hi)/2, (s_hi - s_lo)/2 + 1e-12)
            xs_full = BC.chain_bounds_multi(sig0, T0, T1, M0, M1, k + jmax + 1, d_lo, d, a_lo, a_hi)
            xs = xs_full[k:]          # xs[0] >= |S^(alpha + 1)|, xs[j+1] >= |S^(alpha - j)|
            _, g, _ = BC.state_box(a_lo, a_hi, d_lo, d)
            dt = float(g.upper())*d*d if dress else 0.0
            v = jsum_general(box, xs, d, dt, e0 - y0, jmax, 1500, abs(complex(e0 - y0, T1)) + 1)
            best = max(best, v)
        tot += w*best*(y1 - y0)
    return tot


_PPS = {}
def PP_box_sig(s_lo, s_hi, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, dress):
    """both channels: bound of |S^(sigma + i tau)| for sigma in [0, s_hi] (E', Lemma X: monotone in sigma), and Lemma F with the
    sigma-box [s_lo, s_hi] (frames when tau/a_lo > 0.8); b = 4 by the helicity bound (rt2_kernel.PP_box formula)."""
    key = (s_lo, s_hi, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi)
    if key in _PPS:
        return _PPS[key]
    def p0(sl, sh):
        c = []
        if t_hi <= 0.95*a_lo:
            c.append(PPm.thmE(sh, t_lo, t_hi, a=a_lo))
        else:
            c.append(PPm.PP0(sh, t_lo, t_hi))
        if t_hi/a_lo > 0.8 and t_lo > 0.5*(a_lo + 200) and sl >= 0:
            il, ih = RK.iota_range(a_lo, a_hi, d_lo, d_hi)
            r = TF.lemmaF3_log(sl, sh, t_lo, t_hi, a_lo, a_hi, il, ih, 200, i_s_list=[200, 190, 180, 170, 155])
            if r is not None:
                c.append(math.exp(r[0])*(1 + 1e-9))
        return min(c)
    q0 = p0(s_lo, s_hi); q1 = p0(s_lo + 1, s_hi + 1)
    d0 = d_hi*(1 + 1e-4)
    m = max(q0, q1)
    corr = sum(math.comb(4, i)*d0**i*2.0**i*m for i in range(1, 5))
    p4 = (1 + 1e-5)*(q0 + corr)/((1 - d0)**4*(1 - 1e-6))
    _PPS[key] = v = max(q0, p4)
    return v


def check_near_real(La, J=1500, Y=8.0):
    """REC hypotheses of rt_cert.jsum_near's geometric tail: (A+B)/den <= 2 and F/den <= 1 for |beta| <= J + 3 + Y + 1."""
    d = math.exp(-La)
    A, B, den, F = PPm.rec_step(J + 4 + Y + 1, 9, d)
    assert (A + B)/den <= 2.0 and F/den <= 1.0, (A, B, den, F)
    return (A + B)/den, F/den


if __name__ == '__main__':
    La, Lb = float(sys.argv[1]), float(sys.argv[2]); k0 = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    which = sys.argv[4].split(',') if len(sys.argv) > 4 else ['bottom', 'left', 'top', 'right']
    d = math.exp(-La); t0 = time.time()
    F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    print('near-real REC check (A+B)/den, F/den =', check_near_real(La), flush=True)
    res = {}
    for (name, p, Q, Qp, dress) in (('psi2', 3, RC.Q2, RC.Q2p, True), ('psi3', 4, RC.Q3, RC.Q3p, False)):
        fac = d**3*(2/RC.LAMf)*2 if name == 'psi2' else d**3*(d*d/RC.TAU_COEF)*RC.LAMup*2
        Y = (p - 0.5) - 0.25 + k0
        r = dict(Y=Y)
        if 'bottom' in which: r['bottom'] = fac*RC.piece_bottom(La, p, Q, Qp, dress)
        if 'left' in which: r['left'] = fac*RC.piece_left(La, p, Q, dress, Y)
        if 'top' in which: r['top'] = fac*piece_top(La, p, Q, dress, k0)[0]
        if 'right' in which: r['right'] = fac*piece_right(La, Lb, p, Q, dress, Y)
        res[name] = r
        print(name, {k: (v/F if k != 'Y' else v) for k, v in r.items()}, f'[{time.time()-t0:.0f}s]', flush=True)
    res['F'] = F; res['k0'] = k0
    json.dump(res, open(os.path.join(HERE, 'cache', f'contour3_{La}_{Lb}_k{k0}.json'), 'w'))
    E2 = sum(v for k, v in res['psi2'].items() if k != 'Y'); E3 = sum(v for k, v in res['psi3'].items() if k != 'Y')
    print('F', F, 'E2/F', E2/F, 'E3/F', E3/F)

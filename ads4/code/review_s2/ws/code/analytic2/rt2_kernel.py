"""
rt2_kernel.py: certified per-cell bound of |K_J(nu)| (real nu in [T0, T1]) uniformly over all states with L = -log d in [La, Lb],
a >= 998, b in {0, 4}, using Lemma F (tpf3, turning point), Theorem E' (rt_pp.thmE, decreasing in a), Lemma X (rt_pp.PP0, a = 998
uniform) for the base powers, the helicity bound for b = 4, and Lemma B (bchain) for the beta-chain.  ANALYTIC_SLIVER2 §4.

  |K_J(nu)| <= |B(nu)| sum_{j<J0} |c_j| d^{j-5/2} Mdiff_j,   Mdiff_j = |S^_j|(1 + dt|A0|) + dt |A1| |S^_{j-1}|/d,
  S^_j := S^(5/2 - j + i nu) (S^_{-1} := S^(7/2 + i nu)),  dt = 1/T <= 1/D = g d^2  (g = 1/(d^2 D), bchain.state_box),
  A1, A0 = TAIL_CERT E2 dressing coefficients (rt_cert.Abounds), explicit j < JMAX, geometric tail beyond (rt_cert.jsum_osc argument).
All d-powers at d = d_hi = e^{-La} (the caller multiplies by d^3: d^{j+1/2} is increasing in d).
"""
import math, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common2
import rt_pp as PPm
import rt_cert as RC
import bchain as BC
import tpf3 as TF
from flint import arb, acb, ctx
JMAX = 60
J0 = 1500
_PPC = {}


def t_range(d_lo, d_hi, a_lo):
    """t = a/N over the states (s = 1 or 9 and the helicity-shifted s = 1 state): t/d in [lo, hi] (Arb)."""
    vals = []
    for dd in (d_lo*(1 - 1e-3), d_hi*(1 + 1e-3)):
        for e in (0.0, 9*0.6*d_hi/a_lo):
            ab = arb(dd); ee = arb(e)
            t = ((2 + ee) - ((2 + ee)**2 - 4*ab*(1 + ee)).sqrt())/2
            vals.append((float(t.lower()), float(t.upper())))
    return min(v[0] for v in vals), max(v[1] for v in vals)


def iota_range(a_lo, a_hi, d_lo, d_hi):
    tl, th = t_range(d_lo, d_hi, a_lo)
    return tl/a_hi*(1 - 1e-4), th/a_lo


def PP0_box(sig, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, use_frames=True, K=200):
    """b = 0 bound of |S^(sig + i tau)| over tau in [t_lo, t_hi], a in [a_lo, a_hi], d in [d_lo, d_hi]."""
    key = (sig, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, use_frames, K)
    if key in _PPC:
        return _PPC[key]
    _PPC[key] = v = _PP0_box(sig, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, use_frames, K)
    return v


def _PP0_box(sig, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, use_frames, K):
    c = []
    if t_hi <= 0.95*a_lo:
        c.append(PPm.thmE(sig, t_lo, t_hi, a=a_lo))
    else:
        c.append(PPm.PP0(sig, t_lo, t_hi))      # Lemma X (uniform a >= 998) / E' at 998
    if c[0] < math.exp(3.0) and t_hi <= 0.95*a_lo:
        return c[0]
    if use_frames and a_hi < math.inf and t_hi/a_lo > 0.8 and t_lo > 0.5*(a_lo + K):
        il, ih = iota_range(a_lo, a_hi, d_lo, d_hi)
        r = TF.lemmaF3_log(sig, sig, t_lo, t_hi, a_lo, a_hi, il, ih, K, i_s_list=[K, K - 10, K - 20, K - 30, K - 45])
        if r is not None:
            c.append(math.exp(r[0])*(1 + 1e-9))
    return min(c)


def PP_box(sig, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, **kw):
    """both channels (helicity identity for b = 4, as rt_pp.PP), d0 = d_hi (1 + 1e-4)."""
    p0 = PP0_box(sig, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, **kw)
    p1 = PP0_box(sig + 1, t_lo, t_hi, a_lo, a_hi, d_lo, d_hi, **kw)
    d0 = d_hi*(1 + 1e-4)
    m = max(p0, p1)
    corr = sum(math.comb(4, i)*d0**i*2.0**i*m for i in range(1, 5))
    p4 = (1 + 1e-5)*(p0 + corr)/((1 - d0)**4*(1 - 1e-6))
    return max(p0, p4)


def jsum(nu_box, xs, d, g_hi, jmax=None, J0x=None):
    """sum_{j<J0} |c_j| d^{j-5/2} Mdiff_j with xs[j+1] >= |S^_j| (j = -1..jmax-1), explicit for j < jmax, and for j >= jmax:
    |S^_{j}| <= 2^{j-jmax+1} max(xs[jmax], xs[jmax-1]) (REC: (A+B)/den <= 1.99 up to |beta| = hypot(J0+3, 1151)+2, asserted by the
    caller; G := max(2, (A+B)/den at |beta|max) is used), Mdiff_j <= (1 + dt(A0max + A1max/d)) max(|S^_j|, |S^_{j-1}|), |c_{j+1}/c_j| <= ratio_bound_box (rt_cert)."""
    jmax = jmax or JMAX
    J0x = J0x or J0
    cs = RC.cj_list(nu_box, jmax)
    dt = g_hi*d*d
    tot = 0.0
    nu_hi = float(abs(nu_box).upper())
    for j in range(jmax):
        bmod = math.hypot(2.5 - j, nu_hi) + 1e-9
        A1, A0 = RC.Abounds(bmod)
        md = xs[j + 1]*(1 + dt*A0) + dt*A1*xs[j]/d
        tot += cs[j]*d**(j - 2.5)*md
    # REC growth constant G >= (A+B)/den for |beta| <= hypot(J0+3, nu_hi)+3 (A, B, den increasing in |beta|)
    A_, B_, den_, _ = PPm.rec_step(math.hypot(J0x + 3, nu_hi) + 3, 9, d)
    G = max(2.0, (A_ + B_)/den_)
    rb = RC.ratio_bound_box(nu_box, jmax)*d*G
    assert rb < 0.9, rb
    A1m, A0m = RC.Abounds(math.hypot(J0x + 3, nu_hi) + 3)
    Cd = 1 + dt*(A0m + A1m/d)
    Mx = G*max(xs[jmax], xs[jmax - 1])
    tot += cs[jmax]*d**(jmax - 2.5)*Cd*Mx/(1 - rb)
    return tot*(1 + 1e-10)


def a_ranges(T1):
    """cover [998, inf): unit boxes while the turning region matters, then geometric."""
    A1 = max(1000.0, math.ceil(T1/0.85))
    out = []
    a = 998.0
    while a < A1:
        w = 1.0 if a < T1 + 30 else 2.0
        out.append((a, min(a + w, A1))); a = min(a + w, A1)
    while a < 2.0e5:
        out.append((a, a*1.15)); a = a*1.15
    out.append((a, math.inf))
    return out


def K_cell(T0, T1, La, Lb, verbose=False, use_frames=True, ranges=None, jmax=None, J0x=None):
    jmax = jmax or JMAX
    d_hi = math.exp(-La); d_lo = math.exp(-Lb)
    nu_box = acb(arb((T0 + T1)/2, (T1 - T0)/2))
    Bv = RC.B_abs_real(T0)
    best = (0.0, None)
    for (a_lo, a_hi) in (ranges or a_ranges(T1)):
        M0 = PP_box(2.5, T0, T1, a_lo, a_hi, d_lo, d_hi, use_frames=use_frames)
        M1 = PP_box(3.5, T0, T1, a_lo, a_hi, d_lo, d_hi, use_frames=use_frames)
        xs = BC.chain_bounds_multi(2.5, T0, T1, M0, M1, jmax, d_lo, d_hi, a_lo, a_hi)
        _, g, _ = BC.state_box(a_lo, a_hi, d_lo, d_hi)
        v = Bv*jsum(nu_box, xs, d_hi, float(g.upper()), jmax, J0x)
        if verbose:
            print(f'    a in [{a_lo:.1f},{a_hi:.1f}]: M0={M0:.3e} M1={M1:.3e} |K| <= {v:.3e}', flush=True)
        if v > best[0]:
            best = (v, (a_lo, a_hi))
    return best


if __name__ == '__main__':
    import time
    T0, T1, La, Lb = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    t0 = time.time()
    r = K_cell(T0, T1, La, Lb, verbose=True)
    print('K_cell', T0, T1, La, Lb, r, f'{time.time()-t0:.0f}s')
    # comparison: old K_sup_real
    print('old K_sup_real', RC.K_sup_real(T0, T1, La, 2))

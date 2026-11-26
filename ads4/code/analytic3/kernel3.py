"""
kernel3.py: per-cell bound of |K_J(nu)| on the real axis (ANALYTIC_SLIVER3 §4), = rt2_kernel.K_cell (ANALYTIC_SLIVER2 §4: Theorem E',
Lemma X, Lemma F base powers, helicity for b = 4, Lemma B chain, REC tail with general growth G) with two changes:
  * the twist range is restricted to 998 <= a <= A_STAR = 4500: every state with a >= a*(L) is covered by the high-twist cell
    certificate (ANGULAR_ROUTE §73: a* = 4500 for beta_c in [3.04, 8.24], i.e. L in about [2.3, 7.3]; an imported input);
  * a finer, adaptive partition of the a-range (unit boxes near the turning point, then relative width `ratio`).
All bounds are those of rt2_kernel (imported read-only); only the cover of the a-range differs, and a cover is all the max needs.
"""
import math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'analytic2')); sys.path.insert(0, os.path.join(HERE, '..', 'analytic'))
import common2
import rt2_kernel as RK
A_STAR = 4500.0


def a_ranges(T1, ratio=1.03, a_max=A_STAR):
    out = []; a = 998.0
    while a < a_max:
        w = 1.0 if a < T1 + 30 else (2.0 if a < T1/0.85 else max(2.0, a*(ratio - 1)))
        b = min(a + w, a_max)
        out.append((a, b)); a = b
    return out


def K_cell(T0, T1, La, Lb, ratio=1.03, use_frames=True, jmax=110, J0x=1500, verbose=False):
    return RK.K_cell(T0, T1, La, Lb, verbose=verbose, use_frames=use_frames, ranges=a_ranges(T1, ratio), jmax=jmax, J0x=J0x)


if __name__ == '__main__':
    import time
    T0, La, Lb = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    for ratio in (1.15, 1.03):
        t0 = time.time()
        print(ratio, K_cell(T0, T0 + 1.1, La, Lb, ratio=ratio), f'{time.time()-t0:.0f}s', flush=True)


# ------------------------------------------------------------------ refinements (ANALYTIC_SLIVER3 §4)
import numpy as np
import bchain as BC
import rt_cert as RC
from flint import arb, acb, acb_mat


def solution_frames(sig0, tau, d, g, s, J):
    """point frames T_j = [[P_j, Q_j], [P_{j-1}, Q_{j-1}]] from the two midpoint solutions of x_{j+1} = A_j x_j + B_j x_{j-1}
    with (P_0, P_{-1}) = (1, 0), (Q_0, Q_{-1}) = (0, 1), columns rescaled to unit norm (any invertible point matrices are admissible
    in Lemma B; these make T_{j+1}^{-1} M_j(p0) T_j diagonal at the midpoint)."""
    P = [0j, 1 + 0j]; Q = [1 + 0j, 0j]       # index 0 <-> j = -1
    Ts = []
    sp = 1.0; sq = 1.0
    T0 = np.array([[1, 0], [0, 1]], dtype=complex)
    Ts.append(T0)
    for j in range(J):
        beta = complex(sig0 - j, tau)
        den = 1 - beta*(beta - s)*g*d*d
        A = (1 + beta*(s - 2*beta - 1)*g*d)/den
        B = beta*(beta + 1)*g/den
        P.append(A*P[-1] + B*P[-2]); Q.append(A*Q[-1] + B*Q[-2])
        cp = np.array([P[-1], P[-2]]); cq = np.array([Q[-1], Q[-2]])
        np_, nq = np.linalg.norm(cp), np.linalg.norm(cq)
        P[-1] /= np_; P[-2] /= np_; Q[-1] /= nq; Q[-2] /= nq      # rescale the running solutions (keeps the recurrence linear)
        Ts.append(np.column_stack([np.array([P[-1], P[-2]]), np.array([Q[-1], Q[-2]])]))
    return Ts


def chain_best(sig0, t_lo, t_hi, M0, M1, J, d_lo, d_hi, a_lo, a_hi, s_list=(1, 9)):
    """elementwise min of Lemma B over bchain's frame families and the midpoint solution frames."""
    b = BC.chain_bounds_multi(sig0, t_lo, t_hi, M0, M1, J, d_lo, d_hi, a_lo, a_hi, s_list)
    d, g, s = BC.state_box(a_lo, a_hi, d_lo, d_hi, s_list)
    tau = BC.arb_iv(t_lo, t_hi)
    Ms = []
    for j in range(J):
        beta = acb(arb(sig0 - j) if not isinstance(sig0, arb) else sig0 - j, tau)
        A, B = BC.coef(beta, d, g, s)
        Ms.append(acb_mat([[A, B], [acb(1), acb(0)]]))
    sm = float(arb(sig0).mid()) if isinstance(sig0, arb) else sig0
    Ts = solution_frames(sm, (t_lo + t_hi)/2, float(d.mid()), float(g.mid()), float(s.mid()), J)
    try:
        b2 = BC.chain_core(Ms, Ts, M0, M1, J)
        b = [min(x, y) for x, y in zip(b, b2)]
    except (ZeroDivisionError, ValueError):
        pass
    return b


def a_ranges_fine(T1, ratio=1.03, a_max=A_STAR, w_turn=0.25, turn_hi=None):
    """as a_ranges, with boxes of width w_turn on [998, turn_hi) (turning region)."""
    turn_hi = turn_hi if turn_hi is not None else T1 + 10
    out = []; a = 998.0
    while a < a_max:
        if a < turn_hi:
            w = w_turn
        else:
            w = 1.0 if a < T1 + 30 else (2.0 if a < T1/0.85 else max(2.0, a*(ratio - 1)))
        b = min(a + w, a_max)
        out.append((a, b)); a = b
    return out


def K_cell_refined(T0, T1, La, Lb, ratio=1.03, nsub_nu=1, nsub_L=1, w_turn=1.0, use_frames=True, jmax=110, J0x=1500, solframes=True):
    """sup |K_J| over nu in [T0, T1], L in [La, Lb], a in [998, 4500]: max over nu-subcells x L-subcells x a-boxes of rt2_kernel's
    per-box bound (B(nu) at the subcell's lower end, Lemma F/E'/X base, Lemma B chain incl. solution frames, REC tail)."""
    best = (0.0, None)
    for iL in range(nsub_L):
        la = La + (Lb - La)*iL/nsub_L; lb = La + (Lb - La)*(iL + 1)/nsub_L
        d_hi = math.exp(-la); d_lo = math.exp(-lb)
        for iv in range(nsub_nu):
            t0 = T0 + (T1 - T0)*iv/nsub_nu; t1 = T0 + (T1 - T0)*(iv + 1)/nsub_nu
            nu_box = acb(arb((t0 + t1)/2, (t1 - t0)/2))
            Bv = RC.B_abs_real(t0)
            ranges = a_ranges_fine(t1, ratio, w_turn=w_turn) if w_turn < 1 else a_ranges(t1, ratio)
            for (a_lo, a_hi) in ranges:
                M0 = RK.PP_box(2.5, t0, t1, a_lo, a_hi, d_lo, d_hi, use_frames=use_frames)
                M1 = RK.PP_box(3.5, t0, t1, a_lo, a_hi, d_lo, d_hi, use_frames=use_frames)
                if solframes:
                    xs = chain_best(2.5, t0, t1, M0, M1, jmax, d_lo, d_hi, a_lo, a_hi)
                else:
                    xs = BC.chain_bounds_multi(2.5, t0, t1, M0, M1, jmax, d_lo, d_hi, a_lo, a_hi)
                _, g, _ = BC.state_box(a_lo, a_hi, d_lo, d_hi)
                v = Bv*RK.jsum(nu_box, xs, d_hi, float(g.upper()), jmax, J0x)
                if v > best[0]:
                    best = (v, (a_lo, a_hi, la, t0))
    return best

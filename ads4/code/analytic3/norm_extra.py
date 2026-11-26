#!/usr/bin/env python3
"""
norm_extra.py: pieces needed by norm_cell.py below L = 3.05 (ANALYTIC_SLIVER3 §10).

(1) Exact state constants for the near-real contour pieces (rt_cert.piece_bottom / piece_left, Lemma NB / Lemma V).
    Their code uses g >= h/d^2 through `rt_pp.NB(... N_over_a=40)`, the REC constant 3.9 a^2 <= d^2 D in `rt_pp.rec_step`, and
    tau_T <= d^2/(3.9 a^2) through `rt_cert.TAU_COEF`; all three fail near L = 3.0.  The exact identity (norming README §1)
        h/d^2 = 1/[a^2 (1 + c + sqrt(4c + s^2 d^2/a^2))] <= 1/[a^2 (1 + sqrt(1 - d))^2]     (c = 1 - d)
    gives  g := 1/(998^2 (1 + sqrt(1 - d_hi))^2) >= h/d^2  and  d^2 D >= a^2 (1 + sqrt(1 - d_hi))^2  for every state of the cell
    (the right sides increase with d).  `patch_state` installs these constants (monkey-patching the imported modules in this
    process only; no file is modified).  The patched NB is the same formula with g replaced by this exact bound.
(2) The right edge nu = LAM_R + i y (0 <= y <= Y, LAM_R = 998.4998 < Lam) by complex norming boxes along the imaginary axis:
        |kappa_J - kappa_inf| <= |B(nu)| d^{-Re alpha} [(P M0 + Q M1) + (P + Q)]      (kappa_inf: the h = 0 face with X = 1),
    |B| by rt_cert.B_abs_cell (Arb), sup|psi| by assemble_top.psi_sup (disc covering the segment), M0, M1 by base_bounds.base_bound with
    the y-range; plus the real segment [LAM_R, Lam] (length < 1.3e-4) bounded with the [950, 1040.1] residual box.
"""
import sys, math, json
from fractions import Fraction as F
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str((HERE / '..' / 'norming').resolve())); sys.path.insert(1, str(HERE))
LAM_R = F('998.4998')


def patch_state(La):
    import rt_pp as PPm, rt_cert as RC
    from flint import arb
    d_hi = math.exp(-La)*(1 + 1e-12)
    cst = float((arb(1) + (1 - arb(d_hi)).sqrt())**2*(1 - arb('1e-12')))     # (1 + sqrt(1 - d_hi))^2, rounded down
    def NB(sig_lo, sig_hi, tau, s, dmax, a=998, N_over_a=None):
        assert sig_lo > 0 and dmax <= d_hi*(1 + 1e-9)
        sl = arb(sig_lo); sh = arb(sig_hi); t = arb(tau); dm = arb(dmax)
        g = 1/(arb(a)**2*arb(cst))
        Ct = (t*t*(1/(sl*sl) + 1/sl)/2).exp()
        v = Ct*(g*PPm.poch(sh, 3)/3 + g*(1 + s)*PPm.poch(sh, 2)*dm/2 + arb(4)/15*g*g*PPm.poch(sh, 5)) + arb('1e-30')
        return PPm.up(v)
    def rec_step(bmod, s, d, a=998):
        assert d <= d_hi*(1 + 1e-9)
        b = arb(bmod); A2 = arb(cst)*arb(a)**2; dd = arb(d)
        A = 1 + b*(2*b + s + 1)*dd/A2
        B = b*(b + 1)/A2
        den = 1 - b*(b + s)*dd*dd/A2
        assert den > 0.5
        Fi = b*(b + 1 + s)/A2*arb('1.000001')
        return PPm.up(A), PPm.up(B), float(den.lower()), PPm.up(Fi)
    PPm.NB = NB; PPm.rec_step = rec_step
    RC.TAU_COEF = cst*998**2
    return cst


def right_inputs(box, Y, p, width=F('0.25')):
    """norming configs along Im nu at Re nu = LAM_R."""
    from norm_cell import RHOS
    out = []; y = F(0)
    while y < Y:
        y1 = min(F(Y), y + width)
        out.append(dict(name=f'right{p}_{float(y):.2f}', **box, nu0='998.4998', dnu=str((y1 - y)/2) if (y1 - y)/2 != 0 else '0.125',
                        nuImag=str((y + y1)/2), nuAxis='imag', fullTwist=True, digits=65, K=90, channels=[[p, 0], [p, 4]],
                        epsilon='0.00001', rhos=RHOS))
        y = y1
    for c in out:          # exact decimal strings
        for k in ('dnu', 'nuImag'):
            v = F(c[k]); c[k] = format(float(v), '.6f'); assert F(c[k]) == v
    return out


def right_bound(rows, p, b, DLO, DHI, Y):
    """sum over the y-boxes of (y1 - y0) sup|psi| |B| d^{3 - Re alpha} [(P M0 + Q M1) + (P + Q)] times the contour prefactor."""
    import check_norm as CN, assemble_top as AT
    from base_bounds import base_bound, hi, hcap
    sys.path.insert(0, str(HERE / '..' / 'analytic')); import rt_cert as RC
    e0 = F(2*p - 1, 2)
    cands = []
    for r in rows:
        c = r['config']; d, dd, h, hh, nu, dnu, im = CN.parameters(c)
        if r['p'] == p and r['b'] == b and r['status'] in ('PASS', 'ENCLOSED') and c.get('nuAxis') == 'imag' and nu == LAM_R \
                and d - dd <= DLO and d + dd >= DHI:
            cands.append((im - dnu, im + dnu, r))
    pos = F(0); tot = F(0)
    pref = (4/AT.LAML if p == 3 else 2*AT.LAMU*hcap(DHI, F(998)))
    while pos < Y:
        el = [(l, u, r) for l, u, r in cands if l <= pos < u]
        if not el: return None
        l, u, r = max(el, key=lambda z: z[1]); end = min(F(Y), u)
        P, Q = F(r['P_bound']), F(r['Q_bound'])
        M0, M1 = base_bound(p, b, DLO, DHI, LAM_R, LAM_R, pos, end)
        Bb = F(repr(RC.B_abs_cell(complex(float(LAM_R), float((pos + end)/2)), float((end - pos)/2)*1.0001 + 1e-9)))
        ps = AT.psi_sup(p, LAM_R, F(999), end)
        dpow = AT.pow_upper(DHI, 3 - (e0 - pos))
        tot += (end - pos)*ps*Bb*dpow*((P*M0 + Q*M1) + (P + Q))
        pos = end
    return hi(pref*tot)

"""
loc4.py: Lemma Loc'' -- Lemma Loc' (ANALYTIC_SLIVER §7.3, rt_loc3.py) with a sharp tail ratio.  ANALYTIC_SLIVER2 §5.

Only change against rt_loc3.log_loc3: on real frequency cells [nu_lo, nu_hi] (y = 0) the tail ratio
    sup_{j>=J} |c_{j+1}/c_j| <= 1   whenever   J^2 (4 nu_lo^2 - 2 nu_hi^2 - 3/2) >= (1/4 + nu_hi^2)^2  and 4 nu_lo^2 - 2 nu_hi^2 - 3/2 > 0.
Proof: |c_{j+1}/c_j| = ((j+1/2)^2 + nu^2)/(sqrt((j+1)^2 + 4 nu^2)(j+1)) <= ((x - 1/2)^2 + nu_hi^2)/(x sqrt(x^2 + 4 nu_lo^2)), x = j+1;
it is < 1 iff P(x) := -2x^3 + x^2 (3/2 + 2nu_hi^2 - 4nu_lo^2) - 2x(1/4 + nu_hi^2) + (1/4 + nu_hi^2)^2 < 0, and for x >= J the
hypothesis makes the x^2 term dominate the constant term while the remaining terms are negative.  Hence the tail series obeys
sum_{j>=J} |c_j| r^{j-e0} <= |c_J| r^{J-e0}/(1 - r) for every r < 1, which lets the G-regions extend to |1 - xi| close to 1.
Larger (b1, th0, th1) are added to the parameter search.  Everything else (regions, weights, dressing 41, K0 <= 1.6) as in Loc'.
"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common2
import rt_loc3 as R3
from rt_loc3 import lcs, lsum_pow, lfrac, logsum

SHARP = {}


def log_loc4(p, J, L, lc, lpref, lexp, y, b1, th0, th1, dress=True, L_hi=None, nu=None):
    d = math.exp(-L); d_lo = math.exp(-(L if L_hi is None else L_hi)); om = 0.51*d; e0 = p - 0.5
    dl = 0.255*d if dress else 0.0
    lo_all = 0.49*d_lo - dl
    lyf = y*abs(math.log(lo_all))
    x_eff = (2*lexp/math.pi)
    Rt = (1 + (x_eff + y)/(J + 1))**2/(1 - 2*y/(J + 1))
    if y == 0 and nu is not None:
        nl, nh = nu
        cc = 4*nl*nl - 2*nh*nh - 1.5
        if cc > 0 and J*J*cc >= (0.25 + nh*nh)**2:
            Rt = min(Rt, 1.0)
    r0 = om + th0 + b1 + dl; r1 = om + th1 + b1 + dl
    if r1*Rt >= 0.995:
        return None
    lF = lambda lo: math.log(math.pi/2) + lexp - p*math.log(lo)
    base = lpref + lexp + lyf
    lG1 = base + lc[J] + (J - e0)*math.log(r0) - math.log(1 - r0*Rt)
    lG2 = base + lc[J] + (J - e0)*math.log(r1) - math.log(1 - r1*Rt) + lfrac(th0, d, d_lo=d_lo)
    bs = [b1]
    while bs[-1] < 0.95:
        bs.append(min(0.95, bs[-1]*1.6))
    tiers = []
    for k in range(len(bs)):
        lo = bs[k] - dl
        hi = (om + th1 + bs[k + 1] + dl) if k + 1 < len(bs) else 2 + dl
        lw = 998*math.log(1.02*d/bs[k])
        lT = base + lsum_pow(lc, J, e0, lo, hi)
        tiers.append(lw + logsum([lF(lo), lT]))
    lT = base + lsum_pow(lc, J, e0, lo_all, 2 + dl)
    lTh = logsum([lF(lo_all), lT]) + lfrac(th1, d, d_lo=d_lo)
    parts = dict(G1=lG1, G2=lG2, tiers=logsum(tiers), Th=lTh)
    tot = logsum(list(parts.values()))
    return math.log(1.6) + (math.log(41) if dress else 0.0) + tot + 1e-9, parts


PARAMS = [(b1, th0, th1) for b1 in (0.04, 0.06, 0.08, 0.11, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
          for th0 in (0.02, 0.04, 0.08) for th1 in (0.1, 0.15, 0.2, 0.3, 0.4, 0.5)]


def best_bound(p, J, L, lc, lpref, lexp, y, dress, L_hi=None, nu=None):
    best = None
    for (b1, th0, th1) in PARAMS:
        if th0 >= th1: continue
        r = log_loc4(p, J, L, lc, lpref, lexp, y, b1, th0, th1, dress, L_hi, nu)
        if r is not None and (best is None or r[0] < best[0]):
            best = (r[0], (b1, th0, th1), r[1])
    return best


def cell_bound(p, J, L, nu_lo, nu_hi, dress, L_hi=None):
    if nu_hi <= 0.4 + 1e-12:
        return R3.cell_bound(p, J, L, nu_lo, nu_hi, dress, L_hi)
    import rt_cert as RC
    lc = lcs(nu_lo, nu_hi, 0.0, J + 1)
    lpref = math.log(2*RC.B_abs_real(nu_lo))
    lexp = math.pi*nu_hi/2
    b4 = best_bound(p, J, L, lc, lpref, lexp, 0.0, dress, L_hi, (nu_lo, nu_hi))
    b3 = R3.best_bound(p, J, L, lc, lpref, lexp, 0.0, dress, L_hi)
    c = [b for b in (b3, b4) if b is not None]
    return min(c, key=lambda b: b[0])


if __name__ == '__main__':
    L = float(sys.argv[1]); J = int(sys.argv[2])
    for nu in (1.0, 500.0, 900.0, 1000.0, 1029.0, 1040.0, 1146.0, 1300.0, 1600.0):
        b = cell_bound(3, J, L, nu, nu + 1, True, L + 0.5)
        print(f'L={L} J={J} nu={nu}: log|V_b R| <= {b[0]:.1f} params {b[1]} parts ' + ', '.join(f'{k}: {v:.0f}' for k, v in b[2].items()), flush=True)

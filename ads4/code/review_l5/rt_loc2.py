"""
rt_loc2.py: Lemma Loc' (ANALYTIC_SLIVER §7.2) -- localised bound for the j-truncation remainder of the connection series.

R = R_J(.; nu) = F_{p,nu} - sum_{+-} B(+-nu) sum_{j<J} c_j(+-nu) (1-xi)^{j - e0 -+ i nu}   (e0 = p - 1/2), analytic in |xi| < 1.
For every state (TAIL_CERT §4.1):  L(R) = E_B[ (1/2 pi i) oint_{|z|=rho} (1-z)^{-a-1} R(Bz) z^{-n-1} dz ] / C(N, a),  B ~ Beta(n+s, a).
With rho = n/N, W(theta) = |1 - rho e^{i theta}|^{-a-1}, I_W = (1/2pi) int W, K0 = rho^{-n} I_W / C(N, a) <= 1.6 (Stirling/Robbins for C(N,a),
TAIL's bound on I_W, Gautschi), and the regions (xi = B rho e^{i theta}):
   G1 = {B >= 1-b1, |theta| <= th0}:        |1 - xi| <= r0 := (1-rho) + th0 + b1
   G2 = {B >= 1-b1, th0 < |theta| <= th1}:   |1 - xi| <= r1 := (1-rho) + th1 + b1,   angular fraction <= frac(th0)
   Bd = {B < 1-b1, |theta| <= th1}:          b1 <= |1 - xi| <= 1 (th1 <= 1),            probability <= P_B = (1.02 d/b1)^{998}
   Th = {|theta| > th1}:                     1-rho <= |1 - xi| <= 2,                    angular fraction <= frac(th1)
  frac(th) <= pi sqrt(a) e^{0.51} (1-rho)^{-1} ((1-rho)/q(th))^{a+1},  q(th) = sqrt((1-rho)^2 + 4 rho sin^2(th/2))   (I_W >= e^{-0.51}(1-rho)^{-a}/(pi sqrt a))
  |arg(1 - xi)| <= pi/2 everywhere, so |(1-xi)^{-i nu}| <= e^{pi |nu|/2}.
  On G1, G2 the tail representation R = sum_{+-} B sum_{j>=J} c_j (1-xi)^{j-alpha} converges (r * sup_{j>=J}|c_{j+1}/c_j| < 1).
  On Bd, Th: |R| <= |F| + |T_J|,  |F_{p,nu}(xi)| <= (pi/2) e^{pi|nu|/2} |1-xi|^{-p} (FULLPOS3 Lemma B), |T_J| <= 2|B| e^{pi|nu|/2} sum_{j<J} |c_j| |1-xi|^{j-e0}.
  |L(R)| <= K0 [ sup_G1 |R| + frac(th0) sup_G2 |R| + P_B sup_Bd |R| + frac(th1) sup_Th |R| ].
The dressing tau D_b R (tau <= d^2/(3.9 a^2)) is bounded by Cauchy on circles of radius 0.01 (factor <= 1 + 1e-6 here).
Worst case a = 998 (all fractions decrease in a; 1 - rho in [0.49 d, 0.51 d] for N >= 40 a).
"""
import math
from flint import arb, acb, ctx
ctx.prec = 128
I = acb(0, 1)


def logsum(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def cj_abs_list(nu_box, J):
    """upper bounds of |c_j(nu)|, j = 0..J+1, uniform over real nu in [nu_lo, nu_hi] (the box's real range; Im must be 0 for this form):
    |c_j| = prod_{i<j} |a0+i|^2/(|b0+i|(i+1)),  |a0+i|^2 = (i+1/2)^2 + nu^2 <= (i+1/2)^2 + nu_hi^2,  |b0+i| = sqrt((i+1)^2 + 4nu^2) >= sqrt((i+1)^2 + 4nu_lo^2)."""
    lo = max(0.0, float(nu_box.real.lower())); hi = float(nu_box.real.upper())
    assert float(abs(nu_box.imag).upper()) == 0.0
    out = [arb(1)]
    c = arb(1); h2 = arb(hi)**2; l2 = arb(lo)**2
    for i in range(J + 1):
        c = c*((arb(i) + arb(0.5))**2 + h2)/(((arb(i + 1))**2 + 4*l2).sqrt()*(i + 1))
        out.append(arb(c.upper()))
    return out


def log_loc2(nu_box, p, J, L, b1, th0, th1, Babs):
    """log of the Lemma Loc' bound of |V_b(R_J)| for nu in the (real or complex) box, every state with -log d >= L, a >= 998."""
    d = math.exp(-L)
    om = 0.51*d                       # 1 - rho <= 0.51 d
    a = 998
    e0 = p - 0.5
    nmax = float(abs(nu_box).upper())
    lexp = math.pi*nmax/2
    cs = cj_abs_list(nu_box, J)
    lB = math.log(Babs)
    # tail ratio sup_{j >= J} |c_{j+1}/c_j| <= (1 + |nu|/(J+1))^2 (j >= J >= |nu|), with |b0 + j| >= j + 1 (Im nu >= -0.41 handled by the margin)
    Rt = (1 + (nmax + 0.5)/(J + 0.2))**2
    r0 = om + th0 + b1; r1 = om + th1 + b1
    assert r1*Rt < 0.95, (r1, Rt)
    lcJ = float(cs[J].log().upper())
    lG1 = math.log(2) + lB + lexp + lcJ + (J - e0)*math.log(r0) - math.log(1 - r0*Rt)
    lG2 = math.log(2) + lB + lexp + lcJ + (J - e0)*math.log(r1) - math.log(1 - r1*Rt)
    def lfrac(th):
        q = math.sqrt(om*om + 4*(1 - om)*math.sin(th/2)**2)
        return math.log(math.pi*math.sqrt(a)) + 0.51 - math.log(om*0.49/0.51) + (a + 1)*math.log(om/q)
    # Bd: |1-xi| in [b1, 1]
    lsumc = float(sum(cs[:J], arb(0)).log().upper())
    lT_Bd = math.log(2) + lB + lexp + lsumc + max(0.0, -e0*math.log(b1))
    lF_Bd = math.log(math.pi/2) + lexp - p*math.log(b1)
    lBd = logsum([lT_Bd, lF_Bd]) + 998*math.log(1.02*d/b1)
    # Th: |1-xi| in [1-rho, 2]
    lsum2 = float(sum((cs[j]*arb(2)**j for j in range(J)), arb(0)).log().upper())
    lT_Th = math.log(2) + lB + lexp + lsum2 + e0*math.log(1/(0.49*d))
    lF_Th = math.log(math.pi/2) + lexp - p*math.log(0.49*d)
    lTh = logsum([lT_Th, lF_Th]) + lfrac(th1)
    tot = logsum([lG1, lG2 + lfrac(th0), lBd, lTh])
    return math.log(1.6) + tot + 1e-6, dict(G1=lG1, G2=lG2 + lfrac(th0), Bd=lBd, Th=lTh)


if __name__ == '__main__':
    import sys
    L = float(sys.argv[1]) if len(sys.argv) > 1 else 6
    for nu in (1.0, 100.0, 500.0, 900.0, 1146.0):
        box = acb(arb(nu, 0.5))
        Bv = float((((arb.pi()*arb(max(nu - 0.5, 0.4))).coth()/(4*arb.pi()*arb(max(nu - 0.5, 0.4)))).sqrt()).upper())
        best = None
        for J in (800, 1200, 1600, 2000):
            for b1 in (0.05, 0.08, 0.12):
                for th1 in (0.1, 0.15, 0.2):
                    try:
                        v, parts = log_loc2(box, 3, J, L, b1, 0.02, th1, Bv)
                    except AssertionError:
                        continue
                    if best is None or v < best[0]:
                        best = (v, J, b1, th1, parts)
        print(f'L={L} nu={nu}: best log bound {best[0]:.1f} (J={best[1]}, b1={best[2]}, th1={best[3]})  parts ' +
              ', '.join(f'{k}: {v:.0f}' for k, v in best[4].items()), flush=True)

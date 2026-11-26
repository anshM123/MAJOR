"""
crossing.py: Lemma X (turning-point crossing) for the b = 0 per-power multiplier, ANALYTIC_SLIVER §5.5.
Recurrence S_{k+1} = 2 w_k S_k + S_{k-1}, w_k = eps_k alpha/2 = u_k + i y_k (alpha = sigma + i tau, sigma >= 0), y_k increasing, y_{n-1} < tau/(a+1).
Phases: I (y <= yA = 1 - del), II (window yA < y < yB = 1 + del'), III (y >= yB).
  I   : eigen-coordinates (Theorem E' argument) up to y = yA: factor 2.01 c0 exp( sigma log(2/(1+sqrt(1-yA'^2))) + q[asin yA' + (q/2) log(1/(1-yA'^2))] + q f(yA') dy_max )
        (yA' = yA + 2tau/a^2, dy_max = tau/a^2, f(y) = (1-y^2)^{-1/2} + q y/(1-y^2); excess via Delta eps_k >= eps_k^2/2), giving |v| at the window entry.
  II  : T_k = i^{-k} S_k, T_{k+1} = 2 mu_k T_k - T_{k-1}, mu = y - i u; with D = T_k - T_{k-1} and eta >= 2 max|mu - 1|, the norm
        |T| + |D|/sqrt(eta) grows by <= 1 + sqrt(eta) + eta per step; W <= tau (1/yA - 1/yB) + 2 steps:  factor 2 (1 + 2/sqrt(eta)) e^{W (sqrt(eta)+eta)}.
  III : eigen-coordinates with |lam+| = e^{acosh y}: ||P^{-1}|| <= (1 + |lam|max)/(2 sqrt(yB^2 - 1 - u^2)); sum acosh(y_k) <= G_a(tau) =
        tau acos(a/tau) - a acosh(tau/a); adiabatic <= q[acosh ymax - acosh yB + (q/2) log((ymax^2-1)/(yB^2-1))]; sigma-terms
        sigma log(tau/a)/sqrt(yB^2-1-u^2); final factor (|lam+|+|lam-|)|lam+| <= 2 lmax^2.
The d^sigma normalisation is used only as d_n^sigma <= prod e^{-sigma eps_k/2} <= 1 on phases II-III (sigma >= 0).
Uniform version (all a >= 998): G at a = 998, ymax = tau/999, all phases present.
"""
import math
def log_cross(sig, tau, a=998.0, c=1.0, uniform=True):
    q = math.hypot(sig, tau)/tau
    dl = c*tau**(-2/3); dlp = dl
    yA = 1 - dl; yB = 1 + dlp
    ymax = tau/(a + 1)
    umax = sig/(a + 1)
    # phase I
    # phase-I excess: sum_k (sigma eps_k/2)((1-y_k^2)^{-1/2} - 1) <= sigma int_0^{yA'} ((1-y^2)^{-1/2} - 1) dy/y = sigma log(2/(1+sqrt(1-yA'^2)))
    # (exact inequality  eps_{k+1} - eps_k >= eps_k^2/2,  i.e.  Delta y_k >= y_k^2/tau;  yA' = yA + 2 tau/a^2 covers the last step)
    yAp = min(yA + 2*tau/a**2, 0.999999)
    fA = (1 - yAp*yAp)**-0.5 + q*yAp/(1 - yAp*yAp)
    lI = math.log(2.01*(1 + 2*math.hypot(sig, tau)/a**2)) + sig*math.log(2/(1 + math.sqrt(1 - yAp*yAp))) \
        + q*(math.asin(yAp) + q/2*math.log(1/(1 - yAp*yAp))) + q*fA*tau/a**2
    if ymax <= yA and not uniform:
        return lI
    # phase II: entry N0 <= (2 + e1)|c| + e1 |c|/sqrt(eta),  e1 = 2(u + 1 - yA + dy + sqrt(1 - (yA-dy)^2 + u^2 + 2u))  (|lam± - i| <= e1/2)
    eta = 2*(max(dl, dlp) + umax)
    dyx = 2*tau/a**2
    e1 = 2*(umax + 1 - yA + dyx + math.sqrt(max(0.0, 1 - (yA - dyx)**2) + umax**2 + 2*umax))
    W = tau*(1/yA - 1/yB) + 3
    lII = math.log((2 + e1) + e1/math.sqrt(eta)) - math.log(2.01) + W*(math.sqrt(eta) + eta)
    if ymax <= yB and not uniform:
        return lI + lII
    # phase III entry: |c''| <= N_W [sqrt(eta)(1 + e2) + e2]/(2 rB),  e2 = u + yB + dy - 1 + sqrt((yB+dy)^2 - 1 + u^2 + 2u(yB+dy))
    ym = max(ymax, yB*1.000001)
    lmax = ym + math.sqrt(ym*ym - 1) + umax*2
    rB = math.sqrt(yB*yB - 1 - umax*umax)
    e2 = umax + yB + dyx - 1 + math.sqrt((yB + dyx)**2 - 1 + umax**2 + 2*umax*(yB + dyx))
    lP = math.log((math.sqrt(eta)*(1 + e2) + e2)/(2*rB))
    G = tau*math.acos(min(1.0, a/tau)) - a*math.acosh(max(1.0, tau/a))
    lad = q*(math.acosh(ym) - math.acosh(yB) + q/2*math.log((ym*ym - 1)/(yB*yB - 1)))
    lsig = sig*math.log(max(tau/a, 1.0) + 1e-12)/rB + sig/(a*rB)
    lfin = math.log(2*lmax*lmax)
    return lI + lII + lP + G + lad + lsig + lfin
if __name__ == '__main__':
    import sys
    from flint import arb, acb, ctx
    ctx.prec = 8000
    for (nf, a) in [(1001, 998), (5000, 998), (5000, 1020), (20000, 998), (3000, 1050)]:
        for sig in (0.25, 2.5):
            for tau in (960, 990, 998, 1005, 1020, 1040, 1080):
                N = arb(nf) + arb(a); al = acb(sig, tau)
                Sm, S = acb(1), acb(1)
                for m in range(nf):
                    eps = 2*(2*m + 1)/((N - m)*(N + m + 1)); Sm, S = S, Sm + eps*al*S
                d = 1 - arb(nf*(nf + 1))/(N*(N + 1))
                ex = float(abs(d**al*S).mid())
                G = tau*math.acos(min(1, a/tau)) - a*math.acosh(max(1, tau/a))
                bs = min(log_cross(sig, tau, a, c) for c in (0.5, 1, 2, 4))
                bu = min(log_cross(sig, tau, 998.0, c) for c in (0.5, 1, 2, 4))
                print(f'n={nf} a={a} alpha={sig}+{tau}i: log|S^| = {math.log(ex):.3f}  G_a = {G:.3f}  log bound (this a) = {bs:.2f}  (uniform a>=998) = {bu:.2f}', flush=True)

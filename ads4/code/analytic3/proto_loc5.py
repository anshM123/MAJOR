"""float prototype (DV only) of the 2D (u = 1-B, theta) cell remainder bound with Chernoff B-weights and phase-aware arg bounds."""
import math, sys, cmath
sys.path.insert(0, '../analytic2'); sys.path.insert(0, '../analytic')
import rt_loc3 as R3
import rt_cert as RC

def logsum(xs):
    xs=[x for x in xs if x>-1e300]; m=max(xs); return m+math.log(sum(math.exp(x-m) for x in xs))

def lchern(u, t, a=998):
    # log P(1-B >= u) <= min_m sum_{i<m} log((a+i)/((N+i) u)), N = a/t ; integral form a*int_0^mu log((1+x)/((K+x)u)) dx, K = 1/t
    if u <= t: return 0.0
    K = 1/t
    xs = (K*u - 1)/(1 - u)  # zero of the integrand
    # int_0^X log(1+x) - log(K+x) - log u
    F = lambda X: ((1+X)*math.log(1+X) - X) - ((K+X)*math.log(K+X) - K*math.log(K) - X) - X*math.log(u)
    return min(0.0, a*F(max(0.0, xs - 1.0/a)))

def lfrac(th, t_hi, t_lo, a=998):
    if th <= 0: return 0.0
    q = math.sqrt(t_hi**2 + 4*(1-t_hi)*math.sin(th/2)**2)
    return min(0.0, math.log(math.pi*math.sqrt(a)) + 0.51 - math.log(t_lo) + (a+1)*math.log(t_hi/q))

def cellinfo(u0, u1, th0, th1, t_lo, t_hi, dl):
    # sample-based (DV) extremes of |1-xi|, |arg(1-xi)|, |Im eta| over the cell, enlarged by dl
    los=[]; his=[]; A=[]; E=[]
    for uu in (u0, u1, (u0+u1)/2):
        for th in (th0, th1, (th0+th1)/2):
            for rr in (1-t_lo, 1-t_hi):
                xi = (1-uu)*rr*cmath.exp(1j*th)
                w = 1 - xi
                los.append(abs(w)); his.append(abs(w))
                A.append(abs(cmath.phase(w)) + math.asin(min(1, dl/abs(w))))
                s = cmath.sqrt(xi); eta = cmath.log((1+s)/(1-s))
                E.append(abs(eta.imag) + 2*dl/abs(w))
    return max(0.49*0, min(los) - dl), max(his) + dl, min(max(A), math.pi/2), min(max(E), math.pi/2)

def bound(p, J, L, nu, dress=True, nu_w=1.0):
    d = math.exp(-L); t_hi = d*(0.5 + d/7); t_lo = d*0.5*(1-1e-3)
    dl = 0.5*t_hi if dress else 0.0
    e0 = p - 0.5
    lc = R3.lcs(nu, nu+nu_w, 0.0, J+1)
    lB = math.log(2*RC.B_abs_real(nu))
    nuh = nu + nu_w
    ubr = [0.0] + [t_hi*(1.25**k) for k in range(0, 200) if t_hi*1.25**k < 0.999] + [1.0]
    thb = [0.0] + [t_hi*(1.25**k) for k in range(0, 200) if t_hi*1.25**k < math.pi] + [math.pi]
    terms = []; worst = (-1e9,)
    for i in range(len(ubr)-1):
        u0, u1 = ubr[i], ubr[i+1]
        wB = lchern(u0, t_hi)
        for k in range(len(thb)-1):
            th0, th1 = thb[k], thb[k+1]
            wT = lfrac(th0, t_hi, t_lo)
            if wB + wT < -5000: continue
            lo, hi, A, Ei = cellinfo(u0, u1, th0, th1, t_lo, t_hi, dl)
            cands = []
            # tail series
            if hi < 0.995:
                cands.append(lB + nuh*A + lc[J] + (J-e0)*math.log(hi) - math.log(1-hi))
            # F + T_J
            lF = math.log(math.pi/2) + nuh*Ei - p*math.log(lo)
            lT = lB + nuh*A + logsum([lc[j] + max((j-e0)*math.log(lo), (j-e0)*math.log(hi)) for j in range(J)])
            cands.append(logsum([lF, lT]))
            v = wB + wT + min(cands)
            terms.append(v)
            if v > worst[0]: worst = (v, u0, u1, th0, th1)
    tot = logsum(terms) + math.log(1.6) + (math.log(41) if dress else 0)
    return tot, worst

if __name__ == '__main__':
    for L in [float(x) for x in sys.argv[1].split(',')]:
        for J in [int(x) for x in sys.argv[2].split(',')]:
            for nu in (300., 600., 900., 990., 1029., 1100.):
                b = bound(3, J, L, nu)
                print(L, J, nu, round(b[0],1), [round(x,3) for x in b[1]], flush=True)

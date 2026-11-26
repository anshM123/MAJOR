"""
rt_pp.py: certified per-power multiplier bounds for Route T' (ANALYTIC_SLIVER §7).  All bounds hold for every state with a >= 998,
spin n with N/a >= 40 (true for L >= 7), both channels b in {0, 4} (s = 1, 9), uniformly; they are evaluated in Arb and returned as
floats rounded upward (.upper()).

Notation: S^(beta) = d^beta S(beta), S(beta) = L((1-z)^{-beta}), d = 1 - C0/D.  dmax = upper bound of d on the L-cell.

(NB)  Lemma NB (TAIL_CERT §2.1, re-derived with explicit d): for Re beta = sigma > 0, |Im beta| <= tau,
      |S^ - 1| <= C_tau [ g (sigma)_3/3 + g (1+s)(sigma)_2 dmax/2 + (4/15) g^2 (sigma)_5 ] + 1e-30,
      g = 0.25 (1 + s/N)/(1 - dmax/2)^2 / a^2 >= h/d^2,  C_tau = exp(tau^2 (1/sigma^2 + 1/sigma)/2)  (|(beta)_k| <= C_tau (sigma)_k).
      (1 - E_k <= (d/c) h S1(k) + (4/3)(d/c^2) h^2 S2(k) for k <= n/2, NB factorial moments, Markov for k > n/2.)
(REC) exact three-term recurrence (TAIL E3) in the d-normalisation:
      [D - beta(beta-s)] S^(beta-1) = [D + beta(s-2beta-1)/d] S^(beta) + (beta(beta+1)/d^2) S^(beta+1),  with d^2 D >= 3.9 a^2 (N >= 40a);
      hence |S^(beta-1)| <= [ (1 + |b|(2|b|+s+1) d/(3.9a^2)) |S^(beta)| + |b|(|b|+1)/(3.9a^2) |S^(beta+1)| ] / (1 - |b|(|b|+s) d^2/(3.9a^2)),
      and for Delta = S^ - 1 the same with the inhomogeneity R(beta)/[D - ...], |R|/D <= |b|(|b|+1+s)/(3.9a^2)  (Lemma V).
(E')  Theorem E' (b = 0), closed form, tau <= 0.95 a:  2 c0 e^{sigma/((a+1)sq)} X^sigma e^{q asin k} (1-k^2)^{-q^2/2}.
(X)   Lemma X (turning-point crossing, b = 0), uniform over a >= 998: crossing.py (here in Arb).
b = 4 via the helicity identity S^(9)_{n,N}(beta) = sum_i C(4,i)(-1)^i S^(1)_{n+4,N+4}(beta-i)/U_4  (same a).
"""
import math
from flint import arb, acb, ctx
ctx.prec = 128
A_MIN = 998
KMAX = 0.95


def up(x):
    return float(arb(x).upper())


def poch(x, k):
    v = arb(1)
    for i in range(k):
        v *= (x + i)
    return v


def NB(sig_lo, sig_hi, tau, s, dmax, a=A_MIN, N_over_a=40):
    """upper bound of |S^(beta) - 1| for Re beta in [sig_lo, sig_hi] (sig_lo > 0), |Im beta| <= tau."""
    assert sig_lo > 0
    sl = arb(sig_lo); sh = arb(sig_hi); t = arb(tau); dm = arb(dmax)
    g = arb('0.25')*(1 + arb(s)/(arb(a)*N_over_a))/(1 - dm/2)**2/arb(a)**2
    Ct = (t*t*(1/(sl*sl) + 1/sl)/2).exp()
    v = Ct*(g*poch(sh, 3)/3 + g*(1 + s)*poch(sh, 2)*dm/2 + arb(4)/15*g*g*poch(sh, 5)) + arb('1e-30')
    return up(v)


def rec_step(bmod, s, d, a=A_MIN):
    """coefficients (A, B, den) of |S^(beta-1)| <= (A|S^(beta)| + B|S^(beta+1)|)/den, |beta| <= bmod."""
    b = arb(bmod); A2 = arb('3.9')*arb(a)**2; dd = arb(d)
    A = 1 + b*(2*b + s + 1)*dd/A2
    B = b*(b + 1)/A2
    den = 1 - b*(b + s)*dd*dd/A2
    assert den > 0.5
    F = b*(b + 1 + s)/A2*arb('1.000001')          # inhomogeneity for Delta (Lemma V)
    return up(A), up(B), float(den.lower()), up(F)


def thmE(sigma, t_lo, t_hi, a=A_MIN):
    """Theorem E' over sigma' <= sigma, tau in [t_lo, t_hi] (b = 0), or None if t_hi > 0.95 a.  Increasing in sigma, kappa, q, |alpha|."""
    if t_hi > KMAX*a:
        return None
    s = arb(sigma); aa = arb(a); k = arb(t_hi)/aa
    al = (s*s + arb(t_hi)**2).sqrt()
    sq = (1 - k*k).sqrt(); X = 2/(1 + sq)
    if t_lo > 1e-3:
        q = al/arb(t_lo)
        e1 = q*(k.asin() + q*(-(1 - k*k).log()/2))
    else:
        e1 = al/aa*(k.asin()/k if t_hi > 0 else arb(1)) + (al/aa)**2/(2*sq*sq)
    e2 = s/((aa + 1)*sq)
    c0 = 1 + 2*al/(aa*aa)
    return up(2*c0*X**s*(e1 + e2).exp())


def cross(sigma, t_lo, t_hi, c=1.0, a=A_MIN):
    """Lemma X (Arb), uniform over all a >= 998: bound of |S^| for Re beta <= sigma (sigma >= 0), tau in [t_lo, t_hi], t_hi <= 1150."""
    s = arb(sigma); tl = arb(t_lo); th = arb(t_hi); aa = arb(a)
    al = (s*s + th*th).sqrt(); q = al/tl
    dl = arb(c)*tl**(arb(-2)/3)
    yA = 1 - dl; yB = 1 + dl
    umax = s/(aa + 1)
    yAp = yA + 2*th/aa**2
    assert yAp < 1
    fA = 1/(1 - yAp*yAp).sqrt() + q*yAp/(1 - yAp*yAp)
    lI = (arb('2.01')*(1 + 2*al/aa**2)).log() + s*(2/(1 + (1 - yAp*yAp).sqrt())).log() \
        + q*(yAp.asin() + q/2*(1/(1 - yAp*yAp)).log()) + q*fA*th/aa**2
    eta = 2*(dl + umax)
    dyx = 2*th/aa**2
    e1 = 2*(umax + 1 - yA + dyx + ((1 - (yA - dyx)**2) + umax*umax + 2*umax).sqrt())
    W = th*(1/yA - 1/yB) + 3
    lII = ((2 + e1) + e1/eta.sqrt()).log() - arb('2.01').log() + W*(eta.sqrt() + eta)
    ymax = th/(aa + 1)
    ym = ymax if ymax > yB*arb('1.000001') else yB*arb('1.000001')
    lmax = ym + (ym*ym - 1).sqrt() + 2*umax
    rB = (yB*yB - 1 - umax*umax).sqrt()
    e2 = umax + yB + dyx - 1 + ((yB + dyx)**2 - 1 + umax*umax + 2*umax*(yB + dyx)).sqrt()
    lP = ((eta.sqrt()*(1 + e2) + e2)/(2*rB)).log()
    r_ = th/aa
    if r_ > 1:
        G = th*(aa/th).acos() - aa*r_.acosh()
    else:
        G = arb(0)
    lad = q*(ym.acosh() - yB.acosh() + q/2*((ym*ym - 1)/(yB*yB - 1)).log())
    lsig = s*(r_.log() if r_ > 1 else arb(0))/rB + s/(aa*rB)
    lfin = (2*lmax*lmax).log()
    return up((lI + lII + lP + G + lad + lsig + lfin).exp())


def PP0(sigma, t_lo, t_hi):
    """b = 0 bound of |S^(sigma' + i tau)| for 0 <= sigma' <= sigma, tau in [t_lo, t_hi] (0 <= t_lo <= t_hi <= 1150)."""
    e = thmE(sigma, t_lo, t_hi)
    if e is not None:
        return e
    lo = max(t_lo, 900.0)
    return min(cross(sigma, lo, t_hi, c) for c in (0.7, 1.0, 1.5))


def PP(sigma, t_lo, t_hi, d0):
    """both channels (b = 0 directly; b = 4 by the helicity identity), 0 <= sigma' <= sigma, tau in [t_lo, t_hi].
    d0 >= d of the shifted s = 1 state (<= dmax (1 + 1e-4))."""
    p0 = PP0(sigma, t_lo, t_hi)
    p1 = PP0(sigma + 1, t_lo, t_hi)
    m = max(p0, p1)
    corr = sum(math.comb(4, i)*d0**i*2.0**i*m for i in range(1, 5))
    p4 = (1 + 1e-5)*(p0 + corr)/((1 - d0)**4*(1 - 1e-6))
    return max(p0, p4)


if __name__ == '__main__':
    for (sg, a_, b_) in ((0.25, 0, 1), (0.25, 300, 301), (2.5, 900, 901), (2.5, 940, 948), (2.5, 990, 998.5), (1.25, 1030, 1035), (2.5, 1100, 1146)):
        print(f'PP(sigma<={sg}, tau in [{a_},{b_}]) <= {PP(sg, a_, b_, math.exp(-7)):.4g}')
    print('NB(2.5, 2.5, 0.25, s=9, d=e^-7) =', NB(2.5, 2.5, 0.25, 9, math.exp(-7)))
    print('rec_step(1000, 9, e^-7) =', rec_step(1000, 9, math.exp(-7)))

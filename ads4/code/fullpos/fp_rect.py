"""
fullpos/fp_rect.py (FULLPOS_CERT.md §3.2): every spin J >= J0 at bounded impact parameter, as ONE bivariate certificate
(ANGULAR_ROUTE §34 variables, own implementation, own moments).
  x = bhat^2/4 = Lam^2 (J+1/2)^2/m^2,  v = (J+1/2)^{-2},  m^2 = T - C,  lam = Lam^2:
  C = 1/v - 1/4,  T = E/(x v) with E = lam + x - x v/4,  1/T = x v/E,  r = C/T = (1 - v/4) x/E,
  U_{b,k} = x^k A_k(v)/(lam^k Pi_k),  A_k = prod_{j<k}(1 - (b+j+1/2)^2 v),  Pi_k = prod_{j<k}(1 + (x/lam)(1 - (b+j+1/2)^2 v))
  (exact rewriting of prod_{j<k}(C-u_j)/(T-u_j) with u_j = (b+j)(b+j+1) = (b+j+1/2)^2 - 1/4).
H_ang,b/(1-r)^3 = sum_k D_k^2 U_k [Z_k(1/(k+3) - r/(k+1)) + (2Y_k + V_k Z_k + sigma lam Q_k)/T]
  = (1/(E Pi_K)) * N(x,v) + tail,   N = sum_{k<=K} D_k^2 x^k A_k (Pi_K/Pi_k) [E Zt_k/(k+3) - (1-v/4) x Zt_k/(k+1) + x v (2 lam Yt_k + V_k Zt_k + sigma lam Qt_k)]
with Zt = Z/lam^k, Yt = Y/lam^{k+1}, Qt = Q/lam^k.  For J >= J0 > K + b every k <= K is a genuine ladder term (A_k > 0, Pi_k >= 1), and
  |tail| <= sum_{k>K} D_k^2 x^k [|Zt_k|(1/(k+3) + x/(lam(k+1))) + (x v/lam)(2 lam |Yt_k| + V |Zt_k| + lam |Qt_k|)]
(0 <= A_k <= 1, Pi_k >= 1, 1/T <= x v/lam, r <= x/lam), with certified moments for K < k <= Kc and, for k > Kc, the crude all-k bound
  |f_m(nu)| <= 3^{m-1} exp((pi/sqrt3)(nu + 1/2))   (three-term recurrence, growth factor <= 3 + s/(m+1)^2; prod(1+a/j^2) <= e^{pi sqrt a}),
  |q_n(nu^2)| <= 3^n exp((pi/sqrt3) sqrt(nu^2 + 17/4)),
integrated against |g| <= (1/Lam)[84 sum|x_k| G_k + 189 |L|] (Gaussian exponential moments in closed form).
Certificate: N(x,v) - TB*E*Pi_K > 0 on [0, xmax] x [0, vmax] via tensor Bernstein coefficients (+ bisection in x).   [certified]
Run: python3 fp_rect.py Dg Kmom J0 xmax K
"""
import sys, os, math, time
sys.dont_write_bytecode = True
from flint import arb, arb_poly, fmpq, ctx
from math import comb
import fp_common as FC
import fp_poly as FP


def bern1(coeffs, a, b):
    p = arb_poly(coeffs)
    n = len(coeffs) - 1
    if n == 0:
        return [coeffs[0]]
    return FP.bernstein(p, a, b) if p.degree() == n else FP.bernstein(arb_poly(list(coeffs)), a, b)[:n + 1] + [arb(0)]*0


def bern1_fixed(coeffs, a, b):
    """Bernstein coefficients of degree len-1 (pads if the top coefficients vanish)."""
    n = len(coeffs) - 1
    q = arb_poly(coeffs)(arb_poly([a, b - a]))
    c = [q[i] if i <= q.degree() else arb(0) for i in range(n + 1)]
    fac = [arb(1)]
    for i in range(1, n + 1):
        fac.append(fac[-1]*i)
    Ap = arb_poly([c[j]*fac[n - j]/fac[n] for j in range(n + 1)])
    Ep = arb_poly([1/fac[m] for m in range(n + 1)])
    pr = Ap*Ep
    return [(pr[i] if i <= pr.degree() else arb(0))*fac[i] for i in range(n + 1)]


class Bi:
    """bivariate polynomial: dict (i, j) -> arb coefficient of x^i v^j."""
    def __init__(self, d=None):
        self.d = d or {}

    def __add__(self, o):
        d = dict(self.d)
        for k, c in o.d.items():
            d[k] = d.get(k, arb(0)) + c
        return Bi(d)

    def __mul__(self, o):
        if not isinstance(o, Bi):
            return Bi({k: c*o for k, c in self.d.items()})
        d = {}
        for (i, j), c in self.d.items():
            for (k, l), e in o.d.items():
                key = (i + k, j + l)
                d[key] = d.get(key, arb(0)) + c*e
        return Bi(d)


def crude_tail_bounds(Dg, lam, Lam):
    """returns (G2, Gq): int |g| e^{(pi/sqrt3)(nu+1/2)} dnu and int_0^1 |psi3| e^{(pi/sqrt3)sqrt(lam p^2+17/4)} dp bounds."""
    W = FC.load_weight(Dg)
    a = arb.pi()/arb(3).sqrt()
    tot = arb(0)
    for (c, s), x in zip(W['cs'], W['x']):
        ca = arb(fmpq(c.numerator, c.denominator)); sa = arb(fmpq(s.numerator, s.denominator))
        xa = abs(arb(fmpq(x.numerator, x.denominator)))
        mult = 1 if c == 0 else 2                        # +c and -c members; int over R of G(nu -+ c) e^{a nu} <= e^{a c + a^2 s^2/2} s sqrt(2pi)
        tot += mult*xa*sa*(2*arb.pi()).sqrt()*(a*ca + a*a*sa*sa/2).exp()
    tot = 84*tot
    # log part: |L| <= sum|y_j| nu^{2j}(nu^2+1) log(1+nu^2) e^{-nu^2/2} <= sum|y_j| nu^{2j}(nu^2+1)^2 e^{-nu^2/2}; times e^{a nu}: crude closed bound
    ys = [abs(arb(fmpq(y.numerator, y.denominator))) for y in W['y']]
    # int_0^oo nu^{2j}(nu^2+1)^2 e^{-nu^2/2 + a nu} <= e^{a^2/2} * int_R |nu|^{2j}(nu^2+1)^2 e^{-(nu-a)^2/2}; bound |nu|^n (nu^2+1)^2 <= (|nu|+1)^{2j+4}
    # E(|a + Z| + 1)^{n} <= E(a + 1 + |Z|)^n <= 2^{n} ((a+1)^n + E|Z|^n), E|Z|^n <= (n-1)!! (n even)
    Lt = arb(0)
    for j, y in enumerate(ys):
        n = 2*j + 4
        df = arb(1)
        for i in range(n - 1, 0, -2):
            df *= i
        Lt += y*2**n*((a + 1)**n + df)*(2*arb.pi()).sqrt()
    Lt = 189*Lt*(a*a/2).exp()
    G2 = (tot + Lt)*(a/2).exp()/Lam
    _, p3 = FC.nf02()
    p3abs = sum((abs(arb(fmpq(c.numerator, c.denominator))) for c in p3), arb(0))   # sup_[0,1] |psi3| <= sum |coeffs|
    Gq = p3abs*(a*(arb(lam) + arb(17)/4).sqrt()).exp()
    return G2, Gq


def main():
    Dg = int(sys.argv[1]); Kmom = int(sys.argv[2]); J0 = int(sys.argv[3]); xmax = arb(sys.argv[4]); K = int(sys.argv[5])
    t0 = time.time()
    Z, Y = FP.load_moments(Dg, Kmom)
    if len(sys.argv) > 6:
        ctx.prec = int(sys.argv[6])
    lam = (Dg - 2)*(Dg - 1); L = arb(lam); Lam = L.sqrt()
    Q = [arb(q) for q in FP.Q_exact(Kmom, lam)]
    vmax = 1/(arb(J0) + arb(1)/2)**2
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr(f'Dg={Dg} J>={J0} (v <= {vmax.str(6)}), x <= {xmax.str(6)} (bhat <= {(2*xmax.sqrt()).str(5)}), K = {K}, Kmom = {Kmom}, prec {ctx.prec}')
    G2, Gq = crude_tail_bounds(Dg, lam, Lam)
    for b in (0, 4):
        assert J0 - b > K, 'need every k <= K to be a ladder term for all J >= J0'
        sig = 1 if b == 0 else -1
        # ---- tail bound TB (constant)
        TB = arb(0)
        for k in range(K + 1, Kmom + 1):
            Zt = abs(Z[k])/L**k; Yt = abs(Y[k])/L**(k + 1); Qt = abs(Q[k])/L**k
            Vk = 4*arb(k)/(k + 3) if b == 0 else 4 + arb(12)/(k + 1)
            TB += comb(k + 3, 3)**2*xmax**k*(Zt*(1/arb(k + 3) + xmax/(L*(k + 1))) + xmax*vmax/L*(2*L*Yt + Vk*Zt + L*Qt))
        # k > Kmom: |z_m|,|y_m| <= (4m^2+40m+100) 3^{m+1} * |f-bound|, m = k+2; |Q_k| <= 3^k Gq
        three = arb(3)
        k = Kmom + 1; crude = arb(0)
        while True:
            m = k + 2
            pz = (4*m*m + 40*m + 100)*three**(m + 1)*G2
            term = comb(k + 3, 3)**2*(xmax/L)**k*(pz*(1 + xmax/L) + xmax*vmax/L*(3*pz + L*three**k*Gq))
            crude += term
            if k > Kmom + 50 and term < crude*arb('1e-30'):
                # geometric remainder: ratio <= 3 xmax/lam * (k+4)^6/(k+3)^6 * 1.01 < 1/2
                crude += 2*term
                break
            k += 1
        TB += crude
        # ---- N(x, v) as a list over v-powers of arb_poly in x; forward recurrence F_{k+1} = F_k fac_k + t_{k+1} A_{k+1}, N = F_K
        def vmul_fac(F, j):                       # * [1 + (x/lam)(1 - c v)],  c = (b+j+1/2)^2
            c = (arb(b + j) + arb(1)/2)**2
            p1 = arb_poly([1, 1/L]); px = arb_poly([0, 1/L])
            out = [F[i]*p1 for i in range(len(F))] + [arb_poly([0])]
            for i in range(len(F)):
                out[i + 1] = out[i + 1] - (F[i]*px)*c
            return out
        def vmul_afac(Av, j):                     # v-polynomial (list of arb) * (1 - c v)
            c = (arb(b + j) + arb(1)/2)**2
            out = list(Av) + [arb(0)]
            for i in range(len(Av)):
                out[i + 1] = out[i + 1] - Av[i]*c
            return out
        def tk(k, Av):
            """t_k A_k as v-list of x-polys:  D_k^2 x^k A_k(v) [E Zt/(k+3) - (1 - v/4) x Zt/(k+1) + x v (2 lam Yt + V Zt + sig lam Qt)]."""
            Zt = Z[k]/L**k; Yt = Y[k]/L**(k + 1); Qt = Q[k]/L**k
            Vk = 4*arb(k)/(k + 3) if b == 0 else 4 + arb(12)/(k + 1)
            D2 = arb(comb(k + 3, 3)**2)
            # bracket as v-list: v^0: (lam + x) Zt/(k+3) - x Zt/(k+1);  v^1: -x Zt/(4(k+3)) + x Zt/(4(k+1)) + x (2 lam Yt + V Zt + sig lam Qt)
            b0 = arb_poly([L*Zt/(k + 3), Zt/(k + 3) - Zt/(k + 1)])
            b1 = arb_poly([0, -Zt/(4*(k + 3)) + Zt/(4*(k + 1)) + 2*L*Yt + Vk*Zt + sig*L*Qt])
            xk = arb_poly([0]*k + [D2])
            out = [arb_poly([0])]*(len(Av) + 1)
            for i, a in enumerate(Av):
                out[i] = out[i] + b0*xk*a
                out[i + 1] = out[i + 1] + b1*xk*a
            return out
        def vadd(F, G):
            n = max(len(F), len(G))
            return [(F[i] if i < len(F) else arb_poly([0])) + (G[i] if i < len(G) else arb_poly([0])) for i in range(n)]
        Av = [arb(1)]
        F = tk(0, Av)
        for k in range(K):
            F = vmul_fac(F, k)
            Av = vmul_afac(Av, k)
            F = vadd(F, tk(k + 1, Av))
        N = Bi()
        for jv, px in enumerate(F):
            for ix in range(px.degree() + 1):
                N.d[(ix, jv)] = px[ix]
        # subtract TB * sup(E Pi_K) (constant, positive): E <= lam + xmax, Pi_K <= (1 + xmax/lam)^K
        Emax = L + xmax; PiMax = (1 + xmax/L)**K
        N = N + Bi({(0, 0): -TB*Emax*PiMax})
        dx = max(i for i, j in N.d); dv = max(j for i, j in N.d)
        # self-check against the direct ladder sum at a physical state inside the box: J = J0 + 7, x = xmax/3
        Act = FP.Action(Dg, Z, Y, FP.Q_exact(Kmom, lam))
        Jc = J0 + 7; vc = 1/(arb(Jc) + arb(1)/2)**2; xc = xmax/3
        m2 = L/(xc*vc); Cc = Jc*(Jc + 1); Tc = m2 + Cc
        if Jc - b <= Kmom:
            Hd = Act.H_ang(Jc, b, Tc); rc = arb(Cc)/Tc
            Ec = L + xc - xc*vc/4
            Pic = arb(1)
            for jj in range(K):
                Pic *= 1 + (xc/L)*(1 - (arb(b + jj) + arb(1)/2)**2*vc)
            Nc = sum((c*xc**i*vc**j for (i, j), c in N.d.items()), arb(0))
            pr(f' b={b}: self-check at J={Jc}, x={xc.str(5)}: N/(E Pi_K) = {(Nc/(Ec*Pic)).str(10)}  vs  H_ang/(1-r)^3 = {(Hd/(1-rc)**3).str(10)}')
        pr(f' b={b}: tail bound TB = {TB.str(4)}, N degree (x {dx}, v {dv}), built {time.time()-t0:.0f}s')
        # ---- tensor Bernstein on [0, xmax] x [0, vmax], bisection in x
        rows = []
        for i in range(dx + 1):
            cv = [N.d.get((i, j), arb(0)) for j in range(dv + 1)]
            rows.append(bern1_fixed(cv, arb(0), vmax))
        def cert(xa, xb, depth):
            # Bernstein in v for each x-power (rows), then in x on [xa, xb]
            mn = None; ok = True
            for j in range(dv + 1):
                cx = [rows[i][j] for i in range(dx + 1)]
                bx = bern1_fixed(cx, xa, xb)
                for c in bx:
                    lo = float(c.lower())
                    mn = lo if mn is None or lo < mn else mn
                    if not c > 0:
                        ok = False
            if ok:
                return True, mn, 1
            if depth >= 12:
                return False, mn, 1
            xm = (xa + xb)/2
            o1, m1, n1 = cert(xa, xm, depth + 1)
            if not o1:
                return False, m1, n1
            o2, m2, n2 = cert(xm, xb, depth + 1)
            return o2, min(m1, m2), n1 + n2
        ok, mn, nc = cert(arb(0), xmax, 0)
        pr(f' b={b}: {"PASS" if ok else "FAIL"}  cells={nc}  min Bernstein lower bound of N = {mn:.4e}  (N(0,0) = lam Zt_0/3 = {(L*Z[0]/3).str(8)})  {time.time()-t0:.0f}s')
    open(os.path.join(FC.CACHE, f'fp_rect_{Dg}_J{J0}_x{sys.argv[4]}_K{K}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()

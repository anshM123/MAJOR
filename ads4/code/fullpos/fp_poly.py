"""
fullpos/fp_poly.py (FULLPOS_CERT.md §3): every-mass certificate at fixed spin J for channel b in {0, 4}.
  P_{J,b}(w) = sum_{k<=n} D_k^2 w^k (A_k + B_k w) prod_{j<k}(C-u_j) prod_{j=k}^{n-1}(1 - u_j w/lam),   n = J - b,
  A_k = Zt_k/(k+3), B_k = 2 Yt_k + ((V_{b,k} - C/(k+1))/lam) Zt_k + sigma_b Qt_k,
  Yt = Y/lam^{k+1}, Zt = Z/lam^k, Qt = Q/lam^k;  w = lam/T in [0, w_max], w_max = lam/((J+Dg-2)(J+Dg-1)).
H_ang,b = (1-r)^3 * P * T^n / prod_{j<n}(T-u_j) (positive factor on the heavy range), so P > 0 <=> H_ang,b > 0.
Certificate: Bernstein coefficients of P on [0, w_max] (Arb), all > 0, with de Casteljau bisection if needed.  [certified]
given the moment balls of fp_moments.py.
"""
import sys, os, json, time
sys.dont_write_bytecode = True
from flint import arb, arb_poly, fmpq, fmpq_poly, ctx
import fp_common as FC
import fp_moments as FM


def load_moments(Dg, K):
    d = json.load(open(os.path.join(FC.CACHE, f'moments_{Dg}_K{K}.json')))
    ctx.prec = max(ctx.prec, d['prec'])
    Z = [arb(s) for s in d['Z']]; Y = [arb(s) for s in d['Y']]
    return Z, Y


_QC = {}


def Q_exact(K, lam):
    """Q_k = int_0^1 psi3(p) q_k(lam p^2) dp, exact (fmpq)."""
    key = (K, lam)
    if key in _QC:
        return _QC[key]
    q = FM.cdh1_fmpq(K)
    _, p3 = FC.nf02()
    mom = {}
    def Mi(i):
        if i not in mom:
            mom[i] = sum((fmpq(c.numerator, c.denominator)/(2*i + n + 1) for n, c in enumerate(p3)), fmpq(0))
        return mom[i]
    out = []
    for k in range(K + 1):
        v = fmpq(0)
        for i, c in enumerate(q[k].coeffs()):
            v += c*fmpq(lam)**i*Mi(i)
        out.append(v)
    _QC[key] = out
    return out


class Action:
    def __init__(self, Dg, Z, Y, Q):
        self.Dg = Dg; self.lam = (Dg - 2)*(Dg - 1); self.Z = Z; self.Y = Y; self.Q = [arb(v) for v in Q]
        self.Lam2 = arb(self.lam)

    def P_poly(self, J, b):
        L = FC.ladder(J, b, T=1)
        n, C, u, D, V, sig = L['n'], L['C'], L['u'], L['D'], L['V'], L['sig']
        lam = self.Lam2
        W = arb_poly([0, 1])
        tot = arb_poly([0])
        # suffix products prod_{j=k}^{n-1}(1 - u_j w/lam)
        suf = [None]*(n + 1); suf[n] = arb_poly([1])
        for k in range(n - 1, -1, -1):
            suf[k] = suf[k + 1]*arb_poly([1, -arb(u[k])/lam])
        pre = arb(1)
        for k in range(n + 1):
            Zt = self.Z[k]/lam**k; Yt = self.Y[k]/lam**(k + 1); Qt = self.Q[k]/lam**k
            Ak = Zt/(k + 3)
            Bk = 2*Yt + ((arb(fmpq(V[k].numerator, V[k].denominator)) - arb(C)/(k + 1))/lam)*Zt + sig*Qt
            term = arb_poly([0]*k + [Ak, Bk])*suf[k]*(D[k]**2*pre)
            tot = tot + term
            if k < n:
                pre = pre*(C - u[k])
        return tot

    def H_ang(self, J, b, T):
        """H_ang,b at (T, J) directly (arb): (1-r)^3 sum_k D_k^2 U_k [Z_k(1/(k+3) - r/(k+1)) + (2Y_k + V Z_k + sig Lam^2 Q_k)/T]."""
        L = FC.ladder(J, b, T=1)
        n, C, u, D, V, sig = L['n'], L['C'], L['u'], L['D'], L['V'], L['sig']
        T = arb(T); r = arb(C)/T
        U = arb(1); tot = arb(0)
        for k in range(n + 1):
            Vk = arb(fmpq(V[k].numerator, V[k].denominator))
            tot += D[k]**2*U*(self.Z[k]*(1/arb(k + 3) - r/(k + 1)) + (2*self.Y[k] + Vk*self.Z[k] + sig*self.Lam2*self.Q[k])/T)
            if k < n:
                U = U*(C - u[k])/(T - u[k])
        return (1 - r)**3*tot


def bernstein(p, a, b):
    """Bernstein coefficients of arb_poly p on [a, b]."""
    n = p.degree()
    # q(t) = p(a + (b-a) t)
    q = p(arb_poly([a, b - a]))
    c = [q[i] if i <= q.degree() else arb(0) for i in range(n + 1)]
    from math import comb
    return [sum((arb(comb(i, j))/comb(n, j)*c[j] for j in range(i + 1)), arb(0)) for i in range(n + 1)]


def certify_pos(p, a, b, depth=0, maxdepth=30):
    """True if p > 0 on [a, b] is proved (Bernstein + bisection); returns (ok, min lower bound of Bernstein coeffs, #cells)."""
    bc = bernstein(p, a, b)
    lo = min(float(x.lower()) for x in bc)
    if all(x > 0 for x in bc):
        return True, lo, 1
    if bc[0] < 0 or bc[-1] < 0 or depth >= maxdepth:          # endpoint value negative -> fails (or give up)
        return False, lo, 1
    m = (a + b)/2
    ok1, l1, n1 = certify_pos(p, a, m, depth + 1, maxdepth)
    if not ok1:
        return False, l1, n1
    ok2, l2, n2 = certify_pos(p, m, b, depth + 1, maxdepth)
    return ok2, min(l1, l2), n1 + n2 + 1


def main():
    Dg = int(sys.argv[1]); K = int(sys.argv[2]); Js = sys.argv[3]
    J0, J1 = (int(v) for v in Js.split('-'))
    Z, Y = load_moments(Dg, K)
    lam = (Dg - 2)*(Dg - 1)
    Q = Q_exact(K, lam)
    A = Action(Dg, Z, Y, Q)
    out = []
    for J in range(J0, J1 + 1):
        for b in (0, 4):
            if (b == 0 and J % 2 == 1) or (b == 4 and J < 4):
                continue                      # b=0: even J (++ and odd structure); b=4: J >= 4 (+-, parity-odd for odd J)
            if J - b > K:
                continue
            t0 = time.time()
            P = A.P_poly(J, b)
            wmax = arb(lam)/((J + Dg - 2)*(J + Dg - 1))
            ok, lo, nc = certify_pos(P, arb(0), wmax)
            # scale: P(0), P(wmax)
            p0 = P(arb(0)); p1 = P(wmax)
            relrad = max(float(P[i].rad()/abs(P[i].mid())) if P[i].mid() != 0 else 0 for i in range(P.degree() + 1))
            line = (f'J={J:3d} b={b} n={J-b:3d}  {"PASS" if ok else "FAIL"} cells={nc:3d}  P(0)={p0.str(6)}  P(wmax)={p1.str(6)}'
                    f'  min Bernstein lower={lo:.3e}  max coeff rel.rad={relrad:.1e}  {time.time()-t0:.1f}s')
            print(line, flush=True); out.append(line)
    open(os.path.join(FC.CACHE, f'fp_poly_{Dg}_{J0}-{J1}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()

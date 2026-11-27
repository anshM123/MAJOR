r"""
m1_remainder/m1b_genfun.py  [D for the identities; NE for the numbers]

Exact generating function, in spin, of the per-state light action of Omega^EE on the twist-2 double-twist family
(Delta = J + 2, r = 1, U_{b,k} = 1), Delta_g = 1000:

   hhat_b(J) = T_J sum_{k<=n} a_k + sum_{k<=n} b_k,  n = J - b,  T_J = J(J+1),
   a_k = D_k^2 Z_k (1/(k+3) - 1/(k+1)),   b_k = D_k^2 (2 Y_k + V_{b,k} Z_k + sigma_b Lam^2 Q_k)

(sec_fullpos.tex Eq. (fp-H) at r = 1).  Define H_b(z) := sum_{J>=b} hhat_b(J) z^J.  Then, exactly,

   H_b(z) = theta(theta+1)[ z^b A(z)/(1-z) ] + z^b B(z)/(1-z),     theta = z d/dz,
   A(z) = int g(nu) A_nu(z) dnu,  B(z) = int g(nu) B_nu(z) dnu + sigma_b Lam^2 G_Q(z),

where A_nu, B_nu are finite linear combinations of q(theta) Psi_nu(z) with polynomial q, and of finitely many f_m(nu),
     Psi_nu(z) = sum_{m>=1} m^2 f_m(nu) z^m = z P_{-1/2+i nu}((1+z)/(1-z)) / (1-z)^2        (CDH generating function, d = 0),
     G_Q(z)    = sum_k D_k^2 Q_k z^k = (1-z)^{-4} int_0^1 psi_3(p) P_{-1/2+i Lam p}((1+z)/(1-z)) dp   (d = 1, identity I3).
(Every coefficient polynomial carries the factor m^2 needed to write it through Psi: checked symbolically below.)

Consequence used in the report: for any spin profile W(J) = int rho^J dmu(rho) (completely monotone, mu >= 0),
   sum_{J >= b, J even} W(J) hhat_b(J) = (1/2) int [ H_b(rho) + H_b(-rho) ] dmu(rho)        [P given (fp-H) at light states]
so the e^{pi Lam/2} enters a Bose-symmetric (even-J) sum only through H_b on (-1, 0], the back-to-back point z = -1.

Numerics: per nu, P and P' from mpmath (Ferrers / type 3), Taylor series in z from the Legendre ODE, theta-operators on
truncated series; nu-quadrature: Gauss-Legendre panels of width 1 (12 nodes) on [0, 1060]; chi from m1_weight (mpmath above 940).
Usage: python3 m1b_genfun.py validate | python3 m1b_genfun.py scan  (writes out/m1b_*.json)
"""
import sys, os, json, time, math
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
sys.path.insert(0, HERE)
from m1_weight import Weight

DPS = 40
NSER = 16


# ------------------------------------------------------------------ symbolic coefficient polynomials
def coefficient_polys():
    """For a_k and b_k (b = 0, 4): dict s -> polynomial q_s(m) (sympy, in m) with c_s(m - s) = m^2 q_s(m)."""
    k, m = sp.symbols('k m')
    D2 = ((k + 1)*(k + 2)*(k + 3)/6)**2
    M = k + 2
    # z_M = 12 f_M - (M-1) f_{M-1} + (M+1) f_{M+1};   y_M as in fp_common
    zc = {0: 12, -1: -(M - 1), 1: (M + 1)}
    yc = {0: sp.Rational(1, 2)*(12 - M**2), -1: sp.Rational(1, 2)*(13 - 7*M), 1: sp.Rational(1, 2)*(13 + 7*M),
          -2: (M - 2)**2/4, 2: (M + 2)**2/4}
    out = {}
    # a_k = D2 z_{k+2} (1/(k+3) - 1/(k+1))
    ca = {}
    for d, c in zc.items():
        ca[2 + d] = sp.simplify(D2*c*(1/(k + 3) - 1/(k + 1)))
    out['a'] = ca
    for b in (0, 4):
        V = 4*k/(k + 3) if b == 0 else 4 + 12/(k + 1)
        cb = {}
        for d, c in yc.items():
            cb[2 + d] = cb.get(2 + d, 0) + D2*2*c
        for d, c in zc.items():
            cb[2 + d] = cb.get(2 + d, 0) + D2*V*c
        out[f'b{b}'] = {s: sp.simplify(v) for s, v in cb.items()}
    res = {}
    for key, cs in out.items():
        res[key] = {}
        for s, c in cs.items():
            cm = sp.factor(sp.together(c.subs(k, m - s)))
            q, r = sp.div(sp.Poly(sp.numer(cm), m), sp.Poly(sp.denom(cm)*m**2, m))
            assert r.is_zero, (key, s, cm)
            res[key][s] = [sp.Rational(x) for x in reversed(q.all_coeffs())]   # q_s(m) = sum_i coeff_i m^i
    return res


# ------------------------------------------------------------------ truncated power series in h (z = z0 + h)
def s_mul(a, b, N):
    c = [mp.mpf(0)]*N
    for i, ai in enumerate(a[:N]):
        if ai == 0:
            continue
        for j in range(min(len(b), N - i)):
            c[i + j] += ai*b[j]
    return c


def s_theta(a, z0):
    """theta = (z0 + h) d/dh on a series a (loses the top order)."""
    N = len(a)
    d = [(i + 1)*a[i + 1] for i in range(N - 1)] + [mp.mpf(0)]
    return [z0*d[i] + (d[i - 1] if i >= 1 else 0) for i in range(N)][:N - 1]


def s_poly_theta(q, a, z0):
    """sum_i q_i theta^i a."""
    out = [mp.mpf(0)]*len(a)
    cur = list(a)
    for i, qi in enumerate(q):
        if i > 0:
            cur = s_theta(cur, z0)
        for j in range(len(cur)):
            out[j] += mp.mpf(qi.p)/qi.q*cur[j]
    L = len(a) - (len(q) - 1)
    return out[:L]


def s_pow_z(z0, e, N):
    """series of z^e around z0 (e integer, may be negative)."""
    return [mp.binomial(e, i)*z0**(e - i) for i in range(N)]


def s_inv_1mz(z0, N):
    """1/(1 - z) = sum_i h^i/(1-z0)^{i+1}."""
    return [1/(1 - z0)**(i + 1) for i in range(N)]


def legendre_taylor(lam, X0, N):
    """P^{(j)}(X0)/j!, j < N, from P, P' (mpmath) and the Legendre ODE."""
    typ = 2 if abs(X0) < 1 else 3
    P = mp.legenp(lam, 0, X0, type=typ)
    Pm = mp.legenp(lam - 1, 0, X0, type=typ)
    dP = lam*(X0*P - Pm)/(X0*X0 - 1)
    if typ == 2:
        P, dP = mp.re(P), mp.re(dP)
    else:
        P, dP = mp.re(P), mp.re(dP)
    L = lam*(lam + 1)
    der = [P, dP]
    for j in range(N - 2):
        der.append((2*(j + 1)*X0*der[j + 1] + (j*(j + 1) - L)*der[j])/(1 - X0*X0))
    der = [mp.re(d) for d in der]
    return [der[j]/mp.factorial(j) for j in range(N)]


def compose_X(z0, N):
    """delta(h) = X(z0+h) - X0 = sum_{j>=1} 2 h^j/(1-z0)^{j+1}."""
    return [mp.mpf(0)] + [2/(1 - z0)**(j + 1) for j in range(1, N)]


def series_P_of_z(lam, z0, N):
    X0 = (1 + z0)/(1 - z0)
    tay = legendre_taylor(lam, X0, N)
    dl = compose_X(z0, N)
    out = [mp.mpf(0)]*N
    pw = [mp.mpf(1)] + [mp.mpf(0)]*(N - 1)
    for j in range(N):
        for i in range(N):
            out[i] += tay[j]*pw[i]
        pw = s_mul(pw, dl, N)
    return out


def f_values(nu, M):
    s = mp.mpf(nu)**2 + mp.mpf(1)/4
    f = [mp.mpf(0), mp.mpf(1)]
    for mm in range(1, M):
        f.append((-(s - 2*mm*mm)*f[mm] - (mm - 1)**2*f[mm - 1])/(mm + 1)**2)
    return f


class GenFun:
    def __init__(self):
        mp.mp.dps = DPS
        self.Q = coefficient_polys()
        self.W = Weight(1000)
        self.Lam = mp.sqrt(mp.mpf(998*999))
        # nu nodes: GL panels width 1 on [0, 1060]
        xg, wg = np.polynomial.legendre.leggauss(12)
        nodes, wts = [], []
        for p in range(1060):
            for x, w in zip(xg, wg):
                nodes.append(p + 0.5 + 0.5*x); wts.append(0.5*w)
        self.nu = nodes; self.wnu = wts
        chi = self.W.chi_np(np.array(nodes))
        th = self.W.theta_np(np.array(nodes))
        self.g = []
        for i, nu in enumerate(nodes):
            if nu > 940:
                c = self.W.chi_mp(nu); t = self.W.theta_mp(nu)
            else:
                c = mp.mpf(float(chi[i])); t = mp.mpf(float(th[i]))
            x = mp.mpf(nu)**2
            self.g.append(-9216*(4*x + 1)/((4*x + 9)*(4*x + 49))*c*t/self.Lam)
        # psi_3
        d = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'nf02_psi.json')))
        K0 = mp.mpf('27.48509817')
        self.psi3 = [mp.mpf(Fr(c).numerator)/Fr(c).denominator/K0 for c in d['psi3_exact']]

    def per_nu(self, nu, z0, N=NSER):
        """A_nu(z), B_nu(z) for b = 0, 4 as truncated series around z0 (length >= 3)."""
        lam = mp.mpf(-0.5) + 1j*mp.mpf(nu)
        Py = series_P_of_z(lam, z0, N)
        # Psi = z P/(1-z)^2
        zs = [z0, mp.mpf(1)] + [mp.mpf(0)]*(N - 2)
        inv = s_inv_1mz(z0, N)
        Psi = s_mul(s_mul(zs, Py, N), s_mul(inv, inv, N), N)
        f = f_values(nu, 8)
        res = {}
        for key, qs in self.Q.items():
            tot = None
            for s, q in qs.items():
                ser = s_poly_theta(q, Psi, z0)
                # subtract excluded low terms m = 1..s-1:  m^2 q(m) f_m z^m
                L = len(ser)
                for mm in range(1, s):
                    cm = sum(mp.mpf(qi.p)/qi.q*mm**i for i, qi in enumerate(q))*mm*mm*f[mm]
                    zm = s_pow_z(z0, mm, L)
                    ser = [ser[i] - cm*zm[i] for i in range(L)]
                ser = s_mul(s_pow_z(z0, -s, L), ser, L)
                tot = ser if tot is None else [tot[i] + ser[i] for i in range(min(len(tot), L))]
            res[key] = tot
        return res

    def GQ(self, z0, N=6, npts=3000):
        X0 = (1 + z0)/(1 - z0)
        xg, wg = np.polynomial.legendre.leggauss(npts)
        acc = [mp.mpf(0)]*N
        dl = compose_X(z0, N)
        for x, w in zip(xg, wg):
            p = mp.mpf(0.5 + 0.5*x)
            ps = sum(c*p**i for i, c in enumerate(self.psi3))
            lam = mp.mpf(-0.5) + 1j*self.Lam*p
            tay = legendre_taylor(lam, X0, N)
            ser = [mp.mpf(0)]*N; pw = [mp.mpf(1)] + [mp.mpf(0)]*(N - 1)
            for j in range(N):
                for i in range(N):
                    ser[i] += tay[j]*pw[i]
                pw = s_mul(pw, dl, N)
            for i in range(N):
                acc[i] += mp.mpf(0.5*w)*ps*ser[i]
        inv = s_inv_1mz(z0, N)
        inv4 = s_mul(s_mul(inv, inv, N), s_mul(inv, inv, N), N)
        return s_mul(acc, inv4, N)

    def H(self, z0, verbose=False):
        """H_0(z0), H_4(z0), and F(z0) = sum_k P_k z^k (check), all as mp numbers."""
        z0 = mp.mpf(z0)
        acc = {}
        t0 = time.time()
        for i, nu in enumerate(self.nu):
            r = self.per_nu(nu, z0)
            gw = self.g[i]*self.wnu[i]
            for key, ser in r.items():
                if key not in acc:
                    acc[key] = [mp.mpf(0)]*len(ser)
                for j in range(min(len(ser), len(acc[key]))):
                    acc[key][j] += gw*ser[j]
            if verbose and i % 2000 == 0:
                print(f'   z0={float(z0):+.4f} nu-node {i}/{len(self.nu)} {time.time()-t0:.0f}s', flush=True)
        GQ = self.GQ(z0)
        lam2 = mp.mpf(998*999)
        out = {}
        Lr = 4
        inv = s_inv_1mz(z0, Lr)
        for b, sig in ((0, 1), (4, -1)):
            A = acc['a'][:Lr]; B = [acc[f'b{b}'][i] + sig*lam2*GQ[i] for i in range(Lr)]
            zb = s_pow_z(z0, b, Lr)
            t1 = s_mul(s_mul(zb, A, Lr), inv, Lr)
            th1 = s_theta(t1, z0)                       # theta t1
            th2 = s_theta(th1, z0)                      # theta^2 t1
            val = th2[0] + th1[0] + s_mul(s_mul(zb, B, Lr), inv, Lr)[0]
            out[f'H{b}'] = val
        # F(z) = sum_k P_k z^k where P_k = D_k^2 L_k: equals A(z) + (1-z) M(z) ... (check only the leading-Regge identity)
        out['A'] = acc['a'][0]
        return out


def moments_hhat(K=1000):
    """hhat_b(J), twist 2, J <= K, from the certified moments (arb -> mpf)."""
    sys.path.insert(0, os.path.join(HERE, '..', 'fullpos'))
    import fp_poly as FP
    from flint import arb, ctx
    Z, Y = FP.load_moments(1000, K)
    ctx.prec = 3000
    lam = 998*999
    Q = [arb(q) for q in FP.Q_exact(K, lam)]
    from math import comb
    D = [comb(k + 3, 3) for k in range(K + 1)]
    res = {}
    for b, sig in ((0, 1), (4, -1)):
        al = arb(0); be = arb(0); vals = {}
        for k in range(K + 1):
            V = arb(4*k)/(k + 3) if b == 0 else 4 + arb(12)/(k + 1)
            al += D[k]**2*Z[k]*(1/arb(k + 3) - 1/arb(k + 1))
            be += D[k]**2*(2*Y[k] + V*Z[k] + sig*lam*Q[k])
            J = k + b
            if J <= K:
                vals[J] = (al*J*(J + 1) + be)
        res[b] = vals
    return res


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'validate'
    G = GenFun()
    if mode == 'validate':
        mh = moments_hhat()
        from flint import arb
        out = []
        for z0 in (sys.argv[2:] or ['-0.3', '0.3']):
            r = G.H(mp.mpf(z0), verbose=True)
            for b in (0, 4):
                s = arb(0)
                for J, v in mh[b].items():
                    s += v*arb(z0)**J
                sm = mp.mpf(s.mid().str(40, radius=False))
                line = (f'z0={z0}: H_{b} closed form = {mp.nstr(r[f"H{b}"], 15)}   moments sum_(J<=1000) = {mp.nstr(sm, 15)}'
                        f'   rel.diff = {mp.nstr(abs(r[f"H{b}"]/sm - 1), 3)}')
                print(line, flush=True); out.append(line)
        open(os.path.join(OUT, f'm1b_validate_{"_".join(sys.argv[2:]) or "default"}.txt'), 'w').write('\n'.join(out) + '\n')
    elif mode == 'point':
        res = {}
        for z0 in sys.argv[2:]:
            r = G.H(mp.mpf(z0))
            res[z0] = {k: mp.nstr(v, 25) for k, v in r.items()}
            print(z0, res[z0], flush=True)
        tag = '_'.join(sys.argv[2:])
        json.dump(res, open(os.path.join(OUT, f'm1b_H_{tag}.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()

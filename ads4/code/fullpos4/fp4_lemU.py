"""
fullpos4/fp4_lemU.py (FULLPOS4.md §6) [DV, midpoint arb arithmetic]: the author's Lemma U (ANGULAR_ROUTE §80) at fixed c.
At fixed c the whole twist interval a >= 998 is one polynomial in y = x/x_*, x = h/d^2 (x_* at a = 998):
   H_M = P(c, y)/D^4,  P = sum_m A_m y^m + d^2 x_* y sum_m (q_b A_m + B_m) y^m,  D = 1 + q_b d^2 x_* y,
   A_m = sum_j R_{m,j}(c, x_*) Y^A_j,  B_m = same with Y^B (G_B = D_b G_A + sigma_b G_D),  Y_j = d^3 (c d)^j G^{(j)}(c)/j!.
Integrated G (IDEAL NF02 weights): G_A = 2 int psi2(p) F_{3,Lam p} dp,  G_D = Lam^2 int psi3(p) F_{4,Lam p} dp (Gauss-Legendre in p,
per-frequency jets from the ODE recurrence of fp4_ref).  D_b in t (z = c + c d t):
   [t^j] D_b g = [(j+2)(j+1) g_{j+2} + ((d-c) j + (7-4b) c + 1)(j+1) g_{j+1} + (-c d j(j-1) + (7-4b) c d j) g_j]/(2cd) + (9/2 - b(b-1)) g_j.
Chebyshev on y in [0,1]: Phat_l = (2 - delta_l0) sum_{m>=l} C(2m, m-l) 4^{-m} p_m;  test Phat_0 - sum_{l>=1}|Phat_l| > 0.
Usage: fp4_lemU.py L M Np prec   (both channels; jets computed once per c)
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from flint import arb, acb, ctx
import fp4_ref as R, fp4_int as I

def mid(x):
    return x.mid() if isinstance(x, arb) else acb(x.real.mid(), x.imag.mid())

def jet_mid(p, nu, c, d, w, K):
    y0, y1 = R.seeds(p, nu, c, d, w)
    nu = acb(nu)
    A = [arb(1), d - 2*c, c*c - 2*c*d, c*c*d]
    q = (1 + 2*p)*c
    Bc = [1 - q, -q*d - c*(1 - q), c*q*d]
    Cc = [c*(nu*nu + acb(0.25) - p + p*p*c), acb(c*p*p*c*d)]
    y = [mid(y0), mid(y1)]
    for k in range(0, K - 1):
        acc = acb(0)
        for i in range(1, 4):
            if k - i + 2 >= 0:
                acc += A[i]*(k - i + 2)*(k - i + 1)*y[k - i + 2]
        for i in range(0, 3):
            if k - i + 1 >= 0:
                acc += Bc[i]*(k - i + 1)*y[k - i + 1]
        for i in range(0, 2):
            if k - i >= 0:
                acc += Cc[i]*y[k - i]
        y.append(mid(-acc/((k + 2)*(k + 1))))
    return [v.real for v in y[:K + 1]]

def xstar(c, s, a=998):
    c = arb(c); d = 1 - c; a = arb(a)
    return 1/(a*a*(1 + c + (4*c + s*s*d*d/(a*a)).sqrt()))

def run(L, M, Np, prec):
    ctx.prec = prec
    t0 = time.time()
    c = 1 - arb(-L).exp() if L < 30 else None
    c = arb(c.mid()); d = 1 - c
    K = 3*M + 3
    w2 = I.Weight(2); w3 = I.Weight(3)
    xg, wg = np.polynomial.legendre.leggauss(Np)
    YA = [arb(0)]*(K + 1); YD = [arb(0)]*(K + 1)
    LAM = arb(998*999).sqrt()
    for xn, wn in zip(xg, wg):
        p = arb((xn + 1)/2); wp = arb(wn/2)
        nu = LAM*p
        psi2 = sum((arb(cf.numerator)/cf.denominator*p**i for i, cf in enumerate(w2.c)), arb(0))
        psi3 = sum((arb(cf.numerator)/cf.denominator*p**i for i, cf in enumerate(w3.c)), arb(0))
        y3 = jet_mid(3, nu, c, d, d, K)            # d^3 [t^j] F_3
        y4 = jet_mid(4, nu, c, d, d, K)            # d^4 [t^j] F_4
        fa = 2*psi2*wp; fd = LAM*LAM*psi3*wp/d
        for j in range(K + 1):
            YA[j] = (YA[j] + fa*y3[j]).mid(); YD[j] = (YD[j] + fd*y4[j]).mid()
    tj = time.time() - t0
    out = []
    for b in (0, 4):
        s = 2*b + 1; q_b = b*(b + 1); sig = 1 if b == 0 else -1; kb = arb(9)/2 - b*(b - 1); e = 7 - 4*b
        YB = []
        for j in range(K - 1):
            v = ((j + 2)*(j + 1)*YA[j + 2] + ((d - c)*j + e*c + 1)*(j + 1)*YA[j + 1] + (-c*d*j*(j - 1) + e*c*d*j)*YA[j])/(2*c*d) + kb*YA[j]
            YB.append((v + sig*YD[j]).mid())
        xs = xstar(c, s)
        # R recursion at x = x_*, accumulating A_m, B_m
        Am = [YA[0]]; Bm = [YB[0]]
        prev = [arb(1)] + [arb(0)]*(3*M + 3)
        bigR = 0.0; bigA = 0.0
        for m in range(1, M + 1):
            top = 3*(m - 1) + 1
            Bv = [(prev[j + 1] - prev[j]/c).mid() for j in range(top + 1)] + [arb(0)]*3
            cur = [arb(0)]*(3*M + 4)
            for j in range(0, 3*m):
                acc = d*d*(j*(j + s))*Bv[j] if j <= top else arb(0)
                if j >= 1 and j - 1 <= top:
                    acc += d*(j*(2*j + s - 1))*Bv[j - 1]
                if j >= 2 and j - 2 <= top:
                    acc += (j*(j - 1))*Bv[j - 2]
                cur[j + 1] = (xs*acc).mid()
            A_ = sum((cur[j]*YA[j] for j in range(1, 3*m + 1)), arb(0)).mid()
            B_ = sum((cur[j]*YB[j] for j in range(1, 3*m + 1)), arb(0)).mid()
            Am.append(A_); Bm.append(B_)
            if m % 25 == 0 or m == M:
                bigR = max(bigR, max(float(abs(v).log())/math.log(10) for v in cur[1:3*m + 1] if v != 0))
            prev = cur
        bigA = max(float(abs(v).log())/math.log(10) for v in Am if v != 0)
        dx = d*d*xs
        pm = [Am[0]] + [(Am[m] + dx*(q_b*Am[m - 1] + Bm[m - 1])).mid() for m in range(1, M + 1)] + [(dx*(q_b*Am[M] + Bm[M])).mid()]
        # Chebyshev conversion on y in [0,1]
        D = len(pm) - 1
        four = arb(4)
        Ph = []
        for l in range(D + 1):
            sacc = arb(0)
            for m in range(l, D + 1):
                sacc += arb.bin_uiui(2*m, m - l)/four**m*pm[m]
            Ph.append(((2 if l else 1)*sacc).mid())
        P0 = Ph[0]; rest = sum((abs(v) for v in Ph[1:]), arb(0))
        P1 = sum(pm, arb(0)); Dmax = 1 + q_b*dx
        dec = {l: float(abs(Ph[l]).mid()) for l in [1, 2, 3, 5, 10, 20, 50, 100, 200, 400, 700, 1000] if l < len(Ph)}
        r = dict(L=L, c=float(c.mid()), b=b, M=M, Np=Np, prec=prec, xstar=float(xs.mid()), A0_dGc=float(Am[0].mid()),
                 P_at_y1=float(P1.mid()), H_at_a998=float((P1/Dmax**4).mid()), Phat0=float(P0.mid()), sum_abs_Phat_ge1=float(rest.mid()),
                 margin_rel=float(((P0 - rest)/abs(P0)).mid()), Phat_decay=dec, log10_max_A=bigA, log10_max_R=bigR, t_jets=tj, secs=time.time() - t0)
        print(json.dumps(r), flush=True)
        open(os.path.join(HERE, 'cache', 'fp4_lemU.jsonl'), 'a').write(json.dumps(r) + '\n')
        out.append(r)
    return out

if __name__ == '__main__':
    L = float(sys.argv[1]); M = int(sys.argv[2]); Np = int(sys.argv[3]); prec = int(sys.argv[4])
    run(L, M, Np, prec)

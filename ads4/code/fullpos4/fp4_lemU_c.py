"""
fullpos4/fp4_lemU_c.py (FULLPOS4.md §6.2) [DV]: Lemma U twist polynomial at COMPLEX impact parameter u = beta_c + i eta (Lemma C test).
c = tanh^2(u/2) complex, d = 1 - c, x_* = 1/(a^2 (1 + c + sqrt(4c + s^2 d^2/a^2))) at a = 998 (analytic continuation), seeds
P_nu(cosh u), d_beta P_nu(cosh u); same jets, D_b contraction, R recursion and P(y) = sum A_m y^m + d^2 x_* y sum (q_b A_m + B_m) y^m,
H_M(u, y) = P/D^4.  Reported: max_{y in [0,1]} |H_M(u, y) - H_M(u, 0)| / F_cert(Re u), |H_M(u, 0)|/F_cert, and the real-u reference.
Usage: fp4_lemU_c.py u0 eta1,eta2,... M Np prec
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'fullpos3'))
import numpy as np
from flint import arb, acb, ctx
import fp4_int as I

def cm(x):
    return acb(x.real.mid(), x.imag.mid())

def seeds(p, nu, c, d, w):
    X = (1 + c)/d
    mu = acb(-0.5) + acb(0, 1)*nu
    P0 = X.legendre_p(mu, 0, 3); Pm = X.legendre_p(mu - 1, 0, 3)
    dP = mu*(X*P0 - Pm)/(X*X - 1)
    sh = 2*c.sqrt()/d
    f = (w/d)**p
    return f*P0, f*(p*c*P0 + c.sqrt()*dP*sh)

def jet(p, nu, c, d, K):
    y0, y1 = seeds(p, acb(nu), c, d, d)
    A = [acb(1), d - 2*c, c*c - 2*c*d, c*c*d]
    q = (1 + 2*p)*c
    Bc = [1 - q, -q*d - c*(1 - q), c*q*d]
    Cc = [c*(acb(nu)*nu + acb(0.25) - p + p*p*c), c*p*p*c*d]
    y = [cm(y0), cm(y1)]
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
        y.append(cm(-acc/((k + 2)*(k + 1))))
    return y[:K + 1]

def run(u, M, Np, Fc):
    t0 = time.time()
    c = (u/2).tanh()**2; d = 1 - c
    K = 3*M + 3
    w2 = I.Weight(2); w3 = I.Weight(3)
    xg, wg = np.polynomial.legendre.leggauss(Np)
    YA = [acb(0)]*(K + 1); YD = [acb(0)]*(K + 1)
    LAM = arb(998*999).sqrt()
    for xn, wn in zip(xg, wg):
        p = arb((xn + 1)/2); wp = arb(wn/2)
        nu = LAM*p
        psi2 = sum((arb(cf.numerator)/cf.denominator*p**i for i, cf in enumerate(w2.c)), arb(0))
        psi3 = sum((arb(cf.numerator)/cf.denominator*p**i for i, cf in enumerate(w3.c)), arb(0))
        y3 = jet(3, nu, c, d, K); y4 = jet(4, nu, c, d, K)
        fa = 2*psi2*wp; fd = LAM*LAM*psi3*wp/d
        for j in range(K + 1):
            YA[j] = cm(YA[j] + fa*y3[j]); YD[j] = cm(YD[j] + fd*y4[j])
    res = []
    for b in (0, 4):
        s = 2*b + 1; q_b = b*(b + 1); sig = 1 if b == 0 else -1; kb = arb(9)/2 - b*(b - 1); e = 7 - 4*b
        YB = [cm(((j + 2)*(j + 1)*YA[j + 2] + ((d - c)*j + e*c + 1)*(j + 1)*YA[j + 1] + (-c*d*j*(j - 1) + e*c*d*j)*YA[j])/(2*c*d)
                 + kb*YA[j] + sig*YD[j]) for j in range(K - 1)]
        a = arb(998)
        xs = 1/(a*a*(1 + c + (4*c + s*s*d*d/(a*a)).sqrt()))
        Am = [YA[0]]; Bm = [YB[0]]
        prev = [acb(1)] + [acb(0)]*(3*M + 3)
        for m in range(1, M + 1):
            top = 3*(m - 1) + 1
            Bv = [cm(prev[j + 1] - prev[j]/c) for j in range(top + 1)] + [acb(0)]*3
            cur = [acb(0)]*(3*M + 4)
            for j in range(0, 3*m):
                acc = d*d*(j*(j + s))*Bv[j] if j <= top else acb(0)
                if j >= 1 and j - 1 <= top:
                    acc += d*(j*(2*j + s - 1))*Bv[j - 1]
                if j >= 2 and j - 2 <= top:
                    acc += (j*(j - 1))*Bv[j - 2]
                cur[j + 1] = cm(xs*acc)
            Am.append(cm(sum((cur[j]*YA[j] for j in range(1, 3*m + 1)), acb(0))))
            Bm.append(cm(sum((cur[j]*YB[j] for j in range(1, 3*m + 1)), acb(0))))
            prev = cur
        dx = d*d*xs
        pm = [Am[0]] + [cm(Am[m] + dx*(q_b*Am[m - 1] + Bm[m - 1])) for m in range(1, M + 1)] + [cm(dx*(q_b*Am[M] + Bm[M]))]
        H0 = pm[0]
        worst = 0.0; wy = 0.0; mabs = 0.0; my_ = 0.0
        for yv in np.linspace(0, 1, 401):
            yb = arb(float(yv)); v = acb(0)
            for coef in reversed(pm):
                v = v*yb + coef
            Hy = v/(1 + q_b*dx*yb)**4
            dlt = float(abs(Hy - H0).mid())
            ah = float(abs(Hy).mid())
            if ah > mabs:
                mabs, my_ = ah, float(yv)
            if dlt > worst:
                worst, wy = dlt, float(yv)
        res.append(dict(u=[float(u.real.mid()), float(u.imag.mid())], b=b, M=M, Np=Np, prec=ctx.prec, Fcert=Fc,
                        H0=[float(H0.real.mid()), float(H0.imag.mid())], absH0_over_F=float(abs(H0).mid())/Fc,
                        maxdiff_over_F=worst/Fc, at_y=wy, maxabs_over_F=mabs/Fc, maxabs_at_y=my_, secs=time.time() - t0))
    return res

if __name__ == '__main__':
    u0 = float(sys.argv[1]); etas = [float(x) for x in sys.argv[2].split(',')]
    M = int(sys.argv[3]); Np = int(sys.argv[4]); ctx.prec = int(sys.argv[5])
    import fp3_htdata as HD
    Fc = None
    for ba, bb, F in HD.fcert_cells():
        if float(ba) <= u0 <= float(bb):
            Fc = float(F); break
    for eta in etas:
        for r in run(acb(u0, eta), M, Np, Fc):
            print(json.dumps(r), flush=True)
            open(os.path.join(HERE, 'cache', 'fp4_lemU_c.jsonl'), 'a').write(json.dumps(r) + '\n')

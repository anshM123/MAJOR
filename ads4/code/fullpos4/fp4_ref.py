"""
fullpos4/fp4_ref.py (FULLPOS4.md §5) [DV + ball arithmetic]: the author's §78 quadrature-free reference.
  K_ref = Q_M(theta) F_{p,nu}(c) = sum_{j<=3M} W_j [t^j] Y(t),   Y(t) = w^p F_{p,nu}(c + c d t),
  c = n(n+s)/(N(N+s)), d = 1 - c, h = 1/(N(N+s)), n = J - b, N = n + a, w = 1 - J(J+1)/((Delta-1)(Delta-2)).
  W_j = sum_{m<=M} R_{m,j};  R_{0,0} = 1;  B_j = R_{m-1,j+1} - R_{m-1,j}/c;
  R_{m,j+1} = (h/d^2)[ j(j-1) B_{j-2} + d j(2j+s-1) B_{j-1} + d^2 j(j+s) B_j ],  R_{m,0} = 0 (m > 0).
  Jet: (1+dt)(1-ct)^2 Y'' + (1-ct)[1-(2p+1)c(1+dt)] Y' + c[p^2 c(1+dt) + nu^2 + 1/4 - p] Y = 0,
       Y(0) = (w/d)^p P_nu(cosh beta_c), Y'(0) = (w/d)^p [p c P_nu + sqrt(c) d_beta P_nu].
Exact comparison: w^p sum_{k<=n} U_k g_k.   Inputs may be arb balls (cells): then the output ball is a NAIVE interval enclosure.
Usage: fp4_ref.py J a b nu1,nu2,... [M] [prec]
"""
import sys, os, math, time, json
sys.dont_write_bytecode = True
from flint import arb, acb, ctx


def params(J, a, b):
    J = arb(J); a = arb(a)
    s = 2*b + 1; n = J - b; N = n + a
    c = n*(n + s)/(N*(N + s)); h = 1/(N*(N + s)); d = 1 - c
    Delta = N + b + 2; w = 1 - J*(J + 1)/((Delta - 1)*(Delta - 2))
    return dict(s=s, n=n, N=N, c=c, h=h, d=d, w=w)


def Wcoef(c, h, d, s, M):
    """W_j = sum_{m<=M} R_{m,j}, j <= 3M (arb); also max |R_{m,j}| (log10)."""
    W = [arb(0)]*(3*M + 2)
    W[0] = arb(1)
    prev = [arb(1)] + [arb(0)]*(3*M + 2)
    hd = h/(d*d); big = 0.0
    for m in range(1, M + 1):
        top = 3*(m - 1) + 1
        Bv = [prev[j + 1] - prev[j]/c for j in range(top + 1)] + [arb(0)]*3
        cur = [arb(0)]*(3*M + 3)
        for j in range(0, 3*m):
            acc = d*d*(j*(j + s))*Bv[j] if j <= top else arb(0)
            if j >= 1 and j - 1 <= top:
                acc += d*(j*(2*j + s - 1))*Bv[j - 1]
            if j >= 2 and j - 2 <= top:
                acc += (j*(j - 1))*Bv[j - 2]
            cur[j + 1] = hd*acc
        for j in range(1, 3*m + 1):
            W[j] += cur[j]
            if m % 50 == 0 or m == M:
                x = abs(cur[j].mid())
                if x > 0:
                    big = max(big, float(x.log())/math.log(10))
        prev = cur
    return W, big


def seeds(p, nu, c, d, w):
    nu = acb(nu)
    X = acb((1 + c)/d)
    mu = acb(-0.5) + acb(0, 1)*nu
    P0 = X.legendre_p(mu, 0, 3); Pm = X.legendre_p(mu - 1, 0, 3)
    dP = mu*(X*P0 - Pm)/(X*X - 1)            # dP/dX
    sh = 2*c.sqrt()/d                        # sinh beta_c
    dbP = dP*sh
    f = (w/d)**p
    return f*P0, f*(p*c*P0 + c.sqrt()*dbP)


def jet(p, nu, c, d, w, K):
    y0, y1 = seeds(p, nu, c, d, w)
    nu = acb(nu)
    A = [arb(1), d - 2*c, c*c - 2*c*d, c*c*d]
    q = (1 + 2*p)*c
    Bc = [1 - q, -q*d - c*(1 - q), c*q*d]
    Cc = [c*(nu*nu + acb(0.25) - p + p*p*c), acb(c*p*p*c*d)]
    y = [y0, y1]
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
        y.append(-acc/((k + 2)*(k + 1)))
    return y[:K + 1]


def kref(p, nu, P, W):
    y = jet(p, nu, P['c'], P['d'], P['w'], len(W) - 1)
    return sum((W[j]*y[j] for j in range(len(W))), acb(0)).real


def kexact(p, nu, J, a, b):
    """w^p sum_{k<=n} U_k g_k (integer J, real a)."""
    s = 2*b + 1; n = J - b; N = arb(n) + arb(a)
    Delta = N + b + 2; w = 1 - arb(J*(J + 1))/((Delta - 1)*(Delta - 2))
    sn = arb(nu)**2 + arb(0.25)
    g0, gm = arb(1), arb(0); U = arb(1); tot = arb(0)
    for k in range(n + 1):
        tot += U*g0
        g1 = ((2*k*k + 2*p*k + p - sn)*g0 - (k - 1 + p)**2*gm)/((k + 1)**2)
        gm, g0 = g0, g1
        U = U*(k - n)*(n + s + k)/((k - N)*(N + s + k))
    return w**p*tot


if __name__ == '__main__':
    J = int(sys.argv[1]); a = float(sys.argv[2]); b = int(sys.argv[3])
    nus = [float(x) for x in sys.argv[4].split(',')]
    M = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
    ctx.prec = int(sys.argv[6]) if len(sys.argv) > 6 else 1600
    p = 3
    t0 = time.time()
    P = params(J, a, b)
    W, big = Wcoef(P['c'], P['h'], P['d'], P['s'], M)
    tW = time.time() - t0
    out = []
    for nu in nus:
        t1 = time.time()
        r = kref(p, nu, P, W)
        old = ctx.prec; ctx.prec = max(old, int(1.45*math.pi*nu/2*1.443) + 400)
        e = kexact(p, nu, J, a, b)
        ctx.prec = old
        line = dict(J=J, a=a, b=b, M=M, prec=ctx.prec, nu=nu, kref=r.str(12, radius=True), kexact=e.str(12, radius=True),
                    abs_diff=float(abs((r - e).mid())), log10_maxR=big, tW=tW, t=time.time() - t1)
        print(json.dumps(line), flush=True)
        open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'fp4_ref.jsonl'), 'a').write(json.dumps(line) + '\n')

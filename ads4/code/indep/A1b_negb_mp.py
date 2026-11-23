"""High-precision check of the negative-b direct series and of the Euler identity Ghat_{-b} = v^b Ghat_b
(the double-precision series for b' < 0 loses digits; this re-runs the same recursion in 50-digit arithmetic)."""
import mpmath as mp, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cblocks
mp.mp.dps = 50
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'A1b_negb_mp.txt')
lines = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

def recursion_mp(nu, b, Nmax):
    l = mp.mpc(-0.5, nu); Delta = 1
    C = Delta * (Delta - 3) + l * (l + 1)
    T = cblocks._poly_terms(b)
    T = [[(k, j, mp.mpc(c.real, c.imag)) for (k, j, c) in Tx] for Tx in T]
    B = [dict() for _ in range(Nmax + 1)]; B[0][0] = mp.mpc(1)
    def eta_mult(g):
        out = {}
        for m, c in g.items():
            mu = l + m
            out[m + 1] = out.get(m + 1, 0) + c * (mu + 1) / (2 * mu + 1)
            out[m - 1] = out.get(m - 1, 0) + c * mu / (2 * mu + 1)
        return out
    def D1(g):
        out = {}
        for m, c in g.items():
            mu = l + m; f = mu * (mu + 1) / (2 * mu + 1)
            out[m + 1] = out.get(m + 1, 0) + c * f; out[m - 1] = out.get(m - 1, 0) - c * f
        return out
    for N in range(1, Nmax + 1):
        src = {}
        for kind, Tx in enumerate(T, 1):
            for (k, j, cf) in Tx:
                if k == 0 or k > N or not B[N - k]: continue
                g = B[N - k]; D = Delta + N - k
                if kind == 1: h = {m: c * (D * D - C + (l + m) * (l + m + 1)) for m, c in g.items()}
                elif kind == 2: h = {m: c * D for m, c in g.items()}
                else: h = D1(g)
                for _ in range(j): h = eta_mult(h)
                for m, c in h.items(): src[m] = src.get(m, 0) + cf * c
        D = Delta + N
        for m, c in src.items():
            E = D * D - 3 * D + (l + m) * (l + m + 1) - C
            B[N][m] = -c / E
    return B

def Ghat_mp(B, nu, r, eta):
    l = mp.mpc(-0.5, nu)
    M = len(B) + 1
    P = {0: mp.legenp(l, 0, eta, type=3), 1: mp.legenp(l + 1, 0, eta, type=3)}
    for m in range(1, M):
        mu = l + m; P[m + 1] = ((2 * mu + 1) * eta * P[m] - mu * P[m - 1]) / (mu + 1)
    for m in range(1, M + 1): P[-m] = mp.conj(P[m])
    s = 0; rn = mp.mpf(1)
    for N in range(len(B)):
        s += rn * sum(c * P[m] for m, c in B[N].items()); rn *= r
    return s

if __name__ == '__main__':
    nu = mp.mpf('1.3'); Nmax = 200
    for (rho, beta) in [(0.55, 1.4), (0.85, 0.6), (0.85, 2.2)]:
        rho, beta = mp.mpf(rho), mp.mpf(beta)
        r = rho * mp.exp(-beta); eta = mp.cosh(beta); rhb = r * mp.exp(-beta)
        z, zb = 4 * rho / (1 + rho)**2, 4 * rhb / (1 + rhb)**2; v = (1 - z) * (1 - zb)
        for b in (2, 4):
            Bp = recursion_mp(nu, b, Nmax); Bm = recursion_mp(nu, -b, Nmax)
            gp = Ghat_mp(Bp, nu, r, eta); gm = Ghat_mp(Bm, nu, r, eta)
            # double precision series for comparison
            gm_d = cblocks.Ghat(float(nu), -b, float(r), float(eta), Nmax)
            say(f'rho={float(rho)} beta={float(beta)} b={b}: |Ghat_-b (50 digits) / (v^b Ghat_b) - 1| = {mp.nstr(abs(gm/(v**b*gp)-1),3)};'
                f'  double-precision Ghat_-b rel. error {mp.nstr(abs(gm_d/gm-1),3)}')
    open(OUT, 'w').write('\n'.join(lines) + '\n')

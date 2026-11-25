"""
fullpos/fp_moments.py (FULLPOS_CERT.md §1): CERTIFIED angular moments Y_k, Z_k (and exact Q_k) of the endpoint-repaired k = 2
weight, in closed form, in Arb ball arithmetic.  Nothing is sampled: no quadrature of the RBF part at all.

  g = -(1/Lam) [ R(x) B(nu) + Rt(x) L(nu) ] theta(nu),   x = nu^2,
      R  = 147456 x / ((4x+1)(4x+9)(4x+49))       ( = 9216(4x+1)/((4x+9)(4x+49)) * Q(nu) ),
      Rt = 9216 (4x+1) / ((4x+9)(4x+49)),
      B  = sum_signed x_k exp(-(nu-c_k)^2/(2 s_k^2)),   L = sum_j y_j nu^{2j} (x+1) log(1+x) e^{-x/2}.
For a polynomial p(x) (exact, flint fmpq_poly):
  * p(x) 147456 x = q(x)(4x+1)(4x+9)(4x+49) + rem(x)  ->  int_0^oo p R B = sum_i q_i S_{2i} + sum_a c_a F_a,
      S_n = int_0^oo nu^n B = (1/2) sum_signed x_k m_n(c_k, s_k)  (exact Gaussian moments, recurrence m_{n+1} = c m_n + n s^2 m_{n-1}),
      F_a = int_0^oo B/(4 nu^2 + a) = (1/2) sum_signed x_k (pi/(4 b)) Re w((i b - c_k)/(s_k sqrt 2)),  b = sqrt(a)/2  (Faddeeva w,
            w(z) = exp(-z^2) erfc(-i z), Arb),  c_a = rem(-a/4) / prod_{a' != a} (a' - a).
  * p(x) 9216 (4x+1) = qt(x)(4x+9)(4x+49) + remt(x) ->  int_0^oo p Rt L = sum_i qt_i Lm_i + sum_a d_a LF_a,
      Lm_i = sum_j y_j (K_{i+j+1} + K_{i+j}),  K_n = int_0^oo nu^{2n} log(1+nu^2) e^{-nu^2/2}:
      K_{n+1} = (2n+1) K_n + 2 E_{n+1},  E_{n+1} = G_n - E_n,  G_n = (2n-1)!! sqrt(pi/2),  E_0 = (pi/2) e^{1/2} erfc(1/sqrt2)
      (integration by parts; exact), K_0 and LF_a by acb.integral (rigorous) on [0, 64] plus an explicit tail bound.
  * theta: |int p g (1 - theta)| bounded explicitly (Gaussian completing-the-square majorant), added to the radius.
Outputs cache/moments_<Dg>_K<K>.json : Yt_k = Y_k/lam^{k+1}, Zt_k = Z_k/lam^k (arb strings), Q_k exact, plus Y_k, Z_k.
Run: python3 fp_moments.py Dg K prec
"""
import sys, os, json, time
sys.dont_write_bytecode = True
import flint
from flint import arb, acb, fmpq, fmpq_poly, ctx
import fp_common as FC


def cdh0_fmpq(M):
    f = [fmpq_poly([0]), fmpq_poly([1])]
    X = fmpq_poly([0, 1])
    for m in range(1, M):
        t = -(X + fmpq(1, 4) - 2*m*m)*f[m] - (m - 1)**2*f[m - 1]
        f.append(t*fmpq(1, (m + 1)**2))
    return f


def cdh1_fmpq(N):
    X = fmpq_poly([0, 1])
    q = [fmpq_poly([1])]
    for n in range(N):
        prev = q[n - 1] if n >= 1 else fmpq_poly([0])
        t = -(X + fmpq(1, 4) + 4 - 2*(n + 2)**2)*q[n] - (n*n)*prev
        q.append(t*fmpq(1, (n + 3)**2))
    return q


def zy_fmpq(K):
    f = cdh0_fmpq(K + 5)
    Z, Y = [], []
    for k in range(K + 1):
        m = k + 2
        Z.append(12*f[m] - (m - 1)*f[m - 1] + (m + 1)*f[m + 1])
        Y.append(fmpq(12 - m*m, 2)*f[m] + fmpq(13 - 7*m, 2)*f[m - 1] + fmpq(13 + 7*m, 2)*f[m + 1]
                 + fmpq((m - 2)**2, 4)*f[m - 2] + fmpq((m + 2)**2, 4)*f[m + 2])
    return Z, Y, f


def A(q):
    from fractions import Fraction
    if isinstance(q, Fraction):
        q = fmpq(q.numerator, q.denominator)
    return arb(q)


class Moments:
    def __init__(self, Dg, nmax, verbose=True):
        W = FC.load_weight(Dg)
        self.Dg = Dg
        self.lam = W['lam']
        self.Lam = arb(int(W['lam'])).sqrt()
        self.nmax = nmax                       # need S_{2i}, i <= nmax ; K_n, n <= nmax + 4
        t0 = time.time()
        # signed centres
        sc = []
        for (c, s), x in zip(W['cs'], W['x']):
            sc.append((c, s, x))
            if c > 0:
                sc.append((-c, s, x))
        self.sc = sc
        # --- S_{2i}
        S = [arb(0)]*(nmax + 1)
        sq2pi = (2*arb.pi()).sqrt()
        for c, s, x in sc:
            ca, sa, xa = A(c), A(s), A(x)
            s2 = sa*sa
            m0 = sa*sq2pi; m1 = ca*m0
            mm = [m0, m1]
            for n in range(1, 2*nmax):
                mm.append(ca*mm[n] + n*s2*mm[n - 1])
            for i in range(nmax + 1):
                S[i] += xa*mm[2*i]
        self.S = [v/2 for v in S]
        # --- F_a
        self.F = {}
        sqrt2 = arb(2).sqrt()
        for a in (1, 9, 49):
            b = arb(a).sqrt()/2
            tot = arb(0)
            for c, s, x in sc:
                zeta = acb(-A(c), b)/(A(s)*sqrt2)
                w = (-(zeta*zeta)).exp()*(acb(0, -1)*zeta).erfc()
                tot += A(x)*w.real
            self.F[a] = tot*arb.pi()/(4*b)/2
        if verbose:
            print(f'  gaussian part: {len(sc)} signed centres, {time.time()-t0:.1f}s', flush=True)
        # --- log part
        t0 = time.time()
        self.y = [A(v) for v in W['y']]
        NK = nmax + 6
        G = [(arb.pi()/2).sqrt()]
        for n in range(1, NK + 2):
            G.append(G[-1]*(2*n - 1))
        E0 = arb.pi()/2*arb(arb(1)/2).exp()*(1/arb(2).sqrt()).erfc()
        E = [E0]
        for n in range(NK + 1):
            E.append(G[n] - E[n])
        K0 = self._integral(lambda z: (1 + z*z).log()*(-(z*z)/2).exp())
        Kn = [K0]
        for n in range(NK):
            Kn.append((2*n + 1)*Kn[n] + 2*E[n + 1])
        self.Kn = Kn
        self.Lm = [sum((self.y[j]*(Kn[i + j + 1] + Kn[i + j]) for j in range(3)), arb(0)) for i in range(nmax + 1)]
        self.LF = {}
        for a in (9, 49):
            ya = self.y
            def fL(z, a=a):
                x = z*z
                pol = ya[0] + ya[1]*x + ya[2]*x*x
                return pol*(x + 1)*(1 + x).log()*(-x/2).exp()/(4*x + a)
            self.LF[a] = self._integral(fL)
        if verbose:
            print(f'  log part: {time.time()-t0:.1f}s', flush=True)
        # theta-bound ingredients
        self._theta_setup()

    def _integral(self, f):
        """int_0^oo f, f even-ish analytic near [0, 64], |f| <= nu^8 e^{-nu^2/2} * 10 beyond 64 (all uses): rigorous."""
        def F(z, analytic):
            return f(z)
        v = acb.integral(F, 0, 64).real
        tail = arb(10)**4*(-arb(64)**2/2 + 40).exp()        # crude, ~ e^{-2000}
        return v + arb(0, tail.upper())

    def _theta_setup(self):
        """|1 - theta(nu)| <= erfc((nu_c - nu)/a) <= exp(-(nu_c-nu)^2/a^2) on 0 <= nu <= nu_c, <= 2 beyond.
        For centre (c, s): G_c(nu) e^{-(nu_c-nu)^2/a^2} = exp(-(nu_c-c)^2/(2s^2+a^2)) exp(-(nu-c*)^2/(2 S^2)),
          1/(2S^2) = 1/(2s^2) + 1/a^2,  c* = (c/(2s^2) + nu_c/a^2) 2S^2  (< nu_c).
        beyond nu_c: G_c <= exp(-(nu_c-c)^2/(4s^2)) exp(-(nu-c)^2/(4s^2)).
        Bound: int_0^oo |p| R |B| (1-theta) <= 84 * sum_k |x_k| [ e^{..} M(|p|; c*_max, S) + 2 e^{..} M(|p|; c_max, s sqrt2) ],
        M(|p|; c, s) = sum_i |p_i| E(c + sZ)^{2i} sqrt(2pi) s  (monotone in c >= 0, so c*_max = max over centres)."""
        nu_c = self.Lam + 100; a = arb(10)
        pre = {}
        for c, s, x in self.sc:
            if c < 0:
                continue                       # only nu >= 0 matters; |x| for -c is the same, and G_{-c} <= G_c on nu >= 0
            ca, sa = A(c), A(s)
            w1 = abs(A(x))*(-(nu_c - ca)**2/(2*sa*sa + a*a)).exp()
            w2 = abs(A(x))*(-(nu_c - ca)**2/(4*sa*sa)).exp()
            S2 = 1/(2*(1/(2*sa*sa) + 1/(a*a)))
            cst = (ca/(2*sa*sa) + nu_c/(a*a))*2*S2
            d = pre.setdefault(s, [arb(0), arb(0), arb(0), ca, S2.sqrt(), sa*arb(2).sqrt()])
            d[0] += 2*w1; d[1] += 2*w2                     # factor 2: the -c twin
            if cst > d[2]:
                d[2] = arb(cst.upper())
            if ca > d[3]:
                d[3] = ca
        self.theta_pre = pre
        # log part: on [0, 550] 1-theta <= exp(-((nu_c-550)/10)^2); beyond 550 the Gaussian e^{-nu^2/2} with nu^{2n} is < e^{-100000}
        self.theta_log = (-((nu_c - 550)/10)**2).exp() + (-arb(100000)).exp()

    @staticmethod
    def gmom_abs(pabs, c, s, nmax):
        """sum_i pabs_i * int_R nu^{2i} e^{-(nu-c)^2/(2s^2)} (all terms >= 0)."""
        m0 = s*(2*arb.pi()).sqrt(); mm = [m0, c*m0]
        for n in range(1, 2*len(pabs)):
            mm.append(c*mm[n] + n*s*s*mm[n - 1])
        return sum((pabs[i]*mm[2*i] for i in range(len(pabs))), arb(0))

    def theta_bound(self, p):
        pabs = [abs(A(v)) for v in p.coeffs()]
        tot = arb(0)
        for s, (w1, w2, cst, cmax, S, s2) in self.theta_pre.items():
            tot += w1*self.gmom_abs(pabs, cst, S, None) + 2*w2*self.gmom_abs(pabs, cmax, s2, None)
        tot = 84*tot
        # log part: int |p| Rt |L| <= 189 sum_i |p_i| sum_j |y_j| (K_{i+j+1} + K_{i+j})
        lt = sum((pabs[i]*sum((abs(self.y[j])*(self.Kn[i + j + 1] + self.Kn[i + j]) for j in range(3)), arb(0))
                  for i in range(len(pabs))), arb(0))
        return (tot + 189*lt*self.theta_log)/self.Lam

    def integral(self, p):
        """int_0^oo p(nu^2) g(nu) dnu  (arb ball, rigorous)."""
        D = fmpq_poly([1, 4])*fmpq_poly([9, 4])*fmpq_poly([49, 4])
        q, rem = divmod(p*fmpq_poly([0, 147456]), D)
        val = arb(0)
        qc = q.coeffs()
        assert len(qc) <= self.nmax + 1, (len(qc), self.nmax)
        for i, ci in enumerate(qc):
            val += A(ci)*self.S[i]
        for a in (1, 9, 49):
            den = 1
            for a2 in (1, 9, 49):
                if a2 != a:
                    den *= (a2 - a)
            val += A(rem(fmpq(-a, 4))/den)*self.F[a]
        Dt = fmpq_poly([9, 4])*fmpq_poly([49, 4])
        qt, remt = divmod(p*fmpq_poly([9216, 4*9216]), Dt)
        for i, ci in enumerate(qt.coeffs()):
            val += A(ci)*self.Lm[i]
        for a in (9, 49):
            a2 = 58 - a
            val += A(remt(fmpq(-a, 4))/(a2 - a))*self.LF[a]
        val = -val/self.Lam
        tb = self.theta_bound(p)
        return val + arb(0, tb.upper())


def main():
    Dg = int(sys.argv[1]); K = int(sys.argv[2]); ctx.prec = int(sys.argv[3])
    t0 = time.time()
    Z, Y, f = zy_fmpq(K)
    print(f'polys K={K}: {time.time()-t0:.1f}s', flush=True)
    M = Moments(Dg, K + 4)
    lam = int(M.lam)
    out = dict(Dg=Dg, K=K, prec=ctx.prec, Y=[], Z=[], Yt=[], Zt=[], F=[])
    lamA = arb(lam)
    for k in range(K + 1):
        Zk = M.integral(Z[k]); Yk = M.integral(Y[k])
        out['Z'].append(Zk.str(40, radius=True)); out['Y'].append(Yk.str(40, radius=True))
        out['Zt'].append((Zk/lamA**k).str(40, radius=True)); out['Yt'].append((Yk/lamA**(k + 1)).str(40, radius=True))
        if k % 10 == 0 or k == K:
            print(f'k={k:4d} Z={Zk.str(12, radius=True)}  Y={Yk.str(12, radius=True)}  Y/Z={(Yk/Zk).str(8)}', flush=True)
    for m in range(1, min(K + 4, 40)):
        out['F'].append(M.integral(f[m]).str(30, radius=True))
    Qs = FC.Q_moments(K, M.lam) if K <= 60 else None
    q = cdh1_fmpq(K)
    _, p3 = FC.nf02()
    mom = {}
    def Mi(i):
        if i not in mom:
            mom[i] = sum((fmpq(c.numerator, c.denominator)/(2*i + n + 1) for n, c in enumerate(p3)), fmpq(0))
        return mom[i]
    Q = []
    for k in range(K + 1):
        v = fmpq(0)
        for i, c in enumerate(q[k].coeffs()):
            v += c*fmpq(lam)**i*Mi(i)
        Q.append(v)
    if Qs is not None:
        assert all(fmpq(a.numerator, a.denominator) == b for a, b in zip(Qs, Q))
    out['Q'] = [str(v) for v in Q]
    out['Qt_float'] = [float(arb(v)/lamA**k) for k, v in enumerate(Q)]
    out['time'] = time.time() - t0
    json.dump(out, open(os.path.join(FC.CACHE, f'moments_{Dg}_K{K}.json'), 'w'))
    print(f'done {time.time()-t0:.1f}s')


if __name__ == '__main__':
    main()

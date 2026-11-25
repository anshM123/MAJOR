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
        q.append(t*fmpq(1, (n + 4)**2))
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
    NGL = int(os.environ.get("FP_NGL", "300"))

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
        t0 = time.time()
        self._theta_setup()
        if verbose:
            print(f'  theta edge quadrature: {time.time()-t0:.1f}s', flush=True)

    def _integral(self, f):
        """int_0^oo f, f even-ish analytic near [0, 64], |f| <= nu^8 e^{-nu^2/2} * 10 beyond 64 (all uses): rigorous."""
        def F(z, analytic):
            return f(z)
        v = acb.integral(F, 0, 64).real
        tail = arb(10)**4*(-arb(64)**2/2 + 40).exp()        # crude, ~ e^{-2000}
        return v + arb(0, tail.upper())

    def _theta_setup(self):
        """theta correction, int_0^oo p R B (1-theta) dnu  (+ the log part, bounded).
        * EDGE centres c in [Lam-100, cmax] (c > 0 member only): Gauss-Legendre on 11 panels of width 30 covering
          [Lam-150, Lam+180], n = NGL nodes each; rigorous error (64/15) M rho^{-2n}/(rho^2-1) (L/2) (Trefethen, ATAP Thm 19.3)
          per monomial x^i, with M = sup |h| on the Bernstein ellipse E_rho (Arb boxes covering its boundary) times Xmax^i.
          Outside the panels: |G_c(nu)| <= e^{-d^2/(4s^2)} e^{-(nu-c)^2/(4s^2)}, d >= 50 (left; also 1-theta <= e^{-625})
          or d >= 150 (right).
        * all other members (twins -c, and 0 <= c < Lam-100): bounded via 1 - theta(nu) <= exp(-(nu_c-nu)^2/a^2) (nu <= nu_c), <= 2
          beyond, completing the square (see _rest_bound).  R <= 84 on the real line; |p| <= sum |p_i| x^i."""
        nu_c = self.Lam + 100; a = arb(10)
        self.nu_c = nu_c
        Lam = self.Lam
        cut = Lam - 100
        edge = []; rest = {}; twins = {}
        for c, s, x in self.sc:
            ca, sa, xa = A(c), A(s), A(x)
            if c > 0 and ca > cut:
                edge.append((ca, sa, xa))
            elif c >= 0:
                d = rest.setdefault(s, [arb(0), arb(0), arb(0), arb(0)])
                S2 = 1/(2*(1/(2*sa*sa) + 1/(a*a)))
                cst = (ca/(2*sa*sa) + nu_c/(a*a))*2*S2
                d[0] += abs(xa)*(-(nu_c - ca)**2/(2*sa*sa + a*a)).exp()
                d[1] += abs(xa)*(-(nu_c - ca)**2/(4*sa*sa)).exp()
                if cst > d[2]: d[2] = arb(cst.upper())
                if ca > d[3]: d[3] = ca
            if c < 0:
                d = twins.setdefault(s, [arb(0), arb(0)])
                d[0] += abs(xa)
                d[1] = arb(max(float(d[1].upper()), float(-c)))
        self.rest = rest; self.twins = twins; self.edge = edge
        self.cmax_edge = max(e[0] for e in edge)
        self.sum_edge = {}
        for ca, sa, xa in edge:
            self.sum_edge[float(sa.mid())] = self.sum_edge.get(float(sa.mid()), arb(0)) + abs(xa)
        # ---- panels
        NGL = self.NGL
        nodes = [arb.legendre_p_root(NGL, k, weight=True) for k in range(NGL)]
        lo = Lam - 150; Lp = arb(30); half = Lp/2
        rho = (arb(4)/3 + (arb(16)/9 + 4).sqrt())/2            # b = (L/2)(rho - 1/rho)/2 = 10
        nmax = self.nmax
        T = [arb(0)]*(nmax + 1); Err = [arb(0)]*(nmax + 1)
        prec0 = ctx.prec
        for P in range(11):
            mid = lo + half*(2*P + 1)
            # sup |h| on the ellipse boundary (low precision boxes)
            ctx.prec = 128
            Mh = arb(0); Xm = arb(0); NB = 256
            for j in range(NB):
                ph = 2*arb.pi()*(j + arb(1)/2)/NB
                e = acb(ph.cos(), ph.sin())
                z0 = mid + half*(rho*e + 1/(rho*e))/2
                rad = half*(rho + 1/rho)/2*(arb.pi()/NB)*arb('1.0001')
                z = acb(arb(z0.real.mid(), rad.upper()), arb(z0.imag.mid(), rad.upper()))
                hv = abs(self._h(z, edge, nu_c, a))
                Mh = Mh if Mh > hv else arb(hv.upper())
                xv = abs(z*z)
                Xm = Xm if Xm > xv else arb(xv.upper())
            ctx.prec = prec0
            fac = arb(64)/15*Mh*rho**(-2*NGL)/(rho*rho - 1)*half
            xp = arb(1)
            for i in range(nmax + 1):
                Err[i] += fac*xp; xp *= Xm
            for (t, w) in nodes:
                nu = mid + half*t
                v = self._h(nu, edge, nu_c, a)*w*half
                x = nu*nu
                for i in range(nmax + 1):
                    T[i] += v; v *= x
        self.T = T; self.TErr = Err
        # log part: on [0, 550] 1-theta <= exp(-((nu_c-550)/10)^2); beyond 550 the Gaussian e^{-nu^2/2} with nu^{2n} is < e^{-100000}
        self.theta_log = (-((nu_c - 550)/10)**2).exp() + (-arb(100000)).exp()

    def _h(self, z, edge, nu_c, a):
        """R(z^2) B_edge(z) (1 - theta(z)), z arb or acb."""
        x = z*z
        R = 147456*x/((4*x + 1)*(4*x + 9)*(4*x + 49))
        B = 0
        for ca, sa, xa in edge:
            B += xa*(-(z - ca)**2/(2*sa*sa)).exp()
        if isinstance(z, acb):
            th1 = ((z + nu_c)/a).erfc() + ((nu_c - z)/a).erfc()
        else:
            th1 = ((z + nu_c)/a).erfc() + ((nu_c - z)/a).erfc()
        return R*B*th1/2

    _GM = {}

    def gmom_abs(self, pabs, c, s, nmax=None):
        """sum_i pabs_i * int_R nu^{2i} e^{-(nu-c)^2/(2s^2)} (all terms >= 0); moment arrays cached up to 2(nmax+1)."""
        key = (c.mid().str(30), s.mid().str(30))
        if key not in self._GM:
            m0 = s*(2*arb.pi()).sqrt(); mm = [m0, c*m0]
            for n in range(1, 2*self.nmax + 4):
                mm.append(c*mm[n] + n*s*s*mm[n - 1])
            self._GM[key] = mm[0::2]
        mm = self._GM[key]
        return sum((pabs[i]*mm[i] for i in range(len(pabs))), arb(0))

    def theta_bound(self, p):
        """rigorous bound on everything in int p g (1-theta) except the edge quadrature value sum p_i T_i (times 1/Lam)."""
        pabs = [abs(A(v)) for v in p.coeffs()]
        tot = arb(0)
        sq2 = arb(2).sqrt()
        for s, (w1, w2, cst, cmax) in self.rest.items():
            sa = A(s); S = (1/(2*(1/(2*sa*sa) + 1/arb(100)))).sqrt()
            tot += 84*(w1*self.gmom_abs(pabs, cst, S) + 2*w2*self.gmom_abs(pabs, cmax, sa*sq2))
        nu_c = self.nu_c
        for s, (sx, cm) in self.twins.items():
            sa = A(s)
            tot += 84*sx*((-arb(2750)).exp()*self.gmom_abs(pabs, cm, sa) +
                          (-(nu_c/2)**2/(4*sa*sa)).exp()*self.gmom_abs(pabs, arb(0), sa*sq2))
        for s, sx in self.sum_edge.items():
            sa = A(s)
            tot += 84*sx*((-arb(625) - arb(2500)/(4*sa*sa)).exp() + 2*(-arb(22500)/(4*sa*sa)).exp())* \
                self.gmom_abs(pabs, self.cmax_edge, sa*sq2)
        tot += sum((pabs[i]*self.TErr[i] for i in range(len(pabs))), arb(0))
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
        pc = p.coeffs()
        val += sum((A(pc[i])*self.T[i] for i in range(len(pc))), arb(0))/self.Lam     # + (1/Lam) int p R B_edge (1-theta)
        tb = self.theta_bound(p)
        return val + arb(0, tb.upper())


def main():
    Dg = int(sys.argv[1]); K = int(sys.argv[2]); ctx.prec = int(sys.argv[3])
    t0 = time.time()
    M = Moments(Dg, K + 4)
    lam = int(M.lam)
    out = dict(Dg=Dg, K=K, prec=ctx.prec, Y=[], Z=[], F=[])
    nd = int(ctx.prec/3.33)
    # rolling window of exact CDH polynomials f_{m-2..m+2}
    X = fmpq_poly([0, 1])
    f = [fmpq_poly([0]), fmpq_poly([1])]
    def nextf():
        m = len(f) - 1
        f.append((-(X + fmpq(1, 4) - 2*m*m)*f[m] - (m - 1)**2*f[m - 1])*fmpq(1, (m + 1)**2))
    while len(f) < 5:
        nextf()
    for m in range(1, 5):
        out['F'].append(M.integral(f[m]).str(60, radius=True, more=True))
    for k in range(K + 1):
        m = k + 2
        while len(f) <= m + 2:
            nextf()
        Zp = 12*f[m] - (m - 1)*f[m - 1] + (m + 1)*f[m + 1]
        Yp = (fmpq(12 - m*m, 2)*f[m] + fmpq(13 - 7*m, 2)*f[m - 1] + fmpq(13 + 7*m, 2)*f[m + 1]
              + fmpq((m - 2)**2, 4)*f[m - 2] + fmpq((m + 2)**2, 4)*f[m + 2])
        f[m - 2] = None                                   # free memory
        Zk = M.integral(Zp); Yk = M.integral(Yp)
        out['Z'].append(Zk.str(nd, radius=True, more=True)); out['Y'].append(Yk.str(nd, radius=True, more=True))
        if k % 50 == 0 or k == K:
            print(f'k={k:4d} Z={Zk.str(12)} rel.rad {float((Zk.rad()/abs(Zk.mid()))):.2e}  Y={Yk.str(12)} rel.rad '
                  f'{float((Yk.rad()/abs(Yk.mid()))):.2e}  Y/Z={(Yk/Zk).str(8)}  t={time.time()-t0:.0f}s', flush=True)
    out['time'] = time.time() - t0
    json.dump(out, open(os.path.join(FC.CACHE, f'moments_{Dg}_K{K}.json'), 'w'))
    print(f'done {time.time()-t0:.1f}s')


if __name__ == '__main__':
    main()

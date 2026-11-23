"""deshift/ab_solve.py -- DESHIFT_CLOSE.md §2: continuum de-shift weights in an ANALYTIC basis, solved in the SHIFTED (difference-equation)
form, so that the de-shift identity W[w](t) = Phi(t) is controlled UNIFORMLY in t >= 0.

Set-up (C4: EE3, j = +-1, +-3; C5: C5_pp, j = 0, +-2).  A weight u(nu) (C4: chi, even and real on R; C5: D := c/Nhat5, Hermitian
D(-nu) = conj D(nu), using conj R_j(x) = R_j(-x)) has de-shift transform
        W[u](t) = 1/2 sum_j int_R u(nu) R_j(nu) cos((nu - i j) t) dnu          (= Stage 4's / Stage 5's W on nu > 0).
Shift each j-term to mu = nu - i j (contour R - i j -> R).  If u R_j is analytic in the crossed strip (u has ZEROS at the crossed poles of
R_j) and u has a pole only at nu = i j_max sgn(j) (the "target pole", mu = 0), then for all real t
        W[u](t) = Res-part(t) + 1/2 int_R E(mu) cos(mu t) dmu,   E = sum_j f_j^reg,  f_j(mu) = u(mu + i j) R_j(mu + i j),
where f_j^reg = f_j - S_j, S_j = (Gaussian-regularised principal part at mu = 0), and Res-part(t) = sum_j (i pi s_j / 2)(a1_j - a3_j t^2/2)
(s_j = -1 for j < 0 [contour above the pole], +1 for j > 0) with a_k the Laurent coefficients of f_j at 0 (b2 := a2 must vanish).
Imposing Res-part = Phi exactly gives
        |W[u](t) - Phi(t)| <= int_0^oo |E_even(mu)| dmu =: eps_W      for ALL t >= 0,
        |image(beta) - target(beta)| <= eps_W P_{-1/2}(cosh beta)      for ALL beta (|P_{-1/2+i mu}| <= P_{-1/2}),
        |m_n - m_n^target| <= int_0^oo |E_even| mu^{2n} dmu            (light-side moments).
Basis: u(nu) = Q(nu) B(nu), Q a fixed rational prefactor carrying the forced zeros/poles, B a Gaussian RBF expansion (centres k h,
width sigma; even (C4) or Hermitian (C5)).  Unknown real coefficients from a weighted least-squares problem in 384-bit ball arithmetic
(python-flint): E-rows on a mu-grid, constraint rows (Res-part = Phi, a2 = 0), and weak GAUGE rows pulling the weight towards the stored
Stage-7 XP weight (so the new weight is the same functional up to the stored weight's own accuracy).
  piece   Q(nu)                                              target pole            Phi
  c4v     nu^2 (nu^2+1)(nu^2+4)/[(nu^2+9)(4nu^2+1)(4nu^2+25)]    simple at +-3i     1
  c4d     nu^2 (nu^2+1)(nu^2+4)/[(nu^2+9)^3(4nu^2+1)(4nu^2+25)]  triple at +-3i     -t^2/2
  times 1/(4nu^2 + (3+2sqrt2)^2) for c4v, c4d.
  (the poles at +-i/2, +-5i/2, +-i(3+2sqrt2)/2 sit at zeros of every crossed R_j, so u R_j stays analytic in every crossed strip; they
   are present in the stored XP weights: AAA rational fits of chi_XP on the real axis locate +-0.5000i, +-2.4988i and a cluster
   around +-2.95i = the pair {2.914i, 3i})
  c5      i nu (nu^2+1)/(nu^2+4)                             simple at +-2i         1
Usage: python3 ab_solve.py <piece> [h sigma K lam]  -> cache/ab_<piece>.json (coefficients, 60-digit strings), cache/ab_solve.txt (log)
"""
import os, sys, time, json
import numpy as np
from flint import acb, arb, arb_mat, ctx
import kern

ctx.prec = int(os.environ.get('AB_PREC', '384'))
LOG = os.path.join(kern.CACHE, 'ab_solve.txt')
PI = arb.pi()
I = acb(0, 1)


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True)
    with open(LOG, 'a') as f:
        f.write(s + "\n")


class Piece:
    def __init__(self, name, h=0.5, sigma=1.2, K=120):
        self.name = name
        self.h, self.sigma, self.K = arb(h), arb(sigma), K
        self.s2 = 2*self.sigma**2
        if name in ('c4v', 'c4d'):
            self.R = kern.load_C4(); self.sym = 'even'; self.jp = 3
            self.p = 1 if name == 'c4v' else 3
            self.Phi = (1, 0) if name == 'c4v' else (0, -0.5)          # (constant, t^2 coefficient)
        else:
            self.R = kern.load_C5(); self.sym = 'herm'; self.jp = 2; self.p = 1; self.Phi = (1, 0)
            self.Nhat = kern.Rat(sum(r.expr for r in self.R.values()))
        self.js = sorted(self.R)
        self.k29 = (3 + 2*arb(2).sqrt())**2          # 4 nu^2 + (3 + 2 sqrt 2)^2 = 0 at nu = +-2.9142i, zeros of R_{-+3}
        # unknown labels
        self.lab = [('g', 0)]
        for k in range(1, K + 1):
            self.lab.append(('re', k))
            if self.sym == 'herm':
                self.lab.append(('im', k))
        self.n = len(self.lab)

    # prefactor
    def Q(self, z):
        if self.sym == 'even':
            z2 = z*z
            return z2*(z2 + 1)*(z2 + 4)/((z2 + 9)**self.p*(4*z2 + 1)*(4*z2 + 25)*(4*z2 + self.k29))
        z2 = z*z
        return I*z*(z2 + 1)/(z2 + 4)

    def basis(self, z):
        """list of the n basis functions B_i(z) (without Q)."""
        g0 = (-(z*z)/self.s2).exp()
        e = (z*self.h/(self.sigma**2)).exp(); ei = 1/e
        out = [g0]
        ep, em = acb(1), acb(1)
        for k in range(1, self.K + 1):
            ep *= e; em *= ei
            dk = (-(self.h*k)**2/self.s2).exp()
            gm = g0*ep*dk          # g(z - nu_k)
            gp = g0*em*dk          # g(z + nu_k)
            out.append(gm + gp)
            if self.sym == 'herm':
                out.append(I*(gm - gp))
        return out

    def u_basis(self, z):
        q = self.Q(z)
        return [q*b for b in self.basis(z)]

    def f_basis(self, j, mu):
        z = mu + I*j
        r = self.R[j](z)
        return [r*b for b in self.u_basis(z)]

    # Laurent coefficients a_1..a_p at mu = 0 of f_j (j = +- jp), per basis function, by a Cauchy circle
    def laurent(self, j, M=160, rad=0.2):
        rad = arb(rad)
        acc = [[acb(0)]*self.p for _ in range(self.n)]
        for m in range(M):
            w = (2*PI*I*m/M).exp()
            mu = rad*w
            fb = self.f_basis(j, mu)
            mup = mu**self.p
            for i in range(self.n):
                gval = fb[i]*mup
                for k in range(1, self.p + 1):       # a_k = Taylor coeff of order p-k of g = mu^p f
                    n_ = self.p - k
                    acc[i][k - 1] += gval/(w**n_)
        out = []
        for i in range(self.n):
            out.append([acc[i][k - 1]/(M*rad**(self.p - k)) for k in range(1, self.p + 1)])
        return out              # out[i][k-1] = a_k

    def S_basis(self, j, mu, La):
        """Gaussian-regularised principal part of f_j at 0 (per basis function), at real mu != 0."""
        e = (-(mu*mu)/2).exp()
        out = []
        for i in range(self.n):
            a = La[i]
            if self.p == 1:
                out.append(a[0]*e/mu)
            else:
                b3, b2, b1 = a[2], a[1], a[0] + a[2]/2
                out.append(e*(b3/mu**3 + b2/mu**2 + b1/mu))
        return out

    def E_basis(self, mu, LA):
        """E(mu) per basis function (sum over j of f_j^reg)."""
        tot = [acb(0)]*self.n
        for j in self.js:
            fb = self.f_basis(j, mu)
            if abs(j) == self.jp:
                sb = self.S_basis(j, mu, LA[j])
                fb = [a - b for a, b in zip(fb, sb)]
            tot = [a + b for a, b in zip(tot, fb)]
        return tot

    def Eeven_basis(self, mu, LA):
        a = self.E_basis(acb(mu), LA); b = self.E_basis(acb(-mu), LA)
        return [(x + y)/2 for x, y in zip(a, b)]

    def res_rows(self, LA):
        """(const, t^2) coefficients of the residue part, and the a2's, per basis function (complex)."""
        const = [acb(0)]*self.n; t2 = [acb(0)]*self.n; a2s = []
        for j in (-self.jp, self.jp):
            s = -1 if j < 0 else 1
            for i in range(self.n):
                a = LA[j][i]
                const[i] += I*PI*s/2*a[0]
                if self.p == 3:
                    t2[i] += I*PI*s/2*(-a[2]/2)
            if self.p == 3:
                a2s.append([LA[j][i][1] for i in range(self.n)])
        return const, t2, a2s

    def weight_basis(self, nu):
        """the STORED-convention weight per basis function at real nu: chi (C4) or c = D Nhat5 (C5)."""
        z = acb(nu)
        ub = self.u_basis(z)
        if self.sym == 'even':
            return ub
        nh = self.Nhat(z)
        return [b*nh for b in ub]


def lstsq(rows, rhs, n):
    """weighted least squares via normal equations in ball arithmetic (midpoints)."""
    A = arb_mat(len(rows), n); b = arb_mat(len(rows), 1)
    for r, (row, y) in enumerate(zip(rows, rhs)):
        for i in range(n):
            A[r, i] = row[i]
        b[r, 0] = y
    At = A.transpose()
    G = (At*A).mid(); g = (At*b).mid()
    x = G.solve(g, algorithm='approx')
    return [x[i, 0].mid() for i in range(n)]


def build_and_solve(P, mu_max=80.0, dmu=0.05, lam=1e-4, nug_max=60.0, Wc=1e10, xp=None):
    t0 = time.time()
    LA = {j: P.laurent(j) for j in (-P.jp, P.jp)}
    pr(f"  Laurent coefficients done [{time.time() - t0:.0f}s]")
    rows, rhs = [], []
    sq = arb(dmu).sqrt()
    mus = (np.arange(int(round(mu_max/dmu))) + 0.5)*dmu
    for mu in mus:
        e = P.Eeven_basis(arb(repr(float(mu))), LA)
        rows.append([x.real*sq for x in e]); rhs.append(arb(0))
        rows.append([x.imag*sq for x in e]); rhs.append(arb(0))
    pr(f"  E rows: {len(rows)} [{time.time() - t0:.0f}s]")
    const, t2, a2s = P.res_rows(LA)
    for vec, tgt in ((const, P.Phi[0]), (t2, P.Phi[1])):
        rows.append([x.real*Wc for x in vec]); rhs.append(arb(tgt)*Wc)
        rows.append([x.imag*Wc for x in vec]); rhs.append(arb(0))
    for vec in a2s:
        rows.append([x.real*Wc for x in vec]); rhs.append(arb(0))
        rows.append([x.imag*Wc for x in vec]); rhs.append(arb(0))
    # gauge rows
    ng = 0
    if xp is not None and lam > 0:
        nug, wg = xp
        sc = float(np.max(np.abs(wg)))
        dg = arb(nug[1] - nug[0]).sqrt()
        for nu, w in zip(nug, wg):
            if nu > nug_max:
                break
            vb = P.weight_basis(arb(repr(float(nu))))
            rows.append([x.real*lam*dg/sc for x in vb]); rhs.append(arb(repr(float(np.real(w))))*lam*dg/sc)
            ng += 1
            if P.sym == 'herm':
                rows.append([x.imag*lam*dg/sc for x in vb]); rhs.append(arb(repr(float(np.imag(w))))*lam*dg/sc)
    pr(f"  gauge rows: {ng} [{time.time() - t0:.0f}s]")
    x = lstsq(rows, rhs, P.n)
    pr(f"  solved ({P.n} unknowns) [{time.time() - t0:.0f}s]")
    # diagnostics on the training grid
    res_c = sum(c*xi for c, xi in zip(const, x)); res_t2 = sum(c*xi for c, xi in zip(t2, x))
    pr(f"  residue part: const = {res_c.real.str(12, radius=False)} {res_c.imag.str(3, radius=False)}i (target {P.Phi[0]}); "
       f"t^2 = {res_t2.real.str(12, radius=False)} {res_t2.imag.str(3, radius=False)}i (target {P.Phi[1]})"
       + ("; a2(-/+) = " + ", ".join(sum(c*xi for c, xi in zip(v, x)).str(3, radius=False) for v in a2s) if a2s else ""))
    return x, LA


def save(P, x, LA, extra):
    out = dict(piece=P.name, h=float(P.h.mid()), sigma=float(P.sigma.mid()), K=P.K, sym=P.sym, p=P.p, labels=P.lab,
               coeffs=[xi.str(60, radius=False) for xi in x], **extra)
    json.dump(out, open(os.path.join(kern.CACHE, f'ab_{P.name}.json'), 'w'), indent=0)


if __name__ == '__main__':
    piece = sys.argv[1]
    h, sigma, K, lam = 0.5, 1.2, 120, 1e-4
    if len(sys.argv) > 2:
        h, sigma, K, lam = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
    mu_max = float(os.environ.get('AB_MUMAX', '80')); dmu = float(os.environ.get('AB_DMU', '0.05'))
    P = Piece(piece, h, sigma, K)
    X = kern.load_xp()
    key = {'c4v': 'C4_value', 'c4d': 'C4_deriv', 'c5': 'C5'}[piece]
    pr(f"\n## ab_solve {piece}: h = {h}, sigma = {sigma}, K = {K} (centres to {h*K}), unknowns {P.n}, lam = {lam}, prec {ctx.prec}, "
       f"mu-grid (0, {mu_max}] step {dmu}")
    x, LA = build_and_solve(P, mu_max=mu_max, dmu=dmu, lam=lam, xp=(X[key]['nu'], X[key]['wf']))
    save(P, x, LA, dict(lam=lam, mu_max=mu_max, dmu=dmu))
    pr(f"  saved cache/ab_{piece}.json")

"""kb_h17.py -- certified (Arb ball arithmetic) H17, K and tau_0 from the computed |g3|^2 light functions (ads4/K_B2A4.md).

Inputs (derived in kb_celestial.py, kb_onepoint.py, kb_k2ratio.py):
  k = 3:  -C_3|_{g3}(nu) = -2piG|g3|^2 (nu/Lam)^4 g3(nu),   g3 = 1 - (35/2)/nu^2 + (569/16)/nu^4        (exact polynomial)
  k = 2:  per-nu |g3|^2 light value of Re X_theta Psi^_nu  =  theta chi_re(nu) P(nu)  (times 2piG|g3|^2/Lam^6, h1 units), with
          P(nu) = -(4nu^2+1) N10(nu) / (8 (4nu^2+9)(4nu^2+49)),
          N10   = 1024nu^10 - 195328nu^8 + 8229504nu^6 - 90780512nu^4 + 251603732nu^2 - 114141375,
          to be compared with the flat/h1 image theta chi_re (R_- + R_+) nu^6  (P/((R_-+R_+) nu^6) = 1 - (581/4)/nu^2 + ...).
So, exactly (K normalised as in K_CERT.md, K = K0 + dK0 + H17):
  H17_k3 = 2 int_0^1 psi3 p^4 (g3(Lam p) - 1) dp = (2/Lam^2) [ -(35/2) N1 + (569/16) N0 / Lam^2 ]                       (exact rational)
  H17_k2 = -(2/Lam^7) int_0^oo theta chi_re [ P(nu) - (R_-+R_+)(nu) nu^6 ] dnu      (continuum weight; acb.integral + rigorous tail)
         = -(2/Lam^7) dnu sum_k Re[X_k/(R_-+R_+)(nu_k)] [ P - (R_-+R_+) nu^6 ](nu_k)   (grid weight; exact finite sum)
Reads (read-only) ads4/code/kcert/kc_common.py, kc_dk0.py, kc_grid.py and the published data they read.
Usage: python3 -B kb_h17.py   -> cache/kb_h17.txt
"""
import sys, os, math, time
from fractions import Fraction as Fr
from flint import arb, acb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kcert'))
import kc_common as kc          # noqa: E402  (read-only import)
import kc_dk0, kc_grid          # noqa: E402
import sympy as sp

ctx.prec = 128
nu = sp.Symbol('nu')
N10 = 1024*nu**10 - 195328*nu**8 + 8229504*nu**6 - 90780512*nu**4 + 251603732*nu**2 - 114141375
PEXPR = -(4*nu**2 + 1)*N10/(8*(4*nu**2 + 9)*(4*nu**2 + 49))
B2, B0 = Fr(-35, 2), Fr(569, 16)


def phi_expr():
    Rm_s, Rp_s = kc.kernels()
    Rm = sp.sympify(Rm_s, locals={'nu': nu, 'I': sp.I}); Rp = sp.sympify(Rp_s, locals={'nu': nu, 'I': sp.I})
    S = sp.factor(sp.simplify(Rm + Rp))
    phi = sp.factor(sp.simplify(PEXPR - S*nu**6))
    # sanity: leading behaviour P/(S nu^6) = 1 - (581/4)/nu^2 + ...
    x = sp.Symbol('x')
    ser = sp.series(sp.simplify((PEXPR/(S*nu**6)).subs(nu, 1/x)), x, 0, 5).removeO()
    return phi, S, sp.expand(ser)


def h17_k3(Dg):
    p2, p3 = kc.nf02()
    N0, N1 = kc.poly_moment(p3, 0), kc.poly_moment(p3, 2)
    L2 = Fr((Dg - 2)*(Dg - 1))
    return 2*(B2*N1 + B0*N0/L2)/L2


def tail_bound_phi(C, L):
    """rigorous bound on int_{Lam+60}^oo |theta chi_re Phi|: |Phi| <= 5000 nu^6 for nu >= 300 (checked below), |theta| <= 1,
    |chi_re| <= |Q B| + |sing| with |Q| <= nu^-2 (same estimates as kc_dk0.tail_bound, there for 36 nu^m)."""
    T = L + 60; c = arb(float(C.c.max()))
    assert float(T - c) >= 29
    sx = arb(float(sum(abs(v) for v in C.x))); sy = arb(float(sum(abs(v) for v in C.y)))
    m = 4                      # 5000 nu^6 * nu^-2
    rate = 2*(T - c)/arb('1.62') - arb(m)/T
    g = T**m*(-(T - c)**2/arb('1.62')).exp()/rate
    hsing = 2*T**(m + 11)*(-(T*T)/2).exp()/(T - arb(m + 11)/T)
    return 5000*(2*sx*g + sy*hsing)


def h17_k2_continuum(Dg, phi_rat):
    C = kc.Chi(Dg); L = kc.Lam(Dg); Lf = float(L)

    def f(z, analytic):
        if z.rad() > 4:
            return kc.NAN
        ch = C(z, analytic)
        if not ch.is_finite():
            return kc.NAN
        return kc.theta(z, L)*ch*phi_rat(z)
    tol = arb(10.0)**(-8)*L**5          # integrand ~ 4648 nu^6 chi; result ~ Lam^5 x O(10): relative tolerance ~1e-9 or better
    tot = acb(0); x = 0.0
    while x < Lf + 60:
        y = min(x + 10.0, Lf + 60)
        tot += acb.integral(f, arb(x), arb(y), abs_tol=tol, eval_limit=10**7, deg_limit=200)
        x = y
    tb = tail_bound_phi(C, L)
    J = tot.real + arb(0, float(tb.upper()))
    return -2*J/L**7, J, tb


def h17_k2_grid(Dg, phi_rat):
    fn, dnu = kc_grid.FILES[Dg]
    L = kc.Lam(Dg)
    Rm_s, Rp_s = kc.kernels(); Rm, Rp = kc.Rat(Rm_s), kc.Rat(Rp_s)
    s = arb(0); n = 0
    for r in kc_grid.rows(fn):
        z = acb(arb(r[0])); X = acb(arb(r[1]), arb(r[2]))
        S = Rm(z) + Rp(z)
        s += ((X/S)*phi_rat(z)).real
        n += 1
    return -2*arb(kc.fq(dnu))*s/L**7, n


def up(x, sig=4):
    u = float(x.upper())*(1 + 1e-15)
    e = math.floor(math.log10(abs(u))) - sig + 1
    return '%.*e' % (sig - 1, math.ceil(u/10**e)*10**e)


def main():
    out = []

    def pr(s_):
        print(s_, flush=True); out.append(s_)
    t0 = time.time()
    phi, S, ser = phi_expr()
    pr('# kb_h17: certified H17 from the computed light functions')
    pr('R_- + R_+ = %s' % S)
    pr('P/((R_-+R_+) nu^6) = %s + O(nu^-6)   [h1 reading of the k=2 light function; a4 = coefficient of nu^-2]' % sp.simplify(ser.subs(sp.Symbol('x'), 1/nu)))
    pr('Phi := P - (R_-+R_+) nu^6 = %s' % phi)
    # check |Phi| <= 5000 nu^6 for nu >= 300 (rational bound: numerator/denominator coefficients)
    rat = sp.simplify(phi/nu**6)
    vals = [abs(float(rat.subs(nu, v))) for v in (300, 400, 1000, 10**4, 10**6)]
    assert max(vals) < 5000, vals
    lead = sp.limit(rat, nu, sp.oo)
    pr('Phi/nu^6 -> %s (check |Phi| <= 5000 nu^6 for nu >= 300: sample max %.2f; monotone rational, see K_B2A4 Sec. 5)' % (lead, max(vals)))
    phi_rat = kc.Rat(str(phi).replace('I', 'I'))
    p2, p3 = kc.nf02(); K0 = arb(kc.fq(kc.K0_exact()))
    rows = {}
    for Dg in (300, 1000):
        _, r = kc_dk0.main(Dg)
        Hk2, J, tb = h17_k2_continuum(Dg, phi_rat)
        rows[Dg] = dict(kind='continuum', dK0=r['dK0'], E=16*r['Lam']*r['V3'], Lam=r['Lam'], Hk2=Hk2, tb=tb, M2=r['M2'])
    _, g = kc_grid.main([300, 1000, 2000, 5000])
    for Dg in (300, 1000, 2000, 5000):
        Hk2g, n = h17_k2_grid(Dg, phi_rat)
        key = Dg if Dg in (2000, 5000) else ('%dg' % Dg)
        rr = g[Dg]
        rows[key] = dict(kind='grid', dK0=rr['dK0'], E=16*rr['Lam']*rr['V3'], Lam=rr['Lam'], Hk2=Hk2g, tb=None)
    pr('')
    pr('| Dg | weight | Lam^2 H17_k2 | Lam^2 H17_k3 | Lam^2 H17 | K = K0+dK0+H17 | E | tau0(K) | tau0 upper |')
    for key, r in rows.items():
        Dg = key if isinstance(key, int) else int(key[:-1])
        L = r['Lam']; L2 = L*L
        Hk3 = arb(kc.fq(h17_k3(Dg)))
        H = r['Hk2'] + Hk3
        K = K0 + r['dK0'] + H
        tau = 240/L**4*(r['E']/K).sqrt()
        r.update(H=H, K=K, tau=tau, Hk3=Hk3)
        pr('| %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (key, r['kind'], (r['Hk2']*L2).str(12, radius=True), (Hk3*L2).str(14),
                                                       (H*L2).str(10, radius=True), K.str(14, radius=True), r['E'].str(12),
                                                       tau.str(10, radius=True), up(tau, 5)))
    pr('')
    pr('leading-order check (K_CERT linear form, a4 = -581/4, b2 = -35/2): Lam^2 H17 ~ -2 M2 a4 + 2 N1 b2 =')
    N1 = kc.poly_moment(p3, 2)
    for Dg in (300, 1000):
        lin = -2*rows[Dg]['M2']*arb(kc.fq(Fr(-581, 4))) + 2*arb(kc.fq(N1))*arb(kc.fq(B2))
        pr('   Dg = %d: %s (exact: %s)' % (Dg, lin.str(8), (rows[Dg]['H']*rows[Dg]['Lam']**2).str(8)))
    pr('tail bounds beyond Lam+60 (continuum): 300: %s ; 1000: %s' % (rows[300]['tb'].str(3), rows[1000]['tb'].str(3)))
    targets = {300: ['3.135e-6', '3.14e-6', '3.15e-6'], 1000: ['2.892e-8', '2.90e-8'], 2000: ['1.929e-9', '1.93e-9'],
               5000: ['5.327e-11', '5.33e-11']}
    pr('')
    pr('margins: tau0 <= T needs H17 >= h*(T); computed H17 - h*(T) in units of 1/Lam^2:')
    for Dg, Ts in targets.items():
        r = rows[Dg]; L = r['Lam']; L2 = L*L
        for T in Ts:
            Tq = arb(T)
            Ks = (240/L**4)**2*r['E']/(Tq*Tq)
            pr('   Dg = %d: tau0 <= %s : K* = %s, K = %s, margin (K - K*) Lam^2 = %s' % (Dg, T, Ks.str(10), r['K'].str(12), ((r['K'] - Ks)*L2).str(6)))
    pr('total %.0fs' % (time.time() - t0))
    return out, rows


if __name__ == '__main__':
    o, _ = main()
    open(os.path.join(HERE, 'cache', 'kb_h17.txt'), 'w').write('\n'.join(o) + '\n')

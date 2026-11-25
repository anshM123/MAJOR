"""t0check/chk_lemma2.py (INDEP_CHECK_T0 §3): independent numerical tests of T0_CERT Lemma 2 (curvature in gap units).
  A. M1: the phi-form (Lemma 2 step 2) vs LRINF's t-form (sqrt2/pi) int_0^beta K c dt.
  B. pointwise: h(x) := [w r(tau) - r(x) - eps g1(x)]/eps^2 vs the majorant hbar(x) (sup r'' on [5x/6, x] by dense sampling),
     on x in (0, 2 Lam], several Lam (the lemma holds for every Lam > 0 with x <= 2 Lam).
  C. whole lemma: Lam^2 C_k computed DIRECTLY as int_0^1 psi(p)[P_{-1/2+i Lam p}(cosh beta) - J0(bhat p)] dp with P = 2F1 (mpmath hyp2f1,
     no Mehler-Dirichlet), vs a1 M1 + F_k/24; reports h_eff := Lam^2 (Lam^2 C - a1 M1 - F/24) against H_k (12124 / 12273).
mpmath, dps 40.  Not a proof; a [DV] test.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp
import t0_common as T

mp.mp.dps = 40
PS = {k: [mp.mpf(c.numerator)/c.denominator for c in T.PSI[k]] for k in (2, 3)}
CH2 = {k: [mp.mpf(c.numerator)/c.denominator for c in T.CHI2[k]] for k in (2, 3)}
A1 = {k: mp.mpf(T.A1[k].numerator)/T.A1[k].denominator for k in (2, 3)}
HK = {2: 12124, 3: 12273}


def trig_int(c, x, n):
    """int_0^1 p^n q(p) cos(px) dp (n even) / sin (n odd) as needed for Phi^(n): returns Phi^(n)(x) for q = poly c."""
    # Phi^(n)(x) = Re[ i^n int p^n q(p) e^{ipx} dp ]
    cn = [mp.mpf(0)]*n + list(c)
    if x < 2:
        s = mp.mpf(0)
        for j in range(0, 80):
            mj = sum(cc/(i + j + 1) for i, cc in enumerate(cn))
            s += mj*(1j*x)**j/mp.factorial(j)
        return mp.re((1j)**n*s)
    # exact IBP: int_0^1 f e^{ipx} = sum_j (-1)^j [f^(j) e^{ipx}/(ix)^{j+1}]_0^1
    f = cn; s = mp.mpc(0); j = 0
    while any(v != 0 for v in f):
        f1 = sum(f); f0 = f[0]
        s += (-1)**j*(f1*mp.expjpi(x/mp.pi) - f0)/(1j*x)**(j + 1)
        f = [i*f[i] for i in range(1, len(f))]; j += 1
    return mp.re((1j)**n*s)


def r_n(k, x, n):
    return trig_int(PS[k], x, n) + A1[k]*[1, -2, 6][n]/x**(2 + n)


def M1_phi(beta):
    S = mp.sinh(beta/2)
    def f(ph):
        sg = 2*S*mp.sin(ph); t = 2*mp.asinh(sg/2)
        return 1/sg**2 - 1/(t**2*mp.cosh(t/2))
    return 2/mp.pi*mp.quad(f, [0, mp.pi/4, mp.pi/2])


def M1_t(beta):
    D = mp.cosh(beta) - 1
    c = lambda t: mp.cosh(t/2)/(4*mp.sinh(t/2)**2) - 1/t**2
    return mp.sqrt(2)/mp.pi*mp.quad(lambda t: c(t)/mp.sqrt(D - (mp.cosh(t) - 1)), [0, beta/2, beta])


def Fk(k, b):
    return mp.quad(lambda p: sum(cc*p**i for i, cc in enumerate(CH2[k]))*mp.besselj(0, b*p), mp.linspace(0, 1, 2 + int(b/3)))


def hbar(k, x):
    r0, r1 = abs(r_n(k, x, 0)), abs(r_n(k, x, 1))
    s2 = max(abs(r_n(k, x*(mp.mpf(5)/6 + mp.mpf(i)/(6*40)), 2)) for i in range(41))
    return (3*x**5/640 + x**5/192)*r1 + 3*x**4/128*r0 + x**6/1152*s2


def h_point(k, x, Lam):
    eps = 1/Lam**2
    tau = 2*Lam*mp.asinh(x/(2*Lam)); w = 1/mp.sqrt(1 + x**2/(4*Lam**2))
    g1 = -(x**2/8)*r_n(k, x, 0) - (x**3/24)*r_n(k, x, 1)
    return (w*r_n(k, tau, 0) - r_n(k, x, 0) - eps*g1)/eps**2


def P_conical(nu, beta):
    return mp.re(mp.hyp2f1(mp.mpf(1)/2 - 1j*nu, mp.mpf(1)/2 + 1j*nu, 1, -mp.sinh(beta/2)**2))


def C_direct(k, bh, Lam):
    beta = 2*mp.asinh(bh/(2*Lam))
    psi = lambda p: sum(cc*p**i for i, cc in enumerate(PS[k]))
    nseg = 4 + int((Lam*beta + bh)/2)
    with mp.workdps(30):
        v = mp.quad(lambda p: psi(p)*(P_conical(Lam*p, beta) - mp.besselj(0, bh*p)), mp.linspace(0, 1, nseg))
    return v


def main():
    print('A. M1 phi-form vs t-form (LRINF definition):')
    for beta in (mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('1'), mp.mpf('1.76'), mp.mpf('4')):
        a, b = M1_phi(beta), M1_t(beta)
        print(f'   beta={float(beta):5.2f}: phi-form {mp.nstr(a, 15)}  t-form {mp.nstr(b, 15)}  diff {mp.nstr(a - b, 3)}')
    print('B. pointwise |h(x)| <= hbar(x):')
    worst = 0
    for k in (2, 3):
        for Lam in (mp.mpf(3), mp.mpf(10), mp.mpf(40), mp.mpf(300)):
            for x in [mp.mpf(v) for v in (0.05, 0.3, 1, 2.5, 5, 8, 13, 20, 35, 60, 100, 180, 300, 450, 600)]:
                if x > 2*Lam:
                    continue
                h = h_point(k, x, Lam); hb = hbar(k, x)
                worst = max(worst, abs(h)/hb)
                if abs(h) > hb:
                    print(f'   VIOLATION k={k} Lam={float(Lam)} x={float(x)}: h={mp.nstr(h, 6)} hbar={mp.nstr(hb, 6)}')
            print(f'   k={k} Lam={float(Lam):6.1f}: max |h|/hbar so far = {float(worst):.4f}')
    print('C. direct Lam^2 C_k (2F1, no Mehler-Dirichlet) vs a1 M1 + F/24; h_eff = Lam^2 (Lam^2 C - a1 M1 - F/24):')
    for Lam in (mp.sqrt(298*299), mp.sqrt(998*999), mp.mpf(40)):
        for bh in (mp.mpf(1), mp.mpf('2.5'), mp.mpf(8), mp.mpf(13), mp.mpf(19.6), mp.mpf(60), 2*Lam*mp.mpf('0.95')):
            if bh > 300 and Lam > 900:
                continue
            beta = 2*mp.asinh(bh/(2*Lam)); M1 = M1_phi(beta)
            row = []
            for k in (2, 3):
                C = C_direct(k, bh, Lam)
                pred = A1[k]*M1 + Fk(k, bh)/24
                heff = Lam**2*(Lam**2*C - pred)
                row.append(f'k={k}: L2C={mp.nstr(Lam**2*C, 10)} pred={mp.nstr(pred, 10)} h_eff={mp.nstr(heff, 5)} (|h_eff|/H={float(abs(heff)/HK[k]):.3f})')
            print(f'   Lam={float(Lam):8.2f} bhat={float(bh):8.2f}: ' + ' ; '.join(row), flush=True)


if __name__ == '__main__':
    main()

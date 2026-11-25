"""
contour/cb_common.py (CONTOUR_BOUNDS.md): Arb (python-flint) machinery for the analytically continued effective density of the
repaired k = 2 weight of Omega^EE at Delta_g = 1000, on the contour pieces of ANGULAR_ROUTE.md §38.

Objects.
  chi = chi_re of prl/data/k2/k2_continuum_D_1000_rep.json (the realised, real-analytic weight; ew_common.Weight reads x_re, y_re),
  theta(z) = (1/2)[erf((z + nu_c)/a) - erf((z - nu_c)/a)] = 1 - erfc((z + nu_c)/a)/2 - erfc((nu_c - z)/a)/2,  nu_c = Lam + 100, a = 10,
  G(nu) := theta(nu') chi(nu') R_-(nu'),  nu' = nu - i,
  f(nu) = G(nu) + conj(G(conj nu))       (= chi_th(nu-i)R_-(nu-i) + chi_th(nu+i)R_+(nu+i), since R_+(z) = conj R_-(conj z)),
  B(nu) = Gamma(2 i nu)/Gamma(1/2 + i nu)^2.
Branch.  In chi(nu'), Log(nu'^2 + 1) = Log(nu(nu - 2i)) is replaced by Log(nu) + Log(nu - 2i): for Re nu > 0 the two agree (both
  arguments in (-pi/2, pi/2), nu'^2 + 1 never on (-oo, 0]), and the right-hand side is continuous up to Re nu = 0 away from nu = 0, 2i.
  So this is the analytic continuation from Re nu > 0, with boundary values from Re nu > 0 on the imaginary axis (as in §38).
  conj G(conj nu) is evaluated at conj nu, which also has Re >= 0, so the same formula serves.
Exact kernel algebra (asserted with sympy in ew_common.check_kernels): Q R_- = -32 nu'(2nu'+3i)(2nu'+5i)(2nu'+11i)/((2nu'-i)(2nu'-7i)
  (2nu'-3i)(2nu'+7i)), so  chi R_- = [Q R_-] Bg + Y l R_-,  Bg = sum_k x_k g_k(nu' - c_k),  Y = sum_j y_j nu'^{2j} e^{-nu'^2/2}.
Singularities of G (in nu): poles at 3i/2, 5i/2, 9i/2, -5i/2 (from Q R_- and R_-); nu = i is removable; branch points 0, 2i.
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'epswcert'))
import ew_common as C                    # read-only import
C.CACHE = os.path.join(HERE, 'cache')    # redirect (nothing is written through C anyway)
from ew_common import horner, QRN, QRD, RMN, RMD, amax, up
from flint import arb, acb, acb_series, acb_poly, arb_poly, ctx

ctx.prec = 192
CACHE = os.path.join(HERE, 'cache')
I = acb(0, 1)
REP = os.path.join(ADS4, 'prl', 'data', 'k2', 'k2_continuum_D_1000_rep.json')
KCUT = 20          # RBF terms with |Re nu' - c_k| <= KCUT sigma_k are summed exactly; the rest are bounded


def load():
    assert C.check_kernels()
    W = C.Weight(1000, REP)
    W.nuc = W.Lam + 100
    W.acut = arb(10)
    W.sigs = sorted(W.sumabs.keys())
    return W


def near_terms(W, re):
    import bisect
    lo = bisect.bisect_left(W.cs_sorted, re - KCUT*0.9 - 1)
    hi = bisect.bisect_right(W.cs_sorted, re + KCUT*0.9 + 1)
    out = []
    for i in range(lo, hi):
        c, s, x = W.sc[i]
        if abs(re - c) <= KCUT*s:
            ca, sa, xa, inv = W.sc_arb[i]
            out.append((ca, xa, inv, c, s, x))
    return out


def far_bound(W, imabs, rho):
    """sup over the disc |nu' - w0| <= rho of sum_{far k} |x_k g_k(nu' - c_k)|, far = |Re w0 - c_k| > KCUT sigma_k, |Im w0| <= imabs.
    |g(z)| = exp(-(Re z^2)/(2s^2)) = exp(((Im z)^2 - (Re z)^2)/(2s^2)) <= exp(((imabs + rho)^2 - (KCUT s - rho)^2)/(2s^2)); the sum of
    |x_k| over all k (both signs of c_k) bounds the coefficient sum."""
    tot = arb(0)
    for s in W.sigs:
        sa = arb(s)
        assert KCUT*s - rho > 0
        tot += W.sumabs[s]*(((arb(imabs) + rho)**2 - (KCUT*sa - rho)**2)/(2*sa*sa)).exp()
    return tot


def theta_series(W, z):
    """theta(z) for acb or acb_series z."""
    return 1 - ((z + W.nuc)/W.acut).erfc()/2 - ((W.nuc - z)/W.acut).erfc()/2


def G_series(W, z0, N, near=None):
    """Taylor coefficients (length N) of G_near(z0 + s) at s = 0; z0 acb, Re z0 > 0 (or = 0 away from 0, 2i).  G_near uses the RBF
    terms near Re z0 (near_terms)."""
    S = acb_series([z0, 1], prec=N)
    nup = S - I
    if near is None:
        near = near_terms(W, float(z0.real.mid()))
    Bg = None
    for ca, xa, inv, *_ in near:
        D = nup - ca
        t = (-(D*D)*inv).exp()*xa
        Bg = t if Bg is None else Bg + t
    if Bg is None:
        Bg = S*0
    QR = horner(QRN, nup)/horner(QRD, nup)
    n2 = nup*nup
    Y = (W.y[0] + W.y[1]*n2 + W.y[2]*n2*n2)*((-n2/2).exp())
    l = (n2 + 1)*(S.log() + (S - 2*I).log())
    Rm = horner(RMN, nup)/horner(RMD, nup)
    G = theta_series(W, nup)*(QR*Bg + Y*l*Rm)
    return [G[n] for n in range(N)]


def ring(z0, rho, nb):
    out = []
    two_pi = 2*arb.pi()
    for k in range(nb):
        th = two_pi*(arb(k) + arb(0.5, 0.5))/nb
        out.append(z0 + acb(rho*th.cos(), rho*th.sin()))
    return out


def check_disc(z0, rho):
    """the closed disc |nu - z0| <= rho must avoid the singularities of G (poles 3i/2, 5i/2, 9i/2, -5i/2; branch cuts of Log nu
    (nu <= 0) and Log(nu - 2i) (Im nu = 2, Re nu <= 0))."""
    x = float(z0.real.mid()); y = float(z0.imag.mid())
    for p in (1.5, 2.5, 4.5, -2.5):
        assert abs(complex(x, y - p)) > rho*1.01 + 1e-9, (x, y, rho, p)
    assert x > rho*1.01 or (abs(y) > rho*1.01 and abs(y - 2) > rho*1.01), (x, y, rho)


def G_disc_sup(W, z0, rho, nb=64):
    """(MG, Mfar): MG >= sup_{|nu - z0| <= rho} |G_near(nu)|, Mfar >= sup of |G - G_near| on the same disc.
    |G_near| <= |theta QR| * sum_near |x_k| sup|g_k| + |theta Y l R_-|; the analytic factors theta QR and theta Y l R_- are bounded
    on the circle by a ring of nb boxes (maximum modulus; check_disc), sum_near term by term."""
    check_disc(z0, rho)
    rho_a = arb(rho)
    tq = arb(0); yl = arb(0)
    for nu in ring(z0, rho_a, nb):
        nup = nu - I
        th = theta_series(W, nup)
        tq = amax(tq, abs(th*horner(QRN, nup)/horner(QRD, nup)).upper())
        n2 = nup*nup
        Y = (W.y[0] + W.y[1]*n2 + W.y[2]*n2*n2)*((-n2/2).exp())
        yl = amax(yl, abs(th*Y*(n2 + 1)*(nu.log() + (nu - 2*I).log())*horner(RMN, nup)/horner(RMD, nup)).upper())
    re = z0.real; imabs = abs((z0 - I).imag) + rho_a
    sb = arb(0)
    for ca, xa, inv, c, s, x in near_terms(W, float(re.mid())):
        d = abs(re - ca) - rho_a
        d2 = d*d if d > 0 else arb(0)
        sb += abs(xa)*((imabs*imabs - d2)*inv).exp()
    MG = tq*sb + yl
    Mfar = tq*far_bound(W, abs((z0 - I).imag), rho)
    return arb(MG.upper()), arb(Mfar.upper())


def B_series(z0, N):
    S = acb_series([z0, 1], prec=N)
    Bs = (2*I*S).gamma()*((I*S + arb(1)/2).rgamma())**2
    return [Bs[n] for n in range(N)]


def B_val(z):
    """B at an acb point or box, via duplication and a shift by 5 (Arb's gamma of wide boxes near the poles returns nan):
    B(nu) = 2^{2 i nu - 1} pi^{-1/2} Gamma(z)/Gamma(z + 1/2), z = i nu,
    Gamma(z)/Gamma(z + 1/2) = Gamma(z + 5) rgamma(z + 11/2) prod_{j<5} (z + 1/2 + j)/(z + j)."""
    zz = I*z
    half = arb(1)/2
    r = (zz + 5).gamma()*(zz + 5 + half).rgamma()
    for j in range(5):
        r = r*(zz + half + j)/(zz + j)
    return (2*zz*arb(2).log()).exp()/2/arb.pi().sqrt()*r


def B_disc_sup(z0, rho, nb=64):
    """sup |B| on |nu - z0| <= rho (poles of B at nu = i m, m = 0, 1, 2, ...; the disc must avoid them)."""
    x = float(z0.real.mid()); y = float(z0.imag.mid())
    for m in range(0, 8):
        assert abs(complex(x, y - m)) > rho*1.01, (x, y, rho, m)
    M = arb(0)
    for nu in ring(z0, arb(rho), nb):
        M = amax(M, abs(B_val(nu)).upper())
    return M


def f_line_coeffs(W, z0, d, N):
    """Taylor coefficients in real s of f(z0 + d s), d = 1 (horizontal) or i (vertical):
       G(z0 + d s) = sum G_n(z0) d^n s^n;  conj G(conj z0 + conj(d) s) = sum conj(G_n(conj z0)) d^n s^n."""
    g1 = G_series(W, z0, N)
    g2 = G_series(W, z0.conjugate(), N)
    out = []
    dn = acb(1)
    for n in range(N):
        out.append(dn*(g1[n] + g2[n].conjugate()))
        dn = dn*d
    return out


def B_line_coeffs(z0, d, N):
    b = B_series(z0, N)
    out = []; dn = acb(1)
    for n in range(N):
        out.append(dn*b[n]); dn = dn*d
    return out


def poly_abs_bounds(p, hs):
    """p acb_poly in real s on |s| <= hs: (sup |p| upper, int |p| upper) as arbs (sum of |coefficients|)."""
    sup = arb(0); itg = arb(0); hp = arb(1); hs = arb(hs)
    for n, c in enumerate(p.coeffs()):
        a = arb(abs(c).upper())
        sup += a*hp
        itg += a*2*hp*hs/(n + 1)
        hp = hp*hs
    return arb(sup.upper()), arb(itg.upper())


def line_panel(W, z0, d, h, rho, N, Q, withB=True, nb=64):
    """One panel z0 + d s, |s| <= h (real s), cut into Q subcells.  Returns a list of Q dicts (s-interval, sup|f|, int|f|, and with
    withB also sup|B|, sup|fB|, int|fB|), every number an upper bound.
    Truncation: |f - p_f| <= e_f := (MG(z0) + MG(conj z0)) x^N/(1 - x) + Mfar(z0) + Mfar(conj z0), x = h/rho (Cauchy on the disc
    of radius rho; the far RBF part is not expanded and is bounded on the disc); |B - p_B| <= e_B := MB x^N/(1 - x)."""
    h = arb(h); x = h/rho
    assert x < 1
    M1, F1 = G_disc_sup(W, z0, rho, nb)
    M2, F2 = G_disc_sup(W, z0.conjugate(), rho, nb)
    ef = (M1 + M2)*x**N/(1 - x) + F1 + F2
    pf = acb_poly(f_line_coeffs(W, z0, d, N))
    if withB:
        MB = B_disc_sup(z0, rho, nb)
        eB = MB*x**N/(1 - x)
        pB = acb_poly(B_line_coeffs(z0, d, N))
        pfB = pf*pB
    hs = h/Q
    out = []
    for j in range(Q):
        cj = -h + (2*j + 1)*hs
        sf = pf(acb_poly([acb(cj), 1]))
        supf, intf = poly_abs_bounds(sf, hs)
        r = dict(s0=cj - hs, s1=cj + hs, supf=supf + ef, intf=intf + 2*hs*ef, ef=ef)
        if withB:
            sB = pB(acb_poly([acb(cj), 1])); sfB = pfB(acb_poly([acb(cj), 1]))
            supB, _ = poly_abs_bounds(sB, hs)
            supfB, intfB = poly_abs_bounds(sfB, hs)
            supB_ = supB + eB
            cross = ef*supB + eB*supf + ef*eB
            r.update(supB=supB_, supfB=supfB + cross, intfB=intfB + 2*hs*cross, eB=eB)
        out.append(r)
    return out


def f_point(W, z):
    """f(z) as an acb ball (series of length 1; near RBF terms only, far part bounded and added as a radius)."""
    g1 = G_series(W, z, 1)[0]
    g2 = G_series(W, z.conjugate(), 1)[0]
    return g1 + g2.conjugate()

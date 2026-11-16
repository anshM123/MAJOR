r"""
AdS side of the dictionary  ghat3  <->  bulk W^3 coefficient  <->  t4 (and alpha_T)  in d = 3.
See ads4/DICTIONARY_T4.md, sections 3-5.  Exact arithmetic (sympy) unless stated.

Checks:
 1. CEMZ (arXiv:1407.5597) AdS shock profile varpi_D(rho) (sec. 3.2, eq. \exactsolution) and its D=4,5
    closed forms (eq. \fourd); limits rho->0 (flat) and rho->1 (Hofman-Maldacena collider).
 2. CEMZ time-delay formula (sec. 3.3, eq. \timedelay): t4(rho), t2(rho).  Limit rho->1:
        t4 -> D (D^2-1)(D-2) alpha4 / L^4,      t2 -> (D-1)(D-2) (alpha2 - (4D-2) alpha4)   (alpha_i in units of L)
    and rho->0 with rho=b/2 reproduces the flat formula (App. B, eq. \energycorr) t4 = (D-4)(D-2)D(D+2) alpha4/b^4.
 3. Gauss-Bonnet t2 in general d (Buchel et al. 0911.4257, t2 formula for GB in general D) from alpha2 = 4 lambda_GB^{CEMZ}.
 4. t4 for quasi-topological gravity (d=4, Myers-Paulos-Sinha 1004.2055 (t2t4): 3780 mu) and for Einsteinian cubic
    gravity (d=3, Bueno-Cano-Ruiperez 1802.00018 (t4ecg): -1260 mu) from the flat 3-point map alpha4 = -2 alpha3,
    I2 = I1/2 on-shell (three_point_vertices.py).
 5. Algebraic identity for the 4D Weyl tensor: W_a^c_b^d W_c^e_d^f W_e^a_f^b = (1/2) W_ab^cd W_cd^ef W_ef^ab (random tensors).
 6. Hofman-Maldacena d=3 normalisation: <|e_ij n_i n_j|^2/e*e>_circle = 1/4 = 2/(d^2-1).
 7. C_T = 3 L^2/(pi^3 G) in d=3, c_T/c_T(free real scalar) = 32 L^2/(pi G); ABJM check.
 8. Conditional translation of the flat bound (numbers quoted in DICTIONARY_T4.md section 6).
"""
import sympy as sp
import numpy as np

rho, D = sp.symbols('rho D', positive=True)


def varpi(Dv):
    # CEMZ eq. \exactsolution without the D-independent prefactor 16 pi G |P_u|:
    # varpi(rho) = (1-rho^2) * (4pi)^{(2-D)/2} Gamma(D/2)/((D-1)(D-2)) (rho^2/(1-rho^2))^{2-D} 2F1(D-2, D/2; D; -(1-rho^2)/rho^2)
    pref = (4 * sp.pi) ** sp.Rational(2 - Dv, 2) * sp.gamma(sp.Rational(Dv, 2)) / ((Dv - 1) * (Dv - 2))
    return (1 - rho ** 2) * pref * (rho ** 2 / (1 - rho ** 2)) ** (2 - Dv) * sp.hyper([Dv - 2, sp.Rational(Dv, 2)], [Dv], -(1 - rho ** 2) / rho ** 2)


def t4_of_rho(w, Dv):
    k = (1 - rho ** 2) ** 2 * (sp.diff(w, rho) - rho * sp.diff(w, rho, 2)) / (4 * rho * w)
    return -Dv * k * (Dv * (1 + rho ** 2) ** 2 + 2 * (1 + rho ** 4)) / (4 * rho ** 2)   # times alpha4


def t2_of_rho(w, Dv, a2, a4):
    k = (1 - rho ** 2) ** 2 * (sp.diff(w, rho) - rho * sp.diff(w, rho, 2)) / (4 * rho * w)
    return k * (-a2 + a4 * (Dv * (1 + rho ** 2) ** 2 - 2 * rho ** 2) / rho ** 2)


def closed_forms():
    w4 = -8 * (1 - rho ** 2 + (1 + rho ** 2) * sp.log(rho))
    w5 = 2 * (1 - rho) ** 4 / rho
    return {4: w4, 5: w5}


def check_varpi():
    out = []
    cf = closed_forms()
    for Dv in (4, 5):
        w = varpi(Dv)
        r1 = [sp.N(w.subs(rho, x) / cf[Dv].subs(rho, x), 15) for x in (sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(4, 5))]
        out.append((f"varpi_{Dv}(hyper)/varpi_{Dv}(closed) at rho=0.2,0.5,0.8 (constant expected)", r1))
    return out


def check_limits():
    out = []
    a2, a4 = sp.symbols('alpha2 alpha4')
    for Dv in (4, 5, 6, 7):
        w = closed_forms()[Dv] if Dv in (4, 5) else varpi(Dv)
        if Dv not in (4, 5):
            w = sp.simplify(sp.hyperexpand(w))
        # rho -> 1
        t4_1 = sp.limit(sp.simplify(t4_of_rho(w, Dv)), rho, 1, '-')
        t2_1 = sp.limit(sp.simplify(t2_of_rho(w, Dv, a2, a4)), rho, 1, '-')
        out.append((f"D={Dv}: t4(rho->1)/alpha4  [expect D(D^2-1)(D-2) = {Dv*(Dv**2-1)*(Dv-2)}]", sp.nsimplify(t4_1)))
        out.append((f"D={Dv}: t2(rho->1)  [expect (D-1)(D-2)(alpha2-(4D-2)alpha4)]", sp.expand(t2_1)))
        # rho -> 0 : compare with flat t4 = (D-4)(D-2)D(D+2) alpha4 / b^4, b = 2 rho (only D>4: in D=4 both sides
        # are normalised by a log-divergent leading term)
        if Dv > 4:
            lead = sp.limit(sp.simplify(t4_of_rho(w, Dv) * 16 * rho ** 4), rho, 0, '+')
            out.append((f"D={Dv}: lim rho->0 of 16 rho^4 t4(rho)/alpha4 [flat: (D-4)(D-2)D(D+2) = {(Dv-4)*(Dv-2)*Dv*(Dv+2)}]", lead))
        else:
            # D=4: t4(rho)*varpi(rho) (absolute size of the helicity-flip term) ~ const * alpha4 / rho^4
            lead = sp.limit(sp.simplify(t4_of_rho(w, 4) * w * rho ** 4), rho, 0, '+')
            out.append(("D=4: lim rho->0 of rho^4 t4(rho) varpi(rho)/alpha4 (finite: helicity-flip term ~ alpha4/b^4)", lead))
    return out


def check_GB():
    # CEMZ GB action (1/16piG)(R + lam_C X4); Buchel et al. (1/16piG)(... + L^2 lam/((D-3)(D-4)) X4)  => lam_C = lam L^2/((D-3)(D-4))
    # flat: t2 = 4 lam_C (D-4)(D-2)/b^2 [CEMZ sec. 2, eq. \totaimeGB] and t2 = (D-2)(D-4) alpha2/b^2 [\energycorr] => alpha2 = 4 lam_C
    lam = sp.symbols('lambda')
    res = []
    for Dv in (5, 6, 7, 8):
        lamC = lam / ((Dv - 3) * (Dv - 4))
        t2 = (Dv - 1) * (Dv - 2) * 4 * lamC
        buchel = 4 * lam * (Dv - 1) * (Dv - 2) / ((Dv - 3) * (Dv - 4))   # 0911.4257, leading order in lambda
        res.append((f"GB D={Dv}: CEMZ rho->1 t2 / Buchel et al. t2 (expect 1)", sp.simplify(t2 / buchel)))
    return res


def check_QTG_ECG():
    mu = sp.symbols('mu')
    res = []
    # flat on-shell 3-point: I2 = I1/2 (three_point_vertices.py); L=(1/16piG)[... -(alpha3/6) I1]; alpha4_CEMZ = -2 alpha3
    # QTG (1004.2055 (action)): + (7/8) mu L^4 Z'_5, Z'_5 = I1 + Ricci terms  => -(alpha3/6) = 7/8 mu
    a3_qtg = -sp.Rational(6) * sp.Rational(7, 8) * mu
    t4_qtg = 5 * (25 - 1) * 3 * (-2 * a3_qtg)
    res.append(("QTG d=4: t4 from dictionary (MPS 1004.2055: 3780 mu)", t4_qtg))
    # Myers-Robinson form (7/4) mu L^4 Z_5, Z_5 = I2 + Ricci terms  => same since I2 = I1/2
    a3_qtg2 = -6 * sp.Rational(7, 4) * sp.Rational(1, 2) * mu
    res.append(("QTG d=4 via Z_5 = I2 + ... (7/4 mu): t4", 360 * (-2 * a3_qtg2)))
    # ECG (1802.00018 (ECG)): -(mu L^4/8) P,  P = 12 I2 + I1 + Ricci terms = 7 I1 on-shell  => -(alpha3/6) = -7mu/8
    a3_ecg = 6 * sp.Rational(7, 8) * mu
    t4_ecg = 4 * (16 - 1) * 2 * (-2 * a3_ecg)
    res.append(("ECG d=3: t4 from dictionary (BCR 1802.00018: -1260 mu)", t4_ecg))
    res.append(("d=3 dictionary: t4 = c * alpha3/L^4 with c", 4 * 15 * 2 * (-2)))
    return res


def random_weyl(dim, rng):
    """Random algebraic Weyl tensor (Euclidean metric): sum of Kulkarni-Nomizu products, then trace-free projection."""
    R = np.zeros((dim,) * 4)
    for _ in range(6):
        h = rng.normal(size=(dim, dim)); h = h + h.T
        k = rng.normal(size=(dim, dim)); k = k + k.T
        KN = (np.einsum('ac,bd->abcd', h, k) + np.einsum('bd,ac->abcd', h, k)
              - np.einsum('ad,bc->abcd', h, k) - np.einsum('bc,ad->abcd', h, k))
        R += KN
    g = np.eye(dim)
    Ric = np.einsum('abcb->ac', R)
    S = np.trace(Ric)
    P = (Ric - S / (2 * (dim - 1)) * g) / (dim - 2)    # Schouten
    KNgP = (np.einsum('ac,bd->abcd', g, P) + np.einsum('bd,ac->abcd', g, P)
            - np.einsum('ad,bc->abcd', g, P) - np.einsum('bc,ad->abcd', g, P))
    W = R - KNgP
    return W


def check_weyl_identity():
    rng = np.random.default_rng(7)
    out = []
    for dim in (4, 5, 6):
        rat = []
        for _ in range(3):
            W = random_weyl(dim, rng)
            assert np.max(np.abs(np.einsum('abcb->ac', W))) < 1e-9
            I1 = np.einsum('abcd,cdef,efab->', W, W, W)
            I2 = np.einsum('acbd,cedf,eafb->', W, W, W)
            rat.append(I2 / I1)
        out.append((f"dim={dim}: I2/I1 for random Weyl tensors (1/2 in dim<=5 expected; generic in 6)", np.round(rat, 10)))
    return out


def check_HM_d3():
    th = np.linspace(0, 2 * np.pi, 4001)[:-1]
    out = []
    for e in (np.array([[1, 0], [0, -1]]) / np.sqrt(2), np.array([[0, 1], [1, 0]]) / np.sqrt(2),
              np.array([[1, 1j], [1j, -1]]) / 2):
        n = np.stack([np.cos(th), np.sin(th)])
        enn = np.einsum('ij,in,jn->n', e, n, n)
        val = np.abs(enn) ** 2 / np.sum(np.abs(e) ** 2)
        out.append(("HM d=3: angular mean of |e_ij n_i n_j|^2/e*e (expect 1/4); min,max", (round(val.mean(), 12), round(val.min(), 6), round(val.max(), 6))))
    return out


def check_CT():
    d = 3
    L, G = sp.symbols('L G', positive=True)
    CT = sp.Rational(d + 1) * sp.gamma(d + 1) / ((d - 1) * sp.pi ** sp.Rational(d, 2) * sp.gamma(sp.Rational(d, 2))) * L ** (d - 1) / (8 * sp.pi * G)
    CT_scalar = sp.Rational(d, d - 1) / (2 * sp.pi ** sp.Rational(d, 2) / sp.gamma(sp.Rational(d, 2))) ** 2
    k, N = sp.symbols('k N', positive=True)
    # ABJM: F = pi L^2/(2G) = pi sqrt(2k) N^{3/2}/3  => L^2/G = 2 sqrt(2k) N^{3/2}/3
    cT_abjm = sp.simplify((CT / CT_scalar).subs(G, L ** 2 / (2 * sp.sqrt(2 * k) * N ** sp.Rational(3, 2) / 3)))
    return [("C_T(d=3)", sp.simplify(CT)), ("C_T^{free real scalar}(d=3)", CT_scalar),
            ("c_T = C_T/C_T^scalar", sp.simplify(CT / CT_scalar)),
            ("ABJM c_T (expect 64 sqrt(2k) N^{3/2}/(3 pi))", cT_abjm)]


def check_regulator():
    """Sharp IR cutoff (flat Theorem: L(eps) = int_eps^1 psi2/p^2) vs massive regulator 1/(p^2+a^2) (CHMRS AdS4 graviton
    propagator 1/(nu^2+(d/2)^2), i.e. a = 3/(2 M L) = 3/(2 Delta_gap)). Difference is O(a)."""
    import mpmath as mp
    psi2 = lambda p: (1 - p) ** 3 * p * (65 + 155 * p + 47 * p ** 2)
    out = []
    for a in (mp.mpf('0.01'), mp.mpf('0.001'), mp.mpf('0.0001')):
        sharp = mp.quad(lambda p: psi2(p) / p ** 2, [a, 1])
        massive = mp.quad(lambda p: psi2(p) / (p ** 2 + a ** 2), [0, a, 1])
        out.append((f"a={a}: int psi2/(p^2+a^2) - int_a^1 psi2/p^2  (O(a) expected)", mp.nstr(massive - sharp, 8)))
    return out


def translation():
    A, B = 37.8387, -45.2059          # FLAT_GRAVITON_BOUND.md Theorem 3 (i), m_IR <= M/100
    Bii = -40.790                      # Theorem 3 (ii), eikonal hypothesis, m_IR <= M/20000
    c = 240.0 ** 2
    lines = []
    lines.append(("(t4^2+alphaT^2) <= 57600 |ghat3|^2 M^8 / Delta^8 ;  57600*A", c * A))
    for name, Bv in (("(i)", B), ("(ii)", Bii)):
        # m_IR -> 3/(2L) heuristic (CHMRS (5.36)): log(M/m_IR) = log(2 Delta_gap/3)
        const = c * (Bv - A * np.log(1.5))
        lines.append((f"{name} constant term with log(M/m_IR) -> log(2 Delta/3): 57600*(B - A log 1.5)", const))
    lines.append(("|c|^2 = (t4^2+alphaT^2)/16 = 3600 |ghat3|^2 M^8/Delta^8; 3600*A", 3600 * A))
    lines.append(("CLPS (4.4) optimal (uncertified here): 57600*24.9, 57600*(-27.6)", (c * 24.9, c * -27.6)))
    # Theorem 3 (i) needs m_IR <= M/100, i.e. (with m_IR = 3/(2L), M = Delta_gap/L) Delta_gap >= 150;
    # Theorem 3 (ii) needs m_IR <= M/20000, i.e. Delta_gap >= 30000.  Smaller Delta_gap: illustration only.
    for Dg in (10, 30, 150, 1000, 30000):
        rhs = c * (A * np.log(2 * Dg / 3) + B) / Dg ** 8
        tag = "in range of Thm 3(i)" if Dg >= 150 else "OUTSIDE range of Thm 3 (illustration)"
        lines.append((f"Delta_gap={Dg}: RHS of t4^2+alphaT^2 bound (i), heuristic uplift [{tag}]; Delta^4*sqrt(RHS)", (rhs, Dg ** 4 * np.sqrt(max(rhs, 0)))))
    return lines


if __name__ == "__main__":
    for block in (check_varpi, check_limits, check_GB, check_QTG_ECG, check_weyl_identity, check_HM_d3, check_CT, check_regulator, translation):
        print(f"--- {block.__name__}")
        for name, val in block():
            print(f"{name}: {val}")

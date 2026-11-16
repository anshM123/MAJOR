"""
Numerical check of the flat-space part of the R^3 dictionary (DICTIONARY_T4.md, section 2).

What is checked (random complex 3-particle kinematics, explicit 4x4 tensors, no spinor shortcuts
in the curvature contractions):

  (a) Einstein structure  A1 = e1.e2 e3.p1 + e1.e3 e2.p3 + e2.e3 e1.p2
      |A1(1-,2-,3+)| = |<12>^3/(<23><31>)| / sqrt(2)              -> sqrt(32 pi G) A1^2 = sqrt(8 pi G)<12>^6/(...)
  (b) R^3 structure       A2 = e1.p2 e2.p3 e3.p1
      |A2(+,+,+)|^2 = |[12][23][31]|^2 / 8
  (c) linearised Riemann  R_i = F_i (x) F_i, F_i = p_i e_i - e_i p_i  (g = eta + kappa h, R = (kappa/2) F(x)F)
      I1 = R1^{ab}_{cd} R2^{cd}_{ef} R3^{ef}_{ab}            = 8 A2^2   (each of the 6 assignments)
      I2 = R1_a^c_b^d R2_c^e_d^f R3_e^a_f^b  (ECG/QTG index order) = 4 A2^2
      I2' = C_{mnrs} C^{n a s b} C_a^m_b^r (2205.01495 index order) = 4 A2^2
      => on-shell three-point: I2 = I2' = I1/2,  C^3 := 3 I1 - 4 I2' = 8 A2^2
  (d) duality: for helicity +, *F = s i F with s = +1 or -1 (convention); the parity-odd density
      R R R~ evaluated on (+,+,+) equals s i x I1  => all-plus vertex proportional to alpha3 + s i alpha3~.
  (e) mixed helicity: I1(+,+,-) = 0 (R^3 gives only all-plus / all-minus).
  (f) gauge invariance: everything independent of reference spinors.

Consequences printed at the end (combinatorics: 3! assignments, (kappa/2)^3, kappa^2 = 32 pi G):
  L = (1/16piG)(-alpha3/6) R^(3)  ->  M(+++) = -4 alpha3 sqrt(8piG) A2^2 = -2 alpha3 sqrt(32piG) A2^2
  so |M(+++)| = (|alpha3|/2) sqrt(8piG) |[12][13][23]|^2   (CLPS 2201.06602 (2.2): |ghat3| = |alpha3|)
  and alpha4^{CEMZ} = -2 alpha3  (CEMZ 1407.5597 eq. \totalr normalisation sqrt(32piG)(A1^2 + alpha4 A2^2)).
"""
import numpy as np

rng = np.random.default_rng(20261116)
eta = np.diag([1.0, -1.0, -1.0, -1.0]).astype(complex)        # mostly minus for the spinor map
sig = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
sig = [s.astype(complex) for s in sig]
eps2 = np.array([[0, 1], [-1, 0]], dtype=complex)


def ang(a, b):   # <ab> = eps^{alpha beta} a_alpha b_beta
    return a @ eps2 @ b


def sq(a, b):    # [ab]
    return a @ eps2 @ b


def vec(l, lt):
    """p^mu with p_{alpha alphadot} = l_alpha lt_alphadot = p^mu sigma_mu (lower index = eta p)."""
    P = np.outer(l, lt)
    # P = p^0 1 - p^i sigma_i  (sigma_mu = (1, sigma_i) with lowered index => p_mu sigma^mu)
    p0 = 0.5 * np.trace(P)
    pi = [-0.5 * np.trace(P @ sig[i]) for i in (1, 2, 3)]
    return np.array([p0, *pi])


def lower(v):
    return eta @ v


def dot(a, b):
    return a @ eta @ b


def pol_plus(l, lt, r):
    # eps+^mu = <r|sigma^mu|p] / (sqrt2 <r p>);  the vector <r|sigma^mu|p] has bispinor 2 r lt, i.e. = 2 vec(r, lt)
    return 2 * vec(r, lt) / (np.sqrt(2) * ang(r, l))


def pol_minus(l, lt, rt):
    # eps-^mu = <p|sigma^mu|r] / (sqrt2 [p r])
    return 2 * vec(l, rt) / (np.sqrt(2) * sq(lt, rt))


def rand_spinor():
    return rng.normal(size=2) + 1j * rng.normal(size=2)


def Fmat(p, e):
    """F^{mu nu} = p^mu e^nu - p^nu e^mu (upper indices)."""
    return np.outer(p, e) - np.outer(e, p)


def low2(F):
    return eta @ F @ eta


# Levi-Civita with eps^{0123} = +1 (upper); lower via eta
LC = np.zeros((4, 4, 4, 4))
import itertools
for perm in itertools.permutations(range(4)):
    LC[perm] = np.linalg.det(np.eye(4)[list(perm)])


def dual(F):
    """(*F)^{mu nu} = 1/2 eps^{mu nu}_{rho sigma} F^{rho sigma}, eps^{0123}=+1."""
    Fl = low2(F)
    return 0.5 * np.einsum('mnrs,rs->mn', LC, Fl) * 1.0  # eps^{mn rs} F_{rs}


def riem(F):
    """R^{ab}{}^{cd} (all upper) = F^{ab} F^{cd} (times kappa/2, dropped)."""
    return np.einsum('ab,cd->abcd', F, F)


def lower_all(R):
    return np.einsum('abcd,ai,bj,ck,dl->ijkl', R, eta, eta, eta, eta)


def Rc_mixed(R):
    # R^{ef}_{ab}: upper first pair, lower second pair
    return np.einsum('efab,ai,bj->efij', R, eta, eta)


def I1_(Ra, Rb, Rc):
    # sum_{indices} Ra^{ab}_{cd} Rb^{cd}_{ef} Rc^{ef}_{ab}
    return np.einsum('abcd,cdef,efab->', Rc_mixed(Ra), Rc_mixed(Rb), Rc_mixed(Rc))


def I2_(Ra, Rb, Rc):
    # Ra_a^c_b^d Rb_c^e_d^f Rc_e^a_f^b ; X_a^c_b^d = eta_{aa'} eta_{bb'} X^{a'c b'd}
    def m(R):
        return np.einsum('pcqd,pa,qb->acbd', R, eta, eta)   # index order (a,c,b,d): a,b lower; c,d upper
    A, B, C = m(Ra), m(Rb), m(Rc)
    return np.einsum('acbd,cedf,eafb->', A, B, C)


def I2p_(Ra, Rb, Rc):
    # C_{m n r s} C^{n a s b} C_a^m_b^r
    Al = lower_all(Ra)
    Bu = Rb
    Cm = np.einsum('pmqr,pa,qb->ambr', Rc, eta, eta)   # C_a^m_b^r : a,b lower, m,r upper
    return np.einsum('mnrs,nasb,ambr->', Al, Bu, Cm)


def three_point_plus(conj_ref=True):
    lam = rand_spinor()
    c = rng.normal(size=3) + 1j * rng.normal(size=3)
    lt1, lt2 = rand_spinor(), rand_spinor()
    lt3 = -(c[0] * lt1 + c[1] * lt2) / c[2]
    ls = [c[0] * lam, c[1] * lam, c[2] * lam]
    lts = [lt1, lt2, lt3]
    ps = [vec(ls[i], lts[i]) for i in range(3)]
    r = rand_spinor()
    es = [pol_plus(ls[i], lts[i], r) for i in range(3)]
    return ls, lts, ps, es


def report():
    out = []
    # momentum conservation / on-shell sanity
    ls, lts, ps, es = three_point_plus()
    tot = sum(ps)
    out.append(("momentum conservation |sum p|", np.max(np.abs(tot))))
    out.append(("max |p_i.p_j|", max(abs(dot(ps[i], ps[j])) for i in range(3) for j in range(3))))
    lA, ltA = rand_spinor(), rand_spinor()
    out.append(("normalisation |e+(p).e-(p)| (expect modulus 1; sign is convention)", dot(pol_plus(lA, ltA, rand_spinor()), pol_minus(lA, ltA, rand_spinor()))))
    out.append(("max |e_i.p_i|, |e_i.e_i|", max(max(abs(dot(es[i], ps[i])), abs(dot(es[i], es[i]))) for i in range(3))))

    A2 = dot(es[0], ps[1]) * dot(es[1], ps[2]) * dot(es[2], ps[0])
    s123 = sq(lts[0], lts[1]) * sq(lts[1], lts[2]) * sq(lts[2], lts[0])
    out.append(("(b) |A2(+++)|^2 / |[12][23][31]|^2   (expect 1/8)", abs(A2) ** 2 / abs(s123) ** 2))

    Fs = [Fmat(ps[i], es[i]) for i in range(3)]
    Rs = [riem(F) for F in Fs]
    i1 = I1_(*Rs)
    i2 = I2_(*Rs)
    i2p = I2p_(*Rs)
    out.append(("(c) I1/A2^2   (expect 8)", i1 / A2 ** 2))
    out.append(("(c) I2/A2^2   (ECG order, expect 4)", i2 / A2 ** 2))
    out.append(("(c) I2'/A2^2  (2205.01495 order, expect 4)", i2p / A2 ** 2))
    # all six assignments equal
    import itertools as it
    vals = [I1_(*[Rs[k] for k in perm]) / i1 for perm in it.permutations(range(3))]
    vals2 = [I2_(*[Rs[k] for k in perm]) / i2 for perm in it.permutations(range(3))]
    out.append(("(c) I1 over the 6 assignments / I1 (expect all 1)", np.round(vals, 12)))
    out.append(("(c) I2 over the 6 assignments / I2 (expect all 1)", np.round(vals2, 12)))
    # tr(F1F2F3) = 2 A2
    trF = np.einsum('ab,bc,ca->', Fs[0] @ eta, Fs[1] @ eta, Fs[2] @ eta) if False else np.trace(Fs[0] @ eta @ Fs[1] @ eta @ Fs[2] @ eta)
    out.append(("(c) tr(F1F2F3)/A2 (expect +-2)", trF / A2))

    # (d) duality
    dF = [dual(F) for F in Fs]
    ratios = [dF[i][np.unravel_index(np.argmax(abs(Fs[i])), (4, 4))] / Fs[i][np.unravel_index(np.argmax(abs(Fs[i])), (4, 4))] for i in range(3)]
    resid = [np.max(np.abs(dF[i] - ratios[i] * Fs[i])) for i in range(3)]
    out.append(("(d) *F_i / F_i for helicity + (expect common +-i)", np.round(ratios, 12)))
    out.append(("(d) max residual |*F - ratio F|", max(resid)))
    # parity-odd density R R R~ with the dual on the first pair of the third Riemann
    Rt3 = np.einsum('ab,cd->abcd', dF[2], Fs[2])
    i1t = I1_(Rs[0], Rs[1], Rt3)
    out.append(("(d) [R R R~](+++)/I1 (expect +-i)", i1t / i1))

    # (e) mixed helicity: (+,+,-) kinematics need lambda-tilde proportional for nonzero, but I1 vanishes for any
    # (++-) kinematics; use general complex kinematics with lt proportional (MHV-bar) and also lambda proportional
    lt = rand_spinor(); c = rng.normal(size=3) + 1j * rng.normal(size=3)
    l1, l2 = rand_spinor(), rand_spinor(); l3 = -(c[0] * l1 + c[1] * l2) / c[2]
    ls2 = [l1, l2, l3]; lts2 = [c[0] * lt, c[1] * lt, c[2] * lt]
    ps2 = [vec(ls2[i], lts2[i]) for i in range(3)]
    r, rt = rand_spinor(), rand_spinor()
    es2 = [pol_minus(ls2[0], lts2[0], rt), pol_minus(ls2[1], lts2[1], rt), pol_plus(ls2[2], lts2[2], r)]
    A1 = dot(es2[0], es2[1]) * dot(es2[2], ps2[0]) + dot(es2[0], es2[2]) * dot(es2[1], ps2[2]) + dot(es2[1], es2[2]) * dot(es2[0], ps2[1])
    mhv = ang(ls2[0], ls2[1]) ** 3 / (ang(ls2[1], ls2[2]) * ang(ls2[2], ls2[0]))
    out.append(("(a) |A1(--+)| * sqrt2 / |<12>^3/(<23><31>)|  (expect 1)", abs(A1) * np.sqrt(2) / abs(mhv)))
    Fs2 = [Fmat(ps2[i], es2[i]) for i in range(3)]
    Rs2 = [riem(F) for F in Fs2]
    out.append(("(e) |I1(--+)| (expect 0)", abs(I1_(*Rs2))))
    # (f) gauge invariance: change reference spinor
    ls, lts, ps, es = ls, lts, ps, es
    r2 = rand_spinor()
    es_b = [pol_plus(ls[i], lts[i], r2) for i in range(3)]
    A2b = dot(es_b[0], ps[1]) * dot(es_b[1], ps[2]) * dot(es_b[2], ps[0])
    out.append(("(f) A2 with new reference / A2 (expect 1)", A2b / A2))
    return out


if __name__ == "__main__":
    for name, val in report():
        print(f"{name}: {val}")
    print()
    print("Consequences (exact rationals, from (b),(c)):")
    print("  L=(1/16piG)(-a3/6) I1, 3! assignments, R=(k/2)FxF, k^2=32piG:")
    print("  M(+++) = -(a3/6)(1/16piG) * 6 * (k/2)^3 * 8 A2^2 = -4 a3 sqrt(8piG) A2^2 = -2 a3 sqrt(32piG) A2^2")
    print("  |M(+++)| = 4|a3| sqrt(8piG) |[12][23][31]|^2/8 = (|a3|/2) sqrt(8piG) |[12][13][23]|^2  => |ghat3| = |a3| (CLPS (2.2),(2.11))")
    print("  CEMZ sqrt(32piG)(A1^2 + alpha4 A2^2)  =>  alpha4_CEMZ = -2 a3 ;  2205.01495 (a4/12)C^3 -> +a4 A2^2 consistent (C^3 = 8 A2^2)")

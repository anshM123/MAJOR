"""
Flat-space / algebraic checks for the quartic (R^4) part of the eta/s map.  [computed]

(1) In D = 4 the Gross-Witten/Gubser-Klebanov-Tseytlin quartic Weyl density
        W_GW = C^{hmnk} C_{pmnq} C_h^{rsp} C^q_{rsk} + 1/2 C^{hkmn} C_{pqmn} C_h^{rsp} C^q_{rsk}
    is a combination of O1 = (C_{abcd}C^{abcd})^2 and O2 = (C_{abcd} C~^{abcd})^2 (C~ = 1/2 eps C).
    We fit the two coefficients on random Lorentzian Weyl tensors.

(2) The 4-graviton contact amplitude of L = (1/16 pi G) gamma W_GW, in the conventions of
    ../dictionary/three_point_vertices.py (g = eta + kappa h, R = (kappa/2) F(x)F, eps^{+-} with 1/sqrt2,
    kappa^2 = 32 pi G; these reproduce CLPS's Einstein three-point vertex):
        M_contact = (1/16 pi G) gamma (kappa/2)^4 sum_{sigma in S4} W_GW(R_sigma1..R_sigma4)
    and compare with g4 <12>^4 [34]^4 (CLPS/CHLPS f_low = 8 pi G/(stu) + ... + g4 + ..., M(1-2-3+4+) =
    <12>^4[34]^4 f). Also compute the all-plus contact amplitude of O1 - O2 type operators.
"""
import itertools
import numpy as np

rng = np.random.default_rng(1)

# ---------------------------------------------------------------- (1) algebraic identity, Lorentzian
eta = np.diag([-1.0, 1.0, 1.0, 1.0])
LC = np.zeros((4, 4, 4, 4))
for perm in itertools.permutations(range(4)):
    LC[perm] = np.linalg.det(np.eye(4)[list(perm)])  # eps_{0123} = +1 (lower)


def raise_all(T):
    return np.einsum('abcd,ai,bj,ck,dl->ijkl', T, eta, eta, eta, eta)  # eta^-1 = eta


def rand_weyl():
    R = np.zeros((4, 4, 4, 4))
    for _ in range(6):
        A = rng.normal(size=(4, 4)); A = A - A.T
        B = rng.normal(size=(4, 4)); B = B - B.T
        R += np.einsum('ab,cd->abcd', A, B) + np.einsum('ab,cd->abcd', B, A)
    # remove totally antisymmetric part (first Bianchi)
    Rasym = np.zeros_like(R)
    for perm in itertools.permutations(range(4)):
        s = np.linalg.det(np.eye(4)[list(perm)])
        Rasym += s * np.transpose(R, perm)
    R = R - Rasym / 24
    # traces with metric g = eta (lower indices)
    g = eta
    Ric = np.einsum('abcd,ac->bd', R, np.linalg.inv(g))
    Rs = np.einsum('bd,bd->', Ric, np.linalg.inv(g))
    W = (R - 0.5 * (np.einsum('ac,bd->abcd', g, Ric) - np.einsum('ad,bc->abcd', g, Ric)
                    - np.einsum('bc,ad->abcd', g, Ric) + np.einsum('bd,ac->abcd', g, Ric))
         + Rs / 6 * (np.einsum('ac,bd->abcd', g, g) - np.einsum('ad,bc->abcd', g, g)))
    return W


def invariants(W):
    Wu = raise_all(W)
    O1 = np.einsum('abcd,abcd->', W, Wu) ** 2
    # C~_{abcd} = 1/2 eps_{ab}^{ef} C_{efcd}; eps_{ab}^{ef} = eps_{abgh} eta^{ge} eta^{hf}
    epsmix = np.einsum('abgh,ge,hf->abef', LC, eta, eta)
    Wt = 0.5 * np.einsum('abef,efcd->abcd', epsmix, W)
    O2 = np.einsum('abcd,abcd->', Wt, Wu) ** 2
    # W_GW with index placement: C^{hmnk} C_{pmnq} C_h^{rsp} C^q_{rsk}
    def mix(W, up):  # up: tuple of booleans for which slots are raised
        T = W
        for slot, u in enumerate(up):
            if u:
                T = np.moveaxis(np.tensordot(eta, T, axes=([1], [slot])), 0, slot)
        return T
    C_uuuu = Wu
    C_llll = W
    C_lulu = mix(W, (False, True, True, True))   # C_h^{rsp}
    C_uluu = mix(W, (True, False, False, True))  # C^q_{rs k}: q up, r s lower, k lower -> (True,False,False,False)
    C_q_rsk = mix(W, (True, False, False, False))
    t1 = np.einsum('hmnk,pmnq,hrsp,qrsk->', C_uuuu, C_llll, C_lulu, C_q_rsk)
    t2 = np.einsum('hkmn,pqmn,hrsp,qrsk->', C_uuuu, C_llll, C_lulu, C_q_rsk)
    return O1, O2, t1 + 0.5 * t2


rows, rhs = [], []
for _ in range(6):
    W = rand_weyl()
    O1, O2, GW = invariants(W)
    rows.append([O1, O2]); rhs.append(GW)
coef, res, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
print("(1) W_GW = k1 (W^2)^2 + k2 (W W~)^2 in D=4 (Lorentzian):  k1, k2 =", coef,
      " max rel residual =", np.max(np.abs(np.array(rows) @ coef - rhs) / np.abs(rhs)))

# ---------------------------------------------------------------- (2) contact amplitudes, spinor helicity
etaM = np.diag([1.0, -1.0, -1.0, -1.0]).astype(complex)   # mostly minus, as in three_point_vertices.py
sig = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
sig = [s.astype(complex) for s in sig]
e2 = np.array([[0, 1], [-1, 0]], dtype=complex)
ang = lambda a, b: a @ e2 @ b
sq = ang


def vec(l, lt):
    P = np.outer(l, lt)
    return np.array([0.5 * np.trace(P)] + [-0.5 * np.trace(P @ sig[i]) for i in (1, 2, 3)])


def pol_plus(l, lt, r):
    return 2 * vec(r, lt) / (np.sqrt(2) * ang(r, l))


def pol_minus(l, lt, rt):
    return 2 * vec(l, rt) / (np.sqrt(2) * sq(lt, rt))


rs = lambda: rng.normal(size=2) + 1j * rng.normal(size=2)


def four_kin():
    ls = [rs() for _ in range(4)]
    lt1, lt2 = rs(), rs()
    # solve l3 lt3 + l4 lt4 = -(l1 lt1 + l2 lt2) for lt3, lt4
    Mtot = -(np.outer(ls[0], lt1) + np.outer(ls[1], lt2))
    Mat = np.array([ls[2], ls[3]]).T  # columns l3, l4
    X = np.linalg.solve(Mat, Mtot)    # rows lt3, lt4
    return ls, [lt1, lt2, X[0], X[1]]


etaL = etaM.real


def riem_upper(p, e):
    F = np.outer(p, e) - np.outer(e, p)
    return np.einsum('ab,cd->abcd', F, F)  # R^{abcd}


def lower_all(T):
    return np.einsum('abcd,ai,bj,ck,dl->ijkl', T, etaL, etaL, etaL, etaL)


LCu = np.zeros((4, 4, 4, 4))
for perm in itertools.permutations(range(4)):
    LCu[perm] = np.linalg.det(np.eye(4)[list(perm)])   # eps^{0123} = +1


def multilinear(Rs):
    """Fully symmetrised (sum over 4! assignments) values of O1, O2 and W_GW on four linearised Riemanns."""
    out = np.zeros(3, dtype=complex)
    Rl = [lower_all(R) for R in Rs]
    def mixed(R, up):
        # R given all-upper; lower the slots with up=False
        T = R
        for slot, u in enumerate(up):
            if not u:
                T = np.moveaxis(np.tensordot(etaL, T, axes=([1], [slot])), 0, slot)
        return T
    # dual on first pair: (R~)^{abcd} = 1/2 eps^{ab}_{ef} R^{efcd}
    epsmix = np.einsum('abgh,ge,hf->abef', LCu, etaL, etaL)
    for perm in itertools.permutations(range(4)):
        A, B, C, D = [Rs[k] for k in perm]
        Al, Bl, Cl, Dl = [Rl[k] for k in perm]
        o1 = np.einsum('abcd,abcd->', A, Bl) * np.einsum('abcd,abcd->', C, Dl)
        At = 0.5 * np.einsum('abef,efcd->abcd', epsmix, A)
        Ct = 0.5 * np.einsum('abef,efcd->abcd', epsmix, C)
        o2 = np.einsum('abcd,abcd->', At, Bl) * np.einsum('abcd,abcd->', Ct, Dl)
        Clulu = mixed(C, (False, True, True, True))
        Dq = mixed(D, (True, False, False, False))
        t1 = np.einsum('hmnk,pmnq,hrsp,qrsk->', A, Bl, Clulu, Dq)
        t2 = np.einsum('hkmn,pqmn,hrsp,qrsk->', A, Bl, Clulu, Dq)
        out += [o1, o2, t1 + 0.5 * t2]
    return out


for trial in range(2):
    ls, lts = four_kin()
    ps = [vec(ls[i], lts[i]) for i in range(4)]
    r, rt = rs(), rs()
    # MHV (1-,2-,3+,4+)
    es = [pol_minus(ls[0], lts[0], rt), pol_minus(ls[1], lts[1], rt),
          pol_plus(ls[2], lts[2], r), pol_plus(ls[3], lts[3], r)]
    vals = multilinear([riem_upper(ps[i], es[i]) for i in range(4)])
    ref = ang(ls[0], ls[1]) ** 4 * sq(lts[2], lts[3]) ** 4
    print(f"(2) trial {trial}: MHV  sum_sigma [O1, O2, W_GW] / (<12>^4[34]^4) =", np.round(vals / ref, 10))
    # all-plus
    es = [pol_plus(ls[i], lts[i], r) for i in range(4)]
    vals = multilinear([riem_upper(ps[i], es[i]) for i in range(4)])
    ref = (sq(lts[0], lts[1]) ** 2 * sq(lts[2], lts[3]) ** 2 + sq(lts[0], lts[2]) ** 2 * sq(lts[1], lts[3]) ** 2
           + sq(lts[0], lts[3]) ** 2 * sq(lts[1], lts[2]) ** 2) ** 2
    print(f"(2) trial {trial}: ++++ sum_sigma [O1, O2, W_GW] / ([12]^2[34]^2+[13]^2[24]^2+[14]^2[23]^2)^2 =",
          np.round(vals / ref, 10))

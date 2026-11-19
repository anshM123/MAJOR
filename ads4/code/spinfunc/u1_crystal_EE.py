"""
u1 (EE_FUNCTIONAL.md §1): the D = 4 scattering-crystal (bulk-point) image of the Regge-flip-free (2,4) structure
    EE := V0^4 (z2'.z4')^2   (s1_offdiag_kernel.STRUCTS['EE'], H7 §18.1 x0 = 0 reduction)
and of EH := V0^4 (z2'.z4') H24^.

Covariant form (derived in EE_FUNCTIONAL.md §1.1 from CLS-D's building blocks, the same definitions as g4/g5):
    at x0 = 0, x4 = -x2:  V_{2,04} V_{4,02} / X24 = -(z2.x2)(z4.x2)/x2^2,   H24/X24 = z2.z4 - 2 (z2.x2)(z4.x2)/x2^2
    =>  e24 := (H24 - 2 V_{2,04} V_{4,02}) / X24  reduces to z2.z4 = z2'.z4' (inversion referral), signature-independent;
    photon factor of EE:  E := V_{0,42}^2 X24 e24 / (X20 X40) = (lam/4) e24   (lam := 4 V042^2 X24/(X20 X40); x0 = 0: V0^2 e24)
    photon factor of q :  q := (2 H20 H40 + V_{0,42}^2 H24)/(X20 X40) = 2 Pc + (lam/4) Hc   (g5)
    graviton structures: q^2, EE = E^2, EH = E*(lam/4)Hc (all products, evaluated on null polarisations).
Checks: (0) the reduction e24 -> z2.z4 at a generic x0 = 0 configuration (random complex polarisations);
        (1) at the crystal (S0 point -1 of g5, where q = -Q/(4 nu^2 (m^2-nu^2)^2) exactly), decompose E, Hc, Pc, VVc on the flat bilinears
            {VV, WW, WV, VW} (CLS-D kernel vectors, g4 polarisation map);
        (2) helicity tables (g1/g2 helicity vectors, CLS-D conjugate pairing): flip content and opp/same channel ratio of EE, EH, and the
            minimal admixtures c1 q^2 + EE (+ c2 EH, + c3 PH-type) that are flip-free, with their preserved-channel shapes vs (1-x^2)^-4;
        (3) S0-point +1 consistency (x3 reflection).
Usage: python3 u1_crystal_EE.py -> cache/u1_crystal_EE.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

src = open(os.path.join(HERE, 'g4_images.py')).read().split("for x in (0.1")[0]
G = {'__file__': os.path.join(HERE, 'g4_images.py')}; exec(compile(src, 'g4h', 'exec'), G)
H, V, Xij, ldot, cart_to_lc, flat_vectors, EX, EZ, edot = (G[k] for k in ('H', 'V', 'Xij', 'ldot', 'cart_to_lc', 'flat_vectors', 'EX', 'EZ', 'edot'))
src1 = open(os.path.join(HERE, 'g1_crystal4d.py')).read().split('pr("# g1:')[0]
G1 = {'__file__': os.path.join(HERE, 'g1_crystal4d.py')}; exec(compile(src1, 'g1h', 'exec'), G1)
crystal, hel_vectors, dot, duality = G1['crystal'], G1['hel_vectors'], G1['dot'], G1['duality']

out = []
def pr(s=""):
    print(s, flush=True); out.append(str(s))

rng = np.random.default_rng(1)
pr("# u1: crystal image of EE = V0^4 (z2'.z4')^2 and EH; covariant e24 = (H24 - 2 V_{2,04} V_{4,02})/X24")
# (0) reduction check at x0 = 0, x4 = -x2, generic (complex) polarisations, Lorentzian lightcone metric of g4
pr("\n(0) x0 = 0 reduction: e24 vs z2.z4, and V_{0,42}^2 X24/(X20 X40) vs (z0.x2hat)^2-type factor")
for trial in range(3):
    x2 = rng.normal(size=3) + 0j
    x4 = -x2; x0 = np.zeros(3, complex)
    z2 = rng.normal(size=3) + 1j*rng.normal(size=3); z4 = rng.normal(size=3) + 1j*rng.normal(size=3)
    z0 = np.array([1.0, 0.0, 0.0], complex)
    e = (H(x2, z2, x4, z4) - 2*V(x2, z2, x0, x4)*V(x4, z4, x0, x2))/Xij(x2, x4)
    Hh = H(x2, z2, x4, z4)/Xij(x2, x4)
    v2v4 = V(x2, z2, x0, x4)*V(x4, z4, x0, x2)/Xij(x2, x4)
    pr(f"   trial {trial}: e24 = {e:.12f}, z2.z4 = {ldot(z2, z4):.12f}, |diff| = {abs(e - ldot(z2, z4)):.1e};  "
       f"V24V42/X24 + (z2.x2)(z4.x2)/x2^2 = {abs(v2v4 + ldot(z2, x2)*ldot(z4, x2)/ldot(x2, x2)):.1e}")


def cft_mats(x, sgn, m=1.0):
    nu, s, (V2, W2, V4, W4), fd, B2, B4 = flat_vectors(x, m)
    a = (m - s)/nu; Xc = np.array([-1j*a, 1j*a, 0.0]); x4 = Xc/(-ldot(Xc, Xc)); x2 = -x4
    z0 = np.array([1.0, 0.0, 0.0], dtype=complex)
    Z = [cart_to_lc(0, -1.0, 0.0), cart_to_lc(0, 0.0, 1.0)]
    x0 = np.array([0.0, 0.0, sgn*a], dtype=complex)
    v = V(x0, z0, x4, x2)
    lam = 4*v**2*Xij(x2, x4)/(Xij(x2, x0)*Xij(x4, x0))
    Pc = np.array([[H(x2, Z[i], x0, z0)*H(x4, Z[j], x0, z0)/(Xij(x2, x0)*Xij(x4, x0)) for j in range(2)] for i in range(2)])
    Hc = np.array([[H(x2, Z[i], x4, Z[j])/Xij(x2, x4) for j in range(2)] for i in range(2)])
    VVc = np.array([[V(x2, Z[i], x0, x4)*V(x4, Z[j], x0, x2)/Xij(x2, x4) for j in range(2)] for i in range(2)])
    flat = {}
    for nmA, A2 in (('V', V2), ('W', W2)):
        for nmB, A4 in (('V', V4), ('W', W4)):
            flat[nmA + nmB] = np.array([[fd(B2[i], A2)*fd(B4[j], A4) for j in range(2)] for i in range(2)])
    return dict(nu=nu, lam=lam, Pc=Pc, Hc=Hc, VVc=VVc, E=(lam/4)*(Hc - 2*VVc), q=2*Pc + (lam/4)*Hc, flat=flat,
                pref=4*nu**2*(m*m - nu*nu)**2)


def decompose(M, flat):
    B = np.stack([flat[k].flatten() for k in ('VV', 'WW', 'WV', 'VW')], 1)
    c, res, rk, sv = np.linalg.lstsq(B, M.flatten(), rcond=None)
    return c, np.linalg.norm(B @ c - M.flatten())/np.linalg.norm(M.flatten())


pr("\n(1) photon-level crystal images on the flat bilinears {VV, WW, WV, VW} (x 4 nu^2 (m^2-nu^2)^2), S0 point -1 (g5: q = -Q exactly)")
XS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
DEC = {}
for x in XS:
    C = cft_mats(x, -1)
    row = []
    for nm in ('q', 'Pc', 'Hc', 'VVc', 'E'):
        c, r = decompose(C[nm]*C['pref'], C['flat'])
        DEC[(x, nm)] = c
        row.append(f"{nm}: [" + ", ".join(f"{v.real:+.6f}{v.imag:+.1e}j" for v in c) + f"] res {r:.0e}")
    pr(f"  x={x}: lam = {C['lam']:.10f}")
    for s_ in row:
        pr("      " + s_)

# (2) helicity tables
pr("\n(2) helicity tables at the crystal (g1 helicity vectors; CLS-D conjugate pairing; 'same'/'opp' = preserved channels, flips = (h,-h))")


def hel_tables(x):
    s, sbar, Pm, (V2, W2, V4, W4) = crystal(x)
    p1, p2, p3, p4 = Pm
    H2 = hel_vectors(p2); H4 = hel_vectors(-p4)
    E2 = {h: H2[h] for h in (1, -1)}; E4 = {h: np.conj(H4[h]) for h in (1, -1)}
    lab1 = duality(p1, sbar); lab2 = {h: duality(p2, E2[h]) for h in (1, -1)}
    hsame = -1 if lab2[1] == lab1 else +1
    T = lambda A2, A4: {(a, b): dot(A2, E2[a])*dot(A4, E4[b]) for a in (1, -1) for b in (1, -1)}
    tabs = {'VV': T(V2, V4), 'WW': T(W2, W4), 'WV': T(W2, V4), 'VW': T(V2, W4)}
    pres = [(hsame, hsame), (-hsame, -hsame)]; flip = [(1, -1), (-1, 1)]
    return tabs, pres, flip


def photon_table(c, tabs):
    keys = ('VV', 'WW', 'WV', 'VW')
    return {k: sum(c[i]*tabs[keys[i]][k] for i in range(4)) for k in tabs['VV']}


SUM = {}
for x in XS:
    tabs, pres, flip = hel_tables(x)
    ph = {nm: photon_table(DEC[(x, nm)], tabs) for nm in ('q', 'Pc', 'Hc', 'VVc', 'E')}
    # graviton structures (products of photon factors on the same helicity labels)
    lamq = ph['q']; E = ph['E']; Hl = {k: ph['Hc'][k]*cft_mats(x, -1)['lam']/4 for k in ph['Hc']}
    gr = {'q2': {k: lamq[k]**2 for k in lamq}, 'EE': {k: E[k]**2 for k in E}, 'EH': {k: E[k]*Hl[k] for k in E}}
    SUM[x] = (ph, gr, pres, flip)
    pr(f"  x={x}:")
    for nm in ('q', 'E'):
        T = ph[nm]
        pr(f"     photon {nm:2s}: same {T[pres[0]]:+.6e}  opp {T[pres[1]]:+.6e}  flips {T[flip[0]]:+.6e} {T[flip[1]]:+.6e}")
    for nm in ('q2', 'EE', 'EH'):
        T = gr[nm]
        fl = max(abs(T[k]) for k in flip); pmax = max(abs(T[k]) for k in pres)
        pr(f"     graviton {nm:3s}: same {T[pres[0]]:+.6e}  opp {T[pres[1]]:+.6e}  |flip|/|pres|max {fl/pmax:.3e};  "
           f"opp/same {(T[pres[1]]/T[pres[0]]).real:+.8f}  ((1-x^2)^-4 = {(1 - x*x)**-4:.8f});  same/same_q2 {(T[pres[0]]/gr['q2'][pres[0]]).real:+.8f}")

# (3) bulk-point content of the Regge-flip-free structures, relative to q^2 (scale-free)
pr("\n(3) bulk-point content relative to q^2: max over channels |S - q2|/|q2_pres| (EE), |EH|/|q2_pres|, |EE flips|/|q2_pres|")
for x in XS:
    ph, gr, pres, flip = SUM[x]
    qs = max(abs(gr['q2'][k]) for k in pres)
    dEE = max(abs(gr['EE'][k] - gr['q2'][k]) for k in gr['q2'])/qs
    dEH = max(abs(gr['EH'][k]) for k in gr['q2'])/qs
    fEE = max(abs(gr['EE'][k]) for k in flip)/qs
    pr(f"   x={x}: |EE - q2|/|q2| = {dEE:.1e};  |EH|/|q2| = {dEH:.1e} (null at the bulk point);  EE flips/|q2| = {fEE:.1e}")

# (5) symbolic: at the crystal V_{2,04} is proportional to H20 (and V_{4,02} to H40) on the whole polarisation plane, with
#     V_{2,04} V_{4,02} X20 X40/(X24 H20 H40) = -4 and lam = 1 for all x  =>  E = (lam/4)(Hc - 2 VVc) = Hc/4 + 2 Pc = q  identically.
pr("\n(5) symbolic (sympy, x symbolic, polarisation z2 = z4 = (0, -u, w) arbitrary in the CLS-D plane):")
try:
    import sympy as sp
    xs = sp.symbols('x', positive=True); uu, ww = sp.symbols('u w')
    sq = sp.sqrt(1 - xs**2); aa = (1 - sq)/xs
    sl = lambda A, B: -(A[0]*B[1] + A[1]*B[0])/2 + A[2]*B[2]
    sEX = lambda p: (1, sl(p, p), p); sEZ = lambda p, z: (0, 2*sl(p, z), z)
    sed = lambda A, B: -(A[0]*B[1] + A[1]*B[0])/2 + sl(A[2], B[2])
    sX = lambda a_, b_: -2*sed(sEX(a_), sEX(b_))
    def sH(xi, zi, xj, zj):
        Xi, Zi, Xj, Zj = sEX(xi), sEZ(xi, zi), sEX(xj), sEZ(xj, zj)
        return -2*(sed(Zi, Zj)*sed(Xi, Xj) - sed(Zi, Xj)*sed(Zj, Xi))
    def sV(xi, zi, xj, xk):
        Xi, Zi, Xj, Xk = sEX(xi), sEZ(xi, zi), sEX(xj), sEX(xk)
        return (sed(Zi, Xj)*sed(Xi, Xk) - sed(Zi, Xk)*sed(Xi, Xj))/sed(Xj, Xk)
    Xc = [-sp.I*aa, sp.I*aa, 0]; nr = -sl(Xc, Xc); x4s = [c/nr for c in Xc]; x2s = [-c for c in x4s]
    z0s = [1, 0, 0]; x0s = [0, 0, -aa]; zs = [uu, -uu, ww]          # cart_to_lc(0, -u, w) = (u, -u, w)
    lam_s = sp.simplify(4*sV(x0s, z0s, x4s, x2s)**2*sX(x2s, x4s)/(sX(x2s, x0s)*sX(x4s, x0s)))
    r2 = sp.simplify(sV(x2s, zs, x0s, x4s)/sH(x2s, zs, x0s, z0s)); r4 = sp.simplify(sV(x4s, zs, x0s, x2s)/sH(x4s, zs, x0s, z0s))
    comb = sp.simplify(r2*r4*sX(x2s, x0s)*sX(x4s, x0s)/sX(x2s, x4s))
    pr(f"   lam = {lam_s};  V_(2,04)/H20 depends on (u, w): {bool(sp.simplify(sp.diff(r2, uu)) != 0 or sp.simplify(sp.diff(r2, ww)) != 0)};"
       f"  V_(4,02)/H40 depends on (u, w): {bool(sp.simplify(sp.diff(r4, uu)) != 0 or sp.simplify(sp.diff(r4, ww)) != 0)};"
       f"  V_(2,04)V_(4,02) X20 X40/(X24 H20 H40) = {comb}")
except Exception as exc:                                          # pragma: no cover
    pr(f"   sympy check failed: {exc}")

# (4) S0 point +1: off-diagonal (n_perp x x3) entries odd, diagonal even
pr("\n(4) S0-point +1 vs -1 (x3 reflection: off-diagonal entries flip sign, diagonal unchanged) for E:")
for x in (0.1, 0.5, 0.9):
    Em, Ep = cft_mats(x, -1)['E'], cft_mats(x, +1)['E']
    R = np.diag([1.0, -1.0])
    pr(f"   x={x}: |E(+1) - R E(-1) R|/|E| = {np.linalg.norm(Ep - R @ Em @ R)/np.linalg.norm(Em):.1e}")

open(os.path.join(HERE, 'cache', 'u1_crystal_EE.txt'), 'w').write("\n".join(out) + "\n")

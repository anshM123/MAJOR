"""
x1 (EE_FOUNDATIONS.md §1): is EE = [V_{0,42}^2 (H24 - 2 V_{2,04} V_{4,02})/(X20 X40)]^2 an admissible j' = 0 dual (2,4) structure of the
CLS-D spinning light-ray sum rule (arXiv:2311.04271 eq. spinningfunctional_definition), on the same footing as q^2?

Embedding-space building blocks exactly as CLS-D eq. CFT_building_blocks_VH (and g4_images.py): X = (1, x^2, x), Z = (0, 2 x.z, z),
X.Y = -(X+ Y- + X- Y+)/2 + x.y,  X_ij = -2 X_i.X_j,  H_ij = -2 (Z_i.Z_j X_i.X_j - Z_i.X_j Z_j.X_i),
V_{i,jk} = (Z_i.X_j X_i.X_k - Z_i.X_k X_i.X_j)/(X_j.X_k).   q = (2 H20 H40 + V0^2 H24)/(X20 X40),  E = V0^2 (H24 - 2 V2 V4)/(X20 X40),
V0 = V_{0,42}, V2 = V_{2,04}, V4 = V_{4,02}.

Checks
 (A) weights: homogeneity degrees of q^2 and EE in each X_i, Z_i (generic embedding vectors, d = 3 and d = 4) and transversality
     Z_i -> Z_i + alpha X_i.  [exact algebra; numerics to round-off]
 (B) span: the parity-even monomials V2^a V4^b V0^c H24^n1 H20^n2 H40^n3 (a+n1+n2 = 2, b+n1+n3 = 2, c+n2+n3 = J; 14 of them) have
     numerical rank 13 in d = 3 (one 3d Gram relation) and 14 in d = 4; q^2 lies in CLS-D's printed j' = 0 span {H20^{2-n} H40^{2-n} H24^n}
     (rank test), EE does not; EE lies in the 13-dim span (trivially: 3 monomials).
 (C) the saddle: CLS-D eq. kernel_identity_saddle, -V_{2,40} V_{0,24} = H20 at the scattering crystal with z_i = (0, z^x, z_perp)
     (time component zero), and the same for 4; hence EE and q^2 coincide on the saddle quotient (u1: EE = q^2 at the crystal).
 (D) singular loci: q, E as x0 -> light cone of x2 (Lorentzian): q ~ X02^-1, E ~ X02^-2 (same locus, higher order); and the geometric
     form e24 := (H24 - 2 V2 V4)/X24 = z2.I(x20).I(x04).z4 (parallel transport 4 -> 0 -> 2 by inversion tensors).
 (E) leading Regge (2,4) side: EE/V0^4 at x0 = 0 is (z2'.z4')^2, independent of x-hat, z0 and of the polarisation referral; in the
     v4_hclosed/u2 exact plane-wave transform the (2,4)-side coefficients (c1, c2) of EE are a CONSTANT multiple of the scalar's, so
     h_-^EE/h_-^sc = c(nu + i) is entirely a (1,3)-side (light-ray) ratio:  c(nu + i) = kappa d1^T/d1^sc.
Usage: python3 x1_admissibility.py  -> cache/x1_admissibility.txt   (~1-2 min, one core)
"""
import os, sys, io, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np

buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


rng = np.random.default_rng(7)


# ---------------------------------------------------------------- embedding space, metric diag(-1 or +1 ..., ) on physical part
class Emb:
    def __init__(self, sig):            # sig: physical metric signature, e.g. (-1, 1, 1) Lorentzian 3d, (1, 1, 1) Euclidean 3d
        self.g = np.array(sig, float); self.n = len(sig)

    def pd(self, a, b):
        return np.sum(self.g*a*b)

    def X(self, x):
        return (1.0 + 0j, self.pd(x, x) + 0j, np.asarray(x, complex))

    def Z(self, x, z):
        return (0.0j, 2*self.pd(x, z) + 0j, np.asarray(z, complex))

    def ed(self, A, B):
        return -(A[0]*B[1] + A[1]*B[0])/2 + self.pd(A[2], B[2])

    # building blocks on embedding vectors
    def Xij(self, Xi, Xj):
        return -2*self.ed(Xi, Xj)

    def H(self, Xi, Zi, Xj, Zj):
        return -2*(self.ed(Zi, Zj)*self.ed(Xi, Xj) - self.ed(Zi, Xj)*self.ed(Zj, Xi))

    def V(self, Xi, Zi, Xj, Xk):
        return (self.ed(Zi, Xj)*self.ed(Xi, Xk) - self.ed(Zi, Xk)*self.ed(Xi, Xj))/self.ed(Xj, Xk)

    def blocks(self, P):
        X0, Z0, X2, Z2, X4, Z4 = P
        V0 = self.V(X0, Z0, X4, X2); V2 = self.V(X2, Z2, X0, X4); V4 = self.V(X4, Z4, X0, X2)
        H24 = self.H(X2, Z2, X4, Z4); H20 = self.H(X2, Z2, X0, Z0); H40 = self.H(X4, Z4, X0, Z0)
        X20 = self.Xij(X2, X0); X40 = self.Xij(X4, X0); X24 = self.Xij(X2, X4)
        return dict(V0=V0, V2=V2, V4=V4, H24=H24, H20=H20, H40=H40, X20=X20, X40=X40, X24=X24)

    def q(self, P):
        b = self.blocks(P)
        return (2*b['H20']*b['H40'] + b['V0']**2*b['H24'])/(b['X20']*b['X40'])

    def E(self, P):
        b = self.blocks(P)
        return b['V0']**2*(b['H24'] - 2*b['V2']*b['V4'])/(b['X20']*b['X40'])

    def null(self):
        """random complex null vector of the physical metric."""
        while True:
            v = rng.normal(size=self.n) + 1j*rng.normal(size=self.n)
            # solve last component for v.v = 0
            rest = np.sum(self.g[:-1]*v[:-1]**2)
            v[-1] = np.sqrt(-rest/self.g[-1] + 0j)
            if abs(self.pd(v, v)) < 1e-10:
                return v

    def config(self):
        xs = [rng.normal(size=self.n) for _ in range(3)]
        zs = [self.null() for _ in range(3)]
        return xs, zs

    def P(self, xs, zs):
        (x0, x2, x4), (z0, z2, z4) = xs, zs
        return (self.X(x0), self.Z(x0, z0), self.X(x2), self.Z(x2, z2), self.X(x4), self.Z(x4, z4))


def scaleP(P, idx, lam):
    P = list(P); P[idx] = tuple(lam*c for c in P[idx]); return tuple(P)


def shiftZ(P, iz, ix, al):
    P = list(P); P[iz] = tuple(a + al*b for a, b in zip(P[iz], P[ix])); return tuple(P)


def main():
    t00 = time.time()
    pr("# x1: admissibility of EE as a CLS-D j' = 0 dual (2,4) structure (EE_FOUNDATIONS.md §1)")
    # ------------------------------------------------------------- (A) weights and transversality
    pr("\n(A) homogeneity degrees (log|S(lam .)/S(.)|/log lam, lam = 1.7) and transversality |S(Z+aX) - S|/|S|")
    lab = ['X0', 'Z0', 'X2', 'Z2', 'X4', 'Z4']
    for sig, nm in (((-1, 1, 1), 'd=3 Lorentzian'), ((1, 1, 1, 1), 'd=4 Euclidean')):
        M = Emb(sig)
        xs, zs = M.config(); P = M.P(xs, zs)
        for snm, S in (('q^2', lambda P: M.q(P)**2), ('EE', lambda P: M.E(P)**2), ('q^2 - EE', lambda P: M.q(P)**2 - M.E(P)**2)):
            s0 = S(P)
            degs = [np.log(abs(S(scaleP(P, k, 1.7))/s0))/np.log(1.7) for k in range(6)]
            tr = [abs(S(shiftZ(P, iz, iz - 1, 0.37 - 0.2j)) - s0)/abs(s0) for iz in (1, 3, 5)]
            pr(f"   {nm:15s} {snm:9s}: degrees " + " ".join(f"{l}:{d:+.6f}" for l, d in zip(lab, degs)) +
               "; transversality Z0, Z2, Z4: " + " ".join(f"{t:.1e}" for t in tr))
    pr("   -> q^2 and EE have identical weights (embedding degrees -Delta in X, J in Z) and are transverse: same dual-structure type [P].")

    # ------------------------------------------------------------- (B) span / rank
    pr("\n(B) parity-even monomial space V2^a V4^b V0^c H24^n1 H20^n2 H40^n3 (common prefactor), J = 4 integer stand-in for J_F")
    J = 4
    monos = []
    for n1 in range(3):
        for n2 in range(3 - n1):
            for n3 in range(3 - n1):
                a = 2 - n1 - n2; b = 2 - n1 - n3; c = J - n2 - n3
                monos.append((a, b, c, n1, n2, n3))
    pr(f"   number of monomials: {len(monos)}")
    def mono_vals(M, P):
        bl = M.blocks(P)
        return np.array([bl['V2']**a*bl['V4']**b*bl['V0']**c*bl['H24']**n1*bl['H20']**n2*bl['H40']**n3 for (a, b, c, n1, n2, n3) in monos])
    printed = [(0, 0, J - 4, 0, 2, 2), (0, 0, J - 2, 1, 1, 1), (0, 0, J, 2, 0, 0)]       # V0^{J-4}H20^2H40^2, V0^{J-2}H20H40H24, V0^J H24^2
    ip = [monos.index(m) for m in printed]
    for sig, nm in (((1, 1, 1), 'd = 3 (Euclidean)'), ((-1, 1, 1), 'd = 3 (Lorentzian)'), ((1, 1, 1, 1), 'd = 4 (Euclidean)')):
        M = Emb(sig)
        rows, rq, rE = [], [], []
        for _ in range(60):
            xs, zs = M.config(); P = M.P(xs, zs)
            v = mono_vals(M, P)
            nrm = np.max(np.abs(v)); rows.append(v/nrm)
            bl = M.blocks(P); pref = (bl['X20']*bl['X40'])**2/bl['V0']**(J - 4)
            rq.append(M.q(P)**2*pref/nrm); rE.append(M.E(P)**2*pref/nrm)
        A = np.array(rows); s = np.linalg.svd(A, compute_uv=False)
        rank = int(np.sum(s > 1e-10*s[0]))
        pr(f"   {nm}: singular values (smallest 3): " + " ".join(f"{v:.1e}" for v in s[-3:]) + f" -> rank {rank}")
        Pm = A[:, ip]
        for snm, col in (('q^2', np.array(rq)), ('EE', np.array(rE))):
            c, *_ = np.linalg.lstsq(Pm, col, rcond=None); rp = np.linalg.norm(Pm @ c - col)/np.linalg.norm(col)
            c2, *_ = np.linalg.lstsq(A, col, rcond=None); ra = np.linalg.norm(A @ c2 - col)/np.linalg.norm(col)
            pr(f"      {snm}: residual of projection on CLS-D printed span {rp:.1e}" + (f" (coefficients {np.round(c.real, 6)})" if rp < 1e-8 else "") +
               f"; on the full monomial span {ra:.1e}")
        if nm.startswith('d = 3 (E'):
            # the 3d relation among the 14 monomials
            _, _, Vh = np.linalg.svd(A); rel = Vh[-1].conj(); rel = rel/rel[np.argmax(np.abs(rel))]
            pr("      3d Gram relation (null vector of the 14 monomials): " +
               ", ".join(f"{m}:{r.real:+.4f}" for m, r in zip(monos, rel) if abs(r) > 1e-8))
    pr("   EE = V0^J (H24^2 - 4 H24 V2 V4 + 4 V2^2 V4^2)/(X20 X40)^2 x V0^{-(J-4)}-normalised; q^2 = (4 H20^2H40^2 V0^{J-4} + 4 H20H40H24 V0^{J-2}"
       " + H24^2 V0^J)/(...)^2  [exact expansion]")

    # ------------------------------------------------------------- (C) the saddle identity (CLS-D kernel_identity_saddle)
    pr("\n(C) CLS-D eq. kernel_identity_saddle at the scattering crystal (g4 conventions: lightcone (u, v, y), x4 = -x2 = x/(-x^2),"
       " x = (-ia, ia, 0), x0 = (0, 0, +-a), z0 = (1, 0, 0)); polarisations with zero time component z = (0, z^x, z_y)")
    src = open(os.path.join(HERE, 'g4_images.py')).read().split("for x in (0.1")[0]
    G = {'__file__': os.path.join(HERE, 'g4_images.py')}; exec(compile(src, 'g4h', 'exec'), G)
    Hg, Vg, Xg, ldot, c2l = G['H'], G['V'], G['Xij'], G['ldot'], G['cart_to_lc']
    for xv in (0.2, 0.5, 0.8):
        m = 1.0; nu = xv*m; s = np.sqrt(m*m - nu*nu); a = (m - s)/nu
        Xc = np.array([-1j*a, 1j*a, 0.0]); x4 = Xc/(-ldot(Xc, Xc)); x2 = -x4
        z0 = np.array([1.0, 0.0, 0.0], complex)
        out = []
        for sgn in (+1, -1):
            x0 = np.array([0.0, 0.0, sgn*a], complex)
            for (zx, zy) in ((1.0, 0.0), (0.0, 1.0), (0.3, -0.7)):
                z2 = c2l(0, zx, zy); z4 = c2l(0, -0.4*zx + 0.2, zy + 0.5)
                lhs2 = -Vg(x2, z2, x4, x0)*Vg(x0, z0, x2, x4); rhs2 = Hg(x2, z2, x0, z0)
                lhs4 = -Vg(x4, z4, x2, x0)*Vg(x0, z0, x4, x2); rhs4 = Hg(x4, z4, x0, z0)
                out.append(max(abs(lhs2 - rhs2)/max(abs(rhs2), 1e-300), abs(lhs4 - rhs4)/max(abs(rhs4), 1e-300)))
        pr(f"   x = {xv}: max rel |(-V_{{2,40}} V_{{0,24}}) - H20|, same for 4 (both S0 points, 3 polarisations): {max(out):.1e}")
    pr("   -> at the saddle V_{2,04}, V_{4,02} are proportional to H20, H40 (CLS-D: 'allows us to just remove all V_{2,40}, V_{4,02} structures"
       " when writing down the basis'); the printed list is a basis of the SADDLE quotient, not of the space of covariant structures.")

    # ------------------------------------------------------------- (D) singular loci and the geometric form of e24
    pr("\n(D) singular loci of the photon factors as x0 -> light cone of x2 (Lorentzian d = 3): fitted power of |q|, |E| vs X02")
    M = Emb((-1, 1, 1))
    x2 = np.array([0.3, -0.2, 0.5]); x4 = np.array([-1.1, 0.4, -0.3])
    n = np.array([1.0, np.cos(0.7), np.sin(0.7)])          # null direction
    w = np.array([0.2, 0.9, -0.4])
    z0, z2, z4 = M.null(), M.null(), M.null()
    eps = np.geomspace(1e-3, 1e-6, 4)
    vq, vE, vX = [], [], []
    for e in eps:
        x0 = x2 + 0.8*n + e*w
        P = M.P([x0, x2, x4], [z0, z2, z4])
        vq.append(abs(M.q(P))); vE.append(abs(M.E(P))); vX.append(abs(M.blocks(P)['X20']))
    pq = np.polyfit(np.log(vX), np.log(vq), 1)[0]; pE = np.polyfit(np.log(vX), np.log(vE), 1)[0]
    pr(f"   |q| ~ X02^{pq:+.4f},  |E| ~ X02^{pE:+.4f}   (so q^2 ~ X02^-2, EE ~ X02^-4: same singular locus, order 2 higher)")
    # e24 geometric identity
    def Imat(g, x):
        return np.diag(g)*0 + np.eye(len(g)) - 2*np.outer(x, g*x)/np.sum(g*x*x)
    for sig, nm in (((1, 1, 1), 'd=3 Euclidean'), ((-1, 1, 1), 'd=3 Lorentzian'), ((1, 1, 1, 1), 'd=4 Euclidean')):
        M = Emb(sig); g = M.g
        errs = []
        for _ in range(5):
            xs, zs = M.config(); P = M.P(xs, zs); bl = M.blocks(P)
            x0, x2, x4 = xs; _, z2, z4 = zs
            e24 = (bl['H24'] - 2*bl['V2']*bl['V4'])/bl['X24']
            # z2 . I(x20) . I(x04) . z4 with lowered indices: (I v)^mu = v^mu - 2 x^mu (x.v)/x^2
            Iv = lambda x, v: v - 2*x*M.pd(x, v)/M.pd(x, x)
            geo = M.pd(Iv(x2 - x0, z2), Iv(x4 - x0, z4))
            Hh = bl['H24']/bl['X24']
            errs.append((abs(e24 - geo)/abs(geo), abs(e24 + geo)/abs(geo), abs(Hh - M.pd(z2, Iv(x2 - x4, z4)))/abs(Hh)))
        errs = np.array(errs)
        pr(f"   {nm}: max |e24 - (I(x02)z2).(I(x04)z4)|/|.| = {errs[:, 0].max():.1e} (with + sign: {errs[:, 1].min():.1e}); "
           f"control H24/X24 = z2.I(x24).z4: {errs[:, 2].max():.1e}")
    pr("   -> e24 = (I(x02) z2).(I(x04) z4): both polarisations parallel-transported (by inversion tensors) to the light-ray point x0 and"
       " contracted there; q instead contains (I(x02)z2.z0)(I(x04)z4.z0) and z2.I(x24).z4.")

    # ------------------------------------------------------------- (E) leading-Regge (2,4) side: A_EE = scalar x constant
    pr("\n(E) leading Regge, exact plane-wave transform (v4_hclosed/u2 method): (2,4)-side covector coefficients (c1, c2)")
    import sympy as sp
    import v4_hclosed as H4
    import u2_hclosed_struct as U2
    t0 = time.time()
    c1s, c2s = H4.side_coeffs(False, 'A'); c1e, c2e = U2.side_coeffs_A('EE')
    r1 = sp.simplify(c1e/c1s); r2 = sp.simplify(c2e/c2s)
    pr(f"   c1^EE/c1^sc = {r1};  c2^EE/c2^sc = {r2}   [{time.time() - t0:.0f}s]")
    d1T, d2T = H4.side_coeffs(True, 'B'); d1s, d2s = H4.side_coeffs(False, 'B')
    dd = H4.d
    sub = {dd: sp.Rational(1, 2) + sp.I*H4.nu_s}
    beta1 = sp.factor(sp.simplify((r1*d1T/d1s).subs(sub)))
    nu = H4.nu_s
    cfun = lambda v: 32*(4*v**2 + 1)*(4*v**2 + 9)/((4*v**2 + 25)*(4*v**2 + 81))
    pr(f"   kappa d1^T/d1^sc (nu) = {beta1}")
    pr(f"   kappa d1^T/d1^sc - c(nu + i) = {sp.simplify(beta1 - cfun(nu + sp.I))}  (c(nu) = 32(4nu^2+1)(4nu^2+9)/((4nu^2+25)(4nu^2+81)))")
    # and the R_+ ratio
    l = dd - 1
    Rp_e = r1*c1s*(d1T*(l + 1)/(2*l + 1) + d2T); Rp_s = c1s*(d1s*(l + 1)/(2*l + 1) + d2s)
    pr(f"   R_+^EE/R_+^sc - c(nu - i) = {sp.simplify(sp.simplify(Rp_e/Rp_s).subs(sub) - cfun(nu - sp.I))}")
    pr("   -> the (2,4) side of EE is the scalar's times a constant; c(nu) comes entirely from the light-ray (1,3) side, which is the SAME"
       " for every j' = 0 (2,4) structure [P given H7 §18.1].")
    pr(f"\n[{time.time() - t00:.0f}s]")
    open(os.path.join(HERE, 'cache', 'x1_admissibility.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()

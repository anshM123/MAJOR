"""
s2_lp (OFFDIAG_KERNEL.md §2): the route-A' 2x2 per-state test and NF02's LP family with the EXACT leading-Regge off-diagonal of C^_2'
(s1 tables, geometric = parallel-transported helicity basis) instead of r4's fitted universal g.
Per state (m/Lam, J) at Delta_g:  beta = arccosh(1 + (2J+1)^2/(2 m^2)),
   F_hh = stripped flat B2 + B3 (exact, r3.flat_basis) + 2 int psi2 (P_l - J0) + u sigma int psi3 (P_l - J0)      (r5's diagonal, = F_hyb)
   F_+- = int dp psi2(p) eo2(x) Re Y(Lam p, beta)            [OFF=psi: no de-shift on the off-diagonal]
        or Re int dnu/Lam X_2(nu) eo2(x) Y(nu, beta)         [OFF=chi: exact C^_2' smearing, X = v4_chi solution for the weight]
          (C^_3: leading-Regge off-diagonal = 0 exactly (s1 (a)); OFF3=r4 adds LARGEJ's fitted g3 for comparison)
Y = s1 table (complex), interpolated: nu < 8 at fixed nu, nu >= 8 at fixed t = nu beta, 4-point Lagrange in beta.
Nodes: Gauss-Legendre in p with n = 96 + 3 b-hat (P_l by the Mehler-Dirichlet integral, vectorised), so large b-hat is resolved.
Modes: check (NF02 eigenvalue map), moment (exact first moments M_off(beta), M_1(beta)), lp (NF02 LP + AdS rows, r3.solve).
Usage: S2_DG=300 S2_BMAX=0.3 S2_OFF=psi python3 s2_lp.py check|lp|moment
"""
import os, sys, io, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'omegau')); sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
sys.path.append(os.path.join(HERE, '..', 'robust'))
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import j0
from common import SIGMA
import nf_lp as NL
import r3_redesign as R3

CACHE = os.path.join(HERE, 'cache')
D = R3.D; RR = R3.R


# ------------------------------------------------------------------------------------------------ s1 tables
def load_table(tags=('C2a', 'C2b')):
    T = {}
    for tg in tags:
        fn = os.path.join(CACHE, f's1_table_{tg}.pkl')
        if os.path.exists(fn):
            T.update(pickle.load(open(fn, 'rb')))
    return dict(sorted(T.items()))


class Yint:
    """Y(nu, beta) from a table {beta: (nus, Y)}; Y/beta^2 interpolated (fixed nu below NU0, fixed t above), Lagrange in beta."""
    NU0 = 8.0

    def __init__(self, T):
        self.T = T; self.B = np.array(sorted(T))

    def _one(self, k, nu, be):
        nus, Y = self.T[self.B[k]]; bk = self.B[k]
        Z = Y/bk**2
        out = np.empty(len(nu), complex)
        lo = nu < self.NU0
        out[lo] = np.interp(nu[lo], nus, Z.real) + 1j*np.interp(nu[lo], nus, Z.imag)
        tk = nus*bk; t = nu[~lo]*be
        zr = np.interp(t, tk, Z.real, right=0.0); zi = np.interp(t, tk, Z.imag, right=0.0)
        # beyond the tabulated nu-range (nu > numax or t > tmax) the integrand is set to 0 (|Y|/beta^2 < 1e-2 there)
        out[~lo] = zr + 1j*zi
        return out

    def __call__(self, nu, be):
        nu = np.asarray(nu, float)
        B = self.B
        i = int(np.clip(np.searchsorted(B, be) - 2, 0, len(B) - 4))
        ks = range(i, i + 4)
        xb = np.log(B[i:i + 4]); x = np.log(be)
        out = 0
        for a, k in enumerate(ks):
            w = np.prod([(x - xb[b])/(xb[a] - xb[b]) for b in range(4) if b != a])
            out = out + w*self._one(k, nu, be)
        return out*be**2


class Ycomb:
    """Y_S for S = q^2 + kappa Nq:  Y_S = (Y_q + kappa W_N)/(1 + kappa r),  W_N = o_N/n_q (s3 table x ratio), r = n_N/n_q (s3 ratio)."""
    def __init__(self, Yq, kappa, ntags=('Nqa',)):
        R = pickle.load(open(os.path.join(CACHE, 's3_ratio_Nq.pkl'), 'rb'))
        self.rnu, self.r = R['nu'], R['r']
        T = {}
        for tg in ntags:
            T.update(pickle.load(open(os.path.join(CACHE, f's3_table_{tg}.pkl'), 'rb')))
        T = {be: (nus, Z*self.rf(nus)) for be, (nus, Z) in T.items()}
        self.WN = Yint(dict(sorted(T.items()))); self.Yq = Yq; self.k = kappa
        self.B = np.array(sorted(set(self.WN.B) & set(Yq.B)))

    def rf(self, nu):
        nu = np.asarray(nu, float)
        return np.interp(nu, self.rnu, self.r.real) + 1j*np.interp(nu, self.rnu, self.r.imag)

    def __call__(self, nu, be):
        nu = np.asarray(nu, float)
        r = self.rf(nu)
        return (self.Yq(nu, be) + self.k*self.WN(nu, be))/(1 + self.k*r)


def conical(nus, beta, n=None):
    """P_{-1/2 + i nu}(cosh beta), Mehler-Dirichlet with t = beta - w^2 (regular integrand), vectorised in nu."""
    n = n or int(60 + 1.2*np.max(nus)*beta)
    w, ww = leggauss(n); sb = np.sqrt(beta)
    w = 0.5*sb*(w + 1); ww = 0.5*sb*ww
    t = beta - w*w
    den = np.sqrt(np.maximum(np.cosh(beta) - np.cosh(t), 1e-300))
    f = np.where(w > 0, 2*w/den, 2/np.sqrt(np.sinh(beta)))
    return (np.sqrt(2)/np.pi)*(np.cos(np.outer(nus, t))*(ww*f)).sum(-1)


def basis_poly_coeffs():
    """monomial coefficients (in p) of the D+1 basis functions (1-p)^R p L_i(2p-1) (exact polynomials, degree R+1+D)."""
    from numpy.polynomial import polynomial as Pl, legendre as Lg
    out = []
    for i in range(D + 1):
        L = Lg.leg2poly(np.eye(D + 1)[i])                                 # L_i(y) in powers of y
        # L_i(2p - 1): compose
        c = np.zeros(1)
        for k, a in enumerate(L):
            c = Pl.polyadd(c, a*Pl.polypow([-1.0, 2.0], k))
        c = Pl.polymul(c, Pl.polymul(Pl.polypow([1.0, -1.0], RR), [0.0, 1.0]))
        out.append(c)
    return out


_XCACHE = {}


def chi_basis(Dg, ufun=None, tag='q2'):
    """v4_chi solutions X_i(nu) for each psi2 basis function (min-norm gauge), on v4's nu-grid for Delta_g."""
    key = (Dg, tag)
    if key in _XCACHE:
        return _XCACHE[key]
    fn = os.path.join(CACHE, f's2_chi_{Dg}_{tag}_d{D}.pkl')
    if os.path.exists(fn):
        _XCACHE[key] = pickle.load(open(fn, 'rb')); return _XCACHE[key]
    import v4_chi as C
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    dnu, N = {300: (0.02, 2**14), 1000: (0.04, 2**15), 3000: (0.08, 2**16), 100: (0.01, 2**14)}[Dg]
    Xs = []; res = []
    for c in basis_poly_coeffs():
        S = C.Solver(c, Lm, dnu, N, ufun or C.u_spin, 11.5)
        X, r, it, _ = S.solve_lsqr()
        Xs.append(X); res.append(r)
    out = dict(nu=S.nu, X=np.array(Xs), res=np.array(res), Lam=Lm)
    pickle.dump(out, open(fn, 'wb'))
    _XCACHE[key] = out
    return out


# ------------------------------------------------------------------------------------------------ states
def state_block(Dg, mfl, J, Yf, off='psi', g3=None, chi=None):
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    u = 1/mfl**2
    eta = 1.0 + ((2*J + 1)/mfl)**2/(2*Lm*Lm); be = np.arccosh(eta); bh = (2*J + 1)/mfl
    nq = int(min(2000, 96 + 3*bh))
    xg, wg = leggauss(nq); pn = 0.5*(xg + 1); pw = 0.5*wg
    V = NL.basis_values(pn, D, RR)
    x = pn*np.sqrt(u); y = x*x
    P = conical(Lm*pn, be); jz = j0((2*J + 1)*x)
    n = D + 1
    hs = ['++', '+-'] if J % 2 == 0 else ['+-']
    F = np.zeros((len(hs), len(hs), 2*n)); Nmap = []
    for k, h in enumerate(hs):
        s = SIGMA[h]
        B2b, B3b = R3.flat_basis(h, u, J)
        F[k, k] = np.concatenate([B2b + V @ (pw*2*(P - jz)), B3b + V @ (pw*u*s*(P - jz))])
        Nmap.append((np.concatenate([B2b + V @ (pw*2*(P - jz)), np.zeros(n)]), np.concatenate([np.zeros(n), B3b + V @ (pw*u*s*(P - jz))])))
    if len(hs) == 2:
        eo2 = (2 - y)/(1 - y)**2
        if off == 'psi':
            F[0, 1, :n] = V @ (pw*eo2*Yf(Lm*pn, be).real)
        elif off == 'chi':
            nu = chi['nu']; sel = nu <= Lm
            xx = nu[sel]/Lm*np.sqrt(u); yy = xx*xx
            Yv = Yf(nu[sel], be)*(2 - yy)/(1 - yy)**2
            dnu = nu[1] - nu[0]
            F[0, 1, :n] = (chi['X'][:, sel]*Yv[None, :]).sum(-1).real*dnu/Lm
        elif off == 'none':
            pass
        if g3 is not None:
            t = Lm*pn*be
            F[0, 1, n:] = V @ (pw*u*eo2/(2 - y)*be*be*g3(t))
        F[1, 0] = F[0, 1]
    return dict(mfl=mfl, J=J, DF='s2', hs=hs, F=F, Nmap=Nmap, cov=True, bh=bh, beta=be)


def grid_states(Dg, bmax, mfls, bmin=10.0, nb=60):
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    bhs = np.unique(np.round(np.geomspace(bmin, 2*Lm, nb), 3))
    out = []
    for mfl in mfls:
        for bh in bhs:
            J = int(round((bh*mfl - 1)/2))
            if J < 4:
                continue
            m = mfl*Lm
            alpha = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*m*m))/2
            if alpha < Dg - 2 - 1e-9:
                continue
            be = np.arccosh(1.0 + ((2*J + 1)/mfl)**2/(2*Lm*Lm))
            if be > bmax:
                continue
            out.append((mfl, J))
    return sorted(set(out))


def main():
    mode = sys.argv[1]
    Dg = int(os.environ.get('S2_DG', '300'))
    bmax = float(os.environ.get('S2_BMAX', '0.3'))
    off = os.environ.get('S2_OFF', 'psi')
    tags = tuple(os.environ.get('S2_TABLE', 'C2a,C2b').split(','))
    tag = os.environ.get('S2_TAG', '')
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    Yf = Yint(load_table(tags))
    kap = os.environ.get('S2_KAPPA')
    if kap is not None:
        Yq = Yint({b: v for b, v in load_table(tags).items() if b in pickle.load(open(os.path.join(CACHE, 's3_table_Nqa.pkl'), 'rb'))})
        Yf = Ycomb(Yq, float(kap))
    g3 = None
    if os.environ.get('S2_OFF3', '') == 'r4':
        import r5_universal as R5
        g3 = R5.G3
    chi = chi_basis(Dg) if off == 'chi' else None
    if mode == 'moment':
        # exact first moments at fixed beta: M_off = int_0^oo nu Re Y dnu (Abel mean of the oscillating tail), M_1 = int 2 nu P_l/(e^{2 pi nu}+1)
        from scipy.integrate import cumulative_trapezoid, quad
        pr(f"# s2 moments from tables {tags}: M_off(beta) = int nu Re Y(nu, beta) dnu  (scaling: int t R dt = 47/4?)  vs  M_1(beta)")
        for be in Yf.B:
            if isinstance(Yf, Ycomb):
                nus = Yf.Yq.T[be][0]; Y = Yf(nus, be)
            else:
                nus, Y = Yf.T[be]
            c = cumulative_trapezoid(nus*Y.real, nus, initial=0)
            k = nus*be > 0.6*nus[-1]*be
            M1 = quad(lambda v: 2*v*conical(np.array([v]), be)[0]/(np.exp(2*np.pi*v) + 1), 0, 12)[0]
            pr(f"   beta={be:6.3f}: M_off = {c[k].mean():+8.4f} (+- {c[k].std():.4f}, tail t > {0.6*nus[-1]*be:.0f});  M_1 = {M1:.5f};  "
               f"24 |M_off| M_1(0)/M_1 ... ratio M_off/M_1 = {c[k].mean()/M1:+.1f}")
        open(os.path.join(CACHE, f's2_moment{tag}.txt'), 'w').write(buf.getvalue()); return
    mfls = [float(v) for v in os.environ.get('S2_MFLS', '1.02,1.05,1.1,1.2,1.3,1.5,2.0,3.0,5.0').split(',')]
    sts = grid_states(Dg, bmax, mfls)
    blocks = [state_block(Dg, mfl, J, Yf, off, g3, chi) for mfl, J in sts]
    pr(f"# s2 [{mode}] Delta_g = {Dg} (Lam = {Lm:.1f}); EXACT leading-Regge off-diagonal (tables {tags}, OFF={off}, C3 off "
       f"{'r4 g3' if g3 else '0'}); {len(blocks)} states, beta <= {bmax}; m/Lam in {mfls}  [{time.time() - t0:.0f}s]")
    v_nf = R3.nf_vector()
    if mode == 'check':
        for mfl in mfls:
            bl = sorted([b for b in blocks if b['mfl'] == mfl and len(b['hs']) == 2], key=lambda b: b['bh'])
            if not bl:
                continue
            ls = [(b['bh'], R3.lam(b, v_nf)[0], abs(R3.evalF(b, v_nf)[0, 1])/np.sqrt(R3.evalN(b, v_nf).prod()), b['beta']) for b in bl]
            cross = next((bh for bh, l, o, be in ls if l < 0), None)
            pr(f"  m/Lam = {mfl}: first negative at b-hat = {cross}; (b-hat, beta, lam/N_s, |X|/sqrt(NN)): " +
               " ".join(f"({bh:.0f},{be:.2f},{l:+.2f},{o:.2f})" for bh, l, o, be in ls[::2]))
    elif mode == 'lp':
        G = NL.NFGrid(D, RR)
        R3.summary(blocks, v_nf, pr, 'NF02, exact off-diagonal')
        res = {}
        for dw in [float(v) for v in os.environ.get('S2_DW', '0.0').split(',')]:
            st, v, info = R3.solve(G, blocks, 0.2, dw, Lm, os.environ.get('S2_OBJ', 'A'), v0=v_nf)
            if st != 0:
                pr(f"  LP + AdS rows (dw = {dw}): status {st}: {info.get('msg')}"); res[dw] = None; continue
            pr(f"  LP + AdS rows (dw = {dw}): A = {info['A']:.2f}, B_AdS = {info['B_AdS']:.2f}, bound 300 / 1000 / 3000: "
               f"{R3.bound(info['A'], info['B_AdS'], 300):.3e} / {R3.bound(info['A'], info['B_AdS'], 1000):.3e} / "
               f"{R3.bound(info['A'], info['B_AdS'], 3000):.3e}")
            R3.summary(blocks, v, pr, f'redesign dw={dw}')
            res[dw] = dict(v=v, **info)
        pickle.dump(res, open(os.path.join(CACHE, f's2_lp_{Dg}_b{bmax}_{off}{tag}.pkl'), 'wb'))
    open(os.path.join(CACHE, f's2_{mode}_{Dg}_b{bmax}_{off}{tag}.txt'), 'w').write(buf.getvalue())
    pr(f"[{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()

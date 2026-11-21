"""
bj_ee (BELOW_FOLD_SPINJ.md §1-2): the EE graviton integrand (bg_ee.ee_pieces: w0_ee_common.make_build_even_EE, k = 2 prefactor x (v'-u'),
Lemma-B R for k = 4) at ANY even J >= 4, with the spin-J O-index contraction in closed form (bj_ring.bk_matrices_xs) instead of 2J+1
null-circle samples.

  pref_J[a][c] = N_J R^J sum_k P_J^(k)(t) X_k[a][c],
  X_k[a][c] = sum_{i,j} D[i][j]/(kappa2_i kappa2_j) sum_{n,n'} F_a[4-i][n] Bk_k[n,n'] G_c[4-j][n']
(F_a, G_c: Fourier modes in (z2, z5) / (z4, z5) of the quartic remainders = bg_ee.tto at "J = 4 powers" V5^{4-b15-b25} with the
conserved coefficients of the true (Delta, J); D = E^2 (bg_ee.e_factor), Fourier in (z2, z4)).  R^J = (V^2 W^2)^{J/2} cancels the
J-dependent parts of the 3-pt exponents exactly (V^2 = -m15 m25/m12, W^2 = -m35 m45/m34; asserted), so the logs are the J = 0 ones;
N_J (J only) is dropped.
J/Delta resummation (bg_scalarJ form): t = 1 + a^2 Dhat,  P_J^(k)(t) = sum_n gam_{k,n} Dhat^n,  gam_{k,n} = c_{n+k}(J) (n+k)!/n! a^{2n}
(exact numbers), c_m(J) = P_J^(m)(1)/m!;   pref_J/(N_J R^J) = sum_{n <= N} Dhat^n Y_n,  Y_n = sum_k gam_{k,n} X_k,
truncated triangularly: the n-th term to relative order r_max(n) = max(0, Kr - floor(n/NTRI)) (checked by varying N, NTRI).
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import factorial
import bg_ser as S, bg_geom as GG, bg_ee as BE
from bg_geom import sub, dot
import bj_ring as BR
import bj_xs as XS_

# bg_ser.Poly prunes coefficients below 1e-15 x the largest one.  With high powers of Dhat (degree ~2n in radial atoms whose values in
# the integration region are O(10)) tiny coefficients carry O(1) values: pruning them is a precision bug (found here: Dhat^5 x 1e-16 added
# to an O(1) polynomial changed its value by 1e-4).  In every bj_* process the prune is VALUE-weighted: a term is dropped only if
# |c x monomial| < 1e-16 x max over the polynomial at each of two reference points (radial atoms p_i = 1 and 4; angular atoms generic).
_LPT = np.log(np.array([[1, 1, 1, 1, .5, .75, .8, .6, .3, .95, .7, 1.3],
                        [4, 4, 4, 4, .5, .75, .8, .6, .3, .95, .7, 1.3]], float))


def value_keep(E, C, rel=1e-16):
    """keep a term if |c x monomial| > rel x max_terms |c x monomial| at EITHER reference point."""
    if C.size == 0:
        return np.zeros(0, bool)
    w = np.log(np.abs(C) + 1e-300)[:, None] + E.astype(float) @ _LPT.T
    return np.any(w > w.max(0)[None, :] + np.log(rel), axis=1)


def _prune_safe(self, rel=None):
    if self.C.size == 0:
        return self
    k = value_keep(self.E, self.C)
    if k.all():
        return self
    return S.Poly(self.E[k], self.C[k])


S.Poly.prune = _prune_safe


def tv(c):
    """coefficient vector of a z-linear form (sum_m z_m c_m) -> lightcone vector A with dot(A, z) = that form."""
    return [c[1].scale(-2.0), c[0].scale(-2.0), c[2]]


def mons4(mons):
    return [(a1, a2, 4 - b15 - b25, b12, b15, b25) for (a1, a2, a5, b12, b15, b25) in mons]


def gam(J, Dl, kmax=8, N=12):
    """gam[k][n] = c_{n+k}(J) (n+k)!/n! a^{2n}  (a = 1/Dl), exact."""
    a2 = 1.0/Dl**2
    lc = [0.0]                                   # log|c_m(J)|, sign
    sg = [1.0]
    for m in range(1, N + kmax + 1):
        f = (J*(J + 1) - m*(m - 1))/(2.0*m*m)
        if f == 0:
            lc.append(-np.inf); sg.append(0.0); continue
        lc.append(lc[-1] + np.log(abs(f))); sg.append(sg[-1]*np.sign(f))
    out = []
    for k in range(kmax + 1):
        row = []
        for n in range(N + 1):
            m = n + k
            if sg[m] == 0 or not np.isfinite(lc[m]):
                row.append(0.0); continue
            row.append(sg[m]*np.exp(lc[m] + np.log(factorial(m)/factorial(n)) + n*np.log(a2)))
        out.append(row)
    return out


def pad_add(A, B, top):
    """A + B keeping absolute orders up to `top` (missing coefficients = 0, not a truncation of the sum)."""
    if A is None:
        A, B = B, None
    if B is None:
        return S.ASer(A.d, [A.c[k] if k < len(A.c) else S.Poly() for k in range(top - A.d + 1)])
    d = min(A.d, B.d)
    c = [S.Poly() for _ in range(top - d + 1)]
    for X in (A, B):
        for k, p in enumerate(X.c):
            i = X.d + k - d
            if 0 <= i < len(c) and len(p):
                c[i] = c[i] + p
    return S.ASer(d, c)


def trunc_abs(A, top):
    return S.ASer(A.d, A.c[:max(1, top - A.d + 1)])


def vw(Gm):
    B = Gm.basis
    V = tv([Gm.V(5, B[m], 1, 2) for m in range(3)])
    W = tv([Gm.V(5, B[m], 3, 4) for m in range(3)])
    return V, W


def tser(V, W):
    """t = V.W/sqrt(V^2 W^2) with the branch used by bj_ring.bk_matrices_xs."""
    VV, WW = dot(V, V), dot(W, W)
    c0, e0, eps = (VV*WW).rel('V^2W^2')
    r = S.series_fn(eps, 'pow', -0.5)
    iR = S.ASer(-(VV*WW).d//2, [p.mul_mono(-e0//2, c0**-0.5) for p in r.c])
    return dot(V, W)*iR


def prune_keys(Bk, t, F, Gt, J, Dl, tol, nup=12, rscale=(1.0, 2.5)):
    """B entries whose contribution |P^(k)(t) B_k[n,n'] F_n G_n'| (sample atom points with radial values scaled by rscale, a = 1/Dl)
    stays below tol x the largest contribution at every point are dropped."""
    from bj_ring import legendre_derivs
    pts0 = S.atom_points(nup, 99)
    contrib = {}
    for rs in rscale:
        pts = pts0.copy(); pts[:, 0:4] *= np.sqrt(rs)
        a = 1.0/Dl
        tv_ = t.ev(pts, a)
        Pd = np.array([legendre_derivs(J, complex(tt), 8) for tt in tv_])          # (pts, 9)
        wF = np.zeros((len(pts), 9)); wG = np.zeros((len(pts), 9))
        for kk in range(9):
            wF[:, kk] = np.max([np.abs(F[a_][j][kk].ev(pts, a)) for a_ in range(len(F)) for j in range(5)], axis=0)
            wG[:, kk] = np.max([np.abs(Gt[c][j][kk].ev(pts, a)) for c in range(len(Gt)) for j in range(5)], axis=0)
        for k in range(9):
            for (n, n2), b in Bk[k].items():
                v = np.abs(Pd[:, k]*b.ev(pts, a))*wF[:, n + 4]*wG[:, n2 + 4]
                contrib[(k, n, n2)] = np.maximum(contrib.get((k, n, n2), 0*v), v)
    mx = np.max(np.array(list(contrib.values())), axis=0)
    keep = {k: set() for k in range(9)}
    for (k, n, n2), v in contrib.items():
        if np.any(v > tol*mx):
            keep[k].add((n, n2))
    return keep


def pieces(depth, Delta, J, verbose=False, conserved=None, prune=1e-10):
    """J-generic EE pieces at geometric depth `depth`: logs (J = 0 3-pt exponents), vmu, R, X_k[a][c] (k = 0..8), Dhat, t."""
    import l5_spin as S5
    t0 = time.time()
    cons = conserved or S5.conserved
    mons, sol, _ = cons(float(Delta), int(J))
    mons_t, solt, _ = cons(3.0 - float(Delta), int(J))
    assert tuple(mons) == tuple(mons_t)
    m4 = mons4(mons)
    logs_sc, polys_sc, G, nm = GG.scalar_pieces(depth)
    Gm = BE.EEGeom(depth, G)
    V, W = vw(Gm)
    pts = S.atom_points(6, 7)
    for (vv, num, den) in ((dot(V, V), nm['m15']*nm['m25'], nm['m12']), (dot(W, W), nm['m35']*nm['m45'], nm['m34'])):
        lhs = (vv*den).ev(pts, 1e-3); rhs = (-num).ev(pts, 1e-3)
        assert np.max(np.abs(lhs/rhs - 1)) < 1e-8, ('V^2 identity', np.max(np.abs(lhs/rhs - 1)))
    flip = np.real(complex(tser(V, W).lead('t0')[0])) < 0
    if flip:
        W = [w.scale(-1.0) for w in W]
    G1 = GG.geometry(depth + 1); V1, W1 = vw(BE.EEGeom(depth + 1, G1))
    if flip:
        W1 = [w.scale(-1.0) for w in W1]
    Bk, t, _ = BR.bk_matrices_xs(J, V1, W1, depth, float(Delta))
    G2 = GG.geometry(depth + 2); V2, W2_ = vw(BE.EEGeom(depth + 2, G2))
    if flip:
        W2_ = [w.scale(-1.0) for w in W2_]
    t = tser(V2, W2_)                                   # to order depth + 2, so that Dhat = (t - 1)/a^2 has depth + 1 orders
    if verbose:
        print(f'   Bk {time.time() - t0:.0f}s (W flipped {flip}); entries per k ' + str([len(b) for b in Bk]), flush=True)
    W2, W4 = BE.null_circle(2), BE.null_circle(4)
    x10, x30 = sub(Gm.x[1], Gm.x[0]), sub(Gm.x[3], Gm.x[0])
    F = BE.tto(Gm, 1, x10, 2, W2, 5, W4, 4, m4, sol)
    Gt = BE.tto(Gm, 3, x30, 4, W2, 5, W4, 4, m4, solt)
    F = [BE.fourier2(f, 2, 4) for f in F]; Gt = [BE.fourier2(g, 2, 4) for g in Gt]
    E, v042 = BE.e_factor(Gm, W2)
    D = [[E[i][j]*E[i][j] for j in range(5)] for i in range(5)]
    D = BE.fourier2(D, 2, 2)
    if verbose:
        sz = lambda A: sum(len(p) for p in A.c)
        print(f'   F, G, D {time.time() - t0:.0f}s; sizes F {sz(F[0][2][4])} G {sz(Gt[0][2][4])} D {sz(D[2][2])}', flush=True)
    k2 = BE.kappa(2)
    na, nc = len(F), len(Gt)
    Em = [[[BE.lin([1.0/k2[i] for i in range(5)], [D[i][j]*Fa[4 - i][kk] for i in range(5)]) for kk in range(9)] for j in range(5)]
          for Fa in F]
    X = [[[None]*nc for _ in range(na)] for _ in range(9)]
    keep = prune_keys(Bk, t, F, Gt, J, float(Delta), prune) if prune else {k: set(Bk[k]) for k in range(9)}
    if verbose:
        print(f'   kept B entries per k: {[len(keep[k]) for k in range(9)]} of {[len(b) for b in Bk]}', flush=True)
    for k in range(min(8, J) + 1):
        if not keep[k]:
            continue
        for c in range(nc):
            T = [[None]*9 for _ in range(5)]                    # T[j][n] = sum_n' Bk[n, n'] G_c[j][n']
            for j in range(5):
                for kk in range(9):
                    T[j][kk] = XS_.macc([(Bk[k][(kk - 4, kk2 - 4)], Gt[c][j][kk2], 1.0) for kk2 in range(9)
                                         if (kk - 4, kk2 - 4) in keep[k]])
            for a in range(na):
                X[k][a][c] = XS_.macc([(Em[a][j][kk], T[4 - j][kk], 1.0/k2[j]) for j in range(5) for kk in range(9)])
        if verbose:
            print(f'   X_{k} {time.time() - t0:.0f}s  d/K ' + str([(X[k][a][c].d, X[k][a][c].K) if X[k][a][c] is not None else None
                                                                for a in range(na) for c in range(nc)]), flush=True)
    dl = (t - 1.0).norm()
    assert dl.d >= 2, dl.d
    Dh = S.ASer(dl.d - 2, dl.c)
    keep = ('A', 'xsq', 'FP1', 'FP2', 'FP3', 'dual', 'cel', 'meas')
    logs = [l for l in logs_sc if l[0] in keep]
    logs += [('3pt', nm['m12'], -10/2, 0.0), ('3pt', nm['m34'], -7/2, 0.0), ('3pt', nm['m35'], -3/2, 0.0),
             ('3pt', nm['m45'], -3/2, 0.0), ('v042', v042, -4.0, 0.0)]
    vmu = nm['vp'] - nm['up']
    Rk4 = (vmu*vmu - 1.0)*(nm['up']*nm['vp']).recip('u\'v\'').scale(0.25)
    return dict(logs=logs, vmu=vmu, R=Rk4, X=X, Dh=Dh, t=t, flip=flip, J=J, Delta=Delta, depth=depth, time=time.time() - t0)


def pref_sum(P, N, Kr, ntri=1):
    """sum_n Dhat^n Y_n for each (a, c); triangular truncation r_max(n) = max(0, Kr - floor(n/ntri)) relative to the n = 0 leading power."""
    X, Dh, J, Dl = P['X'], P['Dh'], P['J'], P['Delta']
    g = gam(J, Dl, 8, N)
    na, nc = len(X[0]), len(X[0][0])
    out = [[None]*nc for _ in range(na)]
    for a in range(na):
        for c in range(nc):
            Y = []
            for n in range(N + 1):
                acc = None
                for k in range(9):
                    if X[k][a][c] is None or g[k][n] == 0:
                        continue
                    v = X[k][a][c].scale(g[k][n])
                    acc = v if acc is None else pad_add(acc, v, min(acc.d + acc.K, v.d + v.K))
                Y.append(acc)
            d0 = Y[0].d
            top0 = d0 + Kr
            tot = pad_add(trunc_abs(Y[0], top0), None, top0)
            Dn = None
            for n in range(1, N + 1):
                if Y[n] is None:
                    continue
                rmax = max(0, Kr - n//ntri)
                Dn = Dh.trunc(rmax) if Dn is None else (Dn.trunc(rmax)*Dh.trunc(rmax))
                term = Dn*trunc_abs(Y[n], d0 + rmax)
                tot = pad_add(tot, term, top0)
            out[a][c] = tot
    return out


NUS = np.round((np.arange(1025) + 0.5)*0.04, 4)


def images(P, N, Kr, ntri=1, nus=NUS, kinds=('EE2', 'EE4'), verbose=False):
    """I_kind[a][c](nu) = Delta^{-d_pref} exp(lconst) sum_j (i nu/2)^j sum_k Delta^{-k} M[(k, j)] -- pref's own a-power kept so the four
    (a, c) entries are mutually normalised (the logs, vmu, R and B-correction are common)."""
    t0 = time.time()
    pr = pref_sum(P, N, Kr, ntri)
    Bf = GG.bfac(Kr)
    out = dict(Delta=P['Delta'], J=P['J'], K=Kr, N=N, ntri=ntri, nus=np.asarray(nus), res={})
    na, nc = len(pr), len(pr[0])
    for a in range(na):
        for c in range(nc):
            for kind in kinds:
                polys = [('vmu', P['vmu']), ('pref', pr[a][c])] + ([('R', P['R'])] if kind == 'EE4' else [])
                A = GG.assemble(P['logs'], polys, Kr)
                M = {}
                for j, cj in A['corr'].items():
                    Pk = S.lmul(cj, Bf.c, Kr)
                    for k in range(Kr + 1):
                        M[(k, j)] = S.integrate_fast(Pk[k].E, Pk[k].C, A['E0'], A['E1'], out['nus']) if len(Pk[k]) else 0*out['nus']
                out['res'][(kind, a, c)] = dict(E0=A['E0'], E1=A['E1'], lconst=A['lconst'], M=M, dpref=pr[a][c].d,
                                                nterms=[len(p) for p in A['corr'][0]])
    if verbose:
        print(f'   images N={N} Kr={Kr} ntri={ntri}: {time.time() - t0:.0f}s', flush=True)
    return out


def value(rec, kind, a=0, c=0, Kt=None):
    """I(nu) incl. Delta^{-d_pref}; per-order pieces (k x nu)."""
    r = rec['res'][(kind, a, c)]
    Kt = rec['K'] if Kt is None else Kt
    Dl = rec['Delta']
    nus = rec['nus']
    pre = np.exp(r['lconst'][0] + 1j*nus*r['lconst'][1] - r['dpref']*np.log(Dl))
    pieces = np.zeros((Kt + 1, len(nus)), complex)
    for (k, j), v in r['M'].items():
        if k <= Kt:
            pieces[k] += Dl**(-k)*(0.5j*nus)**j*v
    return pre*pieces.sum(0), pre[None, :]*pieces

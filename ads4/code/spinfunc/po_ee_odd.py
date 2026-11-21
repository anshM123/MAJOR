"""
po_ee_odd (PARITY_ODD_CLOSE.md sect. 2): the PARITY-ODD <T T O> structure in the EE graviton integrand on the analytic-radial a-series
algebra of bg_ser / bg_geom (BELOW_FOLD2), i.e. the bg_ee construction with the conserved symmetric parity-odd structure of
p1_tto_odd (CHL's X^O_opp at odd J, X^O_same at even J) in place of the parity-even l5_spin set.

Term by term a copy of p3_odd_action.tto_odd_samples / eps_cof (numerators: eps(Z_i, Z_j, X_a, X_b, X_5) X_ij x even monomials, real
coefficient vector P1.conserved_odd) and of w0_ee_common.make_build_odd_EE (kind 'EE': D = E (x) E, f = v' - u', the EE slow exponent
_gs), with the extra odd-structure factors (m12 m15 m25)^{-1/2} (12 side) and (m34 m35 m45)^{-1/2} (34 side, O~ at 3 - Delta) moved into
the log list (as bg_ee moves the even 3-pt powers).  The Lemma-B factor R gives k = 4.  Read-only imports: bg_ser, bg_geom, bg_ee,
bg_eeq, p1_tto_odd.
Pointwise check against the project's exact point-wise odd integrand (p8_odd_pw.P8 + w0 make_build_odd_EE, real section of bf_core):
  python3 po_ee_odd.py check Delta J [depth]
Production (state + reference, as bg_eeq.images, 1 x 1):
  python3 po_ee_odd.py run Delta J K depth   -> cache/po_eeq_D{Delta}_J{J}_K{K}.pkl
"""
import os, sys, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_ser as S
from bg_ser import ASer
import bg_geom as GG
from bg_geom import dot, sub
import bg_ee as BE

C = os.path.join(HERE, 'cache')


def _z(K):
    return ASer.const(0.0, K)


def vec5X(x, K):
    return [ASer.const(1.0, K), dot(x, x), x[0], x[1], x[2]]


def vec5Z(x, z, K):
    return [None, dot(x, z).scale(2.0), z[0], z[1], z[2]]          # None = exact zero (first component)


def _mul(a, b):
    if a is None or b is None:
        return None
    return a*b


def _add(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return a + b


def eps_cof(rows):
    """copy of p3_odd_action.eps_cof on ASer (None = 0): cofactors c_i with det([r0; rows]) = sum_i r0_i c_i."""
    R1, R2, R3, R4 = rows

    def m2(A, B, i, j):
        return _add(_mul(A[i], B[j]), None if _mul(A[j], B[i]) is None else -_mul(A[j], B[i]))
    P12 = {(i, j): m2(R1, R2, i, j) for i in range(5) for j in range(i + 1, 5)}
    P34 = {(i, j): m2(R3, R4, i, j) for i in range(5) for j in range(i + 1, 5)}
    c = []
    for i in range(5):
        cols = [k for k in range(5) if k != i]
        d = None
        for a in range(4):
            for b in range(a + 1, 4):
                rest = [cols[k] for k in range(4) if k not in (a, b)]
                perm = [a, b] + [k for k in range(4) if k not in (a, b)]
                sgn = round(np.linalg.det(np.eye(4)[perm]))
                t = _mul(P12[(cols[a], cols[b])], P34[(rest[0], rest[1])])
                if t is not None:
                    d = _add(d, t if sgn > 0 else -t)
        c.append(None if d is None else (d if i % 2 == 0 else -d))
    return c


def _dotc(vs, c):
    """sum_i vs_i c_i (None-aware)."""
    out = None
    for v, ci in zip(vs, c):
        out = _add(out, _mul(v, ci))
    return out


def tto_odd(Gm, ia, za, ib, zbs, i5, z5s, J, mons, sol):
    """p3_odd_action.tto_odd_samples on ASer: array [nb][n5] of ASer (odd numerators, incl. eps.X_ij; no denominators)."""
    K = Gm.K
    B = Gm.basis
    xa, xb, x5 = Gm.x[ia], Gm.x[ib], Gm.x[i5]
    Xa, Xb, X5 = vec5X(xa, K), vec5X(xb, K), vec5X(x5, K)
    Za = vec5Z(xa, za, K)
    Zb = [vec5Z(xb, B[m], K) for m in range(3)]
    Z5 = [vec5Z(x5, B[k], K) for k in range(3)]
    V1 = Gm.V(ia, za, ib, i5)
    V2b = [Gm.V(ib, B[m], i5, ia) for m in range(3)]
    V5b = [Gm.V(i5, B[m], ia, ib) for m in range(3)]
    H12b = [Gm.H(ia, za, ib, B[m]) for m in range(3)]
    H15b = [Gm.H(ia, za, i5, B[m]) for m in range(3)]
    H25b = [[Gm.H(ib, B[m], i5, B[k]) for k in range(3)] for m in range(3)]
    Xab, Xa5, Xb5 = Gm.sq(ia, ib), Gm.sq(ia, i5), Gm.sq(ib, i5)
    c = eps_cof([Za, Xa, Xb, X5])
    E12b = [-_dotc(Zb[m], c) for m in range(3)]
    E15b = [-_dotc(Z5[k], c) for k in range(3)]
    E25m = []
    for m in range(3):
        cm = eps_cof([Zb[m], Xa, Xb, X5])
        E25m.append([-_dotc(Z5[k], cm) for k in range(3)])
    nb, n5 = zbs.shape[0], z5s.shape[0]
    lin = BE.lin
    V2 = [lin(zbs[i], V2b) for i in range(nb)]
    V5 = [lin(z5s[i], V5b) for i in range(n5)]
    H12 = [lin(zbs[i], H12b) for i in range(nb)]
    H15 = [lin(z5s[i], H15b) for i in range(n5)]
    H25 = [[lin([zbs[i, m]*z5s[j, k] for m in range(3) for k in range(3)], [H25b[m][k] for m in range(3) for k in range(3)])
            for j in range(n5)] for i in range(nb)]
    E12 = [lin(zbs[i], E12b)*Xab for i in range(nb)]
    E15 = [lin(z5s[j], E15b)*Xa5 for j in range(n5)]
    E25 = [[lin([zbs[i, m]*z5s[j, k] for m in range(3) for k in range(3)], [E25m[m][k] for m in range(3) for k in range(3)])*Xb5
            for j in range(n5)] for i in range(nb)]
    one = ASer.const(1.0, K)

    def pw(x, kmax):
        out = [one]
        for _ in range(kmax):
            out.append(out[-1]*x)
        return out
    V1p = pw(V1, 2)
    V2p = [pw(v, 2) for v in V2]; V5p = [pw(v, J) for v in V5]; H12p = [pw(v, 2) for v in H12]; H15p = [pw(v, 2) for v in H15]
    H25p = [[pw(v, 2) for v in row] for row in H25]
    used = [i for i in range(len(mons)) if abs(sol[i, 0]) >= 1e-14]
    res = [[None]*n5 for _ in range(nb)]
    for ib_ in range(nb):
        for i5_ in range(n5):
            acc = None
            for i in used:
                typ, a1, a2, a5, b12, b15, b25 = mons[i]
                t = V1p[a1]*V2p[ib_][a2]*H12p[ib_][b12]*V5p[i5_][a5]*H15p[i5_][b15]*H25p[ib_][i5_][b25]
                t = t*(E12[ib_] if typ == 12 else (E15[i5_] if typ == 15 else E25[ib_][i5_]))
                t = t.scale(float(sol[i, 0]))
                acc = t if acc is None else acc + t
            res[ib_][i5_] = acc
    return res


def odd_pieces(K, Delta, J, verbose=False):
    """EE k = 2 integrand with the parity-odd structure on both sides (1 x 1): logs, vmu, pref, R (as bg_ee.ee_pieces)."""
    import p1_tto_odd as P1
    mo, so, _ = P1.conserved_odd(float(Delta), int(J))
    mt, st, _ = P1.conserved_odd(3.0 - float(Delta), int(J))
    assert so is not None and st is not None, ('no symmetric odd structure', Delta, J)
    logs_sc, polys_sc, G, nm = GG.scalar_pieces(K)
    Gm = BE.EEGeom(K, G)
    W2, WJ = BE.null_circle(2), BE.null_circle(J)
    x10, x30 = sub(Gm.x[1], Gm.x[0]), sub(Gm.x[3], Gm.x[0])
    t = time.time()
    F = tto_odd(Gm, 1, x10, 2, W2, 5, WJ, J, list(mo), so)
    if verbose:
        print(f'   F done {time.time() - t:.0f}s', flush=True)
    Gt = tto_odd(Gm, 3, x30, 4, W2, 5, WJ, J, list(mt), st)
    if verbose:
        print(f'   Gt done {time.time() - t:.0f}s', flush=True)
    F = BE.fourier2(F, 2, J); Gt = BE.fourier2(Gt, 2, J)
    E, v042 = BE.e_factor(Gm, W2)
    D = [[E[i][j]*E[i][j] for j in range(5)] for i in range(5)]
    D = BE.fourier2(D, 2, 2)
    k2, kJ = BE.kappa(2), BE.kappa(J)
    n5 = 2*J + 1
    Em = [[BE.lin([1.0/k2[i] for i in range(5)], [D[i][j]*F[4 - i][k] for i in range(5)]) for k in range(n5)] for j in range(5)]
    terms = []; co = []
    for j in range(5):
        for k in range(n5):
            terms.append(Em[j][k]*Gt[4 - j][n5 - 1 - k]); co.append(1.0/(k2[j]*kJ[k]))
    pref = BE.lin(co, terms)
    if verbose:
        print(f'   pref done {time.time() - t:.0f}s; d = {pref.d}, terms per order', [len(p) for p in pref.c], flush=True)
    keep = ('A', 'xsq', 'FP1', 'FP2', 'FP3', 'dual', 'cel', 'meas')
    logs = [l for l in logs_sc if l[0] in keep]
    logs += [('3pt', nm['m12'], -(10 - J)/2, 0.0), ('3pt', nm['m15'], -J/2, 0.0), ('3pt', nm['m25'], -J/2, 0.0),
             ('3pt', nm['m34'], -(7 - J)/2, 0.0), ('3pt', nm['m35'], -(3 + J)/2, 0.0), ('3pt', nm['m45'], -(3 + J)/2, 0.0),
             ('v042', v042, -4.0, 0.0),
             # parity-odd extra factors (p3: ex12 = (m12 m15 m25)^(-1/2), ex34 = (m34 m35 m45)^(-1/2))
             ('odd', nm['m12'], -0.5, 0.0), ('odd', nm['m15'], -0.5, 0.0), ('odd', nm['m25'], -0.5, 0.0),
             ('odd', nm['m34'], -0.5, 0.0), ('odd', nm['m35'], -0.5, 0.0), ('odd', nm['m45'], -0.5, 0.0)]
    vmu = nm['vp'] - nm['up']
    Rk4 = (vmu*vmu - 1.0)*(nm['up']*nm['vp']).recip('u\'v\'').scale(0.25)
    return dict(logs=logs, vmu=vmu, pref=pref, R=Rk4, G=G, nm=nm, v042=v042)


# ------------------------------------------------------------------------------------------------- pointwise check
def atoms_of(Smp, idx):
    P = np.empty((len(idx), S.NV))
    X1, X2, Y1, Y2 = Smp.X1[idx], Smp.X2[idx], Smp.Y1[idx], Smp.Y2[idx]
    P[:, 0], P[:, 1], P[:, 2], P[:, 3] = np.sqrt(X1), np.sqrt(X2), np.sqrt(Y1), np.sqrt(Y2)
    s = Smp.s[idx]; P[:, 4] = s; P[:, 5] = 1 - s*s
    ph = Smp.phi[idx]; P[:, 6] = np.cos(ph); P[:, 7] = np.sin(ph)
    al = Smp.al[idx]; P[:, 8] = np.cos(al); P[:, 9] = np.sin(al); P[:, 10] = 1 - P[:, 8]; P[:, 11] = 1 + P[:, 8]
    return P


def check(Delta, J, depth=4, npts=8):
    """ASer odd prefactor x odd extra factors vs the exact point-wise integrand prefactor (p8 P8 + w0 odd EE build), n0 = +1 sheet."""
    import bf_core as BC
    import p8_odd_pw as P8M
    import w0_ee_common as W0
    PW = BC.pw()
    a = 1.0/Delta
    Smp = BC.Sample(6, 5)
    Z, lj, psi = Smp.raw(a, 1)
    idx = np.arange(npts)
    t = time.time()
    P = odd_pieces(depth, Delta, J, verbose=True)
    print(f'  pieces depth {depth}: {time.time() - t:.0f}s', flush=True)
    At = atoms_of(Smp, idx)
    # the n0 = +1 sheet of bf_core: psi in (0, pi); bg's alpha = psi + phi  (bg_geom: psi = alpha - phi)
    At[:, 8] = np.cos(psi[idx] + Smp.phi[idx]); At[:, 9] = np.sin(psi[idx] + Smp.phi[idx])
    At[:, 10] = 1 - At[:, 8]; At[:, 11] = 1 + At[:, 8]
    nm = P['nm']
    ex = 1.0
    for k in ('m12', 'm15', 'm25', 'm34', 'm35', 'm45'):
        ex = ex*np.exp(-0.5*np.log(nm[k].ev(At, a).astype(complex)))
    ser = [sum(p.ev(At)*a**(P['pref'].d + k) for k, p in enumerate(P['pref'].c[:kk + 1])) for kk in range(len(P['pref'].c))]
    build = W0.make_build_odd_EE(P8M.P8)
    zz = np.vstack([BC._ref_point(a)[None, :], Z[idx]])
    PW.set_batch(zz.shape[0])
    zl = [zz[:, i].copy() for i in range(7)]
    out = build(zl, 1.0, Delta, J, 'EE', 1.0)
    ex_pw = out[2][-1, -1, 1:]
    lines = [f'# po_ee_odd check: Delta={Delta} J={J}, depth {depth}, {npts} real-section points (n0 = +1 sheet)']
    for kk, sv in enumerate(ser):
        r = sv*ex/ex_pw
        lines.append(f'  order <= a^{kk}: ratio ASer/exact = ' + ' '.join(f'{v.real:+.6f}{v.imag:+.6f}i' for v in r[:4])
                     + f' | max |ratio - 1| = {np.abs(r - 1).max():.2e}, max |ratio/ratio_0 - 1| = {np.abs(r/r[0] - 1).max():.2e}')
    print('\n'.join(lines), flush=True)
    return lines


# ------------------------------------------------------------------------------------------------- images (bg_eeq.images, 1 x 1)
NUS = np.round((np.arange(1025) + 0.5)*0.04, 4)


def images(Dl, J, K, depth, nus=NUS, kinds=('EE4',), verbose=True):
    t = time.time()
    P = odd_pieces(depth, Dl, J, verbose=verbose)
    Bf = GG.bfac(K)
    out = dict(Delta=Dl, J=J, K=K, depth=depth, nus=np.asarray(nus), res={}, struct='odd')
    if verbose:
        print(f'  pieces D={Dl} J={J}: {time.time() - t:.0f}s; pref d/K {P["pref"].d}/{P["pref"].K}', flush=True)
    for kind in kinds:
        polys = [('vmu', P['vmu']), ('pref', P['pref'])] + ([('R', P['R'])] if kind == 'EE4' else [])
        A = GG.assemble(P['logs'], polys, K)
        M = {}
        for j, cj in A['corr'].items():
            Pk = S.lmul(cj, Bf.c, K)
            for k in range(K + 1):
                M[(k, j)] = S.integrate_fast(Pk[k].E, Pk[k].C, A["E0"], A["E1"], out["nus"]) if len(Pk[k]) else 0*out['nus']
        out['res'][(kind, 0, 0)] = dict(E0=A['E0'], E1=A['E1'], lconst=A['lconst'], M=M, nterms=[len(p) for p in A['corr'][0]])
        if verbose:
            print(f'   {kind}: terms {out["res"][(kind, 0, 0)]["nterms"]}, {time.time() - t:.0f}s', flush=True)
    return out


def fn(Dl, J, K):
    return os.path.join(C, 'po_eeq_D%g_J%d_K%d.pkl' % (Dl, J, K))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'check':
        Dl, J = float(sys.argv[2]), int(sys.argv[3]); dep = int(sys.argv[4]) if len(sys.argv) > 4 else 4
        L = check(Dl, J, dep)
        open(os.path.join(C, 'po_ee_odd_check.txt'), 'a').write('\n'.join(L) + '\n')
    else:
        Dl, J, K, depth = float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        t0 = time.time()
        rec = images(Dl, J, K, depth)
        pickle.dump(rec, open(fn(Dl, J, K) + '.tmp', 'wb')); os.replace(fn(Dl, J, K) + '.tmp', fn(Dl, J, K))
        print(f'po_ee_odd D={Dl} J={J} K={K} depth={depth}: done {time.time() - t0:.0f}s', flush=True)

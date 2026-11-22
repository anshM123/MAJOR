"""
lj_rs_norm (HLJ_CLOSE.md §2.2): the ANALYTIC normalisation of the closed-form spin-J EE engine (bj_ee images) against the direct
crystal records (z2 / n3 contour; w0 EE integrand), replacing bj_chain.calibrate's fitted constant.

The engine integrand (bg_geom.assemble of bj_ee.pieces) is the a-expansion of the SAME holomorphic integrand as the direct one, written in
the scaled variables of the real section D with every overall a-power (log factors, v'-u', measure) dropped and its own branch
constants.  Hence, point-wise on D (where both are defined and the a-series converges absolutely at small scaled radii),
    f_eng(p; nu) / f_dir(p; nu) = c(Delta, J, nu),      p-independent (up to the engine's truncation O(a^{Kr+1} X^n)),
and by HRS_CLOSE Prop. RS (one populated cohomology class, c_Gamma common to state and reference) the crystal integrals obey
    V_st / V_ref = [a_st/a_ref]^6 [c_ref/c_st] (I_st / I_ref)          (bf_core drops a^6 from d^7z).
So rho_an = kappa_st/kappa_ref (bj_chain notation, V = kappa I) = (a_st/a_ref)^6 c_ref/c_st, with no fitted parameter.
f_dir is bf_core.grav_parts (n1_pw point-wise w0 EE integrand, sheet n0 = +1, psi = alpha - phi in (0, pi)); f_eng is evaluated from the
cached bj pieces exactly as bj_ee.images integrates it (same pref_sum truncation, assemble, B-factor).
Usage: python3 lj_rs_norm.py J Dstate depth N Kr [patch]     -> cache/lj_rs_norm_J<J>_D<D>_d<depth>.txt (+ .pkl)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0


class Pts:
    """minimal bf_core.Sample stand-in: given scaled points (X1, X2, s, Y1, Y2, phi, psi), psi in (0, pi), sheet n0 = +1."""
    def __init__(self, X1, X2, s, Y1, Y2, phi, psi):
        self.X1, self.X2, self.s, self.Y1, self.Y2, self.phi = (np.asarray(v, float) for v in (X1, X2, s, Y1, Y2, phi))
        self.psi = np.asarray(psi, float)
        self.al = self.psi + self.phi           # bf_core.raw: psi = mod(al - n0 phi, pi) -> psi (psi in (0, pi))
        self.theta = 0.0
        self.n = len(self.X1)

    def raw(self, a, n0):
        import bf_core as B
        return B.Sample.raw(self, a, n0)

    def atoms(self):
        P = np.empty((self.n, 12))
        P[:, 0] = np.sqrt(self.X1); P[:, 1] = np.sqrt(self.X2); P[:, 2] = np.sqrt(self.Y1); P[:, 3] = np.sqrt(self.Y2)
        P[:, 4] = self.s; P[:, 5] = 1 - self.s**2
        P[:, 6] = np.cos(self.phi); P[:, 7] = np.sin(self.phi)
        al = self.psi + self.phi               # engine alpha = psi + phi (n0 = +1)
        al = np.mod(al + np.pi, 2*np.pi) - np.pi
        P[:, 8] = np.cos(al); P[:, 9] = np.sin(al); P[:, 10] = 1 - P[:, 8]; P[:, 11] = 1 + P[:, 8]
        return P


def make_points(n, seed=1, rscale=1.0):
    rng = np.random.default_rng(seed)
    X1 = rng.uniform(0.3, 1.5, n)*rscale; X2 = rng.uniform(0.3, 1.5, n)*rscale
    Y1 = rng.uniform(0.3, 1.5, n)*rscale; Y2 = rng.uniform(0.3, 1.5, n)*rscale
    s = rng.uniform(-0.6, 0.6, n); phi = rng.uniform(-0.9, 0.9, n); psi = rng.uniform(0.4, 2.7, n)
    return Pts(X1, X2, s, Y1, Y2, phi, psi)


def engine_pointwise(P, N, Kr, S_, nus, kind='EE2'):
    """f_eng[a, c, pt, nu] from bj pieces P (bj_ee.pieces output), exactly the integrand of bj_ee.images/value."""
    import bj_ee as BJ, bg_ser as S, bg_geom as GG
    pr = BJ.pref_sum(P, N, Kr)
    Bf = GG.bfac(Kr)
    Dl = P['Delta']; a = 1.0/Dl
    A_ = S_.atoms()
    Ssum = S_.X1 + S_.X2 + S_.Y1 + S_.Y2
    la = np.log(A_.astype(complex))
    na, nc = len(pr), len(pr[0])
    out = np.zeros((na, nc, S_.n, len(nus)), complex)
    for ia in range(na):
        for ic in range(nc):
            polys = [('vmu', P['vmu']), ('pref', pr[ia][ic])] + ([('R', P['R'])] if kind == 'EE4' else [])
            A = GG.assemble(P['logs'], polys, Kr)
            l0, l1 = A['lconst']
            for inu, nu in enumerate(nus):
                tot = np.zeros(S_.n, complex)
                for j, cj in A['corr'].items():
                    Pk = S.lmul(cj, Bf.c, Kr)
                    for k in range(Kr + 1):
                        if len(Pk[k]):
                            tot += (0.5j*nu)**j*a**k*Pk[k].ev(A_)
                lead = (la @ (A['E0'] + 1j*nu*A['E1']))
                out[ia, ic, :, inu] = np.exp(l0 + 1j*nu*l1 - pr[ia][ic].d*np.log(Dl) + lead - Ssum)*tot
    return out


def direct_pointwise(S_, Dl, J, nus):
    """f_dir[a, c, pt, nu] = exp((i nu/2) lA + (Delta/2) lB + G) pref (bf_core.grav_parts, EE k = 2, a^6 dropped)."""
    import bf_core as B
    Pp = B.grav_parts(S_, Dl, J, 'EE')
    base = 0.5*Dl*Pp['lB'] + Pp['G']
    out = np.zeros(Pp['pref'].shape + (len(nus),), complex)
    for inu, nu in enumerate(nus):
        out[..., inu] = Pp['pref']*np.exp(base + 0.5j*nu*Pp['lA'])[None, None, :]
    return out


def direct_pointwise_pw(S_, Dl, J, nus, chunk=1 << 10):
    """NOT USED in production: at Delta_ref = 2e4 it disagrees with the w0 contraction by O(1) (t - 1 ~ 1e-9: float64 cancellation in
    the TPS of t_U), at (336, 28) it agrees to 2e-9 (median).
    As direct_pointwise, but the spin-J contraction is pw_grav's exact closed form (f = P_J, mode 0) instead of the w0 Fourier sum
    over 2J+1 null-circle samples, which carries eps^2 4^J round-off at J >~ 40 (FOLD_WINDOW §1.1; HWIN_THIMBLE §1.2: the two agree to
    1e-13 where the Fourier sum is clean).  Same logs/branches as bf_core.grav_parts (n1_pw batch, continued from bf_core._ref_point)."""
    import types
    import bf_core as B, pw_grav as PG, pw_core as PC, l5_spin as S5
    PW = B.pw(); a = 1.0/Dl
    mons, sol, _ = S5.conserved(Dl, J)
    _, solt, _ = S5.conserved(3.0 - Dl, J)
    Z, lj, psi = S_.raw(a, 1)
    z0 = B._ref_point(a, S_.theta)
    ev = types.SimpleNamespace(kind='EE', nu=1.0, Delta=Dl, J=J, mons=mons, sol=sol, solt=solt, mode=0)
    S = PW.S
    k2 = S.kappa(2)
    out = None
    for c0 in range(0, S_.n, chunk):
        zz = np.vstack([z0[None, :], Z[c0:c0 + chunk]])
        PW.set_batch(zz.shape[0])
        zl = [zz[:, i].copy() for i in range(7)]
        _, _, P = PW.Cr.build(zl, 1.0, Dl, 0, n0=1.0, parts=True)
        lamF, Gs, Dm, geo = PG.kind_parts(ev, zl)
        zw, zzz = PC.zw_geom(zl)
        ts = zw*PW.powr(zzz, -0.5)
        x0, x1, x2, x3, x4, e, W2 = (geo[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'W2'))
        Ta, v = PG.side(S, x1, PW.Cr.sub(x1, x0), x2, e, W2, mons, sol)
        Tc, vp = PG.side(S, x3, PW.Cr.sub(x3, x0), x4, e, W2, mons, solt)
        tV = PG.gmet(v, vp)/np.sqrt(PG.gmet(v, v)*PG.gmet(vp, vp))
        sig = np.sign(np.real(tV/ts))
        P5, ln = PG.pairing(Ta, Tc, v, vp, J, sig*ts, 0)
        Ph = S.fourier(S.fourier(P5, 2, 2), 3, 2)[:, :, ::-1, ::-1]
        Dh = S.fourier(S.fourier(Dm, 0, 2), 1, 2)
        pref = np.einsum('i,j,acijn->acn', 1/k2, 1/k2, Ph*Dh[None, None])
        base = (0.5*Dl*P['lB'] + Gs + ln)[1:] + lj[c0:c0 + chunk]
        blk = np.zeros(pref.shape[:2] + (zz.shape[0] - 1, len(nus)), complex)
        for inu, nu in enumerate(nus):
            blk[..., inu] = pref[..., 1:]*np.exp(base + 0.5j*nu*P['lA'][1:])[None, None, :]
        out = blk if out is None else np.concatenate([out, blk], axis=2)
    return out


def run(J, Dst, depth, N, Kr, patch=False, npts=10, nus=(0.5, 3.0, 8.332091574, 16.664183148), pr=print, rscale=1.0):
    import bj_eeq as BQ
    if patch:
        import l5_spin as S5, fw_frame as FW, bf_core as B
        FW.patch_conserved_clean(S5)
        FW.install_pw(B.pw())                 # fw_points' cure of the eps^2 4^J Fourier round-off of the spin-J contraction (J >~ 40)
        FW.patch_conserved_clean(B.pw().S)
    nus = np.asarray(nus, float)
    S_ = make_points(npts, 1, rscale)
    res = {}
    for Dl, (dd, NN, KK) in ((Dst, (depth, N, Kr)), (DREF, (depth, min(N, 4), Kr))):
        t0 = time.time()
        import lj_rs_eeq as LE
        pp = LE.resolve(Dl, J, dd, patch, pieces=True) or 'missing'
        if not os.path.exists(pp):          # the dropped a-powers / branch constants are depth-independent: any depth's pieces will do
            for d2 in (3, 2, 4, 1):
                q = LE.resolve(Dl, J, d2, patch, pieces=True)
                if q:
                    pp, dd, KK = q, d2, min(KK, d2 - 1)
                    break
        P = pickle.load(open(pp, 'rb'))
        fe = engine_pointwise(P, NN, KK, S_, nus)
        fd = direct_pointwise_pw(S_, Dl, J, nus) if os.environ.get('LJ_DIR', 'w0') == 'pw' else direct_pointwise(S_, Dl, J, nus)
        if os.environ.get('LJ_DIRCHK') == '1':          # cross-check of the two direct contractions (clean only at J <~ 30)
            f2 = direct_pointwise(S_, Dl, J, nus)
            pr(f'#   direct pw vs w0 Fourier contraction: max |pw/w0 - 1| = {np.max(np.abs(fd/f2 - 1)):.2e}')
        r = fe/fd
        res[Dl] = dict(r=r, fe=fe, fd=fd)
        lr = np.log(r)
        c = np.exp(np.median(lr.real, axis=(0, 1, 2)) + 1j*np.median(lr.imag, axis=(0, 1, 2)))
        spread = np.max(np.abs(r/c[None, None, None, :] - 1), axis=(0, 1))     # per point, nu
        pr(f'# Delta={Dl:g} J={J}: pieces {os.path.basename(pp)} (N={NN}, Kr={KK}); {time.time() - t0:.0f}s')
        for inu, nu in enumerate(nus):
            pr(f'   nu={nu:8.4f}: c = {c[inu]:.10e}  max_pts,entries |r/c - 1| = {np.max(spread[:, inu]):.2e}  (median {np.median(spread[:, inu]):.1e})')
        res[Dl]['c'] = c
    a_st, a_ref = 1.0/Dst, 1.0/DREF
    rho = (a_st/a_ref)**6*res[DREF]['c']/res[Dst]['c']
    pr(f'# rho_an = (a_st/a_ref)^6 c_ref/c_st per nu: ' + '  '.join(f'{v:.10e}' for v in rho))
    pr(f'#   nu-spread of rho_an: max |rho/rho[0] - 1| = {np.max(np.abs(rho/rho[0] - 1)):.2e}')
    return dict(rho=rho, nus=nus, res={D: dict(c=v['c'], r=v['r']) for D, v in res.items()})


if __name__ == '__main__':
    J, Dst, depth, N, Kr = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    patch = len(sys.argv) > 6 and sys.argv[6] == 'p'
    lines = []

    def pr(s=''):
        print(s, flush=True); lines.append(s)
    out = run(J, Dst, depth, N, Kr, patch, pr=pr)
    for rs in (0.3, 3.0):
        pr(f'# radial scale x{rs} (truncation check)')
        o2 = run(J, Dst, depth, N, Kr, patch, pr=pr, rscale=rs)
        out['rho_rs%g' % rs] = o2['rho']
    tag = 'lj_rs_norm_J%d_D%g_d%d_N%d_K%d%s%s' % (J, Dst, depth, N, Kr, '_p' if patch else '', os.environ.get('LJ_NTAG', ''))
    open(os.path.join(C, tag + '.txt'), 'w').write('\n'.join(lines) + '\n')
    pickle.dump(out, open(os.path.join(C, tag + '.pkl'), 'wb'))

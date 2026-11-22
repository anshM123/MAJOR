"""
lj_rs_chain (HLJ_CLOSE.md §4): the (336, 28/32) nu < 8.3 window through the UNCHANGED hc_chain z4 protocol (HYPOTHESES_CLOSE §5), with the
closed-form engine's C2 values below nu = 8.3 normalised ANALYTICALLY (lj_rs_norm rho_an) instead of bj_chain.calibrate's fitted constant.
Same two-pass wrapper logic as bj_chain.chain2 / chain_win (imported functions are reused where they write nothing); every file written
here is cache/lj_rs_* except the unchanged chain's own output cache/hc_chain_z4_<HC_TAG>.* (HC_TAG always starts with 'lj_rs').
  pass 1:  python3 lj_rs_chain.py nl  J depth N Kr TAG [rho_shift_c2]   -> cache/lj_rs_recs_J<J>_<TAG>.pkl, cache/lj_rs_nl_J<J>_<TAG>.pkl
  pass 2:  BJ_C3H=keep|drop python3 lj_rs_chain.py win J TAG             -> cache/hc_chain_z4_<HC_TAG>.txt/.pkl   (HC_TAG = lj_rs<J>_<TAG>_<c3h>)
rho_shift_c2 (stress): rho -> rho (1 + shift/m^2), i.e. a uniform c0-type shift of the engine values by `shift` c2.
The reference constant per replica (kr_r) is still taken from bj_chain.calibrate: the reference engine equals the direct reference
records entrywise to 2e-7 ... 5e-5 (BELOW_FOLD_SPINJ §2.1), so that constant carries no engine information; only rho is replaced.
"""
import os, sys, pickle, runpy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')


def nl(J, depth, N, Kr, tag, shift=0.0):
    import bj_chain as BC, z3_nodes as Z3
    import lj_rs_val as LV
    BC.engine = lambda Dl, J_, d_, N_, K_: LV.eng(Dl, J_, d_, N_, K_)     # bj_ or lj_rs_bj_ images (LE.resolve)
    cal = BC.calibrate(J, depth=depth, N=N, Kr=Kr, pr=print)
    rho_fit = cal['rho']
    rho = LV.rho_an(J, depth, N, Kr)*(1.0 + shift/cal['m']**2)
    m2 = cal['m']**2
    print(f'# lj_rs_chain: rho_an = {rho:.8e}; fitted rho (bj_chain, cal nodes) = {rho_fit:.8e}; '
          f'm^2 (rho_fit/rho_an - 1) = {m2*(rho_fit/rho - 1):+.1f} c2', flush=True)
    fext = float(os.environ.get('LJ_FEXT', '0'))         # engine tail stress: I -> I_K + fext (I_K - I_{K-1})  (same depth/N)
    if fext:
        nus1, I1 = LV.eng(BC.DST, J, depth, N, Kr - 1)
        assert np.allclose(nus1, cal['nus'])
        cal['Ist'] = cal['Ist'] + fext*(cal['Ist'] - I1)
        print(f'# lj_rs_chain: engine tail stress fext = {fext:+g} x (order {Kr} - order {Kr - 1})', flush=True)
    sst = os.environ.get('LJ_SSTRUCT', '')                # structure stress: I_ac -> I_ac (1 + d_ac/m^2), d = (d01, d10, d11) c2
    if sst:
        d01, d10, d11 = (float(x) for x in sst.split(','))
        fac = 1.0 + np.array([[0.0, d01], [d10, d11]])/m2
        cal['Ist'] = cal['Ist']*fac[None, :, :]
        print(f'# lj_rs_chain: structure stress (d01, d10, d11) = ({d01:+g}, {d10:+g}, {d11:+g}) c2 on every engine node', flush=True)
    cal['rho_r'] = np.full(cal['R'], rho, complex)
    cal['rho'] = rho
    recs = BC.records(J, cal)
    pickle.dump(recs, open(os.path.join(C, 'lj_rs_recs_J%d_%s.pkl' % (J, tag)), 'wb'))
    orig = Z3.zpts
    nucut = BC.NUCUT

    def zpts1(kind, J_):
        P = dict(orig(kind, J_))
        if J_ != J:
            return P
        for k in [k for k in P if k[1] < nucut and k[0] in (BC.DST, BC.DREF)]:
            P.pop(k)
        P.update(recs)
        return P
    Z3.zpts = zpts1
    os.environ['HC_TAG'] = 'lj_rs%d_%s_pass1' % (J, tag)
    sys.argv = ['hc_chain.py', 'z4', '300']
    runpy.run_path(os.path.join(HERE, 'hc_chain.py'), run_name='__main__')
    Z4 = sys.modules['z4_chain']
    rec = Z4.NLROWS[(J, BC.DST)]
    sel = rec['nus'] < nucut
    store = dict(nus=rec['nus'][sel].copy(), **{h: rec[h][sel].copy() for h in ('++', '+-')},
                 c={h: rec[(h, 'c')][sel].copy() for h in ('++', '+-')}, rho=rho, rho_fit=rho_fit)
    pickle.dump(store, open(os.path.join(C, 'lj_rs_nl_J%d_%s.pkl' % (J, tag)), 'wb'))


def win(J, tag, c3h='keep'):
    import bj_chain as BC, z3_nodes as Z3, w2_chain as W2
    from scipy.interpolate import CubicSpline
    store = pickle.load(open(os.path.join(C, 'lj_rs_nl_J%d_%s.pkl' % (J, tag)), 'rb'))
    recs = pickle.load(open(os.path.join(C, 'lj_rs_recs_J%d_%s.pkl' % (J, tag)), 'rb'))
    nucut = BC.NUCUT
    orig2 = Z3.zpts

    def zpts2(kind, J_):
        P = dict(orig2(kind, J_))
        if J_ != J:
            return P
        low = [k for k in P if k[1] < nucut and k[0] in (BC.DST, BC.DREF)]
        if c3h == 'drop':
            for k in low:
                P.pop(k)
        elif kind == 'EE':
            for k in low:
                P[k] = recs[k]
        return P
    Z3.zpts = zpts2
    _sm = W2.set_mode

    def set_mode(*a, **k):
        if 'z4_chain' not in sys.modules or getattr(sys.modules['z4_chain'], '_lj_patched', False):
            return _sm(*a, **k)
        Z4b = sys.modules['z4_chain']
        Z4b._lj_patched = True
        nlo = Z4b._nl_row
        nn = store['nus']; o = np.argsort(nn)
        sp = {h: (CubicSpline(nn[o], store[h][o].real), CubicSpline(nn[o], store[h][o].imag)) for h in ('++', '+-')}

        def _nl_row(J_, D, h, nu):
            v = nlo(J_, D, h, nu)
            if J_ != J or float(D) != BC.DST:
                return v
            nu = np.asarray(nu, float)
            nc = np.clip(nu, nn.min(), None)
            e = sp[h][0](nc) + 1j*sp[h][1](nc)
            gs = float(os.environ.get('LJ_GSHIFT', '0'))
            if gs:
                hm, hp = Z4b.hhat(nu)
                e = e + gs*(hm + hp)
            return np.where(nu < nucut, e, v)
        Z4b._nl_row = _nl_row
        return _sm(*a, **k)
    W2.set_mode = set_mode
    os.environ['HC_TAG'] = 'lj_rs%d_%s_%s%s' % (J, tag, c3h, os.environ.get('LJ_TAGSFX', ''))
    sys.argv = ['hc_chain.py', 'z4', '300']
    runpy.run_path(os.path.join(HERE, 'hc_chain.py'), run_name='__main__')


if __name__ == '__main__':
    mode, J = sys.argv[1], int(sys.argv[2])
    if mode == 'nl':
        depth, N, Kr = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        tag = sys.argv[6]
        shift = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
        nl(J, depth, N, Kr, tag, shift)
    elif mode == 'win':
        win(J, sys.argv[3], os.environ.get('BJ_C3H', 'keep'))

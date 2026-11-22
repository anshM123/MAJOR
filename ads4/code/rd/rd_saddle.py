"""
rd_saddle: crystal saddles of the scalar flow function Phi_w (= n1_direct.phi_w, J = 0, n0 = +1), re-computed with rd_core's exact
gradient: the symmetric crystal S above the fold, and the sigma-conjugate pair A+, A- below it (the n4_thimble construction:
continue S from nu = 20 down to nu_c = 10 with sigma-symmetrisation, kick, Newton, then continue A+- in nu with steps <= 0.25).

track(Delta, nus) -> {nu: dict(kind, c, w, H, L, y0, grad)} for the requested nus (nu >= 11.6: 'S'; nu <= 10: 'A+', 'A-').
Cache: cache/rd_saddles_D<Delta>.pkl
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '2')
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'spinfunc'))
import numpy as np
import rd_core as C
import l1_tps as T            # read-only: csym_cholesky


def sigma(w):
    return np.array([0.0, np.conj(w[3]), np.conj(w[4]), np.conj(w[1]), np.conj(w[2]), 0.0, 0.0], complex)


def asym(w):
    return abs(w[1] - np.conj(w[3])) + abs(w[2] - np.conj(w[4]))


def y0_of(w, c, ref=None):
    z, _, _ = C.zmap(w[:, None], c)
    return C.geom_plain(z, None if ref is None else np.array([ref]))['y0'][0]


def newton(w, nu, Delta, y0, it=100, cap=0.1, tol=1e-12):
    c, _ = C.crystal(nu, Delta)
    F = C.Flow(nu, Delta, c)
    for k in range(it):
        y0 = y0_of(w, c, y0)
        g, _ = F.grad(w[:, None], np.array([y0]))
        H = F.hess(w, y0)
        dw = np.linalg.solve(H, g[:, 0])
        s = 1.0 if np.abs(dw).max() <= cap else cap/np.abs(dw).max()
        w = w - s*dw
        if not np.all(np.isfinite(w)):
            raise RuntimeError('newton diverged')
        if np.abs(dw).max() < tol:
            break
    y0 = y0_of(w, c, y0)
    g, _ = F.grad(w[:, None], np.array([y0]))
    H = F.hess(w, y0)
    return c, w, H, y0, float(np.abs(g).max())


def wshift(w, c_old, c_new):
    w = w.copy()
    for i in C.LOGV:
        w[i] = w[i] + np.log(c_old[i]/c_new[i])
    return w


def rec(kind, nu, c, w, H, y0, gm):
    L = np.linalg.inv(T.csym_cholesky(-H)).T
    return dict(kind=kind, nu=nu, c=c, w=w, H=H, L=L, y0=y0, grad=gm)


def track(Delta, nus, nu_c=10.0, kick=1e-3, verbose=True):
    fn = os.path.join(HERE, 'cache', f'rd_saddles_D{Delta:g}.pkl')
    store = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    want = [v for v in nus if v not in store]
    if not want:
        return {v: store[v] for v in nus}
    # symmetric branch from 20 down
    top = sorted([v for v in want if v > nu_c + 1.0], reverse=True)
    path = sorted(set(list(np.arange(20.0, nu_c - 1e-9, -1.0)) + top + [nu_c]), reverse=True)
    w = np.zeros(7, complex); cprev = None; y0 = None
    for v in path:
        cn, _ = C.crystal(v, Delta)
        if cprev is not None:
            w = wshift(w, cprev, cn)
        if y0 is None:
            y0 = y0_of(w, cn)
        if v > nu_c:
            c, w, H, y0, gm = newton(w, v, Delta, y0)
            w = 0.5*(w + sigma(w)); cprev = cn
            if v in want:
                store[v] = {'S': rec('S', v, c, w, H, y0, gm)}
                if verbose:
                    print(f'  Delta {Delta:g} nu {v:6.2f} S: |grad| {gm:.1e}', flush=True)
        else:
            wS = w; cprev = cn
    # A+- at nu_c
    low = sorted([v for v in want if v <= nu_c], reverse=True)
    if low:
        out = {}
        for which, sg in (('A+', 1), ('A-', -1)):
            w0 = wS.copy(); w0[1] += sg*kick*(1 + 1j); w0[3] -= sg*kick*(1 - 1j)
            c, w, H, y0w, gm = newton(w0, nu_c, Delta, y0)
            if asym(w) < 1e-4:
                raise RuntimeError('returned to the symmetric saddle')
            out[which] = (c, w, H, y0w, gm)
        if (out['A+'][1][1] - np.conj(out['A+'][1][3])).imag < 0:
            out['A+'], out['A-'] = out['A-'], out['A+']
        for which in ('A+', 'A-'):
            c, w, H, y0w, gm = out[which]
            cur = nu_c
            if nu_c in want:
                store.setdefault(nu_c, {})[which] = rec(which, nu_c, c, w, H, y0w, gm)
            for v in low:
                if v >= nu_c:
                    continue
                nsub = max(1, int(np.ceil((cur - v)/0.25)))
                for k in range(1, nsub + 1):
                    vv = cur + (v - cur)*k/nsub
                    cn, _ = C.crystal(vv, Delta)
                    w = wshift(w, c, cn)
                    c, w, H, y0w, gm = newton(w, vv, Delta, y0w)
                cur = v
                store.setdefault(v, {})[which] = rec(which, v, c, w, H, y0w, gm)
                if verbose:
                    print(f'  Delta {Delta:g} nu {v:6.2f} {which}: |grad| {gm:.1e}, asym {asym(w):.3f}', flush=True)
    pickle.dump(store, open(fn, 'wb'))
    return {v: store[v] for v in nus}


if __name__ == '__main__':
    Dl = float(sys.argv[1])
    nus = [float(v) for v in sys.argv[2].split(',')]
    R = track(Dl, nus)
    for v in nus:
        for k, r in R[v].items():
            z, _, _ = C.zmap(r['w'][:, None], r['c'])
            print(f'nu {v:6.2f} {k}: w = {np.round(r["w"], 4)}; raw z = {np.round(z[:, 0]*Dl, 4)} (x Delta)')

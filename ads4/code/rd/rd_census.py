"""
rd_census (RD_CLOSE.md section 4.3): census of critical points of the scalar flow function Phi_w AWAY from the crystal (raw radius
r >= 0.05), i.e. the saddles whose thimbles could enter a Lefschetz decomposition of the physical cycle besides A+-/S.

Damped Newton (exact gradient, holomorphic FD Hessian) from random starts in the far region (w-coordinates: |raw u, v| in [0.02, 3],
random phases; x~0, x_p, y_p random complex O(1) / O(0.5)), vectorised over starts.  Converged points (|grad| < 1e-9) are de-duplicated
and characterised by: raw radius r, (Delta/2) log|B| (sheet-independent part of Re lamF), |A|, the smallest normalised divisor
distance (which divisor), and log(|Q(sigma)|/Q_eff) with Q_eff of the crystal thimble (rd_flow).
Re Phi itself is defined only up to the sheet of A^{i nu/2} (shifts by pi nu k); the B part is not.
Usage: python3 rd_census.py Delta nu [nstart]  -> cache/rd_census_D<Delta>_nu<nu>.txt
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import rd_core as C
import rd_saddle as SD
np.seterr(all='ignore')


def hess_batch(F, w, y0, h=1e-5):
    om = np.array([1, 1j, -1, -1j]); N = w.shape[1]
    W = np.repeat(w[:, :, None], 28, axis=2)            # (7, N, 28)
    for i in range(7):
        for j in range(4):
            W[i, :, 4*i + j] += h*om[j]*np.maximum(1.0, np.abs(w[i]))
    G, _ = F.grad(W.reshape(7, -1), np.repeat(y0, 28))
    G = G.reshape(7, N, 28)
    H = np.zeros((N, 7, 7), complex)
    for i in range(7):
        hi = h*np.maximum(1.0, np.abs(w[i]))
        H[:, :, i] = (sum(G[:, :, 4*i + j]/om[j] for j in range(4))/(4*hi[None, :])).T
    return 0.5*(H + np.transpose(H, (0, 2, 1)))


def run(Delta, nu, nstart=4000, seed=3, iters=80):
    rec = SD.track(Delta, [nu])[nu]
    c, _ = C.crystal(nu, Delta)
    F = C.Flow(nu, Delta, c)
    rng = np.random.default_rng(seed)
    w = np.zeros((7, nstart), complex)
    for i in C.LOGV:
        rad = np.exp(rng.uniform(np.log(0.02), np.log(3.0), nstart))
        w[i] = np.log(rad/np.abs(c[i])) + 1j*rng.uniform(-np.pi, np.pi, nstart)
    w[0] = 0.7*(rng.normal(size=nstart) + 1j*rng.normal(size=nstart))
    w[5] = 0.4*(rng.normal(size=nstart) + 1j*rng.normal(size=nstart)); w[6] = 0.4*(rng.normal(size=nstart) + 1j*rng.normal(size=nstart))
    z, _, _ = C.zmap(w, c); y0 = C.geom_plain(z)['y0']
    act = np.ones(nstart, bool); gn = np.full(nstart, np.inf)
    for it in range(iters):
        idx = np.where(act)[0]
        if idx.size == 0:
            break
        wa = w[:, idx]; ya = y0[idx]
        g, _ = F.grad(wa, ya)
        H = hess_batch(F, wa, ya)
        try:
            dw = np.linalg.solve(H, g.T[:, :, None])[:, :, 0].T
        except np.linalg.LinAlgError:
            dw = np.zeros_like(g)
        mx = np.max(np.abs(dw), axis=0)
        sc = np.where(mx > 0.5, 0.5/mx, 1.0)
        wn = wa - sc*dw
        z, _, _ = C.zmap(wn, c); yn = C.geom_plain(z, ya)['y0']
        w[:, idx] = wn; y0[idx] = yn
        gn[idx] = np.max(np.abs(g), axis=0)
        dead = ~np.all(np.isfinite(wn), axis=0) | (np.max(np.abs(wn.real), axis=0) > 40)
        conv = (mx < 1e-11)
        act[idx[dead | conv]] = False
    g, vals = F.grad(w, y0)
    gn = np.max(np.abs(g), axis=0)
    ok = np.isfinite(gn) & (gn < 1e-9*np.maximum(1, np.max(np.abs(w), axis=0)))
    W = w[:, ok]; Y = y0[ok]
    # dedupe in raw coordinates
    z, _, _ = C.zmap(W, c)
    keys = []; uniq = []
    for j in range(W.shape[1]):
        zz = z[:, j]
        if not any(np.max(np.abs(zz - z[:, u])/(np.abs(z[:, u]) + 1e-3)) < 1e-6 for u in uniq):
            uniq.append(j)
    z = z[:, uniq]; Y = Y[uniq]; W = W[:, uniq]
    gv = C.geom_plain(z, Y)
    vv = {k: v for k, v in gv.items() if not isinstance(v, list)}
    lB = 0.5*Delta*np.real(F.logB(vv))
    _, ljW, _ = C.zmap(W, c)
    Hs = F.logabs_parts(vv) + np.real(ljW)                       # sheet-independent part of Re Phi (arg A term excluded)
    argA = np.imag(F.logA(vv))
    lAabs = np.real(F.logA(vv))
    r = np.max(np.abs(z[1:5]), axis=0)
    # crystal reference
    zc, _, _ = C.zmap(np.stack([rr['w'] for rr in rec.values()], axis=1), c)
    gc = C.geom_plain(zc, np.array([rr['y0'] for rr in rec.values()]))
    vc = {k: v for k, v in gc.items() if not isinstance(v, list)}
    lBc = 0.5*Delta*np.real(F.logB(vc))
    Wc = np.stack([rr['w'] for rr in rec.values()], axis=1)
    _, ljc, _ = C.zmap(Wc, c)
    Hc = F.logabs_parts(vc) + np.real(ljc)
    argAc = np.imag(F.logA(vc))
    Qc = np.abs(C.Qratio(gc, Delta))
    Q = np.abs(C.Qratio(gv, Delta))
    fl = os.path.join(HERE, 'cache', f'rd_flow_D{Delta:g}_nu{nu:g}.pkl')
    qeff = None
    if os.path.exists(fl):
        R = pickle.load(open(fl, 'rb'))
        qeff = np.mean([R[k]['qnear']*R[k]['Qs'] for k in R])
    else:
        qeff = float(np.mean(Qc))
    a = 1.0/Delta
    dv = np.stack([np.abs(1 - (z[1] + z[2]) + (z[1]*z[2] - z[5]**2))/(1 + np.abs(z[1]) + np.abs(z[2])),
                   np.abs(1 - (z[3] + z[4]) + (z[3]*z[4] - z[6]**2))/(1 + np.abs(z[3]) + np.abs(z[4])),
                   np.abs(1 + (z[1] + z[2]) + (z[1]*z[2] - z[5]**2))/(1 + np.abs(z[1]) + np.abs(z[2])),
                   np.abs(1 + (z[3] + z[4]) + (z[3]*z[4] - z[6]**2))/(1 + np.abs(z[3]) + np.abs(z[4])),
                   np.abs(gv['FP1']), np.abs(gv['FP2']),
                   np.abs(1.0 - z[5]**2/(z[1]*z[2])), np.abs(1.0 - z[6]**2/(z[3]*z[4])),
                   np.abs(z[1])/a, np.abs(z[2])/a, np.abs(z[3])/a, np.abs(z[4])/a])
    lab = ['m45=0', 'm35=0', 'm25=0', 'm15=0', 'FP1=0', 'FP2=0', 'x^2=0', 'y^2=0', 'X1=0', 'X2=0', 'Y1=0', 'Y2=0']
    out = dict(Hs=Hs, Hc=Hc, argA=argA, argAc=argAc, z=z, W=W, lB=lB, lBc=lBc, r=r, lAabs=lAabs, lQ=np.log(Q/qeff), dv=dv, nconv=int(ok.sum()), nstart=nstart)
    L = [f'# rd_census: Delta = {Delta:g}, nu = {nu:g}: {nstart} Newton starts, {int(ok.sum())} converged, {len(uniq)} distinct critical points',
         f'#   crystal saddle(s): (Delta/2) log|B| = {np.round(lBc, 3)}, raw radius {np.round(np.max(np.abs(zc[1:5]), axis=0), 4)}',
         f'#   dH = [sum Re(e_k) log|f_k| + Re logJac](sigma) - same at the crystal (the height without the sheet term -(nu/2)(arg A + 2 pi k));',
         f'#   sheet term range: -(nu/2) x [principal arg A difference] shown; each further sheet shifts by pi nu = {np.pi*nu:.2f}',
         '#   r_raw   (Delta/2)log|B|   dH       -(nu/2)dargA   log|A|   log(|Q|/Q_eff)  nearest divisor (normalised distance)']
    order = np.argsort(-(Hs - Hc[0]))
    for j in order:
        k = int(np.argmin(dv[:, j]))
        L.append(f'   {r[j]:9.4f}   {lB[j]:10.3f}   {Hs[j] - Hc[0]:9.2f}   {-(nu/2)*(argA[j] - argAc[0]):+8.2f}   {lAabs[j]:7.2f}   {out["lQ"][j]:+10.2f}     {lab[k]} ({dv[k, j]:.1e})')
    return out, L


if __name__ == '__main__':
    Dl = float(sys.argv[1]); nu = float(sys.argv[2]); ns = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    t = time.time()
    out, L = run(Dl, nu, ns)
    L.append(f'# {time.time() - t:.0f}s')
    print('\n'.join(L))
    tag = f'D{Dl:g}_nu{nu:g}'
    open(os.path.join(HERE, 'cache', f'rd_census_{tag}.txt'), 'w').write('\n'.join(L) + '\n')
    pickle.dump(out, open(os.path.join(HERE, 'cache', f'rd_census_{tag}.pkl'), 'wb'))

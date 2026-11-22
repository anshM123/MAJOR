"""
rd_flow (RD_CLOSE.md section 3): the far-region weight ON THE CRYSTAL DESCENT CYCLE, measured along its own flow lines.

The crystal cycle used by the project below the fold is n4's  n+ J(A+) + n- J(A-)  (Lefschetz thimbles of the scalar flow function
Phi_w, n0 = +1); above the fold the thimble J(S) of the symmetric crystal (n3).  A thimble is the union of the downward flow lines
d xi/dt = -conj(dPhi/dxi) (xi: the saddle's normal coordinates, w = w_s + L xi) leaving the saddle; along each line Im Phi is constant
and h = Re Phi decreases.  Parametrising by the action drop  s = h_s - h  (dw/ds = -L conj(L^T g)/|L^T g|^2, RK4), the scalar weight
of a thimble point relative to the saddle is exactly e^{-s}; the EE weight is e^{-s} |Q(z)|/Q_eff.

Output per (Delta, nu, saddle): for raw radii r (r = max(|u_x|,|v_x|,|u_y|,|v_y|), the light-cone components of x and y):
   c(r)    = min over the sampled thimble points with radius >= r of s/Delta            (scalar far-region exponent)
   cEE(r)  = min over the same points of [s - log(|Q|/Q_eff)]/Delta                      (EE far-region exponent)
   and the extreme values of log(|Q|/Q_eff), |w_x| = |1 - x_p^2/(u_x v_x)| (s -> +-1), |m15|, |m25|, |FP1| (m12 = 0) met.
Q_eff = <|Q|> over the near thimble with the Gaussian weight e^{-s} s^{5/2} ds (7-dim radial measure).
Directions: N random unit vectors in the 7-dim real tangent space, then local refinement of the best lines at each r.

Usage: python3 rd_flow.py Delta nu [N] [smax_over_Delta]  -> cache/rd_flow_D<Delta>_nu<nu>.{txt,pkl}
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import rd_core as C
import rd_saddle as SD
np.seterr(all='ignore')

RADII = (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
R0 = 0.05
DWMAX = float(os.environ.get('RD_DWMAX', '0.03'))
NSTEP = int(os.environ.get('RD_NSTEP', '2500'))
WX_TRUST = 0.03


def _subset(v, m):
    if isinstance(v, list):
        return [_subset(u, m) for u in v]
    if isinstance(v, np.ndarray) and v.ndim == 1 and v.size == m.size:
        return v[m]
    return v


def flow(rec, Delta, dirs, smax, rstop=1.3, tol=2e-3):
    """integrate flow lines from the saddle record rec along unit directions dirs (N, 7). Returns per-line summaries.
    Step control: an RK4 step in the action parametrisation is accepted only if the MEASURED change of Phi (sum_k e_k Log(f_k'/f_k)
    + dlogJac, principal logs of small ratios) satisfies |dPhi + ds| <= tol ds and |Im dPhi| <= tol ds; s accumulates the measured
    -Re dPhi, so s is the exact action drop (h_s - h) up to the accumulated tolerance."""
    nu = rec['nu']; c = rec['c']; L = rec['L']; ws = rec['w']
    F = C.Flow(nu, Delta, c)
    N = dirs.shape[0]
    w = ws[:, None] + L @ (R0*dirs.T).astype(complex)
    y0 = np.full(N, rec['y0'], complex)
    z, _, _ = C.zmap(w, c); y0 = C.geom_plain(z, y0)['y0']
    s = np.full(N, 0.5*R0*R0)
    zs, _, _ = C.zmap(ws[:, None], c)
    gs = C.geom_plain(zs, np.array([rec['y0']]))
    Qs = abs(C.Qratio(gs, Delta)[0])
    nr = len(RADII)
    sfirst = np.full((N, nr), np.inf); cee = np.full((N, nr), np.inf); qmax = np.full((N, nr), -np.inf)
    minwx = np.full(N, np.inf); min15 = np.full(N, np.inf); min25 = np.full(N, np.inf); minFP1 = np.full(N, np.inf)
    rmax = np.zeros(N); sfin = np.zeros(N)
    qw_num = np.zeros(N); qw_den = np.zeros(N)
    last = np.full((9, N), np.nan)
    end = np.zeros(N, int)                     # 1 smax, 2 rstop, 3 step cap, 4 non-finite, 5 step underflow
    keys = [k for k, _ in F.fac]; ex = np.array([e for _, e in F.fac])
    imacc = np.zeros(N); nrej = np.zeros(N, int)
    lqlast = np.zeros(N); nfroz = np.zeros(N, int)
    SCK = 10.0*2.0**np.arange(9)                      # action checkpoints 10, 20, ..., 2560
    sck = np.full((N, SCK.size), np.nan); lqck = np.full((N, SCK.size), np.nan)
    hstep = np.full(N, 1e-3)

    def rhs(wa, ya):
        g, _ = F.grad(wa, ya)
        gx = L.T @ g
        nrm = np.sum(np.abs(gx)**2, axis=0)
        return -(L @ np.conj(gx))/nrm

    def measure(wa, ya):
        z, lj, _ = C.zmap(wa, c)
        gv = C.geom_plain(z, ya)
        return z, lj, gv, np.stack([gv[k] for k in keys])

    def record(idx, z, gv, sa, ds):
        Q = np.abs(C.Qratio(gv, Delta))/Qs
        ok = np.isfinite(Q) & (Q > 0) & np.all(np.isfinite(z), axis=0)
        # double-precision contraction is unreliable for |w_x| < WX_TRUST (rd_check_q): freeze log|Q| at the line's last trusted value
        wxa = np.abs(1.0 - z[5]**2/(z[1]*z[2]))
        trust = ok & (wxa > WX_TRUST)
        lq = np.where(trust, np.log(np.where(ok, Q, 1.0)), lqlast[idx])
        lqlast[idx] = lq
        nfroz[idx] += (~trust & ok)
        Q = np.exp(lq); z = np.where(ok[None, :], z, 0.0)
        for k, sk in enumerate(SCK):
            hitk = (sa >= sk) & np.isnan(sck[idx, k]) & ok
            sck[idx[hitk], k] = sa[hitk]; lqck[idx[hitk], k] = lq[hitk]
        r = np.maximum(np.maximum(np.abs(z[1]), np.abs(z[2])), np.maximum(np.abs(z[3]), np.abs(z[4])))
        r = np.where(ok, r, 0.0)
        rmax[idx] = np.maximum(rmax[idx], r)
        for k, rk in enumerate(RADII):
            beyond = (r >= rk) & ok
            hit = beyond & ~np.isfinite(sfirst[idx, k])
            sfirst[idx[hit], k] = sa[hit]
            cee[idx[beyond], k] = np.minimum(cee[idx[beyond], k], sa[beyond] - lq[beyond])
            qmax[idx[beyond], k] = np.maximum(qmax[idx[beyond], k], lq[beyond])
        wx = np.where(ok, np.abs(1.0 - z[5]**2/(z[1]*z[2])), np.inf)
        minwx[idx] = np.minimum(minwx[idx], wx)
        min15[idx] = np.minimum(min15[idx], np.where(ok, np.abs(gv['m15']), np.inf))
        min25[idx] = np.minimum(min25[idx], np.where(ok, np.abs(gv['m25']), np.inf))
        minFP1[idx] = np.minimum(minFP1[idx], np.where(ok, np.abs(gv['FP1']), np.inf))
        a = 1.0/Delta
        dv = np.stack([np.abs(1 - (z[1] + z[2]) + (z[1]*z[2] - z[5]**2))/(1 + np.abs(z[1]) + np.abs(z[2])),      # m45 zero (B_x = 0)
                       np.abs(1 - (z[3] + z[4]) + (z[3]*z[4] - z[6]**2))/(1 + np.abs(z[3]) + np.abs(z[4])),      # m35 zero (B_y = 0)
                       np.abs(z[1])/a, np.abs(z[2])/a, np.abs(z[3])/a, np.abs(z[4])/a,                         # X1, X2, Y1, Y2 -> 0
                       wx, np.abs(1.0 - z[6]**2/(z[3]*z[4])),                                                  # s -> +-1 ; cos(phi) -> 0
                       np.abs(gv['y0'])/(a*np.sqrt(np.abs(z[3]*z[4])) + 1e-300)])                              # sin(psi) -> 0
        last[:, idx[ok]] = dv[:, ok]
        wt = np.exp(-sa)*sa**2.5*ds
        qw_num[idx] += np.where(ok, wt*Q, 0.0); qw_den[idx] += np.where(ok, wt, 0.0)
        return r

    z, lj, gv, fv = measure(w, y0)
    record(np.arange(N), z, gv, s, np.full(N, 0.5*R0*R0))
    F_prev = fv; lj_prev = lj
    act = np.arange(N)
    for step in range(NSTEP):
        if act.size == 0:
            break
        wa, ya, sa, ds = w[:, act], y0[act], s[act], hstep[act]
        k1 = rhs(wa, ya)
        k2 = rhs(wa + 0.5*ds*k1, ya)
        k3 = rhs(wa + 0.5*ds*k2, ya)
        k4 = rhs(wa + ds*k3, ya)
        wn = wa + ds*(k1 + 2*k2 + 2*k3 + k4)/6
        z, lj, gv, fv = measure(wn, ya)
        dP = ex @ np.log(fv/F_prev[:, act]) + (lj - lj_prev[act])
        err = np.maximum(np.abs(dP.real + ds), np.abs(dP.imag))
        fin = np.isfinite(dP) & np.all(np.isfinite(wn), axis=0)
        acc = fin & (err <= tol*ds)
        # rejected: shrink
        hstep[act[~acc]] *= 0.25
        ia = act[acc]
        if ia.size:
            w[:, ia] = wn[:, acc]; y0[ia] = gv['y0'][acc]
            s[ia] = sa[acc] - dP[acc].real; imacc[ia] += dP[acc].imag
            F_prev[:, ia] = fv[:, acc]; lj_prev[ia] = lj[acc]
            sub = {k: _subset(v, acc) for k, v in gv.items()}
            r = record(ia, z[:, acc], sub, s[ia], ds[acc])
            hstep[ia] = np.minimum(ds[acc]*1.6, np.maximum(0.05, 0.05*s[ia]))
            e = np.zeros(ia.size, int)
            e[s[ia] > smax] = 1; e[r > rstop] = 2
            end[ia[e > 0]] = e[e > 0]; sfin[ia[e > 0]] = s[ia[e > 0]]
        nrej[act[~acc]] += 1
        e2 = (hstep[act] < 1e-12)
        end[act[e2]] = 5; sfin[act[e2]] = s[act[e2]]
        act = act[(end[act] == 0)]
    end[act] = 3; sfin[act] = s[act]
    return dict(sfirst=sfirst, cee=cee, qmax=qmax, minwx=minwx, min15=min15, min25=min25, minFP1=minFP1, rmax=rmax, sfin=sfin,
                end=end, qnear=(np.nansum(qw_num)/np.nansum(qw_den)), Qs=Qs, last=last, imerr=np.abs(imacc), nrej=nrej, nfroz=nfroz, sck=sck, lqck=lqck)


def merge(a, b):
    out = {}
    for k in a:
        if isinstance(a[k], np.ndarray) and a[k].ndim >= 1:
            out[k] = np.concatenate([a[k], b[k]], axis=(1 if k == 'last' else 0))
        else:
            out[k] = a[k]
    return out


def run(Delta, nu, N=1024, smax_fac=2.0, seed=7, refine=2):
    R = SD.track(Delta, [nu])[nu]
    rng = np.random.default_rng(seed)
    res = {}
    for kind, rec in R.items():
        t = time.time()
        d = rng.normal(size=(N, 7)); d /= np.linalg.norm(d, axis=1)[:, None]
        out = flow(rec, Delta, d, smax_fac*Delta); dirs = d
        for it in range(refine):                 # local refinement around the best line of every radius (scalar and EE)
            seeds = []
            for k in range(len(RADII)):
                for arr in (out['sfirst'][:, k], out['cee'][:, k]):
                    if np.isfinite(arr).any():
                        seeds.extend(list(np.argsort(arr)[:3]))
            seeds = sorted(set(seeds))
            sig = 0.15/(it + 1)
            nd = []
            for j in seeds:
                p = dirs[j][None, :] + sig*rng.normal(size=(16, 7))
                nd.append(p/np.linalg.norm(p, axis=1)[:, None])
            nd = np.concatenate(nd)
            o2 = flow(rec, Delta, nd, smax_fac*Delta)
            out = merge(out, o2); dirs = np.concatenate([dirs, nd])
        out['dirs'] = dirs; out['time'] = time.time() - t
        res[kind] = out
    return res


def report(Delta, nu, res):
    lines = [f'# rd_flow: Delta = {Delta:g}, nu = {nu:g}; far-region exponents on the crystal thimble(s); radius r = max |light-cone comp.|']
    for kind, o in res.items():
        Qeff = o['qnear']
        lines.append(f'## saddle {kind}: lines {len(o["end"])}, end codes (smax/rstop/cap/bad/underflow) {[int((o["end"] == k).sum()) for k in (1, 2, 3, 4, 5)]}, '
                     f'max raw radius reached {o["rmax"].max():.3f}; Q_eff/Q_saddle = {Qeff:.3g}; {o["time"]:.0f}s')
        lines.append('   r      c(r)=min s/Delta   artanh(r)   cEE(r)=min[s-log(|Q|/Qeff)]/Delta   max log(|Q|/Qeff) beyond r   lines reaching r')
        for k, rk in enumerate(RADII):
            sf = o['sfirst'][:, k]; ce = o['cee'][:, k] + np.log(Qeff); qm = o['qmax'][:, k] - np.log(Qeff)
            fin = np.isfinite(sf)
            if fin.any():
                lines.append(f'  {rk:5.2f}   {sf[fin].min()/Delta:10.4f}       {np.arctanh(min(rk, 0.999)):8.4f}      {ce[fin].min()/Delta:10.4f}'
                             f'                         {qm[fin].max():+9.2f}                  {fin.sum()}')
            else:
                lines.append(f'  {rk:5.2f}   not reached (max s along lines {o["sfin"].max()/Delta:.3f} Delta)')
        lines.append(f'   flow accuracy: accumulated |Im dPhi| (should stay 0) max {np.nanmax(o["imerr"]):.1e}; rejected steps (median/max per line) '
                     f'{np.median(o["nrej"]):.0f}/{o["nrej"].max()}')
        if 'sck' in o:
            ds_ = np.diff(o['sck'], axis=1); dq_ = np.diff(o['lqck'], axis=1)
            rat = (ds_ - dq_)/ds_
            fin = np.isfinite(rat)
            if fin.any():
                lines.append(f'   EE decay along the lines beyond s = 10 (checkpoints 10, 20, 40, ...): min over lines/segments of '
                             f'(ds - dlog|Q|)/ds = {rat[fin].min():+.3f} (5% quantile {np.quantile(rat[fin], 0.05):+.3f}); '
                             f'> 0 means the EE weight still decreases monotonically toward the line\'s terminal divisor')
        if 'last' in o:
            lab = ['B_x=0 (m45)', 'B_y=0 (m35)', 'X1->0', 'X2->0', 'Y1->0', 'Y2->0', 's->+-1', 'cos phi->0', 'sin psi->0']
            lst = o['last']; okl = np.all(np.isfinite(lst), axis=0)
            amin = np.argmin(np.log(lst[:, okl] + 1e-300) - np.log(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1.0]))[:, None], axis=0)
            dmin = lst[:, okl][amin, np.arange(okl.sum())]
            cnt = {lab[k]: int((amin == k).sum()) for k in range(9) if (amin == k).any()}
            lines.append(f'   terminal divisor of the lines (smallest normalised distance at the last point): {cnt}; '
                         f'median of that distance {np.median(dmin):.1e}')
        lines.append(f'   min over lines of |w_x| = |1 - xp^2/(ux vx)|: {o["minwx"].min():.2e};  min |m15| {o["min15"].min():.3e}, '
                     f'min |m25| {o["min25"].min():.3e}, min |FP1| {o["minFP1"].min():.3e}')
    return lines


if __name__ == '__main__':
    Dl = float(sys.argv[1]); nu = float(sys.argv[2])
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 1024
    sf = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
    res = run(Dl, nu, N, sf)
    L = report(Dl, nu, res)
    print('\n'.join(L), flush=True)
    tag = f'D{Dl:g}_nu{nu:g}'
    open(os.path.join(HERE, 'cache', f'rd_flow_{tag}.txt'), 'w').write('\n'.join(L) + '\n')
    pickle.dump({k: {kk: vv for kk, vv in o.items()} for k, o in res.items()}, open(os.path.join(HERE, 'cache', f'rd_flow_{tag}.pkl'), 'wb'))

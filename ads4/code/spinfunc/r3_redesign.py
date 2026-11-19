"""
r3_redesign (LARGEJ_OFFDIAG.md §4): can a route-A' functional with RE-OPTIMISED flat weights (psi2, psi3) restore per-state positivity
(n = 1) on the computed near-gap states, given the computed graviton corrections (diagonal + real helicity off-diagonal)?

Family: exactly NF02's LP family (robust/nf_lp.py): psi_k = (1-p)^r p sum_i v_i L_i(2p-1), r = 5, d = 8; k(0) = 1; the flat constraints of
NF02 (R1 flat heavy states with margin delta on the full (u, h, J) grid by cutting planes; R2 scaling limit; R3 b-hat -> oo; scalar matter;
Lemma-M proxy BA = 120), i.e. NF02 itself is the solution at delta = 0.2 with no AdS rows (reproduced: A = 3022.34).
New rows (AdS, at a given Delta_g): for every near-gap state of the q3 model (exact-skeleton states e1/e7 + optional hybrid states), the
per-state action matrix is LINEAR in v:
    F_hh(v) = flat B2 + B3 + curvature + scalar remainder + graviton diagonal correction,   F_+-(v) = graviton off-diagonal,
built node by node exactly as q3_analysis.neargap (same grids, same model M), with psi -> basis functions.  2x2 PSD is imposed by
angle cuts  c^2 F_++ + s^2 F_+- + 2cs F_+- >= dw (c^2 N_++ + s^2 N_+-)  (N from NF02, fixed), added at the min-eigenvector angle of each
violating state (exact PSD in the limit).  1x1 states: F >= dw N.
Objective: NF02's (16 x0 + 0.016 Lfin, i.e. ~A), or (OBJ=bound) A log(Lam/1.5) + B, i.e. the squared bound at that Delta_g.
Usage: R3_DG=300 R3_DW=0.0 R3_OBJ=A python3 r3_redesign.py   -> cache/r3_redesign_<Dg>_<obj>_dw<dw>.txt
env: R3_JEVEN (default 4,...,20,24,28,32), R3_EXTRA_M (hybrid states, default '1.4,1.5,1.75,2,2.5'), R3_DELTA (flat margin, 0.2)
"""
import os, sys, io, json, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROB = os.path.join(HERE, '..', 'robust')
import numpy as np
from scipy.optimize import linprog
from scipy.special import j0

os.environ.setdefault('Q3_JEVEN', os.environ.get('R3_JEVEN', '0,2,4,6,8,10,12,14,16,18,20,24,28,32'))
os.environ.setdefault('Q3_EXTRA_M', os.environ.get('R3_EXTRA_M', '1.4,1.5,1.75,2,2.5'))
import r2_ladder                                  # noqa: merges the r1 caches into the q3 loader
import q3_analysis as Q
from spf_common import NearGap, load_states
from common import gl01, SIGMA
from nf_common import flat_parts_nf, envk
import remk
A9 = Q.A

sys.path.append(ROB)
import nf_lp as NL                                # robust/nf_lp.py (read-only)

R, D = int(os.environ.get("R3_R", "5")), int(os.environ.get("R3_D", "8"))
CACHE = os.path.join(HERE, 'cache')


def basis_fn(p):
    return NL.basis_values(np.asarray(p, float), D, R)                 # (D+1, len(p))


def build_models():
    emods, omods = {}, {}
    for J in Q.JEVEN:
        try:
            g = Q.even_grid(J)
        except Exception:
            continue
        if not g:
            continue
        if Q.OFFSCALE and J >= 4:
            g = {Dl: {**r, **{(kd, 'off'): r[(kd, 'off')]*r['m'] for kd in ('C2', 'C3H')}} for Dl, r in g.items()}
        if len(g) == 1:
            Dl = list(g)[0]; g = {Dl: g[Dl], Dl + 1e-6: dict(g[Dl], m=g[Dl]['m']*(1 + 1e-9))}
        emods[J] = Q.NG.model(g, J, ('C2', 'C3H'))
    for J in Q.JODD:
        try:
            g = Q.odd_grid(J)
        except Exception:
            continue
        if not g:
            continue
        if len(g) == 1:
            Dl = list(g)[0]; g = {Dl: g[Dl], Dl + 1e-6: dict(g[Dl], m=g[Dl]['m']*(1 + 1e-9))}
        omods[J] = Q.odd_model(g)
    return emods, omods


_FLAT = {}


def flat_basis(h, u, J):
    """per-basis stripped flat B2, B3 at (h, u, J) (factored basis values; same integrals as flat_parts_nf)."""
    key = (h, round(u, 12))
    if key not in _FLAT or _FLAT[key][0][-1] < J:
        Js, B2, B3 = NL.parts_factored(h, u, max(J, 40), D, R)
        _FLAT[key] = (Js, B2, B3)
    Js, B2, B3 = _FLAT[key]
    i = int(np.searchsorted(Js, J))
    return B2[i], B3[i]


def state_rows(Dg, emods, omods):
    """list of blocks: dict(J, mfl, DF, hs, F (len(hs) x len(hs) x 2(D+1)) linear maps v -> F_ab)."""
    NGp = NearGap(Dg); pn, pw, Lm = NGp.pn, NGp.pw, NGp.Lam
    V = basis_fn(pn)                                                       # (n, nodes)
    states = list(load_states(Dg))
    for mfl in [float(v) for v in os.environ.get('Q3_EXTRA_M', '').split(',') if v]:
        for J in sorted(set(emods) | set(omods)):
            if J < 11 or (mfl <= 1.3 + 1e-9 and J <= 20):
                continue
            mm = mfl*Lm
            alpha = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*mm*mm))/2
            if alpha < Dg - 2 - 1e-9:
                continue
            r0 = dict(mfl=mfl, J=J, m=mfl*Lm, eta=1.0 + ((2*J + 1)/mfl)**2/(2*Lm*Lm), rho={k: np.zeros_like(pn) for k in range(2, 8)})
            for h in ('++', '+-'):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                states.append((h, 'hyb', r0, None))
    n = D + 1
    blocks = {}
    for h, DF, r0, ref in states:
        mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]
        x = pn*np.sqrt(u); y = x*x
        P = remk.Pders_vec(Lm*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
        B2b, B3b = flat_basis(h, u, J)
        fac = (1 - y)**2 if h == '++' else (1 - y)**-2
        k2 = pw*(2*(P - jz) + envk(2, x)*fac*r0['rho'][2])
        k3 = pw*(u*s*(P - jz) + u*s*envk(3, x)*fac*r0['rho'][3])
        base = np.concatenate([B2b + V @ k2, B3b + V @ k3])
        blocks.setdefault((round(mfl, 9), J, DF), {})[h] = dict(base=base, r0=r0, x=x, y=y, u=u, m=r0['m'])
    out = []
    for (mfl, J, DF), b in blocks.items():
        hs = [h for h in ('++', '+-') if h in b]
        F = np.zeros((len(hs), len(hs), 2*n))
        for kk, h in enumerate(hs):
            F[kk, kk] = b[h]['base']
        e = b[hs[0]]; x, y, u, m = e['x'], e['y'], e['u'], e['m']
        rho2, rho3 = e['r0']['rho'][2], e['r0']['rho'][3]
        r20 = np.polyfit(x[:3], rho2[:3], 2)[-1] if np.any(rho2) else 0.0
        cov = False
        if J % 2 == 0 and J in emods:
            cov = True
            G = emods[J][0](x, m)
            if Q.OFFSCALE and len(hs) == 2:
                G = dict(G); G[('C2', 'off')] = G[('C2', 'off')]/m; G[('C3H', 'off')] = G[('C3H', 'off')]/m
            for kk, h in enumerate(hs):
                Dh = A9.Dsh(h, J, x)
                F[kk, kk, :n] += V @ (pw*A9.env(h, x, 2)*(G[('C2', h)]/m**2 - (rho2 - Dh*r20)))
                F[kk, kk, n:] += V @ (pw*u*A9.env(h, x, 3)*(G[('C3H', h)]/m**2 - A9.SIG[h]*(rho3 - Dh*r20)))
            if len(hs) == 2:
                eo2 = np.sqrt(A9.env('++', x, 2)*A9.env('+-', x, 2))
                F[0, 1, :n] = F[1, 0, :n] = V @ (pw*eo2*G[('C2', 'off')])/m
                F[0, 1, n:] = F[1, 0, n:] = V @ (pw*u*eo2/(2 - y)*G[('C3H', 'off')])/m
        elif J % 2 == 1 and J in omods:
            cov = True
            G = omods[J][0](x, m)
            Dh = Q.P4.Dc('+-', J, x)
            F[0, 0, :n] += V @ (pw*Q.P4.env('+-', x, 2)*(G['C2']/m**2 - (rho2 - Dh*r20)))
            F[0, 0, n:] += V @ (pw*u*Q.P4.env('+-', x, 3)*(G['C3H']/m**2 - Q.P4.SIG['+-']*(rho3 - Dh*r20)))
        # N_s pieces (|B2 + cv2| + |B3 + cv3|) as linear maps (sign fixed later at the reference functional)
        Nmap = []
        for h in hs:
            B2b, B3b = flat_basis(h, b[h]['u'], J)
            xx = b[h]['x']; P = remk.Pders_vec(Lm*pn, b[h]['r0']['eta'], nmax=0)[0]; jz = j0((2*J + 1)*xx)
            Nmap.append((np.concatenate([B2b + V @ (pw*2*(P - jz)), np.zeros(n)]),
                         np.concatenate([np.zeros(n), B3b + V @ (pw*b[h]['u']*SIGMA[h]*(P - jz))])))
        out.append(dict(mfl=mfl, J=J, DF=DF, hs=hs, F=F, Nmap=Nmap, cov=cov))
    return out, Lm


def evalF(blk, v):
    return blk['F'] @ v


def evalN(blk, v):
    return np.array([abs(a @ v) + abs(b @ v) for a, b in blk['Nmap']])


def lam(blk, v, t=1.0):
    Fm = evalF(blk, v)
    if t != 1.0:                                           # scale only the graviton part: off-diagonal (and nothing else) for tests
        Fm = Fm.copy(); Fm[np.triu_indices(len(Fm), 1)] *= t; Fm[np.tril_indices(len(Fm), -1)] *= t
    nn = 1/np.sqrt(evalN(blk, v))
    w, U = np.linalg.eigh(nn[:, None]*Fm*nn[None, :])
    return w[0], U[:, 0]*nn


def summary(blocks, v, pr, tag):
    ls = [(lam(b, v)[0], b) for b in blocks]
    w = min(ls, key=lambda t: t[0])
    pr(f"  [{tag}] min lambda/N_s = {w[0]:+.4f} at J={w[1]['J']} m/Lam={w[1]['mfl']} Df={w[1]['DF']} {w[1]['hs']}; "
       f"#negative = {sum(1 for l, b in ls if l <= 0)} of {len(ls)}")
    for J in sorted(set(b['J'] for b in blocks)):
        lj = [(l, b) for l, b in ls if b['J'] == J]
        lmin = min(lj, key=lambda t: t[0])
        pr(f"      J={J:2d}: min {lmin[0]:+.4f} (m/Lam={lmin[1]['mfl']}, {lmin[1]['DF']})")
    return w[0]


def solve(G, blocks, delta, dw, Lm, obj='A', maxit=80, mproxy=dict(BA=120.0), v0=None, pr=print):
    n, nv = G.n, G.nv
    Aeq = np.zeros((1, nv)); Aeq[0, :n] = -2*G.mu2_6; Aeq[0, n:] = 2*G.mu3_4
    c = np.zeros(nv)
    x0v = np.array([cc[1] for cc in G.b])
    if obj == 'A':
        c[:n] = 16*G.Lfin*1e-3 + 16*x0v
    else:                                                  # A log(Lam/1.5) + B  (= A log Lam + B_AdS)
        c[:n] = 16*G.Lfin + 16*x0v*np.log(Lm/1.5)
    fixed, bfix = [], []
    fixed.append(-np.hstack([G.Bsc, np.zeros((len(G.bsc), n))])); bfix.append(-delta*NL.rho(G.bsc))
    for uu in np.linspace(0, 1, 41):
        for sg in (1, -1):
            r = np.zeros(nv); r[:n] = (1 - delta)*G.pp0; r[n:] = sg*(1 + delta)*uu*G.pp0/2
            fixed.append(r[None, :]); bfix.append(np.zeros(1))
    M = G.matter_rows(); fixed.append(M); bfix.append(np.zeros(len(M)))
    Afix = np.concatenate(fixed); bfix = np.concatenate(bfix)
    sc = np.maximum(np.abs(Afix).max(axis=1), 1e-300); Afix = Afix/sc[:, None]; bfix = bfix/sc
    init = np.unique(np.concatenate([np.arange(0, len(G.u), max(1, len(G.u)//3000)), np.where(G.J <= 12)[0]]))
    cuts = []
    for s2 in (1, -1):
        for s3 in (1, -1):
            Rr = -G.rows(init, np.full(len(init), s2), np.full(len(init), s3), delta)
            cuts.append(Rr/np.maximum(np.abs(Rr).max(axis=1), 1e-300)[:, None])
    Mr, naux = NL.m_proxy_nf(G, delta, **mproxy)
    pad = lambda X: np.hstack([X, np.zeros((X.shape[0], naux))])
    Afix = np.concatenate([pad(Afix), Mr/np.maximum(np.abs(Mr).max(axis=1), 1e-300)[:, None]]); bfix = np.concatenate([bfix, np.zeros(len(Mr))])
    cuts = [pad(X) for X in cuts]; Aeq = pad(Aeq); c = np.concatenate([c, np.zeros(naux)])
    bounds = [(-1e6, 1e6)]*nv + [(0, None)]*naux
    # AdS rows: N fixed at the reference functional (NF02) for the margin dw
    Nref = [evalN(b, v0) for b in blocks] if v0 is not None else None
    ads = []
    def ads_cut(bi, ang):
        b = blocks[bi]
        if len(b['hs']) == 1:
            row = -(b['F'][0, 0]); rhs = -dw*Nref[bi][0]
        else:
            cc, ss = np.cos(ang), np.sin(ang)
            row = -(cc*cc*b['F'][0, 0] + ss*ss*b['F'][1, 1] + 2*cc*ss*b['F'][0, 1]); rhs = -dw*(cc*cc*Nref[bi][0] + ss*ss*Nref[bi][1])
        nr = np.abs(row).max()
        return np.concatenate([row/nr, np.zeros(naux)]), rhs/nr
    if dw is not None:
        for bi, b in enumerate(blocks):
            for ang in ((0.0,) if len(b['hs']) == 1 else (0.0, np.pi/2, np.pi/4, -np.pi/4)):
                ads.append(ads_cut(bi, ang))
    for it in range(maxit):
        Aub = np.concatenate([Afix] + cuts + ([np.array([a[0] for a in ads])] if ads else []))
        bub = np.concatenate([bfix, np.zeros(sum(len(X) for X in cuts))] + ([np.array([a[1] for a in ads])] if ads else []))
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1.0], bounds=bounds, method='highs')
        if res.status != 0:
            return res.status, None, {'msg': res.message}
        v = res.x[:nv]
        B2, B3 = G.parts(v); FF = B2 + B3; N = np.abs(B2) + np.abs(B3)
        g = FF - delta*N
        viol = np.where(g < -1e-6*np.maximum(N, 1e-300))[0]
        nadd = 0
        if dw is not None:
            for bi, b in enumerate(blocks):
                Fm = evalF(b, v)
                if len(b['hs']) == 1:
                    continue
                nn = 1/np.sqrt(Nref[bi])
                wv, U = np.linalg.eigh(nn[:, None]*Fm*nn[None, :])
                if wv[0] < dw - 1e-7:
                    vec = U[:, 0]*nn; ang = np.arctan2(vec[1], vec[0])
                    ads.append(ads_cut(bi, ang)); nadd += 1
        if len(viol) == 0 and nadd == 0:
            break
        if len(viol):
            order = viol[np.argsort(g[viol]/N[viol])][:4000]
            Rr = -G.rows(order, np.sign(B2[order]), np.sign(B3[order]), delta)/N[order][:, None]
            cuts.append(pad(Rr))
    else:
        return 9, v, {'msg': 'no convergence'}
    F = NL.to_functional(G, v)
    A, B, BA = F.AB()
    return 0, v, {'A': A, 'B': B, 'B_AdS': BA, 'it': it}


def nf_vector():
    """NF02's LP vector (r = 5, d = 8, delta = 0.2, proxy 120), re-expressed in the (R, D) basis by least squares (exact if the span
    contains it, e.g. R = 5, D >= 8)."""
    nf = json.load(open(os.path.join(ROB, 'cache', 'nf', 'pareto_r5_d8_proxy120.json')))
    v8 = np.array(next(e for e in nf if abs(e['delta'] - 0.2) < 1e-12)['v'])
    if (R, D) == (5, 8):
        return v8
    p = np.linspace(0.001, 0.999, 400)
    B8 = NL.basis_values(p, 8, 5); BB = basis_fn(p)
    out = []
    for k in range(2):
        y = v8[9*k:9*(k + 1)] @ B8
        out.append(np.linalg.lstsq(BB.T, y, rcond=None)[0])
    return np.concatenate(out)


def bound(A, BA, Dg):
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    val = A*np.log(Lm) + BA
    return 240*np.sqrt(val)/Lm**4 if val > 0 else np.nan


def main():
    Dg = int(os.environ.get('R3_DG', '300'))
    dws = [float(v) for v in os.environ.get('R3_DW', '0.0').split(',')]
    obj = os.environ.get('R3_OBJ', 'A')
    delta = float(os.environ.get('R3_DELTA', '0.2'))
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    emods, omods = build_models()
    blocks, Lm = state_rows(Dg, emods, omods)
    pr(f"# r3 redesign at Delta_g = {Dg} (Lam = {Lm:.2f}); {len(blocks)} state blocks; graviton-corrected even J {sorted(emods)}, odd J "
       f"{sorted(omods)}; hybrid m/Lam {os.environ.get('Q3_EXTRA_M')}; flat margin delta = {delta}; objective {obj}  [{time.time() - t0:.0f}s]")
    G = NL.NFGrid(D, R)
    v_nf = nf_vector()
    F0 = NL.to_functional(G, v_nf); A0, B0, BA0 = F0.AB()
    pr(f"  NF02 (LP vector, delta = 0.2, proxy 120): A = {A0:.2f}, B_AdS = {BA0:.2f}; bound at 300 / 1000: {bound(A0, BA0, 300):.3e} / {bound(A0, BA0, 1000):.3e}")
    summary(blocks, v_nf, pr, 'NF02 with the computed corrections')
    st, v1, info = solve(G, blocks, delta, None, Lm, obj, v0=v_nf)
    pr(f"  LP without AdS rows (check: must reproduce NF02 for obj=A): status {st}, A = {info.get('A', np.nan):.2f}, B_AdS = {info.get('B_AdS', np.nan):.2f}")
    res = {}
    for dw in dws:
        st, v, info = solve(G, blocks, delta, dw, Lm, obj, v0=v_nf)
        if st != 0:
            pr(f"  AdS rows with margin dw = {dw}: INFEASIBLE / status {st} ({info.get('msg')})"); res[dw] = None; continue
        pr(f"  AdS rows with margin dw = {dw}: A = {info['A']:.2f}, B_AdS = {info['B_AdS']:.2f}, bound at 300 / 1000: "
           f"{bound(info['A'], info['B_AdS'], 300):.3e} / {bound(info['A'], info['B_AdS'], 1000):.3e} (iterations {info['it']})")
        summary(blocks, v, pr, f'redesign dw={dw}')
        res[dw] = dict(v=v, **info)
    pickle.dump(dict(res=res, v_nf=v_nf), open(os.path.join(CACHE, f'r3_redesign_{Dg}_{obj}{os.environ.get("R3_TAG", "")}.pkl'), 'wb'))
    open(os.path.join(CACHE, f'r3_redesign_{Dg}_{obj}{os.environ.get("R3_TAG", "")}.txt'), 'w').write(buf.getvalue())
    pr(f"[{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()

"""
r5_universal (LARGEJ_OFFDIAG.md §3-4): the near-gap / large-b-hat 2x2 test with the UNIVERSAL (leading-Regge) off-diagonal of r4:
    F_+- density = eo2 beta^2 g2(nu beta) psi2  +  u eo2/(2-y) beta^2 g3(nu beta) psi3        (g2, g3: r4 medians, t in [1, 12]; [M] outside)
    F_hh = stripped flat B2 + B3 (exact, at the state's u, J) + leading-Regge curvature 2 int psi2 (P - J0) + u s int psi3 (P - J0)
           (= ROUTEA_LARGE_B's F_hyb; scalar remainder and graviton DIAGONAL corrections dropped [M]: <= 0.2 N_s at the binding states)
for an arbitrary weight pair in the NF02 LP family (robust/nf_lp basis r = 5, d = 8).  Modes:
    check  : NF02 on a (m/Lam, b-hat) grid at Delta_g -> min eigenvalue map, crossing b-hat per m/Lam
    lp     : NF02's LP (flat constraints, delta = 0.2, proxy 120) + these AdS rows on the grid with beta <= R5_BMAX: feasibility, A, B_AdS
Usage: R5_DG=300 python3 r5_universal.py check|lp   -> cache/r5_<mode>_<Dg>.txt
"""
import os, sys, io, json, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
from scipy.special import j0
from scipy.interpolate import CubicSpline
from common import gl01, SIGMA
import remk
sys.path.append(os.path.join(HERE, '..', 'robust'))
import nf_lp as NL
import r3_redesign as R3

CACHE = os.path.join(HERE, 'cache')
R, D = int(os.environ.get("R3_R", "5")), int(os.environ.get("R3_D", "8"))
GT = pickle.load(open(os.path.join(CACHE, 'r4_g.pkl'), 'rb'))
SCALE_G = float(os.environ.get('R5_GSCALE', '1.0'))      # common factor on g (tolerance scans)


def _gmed(key):
    rows = [r for r in GT['rows'] if r['J2m'] <= 1.0 and r['D'] >= 400]
    tg = GT['tg']; vals = []
    for r in rows:
        ok = (tg >= r['t'][0]) & (tg <= r['t'][-1]); v = np.full(len(tg), np.nan)
        v[ok] = CubicSpline(r['t'], r[key])(tg[ok]); vals.append(v)
    g = np.nanmedian(np.array(vals), axis=0)
    g = np.where(np.isfinite(g), g, 0.0)
    cs = CubicSpline(tg, g)
    def f(t):
        t = np.asarray(t, float)
        return SCALE_G*np.where(t < tg[0], g[0], np.where(t > tg[-1], 0.0, cs(np.clip(t, tg[0], tg[-1]))))
    return f


G2, G3 = _gmed('f2'), _gmed('f3')


def state(Dg, mfl, J, h, V, pn, pw, Lm):
    """linear maps (2n) for F_hh (flat + curvature) and F_+- (universal) at (m/Lam, J, h)."""
    u = 1/mfl**2; s = SIGMA[h]; x = pn*np.sqrt(u); y = x*x
    m = mfl*Lm
    eta = 1.0 + ((2*J + 1)/mfl)**2/(2*Lm*Lm)
    P = remk.Pders_vec(Lm*pn, eta, nmax=0)[0]; jz = j0((2*J + 1)*x)
    B2b, B3b = R3.flat_basis(h, u, J)
    n = D + 1
    base = np.concatenate([B2b + V @ (pw*2*(P - jz)), B3b + V @ (pw*u*s*(P - jz))])
    Nm = (np.concatenate([B2b + V @ (pw*2*(P - jz)), np.zeros(n)]), np.concatenate([np.zeros(n), B3b + V @ (pw*u*s*(P - jz))]))
    be = np.arccosh(eta); t = Lm*pn*be
    eo2 = (2 - y)/(1 - y)**2
    off = np.concatenate([V @ (pw*eo2*be*be*G2(t)), V @ (pw*u*eo2/(2 - y)*be*be*G3(t))])
    return base, off, Nm, be


def grid_states(Dg, Lm, bmax_beta, mfls, bhs):
    out = []
    for mfl in mfls:
        for bh in bhs:
            J = int(round((bh*mfl - 1)/2))
            if J < 4:
                continue
            m = mfl*Lm
            alpha = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*m*m))/2
            if alpha < Dg - 2 - 1e-9:                     # twist gap
                continue
            be = np.arccosh(1.0 + ((2*J + 1)/mfl)**2/(2*Lm*Lm))
            if be > bmax_beta:
                continue
            out.append((mfl, J))
    return sorted(set(out))


def blocks_for(Dg, sts):
    NGp = R3.NearGap(Dg); pn, pw, Lm = NGp.pn, NGp.pw, NGp.Lam
    V = R3.basis_fn(pn)
    blocks = []
    for mfl, J in sts:
        hs = ['++', '+-'] if J % 2 == 0 else ['+-']
        n = D + 1
        F = np.zeros((len(hs), len(hs), 2*n)); Nmap = []
        for k, h in enumerate(hs):
            base, off, Nm, be = state(Dg, mfl, J, h, V, pn, pw, Lm)
            F[k, k] = base; Nmap.append(Nm)
        if len(hs) == 2:
            F[0, 1] = F[1, 0] = off
        blocks.append(dict(mfl=mfl, J=J, DF='uni', hs=hs, F=F, Nmap=Nmap, cov=True, bh=(2*J + 1)/mfl, beta=be))
    return blocks, Lm


def main():
    mode = sys.argv[1]
    Dg = int(os.environ.get('R5_DG', '300'))
    bmax = float(os.environ.get('R5_BMAX', '0.12'))
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    G = NL.NFGrid(D, R)
    v_nf = R3.nf_vector()
    mfls = [float(v) for v in os.environ.get('R5_MFLS', '1.02,1.05,1.1,1.2,1.3,1.5,2.0,3.0,5.0').split(',')]
    bhs = np.unique(np.round(np.geomspace(10, 2*Lm, 60), 3))
    sts = grid_states(Dg, Lm, bmax, mfls, bhs)
    blocks, _ = blocks_for(Dg, sts)
    pr(f"# r5 [{mode}] Delta_g = {Dg} (Lam = {Lm:.1f}); universal off-diagonal (r4 medians, scale {SCALE_G}); {len(blocks)} states, "
       f"beta <= {bmax}; m/Lam in {mfls}  [{time.time() - t0:.0f}s]")
    if mode == 'eval':                                     # another weight vector (pickle path, key) on the same grid
        import pickle as pk
        fn, key = os.environ['R5_VFILE'], os.environ.get('R5_VKEY', '0.1')
        d = pk.load(open(fn, 'rb'))
        v_nf = d['res'][float(key)]['v'] if 'res' in d else d['v']
        pr(f"  evaluating weights from {os.path.basename(fn)} [{key}]")
    if mode == 'states':                                   # explicit states "mfl:J,..." -> ratio table
        sts2 = [(float(a.split(':')[0]), int(a.split(':')[1])) for a in os.environ['R5_STATES'].split(',')]
        bl2, _ = blocks_for(Dg, sts2)
        for b in bl2:
            Fm = R3.evalF(b, v_nf); Nn = R3.evalN(b, v_nf)
            pr(f"  m/Lam={b['mfl']:.4f} J={b['J']:3d} b-hat={b['bh']:.2f} beta={b['beta']:.4f}: F_++/N {Fm[0,0]/Nn[0]:+.3f}, F_+-/N {Fm[1,1]/Nn[1]:+.3f}, "
               f"N_s = {Nn[0]:.4e} / {Nn[1]:.4e}, |X| = {abs(Fm[0,1]):.4e} (|X| Lam^2 = {abs(Fm[0,1])*Lm*Lm:.0f}), |X|/sqrt(NN) = "
               f"{abs(Fm[0,1])/np.sqrt(Nn.prod()):.3f}, lambda/N_s = {R3.lam(b, v_nf)[0]:+.3f}")
    elif mode in ('check', 'eval'):
        for mfl in mfls:
            bl = sorted([b for b in blocks if b['mfl'] == mfl and len(b['hs']) == 2], key=lambda b: b['bh'])
            if not bl:
                continue
            ls = [(b['bh'], R3.lam(b, v_nf)[0], abs(R3.evalF(b, v_nf)[0, 1])/np.sqrt(R3.evalN(b, v_nf).prod())) for b in bl]
            cross = next((bh for bh, l, o in ls if l < 0), None)
            pr(f"  m/Lam = {mfl}: first negative at b-hat = {cross}; (b-hat, lam/N_s, |X|/sqrt(NN)): " +
               " ".join(f"({bh:.0f},{l:+.2f},{o:.2f})" for bh, l, o in ls[::3]))
    else:
        dws = [float(v) for v in os.environ.get('R5_DW', '0.0').split(',')]
        R3.summary(blocks, v_nf, pr, 'NF02, universal off-diagonal')
        for dw in dws:
            st, v, info = R3.solve(G, blocks, 0.2, dw, Lm, os.environ.get('R5_OBJ', 'A'), v0=v_nf)
            if st != 0:
                pr(f"  LP + AdS rows (dw = {dw}): status {st}: {info.get('msg')}")
                continue
            pr(f"  LP + AdS rows (dw = {dw}): A = {info['A']:.2f}, B_AdS = {info['B_AdS']:.2f}, bound 300 / 1000: "
               f"{R3.bound(info['A'], info['B_AdS'], 300):.3e} / {R3.bound(info['A'], info['B_AdS'], 1000):.3e}")
            R3.summary(blocks, v, pr, f'redesign dw={dw}')
            pickle.dump(dict(v=v, info=info), open(os.path.join(CACHE, f'r5_lp_{Dg}_b{bmax}_dw{dw}.pkl'), 'wb'))
    open(os.path.join(CACHE, f'r5_{mode}_{Dg}_b{bmax}_r{R}d{D}{os.environ.get("R5_TAG", "")}.txt'), 'w').write(buf.getvalue())
    pr(f"[{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()

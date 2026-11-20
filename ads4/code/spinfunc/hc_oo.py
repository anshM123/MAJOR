"""
hc_oo (HYPOTHESES_CLOSE.md, H_OO): the even-J parity-odd (OO) block of Omega^EE beyond leading order.
(1) Gamma_OO rows: from the hc_oo_rows K2 records (cache/hc_p3_{EE,C3H}_J<J>.pkl; 3x3 [even(2), odd(1)] action, OO = entry [-1, -1]),
    in exactly the p4/l9 linear-chain conventions with the ++ shapes (X^O_same couples like ++ at leading order: OO = N_same x ++ density,
    N = tan^2(pi tau/2) >= 0, PARITY_ODD_ROUTEA §4/§12): heavy reference Delta_ref = 2e4 order by order (K2), logK - i pi (Delta + J),
    per-row constant a_O (absorbs N) fitted on x <= 0.25, nu >= 15, spinning local de-shift (1 +- u_EE)/2 for k = 2 (Jacobi-fixed Dc),
    none for k = 3, sigma_++ = +1;  Gamma_OO = m^2 [G/(a_O env_++) - D_++];  then the curvature image is subtracted (sm_curvfix, 'row').
(2) Fold into the chain: the EE_NEARGAP main chain (w2/q3, ee1000 loaders, SM_CURV = row) is run with the even rows of J in HC_OOJ
    REPLACED by the OO rows (++ keys = Gamma_OO, +- and off-diagonal keys = 0), so that for every near-gap state the chain's ++ entry is
    base_++ + dF_OO;  the OO block's per-state value is lambda_OO/N_s = (base_++ + dF_OO)/N_s,++ (1x1; N > 0 drops out of the sign).
    Compared with the EE block's own ++ diagonal (base_++ + dF_++)/N_s,++ and the EE 2x2 minimum from the unmodified run.
Usage: HC_OOJ=4,6,8,10 python3 hc_oo.py [Dg ...]   -> cache/hc_oo_<HC_TAG>.txt/.pkl
"""
import os, sys, io, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
os.environ.setdefault('W2_BASIS', 'ee')
os.environ.setdefault('W2_MODES', 'ee')
os.environ.setdefault('SM_CURV', 'row')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
OOJ = [int(v) for v in os.environ.get('HC_OOJ', '4,6,8,10').split(',') if v]
TAG = os.environ.get('HC_TAG', 'x')
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


import sm_curvfix as CF            # noqa: E402  (imports w2_chain; curvature fix wrappers)
import ee1000_chain as CH          # noqa: E402  (ee1000 loaders + Jacobi fix, as sm_chainfix w2)
Q = CH.Q
W2 = CH.W
P4 = Q.P4


def load_oo(kind, J):
    fn = os.path.join(C, 'hc_p3_%s_J%d.pkl' % (kind, J))
    R = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    R = {(k[0], k[1]): v for k, v in R.items() if k[2] == 1.0}
    return Q._strip(R, J)


def getrec(R, D, nu):
    k = (float(D), round(float(nu), 9))
    if k in R:
        return R[k]
    c = [kk for kk in R if kk[0] == float(D) and abs(kk[1] - nu) < 1e-7*nu]
    return R[c[0]] if c else None


def grid_oo(J, K=2, curv=True):
    """per Delta row: m, xs, a_O, Gamma_OO for kinds C2 (= EE crystal) and C3H, h = '++' conventions; LO diagnostics."""
    h = '++'
    RC = {'C2': load_oo('EE', J), 'C3H': load_oo('C3H', J)}
    Ds = sorted(set(k[0] for k in RC['C2'] if k[0] < DREF))
    out = {}
    for D in Ds:
        m = P4.mphys(D, J); mr = P4.mphys(DREF, J)
        xs = sorted(set(round(k[1]/m, 6) for k in RC['C2'] if k[0] == D))
        xs = [x for x in xs if all(getrec(RC[kd], D, x*m) is not None and getrec(RC[kd], DREF, x*m) is not None for kd in RC)]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}; c0 = None; lo = {}
        for kd in ('C2', 'C3H'):
            k = 2 if kd == 'C2' else 3
            G = []; G0 = []
            for x in xs:
                nu = x*m; st = getrec(RC[kd], D, nu); rf = getrec(RC[kd], DREF, nu); xr = nu/mr
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                ords = ('M0', 'M2', 'M4')
                v = np.exp(dl)*P4.pred([st[o][-1, -1] for o in ords], [rf[o][-1, -1] for o in ords], K)*fk
                v0 = np.exp(dl)*st['M0'][-1, -1]/rf['M0'][-1, -1]*fk
                G.append(v*P4.env(h, xr, k)*P4.Dc(h, J, xr)); G0.append(v0)
            rows[kd] = (k, np.array(G))
            lo[kd] = np.array(G0)
        sel = (xs <= 0.25 + 1e-9) & (xs*m >= 15.0)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]

        def fit(g):
            A_ = np.array([P4.Dc(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(A_, (g/P4.env(h, xs, 2))[sel], rcond=None)[0][0]
        a = fit(rows['C2'][1].real); a_raw = a
        sh = {k: P4.shift(J, h, k, xs, m) for k in (2, 3)}
        for _ in range(3):
            a = fit(rows['C2'][1].real - a*sh[2])
        rows = {kd: (k, G - a*sh[k]) for kd, (k, G) in rows.items()}
        rec = dict(m=m, xs=xs, a=a, a_raw=a_raw)
        for kd, (k, G) in rows.items():
            g = m*m*(G.real/(a*P4.env(h, xs, k)) - P4.Dc(h, J, xs))
            rec[(kd, 'raw')] = g.copy()
            if curv:
                g = g - CF.curv_image(kd, h, J, m, xs, Q.A.env, Q.A.SIG)
            rec[kd] = g
            rec[(kd, 'im')] = m*G.imag/(a*P4.env(h, xs, k))
        # LO diagnostic: order-0 C3H/C2 ratio of the OO entry, relative to the ++ expectation (x-independent constant)
        rec['lo_c3c2'] = lo['C3H']/lo['C2']
        out[D] = rec
    return out


def oo_as_even_rows(J, g):
    """OO grid -> the even-row dict format of A.grid (++ keys = Gamma_OO; +- and off keys = 0)."""
    out = {}
    for D, r in g.items():
        z = 0*r['xs']
        out[D] = {'m': r['m'], 'xs': r['xs'], 'a': {'++': r['a'], '+-': r['a']},
                  ('C2', '++'): r['C2'], ('C3H', '++'): r['C3H'], ('C2', '+-'): z, ('C3H', '+-'): z,
                  ('C2', '++', 'im'): r[('C2', 'im')], ('C3H', '++', 'im'): r[('C3H', 'im')], ('C2', '+-', 'im'): z, ('C3H', '+-', 'im'): z,
                  ('C2', 'off'): z, ('C3H', 'off'): z, ('C2', 'off', 'im'): z + 1e-30, ('C3H', 'off', 'im'): z}
    return out


if __name__ == '__main__':
    t0 = time.time()
    Dgs = [int(a) for a in sys.argv[1:]] or [300, 1000]
    W2.set_mode(True)
    G = {}
    pr(f"# hc_oo: OO block (even J parity-odd X^O_same) rows, J = {OOJ}; SM_CURV = {CF.MODE}")
    for J in OOJ:
        g = grid_oo(J)
        G[J] = g
        for D, r in sorted(g.items()):
            pr(f"  J={J} Delta={D:g} m={r['m']:.1f} a_O={r['a']:.6g} (no de-shift {r['a_raw']:.6g}); x = " + " ".join(f"{x:.3f}" for x in r['xs']))
            for kd in ('C2', 'C3H'):
                pr(f"    {kd:3s} Gamma_OO (curv-subtracted): " + " ".join(f"{v:+7.1f}" for v in r[kd]))
                pr(f"    {kd:3s} Gamma_OO (raw)            : " + " ".join(f"{v:+7.1f}" for v in r[(kd, 'raw')]))
            c = r['lo_c3c2']
            pr(f"    LO order-0 C3H/C2 of the OO entry x (2-y) [++ coupling: x-independent]: " +
               " ".join(f"{(cc*(2 - x*x)).real:+.4f}" for x, cc in zip(r['xs'], c)))
    # (2) chain: unmodified run (EE block) and the run with OO rows in place of the even rows of OOJ
    Q.JEVEN[:] = [J for J in OOJ]
    Q.JODD[:] = []
    res0 = Q.run(Dgs, lambda s: None)[1]
    _grid = Q.A.grid
    Q.A.grid = lambda J, K, kinds=('C2', 'C3H'), *a, **k: (oo_as_even_rows(J, G[J]) if J in G and G[J] else _grid(J, K, kinds, *a, **k))
    pr("# rows used for the OO fold: " + "; ".join(f"J={J}: " + ",".join(f"{D:g}" for D in sorted(G[J])) for J in OOJ))
    res1 = Q.run(Dgs, lambda s: None)[1]
    Q.A.grid = _grid
    store = {}
    for Dg in Dgs:
        pr(f"## Delta_g = {Dg}")
        r0 = {(round(r['mfl'], 9), r['J'], r['DF']): r for r in res0[Dg]['recs']}
        r1 = {(round(r['mfl'], 9), r['J'], r['DF']): r for r in res1[Dg]['recs']}
        for J in OOJ:
            ks = sorted(k for k in r1 if k[1] == J and len(r1[k]['hs']) == 2)
            rows = []
            for k in ks:
                a, b = r0[k], r1[k]
                lam_ee = Q.lam(a, 1.0)
                pp_ee = (a['base'][0] + a['dF'][0, 0])/a['Ns'][0]
                oo = (b['base'][0] + b['dF'][0, 0])/b['Ns'][0]
                oo0 = b['base'][0]/b['Ns'][0]
                rows.append((k[0], k[2], lam_ee, pp_ee, oo, oo0, b['dF'][0, 0]/b['Ns'][0], a['dF'][0, 0]/a['Ns'][0]))
            store[(Dg, J)] = rows
            w = min(rows, key=lambda t: t[4])
            pr(f"  J={J:2d}: min lambda_OO/N_s = {w[4]:+.4f} at m/Lam={w[0]:.4f} ({w[1]}) [LO {w[5]:+.4f}]; max |dF_OO|/N_++ "
               f"{max(abs(t[6]) for t in rows):.4f} (EE ++: {max(abs(t[7]) for t in rows):.4f}); EE-block min lam {min(t[2] for t in rows):+.4f}")
            sel_ = [t for t in rows if t[0] <= 1.3 + 1e-9 and abs(round(t[0]*20) - t[0]*20) < 1e-6] or rows[:6]
            pr("      near gap (m/Lam: OO [dF_OO] | EE ++ diag): " + "  ".join(f"{t[0]:.2f}{'h' if t[1] == 'hyb' else ''}: {t[4]:+.4f} [{t[6]:+.4f}] | {t[3]:+.4f}" for t in sel_))
    pr(f"# despiked (even rows, J, Delta, x): {W2.DESPIKE_LOG}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'hc_oo_%s.txt' % TAG), 'w').write(buf.getvalue())
    pickle.dump(dict(grids=G, states=store), open(os.path.join(C, 'hc_oo_%s.pkl' % TAG), 'wb'))

"""
n3_analysis (doc 25): the l9 linear-form chain with DIRECT values on x < 0.2 from the small-nu contour method (n3_run pickles merged
into n1_analysis's direct database), and
  (i)   Gamma per node at the Delta = 100 row: direct vs K2 (n1_analysis.compare);
  (ii)  t(x) per node: t = (Gamma_dir - c2s)/(Gamma_K2 - c2s), c2s = m^2 (rho_2(x) - D_h rho_2(0)) of a representative exact-skeleton
        state (m/Lam = 1.1, same J) -- the graviton-specific correction relative to the base model, as scaled by n2_neargap_split;
  (iii) n at Delta_g = 100 with the direct values wherever they exist (x >= 0.1 at the Delta = 100 row), and the tolerance window of
        the remaining uncontrolled nodes x < N3_XC (default 0.075, i.e. the x = 0.05 node) via n2_neargap_split.
Usage: python3 n3_analysis.py  -> cache/n3_analysis.txt
"""
import os, sys, io, pickle, contextlib
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import n1_analysis as N
import l9_analysis as A

_orig = N.dload
USE_N3 = os.environ.get('N3_USE', '1') == '1'
USE_N4 = os.environ.get('N4_USE', '0') == '1'      # NEARGAP100_SMALLNU.md: also the n4 thimble points at x = 0.05 (default off: old behaviour)
CV = os.environ.get('N3_CV', '0') == '1'


def dload(kind, J):
    d = dict(_orig(kind, J))
    f = os.path.join(HERE, 'cache', 'n3_pts_%s_J%d.pkl' % (kind, J))
    if USE_N3 and os.path.exists(f):
        for k, v in pickle.load(open(f, 'rb')).items():
            d[k] = v
    f4 = os.path.join(HERE, 'cache', 'n4_pts_%s_J%d.pkl' % (kind, J))
    if USE_N4 and os.path.exists(f4):
        for k, v in pickle.load(open(f4, 'rb')).items():
            d[k] = v
    return d


N.dload = dload

SPIN = os.environ.get('N3_SPIN', '0') == '1'
if SPIN:                     # C^_2' (H7_I14_LIGHTSIDE.md sect. 18-22): spinning de-shift for k = 2 (read-only import of h7i14/v3_heavy_impact)
    sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
    import v3_heavy_impact as HV
    A.shift = HV.shift_spin


def tnodes(pr, Js=(4, 6, 8, 10), kinds=('C2', 'C3H')):
    import n3_exact as NX
    e1 = NX.db(100)
    for J in Js:
        N.DIRECT = True
        gd = N.grid_direct(J, 2, kinds, Ds=[100.0])
        N.DIRECT = False
        g2 = N.grid_direct(J, 2, kinds, Ds=[100.0])
        N.DIRECT = True
        if 100.0 not in gd:
            continue
        r, r2 = gd[100.0], g2[100.0]
        st = e1[(1.1, J)]
        xs = r['xs']; m = r['m']
        hels = ('++', '+-') if J >= 4 else ('++',)
        rho_x = np.interp(xs, st['x'], st['rho'][2]); r20 = np.polyfit(st['x'][:3], st['rho'][2][:3], 2)[-1]
        pr(f"## J = {J}, Delta = 100 row (m = {m:.2f}): t(x) = (Gamma_dir - c2s)/(Gamma_K2 - c2s), c2s from the exact state m/Lam = 1.1; "
           f"K2 row re-normalised to the direct a_h (a_K2/a_dir = " + ", ".join(f"{hh}: {r2['a'][hh]/r['a'][hh]:.4f}" for hh in hels) + ")")
        pr("   x       : " + " ".join(f"{x:7.2f}" for x in xs))
        for kd in kinds:
            for h in hels:
                sg = A.SIG[h] if kd != 'C2' else 1.0
                c2s = st['m']**2*(rho_x - A.Dsh(h, J, xs)*r20)*sg
                gdir, gk2, ge = r[(kd, h)], r2[(kd, h)], r[(kd, h, 'err')]
                # common normalisation: the K2 row re-normalised with the direct chain's a_h (a_h is a per-state congruence; it is fit on
                # the smallest-x nodes, so it moves when they become direct)
                DD = sg*A.Dsh(h, J, xs)*m*m
                gk2 = (r2['a'][h]/r['a'][h])*(gk2 + DD) - DD
                t = (gdir - c2s)/(gk2 - c2s)
                pr(f"   {kd:3s} {h:2s} dir: " + " ".join(f"{v:7.0f}" for v in gdir) + "\n            K2 : " + " ".join(f"{v:7.0f}" for v in gk2)
                   + "\n            err: " + " ".join(f"{v:7.0f}" if np.isfinite(v) else "      -" for v in ge)
                   + "\n            c2s: " + " ".join(f"{v:7.0f}" for v in c2s)
                   + "\n            t  : " + " ".join(f"{v:7.2f}" for v in t))
            if len(hels) == 2:
                pr(f"   {kd:3s} off dir: " + " ".join(f"{v:7.0f}" for v in r[(kd, 'off')]) + "\n            K2 : " + " ".join(f"{v:7.0f}" for v in r2[(kd, 'off')]))


def neargap(pr, xc=0.075):
    import n2_neargap_split as NG
    tmp = os.path.join(HERE, 'cache', 'n1_tmp'); os.makedirs(os.path.join(tmp, 'cache'), exist_ok=True)
    NG.HERE = tmp
    NG.A.grid = lambda J, K, kinds=('C2', 'C3H'), **kw: N.grid_direct(J, K, kinds, **kw)
    def run(t, XC):
        NG.TLOW[0] = t; NG.XC = XC
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            NG.main('C23H_K2', [100])
        out = f.getvalue()
        l = [s for s in out.split('\n') if 'Delta_g = 100' in s][0]
        n = float(l.split('n = ')[1].split(' ')[0]) if 'n = None' not in l else None
        lam = float(l.split('min lam/N_s = ')[1].split(' ')[0])
        return n, lam, out
    N.DIRECT = True
    n, lam, out = run(1.0, xc)
    pr("# near-gap 2x2 [C23H_K2], Delta_g = 100, direct values wherever computed (n1: x >= 0.2; n3: x = 0.1, 0.15 [and 0.05 if run])")
    pr(out)
    for XC, lab in ((xc, f'nodes x < {xc}'), (0.2, 'all nodes x < 0.2')):
        def ok(t):
            return run(t, XC)[0] == 1.0
        edges = []
        for sgn in (1, -1):
            lo, hi = 1.0, 1.0 + sgn
            while ok(hi) and abs(hi) < 50:
                lo, hi = hi, hi + sgn*abs(hi - 1.0 + sgn)
            for _ in range(12):
                mid = 0.5*(lo + hi); lo, hi = (mid, hi) if ok(mid) else (lo, mid)
            edges.append(lo)
        pr(f"   tolerance on {lab} (correction there x t, direct where available): n = 1 for t in [{edges[1]:+.2f}, {edges[0]:+.2f}]")
        for t in (-1.0, 0.0, 2.0, 3.0):
            n_, lam_, _ = run(t, XC)
            pr(f"      t = {t:+.1f}: n = {n_}, min lam/N_s = {lam_:+.4f}")


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    tnodes(pr)
    neargap(pr)
    open(os.path.join(HERE, 'cache', 'n3_analysis%s%s.txt' % ('_spin' if SPIN else '', '_n4' if USE_N4 else '')), 'w').write("\n".join(lines) + "\n")

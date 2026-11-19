"""v4 (H7_I14_LIGHTSIDE.md §23.4): parity-odd sector (odd J = 5, 7, 9, +-) at Delta_g = 100 with C^_2', using
  * direct (p8, RQMC on the real contour) values on x >= 0.2 of the Delta = 100 row,
  * the n3 Picard-Lefschetz contour points (H0 §25.7; spinfunc/cache/n3_odd_pts_{C2,C3H}_J{5,7,9}.pkl, p8 format) at x = 0.1, 0.15,
  * the K2 Laplace chain at x = 0.05 (and on the rows Delta = 200 ... 3200),
and the combined even + odd n.  Read-only imports of spinfunc/p8_analysis, p4_odd_analysis, p6_odd_neargap (nothing written there).

Chain: p8_analysis.grid (PARITY_ODD_ROUTEA §11 main protocol: P8_XROWS=0, i.e. the §6 rows, direct on the Delta = 100 row), with
  (i) the n3 contour points merged into p8's direct database (same format; accepted by p8's own rule: rel. RQMC error <= 5 %, else
      abs. error <= 60 c2);
  (ii) the beta phase stripped as usual: p8's P8_PHASE=1 multiplies the state/reference ratio by e^{+i pi (Delta_s - Delta_ref)}, which is
      exactly logK -> logK - i pi (Delta + J) on state and reference (PARITY_ODD_ROUTEA §12.3; J cancels at equal J).  Asserted here on
      the shadow entries used: |Im logK - pi (Delta + J)| < 1e-8;
  (iii) the k = 2 de-shift with the exact spinning u_T (C^_2', v4_hclosed) instead of the scalar g (k = 3 keeps g), as v4_odd_impact.
Near gap: p6_odd_neargap conventions; the correction on GL nodes x < XC (state-dependent x = p/m) can be scaled by t (tolerance of the
unresolved x = 0.05 node region, XC = 0.075 as n3_analysis for the even sector), and the whole correction by t* as p6.
Usage: python3 v4_odd100.py   -> cache/v4_odd100.txt
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ['P8_XROWS'] = os.environ.get('P8_XROWS', '0')
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
sys.dont_write_bytecode = True
import numpy as np
from scipy.special import j0
import p4_odd_analysis as A
import p8_analysis as P8
import p2_shadow_odd as P2
import p6_odd_neargap as P6
from spf_common import NearGap, load_states
from common import gl01, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk
from v4_odd_impact import make_shift, _orig_shift

USE_N3 = [True]
_dload = P8.dload


def dload_merged(kind, J):
    db = dict(_dload(kind, J))
    if USE_N3[0]:
        fn = os.path.join(SF, 'cache', 'n3_odd_pts_%s_J%d.pkl' % (kind, J))
        if os.path.exists(fn):
            for k, v in pickle.load(open(fn, 'rb')).items():
                db[(float(k[0]), float(k[1]))] = v
    return db


P8.dload = dload_merged


def check_phase(J, Ds):
    out = []
    for D in Ds:
        b, logK, _ = P2.shadow_odd(float(D), J)
        out.append(abs(np.imag(logK) - np.pi*(D + J)))
    return max(out)


def grids(shift, use_n3=True):
    USE_N3[0] = use_n3
    A.shift = shift
    try:
        return {J: P8.grid(J) for J in (5, 7, 9)}
    finally:
        A.shift = _orig_shift
        USE_N3[0] = True


def state_parts(G, Dg=100, XC=0.075):
    """p6_odd_neargap.main logic for one Delta_g: per odd-J state base, N_s and the correction split into GL nodes x < XC (dlo) and
    x >= XC (dhi); every state of the exact-skeleton set (odd J >= 11: base only)."""
    P6.GRID_OVERRIDE.clear(); P6.GRID_OVERRIDE.update(G)
    mods = {J: P6.model(J) for J in G}
    P6.GRID_OVERRIDE.clear()
    NG = NearGap(Dg); F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
    recs = []
    for h, DF, r0, ref in load_states(Dg):
        mfl, J = r0['mfl'], r0['J']
        if J % 2 == 0:
            continue
        u = 1/mfl**2; s = SIGMA[h]
        x = pn*np.sqrt(u)
        P = remk.Pders_vec(Lm*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
        pf2, wf2 = gl01(max(400, 2*J + 100))
        Js_, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js_, J)); B2, B3 = B2[i], B3[i]
        p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
        cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
        fac = (1 - x*x)**2 if h == '++' else (1 - x*x)**-2
        rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
        Ns = abs(B2 + cv2) + abs(B3 + cv3)
        base = B2 + B3 + cv2 + cv3 + rem
        dlo = dhi = 0.0
        if J in mods:
            m = r0['m']; Gm = mods[J][0](x, m)
            rho2, rho3 = r0['rho'][2], r0['rho'][3]
            r20 = np.polyfit(x[:3], rho2[:3], 2)[-1]
            Dh = A.Dc('+-', J, x)
            v = (p2*A.env('+-', x, 2)*(Gm['C2']/m**2 - (rho2 - Dh*r20))
                 + p3*u*A.env('+-', x, 3)*(Gm['C3H']/m**2 - A.SIG['+-']*(rho3 - Dh*r20)))
            lo = x < XC
            dlo, dhi = v[lo].sum(), v[~lo].sum()
        recs.append(dict(mfl=mfl, J=J, base=base, Ns=Ns, dlo=dlo, dhi=dhi))
    return recs


def neargap(recs, t=1.0, tlow=1.0):
    lam = [(r['base'] + t*(tlow*r['dlo'] + r['dhi']))/r['Ns'] for r in recs]
    i = int(np.argmin(lam))
    perJ = {J: min(l for l, r in zip(lam, recs) if r['J'] == J) for J in sorted(set(r['J'] for r in recs))}
    argJ = {J: min(((l, r['mfl']) for l, r in zip(lam, recs) if r['J'] == J))[1] for J in perJ}
    return lam[i], recs[i]['J'], recs[i]['mfl'], perJ, argJ


def window(recs, which):
    """exact interval of t (whole correction: 'all'; GL nodes x < 0.075: 'low') keeping every odd-J state positive (linear in t)."""
    lo, hi = -np.inf, np.inf
    for r in recs:
        a = r['base'] + (r['dhi'] if which == 'low' else 0.0)
        b = r['dlo'] if which == 'low' else r['dlo'] + r['dhi']
        if b > 0:
            lo = max(lo, -a/b)
        elif b < 0:
            hi = min(hi, -a/b)
    return lo, hi


MAIN = "C^_2' (exact u_T) + n3 contour points  [MAIN]"

if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(str(s) + "\n")
    t0 = time.time()
    pr(f"# v4_odd100: odd sector at Delta_g = 100; P8_XROWS = {os.environ['P8_XROWS']}, P8_PHASE = {P8.PHASE}, RELMAX = {P8.RELMAX}")
    pr(f"  shadow phase check |Im logK - pi(Delta + J)| over Delta in (100, 2e4), J = 5, 7, 9: "
       f"{max(check_phase(J, (100.0, 20000.0)) for J in (5, 7, 9)):.1e}")
    out = {}
    for lab, sh, n3 in (("baseline: scalar g, no n3 points (= PARITY_ODD §11.4 main protocol)", _orig_shift, False),
                        ("scalar g + n3 contour points at x = 0.1, 0.15", _orig_shift, True),
                        ("C^_2' (exact u_T), no n3 points", make_shift(0.0), False),
                        (MAIN, make_shift(0.0), True)):
        G = grids(sh, n3)
        R = state_parts(G)
        out[lab] = (G, R)
        lm, J, mfl, perJ, argJ = neargap(R)
        tw = window(R, 'all'); tl = window(R, 'low')
        pr(f"\n## {lab}")
        pr(f"   odd-sector min F/N_s = {lm:+.4f} (J = {J}, m/Lam = {mfl}); per J: " +
           ", ".join(f"J={j}: {v:+.4f} ({argJ[j]})" for j, v in perJ.items()))
        pr(f"   t* (whole correction) = ({tw[0]:+.2f}, {tw[1]:+.2f}); correction on GL nodes x < 0.075 scaled by t: positive for "
           f"t in [{tl[0]:+.2f}, {tl[1]:+.2f}]")
        pr(f"   Delta = 100 row: x = " + " ".join(f"{x:.2f}" for x in G[5][100.0]['xs']) + "   (* = direct value)")
        for J in (5, 7, 9):
            g = G[J][100.0]
            for kd in ('C2', 'C3H'):
                d = g[(kd, 'direct')]
                pr(f"     J={J} {kd:3s} Gamma : " + " ".join(f"{v:+7.0f}{'*' if dd else ' '}" for v, dd in zip(g[kd], d)) + f"   (a = {g['a']:.5g})")
                if n3:
                    K2 = g[(kd, 'K2same_a')]
                    sel = [i for i, x in enumerate(g['xs']) if abs(x - 0.1) < 1e-6 or abs(x - 0.15) < 1e-6]
                    pr(f"           x = 0.10/0.15: direct {g[kd][sel[0]]:+.0f} +- {g[(kd, 'err')][sel[0]]:.0f} / {g[kd][sel[1]]:+.0f} +- "
                       f"{g[(kd, 'err')][sel[1]]:.0f};  K2 {K2[sel[0]]:+.0f} / {K2[sel[1]]:+.0f};  direct/K2 "
                       f"{g[kd][sel[0]]/K2[sel[0]]:.2f} / {g[kd][sel[1]]/K2[sel[1]]:.2f}")
    R = out[MAIN][1]
    for tt in (-1.0, 0.0, 2.0, 3.0):
        lm, J, mfl, _, _ = neargap(R, tlow=tt)
        pr(f"   MAIN, t_low = {tt:+.1f}: min F/N_s = {lm:+.4f} (J = {J}, m/Lam = {mfl})")
    pickle.dump({k: v[1] for k, v in out.items()}, open(os.path.join(HERE, 'cache', 'v4_odd100%s.pkl' % ('' if os.environ['P8_XROWS'] == '0' else '_xrows')), 'wb'))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v4_odd100%s.txt' % ('' if os.environ['P8_XROWS'] == '0' else '_xrows')), 'w').write(buf.getvalue())

"""
gaplarge/gl6_p4w.py (GAP_LARGESPIN.md §2.2): Lemma P4W on Theorem L''s range with the UNIFORM rate (AG_U), i.e. AG_PROOF §6.2's family
'U' at theta = 0.97 (the form AG_PROOF says Theorem L' needs), as a function of Delta_g at fixed gap-unit coordinates (beta, mu).  [DV]

Re-typed from agproof/ag_p4w.run_L + its 'U' contour (edge and top with the uniform rate theta Lam beta, rectangle height
25/(theta Lam beta); boundary side = p4win/p4w_contour.contour_bound 'local', unchanged), on the 16 Theorem L' states of
p4win/p4w_L.ROWS mapped to the same (beta, mu) at Dg (m = mu Lam, J = round(m sinh(beta/2) - 1/2)).  Lam enters ONLY through the module
attribute p4w_common.L, which this process sets per gap (read-only imports; nothing of theirs is written).
Premise forms transported with Dg-INDEPENDENT constants in AdS units [M]: Lam_amp = 3 x the scalar CHMRS O(1/m^2) beyond-hybrid term,
the scalar-calibrated edge envelope e_sc (+- transferred, hel = -2) with its AdS-unit coefficients (a0, a_bulk, a_edge) = (10, 4, 5).
Budget: half floor F/2 with F = (i) LRINF_TIGHT's F_cert (Dg = 1000 only: reproduces AG_PROOF's 0.133 / 18.9), (ii) the certified
Dg-uniform floor f_U(beta)/Lam^2 (gl2).
Output: cache/gl6_p4w.{txt,json}
"""
import os, sys, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
CODE = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(CODE, 'p4win'))
from fractions import Fraction as Fr
import numpy as np
from scipy.special import j0
import p4w_common as C          # read-only (C.L patched in this process)
import p4w_contour as K         # read-only
import p4w_L as PL              # read-only (ROWS, mfl_exact)

OUT = []
def pr(s=''):
    print(s, flush=True); OUT.append(s)

G2 = json.load(open(os.path.join(HERE, 'cache', 'gl2_allgaps.json')))
FL = [(float(Fr(f['a'])), float(Fr(f['b'])), f['fU']) for f in G2['floor']]
TH = 0.97


def fU(beta):
    for a, b, v in FL:
        if a <= beta <= b:
            return v
    return None


def _yrule(kl, T):      # agproof/ag_p4w._yrule (copied)
    br = sorted(set([0.0] + [min(T, c/kl) for c in (1, 2, 4, 6, 9, 13, 20, 30, 45, 70)] + [T]))
    xs, ws = np.polynomial.legendre.leggauss(20)
    ys, wy = [], []
    for a, b in zip(br[:-1], br[1:]):
        if b > a:
            ys += list(a + (b - a)*(xs + 1)/2); wy += list(ws*(b - a)/2)
    return np.array(ys), np.array(wy)


def contour_U(beta, mfl, pw, Mw, Mfun, k=2, u=1.0, Aunif=0.0):
    """ag_p4w.contour_bound for family 'U' (theta = TH), mode 'local'."""
    L = C.L
    bd, _, _ = K.contour_bound(beta, mfl, pw, Mw, Mfun, 'local', k=k, u=u, Aunif=Aunif)
    w = 2.0 if k == 2 else u
    k0 = 0.97*beta*L; kl = k0/np.sqrt(1 - (1/mfl)**2)
    T = 25.0/(TH*beta*L)
    ys, wy = _yrule(kl, T)
    ps = np.abs(K.psi_c(k, 1 + 1j*ys))
    Vp = float((ps*np.exp(-TH*beta*L*ys)*wy).sum())
    VA = K.V(k, 1.0, 0.97*beta*L, T)
    M1 = float(Mfun(np.array([1.0]))[0])
    ed = 2*w*(Vp*(M1 + Aunif) + VA*Aunif)
    x, wq = np.polynomial.legendre.leggauss(400); p = pw + (1 - pw)*(x + 1)/2; wq = wq*(1 - pw)/2
    et = np.exp(-TH*beta*L*T)*np.ones_like(p)
    tp = 2*w*float((np.abs(K.psi_c(k, p + 1j*T))*Mfun(p)*et*wq).sum())
    return bd, ed, tp


def state_cost(J, mfl, La=3.0):
    L = C.L
    m = mfl*L; beta = C.beta_of(J, mfl); u = 1/mfl**2
    nuw = 3.0/beta; pw = nuw/L
    nu = np.linspace(0, nuw, 4001)[1:]; dn = nu[1] - nu[0]
    P, dPb = C.Pd(beta, nu); P = P.real
    bb = (2*J + 1)/mfl*nu/L
    B2 = (C.T(1, beta, nu, P, dPb)[0] - nu**2*C.gflat(1, bb))/m**2
    S4 = (C.T(2, beta, nu, P, dPb)[0] - nu**4*(C.gflat(2, bb) + (3/16)*j0(bb)))/m**4
    ps2 = C.psi2(nu/L); ps3 = C.psi3(nu/L)
    act2 = lambda f: 2*float((ps2*f).sum())*dn/L
    w2 = abs(act2(B2)) + u*float((np.abs(ps3)*np.abs(B2)).sum())*dn/L
    w4 = abs(act2(S4)) + u*float((np.abs(ps3)*np.abs(S4)).sum())*dn/L
    nn = np.linspace(0.7*nuw, 1.3*nuw, 41); Pn, dPn = C.Pd(beta, nn); bn = (2*J + 1)/mfl*nn/L
    B2n = (C.T(1, beta, nn, Pn.real, dPn)[0] - nn**2*C.gflat(1, bn))/m**2
    S4n = (C.T(2, beta, nn, Pn.real, dPn)[0] - nn**4*(C.gflat(2, bn) + (3/16)*j0(bn)))/m**4
    Mw = 2*float(np.max(np.abs(La*B2n + S4n)))
    env = K.e_sc(beta, mfl, hel=-2.0)
    Au = K.A_curv(beta, J, mfl)
    bd2, ed2, tp2 = contour_U(beta, mfl, pw, Mw, env, k=2, Aunif=Au)
    bd3, ed3, tp3 = contour_U(beta, mfl, pw, Mw, env, k=3, u=u, Aunif=Au)
    win = La*w2 + w4
    edge = ed2 + ed3 + tp2 + tp3
    return dict(beta=beta, win=win, bd=bd2 + bd3, edge=edge, d=win + bd2 + bd3 + edge)


def main():
    pr(__doc__.split('Output')[0])
    L0 = np.sqrt(998.0*999.0)
    base = []
    for lab, J, mtab, Fn, Fc in PL.ROWS:
        mfl = PL.mfl_exact(lab, J, mtab)
        base.append((lab, J, 2*np.arcsinh((J + 0.5)/(mfl*L0)), mfl, Fc))
    res = {}
    for Dg in (1000, 2000, 5000):
        C.L = np.sqrt((Dg - 2.0)*(Dg - 1.0))
        rows = []
        pr(f'\n## Dg = {Dg}, Lam = {C.L:.3f}   (state | J | m/Lam | beta | d/(F/2) [window, boundary, edge+top] | edge tolerance '
           f'(multiple of the transferred scalar envelope) ; unif-cert floor' + (' ; and with LRINF F_cert' if Dg == 1000 else '') + ')')
        for lab, J0, beta0, mfl, Fc in base:
            if Dg != 1000 and beta0 > 3.1:
                continue                       # Mehler node count ~ 1.5 Lam beta: too costly / memory (J1000 = 1e5, 1e6 at 1000 only)
            J = int(round(mfl*C.L*np.sinh(beta0/2) - 0.5))
            c = state_cost(J, mfl)
            Fu = fU(c['beta'])/C.L**2
            d_u = c['d']/(Fu/2); tol_u = (Fu/2 - (c['d'] - c['edge']))/c['edge']
            row = dict(label=lab, J0=J0, J=J, mfl=mfl, beta=c['beta'], Fu=Fu, d_u=d_u, tol_u=tol_u,
                       parts=[c['win']/(Fu/2), c['bd']/(Fu/2), c['edge']/(Fu/2)])
            s = (f'   {lab:10s} {J:8d} {mfl:6.3f} {c["beta"]:6.3f} | d/(F_U/2) = {d_u:.3e} [{row["parts"][0]:.1e}, {row["parts"][1]:.1e}, '
                 f'{row["parts"][2]:.1e}] | tol {tol_u:.3g}')
            if Dg == 1000:
                row['d_c'] = c['d']/(Fc/2); row['tol_c'] = (Fc/2 - (c['d'] - c['edge']))/c['edge']
                s += f'   || F_cert(1000): d = {row["d_c"]:.3f}, tol {row["tol_c"]:.3g}'
            pr(s)
            rows.append(row)
        res[Dg] = rows
        pr(f'   max d/(F_U/2) = {max(r["d_u"] for r in rows):.4f};  min edge tolerance = {min(r["tol_u"] for r in rows):.3g}'
           + (f';  with LRINF F_cert: max d = {max(r["d_c"] for r in rows):.3f}, min tol = {min(r["tol_c"] for r in rows):.3g} '
              f'(AG_PROOF §6.2: 0.133 / 18.9)' if Dg == 1000 else ''))
    def get(Dg, r0):
        return next((r for r in res[Dg] if r['label'] == r0['label'] and r['J0'] == r0['J0']), None)
    pr('\n## ratio d(1000)/d(Dg) at the same (beta, mu) [unif floor]: ' + '; '.join(
        f'{r0["label"]} J1000={r0["J"]}: ' + ', '.join(f'{Dg}: {r0["d_u"]/get(Dg, r0)["d_u"]:.2f}' for Dg in (2000, 5000))
        for r0 in res[1000] if get(2000, r0) is not None))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE, 'cache', 'gl6_p4w.json'), 'w'), default=float, indent=0)
    open(os.path.join(HERE, 'cache', 'gl6_p4w.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()

"""v5 (H7_I14_LIGHTSIDE.md §24.3-§24.4): size of the null directions X_h of C^_2' on the near-gap states, with the existing chain data.

By Lemma N (§24.1, v5_null.py) a null direction changes the heavy action of a state only through the part of its forward-normalised
Psi^_nu[state] that is NOT of the form A^_h[Y] with Y even and analytic on |Im nu| <= 1 (the leading-Regge kernel K_nu(cosh beta), the
graviton pole, any entire/strip-analytic C-level density drop out exactly).  Whether the beyond-leading-Regge part of Psi^ is of that form
is [I] (§24.2).  Here we bound the effect by treating the beyond-leading-Regge part as GENERIC (a worst case: nothing of it in the range):
  M1: delta(nu) = env_{2,h}(x) [Gamma_{2,h}(x; m)/m^2 - (rho_2 - D_h rho_2(0))] + i env_{2,h} Gamma^im_{2,h}(x; m)/m   (the chain's
      graviton-specific correction, k = 2, C^_2' local de-shift, incl. its imaginary O(1/m) part, which Re no longer removes when the
      weight is complex), off-diagonal env_off (Gamma_off + i Gamma_off^im)/m;
  M2: M1 + [env D_h(x) - 2 J0((2J+1)x)]  (also the crystal-vs-Regge leading-order difference treated as generic: envelope);
  B : the bound int |X_h| |delta_M1|.
Extrapolation of Gamma below the smallest row node x0 (0.05 on the Delta = 100 row): 'chain' (the chain's own natural-spline extrapolation),
'hold' (Gamma(x0)), 'lin0' (Gamma(x0) x/x0).
Pairing (per state, X_h normalised to max |X_h| = max_p psi_2(p)):  P_h = int_0^1 dp X_h(Lam p) delta(p);  Delta F(c) = Re(c P_h).
Reported: max over states of the spectral norm of N_s^{-1/2} Re(e^{i th} P_h) N_s^{-1/2} (max over th) -- the largest possible change of
lam_min/N_s per unit coefficient |c| = 1 -- and the same for the gauge differences between natural conventions (min-norm in the X-, chi- and
psi-relative normalisations; v5 §24.3) and for the low-nu weight-minimising gauge of §24.4.
Chains (read-only imports, nothing written in spinfunc/):
  even: l9 rows (K2) with direct n1 rows (x >= 0.2) and the n3 contour points (x = 0.1, 0.15) at Delta = 100 (as n3_analysis with
        N3_SPIN=1, but with the EXACT u_T of v4_hclosed); Delta_g = 300, 1000: l9 Laplace rows (as v3/v4).
  odd : Delta_g = 100: p8 grid + n3 odd points (v4_odd100 MAIN; P8_XROWS = 0 or 1); 300, 1000: p4 grid (v4_odd_impact).
Usage: python3 v5_null_heavy.py even            -> cache/v5_null_heavy_even.txt
       P8_XROWS=0|1 python3 v5_null_heavy.py odd -> cache/v5_null_heavy_odd_x<0|1>.txt
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
SECTOR = sys.argv[1] if len(sys.argv) > 1 else 'even'
if SECTOR == 'odd':
    os.environ['P8_XROWS'] = os.environ.get('P8_XROWS', '0')
    os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
sys.dont_write_bytecode = True
import numpy as np
from scipy.special import j0
from scipy.interpolate import CubicSpline
from v4_chi import u_spin
from spf_common import NearGap, load_states
from common import gl01, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk

buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


NUL = pickle.load(open(os.path.join(HERE, 'cache', 'v5_null.pkl'), 'rb'))
NUMAX = 60.0
selN = NUL['nu'] <= NUMAX
NU, DNU, VH = NUL['nu'][selN], NUL['dnu'], NUL['V'][:, selN]
UT = u_spin(NU)
THETAS = np.linspace(0, np.pi, 13)[:-1]
EXTRAP = ('chain', 'hold', 'lin0')
MODELS = ('M1', 'M1r', 'Mc', 'B')


def legP_series(lam, eta, kmax=160):
    """P_lam(eta) = 2F1(-lam, lam+1; 1; (1-eta)/2), vectorised in complex lam, eta >= 1 close to 1 (|z| << 1)."""
    z = (1 - eta)/2
    t = np.ones_like(lam, dtype=complex); s = t.copy()
    for k in range(kmax):
        t = t*(-lam + k)*(lam + 1 + k)/((k + 1)**2)*z
        s = s + t
        if np.max(np.abs(t)) < 1e-17*np.max(np.abs(s)):
            break
    return s


def exact_shift(A, orig):
    """l9/p4 'shift' with the exact spinning u_T for k = 2 (k = 3 keeps the scalar g), as v4_odd_impact.make_shift(0)."""
    def sh(J, h, k, xs, m):
        if k != 2:
            return orig(J, h, k, xs, m)
        f = lambda t: A.env(h, t, k)*A.Dc(h, J, t)
        return np.array([(0.5*((1 + complex(u_spin(x*m)))*f(x + 1j/m) + (1 - complex(u_spin(x*m)))*f(x - 1j/m)) - f(x)).real
                         for x in xs])
    return sh


def rowmodel(rows, keys):
    """f(x, m, mode) -> {key: values}: natural cubic spline in x per Delta row (held beyond the last node), linear in 1/m between rows
    (as l9_neargap.model / p6.model); below the first node x0: mode 'chain' (spline extrapolation, = the chains), 'hold', 'lin0'."""
    Ds = sorted(rows); ms = np.array([rows[D]['m'] for D in Ds])
    spl = {D: {k: CubicSpline(rows[D]['xs'], rows[D][k], bc_type='natural') for k in keys if k in rows[D]} for D in Ds}
    x0 = {D: rows[D]['xs'][0] for D in Ds}; x1 = {D: rows[D]['xs'][-1] for D in Ds}
    def one(D, k, x, mode):
        s = spl[D][k]; v = s(np.minimum(x, x1[D]))
        lo = x < x0[D]
        if mode == 'hold':
            v = np.where(lo, s(x0[D]), v)
        elif mode == 'lin0':
            v = np.where(lo, s(x0[D])*x/x0[D], v)
        return v
    def f(x, m, mode='chain'):
        if m <= ms[0]:
            i0 = i1 = 0; t = 0.0
        elif m >= ms[-1]:
            i0 = i1 = len(Ds) - 1; t = 0.0
        else:
            i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
            t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
        return {k: (1 - t)*one(Ds[i0], k, x, mode) + t*one(Ds[i1], k, x, mode) for k in keys if k in spl[Ds[i0]]}
    return f


def gauges(Dg):
    """coefficients c (complex 3-vector; X_h normalised to max|X_h| = max psi_2) of natural gauge changes relative to the min-norm gauge X_0:
    'chi' : min || w (X - psi_2) ||, w = (|h_-| + |h_+|)/|h_- + h_+|   (min-norm in the chi-normalisation);
    'rel' : min || (X - psi_2)/(psi_2 + 0.05 max psi) ||               (relative deviation);
    'low9': min int_{nu < 9} |X|^2                                    (least weight below the validated contour range);
    'n0r' : min int_{4.5 < nu < 7} |X|^2 + int_{elsewhere} |X - X_0|^2  (least weight at the x = 0.05 node, nu = 4.9-6.4, without
            blowing the weight up elsewhere);
    'n0'  : min int_{4.5 < nu < 7} |X|^2 alone (unconstrained: needs |c| ~ 1e3-1e6, shown only to exclude it)."""
    d = pickle.load(open(os.path.join(HERE, 'cache', 'v4_chi_%d.pkl' % Dg), 'rb'))
    nu, X0, ph = d['nu'], d['X'], d['phi']
    s = nu <= NUMAX; nu, X0, ph = nu[s], X0[s], ph[s]
    mx = d['phi'].max()
    Vg = np.array([np.interp(nu, NU, v.real) + 1j*np.interp(nu, NU, v.imag) for v in VH])*mx
    import v5_null as V5
    hm, hp = V5.F_['RmT'](nu), V5.F_['RpT'](nu)
    w = (np.abs(hm) + np.abs(hp))/np.abs(hm + hp)
    win = (nu > 4.5) & (nu < 7)
    out = {}
    for lab in ('chi', 'rel', 'low9', 'n0r', 'n0'):
        if lab == 'chi':
            r0, M = w*(X0 - ph), (w[None, :]*Vg).T
        elif lab == 'rel':
            wt = 1/(ph + 0.05*mx); r0, M = wt*(X0 - ph), (wt[None, :]*Vg).T
        elif lab == 'low9':
            r0, M = X0[nu < 9], Vg[:, nu < 9].T
        elif lab == 'n0':
            r0, M = X0[win], Vg[:, win].T
        else:
            r0 = np.where(win, X0, 0); M = Vg.T
        c, *_ = np.linalg.lstsq(M, -r0, rcond=None)
        red = np.linalg.norm(r0 + M @ c)/np.linalg.norm(r0)
        if lab == 'n0r':
            red = np.linalg.norm((X0 + c @ Vg)[win])/np.linalg.norm(X0[win])
        out[lab] = (c, red)
    prof = {}
    for lab in ('low9', 'n0r'):
        Xg = X0 + out[lab][0] @ Vg
        prof[lab] = [(v, X0[np.argmin(abs(nu - v))]/mx, Xg[np.argmin(abs(nu - v))]/mx, ph[np.argmin(abs(nu - v))]/mx)
                     for v in (1.0, 2.0, 3.0, 4.4, 5.0, 5.5, 6.0, 6.4, 7.0, 9.0, 12.0, 15.0, 20.0)]
    return out, prof


def pairings(r, Lam, G, extra, pn):
    """complex pairings P[(k, e, model, mode)] = int dp X_k(Lam p) delta_e(p) of null vector k (max|X_k| = 1) with:
       M1 : env [Gamma/m^2 - (rho2 - D_h rho2(0))] + i [Im G/a - 2 Re(u_T) Im P_{l-1}(cosh beta)]  (generic graviton-specific part, with
            the imaginary part MINUS its exact leading-Regge value, which is in the range)
       M1r: real part of M1 only;  Mc: M1 restricted to the controlled data x >= 0.075 (delta = 0 below);  B: int |X_k| |delta_M1|."""
    m, mfl, J, hs = r['m'], r['mfl'], r['J'], r['hs']
    p = NU/Lam; x = p/mfl
    rho = np.interp(p, pn, r['rho2'])
    imLO = 2*UT.real*legP_series(-1.5 + 1j*NU, r['eta']).imag
    out = {}
    for mode in EXTRAP:
        g = G(x, m, mode)
        for e in hs:
            env = extra['env'](e, x); Dh = extra['Dsh'](e, J, x)
            re_ = env*(g[extra['key'](e)]/m**2 - (rho - Dh*r['r20']))
            im_ = env*g[extra['key'](e, 'im')]/m - imLO
            out[(e, 'M1', mode)] = re_ + 1j*im_
            out[(e, 'M1r', mode)] = re_ + 0j
            out[(e, 'Mc', mode)] = np.where(x >= 0.075, re_ + 1j*im_, 0)
        if extra.get('off') and len(hs) == 2:
            eo = np.sqrt(extra['env']('++', x)*extra['env']('+-', x))
            d = eo*(g[extra['key']('off')] + 1j*g[extra['key']('off', 'im')])/m
            out[('off', 'M1', mode)] = d; out[('off', 'M1r', mode)] = d.real + 0j; out[('off', 'Mc', mode)] = np.where(x >= 0.075, d, 0)
    P = {}
    for k, v in enumerate(VH):
        for key, d in out.items():
            P[(k,) + key] = np.sum(v*d)*DNU/Lam
            if key[1] == 'M1':
                P[(k, key[0], 'B', key[2])] = np.sum(np.abs(v*d))*DNU/Lam
    return P


def lam_min(F0, Ns, dF=None):
    n = 1/np.sqrt(np.asarray(Ns))
    M = F0 if dF is None else F0 + dF
    return np.linalg.eigvalsh(n[:, None]*M*n[None, :])[0]


def dFmat(P, k_c, hs, model, mode):
    """Re( sum_k c_k P_k ) as a matrix (diag hs, off if 2x2); k_c: complex 3-vector (units: max|X_h| = max psi)."""
    M = np.zeros((len(hs), len(hs)))
    for i, e in enumerate(hs):
        M[i, i] = sum(k_c[k]*P[(k, e, model, mode)] for k in range(3)).real
    if len(hs) == 2:
        M[0, 1] = M[1, 0] = sum(k_c[k]*P[(k, 'off', model, mode)] for k in range(3)).real
    return M


def run_sector(states_fn, Gmods, extra, Dg, label):
    NG = NearGap(Dg); F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
    mxpsi = np.max(F.psi2(np.linspace(0, 1, 20001)))
    recs = [r for r in states_fn(Dg, NG, Gmods)]
    gs, prof = gauges(Dg)
    res = []
    for r in recs:
        if r['J'] in Gmods:
            P = pairings(r, Lm, Gmods[r['J']], extra, pn)
            res.append((r, {k: v*mxpsi for k, v in P.items()}))
    lam0 = [lam_min(r['F0'], r['Ns']) for r in recs]
    i0 = int(np.argmin(lam0))
    pr(f"\n## {label}, Delta_g = {Dg} (Lam = {Lm:.2f}): {len(recs)} states ({len(res)} with chain corrections); X_h normalised to max|X_h| = "
       f"max psi_2 = {mxpsi:.3f}")
    pr(f"   chain baseline (local de-shift, exact u_T): min lam/N_s = {lam0[i0]:+.4f} [J = {recs[i0]['J']}, m/Lam = {recs[i0]['mfl']}]")
    pr("   (a) per unit coefficient |c| = 1 (max over phase): max over states of |Delta lam|/N_s [state]; extrapolation chain / hold / lin0")
    for k in range(3):
        for model in MODELS:
            line = []
            for mode in EXTRAP:
                vals = []
                for r, P in res:
                    if model == 'B':
                        v = (max(abs(P[(k, e, 'B', mode)])/n_ for e, n_ in zip(r['hs'], r['Ns'])) +
                             (abs(P[(k, 'off', 'M1', mode)])/np.sqrt(np.prod(r['Ns'])) if len(r['hs']) == 2 else 0.0))
                    else:
                        n = 1/np.sqrt(np.asarray(r['Ns'])); v = 0.0
                        for th in THETAS:
                            c = np.zeros(3, complex); c[k] = np.exp(1j*th)
                            M = dFmat(P, c, r['hs'], model, mode)
                            v = max(v, np.max(np.abs(np.linalg.eigvalsh(n[:, None]*M*n[None, :]))))
                    vals.append((v, r['J'], r['mfl']))
                vm = max(vals)
                line.append(f"{vm[0]:.2e} [J={vm[1]}, {vm[2]}]")
            pr(f"     X_{k + 1} {model:3s}: " + " / ".join(line))
    pr("   (b) a-priori gauge choices (state-independent) vs the min-norm gauge: |c_k|, objective ratio; new min lam/N_s over all states "
       "(M1 / M1r / Mc; chain extrap.)")
    for lab, (c, red) in gs.items():
        outs = []
        for model in ('M1', 'M1r', 'Mc'):
            lam = list(lam0)
            for j, r in enumerate(recs):
                hit = [P for rr, P in res if rr is r]
                if hit:
                    lam[j] = lam_min(r['F0'], r['Ns'], dFmat(hit[0], c, r['hs'], model, 'chain'))
            jm = int(np.argmin(lam))
            outs.append(f"{lam[jm]:+.4f} [J={recs[jm]['J']}, {recs[jm]['mfl']}]" + ("" if lam[jm] > 0 else " n!=1"))
        pr(f"     {lab:5s}: |c| = " + " ".join(f"{abs(v):.3g}" for v in c) + f"; ratio {red:.3f}; min lam/N_s: " + " / ".join(outs))
    for lab, pf in prof.items():
        pr(f"   gauge '{lab}' profile, nu: Re X_0 / Re X_gauge / psi_2 (units max psi): " +
           "  ".join(f"{v:.1f}: {a.real:+.3f}/{b.real:+.3f}/{c_:.3f}" for v, a, b, c_ in pf))
    # (c) leverage diagnostic (NOT a legitimate choice: tuned on the states): best min lam over |c_k| <= 0.25 (random search, M1r)
    rng = np.random.default_rng(1)
    best = (min(lam0), None)
    for _ in range(400):
        c = (rng.uniform(-1, 1, 3) + 1j*rng.uniform(-1, 1, 3))*0.25/np.sqrt(2)
        lm = min(lam_min(r['F0'], r['Ns'], dFmat(P, c, r['hs'], 'M1r', 'chain')) for r, P in res)
        lm = min([lm] + [l for l, r in zip(lam0, recs) if r['J'] not in Gmods])
        if lm > best[0]:
            best = (lm, c)
    pr(f"   (c) leverage (diagnostic only; tuned on the states, not admissible): max over |c_k| <= 0.25 of min lam/N_s (M1r) = {best[0]:+.4f}"
       f" (baseline {min(lam0):+.4f})")
    return res


def base_parts(Dg, NG, r0, h):
    F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
    mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]; x = pn*np.sqrt(u)
    P = remk.Pders_vec(Lm*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
    pf2, wf2 = gl01(max(400, 2*J + 100))
    Js_, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js_, J)); B2, B3 = B2[i], B3[i]
    p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
    cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
    fac = (1 - x*x)**2 if h == '++' else (1 - x*x)**-2
    rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
    return dict(base=B2 + B3 + cv2 + cv3 + rem, Ns=abs(B2 + cv2) + abs(B3 + cv3), x=x, p2=p2, p3=p3, u=u,
                r20=np.polyfit(x[:3], r0['rho'][2][:3], 2)[-1])


def even_setup():
    import l9_analysis as A
    import l9_neargap as NGm
    import n1_analysis as N1
    _orig = N1.dload
    def dload(kind, J):
        d = dict(_orig(kind, J))
        f = os.path.join(SF, 'cache', 'n3_pts_%s_J%d.pkl' % (kind, J))
        if os.path.exists(f):
            d.update(pickle.load(open(f, 'rb')))
        return d
    N1.dload = dload
    A.shift = exact_shift(A, A.shift)
    origgrid = A.grid
    kinds = ('C2', 'C3H')
    mods = {}
    for Dg in (100, 300):
        if Dg == 100:
            N1.DIRECT = True
            NGm.A.grid = lambda J, K, kinds=('C2', 'C3H'), **kw: N1.grid_direct(J, K, kinds, **kw)
        else:
            NGm.A.grid = origgrid
        mm = {}
        for J in range(0, 11, 2):
            rows = NGm.rows_for(J, 2, kinds)
            keys = [(kd, e) for kd in kinds for e in (('++', '+-', 'off') if J >= 4 else ('++',))]
            keys += [k + ('im',) for k in keys]
            mm[J] = rowmodel(rows, keys)
        mods[Dg] = mm
        NGm.A.grid = origgrid
    mods[1000] = mods[300]
    def states(Dg, NG, Gm):
        blocks = {}
        for h, DF, r0, ref in load_states(Dg):
            b = blocks.setdefault((round(r0['mfl'], 9), r0['J'], DF), dict(mfl=r0['mfl'], J=r0['J'], m=r0['m'], eta=r0['eta'],
                                                                           rho2=r0['rho'][2], rho3=r0['rho'][3], parts={}))
            b['parts'][h] = base_parts(Dg, NG, r0, h)
        out = []
        for b in blocks.values():
            hs = [h for h in ('++', '+-') if h in b['parts']]
            e = b['parts'][hs[0]]; J = b['J']
            b['hs'] = hs; b['Ns'] = [b['parts'][h]['Ns'] for h in hs]; b['r20'] = e['r20']
            F0 = np.diag([b['parts'][h]['base'] for h in hs])
            if J in Gm and J % 2 == 0:
                x, p2, p3, u, m = e['x'], e['p2'], e['p3'], e['u'], b['m']
                G = Gm[J](x, m, 'chain')
                for kk, h in enumerate(hs):
                    Dh = A.Dsh(h, J, x)
                    F0[kk, kk] += (p2*A.env(h, x, 2)*(G[('C2', h)]/m**2 - (b['rho2'] - Dh*e['r20']))).sum()
                    F0[kk, kk] += (p3*u*A.env(h, x, 3)*(G[('C3H', h)]/m**2 - A.SIG[h]*(b['rho3'] - Dh*e['r20']))).sum()
                if len(hs) == 2:
                    eo2 = np.sqrt(A.env('++', x, 2)*A.env('+-', x, 2))
                    F0[0, 1] = F0[1, 0] = ((p2*eo2*G[('C2', 'off')]).sum() + (p3*u*eo2/(2 - x*x)*G[('C3H', 'off')]).sum())/m
            b['F0'] = F0
            out.append(b)
        return out
    extra = dict(env=lambda e, x: A.env(e, x, 2), Dsh=A.Dsh, key=lambda e, *im: ('C2', e) + im, off=True)
    return mods, states, extra


def odd_setup():
    import v4_odd100 as O
    import p4_odd_analysis as A4
    from v4_odd_impact import make_shift, _orig_shift
    G100 = O.grids(make_shift(0.0), True)
    A4.shift = make_shift(0.0)
    G34 = {J: A4.grid(J, K=2, kinds=('C2', 'C3H')) for J in (5, 7, 9)}
    A4.shift = _orig_shift
    keys = ['C2', 'C3H', ('C2', 'im')]
    mods = {100: {J: rowmodel({D: r for D, r in G100[J].items() if D <= 3200.0}, keys) for J in G100},
            300: {J: rowmodel({D: r for D, r in G34[J].items() if D <= 3200.0}, keys) for J in G34}}
    mods[1000] = mods[300]
    def states(Dg, NG, Gm):
        out = []
        for h, DF, r0, ref in load_states(Dg):
            J = r0['J']
            if J % 2 == 0:
                continue
            e = base_parts(Dg, NG, r0, h)
            b = dict(mfl=r0['mfl'], J=J, m=r0['m'], eta=r0['eta'], rho2=r0['rho'][2], hs=['+-'], Ns=[e['Ns']], r20=e['r20'])
            F0 = e['base']
            if J in Gm:
                x, p2, p3, u, m = e['x'], e['p2'], e['p3'], e['u'], r0['m']
                G = Gm[J](x, m, 'chain'); Dh = A4.Dc('+-', J, x)
                F0 += (p2*A4.env('+-', x, 2)*(G['C2']/m**2 - (r0['rho'][2] - Dh*e['r20']))
                       + p3*u*A4.env('+-', x, 3)*(G['C3H']/m**2 - A4.SIG['+-']*(r0['rho'][3] - Dh*e['r20']))).sum()
            b['F0'] = np.array([[F0]])
            out.append(b)
        return out
    extra = dict(env=lambda e, x: A4.env('+-', x, 2), Dsh=lambda e, J, x: A4.Dc('+-', J, x),
                 key=lambda e, *im: 'C2' if not im else ('C2', 'im'))
    return mods, states, extra


if __name__ == '__main__':
    t0 = time.time()
    pr(f"# v5_null_heavy [{SECTOR}{'' if SECTOR == 'even' else ', P8_XROWS = ' + os.environ['P8_XROWS']}]: null directions of C^_2' on the "
       f"near-gap states (§24.3); generic-remainder models M1/M1r/Mc and bound B (see module doc)")
    mods, states, extra = even_setup() if SECTOR == 'even' else odd_setup()
    pr(f"  chains built [{time.time() - t0:.0f}s]")
    for Dg in (100, 300, 1000):
        run_sector(states, mods[Dg], extra, Dg, f"{SECTOR} sector")
    pr(f"[{time.time() - t0:.0f}s]")
    tag = SECTOR if SECTOR == 'even' else 'odd_x' + os.environ['P8_XROWS']
    open(os.path.join(HERE, 'cache', 'v5_null_heavy_%s.txt' % tag), 'w').write(buf.getvalue())

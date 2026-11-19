"""
y2 (EE_DG100.md): per-state positivity n of Omega^EE = int psi2 C'_EE + int psi3 C3 (NF02 weights) at Delta_g = 100 -- the near-gap chain
of S1_S4 (q3_analysis: K2 linear chain, rows linear in 1/m, aligned off-diagonals, c_off = m Gamma_off model, C3 without de-shift,
phase strip logK -> logK - i pi (Delta + J), combined even + odd pass over all exact-skeleton states), with
  * the C^_2 crystal rows of EE (y1 + w1/u5 caches, read-only) in place of q^2's; C3H (o_H q) rows unchanged;
  * DIRECT values on the Delta = 100 row wherever computed: EE (y4: n1 x >= 0.2, n3 x = 0.1, 0.15; y5: n4 thimbles x = 0.05), C3H (n1 /
    n3 project caches; y5 at x = 0.05 even; n4 odd); odd: p8 (x >= 0.2), n3_odd, n4_odd (C3H) + y4/y5 (EE);
  * the k = 2 term evaluated in two ways:
      LOCAL : psi2 x [G/a_h - (A_h[f_LO] - f_LO)] (the chains' local de-shift (1 +- u_EE)/2; singular at the real zeros 1.17, 6.64 of
              h_- + h_+, reported for comparison only where nodes sit near them),
      NONLOCAL (default): Re int X [G/a_h - A_h[f_LO]] + int psi2 f_LO, with X = a Y the Abel-form weight of u3 (min-norm gauge) or of
              x2 (canonical gauge), evaluated in the chi-normalisation  X A_h[f] = Y (h_-^ f(nu+i) + h_+^ f(nu-i)),  h^ = R/n, so nothing is
              divided by h_- + h_+ (Gamma_NL := m^2 a(nu) (G/a_h - A_h[f])/env is smooth through the zeros; it is what is interpolated).
    Im(G) at K2 nodes with nu < Y2_NUIM (default 20) is replaced by the exact leading-Regge value (H7_I14 §24.3 side finding: K2 Im is
    4-7x too large there); direct nodes keep their Im.
  * the q^2 chain through the SAME code (kind C2q: q^2 rows, its own direct values incl. y5's n4 C2 at x = 0.05, local de-shift with the
    exact u_T, and nonlocal with u3's q^2 min-norm weight) for comparison.
Variants via env: Y2_GAUGE = min|can, Y2_XROWS = 0|1 (odd extra rows m = 1.05 ... 1.3 Lam, p8_k2 keys + their direct values),
Y2_EXTRAP = spline|hold|zero (below the first node x = 0.05), Y2_JHI = 1 (even J = 12 ... 28 / odd 11, 13 from the 300/1000 sweep rows,
held below the lowest row [M]), Y2_DIRECT = 0 (K2 only), Y2_MC = N (Gaussian redraws of the direct nodes, replica errors).
Usage: python3 y2_chain100.py [tag]  -> cache/y2_chain100_<tag>.txt, cache/y2_recs_<tag>.pkl
"""
import os, sys, io, glob, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.special import j0
from scipy.interpolate import CubicSpline
import l9_analysis as A
import l8_analysis as L8
import p4_odd_analysis as P4
import p2_shadow_odd as P2
from spf_common import NearGap, load_states
from common import gl01, SIGMA
from nf_common import flat_parts_nf
from adsaction import envk
import remk
import u2_hclosed_struct as U2
from v4_chi import u_spin as u_T

C = os.path.join(HERE, 'cache')
DREF = 20000.0
DG = 100
GAUGE = os.environ.get('Y2_GAUGE', 'min')
XROWS = os.environ.get('Y2_XROWS', '0') == '1'
EXTRAP = os.environ.get('Y2_EXTRAP', 'spline')
JHI = os.environ.get('Y2_JHI', '0') == '1'
DIRECT = os.environ.get('Y2_DIRECT', '1') == '1'
XC = float(os.environ.get('Y2_XC', '0.075'))
TLOW = [1.0]
BELOW = os.environ.get('Y2_BELOW', 'hold')
NUIM = float(os.environ.get('Y2_NUIM', '20'))
GERR = float(os.environ.get('Y2_GERR', '600'))       # direct node accepted if its Gamma-unit (c2) error <= GERR ...
RELACC = float(os.environ.get('Y2_RELACC', '0.5'))   # ... or its relative ratio error <= RELACC
NUFIT = 15.0
XFIT = 0.25
SIG = {'++': 1.0, '+-': -1.0}
LOG = []
_F = U2.numfuncs()


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


# ------------------------------------------------------------------------------------------------ magic-block coefficients
def hhat(struct, nu):
    """(R_-/n, R_+/n) in u3's chi-normalisation n = (|R_-| + |R_+|) e^{i arg(R_- + R_+)_{q^2}}; a = sum."""
    nu = np.asarray(nu, complex)
    Rm, Rp = _F[struct][0](nu), _F[struct][1](nu)
    qm, qp = _F['q2'][0](nu), _F['q2'][1](nu)
    nq = qm + qp
    n = (np.abs(Rm) + np.abs(Rp))*nq/np.abs(nq)
    return Rm/n, Rp/n


def u_of(struct, nu):
    if struct == 'q2':
        return complex(u_T(nu))
    Rm, Rp = _F[struct][0](complex(nu)), _F[struct][1](complex(nu))
    return (Rm - Rp)/(Rm + Rp)


def weight(struct, gauge):
    """(nu grid, Y, a) of the Abel-form weight (X = a Y)."""
    if struct == 'EE' and gauge == 'can':
        d = pickle.load(open(os.path.join(C, 'x2_chi_EE_can_100.pkl'), 'rb'))
    else:
        d = pickle.load(open(os.path.join(C, 'u3_chi_%s_100.pkl' % struct), 'rb'))
    hm, hp = hhat(struct, d['nu'])
    a = hm + hp
    Y = d['X']/np.where(np.abs(a) > 1e-12, a, 1e-12)
    if 'Y' in d:
        dev = np.max(np.abs(d['X'] - a*d['Y']))/np.max(np.abs(d['X']))
        LOG.append(f"weight {struct}/{gauge}: max|X - a Y|/max|X| = {dev:.1e} (a recomputed from u2)")
        Y = d['Y']
    return d['nu'], Y, d['X'], d.get('phi')


# ------------------------------------------------------------------------------------------------ records
def _pk(pats, odd=False):
    R = {}
    for p in pats:
        for fn in sorted(glob.glob(os.path.join(C, p))):
            d = pickle.load(open(fn, 'rb'))
            if odd:
                d = {(k[0], k[1]): v for k, v in d.items() if (len(k) < 3 or k[2] == 1.0)}
            R.update(d)
    return R


def strip(R, J):
    out = {}
    for k, r in R.items():
        if not all(np.all(np.isfinite(r[o])) for o in ('M0', 'M2', 'M4', 'logC', 'logK') if o in r):
            continue
        dev = (r['logK'].imag - np.pi*(k[0] + J)) % (2*np.pi); dev = min(dev, 2*np.pi - dev)
        assert dev < 1e-6, (k, J)
        out[k] = dict(r, logK=r['logK'] - 1j*np.pi*(k[0] + J))
    return out


def get(R, D, nu):
    c = [k for k in R if abs(k[0] - D) < 1e-6 and abs(k[1] - nu) < 1e-7*nu]
    return R[c[0]] if c else None


def even_k2(kind, J):
    if kind == 'EE':
        pats = ['y1_pts_EE_J%d_k2.pkl', 'u5_pts_EE_J%d_k2*.pkl', 'w1_pts_EE_J%d_k2*.pkl']
    elif kind == 'C3H':
        pats = ['l9_pts_C3H_J%d_k2.pkl', 'q2_pts_C3H_J%d_k2*.pkl', 'w1_pts_C3H_J%d_k2*.pkl']
    else:   # C2q (q^2)
        pats = ['l9_pts_C2_J%d_k2.pkl', 'l9_pts_C2_J%d.pkl', 'q2_pts_C2_J%d_k2*.pkl', 'w1_pts_C2_J%d_k2*.pkl']
    pats = pats + ['y1_pts_%s_J%%d_k2.pkl' % ('C2' if kind == 'C2q' else kind)] if kind != 'EE' else pats
    R = _pk([p % J for p in pats])
    if not XROWS:
        R = {k: v for k, v in R.items() if abs(k[0] - round(k[0])) < 1e-9}
    if kind == 'C2q' and J <= 10:
        R0 = dict(L8.records(J)); R0 = {k: dict(v, M6=0*v['M0']) if 'M6' not in v else v for k, v in R0.items()}
        R0.update(R); R = R0
    return strip(R, J)


def odd_k2(kind, J):
    src = {'EE': ['p3_pts_EE_never'], 'C3H': ['p3_pts_C3H_J%d.pkl', 'q2_p3_C3H_J%d_*.pkl', 'w1_p3_C3H_J%d*.pkl'],
           'C2q': ['p3_pts_C2_J%d.pkl', 'q2_p3_C2_J%d_*.pkl', 'w1_p3_C2_J%d*.pkl']}[kind]
    if kind == 'EE':
        src = ['y1_p3_EE_J%d.pkl', 'w1_p3_EE_J%d*.pkl']
    R = _pk([p % J if '%d' in p else p for p in src], odd=True)
    if XROWS and kind != 'EE':
        R.update(_pk(['p8_k2_%s_J%d.pkl' % ('C2' if kind == 'C2q' else kind, J)], odd=True))
    if not XROWS:
        R = {k: v for k, v in R.items() if abs(k[0] - round(k[0])) < 1e-9 and k[0] in (100.0, 200.0, 400.0, 800.0, 1600.0, 3200.0, DREF)}
    return strip(R, J)


def even_direct(kind, J):
    """(Delta, nu) -> (record, source) with n1-format E matrices."""
    if not DIRECT:
        return {}
    src = {'EE': [('y4_n1_pts_EE_J%d.pkl', 'n1'), ('y4_n3_pts_EE_J%d.pkl', 'n3'), ('y5_pts_EE_J%d.pkl', 'n4')],
           'C3H': [('n1_pts_C3H_J%d.pkl', 'n1'), ('n3_pts_C3H_J%d.pkl', 'n3'), ('y5_pts_C3H_J%d.pkl', 'n4')],
           'C2q': [('n1_pts_C2_J%d.pkl', 'n1'), ('n3_pts_C2_J%d.pkl', 'n3'), ('y5_pts_C2_J%d.pkl', 'n4')]}[kind]
    out = {}
    for p, s in src:
        f = os.path.join(C, p % J)
        if os.path.exists(f):
            for k, v in pickle.load(open(f, 'rb')).items():
                out[(float(k[0]), k[1])] = (v, s)
    return out


def odd_direct(kind, J):
    if not DIRECT:
        return {}
    src = {'EE': [('y4_p8_pts_EE_J%d.pkl', 'p8'), ('y4_n3_odd_pts_EE_J%d.pkl', 'n3'), ('y5_odd_pts_EE_J%d.pkl', 'n4')],
           'C3H': [('p8_pts_C3H_J%d.pkl', 'p8'), ('n3_odd_pts_C3H_J%d.pkl', 'n3'), ('n4_odd_pts_C3H_J%d.pkl', 'n4')],
           'C2q': [('p8_pts_C2_J%d.pkl', 'p8'), ('n3_odd_pts_C2_J%d.pkl', 'n3'), ('n4_odd_pts_C2_J%d.pkl', 'n4')]}[kind]
    out = {}
    for p, s in src:
        f = os.path.join(C, p % J)
        if os.path.exists(f):
            for k, v in pickle.load(open(f, 'rb')).items():
                if not XROWS and abs(k[0] - 100.0) > 1e-6 and k[0] < 1e4:
                    continue
                out[(float(k[0]), k[1])] = (v, s)
    return out


SHK2 = {}
NUCUR = [None]
for _f in ['l9_shadow_k2.pkl'] + [os.path.basename(f) for f in sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl')))] + ['y1_shadow_k2.pkl']:
    if os.path.exists(os.path.join(C, _f)):
        SHK2.update(pickle.load(open(os.path.join(C, _f), 'rb')))


def ratio_even(ds, dr, B, Bref, D, J, c0, fk, sq):
    """direct state/reference values per helicity element (++, +-, off) as complex numbers with replica errors (paired = CRN when both
    records come from the same CRN run (n4), else state and reference spreads added)."""
    def H(r, Bm, DD):
        b, logK = SHK2[(float(DD), int(J))][0], SHK2[(float(DD), int(J))][1]
        logK = logK - 1j*np.pi*(DD + J)
        binv = np.linalg.inv(b); na = r['E'].shape[0]
        Hm = Bm['Ti'] @ (r['E'][:na] @ binv.T) @ Bm['Ti'].T
        Hr = np.array([Bm['Ti'] @ (rep[:na] @ binv.T) @ Bm['Ti'].T for rep in r['reps']])
        lN = r['logN']
        if not r.get('ljc', False):            # n1 points stored before the log-Jacobian constant fix (n1_analysis.Hrep)
            import l1_crystal as Cr
            c, _ = Cr.saddle(r['nu'] if 'nu' in r else NUCUR[0], float(DD))
            lN = lN + np.sum(np.log([c[i] for i in (1, 2, 3, 4)]))
        return Hm, Hr, lN - logK
    (rs, ss), (rr, sr) = ds, dr
    Hs, Hsr, ls = H(rs, B, D); Hr, Hrr, lr = H(rr, Bref, DREF)
    e = np.exp(ls - lr - c0)
    n = min(len(Hsr), len(Hrr))
    crn = ss == 'n4' and sr == 'n4' and len(Hsr) == len(Hrr)
    out = {}
    for j, h in enumerate(('++', '+-')[:Hs.shape[0]]):
        v = e*Hs[j, j]/Hr[j, j]*fk
        if crn:
            reps = e*Hsr[:, j, j]/Hrr[:, j, j]*fk
            jk = np.array([e*np.delete(Hsr[:, j, j], i).mean()/np.delete(Hrr[:, j, j], i).mean()*fk for i in range(n)])
            err = np.sqrt((n - 1)/n*np.sum(np.abs(jk - jk.mean())**2))
        else:
            err = np.sqrt(np.var(e*Hsr[:n, j, j]/Hr[j, j]*fk)/n + np.var(e*Hs[j, j]/Hrr[:n, j, j]*fk)/n)
        out[h] = (v, err)
    if Hs.shape[0] == 2:
        v = e*0.5*(Hs[0, 1] + Hs[1, 0])/Hr[0, 0]*fk*sq
        so = 0.5*(Hsr[:n, 0, 1] + Hsr[:n, 1, 0])
        if crn:
            jk = np.array([e*np.delete(so, i).mean()/np.delete(Hrr[:, 0, 0], i).mean()*fk*sq for i in range(n)])
            err = np.sqrt((n - 1)/n*np.sum(np.abs(jk - jk.mean())**2))
        else:
            err = np.sqrt(np.var(e*so/Hr[0, 0]*fk*sq)/n + np.var(e*0.5*(Hs[0, 1] + Hs[1, 0])/Hrr[:n, 0, 0]*fk*sq)/n)
        out['off'] = (v, err)
    return out


def ratio_odd(ds, dr, D, J, c0, fk):
    def H(r, DD):
        b, logK, _ = P2.shadow_odd(float(DD), J)
        bi = 1.0/b[-1, -1]
        logK = logK - 1j*np.pi*(DD + J)
        return r['logN'] - logK, r['E'][0, 0]*bi, np.array([rep[0, 0]*bi for rep in r['reps']])
    (rs, ss), (rr, sr) = ds, dr
    ls, vs, rps = H(rs, D); lr, vr, rpr = H(rr, DREF)
    e = np.exp(ls - lr - c0)
    v = e*vs/vr*fk
    n = min(len(rps), len(rpr))
    if ss == 'n4' and sr == 'n4' and len(rps) == len(rpr):
        jk = np.array([e*np.delete(rps, i).mean()/np.delete(rpr, i).mean()*fk for i in range(n)])
        err = np.sqrt((n - 1)/n*np.sum(np.abs(jk - jk.mean())**2))
    else:
        err = np.sqrt(np.var(e*rps[:n]/vr*fk)/n + np.var(e*vs/rpr[:n]*fk)/n)
    return {'+-': (v, err)}


# ------------------------------------------------------------------------------------------------ grids
def Dsh(h, J, x):
    return A.Dsh(h, J, x) if J % 2 == 0 else P4.Dc(h, J, x)


def fLO(h, J, k, x):
    sg = SIG[h] if k == 3 else 1.0
    return sg*A.env(h, x, k)*(A.Dc(h, J, x) if J % 2 == 0 else P4.Dc(h, J, x))


def build_grid(J, kinds, struct2):
    """rows {Delta: rec}; kinds[0] is the k = 2 kind ('EE' or 'C2q'), kinds[1] = 'C3H'.  Per row: m, xs, a_h, per (kind, element):
    'loc' (Gamma, local de-shift / none for k = 3), 'nl' (k = 2: complex chi-normalised Gamma_NL), 'err' (Gamma units), 'src'."""
    odd = J % 2 == 1
    hels = ('+-',) if odd else (('++', '+-') if J >= 4 else ('++',))
    RK = {kd: (odd_k2(kd, J) if odd else even_k2(kd, J)) for kd in set(kinds) | {'C2q'}}
    DB = {kd: (odd_direct(kd, J) if odd else even_direct(kd, J)) for kd in kinds}
    k2kind = kinds[0]
    if not odd and J >= 4:
        Rq = dict(RK['C2q'])
        Bref = L8.basis(J, DREF, Rq)
        if Bref is None:
            Bref = L8.basis(J, DREF, RK[k2kind])
    Ds = sorted(set(k[0] for k in RK[k2kind] if k[0] < DREF and k[0] <= 3200.0))
    out = {}
    for D in Ds:
        m = mphys(D, J); mr = mphys(DREF, J)
        if not odd and J >= 4:
            B = L8.basis(J, D, RK['C2q'])
            if B is None:
                B = L8.basis(J, D, RK[k2kind])
            if B is None:
                continue
        else:
            B = Bref = dict(Ti=np.eye(1))
        xs = sorted(set(round(k[1]/m, 6) for k in RK[k2kind] if abs(k[0] - D) < 1e-6))
        xs = [x for x in xs if all(get(RK[kd], D, x*m) is not None and get(RK[kd], DREF, x*m) is not None for kd in kinds)]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        c0 = None
        rows = {}
        for kd in kinds + (('C2q',) if (not odd and J >= 4 and k2kind != 'C2q') else ()):
            k = 3 if kd == 'C3H' else 2
            G = {h: [] for h in hels}; Ge = {h: [] for h in hels}; src = []
            Go, Goe = [], []
            for x in xs:
                nu = x*m; st = get(RK[kd], D, nu); rf = get(RK[kd], DREF, nu); xr = nu/mr
                if st is None or rf is None:
                    for h in hels:
                        G[h].append(np.nan); Ge[h].append(np.nan)
                    Go.append(np.nan); Goe.append(np.nan); src.append('-'); continue
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                if odd:
                    vK = {'+-': np.exp(dl)*A.pred([st[o][-1, -1] for o in ('M0', 'M2', 'M4', 'M6') if o in st] + [0]*(4 - sum(o in st for o in ('M0', 'M2', 'M4', 'M6'))),
                                                   [rf[o][-1, -1] for o in ('M0', 'M2', 'M4', 'M6') if o in rf] + [0]*(4 - sum(o in rf for o in ('M0', 'M2', 'M4', 'M6'))), 2)*fk}
                    sq = 1.0
                else:
                    Ms, Mr = A.helmat(st, B, 2), A.helmat(rf, Bref, 2)
                    vK = {h: np.exp(dl)*A.pred([Mo[j, j] for Mo in Ms], [Mo[j, j] for Mo in Mr], 2)*fk for j, h in enumerate(hels)}
                    sq = 1.0
                    if len(hels) == 2:
                        r2 = Mr if kd != 'C3H' else A.helmat(get(RK[k2kind], DREF, nu), Bref, 2)
                        sq = np.sqrt(abs((r2[0][0, 0]/r2[0][1, 1]).real))
                        off = 0.5*np.array([Mo[0, 1] + Mo[1, 0] for Mo in Ms])
                        vK['off'] = np.exp(dl)*A.pred(off, [Mo[0, 0] for Mo in Mr], 2)*fk*sq
                v, s_ = vK, 'K2'
                verr = {e_: np.nan for e_ in vK}
                if kd in DB:
                    ds, dr = get(DB[kd], D, nu), get(DB[kd], DREF, nu)
                    NUCUR[0] = nu
                    if ds is not None and dr is not None and (J >= 4 or odd):
                        vd = ratio_odd(ds, dr, D, J, c0, fk) if odd else ratio_even(ds, dr, B, Bref, D, J, c0, fk, sq)
                        eG = max(m*m*vd[h][1]/max(abs(vd[h][0]), 1e-300)*abs(Dsh(h, J, x)) for h in hels)
                        if kd == 'EE':          # what enters the nonlocal functional is a(nu) x (forward-normalised value)
                            hm_, hp_ = hhat('EE', np.array([nu])); eG = eG*float(abs(hm_ + hp_)[0])
                        LOG.append(f"direct {kd} J={J} D={D:.4f} x={x:.3f} nu={nu:.3f} [{ds[1]}]: Gamma-unit err ~ {eG:.0f} c2 -> "
                                   f"{'used' if eG <= GERR else 'REJECTED (K2 kept)'}")
                        if eG <= GERR:
                            v = {e_: vd[e_][0] for e_ in vd}; verr = {e_: vd[e_][1] for e_ in vd}; s_ = ds[1]
                src.append(s_)
                for h in hels:
                    sc = fLO(h, J, k, xr)
                    G[h].append(sc*v[h]); Ge[h].append(abs(sc)*verr[h])
                if len(hels) == 2:
                    sh_ = np.sqrt(A.env('++', xr, k)*A.env('+-', xr, k)*A.Dsh('++', J, xr)*A.Dsh('+-', J, xr))
                    Go.append(v['off']*sh_); Goe.append(verr['off']*abs(sh_))
            rows[kd] = dict(k=k, G={h: np.array(G[h]) for h in hels}, Ge={h: np.array(Ge[h]) for h in hels}, Go=np.array(Go),
                            Goe=np.array(Goe), src=src)
        # a_h from the k = 2 kind, nodes x <= 0.25 with nu >= 15 (local de-shift regular there), iterated with the local de-shift
        # q^2: l9's rule (x <= 0.25 & nu >= 15, else the 3 smallest x); EE: same, but the fallback avoids nu < Y2_AFITNU (default 9) where
        # the local u_EE de-shift used inside the fit is near its pole (6.64)
        fin = np.isfinite(rows[k2kind]['G'][hels[0]])
        sel = fin & (xs <= XFIT + 1e-9) & (xs*m >= NUFIT)
        if sel.sum() < 3:
            numin = 0.0 if k2kind == 'C2q' else float(os.environ.get('Y2_AFITNU', '9'))
            cand = np.where(fin & (xs*m >= numin))[0]
            sel = np.zeros(len(xs), bool); sel[cand[:3]] = True
        ush = np.array([u_of('q2' if k2kind == 'C2q' else 'EE', x*m) for x in xs])
        def shloc(h, k):
            if k == 3:
                return np.zeros(len(xs))
            return np.array([((1 + uu)/2*fLO(h, J, 2, x + 1j/m) + (1 - uu)/2*fLO(h, J, 2, x - 1j/m) - fLO(h, J, 2, x)).real
                             for x, uu in zip(xs, ush)])
        def fit(h, g):
            Am = np.array([Dsh(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(Am, (g/A.env(h, xs, 2))[sel], rcond=None)[0][0]
        a = {}
        for h in hels:
            a[h] = fit(h, rows[k2kind]['G'][h].real)
            for _ in range(3):
                a[h] = fit(h, rows[k2kind]['G'][h].real - a[h]*shloc(h, 2))
        rec = dict(m=m, xs=xs, a=a, fitx=xs[sel], c0=c0)
        hm_q, hp_q = hhat('q2', xs*m); hm_e, hp_e = hhat('EE', xs*m)
        for kd, rw in rows.items():
            k = rw['k']
            st_ = 'q2' if kd == 'C2q' else 'EE'
            hm, hp = (hm_q, hp_q) if st_ == 'q2' else (hm_e, hp_e)
            an = hm + hp
            isK2 = np.array([s == 'K2' for s in rw['src']])
            for h in hels:
                sg = SIG[h] if k == 3 else 1.0
                env_k = A.env(h, xs, k)
                G = rw['G'][h]
                rec[(kd, h)] = m*m*((G.real - a[h]*shloc(h, k))/(a[h]*env_k) - sg*Dsh(h, J, xs))
                rec[(kd, h, 'err')] = m*m*rw['Ge'][h]/abs(a[h]*env_k)
                if k == 2:
                    fp = np.array([fLO(h, J, 2, x + 1j/m) for x in xs]); fm = np.array([fLO(h, J, 2, x - 1j/m) for x in xs])
                    Gn = G/a[h]
                    Ahf_chi = hm*fp + hp*fm                                  # a(nu) A_h[f]  (finite)
                    ulo = (hm - hp)/an
                    Ahf_im = (0.5*(1 + ulo)*fp + 0.5*(1 - ulo)*fm).imag      # exact LO Im of A_h[f] (singular where an -> 0; only used x an)
                    Gn = np.where(isK2 & (xs*m < NUIM), Gn.real + 1j*Ahf_im, Gn)
                    rec[(kd, h, 'nl')] = m*m*(an*Gn - Ahf_chi)/env_k
                    rec[(kd, h, 'nlerr')] = m*m*np.abs(an)*rw['Ge'][h]/abs(a[h]*env_k)
                rec[(kd, h, 'src')] = rw['src']
            if len(hels) == 2:
                eo = np.sqrt(A.env('++', xs, k)*A.env('+-', xs, k))
                Go = rw['Go']/np.sqrt(a['++']*a['+-'])/eo
                rec[(kd, 'off')] = m*Go.real
                rec[(kd, 'off', 'im')] = m*Go.imag
                rec[(kd, 'off', 'err')] = m*rw['Goe']/np.sqrt(abs(a['++']*a['+-']))/eo
                if k == 2:
                    Gq = np.where(isK2 & (xs*m < NUIM), Go.real + 0j, Go)
                    rec[(kd, 'off', 'nl')] = m*an*Gq
        rec['src'] = rows[k2kind]['src']
        out[D] = rec
    return out


def align(g, J, k2kind):
    """sign per row from the q^2 C2 imaginary O(1/m) mixing m Im Gamma_off (S1_S4 §2.2(a)); applied to all kinds."""
    Ds = sorted(g)
    if J < 4 or J % 2 or len(Ds) < 2:
        return g
    ref = g[Ds[-1]]; xr = ref['xs']; ir = ref[('C2q', 'off', 'im')]
    for D in Ds:
        r = g[D]
        xi = np.intersect1d(np.round(r['xs'], 6), np.round(xr, 6))
        if len(xi) < 3:
            LOG.append(f"align J={J} D={D}: <3 common nodes, not aligned"); continue
        aa = np.array([r[('C2q', 'off', 'im')][np.argmin(np.abs(r['xs'] - x))] for x in xi])
        bb = np.array([ir[np.argmin(np.abs(xr - x))] for x in xi])
        ok = np.isfinite(aa) & np.isfinite(bb)
        ov = float(aa[ok] @ bb[ok]/np.sqrt((aa[ok] @ aa[ok])*(bb[ok] @ bb[ok])))
        LOG.append(f"align J={J} D={D:.0f}: overlap {ov:+.2f}")
        if ov < 0:
            for key in list(r):
                if isinstance(key, tuple) and len(key) >= 2 and key[1] == 'off' and (len(key) == 2 or key[2] in ('im', 'nl')):
                    r[key] = -r[key]
    return g


# ------------------------------------------------------------------------------------------------ model of the grid (x spline, 1/m linear)
def model(g, keys):
    Ds = sorted(g); ms = np.array([g[D]['m'] for D in Ds])
    spl = {}
    for D in Ds:
        r = g[D]
        for key in keys:
            if key not in r:
                continue
            v = np.asarray(r[key])
            ok = np.isfinite(v) if not np.iscomplexobj(v) else np.isfinite(v.real) & np.isfinite(v.imag)
            if ok.sum() < 4:
                continue
            xs = r['xs'][ok]
            spl[(D, key)] = (CubicSpline(xs, v[ok], bc_type='natural'), xs[0], xs[-1])
    def f(key, x, m):
        x = np.asarray(x, float)
        sc_below = 1.0
        if m <= ms[0]:
            i0 = i1 = 0; t = 0.0
            if BELOW == 'scale' and 'off' not in key:      # [M] variant: Gamma ~ m J^2 below the lowest row (relative correction ~ J^2/m)
                sc_below = m/ms[0]
        elif m >= ms[-1]:
            i0 = i1 = len(Ds) - 1; t = 0.0
        else:
            i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
            t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
        def ev(D):
            s, x0, x1 = spl[(D, key)]
            xx = np.minimum(x, x1)
            if EXTRAP == 'hold':
                xx = np.maximum(xx, x0)
            val = s(xx)
            if EXTRAP == 'zero':
                val = np.where(x < x0, 0*val, val)
            return val
        return sc_below*((1 - t)*ev(Ds[i0]) + t*ev(Ds[i1]))
    return f


# ------------------------------------------------------------------------------------------------ near gap
def states_base():
    NGp = NearGap(DG); F = NGp.F; pn, pw, Lm = NGp.pn, NGp.pw, NGp.Lam
    blocks = {}
    for h, DF, r0, ref in load_states(DG):
        mfl, J = r0['mfl'], r0['J']; u = 1/mfl**2; s = SIGMA[h]
        x = pn*np.sqrt(u); y = x*x
        P = remk.Pders_vec(Lm*pn, r0['eta'], nmax=0)[0]; jz = j0((2*J + 1)*x)
        pf2, wf2 = gl01(max(400, 2*J + 100))
        Js, B2, B3, _ = flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
        p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
        cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
        fac = (1 - y)**2 if h == '++' else (1 - y)**-2
        rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum() + (p3*u*s*envk(3, x)*fac*r0['rho'][3]).sum()
        Ns = abs(B2 + cv2) + abs(B3 + cv3)
        blocks.setdefault((round(mfl, 9), J, DF), {})[h] = dict(base=B2 + B3 + cv2 + cv3 + rem, Ns=Ns, r0=r0, x=x, y=y, p2=p2, p3=p3, u=u,
                                                                 m=r0['m'])
    return NGp, blocks


class Weight:
    def __init__(self, struct, gauge, Lam, nq=1200):
        nu, Y, X, phi = weight(struct, gauge)
        self.nu, self.Y, self.X = nu, Y, X
        self.Lam = Lam
        # fine quadrature in p on [0, 1] (graded towards p = 0 where X has structure)
        t, w = np.polynomial.legendre.leggauss(nq)
        s = 0.5*(t + 1); ws = 0.5*w
        self.p = s**2; self.wp = 2*s*ws                     # p = s^2
        nuq = Lam*self.p
        self.Yq = np.interp(nuq, nu, Y.real) + 1j*np.interp(nuq, nu, Y.imag)
        self.nuq = nuq
        self.phi_ok = phi
        self.psi2q = None


def neargap(grids, blocks, NGp, W, k2kind, mode):
    """records per state; mode 'loc' or 'nl' for the k = 2 term."""
    Lm = NGp.Lam
    mods = {J: model(g, sorted(set(kk for D in g for kk in g[D] if isinstance(kk, tuple) and kk[-1] != 'src'), key=str)) for J, g in grids.items()}
    recs = []
    for (mfl, J, DF), b in blocks.items():
        hs = [h for h in ('++', '+-') if h in b]
        base = np.array([b[h]['base'] for h in hs]); Ns = np.array([b[h]['Ns'] for h in hs])
        dF = np.zeros((len(hs), len(hs)))
        e = b[hs[0]]; x, y, p2, p3, u, m = e['x'], e['y'], e['p2'], e['p3'], e['u'], e['m']
        rho2, rho3 = e['r0']['rho'][2], e['r0']['rho'][3]
        r20 = np.polyfit(x[:3], rho2[:3], 2)[-1] if np.any(rho2) else 0.0
        cov = J in mods
        if cov:
            f = mods[J]
            xq = W.nuq/m
            wl = np.where(x < XC, TLOW[0], 1.0); wq = np.where(xq < XC, TLOW[0], 1.0)      # node-region scaling (tolerance scans)
            for kk, h in enumerate(hs):
                Dh = Dsh(h, J, x)
                envd = A.env(h, x, 2) if J % 2 == 0 else P4.env(h, x, 2)
                envt = A.env(h, x, 3) if J % 2 == 0 else P4.env(h, x, 3)
                if mode == 'loc':
                    t2 = (wl*p2*envd*(f((k2kind, h), x, m)/m**2)).sum()
                elif mode == 'locq':
                    envq = A.env(h, xq, 2)
                    t2 = (wq*W.wp*W.psi2q*envq*f((k2kind, h), xq, m)/m**2).sum()
                else:
                    envq = A.env(h, xq, 2)
                    t2 = (wq*W.wp*W.Yq*envq*f((k2kind, h, 'nl'), xq, m)/m**2).sum().real
                dF[kk, kk] = t2 - (wl*p2*envd*(rho2 - Dh*r20)).sum()
                dF[kk, kk] += (wl*p3*u*envt*(f(('C3H', h), x, m)/m**2 - SIG[h]*(rho3 - Dh*r20))).sum()
            if len(hs) == 2:
                eo2 = np.sqrt(A.env('++', x, 2)*A.env('+-', x, 2))
                c3 = f(('C3H', 'off'), x, m)
                if mode == 'loc':
                    o2 = (wl*p2*eo2*f((k2kind, 'off'), x, m)).sum()/m
                elif mode == 'locq':
                    eoq = np.sqrt(A.env('++', xq, 2)*A.env('+-', xq, 2))
                    o2 = (wq*W.wp*W.psi2q*eoq*f((k2kind, 'off'), xq, m)).sum()/m
                else:
                    eoq = np.sqrt(A.env('++', xq, 2)*A.env('+-', xq, 2))
                    o2 = (wq*W.wp*W.Yq*eoq*f((k2kind, 'off', 'nl'), xq, m)).sum().real/m
                o = o2 + (wl*p3*u*eo2/(2 - y)*c3).sum()/m
                dF[0, 1] = dF[1, 0] = o/m      # c_off = m Gamma_off model: grids hold m Gamma_off (see run())
        recs.append(dict(mfl=mfl, J=J, DF=DF, hs=hs, base=base, Ns=Ns, dF=dF, cov=cov))
    return recs


def lam(r, t):
    Fm = np.diag(r['base']) + t*r['dF']; n = 1/np.sqrt(r['Ns'])
    return np.linalg.eigvalsh(n[:, None]*Fm*n[None, :])[0]


def tstar(recs, cond):
    out = []
    for sgn in (1, -1):
        lo, hi = 0.0, 1.0
        while cond(sgn*hi) and hi < 1e4:
            lo, hi = hi, 2*hi
        if hi >= 1e4:
            out.append(sgn*np.inf); continue
        for _ in range(40):
            mid = 0.5*(lo + hi)
            lo, hi = (mid, hi) if cond(sgn*mid) else (lo, mid)
        out.append(sgn*lo)
    return out[1], out[0]


def summarize(recs, pr, label):
    grid = sorted(set(r['mfl'] for r in recs))
    def nval(t, rr=recs):
        return next((n0 for n0 in grid if all(lam(r, t) > 0 for r in rr if r['mfl'] >= n0 - 1e-9)), None)
    w1 = min(recs, key=lambda r: lam(r, 1.0)); w0 = min(recs, key=lambda r: lam(r, 0.0))
    ts = tstar(recs, lambda t: nval(t) == grid[0])
    res = dict(n=nval(1.0), min=lam(w1, 1.0), w=(w1['mfl'], w1['J'], w1['DF'], tuple(w1['hs'])), tstar=ts, perJ={})
    for sect, cond in (('even', lambda r: r['J'] % 2 == 0), ('odd', lambda r: r['J'] % 2 == 1)):
        rr = [r for r in recs if cond(r)]
        ws = min(rr, key=lambda r: lam(r, 1.0))
        res[sect] = dict(n=nval(1.0, rr), min=lam(ws, 1.0), w=(ws['mfl'], ws['J'], tuple(ws['hs'])),
                         tstar=tstar(rr, lambda t: nval(t, rr) == grid[0]))
    pr(f"  [{label}] n = {res['n']} (uncorrected {nval(0.0)}); min lam/N_s = {res['min']:+.4f} [m/Lam={w1['mfl']} J={w1['J']} Df={w1['DF']} "
       f"{w1['hs']}] (uncorr. {lam(w0, 0.0):+.4f} [J={w0['J']}]); t* = ({ts[0]:+.2f}, {ts[1]:+.2f})")
    for sect in ('even', 'odd'):
        s = res[sect]
        pr(f"      {sect}: n = {s['n']}, min {s['min']:+.4f} at (m/Lam, J, h) = {s['w']}, t* = ({s['tstar'][0]:+.2f}, {s['tstar'][1]:+.2f})")
    for J in sorted(set(r['J'] for r in recs)):
        rj = [r for r in recs if r['J'] == J]
        wj = min(rj, key=lambda r: lam(r, 1.0))
        dg = max(abs(r['dF'][k, k])/r['Ns'][k] for r in rj for k in range(len(r['hs'])))
        of = max((abs(r['dF'][0, 1])/np.sqrt(r['Ns'].prod()) if len(r['hs']) == 2 else 0) for r in rj)
        tj = tstar(rj, lambda t: all(lam(r, t) > 0 for r in rj))
        res['perJ'][J] = dict(min=lam(wj, 1.0), unc=min(lam(r, 0.0) for r in rj), at=wj['mfl'], diag=dg, off=of, tJ=tj, cov=wj['cov'])
        pr(f"      J={J:2d}{'' if wj['cov'] else ' (uncovered)'}: min {lam(wj, 1.0):+.4f} (uncorr. {min(lam(r, 0.0) for r in rj):+.4f}) at "
           f"m/Lam={wj['mfl']}; max|dF_hh|/N_h {dg:.4f}; max|dF_+-|/sqrt(NN) {of:.4f}; t*_J = ({tj[0]:+.2f}, {tj[1]:+.2f})")
    return res


def build_all(k2kind, pr):
    jev = [0, 2, 4, 6, 8, 10] + ([12, 14, 16, 18, 20, 24, 28] if JHI else [])
    jod = [5, 7, 9] + ([11, 13] if JHI else [])
    grids = {}
    for J in jev + jod:
        try:
            g = build_grid(J, (k2kind, 'C3H'), 'EE' if k2kind == 'EE' else 'q2')
        except Exception as ex:
            pr(f"  J={J}: no grid ({type(ex).__name__}: {ex})"); continue
        if not g:
            pr(f"  J={J}: no rows"); continue
        g = align(g, J, k2kind)
        if J % 2 == 0 and J >= 4:      # c_off = m Gamma_off (S1_S4 §2.2(b)); neargap divides by the state's m
            for D, r in g.items():
                for key in list(r):
                    if isinstance(key, tuple) and len(key) >= 2 and key[1] == 'off' and (len(key) == 2 or key[2] in ('im', 'nl')):
                        r[key] = r[key]*r['m']
        if len(g) == 1:
            D = list(g)[0]; g[D + 1e-6] = dict(g[D], m=g[D]['m']*(1 + 1e-9))
        grids[J] = g
    return grids


def node_table(grids, pr, k2kind):
    pr(f"# node table ({k2kind}): Delta = 100 row (and extra rows), Gamma (c2 units x LO shape): local / nonlocal(Re, Im) / err, source")
    for J, g in sorted(grids.items()):
        for D, r in sorted(g.items()):
            if D > 150.0:
                continue
            for h in (('+-',) if J % 2 else (('++', '+-') if J >= 4 else ('++',))):
                pr(f"  J={J:2d} D={D:8.3f} m={r['m']:.2f} a_{h}={r['a'][h]:.4g} fit x={list(np.round(r['fitx'], 3))}")
                pr("     x   : " + " ".join(f"{x:8.3f}" for x in r['xs']))
                pr("     src : " + " ".join(f"{s:>8s}" for s in r[(k2kind, h, 'src')]))
                pr("     loc : " + " ".join(f"{v:8.0f}" for v in r[(k2kind, h)]))
                pr("     NLre: " + " ".join(f"{v.real:8.0f}" for v in r[(k2kind, h, 'nl')]))
                pr("     NLim: " + " ".join(f"{v.imag:8.0f}" for v in r[(k2kind, h, 'nl')]))
                pr("     err : " + " ".join(f"{v:8.0f}" if np.isfinite(v) else "       -" for v in r[(k2kind, h, 'err')]))
                pr("     C3H : " + " ".join(f"{v:8.0f}" for v in r[('C3H', h)]) + "   (src " + ",".join(r[('C3H', h, 'src')][:4]) + ")")
                pr("     C3err:" + " ".join(f"{v:8.0f}" if np.isfinite(v) else "       -" for v in r[('C3H', h, 'err')]))


def mc(grids, blocks, NGp, W, k2kind, mode, nmc, seed=1):
    """Gaussian redraws of every direct node (k = 2 kind and C3H; errors from replicas) -> spread of the per-sector minima."""
    rng = np.random.default_rng(seed)
    mins = []
    for it in range(nmc):
        gg = {}
        for J, g in grids.items():
            gg[J] = {}
            for D, r in g.items():
                r2 = dict(r)
                for key in list(r):
                    if not (isinstance(key, tuple) and len(key) == 2 and key[1] in ('++', '+-')):
                        continue
                    kd, h = key
                    e = np.where(np.isfinite(r[(kd, h, 'err')]), r[(kd, h, 'err')], 0.0)
                    z = rng.standard_normal(len(e))
                    r2[key] = r[key] + z*e
                    if (kd, h, 'nl') in r:
                        en = np.where(np.isfinite(r[(kd, h, 'nlerr')]), r[(kd, h, 'nlerr')], 0.0)
                        r2[(kd, h, 'nl')] = r[(kd, h, 'nl')] + en*(z + 1j*rng.standard_normal(len(e)))/np.sqrt(1.0)
                gg[J][D] = r2
        recs = neargap(gg, blocks, NGp, W, k2kind, mode)
        ev = [r for r in recs if r['J'] % 2 == 0]; od = [r for r in recs if r['J'] % 2 == 1]
        mins.append((min(lam(r, 1.0) for r in ev), min(lam(r, 1.0) for r in od)))
    return np.array(mins)


if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else 'main'
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    NGp, blocks = states_base()
    pr(f"# y2_chain100 [{tag}]: Delta_g = {DG}, Lam = {NGp.Lam:.4f}; gauge {GAUGE}; extra rows {XROWS}; extrapolation below x_min: {EXTRAP}; "
       f"J_hi {JHI}; direct {DIRECT}; Im rule: K2 Im replaced by exact LO at nu < {NUIM}")
    store = {}
    for k2kind in os.environ.get('Y2_KINDS', 'EE,C2q').split(','):
        grids = build_all(k2kind, pr)
        W = Weight('EE' if k2kind == 'EE' else 'q2', GAUGE if k2kind == 'EE' else 'min', NGp.Lam)
        W.psi2q = NGp.F.psi2(W.p)
        # weight check: psi2 of the chain vs phi of u3 (same function of p)
        pc = np.linspace(0.02, 0.98, 25)
        dphi = np.max(np.abs(NGp.F.psi2(pc) - np.interp(NGp.Lam*pc, W.nu, W.phi_ok)))/np.max(np.abs(NGp.F.psi2(pc))) if W.phi_ok is not None else np.nan
        pr(f"## k = 2 kind {k2kind}: rows per J: " + "; ".join(f"{J}: " + ",".join(f"{D:.4g}" for D in sorted(g) if D < 1e4) for J, g in sorted(grids.items()))
           + f"   [psi2 (chain) vs u3 phi: max rel dev {dphi:.1e}]")
        node_table(grids, pr, k2kind)
        res = {}
        for mode in os.environ.get('Y2_MODES', 'loc,locq,nl').split(','):
            recs = neargap(grids, blocks, NGp, W, k2kind, mode)
            res[mode] = summarize(recs, pr, f"{k2kind} {mode} {GAUGE if (mode == 'nl' and k2kind == 'EE') else ''}")
            res[mode]['recs'] = [dict(mfl=r['mfl'], J=r['J'], DF=r['DF'], hs=r['hs'], base=r['base'], Ns=r['Ns'], dF=r['dF'], cov=r['cov'])
                                 for r in recs]
        if os.environ.get('Y2_TOL', '0') == '1':
            def okt(t, md='nl'):
                TLOW[0] = t; rr = neargap(grids, blocks, NGp, W, k2kind, md); TLOW[0] = 1.0
                return min(lam(r, 1.0) for r in rr) > 0
            edges = []
            for sgn in (1, -1):
                lo, hi = 1.0, 1.0 + sgn
                while okt(hi) and abs(hi) < 64:
                    lo, hi = hi, hi + sgn*abs(hi - 1.0)
                for _ in range(10):
                    mid = 0.5*(lo + hi); lo, hi = (mid, hi) if okt(mid) else (lo, mid)
                edges.append(lo)
            pr(f"  [{k2kind} nl] tolerance window of the region x < {XC} (its whole correction x t): all states positive for t in "
               f"[{edges[1]:+.2f}, {edges[0]:+.2f}]")
        nmc = int(os.environ.get('Y2_MC', '0'))
        if nmc:
            for mode in [m_ for m_ in ('nl', 'loc') if m_ in res]:
                mm = mc(grids, blocks, NGp, W, k2kind, mode, nmc)
                pr(f"  [{k2kind} {mode}] MC over direct-node errors ({nmc} draws): even min {mm[:, 0].mean():+.4f} +- {mm[:, 0].std():.4f} "
                   f"(range {mm[:, 0].min():+.4f} .. {mm[:, 0].max():+.4f}); odd min {mm[:, 1].mean():+.4f} +- {mm[:, 1].std():.4f} "
                   f"(range {mm[:, 1].min():+.4f} .. {mm[:, 1].max():+.4f})")
                res[mode]['mc'] = mm
        store[k2kind] = dict(res=res, grids=grids)
    for s in LOG:
        pr("# " + s)
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'y2_chain100_%s.txt' % tag), 'w').write(buf.getvalue())
    pickle.dump(store, open(os.path.join(C, 'y2_recs_%s.pkl' % tag), 'wb'))

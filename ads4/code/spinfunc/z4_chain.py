"""
z4 (EE_DG300_SMALLNU.md): the EE_NEARGAP w2 chain (imported read-only; all its conventions: q3/S1_S4 K2 linear chain, EE rows, EE's own
helicity basis, aligned off-diagonals, c_off = m Gamma_off, C^_3 without de-shift, phase strip, Hampel despike, combined near-gap +
hybrid states) with
  (1) the z2 DIRECT crystal points (n3 contour above the fold, n4 thimbles below; state Delta = 336 and reference 2e4 on common random
      numbers; z3 records) injected as extra x-nodes of the Delta = 336 rows of J in Z4_JS (default 24, 28) -- and, with Z4_OVERRIDE = 1
      (default), replacing the K2 nodes at the same keys (nu = 16.68, 25.0, 33.4);   estimator Z4_EST = dir | cv (scalar control
      variate), Z4_JK = i (replica i left out; error propagation);
  (2) a NON-LOCAL k = 2 term at nu < Z4_NUNL (default 25) on every row that has direct nodes: the chain's local de-shifted density
      psi2 env Gamma_loc is replaced by  Re[ (Y/phi)(nu) env Gamma_NL ]  with  Gamma_NL := m^2 a_chi(nu) (G/a_h - A_h[f_LO])/env
      (complex; a_chi = h^_- + h^_+ in u3's chi-normalisation, so Re[(Y/phi) Gamma_NL] = Re[(X/phi) Gamma_c], Gamma_c the complex
      forward-normalised correction; finite through the real zero nu = 6.643 of h_- + h_+);  X = a Y the Abel-form weight of u3
      (Z4_GAUGE = min, min-norm) or x2 (Z4_GAUGE = can, canonical);  phi = psi2.  Gamma_NL is interpolated in NU (not x) on the row and
      held/zeroed/linearly extrapolated below the lowest node (Z4_EXTRAP = hold | zero | lin);  for a bracketing row WITHOUT direct
      nodes the density is Re(X/phi) Gamma_loc (= w2's W2_GAUGE = min convention there).  Rows are interpolated in 1/m as in the chain.
      Z4_NL = 0: local de-shift everywhere (w2 convention; singular next to nu = 6.64 -- comparison only).
  The scalar remainder (rho_2, which cancels between base and dF) keeps the psi2 weight; W2_GAUGE is forced to none.
Usage: Z4_TAG=x [Z4_...] Q3_JEVEN=24,28 Q3_JODD= W2_ROWS=336,500,1600 python3 z4_chain.py 300  -> cache/z4_chain_<tag>.txt/.pkl
"""
import os, sys, io, copy, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
os.environ['W2_GAUGE'] = 'none'
os.environ.setdefault('W2_BASIS', 'ee')
os.environ.setdefault('W2_MODES', 'ee')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.interpolate import CubicSpline
import w2_chain as W2
import z3_nodes as Z3

Q = W2.Q
A = Q.A
C = os.path.join(HERE, 'cache')
JS = [int(v) for v in os.environ.get('Z4_JS', '24,28').split(',') if v]
DIRECT = os.environ.get('Z4_DIRECT', '1') == '1'
EST = os.environ.get('Z4_EST', 'dir')
JK = None if os.environ.get('Z4_JK', '') == '' else int(os.environ['Z4_JK'])
OVERRIDE = os.environ.get('Z4_OVERRIDE', '1') == '1'
NL = os.environ.get('Z4_NL', '1') == '1'
NUNL = float(os.environ.get('Z4_NUNL', '25'))
NUFLAT = float(os.environ.get('Z4_NUFLAT', '14'))     # C2 local values at direct nodes below this are flattened (local de-shift junk)
GAUGE = os.environ.get('Z4_GAUGE', 'min')
EXTRAP = os.environ.get('Z4_EXTRAP', 'hold')
NUMIN = float(os.environ.get('Z4_NUMIN', '0'))        # use only direct nodes with nu >= NUMIN (test of the low-nu nodes)
NUMAX = float(os.environ.get('Z4_NUMAX', '1e9'))
TAG = os.environ.get('Z4_TAG', '')
DG = 300
LOG = []

# ------------------------------------------------------------------------------------------------ record injection
_ee_orig = W2.ee_even
_q2e_orig = Q._q2_even


def _inject(R, kind, J):
    if not DIRECT or J not in JS:
        return R
    add = Z3.records(kind, J, EST, JK, numin=NUMIN, numax=NUMAX)
    R = dict(R)
    nus = sorted(set(k[1] for k in add))
    for nu in nus:
        ks = [k for k in add if k[1] == nu]
        if not OVERRIDE and any(k in R for k in ks):
            continue
        for k in ks:
            R[k] = add[k]
    return R


def ee_even(J):
    return _inject(_ee_orig(J), 'EE', J)


def q2_even(kind, J):
    R = _q2e_orig(kind, J)
    return _inject(R, kind, J) if kind == 'C3H' else R


W2.ee_even = ee_even
Q._q2_even = q2_even


# ------------------------------------------------------------------------------------------------ weights, magic-block coefficients
_F = W2._F


def hhat(nu):
    nu = np.asarray(nu, complex)
    Rm, Rp = _F['EE'][0](nu), _F['EE'][1](nu)
    qm, qp = _F['q2'][0](nu), _F['q2'][1](nu)
    nq = qm + qp
    n = (np.abs(Rm) + np.abs(Rp))*nq/np.abs(nq)
    return Rm/n, Rp/n


def _weight():
    b = pickle.load(open(os.path.join(C, 'u3_chi_EE_%d.pkl' % DG), 'rb'))
    nu = b['nu']; phi = b['phi'].real
    hm, hp = hhat(nu); a = hm + hp
    if GAUGE == 'can':
        c = pickle.load(open(os.path.join(C, 'x2_chi_EE_can_%d.pkl' % DG), 'rb'))
        assert np.allclose(c['nu'], nu)
        X = c['X']; Y = c['Y'] if 'Y' in c else X/np.where(np.abs(a) > 1e-12, a, 1e-12)
    else:
        X = b['X']; Y = b['Y']
    dev = np.max(np.abs(X - a*Y))/np.max(np.abs(X))
    LOG.append(f"weight EE/{GAUGE} (Delta_g = {DG}): max|X - a Y|/max|X| = {dev:.1e}")
    return nu, X/np.where(np.abs(phi) > 0, phi, 1.0), Y/np.where(np.abs(phi) > 0, phi, 1.0)


WNU, XPHI, YPHI = _weight()


def xphi(nu):
    return np.interp(nu, WNU, XPHI.real) + 1j*np.interp(nu, WNU, XPHI.imag)


def yphi(nu):
    return np.interp(nu, WNU, YPHI.real) + 1j*np.interp(nu, WNU, YPHI.imag)


# ------------------------------------------------------------------------------------------------ grid wrapper: NL rows
NLROWS = {}
_A_grid_orig = W2._A_grid


def _ahat_im(J, h, xs, m):
    """Im of A_h[f_LO] (forward-normalised local shift, u = u_EE), f = env_{2,h} D_h."""
    f = lambda t: A.env(h, t, 2)*A.Dc(h, J, t)
    out = []
    for x in xs:
        u = complex(Q.u_spin(x*m)); gp, gm = 0.5*(1 + u), 0.5*(1 - u)
        out.append((gp*f(x + 1j/m) + gm*f(x - 1j/m)).imag)
    return np.array(out)


def _A_grid(J, K, kinds=('C2', 'C3H'), *a, **k):
    g = _A_grid_orig(J, K, kinds, *a, **k)
    if J not in JS or not DIRECT:
        return g
    for D, r in g.items():
        if 'C2' not in kinds:
            continue
        m = r['m']; xs = r['xs']; nus = xs*m
        Rz = Z3.zpts('EE', J)
        dk = np.array([k[1] for k in Rz if k[0] == float(D)])
        dnodes = np.array([dk.size > 0 and np.min(np.abs(dk - v)) < 1e-6*v for v in nus])
        if not dnodes.any():
            continue
        hels = [h for h in ('++', '+-') if ('C2', h) in r]
        sel = nus < 40.0
        rec = dict(m=m, nus=nus[sel], direct=dnodes[sel])
        for h in hels:
            gre = r[('C2', h)][sel]
            gim = m*r[('C2', h, 'im')][sel] - m*m*_ahat_im(J, h, xs[sel], m)/A.env(h, xs[sel], 2)
            gc = gre + 1j*gim
            hm, hp = hhat(nus[sel])
            rec[h] = gc*(hm + hp)                                   # Gamma_NL (finite through the zero)
            rec[(h, 'c')] = gc
            rec[(h, 'C3H')] = r[('C3H', h)][sel].copy()
        NLROWS[(J, float(D))] = rec
        # flatten the local C2 values at direct nodes nu < NUFLAT (local de-shift singular near nu = 6.64; replaced by the NL term)
        lo = (nus < NUFLAT) & dnodes
        if lo.any():
            i1 = int(np.argmax(~lo & (nus >= NUFLAT)))
            for h in hels:
                r[('C2', h)] = np.where(lo, r[('C2', h)][i1], r[('C2', h)])
            LOG.append(f"J={J} D={D:g}: C2 local values flattened at nu = {np.round(nus[lo], 2).tolist()} (NL term used below {NUNL})")
    return g


W2._A_grid = _A_grid


def _nl_row(J, D, h, nu):
    rec = NLROWS[(J, float(D))]
    nn = rec['nus']; gv = rec[h]
    o = np.argsort(nn); nn = nn[o]; gv = gv[o]
    sr = CubicSpline(nn, gv.real, bc_type='natural'); si = CubicSpline(nn, gv.imag, bc_type='natural')
    v = sr(np.clip(nu, nn[0], None)) + 1j*si(np.clip(nu, nn[0], None))
    if EXTRAP == 'zero':
        v = np.where(nu < nn[0], 0.0, v)
    elif EXTRAP == 'lin':
        s = (gv[1] - gv[0])/(nn[1] - nn[0])
        v = np.where(nu < nn[0], gv[0] + s*(nu - nn[0]), v)
    elif EXTRAP == 'const':                 # stress test: Gamma_c := Z4_CONST (real, c2 units) below the lowest node, i.e. Gamma_NL = a_chi C
        hm, hp = hhat(nu)
        v = np.where(nu < nn[0], float(os.environ.get('Z4_CONST', '0'))*(hm + hp), v)
    elif EXTRAP == 'spline':
        v = np.where(nu < nn[0], sr(nu) + 1j*si(nu), v)
    return v


def K2D(J, h, x, m, nu, Gloc, emod):
    """the k = 2 density (Gamma units) used in dF: local Gloc above NUNL; below, the NL density per bracketing row."""
    if not NL or J not in JS or not any(k[0] == J for k in NLROWS):
        return Gloc
    f, Ds = emod
    ms = np.array([Z3.mphys(D, J) for D in Ds])
    if m <= ms[0]:
        i0 = i1 = 0; t = 0.0
    elif m >= ms[-1]:
        i0 = i1 = len(Ds) - 1; t = 0.0
    else:
        i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
        t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
    lo = nu < NUNL
    dens = np.zeros_like(Gloc)
    for i, w in ((i0, 1 - t), (i1, t)):
        if w == 0.0:
            continue
        D = Ds[i]
        if (J, float(D)) in NLROWS:
            d = (yphi(nu)*_nl_row(J, D, h, nu)).real
        else:
            d = xphi(nu).real*f(x, ms[i] if i == i0 else ms[i])[('C2', h)]
        dens = dens + w*d
    return np.where(lo, dens, Gloc)


Q._K2D = K2D
_old = "                dF[kk, kk] = (p2*A.env(h, x, 2)*(G[('C2', h)]/m**2 - (rho2 - Dh*r20))).sum()\n"
_new = "                dF[kk, kk] = (p2*A.env(h, x, 2)*(_K2D(J, h, x, m, Lm*pn, G[('C2', h)], emods[J])/m**2 - (rho2 - Dh*r20))).sum()\n"
assert _old in W2._src
exec(compile(W2._src.replace(_old, _new), 'z4_neargap', 'exec'), Q.__dict__)


# ------------------------------------------------------------------------------------------------ node table
def node_table(pr):
    pr("# NL rows (Delta = 336 rows with direct nodes): per node nu, Gamma_loc (Re Gamma_c, local de-shift), Gamma_c (complex), "
       "Gamma_NL = a_chi Gamma_c, Re[(Y/phi) Gamma_NL] (= the NL density), Re X/phi;  C3H (local, no de-shift)")
    for (J, D), rec in sorted(NLROWS.items()):
        pr(f"  J={J} Delta={D:g} (m = {rec['m']:.2f}):")
        for i, v in enumerate(rec['nus']):
            s = f"    nu={v:7.3f} x={v/rec['m']:.4f} {'dir' if rec['direct'][i] else 'K2 '}"
            for h in ('++', '+-'):
                if h in rec:
                    gc = rec[(h, 'c')][i]; gn = rec[h][i]
                    s += f" | {h}: c {gc.real:+9.0f}{gc.imag:+8.0f}i NL {gn.real:+8.0f}{gn.imag:+8.0f}i dens {float((yphi(v)*gn).real):+8.0f} C3H {rec[(h, 'C3H')][i]:+7.0f}"
            s += f" | ReX/phi {xphi(v).real:+.3f}"
            pr(s)


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    Dgs = [int(a) for a in sys.argv[1:]] or [300]
    t0 = time.time()
    W2.set_mode(True)
    pr(f"######## z4 chain: EE rows + direct z2 nodes (J in {JS}; direct {DIRECT}, est {EST}, jk {JK}, override {OVERRIDE}, nu in "
       f"[{NUMIN}, {NUMAX}]); NL k=2 below nu = {NUNL}: {NL} (gauge {GAUGE}, extrap {EXTRAP}); JEVEN {Q.JEVEN} JODD {Q.JODD}; "
       f"rows {os.environ.get('W2_ROWS', 'all')}; extra m {os.environ.get('Q3_EXTRA_M', '')}")
    grids, out = Q.run(Dgs, pr)
    node_table(pr)
    for s in LOG:
        pr("# " + s)
    pr(f"# despiked nodes (J, Delta, x): {W2.DESPIKE_LOG}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'z4_chain_%s.txt' % TAG), 'w').write(buf.getvalue())
    pickle.dump(dict(out={Dg: {k: v for k, v in o.items() if k != 'recs'} for Dg, o in out.items()},
                     recs={Dg: W2.perstate(out[Dg]['recs']) for Dg in Dgs}, nlrows=NLROWS, grids=grids),
                open(os.path.join(C, 'z4_chain_%s.pkl' % TAG), 'wb'))

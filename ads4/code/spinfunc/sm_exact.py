"""
sm_exact (SMOOTH_BOUND.md §2): EXACT scalar-skeleton remainder beyond the hybrid (flat + leading-Regge curvature) at large J near the
twist gap, for the route-A' weights psi2, psi3 (k = 2, 3 only), at Delta_g = 300 / 1000.

Why not exactmom/skelexact (Mellin, sr_scalcmp.py): BlockMellin's Chebyshev Mack form converts the mp coefficients c_i R^i to doubles;
at (1040, 40) they exceed the double range (see sm_nancheck.py) -> inf/NaN, and the mp set-up is O(N_n J^2) at 60 digits (> 5 min).
Here the Arb Laplace/Barnes route of exactskel/h13_exact.py (Laplace3H, Arb tail weights, FLINT Mack polynomial; validated in H13'
to J = 150) is used unchanged: rho_k = (ex_k - E_k P_J - 2[P_{-1/2+i nu}(cosh b) - J0((2J+1)x)])/E_k (e8 definition), even k exact,
odd k (k = 3) by quartic interpolation in k of ex_k/E_k [M, as e8/h13].  n-sum error = spread over h13's NCONF configurations.

Per state (Delta = Dg + J + dDelta, i.e. at/near the gap) and helicity h (as q3_analysis.neargap's exact-skeleton term):
   rem_h = sum p2 envk(2,x) fac_h rho_2 + sum p3 u s_h envk(3,x) fac_h rho_3,   fac = (1-x^2)^{+2} (++), (1-x^2)^{-2} (+-)
   N_s  = |B2 + cv2| + |B3 + cv3|  (h13_common.hybrid; = ee1000_common.hyb's N_s)
   C_eq = rem_h / [sum p2 envk(2) fac /m^2]   (uniform-c2 equivalent of the remainder; compare SUBLEADING_REGGE §0.1 tolerance)
Also windows nu < 20, nu < 50, and a dense small-nu probe (fold region) of m^2 rho_2 (c2 units).
Usage: python3 sm_exact.py Dg J1,J2,..  [dDelta]   -> cache/sm_exact_<Dg>.pkl (incremental), stdout
Env: SM_XPROBE "a:b:step" extra x nodes (stored after the nu probes; key gets the tag SM_KEYTAG), SM_NPF (p-node factor, default 1.0: n = max(96, SM_NPF*bhat + 64)), SM_DFS (default "3").
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, '..', 'exactskel'); sys.path.insert(0, EX)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import h13_exact as H
from h13_common import F, SIGMA, Lam_of, flat_rows, hybrid
from common import gl01
from nf_common import envk

CACHE = os.path.join(HERE, 'cache')
if os.environ.get('SM_NCONF'):          # e.g. "80:960:1,120:960:1,160:960:1" (N0:nq:Lfactor; the last one is the reported value)
    H.NCONF = tuple(tuple(float(v) if i == 2 else int(v) for i, v in enumerate(c.split(':'))) for c in os.environ['SM_NCONF'].split(','))
NUPROBE = np.array([1.5, 3, 4.5, 6, 7.5, 9, 10, 11, 12, 13, 14, 16, 18, 21, 25, 30, 36, 43, 50])


def near_gap_state(Dg, J, dD=0.0):
    Delta = Dg + J + dD
    alpha = Delta - J - 2.0; abar = alpha + 2*J + 1
    m2 = alpha*abar; m = np.sqrt(m2)
    eta = (alpha**2 + abar**2)/(2*alpha*abar)
    return dict(m=m, m2=m2, eta=eta, beta=float(np.arccosh(eta)), Delta=Delta, mfl=m/Lam_of(Dg), J=J)


def run(Dg, J, dD=0.0, Df=3.0):
    pk = os.path.join(CACHE, 'sm_exact_%d.pkl' % Dg)
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    key = (J, float(dD), Df) + ((os.environ['SM_KEYTAG'],) if os.environ.get('SM_KEYTAG') else ())
    if key in db and not os.environ.get('SM_REDO'):
        return db[key]
    st = near_gap_state(Dg, J, dD)
    Lam = Lam_of(Dg); mfl = st['mfl']; bh = (2*J + 1)/mfl
    npn = max(96, int(float(os.environ.get('SM_NPF', '1.0'))*bh) + 64)
    pn, pw = gl01(npn)
    xpr = os.environ.get('SM_XPROBE', '')                      # extra x-probes "a:b:step" (e.g. the direct-row nodes)
    xp = np.arange(*[float(v) for v in xpr.split(':')]) if xpr else np.zeros(0)
    pp = np.concatenate([pn, NUPROBE/Lam, xp*mfl])
    t0 = time.time()
    H.mack.mack_P_fast(st['Delta'], J)
    r = H.rho_state(Dg, Df, st, int(os.environ.get('H13_PREC', '400')), pn=pp)
    rec = dict(Dg=Dg, J=J, dD=dD, Df=Df, st=st, npn=npn, pn=pn, pw=pw, nus=pp*Lam, x=pp/mfl, rho=r['rho'], rho3=r['rho3'],
               rho_alts=r['rho_alts'], ex=r['ex'], rad=r['rad'], prec=r['prec'], time=None, hel={})
    rec['time'] = time.time() - t0
    rec['hel'] = summarize(rec)
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    db[key] = rec
    pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    return rec


def summarize(rec, rho=None):
    Dg, J, st = rec['Dg'], rec['J'], rec['st']
    Lam = Lam_of(Dg); mfl = st['mfl']; m = st['m']; u = 1/mfl**2
    n = rec['npn']; pn, pw = rec['pn'], rec['pw']; x = pn/mfl; y = x*x
    rho = rec['rho'] if rho is None else rho
    r2, r3 = rho[2][:n], rho[3][:n]
    p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
    out = {}
    for h in ('++', '+-'):
        if h == '+-' and J < 4:
            continue
        s = SIGMA[h]
        fac = (1 - y)**2 if h == '++' else (1 - y)**-2
        B2, B3 = flat_rows(h, u, [J])[J]
        hy = hybrid(h, Dg, mfl, J, B2, B3, beta=st['beta'])
        w2 = p2*envk(2, x)*fac; w3 = p3*u*s*envk(3, x)*fac
        rem2 = float((w2*r2).sum()); rem3 = float((w3*r3).sum())
        sens = float(w2.sum())/m**2
        nu = pn*Lam
        o = dict(Ns=hy['Ns'], Fh=hy['Fh'], rem2=rem2, rem3=rem3, sens=sens, Ceq=(rem2 + rem3)/sens,
                 rem_nu20=float((w2*r2*(nu < 20)).sum() + (w3*r3*(nu < 20)).sum()),
                 rem_nu50=float((w2*r2*(nu < 50)).sum() + (w3*r3*(nu < 50)).sum()),
                 sens50=float((w2*(nu < 50)).sum())/m**2)
        alts = []
        for ra in rec['rho_alts']:
            alts.append(float((w2*ra[2][:n]).sum() + (w3*ra[3][:n]).sum()))
        o['nsum_err'] = max([abs(a - rem2 - rem3) for a in alts] + [0.0])
        o['oddk_err'] = abs(float((w3*rec['rho3'][3][:n]).sum()) - rem3)
        out[h] = o
    return out


def show(rec):
    st = rec['st']; Lam = Lam_of(rec['Dg']); n = rec['npn']
    m = st['m']
    s = (f"Dg={rec['Dg']} J={rec['J']:3d} Delta={st['Delta']:.3f} m/L={st['mfl']:.4f} bhat={(2*rec['J']+1)/st['mfl']:6.1f} "
         f"beta={st['beta']:.4f} npn={n} [{rec['time']:.0f}s, rad {rec['rad']:.0e}]\n")
    for h, o in rec['hel'].items():
        s += (f"   {h}: rem/Ns {(o['rem2']+o['rem3'])/o['Ns']:+.2e} (k2 {o['rem2']/o['Ns']:+.2e}, k3 {o['rem3']/o['Ns']:+.2e}; "
              f"n-sum err {o['nsum_err']/o['Ns']:.0e}, odd-k {o['oddk_err']/o['Ns']:.0e})  Ceq {o['Ceq']:+.2f} c2  "
              f"nu<20 {o['rem_nu20']/o['Ns']:+.1e} nu<50 {o['rem_nu50']/o['Ns']:+.1e}  Ns {o['Ns']:.3e} Fh/Ns {o['Fh']/o['Ns']:+.4f}\n")
    npr = len(NUPROBE)
    G = m*m*rec['rho'][2][n:n + npr]
    s += "   probe nu: " + " ".join(f"{v:6.1f}" for v in rec['nus'][n:n + npr]) + "\n"
    s += "   m^2rho2 : " + " ".join(f"{v:6.1f}" for v in G) + "\n"
    xs = rec['x'][:n]; G2 = m*m*rec['rho'][2][:n]
    s += f"   grid m^2 rho2: rms {np.sqrt(np.mean(G2**2)):.1f}, max|.| {np.max(np.abs(G2)):.1f} at x={xs[np.argmax(np.abs(G2))]:.3f}\n"
    return s


if __name__ == '__main__':
    Dg = int(sys.argv[1]); Js = [int(v) for v in sys.argv[2].split(',')]
    dD = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    for J in Js:
        rec = run(Dg, J, dD)
        print(show(rec), flush=True)

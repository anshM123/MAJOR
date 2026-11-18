"""
e1: exact scalar-skeleton moments on the near-gap heavy states of the Omega^U scan (OMEGA_U.md, u4 grid), Delta_g given.

States: m/Lam on the u4 grid in [1, MFLMAX] (default 2), J <= JMAX (default 10), J in {0,2,..} for ++ and 4..JMAX for +-  (union scanned),
tau >= Delta_g.  Nodes: the 64 Gauss nodes p of u4, nu = Lam p, x = nu/m.
Per state and k = 2..7 (odd k: analytic continuation in k):
   ex_k   = m^{2k} C_{k,nu}[Delta,J]  exact (Mellin, skelexact.py),
   rho_k  = (ex_k - hyb_k)/E_k,  hyb_k = E_k P_J(1-2x^2) + 2[P_{-1/2+i nu}(cosh beta) - J_0((2J+1)x)]   (definition of remk.py),
   also the truncated-series rho, its error estimate and control flag (remk.py, same skeleton), for comparison.
Forward:  R4 = -m^2 d_{nu^2}(m^8 C_4)|_0 - (2 m^2 W - 11/2),  R5 = m^10 C_{5,0} - 2 P_{-1/2}(cosh beta)  exact (nu -> 1e-4, h-step as remk).
Usage: python3 e1_states.py Dg Df [MFLMAX] [JMAX]   -> cache/e1_<Dg>_Df<Df>.pkl (incremental)
Optional env: E1_JLIST (comma list of J, overrides JMAX), E1_MMIN (lowest m/Lam), E1_MSET (comma list of m/Lam, overrides the grid),
              E1_PKL (alternative cache file tag).
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
Dg = int(sys.argv[1]); DF = float(sys.argv[2])
MFLMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
JMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 10
os.environ['OMEGAU_DF'] = '%g' % DF
OMU = os.path.join(HERE, '..', 'omegau'); sys.path.insert(0, OMU)
import numpy as np
from scipy.special import j0, eval_legendre
import remk                                       # read-only import (series remainder, P derivatives)
from common import state, Transforms, KS, gl01   # read-only import
from u4_scan import mfl_grid
from skelexact import ExactState

NP = 64
pn, pw = gl01(NP)


def Ek(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2)


def run():
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    tag = os.environ.get('E1_PKL', '')
    pk = os.path.join(HERE, 'cache', 'e1_%d_Df%g%s.pkl' % (Dg, DF, tag))
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    W0 = Transforms({k: np.array([1.0]) for k in KS})
    Js = sorted(set(range(0, JMAX + 1, 2)) | set(range(4, JMAX + 1)))
    if os.environ.get('E1_JLIST'):
        Js = [int(a) for a in os.environ['E1_JLIST'].split(',')]
    MMIN = float(os.environ.get('E1_MMIN', '0'))
    grid = mfl_grid()
    if os.environ.get('E1_MSET'):
        grid = [float(a) for a in os.environ['E1_MSET'].split(',')]
    t0 = time.time()
    for mfl in grid:
        if mfl > MFLMAX + 1e-9 or mfl < MMIN - 1e-9:
            continue
        for J in Js:
            key = (round(float(mfl), 6), J)
            if key in db:
                continue
            st = state(Dg, Lam, mfl, J)
            if st is None:
                continue
            ts = time.time()
            m, m2, eta, beta = st['m'], st['m2'], st['eta'], st['beta']
            u = 1/mfl**2
            x = pn*np.sqrt(u); nus = Lam*pn
            es = ExactState(st['Delta'], J, DF)
            ex = es.values(nus)
            canc = {k: es.diag[k][0] for k in KS}
            Pd = remk.Pders_vec(nus, eta, nmax=0)[0]
            rec = dict(mfl=float(mfl), J=J, m=m, m2=m2, eta=eta, beta=beta, Delta=st['Delta'], x=x, nus=nus, ex=ex, canc=canc)
            rho = {}
            for k in KS:
                E = Ek(k, x)
                hyb = E*eval_legendre(J, 1 - 2*x*x) + 2*(Pd - j0((2*J + 1)*x))
                rho[k] = (ex[k] - hyb)/E
            rec['rho'] = rho
            _, W, P0 = W0.regge(Lam, beta)
            rec['W'], rec['P0'] = W, P0
            # series (OMEGA_U) remainder on the same nodes
            rec['ser'] = remk.skeleton_state(st, J, x, W, P0)
            # exact forward remainders
            hstep = min(max(0.05, 0.1/beta), 0.02*np.sqrt(m2))
            fn = np.array([1e-4, hstep, np.sqrt(2.0)*hstep])
            es2 = ExactState(st['Delta'], J, DF, tmax_extra=40.0)
            fv = es2.values(fn, ks=(4, 5))
            c4 = fv[4]
            S4 = -m2*(-3*c4[0] + 4*c4[1] - c4[2])/(2*hstep*hstep)
            R4 = S4 - (2*m2*W - 5.5)
            R5 = fv[5][0] - 2*P0
            rec['fw_exact'] = (R4, R5)
            rec['time'] = time.time() - ts
            db[key] = rec
            print(f"Dg={Dg} Df={DF:g} m/Lam={mfl:.2f} J={J:2d} m={m:.2f}  max canc {max(c.max() for c in canc.values()):.1e}"
                  f"  R4 {R4:+.3e} (ser {rec['ser']['fw'][0]:+.3e})  R5 {R5:+.3e} (ser {rec['ser']['fw'][2]:+.3e})  {rec['time']:.1f}s", flush=True)
            pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    print("# done, %d states, %.0f s" % (len(db), time.time() - t0))


if __name__ == '__main__':
    run()

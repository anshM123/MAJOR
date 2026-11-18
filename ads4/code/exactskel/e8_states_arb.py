"""
e8: exact scalar-skeleton moments on near-gap states for any Delta_g with the Arb Laplace/Barnes method (skelexact3.py; no T-contour, so
no e^{-pi nu} cancellation).  Even k = 2,4,6,8,10 exact; odd k = 3,5,7 by quartic interpolation in k of m^{2k}C_k/E_k from the even
values (calibrated at Delta_g = 100 against the exact Gamma-ratio continuation: <= 2e-3 E_k at m <= 1.1 Lam, <= 3e-5 E_k at m >= 1.3 Lam;
t_kinterp.py) [M].  Forward remainders R4, R5 from the CHMRS series (controlled at nu = 0; agree with exact to ~1e-4 at Delta_g = 100).
Output records have the same layout as e1_states.py (rho, ser, fw_exact, ...), so e3_action.py / e7_nf.py work unchanged.
Usage: python3 e8_states_arb.py Dg Df [MFLMAX] [JMAX]   (env E1_MSET, E1_JLIST as in e1)  -> cache/e1_<Dg>_Df<Df>.pkl (incremental)
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
import remk
from common import state, Transforms, KS, gl01
from u4_scan import mfl_grid
from skelexact3 import Laplace3

NP = 64
pn, pw = gl01(NP)
KE = (2, 4, 6, 8, 10)


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
    grid = mfl_grid()
    if os.environ.get('E1_MSET'):
        grid = [float(a) for a in os.environ['E1_MSET'].split(',')]
    t0 = time.time()
    for mfl in grid:
        if mfl > MFLMAX + 1e-9:
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
            L3 = Laplace3(st['Delta'], J, DF, ks=KE, numax=float(nus.max()))
            exE = {}
            rad = 0.0
            for k in KE:
                vals = [L3.value(k, nu) for nu in nus]
                exE[k] = np.array([float(v.mid()) for v in vals])
                rad = max(rad, max(float(v.rad())/max(abs(float(v.mid())), 1e-300) for v in vals))
            ex = {k: exE[k] for k in (2, 4, 6)}
            for ko in (3, 5, 7):
                r = sum(exE[ka]/Ek(ka, x)*np.prod([(ko - kb)/(ka - kb) for kb in KE if kb != ka]) for ka in KE)
                ex[ko] = r*Ek(ko, x)
            Pd = remk.Pders_vec(nus, eta, nmax=0)[0]
            rho = {}
            for k in KS:
                E = Ek(k, x)
                rho[k] = (ex[k] - E*eval_legendre(J, 1 - 2*x*x) - 2*(Pd - j0((2*J + 1)*x)))/E
            _, W, P0 = W0.regge(Lam, beta)
            ser = remk.skeleton_state(st, J, x, W, P0)
            rec = dict(mfl=float(mfl), J=J, m=m, m2=m2, eta=eta, beta=beta, Delta=st['Delta'], x=x, nus=nus, ex=ex, rho=rho,
                       canc={k: np.zeros(NP) for k in KS}, arb_relrad=rad, W=W, P0=P0, ser=ser,
                       fw_exact=(ser['fw'][0], ser['fw'][2]), method='arb-laplace (odd k interpolated)', M=L3.M)
            rec['time'] = time.time() - ts
            db[key] = rec
            print(f"Dg={Dg} Df={DF:g} m/Lam={mfl:.2f} J={J:2d} m={m:.2f}  M={L3.M} arb rel.rad {rad:.0e}  "
                  f"R4 {ser['fw'][0]:+.3e} R5 {ser['fw'][2]:+.3e} (series; ctrl {ser['fw'][4]})  {rec['time']:.1f}s", flush=True)
            pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    print("# done, %d states, %.0f s" % (len(db), time.time() - t0))


if __name__ == '__main__':
    run()

"""
dgscan/tol.py (DG_SCAN.md §2): Delta_g dependence of the per-state positivity MARGINS of Omega^EE in the hybrid near-gap model, i.e. of
everything that does not need new crystal data:
  (a) the leading-Regge + curvature baseline lambda0 = base/N_s (q3/ee1000 'hyb' states: flat NF02 + leading-Regge curvature, scalar
      remainder 0), minimised over helicity, J and m/Lam  -> the band/hybrid margin;
  (b) the failure threshold of a UNIFORM beyond-Regge graviton correction confined to nu < NU_C (ee1000_smallnu.py's C*, same C for the
      C^_2 and C^_3 parts):  C*(J, m) = -lambda0 / W_<,  W_< = sum_{nu < NU_C} [psi2 w env2 + psi3 w u env3] / (m^2 N_h);
  (c) the same for a uniform correction over ALL nu (SMOOTH_BOUND 'uniform', psi2 part only, 0.9 lambda0) and for the bulk x > X_B;
  (d) the H_LJ Regge-window tolerance T_B(c) of HLJ_ANALYTIC §2 (formula of spinfunc/lj_contour.run, re-typed here without the W_X
      term, i.e. HLJ's "T_B C2 without W_X"; lj_contour itself is not imported).
All per-state pieces come from spinfunc/ee1000_common.hyb (read-only import).  Grades: (a) [M] (model states), (b)-(d) [D] as statements
about explicit weights; their use as tolerances needs a correction model (DG_SCAN §2).
Usage: python3 tol.py [Dg,Dg,...] -> cache/tol_<Dg>.json, cache/tol.txt (appended)
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import ee1000_common as E
Q = E.Q
BSTAR = 0.97
NUCS = (8.3, 25.0, 50.0)
MF = (1.0, 1.02, 1.05, 1.1, 1.15, 1.2, 1.3)          # m/Lam grid (plus m_gap); 1.0 means "gap mass"


def jlist(Dg):
    jcap = int(0.345*(Dg - 1))                           # m_gap(J) <= 1.3 Lam
    base = list(range(0, 41)) + [44, 48, 56, 64, 80, 96, 100, 128, 150, 200, 256, 300, 344, 400, 500, 700, 1000, 1400, 1725]
    return [J for J in base if J <= jcap]


def npn_of(J):
    return 256 if J < 100 else (768 if J < 400 else 2048)


def state(Dg, J, mfl, h):
    st = E.hyb(Dg, J, mfl, h, npn_of(J))
    Lm = E.lam_of(Dg)
    x, m, u = st['x'], st['m'], st['u']
    env2 = Q.A.env(h, x, 2); env3 = Q.A.env(h, x, 3)
    w2 = st['p2']*env2/m**2/st['Ns']
    w3 = st['p3']*u*env3/m**2/st['Ns']
    nu = st['pn']*Lm
    lam0 = st['base']/st['Ns']
    r = dict(lam0=float(lam0), Ns=float(st['Ns']))
    for c in NUCS:
        W = float((w2 + w3)[nu < c].sum())
        r['Cs%g' % c] = float(-lam0/W) if W > 0 else float('inf')
    W2all = float(w2.sum())
    r['tol_all'] = float(0.9*lam0/W2all) if W2all > 0 else float('inf')
    # H_LJ window tolerance (HLJ_ANALYTIC §2), C2 part, without W_X
    F = Q.NearGap(Dg).F
    X1 = float(np.polynomial.polynomial.polyval(1.0, F.X))
    beta = 2*np.arcsinh((J + .5)/m); bs = BSTAR*beta
    xm = 1/mfl
    envm = float(Q.A.env(h, np.array([xm*0.999999]), 2)[0])
    php = Lm*bs/np.sqrt(1 - xm**2)
    Eedge = 2*120*abs(X1)*envm/(m**2*php**6*st['Ns'])
    yT = 30/bs
    pt = np.linspace(0, 1, 401) + 1j*yT/Lm
    Etop = 2*np.exp(-bs*yT)*np.max(np.abs(F.psi2(pt)*Q.A.env(h, pt/mfl, 2)))/(m**2*st['Ns'])
    for c in (3.0, 5.0):
        nuB = c*m/(2*J + 1); pB = nuB/Lm; xB = pB/mfl
        WB = float(w2[nu < nuB].sum())
        envB = float(Q.A.env(h, np.array([xB]), 2)[0])
        Ebd = 2*abs(F.psi2(np.array([pB]))[0])*envB/(m**2*Lm*bs*st['Ns'])*(1 + 1/(nuB*bs))
        r['TB%g' % c] = float(0.9*lam0/(WB + Ebd + Etop + Eedge)) if pB < 1 else float('nan')
        r['nuB%g' % c] = float(nuB)
    return r


def run(Dg):
    t0 = time.time()
    Lm = E.lam_of(Dg)
    out = dict(Dg=Dg, Lam=float(Lm), states=[])
    for J in jlist(Dg):
        if J % 2 and J < 5:
            continue                                     # no +- (helicity-4) structure below J = 4
        hs = ('+-',) if J % 2 else (('++',) if J < 4 else ('++', '+-'))
        mg = E.mgap(Dg, J)
        for mf in MF:
            mfl = mg if mf == 1.0 else mf
            if mfl < mg - 1e-12:
                continue
            for h in hs:
                r = state(Dg, J, mfl, h)
                r.update(J=J, mfl=float(mfl), h=h, gap=(mf == 1.0))
                out['states'].append(r)
    out['time'] = time.time() - t0
    json.dump(out, open(os.path.join(HERE, 'cache', 'tol_%d.json' % Dg), 'w'))
    return out


if __name__ == '__main__':
    Dgs = [int(v) for v in sys.argv[1].split(',')] if len(sys.argv) > 1 else [100, 150, 200, 250, 300, 500, 1000, 2000, 5000]
    for Dg in Dgs:
        o = run(Dg)
        S = o['states']
        mn = min(S, key=lambda r: r['lam0'])
        msg = (f"Dg={Dg} Lam={o['Lam']:.2f}: {len(S)} states, {o['time']:.0f}s; baseline min lam0 = {mn['lam0']:+.4f} "
               f"(J={mn['J']}, {mn['h']}, m/Lam={mn['mfl']:.4f})")
        print(msg, flush=True)
        open(os.path.join(HERE, 'cache', 'tol.txt'), 'a').write(msg + "\n")

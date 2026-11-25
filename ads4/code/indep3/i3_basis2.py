"""INDEP_CHECK3, claim 2: even-J structure-basis repair, constrained (identifiable) model, my own implementation.

A general A^T E B fit with A, B in GL(2) is NOT identifiable from one anchor node (i3_basis.py: residual 9e-8 but the solution
family spreads by O(1)).  The physical model: conserved(Delta, J) returns U/colmax|U| with U a real orthonormal basis of a 2-dim
null space, fixed only up to Q in O(2).  So sol_ext = U_loc Q / colmax|U_loc Q| and
   E_ext = M(Q_G)^T E_loc M(Q_H),  M(Q) = pinv(sol_loc) sol_ext(Q),
one angle + one reflection bit for the Delta basis (G) and for the 3 - Delta basis (H); separately for the state Delta and Delta_ref.
sol_loc from l5_spin.conserved (project physics library, read-only).  Fit on the anchor node (8 replicas + mean, EE),
validate on the C3H kind at the anchor and on an unfitted node (both kinds), then map all raw N5 records and compare with
BATCH5_FOLD's mapped points (data only), and recompute the failed-node set.
Usage: python3 i3_basis2.py J LOC_ANCHOR [LOC_SECOND]
"""
import sys, os, glob, pickle, json
import numpy as np
from scipy.optimize import minimize
SF = '/home/user/MAJOR/ads4/code/spinfunc'
sys.dont_write_bytecode = True; sys.path.insert(0, SF)
import warnings; warnings.filterwarnings('ignore')
cwd = os.getcwd(); os.chdir(SF)
import l5_spin as L5
os.chdir(cwd)

J = int(sys.argv[1]); D = 1000.0 + J; REF = 20000.0
RAW = '/home/user/MAJOR/ads4/ext/unpacked/out'
THEIRS = os.environ.get('I3_THEIRS', '/home/user/MAJOR/ads4/code/b5fold/cache/pts/fix')
out = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
def load(p):
    with open(p, 'rb') as fh: return pickle.load(fh)
def keyed(d): return {(round(k[0], 6), round(k[1], 6)): v for k, v in d.items()}
def raw_records(kind):
    recs = {}
    for f in glob.glob(f'{RAW}/N5_J{J}_*/n1_J{J}_*/hc_pts_{kind}_J{J}_D{D:g}.pkl'):
        recs.update(keyed(load(f)))
    return recs
def local(dirn, kind): return keyed(load(f'{dirn}/hc_pts_{kind}_J{J}_D{D:g}.pkl'))
def stack(rec): return np.concatenate([rec['reps'], rec['E'][None]], axis=0)
def rel(a, b): return np.linalg.norm(a - b)/np.linalg.norm(b)

def Mfun(sol, th, s):
    U = sol/np.linalg.norm(sol, axis=0, keepdims=True)
    Q = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]) @ np.diag([1.0, -1.0 if s else 1.0])
    Ue = U @ Q; se = Ue/np.abs(Ue).max(axis=0, keepdims=True)
    return np.linalg.pinv(sol) @ se

def fit(Eloc, Eext, solG, solH):
    def cost(p, sg, sh):
        A, B = Mfun(solG, p[0], sg), Mfun(solH, p[1], sh)
        return sum(np.linalg.norm(A.T @ Eloc[i] @ B - Eext[i])**2/np.linalg.norm(Eext[i])**2 for i in range(len(Eloc)))
    best = None; ths = np.linspace(0, 2*np.pi, 73)[:-1]
    for sg in (0, 1):
        for sh in (0, 1):
            grid = [(cost((a, b), sg, sh), a, b) for a in ths for b in ths]
            grid.sort()
            for c0, a, b in grid[:4]:
                r = minimize(lambda p: cost(p, sg, sh), [a, b], method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-30, maxiter=4000))
                cand = (r.fun, r.x[0] % (2*np.pi), r.x[1] % (2*np.pi), sg, sh)
                if best is None or cand[0] < best[0]: best = cand
            # second-best distinct solution (identifiability)
    A, B = Mfun(solG, best[1], best[3]), Mfun(solH, best[2], best[4])
    res = max(rel(A.T @ Eloc[i] @ B, Eext[i]) for i in range(len(Eloc)))
    return A, B, best, res

def main():
    loc3 = sys.argv[2]; loc6 = sys.argv[3] if len(sys.argv) > 3 else None
    maps = {}
    for Dl in (D, REF):
        solG = L5.conserved(Dl, J)[1]; solH = L5.conserved(3 - Dl, J)[1]
        raw = raw_records('EE'); L3 = local(loc3, 'EE')
        key = [k for k in L3 if abs(k[0] - Dl) < 1e-6][0]
        Eloc, Eext = stack(L3[key]), stack(raw[key])
        A, B, best, res = fit(Eloc, Eext, solG, solH)
        maps[Dl] = (A, B)
        log(f"J={J} Delta={Dl:g} anchor nu={key[1]}: unmapped |E_ext-E_loc|/|E| {rel(Eext[-1], Eloc[-1]):.2e}; fitted "
            f"theta_G={best[1]:.6f} refl {best[3]}, theta_H={best[2]:.6f} refl {best[4]}; max resid (8 reps + mean) {res:.1e}")
        log(f"   A = {np.round(A, 6).tolist()}, B = {np.round(B, 6).tolist()}")
        # other kind at the anchor, same map (not fitted)
        Lc = local(loc3, 'C3H'); rc = raw_records('C3H')
        log(f"   C3H at anchor with the EE map (unfitted kind): {rel(A.T @ Lc[key]['E'] @ B, rc[key]['E']):.1e}")
        if loc6:
            for kind in ('EE', 'C3H'):
                L6 = local(loc6, kind); r6 = raw_records(kind)
                k6 = [k for k in L6 if abs(k[0] - Dl) < 1e-6][0]
                El, Ee = L6[k6]['E'], r6[k6]['E']
                inv = np.linalg.inv(A.T) @ Ee @ np.linalg.inv(B)
                log(f"   unfitted node {kind} nu={k6[1]}: |mapped - local|/|E| = {rel(inv, El):.1e} (unmapped {rel(Ee, El):.1e};"
                    f" local MC err/|E| {np.linalg.norm(L6[k6]['err'])/np.linalg.norm(El):.1e})")
    summary = {}
    for kind in ('EE', 'C3H'):
        raw = raw_records(kind)
        tf = f'{THEIRS}/hc_pts_{kind}_J{J}_D{D:g}.pkl'
        theirs = keyed(load(tf)) if os.path.exists(tf) else {}
        wE = wR = wErr = 0; n = 0; fails = []; failsT = []
        for key, rec in raw.items():
            A, B = maps[key[0]]
            iA, iB = np.linalg.inv(A.T), np.linalg.inv(B)
            Em = iA @ rec['E'] @ iB; Rm = np.array([iA @ r @ iB for r in rec['reps']])
            err_m = np.sqrt((np.abs(Rm - Rm.mean(0))**2).sum(0)/((len(Rm) - 1)*len(Rm)))
            if key[0] == D and np.max(err_m/np.abs(Em)) > 0.2: fails.append(round(rec['nu']/rec['m'], 4))
            if key in theirs:
                t = theirs[key]; n += 1
                wE = max(wE, rel(Em, t['E'])); wR = max(wR, max(rel(Rm[i], t['reps'][i]) for i in range(len(Rm))))
                wErr = max(wErr, float(np.max(np.abs(t['err'] - err_m)/err_m)))
                if key[0] == D and np.max(t['err']/np.abs(t['E'])) > 0.2: failsT.append(round(rec['nu']/rec['m'], 4))
        log(f"{kind}: {len(raw)} raw records, {n} matched in BATCH5's mapped set: max rel diff E {wE:.1e}, replicas {wR:.1e}, err {wErr:.1e}")
        log(f"   failed (>20% elementwise) state nodes, my mapping: {sorted(fails)}; same rule on BATCH5's mapped records: {sorted(failsT)}")
        summary[kind] = dict(n=n, dE=wE, dR=wR, fails=sorted(fails))
    json.dump(summary, open(f'/home/user/MAJOR/ads4/code/indep3/cache/i3_basis2_J{J}.json', 'w'))
    open(f'/home/user/MAJOR/ads4/code/indep3/cache/i3_basis2_J{J}.txt', 'w').write("\n".join(out) + "\n")

if __name__ == '__main__':
    main()

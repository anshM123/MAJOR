"""INDEP_CHECK3, claim 2: independent re-derivation of the even-J structure-basis repair (BATCH5_FOLD §1.2).

Model (re-derived, not read from b5_basis.py): an even-J record is a 2x2 matrix E[a, c], linear in a 2-dim conserved-structure basis
for Delta (index a) and for 3 - Delta (index c).  Two machines may return different bases of the same null spaces, so
   E_ext = A^T E_loc B,   A, B in GL(2, R),   one (A, B) per (J, Delta), common to all nu and to both kinds (EE, C3H).
Fit: vec(E_ext) = K vec(E_loc) with K = B^T (x) A^T a general real 4x4 (linear least squares on the 8 replicas + mean,
complex data), then the nearest Kronecker product (Van Loan-Pitsianis), then refit residual.  Validation on an unfitted node
and on the other kind.  Then map every raw N5 record and compare with BATCH5_FOLD's mapped points (data only).
Usage: python3 i3_basis.py J LOCDIR_x03 [LOCDIR_x06]
"""
import sys, os, glob, pickle
import numpy as np

J = int(sys.argv[1]); D = 1000.0 + J; REF = 20000.0
RAW = '/home/user/MAJOR/ads4/ext/unpacked/out'
THEIRS = '/home/user/MAJOR/ads4/code/b5fold/cache/pts/all'
out = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

def load(path):
    with open(path, 'rb') as fh:
        return pickle.load(fh)

def raw_records(kind):
    recs = {}
    for f in glob.glob(f'{RAW}/N5_J{J}_*/n1_J{J}_*/hc_pts_{kind}_J{J}_D{D:g}.pkl'):
        for k, v in load(f).items():
            recs[(round(k[0], 6), round(k[1], 6))] = v
    return recs

def local(dirn, kind):
    f = glob.glob(f'{dirn}/hc_pts_{kind}_J{J}_D{D:g}.pkl')[0]
    return {(round(k[0], 6), round(k[1], 6)): v for k, v in load(f).items()}

def stack(rec):
    return np.concatenate([rec['reps'], rec['E'][None]], axis=0)   # (9, 2, 2)

def vecs(M):   # column-major vec of each 2x2
    return np.stack([m.T.reshape(-1) for m in M])

def fitK(Eloc, Eext):
    X, Y = vecs(Eloc), vecs(Eext)                      # rows: samples, complex
    Xr = np.concatenate([X.real, X.imag]); Yr = np.concatenate([Y.real, Y.imag])
    K, *_ = np.linalg.lstsq(Xr, Yr, rcond=None)        # Y = X K^T  ->  K^T
    return K.T

def nearest_kron(K):
    # K = B^T (x) A^T  ;  rearrange R[(i,j),(k,l)] = K[2i+k, 2j+l]
    R = np.zeros((4, 4))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    R[2*i+j, 2*k+l] = K[2*i+k, 2*j+l]
    U, s, Vt = np.linalg.svd(R)
    Bt = (np.sqrt(s[0])*U[:, 0]).reshape(2, 2); At = (np.sqrt(s[0])*Vt[0]).reshape(2, 2)
    return At.T, Bt.T, s

def fitAB(Eloc, Eext, nstart=12):
    from scipy.optimize import least_squares
    rng = np.random.default_rng(1)
    def resid(p):
        A = p[:4].reshape(2, 2); B = p[4:].reshape(2, 2)
        r = np.array([(A.T @ Eloc[i] @ B - Eext[i])/np.linalg.norm(Eext[i]) for i in range(len(Eloc))])
        return np.concatenate([r.real.ravel(), r.imag.ravel()])
    sols = []
    for t in range(nstart):
        p0 = rng.normal(size=8)
        sol = least_squares(resid, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        sols.append((np.max(np.abs(sol.fun)), sol.x))
    sols.sort(key=lambda z: z[0])
    best = sols[0][1]; A = best[:4].reshape(2, 2); B = best[4:].reshape(2, 2)
    Kb = np.kron(B.T, A.T)
    good = [z for z in sols if z[0] < 10*sols[0][0] + 1e-12]
    spread = max(np.linalg.norm(np.kron(z[1][4:].reshape(2, 2).T, z[1][:4].reshape(2, 2).T) - Kb)/np.linalg.norm(Kb) for z in good)
    return A, B, (len(good), spread)

def apply(A, B, E):          # E_ext = A^T E_loc B
    return A.T @ E @ B
def unmap(A, B, E):          # E_loc = A^{-T} E_ext B^{-1}
    return np.linalg.inv(A.T) @ E @ np.linalg.inv(B)

def rel(a, b):
    return np.linalg.norm(a - b)/np.linalg.norm(b)

def main():
    loc3 = sys.argv[2]; loc6 = sys.argv[3] if len(sys.argv) > 3 else None
    maps = {}
    for kind in ('EE', 'C3H'):
        raw = raw_records(kind)
        L3 = local(loc3, kind)
        for Dl in (D, REF):
            key = [k for k in L3 if abs(k[0] - Dl) < 1e-6][0]
            Eloc, Eext = stack(L3[key]), stack(raw[key])
            K = fitK(Eloc, Eext)
            A, B, s = fitAB(Eloc, Eext)
            res = max(rel(apply(A, B, Eloc[r]), Eext[r]) for r in range(9))
            resK = max(rel((K @ Eloc[r].T.reshape(-1)).reshape(2, 2).T, Eext[r]) for r in range(9))
            pre = rel(Eext[-1], Eloc[-1])
            log(f"J={J} {kind} Delta={Dl:g} nu={key[1]}: |E_ext-E_loc|/|E| unmapped {pre:.2e}; (unstructured 4x4 resid {resK:.1e},"
                f" underdetermined); structured A^T E B fit: max resid over 8 replicas + mean {res:.1e}; spread of K over {s[0]} starts {s[1]:.1e}")
            maps[(kind, Dl)] = (A, B)
    # kinds share the map?  compare K = B^T (x) A^T, scale-free
    for Dl in (D, REF):
        Ke = np.kron(maps[('EE', Dl)][1].T, maps[('EE', Dl)][0].T); Kc = np.kron(maps[('C3H', Dl)][1].T, maps[('C3H', Dl)][0].T)
        log(f"J={J} Delta={Dl:g}: |K_EE - K_C3H|/|K| = {rel(Kc, Ke):.1e}  (map shared by both kinds)")
        A, B = maps[('EE', Dl)]
        log(f"   A (state index) = {np.round(A, 6).tolist()}, B (shadow index) = {np.round(B, 6).tolist()}, "
            f"|det A|={abs(np.linalg.det(A)):.4f} |det B|={abs(np.linalg.det(B)):.4f}; A^T A = {np.round(A.T@A,4).tolist()}")
    # unfitted second node, both kinds, using EE-fitted map
    if loc6:
        for kind in ('EE', 'C3H'):
            raw = raw_records(kind); L6 = local(loc6, kind)
            for Dl in (D, REF):
                key = [k for k in L6 if abs(k[0] - Dl) < 1e-6][0]
                A, B = maps[('EE', Dl)]
                Eloc, Eext = L6[key]['E'], raw[key]['E']
                log(f"unfitted node {kind} Delta={Dl:g} nu={key[1]}: |A^-T E_ext B^-1 - E_loc|/|E| = {rel(unmap(A, B, Eext), Eloc):.1e}"
                    f" (unmapped {rel(Eext, Eloc):.1e}; MC err/|E| {np.linalg.norm(L6[key]['err'])/np.linalg.norm(Eloc):.1e})")
    # map all records with the EE map; compare to BATCH5_FOLD's mapped points (data only)
    for kind in ('EE', 'C3H'):
        raw = raw_records(kind)
        theirs = {}
        tf = f'{THEIRS}/hc_pts_{kind}_J{J}_D{D:g}.pkl'
        if os.path.exists(tf):
            theirs = {(round(k[0], 6), round(k[1], 6)): v for k, v in load(tf).items()}
        worst = 0; worst_rep = 0; nmatch = 0; worst_err = 0; fails = []
        for key, rec in raw.items():
            A, B = maps[('EE', key[0])]
            Em = unmap(A, B, rec['E']); Rm = np.array([unmap(A, B, r) for r in rec['reps']])
            err_m = np.sqrt((np.abs(Rm - Rm.mean(0))**2).sum(0)/((len(Rm) - 1)*len(Rm)))   # = hc_points' err (checked on raw)
            if key[0] == D and np.max(err_m/np.abs(Em)) > 0.2:
                fails.append(round(rec['nu']/rec['m'], 4))
            if key in theirs:
                nmatch += 1
                worst = max(worst, rel(Em, theirs[key]['E']))
                worst_rep = max(worst_rep, max(rel(Rm[r], theirs[key]['reps'][r]) for r in range(len(Rm))))
                worst_err = max(worst_err, np.max(np.abs(theirs[key]['err'] - err_m)/np.abs(err_m)))
        log(f"{kind}: {len(raw)} raw records; {nmatch} found in BATCH5 mapped set; max rel diff of mapped E vs theirs {worst:.1e};"
            f" replicas {worst_rep:.1e}; err {worst_err:.1e}")
        log(f"   {kind}: state nodes failing the elementwise 20% replica-error rule after MY mapping: {sorted(fails)}")
    np.save(f'/home/user/MAJOR/ads4/code/indep3/cache/i3_basis_J{J}.npy',
            {f'{k[0]}_{k[1]:g}': np.stack(v) for k, v in maps.items()}, allow_pickle=True)
    open(f'/home/user/MAJOR/ads4/code/indep3/cache/i3_basis_J{J}.txt', 'w').write("\n".join(out) + "\n")

def _err_like(Rm, th):
    return 0.0 * th['err']

if __name__ == '__main__':
    main()

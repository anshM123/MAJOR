"""
b4_basis (BATCH4_FOLD.md §1.2; a copy of b5fold/b5_basis.py with the rows moved to Delta = 2000 + J): the even-J structure basis of the external batch-4 records, determined EXACTLY from one local
re-computation per row (b4_points.py), and the external records mapped into the local basis.

hc_points' even record is E = sum over the integrand of F[a] Gt[c], with F linear in the columns of sol = l5_spin.conserved(Delta, J)
and Gt linear in solt = conserved(3 - Delta, J).  conserved() returns sol = U[:, :2]/colmax, U an orthonormal (SVD) basis of the real
2-dim null space -- unique only up to Q in O(2) (the two null singular values are degenerate; LAPACK-dependent).  Hence, exactly,
    sol_ext = sol_loc G,  G = diag(c_loc) Q diag(1/c_ext),  c = column max |U Q|;  likewise solt_ext = solt_loc H (Q' in O(2));
    E_ext = G^T E_loc H     (state Delta and reference Delta_ref = 2e4 separately; the same G, H for both kinds EE and C3H).
Q, Q' (one angle each + reflection bit) are fitted on the local point (8 replicas + mean; the n1 integral is deterministic, so the
residual must be at rounding level), then checked: (i) a general 2x2 G, H fit, (ii) a second local node (J = 16) predicted with
the fitted G, H.  The mapped records E_loc = G^-T E_ext H^-1 (replicas mapped one by one; 'err' recomputed from the mapped replicas as
hc_points does, std/sqrt(R)) are written to cache/pts/fix/ (all non-even files copied unchanged) -> b4_core SET 'fix'.
Usage: python3 b4_basis.py   -> cache/pts/fix/*, cache/b4_basis.txt / .json
"""
import os, sys, glob, json, pickle, shutil
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_guard  # noqa: F401
import numpy as np
from scipy.optimize import least_squares
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
sys.path.insert(0, SF)
cwd = os.getcwd(); os.chdir(SF)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import l5_spin as S
os.chdir(cwd)
B4C = os.path.join(HERE, 'cache'); PTS = os.path.join(B4C, 'pts'); LOC = os.path.join(B4C, 'loc')
DREF = 20000.0
EVEN = (12, 16, 20)
lines, res = [], {}


def pr(s=""):
    print(s, flush=True); lines.append(s)


def ubasis(Delta, J):
    """U (orthonormal real null basis, before the unit-max normalisation) and the local sol = U/colmax (conserved() replicated)."""
    L = S._l4r()
    r = L.conserved_even(float(Delta), int(J))
    sol = r['sol']
    U, sv_, _ = np.linalg.svd(np.hstack([sol.real, sol.imag]), full_matrices=False)
    U = U[:, :sol.shape[1]]
    loc = S.conserved(Delta, J)[1]
    assert np.allclose(U/np.abs(U).max(axis=0, keepdims=True), loc, atol=1e-12)
    return U


def Gof(U, th, refl):
    c, s = np.cos(th), np.sin(th)
    Q = np.array([[c, -s], [s, c]]) if not refl else np.array([[c, s], [s, -c]])
    cl = np.abs(U).max(axis=0); ce = np.abs(U @ Q).max(axis=0)
    return np.diag(cl) @ Q @ np.diag(1/ce)


def fit_GH(Ul, Ut, El, Ee):
    """E_ext = G^T E_loc H over the replica set; returns best (resid, G, H, params)."""
    sc = np.abs(Ee).max()
    best = None
    for rG in (0, 1):
        for rH in (0, 1):
            def res(p):
                G = Gof(Ul, p[0], rG); H = Gof(Ut, p[1], rH)
                d = np.array([G.T @ a @ H - b for a, b in zip(El, Ee)])/sc
                return np.r_[d.real.ravel(), d.imag.ravel()]
            for t0 in np.linspace(0, 2*np.pi, 13)[:-1]:
                for t1 in np.linspace(0, 2*np.pi, 13)[:-1]:
                    r0 = res([t0, t1])
                    if best is not None and np.sqrt(np.mean(r0**2)) > 50*best[0] + 0.5:
                        continue
                    s = least_squares(res, [t0, t1], xtol=1e-15, ftol=1e-15, gtol=1e-15)
                    rr = float(np.max(np.abs(s.fun)))
                    if best is None or rr < best[0]:
                        best = (rr, Gof(Ul, s.x[0], rG), Gof(Ut, s.x[1], rH), (float(s.x[0] % (2*np.pi)), rG, float(s.x[1] % (2*np.pi)), rH))
    return best


def fit_general(El, Ee, G0, H0):
    sc = np.abs(Ee).max()

    def res(p):
        G = p[:4].reshape(2, 2); H = p[4:].reshape(2, 2)
        d = np.array([G.T @ a @ H - b for a, b in zip(El, Ee)])/sc
        return np.r_[d.real.ravel(), d.imag.ravel(), [np.linalg.det(G) - np.linalg.det(G0)]]
    s = least_squares(res, np.r_[G0.ravel(), H0.ravel()], xtol=1e-15, ftol=1e-15, gtol=1e-15)
    return float(np.max(np.abs(s.fun)))


def recs(P, key):
    r = P[key]
    return [np.asarray(q) for q in r['reps']] + [np.asarray(r['E'])]


if __name__ == '__main__':
    os.makedirs(os.path.join(PTS, 'fix'), exist_ok=True)
    for f in glob.glob(os.path.join(PTS, 'all', '*.pkl')):
        shutil.copyfile(f, os.path.join(PTS, 'fix', os.path.basename(f)))
    pr("# b4_basis: E_ext = G^T E_loc H (G, H from O(2) rotations of the orthonormal null bases + unit-max normalisation)")
    for J in EVEN:
        D = 2000.0 + J
        fl = os.path.join(LOC, 'hc_pts_EE_J%d_D%g.pkl' % (J, D))
        if not os.path.exists(fl):
            pr(f"# J={J}: no local point -- skipped"); continue
        PL = pickle.load(open(fl, 'rb'))
        PE = {k: pickle.load(open(os.path.join(PTS, 'all', 'hc_pts_%s_J%d_D%g.pkl' % (k, J, D)), 'rb')) for k in ('EE', 'C3H')}
        nus = sorted(k[1] for k in PL if k[0] == D and (DREF, k[1]) in PL)
        T = {}
        for DD in (D, DREF):
            Ul, Ut = ubasis(DD, J), ubasis(3.0 - DD, J)
            k0 = (DD, nus[0])
            rr, G, H, par = fit_GH(Ul, Ut, recs(PL, k0), recs(PE['EE'], k0))
            rg = fit_general(recs(PL, k0), recs(PE['EE'], k0), G, H)
            T[DD] = (G, H)
            chk = []
            for nu in nus[1:]:
                k1 = (DD, nu)
                d = max(np.abs(G.T @ a @ H - b).max() for a, b in zip(recs(PL, k1), recs(PE['EE'], k1)))/np.abs(PE['EE'][k1]['E']).max()
                chk.append((round(nu, 3), float(d)))
            dE = float(np.abs(np.asarray(PL[k0]['E']) - np.asarray(PE['EE'][k0]['E'])).max()/np.abs(PE['EE'][k0]['E']).max())
            res[f'{J}:{DD:g}'] = dict(resid=rr, resid_general=rg, params=par, G=G.tolist(), H=H.tolist(), check=chk, raw_mismatch=dE)
            pr(f"# J={J} Delta={DD:g}: fitted at nu={nus[0]:.4f}: max resid (rel) {rr:.2e} (general 2x2 G, H: {rg:.2e}); theta_G={par[0]:.6f} "
               f"refl {par[1]}, theta_H={par[2]:.6f} refl {par[3]}; |E_loc - E_ext|/|E| before mapping {dE:.2e}; G = {np.round(G, 6).tolist()}, "
               f"H = {np.round(H, 6).tolist()}; second-node check {chk}")
        # map every external record of the row (state and reference, both kinds) into the local basis
        nmap = 0
        for kind in ('EE', 'C3H'):
            P = dict(PE[kind])
            for k, r in list(P.items()):
                if 'reps' not in r or np.asarray(r['E']).shape != (2, 2):
                    continue
                G, H = T[k[0]] if k[0] in T else (None, None)
                if G is None:
                    continue
                Gi, Hi = np.linalg.inv(G).T, np.linalg.inv(H)
                r = dict(r)
                reps = np.array([Gi @ np.asarray(q) @ Hi for q in r['reps']])
                r['reps'] = reps; r['E'] = reps.mean(0)
                r['err'] = reps.std(0, ddof=1)/np.sqrt(len(reps))
                r['method'] = str(r.get('method', '')) + ' [b4_basis: mapped to the local structure basis]'
                P[k] = r; nmap += 1
            pickle.dump(P, open(os.path.join(PTS, 'fix', 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)), 'wb'))
        # check: mapped vs local at the local nodes
        PF = pickle.load(open(os.path.join(PTS, 'fix', 'hc_pts_EE_J%d_D%g.pkl' % (J, D)), 'rb'))
        dd = max(float(np.abs(np.asarray(PF[k]['E']) - np.asarray(PL[k]['E'])).max()/np.abs(PL[k]['E']).max()) for k in PL if k in PF)
        pr(f"   mapped {nmap} records; max |E_mapped - E_local|/|E| at the local nodes: {dd:.2e}")
        res[f'{J}:mapped'] = dict(n=nmap, chk=dd)
    open(os.path.join(B4C, 'b4_basis.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(B4C, 'b4_basis.json'), 'w'), indent=1)

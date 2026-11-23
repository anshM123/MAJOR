"""
pcyc_cv (PCYC_CLOSE.md section 3.2): zero-variance control variates for the crystal-cycle estimate of <Q>_crys (pcyc_crys runs with
stored points).  On ANY closed twisted cycle Gamma (crystal thimbles included; nothing about P_cyc is assumed)
    int_Gamma d^7w  d/dw_j ( psi e^Phi ) = 0     =>     CV_{j,psi} := d psi/dw_j + psi dPhi/dw_j   has  int_Gamma e^Phi CV = 0,
Phi = the integrand's exponent in the w-chart (rd_core.Flow.grad: exact holomorphic gradient incl. the chart Jacobian), psi = monomials
of degree <= DEG in (z_j/c_j - 1 for the four light-cone variables, x~0, x_p/a, y_p/a).  The residual Q - q - sum beta_k CV_k is fitted
by weighted complex least squares (weights |w F|) on the OTHER replicas (cross-fitting: no fit bias), and
    <Q>_crys = sum_th e^{logN} sum w F (Q - beta.CV)  /  sum_th e^{logN} sum w F (1 - gamma.CV)
with replica jackknife errors (gamma: the same fit for the constant function, reducing the denominator's variance).
Usage: python3 pcyc_cv.py cache/cr_....pkl [DEG] -> printed + appended to cache/pcyc_cv.txt
"""
import os, sys, pickle, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import rd_core as C
np.seterr(all='ignore')


def basis(W, c, Delta, deg):
    """psi monomials and their w-derivatives: returns list of (psi (N,), dpsi (7, N))."""
    a = 1.0/Delta
    v = np.zeros_like(W); dv = np.zeros_like(W)
    for j in (1, 2, 3, 4):
        e = np.exp(W[j]); v[j] = e - 1.0; dv[j] = e
    t = np.tanh(W[0]); v[0] = t; dv[0] = 1 - t*t
    for j in (5, 6):
        v[j] = W[j]/a; dv[j] = 1.0/a
    out = []
    for d in range(deg + 1):
        for combo in itertools.combinations_with_replacement(range(7), d):
            n = np.bincount(np.array(combo, int), minlength=7) if d else np.zeros(7, int)
            psi = np.ones(W.shape[1], complex)
            for j in range(7):
                if n[j]:
                    psi = psi*v[j]**n[j]
            dpsi = np.zeros_like(W)
            for j in range(7):
                if n[j]:
                    rest = np.ones(W.shape[1], complex)
                    for k in range(7):
                        if k != j and n[k]:
                            rest = rest*v[k]**n[k]
                    dpsi[j] = n[j]*v[j]**(n[j] - 1)*dv[j]*rest
            out.append((psi, dpsi))
    return out


def cvs(th, nu, Delta, deg):
    P = th['pts']
    W = np.concatenate([p['W'] for p in P], axis=0).T              # (7, N)
    y0 = np.concatenate([p['y0'] for p in P])
    F = np.concatenate([p['F'] for p in P]); Q = np.concatenate([p['Q'] for p in P])
    wg = np.concatenate([p['wgt'] for p in P]); rr = np.concatenate([np.full(p['F'].size, p['r']) for p in P])
    ok = np.concatenate([p['ok'] for p in P])
    Fl = C.Flow(nu, Delta, th['c'])
    G = np.zeros_like(W)
    for c0 in range(0, W.shape[1], 4096):
        g, _ = Fl.grad(W[:, c0:c0 + 4096], y0[c0:c0 + 4096])
        G[:, c0:c0 + 4096] = g
    cols = []
    for psi, dpsi in basis(W, th['c'], Delta, deg):
        for j in range(7):
            cols.append(dpsi[j] + psi*G[j])
    X = np.array(cols).T                                            # (N, K)
    X[~ok] = 0.0
    good = np.all(np.isfinite(X), axis=1) & ok
    X[~good] = 0.0
    return dict(F=np.where(good, F, 0), Q=np.where(good, Q, 0), wg=wg, r=rr, X=X, good=good)


def fit(D, target, train):
    """complex weighted LS: target ~ q + X beta on rows train (weights |wF|); returns beta (without q)."""
    a = np.abs(D['wg']*D['F'])[train]
    A = np.concatenate([np.ones((train.sum(), 1)), D['X'][train]], axis=1)*a[:, None]
    b = target[train]*a
    # column scaling + truncated SVD (the CV columns are nearly collinear)
    s = np.sqrt(np.sum(np.abs(A)**2, axis=0)); s[s == 0] = 1.0
    U, sv, Vh = np.linalg.svd(A/s[None, :], full_matrices=False)
    keep = sv > 1e-10*sv[0]
    coef = (Vh[keep].conj().T @ ((U[:, keep].conj().T @ b)/sv[keep]))/s
    return coef[1:]


def estimate(res, deg=2, verbose=True):
    nu, Delta, R = res['nu'], res['Delta'], res['R']
    l0 = res['logN0']
    Ds = {k: cvs(t, nu, Delta, deg) for k, t in res['th'].items()}
    num = np.zeros(R, complex); den = np.zeros(R, complex); num0 = np.zeros(R, complex); den0 = np.zeros(R, complex)
    zero_chk = []
    for k, D in Ds.items():
        fac = np.exp(res['th'][k]['logN'] - l0)
        for r in range(R):
            tr = (D['r'] != r) & D['good']; te = D['r'] == r
            bQ = fit(D, D['Q'], tr)
            b1 = fit(D, np.ones_like(D['Q']), tr)
            wF = (D['wg']*D['F'])[te]; Xt = D['X'][te]
            Nr = te.sum()
            num[r] += fac*np.sum(wF*(D['Q'][te] - Xt @ bQ))/Nr
            den[r] += fac*np.sum(wF*(1.0 - Xt @ b1))/Nr
            num0[r] += fac*np.sum(wF*D['Q'][te])/Nr
            den0[r] += fac*np.sum(wF)/Nr
        # identity check: the CVs integrate to ~0 (relative to the scale of the integrand)
        wF = D['wg']*D['F']
        zc = np.abs(np.sum(wF[:, None]*D['X'], axis=0))/np.sqrt(np.sum(np.abs(wF[:, None]*D['X'])**2, axis=0) + 1e-300)
        zero_chk.append(float(np.median(zc)))

    def jk(n, d):
        v = np.array([(n.sum() - n[i])/(d.sum() - d[i]) for i in range(R)])
        return n.sum()/d.sum(), float(np.sqrt((R - 1)/R*np.sum(np.abs(v - v.mean())**2)))
    q0, e0 = jk(num0, den0)
    q1, e1 = jk(num, den)
    if verbose:
        print(f'D={Delta:g} nu={nu:g} deg={deg}: plain <Q> = {q0:.6e} +- {e0:.1e} (rel {e0/abs(q0):.1e});  CV <Q> = {q1:.6e} +- {e1:.1e} '
              f'(rel {e1/abs(q1):.1e});  CV sum check (median |sum|/sqrt(sum^2), ~1 for pure noise) {zero_chk}', flush=True)
    return dict(q0=q0, e0=e0, q1=q1, e1=e1, deg=deg)


if __name__ == '__main__':
    res = pickle.load(open(sys.argv[1], 'rb'))
    deg = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    o = estimate(res, deg)
    with open(os.path.join(HERE, 'cache', 'pcyc_cv.txt'), 'a') as f:
        f.write(f"{os.path.basename(sys.argv[1])} deg={deg}: plain {o['q0']:.8e} +- {o['e0']:.2e}; CV {o['q1']:.8e} +- {o['e1']:.2e}\n")

"""Second-order (infinitesimal) form of the CMI Petz conjecture with classical A.

For sigma>0 on B(x)C, N=tr_C, P = Petz map of (sigma, N), delta Hermitian:
  Q1(delta) = g^BKM_sigma(delta,delta) - g^BKM_{sigma_B}(delta_B,delta_B)
  Q2(delta) = (1/2) g^SLD_sigma(eta,eta),  eta = delta - P(delta_B)
The state rho_ABC = sum_{a=+-} (1/2)|a><a| (x) (sigma + a*eps*delta) has
  I(A;C|B) = (eps^2/2) Q1 + O(eps^3),   -log F(rho,Petz(rho_AB)) = (eps^2/2) Q2 + O(eps^3),
so Q1 >= Q2 is necessary for the Berta-Seshadreesan-Wilde conjecture.
We compute the smallest generalised eigenvalue of (Q1,Q2) and minimise it over sigma."""
import numpy as np, sys
from scipy.optimize import minimize
from scipy.linalg import eigh

def herm_basis(D):
    B = []
    for i in range(D):
        E = np.zeros((D, D), complex); E[i, i] = 1; B.append(E)
    for i in range(D):
        for j in range(i+1, D):
            E = np.zeros((D, D), complex); E[i, j] = E[j, i] = 1/np.sqrt(2); B.append(E)
            E = np.zeros((D, D), complex); E[i, j] = -1j/np.sqrt(2); E[j, i] = 1j/np.sqrt(2); B.append(E)
    return B

def metric_kernel(lam, kind):
    x = lam[:, None]; y = lam[None, :]
    if kind == "sld": return 2/(x+y)
    with np.errstate(divide='ignore', invalid='ignore'):
        k = (np.log(x)-np.log(y))/(x-y)
    k[np.abs(x-y) < 1e-12*(x+y)] = (2/(x+y))[np.abs(x-y) < 1e-12*(x+y)]
    return k

def gram(sig, vecs, kind):
    lam, V = np.linalg.eigh(sig); K = metric_kernel(lam, kind)
    T = [V.conj().T@X@V for X in vecs]
    n = len(T); G = np.zeros((n, n))
    W = [K*t for t in T]
    for a in range(n):
        for b in range(a, n):
            G[a, b] = G[b, a] = np.real(np.vdot(T[a], W[b]))
    return G

def ptrC(X, dB, dC): return np.einsum('icjc->ij', X.reshape(dB, dC, dB, dC))

def forms(sig, dB, dC):
    D = dB*dC; Bs = herm_basis(D)
    sB = ptrC(sig, dB, dC)
    w, V = np.linalg.eigh(sig); sq = (V*np.sqrt(w))@V.conj().T
    w2, V2 = np.linalg.eigh(sB); isq = (V2/np.sqrt(w2))@V2.conj().T
    K = sq@np.kron(isq, np.eye(dC))
    dBs = [ptrC(X, dB, dC) for X in Bs]
    eta = [X - K@np.kron(Y, np.eye(dC))@K.conj().T for X, Y in zip(Bs, dBs)]
    Q1 = gram(sig, Bs, "bkm") - gram(sB, dBs, "bkm")
    Q2 = 0.5*gram(sig, eta, "sld")
    return Q1, Q2, Bs

def min_ratio(sig, dB, dC, return_vec=False):
    Q1, Q2, Bs = forms(sig, dB, dC)
    w, U = np.linalg.eigh(Q2); keep = w > 1e-9*w.max()
    Ur = U[:, keep]/np.sqrt(w[keep])
    M = Ur.T@Q1@Ur; ev, evec = np.linalg.eigh((M+M.T)/2)
    if return_vec:
        c = Ur@evec[:, 0]; return ev[0], sum(ci*X for ci, X in zip(c, Bs))
    return ev[0]

def sig_from(x, D):
    L = (x[:D*D]+1j*x[D*D:]).reshape(D, D); s = L@L.conj().T; return s/np.trace(s).real

if __name__ == "__main__":
    dB, dC = int(sys.argv[1]), int(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4]); D = dB*dC
    rng = np.random.default_rng(seed)
    vals = [min_ratio(sig_from(rng.normal(size=2*D*D), D), dB, dC) for _ in range(300)]
    print(f"dB={dB} dC={dC}: random sigma: min generalised eigenvalue over 300 samples = {min(vals):.5f} (median {np.median(vals):.4f})", flush=True)
    best = (np.inf, None)
    for t in range(trials):
        res = minimize(lambda x: min_ratio(sig_from(x, D), dB, dC), rng.normal(size=2*D*D)*rng.choice([0.5, 1, 2]),
                       method="Nelder-Mead" if D <= 2 else "L-BFGS-B", options=dict(maxiter=300))
        if res.fun < best[0]: best = (res.fun, res.x.copy())
        print(f"   trial {t}: {res.fun:.6f}", flush=True)
    s = sig_from(best[1], D)
    print(f"dB={dB} dC={dC}: optimised min generalised eigenvalue = {best[0]:.6f}; eig(sigma) = {np.linalg.eigvalsh(s)}")
    np.save(f"/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/inf_{dB}{dC}_{seed}.npy", best[1])

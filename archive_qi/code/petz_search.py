"""Adversarial search for counterexamples to the (unrotated) Petz recovery conjectures.

(CMI)  I(A;C|B)_rho >= -log F(rho_ABC, P(rho_AB)),  P(X) = rho_BC^{1/2} rho_B^{-1/2} X rho_B^{-1/2} rho_BC^{1/2}
(GEN)  D(rho||sigma) - D(tr_E rho || tr_E sigma) >= -log F(rho, P_sigma(tr_E rho)),
       P_sigma(X) = sigma^{1/2} (sigma_B^{-1/2} X sigma_B^{-1/2} (x) 1_E) sigma^{1/2}
F(r,s) = ||sqrt r sqrt s||_1^2.  We minimise  gap = LHS - RHS  (natural log) and the ratio LHS/RHS.
"""
import numpy as np, sys
from scipy.optimize import minimize

def herm_fn(M, f):
    w, V = np.linalg.eigh((M + M.conj().T) / 2)
    return (V * f(w)) @ V.conj().T

def psqrt(M):  return herm_fn(M, lambda w: np.sqrt(np.clip(w, 0, None)))
def pisqrt(M): return herm_fn(M, lambda w: 1 / np.sqrt(np.clip(w, 1e-300, None)))
def plog(M):   return herm_fn(M, lambda w: np.log(np.clip(w, 1e-300, None)))

def ent(M):
    w = np.linalg.eigvalsh((M + M.conj().T) / 2); w = w[w > 1e-300]
    return float(-(w * np.log(w)).sum())

def fid(r, s):
    sr = psqrt(r)
    w = np.linalg.eigvalsh(sr @ s @ sr)
    return float(np.sqrt(np.clip(w, 0, None)).sum() ** 2)

def ptr(M, dims, keep):
    n = len(dims); M = M.reshape(dims + dims)
    idx = list(range(n)); cur = n
    for i in sorted(set(idx) - set(keep), reverse=True):
        M = np.trace(M, axis1=i, axis2=i + cur); cur -= 1
    d = int(np.prod([dims[k] for k in keep]))
    return M.reshape(d, d)

def state_from(x, D, rank=None):
    r = D if rank is None else rank
    L = (x[:D * r] + 1j * x[D * r:2 * D * r]).reshape(D, r)
    rho = L @ L.conj().T
    return rho / np.trace(rho).real

def cmi_gap(rho, dims):
    dA, dB, dC = dims
    rAB = ptr(rho, dims, [0, 1]); rBC = ptr(rho, dims, [1, 2]); rB = ptr(rho, dims, [1])
    I = ent(rAB) + ent(rBC) - ent(rB) - ent(rho)
    K = np.kron(np.eye(dA), psqrt(rBC) @ np.kron(pisqrt(rB), np.eye(dC)))
    rec = K @ np.kron(rAB, np.eye(dC)) @ K.conj().T
    F = fid(rho, rec)
    return I, -np.log(F)

def gen_gap(rho, sig, dims):
    dB, dE = dims
    rB = ptr(rho, dims, [0]); sB = ptr(sig, dims, [0])
    def D(r, s): return float(np.trace(r @ (plog(r) - plog(s))).real)
    lhs = D(rho, sig) - D(rB, sB)
    K = psqrt(sig) @ np.kron(pisqrt(sB), np.eye(dE))
    rec = K @ np.kron(rB, np.eye(dE)) @ K.conj().T
    return lhs, -np.log(fid(rho, rec))

def run_cmi(dims, trials, seed, rank=None, mode="gap"):
    rng = np.random.default_rng(seed); D = int(np.prod(dims)); r = D if rank is None else rank
    best = (np.inf, None)
    def obj(x):
        I, mlF = cmi_gap(state_from(x, D, rank), dims)
        return (I - mlF) if mode == "gap" else (I / mlF if mlF > 1e-13 else 10.0)
    for t in range(trials):
        x0 = rng.normal(size=2 * D * r) * rng.choice([1, 0.3, 3])
        res = minimize(obj, x0, method="L-BFGS-B", options=dict(maxiter=400))
        if res.fun < best[0]: best = (res.fun, res.x.copy())
    I, mlF = cmi_gap(state_from(best[1], D, rank), dims)
    return best[0], I, mlF, best[1]

def run_gen(dims, trials, seed, mode="gap"):
    rng = np.random.default_rng(seed); D = int(np.prod(dims)); n = 2 * D * D
    best = (np.inf, None)
    def obj(x):
        lhs, mlF = gen_gap(state_from(x[:n], D), state_from(x[n:], D), dims)
        return (lhs - mlF) if mode == "gap" else (lhs / mlF if mlF > 1e-13 else 10.0)
    for t in range(trials):
        x0 = rng.normal(size=2 * n) * rng.choice([1, 0.3, 3])
        res = minimize(obj, x0, method="L-BFGS-B", options=dict(maxiter=400))
        if res.fun < best[0]: best = (res.fun, res.x.copy())
    lhs, mlF = gen_gap(state_from(best[1][:n], D), state_from(best[1][n:], D), dims)
    return best[0], lhs, mlF, best[1]

if __name__ == "__main__":
    kind = sys.argv[1]; dims = tuple(int(c) for c in sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4])
    mode = sys.argv[5] if len(sys.argv) > 5 else "gap"
    rank = int(sys.argv[6]) if len(sys.argv) > 6 else None
    if kind == "cmi": val, a, b, x = run_cmi(dims, trials, seed, rank, mode)
    else:             val, a, b, x = run_gen(dims, trials, seed, mode)
    print(f"{kind} dims={dims} mode={mode} rank={rank} best={val:.6e} LHS={a:.6e} -logF={b:.6e}", flush=True)
    np.save(f"/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/petz_{kind}_{sys.argv[2]}_{mode}_{rank}_{seed}.npy", x)

import numpy as np, sys
from scipy.optimize import minimize
from petz_infinitesimal import *
from petz_search import cmi_gap
dB, dC, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = dB*dC
rng = np.random.default_rng(seed)
def rand_sig():
    # random spectrum with controllable spread, random unitary
    G = rng.normal(size=(D, D))+1j*rng.normal(size=(D, D)); U, _ = np.linalg.qr(G)
    lam = np.exp(rng.normal(size=D)*rng.choice([0.5, 1.5, 3.0])); lam /= lam.sum()
    return (U*lam)@U.conj().T
pool = []
for _ in range(3000):
    s = rand_sig(); pool.append((min_ratio(s, dB, dC), s))
pool.sort(key=lambda t: t[0])
print(f"({dB},{dC}) sampled: lowest 5 =", [round(p[0], 4) for p in pool[:5]], flush=True)
def param(s):
    L = np.linalg.cholesky(s+1e-15*np.eye(D)); return np.concatenate([L.real.ravel(), L.imag.ravel()])
MINEIG = 1e-4   # keep sigma well-conditioned so that every number below is numerically trustworthy
def obj(x):
    s = sig_from(x, D); s = (1-D*MINEIG)*s + MINEIG*np.eye(D)
    return min_ratio(s, dB, dC)
best = (np.inf, None)
for val, s in pool[:4]:
    res = minimize(obj, param(s), method="L-BFGS-B", options=dict(maxiter=500))
    print(f"   start {val:.4f} -> {res.fun:.6f}", flush=True)
    if res.fun < best[0]: best = (res.fun, res.x.copy())
s = sig_from(best[1], D); s = (1-D*MINEIG)*s + MINEIG*np.eye(D)
r, delta = min_ratio(s, dB, dC, return_vec=True)
print(f"best ratio {r:.6f}; eig(sigma)={np.linalg.eigvalsh(s)}")
# direct finite-epsilon check with the actual cq state (float64)
delta = delta - np.trace(delta)*s; delta /= np.linalg.norm(delta)
for eps in [1e-2, 3e-3, 1e-3]:
    rho = np.zeros((2*D, 2*D), complex)
    rho[:D, :D] = 0.5*(s+eps*delta); rho[D:, D:] = 0.5*(s-eps*delta)
    if np.linalg.eigvalsh(s+eps*delta)[0] < 0 or np.linalg.eigvalsh(s-eps*delta)[0] < 0: print("  eps too large", eps); continue
    I, mlF = cmi_gap(rho, (2, dB, dC)); print(f"  eps={eps:g}: I={I:.6e}  -logF={mlF:.6e}  ratio={I/mlF:.5f}")
np.save(f"/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/infx_{dB}{dC}_{seed}.npy", np.concatenate([s.ravel(), delta.ravel()]))

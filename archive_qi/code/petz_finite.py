"""Finite (non-perturbative) optimisation of the CMI-form violation for cq states:
rho_ABC = sum_a p_a |a><a| (x) rho^a_BC,  minimise ratio I/(-log F) or gap I - (-log F)."""
import numpy as np, sys
from scipy.optimize import minimize
from petz_search import cmi_gap
S="/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/"
def build(x, nA, D, floor):
    blocks = []; 
    for a in range(nA):
        L = (x[a*2*D*D:a*2*D*D+D*D] + 1j*x[a*2*D*D+D*D:(a+1)*2*D*D]).reshape(D, D)
        blocks.append(L@L.conj().T)
    tot = sum(np.trace(b).real for b in blocks)
    rho = np.zeros((nA*D, nA*D), complex)
    for a in range(nA): rho[a*D:(a+1)*D, a*D:(a+1)*D] = (1-floor*nA*D)*blocks[a]/tot + floor*np.eye(D)
    return rho
def run(nA, dB, dC, x0, mode, floor):
    D = dB*dC
    def obj(x):
        I, mlF = cmi_gap(build(x, nA, D, floor), (nA, dB, dC))
        return I/mlF if mode == "ratio" else (I-mlF)
    res = minimize(obj, x0, method="L-BFGS-B", options=dict(maxiter=3000, maxfun=400000, ftol=1e-15, gtol=1e-12))
    rho = build(res.x, nA, D, floor); I, mlF = cmi_gap(rho, (nA, dB, dC))
    return res.x, rho, I, mlF
if __name__ == "__main__":
    mode = sys.argv[1]; floor = float(sys.argv[2]); eps = float(sys.argv[3])
    v = np.load(S+"infx_22_0.npy"); s = v[:16].reshape(4, 4); d = v[16:].reshape(4, 4)
    x0 = []
    for sign in (+1, -1):
        L = np.linalg.cholesky(0.5*(s+sign*eps*d)); x0 += [L.real.ravel(), L.imag.ravel()]
    x, rho, I, mlF = run(2, 2, 2, np.concatenate(x0), mode, floor)
    print(f"mode={mode} floor={floor:g} eps0={eps:g}: I={I:.6e} -logF={mlF:.6e} gap={I-mlF:.4e} ratio={I/mlF:.5f} mineig={np.linalg.eigvalsh(rho)[0]:.3e}", flush=True)
    np.save(S+f"fin_{mode}_{floor:g}_{eps:g}.npy", rho)

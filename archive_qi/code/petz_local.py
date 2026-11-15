"""Targeted search near equality points of the CMI-form Petz conjecture.
gap(rho) = I(A;C|B) + log F(rho, Petz(rho_AB)) ; conjecture: gap >= 0.
Equality holds at Phi_AC (x) tau_B (ratio 1). We perturb around such points (and around
partially entangled / embedded variants) and minimise gap directly."""
import numpy as np, sys
from scipy.optimize import minimize
from petz_search import cmi_gap, state_from

def seed_state(dims, rng, kind):
    dA, dB, dC = dims; D = dA*dB*dC
    k = min(dA, dC)
    lam = np.ones(k) if kind == 0 else rng.dirichlet(np.ones(k))*k
    phi = np.zeros((dA, dC), complex)
    for i in range(k): phi[i, i] = np.sqrt(lam[i]/k)
    b = rng.normal(size=dB) + 1j*rng.normal(size=dB); b /= np.linalg.norm(b)
    psi = np.einsum('ac,b->abc', phi, b).reshape(D)
    if kind == 2:  # entangle B with AC a bit through a random unitary-ish mixing
        G = rng.normal(size=(D, D)) + 1j*rng.normal(size=(D, D)); H = G + G.conj().T
        w, V = np.linalg.eigh(H); psi = (V*np.exp(1j*0.3*w/np.abs(w).max()))@V.conj().T@psi
    return psi

def run(dims, trials, seed, rank, eps_list=(1e-1, 3e-2, 1e-2, 3e-1)):
    rng = np.random.default_rng(seed); D = int(np.prod(dims)); best = (np.inf, None, None)
    def obj(x):
        I, mlF = cmi_gap(state_from(x, D, rank), dims); return I - mlF
    for t in range(trials):
        psi = seed_state(dims, rng, t % 3); eps = eps_list[t % len(eps_list)]
        L = np.zeros((D, rank), complex); L[:, 0] = psi
        L += eps*(rng.normal(size=(D, rank)) + 1j*rng.normal(size=(D, rank)))/np.sqrt(D)
        x0 = np.concatenate([L.real.ravel(), L.imag.ravel()])
        res = minimize(obj, x0, method="L-BFGS-B", options=dict(maxiter=300, ftol=1e-15, gtol=1e-10))
        I, mlF = cmi_gap(state_from(res.x, D, rank), dims)
        if res.fun < best[0]: best = (res.fun, (I, mlF), res.x.copy())
        if res.fun < -1e-7:
            print("  NEGATIVE GAP", dims, rank, res.fun, I, mlF, flush=True)
    return best

if __name__ == "__main__":
    dims = tuple(int(c) for c in sys.argv[1]); rank = int(sys.argv[2]); trials = int(sys.argv[3]); seed = int(sys.argv[4])
    val, (I, mlF), x = run(dims, trials, seed, rank)
    print(f"local dims={dims} rank={rank} min gap={val:.3e} at I={I:.6f} -logF={mlF:.6f}", flush=True)
    np.save(f"/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/petzloc_{sys.argv[1]}_{rank}_{seed}.npy", x)

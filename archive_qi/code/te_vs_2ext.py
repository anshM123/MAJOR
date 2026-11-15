"""Compare two convex sets of bipartite states rho on A (x) E:
  TE  : rho = omega_{AE} for some state omega_{AEE'} with omega_{AE'} = (omega_{AE})^{T_E}   ("transpose-extendible")
  EXT2: rho = sigma_{AE} for some state sigma_{AEE'} with sigma_{AE'} = sigma_{AE}            (symmetric 2-extendible)
A transpose-degradable channel that is not degradable exists iff TE is not contained in EXT2
(rho = Choi state of the complementary channel; see notes).  We compare support functions
h_K(W) = max_{rho in K} tr(W rho) for random Hermitian W:  h_TE(W) > h_EXT2(W) certifies TE not in EXT2."""
import numpy as np, cvxpy as cp, sys

def pt_E(M, dA, dE):   # partial transpose on second factor of a cvxpy/numpy dA*dE matrix
    return cp.partial_transpose(M, (dA, dE), axis=1)

def swap_EEp(dA, dE):
    D = dA*dE*dE; P = np.zeros((D, D))
    for a in range(dA):
        for e in range(dE):
            for f in range(dE):
                P[a*dE*dE + f*dE + e, a*dE*dE + e*dE + f] = 1
    return P

def h_TE(W, dA, dE, solver="CLARABEL"):
    D = dA*dE*dE; om = cp.Variable((D, D), hermitian=True)
    rAE = cp.partial_trace(om, (dA*dE, dE), axis=1)
    P = swap_EEp(dA, dE); rAEp = cp.partial_trace(P@om@P.T, (dA*dE, dE), axis=1)
    cons = [om >> 0, cp.real(cp.trace(om)) == 1, rAEp == pt_E(rAE, dA, dE)]
    prob = cp.Problem(cp.Maximize(cp.real(cp.trace(W@rAE))), cons); prob.solve(solver=solver)
    return prob.value, om.value

def h_EXT2(W, dA, dE, solver="CLARABEL"):
    D = dA*dE*dE; om = cp.Variable((D, D), hermitian=True)
    rAE = cp.partial_trace(om, (dA*dE, dE), axis=1)
    P = swap_EEp(dA, dE); rAEp = cp.partial_trace(P@om@P.T, (dA*dE, dE), axis=1)
    cons = [om >> 0, cp.real(cp.trace(om)) == 1, rAEp == rAE]
    prob = cp.Problem(cp.Maximize(cp.real(cp.trace(W@rAE))), cons); prob.solve(solver=solver)
    return prob.value, om.value

def dist_to_EXT2(rho, dA, dE, solver="CLARABEL"):
    """min ||rho - sigma_AE||_F over symmetric extensions; returns value and the optimal separating direction"""
    D = dA*dE*dE; om = cp.Variable((D, D), hermitian=True)
    rAE = cp.partial_trace(om, (dA*dE, dE), axis=1)
    P = swap_EEp(dA, dE); rAEp = cp.partial_trace(P@om@P.T, (dA*dE, dE), axis=1)
    cons = [om >> 0, cp.real(cp.trace(om)) == 1, rAEp == rAE]
    prob = cp.Problem(cp.Minimize(cp.norm(rho - rAE, "fro")), cons); prob.solve(solver=solver)
    sig = cp.partial_trace(om, (dA*dE, dE), axis=1).value
    return prob.value, rho - sig

if __name__ == "__main__":
    dA, dE, n, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = np.random.default_rng(seed); d = dA*dE; best = (-1, None)
    for t in range(n):
        G = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)); W = (G + G.conj().T)/2
        # see-saw: improve the direction using the separating hyperplane of the TE-optimal point
        for it in range(6):
            v1, om = h_TE(W, dA, dE); rho = np.einsum('ijkj->ik', om.reshape(d, dE, d, dE))
            dist, Wn = dist_to_EXT2(rho, dA, dE)
            if dist < 1e-7: break
            W = (Wn + Wn.conj().T)/2/np.linalg.norm(Wn)
        v1, om = h_TE(W, dA, dE); v2, _ = h_EXT2(W, dA, dE)
        rho = np.einsum('ijkj->ik', om.reshape(d, dE, d, dE)); dist, _ = dist_to_EXT2(rho, dA, dE)
        print(f"dA={dA} dE={dE} trial {t}: h_TE - h_EXT2 = {v1-v2:+.3e}; distance of TE-optimal state to EXT2 = {dist:.3e}" + ("   <<< TE not in EXT2" if dist > 1e-5 else ""), flush=True)
        if dist > best[0]: best = (dist, om)
    np.save(f"/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/te_{dA}{dE}_{seed}.npy", best[1])

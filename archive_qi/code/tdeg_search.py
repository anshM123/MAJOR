"""Search for a transpose-degradable (a.k.a. conjugate-degradable) channel that is not degradable.

WLOG (Stinespring-dilating the degrading map) a transpose-degradable channel is given by an isometry
U: A -> E' (x) F (x) E  with  tr_{FE}(U X U^dag) = [tr_{E'F}(U X U^dag)]^T  for all X,   (TD)
channel Phi = tr_E (output B = E'F), complement Phi^c = tr_{E'F} (output E).
Step 1: find isometries satisfying (TD) by minimising the squared residual.
Step 2: test ordinary degradability of Phi by an SDP (least-squares defect over CPTP D: B->E)."""
import numpy as np, sys, cvxpy as cp
from scipy.optimize import minimize

def iso_from(x, dA, dout):
    M = (x[:dout*dA] + 1j*x[dout*dA:]).reshape(dout, dA)
    Q, R = np.linalg.qr(M); return Q*(np.diag(R)/np.abs(np.diag(R)))
def marginals(U, dA, dE, dF):
    T = U.reshape(dE, dF, dE, dA)                       # (e', f, e, a)
    # Choi-like tensors: R1[a,a',e1,e2] = sum_{f,e} T[e1,f,e,a] conj T[e2,f,e,a'] ; R2 over e',f
    R1 = np.einsum('xfea,yfeb->abxy', T, T.conj())
    R2 = np.einsum('xfea,xfgb->abeg', T, T.conj())
    return R1, R2
def residual(x, dA, dE, dF):
    U = iso_from(x, dA, dE*dF*dE); R1, R2 = marginals(U, dA, dE, dF)
    return np.linalg.norm(R1 - R2.transpose(0, 1, 3, 2))**2
def choi(K_out_from_iso):  # helper not used
    pass
def degradability_defect(U, dA, dE, dF, transpose=False, solver="CLARABEL"):
    """min over CPTP D: B->E of || D o Phi - Phi^c (or T o Phi^c) ||_F  on the Choi level"""
    dB = dE*dF; T = U.reshape(dB, dE, dA)               # (b, e, a)
    # Phi(|a><a'|) = sum_e T[:,e,a] T[:,e,a']^dag ; Phi^c(|a><a'|) = sum_b T[b,:,a] T[b,:,a']^dag
    PhiO = np.einsum('bea,cef->afbc', T, T.conj())      # [a,a',b,b']
    PhiC = np.einsum('bea,bgf->afeg', T, T.conj())      # [a,a',e,e']
    if transpose: PhiC = PhiC.transpose(0, 1, 3, 2)
    J = cp.Variable((dB*dE, dB*dE), hermitian=True)     # Choi of D: J = sum |b><b'| (x) D(|b><b'|)
    cons = [J >> 0, cp.partial_trace(J, (dB, dE), axis=1) == np.eye(dB)]
    # D(Y) = tr_B[(Y^T (x) 1) J]
    expr = []
    for a in range(dA):
        for ap in range(dA):
            Y = PhiO[a, ap]
            out = cp.partial_trace(cp.kron(Y.T, np.eye(dE)) @ J, (dB, dE), axis=0)
            expr.append(cp.vec(out - PhiC[a, ap], order="F"))
    obj = cp.norm(cp.hstack(expr), 2)
    prob = cp.Problem(cp.Minimize(obj), cons); prob.solve(solver=solver)
    return prob.value
if __name__ == "__main__":
    dA, dE, dF = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); seed = int(sys.argv[4]); ntry = int(sys.argv[5])
    rng = np.random.default_rng(seed); dout = dE*dF*dE
    for t in range(ntry):
        res = minimize(residual, rng.normal(size=2*dout*dA), args=(dA, dE, dF), method="L-BFGS-B", options=dict(maxiter=5000, ftol=1e-30, gtol=1e-14))
        U = iso_from(res.x, dA, dout)
        if res.fun > 1e-16: print(f"try {t}: residual {res.fun:.2e} (no exact TD solution reached)", flush=True); continue
        dT = degradability_defect(U, dA, dE, dF, transpose=True); dD = degradability_defect(U, dA, dE, dF, transpose=False)
        print(f"try {t}: TD residual {res.fun:.1e}; SDP transpose-degradability defect {dT:.2e}; ordinary degradability defect {dD:.3e}" + ("   <<< NOT DEGRADABLE?" if dD > 1e-4 else ""), flush=True)
        if dD > 1e-4: np.save(f"/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/tdeg_{dA}{dE}{dF}_{seed}_{t}.npy", U)

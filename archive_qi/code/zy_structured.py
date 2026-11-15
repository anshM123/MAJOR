"""Structured search for Zhang-Yeung violations near saturating states.

For a state rho0 with gap(rho0)=0 and marginal support projectors Pi_R, and any state X,
  gap((1-e) rho0 + e X) = e*log(1/e) * tr[X W] + O(e),   W := sum_R c_R (1 - Pi_R) (x) 1 = -sum_R c_R Pi_R (x) 1
(the gap functional is balanced, sum_R c_R = 0).  A negative eigenvalue of W at a saturating state would give a
violation; if W >= 0 the O(e) term decides, which we probe by direct minimisation over X at fixed small e.
Saturating states are drawn from stabilizer (mixed) states on m qubits per party."""
import numpy as np, sys, itertools
from entropy_opt import *

PAULI = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
def pauli_mat(v, n):          # v = (x|z) in F_2^{2n}
    M = np.array([[1.0+0j]])
    for i in range(n):
        x, z = v[i], v[n+i]; M = np.kron(M, PAULI[{(0,0):0,(1,0):1,(1,1):2,(0,1):3}[(x, z)]])
    return M
def symp(u, v, n): return (u[:n]@v[n:] + u[n:]@v[:n]) % 2
def rank2(M):
    M = [int("".join(map(str, r)), 2) for r in M]; rk = 0
    for bit in reversed(range(64)):
        piv = [i for i, r in enumerate(M) if (r >> bit) & 1]
        if not piv: continue
        p = M[piv[0]]; M = [r ^ p if ((r >> bit) & 1) else r for i, r in enumerate(M) if i != piv[0]]; rk += 1
    return rk
def random_stab(n, k, rng):
    gens = []
    for _ in range(2000):
        if len(gens) == k: break
        v = rng.integers(0, 2, size=2*n)
        if not v.any(): continue
        if all(symp(v, g, n) == 0 for g in gens) and rank2(gens + [v]) == len(gens) + 1: gens.append(v)
    return gens if len(gens) == k else None
def span(gens):
    out = []
    for bits in itertools.product([0, 1], repeat=len(gens)):
        v = np.zeros_like(gens[0])
        for b, g in zip(bits, gens):
            if b: v = (v + g) % 2
        out.append(v)
    return out
def stab_entropy(gens, n, qubits_R):
    """S(R) in bits for the stabilizer state rho = prod (1+g_i)/2 / 2^{n-k}"""
    G = [v for v in span(gens) if v.any() and all(v[i] == 0 and v[n+i] == 0 for i in range(n) if i not in qubits_R)]
    return len(qubits_R) - (rank2(G) if G else 0)
def stab_rho(gens, n):
    """density matrix; signs chosen +1 (other sign choices are local-Pauli equivalent)"""
    P = np.eye(2**n, dtype=complex)
    for g in gens:
        M = pauli_mat(g, n); P = P@(np.eye(2**n) + M)/2
    # generators commute as symplectic vectors but the Pauli matrices might anticommute-free only up to phase: fix
    P = (P + P.conj().T)/2
    w, V = np.linalg.eigh(P); 
    return P/np.trace(P).real

if __name__ == "__main__":
    m = int(sys.argv[1]); seed = int(sys.argv[2]); ntry = int(sys.argv[3]); n = 4*m
    rng = np.random.default_rng(seed)
    parties = [list(range(p*m, (p+1)*m)) for p in range(4)]
    lhs, rhs = zy_parts(); gapc = comb((1, rhs), (-1, lhs))
    dims = (2**m,)*4; D = 2**n
    gapF = EntropyFunctional(dims, gapc)
    seen = set(); found = 0
    for t in range(ntry):
        k = rng.integers(1, n+1); gens = random_stab(n, k, rng)
        if gens is None: continue
        S = {R: stab_entropy(gens, n, sum((parties[p] for p in R), [])) for R in gapc}
        g = sum(c*S[R] for R, c in gapc.items()); l = sum(c*stab_entropy(gens, n, sum((parties[p] for p in R), [])) for R, c in lhs.items())
        if g != 0 or l == 0: continue
        key = tuple(sorted(S.items()))
        if key in seen: continue
        seen.add(key); found += 1
        rho0 = stab_rho(gens, n)
        if abs(gapF.value(rho0)) > 1e-9: print("  (phase issue, skipped)"); continue
        # support operator
        W = np.zeros((D, D), complex)
        for R, c in gapc.items():
            mR = gapF.marg(rho0, R); w, V = np.linalg.eigh(mR); Pi = (V[:, w > 1e-9])@(V[:, w > 1e-9]).conj().T
            W -= c*gapF.embed(Pi, R)
        wmin = np.linalg.eigvalsh((W+W.conj().T)/2)[0]
        # direct minimisation over X at fixed eps
        best = np.inf
        for eps in (1e-2, 1e-4):
            for rep in range(3):
                r = min(D, 4)
                def obj(x):
                    X, L, tt = rho_from_x(x, D, r); rho = (1-eps)*rho0 + eps*X
                    f, G = gapF.value_grad(rho); G = eps*G
                    Gp = G - np.trace(G@X).real*np.eye(D); gL = 2*(Gp@L)/tt
                    return f, np.concatenate([gL.real.ravel(), gL.imag.ravel()])
                res = minimize(obj, rng.normal(size=2*D*r), jac=True, method="L-BFGS-B", options=dict(maxiter=500))
                best = min(best, res.fun/eps)
        flag = "  <<<<< VIOLATION" if (wmin < -1e-9 or best < -1e-7) else ""
        print(f"saturating stabilizer state #{found} (k={k}, 2I(C:D)={l} bits): min eig W = {wmin:+.4f}; min gap/eps = {best:+.3e}{flag}", flush=True)
    print("distinct saturating entropy vectors examined:", found)

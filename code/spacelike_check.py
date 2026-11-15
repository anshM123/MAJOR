"""Check of Lemma 6.3' (pairwise spacelike separation of geodesics with unlinked endpoints on an acausal circle).
Boundary points p_i = a_i b_i^T with a_i=(cos th_i, sin th_i), b_i=(cos ph_i, sin ph_i), th and ph increasing in (0,pi)
(same cyclic order in both chiral factors).  N_ij = w(a_i,a_j) w(b_i,b_j) > 0.
Geodesics x(t)=(e^t p1+e^-t p2)/sqrt(N12), y(s)=(e^s p3+e^-s p4)/sqrt(N34).
Claim: min_{t,s} tr(x^{-1}y)/2 = [sqrt(N13 N24)+sqrt(N14 N23)]/sqrt(N12 N34) and this is > 1 iff the pairs are unlinked."""
import numpy as np
rng = np.random.default_rng(0)
J = np.array([[0., 1.], [-1., 0.]])
def pt(th, ph): return np.outer([np.cos(th), np.sin(th)], [np.cos(ph), np.sin(ph)])
def N(p, q): return np.linalg.det(p+q)                # det of sum of two rank-one matrices = w(a,c) w(b,d)
def adj(x): return np.array([[x[1, 1], -x[0, 1]], [-x[1, 0], x[0, 0]]])
bad_formula = 0; unlinked_fail = 0; linked_fail = 0; shared = 0
for trial in range(20000):
    th = np.sort(rng.uniform(0.05, np.pi-0.05, 4)); ph = np.sort(rng.uniform(0.05, np.pi-0.05, 4))
    P = [pt(th[i], ph[i]) for i in range(4)]
    for pairing, linked in (((0, 1, 2, 3), False), ((0, 3, 1, 2), False), ((0, 2, 1, 3), True)):
        i, j, k, l = pairing
        Nij, Nkl = N(P[i], P[j]), N(P[k], P[l])
        f = (np.sqrt(N(P[i], P[k])*N(P[j], P[l])) + np.sqrt(N(P[i], P[l])*N(P[j], P[k])))/np.sqrt(Nij*Nkl)
        # brute-force minimum of tr(x^-1 y)/2 over a grid around the predicted minimiser
        t0 = 0.25*np.log(N(P[j], P[l])/N(P[i], P[k])) + 0.25*np.log(N(P[j], P[k])/N(P[i], P[l]))
        s0 = 0.25*np.log(N(P[j], P[l])/N(P[i], P[k])) - 0.25*np.log(N(P[j], P[k])/N(P[i], P[l]))
        x = (np.exp(t0)*P[i] + np.exp(-t0)*P[j])/np.sqrt(Nij); y = (np.exp(s0)*P[k] + np.exp(-s0)*P[l])/np.sqrt(Nkl)
        val = np.trace(adj(x) @ y)/2
        bad_formula += abs(val - f) > 1e-9*max(1, abs(f))
        if linked: linked_fail += not (f <= 1 + 1e-12)
        else: unlinked_fail += not (f > 1)
print("formula mismatches:", bad_formula, "| unlinked pairs with min <= 1 (would be causally related):", unlinked_fail,
      "| linked pairs with min > 1:", linked_fail, " (out of 20000 configurations each)")

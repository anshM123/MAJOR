"""Triage: is p_{n,m}(A,B) >= Tr exp(n log A + m log B) for positive definite A,B ?
p_{n,m} = average of Tr W over all words with n letters A and m letters B."""
import numpy as np, sys
from math import comb
from scipy.optimize import minimize
from scipy.linalg import expm, logm
def word_sum(A, B, n, m):
    d = A.shape[0]; M = [[None]*(m+1) for _ in range(n+1)]
    for i in range(n+1):
        for j in range(m+1):
            if i == 0 and j == 0: M[i][j] = np.eye(d, dtype=complex); continue
            X = np.zeros((d, d), complex)
            if i > 0: X += M[i-1][j]@A
            if j > 0: X += M[i][j-1]@B
            M[i][j] = X
    return np.trace(M[n][m]).real
def ratio(x, d, n, m):
    a = np.exp(x[:d]); A = np.diag(a).astype(complex)
    L = (x[d:d+d*d] + 1j*x[d+d*d:]).reshape(d, d); B = L@L.conj().T + 1e-12*np.eye(d)
    w, V = np.linalg.eigh(B); logB = (V*np.log(w))@V.conj().T
    ex = np.trace(expm(n*np.diag(np.log(a)) + m*logB)).real
    return word_sum(A, B, n, m)/comb(n+m, n)/ex
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for d in (2, 3):
        for n in range(1, 7):
            for m in range(n, 7):
                if n+m < 3: continue
                best = np.inf
                for t in range(6):
                    res = minimize(ratio, rng.normal(size=d+2*d*d)*rng.choice([0.5, 1.5]), args=(d, n, m), method="L-BFGS-B", options=dict(maxiter=300))
                    best = min(best, res.fun)
                print(f"d={d} n={n} m={m}: min ratio p/Trexp = {best:.6f}" + ("   <<< BELOW 1" if best < 1-1e-7 else ""), flush=True)

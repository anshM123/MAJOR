"""Independent high-precision (mpmath) evaluation of I(A;C|B) and -log F(rho_ABC, Petz(rho_AB)).
No code shared with petz_search.py. Natural logarithms. F(r,s) = (tr sqrt(sqrt r s sqrt r))^2."""
import mpmath as mp

def hfun(A, f):
    A = (A + A.H)/2
    w, V = mp.eigh(A)
    return V*mp.diag([f(w[i]) for i in range(len(w))])*V.H

def entropy(A):
    w = mp.eigh((A + A.H)/2, eigvals_only=True)
    return -sum(x*mp.log(x) for x in w if x > 0)

def kron(A, B):
    C = mp.zeros(A.rows*B.rows, A.cols*B.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] != 0:
                for k in range(B.rows):
                    for l in range(B.cols):
                        C[i*B.rows+k, j*B.cols+l] = A[i, j]*B[k, l]
    return C

def ptrace(rho, dims, keep):
    """partial trace keeping the subsystems in `keep` (list of indices), dims = (dA,dB,dC)"""
    n = len(dims); D = rho.rows
    def unravel(i):
        idx = []
        for d in reversed(dims): idx.append(i % d); i //= d
        return idx[::-1]
    dk = [dims[k] for k in keep]; Dk = 1
    for d in dk: Dk *= d
    def rav(idx):
        r = 0
        for k in keep: r = r*dims[k] + idx[k]
        return r
    out = mp.zeros(Dk, Dk)
    idxs = [unravel(i) for i in range(D)]
    for i in range(D):
        for j in range(D):
            if all(idxs[i][t] == idxs[j][t] for t in range(n) if t not in keep):
                out[rav(idxs[i]), rav(idxs[j])] += rho[i, j]
    return out

def cmi_and_petz(rho, dims):
    dA, dB, dC = dims
    rAB = ptrace(rho, dims, [0, 1]); rBC = ptrace(rho, dims, [1, 2]); rB = ptrace(rho, dims, [1])
    I = entropy(rAB) + entropy(rBC) - entropy(rB) - entropy(rho)
    IA = mp.eye(dA); IC = mp.eye(dC)
    K = kron(IA, hfun(rBC, mp.sqrt)*kron(hfun(rB, lambda x: 1/mp.sqrt(x)), IC))
    rec = K*kron(rAB, IC)*K.H
    sr = hfun(rho, lambda x: mp.sqrt(max(x, 0)))
    M = sr*rec*sr
    w = mp.eigh((M + M.H)/2, eigvals_only=True)
    rootF = sum(mp.sqrt(x) for x in w if x > 0)
    tr_rec = sum(rec[i, i] for i in range(rec.rows))
    return I, -2*mp.log(rootF), mp.re(tr_rec)

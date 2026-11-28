"""Tripartite information of three spin blocks in the ground state of critical spin chains (exact diagonalisation).
Critical transverse-field Ising (c=1/2, lightest primary sigma: Z2-odd scalar, C_sss = 0) versus the XX chain
(c=1 Dirac fermion, I3 = 0 in the continuum for fermion blocks; spin blocks differ only by JW strings).
Expectation from Agon-Bueno-Casini 2109.09179 eq.(3)/(34) (C_OOO = 0): I3 > 0 at long distances for Ising.
If I3 > 0, MMI (balanced, hence blind to Bell pairs) is violated, so G = S - sum t_p b_p is NOT in the RT cone."""
import itertools, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
def op_chain(L, terms):
    H = sp.csr_matrix((2**L, 2**L))
    for coef, ops in terms:
        mats = [sp.identity(2, format='csr')]*L
        for site, m in ops: mats[site] = sp.csr_matrix(m)
        M = mats[0]
        for m in mats[1:]: M = sp.kron(M, m, format='csr')
        H = H + coef*M
    return H
X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1.])
def ising(L): return op_chain(L, [(-1, [(i, X), ((i+1) % L, X)]) for i in range(L)] + [(-1, [(i, Z)]) for i in range(L)])
def xx(L):    return op_chain(L, [(-0.5, [(i, X), ((i+1) % L, X)]) for i in range(L)] + [(-0.5, [(i, Y), ((i+1) % L, Y)]) for i in range(L)])
def entropy(psi, L, sites):
    sites = sorted(sites); rest = [s for s in range(L) if s not in sites]
    T = psi.reshape([2]*L).transpose(sites + rest).reshape(2**len(sites), -1)
    rho = T @ T.conj().T; w = np.linalg.eigvalsh(rho); w = w[w > 1e-14]
    return float(-(w*np.log(w)).sum())
def i3(psi, L, A, B, C):
    S = lambda *Rs: entropy(psi, L, [s for R in Rs for s in R])
    return S(A)+S(B)+S(C)-S(A, B)-S(A, C)-S(B, C)+S(A, B, C)
if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    for name, Hf in (("Ising", ising), ("XX", xx)):
        H = Hf(L)
        e, v = sla.eigsh(H.astype(complex) if name == "XX" else H, k=1, which='SA'); psi = v[:, 0]
        out = []
        for ell in (1, 2):
            for d in range(1, (L//3) - ell + 1):
                A = list(range(0, ell)); B = list(range(ell+d, 2*ell+d)); C = list(range(2*ell+2*d, 3*ell+2*d))
                if C[-1] >= L - d: continue
                out.append((ell, d, i3(psi, L, A, B, C)))
        print(f"{name} chain L={L} (periodic), E0={e[0]:.6f}:")
        for ell, d, v3 in out: print(f"   blocks of {ell} sites, gaps {d}: I3 = {v3:+.3e}")

def fermion_lattice_i3(ell, d):
    """Free hopping fermion (half filling, infinite chain): entropies of FERMIONIC blocks from the correlation matrix."""
    def S(sites):
        s = np.array(sites); D = s[:, None] - s[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            C = np.where(D == 0, 0.5, np.sin(np.pi*D/2)/(np.pi*D))
        w = np.linalg.eigvalsh(C); w = w[(w > 1e-14) & (w < 1-1e-14)]
        return float(-(w*np.log(w) + (1-w)*np.log(1-w)).sum())
    A = list(range(0, ell)); B = list(range(ell+d, 2*ell+d)); C = list(range(2*ell+2*d, 3*ell+2*d))
    return S(A)+S(B)+S(C)-S(A+B)-S(A+C)-S(B+C)+S(A+B+C)

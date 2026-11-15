"""Second-order test of Zhang-Yeung around reference states where every term vanishes.

Symmetric form:  I(C:D) <= I(A:B) + I(C:D|B) + 2 I(C:D|A) + I(A:C|D) + I(A:D|C)   (*)
(equivalent to 2I(C:D) <= I(A:B)+I(A:CD)+3I(C:D|A)+I(C:D|B); every term on the right is >= 0 by SSA).
At rho0 = rho_A (x) rho_{B1 C} (x) rho_{B2 D}, B=B1B2, all six quantities vanish, so for rho0 + eps*delta
each is (eps^2/2) * (a BKM quadratic form in delta) + O(eps^3):
   I(X:Y|Z)  ->  q = g_{XYZ} + g_Z - g_{XZ} - g_{YZ},   g_R(delta) = BKM metric of delta_R at rho0_R.
(*) requires QL <= QR as quadratic forms; a generalised eigenvalue > 1 gives a violating state."""
import numpy as np, sys, itertools
from scipy.sparse import coo_matrix

def herm_basis(D):
    B = np.zeros((D*D, D*D), complex); k = 0
    for i in range(D):
        B[i*D+i, k] = 1; k += 1
    for i in range(D):
        for j in range(i+1, D):
            B[i*D+j, k] = B[j*D+i, k] = 1/np.sqrt(2); k += 1
            B[i*D+j, k] = -1j/np.sqrt(2); B[j*D+i, k] = 1j/np.sqrt(2); k += 1
    return B

class SecondOrder:
    def __init__(self, qdims, parties):
        """qdims: dims of elementary subsystems; parties: dict name -> list of elementary indices"""
        self.q = tuple(qdims); self.n = len(qdims); self.D = int(np.prod(qdims)); self.parties = parties
        self.B = herm_basis(self.D); self._M = {}
    def idx(self, R): return sorted(sum((self.parties[p] for p in R), []))
    def ptrace_map(self, R):
        key = tuple(self.idx(R))
        if key in self._M: return self._M[key]
        keep = list(key); D = self.D; dk = [self.q[i] for i in keep]; dR = int(np.prod(dk))
        multi = np.array(list(itertools.product(*[range(d) for d in self.q])))
        rk = np.zeros(D, int)
        for i in keep: rk = rk*self.q[i] + multi[:, i]
        drop = [i for i in range(self.n) if i not in keep]
        dkey = np.zeros(D, int)
        for i in drop: dkey = dkey*self.q[i] + multi[:, i]
        I, J = np.meshgrid(np.arange(D), np.arange(D), indexing="ij")
        mask = dkey[I] == dkey[J]
        rows = (rk[I]*dR + rk[J])[mask]; cols = (I*D + J)[mask]
        M = coo_matrix((np.ones(rows.size), (rows, cols)), shape=(dR*dR, D*D)).tocsr()
        self._M[key] = (M, dR); return M, dR
    def marg(self, rho, R):
        M, dR = self.ptrace_map(R); return (M@rho.reshape(-1)).reshape(dR, dR)
    def g(self, rho0, R):
        """real Gram matrix of the BKM form of delta_R at rho0_R in the Hermitian basis"""
        M, dR = self.ptrace_map(R); m = self.marg(rho0, R)
        w, V = np.linalg.eigh((m+m.conj().T)/2); x = w[:, None]; y = w[None, :]
        with np.errstate(divide="ignore", invalid="ignore"): K = (np.log(x)-np.log(y))/(x-y)
        eq = np.abs(x-y) < 1e-12*(x+y); K[eq] = (2/(x+y))[eq]
        Y = (M@self.B).reshape(dR, dR, -1)
        Y = np.einsum("ia,ijk,jb->abk", V.conj(), Y, V)*np.sqrt(K)[:, :, None]
        Y = Y.reshape(dR*dR, -1)
        return (Y.conj().T@Y).real
    def cmi_form(self, rho0, X, Y, Z=()):
        X, Y, Z = list(X), list(Y), list(Z)
        q = self.g(rho0, X+Y+Z) - self.g(rho0, X+Z) - self.g(rho0, Y+Z)
        if Z: q = q + self.g(rho0, Z)
        return q

def zy_forms(so, rho0):
    QL = so.cmi_form(rho0, "C", "D")
    QR = (so.cmi_form(rho0, "A", "B") + so.cmi_form(rho0, "C", "D", "B") + 2*so.cmi_form(rho0, "C", "D", "A")
          + so.cmi_form(rho0, "A", "C", "D") + so.cmi_form(rho0, "A", "D", "C"))
    return QL, QR

def max_ratio(QL, QR, tol=1e-9, vec=False):
    w, U = np.linalg.eigh((QR+QR.T)/2); keep = w > tol*w.max()
    null = U[:, ~keep]
    leak = np.linalg.eigvalsh(null.T@QL@null).max() if null.shape[1] else 0.0
    Ur = U[:, keep]/np.sqrt(w[keep]); Mx = Ur.T@QL@Ur; ev, evec = np.linalg.eigh((Mx+Mx.T)/2)
    if vec: return ev[-1], leak, Ur@evec[:, -1]
    return ev[-1], leak

def rand_state(d, rng, spread):
    G = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)); U, _ = np.linalg.qr(G)
    lam = np.exp(rng.normal(size=d)*spread); lam /= lam.sum(); return (U*lam)@U.conj().T

if __name__ == "__main__":
    seed = int(sys.argv[1]); nsamp = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    # elementary systems: A, B1, B2, C, D  (all qubits)
    so = SecondOrder((2, 2, 2, 2, 2), {"A": [0], "B": [1, 2], "C": [3], "D": [4]})
    def ref(rA, rB1C, rB2D):
        # order A,B1,B2,C,D ; build rho_A (x) rho_{B1C} (x) rho_{B2D} then permute (A,B1,C,B2,D)->(A,B1,B2,C,D)
        T = np.kron(np.kron(rA, rB1C), rB2D).reshape((2,)*10)
        T = T.transpose(0, 1, 3, 2, 4, 5, 6, 8, 7, 9); return T.reshape(32, 32)
    best = 0
    for t in range(nsamp):
        sp = rng.choice([0.3, 1.0, 2.5])
        rho0 = ref(rand_state(2, rng, sp), rand_state(4, rng, sp), rand_state(4, rng, sp))
        QL, QR = zy_forms(so, rho0); lam, leak = max_ratio(QL, QR)
        best = max(best, lam)
        print(f"sample {t} (spread {sp}): max QL/QR = {lam:.6f}   null-space leak = {leak:.2e}", flush=True)
    print("best", best)

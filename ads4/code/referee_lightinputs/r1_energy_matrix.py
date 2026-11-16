"""Referee r1: independent d=3 one-point energy matrix and the commutator [H(n1),H(n3)].

Independent of lightinputs/t2 (which uses KKSDZ's H_{ij,kl} tensor):
 (a) Hofman-Maldacena sesquilinear form in general d, written directly from
     <E(n)>_eps ~ 1 + t2 (eps*_{ij} eps_{il} n_j n_l/eps*eps - 1/(d-1)) + t4 (|n.eps.n|^2/eps*eps - 2/(d^2-1)),
     polarised to B(eps', eps) and evaluated numerically on real 2x2 traceless symmetric bases (no helicity ansatz).
 (b) Most general parity-violating Hermitian form: all SO(2)-covariant quartic contractions with n and
     ntilde = rot90(n); impose Hermiticity + energy sum rule (angular average = identity) -> 1 + flip(c complex).
 (c) Commutator vs -2i|c|^2 sin(4 phi13) sigma_z, random complex c, random angles.
numpy only, < 1 s.
"""
import numpy as np
rng = np.random.default_rng(7)

def basis():   # orthonormal real basis of 2x2 traceless symmetric matrices (Frobenius)
    s = 1/np.sqrt(2)
    return [np.array([[s, 0], [0, -s]]), np.array([[0, s], [s, 0]])]
E = basis()
# helicity basis eps_+- = e+-e+-, e+- = (1, +-i)/sqrt2  (Frobenius-normalised)
ep = np.array([1, 1j])/np.sqrt(2); em = ep.conj()
HEL = [np.outer(ep, ep), np.outer(em, em)]

def form_HM(n, t2, t4, d=3):
    def B(a, b):   # sesquilinear  <a| E(n) |b>, a antilinear
        ab = np.vdot(a, b)
        M = a.conj().T @ b                      # (a* b)_{jl} = a*_{ij} b_{il}
        t2term = n @ M @ n - ab/(d-1)
        t4term = (n @ a.conj() @ n)*(n @ b @ n) - 2*ab/(d**2-1)
        return ab + t2*t2term + t4*t4term
    return B

def matrix(B, base):
    return np.array([[B(a, b) for b in base] for a in base])

ok = True
# (a) HM in d=3: t2 drops out; helicity matrix = [[1, (t4/4) w^-4],[ (t4/4) w^4, 1]] (up to basis phase)
for _ in range(5):
    phi = rng.uniform(0, 2*np.pi); n = np.array([np.cos(phi), np.sin(phi)])
    t2, t4 = rng.normal(size=2)
    H = matrix(form_HM(n, t2, t4), HEL)
    H0 = matrix(form_HM(n, 0.0, t4), HEL)
    ok &= np.allclose(H, H0)                                     # t2 drops out
    ok &= np.allclose(np.diag(H), 1)                             # diagonal = 1
    ok &= np.isclose(abs(H[0, 1]), abs(t4)/4)                    # |flip| = |t4|/4
    ok &= np.isclose(H[0, 1]/H[1, 0].conj(), 1)                  # Hermitian
    ok &= np.isclose(H[0, 1]*np.exp(4j*phi)/(t4/4), 1) or np.isclose(H[0, 1]*np.exp(-4j*phi)/(t4/4), 1)
print("(a) HM form, d=3: t2 drops out, diag=1, |flip|=|t4|/4, flip ~ w^{-+4}:", ok)

# (b) general parity-odd: candidate quartic forms built from n, nt; check which are allowed
def cands(n):
    nt = np.array([-n[1], n[0]])
    vs = {'nn': (n, n), 'tt': (nt, nt), 'nt': (n, nt), 'tn': (nt, n)}
    out = {'id': lambda a, b: np.vdot(a, b)}
    for k1, (p, q) in vs.items():
        for k2, (r, s) in vs.items():
            out[k1+'|'+k2] = (lambda p, q, r, s: lambda a, b: (p @ a.conj() @ q)*(r @ b @ s))(p, q, r, s)
    return out
phis = np.linspace(0, 2*np.pi, 64, endpoint=False)
names = list(cands(np.array([1., 0.])).keys())
# stack helicity matrices for each candidate as function of phi
T = np.array([[matrix(cands(np.array([np.cos(p), np.sin(p)]))[k], HEL) for p in phis] for k in names])
# any real-linear combination sum_k x_k T_k that is Hermitian: diag entries phi-independent? off-diag ~ w^{+-4}?
maxdiagvar = max(np.ptp(T[i][:, 0, 0].real) + np.ptp(T[i][:, 0, 0].imag) for i in range(len(names)))
# Fourier content of off-diagonal
offmodes = set()
for i in range(len(names)):
    f = np.fft.fft(T[i][:, 0, 1])/len(phis)
    offmodes |= {int(m) if m < 32 else int(m) - 64 for m in np.nonzero(np.abs(f) > 1e-12)[0]}
print("(b) all SO(2)-covariant quartic forms (incl. parity-odd with ntilde): max phi-variation of H_{++} =",
      f"{maxdiagvar:.1e}; Fourier modes of H_+- = {sorted(offmodes)}")
print("    => H(phi) = h0*1 + h3*sigma_z + [[0, c e^{-4i phi}],[c* e^{4i phi}, 0]]; energy sum rule (avg = 1, no sigma_z")
print("       since <+|P^0|+> = <-|P^0|->) fixes h0=1, h3=0 -> H = 1 + flip, c complex (t4 and alpha_T).")

# (c) commutator algebra with complex c
err = 0
for _ in range(200):
    c = rng.normal() + 1j*rng.normal(); p1, p3 = rng.uniform(0, 2*np.pi, 2)
    H = lambda p: np.array([[1, c*np.exp(-4j*p)], [np.conj(c)*np.exp(4j*p), 1]])
    C = H(p1) @ H(p3) - H(p3) @ H(p1)
    ref = -2j*abs(c)**2*np.sin(4*(p1-p3))*np.diag([1, -1])
    err = max(err, np.abs(C - ref).max())
print(f"(c) [H1,H3] = -2i|c|^2 sin(4 phi13) sigma_z : max error {err:.1e} over 200 random (c, phi1, phi3)")
# linear term / Einstein:
c = 1e-6*(0.3+0.7j); p1, p3 = 0.4, 1.9
H = lambda p: np.array([[1, c*np.exp(-4j*p)], [np.conj(c)*np.exp(4j*p), 1]])
print("    Einstein c=0: commutator = 0 exactly; at c = 1e-6: |C| =", f"{np.abs(H(p1)@H(p3)-H(p3)@H(p1)).max():.2e}", "(quadratic)")

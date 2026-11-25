"""
np3_tstate.py -- the regulated null-plane identity in the stress-tensor state psi_T = T(F)|0> of the d = 3 free scalar
(a two-particle state; F^(P) = exp(-|P - P0|^2/(2 sig^2)), P0 = (1,0,0): a T wavepacket near rest, helicity +2).

    psi(p1,p2) = F^(p1+p2) * Pi(p1,p2),   Pi = -(3/4)(z.p1)(z.p2) + (1/8)((z.p1)^2 + (z.p2)^2),  z = (0,1,i)/2
    (improved stress tensor, xi = 1/8; same polynomial as code/kernels/fock_free.py).
    LF wavefunction chi(K1,q1,K2,q2) = psi/sqrt(4 K1 K2).

Light-front subtlety handled here: psi_T populates K -> 0 (particles moving along the null generators, i.e. the forward
direction +x).  Slice masses behave like K^{-1/2} and the transverse spread like K^{-1/2}; with f ~ y^2 every term of the
identity has an integrable K^{-1/2} density.  We use K = s^2 (Gauss-Legendre in s), per-K transverse grids.

Static two-body pieces: Schmidt decomposition of each (K1,K2) block.  I_eta: one-body reduced density matrix per K1
(exact, since [A, A_t] is a one-body operator for a free Hamiltonian).
Output: cache/np3_tstate.txt
"""
import numpy as np, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
from lfcore import Slice, gauss_legendre, fm, profile, I_eta_from_G

out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

E0, SIG = 1.0, 0.25
HMAX = E0 + 7*SIG

def psiT(K1, q1, K2, q2):
    h1, h2 = q1**2/(4*K1), q2**2/(4*K2)
    p1 = (h1 + K1, h1 - K1, q1); p2 = (h2 + K2, h2 - K2, q2)
    P0, Px, Py = p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]
    Fh = np.exp(-((P0 - E0)**2 + Px**2 + Py**2)/(2*SIG**2))
    z1 = (p1[1] + 1j*p1[2])/2; z2 = (p2[1] + 1j*p2[2])/2
    Pi = -0.75*z1*z2 + 0.125*(z1**2 + z2**2)
    return Fh*Pi/np.sqrt(4*K1*K2)

class KNode(Slice):
    def __init__(self, K, w, L, dy, Nmax=2**16):
        if 2*L/dy > Nmax:
            dy = 2*L/Nmax
        super().__init__(K, w, L, dy)
        qcut = 2*np.sqrt(K*HMAX)*1.05 + 1e-12
        self.S = np.where(np.abs(self.q) <= qcut)[0]
        self.ph = np.exp(1j*self.q[self.S]*self.y[0])
    def to_y(self, cS):
        """coefficients on the support (l2-normalised as slice norm) -> y-vector with sum |x|^2 dy = sum |c|^2."""
        full = np.zeros(self.N, complex); full[self.S] = cS*self.ph
        return np.fft.ifft(full)*np.sqrt(self.N/self.dy)

def build(ns=40, dy=0.04, cL=45.0):
    s, ws = gauss_legendre(0.0, np.sqrt(HMAX), ns)
    Ks, wK = s**2, 2*s*ws
    nodes = [KNode(K, w, L=8.0 + cL/np.sqrt(K), dy=dy) for K, w in zip(Ks, wK)]
    return nodes

def block(n1, n2):
    Q1, Q2 = np.meshgrid(n1.q[n1.S], n2.q[n2.S], indexing='ij')
    # l2 norm of the coefficients = continuum norm int dq1 dq2 |chi|^2 (up to 2pi's), uniformly across K-nodes
    return psiT(n1.K, Q1, n2.K, Q2)/np.sqrt(n1.N*n1.dy*n2.N*n2.dy)

def tail_check(nodes):
    """fraction of |chi|^2 in |y1| > 0.75 L(K1) for the heaviest K2 partner (y-grid adequacy)."""
    worst = 0.0
    for n1 in nodes:
        mass = [np.sum(np.abs(block(n1, n2))**2) for n2 in nodes]
        n2 = nodes[int(np.argmax(mass))]
        B = block(n1, n2)
        for col in range(0, B.shape[1], max(1, B.shape[1]//6)):
            x = n1.to_y(B[:, col])
            if np.sum(np.abs(x)**2) > 0:
                worst = max(worst, np.sum(np.abs(x[np.abs(n1.y) > 0.75*n1.y[-1]])**2)/np.sum(np.abs(x)**2))
    return worst

def profiles():
    pr = []
    for m in (0, 1, 2, 4):
        for part in ('re', 'im'):
            if not (m == 0 and part == 'im'):
                pr.append((f'{part} f_{m}', (lambda mm, pp: (lambda y: profile(y, mm, pp)))(m, part)))
    return pr

if __name__ == '__main__':
    t0 = time.time()
    nodes = build(ns=28)
    P(f"psi_T: E0={E0}, sig={SIG}; {len(nodes)} K-nodes (K = s^2, Gauss-Legendre in s), K in [{nodes[0].K:.2e}, {nodes[-1].K:.3f}],"
      f" N_y in [{min(n.N for n in nodes)}, {max(n.N for n in nodes)}], support sizes [{min(len(n.S) for n in nodes)}, {max(len(n.S) for n in nodes)}]")
    P(f"  y-grid tail check (worst fraction of |chi|^2 beyond 0.75 L): {tail_check(nodes):.1e}")
    prs = profiles()
    nrm = 0.0
    # ---------------- one-body reduced density matrices rho_1(K1) on the q-support
    rho = []
    for n1 in nodes:
        R = np.zeros((len(n1.S), len(n1.S)), complex)
        for n2 in nodes:
            B = block(n1, n2)
            R += 2*n2.w*(B @ B.conj().T)
        rho.append(R)
        nrm += n1.w*np.real(np.trace(R))/2
    P(f"  <psi|psi> (discrete) = {nrm:.6e};  Tr rho_1 / <psi|psi> = 2 (two particles)")
    eig = []
    lam_max = max(np.max(np.linalg.eigvalsh(R)) for R in rho)
    for n1, R in zip(nodes, rho):
        lam, V = np.linalg.eigh(R)
        keep = lam > 1e-13*lam_max
        eig.append((lam[keep]/nrm, V[:, keep]))
    P(f"  rho_1 ranks kept per K: min {min(len(e[0]) for e in eig)}, max {max(len(e[0]) for e in eig)}")
    def Xy(n, V):
        return [n.to_y(V[:, j]) for j in range(V.shape[1])]
    PK = sum(n.w*np.sum(l*np.array([n.K*np.real(np.vdot(x, x))*n.dy for x in Xy(n, X)])) for n, (l, X) in zip(nodes, eig))
    Ay2 = sum(n.w*np.sum(l*np.array([n.K*np.sum(n.y**2*np.abs(x)**2)*n.dy for x in Xy(n, X)])) for n, (l, X) in zip(nodes, eig))
    P(f"  <A(1)> = <P_K> = {PK:.8f}  (= <(P^0 - P^x)/2> = E0/2 for the rest-frame packet up to O(sig^2)),  <A(y^2)> = {Ay2:.8f}")
    # K^{-1/2} law near K = 0: density of <A(y^2)> in K
    dens = [(n.K, n.w*np.sum(l*np.array([n.K*np.sum(n.y**2*np.abs(x)**2)*n.dy for x in Xy(n, X)]))/n.w) for n, (l, X) in zip(nodes, eig)]
    P("  small-K density of <A(y^2)> (per unit K) times sqrt(K):  " +
      "  ".join(f"K={k:.1e}:{d*np.sqrt(k):.4f}" for k, d in dens[:6]) + "   -> finite limit: integrable K^{-1/2} density")

    # ---------------- one-body: W and I_eta
    etas = (0.1, 0.05, 0.025)
    oneb = {}
    for name, Ff in prs:
        W = 0.0
        for n, (l, X) in zip(nodes, eig):
            F, Fp = Ff(n.y)
            W += n.w*np.sum(l*np.array([0.5*n.K*np.sum(Fp**2*np.abs(x)**2)*n.dy for x in Xy(n, X)]))
        def G(ts, Ff=Ff):
            tot = np.zeros(len(ts))
            for n, (l, X) in zip(nodes, eig):
                F, _ = Ff(n.y); a = n.K*F
                for lj, x in zip(l, Xy(n, X)):
                    Uc = n.evolve(x, ts); Uac = n.evolve(a*x, ts)
                    tot += n.w*lj*2*np.imag(np.sum(np.conj(Uac)*a[None, :]*Uc, axis=1))*n.dy
            return tot
        Is = {e: I_eta_from_G(G, e, n=(8 if name.endswith(('f_0', 'f_1')) else 30)) for e in etas}
        oneb[name] = (W, Is)
        P(f"  {name}:  W = (1/2)<A(F'^2)> = {W:.8f};  I_eta (eta={etas}) = " + "  ".join(f"{Is[e]:.8f}" for e in etas) +
          f";  I(0.025)-W = {Is[0.025] - W:.1e}   [{time.time() - t0:.0f}s]")

    # ---------------- two-body static pieces: support-basis operator matrices + Schmidt blocks
    # E = columns e_k(y) = to_y(delta_k), k in support S (h-eigenvectors).  For a one-body multiplication a(y):
    #   Ad = E^H a E,  AA = (aE)^H (aE),  AhA = (aE)^H h (aE)  (exact grid sums via FFT), and E^H h = diag(h_S) E^H.
    # <chi, (X1 (x) X2) chi> = Tr(Phi^H X1 Phi X2^T) = Tr(S P1 S P2),  Phi = U S W,  P1 = U^H X1 U,  P2 = W X2^T W^H.
    blocks = {}
    for i1, n1 in enumerate(nodes):
        for i2, n2 in enumerate(nodes):
            if i2 < i1:
                continue
            B = block(n1, n2)/np.sqrt(nrm)
            U, S, Wm = np.linalg.svd(B, full_matrices=False)
            r = int(np.sum(S > 1e-8*S[0])) if S[0] > 0 else 0
            if r:
                blocks[(i1, i2)] = (U[:, :r], S[:r], Wm[:r, :])
    P(f"  Schmidt blocks stored: {len(blocks)}, max rank {max(b[1].size for b in blocks.values())}   [{time.time() - t0:.0f}s]")
    def make_E(n):     # N x |S| matrix of support plane waves in y-space
        full = np.zeros((n.N, len(n.S)), complex); full[n.S, np.arange(len(n.S))] = n.ph
        return np.fft.ifft(full, axis=0)*np.sqrt(n.N/n.dy)
    def node_mats(n, E, a, wfun=None):
        aE = a[:, None]*E
        Ad = E.conj().T @ aE*n.dy
        AA = aE.conj().T @ aE*n.dy
        F = np.fft.fft(aE, axis=0)
        AhA = F.conj().T @ (n.h[:, None]*F)*n.dy/n.N
        Wd = None if wfun is None else E.conj().T @ (wfun[:, None]*E)*n.dy
        return Ad, AA, AhA, np.diag(n.h[n.S]).astype(complex), Wd
    def T(blk, X1, X2):
        U, S, Wm = blk
        P1 = U.conj().T @ X1 @ U; P2 = Wm @ X2.T @ Wm.conj().T
        return np.trace((S[:, None]*P1)*1.0 @ (S[:, None]*P2))
    def two_body(mats, hermitian):
        res = np.zeros(4)
        for (i1, i2), blk in blocks.items():
            ww = nodes[i1].w*nodes[i2].w*(1 if i1 == i2 else 2)
            A1, AA1, AhA1, H1, W1 = mats[i1]; A2, AA2, AhA2, H2, W2 = mats[i2]
            I1, I2 = np.eye(len(H1)), np.eye(len(H2))
            Hn = (T(blk, AhA1, I2) + T(blk, A1.conj().T @ H1, A2) + T(blk, H1 @ A1, A2.conj().T) + T(blk, H1, AA2)
                  + T(blk, AA1, H2) + T(blk, A1.conj().T, H2 @ A2) + T(blk, A1, A2.conj().T @ H2) + T(blk, I1, AhA2))
            res[0] += ww*2*np.real(Hn)
            if hermitian:
                Rc = (T(blk, AA1 @ H1, I2) + T(blk, AA1, H2) + 2*T(blk, A1 @ H1, A2) + 2*T(blk, A1, A2 @ H2)
                      + T(blk, H1, AA2) + T(blk, I1, AA2 @ H2))
                res[1] += ww*2*np.real(Rc)
                res[2] += ww*0.5*np.real(T(blk, W1, I2) + T(blk, I1, W2))
                res[3] += ww*np.real(T(blk, I1, I2))
        return res
    acc = {}
    for name, Ff in prs:
        mats = []
        for n in nodes:
            F, Fp = Ff(n.y)
            mats.append(node_mats(n, make_E(n), n.K*F, n.K*Fp**2))
        acc[name] = two_body(mats, True)
        P(f"    two-body {name}: done [{time.time() - t0:.0f}s]  (block-norm check {acc[name][3]:.10f})")
    mlist = (0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32)
    accm = {}
    for m in mlist:
        mats = [node_mats(n, make_E(n), n.K*fm(n.y, m)[0]) for n in nodes]
        accm[m] = two_body(mats, False)[0]/2
    P(f"    m-scan done [{time.time() - t0:.0f}s]")
    P("\n  profile     2||H^1/2 A psi||^2   2Re<Apsi,AHpsi>    difference        W (two-body)     W (rho_1)     lim I_eta")
    for name, _ in prs:
        Hn, Rc, Wb, _ = acc[name]
        P(f"  {name:8s}  {Hn:16.8f}  {Rc:16.8f}  {Hn - Rc:16.8f}  {Wb:14.8f}  {oneb[name][0]:14.8f}  {oneb[name][1][0.025]:14.8f}")
    P("\n  m-scaling (complex A_m, psi_T fixed):   m   ||H^1/2 A_m psi_T||^2     m^2 <A(1)>      ratio")
    for m in mlist:
        P(f"      {m:3d}   {accm[m]:16.8f}   {m*m*PK:14.8f}   {accm[m]/max(m*m*PK, 1e-300) if m else float('nan'):8.5f}")
    ms = np.array(mlist, float); Hs = np.array([accm[m] for m in mlist])
    cf = np.polyfit(ms[2:], Hs[2:], 2)
    P(f"  quadratic fit (m>=2): {cf[0]:.6f} m^2 + {cf[1]:.6f} m + {cf[2]:.6f};  <A(1)> = {PK:.6f}")
    P(f"\n[time {time.time() - t0:.1f}s]")
    open(os.path.join(os.path.dirname(__file__), 'cache', 'np3_tstate.txt'), 'w').write('\n'.join(out) + '\n')

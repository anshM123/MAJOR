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

def build(ns=40, dy=0.04, cL=60.0):
    s, ws = gauss_legendre(0.0, np.sqrt(HMAX), ns)
    Ks, wK = s**2, 2*s*ws
    nodes = [KNode(K, w, L=8.0 + cL/np.sqrt(K), dy=dy) for K, w in zip(Ks, wK)]
    return nodes

def block(n1, n2):
    Q1, Q2 = np.meshgrid(n1.q[n1.S], n2.q[n2.S], indexing='ij')
    return psiT(n1.K, Q1, n2.K, Q2)*np.sqrt(n1.dy*n2.dy/(n1.N*n2.N))   # l2 norm = slice norm

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
    nodes = build()
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
        keep = lam > 1e-15*lam_max
        eig.append((lam[keep]/nrm, [n1.to_y(V[:, j]) for j in np.where(keep)[0]]))
    P(f"  rho_1 ranks kept per K: min {min(len(e[0]) for e in eig)}, max {max(len(e[0]) for e in eig)}")
    PK = sum(n.w*np.sum(l*np.array([n.K*np.real(np.vdot(x, x))*n.dy for x in X])) for n, (l, X) in zip(nodes, eig))
    Ay2 = sum(n.w*np.sum(l*np.array([n.K*np.sum(n.y**2*np.abs(x)**2)*n.dy for x in X])) for n, (l, X) in zip(nodes, eig))
    P(f"  <A(1)> = <P_K> = {PK:.8f}  (= <(P^0 - P^x)/2> = E0/2 for the rest-frame packet up to O(sig^2)),  <A(y^2)> = {Ay2:.8f}")
    # K^{-1/2} law near K = 0: density of <A(y^2)> in K
    dens = [(n.K, n.w*np.sum(l*np.array([n.K*np.sum(n.y**2*np.abs(x)**2)*n.dy for x in X]))/n.w) for n, (l, X) in zip(nodes, eig)]
    P("  small-K density of <A(y^2)> (per unit K) times sqrt(K):  " +
      "  ".join(f"K={k:.1e}:{d*np.sqrt(k):.4f}" for k, d in dens[:6]) + "   -> finite limit: integrable K^{-1/2} density")

    # ---------------- one-body: W and I_eta
    etas = (0.2, 0.1, 0.05, 0.025)
    oneb = {}
    for name, Ff in prs:
        W = 0.0
        for n, (l, X) in zip(nodes, eig):
            F, Fp = Ff(n.y)
            W += n.w*np.sum(l*np.array([0.5*n.K*np.sum(Fp**2*np.abs(x)**2)*n.dy for x in X]))
        def G(ts, Ff=Ff):
            tot = np.zeros(len(ts))
            for n, (l, X) in zip(nodes, eig):
                F, _ = Ff(n.y); a = n.K*F
                for lj, x in zip(l, X):
                    Uc = n.evolve(x, ts); Uac = n.evolve(a*x, ts)
                    tot += n.w*lj*2*np.imag(np.sum(np.conj(Uac)*a[None, :]*Uc, axis=1))*n.dy
            return tot
        Is = {e: I_eta_from_G(G, e, n=40) for e in etas}
        oneb[name] = (W, Is)
        P(f"  {name}:  W = (1/2)<A(F'^2)> = {W:.8f};  I_eta (eta={etas}) = " + "  ".join(f"{Is[e]:.8f}" for e in etas) +
          f";  I(0.025)-W = {Is[0.025] - W:.1e}   [{time.time() - t0:.0f}s]")

    # ---------------- two-body static pieces via Schmidt blocks
    acc = {name: np.zeros(3) for name, _ in prs}       # (2||H^1/2 A psi||^2, 2Re<Apsi,AHpsi>, W)
    mlist = (0, 1, 2, 3, 4, 6, 8, 12, 16)
    accm = {m: 0.0 for m in mlist}
    for i1, n1 in enumerate(nodes):
        for n2 in nodes:
            B = block(n1, n2)
            U, S, Vh = np.linalg.svd(B, full_matrices=False)
            r = int(np.sum(S > 1e-9*S[0])) if S[0] > 0 else 0
            if r == 0:
                continue
            u = np.array([n1.to_y(U[:, k]) for k in range(r)]).T          # N1 x r
            v = np.array([n2.to_y(Vh[k, :]) for k in range(r)]).T         # (B = U S Vh: chi = sum S_k u_k (x) row_k(Vh))
            Sr = S[:r]/np.sqrt(nrm); ww = n1.w*n2.w
            def gram(X, Y, dy):
                return X.conj().T @ Y*dy
            def hgram(n, X, Y):
                FX, FY = np.fft.fft(X, axis=0), np.fft.fft(Y, axis=0)
                return FX.conj().T @ (n.h[:, None]*FY)*n.dy/n.N
            Hu = np.fft.ifft(n1.h[:, None]*np.fft.fft(u, axis=0), axis=0)
            Hv = np.fft.ifft(n2.h[:, None]*np.fft.fft(v, axis=0), axis=0)
            Guu, Gvv = gram(u, u, n1.dy), gram(v, v, n2.dy)
            hGuu, hGvv = hgram(n1, u, u), hgram(n2, v, v)
            def two_body(a1, a2):
                au, av = a1[:, None]*u, a2[:, None]*v
                # A chi = sum_k S_k (au_k (x) v_k + u_k (x) av_k)
                X1 = np.hstack([au, u]); X2 = np.hstack([v, av])            # term j: X1[:,j] (x) X2[:,j], coeff S
                c = np.concatenate([Sr, Sr])
                G1, G2 = gram(X1, X1, n1.dy), gram(X2, X2, n2.dy)
                hG1, hG2 = hgram(n1, X1, X1), hgram(n2, X2, X2)
                Hn = np.real(c @ ((hG1*G2 + G1*hG2) @ c))                  # <A chi, H A chi>
                # A H chi = sum_k S_k [a Hu (x) v + Hu (x) av + au (x) Hv + u (x) a Hv]
                Y1 = np.hstack([a1[:, None]*Hu, Hu, au, u]); Y2 = np.hstack([v, av, Hv, a2[:, None]*Hv])
                cY = np.concatenate([Sr, Sr, Sr, Sr])
                M = gram(X1, Y1, n1.dy)*gram(X2, Y2, n2.dy)
                Rc = np.real(c @ (M @ cY))
                return 2*Hn, 2*Rc
            for name, Ff in prs:
                F1, Fp1 = Ff(n1.y); F2, Fp2 = Ff(n2.y)
                Hn, Rc = two_body(n1.K*F1, n2.K*F2)
                Wb = 0.5*np.real(Sr @ ((gram(u, (n1.K*Fp1**2)[:, None]*u, n1.dy)*Gvv + Guu*gram(v, (n2.K*Fp2**2)[:, None]*v, n2.dy)) @ Sr))
                acc[name] += ww*np.array([Hn, Rc, Wb])
            for m in mlist:
                Hn, _ = two_body(n1.K*fm(n1.y, m)[0], n2.K*fm(n2.y, m)[0])
                accm[m] += ww*Hn/2
        if i1 % 10 == 0:
            P(f"    two-body blocks: row {i1}/{len(nodes)} done [{time.time() - t0:.0f}s]")
    P("\n  profile     2||H^1/2 A psi||^2   2Re<Apsi,AHpsi>    difference        W (two-body)     W (rho_1)     lim I_eta")
    for name, _ in prs:
        Hn, Rc, Wb = acc[name]
        P(f"  {name:8s}  {Hn:16.8f}  {Rc:16.8f}  {Hn - Rc:16.8f}  {Wb:14.8f}  {oneb[name][0]:14.8f}  {oneb[name][1][0.025]:14.8f}")
    P("\n  m-scaling (complex A_m, psi_T fixed):   m   ||H^1/2 A_m psi_T||^2     m^2 <A(1)>      ratio")
    for m in mlist:
        P(f"      {m:3d}   {accm[m]:16.8f}   {m*m*PK:14.8f}   {accm[m]/max(m*m*PK, 1e-300) if m else float('nan'):8.5f}")
    ms = np.array(mlist, float); Hs = np.array([accm[m] for m in mlist])
    cf = np.polyfit(ms[2:], Hs[2:], 2)
    P(f"  quadratic fit (m>=2): {cf[0]:.6f} m^2 + {cf[1]:.6f} m + {cf[2]:.6f};  <A(1)> = {PK:.6f}")
    P(f"\n[time {time.time() - t0:.1f}s]")
    open(os.path.join(os.path.dirname(__file__), 'cache', 'np3_tstate.txt'), 'w').write('\n'.join(out) + '\n')

"""
np3b_small_eta.py -- psi_T (np3_tstate.py), m = 2 and 4 profiles: I_eta down to eta = 1/320 (rho_1 eigenvector route).
Output: cache/np3b_small_eta.txt
"""
import numpy as np, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
from np3_tstate import build, block, profile, I_eta_from_G
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)
t0 = time.time()
nodes = build(ns=28)
rho, nrm = [], 0.0
for n1 in nodes:
    R = np.zeros((len(n1.S), len(n1.S)), complex)
    for n2 in nodes:
        B = block(n1, n2); R += 2*n2.w*(B @ B.conj().T)
    rho.append(R); nrm += n1.w*np.real(np.trace(R))/2
lam_max = max(np.max(np.linalg.eigvalsh(R)) for R in rho)
eig = []
for n1, R in zip(nodes, rho):
    lam, V = np.linalg.eigh(R); keep = lam > 1e-13*lam_max
    eig.append((lam[keep]/nrm, [n1.to_y(V[:, j]) for j in np.where(keep)[0]]))
etas = (0.0125, 0.00625, 0.003125)
for m, part in ((2, 're'), (4, 're'), (4, 'im')):
    Ff = lambda y: profile(y, m, part)
    W = sum(n.w*np.sum(l*np.array([0.5*n.K*np.sum(Ff(n.y)[1]**2*np.abs(x)**2)*n.dy for x in X])) for n, (l, X) in zip(nodes, eig))
    def G(ts):
        tot = np.zeros(len(ts))
        for n, (l, X) in zip(nodes, eig):
            a = n.K*Ff(n.y)[0]
            for lj, x in zip(l, X):
                Uc = n.evolve(x, ts); Uac = n.evolve(a*x, ts)
                tot += n.w*lj*2*np.imag(np.sum(np.conj(Uac)*a[None, :]*Uc, axis=1))*n.dy
        return tot
    P(f"psi_T {part} f_{m}: W = {W:.8f}")
    for e in etas:
        I = I_eta_from_G(G, e, n=30)
        P(f"   eta = {e:.6f}: I_eta = {I:.8f}   I - W = {I - W: .3e}   (I-W)/eta^2 = {(I - W)/e**2: .1f}   [{time.time() - t0:.0f}s]")
open(os.path.join(os.path.dirname(__file__), 'cache', 'np3b_small_eta.txt'), 'w').write('\n'.join(out) + '\n')

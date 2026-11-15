"""How small can the second-order ratio Q1/Q2 get as sigma degenerates?  (dB=dC=2, A classical bit)
Continuation in the eigenvalue floor kappa: sigma = (1-4k) s + k*1."""
import numpy as np, sys
from scipy.optimize import minimize
from petz_infinitesimal import min_ratio, sig_from
S="/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/"
D=4; v=np.load(S+"infx_22_0.npy"); s=v[:16].reshape(4,4); s=(s+s.conj().T)/2
w,V=np.linalg.eigh(s); w=np.clip(w-1e-4,1e-12,None); s0=(V*w)@V.conj().T; s0/=np.trace(s0).real
L=np.linalg.cholesky(s0+1e-13*np.eye(4)); x=np.concatenate([L.real.ravel(),L.imag.ravel()])
for k in [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]:
    f=lambda y: min_ratio((1-D*k)*sig_from(y,D)+k*np.eye(D),2,2)
    best=(f(x),x)
    for rep in range(3):
        r=minimize(f,best[1]+ (0 if rep==0 else 1e-3*np.random.default_rng(rep).normal(size=x.size)),method="L-BFGS-B",options=dict(maxiter=400))
        if r.fun<best[0]: best=(r.fun,r.x)
    x=best[1]; sg=(1-D*k)*sig_from(x,D)+k*np.eye(D)
    print(f"floor={k:g}  min Q1/Q2 = {best[0]:.6f}   eig(sigma)={np.array2string(np.linalg.eigvalsh(sg),precision=3)}",flush=True)
    np.save(S+f"floor_{k:g}.npy",sg)

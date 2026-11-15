import numpy as np
from scipy.integrate import quad
from petz_search import *
SP="/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/"
x = np.load(SP+"petz_gen_22_ratio_None_0.npy"); D=4; n=2*D*D; dims=(2,2)
rho = state_from(x[:n], D); sig = state_from(x[n:], D)
rB = ptr(rho,dims,[0]); sB = ptr(sig,dims,[0])
def cpow(M, z):
    w,V = np.linalg.eigh(M); return (V*np.exp(z*np.log(w)))@V.conj().T
def rotF(t):
    # R^{t}(X) = sig^{(1+it)/2} (sB^{-(1+it)/2} X sB^{-(1-it)/2} (x) 1) sig^{(1-it)/2}
    K = cpow(sig,(1+1j*t)/2) @ np.kron(cpow(sB,-(1+1j*t)/2), np.eye(2))
    return fid(rho, K@np.kron(rB,np.eye(2))@K.conj().T)
lhs, mlF0 = gen_gap(rho,sig,dims)
print("LHS", lhs, " -logF(t=0)", mlF0, " check rotF(0)", -np.log(rotF(0)))
beta = lambda t: np.pi/2/(np.cosh(np.pi*t)+1)
avg = quad(lambda t: -beta(t)*np.log(rotF(t)), -40, 40, limit=400)[0]
print("JRSWW rotated average  -int beta0(t) log F(rho,R^{t}(N rho)) =", avg, "  (theorem: <= LHS) ->", avg <= lhs)
for t in [0,0.25,0.5,1,2,4]: print("  t=",t," -logF=", -np.log(rotF(t)))
# classical sanity: commuting states must satisfy the bound
rng=np.random.default_rng(1); worst=1e9
for _ in range(20000):
    p=rng.dirichlet(np.ones(4)*rng.choice([0.2,1,5])); q=rng.dirichlet(np.ones(4)*rng.choice([0.2,1,5]))
    a,b=gen_gap(np.diag(p).astype(complex),np.diag(q).astype(complex),dims); worst=min(worst,a-b)
print("classical worst gap (must be >=0):", worst)
# full-rank perturbation of the counterexample
rho2 = 0.999*rho+0.001*np.eye(4)/4
print("full-rank rho (min eig %.2e):"%np.linalg.eigvalsh(rho2)[0], gen_gap(rho2,sig,dims))
# random (non-optimised) states: how rare is a violation?
cnt=0; N=20000; mn=1e9
for _ in range(N):
    G=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)); r=G@G.conj().T; r/=np.trace(r).real
    G=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)); s=G@G.conj().T; s/=np.trace(s).real
    a,b=gen_gap(r,s,dims); cnt+= (a<b); mn=min(mn,a/b)
print("random Hilbert-Schmidt pairs: violations %d/%d, min ratio %.4f"%(cnt,N,mn))

import numpy as np, mpmath as mp, sys
from petz_search import state_from, gen_gap
mp.mp.dps = 60
SP="/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/"
x = np.load(SP+"petz_gen_22_ratio_None_0.npy"); D=4; n=2*D*D
rho = state_from(x[:n], D); sig = state_from(x[n:], D)
print("float64:", gen_gap(rho, sig, (2,2)))
print("eig rho", np.linalg.eigvalsh(rho)); print("eig sig", np.linalg.eigvalsh(sig))
def M(a): return mp.matrix(a.tolist())
def hfun(A, f):
    A = (A + A.H)/2; w, V = mp.eigh(A)
    return V*mp.diag([f(w[i]) for i in range(len(w))])*V.H
def ptrE(A):  # dims (2,2): trace out second
    return mp.matrix(2,2, ) if False else mp.matrix([[A[2*i,2*j]+A[2*i+1,2*j+1] for j in range(2)] for i in range(2)])
def kronI(A): # A (x) 1_2
    K = mp.zeros(4,4)
    for i in range(2):
        for j in range(2):
            for e in range(2): K[2*i+e,2*j+e]=A[i,j]
    return K
def Drel(r,s): return mp.re(mp.nsum(lambda *a:0,[0,0]) if False else (r*(hfun(r,mp.log)-hfun(s,mp.log))).__class__ and sum((r*(hfun(r,mp.log)-hfun(s,mp.log)))[i,i] for i in range(r.rows)))
r=M(rho); s=M(sig); r=(r+r.H)/2; s=(s+s.H)/2
rB=ptrE(r); sB=ptrE(s)
lhs = Drel(r,s)-Drel(rB,sB)
K = hfun(s,mp.sqrt)*kronI(hfun(sB,lambda w:1/mp.sqrt(w)))
rec = K*kronI(rB)*K.H
sr = hfun(r,mp.sqrt)
w,_ = mp.eigh((sr*rec*sr+ (sr*rec*sr).H)/2)
F = sum(mp.sqrt(max(w[i],0)) for i in range(4))**2
print("mp: LHS =", mp.nstr(lhs,20), " -logF =", mp.nstr(-mp.log(F),20), " gap =", mp.nstr(lhs+mp.log(F),20))
print("tr rec =", mp.nstr(sum(rec[i,i] for i in range(4)),20))
np.set_printoptions(precision=5, suppress=True, linewidth=200)
print("rho=\n",rho); print("sigma=\n",sig)

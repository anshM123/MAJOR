import numpy as np, mpmath as mp, sys
from mp_cmi import *
mp.mp.dps = 80
S="/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/"
v = np.load(S+"infx_22_0.npy"); D=4; s = v[:16].reshape(4,4); d = v[16:].reshape(4,4)
s=(s+s.conj().T)/2; d=(d+d.conj().T)/2
np.set_printoptions(precision=6, suppress=True, linewidth=200)
print("sigma=\n", s); print("delta=\n", d)
sm = mp.matrix(s.tolist()); dm = mp.matrix(d.tolist())
for eps in ['0.03','0.01','0.003','0.001']:
    e = mp.mpf(eps)
    rho = mp.zeros(8,8)
    for i in range(4):
        for j in range(4):
            rho[i,j] = (sm[i,j]+e*dm[i,j])/2; rho[4+i,4+j] = (sm[i,j]-e*dm[i,j])/2
    wmin = min(mp.eigh((rho+rho.H)/2, eigvals_only=True))
    I, mlF, tr = cmi_and_petz(rho, (2,2,2))
    print(f"eps={eps}: min eig rho={mp.nstr(wmin,5)}  I={mp.nstr(I,15)}  -logF={mp.nstr(mlF,15)}  gap={mp.nstr(I-mlF,8)}  ratio={mp.nstr(I/mlF,8)}  tr(rec)-1={mp.nstr(tr-1,3)}")

import sys; sys.dont_write_bytecode=True
import numpy as np, ls_core as C, json, os
L=C.lam(1000)
d=np.loadtxt('../../prl/data/d_abel_X_1000.txt',comments='#'); nu,X,ps=d[:,0],d[:,1]+1j*d[:,2],d[:,5]
s=nu<60; nu,X,ps=nu[s],X[s],ps[s]
ph,w=C.gl01(400); ph*=np.pi/2; w*=np.pi/2
for beta in (0.05,0.1,0.18,0.3,0.5):
    ss=2*np.arcsinh(np.sinh(beta/2)*np.sin(ph)); P=(2/np.pi)*(np.cos(np.outer(nu,ss))*(w/np.cosh(ss/2))).sum(1)
    sgn=(((X-ps)*P).sum()*0.04/L); ab=(np.abs(X-ps)*np.abs(P)).sum()*0.04/L
    print(beta, 'signed int (X-psi2)P dnu/L =',sgn, ' |.|-version', ab, ' OmLR2(beta)=',C.OmLR(2,[beta])[0])

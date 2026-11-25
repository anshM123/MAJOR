import sys
sys.dont_write_bytecode=True
sys.path.insert(0,'../omegau'); sys.path.insert(0,'../spinfunc')
import numpy as np
import nf_common as NC
F=NC.NF()
p,w=np.polynomial.legendre.leggauss(600); p=(p+1)/2; w=w/2
res={}
for h in ['++','+-']:
    for c in [0.19,0.195,0.199]:
        res[(h,c)]=(9,None)
    mF=(9,None)
    Jlist=[J for J in range(0,400) if (h=='++' and J%2==0) or (h=='+-' and J>=4)]
    for J in Jlist:
        umax=min(1.0,(14.0/(2*J+1))**2)
        for u in np.linspace(umax*1e-3,umax,60):
            Js,B2,B3,_=NC.flat_parts_nf(F,h,u,J,p,w); i=list(Js).index(J); B2,B3=B2[i],B3[i]
            Fv=B2+B3; N=abs(B2)+abs(B3)
            bh=(2*J+1)*np.sqrt(u)
            if Fv<mF[0]: mF=(Fv,(J,u,bh,N))
            for c in [0.19,0.195,0.199]:
                v=Fv-c*N
                if v<res[(h,c)][0]: res[(h,c)]=(v,(J,u,bh))
    print(h,'min F on bhat<=14:',mF)
    for c in [0.19,0.195,0.199]: print('  c',c,'min F-cN',res[(h,c)])

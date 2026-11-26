# Independent exact check of the mixed-profile theorem, from the paper's fp-U, fp-V and the k=3 term sigma_b L(G)/T of fp-H.
from fractions import Fraction as Fr
import math, random
exec(open('../rational_profile/rp_check.py').read().split('# 1. positive basis lemma')[0])   # reuse U, Lb, mul, binpoly, Dop, Vb (exact)
M=998; D=M+6; eps=Fr(1,(M+1)**8); A=Fr(7,10)*M*M
rho=[Fr(2*j+1,200*(D+1)**2) for j in range(1,D+1)]          # t=1/2 member of the family; R has nonneg coeffs, R(1)<=100/99
def Rser(N):
    R=[Fr(1)]+[Fr(0)]*N
    for r in rho:
        for i in range(1,N+1): R[i]+=r*R[i-1]
    return R
N=110; A0=mul(binpoly(M)[:N+1],Rser(N))[:N+1]
F=mul([Fr(0),eps,Fr(1)],A0)[:N+1]; G=[A*c for c in mul([Fr(0),Fr(0),Fr(1)],A0)[:N+1]]
worst=None; cnt=0; bad=[]
for b in (0,4):
    for J in list(range(b,b+14))+[25,40,64,99]:
        for a in (Fr(998),Fr(9983,10),Fr(1500),Fr(10**4),Fr(10**6)):
            Delta=J+a+2; T=(Delta-1)*(Delta-2); sig=1 if b==0 else -1
            H=Vb(F,b,J,Delta)+sig*Lb(G,b,J,Delta)/T          # H_b / w^3
            LF=Lb(F,b,J,Delta); cnt+=1
            if H<=0: bad.append((b,J,a))
            if J-b>=1:
                ratio=H/LF
                if worst is None or ratio<worst[0]: worst=(ratio,b,J,a)
print('mixed states',cnt,'nonpositive',bad)
print('min H/(w^3 L(F)) over n>=1:',float(worst[0]),'at b,J,a=',worst[1],worst[2],float(worst[3]),' (theorem claims >= 1/25 = 0.04)')
# K Lambda^6 = 100 F2 + 36 F3 + 8 G2 >= (28/5)M^2 - 36M + 100(1-eps M) > (11/2)M^2
KL6=100*F[2]+36*F[3]+8*G[2]; print('K Lambda^6 / M^2 =',float(KL6/M**2),' > 5.5:',KL6>Fr(11,2)*M*M)
for Dg in (1000,2*10**4,10**5,10**6): print('tau_mixed(%g)='%Dg, 420*math.sqrt(math.log(Dg))/(Dg-2)**2.5)

"""rv4_light.py -- Sec.2 of the note: exact K(t), K(1/997002), Einstein light coefficient 16 L_A vs the stated bound."""
from fractions import Fraction as Fr
import mpmath as mp
from rv_common import *
mp.mp.dps = 40
c2, c3 = PSI2, PSI3
m2 = lambda j: mom(c2, j); m3 = lambda j: mom(c3, j)
# K(t) = -2[mu26 - 245/4 t mu24 + 9619/16 t^2 mu22 - 51575/64 t^3 mu20] + 2[mu34 - 35/2 t mu32 + 569/16 t^2 mu30]
K = [-2*m2(6) + 2*m3(4),
     -2*(-Fr(245, 4))*m2(4) + 2*(-Fr(35, 2))*m3(2),
     -2*Fr(9619, 16)*m2(2) + 2*Fr(569, 16)*m3(0),
     -2*(-Fr(51575, 64))*m2(0)]
claimed = [Fr('888656439590000/888656439561891'), Fr('-2176144896481250/67209310555101'),
           Fr('-1537157051698780625/1683770096012004'), Fr('78949596683828125/4716442845972')]
print("K coefficients (mine):", [str(k) for k in K])
print("equal to the note's four exact coefficients:", [a == b for a, b in zip(K, claimed)])
D = 997002
KD = sum(k/Fr(D)**i for i, k in enumerate(K))
print("K(1/997002) =", mp.nstr(mp.mpf(KD.numerator)/KD.denominator, 16), " paper: 0.99996752 (K_B2A4: 0.999967523128)")
print("K0-1 =", float(K[0]-1), " paper: 3.2e-11;  t-coefficient", float(K[1]), " paper K-1 ~ -32.38/Lambda^2")
# min of K on [0,1/D]: K is cubic in t; check monotone/minimum
import numpy as np
ts = [Fr(i, 1000*D) for i in range(1001)]
Kmin = min(sum(k*t**i for i, k in enumerate(K)) for t in ts)
print("min K on grid of [0,1/D]:", float(Kmin), "  (> 0.99996?)", Kmin > Fr('0.99996'))
# H17 decomposition check vs K_B2A4: Lambda^2 H17 = k2 part + k3 part ; delta K0 = -68.10
lam2 = Fr(D)
k3 = 2*(-Fr(35, 2))*m3(2) + 2*Fr(569, 16)*m3(0)/lam2
print("Lambda^2 H17_{k=3} (mine) =", float(k3), " K_B2A4: -82.0165")
# Einstein: 16 L_A, L_A = int_0^1 psi2/(p^2 + 9/(4L^2))
def LA(L):
    L = mp.mpf(L)
    e2 = mp.mpf(9)/(4*L*L)
    return mp.quad(lambda p: psi_mp(c2, p)/(p*p+e2), [0, 3/L, 30/L, 0.1, 1])
A = c2[1]
print("A = psi2'(0) =", float(A), " (paper 188.896)")
g = c2[2:]
print("int g = ", float(mom(g, 0)), " (note: -332.95059691167796)")
print("Lambda | 16 L_A | 3023 log(2L/3)-5307 | slack | 16[A log(2L/3)+int g] ")
for L in [mp.sqrt(997002), 1500, 3000, 1e4, 1e5, 1e7]:
    v = 16*LA(L)
    bnd = 3023*mp.log(2*mp.mpf(L)/3) - 5307
    asy = 16*(float(A)*mp.log(2*mp.mpf(L)/3) + float(mom(g, 0)))
    print(mp.nstr(L, 8), mp.nstr(v, 10), mp.nstr(bnd, 10), mp.nstr(bnd-v, 5), mp.nstr(asy, 10))
print("paper 16 L_A at Delta_g=1000: 14339.59")
# Bernstein range of g claimed -511.773 < g <= -188.896: check by dense sampling
gs = [psi_mp(g, mp.mpf(i)/20000) for i in range(20001)]
print("g range on dense grid:", mp.nstr(min(gs), 8), mp.nstr(max(gs), 8))

"""Exact (rational) 3-state primal configuration for the LP {endpoint rule, ANEC harmonic rules R_k for all k}.
States (all at the gap, twist exactly g):  A=(b=0,J=2,Delta=g+2), B=(b=4,J=4,Delta=g+4), C=(b=4,J=6,Delta=g+6).
They satisfy R_k exactly for every k != 2 (A has ladder entries k<=2, B k=0, C k<=2), the endpoint rule with total 3,
and give |c|^2 = rho_C r_2^(4)(C) - rho_A r_2^(0)(A) > 0.  Hence NO bound derived from these rules can be below
tau_LP(g) = 4 sqrt(|c|^2).  Exact rationals (fractions)."""
from fractions import Fraction as Fr
import math, sys

def T(D): return (D-1)*(D-2)
def ladder(b, J, D, kmax):
    C = J*(J+1); TT = T(D); r = [Fr(1)]
    for k in range(kmax):
        s = b+k
        if k+1 > J-b: r.append(Fr(0)); continue
        r.append(r[-1]*Fr((k+4)**2, (k+1)**2)*Fr(C - s*(s+1), 1)/Fr(TT - s*(s+1), 1))
    return r
def W(b, D): return Fr(T(D) - b*(b-1), 3) + Fr(3, 2)

def primal(g):
    g = Fr(g)
    DA, DB, DC = g+2, g+4, g+6
    rA = ladder(0, 2, DA, 3); rC = ladder(4, 6, DC, 3)
    # unknowns rhoA, rhoB, rhoC: R0: rhoB+rhoC=rhoA ; R1: rhoC rC1 = rhoA rA1 ; budget sum rho W = 3
    # set rhoC=1, scale
    rhoC = Fr(1); rhoA = rhoC*rC[1]/rA[1]; rhoB = rhoA - rhoC
    assert rhoB >= 0
    budget = rhoA*W(0, DA) + rhoB*W(4, DB) + rhoC*W(4, DC)
    s = Fr(3)/budget
    rhoA, rhoB, rhoC = rhoA*s, rhoB*s, rhoC*s
    c2 = rhoC*rC[2] - rhoA*rA[2]
    return c2, (rhoA, rhoB, rhoC)

if __name__ == '__main__':
    for g in [10, 20, 50, 100, 300, 1000, 10**4, 10**6]:
        c2, rho = primal(g)
        tau = 4*math.sqrt(c2)
        print(f"g={g}: |c|^2={float(c2):.6e}  tau_LP={tau:.6e}  tau*g^3={tau*g**3:.3f}  rho={[float(x) for x in rho]}")
    print('asymptotic: |c|^2 g^6 -> 21600, tau g^3 -> 4*sqrt(21600) =', 4*math.sqrt(21600))

"""(Delta, J) pairs at (approximately) fixed impact parameter beta and prescribed m:
   alpha = Delta-J-2 = m e^{-beta/2}, alphabar = Delta+J-1 = m e^{beta/2};  J integer  =>  m = (2J+1)/(2 sinh(beta/2))."""
import numpy as np, mpmath as mp

def case(beta, J):
    m = (2*J + 1)/(2*np.sinh(beta/2))
    al, alb = m*np.exp(-beta/2), m*np.exp(beta/2)
    Delta = float((al + alb + 3)/2)
    return Delta, J

def info(Delta, J):
    D = mp.mpf(Delta)
    al, alb = D - J - 2, D + J - 1
    m = mp.sqrt(al*alb)
    return m, mp.log(alb/al), (al*al + alb*alb)/(2*al*alb)

def cases_fixed_beta(beta, mlist):
    out = []
    for m in mlist:
        J = int(round((2*m*np.sinh(beta/2) - 1)/2))
        J = max(J, 0)
        out.append(case(beta, J))
    return out

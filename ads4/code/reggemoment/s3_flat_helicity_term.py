"""
Stage 3, flat-space benchmark for the helicity-dependent 1/m^2 term in the helicity-preserving components.

Flat partial waves (arXiv:2201.06602 (A.1)-(A.3)):  dtilde^J_{lam,lam}(x) = 2F1(lam - J, J + lam + 1; 1; (1-x)/2),
x = 1 - 2 p^2/m^2;  lam = h_1 - h_2 = 0 for ++ and 4 for +-.   Derived here (series in p^2/m^2 at fixed L^2 = J(J+1) = b^2 m^2/4):
   dtilde^J_{lam,lam} - dtilde^J_{0,0} = (1/m^2) [ 2 lam^2 (p/b) J_1(b p) + lam p^2 J_0(b p) ] + O(1/m^4),
   d^J_{lam,lam} - d^J_{0,0}           = (1/m^2) [ 2 lam^2 (p/b) J_1(b p) ] + O(1/m^4)          (unstripped Wigner functions).
"""
import mpmath as mp
mp.mp.dps = 30
def dt(J, lam, x):
    return mp.hyp2f1(lam - J, J + lam + 1, 1, (1 - x)/2)
print("#    b     p     m      m^2 (dt_44 - dt_00)    prediction    m^2 (d_44 - d_00)   prediction")
for (b, p) in ((3, mp.mpf('0.7')), (1, mp.mpf('0.9')), (6, mp.mpf('0.3')), (2, mp.mpf('2.5'))):
    for J in (40, 160, 640, 2560):
        m = 2*mp.sqrt(mp.mpf(J)*(J+1))/b
        x = 1 - 2*p*p/m**2
        lam = 4
        diff = m**2*(dt(J, lam, x) - dt(J, 0, x))
        pred = 2*lam**2*(p/b)*mp.besselj(1, b*p) + lam*p*p*mp.besselj(0, b*p)
        diff2 = m**2*(((1+x)/2)**lam*dt(J, lam, x) - dt(J, 0, x))
        pred2 = 2*lam**2*(p/b)*mp.besselj(1, b*p)
        print(f"  {b:4d}  {mp.nstr(p,3):>4s}  {mp.nstr(m,6):>8s}   {mp.nstr(diff, 10):>14s}   {mp.nstr(pred, 10):>14s}   {mp.nstr(diff2, 10):>14s}   {mp.nstr(pred2, 10):>14s}")

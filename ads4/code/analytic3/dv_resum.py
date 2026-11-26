"""DV (design check for the hump resummation, ANALYTIC_SLIVER3 §6): exact S = sum_j c_j d^j x_j vs the frozen-matrix resummation
S0 = e1^T H(d Mbar) V0 (H = 2F1(a0,a0;b0;.), Mbar = M_{j*}) and the size of the termwise remainder sum |c_j| d^j |x_j - xbar_j|."""
import sys, math
sys.path.insert(0, '../analytic2'); sys.path.insert(0, '../analytic')
import mpmath as mp
from common2 import n_of_L
from fastS import S_vec
mp.mp.dps = 60
L, nu, a = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]); jstar = int(sys.argv[4]) if len(sys.argv) > 4 else 0
s = 1; n = n_of_L(L, a, 0); N = mp.mpf(n) + a; D = N*(N + s); d = 1 - mp.mpf(n*(n + s))/D; g = 1/(d*d*D)
al = mp.mpc(2.5, nu)
S0v, lg = S_vec(n, a, s, [complex(2.5, nu), complex(3.5, nu)])
x = {-1: mp.mpc(complex(S0v[1]))*mp.e**lg[1]*d**(al + 1), 0: mp.mpc(complex(S0v[0]))*mp.e**lg[0]*d**al}
def AB(beta):
    den = 1 - beta*(beta - s)*g*d*d
    return (1 + beta*(s - 2*beta - 1)*g*d)/den, beta*(beta + 1)*g/den
J = 140
for j in range(J):
    A, B = AB(al - j); x[j + 1] = A*x[j] + B*x[j - 1]
a0 = mp.mpf(0.5) - 1j*nu; b0 = 1 - 2j*nu
c = [mp.mpc(1)]
for j in range(J): c.append(c[-1]*(a0 + j)**2/((b0 + j)*(j + 1)))
S = sum(c[j]*d**j*x[j] for j in range(J))
term = sum(abs(c[j])*d**j*abs(x[j]) for j in range(J))
A, B = AB(al - jstar); Mb = mp.matrix([[A, B], [1, 0]])
xb = {-1: x[-1], 0: x[0]}
for j in range(J): xb[j + 1] = A*xb[j] + B*xb[j - 1]
Sb = sum(c[j]*d**j*xb[j] for j in range(J))
err = sum(abs(c[j])*d**j*abs(x[j] - xb[j]) for j in range(J))
lam = [(A + mp.sqrt(A*A + 4*B))/2, (A - mp.sqrt(A*A + 4*B))/2]
H = lambda y: mp.hyp2f1(a0, a0, b0, y)
# e1^T H(d Mbar) V0 via eigen-decomposition
V0 = mp.matrix([x[0], x[-1]])
l1, l2 = lam
# x_bar_j = C1 l1^j + C2 l2^j with C1 + C2 = x0, C1/l1 + C2/l2 = x_{-1}
C1 = (x[-1] - x[0]/l2)/(1/l1 - 1/l2); C2 = x[0] - C1
Sres = C1*H(d*l1) + C2*H(d*l2)
print(f'L={L} nu={nu} a={a} kappa={nu/a:.3f} j*={jstar}: |S| = {mp.nstr(abs(S),5)}  termwise = {mp.nstr(term,5)}  |S_frozen| = {mp.nstr(abs(Sb),5)}'
      f'  |closed form C1 H(d l1) + C2 H(d l2)| = {mp.nstr(abs(Sres),5)} (diff {mp.nstr(abs(Sres - Sb),3)})  termwise remainder = {mp.nstr(err,5)}')
print(f'   xi = nu^3 d^2/(4 a^2) = {nu**3*float(d)**2/(4*a*a):.3f};  |C1| = {mp.nstr(abs(C1),4)}, |C2| = {mp.nstr(abs(C2),4)}, |H(d l1)| = {mp.nstr(abs(H(d*l1)),4)}, |H(d l2)| = {mp.nstr(abs(H(d*l2)),4)}')

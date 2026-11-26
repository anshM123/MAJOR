# check branch formula F_{3,nu}(c) = 2 Re[B(nu) sum_j c_j(nu) d^{j-alpha}], and its continuation to complex nu
from flint import arb, acb, ctx
ctx.prec = 300
def Fplus(nu, c, J=4000):
    nu = acb(nu); d = 1 - arb(c)
    i = acb(0,1)
    B = (2*i*nu).gamma()/((acb(0.5)+i*nu).gamma())**2
    al = acb(2.5) + i*nu
    a = acb(0.5) - i*nu; bb = acb(1) - 2*i*nu
    cj = acb(1); tot = acb(0); dl = acb(d)
    pw = dl**(-al)
    for j in range(J):
        tot += cj*pw
        cj = cj*(a+j)*(a+j)/((bb+j)*(j+1)); pw *= dl
    return B*tot
c = arb('0.445')
for nu in [0.7, 5, 30]:
    fp = Fplus(nu, c)
    X = acb((1+c)/(1-c)); Fex = ((1-c)**(-3)*X.legendre_p(acb(-0.5, nu), 0, 3)).real
    print(nu, 2*fp.real, Fex)
# outgoing via Legendre Q: P+ = -(i/pi) coth(pi nu) Q_{-1/2-i nu}(cosh beta)
import mpmath as mp
mp.mp.dps=40
for nu in [0.7, 5, 30, mp.mpc(5,3)]:
    cc=mp.mpf('0.445'); X=(1+cc)/(1-cc)
    q = -1j/mp.pi*mp.coth(mp.pi*nu)*mp.legenq(-0.5-1j*nu, 0, X, type=3)*(1-cc)**-3
    fp = Fplus(complex(nu), c)
    print(nu, q, fp)

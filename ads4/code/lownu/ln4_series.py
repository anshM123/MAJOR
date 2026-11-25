"""Check of the reviewer's analytic singular de-shift term (k = 2, EE): R_-(-i), R_-'/R_- at -i, and the small-nu structure of
E_sing(nu) = 2 Re[R_-(nu-i) chi_sing(nu-i)],  chi_sing(z) = c w log w (1 + (37/36) w) e^{-z^2/2},  w = z^2+1,  c = -25 kappa/(36 pi sqrt e).
Exact (sympy)."""
import sympy as sp
from ln_common import guard
nu, k = sp.symbols('nu kappa', positive=True)
z = sp.symbols('z')
Rm = -2*(2*z - sp.I)*(2*z + sp.I)**2*(2*z + 3*sp.I)*(2*z + 5*sp.I)*(2*z + 11*sp.I)/(z*(2*z - 7*sp.I)*(2*z - 3*sp.I)*(2*z + 7*sp.I))
out = []
out.append('R_-(-i) = %s' % sp.nsimplify(sp.simplify(Rm.subs(z, -sp.I))))
out.append("R_-'(-i)/R_-(-i) = %s" % sp.simplify((sp.diff(Rm, z)/Rm).subs(z, -sp.I)))
c = -25*k/(36*sp.pi*sp.sqrt(sp.E))
A = sp.symbols('A')   # coefficient replacing 37/36
zz = nu - sp.I; w = sp.expand(zz**2 + 1)        # = nu^2 - 2 i nu
h = c*Rm.subs(z, zz)*w*(1 + A*w)*sp.exp(-zz**2/2)  # multiplies log w = log(2 nu) - i pi/2 + log(1 + i nu/2)
L = sp.symbols('L', real=True)                  # L stands for log(2 nu)
hs = sp.series(h, nu, 0, 6).removeO()
lw = L - sp.I*sp.pi/2 + sp.series(sp.log(1 + sp.I*nu/2), nu, 0, 6).removeO()
Eser = sp.expand(2*sp.re(sp.expand(hs*lw)).subs({sp.re(L): L, sp.im(L): 0}))
Eser = sp.expand(sp.simplify(Eser))
for p in range(0, 5):
    cL = sp.simplify(Eser.coeff(L).coeff(nu, p)); c0 = sp.simplify(Eser.subs(L, 0).coeff(nu, p))
    out.append('nu^%d: log coeff = %s ; plain coeff = %s' % (p, cL, c0))
sol = sp.solve(sp.Eq(sp.simplify(Eser.coeff(L).coeff(nu, 2)), 0), A)
out.append('A solving the nu^2 log nu cancellation: %s' % sol)
txt = '\n'.join(out); print(txt)
open(guard('cache/ln4_series.txt'), 'w').write(txt + '\n')

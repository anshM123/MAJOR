"""Referee r3: c2's forward flip estimate uses a one-sided difference in nu^2 of g(nu) = rho_sc(nu) K_{4,0}(nu,beta), which is not even in nu.
Compare with d/d(nu^2) at 0 of the even part (g(nu)+g(-nu))/2 and size of the nu-odd part."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'superconv'))
from c1_flip_contamination import K
def rho_sc(nu, k=4):
    c = {4: (210, 216, 263), 6: (1890, 536, 2023)}[k]
    return c[0]*(2*nu - 5j)*(2*nu - 1j)/((2*nu + 3j)*(2*nu + 7j)*(16*nu**4 + c[1]*nu**2 + c[2]))
out = []
for b in (0.0363, 0.0761, 0.1437):
    g = lambda nu: rho_sc(nu)*K(nu, b, 4, 0)
    h = 0.05
    dcode = (g(np.sqrt(2*h)) - g(np.sqrt(h)))/h
    ge = lambda nu: (g(nu) + g(-nu))/2
    eps = 1e-2
    deven = (ge(eps) - ge(0.0))/eps**2
    godd = (g(eps) - g(-eps))/(2*eps)
    s = f"beta={b}: |code d|={abs(dcode):.3e}  |d_nu2 even part|={abs(deven):.3e}  |g(0)|={abs(g(0.0)):.3e}  |odd slope dg/dnu|={abs(godd):.3e}"
    print(s); out.append(s)
open(sys.argv[1], 'w').write("\n".join(out) + "\n")

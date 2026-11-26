"""probe1: exact |S^(alpha - j)| at the binding state (a = 998, L ~ 4.5, b = 0), nu = 1031.5, and the true kernel vs termwise."""
import sys, math, time
sys.path.insert(0, '.')
from common2 import *
from flint import arb, acb, ctx
ctx.prec = 200
L = float(sys.argv[1]); nu = float(sys.argv[2]); a = float(sys.argv[3]) if len(sys.argv) > 3 else 998.0
b = int(sys.argv[4]) if len(sys.argv) > 4 else 0; s = 2*b + 1
JJ = int(sys.argv[5]) if len(sys.argv) > 5 else 40
n = n_of_L(L, a, b)
d = dstate(n, a, s)
print('n', n, 'L', float(-d.log()), flush=True)
I = acb(0, 1)
al = acb(2.5, nu)
B = (2*I*nu).gamma()*((acb(0.5) + I*nu).rgamma())**2
c = acb(1); a0 = acb(0.5) - I*nu; b0 = acb(1) - 2*I*nu
tot = acb(0); absum = 0.0; t0 = time.time()
for j in range(JJ):
    S = S_rec(n, a, s, al - j, 200)
    Sh = S*(d**(al - j))
    term = c*S
    tot += term
    absum += float(abs(c).mid())*float((d**(j - 2.5)).mid())*float(abs(Sh).mid())
    print(j, f'|S^|={float(abs(Sh).mid()):.4e}  |c_j| d^j={float(abs(c).mid())*float((d**j).mid()):.3e}  partial |K|/|B|d^-5/2 = {float(abs(tot*d**al).mid()):.4e}   [{time.time()-t0:.0f}s]', flush=True)
    c = c*(a0 + j)*(a0 + j)/((b0 + j)*(j + 1))
K = 2*(B*tot).real
print('K_true ~', K, ' |B| =', abs(B), ' d^-5/2 =', d**(-2.5))
print('termwise |B| sum |c_j| d^{j-5/2}|S^_j| =', float(abs(B).mid())*absum)

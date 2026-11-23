#!/usr/bin/env python3
"""
c1val / v1: the u-channel (shock / CRT) route in d = 3, in the exact Costa-Hansen-Penedones basis.

CHP (arXiv:1707.07689, eq. DH3_grav) write the impact-parameter amplitude of <T T O O> on H_{d-1} as
    sum_k beta_k(nu) D^_k Omega_{i nu}(L),   k = 1 (Einstein), 2 (t2-type), 3 (t4-type),
with
    D^_2 = delta^{{mu}_{{rho}} nabla^{nu}} nabla_{sigma}} + (nu^2 + ((d-2)/2)^2)/(d-1) * delta delta
    D^_3 = 1/2 ( nabla^{mu} nabla_{rho} nabla^{nu} nabla_{sigma} + nabla_{rho} nabla^{mu} nabla_{sigma} nabla^{nu} )
           + (4 nu^2 + d^2 - 2d - 2)/(2(d+1)) * delta^{{mu}_{{rho}} nabla^{nu}} nabla_{sigma}}
contracted with eps_mu eps_nu eps*^rho eps*^sigma (CHP eq. `cases`: z1 = z2^* = eps, i.e. the helicity-PRESERVING element).

Checks (exact sympy on H_2, upper half plane, arbitrary scalar f(x, y)):
 (A) T(e,e*,e,e*) := e^a e*^b e^c e*^d nabla_a nabla_b nabla_c nabla_d f  =  (1/4) Box^2 f   exactly (curvature terms cancel)
     and the same for T(e*,e,e*,e).
 (B) e^a e*^b nabla_a nabla_b f = (1/2) Box f exactly.
 (C) Hence diag D^_3 = (1/4) Box (Box + nu^2 + 1/4) f and diag D^_2 = (1/2)(Box + nu^2 + 1/4) f in d = 3:
     both vanish identically on every H_2 eigenfunction with Box = -(nu^2 + 1/4) (Omega_{i nu}, Pi_{i nu}, any L, any nu).
 (D) Control: the flip element e^a e^b e^c e^d nabla^4 f is NOT zero (a helicity-4 operator).
 (E) CHP's large-L collider form (eq. tensors_n) in d = 3 with circular eps: both the t2 and t4 structures vanish.
Writes cache/v1_chp_d3.txt.
"""
import os, sys
import sympy as sp

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'v1_chp_d3.txt')
lines = []
def log(s=''):
    print(s); lines.append(str(s))

x, y = sp.symbols('x y', real=True)
nu = sp.symbols('nu', real=True)
X = [x, y]
f = sp.Function('f')(x, y)
g = sp.diag(1/y**2, 1/y**2)
gi = g.inv()
n = 2
Gam = [[[sp.simplify(sum(gi[a, e]*(sp.diff(g[e, b], X[c]) + sp.diff(g[e, c], X[b]) - sp.diff(g[b, c], X[e]))
                         for e in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]

def cov(T, rank):
    """T: dict tuple->expr of a covariant rank-`rank` tensor; returns rank+1 tensor, new index FIRST."""
    import itertools
    R = {}
    for idx in itertools.product(range(n), repeat=rank + 1):
        b, rest = idx[0], idx[1:]
        val = sp.diff(T[rest], X[b])
        for i, ai in enumerate(rest):
            for c in range(n):
                r2 = list(rest); r2[i] = c
                val -= Gam[c][b][ai]*T[tuple(r2)]
        R[idx] = val
    return R

T0 = {(): f}
T1 = cov(T0, 0)
T2 = cov(T1, 1)
T3 = cov(T2, 2)
T4 = cov(T3, 3)

box = lambda F: sp.simplify(sum(gi[a, b]*cov(cov({(): F}, 0), 1)[(a, b)] for a in range(n) for b in range(n)))
Boxf = box(f)
Box2f = box(Boxf)

e = [y/sp.sqrt(2), sp.I*y/sp.sqrt(2)]          # contravariant comps of (e_x + i e_y)/sqrt2, orthonormal frame
ec = [y/sp.sqrt(2), -sp.I*y/sp.sqrt(2)]
def contract(T, vecs):
    import itertools
    return sp.expand(sum(T[idx]*sp.prod([vecs[k][idx[k]] for k in range(len(vecs))])
                         for idx in itertools.product(range(n), repeat=len(vecs))))
ee_c = sp.simplify(sum(g[a, b]*e[a]*ec[b] for a in range(n) for b in range(n)))
log('v1: CHP d=3 impact-parameter operators on H_2 (exact sympy)')
log(f'  e.e* = {ee_c}   e.e = {sp.simplify(sum(g[a,b]*e[a]*e[b] for a in range(n) for b in range(n)))}')

A1 = sp.simplify(contract(T4, [e, ec, e, ec]) - Box2f/4)
A2 = sp.simplify(contract(T4, [ec, e, ec, e]) - Box2f/4)
log(f'(A) T(e,e*,e,e*) - Box^2 f/4 = {A1};   T(e*,e,e*,e) - Box^2 f/4 = {A2}')
B = sp.simplify(contract(T2, [e, ec]) - Boxf/2)
log(f'(B) e^a e*^b nabla_a nabla_b f - Box f/2 = {B}')

d = 3
c3 = sp.Rational(1, 2*(d+1))*(4*nu**2 + d**2 - 2*d - 2)
c2 = (nu**2 + sp.Rational(d-2, 2)**2)/(d-1)
diag3 = sp.Rational(1, 2)*(contract(T4, [e, ec, e, ec]) + contract(T4, [ec, e, ec, e])) + c3*ee_c*contract(T2, [e, ec])
diag2 = ee_c*contract(T2, [e, ec]) + c2*ee_c**2*f
lam = sp.symbols('lam')
# reduce on eigenfunctions: replace Box f -> lam f, Box^2 f -> lam^2 f
r3 = sp.simplify(diag3 - (Box2f/4 + c3*Boxf/2))
r2 = sp.simplify(diag2 - (Boxf/2 + c2*f))
log(f'(C) diag D3 - [Box^2/4 + c3(nu) Box/2] f = {r3},  c3(nu) = {sp.simplify(c3)}')
log(f'    diag D2 - [Box/2 + c2(nu)] f = {r2},  c2(nu) = {sp.simplify(c2)}')
ev = -(nu**2 + sp.Rational(1, 4))
d3e = sp.simplify(ev**2/4 + c3*ev/2); d2e = sp.simplify(ev/2 + c2)
log(f'    on Box = -(nu^2+1/4):  diag D3 -> {d3e},  diag D2 -> {d2e}   (for every nu, every L)')

# explicit eigenfunctions (plane waves on H_2): y^(1/2+i nu) and its image based at the boundary point x = 1
fe = y**(sp.Rational(1, 2) + sp.I*nu)
chk_e = sp.simplify(sp.simplify(box(fe)/fe))
sub3 = sp.simplify((diag3.subs(f, fe).doit()))
log(f'    explicit eigenfunction y^(1/2+i nu): Box f/f = {chk_e};  diag D3[f] = {sp.simplify(sub3)}')
sub2 = sp.simplify(diag2.subs(f, fe).doit())
log(f'                                        diag D2[f] = {sub2}')
fe2 = (x**2 + y**2)**0 * (y/((x - 1)**2 + y**2))**(sp.Rational(1, 2) + sp.I*nu)   # boundary-point eigenfunction
s3b = sp.simplify(diag3.subs(f, fe2).doit())
log(f'    eigenfunction (y/((x-1)^2+y^2))^(1/2+i nu):  diag D3[f] = {s3b}')

flip = contract(T4, [e, e, e, e])
fl = sp.simplify(flip.subs(f, fe2).doit())
log(f'(D) control: flip element e e e e nabla^4 f on the same eigenfunction is nonzero: {sp.simplify(fl) != 0}')
num = complex(sp.N(fl.subs({x: sp.Rational(3, 10), y: sp.Rational(7, 10), nu: sp.Rational(13, 10)})))
log(f'    value at (x,y,nu) = (0.3,0.7,1.3): {num:.6e}')

# (E) CHP collider form, d = 3, circular eps
phi = sp.symbols('phi', real=True)
epsv = sp.Matrix([1, sp.I])/sp.sqrt(2); nvec = sp.Matrix([sp.cos(phi), sp.sin(phi)])
ne2 = sp.simplify(sp.Abs((nvec.T*epsv)[0])**2)
t2s = sp.simplify(ne2 - sp.Rational(1, d-1)); t4s = sp.simplify(ne2**2 - sp.Rational(2, (d+1)*(d-1)))
log(f'(E) |n.eps|^2 = {ne2}: t2 structure {t2s}, t4 structure {t4s}  (CHP eq. tensors_n, d = 3)')

ok = (A1 == 0 and A2 == 0 and B == 0 and r3 == 0 and r2 == 0 and d3e == 0 and d2e == 0
      and sp.simplify(sub3) == 0 and sp.simplify(sub2) == 0 and s3b == 0 and t2s == 0 and t4s == 0)
log(f'ALL CHECKS PASS: {ok}')
with open(OUT, 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
sys.exit(0 if ok else 1)

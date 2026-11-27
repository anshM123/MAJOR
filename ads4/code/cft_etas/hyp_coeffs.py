"""General first-order response on the hyperbolic branch and the Taylor coefficients of V at n = 1 (T = T0).
Inputs [DV, hyp_num.py at r_h = 1, 5/4, 3/2, 2, 3, 5, 10; exact rational agreement]:
   e1[theta I] = 4 t0 - 2 f'(r_h) t1,    e1[theta J] = x(-48 - 54 x) t0 + 6 x f'(r_h) t1,   f'(r_h) = 3 r_h - 1/r_h = 4 pi T.
Wald [D, type-D bivector projection, calibrated on the planar engine]: sigma1[theta I] = 4 x t0, sigma1[theta J] = 6 x^2 t0.
Background at the horizon: I = 12 x^2, J = 12 x^3, r_h I' = -72 x^2, r_h J' = -108 x^3, and f'(r_h)/r_h = 2 + x.
So delta = e1 - sigma1 = 4(1-x) th_I - 2(2+x) rh th_I' + (-48x - 60x^2) th_J + 6x(2+x) rh th_J'   (th = F_y at the horizon).
Planar check x = 1: -6 th_I' - 108 th_J + 18 th_J' (paper).  V = delta + (1-x) delta'(x)."""
import sympy as sp
x, n, lam = sp.symbols('x n lambda')
Is, Js = sp.symbols('I J')
def delta_of(F):
    FI, FJ = sp.diff(F, Is), sp.diff(F, Js)
    sub = {Is: 12*x**2, Js: 12*x**3}
    thI, thJ = FI.subs(sub), FJ.subs(sub)
    dthI = (sp.diff(FI, Is)*(-72*x**2) + sp.diff(FI, Js)*(-108*x**3)).subs(sub)
    dthJ = (sp.diff(FJ, Is)*(-72*x**2) + sp.diff(FJ, Js)*(-108*x**3)).subs(sub)
    return sp.expand(4*(1 - x)*thI - 2*(2 + x)*dthI + (-48*x - 60*x**2)*thJ + 6*x*(2 + x)*dthJ)
xn = 1 - 9*n**2/(1 + sp.sqrt(1 + 3*n**2))**2
print('x(n) about n=1:', sp.series(xn, n, 1, 4))
ops = {'J': lam*Js, 'I^2': lam*Is**2, 'IJ': lam*Is*Js, 'J^2': lam*Js**2, 'I': lam*Is, 'I^3': lam*Is**3}
planar_ref = {'J': -108, 'I^2': 864, 'IJ': -1944, 'J^2': -6480, 'I': 0, 'I^3': 31104}
for name, F in ops.items():
    d = delta_of(F)
    V = sp.expand(d + (1 - x)*sp.diff(d, x))
    Vn = sp.series(V.subs(x, xn), n, 1, 4).removeO()
    co = [sp.nsimplify(sp.expand(Vn).coeff(n - 1, k)) if k else sp.nsimplify(Vn.subs(n, 1)) for k in range(4)]
    co = [sp.simplify(sp.diff(V.subs(x, xn), n, k).subs(n, 1)/sp.factorial(k)) for k in range(4)]
    print('%-4s delta(x) = %-28s V(x) = %-34s  V(1)=%s (planar table %s)   Taylor a_k at n=1 (k=0..3): %s' % (
        name, sp.factor(d), sp.expand(V), sp.simplify(V.subs(x, 1)/lam), planar_ref[name], [sp.nsimplify(c/lam) for c in co]))

# Independent low side of NF02 from X, Y in functionals_noforward.json (sympy exact).
import json, sympy as sp
p,m=sp.symbols('p m',positive=True)
D=json.load(open('/home/user/MAJOR/ads4/code/robust/functionals_noforward.json'))['NF02']
r=D['r']; X=sum(sp.Rational(c)*p**i for i,c in enumerate(D['X'])); Y=sum(sp.Rational(c)*p**i for i,c in enumerate(D['Y']))
psi2=sp.expand((1-p)**r*p*X); psi3=sp.expand((1-p)**r*p*Y)
print('psi2 matches json:',sp.Poly(psi2,p).all_coeffs()[::-1]==[sp.Rational(c) for c in D['psi2_exact']])
print('psi3 matches json:',sp.Poly(psi3,p).all_coeffs()[::-1]==[sp.Rational(c) for c in D['psi3_exact']])
mu=lambda f,j: sp.integrate(f*p**j,(p,0,1))
k=-2*mu(psi2,6)+2*mu(psi3,4)
X0=X.subs(p,0)
# L_fin = lim (int_eps^1 psi2/p^2 - X0 log(1/eps))
g=sp.expand(sp.cancel((psi2/p**2 - X0/p)))   # polynomial
Lfin=sp.integrate(g,(p,0,1))   # since int_eps^1 X0/p = X0 log(1/eps)
A=16*X0/k; B=16*Lfin/k; BAdS=B-A*sp.log(sp.Rational(3,2))
print('k =',k,float(k)); print('A =',float(A),' B =',float(B),' B_AdS =',float(BAdS))
print('psi3(0)=',psi3.subs(p,0),' psi3 deriv at 0:',sp.diff(psi3,p).subs(p,0))
print('large-b condition -2psi2\'\'(0) vs |psi3\'\'(0)|:',float(-sp.diff(psi2,p,2).subs(p,0)),float(abs(sp.diff(psi3,p,2).subs(p,0))))
C0=mu(psi2,2)-2*m**2*mu(psi2,0)-mu(psi3,0)
C2=sp.expand((sp.integrate(psi2*(p**2-2*m**2)*(m**4-6*m**2*p**2+6*p**4),(p,0,1))-sp.integrate(psi3*(m**4-6*m**2*p**2+6*p**4),(p,0,1)))/m**4)
print('C0 =',C0); print('C0 at m->0:',float(C0.subs(m,0)))
num2=sp.Poly(sp.expand(C2*m**4),m)
print('m^4 C2 =',num2.as_expr())
import numpy as np
mm=np.linspace(1e-4,1,20001)
f0=sp.lambdify(m,C0); f2=sp.lambdify(m,num2.as_expr())
print('max C0 on (0,1]:',max(f0(mm)),' max m^4C2:',max(f2(mm)))
print('sturm/real roots of m^4C2 in (0,1]:',[rt for rt in sp.Poly(num2.as_expr(),m).real_roots() if 0<rt<=1])

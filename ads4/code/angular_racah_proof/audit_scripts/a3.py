import sympy as S
y,d1,pm,p0,pp,q0,q1,h0=S.symbols('y d1 pm p0 pp q0 q1 h0')
h1=2/p0-h0
dj=(p0*d1+q0*y)/pm
y1=y+d1
dj2=(p0*d1-q1*y1)/pp
Ej=y**2+h0*pm*p0/q0*dj*d1
Ej1=y1**2+h1*p0*pp/q1*d1*dj2
br=p0**2*(h1/q1-h0/q0)+p0/2*(h0-h1)
print('decrement', S.simplify(Ej1-Ej-br*d1**2))
# square form: E_j in terms of y and delta_j (eliminate delta_{j+1})
dj_=S.Symbol('dj')
d1_=(pm*dj_-q0*y)/p0
Ej_b=y**2+h0*pm*p0/q0*dj_*d1_
print('square', S.simplify(Ej_b-((1-h0*q0/4)*y**2+h0*pm**2/q0*(dj_-q0*y/(2*pm))**2)))
# h_j identity and closed form
t,j=S.symbols('t j',positive=True)
for jj in range(0,8):
    hj=2*S.integrate(t**jj*(1+t**8)/(1+t),(t,0,1))
    hj1=2*S.integrate(t**(jj+1)*(1+t**8)/(1+t),(t,0,1))
    pj=S.Rational((jj+1)*(jj+9),2*jj+10)
    assert S.simplify(hj+hj1-2/pj)==0
    Hj=S.Rational(2,jj+1)+S.Rational(2,jj+9)-S.Rational(1,jj+2)-S.Rational(1,jj+10)
    assert S.N(hj)<S.N(Hj) and S.N(hj1)<S.N(hj)
    if jj==5: print('h5 =',S.nsimplify(hj), float(hj), 'H5=',Hj)
print('h identities OK')
pj=(j+1)*(j+9)/(2*j+10)
print('pbound', S.factor(S.simplify(pj/(2*j+9)-S.Rational(1,8))))
Hj=2/(j+1)+2/(j+9)-1/(j+2)-1/(j+10)
num,den=S.fraction(S.factor(S.simplify(S.Rational(4,3)-(2*j+9)/4*Hj)))
print('coerc', S.factor(num), '/', den)
vv=S.Symbol('v')
print('  shifted', S.expand(S.expand(num).subs(j,5+vv)), ' vs 2v^4+81v^3+1075v^2+4806v+3528 x?')
print('  ratio', S.simplify(S.expand(num.subs(j,5+vv))/(2*vv**4+81*vv**3+1075*vv**2+4806*vv+3528)))
# also j=1..4 sign (not needed)
print([S.Rational(4,3)-(2*jj+9)*Hj.subs(j,jj)/4 for jj in range(1,6)])
# zbound
n=S.Symbol('n')
print('zb', S.expand((2*n-2)*(2*j+2*n+7)-4*n**2))
# decrement sign given hypotheses: bracket <= (h0-h1) p0 (1/2 - p0/q0)

"""Referee: independent parse of seriesCfixedJ (CHMRS anc) and check leading-nu part vs bulk-point formula, d=3."""
import sys, re, sympy as sp, pickle, os
src = sys.argv[1]
raw = open(src).read()
i = raw.index('seriesCfixedJ ='); txt = raw[i+len('seriesCfixedJ ='):]
txt = txt.strip().rstrip(';')
txt = txt.replace('\\[CapitalDelta]\\[Phi]', 'Df').replace('\\[Nu]', 'nu').replace('^', '**')
assert '[' not in txt, txt[txt.index('['):txt.index('[')+80]
d,k,Df,nu,m2,J2,x,w = sp.symbols('d k Df nu m2 J2 x w')
e = sp.sympify(txt, locals=dict(d=d,k=k,Df=Df,nu=nu,m2=m2,J2=J2))
e = sp.expand(e.subs(d,3).subs(m2, 1/w))
# BP: (2-x^2)(1-x^2)^(-1-k/2) * P_J(1-2x^2); P_J(1-2x^2)=2F1(-J,J+1;1;x^2); write J via J2: J = (-1+sqrt(1+4J2))/2
# Independent route: generate P_J series coefficients from Legendre ODE recursion in t=x^2:
# t(1-t)y'' + (1-2t)y' + J2 y = 0  -> c_{n+1} = (n(n+1)-J2)/(n+1)^2 c_n
N=5
c=[sp.Integer(1)]
for n in range(N): c.append(sp.expand((n*(n+1)-J2)/sp.Integer((n+1)**2)*c[-1]))
t=sp.Symbol('t')
PJ=sum(c[n]*t**n for n in range(N+1))
Ek=sp.expand(sp.series((2-t)*(1-t)**(-1-k/2), t, 0, N+1).removeO())
BP=sp.expand(sp.series(sp.expand(Ek*PJ), t, 0, N).removeO())
for j in range(5):
    cj = sp.expand(e.coeff(w, j))
    lead = sp.Poly(cj, nu).coeff_monomial(nu**(2*j)) if j>0 else cj
    print(j, 'deg_nu', sp.Poly(cj,nu).degree(), ' lead - BP =', sp.simplify(sp.expand(lead - BP.coeff(t,j))))
# NLO at Df=3,k=2,J=0
R1=[sp.Poly(sp.expand(e.coeff(w,j)),nu).coeff_monomial(nu**(2*j-2)).subs({Df:3,k:2,J2:0}) for j in range(1,5)]
print('R1 (Df=3,k=2,J=0):', R1)
R3=[sp.Poly(sp.expand(e.coeff(w,j)),nu).coeff_monomial(nu**(2*j-6)) for j in range(3,5)]
pickle.dump({'e':sp.srepr(e)}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'cache','fixedJ.pkl'),'wb'))
for kk in (2,3):
  for JJ in (0,2,4):
    print('k',kk,'J',JJ,'R3 (x^0,x^2) =',[float(r.subs({Df:3,k:kk,J2:JJ*(JJ+1)})) for r in R3])

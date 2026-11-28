import sympy as S
T,nu,g,n=S.symbols('T nu g n',positive=True)
al=3*(S.sqrt(2)-1); A=al*nu
F2=1-6*nu/T+6*nu*(nu-1)/(T*(T-2))
P=S.expand(3*T*(T-2)*(T/3+S.Rational(3,2)+A*F2))
Pas=(T-(6-3*S.sqrt(2))*nu)**2*(T+(3+3*S.sqrt(2))*nu)
d=S.expand(S.cancel(P-Pas))
print('P - P_asym =', S.collect(S.simplify(d),T))
disc=S.expand(S.Poly(d,T).discriminant())
print('discriminant in T:', S.collect(S.simplify(disc),nu), ' at nu=4:', S.N(disc.subs(nu,4)))
print('roots of disc in nu:', [S.N(r) for r in S.solve(disc,nu)])
print('A<=W_B with g>=2n:', S.factor(S.expand((4*n**2+10*n-S.Rational(3,2))/3-al*n**2)), S.N((4-9*(S.sqrt(2)-1))))

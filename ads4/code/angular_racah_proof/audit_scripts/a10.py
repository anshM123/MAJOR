import sympy as S
T=S.Symbol('T')
a=S.Symbol('a',positive=True)
al=3*(S.sqrt(2)-1)
f=a**3+3*al*(a**2-6*a+6)   # 3*a^2*(a/3+al(1-6/a+6/a^2)) asymptotic numerator
print('asymptotic min:', [ (S.nsimplify(r), S.N(f.subs(a,r))) for r in S.solve(S.diff(f,a),a)])
print('factor f:', S.factor(S.expand(f),extension=S.sqrt(2)))
for n in [2,3,5,10,50,100,1000,10**4,10**6]:
    nu=n*n; A=al*nu
    M=T/3+S.Rational(3,2)+A*(1-6*nu/T+6*nu*(nu-1)/(T*(T-2)))
    num=S.numer(S.together(S.diff(M,T)))
    crit=[r for r in S.Poly(S.expand(num),T).nroots(n=40) if r.is_real and r>2]
    g_hi=2*n+2
    vals=[(S.N(M.subs(T,c)/(c/3+S.Rational(3,2)),20),S.N(c/nu,8)) for c in crit if c< g_hi*(g_hi+1)]
    print(n, vals)

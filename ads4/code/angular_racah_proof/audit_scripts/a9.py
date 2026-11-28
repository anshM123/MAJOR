import sympy as S, math
T,D=S.symbols('T Delta',positive=True)
def light(g):
    n=math.floor(g/2); nu=n*n
    WB=S.Rational(g**2+5*g,3)-S.Rational(1,2) if isinstance(g,int) else None
    F2=1-6*nu/T+6*nu*(nu-1)/(T*(T-2))
    M=T/3+S.Rational(3,2)+WB*F2      # margin as function of T, light range T in (2, g(g+1))  (Delta in (3,g+2))
    num=S.numer(S.together(M))
    crit=[r for r in S.Poly(S.numer(S.together(S.diff(M,T))),T).nroots() if r.is_real and 2<r<g*(g+1)]
    vals=[(float(M.subs(T,c)/ (c/3+1.5)),float(c)) for c in crit]
    return n,vals
for g in [4,6,10,20,40,60,80,90,100,150,200,400,1000]:
    n,vals=light(g)
    mn=min(vals) if vals else None
    if mn:
        Tc=mn[1]; Dl=1.5+math.sqrt(0.25+Tc)
        print(g,n,'min rel margin %.5f at T=%.1f (Delta=%.3f, Delta/g=%.3f)'%(mn[0],Tc,Dl,Dl/g))

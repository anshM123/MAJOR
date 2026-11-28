import sympy as S
n,Y,nu,h,H,u,t=S.symbols('n Y nu h H u t')
def yx(x,s=9):                      # finite 4F3, symbolic nu=n^2 and Y=N(N+9)
    tot,c=0,S.Integer(1)
    for k in range(1,x+1):
        i=k-1; c*=((i**2-nu)/S.Integer(k*k))*(x*(x+s)-i*(i+s))/(Y-i*(i+s)); tot+=c
    return S.cancel(1+tot)
y4,y5,y6=yx(4),yx(5),yx(6)
p=lambda j:S.Rational((j+1)*(j+9),2*j+10); q=lambda j:(2*j+9)*nu/(Y-j*(j+9))
E5=y5**2+h*p(4)*p(5)/q(5)*(y5-y4)*(y6-y5)             # Sonin energy from its definition
D=Y*(Y-10)*(Y-22)*(Y-36)*(Y-52); w=4*n**2+18*n+S.Rational(25,2)
F=S.cancel(((S.Rational(2,3)*((Y+S.Rational(25,2))/w)**2-E5.subs(h,H/15)).subs(nu,n**2))*1140*D**2*w**2)
for HH in (0,4):
    P=S.Poly(S.expand(F.subs(H,HH).subs(Y,4*n**2+30*n+36+t).subs(n,u+1)),u,t)
    print(HH,len(P.terms()),min(P.coeffs()))   # -> 169 terms, all >0, min 760

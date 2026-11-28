# Independent audit: recurrence, P_k identity, energies
import sympy as S
x,s,k,X,Y,nu,N = S.symbols('x s k X Y nu N')
def Pk(kk,xx,ss):  # product form
    return S.prod([(xx-i)*(xx+ss+i) for i in range(kk)])
def p(xx,ss): xx=S.sympify(xx); return (xx+1)*(xx+ss)/(2*xx+ss+1)
# 1a: P_k identity for explicit k, symbolic x,s
for kk in range(1,9):
    lhs = p(x,s)*(Pk(kk,x+1,s)-Pk(kk,x,s)) - p(x-1,s)*(Pk(kk,x,s)-Pk(kk,x-1,s))
    rhs = kk**2*(2*x+s)*Pk(kk-1,x,s)
    assert S.expand(S.cancel(lhs-rhs))==0, kk
print('P_k identity explicit k=1..8 symbolic s,x: OK')
# 1b: all-k via the ratios (symbolic k) -- confirm the ratios themselves with gamma functions
r_plus = (x+1)*(x+s+k)/((x-k+1)*(x+s))
r_minus = (x-k)*(x+s-1)/(x*(x+s+k-1))
r_low = 1/((x-k+1)*(x+s+k-1))
expr = p(x,s)*(r_plus-1) - p(x-1,s)*(1-r_minus) - k**2*(2*x+s)*r_low
print('all-k ratio identity:', S.simplify(S.cancel(expr)))
# ratios: P_k(X_x) = Gamma(x+1)/Gamma(x-k+1) * Gamma(x+s+k)/Gamma(x+s)
G=S.gamma
PkG = lambda xx: G(xx+1)/G(xx-k+1)*G(xx+s+k)/G(xx+s)
for name,num,den,val in [('plus',PkG(x+1),PkG(x),r_plus),('minus',PkG(x-1),PkG(x),r_minus)]:
    print(name, S.simplify(S.gammasimp(num/den)-val))
PkG1 = G(x+1)/G(x-k+2)*G(x+s+k-1)/G(x+s)
print('low', S.simplify(S.gammasimp(PkG1/PkG(x))-r_low))
# boundary x=0: p_{-1} -> backward term: P_k(X_{-1}) with X_{-1}=1-s ; P_k(X_0)=0 for k>=1
for kk in range(1,6):
    lhs = p(0,s)*(Pk(kk,1,s)-Pk(kk,0,s))
    rhs = kk**2*s*Pk(kk-1,0,s)
    assert S.expand(lhs-rhs)==0
print('x=0 boundary OK')
# 2: recurrence for y_x = 4F3 with symbolic nu (all n), symbolic Y, symbolic s
def yfun(xx,ss,nuv=nu,YY=Y):
    tot=0; c=S.Integer(1)
    for kk in range(0,xx+1):
        if kk>0:
            i=kk-1
            c = c*((i**2-nuv)/S.Integer(kk**2))*(xx*(xx+ss)-i*(i+ss))/(YY-i*(i+ss))
        tot+=c
    return tot
for ss in [1,9,S.Symbol('s')]:
    for xx in range(0,9):
        Xx=xx*(xx+ss)
        q=(2*S.Integer(xx)+ss)*nu/(Y-Xx)
        lhs = p(xx,ss)*(yfun(xx+1,ss)-yfun(xx,ss)) + q*yfun(xx,ss)
        if xx>0: lhs -= p(xx-1,ss)*(yfun(xx,ss)-yfun(xx-1,ss))
        assert S.cancel(S.together(lhs))==0,(ss,xx)
print('spin recurrence x=0..8 for s=1,9,symbolic s, symbolic nu,Y: OK')
# check sum to min(n,x) equals formal sum for integer n
for nn in [2,3,5]:
  for xx in range(0,9):
    full = yfun(xx,9,nuv=nn**2)
    trunc=0;c=S.Integer(1)
    for kk in range(0,min(nn,xx)+1):
        if kk>0:
            i=kk-1; c=c*((i-nn)*(nn+i)/S.Integer(kk**2))*(xx*(xx+9)-i*(i+9))/(Y-i*(i+9))
        trunc+=c
    assert S.cancel(full-trunc)==0
print('truncation consistency OK')
# check against source ladder definition: U_k^(b) with C,T and b=4 -> X,Y
J,D=S.symbols('J Delta')
C=J*(J+1); T=(D-1)*(D-2)
for b in [0,4]:
    Nn=D-b-2; ss=2*b+1
    print(b, S.expand(C-b*(b+1)-(J-b)*(J-b+ss)), S.expand(T-b*(b+1)-Nn*(Nn+ss)))
    for i in range(4):
        print('  ', S.expand((C-(b+i)*(b+i+1)) - ((J-b)*(J-b+ss)-i*(i+ss))), S.expand((T-(b+i)*(b+i+1))-(Nn*(Nn+ss)-i*(i+ss))))

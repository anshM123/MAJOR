import sympy as S
n,Y,nu,T,m,t=S.symbols('n Y nu T m t')
def yx(xx,ss=9,YY=Y):
    tot=0;c=S.Integer(1)
    for k in range(0,xx+1):
        if k>0:
            i=k-1; c=c*((i**2-nu)/S.Integer(k*k))*(xx*(xx+ss)-i*(i+ss))/(YY-i*(i+ss))
        tot+=c
    return S.cancel(tot)
yT=lambda xx: yx(xx,YY=T-20)
print('x1', S.cancel(1-yT(1)-10*nu/(T-20)))
print('x2', S.cancel(1-yT(2)-22*nu*(T-3*nu-27)/((T-20)*(T-30))))
B3=36*(T-30)*(T-42)-234*(nu-1)*(T-42)+364*(nu-1)*(nu-4)
print('x3', S.cancel(1-yT(3)-nu*B3/((T-20)*(T-30)*(T-42))))
T0=(2*n+6)*(2*n+5)
print('B3(T0)', S.expand(B3.subs(nu,n**2).subs(T,T0)-4*(n**4+297*n**3+4405*n**2-1089*n-338)))
print('dB3', S.expand(S.diff(B3,T).subs(nu,n**2).subs(T,T0)-18*(3*n**2+88*n-11)), S.diff(B3,T,2))
Q=(2*Y**3-21*Y**2*nu-115*Y**2+70*Y*nu**2+868*Y*nu+1806*Y-70*nu**3-1540*nu**2-7462*nu-6768)
D4=Y*(Y-10)*(Y-22)*(Y-36)
print('x4', S.cancel(1-yx(4)-26*nu*Q/D4))
# SOS
nn=m+2; N=2*nn+2+t
Qs=S.expand(Q.subs(nu,nn**2).subs(Y,N*(N+9)))
SOS=(2*m**4*(m-2*t-7)**2+m**4*(38*t**2+462*t+1578)+m**3*(152*t**3+3332*t**2+23036*t+52080)
 +m**2*(99*t**4+3486*t**3+40027*t**2+182448*t+281098)
 +m*(24*t**5+1176*t**4+20368*t**3+150612*t**2+480668*t+546588)
 +2*t**6+126*t**5+2987*t**4+32844*t**3+169559*t**2+402738*t+356184)
print('SOS identity', S.expand(Qs-SOS)==0)
neg=[(mon,c) for mon,c in S.Poly(Qs,m,t).terms() if c<0]
print('negative monomials in direct expansion', neg)
# domains: x2: T>=(2n+5)(2n+4); x3 ok; x4 N>=2n+2 (Delta>=8+g>=8+2n -> N=Delta-6>=2n+2)
print('x2 dom', S.expand((2*n+5)*(2*n+4)-3*n**2-27))
# B3 positivity for n>=1: 4405n^2-1089n-338 >= 2978 n^2 <=> 1427n^2-1089n-338>=0 at n=1: 
print(1427-1089-338)
# Remark rho_max
r=S.Symbol('r',positive=True)
print(S.solve(364*r**4-234*r**2+36,r), float(S.sqrt((117-S.sqrt(585))/364)))

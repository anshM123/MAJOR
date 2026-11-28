import sympy as S, json
n,Y,nu,h,H,u,t=S.symbols('n Y nu h H u t')
def yx(xx,ss=9):
    tot=0;c=S.Integer(1)
    for k in range(0,xx+1):
        if k>0:
            i=k-1; c=c*((i**2-nu)/S.Integer(k*k))*(xx*(xx+ss)-i*(i+ss))/(Y-i*(i+ss))
        tot+=c
    return S.cancel(tot)
y4,y5,y6=yx(4),yx(5),yx(6)
pf=lambda j:S.Rational((j+1)*(j+9),2*j+10)
qf=lambda j:(2*j+9)*nu/(Y-j*(j+9))
# E5 from definition, using actual y6 (no recurrence elimination)
E5def = y5**2 + h*pf(4)*pf(5)/qf(5)*(y5-y4)*(y6-y5)
E5dir = y5**2 + h*S.Rational(4225,6156)*(Y-70)/nu*(y5-y4)**2 - h*S.Rational(65,18)*y5*(y5-y4)
print('E5 def == direct:', S.cancel(E5def-E5dir)==0)
# paper's A,P
D=Y*(Y-10)*(Y-22)*(Y-36)*(Y-52)
A=S.cancel(D*y5); P=S.cancel(D*(y4-y5)/(18*nu))
print('A,P polynomial:', S.fraction(A)[1], S.fraction(P)[1])
Aex=(D-70*nu*(Y-10)*(Y-22)*(Y-36)*(Y-52)+1050*nu*(nu-1)*(Y-22)*(Y-36)*(Y-52)
     -5600*nu*(nu-1)*(nu-4)*(Y-36)*(Y-52)+11900*nu*(nu-1)*(nu-4)*(nu-9)*(Y-52)
     -8568*nu*(nu-1)*(nu-4)*(nu-9)*(nu-16))
Pex=(Y**4-28*Y**3*nu-92*Y**3+210*Y**2*nu**2+2030*Y**2*nu+2668*Y**2
     -560*Y*nu**3-10640*Y*nu**2-41664*Y*nu-26400*Y
     +476*nu**4+14840*nu**3+115388*nu**2+224112*nu+57024)
print('A matches', S.expand(A-Aex)==0, ' P matches', S.expand(P-Pex)==0)
B=4225*nu*(Y-70)*P**2+1235*nu*A*P
print('E5 rational form', S.cancel(E5dir-(19*A**2+h*B)/(19*D**2))==0)
# my own F_H: numerator of 2/3 ((Y+25/2)/(4n^2+18n+25/2))^2 - E5(H/15), computed from E5def directly
target=S.Rational(2,3)*((Y+S.Rational(25,2))/(4*n**2+18*n+S.Rational(25,2)))**2
diff=(target-E5def.subs(h,H/15)).subs(nu,n**2)
# multiply by positive factor 60*19*D^2*(4n^2+18n+25/2)^2
FH=S.cancel(diff*60*19*D**2*(4*n**2+18*n+S.Rational(25,2))**2)
num,den=S.fraction(FH)
print('denominator after scaling:',den)
FHpaper=190*(2*Y+25)**2*D**2-(8*n**2+36*n+25)**2*(285*A**2+H*B).subs(nu,n**2)
print('my F_H == paper F_H:', S.expand(FH-FHpaper)==0)
res={}
for HH in [0,4]:
    F=S.expand(FHpaper.subs(H,HH).subs(nu,n**2).subs(Y,4*n**2+30*n+36+t).subs(n,u+1))
    Pl=S.Poly(F,u,t)
    terms=Pl.terms()
    cs=[c for _,c in terms]
    sup=sorted(m for m,_ in terms)
    expected=sorted((a,b) for b in range(13) for a in range(25-2*b))
    print('H=',HH,'nterms',len(terms),'all>0',all(c>0 for c in cs),'min',min(cs),'support exact',sup==expected)
    rows=[]
    for b in range(13):
        row={a:c for (a,bb),c in terms if bb==b}
        m=min(row.values()); am=[a for a,c in row.items() if c==m]
        rows.append((b,m,am))
    print('  row minima',rows)
    res[HH]={f'{a},{b}':int(c) for (a,b),c in terms}
    # affine in H check
F_H=S.expand(FHpaper.subs(nu,n**2))
print('degree in H:', S.Poly(F_H,H).degree())
json.dump({str(k):v for k,v in res.items()},open('my_coeffs.json','w'))
# x=5 E5 : what about Y between 52 and 70 ... not needed. Also sanity: H range covers h5 in (0,4/15] -> H=15h in (0,4]

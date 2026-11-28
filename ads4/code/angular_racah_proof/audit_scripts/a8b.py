import mpmath as mp, random
from fractions import Fraction as Fr
mp.mp.dps=60
F=lambda f: mp.mpf(f.numerator)/f.denominator
random.seed(5)
hc={}
def h(j):
    if j not in hc: hc[j]=2*mp.quad(lambda t: t**j*(1+t**8)/(1+t),[0,1])
    return hc[j]
def y(x,n,N):
    Y=N*(N+9); X=x*(x+9); tot=Fr(1); c=Fr(1)
    for k in range(min(n,x)):
        c=c*Fr((k-n)*(n+k),(k+1)**2)*(X-k*(k+9))/(Y-k*(k+9)); tot+=c
    return tot
wm=-1e9;wc=1e9;wcert=1e9;info=None
for trial in range(150):
    g=Fr(random.randint(400,20000),100); n=int(g//2); x=random.randint(5,int(2*g)+30)
    N=x+g-2+random.choice([Fr(0),Fr(0),Fr(random.randint(0,300),100),Fr(random.randint(0,300*int(g)),100)])
    Y=N*(N+9); nu=n*n
    p=lambda j: Fr((j+1)*(j+9),2*j+10); q=lambda j:Fr((2*j+9)*nu)/(Y-j*(j+9))
    ys=[y(j,n,N) for j in range(0,x+2)]
    E=lambda j: F(ys[j]**2)+h(j)*F(p(j-1)*p(j)/q(j)*(ys[j]-ys[j-1])*(ys[j+1]-ys[j]))
    Es=[E(j) for j in range(5,x+1)]
    for i,(a,b) in enumerate(zip(Es,Es[1:])):
        if b-a>wm: wm=b-a; info=(float(g),n,x,float(N),i+5)
    for j in range(5,x+1): wc=min(wc, E(j)-mp.mpf(2)/3*F(ys[j]**2))
    wcert=min(wcert, F(Fr(2,3)*((Y+Fr(25,2))/(4*n*n+18*n+Fr(25,2)))**2)-Es[0])
print('max E_{j+1}-E_j:',mp.nstr(wm,5),info); print('min E_j-2/3y^2:',mp.nstr(wc,5)); print('min cert slack:',mp.nstr(wcert,5))

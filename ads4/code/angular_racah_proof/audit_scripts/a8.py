import mpmath as mp, random, math
mp.mp.dps=50
random.seed(5)
hc={}
def h(j):
    if j not in hc: hc[j]=2*mp.quad(lambda t: t**j*(1+t**8)/(1+t),[0,1])
    return hc[j]
def y(x,n,N):
    Y=N*(N+9); X=x*(x+9); tot=mp.mpf(1); c=mp.mpf(1)
    for k in range(min(n,x)):
        c=c*mp.mpf((k-n)*(n+k))/(k+1)**2*(X-k*(k+9))/(Y-k*(k+9)); tot+=c
    return tot
worst_mono=-1e9; worst_coer=1e9; worst_cert=1e9
for trial in range(150):
    g=mp.mpf(random.uniform(4,200)); n=int(mp.floor(g/2)); x=random.randint(5,int(2*g)+30)
    N=x+g-2+mp.mpf(random.choice([0,0,random.uniform(0,3),random.uniform(0,3*float(g))]))
    Y=N*(N+9); nu=n*n
    p=lambda j: mp.mpf((j+1)*(j+9))/(2*j+10); q=lambda j:(2*j+9)*nu/(Y-j*(j+9))
    ys=[y(j,n,N) for j in range(0,x+2)]
    E=lambda j: ys[j]**2+h(j)*p(j-1)*p(j)/q(j)*(ys[j]-ys[j-1])*(ys[j+1]-ys[j])
    Es=[E(j) for j in range(5,x+1)]
    for a,b in zip(Es,Es[1:]): worst_mono=max(worst_mono,(b-a)/abs(a))
    for j in range(5,x+1): worst_coer=min(worst_coer, E(j)-mp.mpf(2)/3*ys[j]**2)
    worst_cert=min(worst_cert, mp.mpf(2)/3*((Y+12.5)/(4*n*n+18*n+12.5))**2-Es[0])
print('max relative increase E_{j+1}-E_j (should be <=0):',mp.nstr(worst_mono,5))
print('min E_j-2/3 y_j^2 (>=0):',mp.nstr(worst_coer,5)); print('min cert slack:',mp.nstr(worst_cert,5))

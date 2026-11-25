import mpmath as mp, sys
sys.path.insert(0,'/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/kern')
def F3(nu,z):
    return (1-z)**-3*mp.legenp(-mp.mpf(1)/2+1j*nu,0,(1+z)/(1-z),type=3)
def R1_sum(J,D,b,nu,p=3):
    n=J-b; C=J*(J+1); T=(D-1)*(D-2); u0=b*(b+1); X=mp.mpf(T-u0); h=1/X; c=mp.mpf(C-u0)/X
    s=mp.mpf(nu)**2+mp.mpf(1)/4; K=int(200/(-mp.log(c)))+50
    g=[mp.mpf(1)]
    for k in range(K): g.append(((2*k*k+2*p*k+p-s)*g[k]-(k+p-1)**2*(g[k-1] if k else 0))/(k+1)**2)
    tot=0; ck=mp.mpf(1)
    for k in range(K+1):
        L1=(1/c-1)*((k-1)*k*(2*k-1)/mp.mpf(6)+(2*b+1)*k*(k-1)/mp.mpf(2))
        tot+=ck*mp.e**(-h*L1)*g[k]; ck*=c
    return tot,c,h
def R1_airy(J,D,b,nu):
    C=J*(J+1); u0=b*(b+1); T=(D-1)*(D-2); c=mp.mpf(C-u0)/(T-u0)
    eta=(1-c)/(C-u0); lam=mp.cbrt(eta); Db=b*b+b+mp.mpf(1)/3; rs=c*mp.e**(eta*Db)
    w=mp.e**(1j*mp.pi/6)
    f1=lambda u: mp.airyai(u)*mp.e**(-b*lam*u)*F3(nu,rs*mp.e**(-lam*u))
    f2=lambda u: 2*mp.re(w*mp.airyai(w*u)*mp.e**(1j*b*lam*u)*F3(nu,rs*mp.e**(1j*lam*u)))
    I1=mp.quad(f1,[0,2,5,10,20]); I2=mp.quad(f2,[0,2,5,10,20])
    A=mp.quad(lambda u: abs(f1(u))+abs(f2(u)),[0,2,5,10,20])
    return mp.e**(eta*b**3/3)*(I1+I2), A
for (J,D,b,nu,dps) in [(800,1800,0,5,30),(800,1800,0,50,40),(800,1800,4,50,40),(800,1800,0,300,80)]:
    mp.mp.dps=dps
    s,c,h=R1_sum(J,D,b,nu); mp.mp.dps=dps
    a,A=R1_airy(J,D,b,nu)
    print(J,D,b,nu,'sum',mp.nstr(s,12),'airy',mp.nstr(a,12),'rel diff',mp.nstr(abs(a-s)/abs(s),3),'amplification int|.|/|R1|',mp.nstr(A/abs(s),3),flush=True)

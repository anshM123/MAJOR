import mpmath as mp, sys
def test(J,Delta,b,nu,p=3,orders=(1,2)):
    mp.mp.dps=int(0.75*nu+120)
    n=J-b; C=J*(J+1); T=(Delta-1)*(Delta-2); u0=b*(b+1); X=mp.mpf(T-u0); h=1/X; c=mp.mpf(C-u0)/X
    s=mp.mpf(nu)**2+mp.mpf(1)/4
    Kmax=int(n+ 60/(-mp.log(c)) + 50) if c<1 else n+50
    Kmax=max(Kmax, n+50)
    g=[mp.mpf(1)]; gm=mp.mpf(0)
    for k in range(Kmax):
        nxt=((2*k*k+2*p*k+p-s)*g[k]-(k+p-1)**2*(g[k-1] if k>=1 else 0))/(k+1)**2
        g.append(nxt)
    # exact
    U=mp.mpf(1); Kex=mp.mpf(0)
    for k in range(n+1):
        Kex+=U*g[k]; u=(b+k)*(b+k+1); U*=mp.mpf(C-u)/(T-u)
    out={'exact':Kex}
    for M in orders:
        S=[mp.mpf(0)]*(M+1); tot=mp.mpf(0); ck=mp.mpf(1)
        for k in range(Kmax+1):
            expo=sum(h**l*(c**(-l)-1)/l*S[l] for l in range(1,M+1))
            tot+=ck*mp.e**(-expo)*g[k]
            v=k*(k+2*b+1)
            for l in range(1,M+1): S[l]+=mp.mpf(v)**l
            ck*=c
        out[M]=tot
    # plain Legendre at c
    out['c']=sum(c**k*g[k] for k in range(Kmax+1))
    return c,out
for (J,D,nu) in [(800,1800,50),(800,1800,300),(800,1800,700),(2000,3000,300),(2000,3000,700),(800,3800,700)]:
    c,o=test(J,D,0,nu)
    e=o['exact']
    print(f'J={J} D={D} c={float(c):.4f} nu={nu}: exact={mp.nstr(e,8)}  relerr: Legendre={mp.nstr(abs(o["c"]-e)/abs(e),3)}  exp(L1)={mp.nstr(abs(o[1]-e)/abs(e),3)}  exp(L1+L2)={mp.nstr(abs(o[2]-e)/abs(e),3)}', flush=True)

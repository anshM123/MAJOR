import mpmath as mp
mp.mp.dps=120
def m_seq(b,J,D):
    C=J*(J+1); T=(D-1)*(D-2); m=[mp.mpf(1)]
    for k in range(J-b):
        y=(b+k)*(b+k+1); m.append(m[-1]*(C-y)/(T-y))
    return m
def Ydist(m,M):
    Jb=len(m)-1; P=[]
    for y in range(Jb+1):
        s=mp.mpf(0)
        for k in range(y,Jb+1):
            s+=(-1)**(k-y)*mp.binomial(M-y,k-y)*m[k]
        P.append(mp.binomial(M,y)*s)
    return P
for J,h in [(80,10.5),(100,10.5),(100,12.3),(150,12.3),(150,20.7),(300,20.7),(600,30.5)]:
    for b in (0,4):
        m=m_seq(b,J,mp.mpf(J+h)); Jb=J-b
        res=[]
        for d in [0,1,2,5,10,20,50,100,200]:
            P=Ydist(m,Jb+d); res.append((d,mp.nstr(min(P),3)))
        print(J,b,h,res,flush=True)

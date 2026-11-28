from fractions import Fraction as Fr
import random, math, sys
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 7)
def R(b,J,D,n):
    C=J*(J+1); T=(D-1)*(D-2)
    tot=Fr(1); c=Fr(1)
    for k in range(min(n,J-b)):
        s=(b+k)*(b+k+1)
        c=c*Fr((k-n)*(n+k),(k+1)**2)*(C-s)/(T-s)
        tot+=c
    return tot
def W(b,D):
    T=(D-1)*(D-2); return (T-b*(b-1))/3+Fr(3,2)
worst={0:(9,None),4:(9,None)}
cnt=0
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 3000):
    g=Fr(random.randint(400,20000),100)
    if random.random()<0.2: g=Fr(random.choice([4,5,6,7,8,9,10,20,50,100,150,200]))-Fr(random.choice([0,0,1]),1000) if random.random()<.5 else Fr(random.choice([4,6,10,50,100,200]))
    if g<4: g=Fr(4)
    n=math.floor(g/2)
    WB=W(4,g+4)
    b=random.choice([0,4])
    r=random.random()
    if r<0.4: J=random.randint(0,12)
    elif r<0.8: J=random.randint(0,int(3*g)+10)
    else: J=random.randint(0,2000)
    if b==0: J=2*(J//2)
    if b==4 and J<4: J=4+J
    tt=random.random()
    if tt<0.4: D=J+g
    elif tt<0.8: D=J+g+Fr(random.randint(0,10**6),10**6)*random.choice([1,5,g,10*g])
    else: D=J+g+Fr(random.randint(0,10**6),10**3)
    Rv=R(b,J,D,n)
    m=(W(0,D)+WB*Rv) if b==0 else (W(4,D)-WB*Rv)
    rel=m/W(b,D)
    cnt+=1
    if rel<worst[b][0]: worst[b]=(rel,(float(g),n,J,float(D-J-g)))
    if (b==0 and m<=0) or (b==4 and m<0): print('VIOLATION',b,g,J,D,float(rel)); 
print(cnt,'points; worst relative margins:',{k:(float(v[0]),v[1]) for k,v in worst.items()})

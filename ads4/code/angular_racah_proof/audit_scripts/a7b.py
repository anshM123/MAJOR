from fractions import Fraction as Fr
import math
exec(open('a7.py').read().split('worst=')[0].split('random.seed')[0]+open('a7.py').read().split("random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 7)\n")[1].split('worst=')[0])
glist=[Fr(x,4) for x in range(16,81)]+[Fr(g) for g in range(21,201,7)]+[Fr(199999,1000),Fr(4001,20),Fr(7,1)-Fr(1,10**6)]
res={0:(9,None),4:(9,None)}
for g in glist:
    n=math.floor(g/2); WB=W(4,g+4)
    for b in (0,4):
        Js=list(range(0 if b==0 else 5, int(2*g)+40, 2 if b==0 else 1))+[int(5*g)+1, int(20*g)+3, 1001]
        for J in Js:
            if b==0 and J%2: J+=1
            for off in [Fr(0),Fr(1,3),Fr(1),g/3,g,4*g]:
                D=J+g+off
                Rv=R(b,J,D,n)
                m=(W(0,D)+WB*Rv) if b==0 else (W(4,D)-WB*Rv)
                rel=m/W(b,D)
                assert m>0, (b,g,J,D)
                if rel<res[b][0]: res[b]=(rel,(float(g),n,J,float(off)))
print({k:(float(v[0]),v[1]) for k,v in res.items()})

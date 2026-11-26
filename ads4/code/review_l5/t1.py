import sys, time
from an_common import *
ctx.prec = 4000
J,a,b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
st = state(J,a,b); U = ladder(st); c = st['c']
print('c', c.mid(), 'beta', float((2*(c.sqrt()).atanh()).mid()))
t0=time.time()
for nu in [1,5,20,100,300,500,700,900,950,990,1000,1020,1040]:
    K = Kladder(3, arb(nu), U)
    Fi = Finf(3, nu, c)
    print(nu, float(K.mid()), float(Fi.mid()), float((K-Fi).mid()), float(K.rad()))
print(time.time()-t0)

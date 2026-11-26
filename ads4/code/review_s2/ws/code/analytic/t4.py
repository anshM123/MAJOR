import sys, time
from an_common import *
from kplus import *
ctx.prec = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
J,a,b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
st = state(J,a,b); U = ladder(st)
print('c', float(st['c'].mid()), 'beta', float((2*(st['c'].sqrt()).atanh()).mid()))
# check recurrence vs direct
al = acb(2.5, 30)
d0 = S_direct(al-5, U)
t0=time.time()
for nu in [5, 100, 500, 990]:
    Kp, jj = Kplus(nu, st, U)
    K = Kladder(3, arb(nu), U)
    print(nu, jj, float((2*Kp.real).mid()), float(K.mid()), float(abs(2*Kp.real-K).upper()), f'{time.time()-t0:.1f}s')

"""lemU_cert/lu_explore.py [DV, exploratory]: cancellation factors of the Lemma U contraction.
Usage: lu_explore.py L Np prec nus(comma)"""
import sys, os, math, time, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from flint import arb, acb, ctx
import lu_core as C

L = float(sys.argv[1]); Np = int(sys.argv[2]); prec = int(sys.argv[3]); nus = [float(v) for v in sys.argv[4].split(',')]
ctx.prec = prec
M = 1000; K = 3*M + 3
c = arb((1 - arb(-L).exp()).mid()); d = 1 - c
t0 = time.time()
YA = C.integrated_jets_mid(c, K, Np, which=('A',))['A']
print('jets', time.time() - t0, flush=True)
vecs = [YA] + [C.jets(3, arb(nu), c, K, mid=True) for nu in nus]
s = 1; xs = C.xstar(c, s)
res = C.contract(c, xs, s, vecs, M, absnorm=True)
for i, (Am, Nm) in enumerate(res):
    pm = Am + [arb(0)]       # b = 0 dressing ignored here (diagnostic)
    Ph = C.cheb(pm)
    sabs = sum((abs(v) for v in Ph), arb(0))
    tot = sum(Nm, arb(0))
    lab = 'G' if i == 0 else 'nu=%g' % nus[i - 1]
    print(lab, 'P0=%.4e sum|Ph|=%.4e  sum_m N_m=%.3e  log10 kappa=%.1f  max_m log10 N_m=%.1f  max log10|A_m|=%.1f' % (
        float(Ph[0].mid()), float(sabs.mid()), float(tot.mid()), float((tot/sabs).log().mid())/math.log(10),
        max(float(v.log().mid()) for v in Nm if v > 0)/math.log(10), max(float(abs(v).log().mid()) for v in Am if v != 0)/math.log(10)), flush=True)
print('total', time.time() - t0)

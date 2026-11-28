import time, numpy as np, combo, lp2
g = 20; th = g ** 4
o = lp2.solve(g, 24, 0, th, Jlin=6*g, nJlog=60, nt=24, nlight=24)
print('lp2 baseline tau', o['tau'], 'c*', o['tau']/4)
t0 = time.time()
for c0 in [o['tau']/4*0.9, o['tau']/4*1.1]:
    v, st = combo.solve_c0(g, c0, 24, 0, 0, th, use_flip=False)
    print('c0', c0, 'Phi', v, st, time.time()-t0)

import sys, json, os
import numpy as np
from lp import Grid, solve, CACHE
G = Grid(4, 4, verbose=False)
out = []
for delta in [0.1, 0.2]:
    for BA in [120.0, 160.0, 200.0]:
        for Y in [2.0, 4.0]:
            st, v, info = solve(G, delta, mreg=dict(BA=BA, u1=0.01, Y=Y, amin=0.14, safety=1.4))
            if st == 0:
                print('delta=%.2f BA=%g Y=%g: A=%.4f B=%.4f B_AdS=%.4f a=%.4f gamma0=%.3f' % (delta, BA, Y, info['A'], info['B'], info['B_AdS'], info['a'], -2*(v[1]-3*v[0])), flush=True)
                out.append(dict(d2=4, d3=4, delta=delta, BA=BA, Y=Y, A=info['A'], B=info['B'], B_AdS=info['B_AdS'], a=info['a'], v=list(map(float, v))))
            else:
                print('delta=%.2f BA=%g Y=%g: status %d' % (delta, BA, Y, st), flush=True)
json.dump(out, open(os.path.join(CACHE, 'mtest.json'), 'w'))

"""agproof/ag_growth.py (AG_PROOF.md §3.2): growth rate of the exact scalar correction g = exact - flat (and of flat) along vertical
lines nu_r + i y, averaged over y in [5, 15]/kappa_loc(nu_r) (where the growing wave dominates by >= e^10), in units of the uniform rate
kappa_0 = 0.97 beta and of the local rate kappa_loc = kappa_0/sqrt(1 - x_r^2).  For a real-analytic two-wave function this growth equals
the decay rate of the + wave up to the amplitude factor |a(nu+iy) a(nu-iy)| (slowly varying).  Input: cache/exact_lines_<J>.pkl.
-> cache/growth.txt
"""
import os, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
L = float(np.sqrt(998*999))
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))
pr(__doc__.split('->')[0])
pr(' J   x1    | line p: growth(g)/kappa_0, growth(g)/kappa_loc, growth(flat)/kappa_loc')
for J in (40, 48, 64, 100):
    a = pickle.load(open(os.path.join(HERE, 'cache', 'exact_lines_%d.pkl' % J), 'rb'))
    st = a['st']; m = st['m']; k0 = 0.97*st['beta']
    row = []
    for p, l in sorted(a['lines'].items(), reverse=True):
        xr = p*L/m; kl = k0/np.sqrt(1 - xr*xr)
        ys = l['ys']; g = l['ex']/2 - l['flat']/2; fl = l['flat']/2
        y1, y2 = 5/kl, min(15/kl, ys[-1])
        if y2 <= y1:
            continue
        sg = (np.interp(y2, ys, np.log(np.abs(g))) - np.interp(y1, ys, np.log(np.abs(g))))/(y2 - y1)
        sf = (np.interp(y2, ys, np.log(np.abs(fl))) - np.interp(y1, ys, np.log(np.abs(fl))))/(y2 - y1)
        row.append(f'p={p:.2f}: {sg/k0:.2f}, {sg/kl:.2f}, {sf/kl:.2f}')
    pr(f" {J:3d} {L/m:.3f} | " + ' | '.join(row))
open(os.path.join(HERE, 'cache', 'growth.txt'), 'w').write('\n'.join(OUT) + '\n')

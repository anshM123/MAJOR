# Referee: double vs quad of the stripped tensor at moment-dominant points for the largest cases
import sys, os, time, numpy as np, mpmath as mp
sys.path.insert(0, sys.argv[1]); os.chdir(sys.argv[1])
sys.argv = sys.argv[:1] + sys.argv[2:]
import importlib
src = open('t2_spinning_blocks.py').read().split("y0 = [mp.mpf('0.5')")[0]
exec(src)
from moments import frame_data, frame_positions
from formulas import alpha_pair
l = int(sys.argv[1]); D = float(sys.argv[2]); Jv = int(sys.argv[3]); N = int(sys.argv[4])
al, alb = alpha_pair(D, Jv)
th = np.array([0.3, 1.1, 2.9]); zz = np.stack([np.ones(3), np.cos(th), np.sin(th)], axis=1)
sbd = SpinBlock(l, D, Jv, 1.75, N, prec='d'); t0 = time.time(); sbq = SpinBlock(l, D, Jv, 1.75, N, prec='q'); print("quad setup", time.time()-t0, flush=True)
for (u, t) in ((0.3, 6.0), (0.05, 6.0), (0.01, 10.0), (0.6, 3.0)):
    s = t/(2*(al + alb*u))
    fd = frame_data(np.array([u]), np.array([s]))
    x1, x2, x3, x4 = frame_positions(fd['ru'], fd['xi'])
    Wd = wfull(sbd, x1[0], x2[0], x3[0], x4[0], zz, zz)
    Wq = wfull(sbq, [mp.mpf(float(v)) for v in x1[0]], [mp.mpf(float(v)) for v in x2[0]], [mp.mpf(float(v)) for v in x3[0]], [mp.mpf(float(v)) for v in x4[0]],
               [[mp.mpf(float(c)) for c in row] for row in zz], [[mp.mpf(float(c)) for c in row] for row in zz])
    rel = max(abs(complex(Wd[i, j]/Wq[i, j]) - 1) for i in range(3) for j in range(3))
    print(f"l={l} (Delta,J)=({D},{Jv}) u={u} t={t} s={s:.4g}: max|double/quad-1| = {rel:.2e}", flush=True)

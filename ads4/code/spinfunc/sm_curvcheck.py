"""
sm_curvcheck (SMOOTH_BOUND.md §4): does a crystal row's Gamma contain the leading-Regge AdS curvature that the hybrid base already
contains (cv2 + cv3)?  For each row (Delta, J) of the q3/w2/ee1000 grids (unpatched), compare Gamma_{k,h}(x) with the two curvature
images of sm_curvfix.curv_image: 'bessel' (the hybrid's own 2(P - J0) term) and 'kappa' (exact-in-beta leading-Regge uniform action).
Per x band: rms Gamma, rms(Gamma - image), and the least-squares slope of Gamma on the image (slope ~ 1: the row contains it).
Usage: python3 sm_curvcheck.py "J:D,J:D,..."   (even J: 2x2 rows; odd J: +- only)  -> stdout (cache/sm_curvcheck.txt via tee)
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ['SM_CURV'] = 'none'
import numpy as np
import ee1000_common as E
import sm_curvfix as CF
Q = E.Q
BANDS = ((0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.81))

if __name__ == '__main__':
    for spec in sys.argv[1].split(','):
        J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
        g = Q.even_grid(J) if J % 2 == 0 else Q.odd_grid(J)
        if D not in g:
            print(f"J={J} D={D}: row missing; rows {sorted(g)}"); continue
        r = g[D]; xs, m = r['xs'], r['m']
        print(f"J={J} Delta={D:.0f} m={m:.1f} beta={2*np.arcsinh((J+.5)/m):.4f} nodes {len(xs)} (x {xs[0]:.3f}..{xs[-1]:.3f})")
        keys = [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-')] if J % 2 == 0 else [('C2', '+-'), ('C3H', '+-')]
        for kd, h in keys:
            G = np.asarray(r[(kd, h)] if J % 2 == 0 else r[kd], float)
            env, SIG = (Q.A.env, Q.A.SIG) if J % 2 == 0 else (Q.P4.env, Q.P4.SIG)
            out = []
            for img in ('bessel', 'kappa'):
                CF.IMG = img
                c = CF.curv_image(kd, h, J, m, xs, env, SIG)
                for a, b in BANDS:
                    s = (xs >= a) & (xs < b)
                    if s.sum() < 2:
                        continue
                    sl = float(G[s] @ c[s]/(c[s] @ c[s]))
                    out.append(f"{img[0]}[{a:.1f},{b:.1f}) slope {sl:+.2f} rms G {np.sqrt(np.mean(G[s]**2)):5.0f} img {np.sqrt(np.mean(c[s]**2)):5.0f} res {np.sqrt(np.mean((G[s]-c[s])**2)):5.0f}")
            print(f"   {kd:3s} {h}: " + "\n            ".join(out))

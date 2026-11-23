"""
pwimg/pw1b_e2e.py (PW_IMAGES.md §2.3): end-to-end check, refined.  Same as pw1 (c) at the worst cases (nu = 6, both polarisation
sets) with the direct quadrature refined (nth 128 -> 256, ds 0.02 -> 0.01) and the block series N = 140 -> 200, to show that the
residual of pw1 (c) is the quadrature/series error of the representation K = N r^2 E^1 (EESTEP §1.4), not of the de-shift identity;
plus the identity r^2 E^1_nu = (r/4)/f [alpha E^0_{nu+i} + beta E^0_{nu-i}] itself at the same points (block sums only).
Output: cache/pw1b_e2e.txt
"""
import os, sys, io
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pw_common as PC
OUT = io.StringIO()
pr = lambda s: (print(s, flush=True), OUT.write(s + '\n'))
pr(__doc__.split('Output')[0])
rng = np.random.default_rng(23)
pts = [(0.25, 0.6), (0.55, 1.4), (0.85, 2.2), (0.4, 0.3)]
for ip in range(2):
    pa = tuple(rng.uniform(0, 2*np.pi, 4))
    for nu in (2.5, 6.0):
        for (nth, ds, NB) in ((128, 0.02, 140), (256, 0.01, 200)):
            rat, idn = [], []
            for (rho, beta) in pts:
                r = rho*np.exp(-beta); pol = PC.Q.pols(beta, *pa)
                parts = PC.Q.kernel_parts(r, beta, nth=nth, ds=ds)
                K = PC.Q.integrate(parts, [nu], PC.Q.dressing(parts['th'], r, beta, pol))[0]
                f = np.cosh(beta)*(1 + r*r)/(2*r)
                E1, _ = PC.E_val(nu, r, beta, pol, N=NB, Dp=1.0)
                Ep, _ = PC.E_val(nu + 1j, r, beta, pol, N=NB); Em, _ = PC.E_val(nu - 1j, r, beta, pol, N=NB)
                pred = (r/4)/f*((1 + 0.5j/nu)*Ep + (1 - 0.5j/nu)*Em)
                rat.append(K/pred); idn.append(abs(r*r*E1/pred - 1))
            rat = np.array(rat)
            pr(f'pol set {ip}, nu = {nu}, (nth, ds, N) = ({nth}, {ds}, {NB}): K/pred spread {np.max(np.abs(rat/rat[0] - 1)):.1e};'
               f'  |r^2 E1/pred - 1| max {max(idn):.1e}')
open(os.path.join(PC.CACHE, 'pw1b_e2e.txt'), 'w').write(OUT.getvalue())

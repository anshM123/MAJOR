"""dqsum: print D' quadrature results (pcyc_dquad caches) side by side: relative change of I(nu) vs a base grid."""
import os, sys, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
D = sys.argv[1]; K = sys.argv[2]; base = sys.argv[3]; SUB = sys.argv[4] if len(sys.argv) > 4 else ""
def ld(g):
    d = np.load(os.path.join(HERE, 'cache', SUB, f'dq_D{D}_{K}_{g}.npz')); return d['nus'], d['I']*np.exp(float(d['shift']))
nus, I0 = ld(base)
print(f'# D={D} {K}: rel. change vs {base} at nu = {list(nus)}')
for f in sorted(glob.glob(os.path.join(HERE, 'cache', SUB, f'dq_D{D}_{K}_*.npz'))):
    g = os.path.basename(f)[len(f'dq_D{D}_{K}_'):-4]
    n, I = ld(g)
    print(f'  {g:24s} ' + ' '.join(f'{abs(I[i]/I0[i] - 1):.1e}' for i in range(len(nus))))

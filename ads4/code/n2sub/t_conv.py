"""t_conv.py: checks of n2_common against the project's own objects (read-only imports).
 (1) kappa image vs pw_gamma.pkl 'img' at (1100,100) nodes; (2) kernels.json R_-+ vs u2_hclosed_struct EE (the chain's u_EE)."""
import os, sys, pickle
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N
sys.path.insert(0, os.path.join(N.ADS4, 'code', 'spinfunc'))
d = pickle.load(open(os.path.join(N.ADS4, 'code', 'spinfunc', 'cache', 'pw_gamma.pkl'), 'rb'))
r = d[(100, 1100.0, 12.0, 34.0, 'cv', ('_w',))]
m, nus = r['m'], r['nus']
for h in ('++', '+-'):
    k = N.kappa('C2', h, 100, m, nus/m).real
    print(h, 'kappa mine', np.round(k, 1)); print(h, 'img stored', np.round(r[('C2', h)]['img'], 1))
import u2_hclosed_struct as U2
F = U2.numfuncs()['EE']
Rm, Rp = N.Rk()
v = np.array([0.3, 1.17, 3.0, 6.6, 15.0, 40.0, 2+0.7j, 5-0.4j])
print('R- rel diff', np.max(np.abs(F[0](v)/Rm(v) - 1)), ' R+ rel diff', np.max(np.abs(F[1](v)/Rp(v) - 1)))
# real zeros of R_- + R_+
nv = np.linspace(0.01, 30, 300001); s = (Rm(nv) + Rp(nv))
print('R-+R+ real?', np.max(np.abs(s.imag)/np.abs(s)))
z = nv[:-1][np.sign(s.real[:-1]) != np.sign(s.real[1:])]
print('real zeros of R_- + R_+ on (0,30):', z)

"""md_probe3.py -- boundary exponents on D' of the k = 3 (C3H) and k = 2 (E2) integrands relative to the scalar, J = 0 and 2."""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'pcyc')); sys.path.insert(0, os.path.join(HERE, '..', 'rd'))
import numpy as np
import pcyc_dlib as P
import md_spinJ as SJ
np.seterr(all='ignore')
Delta = float(sys.argv[1]); J = int(sys.argv[2])
base = dict(tx=1.0, rx=0.4, s=0.3, ty=1.2, ry=0.55, ph=0.3, al=1.1)
key = {'tx': 'tx', 'rx0': 'rx', 'rx1': 'rx', 's1': 's', 'ty': 'ty', 'ry0': 'ry', 'ry1': 'ry', 'ph': 'ph', 'al0': 'al', 'alpi': 'al'}
lim = {'tx': 1e-4, 'rx0': 1e-4, 'rx1': 1 - 1e-4, 's1': 1 - 1e-4, 'ty': 1e-4, 'ry0': 1e-4, 'ry1': 1 - 1e-4, 'ph': np.pi/2 - 1e-4,
       'al0': 1e-4, 'alpi': np.pi - 1e-4}
def ev(pt):
    o = P.point(*[np.array([pt[k]]) for k in ('tx', 'rx', 's', 'ty', 'ry', 'ph', 'al')], Delta)
    return o['L'][0], np.log(np.max(np.abs(SJ.ratio(o['g'], Delta, J)))), np.log(np.max(np.abs(SJ.ratio3(o['g'], Delta, J))))
for nm, v in lim.items():
    rows = []
    for f in (1.0, 10.0):
        pt = dict(base); vv = v*f if nm in ('tx', 'rx0', 'ty', 'ry0', 'al0') else (1 - (1 - v)*f if nm in ('rx1', 's1', 'ry1') else
             (np.pi/2 - (np.pi/2 - v)*f if nm == 'ph' else np.pi - (np.pi - v)*f))
        pt[key[nm]] = vv
        rows.append(ev(pt))
    sl = (np.array(rows[0]) - np.array(rows[1]))/np.log(0.1)
    print(f'J={J} {nm:5s}: scalar {sl[0]:+.2f}  E2/scal {sl[1]:+.2f}  C3H/scal {sl[2]:+.2f}', flush=True)

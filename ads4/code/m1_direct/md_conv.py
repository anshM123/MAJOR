"""md_conv.py -- grid-convergence table of the direct D' runs: max_nu |I_g - I_ref|/rms(I_ref) on nu <= numax, per Delta."""
import os, sys, glob
import numpy as np
numax = float(sys.argv[1]) if len(sys.argv) > 1 else 6.0
pat = sys.argv[2] if len(sys.argv) > 2 else 'dq_D*_E2_*.npz'
runs = {}
for f in glob.glob(os.path.join('out', pat)):
    d = np.load(f)
    key = float(d['Delta']) if 'J' not in d.files else (float(d['Delta']), int(d['J']))
    I = d['I']*np.exp(float(d['shift']))
    runs.setdefault(key, []).append(('-'.join(map(str, d['ns'])), d['nus'], I))
for key in sorted(runs):
    rs = sorted(runs[key], key=lambda t: np.prod([int(v) for v in t[0].split('-')]))
    g0, nu0, I0 = rs[-1]
    m = nu0 <= numax
    sh = (lambda I: I.reshape(-1, I.shape[-1])) if I0.ndim > 1 else (lambda I: I[None])
    rms = np.sqrt(np.mean(np.abs(sh(I0)[:, m])**2, axis=1))
    out = []
    for g, nu, I in rs[:-1]:
        mm = np.isin(np.round(nu, 6), np.round(nu0[m], 6))
        d_ = np.max(np.abs(sh(I)[:, mm] - sh(I0)[:, m])/rms[:, None])
        out.append(f'{g}: {d_:.1e}')
    print(f'{key}: ref {g0};  ' + ';  '.join(out))

"""dc_scalseeds (DIRCV.md §3.1): per-wave SCALAR deviation from the exact ratio at (1100,100), all independent seeds pooled."""
import dc_guard, os, sys, pickle
import numpy as np
import z3_nodes as Z3, fw_scalar as FS, pw_gamma as PG
J, D = 100, 1100.0
FS.install(J, D)
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
res = {}
for t, base in (('_w', 'SPIN'), ('_s2', 'SPIN'), ('_v', 'SPIN'), ('_x1', 'CACHE'), ('_x2', 'CACHE')):
    d = pickle.load(open(os.path.join(dc_guard.SPIN + '/cache' if base == 'SPIN' else dc_guard.CACHE, 'pw_pts_J%d_D%g%s.pkl' % (J, D, t)), 'rb'))
    for (DD, nu), r in d.items():
        if DD != D or (Z3.DREF, nu) not in d or r['modes'] != [1, -1]:
            continue
        o = Z3.scalar_check(J, nu, {D: PG.to_z2(r, 'EE') if 'EE' in r['g'] else to_s(r), Z3.DREF: PG.to_z2(d[(Z3.DREF, nu)], 'EE') if 'EE' in d[(Z3.DREF, nu)]['g'] else to_s(d[(Z3.DREF, nu)])}, D) if False else None
        def to_s(rec):
            pm = rec['per_mode']; ref = pm[0]['logN']
            reps = sum(np.asarray(x['reps'])*np.exp(x['logN'] - ref) for x in pm)
            return dict(repsS=reps, logNs=ref, method='pw')
        o = Z3.scalar_check(J, nu, {D: to_s(r), Z3.DREF: to_s(d[(Z3.DREF, nu)])}, D)
        res.setdefault(nu, []).append((t, o['d'].real, o['err'].real, o['dc2'], o['ec2']))
for nu in sorted(res):
    v = np.array([x[1] for x in res[nu]]); e = np.array([x[2] for x in res[nu]])
    w = 1/e**2; mu = np.sum(w*v)/np.sum(w); em = 1/np.sqrt(np.sum(w)); chi = np.sum(((v - mu)/e)**2)
    fac = res[nu][0][3]/res[nu][0][1]
    pr(f"nu={nu:5.1f}: " + "; ".join(f"{t} {x:+.5f}+-{ee:.5f}" for t, x, ee, _, _ in res[nu]) +
       f"  ->  pooled d = {mu:+.5f} +- {em:.5f} ({mu/em:+.1f} sigma; chi2 {chi:.1f}/{len(v)-1});  scalar Dc2 = {fac*mu:+.0f} +- {abs(fac)*em:.0f}")
open(os.path.join(dc_guard.CACHE, 'dc_scalseeds.txt'), 'w').write('\n'.join(lines) + '\n')

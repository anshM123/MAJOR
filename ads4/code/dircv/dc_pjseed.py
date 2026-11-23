"""dc_pjseed (DIRCV.md §7): seed-by-seed consistency of the P_J-cv state quantity at nu = 13.3:
q = log(G_P/S_P) at the state (EE helicity trace; what P_J-cv uses beyond the exact mean), and per-replica outlier scan."""
import dc_guard, os, sys, pickle
import numpy as np
import pw_gamma as PG
J, D, nu = 100, 1100.0, float(sys.argv[1]) if len(sys.argv) > 1 else 13.3
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
tr = lambda M: M[..., 0, 0] + M[..., 1, 1]
vals = []
for t in ('_w', '_x1', '_x2', '_y1', '_y2'):
    for base in (os.path.join(dc_guard.SPIN, 'cache'), dc_guard.CACHE):
        f = os.path.join(base, 'pw_pts_J%d_D%g%sP.pkl' % (J, D, t))
        if os.path.exists(f):
            break
    d = pickle.load(open(f, 'rb'))
    for kind in ('EE', 'C3H'):
        z = PG.to_z2(d[(D, nu)], kind); zr = PG.to_z2(d[(20000.0, nu)], kind)
        G = tr(z['reps']); S = z['repsS']; R = len(S)
        q = lambda s: (np.log(G[s].mean()/S[s].mean()) + z['logN'] - z['logNs']).real
        jk = np.array([q(np.arange(R) != i) for i in range(R)]); e = np.sqrt((R - 1)/R*np.sum((jk - jk.mean())**2))
        # per-replica influence: leave-one-out shift in units of the jackknife error
        infl = (jk - jk.mean())*(R - 1)/e
        # reference graviton/scalar
        qr = (np.log(tr(zr['reps']).mean()/zr['repsS'].mean()) + zr['logN'] - zr['logNs']).real
        diag = d[(D, nu)]['per_mode'][0]['diag']
        pr(f"{t:4s} {kind:3s} R={R:3d}: q_state = {q(np.ones(R, bool)):+.5f} +- {e:.5f};  max |replica influence| = {np.max(np.abs(infl)):.1f} sigma (rep {int(np.argmax(np.abs(infl)))});"
           f"  q_ref = {qr:+.6f};  diag nbad {diag['nbad']} wind {diag['wind']:.2f}")
        if kind == 'EE':
            vals.append((t, q(np.ones(R, bool)), e))
v = np.array([x[1] for x in vals]); e = np.array([x[2] for x in vals]); w = 1/e**2
mu = np.sum(w*v)/np.sum(w)
pr(f"EE: weighted mean {mu:+.5f} +- {1/np.sqrt(np.sum(w)):.5f}; chi2 = {np.sum(((v - mu)/e)**2):.1f}/{len(v) - 1};  pulls " + ", ".join(f"{t} {(x - mu)/ee:+.1f}" for t, x, ee in vals))
open(os.path.join(dc_guard.CACHE, 'dc_pjseed.txt'), 'w').write('\n'.join(lines) + '\n')

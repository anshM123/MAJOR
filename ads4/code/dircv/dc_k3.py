"""
dc_k3 (DIRCV.md §3.5): the k = 3 cross-check (C3H shares a_h, no de-shift: its D_h coefficient must be sigma_h c0(C2)) with
JACKKNIFE errors of the difference (paired on the same replicas), for every estimator / record set.  C2: F1 (A-columns, kappa
de-shifted); C3H: stripped columns {D, x D'/(2J+1), 1} (n2_refit convention).
"""
import dc_guard, os, sys
import numpy as np
import dc_pipe as DP, dc_fit as F
import n2_common as N
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
J, D = 100, 1100.0


def c3fit(o, h, nus_sel):
    idx = F.select(o, nus_sel); m = o['m']; nu = np.asarray(o['nus'])[idx]; x = nu/m
    Dv = N.Dsh(h, J, x).real; Dp = (F.dprime(h, J, x)*x/(2*J + 1)).real
    A = np.array([Dv, Dp, np.ones_like(x)]).T
    y = o[('C3H', h)]['y'][idx]; yj = o[('C3H', h)]['yj'][:, idx]
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    cj = np.array([np.linalg.lstsq(A, yy, rcond=None)[0] for yy in yj])
    return c, cj


if __name__ == '__main__':
    for (est, tag, kf, nus) in (('dir', '_w', False, [13.3, 16.3, 20.8, 25.0]), ('cv', '_w', False, [13.3, 16.3, 20.8, 25.0]),
                                ('cv', '_w', True, [13.3, 16.3, 16.7, 20.8, 25.0, 33.4]), ('dir', '_w', True, [13.3, 16.3, 16.7, 20.8, 25.0, 33.4]),
                                ('cv', '_wP', False, [13.3, 16.3, 20.8, 25.0]), ('dir', '_wP', False, [13.3, 16.3, 20.8, 25.0]),
                                ('cv', '_none', True, [13.3, 16.7, 20.8, 25.0, 33.4]), ('dir', '_none', True, [13.3, 16.7, 20.8, 25.0, 33.4])):
        o = DP.get(J, D, est, [tag], kf)
        s = []
        for h, sg in (('++', 1.0), ('+-', -1.0)):
            f2 = F.fit(o, h, nus, 'F1')
            c3, c3j = c3fit(o, h, nus)
            R = c3j.shape[0]
            dj = c3j[:, 0] - sg*f2['cj'][:, 0]
            e = np.sqrt((R - 1)/R*np.sum((dj - dj.mean())**2))
            e3 = np.sqrt((R - 1)/R*np.sum((c3j[:, 0] - c3j[:, 0].mean())**2))
            d = c3[0] - sg*f2['c'][0]
            s.append(f"{h}: c0_3 = {c3[0]:+6.0f} +- {e3:5.0f} vs sigma c0_2 = {sg*f2['c'][0]:+6.0f};  diff {d:+6.0f} +- {e:5.0f} ({abs(d)/e:.1f} sig)")
        pr(f"  {est:3s} {tag:5s} keepfw={int(kf)} nodes {nus}: " + " | ".join(s))
    open(os.path.join(dc_guard.CACHE, 'dc_k3.txt'), 'w').write('\n'.join(lines) + '\n')

"""LF2: dense float check of an LP functional (cache/lp_*.json) on u' in (0, a^2], all J up to a large cap, both channels;
plus the large-bhat asymptotic coefficient and the scaling-limit B2.  [NE] (floating point, factored weights)."""
import sys, os, json
import numpy as np
from ln_common import guard, CACHE
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'robust'))
from rcommon import flat_integrals, nodes, SIGMA, scaling_B2
from nf_lp import basis_values, basis

def check(fn, Jcap=6000, nu=160):
    I = json.load(open(fn)); a, r, d = I['a'], I['r'], I['d']; v = np.array(I['v']); n = d + 1
    worst = (9, None); a2 = a * a
    us = np.unique(np.concatenate([np.logspace(-6, np.log10(a2), nu), np.linspace(0.9 * a2, a2, 21)]))
    per_u = []
    for u in us:
        Jm = int(min(Jcap, max(60, (600 / np.sqrt(u) - 1) / 2)))
        nq = Jm + 16; p, w = nodes(nq)
        V = basis_values(p, d, r)
        psi2 = v[:n] @ V; psi3 = v[n:] @ V
        for h in SIGMA:
            W = np.vstack([w * psi2 * (2 - u * p * p), w * SIGMA[h] * u * psi3]).T
            Js, R = flat_integrals(h, u, Jm, W, p)
            if h == '++':
                sel = Js % 2 == 0; Js, R = Js[sel], R[sel]
            B2, B3 = R[:, 0], R[:, 1]; N = np.abs(B2) + np.abs(B3); rat = (B2 + B3) / N
            i = np.argmin(rat); per_u.append((u, h, float(rat[i]), int(Js[i])))
            if rat[i] < worst[0]: worst = (float(rat[i]), (float(u), h, int(Js[i]), int(Jm)))
    # scaling limit and 1/bhat^3 coefficient
    b = basis(d, r); bsc = np.linspace(0, 2000, 40001)
    S = scaling_B2(b, bsc) @ v[:n]
    pp2 = sum(v[i] * 2 * b[i][2] for i in range(n)); pp3 = sum(v[n + i] * 2 * b[i][2] for i in range(n))
    asym = (-pp2 - a2 * abs(pp3) / 2)
    return dict(file=os.path.basename(fn), a=a, delta=I['delta'], tau=I['tau'], worst=worst, minScaling=float(S.min()),
                argminScaling=float(bsc[np.argmin(S)]), psi2pp0=float(pp2), psi3pp0=float(pp3), asym_coef=float(asym),
                worst_u_list=sorted(per_u, key=lambda x: x[2])[:6])

if __name__ == '__main__':
    out = [check(f) for f in sys.argv[1:]]
    for o in out:
        print(json.dumps(o))
    json.dump(out, open(guard(os.path.join(CACHE, 'ln3_check_%s.json' % os.path.basename(sys.argv[1])[:-5])), 'w'))

"""
edger2/e2_merge.py (EDGE_R2.md §2): the batch-3 partial external results (ads4/ext/results/ext_results_batch3_partial.json, read-only)
decoded with ext/unpack.dec (imported, read-only; no external pickle is loaded) and merged per state into hc_points cache files
edger2/cache/pts/hc_pts_<kind>_J<J>_D<D>.pkl ({(Delta, nu): record}, state and Delta_ref = 2e4 reference), plus a node inventory
against ext/jobs_batch3.txt (which nodes of each state's job list are present, n1 vs n3).
Usage: python3 e2_merge.py   -> cache/pts/*.pkl, cache/merge.txt
"""
import e2_guard  # noqa: F401
import os, sys, json, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.abspath(os.path.join(HERE, '..', '..', 'ext'))
sys.path.insert(0, EXT)
from unpack import dec
OUT = os.path.join(HERE, 'cache', 'pts'); os.makedirs(OUT, exist_ok=True)
TAGS = ['P1_J26_m1.026', 'P1_J30_m1.030', 'P1_J34_m1.033', 'P2_J33_m1.033', 'P2_J38_m1.037', 'P2_J42_m1.041', 'P2_J46_m1.045',
        'P2_J50_m1.049']


def jobs():
    out = {}
    for ln in open(os.path.join(EXT, 'jobs_batch3.txt')):
        if ln.startswith('#') or not ln.strip():
            continue
        t, meth, J, D, seed, nus = ln.split()
        out[t] = dict(method=meth, J=int(J), D=float(D), nus=[float(v) for v in nus.split(',')])
    return out


def main():
    lines = []
    data = json.load(open(os.path.join(EXT, 'results', 'ext_results_batch3_partial.json')))
    merged = {}
    for r in data:
        tag = r['dir'].split('/')[1]
        if not tag.startswith(('P1_', 'P2_')):
            continue
        merged.setdefault(r['file'], {}).update(dec(r['data']))
    for f, d in merged.items():
        pickle.dump(d, open(os.path.join(OUT, f), 'wb'))
    jb = jobs()
    for t in TAGS:
        j1 = jb[t + '_n1']; j3 = jb.get(t + '_n3a', dict(nus=[]))
        J, D = j1['J'], j1['D']
        d = merged.get('hc_pts_EE_J%d_D%g.pkl' % (J, D), {})
        c = merged.get('hc_pts_C3H_J%d_D%g.pkl' % (J, D), {})
        have = sorted(nu for (DD, nu) in d if DD == D and (20000.0, nu) in d and (D, nu) in c and (20000.0, nu) in c)
        miss1 = [nu for nu in j1['nus'] if not any(abs(nu - h) < 1e-5 for h in have)]
        al = D - J - 2; m = np.sqrt(al*(al + 2*J + 1))
        n3only = [nu for nu in j3['nus'] if not any(abs(nu - v) < 1e-5 for v in j1['nus'])]
        lines.append(f"{t}: J = {J}, Delta = {D:g}, m = {m:.4f}; n1 nodes present {len(have)}/{len(j1['nus'])} (missing {miss1}); "
                     f"x range {have[0]/m:.4f} .. {have[-1]/m:.4f} (nu/Delta <= {have[-1]/D:.4f}); n3-only nodes NOT run: {len(n3only)} "
                     f"(x = {', '.join(f'{v/m:.4f}' for v in sorted(n3only))})")
    open(os.path.join(HERE, 'cache', 'merge.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()

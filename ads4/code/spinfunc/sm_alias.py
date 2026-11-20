"""
sm_alias (SMOOTH_BOUND.md §3-4): representation (aliasing) error of a chain on a KNOWN function with the oscillation of Gamma.
Input: a sm_chainfix run with SM_CURV=only (every row's diagonal Gamma replaced by its own curvature image), i.e. per hybrid state the
chain's dF_hh = model integral of the curvature.  Compared with the exact state curvature cv2 + cv3 (ee1000_common.hyb, same p-rule as the
chain: 64 nodes for w2/q3, EE1000_NPN for ee1000).  err_h = (dF_hh - cv_h)/N_h is the chain's representation error for an oscillating
correction of the curvature's size; it is also the difference between the two double-count fixes ('row' vs 'state').
Usage: python3 sm_alias.py w2 <tag> [npn]  |  python3 sm_alias.py ee1000 <tag> [npn]
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import ee1000_common as E

mode, tag = sys.argv[1], sys.argv[2]
npn = int(sys.argv[3]) if len(sys.argv) > 3 else (64 if mode == 'w2' else 256)
d = pickle.load(open(os.path.join(E.C, 'sm_chainfix_%s_%s.pkl' % (mode, tag)), 'rb'))
items = []
if mode == 'w2':
    for Dg, recs in d.items():
        for r in recs:
            if r['DF'] == 'hyb' and r['J'] >= 16:
                items.append((Dg, r))
else:
    Dg = int(os.environ.get('EE1000_DG', '1000'))
    items = [(Dg, r) for r in d['recs']]
print(f"# sm_alias {mode} {tag} (npn {npn}): per hybrid state, chain integral of the curvature image vs exact cv, in N_h units")
rows = []
for Dg, r in items:
    J, mfl = r['J'], r['mfl']
    errs = []
    for i, h in enumerate(r['hs']):
        st = E.hyb(Dg, J, mfl, h, npn)
        cv = st['cv2'] + st['cv3']
        dFh = r['dF'][i][i] if np.ndim(r['dF']) == 2 else r['diag'][i]*r['Ns'][i]
        errs.append(((dFh - cv)/st['Ns'], cv/st['Ns']))
    rows.append((Dg, J, mfl, errs))
for Dg in sorted(set(x[0] for x in rows)):
    for J in sorted(set(x[1] for x in rows if x[0] == Dg)):
        R = [x for x in rows if x[0] == Dg and x[1] == J]
        s = "  ".join(f"{mfl:.3f}: " + "/".join(f"{e:+.3f}({c:+.3f})" for e, c in errs) for _, _, mfl, errs in R if mfl <= 1.5)
        print(f"Dg={Dg} J={J:3d}  m/L: (model - cv)/N_h (cv/N_h) per helicity:  {s}")

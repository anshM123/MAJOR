"""Referee: check the CHL closed forms (t1b) at new exact points, incl. non-integer J and reflected labels, different random kinematics."""
import sys, os, pickle, sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fourpt'))
import chl_map
from t1b_helicity_l1 import closed
from embed2 import CACHE
ws, vecs = pickle.load(open(os.path.join(CACHE, 'consbasis_l1.pkl'), 'rb'))
Qs, BV, syms = chl_map.frame_Q()
orig_points = chl_map.points
chl_map.points = lambda s, n=7, seed=11: orig_points(s, n, seed=97)
R = sp.Rational
for (D, J) in [(R(171, 13), 7), (R(29, 3), R(5, 2)), (3 - R(171, 13), -8), (R(83, 7), R(13, 3)), (3 - R(83, 7), -1 - R(13, 3))]:
    try:
        r = chl_map.identify(D, J, Qs, BV, syms, ws, vecs)
    except AssertionError as e:
        print(D, J, 'degenerate', e); continue
    c = closed(D, J)
    ok = all(sp.simplify(r[h][0] - c[h][0]) == 0 and sp.simplify(r[h][1]/r[h][0] - c[h][1]) == 0 for h in ('opp', 'same'))
    print(f"(Delta,J)=({D},{J}): closed forms exact: {ok}", flush=True)

"""Build (and cache) the 14 differential-basis words for <T T O> (l = 2): Qrel_w relative to M_ref (no simplification)."""
import time, pickle, os
import sympy as sp
import embed2
from embed2 import word3, shifts, D1, D2, CACHE
from structures2 import words
t = time.time()
out = {}
for w in words(2):
    e, S = embed2._mk_engine(embed2.PAIRS3)
    v5 = (S['z51']*S['s25'] - S['z52']*S['s15'])/2
    s1, s2 = shifts(w)
    d1, d2 = 3 + s1, 3 + s2
    al = (d1 + d2 - embed2.Delta + embed2.J)/2; be = (embed2.Delta + embed2.J + d1 - d2)/2; ga = (embed2.Delta + embed2.J - d1 + d2)/2
    e.logK = embed2.J*sp.log(v5) - al*sp.log(S['s12']) - be*sp.log(S['s15']) - ga*sp.log(S['s25'])
    Phi = sp.Integer(1)
    for kind in reversed(w):
        Phi = embed2.op(e, Phi, kind, 1, 2)
    h = sp.Rational(s1 - s2, 2)
    out[w] = sp.expand(Phi*S['s15']**(-h)*S['s25']**h)
    print(w, len(sp.Add.make_args(out[w])), "%.0f s" % (time.time() - t), flush=True)
pickle.dump((out, S, v5), open(os.path.join(CACHE, 'rel_l2.pkl'), 'wb'))
print("done", flush=True)

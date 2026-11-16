"""
<T T O> exact helicity map at a rational point and the Weyl parity of its MFT normalisation (coordinator goal 2).
 1. conserved space at (Delta, J) in the frame at infinity (infcons2, exact), d = 3 relation among the words;
 2. exact helicity lines by the momentum-space projection calibrated on l = 1 (ft_helicity.helicity_components_exact):
    opposite <=> T(+,-) = 0, same <=> T(+,+) = 0  (modulo the relation);
 3. Weyl-covariant family: e_h(Delta~) := T e_h(Delta)  (word-wise b^(0)/b_w); check that it lies on the helicity line at (3-Delta, -1-J);
 4. CHL's T^{E,h} = N_h e_h (chl_momentum; exact up to 30 digits) at both points;  kappa_h := N_h(Delta)/N_h(Delta~).
    MFT normalisation of e_h is Weyl-even  <=>  kappa_h A_h is helicity-independent and equal to its l = 1 value (same (Delta,J)), with
    A_opp = (-J)_{2l}/(J+1)_{2l}, A_same = (Delta-1)_{2l}/(2-Delta)_{2l} (CHL's MFT helicity factors; for l = 1 kappa_h A_h = 23170.475..., t5 l=1 run).
Output: out_t5_helicity_l2.txt
"""
import os, sys, time, pickle
import sympy as sp
from infcons2 import conserved_space
from ft_helicity import helicity_components_exact, R_inv
from chl_momentum import chl_components, our_components
from infframe import words, CACHE
from embed2 import shifts
from t1_structures_l1 import b_ratio
from embed2 import Delta as sD, J as sJ
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
R = sp.Rational
import sys as _s
Dv, Jv = (R(97, 10), R(13, 3)) if len(_s.argv) < 3 else (R(_s.argv[1]), R(_s.argv[2]))
Dt, Jt = 3 - Dv, -1 - Jv
l = 2
t = time.time()
fn = os.path.join(CACHE, 'cons_l2_%s_%s.pkl' % (str(Dv).replace('/', '_'), str(Jv).replace('/', '_')))
if os.path.exists(fn):
    ws, ns, rel = pickle.load(open(fn, 'rb'))
else:
    ws, ns, rel = conserved_space(l, Dv, Jv)
    pickle.dump((ws, ns, rel), open(fn, 'wb'))
pr(f"# <TTO> at (Delta,J) = ({Dv},{Jv}): conservation null space {len(ns)}, word relations {len(rel)}  ({time.time()-t:.0f} s)")
vecs = [{w: v[i] for i, w in enumerate(ws)} for v in ns]
T = [helicity_components_exact(l, c, Dv, Jv) for c in vecs]
pr("  helicity components of the null vectors (exact):", [{k: complex(sp.N(x)) for k, x in Ti.items()} for Ti in T])
lines = {}
for name, cond in (('opp', (1, -1)), ('same', (1, 1))):
    cs = sp.symbols('c0:%d' % len(vecs))
    eqs = [sum(cs[i]*T[i][cond] for i in range(len(vecs))), sum(cs[i]*T[i][(-cond[0], -cond[1])] for i in range(len(vecs)))]
    sol = sp.linsolve(eqs, cs)
    pr(f"  {name}: solution set {sol}")
    sol = list(sol)[0]
    free = [s for s in cs if any(s in e.free_symbols for e in sol)]
    # pick a representative with nonzero helicity content
    for choice in ([1] + [0]*(len(free) - 1), [0, 1] + [0]*(len(free) - 2) if len(free) > 1 else [1]):
        vals = [e.subs(dict(zip(free, choice))) for e in sol]
        e_h = {w: sp.nsimplify(sum(vals[i]*vecs[i][w] for i in range(len(vecs)))) for w in ws}
        Th = helicity_components_exact(l, e_h, Dv, Jv)
        if any(x != 0 for x in Th.values()): break
    lines[name] = e_h
    pr(f"  {name} line representative: helicity components {[complex(sp.N(x)) for x in Th.values()]}")
# Weyl-covariant images and checks at the reflected labels
res = {}
for name, idx in (('same', 0), ('opp', 1)):
    e = lines[name]
    f = {w: sp.nsimplify(e[w]/b_ratio(shifts(w)[0] - shifts(w)[1]).subs({sD: Dv, sJ: Jv})) for w in ws}
    Tf = helicity_components_exact(l, f, Dt, Jt)
    pr(f"  {name}: T e_h at reflected labels, helicity components {[complex(sp.N(x)) for x in Tf.values()]}  (the 'wrong' pair must vanish)")
    C1 = chl_components(l, Dv, Jv); C2 = chl_components(l, Dt, Jt)
    O1 = our_components(l, e, Dv, Jv); O2 = our_components(l, f, Dt, Jt)
    r1 = [C1[h][idx]/O1[h] for h in O1 if abs(O1[h]) > 1e-25*max(abs(x) for x in O1.values())]
    r2 = [C2[h][idx]/O2[h] for h in O2 if abs(O2[h]) > 1e-25*max(abs(x) for x in O2.values())]
    kappa = r1[0]/r2[0]
    A = (sp.rf(-Jv, 2*l)/sp.rf(Jv + 1, 2*l)) if name == 'opp' else (sp.rf(Dv - 1, 2*l)/sp.rf(2 - Dv, 2*l))
    res[name] = kappa*float(A)
    pr(f"  {name}: N_h(Delta) on the two components {['%.12g%+.3gj' % (z.real, z.imag) for z in r1]};  kappa_h = {kappa:.12g};  kappa_h A_h = {res[name]:.12g}")
pr(f"  RESULT: kappa_opp A_opp / kappa_same A_same = {res['opp']/res['same']}   (Weyl-even MFT normalisation requires 1; l = 1 value of kappa A: 23170.475005920...)")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t5_helicity_l2_%s_%s.txt' % (str(Dv).replace('/', '_'), str(Jv).replace('/', '_'))), 'w').write("\n".join(out) + "\n")

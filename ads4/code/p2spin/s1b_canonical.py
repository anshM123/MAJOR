"""
p2spin/s1b_canonical.py (P2_SPINNING.md §1.2): canonical exact bases of the conserved structures (relations quotiented), printed with
factored coefficients.  Reads cache/s1_structures.pkl.  -> cache/s1b_canonical.txt, cache/s1b_canonical.pkl
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import p2_core as C
import s1_structures as S

OUT = []
pr = lambda *a: (print(' '.join(str(x) for x in a), flush=True), OUT.append(' '.join(str(x) for x in a)))
name_e = lambda t: f'V1^{2-t[0]-t[1]} V2^{2-t[0]-t[2]} V5^(J-{t[1]+t[2]}) H12^{t[0]} H15^{t[1]} H25^{t[2]}'
name_o = lambda t: f'e{t[0]} * (a,b,c)=({t[1]},{t[2]},{t[3]})'
res = {}
EV_ORDER = [C.EVEN.index(t) for t in [(2, 0, 0), (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1),
                                       (0, 2, 0), (0, 0, 2), (1, 1, 1), (0, 2, 1), (0, 1, 2), (0, 2, 2)]]
OD_ORDER = list(range(len(C.ODD)))
for sector, order, nm in (('even', EV_ORDER, name_e), ('odd', OD_ORDER, name_o)):
    for sigma in (1, -1):
        ns, rel = S.load_sol(sector, sigma)
        if not ns or len(ns) == len(rel):
            pr(f'{sector} sigma={sigma}: no physical structure'); continue
        if rel:
            pr(f'{sector} sigma={sigma}: relations:')
            for r in rel:
                pr('   ', ' + '.join(f'({x})*[{nm((C.EVEN if sector=="even" else C.ODD)[i])}]' for i, x in enumerate(r) if x != 0))
        B = S.canonical(ns, rel, order)
        res[(sector, sigma)] = B
        pr(f'{sector} sigma={sigma}: {len(B)} canonical basis vector(s):')
        for k, v in enumerate(B):
            pr(f'  basis {k}:')
            for i, x in enumerate(v):
                if x != 0:
                    pr(f'     {nm((C.EVEN if sector=="even" else C.ODD)[i]):48s}  {x}')
pickle.dump({k: [[sp.srepr(x) for x in v] for v in B] for k, B in res.items()}, open(os.path.join(HERE, 'cache', 's1b_canonical.pkl'), 'wb'))
open(os.path.join(HERE, 'cache', 's1b_canonical.txt'), 'w').write('\n'.join(OUT) + '\n')

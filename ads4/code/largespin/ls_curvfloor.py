"""
largespin/ls_curvfloor.py (LARGESPIN_ANALYTIC.md §3.3): at large b-hat the hybrid is curvature-dominated, so the flat-only floor N_flat/5
under-states the margin.  For the states of cache/struct_<Dg>.json with b-hat >= 300, recompute the floor as N_flat/5 + curv (curv = exact
leading-Regge curvature, Mehler-Dirichlet, [NE]; certified positive on b-hat in [91.2, 2 Lam] for m >= 1.5 Lam by ROUTEA_LARGE_B) and
print the ratio (new floor)/(old floor), by which every tolerance of ls_struct scales.  Output: cache/curvfloor_<Dg>.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
import ls_core as C
Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
rows = json.load(open(os.path.join(C.CACHE, 'struct_%d.json' % Dg)))
out = ['%5s %7s %8s %9s %9s %9s %s' % ('J', 'm/Lam', 'bhat', 'old/Ns', 'new/Ns', 'ratio', 'a_tol[F], a_tol[D], a_tol[A] rescaled')]
for r in rows:
    if r['bhat'] < 300:
        continue
    st = C.hyb_fast(Dg, r['J'], r['mfl'], r['h'])
    Nf = abs(st['B2']) + abs(st['B3']); cv = st['cv2'] + st['cv3']
    new = (Nf/5 + cv)/st['Ns']; old = r['floor_rel']; q = new/old
    out.append('%5d %7.3f %8.1f %9.4f %9.4f %9.2f   %.3g  %.3g  %.4g' % (r['J'], r['mfl'], r['bhat'], old, new, q, q*r['atol']['F'], q*r['atol']['D'], q*r['atol']['A']))
txt = '\n'.join(out); print(txt)
open(os.path.join(C.CACHE, 'curvfloor_%d.txt' % Dg), 'w').write(txt + '\n')

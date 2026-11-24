"""
gc_table (GAP_C3.md §4): Delta_1(J, xi) from the certified node data (cache/gc_cert_J*.json, cache/gc_int_J*.json) and Lemma L3.

Target (Lemma T form, GAP_UNIFORM §2-3): positivity at gap Lambda needs  m^2 Pi rel_env + sup|Gamma_K2| <= G_cert(Lambda0) (Lambda/Lambda0)^2;
for Lambda >> Lambda0 this is  rel_env <= rel* = 0.9 G_cert / (Pi mu^2 Lambda0^2),  mu = 1.3 (largest tabulated m/Lambda), Pi = 3
(state + reference + a_h fit, GAP_C3 §4.2 [M]).  Delta_1 = smallest Delta (log grid) at which the certified bound is <= rel*; the bound is
decreasing in Delta, so it then holds for every Delta >= Delta_1.  Also reported: C3 := sup_{Delta >= Delta_1} rel Delta^3 (the explicit
constant), and the order-2 (remainder after the leading term) Delta_1 for comparison.
Usage: python3 gc_table.py  -> cache/gc_table.txt, cache/gc_table.json
"""
import os, sys, json, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_guard  # noqa: F401
import numpy as np
import gc_lemma as LM

GCERT = {0: 1.99e5, 2: 6.34e4, 4: 1.10e4, 5: 3.01e4, 6: 3.60e3, 7: 2.34e3, 8: 316.0, 10: 153.0, 12: 30.5}   # GAP_UNIFORM §3.2
LAM0 = 998.50
MU = 1.3
PI = 3.0


def target(J):
    return 0.9*GCERT[J]/(PI*MU*MU*LAM0*LAM0)


def main():
    out = {}
    lines = ['GAP_C3 Delta_1 table (certified node data + Lemma L3; rel* = 0.9 G_cert/(Pi mu^2 Lambda0^2), Pi = %g, mu = %g)' % (PI, MU),
             '%-4s %-16s %-6s %-4s %-9s %-10s %-10s %-10s %-10s' % ('J', 'xi', 'rho', 'deg', 'rel*', 'Delta_1', 'C3', 'tau', 'Delta_1(o2)')]
    # cache/prefix/: records made before the exact-xi phase fix (|grad Phi(z*)| <= 5e-17, GAP_C3 §3); records in cache/ override them
    merged = {}
    for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'prefix', 'gc_cert_J*.json'))) + sorted(glob.glob(os.path.join(HERE, 'cache', 'gc_cert_J*.json'))):
        merged.setdefault(os.path.basename(fn), {}).update(json.load(open(fn)))
    for name, db in sorted(merged.items()):
        best = {}
        for k, d in db.items():
            if d['bad'] or not np.isfinite(d['Mg']):
                continue
            J = d['J']
            if J not in GCERT:
                continue
            xk = tuple(d['xint']) if d.get('xint') else (round(d['xi'], 4),)
            D1, r = LM.delta1(d, target(J))
            if D1 is None:
                continue
            if (J, xk) not in best or D1 < best[(J, xk)][0]:
                c3 = r[0]*D1**3
                D2, _ = LM.delta1(d, target(J), order=2)
                best[(J, xk)] = (D1, d['rho'], d['deg'], c3, r[1], D2)
        for (J, xk), (D1, rho, deg, c3, tau, D2) in sorted(best.items()):
            xs = ('%.4f' % xk[0]) if len(xk) == 1 else ('[%.4f,%.4f]' % xk)
            lines.append('%-4d %-16s %-6g %-4d %-9.3g %-10.3g %-10.3g %-10.2f %-10.3g' % (J, xs, rho, deg, target(J), D1, c3, tau,
                                                                                    D2 if D2 else float('nan')))
            out['%d_%s' % (J, xs)] = dict(J=J, xi=xk, Delta1=D1, rho=rho, deg=deg, C3=c3, tau=tau, Delta1_order2=D2, target=target(J))
    open(os.path.join(HERE, 'cache', 'gc_table.txt'), 'w').write('\n'.join(lines) + '\n')
    json.dump(out, open(os.path.join(HERE, 'cache', 'gc_table.json'), 'w'), indent=1)
    print('\n'.join(lines))


if __name__ == '__main__':
    main()

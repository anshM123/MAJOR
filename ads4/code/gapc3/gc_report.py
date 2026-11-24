"""
gc_report (GAP_C3.md §4-5): Delta_1 for the spinning-amplitude certificates (cache/gc_spin.json) and the [NE] ideal sup-norm probes
(cache/gc_ideal.json), with the same Lemma L3 and targets as gc_table.  Prints a text report.
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_guard  # noqa: F401
import numpy as np
from scipy.special import eval_legendre
import gc_lemma as LM
import gc_table as TB


def spin_report():
    fn = os.path.join(HERE, 'cache', 'gc_spin.json')
    if not os.path.exists(fn):
        return ['(no gc_spin.json)']
    db = json.load(open(fn))
    best = {}
    for k, d in db.items():
        if d['bad'] or not np.isfinite(d['Mg']):
            continue
        D1, r = LM.delta1(d, TB.target(d['J']))
        if D1 is None:
            continue
        key = (d['kind'], d['J'], round(d['xi'], 4))
        if key not in best or D1 < best[key][0]:
            best[key] = (D1, d['rho'], d['deg'], r[0]*D1**3, d['Mg'], d['h0'], d['S3'])
    lines = ['spinning amplitudes (frozen <TTO> vectors at Delta = 1002), certified; phase from the log-factored integrand',
             '%-4s %-3s %-7s %-6s %-4s %-10s %-10s %-9s %-8s %-6s' % ('kind', 'J', 'xi', 'rho', 'deg', 'Delta_1', 'C3', 'Mh_env', 'D_h', 'S3')]
    for (kd, J, xi), (D1, rho, deg, c3, Mh, h0, S3) in sorted(best.items()):
        lines.append('%-4s %-3d %-7.4f %-6g %-4d %-10.3g %-10.3g %-9.3g %-8.3g %-6.3g' % (kd, J, xi, rho, deg, D1, c3, Mh, h0, S3))
    return lines


def ideal_report():
    fn = os.path.join(HERE, 'cache', 'gc_ideal.json')
    if not os.path.exists(fn):
        return ['(no gc_ideal.json)']
    db = json.load(open(fn))
    lines = ['[NE] ideal (sampled) sup norms -> Lemma L3 (Schwarz form, D = 3): Delta_1^ideal; rel bound at Delta = 1000 / 1e4',
             '%-4s %-3s %-6s %-6s %-9s %-9s %-10s %-10s %-10s' % ('kind', 'J', 'xi', 'rho', 'M3/rho^3', 'Mh_env', 'Delta_1', 'rel(1e3)', 'rel(1e4)')]
    for k, r in sorted(db.items()):
        J, xi = r['J'], r['xi']
        Dh = abs(eval_legendre(J, 1 - 2*xi*xi)) if J > 0 else 1.0
        Mh = r['Mh_rel']*Dh
        d = dict(rho=r['rho'], deg=3, Bphi=[0, 0, 0, r['M3']], remPhi=0.0, Bh=[Mh, 0, 0, 0], remh=0.0, E=0.0)
        D1, _ = LM.delta1(d, TB.target(J))
        b3 = LM.best(d, 1e3); b4 = LM.best(d, 1e4)
        lines.append('%-4s %-3d %-6.3f %-6.3f %-9.3g %-9.3g %-10.3g %-10s %-10s' % (
            r['kind'], J, xi, r['rho'], r['M3']/r['rho']**3, Mh, D1 if D1 else float('nan'),
            ('%.3g' % b3[0]) if b3 else 'n/a', ('%.3g' % b4[0]) if b4 else 'n/a'))
    return lines


if __name__ == '__main__':
    print('\n'.join(spin_report() + [''] + ideal_report()))

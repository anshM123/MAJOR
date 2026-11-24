"""
gm_report (GAP_MORSE.md §4-5): Delta_1(J, xi) tables (scalar, from cache/gm_table.json; spinning graviton amplitudes from
cache/gm_spin.json with the J-independent phase records of cache/gm_cert_J0.json), side by side with GAP_C3's certified Delta_1
(gapc3/cache/gc_table.json and the gc_spin records, read only).  Usage: python3 gm_report.py -> cache/gm_report.txt, gm_spin_table.json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gm_guard  # noqa: F401
import gm_table as TB

C3 = os.path.join(HERE, '..', 'gapc3', 'cache')
XIS = (0.1, 0.3, 0.5, 0.7, 0.85)
C3SPIN = {('C3H', 0, 0.5): 1.41e6, ('C3H', 2, 0.5): 1.5e6, ('C3H', 2, 0.7): 1.78e6, ('EE', 0, 0.5): 1.0e6, ('EE', 2, 0.3): 1.68e6,
          ('EE', 2, 0.5): 1.0e6, ('EE', 2, 0.7): 1.19e6}          # GAP_C3 §4.3 / gc_report.txt


def scalar_lines():
    T = json.load(open(os.path.join(HERE, 'cache', 'gm_table.json')))
    old = json.load(open(os.path.join(C3, 'gc_table.json')))
    lines = ['scalar k = 2 with Q_J: certified Delta_1 (GAP_MORSE, Lemma M)  [GAP_C3 Lemma L3 value]',
             '%-3s %-9s ' % ('J', 'rel*') + ' '.join('%-22s' % ('xi = %g' % x) for x in XIS) + ' max over xi']
    for J in sorted(TB.GCERT):
        row, mx = [], 0.0
        for x in XIS:
            o = T.get('%d_%.4f' % (J, x))
            c = old.get('%d_%.4f' % (J, x), {}).get('Delta1')
            v = o.get('Delta1') if o else None
            if v:
                mx = max(mx, v)
            else:
                mx = float('inf')
            row.append('%-22s' % (('%.2e' % v if v else ('none' if o else '-')) + ' [%s]' % ('%.1e' % c if c else '-')))
        lines.append('%-3d %-9.2g ' % (J, TB.target(J)) + ' '.join(row) + ' %.2e' % mx)
    return lines


def spin_table():
    fn = os.path.join(HERE, 'cache', 'gm_spin.json')
    if not os.path.exists(fn):
        return ['(no gm_spin.json)'], {}
    sp = json.load(open(fn))
    ph = list(json.load(open(os.path.join(HERE, 'cache', 'gm_cert_J0.json'))).values())
    out = {}
    keys = sorted(set((o['kind'], o['J'], round(o['xi'], 4)) for o in sp.values()))
    for (kind, J, xi) in keys:
        amps = [o for o in sp.values() if (o['kind'], o['J'], round(o['xi'], 4)) == (kind, J, xi)]
        recs = [o for o in ph if round(o['xi'], 4) == xi]
        b = TB.node_best(recs, J, amps=amps, mixed=True)
        if b is None:
            out['%s_%d_%.4f' % (kind, J, xi)] = dict(kind=kind, J=J, xi=xi, Delta1=None)
            continue
        D1, info, rp, rh, deg, degh = b
        out['%s_%d_%.4f' % (kind, J, xi)] = dict(kind=kind, J=J, xi=xi, Delta1=D1, rho_Phi=rp, rho_h=rh, deg=deg, deg_h=degh, **info)
    lines = ['spinning graviton amplitudes (frozen <TTO> at Delta = 1002), certified Delta_1 (Lemma M)  [GAP_C3]',
             '%-4s %-3s %-6s %-9s %-7s %-7s %-6s %-7s %-8s %-9s' % ('kind', 'J', 'xi', 'Delta_1', 'rhoPhi', 'rho_h', 'deg_h', 'r', 'M_H',
                                                                  'GAP_C3')]
    for k, o in sorted(out.items()):
        c = C3SPIN.get((o['kind'], o['J'], o['xi']))
        if o['Delta1'] is None:
            lines.append('%-4s %-3d %-6.3g none' % (o['kind'], o['J'], o['xi']))
            continue
        lines.append('%-4s %-3d %-6.3g %-9.3g %-7.3g %-7.3g %-6d %-7.3g %-8.3g %-9s' % (
            o['kind'], o['J'], o['xi'], o['Delta1'], o['rho_Phi'], o['rho_h'], o['deg_h'], o['st']['r'], o['st']['MH'],
            '%.2e' % c if c else '-'))
    json.dump(out, open(os.path.join(HERE, 'cache', 'gm_spin_table.json'), 'w'), indent=1)
    return lines, out


if __name__ == '__main__':
    L = scalar_lines() + [''] + spin_table()[0]
    open(os.path.join(HERE, 'cache', 'gm_report.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))

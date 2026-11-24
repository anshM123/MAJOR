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


def sci(v):
    if v is None:
        return '—'
    e = int(f'{v:.1e}'.split('e')[1]); m = v/10**e
    sup = str(e).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
    return f'{m:.1f}·10{sup}'


def markdown_scalar():
    T = json.load(open(os.path.join(HERE, 'cache', 'gm_table.json')))
    old = json.load(open(os.path.join(C3, 'gc_table.json')))
    L = [r'| J | rel*_J | ' + ' | '.join('ξ = %g' % x for x in XIS) + r' | **Δ₁\*(J) = max** | GAP_C3 max |',
         '|---|---|' + '---|'*len(XIS) + '---|---|']
    for J in sorted(TB.GCERT):
        vals = [T.get('%d_%.4f' % (J, x), {}).get('Delta1') for x in XIS]
        olds = [old.get('%d_%.4f' % (J, x), {}).get('Delta1') for x in XIS]
        mx = None if any(v is None for v in vals) else max(vals)
        omx = None if any(v is None for v in olds) else max(olds)
        L.append('| %d | %s | ' % (J, sci(TB.target(J))) + ' | '.join(sci(v) for v in vals) + ' | **%s** | %s |' % (sci(mx), sci(omx)))
    return L


def breakdown_lines():
    T = json.load(open(os.path.join(HERE, 'cache', 'gm_table.json')))
    lines = ['binding set-ups (scalar): radii, chart quantities and the four remainder terms at Delta_1',
             '%-3s %-5s %-8s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-8s %-8s %-8s %-8s' % (
                 'J', 'xi', 'Delta1', 'rhoPhi', 'rho_h', 'rho2', 'r', 'rc', 'eta2', 'lamJ2', 'M_H', 'kappa', 't1', 't3', 't4')]
    for k in sorted(T, key=lambda k: (T[k]['J'], T[k]['xi'])):
        o = T[k]
        if not o.get('Delta1'):
            continue
        st, p = o['st'], o['parts']
        lines.append('%-3d %-5.2f %-8.3g %-6.3g %-6.3g %-6.3g %-6.3g %-6.3g %-6.3g %-6.3g %-6.3g %-8.2g %-8.2g %-8.2g %-8.2g' % (
            o['J'], o['xi'], o['Delta1'], o['rho_Phi'], o['rho_h'], st['rho2'], st['r'], st['rc'], st['eta2'], st['lamJ2'], st['MH'],
            st['kap'], p['t1'], p['t3'], p['t4']))
    return lines


def markdown_breakdown(xis=(0.1, 0.5)):
    T = json.load(open(os.path.join(HERE, 'cache', 'gm_table.json')))
    L = ['| J | ξ | Δ₁ | ρ_Φ (deg) | ρ_h | ρ₂ | r_c | η(ρ₂) | λ_J(ρ₂) | M_H | κ | t₁ | t₃ | t₄ |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for J in sorted(TB.GCERT):
        for x in xis:
            o = T.get('%d_%.4f' % (J, x))
            if not o or not o.get('Delta1'):
                continue
            st, p = o['st'], o['parts']
            L.append('| %d | %g | %s | %.3g (%d) | %.3g | %.3g | %.3g | %.2f | %.2f | %.2f | %.3f | %.1e | %.1e | %.1e |' % (
                J, x, sci(o['Delta1']), o['rho_Phi'], o['deg'], o['rho_h'], st['rho2'], st['rc'], st['eta2'], st['lamJ2'], st['MH'],
                st['kap'], p['t1'], p['t3'], p['t4']))
    return L


def markdown_spin():
    fn = os.path.join(HERE, 'cache', 'gm_spin_table.json')
    if not os.path.exists(fn):
        return ['(none)']
    T = json.load(open(fn))
    L = ['| kind | J | ξ | rel*_J | Δ₁ (certified) | ρ_h (deg 6) | r_c | M_H | GAP_C3 Δ₁ |', '|---|---|---|---|---|---|---|---|---|']
    for k in sorted(T, key=lambda k: (T[k]['kind'], T[k]['J'], T[k]['xi'])):
        o = T[k]
        c = C3SPIN.get((o['kind'], o['J'], o['xi']))
        if not o.get('Delta1'):
            L.append('| %s | %d | %g | %s | none | — | — | — | %s |' % (o['kind'], o['J'], o['xi'], sci(TB.target(o['J'])), sci(c)))
            continue
        L.append('| %s | %d | %g | %s | %s | %.3g | %.3g | %.2f | %s |' % (o['kind'], o['J'], o['xi'], sci(TB.target(o['J'])),
                                                                     sci(o['Delta1']), o['rho_h'], o['st']['rc'], o['st']['MH'], sci(c)))
    return L


def markdown_intervals():
    fn = os.path.join(HERE, 'cache', 'gm_int_table.json')
    T = json.load(open(fn)) if os.path.exists(fn) else {}
    N = json.load(open(os.path.join(HERE, 'cache', 'gm_table.json')))
    L = ['| J | ξ-interval | ρ | ‖E′‖_F (uniform) | λ_min | Δ₁ (whole interval) | node Δ₁ at ξ₀ |', '|---|---|---|---|---|---|---|']
    best = {}
    for k, o in T.items():
        if not o.get('Delta1'):
            continue
        key = (o['J'], o['xi0'], o['rx'])
        if key not in best or o['Delta1'] < best[key]['Delta1']:
            best[key] = o
    for (J, x0, rx), o in sorted(best.items()):
        nd = N.get('%d_%.4f' % (J, x0), {}).get('Delta1')
        L.append('| %d | [%.4f, %.4f] | %.3g | %.2g | %.3f | %s | %s |' % (J, x0 - rx, x0 + rx, o['rho'], o['EF'], o['lam_min'],
                                                                         sci(o['Delta1']), sci(nd)))
    return L


if __name__ == '__main__':
    L = scalar_lines() + [''] + breakdown_lines() + [''] + spin_table()[0]
    L += [''] + markdown_scalar() + [''] + markdown_breakdown() + [''] + markdown_spin() + [''] + markdown_intervals()
    open(os.path.join(HERE, 'cache', 'gm_report.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))

"""
gm_table (GAP_MORSE.md §4): certified Delta_1(J, xi) from the Taylor-model records (cache/gm_cert_J*.json) and Lemma M (gm_chart).
For each node all pairs (phase record at rho_Phi, amplitude record at rho_h <= rho_Phi, same degree) are tried; the polynomial parts are
radius-independent (exact Taylor coefficients), the radii enter only through the remainders.  Targets as GAP_C3 §4.1:
rel* = 0.9 G_cert(J)/(Pi mu^2 Lambda0^2), Pi = 3, mu = 1.3 (GAP_UNIFORM §3.2 tolerances).
Usage: python3 gm_table.py [J ...]  -> cache/gm_table.json, cache/gm_table.txt
"""
import os, sys, json, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gm_guard  # noqa: F401
import numpy as np
import gm_chart as G

GCERT = {0: 1.99e5, 2: 6.34e4, 4: 1.10e4, 5: 3.01e4, 6: 3.60e3, 7: 2.34e3, 8: 316.0}
LAM0, MU, PI = 998.50, 1.3, 3.0


def target(J):
    return 0.9*GCERT[J]/(PI*MU*MU*LAM0*LAM0)


def node_best(recs, J, tgt=None, amps=None, mixed=False):
    """recs: phase records (gm_cert); amps: amplitude records (default: the same records)."""
    tgt = tgt or target(J)
    best = None
    ph = [r for r in recs if r.get('badPhi', r.get('bad')) == 0 and np.isfinite(r.get('remPhi', np.nan))]
    am = [r for r in (amps if amps is not None else recs) if r.get('bad') == 0 and np.isfinite(r.get('remh', np.nan))]
    for a in ph:
        for b in am:
            if b['rho'] > a['rho'] or (b['deg'] != a['deg'] and not mixed):
                continue
            d = dict(deg=a['deg'], EF=a['EF'], trE=a['trE'], bPhi=a['bPhi'], bTr=a['bTr'], rho=a['rho'], remPhi=a['remPhi'],
                     rho_h=b['rho'], remh=b['remh'], bh=b['bh'], deg_h=b['deg'])
            D1, info = G.delta1(d, tgt)
            if D1 is not None and (best is None or D1 < best[0]):
                best = (D1, info, a['rho'], b['rho'], a['deg'], b['deg'])
    return best


def main(Js=None):
    out = {}
    fn = os.path.join(HERE, 'cache', 'gm_table.json')
    if os.path.exists(fn):
        out = json.load(open(fn))
    fph = os.path.join(HERE, 'cache', 'gm_phase.json')
    extra = list(json.load(open(fph)).values()) if os.path.exists(fph) else []      # phase-only records (J-independent)
    for f in sorted(glob.glob(os.path.join(HERE, 'cache', 'gm_cert_J*.json'))):
        db = json.load(open(f))
        recs = list(db.values())
        J = recs[0]['J']
        if Js and J not in Js:
            continue
        for xi in sorted(set(round(r['xi'], 4) for r in recs)):
            rr = [r for r in recs if round(r['xi'], 4) == xi]
            ex = [r for r in extra if round(r['xi'], 4) == xi]
            b = node_best(rr + ex, J, amps=rr, mixed=True)
            key = '%d_%.4f' % (J, xi)
            if b is None:
                out[key] = dict(J=J, xi=xi, Delta1=None, target=target(J))
            else:
                D1, info, rp, rh, deg, degh = b
                out[key] = dict(J=J, xi=xi, Delta1=D1, target=target(J), rho_Phi=rp, rho_h=rh, deg=deg, deg_h=degh, **info)
            print(key, out[key].get('Delta1'), flush=True)
    json.dump(out, open(fn + '.tmp', 'w'), indent=1); os.replace(fn + '.tmp', fn)
    lines = ['GAP_MORSE certified Delta_1 (Lemma M; rel* = 0.9 G_cert/(Pi mu^2 Lambda0^2), Pi = 3, mu = 1.3)',
             '%-3s %-6s %-9s %-9s %-7s %-7s %-4s %-7s %-7s %-6s %-7s %-7s %-7s' % ('J', 'xi', 'rel*', 'Delta1', 'rhoPhi', 'rho_h', 'deg',
                                                                             'rho2', 'r', 'tau', 'M_H', 'kappa', 's_out')]
    for k in sorted(out, key=lambda k: (out[k]['J'], out[k]['xi'])):
        o = out[k]
        if o.get('Delta1') is None:
            lines.append('%-3d %-6.3f %-9.3g %s' % (o['J'], o['xi'], o['target'], 'none'))
            continue
        st = o['st']
        lines.append('%-3d %-6.3f %-9.3g %-9.3g %-7.3g %-7.3g %-4d %-7.3g %-7.3g %-6.2g %-7.3g %-7.3g %-7.3g' % (
            o['J'], o['xi'], o['target'], o['Delta1'], o['rho_Phi'], o['rho_h'], o['deg'], st['rho2'], st['r'], st['tau'], st['MH'],
            st['kap'], st['sout']))
    open(os.path.join(HERE, 'cache', 'gm_table.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))




def intervals():
    """xi-interval certificates (cache/gm_int_J*.json, gm_cert8): Delta_1 valid for every xi in [xi0 - rx, xi0 + rx]."""
    out = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'cache', 'gm_int_J*.json'))):
        for k, r in json.load(open(f)).items():
            if r.get('bad') != 0 or not r.get('grad_ok') or not (r.get('lam_min', -1) > 0):
                out[k] = dict(J=r['J'], xi0=r['xi0'], rx=r['rx'], rho=r['rho'], Delta1=None)
                continue
            d = dict(r)
            D1, info = G.delta1(d, target(r['J']))
            out['%d_%s' % (r['J'], k)] = dict(J=r['J'], xi0=r['xi0'], rx=r['rx'], rho=r['rho'], deg=r['deg'], Delta1=D1,
                                             EF=r['EF'], lam_min=r['lam_min'], **(info or {}))
            print(r['J'], r['xi0'], r['rx'], r['rho'], r['deg'], 'EF %.3g' % r['EF'], 'Delta1', D1, flush=True)
    json.dump(out, open(os.path.join(HERE, 'cache', 'gm_int_table.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    if sys.argv[1:2] == ['int']:
        intervals()
    else:
        main([int(a) for a in sys.argv[1:]] or None)

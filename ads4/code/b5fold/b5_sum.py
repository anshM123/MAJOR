"""
b5_sum (BATCH5_FOLD.md §0): the per-state table.  Inputs (this note's cache, all from the scripts above):
  even J (own-row route): state/b5_state_J<J>.json (ng_bind budget with the direct row; E_interp also on lam_c) + state/sn_even_J<J>.json
     (measured window).  Sign-free L = lam_c / (E_trunc + E_interp[lam_c] + E_edge + E_mc + window x3 envelope); lambda-level: E_trunc
     -> 3|lam_K2 - lam_K1| (NEARGAP2 §1.2); also with the 'measured' window (sn_even) and with the protocol E_interp (on lam).
  odd J, hs_odd route (J = 23, 25, 31; 17, 19 diagnostics): odd_<tag>/b5_odd_J<J>.json (L x3 / signed with E_repr on lam_c).
  odd J, own-row route (J = 17, 19; 25 cross-check): state/b5_state_J<J>.json + state/b5_window_J<J>.json.
A state is 'held' when its mass exceeds every row of its J (the model holds the own row in 1/m): as the hs_odd protocol does for every
odd-J mass and NEARGAP_FILL §2.3 for J = 11, 13, 15.
Usage: python3 b5_sum.py   -> cache/b5_sum.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b5_guard  # noqa: F401
import numpy as np
C = os.path.join(HERE, 'cache')
LAM = 998.5   # Lambda(1000) (only for the held flag; exact value from the state files' rows)


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def ld(p):
    p = os.path.join(C, p)
    return json.load(open(p)) if os.path.exists(p) else None


def even(J):
    st = ld('state/b5_state_J%d.json' % J); sn = ld('state/sn_even_J%d.json' % J)
    out = []
    if not st or not sn:
        return out
    for s in st['states']:
        d = s['direct']; mfl = s['mfl']
        w = sn['states'][[k for k in sn['states'] if abs(float(k) - mfl) < 1e-6][0]]
        lamc = d['lam_central']
        ec = s['interp_c']['E_interp_c']
        rest = s['E_tot'] - d['E_low'] - d['E_interp'] + ec
        rest_l = s['E_tot_lam'] - d['E_low'] - d['E_interp'] + ec
        mrow = max(mphys(D, J) for D in d['rows'])
        held = mfl*LAM > mrow*(1 + 1e-4)
        out.append(dict(J=J, mfl=mfl, route='own row (even)', lam=lamc, jk=s['jk_sigma'], pub=(s['k2'] or {}).get('lam_model'),
                        E_trunc=d['E_trunc'], E_trunc_lam=d['E_trunc_lam'], E_interp=d['E_interp'], E_interp_c=ec, E_edge=d['E_edge'],
                        E_mc=s['E_mc'], E_low_conv=d['E_low'], E_env3=w['E_env3'], E_env10=w['E_env10'], E_meas=w['E_meas'],
                        budget=rest + w['E_env3'], budget_lam=rest_l + w['E_env3'], L=lamc/(rest + w['E_env3']),
                        L_lam=lamc/(rest_l + w['E_env3']), L_meas=lamc/(rest + w['E_meas']), L_10=lamc/(rest + w['E_env10']),
                        L_prot=lamc/(s['E_tot'] - d['E_low'] + w['E_env3']), L_conv=s['L'], held=held, rows=d['rows'],
                        xmax=d['xmax'], G_edge=d['G_edge'], failed=s['failed'], rz=s['rz']))
    return out


def odd_hs(J, tag='main'):
    r = ld('odd_%s/b5_odd_J%d.json' % (tag, J))
    out = []
    for k, o in sorted((r or {}).items(), key=lambda t: float(t[0])):
        out.append(dict(J=J, mfl=o['mfl'], route='hs_odd ' + tag, lam=o['lamc'], lam_main=o['main'], jk=o['jk'], pub=o['pub'],
                        E_mc=o['E_mc'], E_quad=o['E_quad'], E_miss=o['E_miss'], E_repr=o['E_repr'], E_repr_c=o['E_repr_c'],
                        E_env3=o['E_env3'], E_signed=o['E_signed'], budget=o['budget3c'], L=o['L3c'], L_signed=o['Lsc'], L_prot=o['L3'],
                        L_prot_signed=o['Ls'], L_10=o['L10c'], probe=o['probe'], d=o['d'], sd=o['sd'], held=False))
    return out


def odd_own(J):
    st = ld('state/b5_state_J%d.json' % J); w = ld('state/b5_window_J%d.json' % J)
    out = []
    if not st or not w:
        return out
    for s in st['states']:
        d = s['direct']; mfl = s['mfl']
        o = w['states'][[k for k in w['states'] if abs(float(k) - mfl) < 1e-6][0]]
        ec = s['interp_c']['E_interp_c']
        rest_l = s['E_tot_lam'] - d['E_low'] - d['E_interp'] + ec
        mrow = max(mphys(D, J) for D in d['rows'])
        out.append(dict(J=J, mfl=mfl, route='own row (odd)', lam=o['lamc'], jk=s['jk_sigma'], pub=(s['k2'] or {}).get('lam_model'),
                        E_trunc=d['E_trunc'], E_trunc_lam=d['E_trunc_lam'], E_interp=d['E_interp'], E_interp_c=ec, E_edge=d['E_edge'],
                        E_mc=s['E_mc'], E_low_conv=d['E_low'], E_env3=o['E_env3'], budget=o['rest_c'] + o['E_env3'],
                        L=o['Lc_env3'], L_lam=o['lamc']/(rest_l + o['E_env3']), L_signed=o['Lc_signed'], L_prot=o['L_env3'],
                        L_10=o['Lc_env10'], held=mfl*LAM > mrow*(1 + 1e-4), rows=d['rows'], xmax=d['xmax'], G_edge=d['G_edge'],
                        failed=s['failed']))
    return out


def grade(r):
    if r['lam'] - r['budget'] > 0 and r['L'] >= 1:
        g = '[NE±]'
    elif r.get('L_lam', 0) >= 1:
        g = '[NE±] (λ-level)'
    elif r.get('L_signed', 0) >= 1:
        g = '[NE±] (signed window, thin)'
    else:
        g = '[NE + M]'
    if r['held'] and g.startswith('[NE±]'):
        g += ', held'
    return g


if __name__ == '__main__':
    rows = []
    for J in (18, 20, 22, 24, 26, 28, 30, 32):
        rows += even(J)
    for J in (23, 25, 31):
        rows += odd_hs(J)
    for J in (17, 19):
        rows += odd_own(J)
    extra = {f'{J}:{t}': odd_hs(J, t) for J, t in ((17, 'main'), (19, 'main'), (17, 'kb5'), (19, 'kb5'), (23, 'old'), (31, 'old'),
                                                    (23, 'x09'), (31, 'x09'))}
    extra['25:own'] = odd_own(25)
    lines = []
    for r in sorted(rows, key=lambda t: (t['J'], t['mfl'])):
        r['grade'] = grade(r)
        alt = (f"λ-level {r['L_lam']:.2f}" if 'L_lam' in r else '') + (f", signed {r['L_signed']:.2f}" if 'L_signed' in r else '')
        lines.append(f"J={r['J']:2d} m/Lam={r['mfl']:.4f} [{r['route']}] lam {r['lam']:+.4f} +- {r['jk']:.4f} (pub {r['pub'] if r['pub'] is None else round(r['pub'], 4)}) "
                     f"budget {r['budget']:.3f} L {r['L']:.2f} ({alt}; protocol {r['L_prot']:.2f}; x10 {r['L_10']:.2f}) {'HELD ' if r['held'] else ''}-> {r['grade']}")
    for k, v in extra.items():
        for r in v:
            lines.append(f"[extra {k}] J={r['J']} m/Lam={r['mfl']:.4f} lam {r['lam']:+.4f} budget {r['budget']:.3f} L {r['L']:.2f}" +
                         (f" signed {r['L_signed']:.2f}" if 'L_signed' in r else '') + f" protocol {r['L_prot']:.2f}")
    print("\n".join(lines))
    open(os.path.join(C, 'b5_sum.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(dict(rows=rows, extra=extra), open(os.path.join(C, 'b5_sum.json'), 'w'), indent=1, default=float)

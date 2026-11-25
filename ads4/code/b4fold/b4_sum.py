"""
b4_sum (BATCH4_FOLD.md §0): the per-state table at Delta_g = 2000.  b5fold/b5_sum.py with the inputs moved to b4fold/cache and
Lambda = Lambda(2000); grading rules unchanged (BATCH5_FOLD §0):
  even J (own-row route, J = 12, 16, 20): state/b4_state_J<J>.json + state/sn_even_J<J>.json.  Sign-free L = lam_c / (E_trunc +
     E_interp[lam_c] + E_edge + E_mc + window x3 envelope); lambda-level: E_trunc -> 3|lam_K2 - lam_K1|.
  odd J, hs_odd route (J = 23 [tag gk9], 27, 29; J = 21 diagnostic): odd_<tag>/b4_odd_J<J>.json (L x3 / signed, E_repr on lam_c).
  odd J, own-row route (J = 21; 23, 27, 29 cross-checks): state/b4_state_J<J>.json + state/b4_window_J<J>.json.
'held': the mass exceeds every row of its J (at Delta_g = 2000 every mass above the 2000 twist edge: no row above 2000 + J exists;
the K2-only entries of GAP_NEARGAP / GAP_SEAM / GAP_SEAM2 held the same own row).  The hs_odd route holds the own row by construction and
is not flagged, as in BATCH5_FOLD / EXT_FOLD.
OLD: the current K2-only entries at 2000 (supplement Table tab:gap2000 = GAP_NEARGAP (J = 12), GAP_SEAM (J = 16, 21; Delta_ref = 5e4),
GAP_SEAM2 (J = 27, 29; 5e4, 1000 grid u edge / exact grid)): value, sign-free L, lambda-level L.
Usage: python3 b4_sum.py   -> cache/b4_sum.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_guard  # noqa: F401
import numpy as np
C = os.path.join(HERE, 'cache')
LAM = float(np.sqrt(1998.0*1999.0))
OLD = {  # (J, m/Lam rounded to 4) -> (value, L sign-free, L lambda-level, grade) at 2000, K2 only
    (12, 1.0060): (0.4094, 0.7, 5.4, '[NE±] λ-level'), (12, 1.0119): (0.8496, 1.1, 8.0, '[NE±] λ-level'),
    (12, 1.0200): (None, 1.3, 11.6, '[NE±] λ-level'), (12, 1.0500): (None, 1.6, 18.2, '[NE±] λ-level'),
    (16, 1.0080): (0.9997, 0.4, 4.1, '[NE±] λ-level'), (16, 1.0159): (0.9399, 0.4, 5.1, '[NE±] λ-level'),
    (16, 1.0200): (None, 0.4, 4.6, '[NE±] λ-level'), (16, 1.0300): (None, 0.5, 5.1, '[NE±] λ-level'), (16, 1.0500): (None, 0.6, 6.4, '[NE±] λ-level'),
    (21, 1.0105): (0.9973, 0.1, 1.7, '[NE±] λ-level, thin'), (21, 1.0208): (0.9979, 0.2, 2.3, '[NE±] λ-level'),
    (21, 1.0300): (None, None, 2.7, '[NE±] λ-level'), (21, 1.0500): (None, None, None, '[NE±] λ-level'),
    (21, 1.0750): (None, None, None, '[NE±] λ-level'), (21, 1.1000): (None, None, 3.9, '[NE±] λ-level'),
    (27, 1.0134): (0.9670, 0.06, 0.46, '[NE + M]'), (27, 1.0175): (0.9098, 0.07, 0.56, '[NE + M]'),
    (27, 1.0267): (0.9983, 0.08, 0.77, '[NE + M]'), (27, 1.0300): (0.9984, 0.08, 0.85, '[NE + M]'),
    (27, 1.0400): (0.9920, 0.10, 1.07, '[NE±] λ-level, thin'), (27, 1.0500): (0.9620, 0.11, 1.30, '[NE±] λ-level, thin'),
    (27, 1.0750): (0.9989, 0.15, 1.62, '[NE±] λ-level, thin'), (27, 1.1000): (0.9740, 0.20, 1.79, '[NE±] λ-level, thin'),
    (29, 1.0144): (0.9143, 0.05, 0.16, '[NE + M]'), (29, 1.0150): (0.9133, 0.05, 0.16, '[NE + M]'),
    (29, 1.0286): (0.9960, 0.06, 0.25, '[NE + M]'), (29, 1.0400): (0.9673, 0.07, 0.37, '[NE + M]'),
    (29, 1.0500): (0.9836, 0.08, 0.52, '[NE + M]'), (29, 1.0600): (0.9983, 0.09, 0.70, '[NE + M]'),
    (29, 1.0750): (0.9867, 0.11, 1.01, '[NE±] λ-level, thin'), (29, 1.1000): (0.9926, 0.14, 1.42, '[NE±] λ-level, thin'),
}


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def ld(p):
    p = os.path.join(C, p)
    return json.load(open(p)) if os.path.exists(p) else None


def old(J, mfl):
    return OLD.get((J, round(mfl, 4)))


def even(J):
    st = ld('state/b4_state_J%d.json' % J); sn = ld('state/sn_even_J%d.json' % J)
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
        k2 = s['k2'] or {}
        out.append(dict(J=J, mfl=mfl, route='own row (even)', lam=lamc, lam_full=d['lam_K2'], jk=s['jk_sigma'], pub=k2.get('lam_model'),
                        k2=k2.get('lam_central'), k2_L=k2.get('L'), k2_Llam=k2.get('L_lam'),
                        E_trunc=d['E_trunc'], E_trunc_lam=d['E_trunc_lam'], E_interp=d['E_interp'], E_interp_c=ec, E_edge=d['E_edge'],
                        E_mc=s['E_mc'], E_low_conv=d['E_low'], E_env3=w['E_env3'], E_env10=w['E_env10'], E_meas=w['E_meas'],
                        budget=rest + w['E_env3'], budget_lam=rest_l + w['E_env3'], L=lamc/(rest + w['E_env3']),
                        L_lam=lamc/(rest_l + w['E_env3']), L_meas=lamc/(rest + w['E_meas']), L_10=lamc/(rest + w['E_env10']),
                        L_prot=lamc/(s['E_tot'] - d['E_low'] + w['E_env3']), L_conv=s['L'], held=held, rows=d['rows'],
                        xmax=d['xmax'], G_edge=d['G_edge'], failed=s['failed'], rz=s['rz']))
    return out


def odd_hs(J, tag='main'):
    r = ld('odd_%s/b4_odd_J%d.json' % (tag, J))
    out = []
    for k, o in sorted((r or {}).items(), key=lambda t: float(t[0])):
        out.append(dict(J=J, mfl=o['mfl'], route='hs_odd ' + tag, lam=o['lamc'], lam_full=o['main'], jk=o['jk'], pub=o['pub'],
                        E_mc=o['E_mc'], E_quad=o['E_quad'], E_miss=o['E_miss'], E_repr=o['E_repr'], E_repr_c=o['E_repr_c'],
                        E_env3=o['E_env3'], E_signed=o['E_signed'], budget=o['budget3c'], L=o['L3c'], L_signed=o['Lsc'], L_prot=o['L3'],
                        L_prot_signed=o['Ls'], L_10=o['L10c'], probe=o['probe'], d=o['d'], sd=o['sd'], held=False))
    return out


def odd_own(J):
    st = ld('state/b4_state_J%d.json' % J); w = ld('state/b4_window_J%d.json' % J)
    out = []
    if not st or not w:
        return out
    for s in st['states']:
        d = s['direct']; mfl = s['mfl']
        o = w['states'][[k for k in w['states'] if abs(float(k) - mfl) < 1e-6][0]]
        ec = s['interp_c']['E_interp_c']
        rest_l = s['E_tot_lam'] - d['E_low'] - d['E_interp'] + ec
        mrow = max(mphys(D, J) for D in d['rows'])
        k2 = s['k2'] or {}
        out.append(dict(J=J, mfl=mfl, route='own row (odd)', lam=o['lamc'], lam_full=d['lam_K2'], jk=s['jk_sigma'], pub=k2.get('lam_model'),
                        k2=k2.get('lam_central'), k2_L=k2.get('L'), k2_Llam=k2.get('L_lam'),
                        E_trunc=d['E_trunc'], E_trunc_lam=d['E_trunc_lam'], E_interp=d['E_interp'], E_interp_c=ec, E_edge=d['E_edge'],
                        E_mc=s['E_mc'], E_low_conv=d['E_low'], E_env3=o['E_env3'], budget=o['rest_c'] + o['E_env3'],
                        L=o['Lc_env3'], L_lam=o['lamc']/(rest_l + o['E_env3']), L_signed=o['Lc_signed'], L_prot=o['L_env3'],
                        L_10=o['Lc_env10'], held=mfl*LAM > mrow*(1 + 1e-4), rows=d['rows'], xmax=d['xmax'], G_edge=d['G_edge'],
                        failed=s['failed'], rz=s['rz']))
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


def fmt(r):
    o = old(r['J'], r['mfl'])
    os_ = ('old ' + ('—' if o[0] is None else f"{o[0]:+.4f}") + f" / L {o[1]} (λ {o[2]}) {o[3]}") if o else 'old: not computed'
    alt = (f"λ-level {r['L_lam']:.2f}" if 'L_lam' in r else '') + (f", signed {r['L_signed']:.2f}" if 'L_signed' in r else '')
    k2 = f" | K2 own row (2e4) {r['k2']:+.4f} L {r['k2_L']:.2f} (λ {r['k2_Llam']:.2f})" if r.get('k2') is not None else ''
    return (f"J={r['J']:2d} m/Lam={r['mfl']:.4f} [{r['route']}] lam_c {r['lam']:+.4f} (full {r['lam_full']:+.4f}) +- {r['jk']:.4f} "
            f"budget {r['budget']:.3f} L {r['L']:.2f} ({alt}; protocol {r['L_prot']:.2f}; x10 {r['L_10']:.2f}) {'HELD ' if r['held'] else ''}"
            f"-> {r['grade']} || {os_}{k2}")


if __name__ == '__main__':
    rows = []
    for J in (12, 16, 20):
        rows += even(J)
    rows += odd_hs(23, 'gk9')                 # J = 23: one GK node failed (x = 0.9139): GK edge with it in E_miss (B5_GKMIN = 9)
    for J in (27, 29):
        rows += odd_hs(J)
    rows += odd_own(21)
    extra = {f'{J}:hs': odd_hs(J) for J in (21,)}
    extra['23:espl (hs_odd default: < 10 valid GK nodes)'] = odd_hs(23)
    extra.update({f'{J}:own': odd_own(J) for J in (23, 27, 29)})
    for t in ('xb05', 'rzodd'):
        extra.update({f'{J}:{t}': odd_hs(J, t) for J in (21, 23, 27, 29)})
    lines = []
    for r in sorted(rows, key=lambda t: (t['J'], t['mfl'])):
        r['grade'] = grade(r)
        lines.append(fmt(r))
    for k, v in extra.items():
        for r in v:
            r['grade'] = grade(r)
            lines.append(f"[extra {k}] " + fmt(r))
    print("\n".join(lines))
    open(os.path.join(C, 'b4_sum.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(dict(rows=rows, extra=extra), open(os.path.join(C, 'b4_sum.json'), 'w'), indent=1, default=float)

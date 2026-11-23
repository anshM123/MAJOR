"""
ef_summary (EXT_FOLD.md §0): collect every result of this note into one table (read: extfold/cache/*.json; HIGHSPIN / SMALLNU result
files read-only for the 'old' columns).  J = 31's corrected budget (ef_odddiag (b): spline-edge tail -> sign-free envelope) is
assembled here.  Usage: python3 ef_summary.py -> cache/ef_summary.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_guard  # noqa: F401
import numpy as np
EFC = os.path.join(HERE, 'cache')
HSC = os.path.join(HERE, '..', 'highspin', 'cache'); SNC = os.path.join(HERE, '..', 'smallnu', 'cache')
ld = lambda *p: json.load(open(os.path.join(*p)))
lines, out = [], {}


def pr(s=""):
    print(s, flush=True); lines.append(s)


# ---- node checks (direct vs K2 in the bulk; join with the SMALLNU small-nu nodes)
pr("# node checks: direct (ext job A) vs K2 on the bulk nodes 0.05 <= x < 0.8 (Gamma - kappa, c2); join to SMALLNU's x < 0.05 n1 nodes")
for J in (23, 27, 31):
    R = ld(EFC, 'hs_odd_J%d.json' % J)
    nd = [n for n in R['nodes'] if 0.05 - 1e-9 <= n['x'] < 0.8 - 1e-9 and n['C2']['K2'] is not None]
    d = np.array([abs(n[k]['dir'] - n[k]['K2']) for n in nd for k in ('C2', 'C3H')])
    z = np.array([abs(n[k]['dir'] - n[k]['K2'])/n[k]['err'] for n in nd for k in ('C2', 'C3H')])
    e = np.array([n[k]['err'] for n in nd for k in ('C2', 'C3H')])
    lo = [n for n in R['nodes'] if n['x'] < 0.05 - 1e-9]; first = [n for n in R['nodes'] if abs(n['x'] - 0.05) < 1e-6][0]
    a = lo[-1]
    jz = [(first[k]['dir'] - a[k]['dir'])/np.hypot(first[k]['err'], a[k]['err']) for k in ('C2', 'C3H')]
    out[f'nodes_J{J}'] = dict(n=len(nd), med=float(np.median(d)), max=float(d.max()), medz=float(np.median(z)), maxz=float(z.max()),
                              sig_med=float(np.median(e)), join=jz)
    pr(f"  J={J}: {len(nd)} nodes x 2 kinds: |dir - K2| median {np.median(d):.0f}, max {d.max():.0f} c2; median {np.median(z):.1f} sigma, "
       f"max {z.max():.1f} sigma (direct sigma median {np.median(e):.0f}); join x={a['x']:.4f} -> 0.05: C2 {jz[0]:+.1f} sigma, C3H {jz[1]:+.1f} sigma")

# ---- odd J = 23, 27, 31: budgets
pr("\n# odd J = 23, 27, 31 (rows 1000 + J): direct bulk row (job A) + GK15 / HC edge + measured window (sn_odd, FULL normalisation)")
old_hs = {J: ld(HSC, 'hs_odd_J%d.json' % J) for J in (23, 27, 31)}
old_sn = {J: ld(SNC, 'sn_odd_J%d.json' % J) for J in (23, 27, 31)}
for J in (23, 27, 31):
    hs = ld(EFC, 'hs_odd_J%d.json' % J); sn = ld(EFC, 'sn_odd_J%d.json' % J); dg = ld(EFC, 'ef_odddiag_J%d.json' % J)
    out[f'J{J}'] = {}
    for key, s in sn['states'].items():
        h = hs['states'][key]; o_h = old_hs[J]['states'].get(key); o_s = old_sn[J]['states'].get(key); g = dg['states'][key]
        mfl = float(key)
        r = dict(mfl=mfl, main=h['main'], jk=h['jk'], K2gk=h.get('K2/gk', h.get('K2/espl')), lamc=s['lamc'], lam_signed=s['lam_signed'],
                 E_mc=h['E_mc'], E_quad=h['E_quad'], E_miss=h['E_miss'], E_repr=h['E_repr'], E_base=s['E_base'], E_env3=s['E_env3'],
                 E_signed=s['E_signed'] + s['E_lo'], L_env3=s['L_env3'], L_signed=s['L_signed'], L_env10=s['L_env10'], L_conv=s['L_conv'],
                 old_main=o_h['main'] if o_h else None, old_L_env3=o_s['L_env3'] if o_s else None, old_Lcal=o_s['Lcal_env3'] if o_s else None,
                 old_Eb=o_s['E_base'] if o_s else None, K7=g['K7'], K7_jk=g['K7_jk'], spline_bulk=g['spline_bulk'])
        # noise-aware representation alternative: degree 4/5, half samples and the model-free spline bulk (no degree-7 term)
        alt = [h['repr'][k] for k in ('K4', 'K5', 'half0', 'half1')] + [g['spline_bulk']]
        r['E_repr_alt'] = 3*max(abs(h['main'] - v) for v in alt)
        if J == 31:
            r['lamc_corr'] = g['lam_espl_cut'] - (h['main'] - s['lamc'])
            r['E_tail'] = g['E_tail']; r['E_edge_interp'] = g['E_edge_interp']
            Eb = h['E_mc'] + g['E_edge_interp'] + g['E_tail'] + h['E_repr'] + h['E_miss']
            r['E_base_corr'] = Eb
            r['L_env3_corr'] = r['lamc_corr']/(Eb + s['E_env3'])
            r['L_signed_corr'] = (r['lamc_corr'] + s['D_mid'])/(Eb + s['E_lo'] + s['E_signed'])
        else:
            Eb_alt = s['E_base'] - h['E_repr'] + r['E_repr_alt']
            r['L_env3_alt'] = s['lamc']/(Eb_alt + s['E_env3'])
        out[f'J{J}'][key] = r
        pr(f"  J={J} m/Lam={mfl:.5f}: value {h['main']:+.4f} +- {h['jk']:.4f} (K2 bulk {r['K2gk']:+.4f}; old {r['old_main'] if r['old_main'] is None else round(r['old_main'], 4)}) "
           f"| lam_c {s['lamc']:+.4f} | budget mc {h['E_mc']:.3f} quad {h['E_quad']:.3f} miss {h['E_miss']:.3f} repr {h['E_repr']:.3f} -> {s['E_base']:.3f}; "
           f"window x3 {s['E_env3']:.3f} | L x3 {s['L_env3']:.2f} (signed {s['L_signed']:.2f}, x10 {s['L_env10']:.2f}, 1e4 {s['L_conv']:.2f}) | "
           + (f"old L x3 {r['old_L_env3']:.2f} (calibrated {r['old_Lcal']:.2f})" if r['old_L_env3'] is not None else "old: mass not evaluated"))
        if J == 31:
            pr(f"        J=31 corrected edge: lam_c {r['lamc_corr']:+.4f}; budget mc {h['E_mc']:.3f} edge-interp {g['E_edge_interp']:.3f} tail {g['E_tail']:.3f} "
               f"repr {h['E_repr']:.3f} -> {Eb:.3f} + window x3 {s['E_env3']:.3f} -> L {r['L_env3_corr']:.2f} (signed {r['L_signed_corr']:.2f})")
        else:
            pr(f"        repr diagnostics: K7 {g['K7']:+.4f} (jk {g['K7_jk']:.4f}), spline bulk {g['spline_bulk']:+.4f}; noise-aware repr {r['E_repr_alt']:.3f} "
               f"-> L x3 {r['L_env3_alt']:.2f}")

# ---- J = 29
pr("\n# J = 29 (row 1029): window with the pooled n3 nodes (31 distinct streams)")
new = ld(EFC, 'sn_odd_J29.json'); old = ld(EFC, 'old_repro', 'sn_odd_J29.json'); env = ld(EFC, 'ef_j29env.json')
out['J29'] = {}
for key, s in new['states'].items():
    o = old['states'][key]; e = env[key]
    out['J29'][key] = dict(lamc=s['lamc'], lam_signed=s['lam_signed'], L_env3_old=o['L_env3'], L_env3=s['L_env3'], L_signed_old=o['L_signed'],
                           L_signed=s['L_signed'], E_signed_old=o['E_signed'], E_signed=s['E_signed'], S3=e['S3 nu=10-17.5 at |Gamma|']['L'],
                           S4=e['S4 S3 + nu=20 at |Gamma|']['L'], S5=e['S5 S1 without the nu=20 n1 node']['L'], E_base=s['E_base'], E_env3=s['E_env3'], E_lo=s['E_lo'])
    pr(f"  m/Lam={float(key):.5f}: sign-free x3 L {o['L_env3']:.2f} -> {s['L_env3']:.2f} (S3 {e['S3 nu=10-17.5 at |Gamma|']['L']:.2f}, "
       f"S4 {e['S4 S3 + nu=20 at |Gamma|']['L']:.2f}, S5 {e['S5 S1 without the nu=20 n1 node']['L']:.2f}) | signed L {o['L_signed']:.2f} -> {s['L_signed']:.2f} "
       f"(window MC/interp error {o['E_signed']:.4f} -> {s['E_signed']:.4f}); lam signed {s['lam_signed']:+.4f}")

# ---- J = 7
pr("\n# J = 7 (own direct row 1007)")
J7 = ld(EFC, 'ef_j7.json'); out['J7'] = J7
for o in J7:
    pr(f"  m/Lam={o['mfl']:.5f}: published K2 chain {o['k2']['lam_K2']:+.5f} (NEARGAP2-style budget on the published rows: L {o['k2']['L']:.1f} / "
       f"{o['k2']['L_lam']:.1f}) -> own direct row {o['direct']['lam_central']:+.5f} +- {o['jk_sigma']:.5f}; budget {o['E_tot']:.4f} "
       f"(edge {o['direct']['E_edge']:.4f}, low {o['direct']['E_low']:.4f}, trunc {o['direct']['E_trunc']:.4f}, mc {o['E_mc']:.5f}); "
       f"L = {o['L']:.1f} sign-free ({o['L_lam']:.1f} lambda-level); certified-floor route {o['direct']['base']['+-']['rig']:+.4f}")
open(os.path.join(EFC, 'ef_summary.txt'), 'w').write("\n".join(lines) + "\n")
json.dump(out, open(os.path.join(EFC, 'ef_summary.json'), 'w'), indent=1, default=float)

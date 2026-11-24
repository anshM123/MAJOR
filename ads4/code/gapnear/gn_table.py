"""
gn_table (GAP_NEARGAP.md §0): the headline table, per state (J, m/Lam):
  Delta_g = 1000: the GRADED value / L of the error-budgeted protocol (NEARGAP2 / HIGHSPIN / SMALLNU_WINDOW / EXT_FOLD; quoted, not
                  recomputed) and the same-protocol K2 value / L of this note's rerun (cache/gn_bind_base1000.json);
  Delta_g = 2000 / 5000: own-row K2 value, sign-free L and lambda-level L (cache/gn_bind_<tags>.json).
Usage: python3 gn_table.py tag[,tag...] > cache/table.txt
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, 'cache')
# graded Delta_g = 1000 values (twist edge unless m/Lam given): (value, L, source)
G1000 = {0: ('+1.0000', '4362', 'NEARGAP2 own row 990'), 2: ('+1.0000', '913', 'NEARGAP2 own row 990'),
         4: ('+0.2041', '22 (λ 80)', 'NEARGAP2'), 5: ('+0.2009', '78 (λ 130)', 'NEARGAP2'), 6: ('+0.2295', '5.0 (λ 28)', 'NEARGAP2'),
         7: ('+0.2683', '6.7', 'EXT_FOLD direct row'), 8: ('+0.4297', '1.2 (λ 4.6)', 'NEARGAP2'),
         12: ('+0.8528', '7.4', 'HIGHSPIN direct row + SMALLNU'), 16: ('+0.9402', '4.4', 'HIGHSPIN direct row + SMALLNU'),
         21: ('+0.998 ± 0.020', '1.20', 'SMALLNU direct bulk'), 27: ('+0.976 ± 0.019', '1.64', 'EXT_FOLD direct bulk'),
         29: ('+0.972 ± 0.018', '1.24 signed / 0.94 ×3', 'HIGHSPIN + SMALLNU + EXT_FOLD')}


def load(tag):
    return json.load(open(os.path.join(C, 'gn_bind_%s.json' % tag)))


base = {(r['J'], round(r['mfl'], 5)): r for r in load('base1000')}
rows = []
for t in sys.argv[1].split(','):
    rows += load(t)
print("| state (J, m/Λ) | Δ_g = 1000 graded value / L | 1000 same-protocol K₂ value / L (λ-level) | 2000 own-row K₂ value / L (λ-level) | 5000 / L (λ-level) |")
print("|---|---|---|---|---|")
for r in sorted(rows, key=lambda r: (r['J'], r['mfl'], r['Dg'])):
    if r['Dg'] == 5000:
        continue
    k = (r['J'], round(r['mfl'], 5)); b = base.get(k)
    g = G1000[r['J']] if (b is not None and r.get('tok') in ('tw1000',) or (r['J'] in (0, 2))) else None
    five = [s for s in rows if s['Dg'] == 5000 and (s['J'], round(s['mfl'], 5)) == k]
    f = f"{five[0]['lam_K2']:+.4f} / {five[0]['L']:.1f} ({five[0]['L_lam']:.1f})" if five else '—'
    print(f"| {r['J']}, {r['mfl']:.5f}{' (2000 twist edge)' if r.get('tok') == 'tw' else ''} | "
          f"{(g[0] + ' / ' + g[1]) if g else '—'} | "
          f"{(format(b['lam_K2'], '+.4f') + ' / ' + format(b['L'], '.1f') + ' (' + format(b['L_lam'], '.1f') + ')') if b else '—'} | "
          f"{r['lam_K2']:+.4f} / {r['L']:.1f} ({r['L_lam']:.1f}) | {f} |")
for s in rows:
    if s['Dg'] == 5000 and not any(r['Dg'] != 5000 and (r['J'], round(r['mfl'], 5)) == (s['J'], round(s['mfl'], 5)) for r in rows):
        print(f"| {s['J']}, {s['mfl']:.5f}{' (5000 twist edge)' if s.get('tok') == 'tw' else ''} | — | — | — | "
              f"{s['lam_K2']:+.4f} / {s['L']:.1f} ({s['L_lam']:.1f}) |")

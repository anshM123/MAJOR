"""
gsm_table (GAP_SEAM.md §0): the seam / near-gap table from gsm_bind JSONs: per state (J, m/Lam) the Delta_g = 1000 same-protocol value /
L (lambda-level L), the 2000 own-row value / L, the grade at 2000, and the chart-A ratios of the budget terms (2000 / 1000, expected
(Lam_1000/Lam_2000)^2 = 0.2496 for c2-type terms, ~0.062 for E_low) where the 1000 run exists at the same m/Lam.
Grades: [NE±] (sign-free budget closes, L > 1), [NE±] lam (only the lambda-level budget closes), [NE + M] otherwise.
OO rows (even J >= 4): the ++-diagonal 1x1 value and budget (Gamma_OO = Gamma_++).
Usage: python3 gsm_table.py BASE1000_TAG TAG2000[,TAG2000...]   -> stdout
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))


def load(tag):
    return json.load(open(os.path.join(HERE, 'cache', 'gsm_bind_%s.json' % tag)))


def grade(L, Ll):
    return '[NE±]' if L > 1 else ('[NE±] λ-level' if Ll > 1 else '[NE + M]')


def key(r):
    return (r['J'], round(r['mfl'], 5))


def fmt(r, oo=False):
    if r is None:
        return '—'
    s = r['OO'] if oo else r
    return f"{s['lam_central']:+.4f} / {s['L']:.1f} ({s['L_lam']:.1f})"


if __name__ == '__main__':
    b = {}
    for t in sys.argv[1].split(','):
        for r in load(t):
            b[key(r)] = r
    for oo in (False, True):
        print(('## OO block (++ diagonal)' if oo else '## EE block') + "\n| J | m/Λ (tok) | 1000: value / L (λ) | 2000: value / L (λ) | grade 2000 | "
              "ratio E_trunc, E_trunc^λ, E_low, E_edge (2000/1000) |\n|---|---|---|---|---|---|")
        for t in sys.argv[2].split(','):
            for r in load(t):
                if oo and 'OO' not in r:
                    continue
                r1 = b.get(key(r))
                if oo and r1 is not None and 'OO' not in r1:
                    r1 = None
                s, s1 = (r['OO'], r1['OO'] if r1 else None) if oo else (r, r1)
                rat = '—'
                if s1 is not None:
                    rat = ", ".join(f"{s[k]/s1[k]:.3f}" if s1[k] > 0 else '—' for k in ('E_trunc', 'E_trunc_lam', 'E_low', 'E_edge'))
                print(f"| {r['J']} | {r['mfl']:.5f} ({r['tok']}) | {fmt(r1, oo)} | {fmt(r, oo)} | {grade(s['L'], s['L_lam'])} | {rat} |")

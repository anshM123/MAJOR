"""
nf_reprsum (NEARGAP_FILL.md §3): the odd-J direct-bulk budgets (J = 23, 27, 31) with the representation term evaluated on the quantity
the budget bounds.  Inputs (this note's cache, all from unchanged project code): repr_ols/hs_odd_J<J>.json (hs_odd lambda-level terms;
reproduces EXT_FOLD), repr_ols/sn_odd_J<J>.json (lam_c = value with the nu < 50 window REMOVED, window envelope x3 / signed window),
nf_reprcut_J<J>.json (every representation probe on lam and on lam_c, with jackknife sigmas).
  protocol : L = lam_c / (E_mc + E_quad + E_miss + E_repr[lam] + window)                 (EXT_FOLD's numbers)
  lam_c    : same, with E_repr[lam_c] = 3 max_probes |lam_c(probe) - lam_c(K6)|  (the probes' polynomial continuation below the first
             bulk node, x < 50/m, is not part of lam_c: that window is covered by the measured envelope / signed window)
Usage: python3 nf_reprsum.py -> cache/nf_reprsum.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_guard  # noqa: F401
C = os.path.join(HERE, 'cache')
PROBES = ('K4', 'K5', 'K7', 'half0', 'half1')
lines, res = [], {}


def pr(s=""):
    print(s, flush=True); lines.append(s)


pr("| J | m/Λ | λ_c | mc / quad / miss | E_repr on λ (protocol) → on λ_c | window ×3 / signed | L ×3: protocol → λ_c | L signed: protocol → λ_c | largest λ_c probe (Δ ± jk) |")
pr("|---|---|---|---|---|---|---|---|---|")
for J in (23, 27, 31):
    hs = json.load(open(os.path.join(C, 'repr_ols', 'hs_odd_J%d.json' % J)))['states']
    sn = json.load(open(os.path.join(C, 'repr_ols', 'sn_odd_J%d.json' % J)))['states']
    rc = json.load(open(os.path.join(C, 'nf_reprcut_J%d.json' % J)))
    for key in sorted(hs, key=float):
        h = hs[key]; s = sn[[k for k in sn if abs(float(k) - float(key)) < 1e-6][0]]
        r = rc[[k for k in rc if abs(float(k) - float(key)) < 1e-6][0]]
        Ec = 3*max(abs(r[p]['dc']) for p in PROBES); pk = max(PROBES, key=lambda p: abs(r[p]['dc']))
        base_c = h['E_tot'] - h['E_repr'] + Ec
        L3, Ls = s['lamc']/(h['E_tot'] + s['E_env3']), s['lam_signed']/(h['E_tot'] + s['E_signed'] + s['E_lo'])
        L3c, Lsc = s['lamc']/(base_c + s['E_env3']), s['lam_signed']/(base_c + s['E_signed'] + s['E_lo'])
        res[f'{J}:{key}'] = dict(J=J, mfl=float(key), lamc=s['lamc'], E_mc=h['E_mc'], E_quad=h['E_quad'], E_miss=h['E_miss'], E_repr=h['E_repr'],
                                E_repr_c=Ec, E_env3=s['E_env3'], E_signed=s['E_signed'] + s['E_lo'], L3=L3, L3c=L3c, Ls=Ls, Lsc=Lsc, probe=pk,
                                d=r[pk]['dc'], sd=r[pk]['sdc'], lower3c=s['lamc'] - base_c - s['E_env3'])
        pr(f"| {J} | {float(key):.4f} | {s['lamc']:+.3f} | {h['E_mc']:.3f} / {h['E_quad']:.3f} / {h['E_miss']:.3f} | {h['E_repr']:.3f} → {Ec:.3f} | "
           f"{s['E_env3']:.3f} / {s['E_signed'] + s['E_lo']:.3f} | {L3:.2f} → **{L3c:.2f}** | {Ls:.2f} → {Lsc:.2f} | {pk} {r[pk]['dc']:+.3f} ± {r[pk]['sdc']:.3f} |")
open(os.path.join(C, 'nf_reprsum.txt'), 'w').write("\n".join(lines) + "\n")
json.dump(res, open(os.path.join(C, 'nf_reprsum.json'), 'w'), indent=1)

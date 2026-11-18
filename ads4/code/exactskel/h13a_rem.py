"""
h13a_rem: exact scalar-skeleton remainder for the ROUTE-A' D-safe functional on all exact large-J states of h13_exact.py
(cache/h13_ex_<Dg>.pkl; rho_k are functional-independent, so nothing is recomputed).  F = F_hyb^{A'} + S^{A'}, S^{A'} from
h13a_common.remA (variants: 'strip' = stripped envelope, f6's headline; 'Bp' = B' envelope).  Reliability: n-sum spread of S^{A'}
between rho and rho_alt (another (N0, nq, L) configuration) <= 1e-3 N_s, and the old-functional spread <= 1e-3 N_s.
Output: cache/h13a_rem.txt, cache/h13a_rem_table.md, cache/h13a_rem_recs.pkl
"""
import os, sys, glob, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from h13_common import flat_rows, Lam_of, CACHE
from common import gl01
from h13a_common import hybridA, remA

out, md, allrecs = [], [], []
def pr(s=''):
    print(s, flush=True); out.append(s)

pr("# h13a_rem: route-A' exact skeleton remainder on the h13_exact states (both skeletons, both channels)")
for Dg in (100, 300, 1000):
    fn = os.path.join(CACHE, 'h13_ex_%d.pkl' % Dg)
    db = pickle.load(open(fn, 'rb'))
    recs = []
    for key, r in db.items():
        npn = r.get('NP', 64); pn, pw = gl01(npn)
        olderr = max(abs(v['S'] - v.get('Salt', v['S']))/v['Ns'] for v in r['hel'].values())
        for h in r['hel']:
            B2, B3 = flat_rows(h, 1/r['mfl']**2, [r['J']])[r['J']]
            hy = hybridA(h, Dg, r['mfl'], r['J'], B2, B3, beta=r['beta'])
            S, Sabs = remA(h, r['mfl'], r['rho'], pn, pw, 'strip')
            Sb, _ = remA(h, r['mfl'], r['rho'], pn, pw, 'Bp')
            Sa, _ = remA(h, r['mfl'], r['rho_alt'], pn, pw, 'strip') if 'rho_alt' in r else (S, 0)
            S3, _ = remA(h, r['mfl'], r['rho3'], pn, pw, 'strip')
            err = abs(S - Sa)/hy['Ns']
            rec = dict(Dg=Dg, Df=r['Df'], h=h, mfl=r['mfl'], J=r['J'], bh=r['bh'], beta=r['beta'], m=r['m'], Fh=hy['Fh'], Ns=hy['Ns'],
                       S=S, SBp=Sb, Sabs=Sabs, err=err, olderr=olderr, odd=abs(S - S3)/hy['Ns'], ok=(err <= 1e-3 and olderr <= 1e-3))
            rec['F'] = hy['Fh'] + S
            recs.append(rec); allrecs.append(rec)
    good = [q for q in recs if q['ok']]
    w = min(good, key=lambda q: q['F']/q['Ns'])
    pr(f"## Dg = {Dg}: {len(recs)} helicity records, reliable {len(good)}; non-positive (reliable): {sum(q['F'] <= 0 for q in good)}")
    pr(f"   min F/N_s = {w['F']/w['Ns']:+.5f} ({w['h']}, Df={w['Df']:g}, m/Lam={w['mfl']:.3f}, J={w['J']}, bhat={w['bh']:.1f});"
       f"  max |S|/F_hyb = {max(abs(q['S']/q['Fh']) for q in good):.2e} (strip), {max(abs(q['SBp']/q['Fh']) for q in good):.2e} (B' env);"
       f"  max |S| m^2/F_hyb = {max(abs(q['S']/q['Fh'])*q['m']**2 for q in good):.1f};  max sum|.|/N_s = {max(q['Sabs']/q['Ns'] for q in good):.1e};"
       f"  max odd-k sens/N_s = {max(q['odd'] for q in good):.1e};  bhat range reliable [{min(q['bh'] for q in good):.1f}, {max(q['bh'] for q in good):.1f}]")
    md.append(f"\n**Δ_g = {Dg}**\n")
    md.append("| `m/Λ` | `J` | `b̂` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/F_hyb` (strip) | max `|S|/F_hyb` (B′ env) | status |")
    md.append("|---|---|---|---|---|---|---|---|")
    for (mflt, J) in sorted(set((round(q['mfl'], 3), q['J']) for q in recs), key=lambda a: (a[1], a[0])):
        v = [q for q in recs if round(q['mfl'], 3) == mflt and q['J'] == J]
        g = [q for q in v if q['ok']]
        cell = lambda h: (f"{min(q['F']/q['Ns'] for q in g if q['h'] == h):+.4f}" if any(q['h'] == h for q in g) else '—')
        dfs = sorted(set(q['Df'] for q in g), reverse=True)
        st = ('ok (Δφ=' + ','.join(f'{d:g}' for d in dfs) + ')') if g else '**unreliable** (n-sum)'
        md.append(f"| {mflt:.2f} | {J} | {v[0]['bh']:.1f} | {cell('++')} | {cell('+-')} | "
                  + (f"{max(abs(q['S']/q['Fh']) for q in g):.1e} | {max(abs(q['SBp']/q['Fh']) for q in g):.1e}" if g else "— | —") + f" | {st} |")
pickle.dump(allrecs, open(os.path.join(CACHE, 'h13a_rem_recs.pkl'), 'wb'))
open(os.path.join(CACHE, 'h13a_rem.txt'), 'w').write('\n'.join(out) + '\n')
open(os.path.join(CACHE, 'h13a_rem_table.md'), 'w').write('\n'.join(md) + '\n')

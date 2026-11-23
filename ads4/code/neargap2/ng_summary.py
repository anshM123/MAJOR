"""ng_summary (NEARGAP2.md): compact table of every ng_bind_*.json and the state lines of every ng_val_*.txt -> cache/ng_summary.txt"""
import os, glob, json
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache')
L = ["tag J m/Lam | hybrid model own K1 | E_trunc E_trunc^lam E_interp E_low E_edge | lower(sf) lower(lam) | L(sf) L(lam)"]
for f in sorted(glob.glob(os.path.join(C, 'ng_bind_*.json'))):
    for r in json.load(open(f)):
        L.append("%-10s %2d %.5f | %+.5f %+.5f %+.5f %+.5f | %.1e %.1e %.1e %.1e %.1e | %+.4f %+.4f | %.1f %.1f" % (
            os.path.basename(f)[8:-5], r['J'], r['mfl'], r['lam_hybrid'], r['lam_model'], r['lam_K2'], r.get('lam_K1', float('nan')),
            r['E_trunc'], r.get('E_trunc_lam', float('nan')), r['E_interp'], r['E_low'], r['E_edge'], r['lower'],
            r.get('lower_lam', float('nan')), r['L'], r.get('L_lam', float('nan'))))
for f in sorted(glob.glob(os.path.join(C, 'ng_val_J*_*.txt'))):
    L += [os.path.basename(f)] + [l for l in open(f) if 'state m/Lam' in l or 'normalisation' in l]
open(os.path.join(C, 'ng_summary.txt'), 'w').write("\n".join(x.rstrip('\n') for x in L) + "\n")
print("\n".join(x.rstrip('\n') for x in L))

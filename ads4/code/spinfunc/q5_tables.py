"""q5 (S1_S4_FIXES.md §2): tables from q3_analysis.run -- per (J, m/Lam) min eigenvalue lambda/N_s with and without the computed
graviton correction, and the off-diagonal ratio |dF_+-|/sqrt(N N); env as q3 (Q3_*).  Usage: python3 q5_tables.py Dg -> stdout + cache"""
import os, sys, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import q3_analysis as Q

if __name__ == '__main__':
    Dg = int(sys.argv[1]); Jmin = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    grids, out = Q.run([Dg], pr)
    recs = out[Dg]['recs']
    ms = sorted(set(r['mfl'] for r in recs if r['J'] >= Jmin))
    pr(f"# per-state min eigenvalue lambda/N_s (corrected | uncorrected | off/sqrt(NN)), Delta_g = {Dg}; Df = 3 or hyb (worst of the two)")
    pr("  m/Lam : " + " ".join(f"{m:>17.2f}" for m in ms))
    for J in sorted(set(r['J'] for r in recs if r['J'] >= Jmin)):
        row = []
        for m in ms:
            rr = [r for r in recs if r['J'] == J and abs(r['mfl'] - m) < 1e-9 and r['DF'] in (3.0, 'hyb')]
            if not rr:
                row.append(f"{'':>17s}"); continue
            r = min(rr, key=lambda r: Q.lam(r, 1.0))
            of = abs(r['dF'][0, 1])/np.sqrt(r['Ns'].prod()) if len(r['hs']) == 2 else 0.0
            row.append(f"{Q.lam(r, 1.0):+.3f}|{Q.lam(r, 0.0):+.2f}|{of:.2f}")
        pr(f"  J={J:2d}  : " + " ".join(f"{c:>17s}" for c in row))
    open(os.path.join(HERE, 'cache', 'q5_tables_%d%s.txt' % (Dg, os.environ.get('Q3_TAG', ''))), 'w').write(buf.getvalue())

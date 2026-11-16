"""Collect the t3 outputs for both skeleton choices into one table (cache/out_summary.txt)."""
import re, os
HERE = os.path.dirname(os.path.abspath(__file__))
L = []
for tag, lab in (('', 'Df=3'), ('_Df2', 'Df=2')):
    for Dg in (20, 50, 100, 300, 1000):
        fn = os.path.join(HERE, 'cache', f'out_t3_{Dg}{tag}.txt')
        if not os.path.exists(fn):
            continue
        t = open(fn).read().replace('np.float64(', '(')
        for n in ('M02', 'M01', 'Fg3'):
            for h in ('++', '+-'):
                m = re.search(rf"  {n}\s+\{h[0]}\{h[1]}\s+states.*?min F_hyb/N = (\S+).*?\n\s+Cs=   0: .*?= (\S+) at.*?kappa_crit = (\S+).*?\n\s+Cs=  40: .*?= (\S+) at.*?kappa_crit = (\S+).*?\n\s+C1_eff.*?= (\S+) at", t, re.S)
                bands = re.findall(rf"  {n} \{h[0]}\{h[1]} m/Lam in \[(\S+),(\S+)\): E/N <= (\S+)\s+U/N <= (\S+)\s+min\(F-E\)/N = (\S+)\s+kappa_crit = (\S+)\s+model: (\S+)", t)
                nm = re.search(rf"  {n} \{h[0]}\{h[1]}: (.*)", t)
                L.append(f"{lab} Dg={Dg:4d} {n} {h}: Fhyb/N={m.group(1)} ctrl-min={m.group(2)} kap={m.group(3)} (Cs40: {m.group(4)}, {m.group(5)}) C1glob={m.group(6)}")
                for b in bands:
                    L.append(f"      band [{b[0]},{b[1]}): E/N<={b[2]} U/N<={b[3]} min(F-E)/N={b[4]} kappa={b[5]} model={b[6]}")
                if nm:
                    L.append("      near-gap model: " + nm.group(1)[:120])
open(os.path.join(HERE, 'cache', 'out_summary.txt'), 'w').write("\n".join(L) + "\n")
print("\n".join(L))

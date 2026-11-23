"""census_table: one-line-per-(Delta, nu) summary of the pcyc_census runs (counts per raw-radius class, max dH, (Delta/2)log|B|)."""
import os, sys, glob, re
HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for f in sorted(glob.glob(os.path.join(HERE, 'cache', 'census_D*_s*.txt'))):
    L = open(f).read().splitlines()
    h = L[0]
    m = re.search(r'Delta=(\S+) nu=(\S+) starts=(\d+) seed=(\d+): (\d+) converged, (\d+) distinct', h)
    if not m:
        continue
    D, nu, st, sd, nc, nd = m.groups()
    cl = []
    for ln in L[1:6]:
        mm = re.search(r':\s+(\d+) points; max dH\s+(\S+); \(Delta/2\)log\|B\| in \[(\S+), (\S+)\]', ln)
        if mm:
            cl.append(f'{mm.group(1):>4s} ({mm.group(2)})')
    rows.append((float(D), float(nu), int(sd), f'{float(D):6g} {float(nu):4g} {st:>5s} {sd:>3s} {nd:>5s} | ' + ' | '.join(cl)))
rows.sort()
print('# Delta  nu  starts seed distinct | Regge r<0.1 | mid 0.1-0.3 | FAR 0.3-1.5 | beyond 1.5-10 | raw inf >=10   (count (max dH))')
for r in rows:
    print(r[3])

import re, sys, math
sys.path.insert(0, '.')
import common2, rt_cert as RC
La, Lb = float(sys.argv[1]), float(sys.argv[2])
d = math.exp(-La); F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb)); fac = d**3*(2/RC.LAMf)*2/F
best = {}
for p in (0, 1):
    for line in open(f'logs/resid_{La}_{Lb}_p{p}.log'):
        m = re.match(r'(\w+) \[([\d.]+),([\d.]+)\] m=([\de.+-]+): sup\|K_J\| <= ([\de.+-]+)', line)
        if m:
            k = float(m.group(2)); v = fac*float(m.group(5))*float(m.group(4))
            best[k] = max(best.get(k, 0), v)
tot = sum(best.values())
print('cells', len(best), 'partial E_r/F', tot)
for k, v in sorted(best.items(), key=lambda x: -x[1])[:8]: print(k, v)

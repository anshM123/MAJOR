"""rv_cell1031.py (review): per-unit-cell residual contributions at L = 5 near the turning point, split into multiplier / j-sum parts,
and the same with the exact per-power multiplier at (n, a) = (295234, 998) for comparison (DV)."""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rt_cert as RC, rt_pp as PPm
from flint import arb, acb
L = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
d = math.exp(-L)
F = RC.fcert_min(RC.beta_of_L(L), RC.beta_of_L(L + 0.5))
cells = json.load(open(os.path.join(RC.HERE, 'cache', 'rt_rmass2.json')))
pref = d**3*(2/RC.LAMf)
tot = 0
for (T0, T1, m) in cells:
    if T1 > 1e8 or T0 < 1020 or T0 > 1045: continue
    K = RC.K_sup_real(float(T0), float(T1), L, 2)
    c = 2*K*m*pref/F
    tot += c
    pp = PPm.PP(2.5, T0, T1, d*1.0001)
    print(f'cell [{T0},{T1}]: mass {m:.3e}  sup|K_J| <= {K:.3e}  PP(2.5) {pp:.3e}  contribution/F = {c:.3e}', flush=True)
print('sum over [1020,1045] cells /F =', tot, ' F =', F)

"""t0check/chk_scan.py (INDEP_CHECK_T0 §5): hostile float scan [DV] of lambda = H/N_s with the independent evaluator chk_eval.py,
at Delta_g = 300, 1000, 5000 (argv[1]); fresh random grid (seed = argv[2]) incl. near-forward (bhat -> 0), the twist edge, m = Lam
(u = 1), huge J (to 1e5) and large bhat (to ~ 2e4, beyond 2 Lam).  Output: cache/chk_scan_<Dg>.txt"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import chk_eval as E

Dg = float(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
rng = np.random.default_rng(seed)
L = E.lam(Dg)
def mgap(J):
    al = Dg - 2.0
    return np.sqrt(al*(al + 2*J + 1))/L

pts = []
# 1. every small J, u = 1 (m = Lam), the twist edge, and random m
for h, Js in (('++', range(0, 61, 2)), ('+-', range(4, 61))):
    for J in Js:
        for mf in [1.0, mgap(J)] + list(np.exp(rng.uniform(0, np.log(200), 6))):
            pts.append((h, J, mf, 'smallJ'))
# 2. near-forward: bhat in [1e-3, 2]
for _ in range(150):
    h = rng.choice(['++', '+-']); J = int(rng.integers(0, 40)) if h == '++' else int(rng.integers(4, 40))
    if h == '++': J -= J % 2
    bh = np.exp(rng.uniform(np.log(1e-3), np.log(2)))
    pts.append((h, J, (2*J + 1)/bh, 'fwd'))
# 3. medium/large J, random bhat (log in [0.05, 2e4]) with m >= Lam
for _ in range(500):
    h = rng.choice(['++', '+-']); J = int(np.exp(rng.uniform(np.log(60), np.log(3000))))
    if h == '++': J -= J % 2
    bh = np.exp(rng.uniform(np.log(0.05), np.log(min(2e4, 2*J + 1))))
    pts.append((h, J, max(1.0, (2*J + 1)/bh), 'midJ'))
# 4. huge J (1e4 .. 1e5), bhat <= 300
for _ in range(60):
    h = rng.choice(['++', '+-']); J = int(np.exp(rng.uniform(np.log(1e4), np.log(1e5))))
    if h == '++': J -= J % 2
    bh = np.exp(rng.uniform(np.log(0.05), np.log(300)))
    pts.append((h, J, (2*J + 1)/bh, 'hugeJ'))
# 5. twist edge, random J to 3000
for _ in range(120):
    h = rng.choice(['++', '+-']); J = int(np.exp(rng.uniform(np.log(4), np.log(3000))))
    if h == '++': J -= J % 2
    pts.append((h, J, mgap(J), 'edge'))

t0 = time.time(); res = []
for h, J, mf, tag in pts:
    d = E.hyb(Dg, J, mf, h)
    lam_ = d['base']/d['Ns']
    res.append((lam_, h, J, mf, d['bhat'], tag, d['B2'] + d['B3'], abs(d['B2']) + abs(d['B3'])))
res.sort()
out = [f'Delta_g = {Dg}, Lam = {L:.3f}, seed {seed}: {len(res)} states, {time.time() - t0:.0f}s']
by = {}
for r in res:
    by.setdefault(r[5], r)
out.append('min H/N_s per family: ' + '; '.join(f'{k}: {v[0]:.9f} ({v[1]}, J={v[2]}, m/Lam={v[3]:.4g}, bhat={v[4]:.4g})' for k, v in by.items()))
out.append('lowest 8 states (lambda, h, J, m/Lam, bhat, family, flat ratio F/N):')
for r in res[:8]:
    out.append(f'   {r[0]:.10f}  {r[1]} J={r[2]} m/Lam={r[3]:.6g} bhat={r[4]:.5g} {r[5]}  flat F/N = {r[6]/r[7]:.10f}')
nb = sum(1 for r in res if r[0] < 0.1999)
out.append(f'states with H/N_s < 0.1999: {nb};  < 0.199: {sum(1 for r in res if r[0] < 0.199)}')
print('\n'.join(out))
open(os.path.join(HERE, 'cache', f'chk_scan_{int(Dg)}_{seed}.txt'), 'w').write('\n'.join(out) + '\n')

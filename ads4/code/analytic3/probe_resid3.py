"""DV estimate of the residual budget per band at L-cell [La, Lb] using rt2_kernel.K_cell on sample cells and the EPSW per-cell masses."""
import sys, os, math, json, time
sys.path.insert(0, '../analytic2'); sys.path.insert(0, '../analytic')
import common2, rt2_kernel as RK, rt_cert as RC, kernel3 as K3
La, Lb = float(sys.argv[1]), float(sys.argv[2]); jm = int(sys.argv[3]) if len(sys.argv) > 3 else 110
d = math.exp(-La)
F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
fac = d**3*(2/RC.LAMf)*2/F
masses = {}
for fn in ('../analytic2/cache/rmass3_1_40.json', '../analytic2/cache/rmass3_40_898.json', '../analytic/cache/rt_rmass2.json'):
    for (T0, T1, m) in json.load(open(fn)):
        if T1 < 1e8: masses[int(T0)] = m
print('L', La, 'F', F, 'fac', fac, flush=True)
tot = 0.0; prevK = None; prevT = None
samples = [1, 5, 20, 40, 100, 200, 300, 400, 500, 600, 700, 800, 850, 900, 950, 990, 1000, 1010, 1020, 1031, 1040, 1050]
Ks = {}
for T in samples:
    t0 = time.time()
    v, ar = K3.K_cell(float(T), T + 1.1, La, Lb, ratio=float(sys.argv[4]) if len(sys.argv) > 4 else 1.03, jmax=jm)
    Ks[T] = v
    print(f'  nu {T}: K <= {v:.3e} (a {ar})  K*m*fac = {v*masses.get(T,0)*fac:.3e}  [{time.time()-t0:.0f}s]', flush=True)
# crude band estimate: piecewise geometric interpolation of K between samples
est = 0.0
for T in range(1, 1059):
    lo = max(s for s in samples if s <= T); hi = min([s for s in samples if s >= T] or [samples[-1]])
    if hi == lo: K = Ks[lo]
    else: K = math.exp(math.log(Ks[lo]) + (T - lo)/(hi - lo)*(math.log(Ks[hi]) - math.log(Ks[lo])))
    est += K*masses.get(T, 0)*fac
print('estimated residual/F (cells >= 1):', est)

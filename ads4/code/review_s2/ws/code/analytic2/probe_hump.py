import sys, math, time
sys.path.insert(0, '.')
import rt2_kernel as R
import rt_cert as RC
La = float(sys.argv[1]); Lb = La + 0.5
d = math.exp(-La)
F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
fac = d**3*(2/RC.LAMf)*2/F
print('L', La, 'F', F, 'per unit (K*mass):', fac, flush=True)
for (T0, T1) in [(40, 42), (200, 202), (400, 402), (600, 602), (800, 802), (896, 898), (950, 951.1), (1000, 1001.1), (1031, 1032.1)]:
    t0 = time.time()
    r = R.K_cell(T0, T1, La, Lb)
    print(f'  nu in [{T0},{T1}]: K <= {r[0]:.3e} at a {r[1]}  -> per 1e-9 mass: {fac*r[0]*1e-9:.3e}   [{time.time()-t0:.0f}s]', flush=True)

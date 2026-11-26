import sys, math, time
sys.path.insert(0, '../analytic2'); sys.path.insert(0, '../analytic')
import rt2_kernel as RK
T0, La, Lb = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
T1 = T0 + 1.1
def ranges(ratio, w0):
    out = []; a = 998.0
    while a < 2.0e5:
        w = max(w0, a*(ratio - 1)) if a > T1 + 30 else 1.0
        out.append((a, a + w)); a += w
    out.append((a, math.inf)); return out
for (ratio, w0) in ((1.15, 2.0), (1.03, 1.0), (1.01, 1.0)):
    t0 = time.time()
    R = ranges(ratio, w0)
    v = RK.K_cell(T0, T1, La, Lb, ranges=R, jmax=110)
    print(ratio, w0, len(R), v, f'{time.time()-t0:.0f}s', flush=True)

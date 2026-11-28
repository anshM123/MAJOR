"""Validation: disc HRT growth in thin-shell Vaidya-AdS4 (M=1, z_h=1).
Liu-Suh linear regime: d(A/2pi)/dt -> v_E * R / z_h^2, v_E = sqrt(3)/2^(4/3) = 0.6874
(for R >> t >> z_h); saturation value: static black-brane disc."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from vaidya4 import Vaidya
from hrt_pieces import Pieces
V = Vaidya(M=1.0, dlt=0.01)
P = Pieces(V, cap_grid=(260, 220))
P.cap_scan(vrange=(-9.0, 9.0), zrange=(0.02, 7.0))
R = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
ts = np.arange(0.0, 2 * R + 1.01, 0.25)
prev = None
for t in ts:
    s = P.caps(R, t)
    A = sorted(r[2] for x, r in s)
    a = A[0] if A else np.nan
    sl = (a - prev) / 0.25 if prev is not None else np.nan
    print(f"t={t:5.2f} nbranch={len(s)} A_min={a:+.6f} all={np.round(A,5)} slope/R={sl / R:.4f}", flush=True)
    prev = a
Vs = Vaidya(M=1.0, dlt=0.01)
r = P.caps(R, 30.0)
print('late time t=30:', [x[1][2] for x in r])
print('v_E =', np.sqrt(3) / 2 ** (4 / 3))

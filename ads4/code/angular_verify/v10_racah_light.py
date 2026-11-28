"""v10: does the claimed Racah dual stay feasible on LIGHT spin-0/2 single traces (b=0 only, twist < g, allowed by
the paper's gap assumption, O(1) amplitudes, endpoint cost W)?  Feasibility there needs W_0 - f_0 >= 0.
usage: python3 v10_racah_light.py g n"""
import sys
from fractions import Fraction as Fr
from common import W
g = int(sys.argv[1]); n = int(sys.argv[2]); s = 0
A = -W(4, Fr(g + 4))
c = [Fr(1)]
for k in range(n): c.append(c[-1] * (-n + k) * (n + s + k) / ((k + 1) * (s + 1 + k)))
worst = (10.0, None)
for i in range(1, 4000):
    D = Fr(3) + Fr(i, 4000) * (g - 1)
    T = (D - 1) * (D - 2)
    r1 = Fr(6) / T; r2 = r1 * 4 / (T - 2)
    f0 = A * (c[0] + c[1] * r1 + c[2] * r2)
    m = float((W(0, D) - f0) / W(0, D))
    if m < worst[0]: worst = (m, float(D))
print(f"g={g} n={n}: light b=0 J=2 states, 3<D<g+2: min relative margin {worst[0]:.4f} at D={worst[1]:.3f}  (J=0: margin>0 trivially)")

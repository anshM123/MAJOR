"""check_F.py: [DV] Lemma F (tpf3) bound vs exact |S^_n| (dual-Hahn recurrence, complex128 with rescaling) at random states."""
import sys, math, random
sys.path.insert(0, '.')
import numpy as np
from common2 import n_of_L
from fastS import S_vec, d_of
import tpf3 as TF
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
worst = 0
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 12):
    a = random.uniform(998, 1100); L = random.uniform(3.5, 5.0); tau = random.uniform(950, 1150); sig = random.choice([0.25, 1.25, 2.5, 3.5, 4.5])
    n = n_of_L(L, a); N = n + a; io = 1/N
    S, lg = S_vec(n, a, 1, [complex(sig, tau)])
    ex = math.log(abs(S[0])) + lg[0] + sig*math.log(d_of(n, a, 1))
    # a box of the size used in the certificate around the state
    da = random.uniform(0, 1.0); dt = random.uniform(0, 1.1)
    r = TF.lemmaF3_log(sig, sig, tau - dt*random.random(), tau + dt, a - da*random.random(), a + da, io*0.85, io*1.15, 200, i_s_list=[200, 190, 180, 170, 155])
    b = r[0] if r else float('inf')
    worst = max(worst, ex - b)
    print(f'a={a:.2f} L={L:.3f} n={n} tau={tau:.2f} sigma={sig}: log|S^| exact {ex:.3f}  bound {b:.3f}  gap {b-ex:.2f}', flush=True)
print('max(exact - bound) =', worst, '(must be < 0)')

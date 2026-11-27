"""P2: check the closed-form optimal (Chebyshev-Markov-Krein / LP-duality) bounds against a direct LP.

Spectrum: positive weights c_i at masses m_i >= m_gap = 1.  Variable x = (m_gap/m)^2 in (0, 1].
Known:  int x^{a_k} drho <= B_k   (a_k = p_k/2 for weight m^{-p_k}).   Target: int x^b drho  (weight m^{-2b}).
Closed forms (README Sec. P2):
  (i)   b <  min a_k              : sup = +infinity
  (ii)  b >= max a_k              : sup = min_k B_k             (atom at m = m_gap)
  (iii) a_1 <= b <= a_2, two data : sup = B_1^th * min(B_1,B_2)^(1-th),  th = (a_2-b)/(a_2-a_1)   (Lyapunov, one atom)
Output: out_p2_cmk.txt
"""
import numpy as np
from scipy.optimize import linprog

x = np.unique(np.concatenate([np.logspace(-8, 0, 4000), np.linspace(0, 1, 2001)[1:]]))


def lp_sup(a, B, b):
    A = np.vstack([x ** ak for ak in a])
    res = linprog(-(x ** b), A_ub=A, b_ub=np.array(B, float), bounds=(0, None), method="highs")
    if res.status == 3:
        return np.inf
    return -res.fun


def closed(a, B, b):
    a, B = list(a), list(B)
    if b < min(a):
        return np.inf
    if b >= max(a):
        return min(B)
    (a1, B1), (a2, B2) = sorted(zip(a, B))
    th = (a2 - b) / (a2 - a1)
    return B1 ** th * min(B1, B2) ** (1 - th)


if __name__ == "__main__":
    cases = [  # (weights known: m^-4, m^-8 ; target m^-6), etc.
        ((2, 4), (1.0, 0.2), 3), ((2, 4), (1.0, 0.9), 3), ((2, 4), (1.0, 2.0), 3), ((2, 4), (1.0, 0.01), 2.5),
        ((2, 3), (1.0, 0.3), 4), ((2, 3), (1.0, 0.3), 1),   # target m^-2 from m^-4, m^-6: unbounded
    ]
    out = []
    for a, B, b in cases:
        lp = lp_sup(a, B, b)
        cf = closed(a, B, b)
        note = "(unbounded: mass escaping to m -> infinity)" if np.isinf(cf) else ""
        out.append(f"known x^{a} <= {B}, target x^{b}:  LP = {lp:.6g}   closed form = {cf:.6g} {note}")
    txt = "\n".join(out)
    print(txt)
    open(__file__.replace("p2_cmk_check.py", "out_p2_cmk.txt"), "w").write(txt + "\n")

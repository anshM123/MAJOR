"""Closed forms of the kernel factors  M_pres (coefficient of P_{-1/2+i nu}) and M_flip (coefficient of F^nu_{2l}) for the
parity-even block tensor  sum_{h1,h3} (z_1.eps_{h1})^l (z_3.eps_{h3})^l,  for selected modes.   d_b := nu^2 + (a-b)^2."""
import sympy as sp
from laplace import M_table, nu, a, g1
out = []
def pr(*x):
    s = " ".join(str(v) for v in x); print(s, flush=True); out.append(s)
for l in (1, 2):
    T = M_table(l)
    pr(f"l = {l}:")
    for (q1, q3) in ([(1, -1), (1, 1), (0, 0), (1, 0)] if l == 1 else [(2, -2), (2, 2), (0, 0), (1, -1), (1, 1), (2, 0)]):
        pres = flip = 0
        for h1 in (1, -1):
            for h3 in (1, -1):
                R, S, mu = T[(q1, q3, h1, h3)]
                if mu == 0: pres += R + g1*S
                else: flip += R + g1*S
        pres = sp.factor(sp.cancel(sp.together(pres))); flip = sp.factor(sp.cancel(sp.together(flip)))
        pr(f"  mode ({q1:+d},{q3:+d}):  M_pres = {pres}")
        pr(f"               M_flip = {flip}")
        # large nu
        lp = sp.limit(pres, nu, sp.oo); lf = sp.limit(flip, nu, sp.oo)
        pr(f"               nu -> oo:  M_pres -> {lp},  M_flip -> {lf}")
open('out_mforms.txt', 'w').write("\n".join(out) + "\n")

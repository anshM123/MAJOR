"""
Independent internal check of the eta-derivative tails (l >= 1) of the LC representation in the deep corner, where no radial
reference exists and the Casimir check of the note uses l = 0 only:  S[k, l] must equal d^l/d eta^l S[k, 0] at fixed r, and
S[k+1, 0] must equal r d/dr S[k, 0] + ... (here only the eta test and the r test S[1,0] = r d/dr (r^Delta Gt)/r^Delta).
The l = 0 and l >= 1 tails are fitted independently, so agreement tests the l >= 1 Darboux ansatz (exponent -2 sigma - l).
Finite differences: 6th-order central stencil in chi with step h = 1e-3 (1 - rho_b).
Usage: python3 fd_check.py Delta J c1 c2
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'lowk')); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np
from lcblock import LCBlock
D, J, c1, c2 = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
lb = LCBlock(D, J, c1, c2, N=10000, I=1500, Q=8000, M=8)
st = np.array([-3, -2, -1, 1, 2, 3]); cf = np.array([-1, 9, -45, 45, -9, 1])/60.0
print(f"# Delta={D} J={J} c=({c1},{c2}); max fit residual l=0..4 (i<=300): " + ", ".join(f"{lb.fitres[:301, l, :].max():.1e}" for l in range(5)))
for r in (0.3, 0.6):
    for om in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
        chi0 = np.log((1 - om)/r)
        h = 1e-3*om
        chis = chi0 + st*h
        S0 = lb.series(np.full(1, r), np.array([chi0]))[0]
        Sp = lb.series(np.full(6, r), chis)
        # d/d eta = (1/sinh chi) d/d chi ;  second derivative: d/deta (1/sinh d/dchi)
        d1 = (cf @ Sp[:, 0, 0])/h/np.sinh(chi0)
        d1k1 = (cf @ Sp[:, 1, 0])/h/np.sinh(chi0)
        rel = abs(d1/S0[0, 1] - 1); rel2 = abs(d1k1/S0[1, 1] - 1)
        # l = 2 from the l = 1 output
        d2 = (cf @ Sp[:, 0, 1])/h/np.sinh(chi0); rel3 = abs(d2/S0[0, 2] - 1)
        d3 = (cf @ Sp[:, 0, 2])/h/np.sinh(chi0); rel4 = abs(d3/S0[0, 3] - 1)
        d4 = (cf @ Sp[:, 0, 3])/h/np.sinh(chi0); rel5 = abs(d4/S0[0, 4] - 1)
        print(f"  r={r} 1-rho_b={om:.0e}: |dS00/deta / S01 - 1| = {rel:.1e};  k=1: {rel2:.1e};  l=2 from l=1: {rel3:.1e};  l=3 from l=2: {rel4:.1e};  l=4 from l=3: {rel5:.1e}", flush=True)

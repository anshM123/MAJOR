"""
sm_nancheck (SMOOTH_BOUND.md §2.1): why exactmom.BlockMellin (Chebyshev Mack form, J > 10) returns NaN at (Delta, J) = (1040, 40).
Evaluates mack.Qpoly_cheb on a few n-nodes of BlockMellin's grid (direct n, EM-tail Gauss nodes) and reports the size of the double
Chebyshev coefficients c_i R^i and of log K_n (the weights), i.e. where double precision overflows.  The FLINT Mack polynomial of
h13_exact replaces mack.mack_P_fast (identical rationals, faster).  Usage: python3 sm_nancheck.py [Delta J]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); EX = os.path.join(HERE, '..', 'exactskel'); sys.path.insert(0, EX)
import numpy as np
import h13_exact as H          # patches mack.mack_P_fast -> FLINT
import mack

D = float(sys.argv[1]) if len(sys.argv) > 1 else 1040.0
J = int(sys.argv[2]) if len(sys.argv) > 2 else 40
Df = 3.0; N0 = 400; nq = 240
L = max(4.0*N0, 0.5*D**2/max(2*Df - 1.5, 0.5))
x, w = np.polynomial.legendre.leggauss(nq); uu = 0.5*(x + 1)
nint = N0 - 0.5 + L*uu/(1 - uu)
ns = np.array([0.0, 50.0, 399.0, nint[0], nint[nq//2], nint[-20], nint[-1]])
import warnings; warnings.filterwarnings(os.environ.get('SM_WARN', 'ignore'))
try:
    R, C = mack.Qpoly_cheb(D, J, ns)
    for n, r, c in zip(ns, R, C):
        print(f"n={n:12.4g} R={r:10.4g}  max|c_i| {np.max(np.abs(c)):.3e}  finite {np.all(np.isfinite(c))}")
except Exception as e:
    print('Qpoly_cheb raised:', type(e).__name__, e)

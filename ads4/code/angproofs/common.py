"""Shared setup for ANGULAR_PROOFS code.  Imports code/kernels read-only; caches redirected to angproofs/cache."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['KERNELS_CACHE'] = os.path.join(HERE, 'cache')
sys.path.insert(0, os.path.join(HERE, '..', 'kernels'))
import sympy as sp
S, J = sp.symbols('s J')
CACHE = os.path.join(HERE, 'cache')

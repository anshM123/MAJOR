# Referee: recompute block tensors from scratch (no cache) for a few cases and compare with cached data_modes and predictions
import sys, os, time, numpy as np, mpmath as mp
sys.path.insert(0, sys.argv[1])
os.chdir(sys.argv[1])
from t4_spinning_e2e import get_moments, moments_of, fname, default_N
from moments import TensorMoments
from predict import prediction
from cases import info
l = int(sys.argv[2]); Delta = float(sys.argv[3]); J = int(sys.argv[4])
t0 = time.time()
N = default_N(l, Delta)
fresh = TensorMoments(l, Delta, J, 1.75, N)
print("fresh compute time", time.time()-t0, "N", N, flush=True)
cached = TensorMoments.load(fname(l, Delta, J))
for k in (4, 6):
    for nu in (0.0, 3.0):
        for (q1, q3) in [(l, -l), (l, l), (0, 0)]:
            v1, s1 = fresh.moment(k, nu, q1, q3); v2, s2 = cached.moment(k, nu, q1, q3)
            ex = (-1)**J*mp.mpc(complex(v1))*mp.exp(mp.mpf(s1))
            pd = prediction(l, Delta, J, 1.75, k, nu, q1, q3)
            print(f"k={k} nu={nu} ({q1},{q3}) fresh/cached-1 = {abs(v1*np.exp(s1-s2)/v2-1):.2e}   fresh/pred = {complex(ex/pd).real:.6f}", flush=True)

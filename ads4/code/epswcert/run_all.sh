#!/bin/bash
# epswcert/run_all.sh (EPSW_CERT.md §7): full reproduction, <= 2 processes, one thread each, ~ 6 CPU-min.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 -B ew_repair.py > cache/ew_repair.log 2>&1                         # repaired weights -> cache/k2_continuum_D_<Dg>_{rep,repk}.json
for W in pub rep repk; do
  python3 -B ew_cert.py 300 re $W > cache/ew_cert_300_re_$W.log 2>&1 &      # eps_W, B, r(0), r'(0+) in ball arithmetic
  python3 -B ew_cert.py 1000 re $W > cache/ew_cert_1000_re_$W.log 2>&1
  wait
done
T0_DGMIN=1000 T0_C=1999/10000 python3 -B t3_flat.py 293/10000000 > cache/t3_1000_c1999_shift.log 2>&1   # byte copy of t0cert/t3_flat.py, larger shift
python3 -B ew_t0.py > cache/ew_t0.log 2>&1                                  # realised T0 corollary (pub: beta caps + certified failure; rep: all beta)
python3 -B ew_light.py > cache/ew_light.log 2>&1                            # light side (K, E, tau_0) pub vs rep, via kcert read-only
grep -h RESULT cache/ew_cert_*.log; grep -hE "^---|PASS|FAIL|beta >=|FAILURE" cache/ew_t0.log

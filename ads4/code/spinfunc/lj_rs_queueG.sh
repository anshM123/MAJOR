#!/bin/bash
# HLJ_CLOSE §5: analytic normalisation (fw patched basis + fw_frame Fourier cure) and fw-convention Gamma - kappa for the high-J rows.
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BJ_PATCH=1
while pgrep -f lj_rs_queueH.sh > /dev/null; do sleep 30; done
python3 lj_rs_norm.py 100 1100 2 4 1 p > cache/lj_rs_norm_J100.log 2>&1
python3 lj_rs_norm.py 64 364 2 4 1 p > cache/lj_rs_norm_J64.log 2>&1
for K in 1 0; do
  LJ_KR=$K BJ_NUCUT=13.0 FW_EST=dir python3 lj_rs_gamma.py 100:1100:2:5:4 0 16.3 > cache/lj_rs_gamma_J100_K$K.log 2>&1
  LJ_KR=$K BJ_NUCUT=9.9 FW_EST=dir python3 lj_rs_gamma.py 64:364:2:5:4 0 8.3 > cache/lj_rs_gamma_J64_K$K.log 2>&1
done
echo done > cache/lj_rs_queueG.done

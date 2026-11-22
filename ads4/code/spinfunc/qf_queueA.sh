#!/bin/bash
# qf queue A (H_OO^Q): parity-odd structure at even J = 0, 2 (OO channel), references and near-gap states.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
run() { python3 qf_eeq.py "$@" >> cache/qf_eeq_$1_D$2_J$3.log 2>&1; }
run odd 20000 0 2 3
run odd 300 0 4 5
run odd 1000 0 4 5
run odd 20000 2 2 3
run odd 314.935 2 4 5
run odd 1049.93 2 4 5
run odd 314.925 0 4 5
# m = 1.3 Lam (Dg = 300): OO and ++ at J = 0, 2 (m-dependence within the channel)
run odd 389.54977773476435 0 4 5
run even 389.54977773476435 0 4 5
run odd 389.5575086246883 2 4 5
run even 389.5575086246883 2 4 5
# helicity bases T0(Delta) for the even J = 6 states
L1_DEG=6 python3 qf_helbasis.py 6 20000 306.0297506031432 441.8348597942253 >> cache/qf_helbasis.log 2>&1

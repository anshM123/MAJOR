#!/bin/bash
# order-0 diagnostic: state <TTO> coefficients at Delta = 300 with the geometry at Delta = 1e6, 2e6 (a -> 0 by Richardson)
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
DCOEF=300 python3 direct_rho.py 1000000 8,8,16,16,6,16,16 f80 4.1 | tail -1 &
DCOEF=300 python3 direct_rho.py 2000000 8,8,16,16,6,16,16 f80 4.1 | tail -1; wait
echo ORDER0_DONE

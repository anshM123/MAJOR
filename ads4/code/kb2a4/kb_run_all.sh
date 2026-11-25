#!/bin/bash
# ads4/K_B2A4.md: all outputs to cache/ ; one process at a time, one thread.  ~8 CPU-min (kb_celestial full ~4 min).
ulimit -v 3000000; export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
python3 -B kb_onepoint.py        > cache/kb_onepoint.log      # one-point traces (k = 2 first moments)
python3 -B kb_k2ratio.py         > cache/kb_k2ratio.log       # k = 2: Einstein check, Theta(nu), a4, top-coefficient check
python3 -B kb_celestial.py full  > cache/kb_celestial.log     # k = 3: S(4)/S(2) = (nu^4 - 35/2 nu^2 + 569/16)/400, b2
python3 -B kb_h17.py             > cache/kb_h17.log           # certified H17, K, tau0 (imports ../kcert read-only, -B: no bytecode)

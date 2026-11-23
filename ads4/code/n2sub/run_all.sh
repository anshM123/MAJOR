#!/bin/bash
# N2_SUBREGGE.md: reproduce everything (~3 CPU-min; one process at a time). Reads only published data and read-only project caches.
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 t_conv.py   > cache/t_conv.txt    # kappa image vs pw_gamma.pkl; kernels.json R_-+ vs the chain's u_EE; real zeros of R_- + R_+
python3 t_deshift.py 1000                 # de-shift identity with the published X (1e-13); raw-read-as-stripped mismatch -> cache/t_deshift_1000.txt
python3 n2_tb.py                          # T_B with psi2 / |X| / |Y|; G bump; D_h class R1/R2 -> cache/n2_tb.{txt,json}
python3 n2_table.py                       # corrected tolerance table, both floors -> cache/n2_table.txt
python3 n2_refit.py                       # (1100,100), (1064,64), (364,64) window data refit in the raw convention -> cache/n2_refit.txt
python3 n2_refit_b.py                     # + the CHMRS ring-derivative class nu^2 D -> cache/n2_refit_b.txt
python3 s_chmrs.py                        # CHMRS (3.20) vs the exact scalar skeleton -> cache/s_chmrs.txt
python3 s_cost.py                         # action of the scalar first subleading term at the H_LJ states -> cache/s_cost.txt

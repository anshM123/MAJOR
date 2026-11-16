#!/bin/sh
# Regenerates the outputs of ads4/code/fourpt/tttt/ (SPINNING_TTTT.md; SPINNING_4PT.md sects. 4.5-4.6).  Caches in cache/ (git-ignored).
# Needs ../cache/rel_l2.pkl (fourpt/build_l2.py) and the <JJJJ> mode caches of fourpt (t4_compute_l1.py) for t4_sign_l1.
set -e
cd "$(dirname "$0")"
python3 t1_mft_l1.py              # exact Fock MFT data of [JJ] in the helicity basis   (~2 min)
python3 t2_mft_heavy_l1.py        # absolute two-pair heavy constant
python3 t3_lightcone_l1.py        # block-level leading-twist test (~40 min)
python3 t4_sign_l1.py             # physical sign of the opposite-sign term
python3 t6_kappa_l1.py            # kappa_h A_h for currents
python3 t5_helicity_l2.py         # <TTO> point A (~40 min)
python3 t5_helicity_l2.py 53/4 11/2   # <TTO> point B (~40 min)
python3 t7_fock_l2.py l1          # control: CHL x-space helicity structures, Fock data for currents (~1 min)
python3 t7_fock_l2.py l2 6        # [TT]_{N,J} exact Fock data in the helicity basis, N <= 6 (~5 min)
python3 t7b_closed_forms.py       # closed forms p[T^{E,h}] = A_h P_FK(1,1)/K_l
rm -f out_t8_kappa_exact.txt
python3 t8_kappa_exact.py 1 97/10 13/3 53/4 11/2                              # exact kappa_h, l = 1
python3 t8_kappa_exact.py 2 97/10 13/3 53/4 11/2 61/7 29/5 127/10 47/7        # exact kappa_h, l = 2 (~20 min)
python3 t10_kappa_symbolic.py 1   # symbolic proof of kappa_h A_h = 2^{3J+3/2}, l = 1
python3 t10_kappa_symbolic.py 2   # symbolic, l = 2 (long)
python3 t9_taylor_validate.py     # Taylor-mode operator application vs symbolic word pairs
gcc -O3 -march=native -fopenmp -shared -fPIC cblock12.c -o libcblock12_d.so -lm
python3 t11_feas_taylor.py 1.0 40 2   # <TTTT> feasibility at a heavy point (~10 min)
# l = 2 heavy vertices, absolute constant, end-to-end <TTTT> (third request; long runs: use nohup, logs in cache/)
python3 t12_vertex_l2.py b0.4 1 2 4 8     # heavy vertices of CHL's X_h along beta = 0.4 (~1.5 h)
python3 t14_vertex_fit_l2.py              # E = m^8 T^hh/lambda_h -> 8; asymmetry of the mixed waves
python3 t15_mft_constant_l2.py            # absolute two-pair MFT constant for l = 2
python3 tmoments_l2.py 1.0 40 float64     # Taylor-mode <TTTT> waves on the moment grid, beta = 1, m ~ 40 (~4 h on a loaded machine)
python3 tmoments_l2.py 1.0 70 float64     # second mass (~4 h)
python3 t13_analysis_l2.py cache/modes_l2_D45.861045_J20_d.npz 8 6 0
python3 t13_analysis_l2.py cache/modes_l2_D80.484300_J36_d.npz 8 6 0

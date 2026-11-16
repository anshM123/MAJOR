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

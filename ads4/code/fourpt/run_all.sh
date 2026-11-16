#!/bin/sh
# Regenerates every output of ads4/code/fourpt/ (SPINNING_4PT.md).  Caches (git-ignored) go to cache/.
# Run time (4 cores, shared): symbolic set-up ~15 min; <JJJJ> block tensors ~30 min per (Delta, J) case (10 cases); analysis ~20 min;
# <TTO> words ~5 min, t6 ~10 min, t7 ~1 h.   'sh run_all.sh fast' skips the block-tensor computation if the cache exists.
set -e
cd "$(dirname "$0")"
python3 t1_structures_l1.py            # <JJO> words, conserved basis (v1, v2), exact shadow map; CHL frame (slow part may be skipped)
python3 t1b_helicity_l1.py             # exact CHL helicity structures in the (v1, v2) basis
python3 gen_l1.py                      # symbolic word pairs of <JJJJ> and generated numerical code (cache/gen4_l1.py)
python3 t2_validate_l1.py              # conservation, pair exchange, scalar (H,H) limit
python3 t3_vertex_l1.py                # heavy vertex closed forms, pointwise; kernel conventions
if [ "$1" != "fast" ]; then
  python3 t4_compute_l1.py 0.4 20 40 70 100 10 &
  python3 t4_compute_l1.py 1.0 20 40 70 100 10 &
  wait
fi
python3 t5_analysis_l1.py              # end-to-end comparison, +- structure, 1/m fits, mutation tests
python3 t8_mft_l1.py                   # MFT normalisation of the helicity components
python3 build_l2.py                    # <TTO> words
python3 t6_structures_l2.py            # <TTO>: 3d relation, conserved space, shadow-map consistency
python3 t7_vertex_l2.py                # <TTTT> heavy shadow integrand with Weyl-covariant families

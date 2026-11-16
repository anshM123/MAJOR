#!/bin/sh
# Regenerates every number quoted in ads4/SPINNING_E2E.md.   Usage:  sh run_all.sh [fast]
#   fast: skips the convergence study t4b_convergence.py.  The exact block tensors on the integration grids are cached in data_modes/
#         (git-ignored) and reused; delete data_modes/ to recompute them (about 2.5 h on 4 cores, dominated by l = 2 at m >= 70).
# Outputs: out_*.txt in this directory.
set -e
cd "$(dirname "$0")"
gcc -O3 -march=native -fopenmp -shared -fPIC cblock.c -o libcblock_d.so -lm
gcc -O2 -march=native -fopenmp -shared -fPIC -DQUAD cblock.c -o libcblock_q.so -lquadmath -lm
python3 casimir_terms.py > out_casimir_terms.txt          # radial Casimir recursion for general external dimensions -> terms.json
python3 structures.py > out_structures.txt                # conserved structures (cached in cache_struct_l*.pkl)
python3 laplace.py > out_laplace_table.txt                # tensor Fourier-Laplace transform, table of M (cached in cache_Mtable_l*.pkl)
python3 mforms.py                                         # closed forms of M_pres, M_flip          -> out_mforms.txt
python3 t1_scalar_blocks.py                               # scalar blocks, general external dims     -> out_t1_scalar_blocks.txt
python3 t1b_mft_scalar.py                                 # mean-field decomposition (scalar seeds)  -> out_t1b_mft_scalar.txt
python3 t2_spinning_blocks.py                             # conservation, OPE limit, double vs quad  -> out_t2_spinning_blocks.txt   (~5 min)
python3 t2b_shadow_consistency.py                         # shadow coefficients vs conservation      -> out_t2b_shadow_consistency.txt
python3 t5_shadow_direct.py                               # direct test of the shadow representation -> out_t5_shadow_direct.txt     (~10 min)
python3 vertex.py                                         # heavy limit of the shadow integrand      -> out_vertex.txt
python3 t6_laplace_check.py                               # M by direct H^2 integration              -> out_t6_laplace_check.txt     (~30 min)
python3 t3_scalar_mixed.py                                # scalar skeleton end-to-end               -> out_t3_scalar_mixed.txt
python3 t4_spinning_e2e.py 1                              # currents:       exact vs predicted       -> out_t4_spinning_e2e_l1.txt
python3 t4_spinning_e2e.py 2                              # stress tensors: exact vs predicted       -> out_t4_spinning_e2e_l2.txt
python3 t4_spinning_e2e.py 1 small                        # small impact parameter cases             -> out_t4_spinning_e2e_l1_small.txt
python3 t4_spinning_e2e.py 2 small                        #                                          -> out_t4_spinning_e2e_l2_small.txt
if [ "$1" != "fast" ]; then python3 t4b_convergence.py; fi   # truncation / grid error controls      -> out_t4b_convergence.txt       (~30 min)
python3 t7_corrections.py 1                               # components, O(1/m), O(1/m^2)             -> out_t7_corrections_l1.txt
python3 t7_corrections.py 2                               #                                          -> out_t7_corrections_l2.txt
python3 t7b_scalar_skeleton.py                            # scalar reference for the corrections     -> out_t7b_scalar_skeleton.txt
python3 t7c_spin_dependence.py                            # spin-dependent part of the 1/m^2 terms   -> out_t7c_spin_dependence.txt
python3 t8_uniformity.py 1                                # nu = x m, small beta                     -> out_t8_uniformity_l1.txt
python3 t8_uniformity.py 2                                #                                          -> out_t8_uniformity_l2.txt
python3 t9_tables.py                                      # tables of SPINNING_E2E.md                -> out_t9_tables.txt
python3 t9b_summary.py                                    # compact ratio tables (sect. 5)           -> out_t9b_summary.txt
python3 clsd_kernel.py                                    # d = 3 Regge reduction of the CLS-D kernel -> out_clsd_kernel.txt

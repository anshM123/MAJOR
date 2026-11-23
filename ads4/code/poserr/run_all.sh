#!/bin/bash
# poserr/run_all.sh (POS_ERRORS.md §6): reproduce every number of POS_ERRORS.md.  <= 2 processes, 3 GB, one thread each (~12 min).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p cache
python3 curv_lemma.py                                              # Lemma C check (mpmath, ~1 min)      -> cache/curv_lemma.txt
python3 eng_orders.py                                              # engine orders, J = 0 and spin J     -> cache/eng_orders{,_bj}.txt
python3 large_m_map.py 1000 1.3,1.5,2,3 400 &                      # D, Delta_g = 1000, part a (~2 min)  -> cache/large_m_1000.json
LM_TAG=b python3 large_m_map.py 1000 5,10,20 400; wait             #                    part b           -> cache/large_m_1000b.json
python3 large_m_map.py 300 1.3,1.5,2,3,5,10,20 400 &               # D, Delta_g = 300 (~3 min)           -> cache/large_m_300.json
python3 large_m_map.py 1000 gap 40; wait                           # near-gap strip, 1000 (~1 min)       -> cache/gap_m_1000.json
python3 q_large_m.py 1000; python3 q_large_m.py 300                # Q flat parts and weights (seconds)  -> cache/q_large_m_*.json
python3 analyze_D.py 1000; python3 analyze_D.py 300                #                                     -> cache/analyze_D_*.txt
python3 analyze_Q.py 1000; python3 analyze_Q.py 300                #                                     -> cache/analyze_Q_*.txt
python3 analyze_gap.py                                             #                                     -> cache/analyze_gap.txt

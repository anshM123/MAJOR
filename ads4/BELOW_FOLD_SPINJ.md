# Below the fold at large spin: the closed-form spin-J EE engine, the `ν < 8.3` window at `(336, 28/32)`, and the high-J window

*21–22 Nov 2026, below-fold-spinJ agent. Unrefereed, not committed. No existing document or code file was edited; nothing of the
`oo_*` agent was touched or run. New code `ads4/code/spinfunc/bj_*.py` (+ `bj_queue*.sh`), caches `spinfunc/cache/bj_*`, `z2bj_*`,
`hc_chain_z4_bj*` (git-ignored; the last are outputs of the unchanged `hc_chain.py`). ≤ 2 processes, `OMP_NUM_THREADS=1`,
`ulimit -v 3 GB` on every job (one early unconstrained ring test was OOM-killed at ~13 GB; all later jobs ran under the limit).
Read-only imports: `bg_ser`, `bg_geom`, `bg_ee`, `bh_spinJ`, `bg_scalarJ`, `l5_spin`, `fw_frame`, `z3_nodes`, `hc_chain`/`z4_chain`/
`w2_chain` (run unchanged, in-process wrappers only), `l9_analysis`, `fw_gamma`.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model;
**[I]** open. `c₂ = m² ×` relative correction. "Engine" = the analytic-radial 1/Δ expansion on the real section D (BELOW_FOLD2); its
EE beyond-LO equivalence with the crystal cycle is assumed there — §2.4 below is the first test of it at large J.

@@SUMMARY@@

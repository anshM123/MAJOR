# `⟨TTTT⟩` beyond one point: H1 (heavy formula at more `k, ν, β, m`) and H2 (reading B for gravitons)

Status: written 18 Nov 2026 by a derivation/numerics agent; **not refereed; not committed**. Code: `ads4/code/fourpt/tttt/`
(`tmoments_l2.py` production, `t16_kernels_ext.py`, `t17_readingB.py`; logs and caches in `tttt/cache/`, git-ignored). Builds on
`SPINNING_TTTT.md` (Taylor-mode waves, `l = 2` vertices `T^{hLhR} = 8√(λ_{hL}λ_{hR})/m⁸`, Fock data), `SPINNING_ODDK.md` (robust extraction)
and `SKELETON_TEST.md` §8 (reading A/B for currents). Targets hypotheses H1 and H2 of `RESULT_D.md`.

Marks: **[D]** derived; **[DVE]** derived and verified end to end on exact blocks; **[NE]** numerical evidence; **[ID]** closed form identified
from exact data; **[I]** inferred; **[ND]** not done.

## 0. Summary

SUMMARY

## 1. Data and method

**Production.** `⟨TTTT⟩` waves `W[X_{hL},X_{hR}]` (CHL's exact helicity structures, all 196 word pairs by Taylor mode, `SPINNING_TTTT.md` §4)
on the moment grid of `e2e/moments.py`, float64, `N = 8000`:

| `β` | `m` | `(Δ, J)` | grid points | wall time (1 core, shared machine) |
|---|---|---|---|---|
| 1 | 39.3 | (45.861, 20) | 3336 | (earlier) |
| 1 | 70.0 | (80.484, 36) | 3475 | (earlier) |
| 1 | 100.7 | (115.108, 52) | 3272 | 8820 s (new) |
| 0.4 | 42.2 | (44.565, 8) | 3463 | B04M40TIME (new) |
| 0.4 | 72.0 | (74.964, 14) | 3397 | 8970 s (new) |

The production code now builds the 25 seed blocks one at a time and precomputes their radial sums for all points. Memory is
0.25–0.9 GB per run, against 3.8 GB before, which allows two runs side by side. The cached modes are reproduced bit for bit.

`β ≈ 0.4` is the coordinator's first choice, with masses chosen so that `A_opp` is not tiny. At fixed `β`, `J ∝ m`, so `J = 8, 14` are
forced. There `A_opp = 0.14, 0.33`, and the Weyl-odd `O(1/m)` term of the integrand is large (`λ_o/λ_s = 8.9, 2.9`; `SPINNING_TTTT.md` §3.6).

**Extraction (robust; `REFEREE_SPINNING_ODDK.md` M1).** No class fit is used. A fixed kernel contraction of the measured modes is compared
with the same contraction of the predicted modes. On the polarisation circles (`SPINNING_4PT.md` §4.5) `D_1 = cos(θ_out − θ_in)`, `D_2 = −1` and
`D_5 = sin(θ_out − θ_in)`. Per pair:

| kernel | definition | modes `(q_in, q_out)`: weight |
|---|---|---|
| E | `D_1²` (the project's even "`D₁⊗D₁`" for `l = 2`) | `(0,0)`: ½; `(±2,∓2)`: ¼ |
| FF | `𝔖_ff = D_1² + D_1D_2 = D_1(z_in·z_out)` (flip-free, `SUPERCONVERGENT_UPLIFT.md`) | E plus `(±1,∓1)`: −½ |
| O | `𝔄 = D_5D_1` | `(2,−2)`: `1/(4i)`; `(−2,2)`: `−1/(4i)` |

These are used as `E⊗E`, `FF⊗FF` and `O⊗O`, together with the leading single mode of each wave. The prediction is that of
`SPINNING_TTTT.md` §4.5, with `a = k + 5` and nothing fitted.

The tiny-`t` (dropped) region is estimated only for the leading mode. The `oddk` estimator fails for modes that are small in the bulk: it
exceeds the `(0,0,0,0)` mode of `W[o,o]` thirtyfold. Kernel ratios are therefore raw.

## 2. H1: results

RESULTS_H1

## 3. H2: reading A / B / B̃ for gravitons at fixed `b̂ = 1.04` (`t17_readingB.py`, `out_t17_readingB.txt`)

RESULTS_H2

## 4. What is established, and what is not

MARKS

## 5. Code

| file | content | output |
|---|---|---|
| `tmoments_l2.py` | Taylor-mode `⟨TTTT⟩` waves (blocks one at a time, precomputed radial sums) | `cache/modes_l2_*.npz`, `cache/prod_*.log` |
| `t16_kernels_ext.py` | kernel moments `E⊗E`, `FF⊗FF`, `O⊗O` and leading modes vs prediction, `k = 2…7`, `ν = 0, 1, 3`; convergence in `m` | `out_t16_kernels_ext.txt` |
| `t17_readingB.py` | reading A / B / B̃ at fixed `b̂`, single and double (opp/same) ratios | `out_t17_readingB.txt` |
| `t7_fock_l2.py l2 8 7`, `t7b` (ext) | Fock data for `[TT]` at `N = 7, 8` and the closed-form check | `out_t7_fock_l2ext.txt`, `out_t7b_closed_forms_ext.txt` |
| `t12_vertex_l2.py b0.9 1 2`, `t14_vertex_fit_l2.py` | `E` on a second ray | `out_t14_vertex_fit_l2.txt` |

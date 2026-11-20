# `Ω^EE` at `Δ_g = 1000`: the large-`J` near-gap region (even `J ≥ 32`, odd `J ≥ 17`)

*20 Nov 2026, EE-1000-large-J agent. Unrefereed, not committed. No existing document or code file was edited. New code
(`ads4/code/spinfunc/`): `ee1000_common.py`, `ee1000_rows.py`, `ee1000_direct.py`, `ee1000_chain.py`, `ee1000_w2rerun.py`,
`ee1000_checks.py`, `ee1000_directcmp.py`, `ee1000_shapecheck.py`, `ee1000_band.py`, `ee1000_smallnu.py`, `ee1000_tables.py`,
`ee1000_ensemble.sh`, `ee1000_queue*.sh`, `ee1000_run.sh`. Caches: `spinfunc/cache/ee1000_*`. Project modules are imported read-only,
and patches are applied inside the process only.*

Grades: **[P]** proved, **[D]** derived, **[NE]** numerical evidence, **[M]** model, **[I]** open — as in `EE_NEARGAP.md`.

**STATUS: IN PROGRESS.** This file is written early and updated as rows finish; the queue state is in §8.

## 0. Summary (provisional; numbers final only where marked)

1. **The "K₂ `+−` breakdown" (`EE_NEARGAP.md` §5, and S1_S4's `J = 32` on 1600) is not a breakdown of the Laplace series. It is a
   numerical bug in the local de-shift** [NE, checked against mpmath at 50 digits; `cache/ee1000_checks_jacobi.txt`].
   * `q3._shift_generic` evaluates the `+−` LO shape `((1+c)/2)⁴ P^{(0,8)}_{J−4}(c)` at complex `c = 1 − 2(x ± i/m)²`. It does this
     through `scipy.special.eval_jacobi`.
   * The complex branch of that function loses all accuracy for `n = J − 4 ≳ 24` at `x ≳ 0.5`. Relative error: `5·10⁻²` at
     `n = 24, x = 0.8`; `1.1` at `n = 28, x = 0.8`; `2·10⁶` at `n = 36, x = 0.8`; `7·10²⁰` at `n = 60`.
   * The raw Laplace orders are fine there.
   * Fix: replace it by the three-term recurrence (error `≤ 7·10⁻¹⁵`; `ee1000_common.Dc_fixed`, applied in-process).
   * Effect on `(1040, 40)`: `Γ_{+−}` was `7·10³ … 4·10¹⁰` on `x = 0.5–0.8`; it is now `|Γ| ≲ 2·10³`. That is the size of `Γ_{++}` and of
     `C3H`, which have no de-shift and were never affected (`cache/ee1000_checks_deshift.txt`).
2. **Re-run of the EE_NEARGAP main chain with the fix** (`ee1000_w2rerun.py`, `cache/ee1000_w2rerun_main_fix.txt`; same rows and states):
   * `Δ_g = 1000`, even `J ≤ 24` and odd `J ≤ 15`: unchanged.
   * `J = 28`: max `|dF_hh|/N_h` falls from 0.46 to 0.047. The 0.46 was the bug.
   * `J = 32`, rows `1040 + 1600`: `+0.9535`. It was `−3.63`, which is why EE_NEARGAP had dropped the 1600 row.
   * The 1040-only `J = 32` value of EE_NEARGAP (`+0.967`) contained a spurious `+33 N_s` on the `+−` diagonal at `1.05Λ`.
   * Overall: **`n = 1`, min `+0.2021`** [NE + M].
   * `Δ_g = 300` (reported for the 300 agent, not pursued here): `J = 28` at `1.1Λ` in convention A goes from `−0.095` to **`−0.310`**.
     The 336 row at `x ≥ 0.7` was affected; `J ≤ 24` is unchanged.
3. Large-`J` rows, direct validation, the representation of an oscillating `Γ`, and the small-`ν` window: §3–§6 (in progress).

## 1. The Jacobi bug (details)

`cache/ee1000_checks_jacobi.txt`. The table gives the relative error of `P^{(0,8)}_n(c)` at `c = 1 − 2(x + i/1040)²`, against mpmath at
50 digits.

| `n` | scipy, `x = 0.5` | scipy, `x = 0.7` | scipy, `x = 0.8` | recurrence (all `x ≤ 0.9`) |
|---|---|---|---|---|
| 20 | `9·10⁻⁹` | `1.5·10⁻⁵` | `5.6·10⁻⁵` | `≤ 2·10⁻¹⁵` |
| 24 | `9·10⁻⁷` | `1.3·10⁻³` | `5·10⁻²` | `≤ 2·10⁻¹⁵` |
| 28 | `6·10⁻⁴` | `0.35` | `1.1` | `≤ 7·10⁻¹⁵` |
| 36 | `7·10⁻²` | `5·10⁴` | `2·10⁶` | `≤ 3·10⁻¹⁵` |
| 60 | `2·10⁸` | `5·10¹⁷` | `7·10²⁰` | `≤ 5·10⁻¹⁵` |

* Real arguments (`Dsh`, the LO shapes) are accurate, as is numpy's Legendre `legval` at complex argument (`++`).
* Only the `+−` de-shift term is affected: even `J ≥ 28` at large `x`, odd `J ≥ 29` (parity-odd `+−`).
* The fix touches neither the crystal records nor the `a_h` fits (`a_{+−}` changes by `< 10⁻¹⁰`).
* It is applied to `Q.A.shift` and to `Q.P4.shift`.

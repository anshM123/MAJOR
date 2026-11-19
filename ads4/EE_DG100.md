# `Ω^EE` at `Δ_g = 100`: per-state positivity near the gap (IN PROGRESS — partial; hand-over at the end)

*19 Nov 2026, EE-Δ_g=100 agent. Unrefereed, not committed. No existing document or project code edited. New code:
`ads4/code/spinfunc/y1_rows.py` (K₂ EE rows), `y1_queue.sh`, `y2_chain100.py` (chain), `y4_direct.py` (n1/n3/p8 EE points),
`y5_n4.py` (n4 thimble EE points); caches `spinfunc/cache/y*`. Project modules (and the 300/1000 sweep's `w0_ee_common`, `w1` caches)
are imported/read only.*

Grades: [P] proved, [D] derived, [NE] numerical evidence, [M] model, [I] open.

## 1. Chain (`y2_chain100.py`)

* S1_S4 chain (q3_analysis conventions: `K₂` linear chain, rows linear in `1/m`, aligned off-diagonals (sign from `q²`'s
  `m Im Γ_off`), `c_off = mΓ_off`, `Ĉ₃` without de-shift, phase strip `logK − iπ(Δ+J)`, combined even + odd pass over all
  exact-skeleton states of `Δ_g = 100`).
* **Validation [NE]:** run on the `q²` rows (`K₂` only) it reproduces `q3_analysis.run([100])` to all printed digits.
  Side result: the `q²` functional in this (S1_S4-corrected) chain at `Δ_g = 100`, `K₂` only, gives **`n = 1.35`** (`J = 10`, 1.21Λ,
  −0.329; `J = 8` −0.054), whereas the older l9 conventions gave `n = 1`. Alignment overlaps reproduce S1_S4 §2.2(a).
* `k = 2` term two ways: LOCAL (`(1 ± u_EE)/2` de-shift, singular at 1.17, 6.64) and NONLOCAL `Re∫X(G/a_h − Â_h[f_LO]) + ∫ψ₂f_LO`
  in the χ-normalisation (`XÂ_h[f] = Y(ĥ₋f(ν+i) + ĥ₊f(ν−i))`, `ĥ = R/n`); the interpolated quantity
  `Γ_NL = m²a(ν)(G/a_h − Â_h[f])/env` is finite through the zeros. `K₂` imaginary parts at `ν < 20` replaced by the exact LO value.

(results to follow)

### Interim snapshots (19 Nov, 20:00–21:10; superseded by §2 when final)
* Main protocol, EE, nonlocal min-norm, with n1/n3 (even) and n4 (x = 0.05: even J = 4, 6; odd J = 5, 7, 9) direct: even min
  +0.240 (J = 4, 1.04Λ); **odd J = 9 at 1.12Λ ≈ −0.005** (driven by the Δ = 100 row's n4 x = 0.05 EE value, Γ error ≈ 580 c2, interpolated in
  1/m towards the Δ = 200 row's `K₂` node). With the extra rows (odd n4 at 1.1Λ used) J = 9 is +0.184 (1.09Λ).
* Local de-shift (u_EE) version: odd J = 9 +0.15; the local form is singular at ν = 6.64 and is reported for comparison only.

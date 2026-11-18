# Direct test of H4 (and partly H1): exact ⟨TTTT⟩ moments at two near-gap binding states, Δ_g = 100

Status: 18 Nov 2026, numerics agent. **Unrefereed. Not committed.** Code: `ads4/code/directgrav/` (`run_all.sh`; caches and logs in
`directgrav/cache/`, git-ignored by `directgrav/.gitignore`). It imports `fourpt/tttt/` (Taylor-mode waves, the t13 prediction), `e2e/`,
`exactskel/` and `omegau/` **read-only**. Nothing outside `directgrav/` and this file was written.

Marks: **[D]** derived; **[DVE]** derived and checked end to end on exact blocks; **[NE]** numerical evidence; **[I]** inferred or assumed;
**[ND]** not done.

## 0. Summary

> **Update (§7, same day): the `+−` excess of §3 was a numerical artefact of the cached waves, not physics.** The t19 cleaning of
> `SPINNING_TTTT_EXT.md` §7 removes it: mask cut 40, noise floor, row tails and wedge. After cleaning:
> * `+−` under `𝔄⊗𝔄` agrees with the heavy formula wherever it is resolved (`k = 5` at small `ν`, `k = 7` to `ν ≈ 6`).
> * At `k = 7` the exact `+−` moment matches the reading-B-transferred skeleton to `≤ 1.6·10⁻³`. Reading A is off by `3–4 %`.
> * `+−` under `K*` and all `+−` moments at `k ≤ 4` (`k ≤ 6` under `K*`) remain **noise-limited** in float64 (§7.3). That is undetermined,
>   not contradicted.
> * The `++` conclusions are unchanged and slightly strengthened.
>
> §3 below is kept as the record of the raw extraction; read it together with §7 (at the end).

**Verdict.**
* **H4 is not settled.** The test reaches only `x = ν/m ≲ 0.05`. Those nodes carry **0.2–0.4 %** of the remainder norm of either
  functional. The near-gap action and both D-strong failures live at `x ≳ 0.5`, which is out of reach.
* **`++` (same helicity), on the nodes that could be resolved [NE]:**
  * the exact graviton moments agree with the exact scalar-skeleton moments up to one overall factor, `c* = 1.013` (`m = 102`) and
    `1.008–1.010` (`m = 128`);
  * after removing `c*`, the shape agrees to `≤ 3·10⁻³` for `k = 3…7` and `x ≤ 0.034`;
  * the graviton-minus-skeleton shift of the D-safe action on those nodes is `|ΔF| ≤ 2·10⁻³ N_s`, against margins of `0.27` and `0.32`.
* **`+−` (opposite helicity): the exact moments cannot be resolved at these low-`J` states [NE, negative].** For `k ≤ 5` the measured
  opposite-helicity kernel moments exceed the leading prediction by factors of `10²–10¹⁰`, under every kernel. So on `+−` this test
  says **nothing for or against** H4, H1 or H2.

**Scope actually run.** Two `(m, J)` states at `Δ_g = 100`. Each state gives both helicity waves.

| state | `m/Λ` | `J` | `m` | `Δ` | AdS `β` | role |
|---|---|---|---|---|---|---|
| A | 1.04 | 4 | 102.44 | 104.037 | 0.088 | D-safe `+−` J=4 (`0.241`) and `++` J=4 (`0.266`); D-strong `+−` pocket (`−1.39 N_s`) |
| C | 1.30 | 6 | 128.05 | 129.713 | 0.101 | D-safe `++` J=6 (`0.316`, near-gap binding `++`) |

Not run: `J = 0, 2`; `m/Λ = 1.0, 1.05, 1.1`. Each state costs about 2.4 h of one core (see §1), and the budget allowed two.

**Two corrections made during the run.**
* **`J = 5` (`+−`) is not a parity-even exchange for identical `T`.** By Bose symmetry the `+−` structure has parity `(−1)^J`. So the
  odd-`J` `+−` states of the flat functional are parity-odd exchanges (H12), which this code does not treat. The `J = 5` run was stopped
  and replaced by state C.
* **Corrected even-kernel weights.** The even kernels use the corrected tensor circle weights of `SPINNING_TTTT_EXT.md` §6.2
  (`n_|q| = 1/36, 1/16, 1`), not the naive weights of `t16_kernels_ext.py`.

## 1. Feasibility (step 1)

**What the functionals need near the gap [D, from `OMEGA_NF.md` / `OMEGA_U.md` / `exactskel/e7`, `e3`].**
`Ω^{NF,SC} = ∫₀¹dp Σ_{k=2}^{7} φ_k(p) ε_kΛ^{2k}Ĉ_{k,Λp}` is evaluated on 64 Gauss nodes `p_i`, with `ν_i = Λp_i ∈ (0, 98.5)` and
`x = ν/m` up to `1/(m/Λ)`. The kernels per `k` are:

| `k` | `Ω^{NF,SC}` (D-safe) | `Ω^U` (M02, D-strong) |
|---|---|---|
| 2 | `𝔖_ff⊗𝔖_ff` | `D₁²⊗D₁²` |
| 4, 6 | `𝔖_ff⊗𝔖_ff`, via `Ψ_k/N^{ff}_k` | `D₁²⊗D₁²` |
| 3, 5, 7 | `𝔄⊗𝔄`, via `Ψ_k/N_k` | `𝔄⊗𝔄` |

`Ω^U` also has forward terms at `ν = 0`.
* `++` needs `W[X_s,X_s]` and `+−` needs `W[X_o,X_o]`. `W[X_s,X_o]` is the interference.
* **Where the weights sit.** The NF02 weights vanish like `p(1−p)⁵` but carry `p^{2(k−2)}` factors (`φ_k` in `nf_common`). So the remainder
  norm `Σ|u^{k−2}φ_k w env_k|` sits at large `p`:

  | computed nodes | share of the norm, NF02 | share of the norm, M02 |
  |---|---|---|
  | `ν ≤ 10.7` | 5–7 % | 1–3 % |
  | resolved set actually used (§3) | 0.4 % | 0.3–1.1 % |

  Almost all of the computed share comes from `k = 2, 3`; for `k ≥ 4` it is `≤ 0.5 %`.

**What the code can compute.**
* `fourpt/tttt/tmoments_l2.TaylorWaves` computes the exact waves at any parity-even `(Δ, J)`. `directgrav/gwaves.py` is a resumable driver
  for it, with the same numerics: float64, `N = 8000`, the e2e moment grid and the same convergence mask.
* **Cost:** `2250 s` setup plus `16 × 370 s` of chunks, i.e. **≈ 2.3–2.4 h per state on one core**, at about 0.25 GB.
* The moments are real-section (position-space) moments. They lose about `e^{−πν}` to cancellation.

**Calibration on the scalar skeleton, where the exact answer is known [NE; `c1_calib.py`, `cache/c1_100_1.040_4_Df3.txt`].** The
position-space moment of the scalar block on the same grid, compared with the exact Mellin/Barnes moment (exactskel):
* `R_pos − R_ex` is `+5·10⁻⁴` to `+8·10⁻⁴`, flat in `ν` (an offset common to graviton and scalar, irrelevant for ratios);
* it breaks down at `ν ≈ 12` (`k = 2`) and `ν ≈ 17` (`k ≥ 3`), i.e. `x ≲ 0.12–0.17`.

So **the full ν-integrated functional is out of reach** with this route: it needs `x` up to 1, which takes 100–150 digits, or a Mellin
representation of spinning blocks, which does not exist here.

**The reduced test, defined.**

| | content |
|---|---|
| state | `(m, J)`, with waves `W[s,s]` and `W[o,o]` |
| nodes | the functional's own Gauss nodes with `ν ≤ 15.5` (17 of 64), plus `ν = 0, 1, 3` |
| moments | `k = 2…7` |
| kernels | `𝔖_ff`, `D₁²`, `K* = (D₁+D₂)(D₁−½D₂)` and `𝔄`, each tensored with itself |

1. **`R_grav(k,ν)`** = exact kernel moment / leading heavy prediction (t13 `design`/`Tpred`, `a = k+5`, nothing fitted).
2. **Normalised graviton moment** `ĉ = 2P_{−½+iν}(cosh β)·R_grav`. This is exactly the functional's normalisation: `Ψ_k/N_k` is defined
   so that its leading heavy action is `2P|c|²/m^{2k}`. An overall per-state constant (vertex normalisation, the `¼`) multiplies `F` and
   `N_s` alike and cannot change a sign.
3. **Comparison with the exactskel model on the same node.**
   * `++`: the model is `ĉ_model = m^{2k}C^{skel}_{exact}`, i.e. `R_skel = ex/(2P)` (Df = 3 and 2).
   * `+−`: the model adds the reading-B term, `E_k(d^J_{44} − P_J)`.
4. **Action shift** `ΔF = Σ_k Σ_{i computed} σ^k u^{k−2}φ_k(p_i)w_i ε_k 2P_i (R_grav/c* − R_skel)`, the clean H4 metric. Also reported:
   the shift in exactskel's model variables (`a2`), which mixes in the odd-`k` B′-envelope convention (I3).
5. **`d_crit`** = `F_ex / Σ_all|w env|`: the uniform pointwise remainder, on all 64 nodes, that would just cancel the skeleton margin.

**What this tests, and what it does not.**
* **It tests H4 pointwise at `x ≲ 0.05` on the resolved moments,** and **it bypasses H1 there.** No heavy formula is used beyond the
  normalisation that defines `Ĉ_k`.
* **It does not test:**
  * `x ≳ 0.1`, where the whole near-gap action and both D-strong failures live;
  * `Ψ_k/N_k` versus the project moment `C_k` beyond leading order (the identification exactskel also makes);
  * forward terms;
  * states other than A and C.

## 2. Results: same helicity (`++`) [NE]

**`R_grav` at state A, `J = 4`, `m = 102.4`** (from `cache/a1_A.log`; skeleton Df = 3 shown):

| `k` | `ν` | `K*` | `𝔄` | `𝔖_ff` | `D₁²` | `R_skel` (exact) |
|---|---|---|---|---|---|---|
| 2 | 0.03 / 1.31 / 2.64 | 1.0028 / 1.0015 / 0.9887 | 1.0144 / 1.0202 / 1.0218 | 1677 / 2393 / 5937 | 1032 / 1844 / 6862 | 1.0006 / 1.0009 / 1.0018 |
| 3 | 0.03 / 2.64 / 5.45 | — | 1.0125 / 1.0135 / 1.0136 | — | — | 1.0001 / 1.0016 / 1.0065 |
| 4 | 0.03 / 2.64 / 5.45 | 1.0121 / 1.0134 / 1.0179 | — | 1.110 / 1.94 / 10.9 | 1.128 / 1.48 / 19.5 | 0.9991 / 1.0009 / 1.0069 |
| 5 | 0.03 / 2.64 / 5.45 | — | 1.0108 / 1.0125 / 1.0180 | (1.011 / 1.062 / 1.59) | | 0.9975 / 0.9997 / 1.0068 |
| 6 | 0.03 / 2.64 / 5.45 | 1.0089 / 1.0111 / 1.0183 | — | 1.0054 / 1.0131 / 1.0825 | 1.0068 / 1.0104 / 1.0370 | 0.9953 / 0.9978 / 1.0060 |
| 7 | 0.03 / 2.64 / 5.45 | — | 1.0062 / 1.0087 / 1.0168 | | | 0.9924 / 0.9952 / 1.0045 |

**What the table shows.**
* **`K*` (even `k`) and `𝔄` (odd `k`) are smooth for `ν ≲ 6–7` (`k ≥ 3`) and `ν ≲ 2` (`k = 2`).** They sit about `1 %` above the
  skeleton.
  * The offset is almost the same for every `k` and `ν`, so it is a common factor `c*`: 1.013 (A) and 1.008–1.010 (C).
  * `m²(c*−1) ≈ 135` (A) and `≈ 140–170` (C). This is the size of the known `O(1/m²)` correction to the `l = 2` vertex,
    `m²(E/8−1) ≈ 120–140` (`SPINNING_TTTT.md` §3.6). It is a rescale and harmless.
* **The isotropic-mode kernels `𝔖_ff` and `D₁²` are not resolved at `k = 2, 3, 4`.**
  * At `k = 2` they are off by `10³` even at `ν = 0`; at `k = 4` by `10–100 %`.
  * At `k = 6` they are off by up to `14 %` at `ν ≤ 6.6`; at `k = 5, 7` they are fine.
  * `d1_tprofile.py` locates the `k = 2` `++` contamination of the `(0,0)` modes in the small-`t` region, as in `SPINNING_TTTT_EXT.md` §6.3.
* **So, for `Ω^{NF,SC}` as defined (`𝔖_ff`), only `k = 3, 5, 7` and `k = 6` at small `ν` are resolved. `K*` is used as the even-`k`
  proxy.**

**Pointwise H4 residual `max|R_grav/(c* R_skel) − 1|` on the resolved nodes** (`a3_restricted.py`, `cache/a3_summary.txt`; `ν ≤ 6.6` for
`k ≥ 3`, `ν ≤ 1.92` for `k = 2`; Df = 3, with Df = 2 the same to `±4·10⁻⁴`):

| state | `k = 2` (`K*`) | 3 | 4 | 5 | 6 | 7 | `d_crit` (NF02) |
|---|---|---|---|---|---|---|---|
| A (`J=4`, `1.04Λ`) | 1.2e−2 | 3.1e−2 (1.5e−3 at `ν ≤ 3.5`) | 3.0e−3 | 2.5e−3 | 1.2e−3 | 1.2e−3 | 1.5e−2 |
| C (`J=6`, `1.3Λ`) | 6.1e−2 | 1.1e−2 (1.0e−3 at `ν ≤ 3.5`) | 1.2e−3 | 1.9e−3 | 1.1e−3 | 1.1e−3 | 1.6e−2 |

* At `k = 3` the residual grows sharply between `ν = 5.5` and `6.6`, and the moment breaks down at `ν = 7.8`. The `ν ≤ 3.5` value is the
  reliable one.
* At `k = 2` (`K*`) the residual at C is already `5 %` at `ν = 0`, and `R` drifts with `ν`. This is the unresolved small-`t` region
  (`SPINNING_TTTT_EXT.md` §6: "`k = 2` not resolved under any even kernel"). It is not counted as H4 evidence.

**Action shift (clean H4 metric, `c*` removed, resolved nodes only):**

| state | functional (even kernel) | `ΔF/N_s` | skeleton margin `F/N_s` | coverage of the remainder norm |
|---|---|---|---|---|
| A `++` | NF02 (`K*`/`𝔄`) | `−7.4·10⁻⁴` (`−2.3·10⁻⁴` at `ν ≤ 3.5`) | +0.266 | 0.39 % |
| C `++` | NF02 (`K*`/`𝔄`) | `+1.9·10⁻³` (`+1.6·10⁻³` at `ν ≤ 3.5`) | +0.316 | 0.41 % |

With `𝔖_ff` or `D₁²` (`Ω^U`) as literally defined, the shifts are `O(10–40 N_s)`: these are the garbage `k = 2, 4` moments, not physics.

**Interference.** `|K(W[s,o])|/√|K(ss)K(oo)|` is `10⁻²–10` (`a1` tables). It is not interpretable while `W[o,o]` is unresolved (§3).

## 3. Results: opposite helicity (`+−`) [NE, negative]

**Measured `|K(W[o,o])|/|K(W[s,s])|` against predicted** (predicted: `82` (A), `16` (C), from `λ_o/λ_s`; `cache/a1_*.log`, scratch
`sz.py` output reproduced in the logs):

| `k` | A, `𝔄` | A, `K*` | C, `𝔄` | C, `K*` |
|---|---|---|---|---|
| 2 | `8·10⁸` | `8·10⁸` | `10¹⁰` | `3·10¹⁰` |
| 4 | `6·10⁴` | `9·10⁴` | `9·10⁵` | `2·10⁶` |
| 5 | `2·10²` | `5·10³` | `5·10³` | `10⁴` |
| 6 | **81** | 510 | 37 | 63 |
| 7 | **81** | 110 | **16** | **16** |

* **Only `𝔄⊗𝔄` at `k = 6, 7` (A) and `k = 7` (C) agree with the heavy formula.** There `R = 0.99–1.00` for `ν ≲ 5` at A, and
  `R = 1.0002–1.02` at C.
* **For `k ≤ 5` the opposite-helicity wave is dominated by a contribution the leading formula does not contain.** The excess falls by
  about `10²` per unit of `k`, the signature of a small-radius piece.
* **`𝔄` kernel (`d1_tprofile.py`, A).** Removing the grid columns `t < 2` changes `R(k=5)` from `−2.8` to `+1.44`, and `t < 4` gives `1.09`.
  So for the odd kernel the excess sits in the small-`t` region.
* **`K*` kernel.** The excess (`R ≈ 60` at `k = 5`, `≈ 6` at `k = 6`) does not move under the same cuts. It is a bulk effect.
* **Plausible origin [I].** At `J = 4, 6` the opposite-helicity vertex factor is tiny: `A_opp = (−J)₄/(J+1)₄ = 0.014` at `J = 4`. The Weyl-odd
  `O(1/m)` terms, which `SPINNING_TTTT.md` §3.6 found to be "large at small `J`", are then not suppressed at `m ≈ 100–130`.
* **Roundoff check of the Taylor stage (`d2_ldcheck.py`, float64 against long double at 20 points; `cache/d2_A.txt`): inconclusive.**
  * The long-double values differ from float64 by factors of `10²–10⁴`, for `W[o,o]` and also for the control `W[s,s]`.
  * This holds even at bulk points (`t = 8–15`), where the float64 `W[s,s]` moments match the heavy formula to `10⁻³`.
  * So the long-double path, as driven here (float64 seed sums, long-double word application), is itself not valid. It is not evidence
    either way.
  * **Whether the `+−` excess is physical or numerical is undecided.**
* **Consequence.** At the binding `+−` states (A: D-safe `0.241`, D-strong `−1.39 N_s`) this test gives no evidence on H4, H1 or H2.
  The exact `+−` moments at low `J` are not known even at small `ν`.

## 4. Grades and consequences for `RESULT_D.md`

**H4, `++`:**
* **[NE], small `x` only.** At two binding states (`J = 4` at `1.04Λ`, `J = 6` at `1.3Λ`), for `x ≤ 0.034–0.065`, `k = 3…7` (`k = 6` and
  `k = 4` via `K*`, not `𝔖_ff`): exact graviton moments equal exact scalar-skeleton moments (Df = 2 or 3) times one common factor, to
  `≤ 3·10⁻³`. The common factor is `1 + O(140/m²)`, consistent with the vertex normalisation.
* The resulting D-safe action shift on these nodes is `≤ 2·10⁻³ N_s`, i.e. `≤ 1 %` of the margin. But these nodes carry `≤ 0.4 %` of the
  remainder norm.
* **Not tested** where the action is decided (`x ≳ 0.1`). H4 stays **[I]** there.

**H4, `+−`:** **untested** (moments unresolved at low `J`).

**H1 at these states.**
* The test bypasses H1 only where it resolves the moments. There:
  * `++` under `K*`, `𝔄`, and `𝔖_ff` at `k = 5, 7`: the heavy formula holds, `R − 1 ≈ +0.6…1.3 %`, i.e. `c₂ ≈ 100–170` [NE];
  * `+−`: only under `𝔄` at `k = 6, 7`.
* **New negative finding [NE]:** at low `J` (`β ≈ 0.09–0.10`) the even kernels with an isotropic mode (`𝔖_ff`, `D₁²`) are not resolved at
  `k ≤ 4` even for **same** helicity. At `β = 1` they were (`SPINNING_TTTT_EXT.md` §6, `k ≥ 4`).
* **`+−` at `k ≤ 5` is far from the heavy formula** under every kernel, `K*` included. The bulk part of this (`K*`) is not explained by
  the small-`t` artefact. If it is physical, the H1 problem at the binding `+−` states is worse than at `β = 1`.

**Nothing here changes the conditional numbers of `RESULT_D.md`.** The D-safe `n = 1` statement still rests on H4 [I] at `x ≳ 0.05` and
on the `+−` moments, on which this test is silent.

## 5. Open items

1. **Decide whether the `+−` low-`J` excess is numerical.**
   * The long-double check `d2` failed its own control (`W[s,s]`) and needs a correct extended-precision path for the whole chain: the
     seed radial sums in quad precision (`libcblock12_q`), and the word application in long double or mpmath at a few points.
   * Also needed: a run at `N = 12000` to move the convergence mask.
2. **Physical diagnosis.** The `+−` heavy formula at fixed small `J`: the Weyl-odd `O(1/(m A_opp))` terms as `J → 4`.
3. **More states.** `J = 0` at `m = Λ` (the D-strong `++` failure), `J = 2`, `m/Λ = 1.05, 1.1`. About 2.4 h per state per core.
4. **Large `x`.** A spinning Mellin/Barnes representation (weight-shifting operators acting on exactskel's scalar Barnes moments) is the
   only route seen to `x ~ 1`. The position-space route stops at `x ≈ 0.15`.
5. **`K*` instead of `𝔖_ff`** in `Ω^{NF,SC}` (re-deriving `N^{K*}`; `SPINNING_TTTT_EXT.md` §6). At these states, `𝔖_ff` is unresolved at
   `k ≤ 4`.

## 6. Files (`ads4/code/directgrav/`)

| file | content | output (`cache/`) |
|---|---|---|
| `gwaves.py` | resumable exact `⟨TTTT⟩` wave driver (read-only import of `tttt.TaylorWaves`) | `gw_100_1.040_4.npz`, `gw_100_1.300_6.npz`, `gw_A.log`, `gw_C.log` |
| `c1_calib.py` | position-space vs exact scalar-skeleton moments (usable `ν` range) | `c1_100_1.040_4_Df3.txt` |
| `a1_moments.py` (+ `a1_moments_ker.py`) | kernel moments with corrected weights (`𝔖_ff`, `D₁²`, `K*`, `𝔄`) vs prediction; exact skeleton at the same nodes | `a1_gw_*.pkl`, `a1_A.log`, `a1_C.log` |
| `a2_action.py` | shift in exactskel's model variables, all computed nodes (dominated by unresolved moments; diagnostic only) | `a2_summary.txt` |
| `a3_restricted.py` | resolved-node shift, clean H4 metric, `c*`, `d_crit` | `a3_summary.txt`, `a3_summary_nu3.5.txt` |
| `d1_tprofile.py` | `t`-profile of kernel moments (small-`t` contamination) | `d1_A.txt` |
| `d2_ldcheck.py` | float64 vs long-double Taylor stage at selected points (control failed: inconclusive) | `d2_A.txt` |
| `r1_clean.py` | §7: cleaned re-extraction with `t19_k2` (cut 40, noise floor, tails, wedge; spread and envelope) | `r1_gw_*.pkl`, `r1_A.log`, `r1_C.log` |
| `a4_h4_both.py` | §7: H4/H1/H2 metrics for both helicities from the cleaned values | `a4_summary.txt`, `a4_summary_sp1e-2.txt` |
| `d3_snr.py` | §7.3: signal-to-noise of the kernel densities | `d3_A.txt`, `d3_C.txt` |

## 7. Cleaned re-extraction (t19 method) and the `+−` states (`r1_clean.py`, `a4_h4_both.py`, `d3_snr.py`)

Added after `SPINNING_TTTT_EXT.md` §7 found two defects of the cached `l = 2` waves:
* the production truncation mask (`X − 10 log X > 25`) was calibrated for `l = 1`. Next to the edge, opposite-helicity nodes are wrong by up
  to `10⁴×`;
* at tiny `t` the waves are noise-dominated. The imaginary parts of the real-by-symmetry modes give the noise floor.

**Method.**
* The waves of states A and C were **re-extracted, with no new production**, using `t19_k2.cleaned`/`analyse` (imported read-only;
  `design` memoised):
  * mask cut 40 (production 25);
  * per-row start where three consecutive nodes exceed `10σ`, with power-law tails before it;
  * light-cone wedge extrapolation for `t ≥ 3`.
* The **model spread** is the maximum deviation over the variants: shifted fit nodes, two-node fits, `κ = 3, 30`, wedge from `t ≥ 5`, cut 50.
  The **envelope** adds constant tails and cut 32.
* Grid: `k = 2…7`; the functional's Gauss nodes with `ν ≤ 10.7` (14), plus `ν = 0, 1, 3`; waves `ss` and `oo`; kernels `K*`, `𝔄`, `𝔖_ff`
  (info) and the leading mode.
* Nodes kept at cut 25/32/40/50: `3165/3098/3020/2930` (A) and `3024/2948/2862/2772` (C).
* **Cost:** about 11 min per state on one core.

### 7.1 `++` after cleaning [NE]

* The table of §2 is reproduced and extended, with the cleaning changing `k ≥ 4` by `≲ 10⁻³`.
* `R_grav` (`K*` even, `𝔄` odd) has model spread `≤ 10⁻⁴` (`k ≥ 4`, `ν ≤ 4.4`) and `≤ 10⁻⁷` at `k = 7`.
* `𝔖_ff` is now also clean at `k ≥ 4`, but not at `k = 2, 3`.

**H4 metrics** (`a4_h4_both.py`, `cache/a4_summary.txt`). A node counts as resolved if its spread is `≤ 3·10⁻³` and its envelope is
`≤ 3·10⁻²` wide.

| state | `c*` (`m²(c*−1)`) | resolved nodes, `k = 2 / 3 / 4 / 5 / 6 / 7` | `max\|R/(c*R_skel) − 1\|` per `k` (3…7) | clean H4 shift `ΔF/N_s` | margin |
|---|---|---|---|---|---|
| A, Df 3 | 1.0127 (134) | 1 / 10 / 9 / 13 / 13 / 14 (`ν` up to 1.3 / 5.5 / 4.4 / 9.2 / 9.2 / 10.6) | 3.5e−3, 2.2e−3, 4.2e−3, 2.8e−3, 4.3e−3 | `−6.9·10⁻⁵` | +0.266 |
| C, Df 3 | 1.0082 (135) | 0 / 9 / 10 / 12 / 14 / 14 | 1.2e−3, 7.4e−4, 2.8e−3, 2.3e−2 (`ν ≤ 10.6`), 3.9e−3 | `−1.6·10⁻⁵` | +0.316 |

* Df = 2 gives the same to `±4·10⁻⁴`, with `c* = 1.015 / 1.010`.
* The C `k = 6` entry of `2.3·10⁻²` comes from the largest nodes, `ν ≈ 9–10.6`. At `ν ≤ 6.6` it is `≤ 3·10⁻³`. The position-space calibration
  (§1) makes `ν ≳ 9` marginal there.
* `m²(c*−1) ≈ 135` at both masses, so the common factor scales like `1/m²`, as a vertex-normalisation correction should.
* With the looser criterion (spread `≤ 10⁻²`; `cache/a4_summary_sp1e-2.txt`) `k = 2` enters at A (`ν ≤ 2.6`, deviation `8.8·10⁻³`), and
  `ΔF/N_s = −5.9·10⁻⁴`.

### 7.2 `+−` after cleaning [NE]

**`R_grav` (cleaned) [model spread] at `ν ≈ 0.03 / 1.3 / 3.5`:**

| `k` | A, `𝔄` (odd `k`) | A, `K*` (even `k`) | C, `𝔄` (odd `k`) | C, `K*` (even `k`) |
|---|---|---|---|---|
| 2 | — | `3·10³ [10⁴]`: noise | — | `12 [10]`: noise |
| 3 | 2.40 [2] / 0.91 [0.1] / 0.31 [1]: noise | — | 1.05 [0.5] / 1.009 [0.04] / 1.11 [0.3] | — |
| 4 | — | 62 [40] / −5.7 [4] / 107 [60]: noise | — | 4.7 [3] / 0.88 [0.5] / 1.85 [1]: noise |
| 5 | **0.988 [0.08] / 1.005 [0.005] / 1.032 [0.02]** | — | **1.002 [0.1] / 1.002 [0.007] / 1.017 [0.03]** | — |
| 6 | — | 1.49 [0.2] / 1.02 [0.04] / 1.37 [0.3] | — | 1.20 [0.2] / 0.997 [0.02] / 0.96 [0.05] |
| 7 | **0.9916 [10⁻³] / 0.9965 [6·10⁻⁵] / 1.0175 [2·10⁻⁴]** | — | **0.9955 [2·10⁻³] / 0.9977 [10⁻⁴] / 1.0111 [3·10⁻⁴]** | — |

* At `ν = 1.0` (not a Gauss node) the even-`k` `𝔄` contractions are also clean, e.g. `k = 4`: `1.007 [3·10⁻³]` (A) and `1.006 [3·10⁻³]` (C).
  Their `ν = 0` and `3` neighbours are not.
* **Where the moments are resolved, they agree with the heavy formula.** At `k = 5, 7` under `𝔄` the agreement is `0.99–1.03`, with
  `R − 1` of order `10⁻²`, like `++`.
* So the raw `10²–10¹⁰` excess of §3 was the mask and noise defects. The leading-mode `+−` ratios are noisy as well (`|lead|` up to
  `10⁶` at C, `k = 2`), which is the noise signature.

**H4 and H2 at the `+−` binding states** (resolved nodes only; `cache/a4_summary.txt`):

| state | resolved | reading B (exactskel model): `c*`, max dev | reading A (`D = P_J`): `c*`, max dev | `ΔF/N_s` | margin |
|---|---|---|---|---|---|
| A, `+−` J=4 | `k = 7`, 10 nodes (`ν ≤ 5.5`) | 1.0008, **1.6·10⁻³** | 1.0050, 4.4·10⁻² | `−2·10⁻¹⁷` | +0.241 (D-safe) |
| C, `+−` J=6 | `k = 7`, 10 nodes; `k = 5`, 2 nodes (`ν ≤ 0.8`) | 1.0006, **8·10⁻⁴** (`k = 5`: 2.8e−3) | 1.0035, 2.8·10⁻² | `−4·10⁻¹⁵` | +1.002 |

* With spread `≤ 10⁻²`, `k = 5` adds 4 nodes (`ν ≤ 1.9`). Reading B then gives deviation `7.7·10⁻³` (A) and `2.5·10⁻³` (C).
* **H2 (reading B) at these near-gap low-`J` states [NE].** At `k = 7` the exact `+−` moment follows the reading-B-transferred skeleton to
  `≤ 1.6·10⁻³`, against `3–4 %` for reading A.
* **H4 for `+−`.** It is verified only at `k = 7` (and two `k = 5` nodes). These carry essentially none of the NF02 weight (`p^{10}`), so the
  action shift is `~10⁻¹⁵ N_s`. **No statement on the `+−` action is possible.**

### 7.3 Numerical vs physical: the unresolved `+−` moments are noise-limited (`d3_snr.py`, `cache/d3_*.txt`) [NE]

`|M|/Σσ|W|` compares the integrated moment with the noise floor integrated over the grid. It is a conservative ratio: the noise floor
summed coherently. For incoherent noise the ratio is about 50× larger.

| `k` | `++` (`K*`, `𝔄`) | `+−`, A: `K*` / `𝔄` | `+−`, C: `K*` / `𝔄` |
|---|---|---|---|
| 2 | 15–80 | 2·10⁻³ / 2·10⁻³ | 3·10⁻⁵ / 10⁻⁶ |
| 4 | 4·10⁴ | 6·10⁻³ / 4·10⁻⁴ | 9·10⁻⁵ / 2·10⁻⁵ |
| 5 | 5·10⁵ | 6·10⁻² / 10⁻² | 10⁻³ / 10⁻³ |
| 6 | 3·10⁶ | 1.4 / 1.4 | 0.12 / 0.13 |
| 7 | 10⁷ | 140 / 150 | 17 / 19 |

* **The opposite-helicity kernel moments at low `J` are a cancellation.** Most nodes are individually above `10σ`, but the integrated
  moment is `10²–10⁵` below the integrated noise floor for `k ≤ 4`, and `≲ 1` at `k = 5, 6`.
* At those `k` the float64 waves cannot determine the `+−` moment. The values in the raw and cleaned tables there are numerical, not
  physical.
* Where `|M|/Σσ ≳ 1` (`k = 6, 7`, and `k = 5` at A with the cleaning), the moments agree with the heavy formula.
* **So the evidence is that the `+−` excess was numerical [NE].** No sign of a physical deviation from H1 was found at these states.
* **Not done:** the requested mpmath end-to-end check at one node. A correct extended-precision path needs the quad-precision seed radial
  sums plus the Taylor-mode word application (196 word pairs) in mpmath. That is beyond the 1 h limit. The signal-to-noise analysis is
  used instead; it settles "noise-limited", not the true value.

### 7.4 Revised grades

**H4:**
* `++`: unchanged, [NE] at `x ≲ 0.03–0.09`. Graviton moments equal skeleton moments up to `c* = 1 + 135/m²`, to `≤ 4·10⁻³` for
  `k = 3…7`. The action shift is `≤ 6·10⁻⁴ N_s`, but coverage is `≤ 0.4 %` of the remainder norm.
* `+−`: [NE] only at `k = 7` (and `k = 5`, `ν ≲ 2`), with reading B. The action is not covered.
* **Still [I] where the action is decided** (`x ≳ 0.1`).

**H1 at the binding states:**
* `++`: holds, `R − 1 = 0.4–2 %`, both `K*` and `𝔄`, `k = 3…7`.
* `+−`: holds where resolved (`𝔄`, `k = 5, 7`).
* `+−` under `K*` at `k ≤ 6`, and all `+−` at `k ≤ 4`: noise-limited in float64, undetermined.

**H2 [NE] at near-gap low `J`:** reading B over reading A, at `k = 7`.

**What would close `+−`:** an extended-precision production (quad seed sums plus long-double or mpmath Taylor stage) at one or two states,
or a formulation that avoids the cancellation, e.g. extracting the opposite-helicity moments from a better-conditioned mode combination.

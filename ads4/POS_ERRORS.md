# POS_ERRORS: which parts of the per-state positivity evidence are models, how large their errors could be, and what can be error-bounded

*23 Nov 2026. Unreviewed, not committed. No existing file was edited. `A1PRIME.md` and `code/a1p/` were not touched.*
* *New code:* `ads4/code/poserr/`: `curv_lemma.py`, `eng_orders.py`, `large_m_map.py`, `q_large_m.py`, `analyze_D.py`, `analyze_Q.py`,
  `analyze_gap.py`, `run_all.sh`.
* *Caches:* `ads4/code/poserr/cache/` (git-ignored by its own `.gitignore`).
* *Read-only imports:* `spinfunc/ee1000_common.py` (hybrid states: flat NF02 with exact Wigner functions plus exact leading-Regge
  curvature), `spinfunc/bg_eeq.py`, `spinfunc/bj_ee.py` (engine value functions only).
* *Read-only data:* `spinfunc/cache/bg_eeq_D{300,1000}_J0_K4.pkl` and `bj_eeq_D*_J*_d*.pkl` (the existing analytic-radial engine
  images), `prl/data/q3_full_delta0.2.json`.
* *Compute:* about 20 CPU-minutes, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.
* No crystal or engine image was recomputed. Every engine number below is read from caches that already existed.

Grades are those of UNIFIED_THEOREM_v2 §0. One label is added:
* **[NE±]**: numerical evidence with an explicit error budget. Every model input has a stated safety factor, and the conclusion
  survives the stated multiple of it. It is weaker than [D] and stronger than [NE + M].

"Leverage" `L` = margin ÷ (plausible size of the model error, or a stated envelope of it). `c₂` = `m² ×` relative correction
(UNIFIED_THEOREM_v2 §2.5). Unless stated otherwise, D means `Ω^EE_θ` and Q means `Ω^Q_up` with the δ = 0.2 weights.

---

> **Update (NEARGAP2.md, unreviewed): Result D's binding state is now error-bounded [NE + M] → [NE±]; light J = 0, 2 states are positive up to just below the gap.**
> * **Binding state** (J = 5, `+−`, twist edge 1.0049925Λ, Δ_g = 1000), on a K₂ row at its own mass (Δ = 1005, 36 nodes): λ/N_s = +0.20088. Error budget, sign-free with ×3 on every model term: 2.6·10⁻³, so **λ/N_s ≥ +0.198, leverage L = 78** (130 with the λ-level truncation). Certified flat floor minus the Lemma-C bound gives ≥ +0.197 with no quadrature input. 9 direct RQMC nodes agree with K₂ (λ moves < 10⁻⁶); the parity-odd analytic-radial engine agrees on ν = 12.5–38. At Δ_g = 300: ≥ +0.146, L = 3.5 [NE±], with no direct check.
> * **Near-gap J = 4, 6** (own rows): [NE±], L ≥ 22 / 5.0 sign-free. **J = 8**: [NE±] on the λ-level budget only (L = 4.3–8.4). **J = 12, 16**: [NE + M]; the budget does not close (L ≈ 1–3), and the next Laplace order is needed. The published chain differs from the own-row values by at most 0.015 N_s.
> * **H_cont, J = 15 twist edge**: computed, +0.6900 [NE + M] (was [M]). **H_edge**: J = 18, 19, 24, 25 shift by ≤ 0.032 [NE + M]. **J = 29 is not resolved:** one K₂ outlier at x = 0.8 drives the plain w2 odd chain to −3.62, while the direct nodes give +1.39 ± 0.03. It stays [M], and is carried into the high-spin task.
> * **Light J = 0, 2 (A1′):** positive at every Δ_L up to Δ_g − 1 in the leading-Regge model (F/N ≥ 0.788 / 0.638), and Lemma-C certified there for Δ_L ≥ 4.5. Own rows at Δ_L = 990 give L = 4362 / 913 [NE±]; Δ_L ≳ 0.7Δ_g gives [NE + M]. Whether the leading-Regge form is the true AdS image of a light block at Δ_L = O(1) remains [D-sketch].
> * **Correction to POS_ERRORS §3.1:** odd J has only the parity-odd `+−` channel, a 1×1 block; the sign-free tolerance at the binding state is 3.6·10⁴ c₂, not 6.6·10³.

> **Update (LARGESPIN_ANALYTIC.md, unreviewed): no large-spin hypothesis is replaced by a proof. H_R′ is upgraded [M] → [D-sketch + M]; H_LJ(1000) and H_bulk keep their grades. A possible convention problem (N2) affects H_LJ's leverage.**
> * **Common missing input:** the beyond-leading-Regge correction inside the Regge window (ν ≲ 3/β). It needs the ⟨TTTT⟩ heavy-block Regge formula at subleading order (the spinning analogue of CHMRS (3.17)–(3.20)), which is not derived.
> * **Lemma I [P]:** for an image ∫P_{−½+iν}(cosh t)σ(dt), the action is ∫Ω_LR(t)σ(dt). This gives explicit tolerance profiles for H_LJ (≥ 0.062 of the ring weight at t = β) and tolerated sizes of candidate correction shapes (ring amplitude a ≥ 80/m, ν² a ≥ 37/m², index-shifted a ≥ 10.8/m, …). The remaining model input R_win (non-ring window correction): exact scalar ≈ 10 c₂ (×38); graviton J = 100 ≈ 550 c₂ (leverage 1.1–1.5); no data at J ≥ 150.
> * **H_R′ [D-sketch + M] for b̂ ≥ b̂_A(m).** Lemma C′ [P + DV], a sign-free curvature bound using the decay of Φ_k, makes the hybrid floor positive for b̂ ≤ 119 at every m/Λ from 1.3 to 20 (the old Lemma C reached b̂ < 10). With ROUTEA, the leading-Regge hybrid is rigorous for all b̂ ≤ 2Λ when m ≥ 1.5Λ. The AG saddle test passes on 5 states. 5 ≤ b̂ < b̂_A stays [NE + M].
> * **H_bulk [M], restated:** the bulk correction may reach ≥ 3.9 % of the flat envelope (9.8 % at 1.3Λ) before any state fails. There is no analytic route.
> * **N2 (open, may weaken H_LJ):** Ω^EE's k = 2 part acts on raw images with the de-shift weight X, not ψ₂, and |X−ψ₂|/ψ₂ = 0.3–0.9 at ν ≤ 10. HLJ's T_B used ψ₂ inside the window. With |X|, the leverage at (1100, 100) drops to ≈ 1.1, and the sign-free ring-amplitude tolerance (3.2) is below the measured excess (2.7–5.9). Which convention the fw/pw Γ data use must be checked. **Until then quote the H_LJ(1000) leverage as 1.1–4.8 (not 2.7–4.8).**
> * **Corrections:** ROUTEA's certified window [91.2, 2Λ] holds only for m ≥ 1.5Λ (UT3 §5.2 and PE §2.4 quote it for m ≥ 1.3Λ). HLJ §5.2's [D-sketch] (subleading terms are P_{l+n}(cosh β) of the same β) does not by itself make window corrections harmless (N3). PE's "no sign-free argument can work" is true for δ but false for the curvature.

## 0. Verdict

| item | result | grade |
|---|---|---|
| **Inventory** (§1) | 13 rows for D at Δ_g = 1000 (10 of them [M], [I] or [D-sketch]) and 9 for Q, each with failure state, margin, plausible error, leverage and cheapest route. By leverage, the weakest are `H_LJ(1000)` (×2.7–4.8, refined window only; the sign-free window is violated at (1100, 100)) and the bulk of `H_R` at large `b̂` (sign-free tolerance below 10 `c₂`). The strongest are `H_OO` (≥ 10³), `H_smear^Q` (≈ 8·10³) and Q's `H_RS` truncation (×85 / ×850). | — |
| **Lemma C** (§2.1) | `|P_{−½+iν}(cosh β) − J₀(2ν sinh(β/2))| ≤ (1 − sech(β/2)) + min(2, ν(2 sinh(β/2) − β))` for all real `ν`, `β`. The flat Bessel argument is exactly `b̂p`, so this bounds the hybrid's AdS curvature with no special-function evaluation. Checked against 30-digit values: worst ratio 0.73. | **[P]** (given the Mehler–Dirichlet formula [Q]) + [DV] |
| **Engine structure at large m** (§2.3) | At J = 0 the 1/Δ engine's order-k coefficients are Δ-independent: 0.04 % (k = 0) to 2.3 % (k = 4) between Δ = 300 and 1000. In `c₂` units at **every** `m`: <br>• order 1 = a ν-flat multiplicative factor `1 + 20.97/m` (EE2) or `1 + 26.96/m` (EE4), plus ≤ 108 `c₂` (EE2) / ≤ 82 (EE4) of ν-dependent remainder; <br>• order 2 ≤ 333 / 445 `c₂`; <br>• orders ≥ 3 fall like `12/Δ` per order. <br>At spin J (cached `bj` images, J = 28–100, Δ = 336–2·10⁴): the ν-flat factor is J-independent (`20.9 m` at Δ = 2·10⁴). The ν-dependent order-1 part is ≤ `450 + 100J` `c₂` on ν < 12. **The D̂ series is converged only for ν < 12–20** (`|N₅ − N₄| ≤ 0.7 c₂` below 12; unconverged above 20). | [D given RD] (expansion, HRS_CLOSE) + [NE] (coefficients) |
| **D: large mass** (§2.4) | **There is no `m_*` above which all J are covered.** Fixed-`J` large `m` is covered, and so are bounded `b̂ = (2J+1)Λ/m`. Bound: certified flat floor `N/5`, minus Lemma C, minus the engine envelope with ×3 safety, all sign-free. **On the strip `S_D = {m ≥ 1.3Λ, b̂ < 5}` (J ≤ 2 at 1.3Λ … J ≤ 49 at 20Λ): `L ≥ 8.7`** (10.3, 12.5, 19.3, 32, 64, 129 at 1.5, 2, 3, 5, 10, 20Λ), growing ∝ m. The one model input (the bulk ν ≥ 20 envelope at J > 0) carries ≤ 15 % of the cost and may be wrong by **×76** before any state fails. L ≥ 3 extends to `b̂ ≤ b̂_A(m) = 8.5, 8.7, 10.5, 13, 15, 15.3, 13.8`, but there the bulk envelope carries 40–95 % of the cost. **At larger `b̂` the sign-free tolerance falls to 10²–10⁻¹ `c₂`, so no sign-free argument can work.** That region stays `H_R` [M]. | `S_D`: **[NE±]** (certified + P + D given RD + NE; bulk model needed to ×76). Rest: unchanged |
| **Q: large mass** (§2.5) | Q's forward term `a·2J(J+1)u³` keeps `N` large at large `b̂`, so sign-free bounds work **at every `b̂ ≤ 400`**. At 1000, the sign-free tolerance on the smeared parts is ≥ 2·10⁴ `c₂` on every computed state with `m ≥ 1.3Λ`. The derivative-piece failure threshold is `s_fail′ ≥ 4.8` (1.3Λ), rising as `u⁻³`, against `|s_d| ≤ 0.6`. `s_d` is measured to be m-independent from the gap to 20Λ (EESTEP §5). On the engine-covered states (`J²/Δ ≤ 1`) the floor 9/50 gives `L ≥ 3.95` for `1.3 ≤ m/Λ ≤ 3` and `L ≥ 1.65` up to 20Λ; the computed flat floor gives `L ≥ 5.0`. **At 300 the same budget gives `L ≥ 1.2`** (`m ≤ 2Λ`) **and `L < 1` at `m ≥ 3Λ`, `b̂ ≥ 5`**: not converted. | 1000, `m ≥ 1.3Λ`, `J²/Δ ≤ 1`: **[NE±]**, `L ≥ 1.65` (≥ 5 with the computed flat floor). 300: unchanged |
| **Near-gap window** (§3) | The sign-free route converts only J = 0, 2 at 1000 (`L = 7.0, 6.0`). At J = 4, 5, 6 it gives `L = 2.7, 1.8, 1.0`, and **0.09 or less from J = 8 on**, with ≥ 63 % of the cost in the bulk. **The binding state (J = 5, `+−`, +0.2009) is not error-bounded by this route.** The window stays [NE + M]. §3.2 lists exactly which states are model-only. | — |
| **Calibration** (§3.3) | The documented beyond-Regge corrections at `m ≥ 1.3Λ` equal **uniform** corrections of only 20–50 `c₂` (300) and 1.6–3.3 `c₂` (1000). The sign-free envelope is 10²–10⁴ times larger. **The true corrections cancel (they oscillate with the partial wave).** That cancellation is what a proof at large `b̂` would need; it is not derived. | [NE] |
| **Large-J tail** (§4) | This hits the `H_LJ` resummation wall. At fixed Λ and growing J, the state approaches the twist edge (`b̂ → m`). The sign-free tolerance falls like `b̂⁻³`, and the engine's structure series is asymptotic once `J²/Δ ≳ 2` (HLJ_CLOSE). **Stopped.** | [I] (unchanged) |

**Bottom line.** Three results:
* **Flat limit:** fully certified (CERT_CLOSE). This note adds nothing there.
* **Small `b̂`, `m ≥ 1.3Λ`:** now [NE±] for D at 1000 and for Q at 1000. The margins are explicit, every model input carries an explicit
  safety factor, and only the 1/Δ engine and the [D given RD] structure behind it are used.
* **Everything else:** still [NE + M] or [M].
  * Near-gap states.
  * Large `b̂` at any mass.
  * Near-gap large J (`H_LJ(1000)`).

  The reason is structural. A per-ν sign-free error bar cannot certify a state whose normalisation `N_s ∝ b̂⁻³` comes from
  cancellations. Doing so needs the oscillating form of the corrections, which is modelled, not derived.

---

## 1. Inventory of the model and open components of the positivity evidence

Columns:
* **what fails** = the state that would flip first if the component were wrong, with its margin (min `λ/N_s` or FF/N);
* **plausible error** = the largest discrepancy measured or estimated anywhere in the record;
* **L** = margin ÷ plausible error;
* **route** = the cheapest way to make it error-bounded.

### 1.1 Result D, Δ_g = 1000

| # | component | grade now | what fails (margin) | plausible error | L | cheapest route to an error-bounded version |
|---|---|---|---|---|---|---|
| D1 | **`H_LJ(1000)`**: near-gap large spin (even `J ≥ 50`, odd `J ≥ 33` to `J ≈ 344`, `m ≤ 1.3Λ`; and even J = 26, 30, 34, 38, 42–46) | [D-sketch + M] | `(1100, 100)` window (`T_B = 1409 c₂`); κ-model state values ≥ +0.90 | measured window remainder ≈ 250 `c₂` (refined) / 450 (full); **sign-free form violated at (1100, 100)** (−2154 ± 43) | **2.7–4.8** (refined, DG P2) | **none cheap.** The engine cannot decide `H_win` at `J²/Δ ≳ 9` (HLJ §5). It needs a resummation of the helicity-structure series, or below-fold thimble contours. |
| D2 | **`H_R`**: beyond-leading-Regge at `m ≥ 1.3Λ` | [M] | at large `b̂`: hybrid ≥ 0.93 (`b̂ ≥ 20`), with `N_s ∝ b̂⁻³`; at small `b̂`: hybrid ≥ 0.32 (`b̂ ∈ [10, 20)`, 1.3Λ) | measured ≤ 0.005 `N_s` (`J ≤ 32`), ≤ 0.02 (`J ≤ 40`); scalar skeleton `|S|/F_hyb ≤ 122/m²` (`b̂ ≤ 200`) | ≥ 15 (`J ≤ 40` measured); ≥ 150 (scalar transfer) | **small `b̂`: done here (§2.4, strip `S_D`).** Large `b̂`: needs the oscillating structure (the κ-model of SUBLEADING_REGGE) derived to an explicit order, or crystal rows at `m = 2–5Λ`, `J ≤ 100` |
| D3 | **`H_M`**: the near-gap chain (`K₂` rows, local/non-local de-shift, value held below ν ≈ 8.3) on the computed states | [M], validated [NE] on rows with direct data | binding `+−` J = 5 at the twist edge, **+0.2009**; J = 4: +0.2041; J = 6: +0.2296 | the J = 5 correction itself is ≤ 10⁻³ (hybrid 0.20090 = chain); edge shift ≤ 0.006; `K₂` vs direct ≤ 0.05 `N_s` (`(1040, 40)`) | ≈ 30 (J = 5); 11 (J = 48: `0.915 ± 0.08`) | **Engine values at the near-gap states J = 4–16** (`J²/Δ ≤ 0.26`, convergent regime; `bj` depth 3, ν < 12, about 25 min per J). This replaces the small-ν part of the chain by [D given RD] + truncation. The bulk (ν ≳ 40, light-cone edge) still needs rows. |
| D4 | **`H_edge`**: `x ∈ (0.8, Λ/m]` on even J = 18–32, odd 19–29, J = 48 | [M] there; [NE + M] where computed | those states, margins 0.84–0.98 | largest computed shift −0.016 `N_s` (near gap) | ≥ 50 | edge rows (`w1_rows`, existing code) at the listed J |
| D5 | **`H_OO`** for J ≥ 12: OO block via `Γ_OO = Γ_++` | [D given the identity]; identity [P] at order 0, [NE] beyond | OO J = 4: **+0.215** (computed); J ≥ 12: the `++` diagonal (≥ 0.67) | identity residual ≤ 0.05 / 0.12 `c₂` | ≥ 10³ | a derivation of `X = 0` beyond order 0 (analytic) |
| D6 | **`H_cont`** between sampled masses; **J = 15 twist edge** not sampled | [NE]; J = 15 [M] | J = 15 at its twist edge: hybrid ≈ **+0.69** | the J = 14 and 16 twist-edge corrections are ≤ 0.05 | ≥ 10 | the J = 15 twist-edge row (CL protocol), minutes |
| D7 | `K₂` band `J²/m = 1.6–2.2` (near-gap J ≈ 41–47) | [I] (inside D1 except J = 48) | J = 48: +0.915 ± 0.08 | `K₂` fails ×2–5 at `J²/m ≥ 2.2` | (inside D1) | direct sub-rows at J = 42–46 |
| D8 | value held below the lowest direct node (ν ≈ 8.3); C3H below 8.3 | [M] | near-gap J ≥ 24 | measured small-ν `|Γ| ≤ 350 c₂` (J = 64) | `C*_25` = 1.2·10⁶ (J = 8) … 1.1·10⁴ (J = 48) `c₂` → **≥ 30** | engine values below 8.3 (as D3) |
| D9 | odd J: local de-shift only (z4's non-local term is even-J only) | [M] | binding J = 5 | J = 5 correction ≤ 10⁻³ | ≥ 100 | implement the odd-J non-local term (code change) |
| D10 | **RD** (Regge dominance of the crystal integral; premise of every engine use) | [D-sketch] (EE ν < 10); [DV] scalar | wherever the engine is used | far region ≲ `e^{−26}` after the volume factor | ≫ 10⁶ | a written far-region bound for the EE kernel |
| D11 | flat limit (NF02) | **[certified]** `≥ N/5`, attained to 10⁻⁸ (CERT_CLOSE) | — | — | — | done |
| D12 | in-sample states (300: J = 18, 20, 24, 28) | disclosure | — | — | — | — (1000 has none) |
| D13 | `σ_hU` | [I], not used | — | — | — | — |

### 1.2 Result Q, Δ_g = 300 and 1000

| # | component | grade now | what fails (margin) | plausible error | L | cheapest route |
|---|---|---|---|---|---|---|
| Q1 | **`H_RS`** at the binding pair (`++`/OO J = 0, gap) | [D given RD] (HRS_CLOSE); truncation [NE] | FF/N **+0.1427 / +0.1428**; `s_fail − s_d = 1.27` | truncation `|Δs_d| ≤ 0.015 / 0.0015` (the a⁴ term, ×24 the a⁵) | **×85 / ×850** | already error-bounded given RD (see D10) |
| Q2 | **HYB**: scalar skeleton + leading Regge on every state other than the measured ones; `β > 30` | [M] (images [DV]) | full-region HYB minimum 0.197 (Stage 6) | HYB − exact skeleton ≤ 6·10⁻⁴ (300) / 5·10⁻⁵ (1000) | ≥ 300 (relative to the scalar skeleton; the transfer to gravitons is the model) | **large `m`: done here at 1000 (§2.5).** Near gap: as Q3 |
| Q3 | **`H_fwd^Q`**: unmeasured forward channels; the 300 J-drift | [NE + M]; drift [I] | the channel with the smallest `s_fail` among the unmeasured ones (≥ 11) | measured `|s_d| ≤ 0.59`; drift −0.07 (favourable) | **≥ 19** | engine `s_d` for `++` J = 10, 12, OO J = 4 (bg engine, minutes each) |
| Q4 | **`H_smear^Q`**: graviton corrections to the smeared k = 2, 3 parts | [NE + M] | binding states, tolerance ≈ 24 % of B₂ | ≤ 3·10⁻⁵ of B₂ | ≈ 8·10³ | done in effect; at large `m` see §2.5 |
| Q5 | **`H_LJ^Q`**: large-J near-gap corrections | [M] | tolerance ≥ 1.1·10⁴ `c₂` (`+−` J = 8, 300), ≥ 1.6·10⁵ (J ≥ 16) | measured 10²–10⁴ `c₂` | ×1–10³ (worst at 300, J = 8) | engine values at `+−` J = 8–12, 300 |
| Q6 | `H_OO^Q` | [D given `H_RS`] | OO J = 0 = `++` J = 0 | ≤ 5·10⁻⁶ in `s_d` | ≫ 10³ | — |
| Q7 | C5 piece beyond leading Regge | inside HYB [M] | derivative combination | tolerance `G* = 1.18·10⁵ c₂` (300) | ≥ 10 for any correction below 10⁴ `c₂` | none cheap (no k = 5 engine) |
| Q8 | flat limit (Q3\*) | [certified] `≥ 9/50` on most regions; [NE] in two (float minima 0.304, 0.997) | — | — | — | Lemma B for `+−`, `m ≤ 2M`, `J ≥ 400` (CERT_CLOSE §2.3) |
| Q9 | de-shift ε-DS | [DV] | — | `ε_W ≤ 4.2·10⁻¹¹` | ≫ 10⁶ | done |

**Ranking by leverage (weakest first):** D1 (2.7–4.8), Q5 at 300 (≈ 1–10), D2 at large `b̂` (sign-free < 1; measured ≥ 15), D6 (≥ 10), Q3 (≥ 19), D3 (≈ 30 at the binding state), D8 (≥ 30), D4 (≥ 50), Q1 (85 / 850), D9 (≥ 100), Q2 (≥ 300), D5 (≥ 10³), Q4 (8·10³), D10.

---

## 2. The large-mass region

### 2.1 Per-state decomposition and Lemma C

For a heavy state `(J, m)` and helicity `h`, the diagonal entry of `M_Ω` is exactly

    diag_h = F_h(flat) + curv_h + δ_h,

with these pieces:
* `F_h` = the stripped NF02 flat action, with exact Wigner functions;
* `curv_h = 2∫ψ₂[P_{−½+iΛp}(cosh β) − J₀(b̂p)] + σ_h u∫ψ₃[P − J₀]` = the leading-Regge AdS curvature (the hybrid's AdS part,
  `ee1000_common.hyb`);
* `δ_h = ∫ψ₂ env₂ Γ₂/m² + u∫ψ₃ env₃ Γ₃/m²` = the beyond-leading-Regge correction (`Γ` in `c₂`).

The flat and curvature off-diagonals vanish. The off-diagonal is then `∫ψ₂ eo₂Γ_off/m² + …`, with the same structure
(`ee1000_common.dF_off`).

So `λ_min ≥ min_h (diag_h/N_h) − |off|/√(N₊N₋)`, and the inputs are graded as follows:
* `F_h ≥ N_flat/5`: **[certified]** for every flat state (CERT_CLOSE §1).
* `curv_h`: an explicit special-function integral. It is evaluated by quadrature [NE] or bounded by Lemma C [P].
* `δ_h`: the only model content.

**Lemma C [P].** For real `ν ≥ 0` and `β > 0`, let `B = 2 sinh(β/2)`. Then

    |P_{−½+iν}(cosh β) − J₀(νB)| ≤ (1 − sech(β/2)) + min(2, ν(B − β)).

*Proof.*
1. Start from Mehler–Dirichlet [Q]: `P_{−½+iν}(cosh β) = (√2/π)∫₀^β cos(νt)dt/√(cosh β − cosh t)`.
2. Substitute `s = 2 sinh(t/2)`. Then `cosh β − cosh t = (B² − s²)/2` and `dt = ds/cosh(t/2)`, so
   `P = (2/π)∫₀^B cos(νt(s)) ds /(cosh(t/2)√(B² − s²))`.
3. Also `J₀(νB) = (2/π)∫₀^B cos(νs) ds/√(B² − s²)`.
4. Pointwise in `s`: `|cos(νt)/cosh(t/2) − cos(νs)| ≤ (1 − sech(t/2)) + min(2, ν(s − t))`.
5. Both terms on the right increase with `t`, because `d(2 sinh(t/2) − t)/dt = cosh(t/2) − 1 ≥ 0`.
6. Finally `(2/π)∫₀^B ds/√(B² − s²) = 1`. ∎

**Why it applies to the hybrid.** `sinh(β/2) = (J + ½)/m` exactly, so `νB = Λp(2J + 1)/m = b̂p`: Lemma C compares exactly the two
functions in `curv_h`. Hence

    |curv_h| ≤ ∫(2|ψ₂| + u|ψ₃|)[(1 − sech(β/2)) + min(2, Λp(B − β))]dp,      β ≈ b̂/Λ.

Check (`curv_lemma.py`, 30-digit mpmath, `β ∈ [10⁻³, 4]`, `ν ∈ [0, 10³]`): worst `|P − J₀|/bound = 0.73` **[DV]**.

What it buys (`large_m_map.py`, column `rig` = `(N_flat/5 − |curv| bound)/N_s`):
* For `b̂ < 10` the certified floor survives the curvature bound almost intact: `rig ≥ 0.1994` at 1000 (≥ 0.1933 at 300).
* For `b̂ ≥ 20` the sign-free bound is useless (`rig < 0`). The true curvature is ~10³ times smaller: it is the oscillating cancellation
  of an `O(1)` integrand.

There I use the computed curvature instead [NE]. It is positive on every record (ROUTEA_LARGE_B), and the floor
`(N_flat/5 + min(0, curv))/N_s` is ≥ 0.110 at every computed state (≥ 0.19 for `b̂ < 160`). ROUTEA_LARGE_B's certified windows
(`F_hyb > 0` on `b̂ ∈ [91.2, 2Λ]`) remain the rigorous statement at large `b̂`.

### 2.2 What a sign-free bound on `δ` can and cannot do

Define `W_h = ∫(|ψ₂|env₂ + u|ψ₃|env₃)/(m²N_h)` (and `W_off` likewise). A uniform sign-free `|Γ| ≤ G` on every ν moves `λ` by at most
`G(W_h + W_off)`. The tolerance `G_tol = floor/(W + W_off)` at Δ_g = 1000 (`analyze_D.py`, min over states in each band, `c₂`):

| `m/Λ` \ `b̂` | [0,5) | [5,10) | [10,20) | [20,40) | [40,80) | [80,160) | [160,400] |
|---|---|---|---|---|---|---|---|
| 1.3 | 1.8·10⁵ | 1.4·10⁴ | 730 | 76 | 8.9 | 1.4 | 0.09 |
| 2 | 2.1·10⁵ | 2.1·10⁴ | 2.0·10³ | 250 | 36 | 4.7 | 0.31 |
| 5 | 1.1·10⁶ | 1.1·10⁵ | 1.4·10⁴ | 2.1·10³ | 280 | 36 | 2.4 |
| 20 | 1.6·10⁷ | 1.8·10⁶ | 2.3·10⁵ | 3.5·10⁴ | 4.6·10³ | 600 | 40 |

(At 300 every entry is about 11× smaller: `G_tol ∝ m²` at fixed `m/Λ`.)

Two facts follow.
* `G_tol ∝ m²` at fixed `b̂` [D: `N_s` depends on `(b̂, u)` and `W ∝ m⁻²`]. So a correction that grows like `m` (H_mono) loses
  eventually **at fixed `b̂`**.
* `G_tol ∝ b̂^{−3…−5}` at fixed `m`. This is SUBLEADING_REGGE's "smooth corrections are amplified by `b̂⁵/Λ²`". Beyond `b̂ ≈ 20` no
  per-ν error bar of any plausible size (the measured `|Γ|` reaches 10²–10⁴ `c₂` locally) can be propagated sign-free. Positivity
  there rests on the corrections oscillating with the partial wave (§3.3). That is a property of the true corrections that is
  modelled ([M]: κ-model, SUBLEADING_REGGE), not derived.

### 2.3 The size of `δ` at large mass: the 1/Δ engine, read order by order

**J = 0** (`eng_orders.py`; existing images `bg_eeq_D{300,1000}_J0_K4` and their Δ = 2·10⁴ reference).

The engine is `I(ν; Δ) = e^{l(ν)} Σ_k Δ^{−k} M_k(ν)` [D given RD, HRS_CLOSE Prop. RS].
* **The coefficients `M_k` are Δ-independent:** `max|M_k(300) − M_k(1000)|/max|M_k|` = 0.04 %, 1.0 %, 1.7 %, 2.1 %, 2.3 % for
  k = 0 … 4 [NE].
* Successive maxima have ratio `|M_{k+1}|/|M_k| ≈ 21, 13, 12.7, 11`, so the series is geometric with ratio ≈ 12/Δ (1/25 at 300,
  as quoted in HRS_CLOSE).

In `c₂` units (`m²|T_k|/env|T₀|`, max over ν windows):

| kind | order-1 ν-flat factor | order-1 remainder [0,2) / [2,10) / [10,41) | order 2 | order 3 (Δ = 300 → 1000) | order 4 |
|---|---|---|---|---|---|
| EE2 (D's k = 2) | **20.895 / 20.969 / 20.998 · m** (Δ = 300 / 10³ / 2·10⁴) | 108 / 29 / 0.85, **identical at all three Δ** | 322–333 / 271 / 40, **identical at all three Δ** | 13.6 → 4.2 (∝ 1/m) | 0.50 → 0.046 (∝ 1/m²) |
| EE4 (Q's k = 4) | 26.86 / 26.96 / 27.00 · m | 81 / 17 / 1.2, identical | 434–445 / 395 / 266, identical | 18.5 → 5.7 | 0.68 → 0.06 |

So, explicitly, at J = 0 and **any** `m ≥ 300`:

    Γ(ν; m) = 21.0 m · (image)/env   [ν-flat, multiplicative]   +  g(ν)  +  O(14·300/m),      |g| ≤ 108 + 333 c₂  (ν < 2),
                                                                                         ≤ 29 + 272 (2 ≤ ν < 10),  ≤ 0.9 + 40 (ν ≥ 10),

with `27.0 m` and `82 + 445` for EE4. The ν-flat term multiplies the whole image, so it shifts `λ/N_s` by at most `21/m` in either
direction (0.016 at 1.3Λ). It is also positivity-neutral to first order when it is common to the k = 2 and k = 3 parts. This turns
`H_mono` at J = 0 ("Γ grows at most like `m`") from [M] into **[D given RD] + [NE]**: the growth is exactly the ν-flat
`21 m` term, and the rest is m-independent.

**Spin J** (`eng_orders.py`, cached `bj` images; worst entry of the 2×2 conserved basis, absolute deviation over a common smooth
envelope).

| `(Δ, J)` | `J²/Δ` | ν-flat factor | order-1 remainder [0,2) / [2,6) / [6,12) / [12,20) (`c₂`) |
|---|---|---|---|
| (2·10⁴, 28) d3 | 0.04 | 20.93 m | 2.8·10³ / 620 / 150 / 760 (order 2: 150 / 25 / 47 / 1.8·10³) |
| (2·10⁴, 64) | 0.20 | 20.95 m | 4.1·10³ / 3.0·10³ / 1.5·10³ / 2.3·10³ |
| (2·10⁴, 100) | 0.50 | 20.91 m | 6.7·10³ / 5.0·10³ / 2.5·10³ / 3.7·10³ |
| (336, 28) d3 | 2.3 | 18.1 m | 2.2·10³ / 430 / 91 / 7.4·10³ (order 2: 580 / 560 / 240 / 590) |
| (1100, 100) | 9.1 | 16.9 m | 8.0·10³ / 2.2·10³ / 260 / 670 |

The D̂ series (N = 4 vs 5; `(1100, 100)`, `(336, 28)`, `(336, 32)`) differs by ≤ 0.7 `c₂` for ν < 12, by ≤ 10 `c₂` for ν < 20, and by
10³–10⁵ `c₂` for ν ≥ 20. **For J > 0 the engine says nothing at ν ≥ 20**, and nothing at all at ν > 41 for any J: that is the
bulk-point regime `ν ~ xm`, outside the Regge expansion. The larger numbers on 12 ≤ ν < 20 in the depth-2 files (2.3·10⁴ at J = 28, 2.6·10⁴ at J = 32)
are depth truncation: the J = 28 depth-3 file gives 760 there.

**Envelope E1′ used below** (×3 on every engine term; the ν-flat factor ×1 because it is measured to 10⁻³):
* ν-flat: relative `27/m` on the whole image, sign-free;
* ν < 20, `J²/Δ ≤ 1`: `|Γ_non-flat| ≤ g(J) = 3(450 + 100J)` `c₂`. This covers every entry above: J = 0 → 450, J = 28 → 2.9·10³,
  J = 100 → 8.0·10³. Grade [NE] coefficients on [D given RD] structure.
* **ν ≥ 20 (the bulk), J > 0:** the same `g(J)`, as a **model** [M]. Its support:
  * the J = 0 engine to ν = 41: ≤ 40 `c₂`;
  * the exact scalar skeleton in the bulk-point regime: `m²(exact/BP − 1) ∈ [−80, −4]` (UNIFORM_REMAINDER t2, [DV]), and
    `|S|/F_hyb ≤ 122/m²` up to `b̂ = 200` (ROUTEA);
  * the calibration of §3.3 (documented EE corrections ≈ uniform 1.6–3.3 `c₂` at 1000).

  The conclusions below state how much of each cost this model carries and by what factor it may fail.
* k = 3 part (no C3H engine): the same envelope [M]. Its share of the cost is ≤ 1 % on `S_D` and ≤ 29 % anywhere.

### 2.4 Result D at Δ_g = 1000: the converted strip

`analyze_D.py 1000` (575 states; `m/Λ ∈ {1.3, 1.5, 2, 3, 5, 10, 20}`; every J up to `b̂ = 400` or the twist cap; both helicities;
off-diagonal included).

Leverage `L = floor/cost` under E1′, on the states with `J²/Δ ≤ 1` (min over the band; "—" = no such state):

| `m/Λ` \ `b̂` | [0,5) | [5,10) | [10,20) | [20,40) | bulk (ν ≥ 20) share of the cost on [0,5) / [5,10) / [10,20) |
|---|---|---|---|---|---|
| 1.3 | **8.7** | 3.3 | 0.15 | 0.009 | 0.09 / 0.39 / 0.95 |
| 1.5 | **10.3** | 3.6 | 0.18 | 0.012 | 0.07 / 0.41 / 0.95 |
| 2 | **12.5** | 3.8 | 0.28 | 0.019 | 0.15 / 0.49 / 0.94 |
| 3 | **19.3** | 5.7 | 0.48 | 0.048 | 0.13 / 0.49 / 0.94 |
| 5 | **32** | 9.7 | 0.86 | 0.28 | 0.13 / 0.48 / 0.93 |
| 10 | **64** | 20 | 1.9 | — | 0.14 / 0.48 / 0.93 |
| 20 | **129** | 42 | 21 | — | 0.13 / 0.46 / 0.92 |

**Statement (D, Δ_g = 1000) [NE±].** Every state with `m ≥ 1.3Λ` and `b̂ = (2J+1)Λ/m < 5` is positive. In terms of `J`:
J ≤ 2 (1.3Λ), 2 (1.5Λ), 4 (2Λ), 6 (3Λ), 11 (5Λ), 24 (10Λ), 49 (20Λ). Ingredients:
* certified flat floor `N/5`, minus Lemma C [P], minus E1′;
* `L ≥ 8.7`, growing ∝ m;
* the only model input (the bulk-ν envelope at J > 0) may be **too small by ×76** (×135, 76, 141, 239, 467, 980 at 1.5 … 20Λ) on
  top of its ×3 before any state fails.

`L ≥ 3` extends to `b̂ < b̂_A(m)` = 8.5, 8.7, 10.5, 13, 15, 15.3, 13.8. On `5 ≤ b̂ < b̂_A` the bulk model carries 40–95 % of the cost,
so that extension is [NE + M], not [NE±].

Above 20Λ:
* the engine-covered set is `J ≤ √Δ`, i.e. `b̂ ≲ 2Λ/√m`;
* there `G_tol ∝ m²` while E1′ grows ∝ m, so `L` grows ∝ m [D scaling];
* fixed `b̂` with `m → ∞` means `J²/Δ → ∞` (the eikonal regime), which the engine does not cover.

**What this does not give: an `m_*`.** At every mass the states with `b̂ ≳ 10–15` (all J beyond `b̂_A(m)m/2Λ`) fail the sign-free
test by orders of magnitude (§2.2). They remain `H_R` [M].
* Hybrid floor there: [certified] on `b̂ ∈ [91.2, 2Λ]` (ROUTEA), [NE] elsewhere.
* Beyond-Regge evidence: the measured `≤ 0.005–0.02 N_s` (J ≤ 40) and the scalar-skeleton remainder `≤ 122/m²` relative ([NE + M],
  margin ≥ 150).

**At Δ_g = 300** (`analyze_D.py 300`): `L ≥ 2.1` on `S_D` (≥ 4.4 at `m ≥ 3Λ`), with bulk share ≤ 0.37 and tolerable bulk factor
≥ 5.6. `L ≥ 3` only for `m ≥ 2Λ`, `b̂ ≤ 2.5–8.9`. This is weaker than at 1000, as expected (`G_tol ∝ m²`, envelope Δ-independent).

### 2.5 Result Q: large mass

The normal form is `FF = B₂ + B₃ + a(2J(J+1) − c₄ + σ)u³ − 2cu²`, with `N = |B₂| + |B₃| + |FW| + 2cu²` (`q_large_m.py`, flat parts
with exact Wigner functions). Costs:
* smeared k = 2, 3 parts: `g·W_q`, with `W_q` as in §2.2 for Q's weights;
* derivative piece: `|s_d|·W_d`, `W_d = 2au³/N` (RESULT_Q_FWD_CLOSE: the correction enters as `a·2u³s_d`);
* ν-flat factor: `27/m`.

The forward term `a·2J(J+1)u³` keeps `N` from collapsing at large `b̂` (unlike D's `N_s ∝ b̂⁻³`). As a result:

| Δ_g = 1000, `m/Λ` | 1.3 | 1.5 | 2 | 3 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|---|
| min flat FF/N (all `b̂ ≤ 400`) | 0.243 | 0.339 | 0.420 | 0.379 | 0.248 | 0.416 | 0.503 |
| min `G_q` (`c₂`, all `b̂ ≤ 400`) | 1.0·10⁵ | 9.1·10⁴ | 7.7·10⁴ | 5.1·10⁴ | 2.8·10⁴ | 3.0·10⁴ | 2.1·10⁴ |
| min `s_fail′ = 0.18/W_d` | 4.8 | 6.0 | 10.9 | 29 | 110 | 1.9·10³ | 2.0·10⁴ |
| min `L` (floor 9/50), `J²/Δ ≤ 1` | **3.95** (J ≤ 36) | **4.26** (≤ 38) | **5.72** (≤ 44) | **6.15** (≤ 54) | 3.31 (≤ 67) | 1.65 (≤ 96) | 1.87 (≤ 137) |
| min `L` (computed flat FF/N), `J²/Δ ≤ 1` | 5.4 | 8.0 | 14.3 | 14.2 | 5.4 | 5.0 | 10.3 |
| max share of the cost from ν ≥ 12 (bulk) | 0.12 | 0.15 | 0.26 | 0.56 | 0.88 | 0.97 | 0.99 |

* The derivative piece is the only measured, m-dependent quantity, and it is controlled by `u³`. `s_d` was measured m-independent
  from the gap to 20Λ (0.5081 / 0.5095 / 0.5096 / 0.5095 at `m/Λ` = 1 / 3.35 / 10 / 20, EESTEP §5) [NE given `H_RS` = D given RD].
  Against `|s_d| ≤ 0.6`, `s_fail′ ≥ 4.8` at 1.3Λ, ×8. Its share of the cost is ≤ 50 % at 1.3Λ, ≤ 13 % at 3Λ and ≤ 2 % from 5Λ on.
* **Statement (Q, Δ_g = 1000) [NE±].** Every state with `1.3Λ ≤ m ≤ 3Λ` and `J²/Δ ≤ 1` has `L ≥ 3.95` against the certified floor
  9/50. There the bulk model carries ≤ 56 % of the cost. Up to 20Λ, `L ≥ 1.65` (certified floor) or ≥ 5.0 (computed flat floor), but
  the bulk share rises to 0.88–0.99 at `m ≥ 5Λ`. So the conversion is clean for `m ∈ [1.3, 3]Λ` and [NE + M] above.
* For J beyond the engine cover, the sign-free tolerance `G_q ≥ 2.1·10⁴ c₂` holds on every computed state with `m ≥ 1.3Λ`,
  `b̂ ≤ 400`. That is a clean, explicit replacement of HYB's "[M]" by the hypothesis `|Γ| ≤ 2·10⁴ c₂` [M]. Every measured value is
  below 10⁴ `c₂`, so `L ≥ 2`.
* **At 300**: `G_q` is ~10× smaller (`≥ 1.9·10³ c₂` at 20Λ); `L ≥ 1.2` on the engine-covered states for `m ≤ 2Λ`, and `L < 1` at
  `m ≥ 3Λ`, `b̂ ≥ 5`. **Not converted at 300.**

### 2.6 What §2 changes

| hypothesis | before | after |
|---|---|---|
| `H_R` (D, 1000) | [M] on all of `m ≥ 1.3Λ` | **strip `m ≥ 1.3Λ`, `b̂ < 5`: [NE±]** (`L ≥ 8.7`; bulk model needed to ×76); `b̂ < b̂_A(m)`: [NE + M], `L ≥ 3`; rest unchanged [M] |
| `H_mono` at fixed J (corollary of UT §8) | [M] | **J = 0: [D given RD] + [NE]**: `Γ = 21.0m·(ν-flat, multiplicative) + g(ν) + O(1/m)`, with `g` m-independent. J = 28–100: the ν-flat part is J-independent (`20.9 m`); the ν-dependent part ≤ 450 + 100J `c₂` on ν < 20, [NE]. At fixed `(β, ν)` in the window: unchanged [M] |
| HYB (Q, 1000) | [M] everywhere | `m ∈ [1.3, 3]Λ`, `J²/Δ ≤ 1`: **[NE±]**, `L ≥ 3.95`; `m > 3Λ`: [NE + M], `L ≥ 1.65`; elsewhere unchanged |
| curvature of the hybrid | [NE] (quadrature) | `b̂ < 10`: **[P]** (Lemma C); `b̂ ∈ [91.2, 2Λ]`: [certified] (ROUTEA, unchanged) |

---

## 3. The near-gap compact window

### 3.1 The sign-free route there (`analyze_gap.py`)

Strip: Δ_g = 1000, J = 0–40, `m/Λ ∈ {m_gap(J), 1.02, 1.05, 1.1, 1.2}`.

| J | 0 | 2 | 4 | 5 | 6 | 7 | 8 | 10 | 16 | 24 | 32 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| min hybrid `λ/N_s` | 1.000 | 1.000 | 0.204 | **0.2009** | 0.230 | 0.269 | 0.429 | 0.588 | 0.767 | 0.964 | 0.979 | 0.994 |
| `G_tol` (`c₂`) | 2.0·10⁵ | 6.3·10⁴ | 1.1·10⁴ | 6.6·10³ | 3.6·10³ | 1.4·10³ | 343 | 205 | 49 | 15 | 6.7 | 3.7 |
| `L` under E1′ | 7.0 | 6.0 | 2.7 | 1.8 | 1.0 | 0.38 | 0.09 | 0.05 | 0.008 | 0.002 | 6·10⁻⁴ | (J²/Δ > 1) |
| bulk share | 0.05 | 0.19 | 0.63 | 0.76 | 0.86 | 0.95 | 0.99 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 |

The near-gap states are converted to [NE±] only at J = 0 and 2 (margins ≈ 1 there anyway). At the binding J = 5 `+−` state:
* the sign-free leverage is 1.8;
* 76 % of the cost is the bulk (`x = ν/m` up to 1, the light-cone edge);
* the 1/Δ engine does not reach the bulk.

**So the binding state is not error-bounded by this route.** Its evidence stays the chain: the hybrid value 0.20090, with corrections
≤ 10⁻³ and an edge shift ≤ 0.006 [NE + M], `L ≈ 30` against the measured discrepancies.

### 3.2 State-by-state status at Δ_g = 1000 (D)

| states | value | evidence | error bar available? | status after this note |
|---|---|---|---|---|
| J = 0, 2 near gap | ≈ +1.0 | EE_NEARGAP chain | — | **[NE±]** (§3.1, `L ≥ 6`) |
| **J = 5 `+−` twist edge (binding)** | **+0.2009** | hybrid + w2 chain + edge | correction ≤ 10⁻³, edge ≤ 0.006 | [NE + M] (H_M, D9) |
| J = 4, 6–14, 16 | +0.204 … +0.72 | twist-edge grid + edge (CL) | edge ≤ 0.016 | [NE + M] (H_M) |
| **J = 15 twist edge** | ≈ +0.69 (hybrid) | grid from 1.01Λ only | none | **[M]** (D6) |
| even J = 18–24, 28 | positive, ≥ binding | EE_NEARGAP chain | — | [NE + M]; **H_edge [M]** (D4) |
| J = 17, 31 | +0.840, +0.988 | rows + edge | edge computed | [NE + M] |
| odd J = 19–29 | +0.89 … +0.98 | rows (`Δ = 1000 + J`) | — | [NE + M]; **H_edge [M]** |
| J = 32, 36, 40 | +0.904, +0.957, +0.980 | own-mass rows (+ edge at 36, 40); J = 40 validated by a direct sub-row | ≤ 0.05 `N_s` (J = 40) | [NE + M] |
| J = 48 | +0.915 ± 0.08 | direct + sub-row | jackknife ± 0.08 | [NE + M]; **H_edge [M]** |
| **even J = 26, 30, 34, 38, 42, 44, 46** | — | interpolation in J | none | **[M]** (inside `H_LJ(1000)`) |
| **even J ≥ 50, odd J ≥ 33, `m ≤ 1.3Λ` (to J ≈ 344)** | κ-model ≥ +0.90 | model | none | **[D-sketch + M]** (`H_LJ(1000)`) |
| OO block J ≥ 12 | = `++` diagonal | identity | ≤ 0.12 `c₂` | [D given identity] (D5) |
| `m ≥ 1.3Λ`, `b̂ < 5` | ≥ 0.2 | §2.4 | `L ≥ 8.7` | **[NE±]** |
| `m ≥ 1.3Λ`, `5 ≤ b̂ < b̂_A(m)` | ≥ 0.32 | §2.4 | `L ≥ 3`, bulk-model-dominated | [NE + M] |
| `m ≥ 1.3Λ`, `b̂ ≥ b̂_A(m)` | ≥ 0.93 (hybrid) | measured ≤ 0.02 (J ≤ 40); scalar 122/m² | — | **[M]** (`H_R`, large `b̂`) |

**Model-only states at Δ_g = 1000 (D), exactly:**
1. the J = 15 twist edge;
2. even J = 26, 30, 34, 38, 42, 44, 46 near the gap;
3. every near-gap state with even J ≥ 50 or odd J ≥ 33 (`m ≤ 1.3Λ`);
4. the edge window (`x > 0.8`) of even J = 18–32, odd J = 19–29 and J = 48;
5. every state with `m ≥ 1.3Λ` and `b̂ ≥ b̂_A(m)`. On `5 ≤ b̂ < b̂_A(m)` the bulk is model-dominated.

All other near-gap states are computed [NE + M]: computed inside the chain, whose accuracy (`H_M`) is validated where direct data
exist but carries no error bar of its own.

**Q, both Δ_g.**
* The binding pair is error-bounded given RD (Q1).
* The measured channels `++` J = 0–8, OO J = 0, 2, `+−` J = 4–9, 11 are [NE + M given `H_RS`]. Their engine truncation is measured
  and the chain is HYB.
* Model-only:
  * the unmeasured channels (`s_fail ≥ 11`);
  * every state not in a measured channel at `m < 1.3Λ`;
  * at 1000, `m > 3Λ` (bulk-dominated, `L ≥ 1.65`) and `J²/Δ > 1`;
  * at 300, everything at `m ≥ 1.3Λ` except `b̂ < 5`.

### 3.3 Calibration: how the documented corrections act

The documented beyond-Regge corrections (EE_NEARGAP §6), read as the uniform sign-free `|Γ|` that would produce them
(`Γ_eff = |dF/N_s|/W`):

| Δ_g, `m/Λ` | J = 12 | 14 | 16 | 18 | 20 | 24 | 28 | 32 | 40 |
|---|---|---|---|---|---|---|---|---|---|
| 300, 1.3 (`dF/N` = −0.04 … −0.49) | 24 | 23 | 30 | 30 | 32 | 32 | 22 | — | — |
| 300, 1.5 / 2 | — | — | — | — | — | 39 / 52 | 24 / 24 | — | — |
| 1000, 1.3 (≤ 0.005 / 0.02) | — | — | — | — | — | — | — | ≤ 1.6 | ≤ 3.3 |

* The actual corrections are 10²–10⁴ `c₂` locally, yet they act like uniform corrections of only 1.6–50 `c₂`, about 1/100–1/1000 of
  E1′. They cancel against the oscillation of the partial wave.
* At fixed `(J, m/Λ)` the effective size falls by ×10–20 from 300 to 1000 (`p ≈ 2–2.5`), consistent with DG_SCAN's `p_emp`.
* This is what makes the large-`b̂` states positive with large margins, and it is precisely what is not derived.

---

## 4. The large-J tail at fixed Λ: the `H_LJ` wall

At fixed Δ_g, increasing J at fixed `m` drives the state towards the twist edge, with `b̂ → m` and `β ≈ 2 ln(2m/Λ)`.
* **Sign-free:** `G_tol` falls like `b̂⁻³…⁻⁵` (§2.2): below 1 `c₂` at `b̂ ≥ 160`, `m ≤ 3Λ`.
* **Engine:** at the twist edge `J²/Δ ≈ J` grows without bound. The helicity-structure series is asymptotic for `J²/Δ ≳ 2`
  (HLJ_CLOSE §1.2), and the D̂ series is unconverged at ν ≥ 20.
* **Leading-twist / large-spin asymptotics** would need the smooth (non-oscillating) part of the correction at small ν, i.e. `H_win`.
  That is exactly the resummation HLJ_CLOSE could not do.

**Stopped here, as instructed.** `H_LJ(1000)` stays [D-sketch + M] and `H_LJ(300)`'s window stays [I].

---

## 5. Final evidence structure

### 5.1 Result D at Δ_g = 1000

| region of `(J, m)` | grade |
|---|---|
| flat limit (all states) | **[certified]** (`𝔉 ≥ N/5`) |
| `m ≥ 1.3Λ`, `b̂ < 5` (J ≤ 2 … 49 for `m/Λ` = 1.3 … 20; above 20Λ, `J ≤ √Δ`) | **[NE±]**: certified floor + Lemma C [P] + engine [D given RD] with ×3; `L ≥ 8.7`; bulk model tolerated to ×76 |
| near gap, J = 0, 2 | **[NE±]** (`L ≥ 6`) |
| near gap, computed J (4–14, 16, 17, 31, 32, 36, 40, 48, and even 0–28, odd 19–29) | [NE + M] (`H_M`; `H_edge` [M] on the listed J) |
| `m ≥ 1.3Λ`, `5 ≤ b̂ < b̂_A(m)` | [NE + M] (`L ≥ 3`, bulk-model-dominated) |
| `m ≥ 1.3Λ`, `b̂ ≥ b̂_A(m)` | [M] (`H_R`, large `b̂`); hybrid [certified] on `b̂ ∈ [91.2, 2Λ]` |
| J = 15 twist edge; even J = 26–46 without rows; near-gap even J ≥ 50, odd J ≥ 33 | [M] / [D-sketch + M] (`H_LJ(1000)`) |

### 5.2 Result Q

| region | 1000 | 300 |
|---|---|---|
| binding pair `++`/OO J = 0 | [D given RD] + truncation [NE] (×850) | same (×85) |
| measured channels near the gap | [NE + M given `H_RS`] | same |
| `1.3Λ ≤ m ≤ 3Λ`, `J²/Δ ≤ 1` | **[NE±]**, `L ≥ 3.95` | [NE + M] (`L ≥ 1.2`) |
| `m > 3Λ`, `J²/Δ ≤ 1` | [NE + M], `L ≥ 1.65` (≥ 5 with the computed flat floor) | [M] |
| `J²/Δ > 1`, `m ≥ 1.3Λ` | [M], sign-free tolerance ≥ 2·10⁴ `c₂` | [M], ≥ 1.9·10³ `c₂` |
| near gap, unmeasured channels | [NE + M] (`H_fwd^Q`, `s_fail ≥ 11`) | same |

### 5.3 Covering hypotheses after this note (replaces UNIFIED_THEOREM_v2 §5.2 for Δ_g = 1000)

| # | hypothesis | domain now | grade | leverage |
|---|---|---|---|---|
| `H_LJ(1000)` | near-gap large spin | unchanged | [D-sketch + M] | 2.7–4.8 (refined window) |
| `H_R′` | beyond-Regge at `m ≥ 1.3Λ` | **only `b̂ ≥ 5`** (the strip `b̂ < 5` is removed) | [M] (`b̂ < b̂_A`: [NE + M]) | measured ≥ 15; scalar ≥ 150 |
| `H_M` | near-gap chain | computed near-gap states with J ≥ 4 | [M], validated [NE] | ≈ 30 at the binding state |
| `H_edge` | edge window | even J = 18–32, odd 19–29, J = 48 | [M] | ≥ 50 |
| `H_OO` | OO = `++` beyond order 0 | J ≥ 12 | [D given identity] | ≥ 10³ |
| `H_cont` | between masses; J = 15 twist edge | unchanged | [NE]; J = 15 [M] | ≥ 10 |
| `H_bulk` (new name for an old input) | the bulk-ν (`ν ≥ 20`, J > 0) beyond-Regge envelope `≤ 3(450 + 100J) c₂` | used on `S_D` and on Q's `[1.3, 3]Λ` | [M] | tolerated ×76 on `S_D`; ≥ ×1 elsewhere |
| RD | Regge dominance of the crystal integral | everywhere the engine is used | [D-sketch] (EE) | ≫ 10⁶ |
| Q: HYB | as before | **minus `m ∈ [1.3, 3]Λ`, `J²/Δ ≤ 1` at 1000** | [M] | ≥ 300 (scalar) |
| Q: `H_fwd^Q`, `H_smear^Q`, `H_LJ^Q`, `H_RS` | unchanged | unchanged | as before | 19; 8·10³; 1–10³; 85/850 |

Net at Δ_g = 1000: the weakest link is still `H_LJ(1000)` [D-sketch + M]. **No covering hypothesis is removed.** `H_R` loses its
small-`b̂` strip, and HYB (Q) loses `[1.3, 3]Λ`.

### 5.4 Letter wording

> *Per-state positivity.* "In the flat limit both functionals are positive on every heavy state, by a computer-assisted
> (interval-arithmetic) proof: `𝔉 ≥ N/5` for `Ω^EE`, and `≥ 0.18N` for `Ω^Q_up` on all but two corners, which are checked in floating
> point. At `Δ_g = 1000`, for heavy states with `m ≥ 1.3Λ` and impact parameter `b̂ < 5`, we bound every AdS correction explicitly.
> The curvature is bounded by an elementary Mehler–Dirichlet estimate. The beyond-leading-Regge terms are bounded by the analytic
> `1/Δ` expansion, whose coefficients we find to be Δ-independent, with a factor-3 safety margin. Positivity there holds with margin
> ≥ 8.7 times the bound (≥ 4 for `Ω^Q_up` at `1.3Λ ≤ m ≤ 3Λ`). Elsewhere, near the gap, at large impact parameter and at large spin,
> positivity is computed on finite sets of states inside the near-gap chain model and assumed in between. We state these covering
> hypotheses explicitly, with the margins by which each would have to fail (the weakest is the near-gap large-spin window, margin
> 2.7–4.8). As in the flat-space literature, per-state positivity is the one hypothesis of our bounds that is supported numerically
> rather than proved."

Do not write "positivity proved for `m ≥ m_*`": no such `m_*` exists for D (§2.4). Do not write "error-bounded" for the binding state:
it is not (§3.1).

---

## 6. Reproduce

    cd ads4/code/poserr; bash run_all.sh      # ~12 min, <= 2 processes, 3 GB, one thread each

| step | output |
|---|---|
| `curv_lemma.py` | `cache/curv_lemma.txt` (Lemma C check) |
| `eng_orders.py` | `cache/eng_orders.txt` (J = 0 per-order table), `cache/eng_orders_bj.txt` (spin J) |
| `large_m_map.py 1000 …`, `… 300 …`, `… 1000 gap 40` | `cache/large_m_*.json`, `cache/gap_m_1000.json` |
| `q_large_m.py 1000 / 300` | `cache/q_large_m_*.json` |
| `analyze_D.py`, `analyze_Q.py`, `analyze_gap.py` | `cache/analyze_*.txt` (every table of §2–§3) |

The D̂-convergence check of §2.3 (N = 4 vs 5) is a ten-line script over the same `bj` caches. It was run from the scratchpad and is
not part of `run_all.sh`; its numbers are quoted in §2.3.

## 7. Caveats

1. Only grid masses `m/Λ ∈ {1.3, 1.5, 2, 3, 5, 10, 20}` were scanned. Between them, `W` and the floor are explicit smooth functions and
   `L` rises monotonically along the grid, but nothing is certified between grid points [NE].
2. The flat floor `N/5` is certified. The curvature is bounded rigorously only for `b̂ < 10` (Lemma C). Elsewhere the computed
   curvature is used [NE].
3. The ×3 safety factor on the engine terms is a choice, not a derived bound. At J = 0 the engine's own truncation (the next order,
   ≤ 14·300/m `c₂`) is below 0.3 % of the envelope at `m ≥ 1.3Λ`. At J > 0 only orders 1–2 exist, and at J = 28 order 2 is
   ≤ 150 `c₂` on ν < 12 but larger than order 1 on 12 ≤ ν < 20 (1.8·10³ `c₂`): the 12–20 window is not convergent there.
4. E1′ uses the worst entry of the conserved 2×2 basis, not the helicity basis, and treats `++`/`+−` alike. This is conservative.
5. The k = 3 (C3H) part uses the k = 2 envelope [M]. Its share of the cost is ≤ 1 % on `S_D`.
6. Q's analysis uses flat `N` and `FF`. AdS curvature of the forward pieces is not included: it is positive for the value piece
   (`P_{−½}(cosh β) < 1`), and Stage 6's HYB minimum is 0.197.

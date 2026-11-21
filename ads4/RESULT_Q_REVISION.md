# Result Q, revision round 1: the smeared graviton corrections, the swapped form of V₄, and the W³W³ residue

*21 Nov 2026, revision agent. Unrefereed, not committed. No existing document or code file was edited; no `bj_*` file was touched or
run. New code: `ads4/code/resultq/rv_*.py`; outputs in `resultq/cache/rv_*` (git-ignored). Read-only imports: `stage6_model`,
`stage6_stress`, `stage6_scan`, `common_q`, `stage4_hj`, `spinfunc/{l9_neargap, l7_dc2_analysis, bf_qan, bg_eeq, u3_chi_EE}`,
`h7i14/v4_chi`, `omegau/nf_common`. Read-only data: Result D's chain grids `spinfunc/cache/w2_recs_final_xc0.pkl` (EE_NEARGAP run A),
BELOW_FOLD2's `bg_eeq_*` EE images, the CHMRS LaTeX source in the scratchpad. Compute: about 15 CPU-minutes, at most 2 cores and under 1 GB.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and numerically verified; **[NE]** numerical evidence;
**[M]** model; **[I]** open. `c₂` = `m² ×` relative correction. Functional: δ = 0.2 (`Q3-full-delta0.2`, `A = 17.241`, `a = 0.1694`,
`c = ½`, `κ = 0.1099`). Gauge for the forward `s_d`: H0 (SB in brackets where it matters).

---

## 0. Answer in brief

1. **H_smear^Q (review §3.1): computed at the binding states. The graviton-specific corrections of the smeared parts are 10⁻⁵–10⁻⁴ of
   `N`, four orders of magnitude below the tolerance.** [NE + M]
   * **Source.** Result D's own EE `K₂` crystal rows (`Γ_{C2}`, `Γ_{C3H}` per helicity, `x = ν/m ∈ [0.05, 0.8]`, rows `Δ = 200…1600`)
     are projected on **Q's** `ψ₂`, `ψ₃` with Result D's chain formula. Below the fold, the closed-form engine's EE `k = 2` images are
     projected with the exact (min-norm) Abel weight of Q's `ψ₂`.
   * **Validation.** The same code with NF02's `ψ`s reproduces Result D's stored `dF` on all 30 hybrid test states to `6·10⁻¹⁵` [DV].
   * At `++` `J = 0`, at the gap: `δB₂/B₂ = +3.1·10⁻⁵` (Δ_g = 300) and `+1.2·10⁻⁶` (1000). The tolerance is ≈ 24 %.
   * **Corrected minimum:** `FF/N = +0.1427 ± 0.0017` (SB `+0.1410`) at both Δ_g, unchanged in the fourth digit (shift `+2·10⁻⁵`).
   * **Not covered:**
     * the edge window `x > 0.8` (held there, as in Result D);
     * the `k = 3` part below the fold (no `o_H·q` engine was built).

     Both have tolerances of `≳ 3·10⁵ c₂` at `++` and `≳ 760 c₂` at the `+−` gap states. The measured `+−` edge values are `≤ 0.2 c₂`.
2. **V₄ in swapped form (review §2.1): partly settled.**
   * **The infinite `[TT]_{0,J}` J-sum is never needed** [D]. Swap only the heavy blocks (margin ≈ Δ_g at the edge). The light side is
     then `Ψ₄` acting on the function `𝒢_rest = 𝒢 − Σ_heavy`.
   * **In the scalar Mellin form the leading family enters J-resummed.** It appears as `M` and `∂M` at the double-twist subtraction
     points. For a spin-2 Regge-channel (u-channel) exchange these individually nonzero `O(G)` light-family terms **cancel identically**
     against the crossed double-twist terms, so `B̂_{4,T} = 0` exactly [DV, sympy].
   * **Not done:** an evaluation of `Ψ₄`/`Ψ₅` on MFT or Einstein `⟨TTTT⟩` in spinning swapped form. The spinning
     position-space ↔ Mellin equivalence remains [I].
   * **What is left.** The open part reduces to one place: integrability at the corner where the light-cone edge `u′ → 0` meets the
     large-`η` Regge corner (`r → 0`, `ρ = 2rη → 1`, with `R ≈ 4ρ²/(1−ρ²)²`). There the `[TT]₀` tail builds the u-channel singularity.
     A surviving residue would still have leverage `10¹⁰–10²²`. **V₄ stays [D-sketch], now with a precise statement (V₄′, §2.4).**
3. **W³W³ (review §2.3): no decaying residue in `∂_{ν²}Ĉ₄(0)` in the scalar-Mellin form** [DV algebra; spinning D-sketch].
   * **u-channel W³W³ diagram:** it is spin 2 in the Regge channel, so it is annihilated by the `k = 4` functional (`B̂₄ ≡ 0`). Its
     `1/(ν² + 9/4)` residues live only in `Ĉ′_EE` (`k = 2`), which is Result D's H17.
   * **s- and t-channel diagrams:** they act like contacts. `B̂_{4,T}` is exactly a polynomial in `T` for every pole position, hence a
     polynomial in `ν²`. The flat-unit argument then gives `−Λ²∂_{ν²}Ĉ₄(0)|_{W³W³} = 4πG|ĝ₃|²` exactly.
   * **The lower slope with the measured `s_d`:** `0.491 / 0.475` (H0/SB; binding state `++` `J = 0` at the gap). If `s_d` is double
     counted (review §3.3, reading (a)) it is `0.255`. Both are ≥ ¼.
   * **Implied closure:** `|ĝ₃|²Λ⁸ ≤ 120–200` (Δ_g = 300) and `154–257` (1000).
   * **Review issue 3:** there is no contradiction with Result D's `1.07·10⁴ / 1.43·10⁴`. Q's closure is a *stronger statement under more
     hypotheses*. Result D's headline is **not** superseded (§3.4).

---

## 1. H_smear^Q: graviton corrections of `∫ψ₂Ĉ′_EE` and `∫ψ₃Ĉ₃` in Ω^Q

### 1.1 Method (`rv_smear.py`)

**Approach.** Ω^Q's smeared parts are Result D's functionals `Ĉ′_EE` (k = 2, EE structure, min-norm de-shift) and `Ĉ₃` (k = 3, `o_H·q`,
no de-shift), with different weights `ψ₂`, `ψ₃`. So Result D's measured beyond-leading-Regge EE corrections apply as they are. Only
the projection weights change.

**Correction of a hybrid state.** This is `q3_analysis.neargap` with the scalar remainder `ρ = 0`, i.e. exactly the HYB base of Q:
* `dF_hh = Σ w ψ₂ env(h,x,2) Γ_{C2,h}(x,m)/m²  +  Σ w ψ₃ u env(h,x,3) Γ_{C3H,h}(x,m)/m²`;
* the off-diagonal term (`J ≥ 4`), as in `q3_analysis`;
* `x = p√u`;
* the `Γ` model (`l9_neargap.model`): cubic spline in `x`, linear in `1/m` between `Δ` rows, held outside.

`env(++) = 2 − y` and `env(+−) = (2 − y)/(1 − y)⁴` match the normalisation of `stage6_model.hyb`. The `B₂` and `B₃` of Q and of NF02
are the same objects (`∫ψ₂(2−y)D_h`, `σu∫ψ₃D_h`).

**Validation [DV].**
* The NF02 `ψ`s through `rv_smear` reproduce Result D's stored `dF` (run `final_xc0`): 30 hybrid states (`J = 0, 2, 4, 5`;
  Δ_g = 300, 1000), diagonal and off-diagonal, max relative deviation `5.7·10⁻¹⁵`.

**Variants** (the parts no row covers):
* **A:** Result D's convention (spline below the first node, held above `x = 0.8`);
* **B:** held below;
* **Z:** zero outside the rows;
* **E:** relative `+−` correction held above 0.8;
* **X:** edge window × 10.

### 1.2 Results at the binding states (`cache/rv_smear.txt`)

`FF/N` columns: HYB (`K = 0`) → with the measured forward `s_d` (BELOW_FOLD2/3) → with `δB₂ + δB₃` (variant A).

| Δ_g | state (`m/Λ`) | `B₂` | `B₃` | `FF/N`: HYB → `+s_d` → `+δB` | `δB₂/B₂` | `δB₃/|B₃|` | `(δB₂+δB₃)/N`: A / Z / X |
|---|---|---|---|---|---|---|---|
| 300 | `++` J=0, 1.000 | 1.805 | −0.203 | 0.2000 → **0.1427 → 0.1427** | `+3.1·10⁻⁵` | `+5.9·10⁻⁵` | `+2.3 / +2.4 / +1.0 ·10⁻⁵` |
| 300 | `++` J=0, 1.05 | 1.826 | −0.184 | 0.2891 → 0.2436 → 0.2437 | `+2.7·10⁻⁵` | `+5.3·10⁻⁵` | `+2.1 / +2.2 / +1.4 ·10⁻⁵` |
| 300 | `++` J=2, gap 1.0067 | 0.388 | −0.421 | 0.2534 → 0.2090 → 0.2089 | `−5.8·10⁻⁵` | `−7.0·10⁻⁴` | `−8.5 / −3.3 / −53 ·10⁻⁵` |
| 300 | `++` J=2, 1.118 | 0.534 | −0.410 | 0.2002 → 0.1665 → 0.1665 | `+5.4·10⁻⁵` | `−4.0·10⁻⁴` | `−5.1 / −2.7 / −26 ·10⁻⁵` |
| 300 | `+−` J=4, gap 1.0133 | 1.811 | +0.198 | 0.2283 → 0.2777 → 0.2777 | `+6.1·10⁻⁵` | `−9.0·10⁻⁵` | `+2.8 / +1.0 / +20 ·10⁻⁵` |
| 300 | `+−` J=5 (odd), gap 1.0166 | −1.776 | +1.869 | 0.2612 → 0.2612 → 0.2611 | `+2.2·10⁻⁴` | `−3.4·10⁻⁴` | `−14 / +2.0 / −155 ·10⁻⁵` |
| 300 | `+−` J=5 (odd), 1.15 | −1.082 | +1.175 | 0.2025 → 0.2025 → 0.2025 | `−1.3·10⁻⁴` | `−1.3·10⁻⁵` | `+3.1 / +3.6 / −0.4 ·10⁻⁵` |
| 1000 | `++` J=0, 1.000 | 1.805 | −0.203 | 0.2000 → **0.1428 → 0.1428** | `+1.2·10⁻⁶` | `+5.8·10⁻⁶` | `+1.1 / +1.3 / −0.5 ·10⁻⁶` |
| 1000 | `++` J=2, gap 1.0020 | 0.382 | −0.419 | 0.2577 → 0.2130 → 0.2130 | `−1.2·10⁻⁵` | `−2.9·10⁻⁵` | `−4.4 / −2.2 / −24 ·10⁻⁶` |
| 1000 | `+−` J=4, gap 1.0040 | 1.807 | +0.202 | 0.2086 → 0.2594 → 0.2594 | `+5.6·10⁻⁶` | `+2.0·10⁻⁴` | `+15 / −0.1 / +152 ·10⁻⁶` |
| 1000 | `+−` J=5 (odd), gap 1.0050 | −1.845 | +1.952 | 0.2696 → 0.2696 → 0.2695 | `−1.0·10⁻⁴` | `−6.6·10⁻⁴` | `−14 / +0.6 / −149 ·10⁻⁵` |

Also computed (in the file): `++` J=4 at the gap, `+−` J=4 at 1.05Λ, `+−` J=5 at 1.2Λ. All have `|δ|/N ≤ 4·10⁻⁵`. Off-diagonal mixing:
`|off|/N ≤ 2.2·10⁻⁵`.

**Reading.**
* At the binding state (`++` J = 0, gap) the relative correction of `B₂` is **`3·10⁻⁵` (300) and `10⁻⁶` (1000)**. The review's
  tolerance is 24 % (`ε ≲ 0.24`), so the margin is about four orders of magnitude.
* The correction is positive there, i.e. it helps.
* The largest |shift| of `FF/N` at any listed binding state is `1.5·10⁻³` (variant X at `+−` J = 5, gap).
* **The full-set minimum is unchanged: `+0.1427 ± 0.0017` (H0), `+0.1410` (SB), Δ_g = 300 and 1000** [NE + M; the smeared part is
  now measured at these states].

### 1.3 Below the fold: the engine's `k = 2` images (`rv_k2eng.py`, `cache/rv_k2eng.log`)

**Chain.** This is the k = 4 `bg_qan` chain with Stage 4's exact EE `k = 2` kernel (`EE1`):
* `Γ_exact = m²(c T_ref I_st/I_ref − T_st)`;
* HYB's own term `γ₂ Σ_j R_j(ν − ij)²`;
* the `a_h` class removed.

**`γ₂` is a prediction, not a fit choice [DV].** The flat-exact/LR_flat ratio with the `(2 − y)` envelope predicts `γ₂ = −¼`. The
engine gives:

| state | fitted `γ₂` |
|---|---|
| J = 0, Δ_g = 300 | −0.2500 |
| J = 0, 1.05Λ | −0.2500 |
| J = 0, Δ_g = 1000 | −0.2500 |
| J = 2 | −0.2509 / −0.2501 |

This is a further check that the engine reproduces the HYB finite-`x` structure for `k = 2`.

**Per-harmonic correction `Re Γ_K/N̂`.**
* J = 0, Δ_g = 300: `|Γ| ≤ 0.015 c₂` on `ν < 15`, rising to −1.8 at `ν = 40`. There the engine's a⁴ step is −4.6, i.e. it is
  unconverged at `x ≈ 0.13`.
* Δ_g = 1000: `≤ 0.19 c₂` on all of `ν < 41`.
* The K₂ chain's spline gives −0.66 … +2.0 c₂ on the same window. The two agree at the `O(1) c₂` level, within the engine's truncation.
* J = 2 at 300: the engine is unconverged above `ν ≈ 20` (a⁴ step up to 40 c₂), as BELOW_FOLD3 §2 found.

**Window contribution** (`ν < 41`, exact Abel weight `X^Q`):
* `δB₂ = −1.0·10⁻⁶` (engine) vs `+2.2·10⁻⁶` (chain, local `ψ₂`) at `++` J = 0, gap, Δ_g = 300;
* `−1.7·10⁻⁵` vs `+2.0·10⁻⁵` at J = 2;
* `≤ 10⁻⁸` at 1000.

All are `≤ 10⁻⁵` of `B₂`. **The fold region is irrelevant for Ω^Q's smeared part**, because `ψ₂ ∝ p` there [NE given D ≡ crystal for EE].

Abel solves for Q's `ψ₂`: LSQR residual `10⁻¹⁵`; the graviton-pole identity (V3) holds to `10⁻¹⁴` [DV].

### 1.4 Uncomputed windows: tolerances (`rv_tol.py`, `cache/rv_tol.txt`)

The table gives the uniform `Γ` (in `c₂`) that would drive the state to `FF = 0`.

| state (Δ_g = 300) | k=2, `x > 0.8` | k=3, `x > 0.8` | k=3, `ν < 15` (no engine) |
|---|---|---|---|
| `++` J = 0, gap | `−7.9·10⁵` | `−3.1·10⁵` | `+2.0·10⁶` |
| `++` J = 2, gap | `−1.6·10⁶` | `−6.1·10⁵` | `+3.6·10⁶` |
| `+−` J = 4, gap | `−3.5·10³` | **`−7.6·10²`** | `+4.4·10⁶` |
| `+−` J = 5, gap | `−9.1·10³` | `−2.0·10³` | `+9.4·10⁶` |

At Δ_g = 1000 every entry is larger (`+−` J = 4 gap, k = 3 edge: `−3.0·10³`).

**Why the `+−` gap states are the exception.** They are the only states with a tolerance below `10⁴ c₂`: `env₊₋ ∝ (1 − y)⁻⁴` and
`y → u_gap ≈ 0.97` at `p → 1`. The measured `Γ_{C2,+−}` and `Γ_{C3H,+−}` at `x = 0.8` are `0.2 c₂` and still falling.

Result D's H_edge concern (`−0.2 … −0.4 N_s` at J ≈ 30–40) is about large-J states, not these.

### 1.5 Grade and what remains

| claim | grade |
|---|---|
| `rv_smear` = Result D's chain on NF02 (`6·10⁻¹⁵`) | [DV] |
| H_smear^Q at the listed binding states (`++` J = 0, 2; `+−` J = 4, 5; Δ_g = 300, 1000): `|δB|/N ≤ 1.5·10⁻⁴` (A), `≤ 1.6·10⁻³` (X) | [NE + M] (= Result D's H_M: K₂ chain + engine) |
| min `FF/N` unchanged, `+0.1427 ± 0.0017` | [NE + M] |
| `γ₂ = −¼` for `k = 2` (engine reproduces HYB's finite-`x` term) | [DV] |
| below-fold `k = 2` window `≤ 10⁻⁵ B₂` | [NE given D ≡ crystal] |
| `k = 3` below the fold; edge window `x > 0.8` | [I], tolerance-bounded (§1.4) |
| Δ_g = 300, J = 0, 2 states below the lowest row (`m = 298 < 398`): `Γ` held in `c₂` units (`m²`-scaling) | [M] |

**Remaining unmeasured.** These are covered only through HYB plus Result D's model hypotheses:
* states between the listed masses (`H_cont`);
* `J ≥ 6` near the gap (Result D's grids exist to J = 28/32, but Q's binding states are at J ≤ 5);
* the OO block (`H_OO^Q` unchanged).

**Conclusion.** H_smear^Q, as the review posed it (graviton corrections of the smeared parts at the binding states vs the 24 %
margin), is answered with a margin of about 10⁴.

---

## 2. V₄ in swapped form and the `[TT]_{0,J}` tail

### 2.1 What the review asks, and where the J-sum enters

**The review's worry.** The swapped (s-channel) form of `Ψ₄` contains the leading family `[TT]_{0,J}` with a single zero, for every
`J`. At `O(G)`:
* the family has `γ_{0,J} ~ 1/J` (T exchange, `τ_T = 1`);
* per-block Gel'fand–Shilov continuation followed by the J-sum need not equal the continuation of the sum.

**This is correct as a statement about the J-sum.** CHMRS say the same in their own setting (sec. light, below `eq:abusenotation`):
"`Π_{k,η}|_light` is *not* equal to a sum over light operators of OPE coefficients times `Π[G]`, since the expansion in Regge moments
does not necessarily converge on light operators."

**Their fix.** CHMRS (App. `omegalight`) never sum light blocks. They define `ω|_light := ∫dDisc[ℋ·𝒢^EFT]`, ω acting on the EFT
correlator as a *function*, and evaluate it by an arc at infinity.

### 2.2 Reorganisation: only heavy blocks are swapped [D]

`Ψ₄[𝒢] = 0` holds in commutator form. Split `𝒢 = Σ_{τ ≥ Δ_g} p G + 𝒢_rest`. Then:
* the heavy swap `Ψ₄[Σ_heavy] = Σ_heavy p Ψ₄[G]` needs only K1′/K2/K3. At the edge `u′ → 0`, heavy blocks are integrable with margin
  ≈ Δ_g (Stage 3 §1.4);
* the light side is `Ψ₄[𝒢_rest]`, with `𝒢_rest = 1 + T + Σ[TT]` **as a function**.

At tree level, `𝒢_rest^{(1)}` is the Einstein + `E×W³` + `W³W³` + contact correlator, plus the s-channel double-trace parts of heavy
exchange diagrams. Lemma A acts on the Einstein and `E×W³` correlators as functions. **No J-sum of the `[TT]₀` family is ever
performed** [D].

What remains is that `Ψ₄[𝒢_rest]` is an absolutely convergent integral, and that its swap with the commutator is legitimate. Near
`u′ → 0` the s-channel light-cone expansion of `𝒢_rest` converges for `v′` in the integration domain:
* in the CHMRS frame `v′ ≥ ¼`, so the crossed pole `1/v′` of `f` is never reached;
* on the edge `u′ = 0`, `v′ = (1+r²)²/(4r²) ∈ [1, ∞)`, i.e. `z̄ ∈ (−∞, 0]` and `|ρ̄| < 1`.

The J-tail therefore converges **pointwise along the edge** and resums into the light-cone coefficient functions. The only place
where it does not converge uniformly is `ρ̄ → −1`, i.e. `r → 0` along the edge. That is the corner where the light-cone edge meets
the large-`η` Regge corner (`ρ = 2rη → 1`, `R ≈ 4ρ²/(1−ρ²)²`). There the `[TT]₀` tail builds the u-channel light-cone singularity.
(The review's phrase "the crossed-channel light-cone singularity, which is exactly where the poles of `f` sit" applies here to the
`u′`-pole at this corner. The `v′`-pole is outside the domain.)

### 2.3 Scalar Mellin form: the light-family terms, J-resummed, cancel identically (`rv_w3w3.py` (2)) [DV, algebra]

**Setup.** In CHMRS's Mellin form (`B̂_{k,T}` at fixed `U′ = 2Δ_φ − T`, the u-channel Regge variable), the s-channel `[φφ]_{0,1}`
families appear only through the kernel poles at `S′ = 2Δ_φ, 2Δ_φ + 2`. Their residues are the Mellin amplitude and its `S′`
derivative at the double-twist points, i.e. the whole J-tail resummed.

**Test.** Take a spin-2 Regge-channel pole term `M = S′²P(U′)/(U′ − c)`, the Einstein-type `O(G)` growth, Δ_φ = 3, `k = 4`. The four
finite residues are individually nonzero and `O(G)`:

`S′ = 6: −48 P(T)(5T² − 58T + 156)/((T−8)²(T−6)²(T+c−6))`, `S′ = 8: 256 P(T)/((T−10)(T−8)(T+c−6))`,
`S′ = T: 8T²P(T)/((T−8)(T−6)²(T+c−6))`, `S′ = T−2: −8(T−2)²P(T)/((T−10)(T−8)²(T+c−6))`

**Result.** They sum to **exactly 0** for every `T` and every pole position `c`. The same term at `k = 2` gives the pole
`4P(T)/(T+c−6)`, which is the origin of `16πG/(ν²+9/4)`.

**Reading.** In the scalar proxy the swapped light side, "identity + T + `[TT]₀` edge terms", vanishes identically when the
`[TT]₀` tail is J-resummed by the Mellin amplitude: **no `O(G)` remnant**. This is the content of CHMRS's "graviton exchange
contributes only to `k = 2`" [Q], made explicit term by term [DV].

### 2.4 Verdict and the precise residual hypothesis

**Not done.**
* No spinning MFT or Einstein `⟨TTTT⟩` swapped-form evaluation of `Ψ₄`/`Ψ₅`. No engine evaluates the spinning `EE`/`f⁺₄` kernel off the
  heavy-block saddle (Stage 3 §1.6), and building one exceeded this round.
* **No numerical bound on the large-J tail is given.** The scalar Mellin cancellation is exact, but it is a proxy.

**V₄′ (replaces "V₄ with the `[TT]₀` tail").**
* **Statement.** On tree-level `⟨TTTT⟩` correlators with `J₀ < 4`, the CLS-D integral `Ψ₄[f⁺₄ EE · 𝒢_rest]`, `Ψ₅[f⁻₅ o_Hq · 𝒢_rest]`
  converges absolutely at the corner `u′ → 0`, `r → 0`, `ρ → 1`, and equals the spinning analogue of CHMRS's residue-at-infinity
  (Mellin / eftformula) form.
* **Grade.** [D-sketch] (scalar: [Q, CHMRS App. omegalight] + [DV] cancellation).

**Leverage if it failed.** Unchanged from Stage 12 (`10¹⁰–10²²` per unit relative `O(G)` residue). V₄′ stays the single load-bearing
analytic input of the upper bound.

**Cheapest decisive test (not done).** The edge integral of Stage 3 §1.6 for MFT (`W_t = ⟨41⟩⟨32⟩` on `u′ = 0`), evaluated with the
spinning kernel's edge Taylor coefficients near `r → 0`. Alternatively, the scalar position-space `C_{4,ν}` on a d = 3 GFF + Witten
exchange, as a check of the corner power counting.

---

## 3. The W³W³ residue in `∂_{ν²}Ĉ₄(0)` and the `|ĝ₃|²` sector

### 3.1 Channel bookkeeping [D]

The `|ĝ₃|²` correlator is the graviton exchange between two W³ vertices in the s, t and u channels.

**u-channel diagram (the Regge channel of all Ĉ_k).**
* It grows like spin 2: the large-`s` dependence comes only from contracting the exchanged graviton's two index pairs. The W³ vertex
  factors depend on the momenta and polarisations at each vertex, and `ε_⊥·p_other = 0` in the Regge frame [D-sketch for the
  light-ray-polarised CFT version].
* So by Lemma A / CHMRS's eftformula (`r^{−1}` only) it contributes **nothing** to `Ĉ₄`, `Ĉ₅`.
* Its `1/(ν² + 9/4)` residues sit in `Ĉ′_EE` (`k = 2`). That is Result D's H17, `O(Λ⁻²)` in `κ`.

**s- and t-channel diagrams.** These carry the flat `4πG|ĝ₃|²` and are "contact-like" in the u-channel Regge limit (ADS4_LIGHT §2.4).

### 3.2 Polynomiality of the s-/t-channel light side (`rv_w3w3.py` (1)) [DV algebra, scalar]

* For `M = S′^p U′^q/(T′ − c)`, where the pole runs with `S′` at fixed `U′`, `B̂_{4,T}` is the contour at `|S′| = ∞`.
* It equals the sum of the finite residues, and **is a polynomial in `T` for every `c`**. Checked for `(p,q) = (5,0), (4,1), (3,2),
  (6,0)`; e.g. `B̂ = 16` and `B̂ = 48T − 16c + 384`.
* The `T`-poles of the individual residues cancel.
* With CHMRS's Barnes integral, a polynomial `B̂` gives a polynomial in `ν²`.

**Consequence.** In AdS units, `−Ĉ₄(ν)|_{W³W³} = G|ĝ₃|²(c₁ν² + c₀ + …)`. The flat-unit argument (Stage 12 §2.2) kills every `ν^{2j}`
term with `j ≥ 2` and fixes `c₁ = 4π`.

**Result:** **`−Λ²∂_{ν²}Ĉ₄(0)|_{W³W³} = 4πG|ĝ₃|²` exactly; there is no `c/(ν²+9/4)ⁿ` term.** The curvature constants (`c`, satellite
sums) enter only `c₀ = b₀`, the value piece, which is `Λ⁻²`-suppressed (Stage 12 leverage row (v)).

**Grade.** [D-sketch] for `⟨TTTT⟩`. The scalar Mellin algebra is [DV]. The spinning step assumes a CLS-D analogue of CHMRS's eftformula
(light side = Regge Laurent coefficient of `𝒢^EFT`), the same input as V₄′.

**The review's mechanism.** "Decaying residues `O(1)·G|ĝ₃|²` in the derivative piece" needs a Regge-channel pole in a diagram with
Regge spin ≥ 4. No W³W³ diagram has both.

### 3.3 `κ` and the lower slope, re-evaluated with the measured `s_d` (`rv_lower.py`, `cache/rv_lower.txt`)

**`κ`.** `κ = 0.1099` (flat moments of the δ = 0.2 functional), plus H17-type `O(Λ⁻²)` from the k = 2, 3 parts, plus `b₀` (`Λ⁻²`).
* The slope term `4a` in `k` is now exact (§3.2).
* Grade: [D-sketch + M], up from [I].

**Lower bound `g₄Λ⁶/(8πG) ≥ |ĝ₃|²Λ⁸/(4c_L)`.** Here `c_L = max (∂B₄ − B₅)/B₄` over Stage 6's state set, with
`∂B₄ → ∂B₄ + 2u³s_d`.

| Δ_g | `K = 0` (Stage 3) | measured channels | universal `σ_h s₀`, odd-J `+−` worse sign (H0 / SB) |
|---|---|---|---|
| 300 | `c_L = 0.981` (`+−` J=4 gap) → slope **0.255** | `c_L = 0.509` (`++` J=0 gap) → **0.491** | 0.491 / **0.475** |
| 1000 | `0.993` → 0.252 | `0.508` → 0.492 | 0.491 / 0.475 |

**Why the binding state moves.** The measured `s_d` makes the `+−` J = 4 gap state gain (`c_L` there drops to 0.47). It moves the
binding state to `++` J = 0, where the flat `∂B₄ − B₅ = (1 − 0 − 1)u³ = 0`. There `c_L = s_d` exactly.

**So the lower slope is 0.475–0.49 if the `σ_hU` correction is real.** It is 0.255 if it is double counted (review §3.3 reading (a)).
**¼ survives in both readings** [M + NE], given §3.2.

### 3.4 Review issue 3: the implied `|ĝ₃|` closure vs Result D's headline

Q's two-sided statement gives `(κ + slope)|ĝ₃|²Λ⁸ ≤ 2L_AdS`:

| Δ_g | `2L_AdS` | slope 0.255 (`K = 0`) | slope 0.475–0.49 (measured) | Result D (`16L_A/k(0)`) |
|---|---|---|---|---|
| 300 | 72.04 | **198** | **120–123** | `1.07·10⁴` |
| 1000 | 92.81 | **257** | **154–159** | `1.43·10⁴` |

If the W³W³ slope were halved or zero (the review's worst case), the closure would be 304 / 656 (300). §3.2 excludes that in the
Mellin form.

**Which is right: both, under different hypotheses.**
* **Not a contradiction.** Adding the forward `k = 4, 5` sum rules to the smeared `k = 2, 3` ones, and eliminating `g₄` between the
  upper and the lower functional, *must* give a stronger `|ĝ₃|` bound than the `k = 2, 3` functional alone. The same happens in flat
  space: Stage 0's `x² ≤ 110` vs NF02's weaker closure; Result D is about 10× weaker than the transplanted CHLPS optimum.
* **Q's closure is valid only given everything Q needs.** That is:
  * V₄′;
  * (F) + (P) spinning, now including §3.2's spinning eftformula;
  * the de-shift existence (C4/C5);
  * H_fwd^Q, H_OO^Q (the lower functional's binding state is `++` J = 0, so H_OO^Q at J = 0 enters here too);
  * the `σ_h U` double-counting question, which moves the closure by ×1.6;
  * the realised-weight limits.

  **Grade:** [D-sketch + NE + M].
* **Result D's headline (`1.07·10⁴ / 1.43·10⁴`) needs none of the forward-functional inputs.** It stands as the bound under D's
  hypotheses and is **not superseded**.
* **Recommendation.** Record the Q-closure as a *conditional strengthening* (`|ĝ₃|²Λ⁸ ≲ 120–260` given Q's hypothesis list). Do not
  quote it as Result D's replacement until V₄′ is settled.

---

## 4. Hypothesis list of RESULT_Q_STATEMENT, deltas

| # | item | was | now |
|---|---|---|---|
| 3 | **V₄** incl. the `[TT]₀` large-J tail | D-sketch, "finitely many light-family terms" false | **V₄′** (§2.4): the J-sum is not needed; the residual is corner integrability (`u′ → 0` ∩ large-`η`) plus spinning eftformula equivalence. Scalar Mellin: light-family terms cancel identically [DV]. Still [D-sketch]; leverage unchanged |
| 7 | **H_smear^Q** | [I] (never computed) | **[NE + M]** at the binding states: `|δB|/N ≤ 1.5·10⁻⁴` against 0.143 (factor ≈ 10³–10⁴). Residual [I]: edge window `x > 0.8` and `k = 3` below the fold, with tolerances `≥ 3·10⁵ c₂` (`++`) and `≥ 7.6·10²` (`+−` gap, 300) |
| — | W³W³ residue in `∂Ĉ₄(0)` (review item 7) | [I] | **zero** (u-channel annihilated; s-/t-channel polynomial): [D-sketch], scalar algebra [DV] |
| — | `κ` | [I] | [D-sketch + M] (`0.1099` + H17 `O(Λ⁻²)` + `b₀Λ⁻²`) |
| — | lower slope | [I] / ¼ [M] | 0.255 (`K = 0` or double counted) to 0.49 (measured `s_d`); ≥ ¼ in both [M + NE] given §3.2 |
| — | `|ĝ₃|` closure | unremarked | `|ĝ₃|²Λ⁸ ≤ 120–260` [conditional, D-sketch + NE + M]; Result D's headline not superseded |
| new | **spinning eftformula** (the CLS-D light side = Regge Laurent coefficient of `𝒢^EFT`; contains (P) for exchanges in non-Regge channels) | — | [D-sketch]; scalar [Q] |
| 1, 2, 4–6, 8–10 | — | — | unchanged (the Δ_g = 100 column stays dropped) |

## 5. Files (`ads4/code/resultq/`)

| file | content | output | time |
|---|---|---|---|
| `rv_smear.py` | Result D's K₂ chain on Q's `ψ₂`, `ψ₃` at the binding states; validation on NF02; variants A/B/Z/E/X | `cache/rv_smear.{txt,pkl}` | 1 min |
| `rv_k2eng.py` | engine EE `k = 2` images below the fold, `γ₂` fit, exact Abel weight of Q's `ψ₂`, window projections | `cache/rv_k2eng.{log,txt}`, `rv_k2eng_{T,Y}_*.pkl` | 5 min |
| `rv_tol.py` | uniform-`Γ` tolerances of the uncovered windows | `cache/rv_tol.txt` | 10 s |
| `rv_lower.py` | lower-bound functional `c_L` with measured `s_d`; implied closure | `cache/rv_lower.txt` | 1 min |
| `rv_w3w3.py` | scalar-Mellin `B̂_{4,T}` of t-channel (polynomial) and u-channel (zero; residue breakdown) pole terms | `cache/rv_w3w3.txt` | 10 s |

Order: `rv_smear.py` → `rv_tol.py`; `rv_k2eng.py` (needs `rv_smear`); `rv_lower.py`; `rv_w3w3.py`.

# H14 (object identity near the gap) for scalar-type heavy blocks

Status: 18 Nov 2026, numerics agent. **Unrefereed. Not committed.** Code: `ads4/code/objid/` (`o1_flat.py`, `o2_ads.py`, `o3_envelope.py`;
outputs in `objid/cache/`, git-ignored). It imports `omegau/`, `exactskel/` (e1/e3/e7 caches, `evenkernel`) and `robust/` **read-only**.

Sources:
* CHMRS LaTeX source (arXiv:2106.10274): `sec_sumrules.tex` (`magicmoment`, `Bkv r eta`), `sec_heavy.tex` (`bulkpointresult`, `saddle`,
  `extrafactors`, `CfromBMellin`), `app_heavy.tex` (`definitionofomegatilde`, `CExpansionInPiHat`), `app_superconvergence.tex`
  (`magicblockidentity`, `c2intermsofsuperconvergence`, `cfoursuper`, `cknusuper`, `largenukernelexpression`);
* the project files `OMEGA_NF.md`, `OMEGA_U.md`, `SUPERCONVERGENT_UPLIFT.md`, `OMEGA_PRIME.md`.

Marks: **[Q]** quoted, **[D]** derived, **[D-sketch]**, **[DV]** derived and checked, **[NE]** numerical, **[M]** model, **[I]** assumed.

## 0. Verdict

1. **The project's `Ψ_k/N_k` are defined only through their leading Regge kernel.**
   * Their normalisation `N_k(ν)` is fixed at leading order, and the light-ray Jacobian `j_k(ν)` is untracked (`OMEGA_PRIME.md` item 5,
     `SUPERCONVERGENT_UPLIFT.md` §2).
   * So their exact near-gap action is **not defined by the project's own definitions**, for gravitons or for a scalar analogue.
   * The skeleton model replaces every moment by CHMRS's dispersive `C_{k,ν}` (exact, exactskel).
   * **H14 is therefore a genuine hypothesis.** It cannot be closed by computation within the current definitions [D].
2. **Can an exact (Mellin/Barnes) value be computed for a scalar analogue?** CHMRS show that `C_{k,ν}` *is* a superconvergence (light-ray)
   sum rule, `C₂ ∝ Φ_{−2,ν}[G/(2⁸u′v′)]`, with a meromorphic weight [Q]. That is why it has a Mellin form. The project's analogues have:
   * **even `k`:** the weight `sgn(v′−u′)r^{k−3}`, which is non-meromorphic. A Mellin form like `C`'s does not exist [D];
   * **odd `k`:** no scalar analogue at all. For identical scalars the odd sum rules vanish; the skeleton uses `k`-continuation [M].

   **So the exact near-gap value is not computable with `skelexact3`** (item 4 of the request applies).
3. **The size of the effect, quantified at the bulk-point level** (§3; [D-sketch] densities, [M] for their use on the exact skeleton). The
   flat-limit density of each candidate object, relative to `C_{k,ν}`, is an explicit function of `x` through the CHMRS saddle
   `r* = (1−s)/(1+s)`, `s = √(1−x²)`:

   | object | density ratio `R(k,x)` relative to `C_{k,ν}` | `x → 1` behaviour |
   |---|---|---|
   | `Π̂_{k,ν}` (Regge moment) | `(1−r*²)^k/(1+r*²)` | density `∼ (1−x²)^{−1}` |
   | scalar bare light-ray commutator, project weight | `(1−r*²)^{k+2}/(1+r*²) × S_k` | finite density for even `k` |

   * For comparison, `C`'s density `E_k` diverges like `(1−x²)^{−1−k/2}`.
   * The differences are `O(x⁴)` at small `x` and `O(1)` at `x ~ 1`.
4. **Effect on the near-gap claims (`Δ_g = 100`, exact skeleton, model of §3):**
   * **D-safe (`Ω^{NF,SC}`, NF02).** **`n = 1` survives for every candidate object and for every power-law family `R = (1−r*²)^α`,
     `α ∈ [−6, 12]` (all `k`, or `k ≥ 3` only).** The minimum margin is `F/N_s = 0.202–0.239` (baseline `0.227`). The flat limit stays
     positive with minimum `0.169–0.272` (baseline `0.200`). The NF02 weights vanish like `(1−p)⁵`, so the `x → 1` region where the objects
     differ carries little weight.
   * **D-strong (`Ω^U`, M02).** **`n` moves within `1.00–1.10`** (baseline `1.08`; `1.06–1.07` for the three candidate objects).
     `√(t₄²+α_T²)` at `Δ_g = 100` correspondingly moves within `3.2·10⁻⁵ … 4.5·10⁻⁵` (baseline `4.2·10⁻⁵`). The failing pocket (`+−`, `J = 4`,
     `m ≈ 1.04Λ`; `++`, `J = 0`, `m = Λ`) is exactly where the objects differ most.
   * **Not covered:**
     * the light side of an object change: contacts, `k(0)`, the `f̃` terms [I];
     * the true spin-2 light-ray radial profile of the graviton `Ψ_k`, if it is not of the power forms tested [I];
     * corrections to the multiplicative bulk-point model at finite `m` [M].

## 1. What each functional applies, term by term [Q from the project files]

Notation: `ν = Λp`, `x = ν/m`, `σ = ±1` for `++`/`+−`, `ε_k = 1` (even) or `½` (odd). The weights `φ_k(p)` are built from the flat
functional's `ψ₂, ψ₃` (and `a`) as in `FLAT_UNSTRIPPED.md` item 5: `φ₂ = ψ₂`, `φ₃ = ψ₃ − 4p²ψ₂`, `φ₄ = 2p⁴ψ₂ − p²ψ₃`,
`φ₅ = −3p⁶ψ₂ + 2p⁴ψ₃`, `φ₆ = ½(p⁸ψ₂ − p⁶ψ₃)`, `φ₇ = ½(−p¹⁰ψ₂ + p⁸ψ₃)` (`nf_common.NF.phi`).

| term | `Ω^U` (D-strong; M02 `ψ`, forward `a`) | `Ω^{NF,SC}` (D-safe; NF02 `ψ`, no forward terms) | object actually meant |
|---|---|---|---|
| `k = 2` | `∫φ₂ Λ⁴ C⁺_{2,Λp}` | `∫φ₂ Λ⁴ Ĉ₂`, `Ĉ₂ = C⁺₂` with `𝔖_ff⊗𝔖_ff` (kept, "project moment") | a spin-2 dispersive functional whose leading Regge moment is `Π̂₂` under the kernel (I10-type assumption); graviton pole from its leading kernel |
| `k = 3, 5, 7` | `∫φ_k ½Λ^{2k} C⁻_{k,Λp}` (project moments) | `∫φ_k ½Λ^{2k} Ψ_k/N_k`; `𝔄⊗𝔄`; weights `A₃ = 1`, `A₅ = (u′+v′)/(32u′v′)`, `A₇ = 1/(256u′v′)`; `j′ = 1` | `Ω^U`: assumed dispersive; `Ω^{NF,SC}`: explicit light-ray commutator sum rules, normalised by the leading-Regge `N_k(ν)` |
| `k = 4, 6` | `∫φ_k Λ^{2k} C⁺_{k,Λp}` | `∫φ_k Λ^{2k} Ψ_k/N^{ff}_k`; `𝔖_ff`; weight `sgn(v′−u′)r^{k−3}` (`= r^{k−3}` on the domain); `j′ = 0` | as above |
| forward | `−aΛ¹⁰∂_{ν²}C⁺_{4}|₀ + (5a/2)Λ¹⁰C⁻_{5,0}` (in `Ω^{SC}`: `Ψ₄/N₄^{ff}`, `Ψ₅/N₅` at `ν = 0`) | none | — |

**Heavy-action model** (both functionals; `exactskel/e3`, `e7`):

`ε_k m^{2k}Ĉ^h_k = e_k(x²)D_h + 2ε_kσ^k[P_{−½+iν}(cosh β) − J₀((2J+1)x)] + σ^k env_k ρ_k`,

where:
* `e_k` is the bulk-point (B′) density. For even `k` it is CHMRS's `E_k = (2−x²)/(1−x²)^{1+k/2}`, which is the density of the
  **dispersive `C_{k,ν}`**. For odd `k` it is the B′ form `σ(1−x²)^{−(k+1)/2}`.
* `ρ_k` is the exact remainder of the **scalar CHMRS `C_{k,ν}`** (exactskel; odd `k` by continuation in `k`).

**So the skeleton model identifies every term, near the gap, with CHMRS's dispersive `C_{k,ν}`.** It is **not** the Regge moment `Π̂_{k,ν}`
(the object `e2e`/`tttt`/`directgrav` compute), and **not** `Ψ_k/N_k` beyond leading order.

## 2. Relations between the objects for scalar-type heavy blocks

1. **`Π̂` and `C`** [Q].
   * Exact kernel: `C_{k,ν} = ∫∫Ω̃_{k,ν}(r,η)dDisc G`, with
     `Ω̃ = (a/γ²) r^{k−2}(1−r⁴)D^{−(k+1)/2}[P_ν + Σ_{n≥1}(…)(−4r²/D)ⁿI_n(η)]` and `D = (1+r²)²−4η²r²`.
   * Hence `(γ²/a)C_k = Π̂_k + O(Π̂_{k+2}, Π̂_{k+2,ν±2i})` (`CExpansionInPiHat`).
   * At large `ν` the kernel is `(γ²/a)Ω̃ → (1−r⁴)r^{k−2}D^{−(k+1)/2}P_ν`.
2. **`C` as light-ray / superconvergence sum rules** [Q].
   * `(γ²_{2Δφ+1}/a)C_{2,ν} = Φ_{−2,ν}[G/(2⁸u′v′)]`, where `Φ_{J,ν}` is the `u`-channel Lorentzian inversion integral.
   * `Φ_{−1,ν}[(v′−u′)G/(u′v′)] = 0` is the subtracted null-integrated commutator `∫dx₁⁺dx₃⁺((v′−u′)/(u′v′))⟨φ[φ,φ]φ⟩ = 0`.
   * The identity `(v′−u′)G̃_{−1,ν} = c₊G̃_{−2,ν+i} − c₋G̃_{−2,ν−i}` links the two.
   * `C₄` and the general `C_k` are finite combinations of `Φ_{−4}`, `Φ_{−2}` (and `Φ_{−k}`) at shifted `ν`, applied to `G/(2⁸u′v′)^{k/2}`.
   * **So `C_{k,ν}` is itself a light-ray commutator sum rule with a specific meromorphic weight.**
3. **A scalar bare commutator with weight `A`.** Define `Ψ^A := Φ_{−1,ν}[A·G]`.
   * Its large-`ν` kernel, from the identity above at large `ν` (`P_{μ+1} + P_{μ−1} ≈ 2ηP_μ`), `v′−u′ = η(1+r²)/(2r)` and `2⁸u′v′ = D²/r⁴`
     [D-sketch]:

     `K_A(r,η) ∝ (1−r²)D^{1/2}r^{−3}·A·P_ν(η)`.

     Check: for `A = (v′−u′)/(2⁸u′v′)` this reproduces `C₂`'s `(1−r⁴)D^{−3/2}` [DV, algebra].
   * The project's spin-2 `Ψ_k` has leading index `3` (`L[T]` commutator, `J₁+J₃−1 = 3`) and weights `A_k = r^{k−3}S_k`. The scalar
     analogue is the same shape on the scalar base: `A = r^{k+1}S_k`, with
     * `S_even = S₃ = 1`;
     * `S₅ = (1+6r²+r⁴)/(1−r²)⁴`;
     * `S₇ = (1−r²)^{−4}` (evaluated at `η = 1`).
   * **Caveat [I].** The graviton's own light-ray kernel is the `J = 3` one. Its radial profile at large `ν` is not derived, so the scalar
     analogue is only indicative.
4. **Computability.**
   * `C` (and CHMRS's alternatives `C̃_k = C₂[G/(2⁸u′v′)^{(k−2)/4}]` for `k ≡ 2 mod 4`, which share the leading moment and the large-`ν` limit)
     are Mellin-computable, because the weights are monomials in `u′, v′` [Q, D].
   * `Ψ^A` with `A = r^{k+1}` (even `k`) is not: `r` is algebraic in `(u′, v′)`.
   * The project's graviton `Ψ_k` is not defined beyond leading order (item 1 of §0).

   **⇒ No exact near-gap value of `Ψ_k/N_k` can be computed with `skelexact3`, for gravitons or scalars [D].**
5. **Bulk-point (flat-limit) densities** [D-sketch, from CHMRS `bulkpointresult` + `extrafactors`].
   * At the saddle `r* = ν²/(m+√(m²−ν²))² = (1−s)/(1+s)`, `η* = 1`, the density of a functional with kernel `r^{k−2}W(r,η)P_ν` is
     `E^Π̂_k(x)·W(r*,1)`, normalised to 2 at `x = 0`. Hence:

     | object | density |
     |---|---|
     | `C_k` | `E_k = (2−x²)/s^{k+2}` |
     | `Π̂_k` | `2^{2k−1}/(s²(1+s)^{2k−2})` |
     | `Ψ^{r^{k+1}}` | `2^{2k+3}/(1+s)^{2k+2}` (finite at `x = 1`) |
     | odd `k` | the last row times `S_k(r*)` |

   * The ratios `R(k,x)` at `x = 0.3 / 0.6 / 0.9 / 0.99`, e.g. `k = 4` (`cache/o1_flat.txt`):

     | object | 0.3 | 0.6 | 0.9 | 0.99 |
     |---|---|---|---|---|
     | `Π̂` | 0.997 | 0.940 | 0.443 | 0.023 |
     | `Ψ` | 0.996 | 0.917 | 0.317 | 0.004 |

## 3. Quantification on the near-gap actions

**Model [M].** Object `X`'s moment = `R_X(k,x)` × (exactskel model moment), node by node (64 Gauss nodes; exact `ρ_k`, Df = 3 and 2):
* the curvature term and the forward terms (`ν = 0`, `R = 1`) are kept;
* this is exact at leading order in `1/m` at fixed `x`;
* with `R = 1` the node sum reproduces `F_ex` of e7 to `≤ 6·10⁻¹⁰N_s` (`o2` validation).

**Variants:**
* `C` (baseline);
* `Π̂`;
* `Ψ` (shapes `S_k`);
* `PsiPure` (`S = 1`);
* `mix` (`k = 2` → `C`, `k ≥ 3` → `Ψ`: the actual `Ω^{NF,SC}` assignment);
* `o3`: the power-law family `R = (1−r*²)^α`, `α ∈ {−6, −4, −2, 2, 4, 8, 12}`, applied to all `k` or to `k ≥ 3` only.

### 3.1 Flat limit, NF02 heavy side (`o1_flat.py`) [D-sketch + NE]

Minimum `F/N_s` over `J ≤ 80` and `m/M ∈ [1, 10]`:

| object | `++` | `+−` |
|---|---|---|
| `C` (= B′, the identity checked to `10⁻¹³`) | 0.2000 (`m = M`, J = 4) | 0.2000 (J = 5) |
| `Π̂` | 0.2006 | 0.1700 (`m = M`, J = 17) |
| `Ψ` | 0.2408 | 0.1776 (J = 17) |
| `PsiPure` | 0.1977 | 0.1687 (J = 17) |
| `mix` | 0.2715 (J = 6) | 0.1782 (J = 17) |

Above `m = 1.5M` all variants agree to `≤ 10⁻³`.

### 3.2 AdS near-gap, `Δ_g = 100` (`o2_ads.py`, `o3_envelope.py`) [M + NE]

**Binding states** (Df = 3; `F/N_s`; `cache/o2_ads.txt` has all of them):

| functional | state | C | mix | Ψ | Π̂ |
|---|---|---|---|---|---|
| NF02 | `+−` 1.04Λ J=4 | 0.241 | 0.244 | 0.240 | 0.242 |
| NF02 | `++` 1.04Λ J=4 | 0.266 | 0.320 | 0.296 | 0.271 |
| NF02 | `++` 1.10Λ J=6 | 0.348 | 0.333 | 0.329 | 0.338 |
| NF02 | `++` 1.30Λ J=6 | 0.316 | 0.323 | 0.313 | 0.314 |
| NF02 | `++` Λ J=0 | 0.998 | 1.029 | 1.001 | 0.988 |
| M02 | `+−` 1.04Λ J=4 | **−1.391** | −0.465 | −0.472 | −0.254 |
| M02 | `++` Λ J=0 | **−0.099** | +5.65 | +5.51 | +0.416 |
| M02 | `++` 1.04Λ J=2 | 0.412 | 0.230 | 0.238 | 0.423 |
| M02 | `++` 1.30Λ J=2 | 0.254 | 0.245 | 0.247 | 0.257 |

**Thresholds over all scanned states** (both Df, `J ≤ 10` plus the high-`J` sample):

| functional | object | `n` | min `F/N_s` | `√(t₄²+α_T²)` at `n` |
|---|---|---|---|---|
| NF02 (D-safe) | C / mix / Ψ / PsiPure / Π̂ | **1 / 1 / 1 / 1 / 1** | 0.227 / 0.218 / 0.234 / 0.227 / 0.225 | 2.21·10⁻⁴ (unchanged) |
| NF02 | `(1−r*²)^α`, `α ∈ [−6, 12]` | **1** for all `α` (all `k` or `k ≥ 3`) | 0.202 – 0.239 | 2.21·10⁻⁴ |
| M02 (D-strong) | C / mix / Ψ / Π̂ | 1.08 / 1.07 / 1.07 / 1.06 | −1.39 / −0.47 / −0.48 / −0.27 | 4.24 / 4.09 / 4.09 / 3.94 ·10⁻⁵ |
| M02 | `(1−r*²)^α`, `α ∈ [−6, 12]` | 1.10 (`α = −6`) … 1.00 (`α = 12`) | −8.1 … +0.15 | 4.5·10⁻⁵ … 3.2·10⁻⁵ |

## 4. Grades and what remains

| statement | grade |
|---|---|
| the skeleton model uses CHMRS dispersive `C_{k,ν}` (even-`k` density `E_k`; exact `C` remainder) for every term | [D] (read off e3/e7, `remk`) |
| `C_{k,ν}` is a superconvergence / light-ray sum rule with a meromorphic weight; it has a Mellin form | [Q] CHMRS app. superconvergence |
| `Ψ_k/N_k` defined only via the leading Regge kernel (`N_k`, P1, untracked `j_k(ν)`); its subleading and near-gap action are undefined | [D] from the project files |
| no Mellin/Barnes form for the scalar analogue of even `Ψ_k` (non-meromorphic weight); no scalar odd analogue | [D] |
| bulk-point densities of `Π̂`, `C`, scalar `Ψ^A` (formulas of §2.5) | [D-sketch] (large-`ν` kernel via the identity of §2.3; CHMRS saddle) |
| D-safe `n = 1` and margin `≥ 0.20` robust to the object (candidates and `α ∈ [−6, 12]`) | [NE + M] |
| D-strong `n ∈ [1.00, 1.10]` depending on the object | [NE + M] |

**What remains:**
* **Light side of an object change [I].**
  * The forward-free contact cancellation of NF02, `k(0)` and the `f̃` coefficients are derived for B′/`C`-type densities.
  * A different object with the same leading Regge moment differs by higher-spin sum rules. These annihilate contacts below their spin but
    not above, so the change is power-suppressed in the EFT counting (H10) [D-sketch].
  * `C⁺₂` (graviton pole) is unchanged in the `mix` assignment.
* **The graviton `Ψ_k`'s own large-`ν` kernel** (spin-2 light-ray, `J = 3`) and `j_k(ν)` [I]. The envelope `α ∈ [−6, 12]` covers power-law
  deviations only.
* **The multiplicative bulk-point model** at finite `m` [M]. At `m ≈ Λ`, `x → 1` it is least reliable, but there the NF02 weights vanish
  like `(1−p)⁵`.

**Recommendation.** Make the object identity a definition rather than a hypothesis by choosing Mellin-computable dispersive functionals:
* `C_{k,ν}` itself, or CHMRS's `C̃_k`, or `Φ_{−2,ν}` with meromorphic weights;
* for the spinning case, the CHMRS-suggested spinning extension via superconvergence `Φ_{−2}`;
* in place of the leading-order-defined `Ψ_k/N_k`.

The exact near-gap action is then computable for scalars now (exactskel), and for gravitons once spinning Mellin amplitudes are available.

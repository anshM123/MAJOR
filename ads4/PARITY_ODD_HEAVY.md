# H12: parity-odd and indefinite-parity heavy exchanges

Status: written 18 Nov 2026 by the `⟨TTTT⟩` derivation/numerics agent; **not refereed; not committed.**
* Code: `ads4/code/exactskel/e11_parity.py` (output `exactskel/cache/e11.txt`). It reads the per-`k` heavy parts saved by `e10_kernel.py`.
* Paper text: arXiv:2201.06602 v1 (PDF fetched 18 Nov 2026; the equation numbers below are v1).
* Builds on `FLAT_GRAVITON_BOUND.md`, `ROBUST_FUNCTIONAL.md`, `NOFORWARD_FUNCTIONAL.md`, `KERNELS_D3.md` (helicity map),
  `REGGE_MOMENT_TTTT.md` §4–5, `SPINNING_4PT.md` §6 (no-`O(1/m)` argument), `SUPERCONVERGENT_UPLIFT.md`, `EXACT_SKELETON_NEARGAP.md` §9.

Marks:
* **[D]** derived;
* **[Q]** quoted from the literature (checked against the source text);
* **[NE]** numerics;
* **[I]** inferred / assumed;
* **[ND]** not done.

## 0. Verdict

| item | grade |
|---|---|
| Flat side: every flat functional used (`F_{g3}`, M01, M02, NF02) is built only from the MHV amplitude `f`. Its heavy side is `Σ_X [|c^{+−}_X|² F_{+−} + |c^{++}_X|² F_{++}]`, whatever the parity of `X`, in parity-violating theories too. There are no cross terms, so the certificates already cover parity-odd and indefinite-parity heavy states. | **[D]** (with the paper's (2.13)–(2.16) [Q]) |
| AdS side, interference: every AdS functional is parity-invariant (kernels `𝔖_ff`, `K*`, `𝔄⊗𝔄`, `j′`-matched superconvergent kernels). So the moment is **exactly even in the parity-odd OPE coefficient**: no `λ_E λ_O` term at any order in `1/m` or `ν/m`. The phase-sensitive flip terms (single flips `c_s c_o*`, double flips `(c^{++})²`) vanish identically for `𝔖_ff`, `K*`, `𝔄⊗𝔄`. They do **not** vanish for the superconvergent `Ψ₄`, `Ψ₆` (the `(1,3)` `j′ = 0` structure has a flip response). There the parity-odd part of the interference, `∝ λ_Eλ_O`, is removed by parity invariance; what remains is the parity-even residual `ρ_k` already bounded in `SUPERCONVERGENT_UPLIFT.md` item 3 (corrected 18 Nov, referee). | **[D]** (given Lemma L flip responses [DV]) |
| AdS side, leading heavy formula: a parity-odd (or mixed) exchange enters as `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)`, with `c_same` the full complex helicity coupling (`∝ r_Eλ_E + i r_Oλ_O`). At even `J` the parity-odd structure feeds only `c_same`. | **[D given REGGE_MOMENT_TTTT (H1)]**; not tested end to end on a parity-odd structure **[ND]** |
| AdS side, no `O(1/m)`: mixed waves `W[E,O]` drop out exactly (above). For `W[O,O]` the reflection argument of `SPINNING_4PT.md` §6.1 goes through: two `ε`-tensors give a squared sign. The step to physical couplings needs the parity-odd MFT normalisation to be Weyl-even. | **[D given Assumption T, Lemma S]** for the structures; **[I]** for the normalisation |
| Margin if an `O(1/m)` term did exist (`e11`): thresholds unchanged for `\|c₁\| ≤ 8.4` (D-strong, `++`, `Δ_g = 100`, `δ = 0.03`), `≤ 19` (D-safe). The mutation scale is `\|c₁\| ≈ 4–10`. | **[NE]** |

**H12 is settled except for one reduced item.** Positivity of D-safe and D-strong extends to parity-odd and indefinite-parity heavy states,
with no interference terms and no change in any number. What remains is **H12′ = H1 for parity-odd `⟨TTO⟩` structures**: the two-pair
heavy formula, and the absence of `O(1/m)` through the Weyl parity of the odd MFT normalisation, are derived but not tested end to end.

Even if H12′ failed by a mutation-size `O(1/m)` term, D-safe would not move: `n = 1` at all `Δ_g` for `|c₁| ≤ 19–150`. D-strong moves
at most from `1.08` to `1.09` (bound `4.2 → 4.4·10⁻⁵`) for `|c₁| ≤ 8` at `Δ_g = 100`.

## 1. Flat side

### 1.1 What the paper does [Q]

arXiv:2201.06602 v1 treats parity-violating theories explicitly.
* "Hatted couplings are complex (the real and imaginary part representing parity-even and parity-odd couplings, respectively)" (below
  (2.7)). So `ĝ₃` (in `f_low` through `|ĝ₃|²`) contains the parity-odd cubic coupling: `α_T` in the AdS dictionary.
* Unitarity (2.13): `i[(a_J^{−h₄,−h₃,−h₂,−h₁})* − a_J^{h₁h₂h₃h₄}] = Σ_X (a_J^{−h₃,−h₄→X})* a_J^{h₁,h₂→X}`, abbreviated `2(c^{−h₄,−h₃})* c^{h₁h₂}`.
  "The right-hand-side is a positive semi-definite matrix."
* Specialised to the MHV amplitude (2.14): `Im a_J^{+−+−} = |c^{+−}_J|²`, `Im a_J^{++−−} = |c^{++}_J|²`, "positive". Here `c^{+−}` is
  real and `c^{++}` complex: "the complex couplings `ĝ₄` and `c^{++}_{J,m²}` combine parity even and odd parts" (below (2.37)).
* (2.16a,b): the `s`-channel cut of `f` carries `|c^{+−}|²`, all `J ≥ 4`; the `t`-channel cut carries `|c^{++}|²`, `J ≥ 0` even.
* The non-MHV amplitudes carry interference (2.17)–(2.19): `Im a^{+++−} = c^{++}c^{+−}`, `Im a^{++++} = (c^{++})²`, with
  `Im g ∝ Σ c^{++}c^{+−} d̃^J_{4,0}` and `Im h ∝ Σ (c^{++})² d̃^J_{0,0}`.
* The heavy three-point amplitudes (A.7), (A.8) are fixed by Poincaré symmetry up to the constants `c^{+−}_{J,m²}`, `c^{++}_{J,m²}`. They make
  no reference to the parity of `X`.

### 1.2 Derivation [D]

* The two cuts of `f` are diagonal elements of the absorptive part, one per channel:
  * the `s` channel is `(1⁺2⁻) → (3⁻4⁺)`, i.e. the physical elastic process `+− → +−`;
  * the `t` channel is `(1⁺4⁺) → (2⁻3⁻)`, i.e. incoming `++` to outgoing `++` (an incoming negative helicity is an outgoing positive one).
* By (2.13) each is `Σ_X |⟨X|T|h₁h₂⟩|² ≥ 0` (generalised optical theorem for an elastic diagonal element). Nothing here uses parity, and
  states of either intrinsic parity, or of none, appear on the same footing.
* Bose symmetry of two identical same-helicity gravitons allows only even `J` in `++`, for either parity of `X`. Both a `0⁺` (`φR²`) and a
  `0⁻` (`φRR̃`) couple to `++`.
* A state that couples to both channels appears in `f` as `|c^{++}_X|²` at `t = m²` and `|c^{+−}_X|²` at `s = m²`. These are two separate
  terms, so its action is `|c^{++}_X|² F_{++}[m,J] + |c^{+−}_X|² F_{+−}[m,J]`. The certificates (Theorem 1 of `FLAT_GRAVITON_BOUND.md`;
  Theorem R of `ROBUST_FUNCTIONAL.md`, `𝔉 ≥ (9/50)N` per channel; NF02) are per channel. They therefore bound every such sum: there is no
  cross term `c^{++}c^{+−*}` or `c^{++}c^{−−*}` to control.
* Example: a scalar with `φ(aR² + bRR̃)` (`a`, `b` real for a Hermitian Lagrangian; parity-violating when `ab ≠ 0`).
  * The self-dual decomposition gives `c^{++} ∝ a + ib` and `c^{−−} ∝ a − ib = (c^{++})*` (CPT for a self-conjugate field).
  * The `t`-channel residue of `f` is `c^{++}c^{−−} ∝ a² + b² ≥ 0`.
  * The all-plus `h` carries `(c^{++})² ∝ a² − b² + 2iab`. Here the parity-odd part enters with the opposite sign, which is the
    interference the MHV functionals avoid.
* **Which functionals.**
  * `F_{g3}` (completed), M01, M02 (`ROBUST_FUNCTIONAL.md` item 1) and NF02 use only `B_k^{(1)}`, i.e. `f(s,−p²) ± f(p²−s,−p²)`.
  * The JHEP variant adds `B₄^{(2)}`: `f` at fixed `t`, whose cuts are `s` and `u`, both `+−`, still diagonal.
  * None uses `g` or `h`.

  **So H12 is not an assumption on the flat side** [D].
* **Light side.** The parity-odd cubic coupling is inside `|ĝ₃|²` (2.7a), and light matter of either parity contributes `a_J = |g_J|²`
  (B.1). That is light-side physics already in Theorem 2, not H12.

### 1.3 What would need more

A functional that used `g` or `h`, like the paper's §4 bounds on `ĝ₄` (Fig. 10), would face `c^{++}c^{+−}` and `(c^{++})²`. It would
need matrix positivity, as in the paper's footnote 6 (`S`-matrix blocks of `|++⟩`, `|−−⟩`, `(|+−⟩+|−+⟩)/√2`). The project's functionals do
not use them.

## 2. AdS side

### 2.1 Helicity couplings of an exchanged operator [Q + D]

`KERNELS_D3.md` §5.1 (all conserved `⟨TTO⟩` structures, `J ≤ 12`, symbolic in `Δ`):
* **Even `J`.** The parity-odd structure has `P_{m≥4} = 0`: it feeds only the same-helicity sector `a^{++}`, "with relative phase `i`".
  The two parity-even structures feed `a^{++}` (combination `S^{++}`) and `a^{+−}` (combination `S^{+−}`). With `h = −2`:
  `â_{m,−2}(λ_E, λ_O) = â_{−m,2}(λ_E, −λ_O)`.
* **Odd `J ≥ 5`.** There is a single structure, which feeds only `a^{+−}`.
* Hence, for real three-point coefficients, `c^{±±} ∝ r_Eλ_E ± i r_Oλ_O`, `c^{±∓} ∝ r′_Eλ′_E`, and `|c^{++}|² = |c^{−−}|² = r_E²λ_E² + r_O²λ_O²`.
  * Reality of the coefficients of a real structure basis is the CFT image of CPT, i.e. CRT of the Lorentzian theory [I, standard]. It is
    what makes `ρ^h = |λ^h|²` the positive diagonal block coefficient of `REGGE_MOMENT_TTTT.md` (H1), "each summed over the two parity
    images".
  * A theory without parity simply has both `λ_E` and `λ_O` nonzero. It never has a different modulus for `++` and `−−`.

### 2.2 All AdS functionals are parity-invariant ⇒ no `λ_Eλ_O` terms, to all orders [D]

The kernels are:
* even `k`: `𝔖_ff = D₁²+D₁D₂` or `K* = (D₁+D₂)(D₁−½D₂)` at both pairs. `D₁`, `D₂` are scalar products, parity-even.
* odd `k`: `𝔄⊗𝔄`, `𝔄 = D₅D₁` with `D₅ = ε(z₁,z₃,ŷ)` parity-odd at each pair, so the product is even.
* superconvergent replacements: `(1,3)` celestial structure of parity `j′ = k mod 2` times the `(2,4)` lift of the same parity
  (`SUPERCONVERGENT_UPLIFT.md`), so the product is even.

The Regge weight `P_{−1/2+iν}(−x̂·ŷ)`, the radial weights and `A_k` are parity-even. So `ω∘𝒫 = ω` for every functional of D-safe and
D-strong.

Under parity a `⟨TTTT⟩` block with coefficients `(λ_E, λ_O)` at each vertex goes to the block with `(λ_E, −λ_O)`, because each odd
structure carries one `ε`. Hence

`ω[G_O] = λ_E² ω[W_{EE}] + λ_O² ω[W_{OO}] + λ_Eλ_O · 0`,

exactly: at every `m`, `ν`, `β`, not only at leading order. Mixed waves cannot produce interference or an `O(1/m)` term of their own.

### 2.3 Leading two-pair formula; flip terms [D given H1 of `REGGE_MOMENT_TTTT.md`]

* The proposition of `REGGE_MOMENT_TTTT.md` §5 is stated in the helicity basis and uses only the helicity content of the vertices.
  * The helicity-preserved part is `𝒩 (π/4)γ² m^{−2a} P_{−1/2+iν}(cosh β){κ₊M_𝔖²(ρ^{++}+ρ^{+−}) + κ₋|M_𝔄|²(ρ^{++}−ρ^{+−})}`.
  * With §2.1, `ρ^{++} = r_E²λ_E² + r_O²λ_O²`.
  * So a parity-odd exchange (or the odd part of an indefinite-parity one) is an ordinary same-helicity contribution. It enters positively in
    `C⁺` and in `C⁻` (`++ ↔ σ = +1`), exactly as in the flat `++` channel.
  * The near-gap scans (`e3`, `e7`, `e10`), which are per `|c^h|²` with `++` at even `J`, already cover it.
* **Flip terms.** These are the only phase-sensitive pieces of a two-pair moment:
  * one pair flipped: `∝ c_s c_o*`;
  * both pairs flipped: `∝ λ^{++}λ^{++} = (c^{++})² ∝ λ_E² − λ_O² + 2iλ_Eλ_O`, where a parity-odd exchange enters with the opposite sign.

  They are multiplied by the pair flip responses. These vanish identically for `𝔖_ff`, `K*` (`exactskel/cache/evenkernel.txt`) and `𝔄`
  (`superconv/cache/s1.log`: flip even = flip odd = 0), so there the flip terms are absent [D].
  * **Correction (referee, 18 Nov).** The earlier sentence "flip responses vanish" is false for the superconvergent replacements `Ψ₄`,
    `Ψ₆`, which use the `(1,3)` `j′ = 0` celestial structure: its flip response is nonzero (`s1.log`). The conclusion still holds, for a
    different reason. The surviving single flip `(1,3) × preserved (2,4)` is `∝ c_s c_o*`. At even `J` only the parity-even structure
    feeds `c_o`, so `c_s c_o* = (r_Eλ_E + i r_Oλ_O) r′_Eλ_E`. Its `λ_Eλ_O` part is excluded by the parity invariance of the functional
    (§2.2). What is left, `∝ λ_E²`, is the parity-even residual `ρ_k` bounded in `SUPERCONVERGENT_UPLIFT.md` item 3
    (`≤ 6.5·10⁻⁴ N_s`).
  * With the old `D₁⊗D₁` kernel they would not be: the double-flip interference flips sign for parity-odd exchanges. That is a second
    reason, besides `SUPERCONVERGENT_UPLIFT.md` item 6, to keep flip-free kernels.
* **Not tested end to end [ND].** The two-pair checks (`⟨JJJJ⟩`, `⟨TTTT⟩`, all orders `k = 2–7`) were run on parity-even exchanged
  structures only. `SPINNING_E2E.md` / `SPINNING_4PT.md` §6.3 (referee M4) list parity-odd structures as not built, because they need the
  `ε`-operator `D̃` at both points.

### 2.4 `O(1/m)` for parity-odd exchanges [D given Assumption T, Lemma S; I for the normalisation]

* **Mixed waves `W[E,O]`.** They drop out exactly (§2.2).
* **`W[O,O]`.** The argument of `SPINNING_4PT.md` §6.1 maps `b𝒯^{ab}(z₁,z₂,z₃,z₄;x,y;A_±)` to `b𝒯^{ba}(z₃,z₄,z₁,z₂)` under
  `(x,y,A_±) → (−x,−y,−A_±)`.
  * For a parity-odd structure each three-point factor picks up the sign of its `ε` under the point reflection (`det = −1` in `d = 3`).
  * In `W[O,O]` the two signs multiply to `+1`, so `t₁^{OO}(z₁…z₄) = −t₁^{OO}(z₃,z₄,z₁,z₂)` holds as for even structures.
  * The two-pair mode symmetry (Lemma S: pair exchange with `v → −v`, a single reflection) likewise acts with a squared sign.
  * So no mode of `W[O,O]` has an `O(1/m)` term in the Weyl-odd variables (given Assumption T, as for parity-even exchanges).
* **Physical couplings.** The step from the structure basis to `|c^{++}|²` needs the MFT normalisation of the parity-odd helicity
  structure to be Weyl-even.
  * For the parity-even structures this is [D from Q] at `l = 1` and [ID] at `l = 2` (`SPINNING_TTTT.md`).
  * For the odd structure the only input is CHL's statement that the odd MFT data are `½c(1,2)^s` times the same helicity matrix
    (`SPINNING_4PT.md` §6.3) [Q, not re-derived]. If so, its Weyl parity is that of the even case.
  * The `(R0)` universality of `KERNELS_D3.md` (all structures with `P₀ ≠ 0`, odd included, have the same helicity profile
    `|â_{−k}|²/|â₀|² = W⁰_k(Δ,J)`, Weyl-even in `Δ`) points the same way. Neither is a proof: **[I]**.

### 2.5 Margin if an `O(1/m)` term did exist [NE] (`e11_parity.py`, `cache/e11.txt`)

The model multiplies every heavy moment of the affected channel by `1 + c₁/m` with unknown sign, and bounds the effect by
`|c₁|/m (Σ_k |P_k| + |FW|)` on the exact-skeleton states (`K*`, both skeletons, high-`J` samples). This is conservative: it ignores the
strong cross-`k` cancellations. Scopes:
* "`++`" is the parity-odd exchanges of §2.1;
* "all" covers the case where the odd-`J` `+−` structure is also affected.

Table: largest `|c₁|` that leaves `n` unchanged (`δ = 0` / `δ = 0.03`):

| `Δ_g` | D-strong `++` | D-strong all | D-safe (both scopes) |
|---|---|---|---|
| 100 | 8.7 / 8.4 | 4.1 / 2.2 | 21.7 / 19.2 |
| 300 | 15.1 / 14.7 | 7.8 / 2.1 | 51.2 / 44.5 |
| 1000 | 9.5 / 8.9 | 9.5 / 8.9 | 152 / 131 |

* Thresholds at mutation-size `|c₁|`, D-strong, `δ = 0.03`, scope "all":
  * `c₁ = 4`: `n = 1.09 / 1.06 / 1.01`; bounds `4.4·10⁻⁵ / 5.5·10⁻⁷ / 4.1·10⁻⁹`, against `4.2·10⁻⁵ / 5.3·10⁻⁷ / 4.1·10⁻⁹`.
  * `c₁ = 16`: `n = 1.15 / 1.06 / 1.02`; bounds `5.4·10⁻⁵ / 5.5·10⁻⁷ / 4.3·10⁻⁹`.
* D-safe stays at `n = 1` for every `|c₁| ≤ 16` at all three `Δ_g`.
* For scale: measured `|c₁|` for parity-even exchanges is `≤ 0.2` (`⟨JJJJ⟩`) and `≤ 0.05` (`⟨TTTT⟩`, `K*`). Non-Weyl-even mutations give
  `|c₁| ≈ 4–10` (`SPINNING_4PT.md` §6.2).

## 3. Conclusion

1. **Flat [D].** Parity-odd and indefinite-parity heavy states enter the MHV sum rules only through `|c^{++}|²` and `|c^{+−}|²`. The
   certified positivity of `F_{g3}`, M01, M02 and NF02 covers them. No new assumption is needed, and no cross terms arise.
2. **AdS, structure [D].** The D-safe and D-strong functionals are parity-invariant and flip-free. So:
   * a parity-odd exchange contributes `λ_O²` times its own diagonal moment, and never `λ_Eλ_O`, at any order;
   * the phase-sensitive flip terms vanish identically for `𝔖_ff`, `K*`, `𝔄`. For `Ψ₄`, `Ψ₆` only their parity-even part survives,
     which is the bounded `ρ_k` (correction of §2.3).
3. **AdS, leading action [D given H1].** It equals that of a same-helicity parity-even state with the same `|c^{++}|²`. The near-gap scans
   (per `|c^h|²`) already cover it. Margin effect: none.
4. **Remaining assumption: H12′.** This is the two-pair heavy formula (H1 of `REGGE_MOMENT_TTTT.md`) for the parity-odd `⟨TTO⟩`
   structure, including the absence of `O(1/m)`. That absence rests on the Weyl parity of the odd MFT normalisation [I]. Even a
   mutation-size violation leaves D-safe unchanged and moves D-strong's `n` by `≤ 0.01` (`|c₁| ≤ 8` at `Δ_g = 100`). Proposed grade for H12
   in `RESULT_D.md`: **"[D] flat and interference-free; reduced to H12′ (H1 on the parity-odd structure) [I], with tolerance `|c₁| ≲ 8`
   (D-strong) / `≲ 19` (D-safe)".**
5. **Shortest test of H12′ [ND].**
   * What to build: the parity-odd conserved `⟨TTO⟩` (or `⟨JJO⟩`) structure with the `ε`-operator `D̃` in the Taylor-mode code, then one
     end-to-end point (`β = 1`, `k = 4–6`, three masses, as in `SPINNING_TTTT_EXT.md`) for `W[O,O]` and `W[E,O]`.
   * What it should show: `W[E,O]` vanishes identically under `K*⊗K*` and `𝔄⊗𝔄`, which is a sharp check of §2.2. `W[O,O]` gives the
     same-helicity formula with `|c₁| ≪ 1`.
   * Cost estimate: a few hours of development, plus production runs of about 2.5 h per mass on one core (as for the parity-even waves).

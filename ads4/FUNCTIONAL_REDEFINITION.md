# Making the object identity a definition: which sum rules should replace `Ψ_k/N_k`?

Status: 18 Nov 2026, numerics agent. **Unrefereed. Not committed.** Follows `OBJECT_IDENTITY.md` (H14).
* Code: `ads4/code/objid/r1_redef.py` (output `objid/cache/r1_redef.txt`). It imports `omegau/` and `exactskel/` (`evenkernel`, e1/e7 caches)
  read-only.
* Sources: CHMRS LaTeX (arXiv:2106.10274; `sec_light.tex`, `app_superconvergence.tex`, `app_heavy.tex`), `LIGHT_INPUTS.md`, `OMEGA_PRIME.md`,
  `SUPERCONVERGENT_UPLIFT.md`.

Marks: **[Q]** quoted, **[D]**, **[D-sketch]**, **[NE]**, **[M]**, **[I]**.

## 0. Result

1. **For scalar-type states the natural definition is CHMRS's dispersive `C_{k,ν}` (even `k`).**
   * It satisfies all five constraints (§1, §2).
   * It is already the object the skeleton model computes exactly. So with `Ĉ_k :≡ C_{k,ν}` the near-gap action of the D-safe functional on the
     scalar skeleton **is** the exactskel result, with no object-identity gap for `k = 2, 4, 6`.
   * Recomputed (`r1_redef.py`; odd `k` bracketed by the B′ and E readings): **`n = 1` at `Δ_g = 100 / 300 / 1000`, min `F/N_s = 0.227 / 0.206 / 0.202`
     (B′) and `0.228 / 0.208 / 0.204` (E). `√(t₄²+α_T²) ≤ 2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸`** [NE; unchanged from `RESULT_D.md`].
   * The flat limit is then **derived** (CHMRS bulk-point density `E_k` = the even B′ density) rather than assumed via the B′ dictionary.
2. **Odd `k` has no scalar functional** (identical scalars, `s↔t` antisymmetry). At scalar level it remains a model (H3): the two readings differ
   by `≤ 0.002 N_s` at the minimum.
3. **For gravitons the redefinition is a trade, not a free gain.**
   * `Ψ_k/N_k` are explicit sum rules, Einstein-free given RB₂ (the reason they replaced the project's `C_k`), but undefined near the gap.
   * A spinning CHMRS-type `C_k` would be well defined and computable near the gap, but it has to be **constructed**. Without a construction one
     is back to the "completion" assumption (I10) that `SUPERCONVERGENT_UPLIFT.md` removed.
   * The construction that keeps both properties is identified in §3:
     * meromorphic-weight light-ray commutators;
     * for even `k`, followed by the CHMRS `η`-inversion.

     **For scalars these reduce exactly to CHMRS's `C₂` / `C̃_k`.**
   * Their exact near-gap evaluation for gravitons needs spinning Mellin amplitudes for heavy `⟨TTTT⟩` exchanges, which do not exist here
     (`SPINNING_MELLIN.md`). **Distance: research-level (weeks or more).**
4. **What is removed, and what remains:**
   * **Removed at scalar-skeleton level:** H14 for `k = 2, 4, 6` (object defined as `C_{k,ν}` and computed exactly), and the B′ dictionary for even `k`
     (flat limit derived).
   * **Remaining:**
     * odd `k` (H3);
     * for gravitons, the existence of the spinning analogue with properties (i)–(v) (call it **H14′**, [I]) plus the graviton-vs-skeleton
       transfer for that object, which **merges into H4/H1**.

## 1. Constraints on a replacement for `Ψ_k/N_k` (`k = 3…7`)

| | constraint | why |
|---|---|---|
| (i) | same leading heavy Regge action `(2P_{−½+iν}(cosh β)/m^{2k})(\|c_s\|² ± \|c_o\|²)` **and** the B′ bulk-point density `e_k(x²)` | NF02's flat limit and certified flat positivity are statements about the bulk-point densities (`OBJECT_IDENTITY.md` §2.5); leading Regge equality alone does not fix them |
| (ii) | no Einstein term given `J₀(G^E) ≤ 2` | the reason `Ψ_k` replaced the project's `C_k`; the latter were defined only by their leading kernel, and `C⁻₃` is not the leading moment of any bare commutator (`LIGHT_INPUTS.md` item 3) |
| (iii) | no forward (`ν = 0`) pieces | D-safe |
| (iv) | exactly computable near the gap: scalar skeletons now, gravitons in principle | to remove H14 |
| (v) | well-defined light side (graviton pole only at `k = 2`; contacts and `[TT]` analysable; `k(0)`, `A`, `B_AdS`) | the bound |

## 2. Candidates for scalar-type states

| candidate | (i) leading / bulk-point | (ii) Einstein-free (`k ≥ 3`) | (iii) | (iv) exact | (v) light side |
|---|---|---|---|---|---|
| **CHMRS `C_{k,ν}`** (even `k`) | `Π̂_k` + `O(Π̂_{k+2})` [Q]; bulk-point `E_k(x)` = even B′ density [Q, `bulkpointbknuagain`] | yes for `k ≥ 4`: `ω\|_light` comes only from EFT Regge growth `n > k` [Q, `sec_light` + app. omega light], so graviton exchange (`n = 2`) and `n < 2` terms vanish | yes | **yes**: `CfromBMellin`; Arb Barnes (`skelexact3`) | yes [Q]: support on the first `k/2` double-trace trajectories; contacts via `B̂_{k,T}[M_{a,b}]` (polynomial in `T`); zero on contacts with zeros at `S = 2Δφ…2Δφ+k−2`; graviton pole in `C₂` only |
| CHMRS `C̃_k` (`= C₂[G/(2⁸u′v′)^{(k−2)/4}]`, `k ≡ 2 mod 4`; `C₄`-based for `k ≡ 0 mod 4`) | same leading moment; same large-`ν` kernel at leading order (`largenukernelrelation`), so the same bulk-point density | yes, as above (spin-`k` superconvergence) | yes | yes in principle: the weight is a monomial, so a shifted Mellin amplitude with double-trace double poles; not implemented | yes [D-sketch] |
| `Φ_{−k}[G^{(k)}]` alone (leading part of `C_k`) | same leading moment; bulk-point density not derived | yes | yes | yes (Mellin monomial weight) | yes |
| project-type bare commutator, weight `r^{k+1}` (scalar analogue of even `Ψ_k`) | leading yes; bulk-point `2^{2k+3}/(1+s)^{2k+2}` ≠ `E_k` (`OBJECT_IDENTITY.md` §2.5) | yes | yes | **no** (non-meromorphic weight) | (P) not analysable (no `dDisc` rewriting) |
| Regge moment `Π̂_{k,ν}` | leading yes; bulk-point differs | **not a sum rule** (`Π̂` does not vanish on physical correlators, CHMRS) | — | no (cancellation; `SPINNING_MELLIN.md` §5.3) | — |

**Conclusion [D from Q].** For even `k`, `C_{k,ν}` satisfies (i)–(v) at scalar level. `C̃_k` does too, with a different subleading structure but the
same bulk-point density.

**Why did the project not use it?** It never rejected CHMRS's *scalar* `C_k`.
* What was rejected (`LIGHT_INPUTS.md`, `SUPERCONVERGENT_UPLIFT.md`) are the project's *graviton* moments `C^±_k`. These are defined only through
  their leading tensor kernel (`D⊗D`, `𝔄⊗𝔄`).
* Their completion to a dispersive functional was assumed (I10), and for `C⁻₃` it was shown **not** to be the leading moment of any bare-commutator
  functional.
* So the problem is the **spinning completion**, not `C_k` itself.

## 3. Gravitons: what a spinning analogue requires

* **Structure [Q + D-sketch].** CHMRS (app. superconvergence) show that `C₂ ∝ Φ_{−2,ν}[G/(2⁸u′v′)]` is the null-integrated commutator
  `Φ_{−1,ν}[(v′−u′)G/(u′v′)]` followed by the `ν`-finite-difference inversion of multiplication by `η`:
  * `(v′−u′)G̃_{−1,ν} = c₊G̃_{−2,ν+i} − c₋G̃_{−2,ν−i}`;
  * the inversion is invertible because `η ≥ 1`.

  This is exactly the device that resolves `SUPERCONVERGENT_UPLIFT.md`'s objection to meromorphic even weights. Any antisymmetric meromorphic
  weight gives a leading kernel `∝ η·P_ν`, and the `η`-inversion restores a pure harmonic.
* **Proposed graviton definitions.**
  * **Odd `k`:** keep the meromorphic-weight bare commutators already in use: `A₃ = 1` (CLS-D's `k = 3`), `A₅ = (u′+v′)/(32u′v′)`,
    `A₇ = 1/(256u′v′)`. Fix the normalisation by the **exact** leading coefficient rather than the `d₀`/`j(ν)` convention.
  * **Even `k`:** replace the non-meromorphic `sgn(v′−u′)r^{k−3}` by an antisymmetric meromorphic weight, for example
    `(v′−u′)(u′+v′)^{b}/(2⁸u′v′)^{c}`, followed by the `η`-inversion. For identical scalars these reduce to `C₂`, `C̃_k` and the `C₄` combination of
    `cfoursuper`.
* **Properties to establish.**
  * (ii) Einstein-free given RB₂: by the counting rule (`LIGHT_INPUTS.md` item 1), *provided* the meromorphic weights do not spoil convergence at the
    single-ray endpoints [I; the same open point as M2 of `SUPERCONVERGENT_UPLIFT.md`].
  * (v) (P): meromorphic weights restore the `dDisc` rewriting and the double-zero pattern. That is an improvement over the current even `Ψ_k`,
    whose (P) content is unanalysable.
  * (i) the bulk-point density of the spinning functional: it must be derived and compared with the B′ density, as for scalars.
* **Exact near-gap evaluation (iv) for gravitons** needs either:
  * spinning Mellin amplitudes / spinning Mack polynomials for heavy `⟨TTTT⟩` exchanges in `d = 3`; or
  * the CPPR reduction in a covariant tensor basis whose structures are monomials in `u′, v′`. The frame factors `r^c`, `η^j` are what break
    `C`'s Mellin form (`SPINNING_MELLIN.md` §5.2).

  Neither exists in the repo; both are research-level.

## 4. Scalar-skeleton results for `Ω^{NF,C}` (`r1_redef.py`) [NE + M for odd `k`]

`Ω^{NF,C} := ∫dp Σ_k φ_k(p)ε_kΛ^{2k}C_{k,Λp}`: NF02 weights, no forward terms, `C_{k,ν}` by definition (even `k` exact). Odd `k` is taken in two
readings:
* **B′:** odd B′ density with the continued scalar remainder (= exactskel e7);
* **E:** `ε_k` × the exact `k`-continued scalar `C`.

All e1 states are included: both skeletons, `J ≤ 10`, and the high-`J` samples.

| `Δ_g` | `n` (B′ / E; `δ = 0` and `0.03`) | min `F/N_s` (B′ / E) | worst state | `√(t₄²+α_T²)` |
|---|---|---|---|---|
| 100 | 1 / 1 | 0.2267 / 0.2281 | `+−`, `1.05Λ`, `J = 5` | 2.21·10⁻⁴ |
| 300 | 1 / 1 | 0.2057 / 0.2077 | `+−`, `1.02Λ`, `J = 5` | 3.13·10⁻⁶ |
| 1000 | 1 / 1 | 0.2022 / 0.2043 | `+−`, `1.01Λ`, `J = 5` | 2.89·10⁻⁸ |

These coincide with `EXACT_SKELETON_NEARGAP.md` §8 and `RESULT_D.md` (by construction for even `k`). **The redefinition changes the status of
these numbers, not their values.** On the scalar skeleton they are now the exact action of the defined functional, rather than a model of an
undefined one.

## 5. Grades and what remains

| statement | grade |
|---|---|
| CHMRS `C_{k,ν}` (even `k`) satisfies (i)–(v) at scalar level | [Q + D] |
| with `Ĉ_k :≡ C_{k,ν}`, the D-safe near-gap action on the scalar skeleton is exact: `n = 1` at `Δ_g = 100/300/1000`, bounds unchanged | [NE] (exactskel numerics; even `k` exact, odd `k` [M] bracketed) |
| the flat limit for even `k` is derived (bulk-point density `E_k`), not assumed | [Q] |
| odd `k` at scalar level | [M] (H3), readings differ `≤ 0.002 N_s` at the minimum |
| a spinning analogue with (i)–(v) exists (meromorphic light-ray commutators + `η`-inversion) | **H14′ [I]**: construction sketched (§3); endpoint convergence, normalisation, spinning bulk-point density, and (P) are to be shown |
| exact near-gap graviton values for the redefined functional | not available; needs spinning Mellin (research-level) |

**Net effect on the hypothesis list.**
* H14 for scalar skeletons and even `k`: **removed by definition.**
* H14 for gravitons becomes H14′ (existence and properties of the spinning `C`-type functionals) plus the transfer question H4 (graviton moments of
  that functional against the scalar `C` skeleton).
* H3 (odd `k`) is unchanged.
* The light side of D-safe changes only if the even-`k` graviton terms are actually switched from `Ψ_k/N^{ff}_k` to the meromorphic construction.
  Then `A`, `B_AdS` are unchanged at leading order (same leading kernels and same flat limit), and the (P) analysis becomes possible [D-sketch].

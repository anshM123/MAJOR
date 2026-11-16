# Referee report on `ads4/LIGHT_INPUTS.md` and `ads4/code/lightinputs/`

Referee: independent, hostile. Date 16 Nov 2026. Scratch code: `ads4/code/referee_lightinputs/` (`r1_energy_matrix.py`,
`r2_clsd_crystal.py`; logs in `cache/`, git-ignored). No file under review was modified. `run_all.sh` runs (3.8 s); all four logs
agree with the text.

## Verdicts

* **Correctness:** correct. I reproduced `t2` independently, with no use of the KKSDZ tensor (`r1`):
  * I polarised the Hofman–Maldacena one-point form in `d = 3`. The `t₂` term drops out identically, the diagonal is 1, and `|flip| = |t₄|/4 ∝ w^{∓4}`.
  * I classified every SO(2)-covariant quartic form, including the parity-odd ones built with `ñ`. Covariance gives a constant diagonal and an off-diagonal `∝ e^{−4iφ}`. The energy sum rule then forces `H = 𝟙 + flip`, with `c` complex. So no other structure exists, and `α_T` can only enter as `Im c`.
  * `[H₁,H₃] = −2i|c|² sin4φ₁₃ σ_z` holds to `10⁻¹⁴`.

  I also re-derived CLS-D's crystal values from scratch (`r2`). Solving the `f⁺` lines gives `u′ = m⁴/ν⁴`, `v′ = (m²−ν²)²/ν⁴`. All four cases of
  (106)/(107) then reproduce (108), and the Regge shift is `n = k−2J_e+1` for both parities. So `J₀ < k` holds for every `J_e`, `k`.
* **Soundness:**
  * Claim 1 is sound, given the CHMRS App. E hypotheses. The CLS-D slip is real, and the cleanest proof of it is internal to CLS-D. Their spinning
    section states that the bare commutator converges for `J₁+J₃−1 > J₀`, and their (107) gives `f⁻_{3;2} = const`. Together these give `3 > J₀`,
    whereas their text gives `k−1 > J₀`, i.e. `2 > J₀`. This argument does not depend on the crystal↔Regge map.
  * Claim 2 is sound. I checked one block analytically: `Ĥ₁₃ = −2(z₁·z₃+2)`, which is independent of `r`.
  * Claim 3 is sound, modulo (I1) (M1).
  * Claim 4 is conditional and honestly labelled.
  * Claim 5 is a sketch and is over-labelled (M3).
* **Labels:** mostly honest. The exceptions are M3 and m1.
* **Applicability:** the `C⁻₃` zero (leverage `~10⁴`) rests entirely on (I1). (I1) conflicts with `ADS4_LIGHT.md` §0 item 6 and §3 (M1).
* **Novelty:** low. Superconvergence killing of `J < k` exchanges is standard (KKSDZ, CLS-D, CHMRS). Claim 3 is KKSDZ's planar sum rule (M2). What is new:
  * the CLS-D off-by-one (an erratum to a conservative statement, which harms nothing in CLS-D);
  * the explicit `d = 3`, parity-odd, two-line commutator;
  * the bookkeeping.

## Issues

**Fatal:** none.

**M1 (major: internal conflict, applicability).** §3 says the `k = 3` functional has double zeros on *all* `[TT]` families. `ADS4_LIGHT.md`
§0 item 6 and §3 say that `⌈k/2⌉ = 2` families enter `C⁻₃` at `O(1/c_T)` (labelled [I] there, extrapolated from the scalar CHMRS `C_k`). Both cannot hold.
If `ADS4_LIGHT` is right about the project's `C⁻₃`, then `C⁻₃` is not the bare commutator, and (I1) fails.

The double-zero claim *is* correct for the bare commutator, and a better argument exists than the CLS-D quote. In MFT `[T,T]` is a c-number.
So `⟨0|[T₁,T₂]|[TT]⟩ = O(γ)`, and the double commutator in `Ψ⁺` is `O(γ²)` for every family and every structure.
*Fix:* state the conflict, use the factorisation argument, and carry out the `C⁻₃` → CLS-D `k = 3` mode map. This is the only remaining link.

**M2 (major: attribution/novelty).** KKSDZ already derive the content of item 3 in general `d`:
* the `T` contribution to `⟨T|[E,E]|T⟩` is `∝ ((4+d)t₂+4t₄)²`, with no Einstein term and no linear term (their App., the `{2,2|2|2,2}` sum rules; (superconvergenceexamples));
* "satisfied if all single traces have minimal couplings";
* double traces are `O(1/N⁴)`.

The `d = 3` matrix is HM/CDP. *Fix:* cite these, and present item 3 as their `d = 3` specialisation with `α_T` included.

**M3 (major: label).** "`(ν²+9/4)⁻²` excluded **[D]**" (§0.5, §4.3, table) assumes two things:
* any IP failure is a re-weighting of the same transverse plane waves (the `M`-table ratios);
* the density is a polynomial differential operator × `Π_{Δ₀=2}` (the MP form).

IP failure means longitudinal pieces, whose tables were never computed. *Fix:* relabel it as [D given MP form + transverse class]. Also relabel the `γ`
cancellation (a residue-in-`k` sketch) as [D-sketch].

**m1.** CLS-D (111) [(3.80)] prints `f⁻_{3;2} = +1`. The value `−1` comes from (107)/(108) [(3.76)/(3.77)]. Quote it correctly and give the arXiv numbering.
**m2.** "No Regge assumption" is not quite true. The `k = 3` sum rule needs `J₀^planar < 3` for the full UV-complete correlator (KKSDZ via the chaos bound). Say so.
**m3.** Item 6: KKSDZ's Einstein–Maxwell `H_ij` governs `⟨J|EE|J⟩ = ⟨JTTJ⟩`, not `⟨JJJJ⟩`. For `J_e = 1`, the `k = 3` functional is subtracted (`n = 2`).
**m4.** Lemma F: KKSDZ's `J₀` frame (`x₂ = −z₂/α₂`, …) is not the CHMRS frame, so one sentence on the `r`-independent map is needed. Also, frame modes are *bounded* by `r^{1−J₀}`, not equal to it.
**m5.** The `T`-block route implies that no contact enters `C⁻₃` at `O(1/c_T)`. Reconcile this with the `Σ_{k′≥k}` tower of §2.2, which as written lets `J ≥ 3` contacts in. This bears on (P).
**m6.** RB₃: KKSDZ only assert `J₀^planar = 2`. Keep it as [I] for the Einstein truncation, or compute it from the explicit AdS4 graviton results.

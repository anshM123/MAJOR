# Double-lightcone estimates in the rigorous lightcone-bootstrap literature, and what they imply for the target of `TARGET.md`

Reader's report, 15 November 2026. Nothing here is refereed.

**Marking convention.** Every statement is one of

* **[Q]** quoted or transcribed from the paper named, with its own numbering (formulas are transcribed from `pdftotext` output of the arXiv PDF and re-typeset; wording inside quotation marks is verbatim);
* **[I]** an immediate inference from quoted statements (the step is given);
* **[A]** my own assessment or my own derivation, done in this session, checked only by me. Sections 5–7 are almost entirely [A].

**Sources actually read (full text, arXiv PDF → `pdftotext -layout`).**

| Paper | Version read | How much |
|---|---|---|
| Pal, Qiao, Rychkov, "Twist accumulation in conformal field theory. A rigorous approach to the lightcone bootstrap", arXiv:2212.04893 (CMP 402 (2023) 2169) — "PQR" | v3, 17 May 2023 | §§1–2, 4, 6, App. A in full; §3, 5, App. B skimmed |
| van Rees, "Theorems for the Lightcone Bootstrap", arXiv:2412.06907 (SciPost Phys. 18 (2025) 207) — "vR" | v1, 9 Dec 2024 (8 pages) | in full |
| Pal, Qiao, van Rees, "Universality of the microcanonical entropy at large spin", arXiv:2505.02897 — "PQvR" | v1 | intro, §3.2, outlook |
| Pal, Qiao, "Lightcone Modular Bootstrap and Tauberian Theory", arXiv:2307.02587 | v2 | abstract, intro, definitions of limits |
| Pal, Rogelberg, "Liberation in the Heavy Limit", arXiv:2611.11594 | v1, 8 Nov 2026 | abstract, intro, §6.2 (relay lemmas) |
| Li, Meltzer, Poland, "Conformal Collider Physics from the Lightcone Bootstrap", arXiv:1511.08025 — "LMP" | — | abstract, §§1, 3.2, 4.2, 5.2, 8 (passages) |
| Hofman, Li, Meltzer, Poland, Rejon-Barrera, arXiv:1603.03771; Caron-Huot, arXiv:1703.00278; Fitzpatrick–Kaplan–Poland–Simmons-Duffin, arXiv:1212.3616 | — | abstract + keyword-located passages only |
| Kundu arXiv:2002.12390; Komargodski–Zhiboedov arXiv:1212.4103; Faulkner et al. arXiv:1605.08072; Hartman–Kundu–Tajdini arXiv:1610.05308 | — | **abstract only** |

**Follow-up search.** INSPIRE `refersto` lists of both papers (23 and 51 citing records on 15 Nov 2026, all titles read), plus arXiv API searches by author (Pal, Qiao, van Rees) and by abstract keywords ("lightcone bootstrap" with rigorous/theorem/Tauberian; "twist accumulation"; Tauberian with "large spin"/"anomalous dimension"/"double-twist"). **Result: I found no paper, 2023–Nov 2026, proving a quantitative correction (rate of approach, anomalous dimension) for the four-point function in d ≥ 3.** The rigorous follow-ups are all two-dimensional modular results (2307.02587, 2410.18174, 2505.02897) or large-Δ results at exponential accuracy (2611.11594). "Not found" means not found by these searches.

**Correction to the brief.** arXiv:1811.05710 is Collier–Gobeil–Maxfield–Perlmutter, "Quantum Regge Trajectories and the Virasoro Analytic Bootstrap" (2D, Virasoro fusion kernel; a derivation), not Mukhametzhanov–Zhiboedov. The Mukhametzhanov–Zhiboedov Tauberian papers are arXiv:1808.03212 ("Analytic Euclidean Bootstrap", cited by vR as his ref. [22]) and arXiv:1904.06359 (from memory, not re-checked); both concern the Euclidean/OPE-in-dimension limit, not the lightcone limit.

Notation of this file: `y := 1 − z̄`, `G̃(z,y) := G(z, 1−y)`, `t := (τ_*+δ)/2`, `ε := z^{τ_*/2}`, `K(z,y) := (z(1−y)/((1−z)y))^Δ`. Crossing reads `G̃(z,y) = K(z,y) G̃(y,z)`.

---

## 1. (a) Hypotheses and main theorems

### 1.1 PQR

**Setting [Q, §2].** Unitary CFT in `d ⩾ 2`; `g(z,z̄)` the four-point function of a Hermitian scalar quasiprimary of dimension `Δ_φ`; crossing (2.1) `g(z,z̄) = (z z̄/((1−z)(1−z̄)))^{Δ_φ} g(1−z,1−z̄)`; expansion (2.2) `g = 1 + Σ p_{τ,ℓ} g_{τ,ℓ}`, `p ⩾ 0`; block expansion (2.3) `g_{τ,ℓ}(z,z̄) = (z z̄)^{τ/2} Σ_{n,m⩾0} a_{n,m} z^n z̄^m`, `a_{n,m} ⩾ 0`; twist gap (2.4) `τ ⩾ τ_gap > 0` for all non-identity blocks. All arguments are for **real** `0 < z, z̄ < 1`. Block normalisation (2.26): `lim_{x→0+} x^{−Δ} g_{τ,ℓ}(x,x) = 1` (this differs from the normalisation of `TARGET.md` by the factor in (2.28) below).

**Theorem 2.1 [Q].** "Let φ be a Hermitian scalar quasiprimary in a unitary CFT in d ⩾ 2. Suppose that • all nontrivial quasiprimary operators in the OPE φ × φ satisfy the twist gap assumption (2.4) with τ_gap > 0; • the spectrum of operators in this OPE is discrete, i.e. there is a finite number of operators with Δ ⩽ Δ_* for any Δ_*. Then: (a) τ_gap ⩽ 2Δ_φ; (b) The OPE φ × φ contains an infinite sequence of quasiprimary operators (O_i) whose twist τ_i → 2Δ_φ while spin ℓ_i → ∞."

No restriction on `Δ_φ` beyond unitarity, no restriction on `d ⩾ 2`, footnote 3: `τ_gap` need not be the maximal gap. §2.3.1 removes the discreteness assumption (Proposition 2.13, Lemma 2.14). §2.3.2 only *sketches* non-identical scalars ("We did not try to work out a proof of (2.60)") and spinning operators.

**Theorem 3.1 [Q]** is the 2D modular analogue (`c > 1`, normalisable vacuum, twist gap in Virasoro primaries ⇒ `τ_gap ⩽ (c−1)/12` and a family with `h → (c−1)/24`, `h̄ → ∞`). Not used below.

**A dependency the authors flag [Q, App. A.2].** The block bound Lemma 2.7 is proved in `d ⩾ 3` from Hogervorst's dimensional-reduction coefficients (A.13), about which: "The complete proof that (A.13) is a solution of Hogervorst's recursion relations is so far missing, but there is hardly any doubt that they are correct." So, strictly, the `d ⩾ 3` low-twist step of Theorem 2.1(b) rests on that closed formula. (The direct-channel estimates of §2 below do **not**.)

### 1.2 van Rees

**Setting [Q, §2].** Functions (1) `G(z,z̄) = z^{−Δ} z̄^{−Δ} Σ_k Σ_{j=0}^{l_k} c_{k,j} z^{(E_k−l_k)/2+j} z̄^{(E_k+l_k)/2−j}`, `c_{k,j} = c_{k,l_k−j} ⩾ 0`, `E_k ⩾ 0`, `Δ > 0`, absolutely convergent for `0 < |z|,|z̄| < 1`; discrete (`E_k → ∞`); identity `c_{0,0} = 1`; twist gap `inf_{k>0}(E_k − l_k) > 0`; crossing (2) `G(z,z̄) = G(1−z,1−z̄)` (the prefactor is absorbed in `G`). No dimension `d` enters the main text; `d` appears only in the block statement of App. B.

**Main results [Q].** With `z = e^w`, `g(w,z̄)` the continuation to `Re w < 0`, `g = Σ_n f_n(w) z̄^{h_n−Δ}` (3) and `F_Λ(w) := Λ^{−Δ}Γ(Δ+1) Σ_{h_n ⩽ Λ} f_n(w)` (5):

* Proposition 3.2: "Pointwise for any real w < 0, lim_{Λ→∞} F_Λ(w) = (1 − e^w)^{−Δ}." (Hardy–Littlewood, Theorem A.1.)
* Theorem 3.3: "For any w with Re(w) < 0, lim_{Λ→∞} F_Λ(w) = (1 − e^w)^{−Δ}, with the limit being uniform on compact subsets." (Vitali.)
* (13): for `φ ∈ D(ℝ)`, `lim Λ^{−Δ} ∫ φ H_Λ = Σ_m Γ(m+Δ)/(m! Γ(Δ)Γ(Δ+1)) φ(Δ+m)`: "as an element of D′(ℝ), the monomial twist density becomes that of the generalized free field to leading order at large Λ."
* (14),(16): `lim_{z̄↗1}(1−z̄)^Δ g(w,z̄) = (1−e^w)^{−Δ}` uniformly on compacts of `Re w < 0`, hence `lim (1−z̄)^Δ dDisc_s g = 0`.
* Theorem B.2 / (51): the same for the primary twist density, with the GFF coefficients; existence of infinitely many primaries with `h ∈ [Δ+k−ε, Δ+k+ε]` "If Δ > d/2 − 1".

All limits are `z̄ → 1` (equivalently `Λ → ∞`) **at fixed `w`** or uniformly on compact sets of `w`. No rate.

---

## 2. (b) The double-lightcone limit they use, and every remainder estimate

### 2.1 Definition

**PQR (2.12) [Q]:** "DLC∗ limit: z → 0, z̄ → 1, z^{τ_gap/2}/(1 − z̄)^{Δ_φ} → 0. The third requirement is what makes DLC∗ different from DLC." In the notation of this file: `y ≫ z^{τ_gap/(2Δ)}`. Figure 1 caption: "DLC*: z = o((1 − z̄)^{2Δ_φ/τ_min})". For non-identical scalars the sketched analogue is (2.59) `z^{τ_gap}/(1−z̄)^{Δ_1+Δ_2} → 0`.

**vR:** uses **no** double-lightcone limit. He says so [Q, §1]: previous approaches "proceeded to series expand around z = 0, but this is not justified because the Tauberian theorem only works pointwise for 0 < z < 1"; his limit is `z̄ → 1` at fixed `z = e^w`.

### 2.2 PQR estimates, complete list for §2 (all [Q]; ranges as printed)

Direct channel.

1. **Lemma 2.2, (2.16).** `g_{τ,ℓ}(z,z̄) ⩽ (z/z_0)^{τ/2} g_{τ,ℓ}(z_0,z̄)` for `0 ⩽ z ⩽ z_0 < 1`, `0 ⩽ z̄ < 1`.
2. **Lemma 2.3, (2.18).** `g(z_0,z̄) ⩽ B/(1−z̄)^{Δ_φ}` for all `0 ⩽ z̄ < 1`; Remark 2.4: `B = max(g(z_0,z_0), (z_0/(1−z_0))^{Δ_φ} g(1−z_0,1−z_0))`.
3. **Proposition 2.5, (2.22).** With `g_{τ⩾τ_*}` the sum of all direct-channel terms of twist `⩾ τ_*`: `g_{τ⩾τ_*}(z,z̄) ⩽ C z^{τ_*/2}/(1−z̄)^{Δ_φ}` for `0 ⩽ z ⩽ z_0 ⩽ z̄ < 1`, `C = B/z_0^{τ_*/2}`.
4. **(2.14).** `1 ⩽ g(z,z̄) ⩽ 1 + C z^{τ_gap/2}/(1−z̄)^{Δ_φ}` on the same range; hence (2.15) `g → 1` when `z^{τ_gap/2}/(1−z̄)^{Δ_φ} → 0`.
5. **Remark 2.6.** "A weaker version of Eq. (2.14), with z^{τ_gap/2}/(1−z̄)^{Δ_φ} replaced by z^{τ_gap/2}/(1−z̄)^{2Δ_φ}, was previously shown in [12], App. F. It would be interesting to know if Eq. (2.14) is optimal or if it can be further improved."

Blocks.

6. **Lemma 2.7, (2.29).** With `F_{τ,ℓ}(z) := ((ν)_ℓ/(2ν)_ℓ) 4^{τ/2} k_{τ+2ℓ}(z)` (2.28), `ν = (d−2)/2`, `ρ̄ = z̄/(1+√(1−z̄))²`: `1 ⩽ g_{τ,ℓ}(z,z̄)/(ρ̄^{τ/2} F_{τ,ℓ}(z)) ⩽ K_d(a,b)(1 + θ(d⩾3)δ_{ℓ0}/(Δ−ν))` for `0 ⩽ z̄ ⩽ a < b ⩽ z < 1`, `K_d(a,b)` "independent of τ, ℓ". Equivalent form Lemma A.1, (A.1), with denominator `((ν)_ℓ/(2ν)_ℓ) k_τ(z̄) k_{τ+2ℓ}(z)`. **Remark A.2:** "it's easy to see that K_d(a,b) = 1 + O(a) for a → 0 and b fixed."
7. **Remark 2.8(e), (2.33).** `F_{τ,ℓ}(z) ∼ ((ν)_ℓ/(2ν)_ℓ) 4^{τ/2} C_{τ+2ℓ} log(1/(1−z))` as `z → 1`; footnote 8: "Surprisingly, we do not know any other general rigorous proof of this result."
8. **Lemma 2.14, (2.51).** `g_{τ,ℓ}(z,z̄)/g_{τ,ℓ}(b,a) ⩽ B z̄^{τ/2} log(1/(1−z))` for `0 ⩽ z̄ ⩽ a < b ⩽ z < 1`, `τ ⩽ τ_max`, `ℓ ⩽ ℓ_max`.

Crossed channel.

9. **Lemma 2.9, (2.34)–(2.37).** Identity: `(z/(1−z̄))^{Δ_φ} ⩽ z^{τ_gap/2}/(1−z̄)^{Δ_φ}`; a single block: `g_{τ,ℓ}(1−z,1−z̄) ⩽ C (1−z̄)^{τ_gap/2} log(1/z)`, so its contribution is `⩽ C [z^{τ_gap/2}/(1−z̄)^{Δ_φ}]^κ z^{Δ_φ−τ_gap κ/2} log(1/z)`, `κ = 1 − τ_gap/(2Δ_φ)`.
10. **High twists, (2.39)–(2.40).** `g_{τ⩾τ_*}(1−z,1−z̄) ⩽ C (1−z̄)^{τ_*/2}/z^{Δ_φ}` for `0 < z ⩽ z_0 ⩽ z̄ < 1`; for `τ_* = 2Δ_φ+ε` the contribution to crossing is `⩽ const·(1−z̄)^{ε/2}`. "(We only need z̄ → 1 but not the full strength of DLC∗.)"
11. **Low twists, (2.42)–(2.49).** From the two-sided bound (2.42): (2.44) `g_{τ,ℓ}(z,z̄) ⩽ C (z̄/z̄′)^{τ_*/2} g_{τ,ℓ}(z,z̄′)` for `0 < z̄′ ⩽ z̄ ⩽ a < b ⩽ z < 1`, `τ ⩽ τ_*`; hence (2.47): for `τ_* = 2Δ_φ − ε` the crossed low-twist contribution is `⩽ C ((1−z̄′)/(1−z̄))^{ε/2} g(z,z̄′)` for `1−a ⩽ z̄ ⩽ z̄′ < 1`; with the choice (2.48) `z̄′ = 1 − z^{τ_gap/(2Δ_φ)}` this is `⩽ C′ (z^{τ_gap/(2Δ_φ)}/(1−z̄))^{ε/2}`. Remark 2.11: "The nontrivial equation was (2.44)."
12. **Proposition 2.13, (2.53)–(2.54).** Blocks with `τ ⩽ τ_max`, `ℓ ⩽ ℓ_max` contribute `⩽ B′ z^{Δ_φ}(1−z̄)^{τ_gap/2−Δ_φ} log(1/z)`.
13. **(2.58), sketched only.** `g_{1122} = 1 + O(z^{τ_gap/2}/(1−z̄)^{(Δ_1+Δ_2)/2})`.

### 2.3 Do any of these subtract image terms?

**No [I].** Every direct-channel estimate (items 1–5) is an *upper bound on a sum of positive terms* by the full correlator at a reference point; the only identity used on the direct side is crossing at `z = z_0` inside Lemma 2.3. The `(z/y)^Δ` behaviour (image of the crossed-channel identity) is exactly what Lemma 2.3 feeds into Proposition 2.5 as a bound; it is never isolated or removed. The crossed-channel items 9–12 bound the identity and low-twist terms *from above*; none compares them with the corresponding part of the direct channel. The words "image", "subtract", "second order" do not occur. vR has no remainder estimate at all in a double limit: his only two-variable statements are (7) (the crossed identity term `(1−z)^{−Δ}(1−z̄)^{−Δ}` plus a correction of relative order `(1−z̄)^{ĥ}`, `ĥ = inf_{n>0} h_n`, which "holds pointwise for any z as long as 0 < |1 − z| < 1") and the uniform-on-compacts limits (9), (14).

PQvR (2D, modular) is the only place where a bound is obtained by using the symmetry *twice*: [Q, intro] "A crucial difference appears in the proof of lemmas 3.4 and 3.5, where we use modular invariance twice to bound the non-universal term"; [Q, after (3.39)] "The trick is now to use modular invariance again". For four-point functions: [Q] "The application of the ideas presented in this work to the four-point function bootstrap will appear elsewhere." I found no such paper. Pal–Rogelberg 2611.11594 iterate a twist inequality and crossing along an orbit ("relay", their Lemmas 6.10–6.11, taken from Dey–Pal–Qiao 2410.18174), but at exponential accuracy in `Δ → ∞` at fixed cross-ratios; their abstract records that for twist cutoffs below `Δ` "composites of two spinning generalized free fields fail to be liberated in a region that fills the band near the double lightcone limit".

---

## 3. (c) Do their estimates give the direct channel to `o(z^{τ_*/2})` in a window `y ⩽ z^b`, `b > τ_*/2`?

**As stated: no, except under `δ > Δτ_*`.**

* PQR's own direct-channel statement (2.14) uses only `τ_gap` and gives `G̃ = 1 + O(z^{τ_gap/2} y^{−Δ})`: the first non-identity term is inside the error.
* [I] Under hypothesis (H) of `TARGET.md`, Proposition 2.5 with `τ_* → τ_*+δ` gives `R(z,y) ⩽ C z^{(τ_*+δ)/2} y^{−Δ}`; the finitely many twist-`τ_*` blocks are `z^{τ_*/2}[k_{τ_*+2ℓ}(1−y) + O(z y^{−Δ})]` (positivity of the `z`-expansion plus Lemma 2.3; or `O(z log(1/y))` using Lemma 2.14). So
  `G̃(z,y) = 1 + z^{τ_*/2} Σ_* p_O k_{τ_*+2ℓ_O}(1−y) + o(z^{τ_*/2})` uniformly for `y ⩾ z^{a}`, any fixed `a < min(δ,2)/(2Δ)`.
  This is exactly the crude bound of `TARGET.md` §4; it meets `y ⩽ z^b`, `b > τ_*/2`, iff `δ > Δτ_*`.
* vR's estimates give nothing here: no uniformity as `w → −∞` jointly with `z̄ → 1`.

**The step that would have to be strengthened** is Lemma 2.3 → Proposition 2.5: the whole high-twist sum at `z = z_0` is bounded by the full correlator, `g(z_0,z̄) ⩽ B(1−z̄)^{−Δ}`. That singularity belongs to direct twists near `2Δ` (the image of the crossed identity), but the proof charges it to *every* twist above the gap, at the rate `(z/z_0)^{τ_*/2}` of the lowest one. Remark 2.6 (item 5) is the authors' acknowledgement that the exponent pair in (2.14) may not be optimal. Section 5.2 below shows [A] that iterating Proposition 2.5 with crossing improves `y^{−Δ}` to `y^{−max(t,Δ−t)−0}`, and that positivity alone cannot do better.

---

## 4. (d) What they say about corrections and open problems (all [Q])

* PQR, Remark 2.6: quoted in item 5 above.
* PQR §6: "We note that our analysis did not take advantage of the analyticity of CFT 4pt functions and their conformal blocks, which is essential for a stronger result of the lightcone bootstrap: the OPE of φ × φ also contains a tower of families of quasiprimaries O_{τ,ℓ} with ℓ → +∞ and τ → 2Δ_φ + 2n." … "It would be nice to provide a mathematical proof of the MFT universality at large spin." … "The lightcone bootstrap also leads to predictions about the OPE spectral density of the operators whose twist asymptotes to 2Δ_φ + 2n [1, 2]. It would be interesting to prove those results rigorously, using e.g. such techniques such as the Tauberian theory [9, 12, 32]." Anomalous dimensions are not mentioned as proved or as a stated goal.
* PQR footnote 16 (a warning relevant to any first-moment statement): "we may not differentiate after estimate."
* vR §5: "It will be interesting to understand the size of the corrections to the results presented here. The best estimates will likely utilize our knowledge of the complex z plane and rely on one of the Tauberian theorems with remainder proven by Subhankulov [21]. We imagine this will allow us to estimate corrections for the behavior of F_Λ(w) and H_Λ(η) as described in theorem A.2 and in equation (13)."
* vR §4: "This simpler proof however only works at the leading order, and for subleading terms the use of the Abelian theorem seems unavoidable." and "In the future we plan to use this result as a stepping stone towards proving several of the claims in the works mentioned in the introduction."
* PQvR: "A rigorous analysis of the off-diagonal or lightcone limit proved to be more stubborn"; "Proving theorems for the lightcone bootstrap at higher twists is more difficult because of the need to perform two inverse Laplace transforms, which naively invalidates the direct use of Tauberian theorems."
* Pal–Qiao 2307.02587, footnote 5: "in the double lightcone limit, the conformal blocks exhibit approximate factorization (see [28], appendix A). Based on this observation, we anticipate that a similar Tauberian theorem to theorem 2.1 can be established".
* Pal–Rogelberg 2611.11594 describe large-spin liberation as "discovered in [1, 2] and partly proven in [3, 4]" ([3,4] = PQR, vR).

---

## 5. (e) Assessment [A]

Everything in this section is my own analysis. I have checked the algebra of §5.4 numerically on a synthetic measure (script in the session scratchpad, not in the repository); I have **not** tested the full statement on an interacting correlator.

### 5.1 Window arithmetic of `TARGET.md` §4, redone for general `d`

Crude route: window `y = z^β` with `τ_*/2 < β < δ/(2Δ)`, i.e. `δ > Δτ_*`; (H) forces `δ ⩽ 2Δ − τ_*`. Hence `τ_* < 2` and `Δ > τ_*/(2−τ_*)` are necessary.

Stress tensor, `τ_* = d−2`; if the double-stress-tensor tower appears in `φ×φ` then `δ ⩽ d−2` (its twist tends to `2(d−2)`, expected from below).

| `d` | crude window non-empty iff | result |
|---|---|---|
| general | `(d−2)Δ < δ ⩽ min(2Δ−(d−2), d−2)` | needs `(d−2)/(4−d) < Δ < 1`, possible only for `d < 3` |
| 3 | `Δ < δ ⩽ min(2Δ−1, 1)` | **empty** (needs `Δ > 1` and `Δ < 1`); exactly marginal |
| ⩾ 4 | `δ > (d−2)Δ ⩾ 2Δ > 2Δ−(d−2) ⩾ δ` | **empty for every `Δ, δ`, even without `[TT]`** |

For `d ⩾ 4` the reason is structural: `y ≪ z^{τ_*/2} ⩽ z` lies beyond the anti-diagonal `y = z`, where `(z/y)^Δ ⩾ 1` and the direct channel is not near 1 at all. So the arithmetic of `TARGET.md` is confirmed, and it is worse in higher `d`. For a scalar exchange with `τ_* < 2` the crude window is non-empty iff `Δτ_* < δ ⩽ 2Δ − τ_*`.

### 5.2 What iterating between the two channels gives, and where it stops

**Proposition (relay improvement of PQR Prop. 2.5) [A].** Assume (H). For every `e > e_∞ := max(t, Δ−t)` there is `C_e` with `R(z,y) ⩽ C_e z^{t} y^{−e}` for `0 < z ⩽ y ⩽ y_0`.

*Proof sketch.* Lemma 2.2 summed over twists `⩾ 2t`, then crossing: for `z ⩽ z′`,
`R(z,y) ⩽ (z/z′)^t G̃(z′,y) = (z/z′)^t (z′/y)^Δ ((1−y)/(1−z′))^Δ G̃(y,z′)`.  (★)
Suppose `R(z,y) ⩽ C_k z^t y^{−e_k}` on `z ⩽ y` (true for `e_0 = Δ` by Prop. 2.5). Put `z′ = y^{t/e_k} ⩾ y`. Then at the point `(y,z′)`: `R(y,z′) ⩽ C_k`, the twist-`τ_*` part is `O(y^{τ_*/2} log(1/z′))` (Lemma 2.14), so `G̃(y,z′) ⩽ C`. (★) gives `R(z,y) ⩽ C′ z^t y^{−e_{k+1}}`, `e_{k+1} = Δ − t(Δ−t)/e_k`. The map has fixed points `t` and `Δ−t`; from `e_0 = Δ` it decreases to the larger one. ∎

Consequences.
* It answers Remark 2.6 partially: (2.14) is not optimal when there is a gap structure; the exponent of `1/(1−z̄)` can be lowered from `Δ` to `max(t, Δ−t)+0`.
* It matches the heuristic "image" picture: in a symmetric double expansion of `(y(1−z))^Δ G̃` in monomials `z^a y^b`, a tower of twist `2t` in both channels gives a term `z^t y^{t−Δ}` in `G̃`, and the identity image restricted by Lemma 2.2 gives `(z/y)^t`.
* Window with this bound: `τ_*/2 < β < δ/max(τ_*+δ, 2Δ−τ_*−δ)`, i.e. `δ > max(τ_*²/(2−τ_*), τ_*(2Δ−τ_*)/(2+τ_*))`. Stress tensor with `δ ⩽ τ_*`: needs `τ_* < 1`, i.e. `d < 3`. **`d = 3` is again exactly marginal, `d ⩾ 4` empty.**
* Ceiling: any argument built from Lemma 2.2 and crossing bounds `R(z,y)` by `(z/y)^t R(y,y)`, and `R(y,y) ⩾ 1 − o(1)` because on the anti-diagonal the identity image equals 1. Going below `(z/y)^t` requires separating `R` by direct-channel twist near the anti-diagonal; monotonicity in `z` cannot do it (the "reverse" inequality (2.44) costs the maximal twist, the forward one gains the minimal twist), and in `y` no inequality uniform in spin exists (`k_{2h̄}(1−y)/k_{2h̄}(1−y′)` is unbounded in `h̄`). **So a genuinely subtracted second-order DLC estimate in a window `β > τ_*/2` is, in my judgement, out of reach of positivity-plus-crossing; it needs analyticity (vR's second sheet, i.e. a rigorous dDisc).**

### 5.3 Is the obstruction of `TARGET.md` §4 real?

**For the route written there, yes. For the target itself, I believe not.** The lower limit `β > τ_*/2` comes from treating three crossed-channel `O(y)` effects as errors: the factor `(1−y)^{−Δ}`, the `O(y)` corrections in the factorisation of crossed blocks, and the `n = 1` tower. All three can be handled *exactly*:

1. Crossing can be written as an **exact Laplace transform of a positive measure**: by (2.3), `G(1−z,y) = Σ_h A_h(z) y^h` with `A_h(z) ⩾ 0`, so for all `y ∈ (0,1)`
   `(1−y)^{−Δ} G̃(z,y) = ∫ y^{x/2} dM_z(x)`,  `M_z := (z/(1−z))^Δ Σ_h A_h(z) δ_{x = 2h−2Δ} ⩾ 0`.
   `M_z` is the crossed-channel *monomial* twist measure (the object of vR's (11)–(12)); its part from leading monomials of primaries is `(1−z)^{−Δ} μ_z`.
2. On the other side `(1−y)^{−Δ}·1 = Σ_m ((Δ)_m/m!) e^{−2ms}` is explicit, with exponents exactly at the GFF tower positions `x = 2m`.
3. Replace the three-point second difference of the Tauberian lemma by a `(2M+1)`-point functional whose test function has **double zeros at `x = 0, 2, …, 2(M−1)`**: `φ(x) = Π_{m<M}(e^{−h(x−2m)} − 1)² = Σ_{j⩽2M} c_j e^{−jhx} ⩾ 0`, and `ψ(x) = −h^{−1}(e^{−hx}−1) Π_{1⩽m<M}(e^{−h(x−2m)}−1)²/(e^{2mh}−1)²`, so `ψ(0)=0`, `ψ′(0)=1`, `ψ, ψ′` vanish at `2m`, `1 ⩽ m < M`. With `Λ_φ[f] := Σ_j c_j f(S+jh)`: `Λ_φ` annihilates `1`, `s`, `e^{−2ms}` (`m < M`); `Λ_ψ[s] = −1`. The towers `n ⩾ 1` then sit **inside** the positive measure at (near-)zeros of the test function and need no estimate. This is a real-variable, finite-order substitute for the double discontinuity, whose kernel `2 sin²(πx/2)` has the same zeros.

Then no lower limit on `β` remains, and the window can be placed deep inside the region where Proposition 2.5 is ample.

### 5.4 The statement that looks provable

**Theorem T (candidate) [A].** Let `d ⩾ 2`, unitary CFT, discrete spectrum, identical Hermitian scalars of dimension `Δ`, blocks normalised as in `TARGET.md` §2 with non-negative double power series (PQR (2.3)). Assume (H): besides the identity, finitely many operators in `φ×φ` have twist `< τ_*+δ`, all of twist `τ_* > 0`; `b_* := Σ p_O Γ(τ_*+2ℓ_O)/Γ(τ_*/2+ℓ_O)²`. Fix `0 < β < min(δ,2)/(2Δ)`, `0 < η ⩽ 1`, an integer `M > τ_*/(2β) + 1`, and put `y = z^β`. Then as `z → 0`

* (i) `∫_{|x|<η} x² y^{x/2} dM_z(x) ⩽ C ρ(z)`;
* (ii) `| ∫_{|x|<η} x y^{x/2} dM_z(x) + 2 b_* z^{τ_*/2} | ⩽ C ρ(z)`;
* (iii) `∫_{|x−2m|<η} (x−2m)² y^{x/2} dM_z ⩽ Cρ(z)` for `1 ⩽ m < M`, and `∫_{dist(x,{0,2,…,2M−2}) ⩾ η} y^{x/2} dM_z ⩽ Cρ(z)`;

with `ρ(z) = z^{τ_*/2}[ z^{δ/2−βΔ} + z^{1−βΔ} + z^{β} log(1/z) ] + z^{Mβ} = O(z^{τ_*/2+κ})`, `κ > 0`. No relation between `δ` and `Δτ_*`, no restriction on `d`, `Δ` or `c_T`.

*Proof sketch.* `L(s) := ∫ e^{−sx} dM_z = (1−y)^{−Δ}[1 + S_*(z,y) + R(z,y)]`, `y = e^{−2s}`, exactly. (1) `(1−y)^{−Δ} = Σ_{m<M} g_m e^{−2ms} + O(y^M)`. (2) `S_* = z^{τ_*/2} Σ_* p k_{τ_*+2ℓ}(1−y) + O(z^{τ_*/2+1} y^{−Δ})` (positive `z`-expansion, Lemma 2.2 with exponent `τ_*/2+1`, Lemma 2.3), and `k_{2h}(1−y) = (Γ(2h)/Γ(h)²)(2s − 2ψ(h) + 2ψ(1)) + O(y log(1/y))`, so `(1−y)^{−Δ}S_* = ε(c + 2b_* s) + O(ε y s) + O(ε z y^{−Δ})`. (3) `R ⩽ C z^t y^{−Δ}` (Prop. 2.5). Hence `L(s) = Σ_{m<M} g_m e^{−2ms} + ε(c+2b_*s) + r(s)`, `|r| ⩽ Cρ` at the `2M+1` points `S + jh`, `e^{−2S} = y`. Apply `Λ_φ`: `∫ e^{−Sx} φ dM_z = Λ_φ[r] ⩽ ‖c‖_1 Cρ`, which gives (i), (iii) since `φ ⩾ c_η x²` near 0, `⩾ c_η (x−2m)²` near `2m`, `⩾ c_η` elsewhere on `x ⩾ −2Δ`. Apply `Λ_ψ`: `∫ e^{−Sx} ψ dM_z = −2b_* ε + Λ_ψ[r]`, and `|ψ(x) − x 1_{|x|<η}| ⩽ C φ(x)` on `x ⩾ −2Δ`, which gives (ii). ∎

Remarks.
* *Primary measure.* `M_z = (1−z)^{−Δ}(μ_z + descendant monomials)`. A descendant at level `m ⩾ 1` lands in `|x| < η` only if its primary has `x ∈ (−2m−η, −2m+η)`, where by (iii) the mass is `⩽ Cρ y^{m}`; with the uniform bound on level-`m` over level-0 coefficients that follows from PQR Lemma A.1, (i)–(ii) hold for `μ_z`, i.e. the target of `TARGET.md` §2 in its stated form. This step (only) inherits the Hogervorst-formula caveat of §1.1 in `d ⩾ 3`.
* *Tilt.* The weight `y^{x/2} = z^{βx/2}` favours negative `x`; (i) shows the tilted measure has mean `−2b_*ε` and second moment `o(ε)` (not `o(ε²)`). The theorem is about the **total** weight in the twist window, averaged over spins `ℓ ~ z^{−1/2}`; it does not isolate a trajectory, and says nothing pointwise in spin.
* *Non-uniformity.* `C` depends on the correlator (through `G` at reference points) and the regime is `z^{δ/2−βΔ} ≪ 1`: with `δ ~ 1/N` (large-`N` O(N), §6) or `δ ≈ 0.02` (3d Ising `⟨σσσσ⟩`, where (H) holds with `τ_* = 1`, the stress tensor, if the spin-4 operator has twist `≈ 1.02`) the asymptotic regime is at astronomically large spin.
* *What must be checked before believing it.* (1) the positivity and leading-term normalisation of the block expansion in the stated normalisation; (2) the `O(y log(1/y))` remainder of `k_{2h}(1−y)`; (3) the elementary inequality `|ψ − x1_{|x|<η}| ⩽ Cφ`; (4) an independent test on a correlator with a genuine exchange (e.g. to first order in `1/c_T` with a tree-level exchange added to a generalised free field, where both sides are computable).

**Coverage.** Stress-tensor exchange in `d = 3` and in every `d ⩾ 3` (τ_* = d−2, `δ ⩽ d−2` is harmless), light-scalar exchange, and several operators of the same minimal twist. Excluded: infinitely many operators at the minimal twist (free and higher-spin-symmetric theories), accumulation of twists at `τ_*` (δ = 0), non-identical or charged externals (no positivity of the direct coefficients).

### 5.5 What remains genuinely out of reach

The estimate asked for in `TARGET.md` §4 (subtracted remainder in a window `β > τ_*/2`), any statement uniform up to the anti-diagonal, any pointwise-in-spin anomalous dimension, and the second sheet. For those the obstruction is real and §5.2 says why.

---

## 6. Free-field falsification test of (H) and of the conditional target [A unless marked]

**Free scalar in `d = 3`, external `φ`** (`Δ = 1/2`): `G = 1 + u^{1/2} + (u/v)^{1/2}`; `φ×φ` contains one conserved current of twist exactly 1 for every even spin, anomalous dimension zero, stress tensor included. (H) **excludes** it twice over: infinitely many operators at the minimal twist `τ_* = 1`, and `τ_* = 2Δ`, so no `δ > 0` with `τ_*+δ ⩽ 2Δ` exists (the minimal-twist operators *are* the double-twist tower). Consistently, `b_*` would be a divergent sum: the direct-channel twist-1 sum is the power `(z/y)^{1/2}`, not a logarithm.

**Every theory with infinitely many twist-1 currents in `φ×φ`** is excluded by "finitely many". [M, not re-checked] By Maldacena–Zhiboedov, a CFT3 with a unique stress tensor and one higher-spin conserved current has infinitely many and free-field current correlators; so for `τ_* = 1` (H) in effect says: the twist-1 operators in `φ×φ` are stress tensor(s) and possibly scalars of dimension exactly 1.

**Other candidates for a counterexample (a theory satisfying (H) with non-negative averaged leading anomalous dimension).**

| Theory, external operator | (H)? | Outcome |
|---|---|---|
| Free scalar, external `φ²` (Δ = 1): `G = 1 + u + u/v + c[u^{1/2} + (u/v)^{1/2} + u v^{−1/2}]` | fails: all even-spin twist-1 currents appear | zero anomalous dimensions, no conflict; the `(u/v)^{1/2}` power is matched by crossed twist-1 currents, not by a `log y` |
| Free scalar, external `φ^n` | fails, same reason | — |
| Generalised free field | holds vacuously with no operator below `2Δ`, `b_* = 0` | prediction `0`; true |
| O(N) at `N = ∞`, external `φ_i` | fails (exact higher-spin currents) | — |
| O(N) at large finite `N` | holds with `τ_* = 1` (stress tensor), `δ = τ_{ℓ=4} − 1 = O(1/N) > 0` | [M] `τ_ℓ = 1 + 2γ_φ(ℓ−2)/(ℓ−1/2)`, `γ_φ = 4/(3π²N)`, so `τ_ℓ − 2Δ_φ = −3γ_φ/(ℓ−1/2) < 0`: negative, `∝ 1/ℓ`. Sign agrees. The *coefficient* at accessible spins is `O(1/N)`, not the stress-tensor value: all nearly conserved currents compete until `log ℓ ≫ N`. This is the non-uniformity in `δ`, not a contradiction |
| 3d Ising, external `σ` | holds with small `δ` (numerical spectrum) | leading tower approaches `2Δ_σ` from below; sign agrees |
| A scalar of dimension `< 1` in the OPE | then `τ_*` is that scalar's dimension, not 1 | prediction `−2 p Γ(τ_*)/Γ(τ_*/2)² z^{τ_*/2} < 0`; the stress tensor is a competitor inside `R` |

I found **no counterexample**, and none can exist if Theorem T is right: for identical scalars `b_* > 0` whenever a minimal-twist operator is exchanged, because `p_O ⩾ 0`. Positive leading anomalous dimensions do occur, but outside the hypotheses: non-identical or charged externals (current exchange between like charges), or spinning externals violating collider-type bounds (§7).

Consistency with the standard formula [Q for the formula, `LITERATURE.md`]: `γ_0 = −2P_m Γ(Δ)²Γ(τ_m+2ℓ_m)/(Γ(Δ−τ_m/2)²Γ(τ_m/2+ℓ_m)²)`; averaging `γ_0 ℓ^{−τ_m}` against the GFF weight at `ℓ ~ z^{−1/2}` gives the factor `Γ(Δ−τ_m/2)²/Γ(Δ)² z^{τ_m/2}`, hence `−2b_* z^{τ_*/2}` [A].

---

## 7. Li–Meltzer–Poland, arXiv:1511.08025, and the collider bounds

**What they derive [Q, abstract].** "We analytically study the lightcone limit of the conformal bootstrap equations for 4-point functions containing global symmetry currents and the stress tensor in 3d CFTs. We show that the contribution of the stress tensor to the anomalous dimensions of large spin double-twist states is negative if and only if the conformal collider physics bounds are satisfied." Correlators: `⟨JJφφ⟩`, `⟨JJJJ⟩`, `⟨TTφφ⟩`, `d = 3` only ("The restriction to 3d is technical, as in this dimension all the relevant conformal blocks are known"). Example (4.18)–(4.19): `γ_{[Jφ]_0} ∝ −(3C_J − 8πλ_{JJT})` and the parity-odd family `∝ −(16πλ_{JJT} − 3C_J)`, non-positive iff `3C_J/(16π) ⩽ λ_{JJT} ⩽ 3C_J/(8π)`.

**Level of rigour: derivation.** Assumptions, as they state them [Q]:
* large spin, with the form (3.14) `γ(n,ℓ) = γ_n/ℓ^δ` imposed;
* non-degeneracy: "we will assume that for each τ_n = Δ_1 + Δ_2 + 2n there exists a unique operators at each spin ℓ, with twist approaching τ_n as ℓ → ∞. We will relax this assumption later when considering 4-point functions of conserved currents." (mixing among `[JJ]` structures is treated in their §5 and App. C);
* isolated minimal twist: "the next-to-minimal twist operators, which we will assume to be conserved spin 1 and 2 operators"; footnote 7: "we are assuming that there are no scalars with 1/2 < Δ < 1"; no higher-spin currents;
* parity: "we will only consider contributions in the s-channel from 3-point structures of even parity".

**Negativity is not derived there.** [Q] "Since gravity is expected to be attractive, we expect the anomalous dimensions to be negative." The paper proves an *equivalence of signs* (sign of the isolated stress-tensor coefficient ⇔ collider bounds), at derivation level, and for `⟨JJJ⟩` only conjectures a bound ("We do not have a rigorous way to forbid this situation, but it is tempting to conjecture that the bound (4.23) holds in all unitary CFTs"). No Nachtmann-type convexity is proved or used as input; convexity appears as the expected consequence.

**Relation to the proven bounds.** Hofman–Maldacena obtained the bounds from positivity of the energy flux. Proofs, all at physics level of rigour: Hofman–Li–Meltzer–Poland–Rejon-Barrera 1603.03771 [abstract Q: "for any unitary parity-preserving conformal field theory (CFT) with a unique stress tensor in spacetime dimensions larger than 2 … relies only on CFT first principles"; text: reflection positivity and analyticity in the lightcone limit, no higher-spin symmetry, light scalars shown not to spoil it]; Faulkner–Leigh–Parrikar–Wang 1605.08072 (ANEC from monotonicity of relative entropy) [abstract only]; Hartman–Kundu–Tajdini 1610.05308 (ANEC from causality) [abstract only]. Combining: collider bounds (derived) + LMP equivalence (derived) ⇒ the isolated stress-tensor contribution to spinning double-twist anomalous dimensions in `d = 3` is non-positive. Nothing in this chain is a theorem in the sense of PQR/vR, and nothing bounds competing exchanges.

---

## 8. Three statements, kept apart

**(i) Negativity of the *total* leading double-twist anomalous dimension at large spin.**
Status: **derivation**, for the minimal-twist family of even spin in the OPE of identical (or reflection-symmetric) scalars: Nachtmann's theorem; Komargodski–Zhiboedov 1212.4103 (via deep inelastic scattering in a gapped phase) [abstract]; Kundu 2002.12390 [abstract Q: "twist must be a monotonically increasing convex function of the spin … completely non-perturbative"]. Convexity plus the limit `2Δ` gives approach from below. Assumptions: Lorentzian analyticity and positivity, spins above a critical value, a well-defined minimal-twist family. Rigorous at theorem level: only `τ_gap ⩽ 2Δ` and accumulation at `2Δ` (PQR), the GFF density (vR); **no sign**. **False** in general for non-identical or charged externals and for non-minimal families.

**(ii) Negativity of the *isolated stress-tensor* contribution.**
Identical scalars: the coefficient is `−2P_T×(positive)` with `P_T = λ²_{φφT} ∝ Δ²/c_T > 0` by the Ward identity; **derivation** (FKPS, KZ) of the fact that this coefficient governs large-spin anomalous dimensions, the sign itself being trivial. Spinning externals: sign depends on `⟨JJT⟩`, `⟨TTT⟩` parameters; non-positive iff collider bounds (LMP, `d = 3`, **derivation**), bounds themselves derived (§7). Without the collider bounds as input the statement would be **false** as a consequence of conformal symmetry alone.

**(iii) A rigorous quantitative estimate, at finite `c_T`, of the stress-tensor contribution with a bound on competing exchanges.**
Status: **not in the literature I could find.** PQR and vR stop at leading order (§4); Caron-Huot 1703.00278 is a derivation whose error control is stated as a programme [Q: "Assuming that individual contributions are sufficiently regular and close to the average, this then provides an asymptotic 1/J expansion of this data … Gaining control over the errors in this expansion will be crucial"; "An interesting possibility is to use the convergent sums to control the errors"]; its convergent sums hold [Q] "for any individual spin J > 1". Kravchuk–Simmons-Duffin and the dispersive papers are derivations (per `LITERATURE.md`; not re-read). The Tauberian papers with remainders are Euclidean (Mukhametzhanov–Zhiboedov) or two-dimensional modular (Pal–Qiao, PQvR).

*Would a second-order DLC remainder bound be new?* Yes relative to PQR and vR: no estimate there retains the first non-identity term in a double limit. Theorem T, if correct, is such a result in the only form I can see a proof of (window `β < δ/(2Δ)`, not the window of `TARGET.md` §4).

*Physical content beyond Nachtmann and the inversion formula?* **None as physics.** The inversion formula, at derivation level, gives more: individual spins `J > 1`, the full `1/J` series, each trajectory. What Theorem T adds is status: theorem level; finite `c_T` and no large-`N` or sparseness input; no regularity-in-spin assumption; no non-degeneracy or mixing assumption (it is a statement about the total measure in a twist window, so degeneracy is irrelevant, and for the same reason it cannot resolve trajectories); an explicit power-saving bound `z^{δ/2−βΔ}` on all competing exchanges. Relative to Nachtmann: it supplies the coefficient of an average, not a sign for each spin.

*Assumptions any such theorem needs, explicitly.* (1) Large spin only as an asymptotic average, `ℓ ~ z^{−1/2}`, `z → 0`; nothing at fixed spin. (2) A spectral condition in twist, (H): the exchanged minimal-twist operators finite in number and isolated by a gap `δ > 0` from all other twists, in particular from `2Δ`. (3) `c_T` arbitrary and finite; constants are not uniform in `c_T`, in `δ`, or in the correlator. (4) Identical Hermitian external scalars (positivity in both channels). (5) Discreteness of the spectrum for the passage from monomials to primaries. (6) No assumption on mixing, degeneracy, parity, or the Regge limit.

# m1_interchange: is the block-wise light side the light side? (referee point M1, remaining question)

27 Nov 2026. Unrefereed. Nothing outside `ads4/code/m1_interchange/` was edited and nothing was committed.

Read-only imports: `m1_direct/md_ladder.py`, `m1_direct/md_k3.py`, `m1_remainder/m1_weight.py`. Data read: `prl/data/kernels.json`
(`D_k2_EE`: R₋, R₊), `prl/data/nf02_psi.json`, `prl/data/k2/k2_continuum_D_1000_rep.json` (through `m1_weight`),
`m1_remainder/out/m1a_hhat_*_b*.json`, `m1_remainder/out/m1c_physical_sum.txt`. Literature: CHMRS (arXiv:2106.10274), read
from the PDF (§4.1 and App. E).

Grades: [P] proved, [D] derived, [D-sketch], [DV] derived and checked numerically, [NE] numerical evidence, [E] estimate,
[A] assumption, [I] open.

## 0. Verdict

**Both alternatives of the brief hold in part. The m1 tension is not removed.**

1. **The block-wise evaluation fails on the light tower.** [P given that (fp-H) is the continuation, row L1 and Lemma R]
   * Test object: the exactly known tree-level Einstein `⟨TTTT⟩`.
   * Its O(1/c_T) light side per ν is `−Ĉ′_{EE,ν}|_E = 16πG/(ν² + 9/4)` (row L1 / Lemma R).
   * Its O(1/c_T) light content in the s-channel is the T block plus the `[TT]` double traces at double zeros.
   * The continued per-block action of the T block, i.e. (fp-H) at (Δ, J) = (3 + ε, 2), has a **simple pole at Δ = 3**. The same
     holds at every twist-1 current, Δ = J + 1.
   * The double traces are regular at every even twist (checked to J ≤ 40). So they contribute 0 × finite = 0.
   * Every continued per-ν image is (fixed rational factor ρ(ν)) × (polynomial in ν²). The honest value is not of this form in any
     of the three normalisations used in the project: it has a pole at ν = ±3i/2 or ±i/2.
   * Hence no regularisation of `λ²_TTT ω^cont[T] + Σ 0` reproduces the light side. This holds per ν and for every weight.
   * Consequence: the term-wise numbers of the twist-2, 3, 4 families are not shares of the remainder. That includes
     m1_remainder's |c_T²R| ≈ 10^{679}–10^{687}, which is twist-2 dominated.
   * The paper's own L₁ is computed in a different way, from Sec. `sec:arc` (light-ray matrix elements, Lemma R or L4). So the
     project has two inequivalent "actions on a light block". L0's continuation is not the one in the sum rule.
2. **The giant values are real where the per-block integral converges.** [DV for the values; D-sketch for the family sums]
   * For twist > 5 CLS-D's per-block integral converges absolutely (m1_direct §4.5), so the per-state actions are honest.
   * They are huge, and they come from ν ≤ Λ: the weight beyond Λ contributes ≤ 4·10⁻⁹. Examples for twist 6 (Δ_g = 1000), by
     independent ν-quadrature, agreeing with the moment balls to all printed digits:

     | J | ĥ |
     |---|---|
     | 2 | 1.8·10¹⁵ |
     | 8 | 4.0·10³⁸ |
     | 20 | 1.8·10⁷⁹ |

   * Family sums at fixed twist are absolutely convergent (Fubini), so they too are honest.
   * With tree-level Einstein data the twist-6 `[TT]_{2,J}` family alone gives |c_T²R_{>5}| ≈ 10^{668}–10^{673} [E].
3. **Resolution.** [D, given the exact sum rule and Σ_heavy ≥ 0]
   * Write `R = R_{>5} + R_{≤5}`. Here R_{≤5} is the light-cone sector: the twist-2, 3, 4, 5 states, other than T. It is defined
     only on the light correlator as a whole, as `Ω[𝒢_{≤5}] − c_T⁻¹L₁`.
   * The exact sum rule fixes the total: `R = −Σ_heavy − c_T⁻¹L₁`.
   * In any theory with Einstein-like light data at c_T ≲ 10^{668}, R_{≤5} must therefore cancel R_{>5} to ~670 digits. No input
     of the sum rule has to fail for that.
   * What fails is the idea that R can be bounded from light OPE data. Light OPE data enter amplified by up to e^{πΛ/2}10^{4.5}.
     A2′ controls only Σ|a|²(1+J)^{5/2}log(2+J), and the light-cone sector has no per-block representation at all.
4. **Honest size of the light remainder.** At O(1/c_T) it is L₁ exactly (Lemma R), i.e. `p_Tω_T = −16L_A n₀G = −0.066/c_T` in the
   paper's units. At O(1/c_T²) it is **not obtainable** [I].
   * Term by term it is meaningless for twist ≤ 5, and ≳ 10^{668}/c_T² for twist > 5.
   * The only sound route is the whole light correlator: CHMRS's Regge-contour formula with loops, which "doesn't necessarily
     localize to a residue". No one has evaluated it for this functional.
   * A flat-space-limit model (one-loop arc ~ GΛ² × tree) would give c_T* ~ 10⁵–10⁸ [M]. It is not usable, because (2) shows that
     the whole must cancel a 10^{670} block-wise piece that the model knows nothing about.
5. **c_T threshold for the cubic bound.** None can be derived. Eq. (5) is a statement about c_T → ∞ at fixed Δ_g.
   * The planar step itself needs repair. "By L0, c_T|R| ≤ … so A2′ gives c_TR → 0" uses the block-wise light side.
   * For twist > 5 it stands [D-sketch].
   * For the light-cone sector it needs an assumption on the light correlator as a whole [A]: "c_T R_{≤5} → 0". This is CHMRS's and
     CLS-D's standing assumption (light data = tree-level EFT up to o(1/c_T)), or a bound on the Sec. `sec:arc` block action of
     twist ≤ 5 states, which is not available.
   * The angular bound (τ_* ≈ 0.0178, c_T ≳ 2·10⁸) is unaffected.

## 1. Ledger

| step | what | file / output | grade |
|---|---|---|---|
| 0 | read m1_direct, m1_remainder, supplement (Order of limits, sec:HK, sec:arc, sec:c1, L0, L1, H_K, Scope), letter l. 136–140 | — | — |
| 1 | CHMRS §4.1 and App. E from the PDF: light side defined on the EFT correlator as a whole; remark after (4.7) | §5 | [Q] |
| 2 | Step 1: R3's per-block integral on light blocks | §2 | see §2 |
| 3 | GFF at O(1): double commutator = s-channel identity only | §3.1 | [P] (Wick) |
| 4 | tree-level Einstein test: ladder at the T block, twist-1 poles, even-twist regularity, per-ν profile | `m1i_tblock.py` → `out/tblock.txt` | [P] (exact sympy) |
| 5 | twist ≥ 6 per-state actions by an independent route; share of ν > Λ | `m1i_tw6.py` → `out/tw6.txt` | [DV] |
| 6 | twist-6 family size with tree-level Einstein data | from `m1_remainder/out/m1c_physical_sum.txt` | [E] |
| 7 | verdict, honest R, threshold, replacement text | §0, §4, §6 | — |

## 2. Step 1: does R3's cross-ratio integral converge on a single low-twist block?

* **On the Regge box: yes, for every block and every sub-family** [D, given R1′, R3 as stated in sec:HK]. Per-block dDisc
  positivity and Tonelli give `∫_{B_c}|K| dDisc F ≤ ∫_{B_c}|K| dDisc 𝒢 ≤ 2B(c)∫_{B_c}|K|` for any sum F of blocks, light or heavy.
  So whatever goes wrong for light blocks happens outside B_c (c < 1), at the light-cone edges.
* **At the light-cone edge: not absolutely, and not analytically in Δ down to twist 1.**
  * CLS-D's per-block D′ integral diverges for τ ≤ 5, as B^{(τ−5)/2−1} (m1_direct §4.5) [D + NE].
  * Suppose R3's cross-ratio integral of the J = 2 block converged absolutely and locally uniformly on τ ∈ [1, 5 + δ). Then it
    would be analytic there and equal to the continuation (fp-H). But (fp-H) has a pole at τ = 1 (§3.2). So it does not.
    [P given (fp-H) = continuation on τ > 5]
  * For a twist-2 block specifically, convergence in the cross-ratio representation was **not determined** [I]. It is not needed
    for the verdict, because the tree-level test (§3.2) decides the bookkeeping question directly.
* **Consequence for H_K.** By R3's own inequality, a divergent per-block light-cone integral of a block with nonzero dDisc
  coefficient (T itself at O(1/c_T), the twist-2, 4 double traces at O(1/c_T²)) forces ∫|K| dDisc 𝒢 = ∞ over the full region.
  * H_K's absolute convergence is therefore a statement about the Regge box plus the heavy blocks' edges (K2). It says nothing
    about the light-cone sector of the light states.
  * This is consistent with sec:HK as written. The text should not be read as covering the light tower.

## 3. Step 2: the interchange on exactly known correlators

### 3.1 GFF at O(1) [P]: not discriminating
By Wick's theorem the crossed terms cancel in `⟨[φ₄,φ₃][φ₁,φ₂]⟩`: ⟨4312⟩ and ⟨4321⟩ have the same (41)(32) and (42)(31)
pairings. So the GFF double commutator is the s-channel identity alone. Each double trace has sin² = 0 and a finite continued
value, so both sides are 0. Nothing is tested.

### 3.2 Tree-level Einstein `⟨TTTT⟩` at O(1/c_T) [P given (fp-H) as continuation, Lemma R, L1] (`out/tblock.txt`)

**Whole.** The light side is Lemma R's residue of the full tree correlator: `−Ĉ′_{EE,ν}|_E = 16πG/(ν²+9/4)` and `Ĉ₃|_E = 0`.

**Block-wise.** The O(1/c_T) s-channel light content is the T block plus `[TT]_{n,J}` at double zeros (L1).
* **Exact results:**
  * **The T block is singular.** (fp-U) has the factor (b + 2 − Δ)_k in the denominator. At (Δ, J) = (3 + ε, 2) the k = 2 term of
    the b = 0 ladder has T − u₁ = ε(ε+3). The per-ν images have simple poles at ε = 0, with polynomial residues:
    * k = 2 residue: `−(s⁵ − 192s⁴ + 8228s³ − 94752s² + 291552s − 178560)/144`, with s = ν² + ¼;
    * k = 3 residue: `s² − 18s + 40`.
  * **Every twist-1 current is singular.** At Δ = J + 1 the denominator is T − u_{J−1} = ε(ε + 2J − 1) (J = 2…12 listed): a simple
    pole.
  * **The double traces are regular.** No pole at any even twist τ = 2…12, J ≤ 40, b = 0, 4; the smallest |T − u_j| on the ladder
    is 2.
  * **The per-ν profiles differ.** A continued block image is ρ(ν)·S(s), with S a polynomial (finite ladder) and ρ one of the
    conventions in use. In each, (16π/(ν²+9/4))/ρ has a pole at finite ν:

    | ρ | (16π/(ν²+9/4))/ρ |
    |---|---|
    | c = (4ν²+1)/((4ν²+9)(4ν²+49)) (identity I2) | 64π(4ν²+49)/(4ν²+1) |
    | c/(R₋+R₊), with R₋ + R₊ = −8(4ν²+1)(16ν⁴−728ν²+969)/((4ν²+9)(4ν²+49)) | −512π(16ν⁴−728ν²+969)/(4ν²+9) |
    | 1 (g-units, m1_n0) | 64π/(4ν²+9) |

* **Conclusion.**
  * Take the residue, the finite part or any other regularisation of λ²_TTT ω^cont[T]. Adding Σ 0 from the double traces never
    gives the light side, at any ν.
  * The light side at O(1/c_T) lies outside the span of the continued images of any finite set of blocks.
  * So the identity "light side = Σ_{light} p_O × (L0 action)" is false on an exact correlator, at the first order where the light
    side is nonzero.
* **What this does and does not say.** The paper evaluates its T term in Sec. sec:arc and sec:c1 (light-ray matrix elements
  ⟨Ψ₄|A₃P_TA₁|Ψ₂⟩; L4) and gets a finite value. That is a different per-block object from L0's continuation. It must agree
  with the ladder on heavy states (the heavy side uses the ladder), but it does not at τ = 1.
  * So L0 (continued) is not the action that enters the sum rule on twist ≤ 5 states.
  * Whether the Sec.-arc action of a twist-2 double trace is finite, and how large, is not known [I].
* **Side remark** [I, not load-bearing]. m1_n0.py writes the T term in g-units as `16π/(ν²+9/4)`, i.e. with ρ = 1.
  * Converting the stated light value `−Ĉ′_{EE,ν}|_E = 16πG/(ν²+9/4)` to g-units costs a factor 1/c (χ-weight convention,
    Ψ̂_ν images = R₋P_{l+1} + R₊P_{l−1}) or (R₋+R₊)/c (X-weight convention), not 1.
  * These differ by (4ν²+9)(4ν²+49)/(4ν²+1), resp. by the quartic −8(16ν⁴−728ν²+969).
  * n₀, and with it every log₁₀ c_T* quoted from m1_remainder, should be re-derived in one stated convention.

### 3.3 Small gap
The brief proposed starting at small gap. That is not needed: the tree-level test is exact and gap-independent (per ν). The
twist-1 pole and the profile mismatch hold for every weight.

## 4. Twist > 5: the giant values are honest (`out/tw6.txt`)

**Independent per-state values.** ν-quadrature of the published Δ_g = 1000 continuum weight against own ladder polynomials
(k = 2), plus exact Q_k (k = 3). They agree with m1a's moment balls to all printed digits.

| state (twist, J, b) | ĥ⁽²⁾ | ĥ⁽³⁾ | ĥ | log₁₀\|ĥ\| [m1a] | share of ĥ⁽²⁾ from ν > Λ |
|---|---|---|---|---|---|
| 6, 2, 0 | −5.02·10¹⁴ | +2.27·10¹⁵ | +1.769·10¹⁵ | 15.248 [15.248] | −3·10⁻¹³ |
| 6, 4, 0 | −9.20·10²² | +3.06·10²³ | +2.136·10²³ | 23.330 [23.330] | −2·10⁻¹² |
| 6, 8, 0 | −2.71·10³⁸ | +6.67·10³⁸ | +3.957·10³⁸ | 38.597 [38.597] | −2·10⁻¹¹ |
| 6, 8, 4 | −2.44·10²³ | −8.12·10²³ | −1.056·10²⁴ | 24.024 [24.024] | −2·10⁻¹² |
| 6, 20, 0 | −1.47·10⁷⁹ | +3.29·10⁷⁹ | +1.819·10⁷⁹ | 79.260 [79.260] | −4·10⁻⁹ |
| 6, 20, 4 | −4.05·10⁶⁶ | −8.52·10⁶⁶ | −1.258·10⁶⁷ | 67.100 [67.100] | −8·10⁻¹⁰ |
| 8, 8, 0 | | | +2.57·10³⁷ | 37.410 | −2·10⁻¹¹ |
| 10, 8, 0 | | | +2.61·10³⁶ | 36.416 | −2·10⁻¹¹ |
| 22, 8, 0 | | | +2.61·10³² | 32.417 | −2·10⁻¹¹ |

* **Not a tail artefact.** The k = 3 piece has compact support ν ≤ Λ, and the k = 2 tail beyond Λ is negligible.
  * The size comes from the per-ν images of light states growing above their mass. This is the scalar-proxy statement of the
    Supplement's Cutoff section, "ln Γ ≈ m G(ν/m), rising to π/2".
  * The weight's reach ν ≈ Λ does the rest.
  * Changing the cutoff or the tail of the weight therefore cannot remove the giant values.
* **Twist 22, J = 8 is in m1_direct's verified range** (direct D′ integral = ladder to 6–9·10⁻⁴, at ν ≤ 6). Its action is
  2.6·10³². The direct checks stop at ν ≤ 6–8, below these states' masses. The large-ν growth rests on the ladder's derivation
  plus analyticity [D], not on a direct large-ν check.
* **Family sums** [D-sketch].
  * On D′ the modulus of the integrand is ν-independent and the smeared weight is absolutely integrable (cutoff θ). Absolute
    convergence of a fixed-twist family therefore reduces to Σ_J p_J ∫_{D′}|Φ_J| < ∞, which is not checked.
  * Given it, Fubini makes Σ_J p_J ω[O_J] the honest family action, with no continuation involved.
* **Size with tree-level Einstein data** [E] (m1c, A_n = κ_n = 1, J ≤ 1000, even J). For the `pp` twist-6 family:
  * signed Σ c_T|a|²ĥ = 9.4·10⁶⁵⁹, abs 4.8·10⁶⁶²;
  * × 64Λ⁴/9 = 7.1·10¹²: |c_T²R_{>5}| ≈ 7·10⁶⁷² (truncated), 3·10⁶⁷⁵ term-wise;
  * a smooth full spin sum gained 10^{4–5} at twist 2 (m1_remainder §3), so ~10^{668} is the plausible order;
  * threshold read naively: c_T* = π|R₂|/(512 n₀L_A) ≈ 10^{669}–10^{674}.

## 5. Literature [Q] (CHMRS, arXiv:2106.10274)

* **§4.1, after (4.7).** "Note that Π_{k,η}|light is not equal to a sum over light operators of OPE coefficients times
  Π_{k,η}[G^s_{Δ,J}] (since the expansion in Regge moments does not necessarily converge on light operators)."
  * Ω^EE is a ν-superposition of such moments (C_{k,ν}-type per-ν functionals).
  * L0's per-state light values are exactly that sum read state by state.
* **§4.1, after (4.4).** "when loops are included, ω|light can still be computed from its Regge moments via a contour integral
  close to the Regge limit; however the integral doesn't necessarily localize to a residue."
  * This is the only route to an honest O(1/c_T²) light side. It is not done for Ω^EE.
* **App. E.** The light side is defined on the EFT correlator as a whole, eq. (E.3). It assumes that dDisc_s G^EFT decays along the
  negative axes so that (E.3) is well defined. For ⟨TTTT⟩ this is Lemma R's R5, "assumed, not checked" (H7_I14 §3).

## 6. Proposed replacement text (not applied)

**Supplement, "Order of limits" (sec:HK), replace "By L0, $c_T|R|\le(2\pi c_TG\Lambda^4/3)\max(C_h,C')\,\mathcal E_{\rm light}(c_T)$
up to $O(1)$ factors, so A2$'$ gives $c_TR\to0$; at each fixed spin the double zeros make the contribution $O(c_T^{-2})$ under
factorisation alone, and only the uniformity in spin is assumed." with:**

> Split $R=R_{>5}+R_{\le5}$ by twist. For twist $>5$ CLS-D's per-block integral converges absolutely, so the family sums are
> absolutely convergent and L0 with A2$'$ gives $c_TR_{>5}\to0$ [D-sketch]; at each fixed spin the double zeros make the
> contribution $O(c_T^{-2})$ under factorisation alone, and only the uniformity in spin is assumed. The light-cone sector $R_{\le5}$
> (twist $2$--$5$, $T$ excluded) has no per-block integral, and the continuation of L0 is not its value (below); $c_TR_{\le5}\to0$
> is an assumption on the light correlator as a whole, the standing assumption of Refs.~\cite{CHMRS,CLSD} that the light side is that
> of a tree-level bulk theory up to $o(1/c_T)$ [A].

**Same paragraph, replace from "The size of the light-side remainder is not a power of $\Delta_g$." to "(its light remainder is
$\le70\,\mathcal E_{\rm light}$ with no exponential factor)." with:**

> The light-side remainder is not controlled at finite $c_T$. For twist $\le5$ the action of $\Omega^{EE}$ on a light block exists only
> as the continuation of Eq.~\eqref{eq:fp-target} in $\Delta$, and that continuation is not the light-side action. On the exact
> tree-level Einstein correlator, whose $O(1/c_T)$ light content is the $T$ block and the $[TT]$ at double zeros (L1), the continued
> action of the $T$ block has a pole at $\Delta=3$ (at every twist-1 current, $\Delta=J+1$), and every continued per-$\nu$ image is a fixed
> rational factor times a polynomial in $\nu^2$. The light side, $16\pi G/(\nu^2+\tfrac94)$, is not of this form in any normalisation
> (\texttt{code/m1\_interchange}; cf.\ the remark after Eq.~(4.7) of Ref.~\cite{CHMRS}). The term-wise numbers of the twist-2, 3 and 4
> families, $|c_T^2R|\approx10^{679}$--$10^{687}$ with tree-level Einstein data, are therefore not shares of the remainder.
>
> For twist $>5$ the per-block actions are honest and large. Above its mass a light state's per-$\nu$ image grows like $e^{\pi\nu/2}$, so
> $[TT]_{2,J}$ has $\hat h=1.8\times10^{15}$, $4.0\times10^{38}$ and $1.8\times10^{79}$ at $J=2$, $8$ and $20$, all from $\nu\le\Lambda$.
> With tree-level Einstein data the twist-6 family alone gives $|c_T^2R_{>5}|\sim10^{668}$--$10^{673}$ [E].
>
> The sum rule fixes the total, $R=-\Sigma_{\rm heavy}-c_T^{-1}L_1$. In a theory with Einstein-like light data at moderate $c_T$ the
> light-cone sector must therefore cancel $R_{>5}$. No bound on $|R|$ at finite $c_T$ follows from light OPE data, which enter amplified
> by up to $e^{\pi\Lambda/2}10^{4.5}$, and none is claimed. Eq.~(5) is a statement about $c_T\to\infty$ at fixed $\Delta_g$, with or without
> a bulk effective theory. The angular bound, Eq.~(4), carries no such constant: its light remainder is $\le70\,\mathcal E_{\rm light}$, with
> no exponential factor.

**Supplement, Scope, replace "and any control of the light remainder at finite $c_T$: block by block it is
$\sim e^{\pi\Lambda/2}c_T^{-2}$ ($c_T^*\approx10^{681}$ at $\Delta_g=1000$, row L0), the per-block actions of the twist-2 and twist-4
double traces exist only as analytic continuations, and a remainder of this size is inconsistent with the sum rule (open)" with:**

> and any control of the light remainder at finite $c_T$, together with the planar step for the light-cone sector (twist $\le5$), which
> is assumed. For twist $\le5$ the continued actions of row L0 are not the light-side action: they fail on the exact tree-level
> correlator. The absolutely convergent twist $>5$ families are $\sim10^{670}/c_T^2$ term by term with tree-level Einstein data. Only
> the total, fixed by the sum rule, is meaningful (open).

**Table inputs, row L0, grade column**: after "for twist $\le5$ the unique continuation in $\Delta$" insert "; this continuation is
not the light-side action (pole at every twist-1 current; fails on the tree-level Einstein correlator, [P given the ladder])", and
replace "the block-wise sum as light side of the finite-$c_T$ sum rule [I]" with "the block-wise sum as light side: twist $>5$
[D-sketch], twist $\le5$ false". **Row $H_K$**: replace "light remainder [D-sketch] through L0" with "light remainder: twist $>5$
[D-sketch] through L0, twist $\le5$ [A]". **Table c1, row (ii)**: replace "$c_TR\to0$ [D-sketch given A2$'$] (L0)" with "$c_TR\to0$:
twist $>5$ [D-sketch given A2$'$] (L0), twist $\le5$ [A]". The "Net grades" sentences (l. 389, 427) then read "with the light-remainder
step [D-sketch] for twist $>5$ through L0 and [A] for the light-cone sector".

**Letter (letter_v2.tex l. 138), replace "The functional's action on these states is exact but grows like $e^{\pi\Lambda/2}$, so
their remainder is not controlled at any realistic $c_T$ (term by term it is $\sim10^{681}/c_T^2$), and Eq.~\eqref{eq:cubic} is a
statement about the planar limit; the angular bound below has no such constant." with:**

> The functional's action on light states grows like $e^{\pi\Lambda/2}$, and on the lowest-twist ones it is not defined state by
> state, so the light remainder is not controlled at any realistic $c_T$ and Eq.~\eqref{eq:cubic} is a statement about the planar
> limit; the angular bound below has no such constant.

## 7. Files
* `m1i_tblock.py`: the tree-level test. Covers the T-block Laurent expansion of (fp-H), the twist-1 poles, regularity at even
  twist, and the per-ν profile in the three normalisations. Exact sympy, about 2 min. Output: `out/tblock.txt`.
* `m1i_tw6.py`: twist ≥ 6 per-state actions by ν-quadrature of the published weight, plus the share from ν > Λ. About 8 min;
  the weight evaluation dominates. Output: `out/tw6.txt`.

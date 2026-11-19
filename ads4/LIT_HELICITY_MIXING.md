# Literature check: helicity mixing, large impact parameter, and t₄ bounds in AdS4/CFT3

*19 Nov 2026, literature agent. Not committed. No existing file was edited. Sources were downloaded to the session scratchpad
(`…/scratchpad/lit/helmix/*.pdf|txt`) and read as text. Equation and section numbers are those of the arXiv PDFs (v1 or the latest
version). Marks: **[R]** I read the passage. **[A]** abstract only. **[I]** my inference, not stated in the paper. **[U]** not verified.*

---

## 0. Bottom line

1. **No published work derives graviton-coupling bounds from spinning CFT sum rules, in any d.** CLS-D (2311.04271) stop at the
   bulk-point dictionary. They explicitly defer the Regge regime and the remaining sum rules [R, §2.2.2 below (2.24); §5]. None of the
   15 INSPIRE citations of CLS-D (checked 19 Nov 2026, titles + abstracts of the relevant ones) does this either. **I found no paper
   that treats a curvature-induced helicity off-diagonal of a CFT functional.** The obstruction in `LARGEJ_OFFDIAG.md` is new to the literature.
2. **The obstruction has a flat-space twin with the same scaling.** In D = 4 the IR regulator makes impact-parameter functionals
   negative on the band `m_IR^{-2/3} ≲ b ≲ m_IR^{-1}` (2201.06602 eq. (3.6), Fig. 5; 2106.10274 eq. (5.59)). With `m_IR → 1/R_AdS` and
   `b̂ = bM`, the band is `b̂ ≳ Λ^{2/3}`, exactly the project's `b̂_c ≈ 0.65 Λ^{2/3}`. For scalars CHMRS showed that the AdS curvature
   *replaces* this negative IR tail with a positive one (Q₁(β) > 0; 2106.10274 §5.5). The project's finding is that for gravitons the same
   curvature order also produces a helicity off-diagonal that undoes this cure. The flat-space literature never cured its own band: it
   accepted the negativity and called the bounds rigorous only "if one assumes additional information about the spectrum at large impact
   parameters" (2106.10274, below (5.60); 2201.06602 §3.1, after (3.2)).
3. **The weight-independence of the obstruction (`|X|/F → 24|G₁|` for every `ψ₂` with `ψ₂′(0) ≠ 0`) matches a literature theorem.**
   In D = 4, positivity of the transverse Fourier transform forces `lim_{p→0} ψ(p)/p > 0` (2106.10274 eq. (5.54); 2201.06602 §3.1:
   "finiteness on the graviton pole requires ψ(p) to vanish faster than p … impossible for the Fourier transform of a positive
   function"). So `ψ₂′(0) > 0` is forced. Within a single k = 2 rule like Ĉ₂ (= B⁽¹⁾₂), nothing can be traded against it [I].
4. **In flat D = 4, the k = 2 rule that Ĉ₂ maps to is exactly helicity-diagonal.** The only k = 2 sum rule with a `(++)(+−)` cross term
   is a *different* rule, B⁽³⁾₂ (fixed-t, the helicity-violating amplitude g). CLS-D's construction does not produce it (§2 below). This
   makes "a third sum rule with its own off-diagonal" concrete: the missing object is the CFT uplift of B⁽³⁾₂, or of the CEMZ helicity
   *matrix* of sum rules (2201.06602 §3.5).
5. **Known t₄-vs-gap scaling in d = 3 is parametric only** (§3 below): `t₄ ≲ Δ_g⁻⁴` (CEMZ, CHP, Meltzer–Perlmutter) or
   `≲ log Δ_g / Δ_g⁴` (Afkhami-Jeddi–Kundu–Tajdini, explicitly d = 3). The sharp-coefficient expectation is
   `|t₄| ≲ c √(log Δ_g)/Δ_g⁴`, from the flat D = 4 bound `|ĝ₃|²M⁸ ≤ 24.9 log(M/m_IR) − 27.6` (2201.06602 eq. (4.4)), which the
   authors say should carry over to AdS4 with `log(M R_AdS)` (§5). This is asserted there, not derived. No paper gives a CFT-derived
   sharp t₄ bound in d = 3.

---

## 1. CLS-D, arXiv:2311.04271 [R]

* **Bounds?** None. The introduction states the goal (eq. (1.3), `|α₂| ≤ #/Δ²_gap`). The discussion (§5) lists what is still needed
  [R, verbatim gist]: (i) the "additional CFT sum rules that give rise to the other flat space sum rules which cannot be written in
  terms of the shock amplitudes"; (ii) "check the positivity of the sum rules in the Regge limit, where Δ and J are both large";
  (iii) a light-state dictionary. Only after these "one can then derive bounds … in terms of the large CFT gap".
* **Positivity on heavy states.** They review it as "the spectral density Im a_J(s) should be positive (or more generally a positive
  semidefinite matrix)" (§2.2.1, text below (2.17)). Their heavy action is a **matrix in exchanged 3-pt structures**,
  `Ψ[Δ,ρ]_{(a′b′)}` (eq. (4.34)), related to flat partial waves `π_{ρ,(a′b′)}` by (4.33). In the forward limit the partial wave is
  "positive (or positive semi-definite for multiple tensor structures) by construction", and the CFT-structure basis preserves
  this (App. C.1, around (C.3)). So the intended criterion is PSD over structures, as in the project's 2×2. They never test it on an
  actual functional.
* **Large impact parameter / Regge / AdS curvature.** Not treated. Verbatim: "We leave an exploration of the Regge regime for
  four-point functions of spinning operators to future work" (§2.2.2, after (2.24)). Their saddle analysis is for `Δ, ν ≫ 1` with
  `ν < m`. All results carry `1 + O(J²/m²)` corrections, i.e. bulk point only.
* **Helicity off-diagonals.** Not encountered: no Regge or curvature analysis is done. The dictionary is computed in generic D
  (29-generator basis, D ≥ 8; eqs. (4.81)–(4.88)) with **no D = 4 reduction** [R: no D = 4 / d = 3 specialisation in the source].
* **Which sum rules.** The CFT rules map only to the **lowest-subtracted, fixed-u, shock-polarised** flat sum rules: 7 even-k
  (4.82) and 3 odd-k (4.83). Dictionary matrices `M^±_{4g}` are in (4.86)–(4.88). Allowed k: `k ≥ 2` for even j′ and `k ≥ 3` for odd
  j′. For `J_e = 2` the k = 2, 3 rules are unsubtracted and vanish on all double traces (§3.5, text after (3.75)–(3.76)). Rules
  with other external polarizations, i.e. derivatives of light-transformed operators, are left open (§5).
  * [I] In D = 4 the shock-type fixed-u rules of 2201.06602 are the MHV rules B⁽¹⁾_k only (2201.06602 (2.26)). This agrees with the
    project's finding that Ĉ₂, Ĉ₃ map to stripped B⁽¹⁾₂, B⁽¹⁾₃. The fixed-t rules B⁽²,³,⁴⁾ (2201.06602 (2.27)) are of the
    "non-shock" kind CLS-D leave open.

## 2. CHLPS 2201.06602 (D = 4) and 2205.01495 (D ≥ 5) [R]

**Helicity structure of the D = 4 sum rules.**
* Unitarity: `Im a_J` is a PSD matrix, written `2(c^{−h₄,−h₃})* c^{h₁h₂}` (eq. (2.13)). Diagonal entries: `Im a^{+−+−}_J = |c^{+−}|²`,
  `Im a^{++−−}_J = |c^{++}|²` (2.14). Off-diagonal: `Im a^{+++−}_J = c^{++}c^{+−}` (2.17).
* **B⁽¹⁾_k (MHV f, fixed-u) is helicity-diagonal.** Heavy side: `⟨|c^{++}|² d̃^J_{00}(x) + |c^{+−}|² d̃^J_{44}(x)⟩` (2.33; odd k: (2.35a)).
  The cross term `c^{++}c^{+−} d̃^J_{40}` appears **only in B⁽³⁾_k** (fixed-t, amplitude g = M(+++−)) (2.27), (2.35b).
* Positivity check in practice: plot `F^{++}` and `F^{+−}` separately against `(m, b = 2J/m)` for J ≤ 300 (§3.2, Figs. 3–4). In the
  scaling limit `F^{++} = F^{+−} = ∫ψ J₀(pb)` (3.5). For the explicit functionals (3.3) built from B⁽¹⁾₂, B⁽¹⁾₃ and forward B⁽¹⁾₄, this
  diagonal check is complete because those rules have no off-diagonal [I from (2.33)].
* Full numerics (§3.4, §4.2, App. C): all improved k = 2 rules B⁽¹⁻⁴⁾₂ including B⁽³⁾₂, polynomial wavefunctions (3.17), and SDPB
  with **n×n PSD heavy densities** per little-group irrep (C.1). Two points matter here:
  * "We require that all sum rules possess the same large-b behavior in the scaling limit … (One can prove that this is a necessary
    condition for obtaining a positive functional.)" The minimal powers are `(n⁽¹⁾,n⁽²⁾,n⁽³⁾,n⁽⁴⁾)_min = (1, 10, 6, 1)` (3.16). So the
    off-diagonal rule B⁽³⁾₂ enters with `ψ ∼ p⁶` and is subleading at large b.
  * Large-J at fixed m is controlled explicitly by a polynomial-in-J matrix positivity ((C.3); "unlike … scalar particles … it is
    important").
* **Large-b, matrix-valued trick (App. C, (C.6)–(C.14)).** At large b the scaling-limit action is `A(b) + B(b)cos φ + C(b)sin φ`
  with `φ = b − π(d−1)/4`. Treat φ as free and require the `2n×2n` block matrix `M_λ(b)` to be PSD. By a minimal-degree-variety theorem
  (Segre embedding), this is *equivalent* to positivity for all φ and all vectors v (C.12)–(C.14). Then expand in 1/b and impose PSD
  of a matrix polynomial for `b ≥ B`.
* **CEMZ as a PSD matrix of sum rules (§3.5, (3.20)–(3.28)).** Scatter a graviton of arbitrary helicity superposition against a
  fixed-helicity target. The forward 2×2 matrix of sum rules is
  `[[½B⁽¹⁾₂, p̂⁴B⁽³⁾₂],[p̂*⁴B⁽³⁾₂*, ½B⁽¹⁾₂]]` (3.22). Its heavy side in b-space is a Gram matrix, hence PSD per state (3.24). Its light
  side is `4G[[log(1/(bm_IR)), 24ĝ₃/b⁴],[c.c., log]]` (3.25). This gives the parametric CEMZ bound `|ĝ₃| ≲ M⁻⁴ log(1/(bm_IR))` (3.28).
  This is the literature's model of a **PSD functional that mixes a helicity-diagonal and a helicity-off-diagonal sum rule**.

**Why D = 4 needs an IR cutoff (2201.06602 §3.1; 2106.10274 §5.5).** Positivity at all b needs `ψ(p) ∼ p` at small p
(2106.10274 (5.54)), so the action on the 1/p² graviton pole diverges logarithmically. The cutoff `p > m_IR` makes the functional
negative at large b: `lim F ≈ 160/(m⁴b³) − 260 m_IR² J(m_IR b)/(m⁴ m_IR b)` (3.6). There is a negative plateau for
`m_IR^{-2/3} ≲ b ≲ m_IR^{-1}` (Fig. 5). The authors "accept that this causes negativity at large impact parameters" (§3.1). Results:
`|ĝ₃|²M⁸ ≤ 24.9 log(M/m_IR) − 27.6` and `g₄M⁶/(8πG) ≤ 12.3 log − 13.5` ((4.4)–(4.5)). Conclusions (§5): "(4.4) translates to a bound on
three-point OPE coefficients of stress tensors in holographic 3d CFTs. In AdS4, infrared logarithms are simply replaced by
log(M R_AdS)." This is stated as an expectation, not derived.

**D ≥ 5 (2205.01495).** There are no IR logs ("D ≥ 5 where the S-matrix does not suffer from infrared ambiguities", abstract).
Unitarity is `0 ≤ Im a_ρ ≤ 2` as a matrix inequality (§II.B, below (16)). The partial waves are full SO(D−1) with multiple
vertices per irrep (Young-tableau construction, §III). Functionals (App. E): improved `B_k^imp` up to k = 5 or 7, 5–6 powers of p per
rule, a factor `p^α(M−p)` so k = 2 rules decay at large b (`∼1/b³` in D = 5), 200–400 trial rules. SDPB enforces positivity on
10⁴–2·10⁵ sampled states (J ≤ 400) plus the m → ∞ scaling limit. AdS5 output: `(a−c)/c ≤ 23/Δ²_gap + O(Δ_gap⁻⁴)` (30). **The paper
says nothing about D = 4 except that it is the IR-problematic case** [R]. Its uplift statement ("flat-space dispersive bounds uplift in
AdS to rigorous bounds on holographic CFTs, up to corrections suppressed by 1/(M R_AdS)", §IV) is asserted by citing CHMRS. It is not
checked for spinning external states [R].

## 3. CHMRS 2106.10274 and gap-scaling results for t₄ in d = 3

**CHMRS (scalar ⟨φφφφ⟩ only) [R].** Large-impact-parameter strategy (§5.3, Fig. 14, (5.28)):
* bulk point: positive by the flat limit;
* fixed m, `β ≫ 1/M`: the *improvement* term (forward k ≥ 4 derivatives, growing like `P′_J(1) ∝ J²`) dominates and is positive. Its
  dominance is kept at all β because `(−½∂_{ν²}P_{(2−d)/2+iν}(cosh β))/P_{(2−d)/2}(cosh β)` is monotone in β ((5.33)–(5.35));
* β fixed, `m → ∞`: the *unimproved* k = 2 part reduces to `2f₀ Q_n(β)/M^{n+1}` (5.42). `Q_n > 0` was checked numerically for
  `0 ≤ n ≤ 1` in d = 3 (below (5.44)).

AdS4 (§5.5): the d = 3 functional (5.55) has `f → 1310p − 1080p²`. The `a₁p` term integrates to 0 in flat space (5.57). In AdS it
gives `a₁Q₁(β)/M²`, everywhere positive (Fig. 18). "This is the crucial distinction between AdS4 and the momentum cutoff in (5.59): the
latter goes negative at large b, but the former stays positive!" Result: `g₂ ≥ −17.64(8πG/M²) log(0.124 M R_AdS)` (5.64). For the
project, the diagonal `A₁M₁(β)/Λ²` term (certified positive) is the spinning analogue of `a₁Q₁`. The off-diagonal has no scalar analogue [I].

**Gap scaling of t₄ (the ⟨TTT⟩ structure dual to W³/Riem³) in d = 3.** All of these are parametric. None has a sharp constant.

| paper | statement | d = 3? | assumptions / status |
|---|---|---|---|
| CEMZ 1407.5597 [R] | `(α₄)^{1/4} ≲ 1/Δ_gap` (5.19), so `t₄ ∼ α₄ ≲ Δ_gap⁻⁴`; in D = 4 the Einstein log time delay "L → R_AdS4 … not a real issue" (§3.5) | yes (bulk D = 4 discussed) | eikonal/time-delay thought experiment; `≲` with unknown constant |
| Costa–Hansen–Penedones 1707.07689 [R] | `C⁽³⁾_TTT ∼ 1/Δ_g⁴` (4.36) (`C⁽²⁾ ∼ Δ_g⁻²`) | general `d = 2h` | assumes OPE data along the leading Regge trajectory is analytic in `D = 1/Δ_g²` with O(1) coefficients (4.32)–(4.34) |
| Meltzer–Perlmutter 1712.04861 [R] | quotes `α_n ≲ Δ_gap⁻ⁿ` (1.3); d = 3 parity-odd `C̃_TTT ∼ Δ_gap⁻⁴` (App. C.4, (C.35)); notes the d = 3 log, `Π ∼ −log L/2π` (3.13), "left understood" (fn. 9) | yes | Regge-limit matrix unitarity `D₁₁D₂₂ ≥ D₁₂²` in impact-parameter space ((3.4)–(3.6), (3.14)) |
| Afkhami-Jeddi–Kundu–Tajdini 1805.07393 [R] | **d = 3:** `|ñ_s − 4ñ_f|/c_T ≲ ln Δ_gap/Δ_gap⁴`, `|ñ_odd|/c_T ≲ ln Δ_gap/Δ_gap⁴` (7.7); bulk `α₄, α̃₄ ∼ ln Δ_gap/Δ_gap⁴` for `W³, W̃W²` (7.8); "logarithmic violation" special to d = 3 (§7, (2.28)) | **yes, explicit** | holographic null energy condition; assumes OPE coefficients in δG small (fn. 13) |
| AHKT 1610.09378 [A] | ⟨TTT⟩ → Einstein, a ≈ c | d = 4 only | chaos bound + assumptions |
| KKSZ 1904.05905 [R, intro + §6.1] | 4d sum rules `(t₄+2t₂)² = Σ|λ_TTφ|²(…) + …` (1.6); argues heavy contributions are gap-suppressed | d = 4 | "non-rigorously" (§1), "no outliers" assumption (§6.1) |
| CHLPS 2201.06602 [R] | flat D = 4: `|ĝ₃|²M⁸ ≤ 24.9 log(M/m_IR) − 27.6` (4.4); claimed to translate to AdS4 with `log(M R_AdS)` (§5) | expected, not derived | S-matrix axioms; IR cutoff |

**What Result D must reproduce or beat:** `|t₄| ≲ C √(A log Δ_g + B)/Δ_g⁴`. This is the √log form, because the light side is quadratic
in ĝ₃ / t₄ (2201.06602 (2.31a): `−B⁽¹⁾₂|_low = 16πG/p² + 2πG|ĝ₃|²p⁶`), as in the project's `240√(A log Λ + B_AdS)/Λ⁴`. That is
already parametrically better than AKT's `ln Δ/Δ⁴`, and C is the new content. [U]: I did not convert the project's `A = 3022` to
CHLPS's 24.9 (normalisations of ĝ₃ vs t₄ and of `M ↔ Λ` differ, and NF02 is not the optimal functional).

## 4. Techniques relevant to a curvature-induced off-diagonal at large b

1. **A second k = 2 rule whose heavy action is purely off-diagonal (B⁽³⁾₂-type).** In flat D = 4 this is the only source of a
   tunable `c^{++}c^{+−}` term (2201.06602 (2.35b)). Its b-profile is `∝ ∫dp p⁵ J₄(bp) B⁽³⁾₂` (3.23). [I] In CFT language it would come
   from ⟨TTTT⟩ with one helicity-flipped external polarization, i.e. a "non-shock" rule. CLS-D §5 propose building it from a conformally
   covariant differential operator `𝒟 L[O₁]L[O₃]`. It is not constructed anywhere. Caveats [I]: (a) its light side contains ĝ₃ *linearly*
   (3.25), which is fine for bounding |t₄| but changes the optimisation; (b) the flat kernel decays at large b (J₄), so it can only cancel
   a b-independent AdS off-diagonal if its own AdS curvature part has the right shape. That requires the same crystal-plus-curvature
   computation as for Ĉ₂.
2. **CEMZ-type PSD matrix of sum rules** (2201.06602 §3.5). Make the heavy side a Gram matrix in external helicity superpositions,
   so per-state PSD is automatic and positivity moves to the light-side 2×2 (3.25). This is the cleanest literature template for "a
   PSD functional mixing several sum rules". In AdS it again needs the B⁽³⁾-type rule (item 1).
3. **Improvement terms (forward k ≥ 4 derivatives).** At fixed m and large β they dominate because `P′_J(1) ∝ J²` grow with J, and
   they stay dominant at all β by the monotonicity (5.35) (2106.10274 §5.3.1). [I] Their action at fixed b̂ scales like
   `b̂²(Λ/m)⁴`-ish relative to k = 2 [U]. So they can cover the band at `m ≈ Λ` but fade at `m ≫ Λ`, where the project finds the failure
   persists. They also bring back the contact/forward terms (project's H3, H9b, H15). CHLPS also add `B⁽¹⁾₄` forward derivatives,
   (3.3), (4.4).
4. **Matrix-valued large-b asymptotics** (2201.06602 App. C, (C.11)–(C.14)). Once a family with a tunable off-diagonal exists, this
   is the standard way to impose `2×2` PSD on the whole tail `b ≥ B` with a free phase, via a matrix-polynomial SDP. Large-J at fixed m
   needs its own polynomial-in-J constraint (C.3).
5. **Exact helicity basis in d = 3** (Caron-Huot–Li 2102.08160 [R]). Helicity `h` (2.19)–(2.26) commutes with conformal maps (shown
   for spin 1 in (2.22); conjectured for higher spin). Helicity-eigen 3-pt structures (2.10), (2.14) make MFT OPE data diagonal. **But**
   "taking discontinuities or commutators do not preserve h eigenstates" (end of §2.2). The light-transform matrix is off-diagonal in
   the Even/Odd basis, `L_E ∝ [[0,1],[1,0]]` (2.57), and "not only non-diagonal, but also dense" in the helicity basis (§2.6). In AdS4
   YM, Einstein-like × Einstein-like exchange gives helicity-diagonal anomalous dimensions to all orders in 1/β ((4.7), (4.9)–(4.11)).
   Mixed minimal × non-minimal vertices are purely off-diagonal (4.12). [I] This is the structural reason a light-ray-commutator
   functional (CLS-D) can mix helicity at finite curvature while being diagonal at the bulk point. A congruence to this basis cannot
   remove a PSD failure, because PSD is basis-invariant. It would make the order-0 off-diagonal *exactly* zero (cf. the project's
   residual ≤ 0.24 near zeros of `D_h`), isolating the genuine mixing. Meltzer–Perlmutter's matrix positivity uses the same logic:
   "every extra derivative in D₁₂, relative to the diagonal terms, must have a coefficient with one more power of Δ_gap⁻¹" ((3.14)).
   [I] That is a possible way to argue that the off-diagonal *couplings* `λ₊₊λ₊₋` of the states populating the band are suppressed, i.e.
   a non-per-state argument.
6. **Conditional fallback (what the D = 4 literature actually does).** State the bound as rigorous modulo an assumption on the
   heavy spectrum in the band `b̂ ∈ [~0.65Λ^{2/3}, O(Λ)]`, as CHLPS/CHMRS do for their `[m_IR^{-2/3}, m_IR^{-1}]` band
   (2106.10274, below (5.60)). [I] Physically the band is where eikonal/Regge physics is expected to dominate (Häring–Zhiboedov
   2410.21499 [A]: "At large impact parameters … the eikonal model … as a universal mechanism that generates the graviton pole"). For
   Einstein gravity the eikonal phase is helicity-diagonal (2201.06602 (3.26): off-diagonal ∝ ĝ₃).

## 5. Most useful for getting around the obstruction (ranked)

1. **The failure is the AdS version of the known D = 4 IR band** (`b ≳ m_IR^{-2/3}` ↔ `b̂ ≳ Λ^{2/3}`; 2201.06602 (3.6), 2106.10274
   (5.59)), and the weight-independence follows from the D = 4 Fourier-positivity theorem `ψ′(0) > 0` (2106.10274 (5.54)). So no Ĉ₂/Ĉ₃
   reweighting can work. This agrees with `LARGEJ_OFFDIAG.md` §3, and the obstruction is not a numerical artefact.
2. **The literature's way out in flat space was to accept the band** and make the bound conditional on large-b spectral data. For Result
   D the honest analogue is to state the bound conditional on the `(++,+−)` couplings in the band (§4.6), quantified by the computed
   `λ/N_s`.
3. **The only literature-backed functional fix is a helicity-off-diagonal k = 2 sum rule (B⁽³⁾₂ uplift) or the CEMZ matrix of sum
   rules** (2201.06602 (2.35b), §3.5). Both need CFT sum rules beyond CLS-D's shock family (CLS-D §5 open problem). That is research,
   not bookkeeping.
4. **Improvement (k = 4 forward) terms** can dominate at fixed m and large b with monotone AdS continuation (2106.10274 (5.33)–(5.35)).
   They are a partial fix near `m ≈ Λ` and reintroduce contacts.
5. **For the SDP once a tunable off-diagonal exists:** use CHLPS's matrix-polynomial large-b PSD with a free phase ((C.11)–(C.14)) and
   explicit large-J control (C.3). Optionally use the Caron-Huot–Li helicity basis to make the order-0 matrix exactly diagonal.

## 6. Not verified / not done

* Whether CLS-D's 7+3 rules reduce exactly to B⁽¹⁾ in D = 4 [I; consistent with project H0 results, not checked in the paper].
* The AdS curvature profile of any B⁽³⁾₂-type CFT rule: no source.
* Body of AHKT 1610.09378 and of "Shockwaves from the OPE" 1709.03597: abstract only (both are d = 4 / general ⟨TTT⟩ consistency, no d = 3
  sharp bound found).
* Citation search: INSPIRE "refersto" for 2311.04271 (15 papers), 2201.06602 (206), 2106.10274 (141) and 2102.08160 (56), filtered by
  title. Relevant-looking recent items (2610.19149 "IR-finite positivity bounds in 4D gravity" (half-space, scalars), 2501.17949
  "Graviton loops and negativity", 2407.02346 "Non-forward UV/IR relations", 2604.15235 "Sampling the graviton pole") were checked at
  abstract level only. None treats spinning AdS4 functionals or helicity mixing.

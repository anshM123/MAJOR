# `c₁`: what the AdS shock / conformal-Regge route gives, what the contact-like part really is, and the conditional value `c₁ = 0`

*23 Nov 2026, research agent. Unreviewed, not committed. No existing file was edited. New code: `ads4/code/c1val/` (outputs in
`code/c1val/cache/`, git-ignored except its `.gitignore`). Literature read from arXiv LaTeX source in the session scratchpad: CHP
1707.07689 (§2.5, eqs. `DTensors_st`, `DH3_grav`, `cases`, `tensors_n`, App. `t2_and_t4`); Meltzer–Perlmutter 1712.04861 (eq.
`TTdiagonal`); Maldacena–Pimentel 1104.2846 ("On graviton non-gaussianities during inflation", eqs. `Wsp`, `ppp`, `newform`); Jain–John–
Mehta–Nizami–Suresh 2106.00016 ("Higher spin 3-point functions in 3d CFT using spinor-helicity variables", JHEP 09 (2021) 041; eqs.
`JsJstwopointfn`, `evenTTTso`, `tttoddsphresult`, `TTTshal`, `tttsph1`, §`contact`). CHP and MP are in `prl/refs_verified.bib`; the
other two are not.*

**Grades.** [P] proved (incl. exact algebra); [D] derived; [D-sketch]; [DV] derived and checked numerically/symbolically; [Q] quoted;
[NE]; [M] model / natural size; [A] assumption; [I] open.

---

## 0. Summary

| item | result | grade |
|---|---|---|
| **Value of `c₁`** | **`c₁ = 0` under premise F₂** (§4.2): the `k = 2` analogue of H7's factorisation F, i.e. the `s`-channel `T`-family term of the `k = 2` functional factorises through light-ray one-point matrices, with the Einstein one-point matrix diagonal and the `t₄` one-point matrix off-diagonal in helicity. **Without F₂ the value is not determined.** | `c₁ = 0`: **[D-sketch given F₂]**; F₂: **[A]**. Its `k = 3` analogue is [D] (H7 §4). Unconditional value: **[I]** |
| Rigorous bound on `|c₁|` | none | [I] |
| **Suggested route: AdS shock / CRT phase shift** (CHP impact-parameter operators `D̂₁, D̂₂, D̂₃` on `H₂`) | In `d = 3` the preserved (`ε ε ε* ε*`) element of the `t₄`-type operator is `¼□(□ + ν² + ¼)`, and that of the `t₂`-type operator is `½(□ + ν² + ¼)`. **Both vanish identically on every `H₂` eigenfunction, for every `ν` and `L`.** So "`ν²P(ν²)/(ν²+9/4)` = polynomial + pole" has **zero** polynomial remainder in the preserved sector. This route gives `c₁^{(U)} = 0`, as C1_COMPUTE §3.3 found, and **cannot produce a non-zero `c₁`.** | identity **[P]** (exact sympy, `v1`); applied to CHP's basis **[Q + P]** |
| Can a shock computation see `c₁` at all? | **No.** A `ν`-independent term in `β(ν)` is a `δ`-function at zero impact parameter on `H₂`. Shock / phase-shift computations at `L > 0` fix `β(ν)` only modulo polynomials in `ν²`. | [D] |
| **ASD-closure lemma (new)** | In the bulk Weyl-type configuration (legs 1, 2 self-dual; legs 3, 4 anti-self-dual), which is the preserved ⊗ preserved one in C1_COMPUTE's dictionary, the whole `O(ĝ₃)` tree contribution (the `s`- and `t`-channel `E × W³` exchanges plus the quartic seagull) collapses **at integrand level** to `3ĝ₃∫√ḡ Tr(Ψ₁Ψ₂Φ₃₄)`. Here `Φ₃₄` is a **free** regular self-dual spin-2 field, generated only by reflection at the AdS boundary. **There is no independent zero-impact-parameter seagull.** The remnant is an energy-channel (`s`/`t`) single-trace-`T` factorised term, `⟨T₁⁺T₂⁺T⁺⟩_{W₊} · [⟨T⁺T₃⁻T₄⁻⟩_E − κ⟨T⁺T₃⁻T₄⁻⟩_{CS}] / ⟨T⁺T⁺⟩`, plus its parity mirror. It has no total-energy singularity, so its flat limit vanishes. | **[D]**. Inputs: SD/ASD of AdS₄ helicity modes [DV, `v2`]; existence of a regular ASD second-order solution [D-sketch]; boundary relation checked against JJ21 [Q + DV, `v3`] |
| Consistency check of the lemma's boundary analysis | It predicts that on two ASD legs the Einstein `⟨TTT⟩` equals `κ ×` (gravitational Chern–Simons contact) **exactly in the all-same-label slot**, and differs in the mixed slot. JJ21's published `⟨T⁻T⁻⟩`, `⟨T⁻T⁻T⁻⟩_nh` share the factor `c_T + i c_T′`. The mixed slot differs by `c_T(X−Y)`, with `X − Y = −2k₃(k₁²+4k₁k₂+k₂²−k₃²)/(k₁²k₂²E²) ≠ 0`. | [Q + DV] (`v3`) |
| `c₁` is universal | It is fixed entirely by `⟨TT⟩` and `⟨TTT⟩` data (`c_T`, `t₄`): no bulk quartic coupling enters. This confirms C1_COMPUTE's universality claim by a different route. | [D] |
| Reconciliation of the two pictures | The shock picture (`u`-channel) contributes exactly 0. C1_COMPUTE's "contact-like part (C)" is not a contact: it is the energy-channel `T` exchange above. That is exactly the object the term-wise functional `Ω_OPE` sees (the `s`-channel `T` block, `E × W` structure). | [D] |
| Flat limit / HM parametrisation in `d = 3` / `EE` projection | Flat limit 0 [P: integrand-level for plane waves; no total-energy pole]. HM: in `d = 3` with circular `ε` both CHP collider structures (`t₂`, `t₄`) vanish, and the `t₂` structure is absent identically in the preserved sector [P, `v1` (C), (E)]. `EE`: the `u`-channel `t₄` part is pure flip (killed). Under F₂ the energy-channel part is a dual-line flip (killed, [P] by C1 §3.2). | as stated |
| **Effect on Result D** (`Δ_g = 300 / 1000`) | `c₁ = 0`: headline unchanged, `√(t₄²+α_T²) ≤ 3.132·10⁻⁶ / 2.891·10⁻⁸`. For any `c₁`: `R₀(x+√(1+x²))`, `x = 2.160·10⁻³|c₁| / 1.670·10⁻⁴|c₁|`. Moves by `< 1 %` iff `|c₁| < 4.61 / 59.6`. | [D] (C1 §4.1 formula) |
| **Effect on Result Q** | `c₁ = 0`: upper end `2L_AdS = 72.04 / 92.81`, unchanged. For any `c₁`: `+6.81·10⁻⁵c₁² / +5.44·10⁻⁷c₁²`. Moves by `< 1 %` iff `|c₁| < 103 / 1306`. | [D given κ > 0] (C1 §4.2) |

**Bottom line.**
* The suggested shock / conformal-Regge route is **provably blind** to `c₁`: the `u`-channel part is identically zero in `d = 3`.
* What is left is an energy-channel single-trace `T` exchange, not a seagull.
* `c₁ = 0` holds if the `k = 2` functional sees that exchange only through light-ray one-point matrices (F₂).
* Otherwise `c₁` stays an uncomputed universal constant. C1_COMPUTE's natural size (`≲ 0.1`, [M]) moves either bound by `≤ 0.022 %`.

---

## 1. What is being computed (recap of C1_COMPUTE §§1–3)

* **Object.** `c₁` is defined by `−Ĉ′_{EE,ν} ⊃ 16πG c₁ t₄`, in H7's prescription:
  * the Lemma-R residue `[𝒦₊]_{r⁻¹}` of the full tree `⟨TTTT⟩`;
  * the D3 plane-wave map;
  * E5 flat normalisation, so that `−Ĉ′_EE|_E = 16πG/(ν²+9/4)`.
* **Normalisations.** `8πG = 256/c_T` and `√(t₄²+α_T²) = 240|ĝ₃|`. `t₄` is the `d = 3` HM parameter.
* **Already shown in C1_COMPUTE:**
  * the `k = 3` part is 0;
  * the `u`-channel `T`-exchange part is 0 (its `t₄` vertex is pure flip and `EE` kills flips);
  * what is left, part (C), is `ν`-independent.
* **The question here:** the value of (C).

---

## 2. The suggested route: AdS shock / conformal-Regge phase shift (CHP basis) — it gives zero

### 2.1 The CHP operators [Q]
CHP 1707.07689, eqs. `calB`, `DH3_grav`, `cases`, write the impact-parameter amplitude of `⟨TT𝒪𝒪⟩` as `Σ_k β_k(ν)D̂_kΩ_{iν}(L)`. Conservation
leaves `k = 1, 2, 3`:
* `D̂₁ = δδ`;
* `D̂₂ = δ^{{μ}}_{{ρ}}∇^{ν}}∇_{σ}} + (ν²+((d−2)/2)²)/(d−1)·δδ`;
* `D̂₃ = ½(∇^μ∇_ρ∇^ν∇_σ + ∇_ρ∇^μ∇_σ∇^ν) + (4ν²+d²−2d−2)/(2(d+1))·δ∇∇`.

The polarisations are contracted as `ε_μ ε_ν ε*^ρ ε*^σ`, with `z₁ = z₂* = ε` (CHP: "only in this configuration … the norm of a state"). This is
the **helicity-preserving** element. Their large-`L` form (eq. `tensors_n`) is the HM structure
`1 + t₂(|n·ε|²/|ε|² − 1/(d−1)) + t₄(|n·ε|⁴/|ε|⁴ − 2/(d²−1))`.

### 2.2 Exact identities on `H₂` [P] (`v1_chp_d3.py`, sympy, arbitrary `f(x,y)` on the upper half plane)
* `ε^aε*^bε^cε*^d∇_a∇_b∇_c∇_d f = ¼□²f` exactly, and the same for the `(ε*,ε,ε*,ε)` ordering.
  * Reason: `ε^aε*^b = ½(g^{ab} + iω^{ab})`, and every `ω` contraction gives `ω^{ab}∇_a∇_b(scalar) = 0`.
  * Curvature terms cancel; they do not merely shift the coefficients.
* `ε^aε*^b∇_a∇_b f = ½□f` exactly.
* Hence, for `d = 3`:
  * `diag D̂₃ = ¼□(□ + ν² + ¼)`;
  * `diag D̂₂ = ½(□ + ν² + ¼)`.
* Both vanish on every eigenfunction with `□ = −(ν²+¼)`: `Ω_{iν}`, `Π_{iν}`, and the plane waves `y^{½+iν}` (checked explicitly), for
  every `ν` and every `L`.
* **Control:** the flip element `εεεε∇⁴f` is non-zero, `−4.17 + 4.52i` at a test point.
* **Collider form (E):** for circular `ε`, `|n·ε|² = ½`, so both CHP collider structures vanish in `d = 3`.

### 2.3 Consequence
* In the preserved sector only `β₁(ν)D̂₁` survives.
* **The "polynomial remainder of the `t₄` operator over the pole" is identically zero.** It is not a constant to be fixed.
* Whether `β₁(ν)` itself acquires an off-pole `t₄` dependence is a separate statement about the bulk vertex. It is excluded by C1_COMPUTE
  §3.1(b) (`Tr W³ = Tr W₊³ + Tr W₋³`; the cubic vertex needs both legs of a line to have the same Weyl type, i.e. flip). The CHP identity
  is consistent with this: the `d = 3` structures that carry `t₄` have no preserved projection at all.
* **So the route returns `c₁^{(U)} = 0` [Q + P + C1's D].**

### 2.4 Why no shock computation can return `c₁` [D]
* `β(ν) ∋ const` is, by harmonic analysis on `H₂`, a `δ`-function at `L = 0`.
* Shock / phase-shift computations are done at `L > 0` (AHKT, KPZ, CHP, KKSDZ §3.3). They determine `β` modulo polynomials in `ν²`.
* `c₁` is by definition such a polynomial: C1_COMPUTE's zero-impact-parameter remnant. It is therefore invisible to the suggested route,
  even in principle.

---

## 3. The contact-like part (C): ASD closure

### 3.1 Setting
* Euclidean AdS₄ (`R = 1`), with action `S_EH,Λ + ĝ₃I₊ + ĝ₃*I₋` and `I_± = ∫√g Tr W_±³`.
* `Tr W³ = Tr W₊³ + Tr W₋³` holds identically in 4D [P]. `t₄ ∝ Re ĝ₃` and `α_T ∝ Im ĝ₃` (DICTIONARY_T4).
* `Ψ := W₊` is the self-dual (SD) Weyl spinor.
* At `O(ĝ₃)` the generating functional changes by `ĝ₃I₊[g_E[σ]] + c.c.`: the first-order change of an on-shell action is the
  perturbation evaluated on the unperturbed (Einstein) solution `g_E[σ]`.
  * This uses that `I_±` needs no boundary counterterm: the Weyl tensor of an asymptotically AdS₄ metric is `O(z³)` in an orthonormal
    frame [D-sketch].
  * It includes every `O(ĝ₃)` tree diagram (`W³` exchanges in all channels and the quartic `W^{(1)}W^{(1)}W^{(2)}` seagull) automatically.

### 3.2 Lemma (ASD closure) [D]
**Statement.** Let the external legs be linearised solutions with `Ψ^{(1)}[h₃] = Ψ^{(1)}[h₄] = 0` (ASD) and `Ψ^{(1)}[h₁], Ψ^{(1)}[h₂] ≠ 0`
(SD). Then the four-leg term of `I₊[g_E]` is exactly

  `I₊^{(4)} = 3∫√ḡ Tr(Ψ^{(1)}[h₁] Ψ^{(1)}[h₂] Φ₃₄)`,   with `Φ₃₄ := Ψ^{(2)}[g_E]|_{h₃h₄}`,

and `Φ₃₄` is a **free, regular, self-dual massless spin-2 field** on AdS₄.

**Proof sketch.**
1. **Only one partition survives.** Every term of `√g TrΨ³` at fourth order has the form `Ψ^{(1)}Ψ^{(1)}Ψ^{(2)}`, or `√g^{(1)}(Ψ^{(1)})³`
   from three distinct legs. With only legs 1 and 2 SD, the only surviving partition is `{1,2}{3,4}`.
   * `Φ₃₄` includes the second-order Einstein field `δh₃₄` (the `s`-channel exchange) and the second-order part of `Ψ` built from
     `h₃h₄` (the seagull).
   * The frame dependence of `Ψ^{(2)}` cancels in the invariant `TrΨ³`.
2. **`Φ₃₄` is free.** The contracted Bianchi identity `∇^A_{A′}Ψ_{ABCD} = 0` holds on every Einstein space. At second order its source
   terms are `∇^{(1)}[h₃]Ψ^{(1)}[h₄] + (3↔4) = 0`. Hence `∇̄Φ₃₄ = 0`: Penrose's spin-2 equation, conformally invariant.
3. **ASD closure.** ASD Einstein metrics with `Λ` form a consistent sub-sector (Plebański / nonlinear graviton with `Λ`). Assume there is a
   regular second-order solution `δh^{ASD}₃₄` with vanishing SD Weyl [D-sketch].
   * Then `χ := δh^W₃₄ − δh^{ASD}₃₄` is a regular solution of the linearised Einstein equations, with boundary source `−b₃₄`. Here `b₃₄` is
     the boundary metric that the ASD solution needs and `δh^W` is the normalisable Witten field.
   * So `Φ₃₄ = Ψ^{(1)}[χ]`.
   * By `v2` (below), the regular modes of 3d label `+` are SD and those of label `−` are ASD. So **`Φ₃₄ = −Σ_K b₊(K)Ψ[regular label-+ mode at
     K]`**, with `K = k₃+k₄`.
   * In flat space the same argument gives `Φ₃₄ ≡ 0` for plane waves: a free massless field cannot carry the non-null momentum `p₃+p₄`. This
     is the integrand-level reason that no `O(ĝ₃)` MHV amplitude exists. In AdS, **`Φ₃₄` exists only through reflection at the boundary**.
4. **Fixing `b₊` (boundary relation).** For an ASD solution the holographic stress tensor equals `κ·Cotton[γ]` (de Haro–Leigh–Petkou type).
   At linear order this is `G = κC` on ASD modes.
   * Linear algebra on regular and irregular modes (`v2`: the irregular `e^{+kz}` mode of each label has the **opposite** duality) gives
     `τ^W − κC^{(2)}[3,4] = (κC − G)b`.
   * This is **zero in the label of the ASD legs** and `−2G b₊` in the other label.
   * Hence `b₊ = −[⟨T⁺_{k₁+k₂}T₃⁻T₄⁻⟩_E − κ⟨T⁺_{k₁+k₂}T₃⁻T₄⁻⟩_{CS}]/(2⟨T⁺T⁺⟩)`.

**Result (momentum space, 3d helicity labels, up to terms analytic in some momenta, i.e. contact terms; overall constant `𝒩` not fixed
here):**

  `⟨T₁⁺T₂⁺T₃⁻T₄⁻⟩|_{O(ĝ₃)} = 𝒩 ĝ₃ ⟨T₁⁺T₂⁺T⁺_{k₃+k₄}⟩_{W₊} · Δ(k₁+k₂; 3,4) / ⟨T⁺T⁺⟩_{k₁+k₂} + (1,2 ↔ 3,4; W₋, conj.)`,

where `Δ := ⟨T⁺T⁻T⁻⟩_E − κ⟨T⁺T⁻T⁻⟩_{CS}`, i.e. Einstein minus the ASD-normalised gravitational Chern–Simons contact. The other preserved
arrangement (1, 4 SD) gives the same with the `t`-channel.

### 3.3 Checks of the lemma's inputs
* **`v2` [DV]** (exact sympy, Euclidean Poincaré AdS₄, TT radial gauge):
  * the regular graviton modes `(1+kz)e^{−kz}` of 3d helicity `±` solve the linearised Einstein equation and have purely SD / ASD Weyl
    tensor;
  * the irregular `(1−kz)e^{+kz}` modes have the opposite duality.

  This is the basis of MP's "`W³` gives only `+++`/`−−−`" and of step 4.
* **`v3` [Q + DV]** (JJ21's published two- and three-point functions):
  * `⟨T⁻T⁻⟩ ∝ (c_T + ic_T′)` and the non-homogeneous `⟨T⁻T⁻T⁻⟩ ∝ (c_T + ic_T′)F`. So subtracting `κ·CS` with the ASD normalisation
    (`c_T′ → c_T′ + ic_T`) removes both. This is exactly the predicted zero in the all-same-label slot.
    * Corollary: the Einstein all-minus component (MP's "non-flat" component) is a pure contact term on this sector.
    * C1_COMPUTE §3.5's analogy ("non-flat `+++` of Einstein `⟨TTT⟩`") therefore points at a piece that drops out.
  * The mixed slot keeps `c_T(X − Y)`, with `X − Y = −2k₃(k₁²+4k₁k₂+k₂²−k₃²)/(k₁²k₂²E²)` (`k₃` the `+` leg). It carries the flat-space
    `E⁻²` pole of `X` (`E²X → 4/k₁ + 4/k₂`), while `Y` has none.
  * So `b₊ ≠ 0`: **the remnant is generically non-zero in AdS₄.**
  * JJ21's displayed `TTTshal` line for the odd `−−+` term is dimensionally inconsistent. The version after their eq. `tttsph1` was used.
* **Flat limit [P].** The remnant is a product of three-point functions with partial-energy singularities only (`E_L⁻⁶` from `W₊`, `E_R⁻²`
  from `Δ`). It has **no total-energy singularity**, so its flat-space limit vanishes, as it must.

### 3.4 What the lemma settles and what it does not
**Settles [D]:**
* there is no independent quartic-seagull or zero-impact-parameter contact contribution;
* the `O(t₄)` preserved ⊗ preserved remnant is a single-trace `T` exchange in the **energy** channels (`s`, `t`);
* it is built only from `⟨TT⟩`, `⟨TTT⟩_E`, `⟨TTT⟩_{W₊}` and the CS contact. So **`c₁` is universal** (C1's claim, now by a different
  route).
* **Reconciliation.**
  * The shock picture sees only the `u` channel, and gives 0 (§2).
  * C1's "contact-like (C)" is this energy-channel exchange.
  * In the term-wise functional `Ω_OPE` (LIGHTSIDE_DIST) it is precisely the `E × W` term of the `s`-channel `T` block. At `k = 2`
    double traces have double zeros (CLS-D via EE_FOUNDATIONS §1.4), so nothing else enters.

**Does not settle:** whether this energy-channel term carries a Regge-spin-2 top coefficient, i.e. a non-zero Lemma-R residue.

---

## 4. The value

### 4.1 Unconditional: [I]
* `c₁` is the `EE`-projected Lemma-R residue of the §3.2 remnant, equivalently its `Ω_OPE` image (LIGHTSIDE_DIST: equal on tree terms
  [D]).
* Neither a closed form nor a rigorous bound was obtained. The natural size `|c₁| ≲ 10⁻³–10⁻¹` of C1_COMPUTE §3.5 [M] is unchanged.
  The remnant has the same `W³ × Einstein / c_T` counting.

### 4.2 Conditional: `c₁ = 0` [D-sketch given F₂]
**Premise F₂ [A].** The `s`-channel `T`-family term of the `k = 2` CLS-D functional (`f = v′ − u′`) factorises as
`Σ_h ⟨T₄|A₃|T,h⟩⟨T,h|A₁|T₂⟩`. Here the `A_i` are the (weighted) light-ray integrals of `T₁`, `T₃`, and the matrix elements in the
light-transform helicity basis are **diagonal for Einstein and purely off-diagonal (`Δh = ±4`) at `O(t₄)`**.
* For `k = 3` (`f = −1`, bare light transforms) this is H7's F [D] together with Q8 (HM one-point matrix `𝟙 + (t₄/4)·flip`, exact [Q]).
* For `k = 2` the weight `f` turns the light transforms into first moments `∫dα α T`. Their matrix elements between momentum eigenstates
  are momentum derivatives of the light-transform ones. The derivative of a helicity basis is helicity-diagonal (a phase connection; the
  longitudinal parts drop for conserved `T`), so the diagonal / off-diagonal pattern survives [D-sketch]. This is not written out.

**Argument.**
1. Given F₂, the `O(t₄)` part of the `T`-family term is `𝟙·flip + flip·𝟙`.
2. That is a **dual-line (`2,4`) flip**: out-state helicity ≠ in-state helicity.
3. `EE` annihilates every dual-line flip exactly, since `(z₂′·z₄′)² = (ε_h·ε_h)² = 0` [P] (C1_COMPUTE §3.2, EE_FOUNDATIONS §2.1).
4. At `O(1/c_T)` the `k = 2` light side of `Ω_OPE` is the `T` family alone (double traces have double zeros), so `c₁^{Ω_OPE} = 0`.
5. LIGHTSIDE_DIST's `Ω_k[X] = D_LR[X] + 𝒜_k[X]` (D_LR = 0 for microcausal terms) then transfers this to H7's Lemma-R prescription, and
   **`c₁ = 0`**.

**Consistency with §3.**
* The Euclidean remnant of §3.2 is non-zero as a function, but this does not contradict `c₁ = 0`.
* `c₁` measures only the remnant's Regge-spin-2 top coefficient. A non-zero, Regge-subcritical energy-channel term has zero light value
  (LIGHTSIDE_DIST, "Lemma A′").
* Independent heuristic pointing the same way [D-sketch]:
  * a non-zero `c₁` is a degree-2 Mellin polynomial in the preserved sector;
  * momentum-space contact diagrams always carry total-energy singularities, and the remnant has none.
  * Gap: pole and polynomial parts could cancel each other's `E_T` singularities. The spinning Mellin ↔ momentum dictionary is not
    established.

**Where it could fail.**
* F₂ is the single load-bearing assumption.
* If the `k = 2` weight produced a helicity-diagonal `O(t₄)` one-point element, `c₁` would equal the `EE` image of that element on the
  `T` block. That is computable from three-point data alone (§3.4), but was not computed here.
* An unresolved bookkeeping point. In the bulk Weyl-type dictionary of C1_COMPUTE §2 (preserved line = one SD, one ASD leg), the §3.2
  remnant sits in preserved ⊗ preserved. In the light-ray state language of F₂, the `T`-family `O(t₄)` term is a dual-line flip.
  * Both statements are compatible only if the remnant is Regge-subcritical, which is what `c₁ = 0` asserts.
  * The relation between Euclidean 3d momentum labels and Lorentzian light-ray helicity for the intermediate `T` state was not settled
    here. It should be checked when F₂ is written out.

---

## 5. Checks requested by the task

| check | outcome | grade |
|---|---|---|
| Flat limit of the answer vanishes | Yes. `c₁^{(U)} = 0` identically. The remnant has no total-energy pole, and `Φ₃₄ ≡ 0` for flat plane waves (integrand level). `c₁ = 0` is trivially consistent. | [P] |
| `d`-dependence / HM normalisation, `t₂ ≡ 0` in `d = 3` | CHP's `ν`-dependent subtraction constants at `d = 3` are `c₂ = c₃ = ν²/2 + 1/8`. With them both preserved elements vanish identically. For circular `ε` the two collider structures are `½ − ½ = 0` and `¼ − ¼ = 0`. So `t₂` is absent from the preserved sector and the diagonal is `t₄`-free, matching the `d = 3` HM matrix `𝟙 + (t₄/4)(e^{−4iφ}σ₊ + h.c.)`. | [P] (`v1`) + [Q] |
| Does `EE` kill or keep it? | `u` channel: pure flip, killed (C1 §3.2–3.3). Energy channel: preserved ⊗ preserved in the bulk dictionary (§3), but under F₂ its `T`-family image is a dual-line flip, killed [P]. | as stated |
| Normalisation-independent statement | `c₁ = 0` (under F₂) needs no normalisation convention. The §3.2 overall constant `𝒩` is not needed for it. | — |

---

## 6. Consequences for both bounds (exact formulas of C1_COMPUTE §4; `code/c1val/v4_bounds.py`)

| `Δ_g` | `|c₁|` | D: `√(t₄²+α_T²) ≤` | D ratio | Q upper end `g̃₄Λ⁶/8πG ≤` |
|---|---|---|---|---|
| 300 | 0 (F₂) | **3.13206·10⁻⁶** | 1 | **72.0406** |
| 300 | 0.1 ([M] upper) | 3.13274·10⁻⁶ | 1.000216 | 72.0406 (+9·10⁻⁹ rel.) |
| 300 | 1 | 3.13883·10⁻⁶ | 1.00216 | 72.0407 |
| 300 | 10 | 3.20044·10⁻⁶ | 1.0218 | 72.0474 |
| 1000 | 0 (F₂) | **2.89126·10⁻⁸** | 1 | **92.8135** |
| 1000 | 0.1 | 2.89131·10⁻⁸ | 1.000017 | 92.8135 |
| 1000 | 1 | 2.89174·10⁻⁸ | 1.000167 | 92.8135 |
| 1000 | 10 | 2.89609·10⁻⁸ | 1.00167 | 92.8136 |

* D: `R₀(x+√(1+x²))`, with `x = 2.160·10⁻³|c₁|` (300) / `1.670·10⁻⁴|c₁|` (1000).
* Q: `2L_AdS + 6.81·10⁻⁵c₁²` / `5.44·10⁻⁷c₁²`, D-free, given `κ > 0`.
* **Tolerances.** The bounds move by `< 1 %` iff:

  | bound | `Δ_g = 300` | `Δ_g = 1000` |
  |---|---|---|
  | D | `|c₁| < 4.61` | `|c₁| < 59.6` |
  | Q | `|c₁| < 103` | `|c₁| < 1306` |

  At 0.1 %: D `|c₁| < 0.46 / 6.0`; Q `|c₁| < 32.5 / 413`.
* **Under F₂ both headlines are exact as printed (`c₁ = 0`).** Without F₂, the [M] natural size moves D by `≤ 0.022 % / 0.0017 %` and Q by
  `< 10⁻⁸`.

---

## 7. Suggested one-sentence wording for the letter

> "The only `O(t₄)` light term is a universal constant `c₁`, fixed by `⟨TTT⟩` data alone. Its shock-wave (`u`-channel) part vanishes
> identically in `d = 3`, and the remainder is an energy-channel stress-tensor exchange that vanishes if the `k = 2` functional sees it only
> through light-ray one-point matrices (as the `k = 3` functional does). Without that input the cubic bound changes by less than 1 % for
> `|c₁| < 4.6` (`Δ_g = 300`; `< 60` at `Δ_g = 1000`), and Result Q changes only by `+6.8·10⁻⁵c₁²`."

A shorter alternative: "… a universal constant `c₁`, which vanishes given a light-ray factorisation established for the companion `k = 3`
functional; the bounds are insensitive to it for `|c₁| ≲ 5`."

---

## 8. What in existing documents should be corrected (no file was edited)

1. **C1_COMPUTE §3.5, "Why §3.1 does not force it to zero".**
   * The quartic seagull is **not** an independent source. In the preserved ⊗ preserved Weyl-type sector, seagull plus non-factorisable
     `s`/`t` exchange cancel identically (ASD closure, §3.2). The remnant is an energy-channel single-trace exchange produced by
     boundary reflection.
   * The Maldacena–Pimentel analogy ("non-flat `+++` components of Einstein `⟨TTT⟩`") should be replaced. The all-same-helicity Einstein
     component equals the CS contact on this sector and drops out. The remnant is driven by the mixed component, Einstein − κCS.
2. **C1_COMPUTE §7 item 1, third bullet** ("factorises through shock vertices … would give `c₁ = 0`"). The factorisation that holds is
   through the **energy** channel, not through `u`-channel shock vertices. The missing step is F₂ (§4.2), not "no seagull".
3. **ADS4_LIGHT §2.3, `t₄`-dependence bullet,** and §0 item 3 (`t₄[V(ν) − V(−3i/2)]/(ν²+9/4)` "is the constant `16πG c₁t₄`"). In `d = 3`
   the `t₄`- and `t₂`-type impact-parameter structures have identically zero preserved element (CHP basis, all `ν`, `L`; §2). This
   polynomial is zero. `c₁` does not originate there (consistent with C1_COMPUTE §3.1(b)).
4. **Task premise** ("the `t₄` structure acts as a polynomial operator; the remainder is a candidate for `c₁`"). This does not apply in the
   `d = 3` preserved sector. It would in `d ≥ 4`, where linear polarisations give a non-zero `t₄` diagonal.
5. **UNIFIED_THEOREM_v2 / letter wording.** "Free parameter" should become "uncomputed universal constant, zero given F₂". The numerical
   bound tables need no change if F₂ is adopted.

---

## 9. Reproduce

All runs use one core with `ulimit -v 3000000; export OMP_NUM_THREADS=1`. Outputs go to `ads4/code/c1val/cache/`. Nothing is imported
from other project code.

```
cd ads4/code/c1val
python3 -B v1_chp_d3.py         # §2: CHP d=3 preserved elements of D2, D3 on H2 (exact)     -> cache/v1_chp_d3.txt        (~6 s)
python3 -B v2_ads4_sd_modes.py  # §3.3: SD/ASD of AdS4 helicity modes, regular vs irregular  -> cache/v2_ads4_sd_modes.txt  (~40 s)
python3 -B v3_jj21_slots.py     # §3.3: JJ21 slot structure under subtraction of kappa*CS    -> cache/v3_jj21_slots.txt    (seconds)
python3 -B v4_bounds.py         # §6: D and Q versus c1, tolerances                          -> cache/v4_bounds.txt        (instant)
```

Every script prints `ALL CHECKS PASS: True` where applicable; all did.

## 10. Grade table

| statement | grade |
|---|---|
| CHP `d = 3`: `diag D̂₃ = ¼□(□+ν²+¼)`, `diag D̂₂ = ½(□+ν²+¼)`; zero on all eigenfunctions | [P] (sympy exact) |
| Shock / CRT route gives `c₁^{(U)} = 0` and cannot see a `ν`-constant | [Q + P]; [D] |
| AdS₄ helicity modes: regular `±` are SD/ASD; irregular have the opposite duality | [DV] |
| ASD-closure lemma: `O(ĝ₃)` preserved ⊗ preserved = `3ĝ₃∫TrΨ₁Ψ₂Φ₃₄`, `Φ₃₄` free; no seagull | [D] (ASD second-order existence [D-sketch]) |
| Remnant = energy-channel `W₊ ⊗ (E − κCS)` factorised term; flat limit 0 | [D]; flat limit [P] |
| Einstein all-same-label slot = κ·CS contact; mixed slot ≠ (JJ21) | [Q + DV] |
| `c₁` universal (three-point data only) | [D] |
| `c₁ = 0` | [D-sketch given F₂]; F₂ [A] (`k = 3` analogue [D]) |
| `c₁` unconditional value; rigorous bound | [I]; none |
| natural size `|c₁| ≲ 0.1` | [M] (C1_COMPUTE, unchanged) |
| D, Q bounds versus `c₁` | [D] (C1_COMPUTE §4 formulas), [D given κ > 0] for Q |

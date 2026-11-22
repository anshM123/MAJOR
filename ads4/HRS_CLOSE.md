# Closing `H_RS` at low ν: the analytic-radial engine and the crystal cycle agree to all orders in 1/Δ (contiguity argument), with the bookkeeping checked and the truncation measured

*22 Nov 2026, H_RS-close agent. Unrefereed, not committed. No existing document or code file was edited. New code:
`ads4/code/spinfunc/rs2_bookkeep.py`, `rs2_uniform.py`, `rs2_parity.py`, `rs2_engine5.py`, `rs2_trunc.py`, `rs2_queue{A,B,C,D}.sh`.
Caches: `spinfunc/cache/rs2_*` only (no existing cache was written). Read-only imports: `bg_ser`, `bg_geom`, `bg_ee`, `bg_eeq`, `bg_qan`,
`bf_qan`, `l5_spin` (via `bg_ee`), `resultq/stage6_*`, `common_q` (via `bg_qan.ffn`). Existing caches read, never written:
`bg_eeq_D{300,1000}_J0_K4.pkl` (production engine images), `bg_qan_T_*.pkl`, `bg_validate.txt`. Compute: ≈ 1.5 CPU-h, ≤ 2 workers,
`ulimit -v 3 GB`, `OMP/OPENBLAS = 1`, ≤ 0.8 GB per process. Nothing in `code/resultq` (`ein_*`) was touched.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]** derivation sketched; **[DV]** derived and checked numerically; **[NE]** numerical
evidence; **[M]** model; **[A]** assumption; **[I]** open.
Unit: **`c₂` per harmonic** = `Γ/Σ_j|R_j(ν)|` (BELOW_FOLD2 §4.2, HRS_CHECK). `s_d = 0.04 Σ_ν χ_d(ν) Re Γ_K(ν)` (Stage 7 §2).
**Engine** = BELOW_FOLD2's analytic-radial `1/Δ` engine on the real section `D` (scaled variables `X₁ = p₁², X₂ = p₂², s, Y₁ = q₁², Y₂ = q₂²,
φ, α`; atoms `p₁, p₂, q₁, q₂, s, w = 1−s², cf, sf, ca, sa, u = 1−cos α, v = 1+cos α`). **Crystal** = the crystal-cycle integral
(n3 contour above the fold, `A₊ + A₋` thimbles below). `a = 1/Δ`.

---

## 0. Summary

| item | result | grade |
|---|---|---|
| **Route** | (a), analytic. The numerical routes (b)–(d) cannot reach 10–20 `c₂` at ν < 2 with the available crystal-node machinery (HRS_CHECK §7; error ∝ Δ, CRN collapse below ν = 2). They are also unnecessary once (a) holds: the cycle question reduces to the engine's own truncation, which is measured here (§4). | — |
| **The argument** (§2) | To every order in `a`, the engine's integrand is `Ψ_λ × (Laurent polynomial in the atoms)`. `Ψ_λ = Π atom^{λ} e^{−S}`, `λ = E₀ + iνE₁`, is a **product of one-variable twisted local systems** (4 radial, `s`, `φ`, `α`). Integration by parts (contiguity) reduces every monomial form to a rational multiple of `Ψ_λ`, with coefficients that do not depend on the cycle. The `s`-odd and `sf`-odd forms are exact. The only other cohomology class (the `α`-odd one) is not populated. Hence **any closed twisted cycle Γ of this local system gives `∫_Γ = c_Γ(ν) × (engine)` order by order, with one constant `c_Γ(ν)` common to every order, state, reference and kind (scalar, EE `k = 2`, EE `k = 4`)**. So the state/reference ratio `ρ = I_st/I_ref`, which is all Result Q uses, is cycle-independent to all orders in `1/Δ`, at every real ν including ν → 0. | **[D]** (Lemmas: [P], classical) |
| **Hypotheses of the argument** (§3) | Checked on the engine's own pieces, order by order to `a⁴` (`a⁵` for the binding state), for J = 0 (Δ = 100, 300, 1000, 2·10⁴) and J = 2 (Δ = 314.935, 2·10⁴):<br>(B1) the ν-part `E₁` of the leading exponents is identical for the scalar, EE2 and EE4, and for state and reference;<br>(B2) all exponent differences lie on the integer lattice: EE − scalar = `(p₁ +16, p₂ +8, w +8)` (even in `p`, i.e. integer in `X`); no `p`/`q`/`cf` parity violations among surviving monomials;<br>(B3) no negative power of `s`, `sf`, `ca`, `sa` (no pole at interior points or at the crystal angles);<br>(B4) the parity rules (`p₁ ≡ p₂ ≡ s`, `q₁ ≡ q₂`, `s+sf+sa` even, `cf+sa` even) hold exactly; the only "violators" are DFT round-off, `|c| ≤ 4.7·10⁻¹³` of the largest coefficient, integrated weight `≤ 4·10⁻¹⁰` even on a single α-arc;<br>(B5) the local system is non-resonant at every divisor for all real ν, ν = 0 included;<br>(B6) **no boundary layer**: every logged or inverted factor (40 of them: all scalar and EE logs, the 12 inverted distances of the EE prefactor, `u′v′`, `m₁₃m₂₄`, `−x²`) is `c·atoms^e·(1 + ε)` with `ε` bounded at every divisor (growth exponent ≥ −0.005 over all factors, orders and divisors). | **[DV]** |
| **Lemma-B `R` and the "new poles"** | On the Regge region `R = 4a⁴C²(1 + O(a⁴))` exactly (`C = X₁Y₂ + X₂Y₁ − 2x̂_pŷ_p`), a polynomial. Its poles `u′v′ = 0` are `m₁₂ = m₃₄ = 0` and `m₁₄ = m₂₃ = 0`, which sit at `a²C ≈ ±1`, i.e. at scaled radial distance `~Δ` (`e^{−O(Δ)}`). They are **not new**: the scalar `k = 2` integrand already has `1/(u′v′)` (its subtraction `f = (v′−u′)/(u′v′)`), with the same double pole at `m₁₄ = 0`. EE4 = scalar × (rational function) with no new branch divisor. | [P] (algebra) + [DV] |
| **Remaining premise RD** (§2.6) | The crystal-cycle integral equals, up to terms beyond all orders in `1/Δ`, the integral over a closed twisted cycle of the Regge-scaled local system. In other words its far region (`|raw x| ~ 1`, where `m₁₂ = 0`, `m₁₄ = 0` and the pole of `B` live) carries only `e^{−O(Δ)}` weight. Status by kind:<br>• scalar, on the same local system (B2): [DV], since engine = exact Mellin skeleton to ≤ 5·10⁻⁴ `c₂` at Δ = 300, ν ≤ 10;<br>• EE, ν ≥ 10: [NE] (HRS_CHECK: crystal = engine to 0.3–1 % of the beyond-LO content);<br>• EE, ν < 10: [D-sketch] (descent cycles avoid `m₁₅ = m₂₅ = 0`; EE/scalar is polynomially bounded in the Regge region). | scalar [DV]; EE [NE] / [D-sketch] |
| **Truncation, the only residual error** (§4) | New `K = 5` images. Binding state `++` J = 0 at the gap:<br>• the `a⁵` term is ≤ 0.0096 `c₂` per harmonic on ν < 4 at Δ_g = 300 (≤ 0.0003 at 1000);<br>• the `a⁴` term is ≤ 0.20 / 0.020;<br>• the tolerances are 16.9 `c₂` on [1,2] and 65–80 on [2,4].<br>In `s_d` the `a⁵` term is **+0.00064** (300) and **+0.00002** (1000); the `a⁴` term is +0.0152 / +0.0015. | [NE] |
| **Propagated** (§6) | `Δs_d(H_RS) = 0` to all orders. Truncation bound (size of the `a⁴` term, ×24 the measured `a⁵`): `|Δs_d| ≤ 0.015` (300) / `0.0015` (1000), i.e. `|ΔFF/N| ≤ 0.0017` (δ = 0.2) / `0.0013` (δ = 0.05).<br>FF/N with `K = 5` (H0): **+0.14265** (300) / **+0.14282** (1000) at δ = 0.2; +0.00608 / +0.00618 at δ = 0.05 (SB: +0.1412 / +0.1411; +0.0050 / +0.0049). The margin against `s_fail − s_d = 1.27` is ×85 (300) / ×850 (1000). | [D] + [NE] |
| **H_RS regrade** (§7) | **Binding `++` J = 0 (Δ_g = 300, 1000) and `++` J = 2 (1.05Λ): from [I] to [D given RD], net [D-sketch]**. That holds in every ν window, ν < 1 and ν < 2 included. The decisive change: the "×3.5 engine error" scenario of REVIEW_UNIFIED issue 3 now requires RD to fail, i.e. a non-perturbative far-region contribution of relative size ≥ 2·10⁻⁴ (300) / 2·10⁻⁵ (1000) that is moreover identical at both Δ_g (the measured `s_d` agrees to 0.002). OO J = 0: via `Γ_OO = Γ_++` (unchanged) or directly by the same argument [D-sketch]. | — |

**Bottom line.**
* The real-section ↔ crystal-cycle question that HRS_CHECK could not decide numerically below ν ≈ 6 has a structural answer. The Regge-scaled integrand is a product of one-variable hypergeometric-type local systems times polynomials. By integration by parts, every order of the `1/Δ` expansion is the leading integral times a cycle-independent rational function of the exponents. The engine's Γ/Beta functions are exactly these relations evaluated on `D`.
* What this needs about the EE `k = 4` integrand (same local system as the scalar, no new divisors in the Regge region, no boundary layers, populated cohomology one-dimensional) has been checked on the engine's pieces.
* What remains is RD, the premise that the crystal-cycle integral is Regge-dominated. It is the same premise on which the scalar engine is verified ([DV], all ν), and it can fail only non-perturbatively.

---

## 1. What is being closed

* H_RS (HRS_CHECK §1, UNIFIED_THEOREM §3.2): the engine's EE `k = 4` image equals the crystal-cycle integral beyond leading order.
* Result Q uses only the state/reference ratio `ρ(ν) = I_st/I_ref`, through `Γ_exact = m²(c T_ref ρ − T_st)` with `c` fixed at ν = 0.02.
  Then `Γ_K = Γ_exact − 0.75 Σ_j R_j(ν−ij)² − const·N̂₄`.
* HRS_CHECK tested `ρ` node by node with direct crystal nodes:
  * [NE] at ν ≥ 10;
  * [NE, weak] at 6–10;
  * uninformative at 2–6 (errors ±40 … ±1800 `c₂`, ×3 inflated);
  * untested below 2, which carries 68 % of the `s_d` weight.
* **Why no numerical route was pursued first.**
  * The crystal-node error per harmonic grows ∝ Δ.
  * The CRN correlation is lost below ν = 2.
  * The `k = 4` errors are heavy-tailed.
  * Reaching 17 `c₂` at ν ∈ [1,2] would need ×100 in precision.
  * Route (d) (lower Δ) gains at most ×6 before the engine's own truncation (∝ Δ^{2−k} in `c₂`) reaches the same size.
  * Route (c) (Gauss rules built from the weight's moments) is, by construction, blind to the cycle: Gauss rules reproduce the moments, and the moments are exactly what §2 shows to be cycle-independent. Applied to the exact integrand it would only re-measure the truncation, which §4 does more cheaply.
  * Route (b) needs an independently known control integral on the crystal cycle; the only candidates are engine values, so it is circular.

  Route (a) was therefore taken first, and it closes.

## 2. The contiguity argument

### 2.1 Structure of the engine's integrand [P from the construction; DV §3]
For a kind `κ ∈ {scalar k = 2, EE2, EE4}`, a state `(Δ, J)` and its reference, `bg_geom.assemble` writes the exact integrand on the
Regge-scaled variables as
```
Φ_κ(ξ; a, ν) = a^{d_κ} e^{l_κ(ν)} Ψ_λ(ξ) [ F_κ^{(K)}(ξ; a, ν) + O(a^{K+1}) ],   Ψ_λ = Π_atoms atom^{λ_atom} e^{−(X₁+X₂+Y₁+Y₂)},
F_κ^{(K)} = Σ_{k≤K} a^k Σ_j (iν/2)^j corr_{κ,j,k}(ξ),     λ = E₀ + iνE₁,
```
* `corr` are Laurent polynomials in the atoms; the `B`-correction is multiplied in.
* `λ` comes only from the factors raised to non-integer or ν-dependent powers (`A^{iν/2}`, measure, dual, celestial, 3-point).
* The polynomial factors (`f = (v′−u′)/(u′v′)` for the scalar; `v′ − u′`, the EE prefactor and `R` for EE) and all relative corrections
  `ε` enter `corr` with integer powers.
* The divisor atoms are `p₁, p₂, q₁, q₂` (`X_i = 0`), `w` (`s = ±1`, the `x` light cone), `cf` (`φ = ±π/2`, the `y` light cone), and
  `u, v` (`α = 0, π`).
* In the coordinates `X_i`, `s`, `t = sin φ`, `z = e^{iα}`, `Ψ_λ` is a product of seven one-variable multivalued functions:
  * radial: `X^{c}e^{−X}` on `ℂ*`;
  * `s`: `(1−s²)^{g}` on `ℂ∖{±1}`;
  * `t`: `(1−t²)^{(h−1)/2}` on `ℂ∖{±1}`;
  * `z`: `(z−1)^{2A}(z+1)^{2B}z^{−A−B−1}` on `ℙ¹∖{0, ±1, ∞}`.

### 2.2 Lemma 1 (one-variable contiguity) [P]
* For each factor, the twisted de Rham cohomology `H¹` (forms `Ψ·(Laurent polynomial)·dξ` modulo `∇`-exact ones) has dimension:
  * 1 for each radial factor and for the `s` and `t` factors;
  * 2 for `z`. It splits under `z → 1/z` (`α → −α`) into an even and an odd line.
* So every monomial form in a given parity sector is cohomologous to `R_m(λ) ×` the sector's generator. `R_m` is rational in `λ`.
  Examples:
  * radial: `X^{c+1}e^{−X}dX = (c+1)X^c e^{−X}dX − d(X^{c+1}e^{−X})`;
  * `s`: `s·w^g ds = −d(w^{g+1})/(2(g+1))`;
  * `t`: `t·(1−t²)^{γ}dt = −d((1−t²)^{γ+1})/(2(γ+1))`.
* For any closed (rapid-decay or regularised) twisted cycle γ the boundary terms vanish, so `∫_γ Ψ m = R_m ∫_γ Ψ`, **with the same
  `R_m` for every γ**. The engine's closed forms are the case γ = the real segments:
  * radial `Γ(1 + T/2)`;
  * `s`: `Γ(½)Γ(G+1)/Γ(G+3/2)`;
  * `φ`: `Γ(½)Γ((M+1)/2)/Γ(M/2+1)`;
  * `α`: `2^{A+B+1}Σ_j C(r,j)(−2)^j B(A+j+½, B+½)`.

  They satisfy the same relations.

### 2.3 Lemma 2 (which classes are populated) [P + DV]
* By Künneth, the 7-form cohomology of the product is the tensor product of the factors'. A product form with an exact factor is exact:
  the other factors are top forms in their own variables, hence closed.
* Therefore all `s`-odd and `sf`-odd (`t`-odd) monomials integrate to zero on **every** closed cycle, not only on `D`.
* The one class that a non-real cycle could weigh differently from `D` is the `α`-odd (`sa`-odd) class with `s`, `sf` even. A one-arc
  cycle `c_U U + c_L L` pairs it with `c_U − c_L` instead of `c_U + c_L`.
* (B4) shows that this class is not populated:
  * by exact symmetries [D]: transverse reflection `(s, φ, α) → (−s, −φ, −α)`, under which the parity-even integrand is invariant, and the
    invariance of the raw point under `(φ, ψ) → (π−φ, π−ψ)`;
  * numerically [DV]: only round-off, of integrated weight ≤ 4·10⁻¹⁰ (§3.3).

  Hence **the populated cohomology is one-dimensional**.

### 2.4 Proposition RS [D given Lemmas 1–2 and §3]
Let Γ be any closed twisted 7-cycle of the local system of `Ψ_λ`:
* its ends are at radial infinity with `Re X → +∞`;
* it stays at fixed scaled distance from the divisors, or is regularised there (Pochhammer or Hankel loops);
* it may cover one or both α-arcs and one or both sheets `n₀`.

Then for every kind, state, reference and order `K`
```
∫_Γ Φ_κ = c_Γ(ν) · I_κ^{engine,(K)}(ν; Δ) + O(a^{K+1}),      c_Γ(ν) = ∫_Γ Ψ_λ / ∫_D Ψ_λ,
```
with `c_Γ` independent of κ, Δ, J and K.

*Proof.*
1. Insert the expansion of §2.1.
2. The remainder is `O(a^{K+1})·(1+|X|)^N` uniformly on Γ. This holds near the divisors because the `ε` are bounded there (B6), and on
   the radial ends by Watson's lemma (`e^{−X}`).
3. For each populated monomial, `∫_Γ Ψ_λ m = R_m ∫_Γ Ψ_λ = c_Γ R_m ∫_D Ψ_λ = c_Γ ∫_D Ψ_λ m` (Lemmas 1, 2).
4. `λ` is common to scalar, EE2 and EE4, and to state and reference, modulo the lattice (B1, B2); `R_m` absorbs integer shifts. ∎

**Corollaries.**
1. **H_RS to all orders** [D]. `ρ_Γ(ν) = ρ_engine(ν)` to all orders in `1/Δ`, for EE2 and EE4 at J = 0 and J = 2. In
   HRS_CHECK's notation this means `d(ν) = log ρ_cr − log ρ_eng` is **exactly ν-independent** order by order. The fitted κ there is
   only the Δ-only constants the engine drops (`β_{Δ,J}`, the overall `a`-power).
2. **ν < 1** [D]. The local system is non-resonant at every divisor for all real ν, including ν = 0 (B5). At ν = 0 the X-exponents of EE4
   are 9.25, 4.75, 1.25, 0.75 (plus integers), the `w` and `t` exponents are 0.75 mod 1, and the `z` exponents have real parts in
   ½ + ℤ. None is an integer. So the proposition and `c_Γ(ν) ≠ 0` hold uniformly down to ν = 0, and nothing distinguishes ν < 1 from the validated ν ≥ 10.
3. **Absolute normalisation.** Scalar and EE share `c_Γ` (B2). The scalar's agreement with the exact Mellin skeleton at all ν (state/ref,
   ≤ 5·10⁻⁴ `c₂` at Δ = 300) is therefore a direct measurement of the same cycle class that the EE images use.
4. **Single-sheet crystal codes are fine.** The graviton builders hard-wire `n₀ = +1` and HRS_CHECK's nodes use one sheet. Lemma 2 makes
   the one-arc class irrelevant.
5. **State and reference must be integrated on the same cycle.** HRS_CHECK's nodes satisfy this: the proposal is adapted on the state and
   re-used for the reference. The physical integral satisfies it trivially.

### 2.5 What the argument does not say
* It concerns the `1/Δ` asymptotic series at fixed ν and J. It says nothing about the uniformity in `J²/Δ`, so not about Result D's
  `(336, 28/32)` window, where the fixed-order expansion fails and the resummed EE engine does not exist.
* It does not say that the *exact* integral over the real section `D` equals the crystal integral at finite Δ. BELOW_FOLD §3.2's J > 0
  real-section QMC failure is a far-region (non-perturbative) effect of `D`, and is not in conflict with it.
* It does not test the crystal-node numerics. It makes them unnecessary for H_RS.

### 2.6 The remaining premise RD
**RD**: the crystal-cycle (equivalently, the physical) integral equals the integral of the exact integrand over a closed twisted cycle of
the Regge-scaled local system, up to terms smaller than any power of `1/Δ`.

* **Scalar: [DV].** The engine (= the `D` class by §2.4) equals the exact Mellin skeleton to ≤ 5·10⁻⁴ `c₂` at Δ = 300, ν ≤ 10, and to
  ≤ 0.001 at Δ = 1000 (`cache/bg_validate.txt`, BELOW_FOLD2 §2). By B2 this is the **same local system** as EE4.
* **EE at ν ≥ 10: [NE]** (HRS_CHECK §3). Direct crystal nodes reproduce the whole beyond-LO `Γ_exact` (up to −209 `c₂`) to ±0.5–1.8
  (Δ = 100, 300), at J = 0 and 2, for `k = 4` and `k = 2`.
* **EE at ν < 10: [D-sketch].** EE4/scalar = `m₁₂⁻² m₃₄⁻² v₀₄₂⁻⁴ pref ((v′−u′)² − 1)/4`, a rational function (code reading of
  `bg_ee.ee_pieces` vs `bg_geom.scalar_pieces`). Its singular loci are:
  * the atom divisors, where it is uniform by B6;
  * `m₁₂ = m₃₄ = 0` and `m₁₄ = m₂₃ = 0` (also singular for the scalar);
  * `m₁₅, m₂₅, m₃₅, m₄₅ = 0`.

  All of them except the atom divisors lie at `|raw x| ~ 1`: scaled radial distance `~Δ`, or `a²C ≈ ±1`. There:
  * `B^{Δ/2}` gives `e^{−O(Δ)}` on `D` (`m₁₅, m₂₅ > 1` on `D`);
  * descent cycles of `λF ∋ (Δ/2) log B` cannot approach the pole of `B` at `m₁₅m₂₅ = 0`.

  A polynomially bounded rational factor cannot lift an `e^{−O(Δ)}` region to the `10⁻⁴` (300) / `10⁻⁵` (1000) relative level that a
  δ = 0.2 failure needs.
* **Moreover a failure would have to be Δ_g-independent in `c₂`.** The binding `s_d` is +0.5096 (300) and +0.5077 (1000); a
  non-perturbative term cannot do that.

## 3. Bookkeeping [DV] (`rs2_bookkeep.py`, `rs2_uniform.py`, `rs2_parity.py`)

### 3.1 Leading exponents (B1, B2, B5)
| kind (J = 0; all Δ) | `E₀` (`p₁ p₂ q₁ q₂ w cf u v`) | `E₁` (same atoms) |
|---|---|---|
| scalar `k = 2` | 2.5, 1.5, 2.5, 1.5, 0.75, 2.5, −0.25, −0.25 | −1, +1, +1, −1, −½, +1, −½, −½ |
| EE2, EE4 | 18.5, 9.5, 2.5, 1.5, 8.75, 2.5, −0.25, −0.25 | identical |
| EE4, J = 2 (Δ = 314.935 and 2·10⁴) | identical to J = 0 | identical |

* `E₀(EE) − E₀(scalar) = (p₁ +16, p₂ +8, w +8)`: even in `p` (integer in `X`), integer in `w`. The local systems coincide.
* `E₁ ≠ 0` at every divisor; at ν = 0 every exponent is non-integer mod the lattice (B5).
* State, reference and Δ = 100, 300, 1000, 2·10⁴ have the same `E₀, E₁` (the Δ-dependence is all in `corr`).

### 3.2 Monomial content (B3, B4), Δ = 300, J = 0, K = 4 (other states: same pattern; `cache/rs2_bookkeep_*.txt`)
| kind | monomials | negative powers of `s, sf, ca, sa` | parity violations (`p₁≡p₂≡s`, `q₁≡q₂`, `cf+sa`) | `s+sf+sa` odd / (`sa` odd, `s, sf` even) | surviving: lattice violations |
|---|---|---|---|---|---|
| scalar | 1 114 | 0 | 0 | 0 / 0 | 0 |
| EE2 | 161 078 | 0 | 0 | 50 058 / 13 246, all `|c| ≤ 4.4·10⁻¹³` of max | 0 |
| EE4 | 294 043 | 0 | 0 | 69 992 / 18 278, all `|c| ≤ 4.7·10⁻¹³` | 0 |
| EE4, J = 2, Δ = 314.935 | 287 781 | 0 | 0 | 69 156 / 18 092, `≤ 2.2·10⁻¹³` | 0 |

* After removing round-off residues (`|c| < 10⁻¹⁰` of the block maximum), every surviving EE4 monomial is absolutely integrable on `D`:
  * real exponents `X₁^{≥2.25}, X₂^{≥1.75}, Y₁^{≥0.25}, Y₂^{≥−0.25}`;
  * `w^{≥1.75}`, `cf^{≥0.5}`.
* The `u^{−1}, v^{−1}` shifts are representation redundancies (`ca = 1 − u = v − 1`). The relevant combinations stay bounded as `α → 0, π`
  (§3.4).
* So at J = 0 the engine's closed forms are ordinary convergent integrals, and `reg D = D`.

### 3.3 The `α`-odd class, integrated (B4′, `rs2_parity.py`)
| order `k` (j = 0, 1) | 0 | 1 | 2 | 3 | 4 | 4 (j = 1) |
|---|---|---|---|---|---|---|
| `|class on one α-arc| / |image|`, max over ν = 0.02 … 8.02 | 8·10⁻¹² | 1.3·10⁻¹¹ | 2.6·10⁻¹¹ | 3.3·10⁻¹¹ | 7.7·10⁻¹¹ | 3.8·10⁻¹⁰ |

That is ≤ 4·10⁻⁵ `c₂` (relative 4·10⁻¹⁰ × `m²`) even for a cycle with `c_L = 0`.

### 3.4 Uniformity at the divisors (B6, `rs2_uniform.py`, `cache/rs2_uniform_K4.txt`)
* Each of the 40 logged or inverted factors was tested as `ε_k` for `k ≤ 4`, at each of the 8 divisors with `t = 10⁻²…10⁻⁶`.
* The growth exponent of `max|ε_k|` is ≥ **−0.005** (fit noise) everywhere: no boundary layer.
* The `u^{−1}, v^{−1}` monomials of `A`'s `ε₄` sum to a bounded function: `|ε₄| = 44.1` at `α = 0.003` vs 63 at 0.1.
* All `ε` are finite polynomials in `a` of degree ≤ 4 (BELOW_FOLD2 §1.2), so the uniformity holds to all orders.

### 3.5 Lemma-B `R` and `v′ − u′` exactly
* `v′ − u′ = C(1+P)/(4a²q_xq_y)` and `R = 4a⁴[C²(1+P)² − 16a⁴q_x²q_y²]/((1+P)² − a⁴C²)²`, with `P = a⁴q_xq_y`.
* The engine's `R.c[0]/(4C²) = 1` to 10⁻¹⁵; the series matches the exact `R` to 10⁻¹⁵ at depth 8.
* **Minor finding.** At the production depth (K + 1 = 5), `v′ − u′` and `R` carry their relative corrections only to `a³`
  (checked: depth 5 gives `K = 3` for both, depth 6 gives `K = 4` including the `a⁴` term). Their `a⁴`
  pieces (`P`, `2a⁴C²`, …) are missing from the production `a⁴` term.
  * Effect: ≤ 0.023 / 0.045 / 0.037 `c₂` per harmonic on ν ∈ [0,2) / [2,4) / [4,8); −0.0002 in `s_d`.
  * Harmless, and corrected in the `K = 5` images (depth 6 is complete through `a⁵`; a depth-9 `R` and `v′ − u′` give identical
    results at Δ = 300 and 1000).

## 4. Truncation: the `a⁵` term (`rs2_engine5.py`, `rs2_trunc.py`) [NE]

* The binding state `++` J = 0 at the gap has new `K = 5` (depth 6) EE4 images, Δ = 300 and 1000. The reference is the production one
  (Δ_ref = 2·10⁴, whose `a^k` terms are `~(5·10⁻⁵)^k`).
* Per harmonic, max over the window; Δ_g = 300 (1000 in brackets):

| ν window | `max|Γ_K|` | `a⁴` term | `a⁵` term | `a⁵` term in `s_d` (H0) |
|---|---|---|---|---|
| [0,1) | 4.26 | 0.127 (0.013) | 0.0058 (0.0002) | +0.00016 (+0.00001) |
| [1,2) | 10.38 | 0.203 (0.020) | 0.0096 (0.0003) | +0.00060 (+0.00002) |
| [2,3) | 10.53 | 0.202 (0.020) | 0.0096 (0.0003) | −0.00004 |
| [3,4) | 7.61 | 0.110 (0.011) | 0.0072 (0.0002) | −0.00010 |
| [4,6) | 2.15 | 0.029 (0.003) | 0.0027 (0.0001) | −0.00001 |
| [6,8) | 1.45 | 0.023 (0.002) | 0.0011 | 0 |
| [8,10) | 0.59 | 0.028 (0.003) | 0.0019 | 0 |
| [10,18) | 2.05 | 0.12 (0.012) | 0.019 (0.0006) | +0.00003 |
| [18,41) | 4.47 | 9.5 (0.92) | 0.83 (0.024) | 0 (0.1 % of the weight) |

| gauge | `s_d` (300): K = 3 / 4 / 5 | `a⁵` term | `s_d` (1000): K = 3 / 4 / 5 | `a⁵` term |
|---|---|---|---|---|
| H0 | +0.4937 / +0.5090 / **+0.5096** | +0.00064 | +0.5062 / +0.5077 / **+0.5077** | +0.00002 |
| H2 | +0.4943 / +0.5092 / +0.5098 | +0.00060 | +0.5068 / +0.5083 / +0.5083 | +0.00002 |
| H4 | +0.4926 / +0.5075 / +0.5081 | +0.00060 | +0.5051 / +0.5066 / +0.5066 | +0.00002 |
| SB | +0.5090 / +0.5218 / +0.5222 | +0.00037 | +0.5215 / +0.5228 / +0.5228 | +0.00001 |
| XP | +0.4993 / +0.5137 / +0.5142 | +0.00054 | +0.5117 / +0.5132 / +0.5132 | +0.00002 |

* The series converges geometrically. The `a⁴`/`a⁵` ratio in `s_d` is 24 (300) and 80 (1000), consistent with one power of `1/Δ` per order.
* BELOW_FOLD2's error convention (size of the last term, ±0.0155 / ±0.0016) is conservative by ×24 / ×80.

## 5. Tables per ν window (binding `++` J = 0 at the gap; Δ_g = 300, 1000 in brackets)

| ν | `s_d` weight (H0) | tolerance per harmonic δ = 0.2 / 0.05 | cycle question | truncation (`a⁴` term, bound) | direct crystal check (HRS_CHECK, 1σ ×3) | status |
|---|---|---|---|---|---|---|
| [0, 1) | 24.5 % | 30.2 / 1.67 | Prop. RS [D] (non-resonant at ν = 0) | 0.13 (0.013) | not reachable | **[D given RD]** |
| [1, 2) | 43.9 % | 16.9 / 0.94 | [D] | 0.20 (0.020) | not reachable (CRN lost) | **[D given RD]** |
| [2, 3) | 9.2 % | 80.2 / 4.4 | [D] | 0.20 (0.020) | ±1800 (300) | **[D given RD]**; direct consistent, uninformative |
| [3, 4) | 11.4 % | 65.2 / 3.6 | [D] | 0.11 (0.011) | ±390 (100), ±520 (300) | as above |
| [4, 6) | 3.7 % | 200 / 11.1 | [D] | 0.03 (0.003) | ±130 (J = 2), ±210 (300) | as above |
| [6, 8) | 2.8 % | 265 / 14.7 | [D] | 0.02 | ±39 (300) | [D given RD] + [NE, weak] |
| [8, 10) | 0.8 % | 979 / 54 | [D] | 0.03 | ±61 (300) | [D given RD] + [NE, weak] |
| [10, 18) | 3.7 % | 200 / 11.1 | [D] | 0.12 (0.012) | ±2–6 | [D] + [NE] (RD tested) |
| [18, 41) | 0.1 % | 8948 / 496 | [D] | 9.5 (0.9) | ±2–5 | [D] + [NE] |

* The tolerances (uniform per-harmonic `|dG|` that would use up `s_fail − s_d`) are recomputed from Stage 7's H0 weights. They agree with
  HRS_CHECK §4 (17 on [1,2], 65–80 on [2,4]).
* The worst ratio of truncation bound to tolerance is 0.20/16.9 = 1.2 % (δ = 0.2, ν ∈ [1,2], 300).
* For δ = 0.05 the same ratio is 22 %, i.e. **δ = 0.05 is also not threatened by truncation** at the binding state. It stays limited by
  `σ_hU` and the `+−` J = 5 state (RU issue 8), not by H_RS.
* `++` J = 2 at 1.05Λ (Δ_g = 300): same bookkeeping (§3.1–3.2), same conclusion. Truncation (BELOW_FOLD2's `a⁴` term) ±0.010 in `s_d`
  (300), ±0.001 (1000). Its δ = 0.2 FF/N is +0.181, a larger margin than at J = 0.

## 6. Propagation: `Δs_d` and `ΔFF/N` (binding `++`/OO J = 0 at the gap)

| source | Δ_g = 300 | Δ_g = 1000 | grade |
|---|---|---|---|
| H_RS (real section vs crystal), all orders | `Δs_d = 0` | `0` | [D given RD] |
| truncation: measured `a⁵` term (H0) | +0.00064 | +0.00002 | [NE] |
| truncation: honest bound (size of the `a⁴` term, ×24 / ×80 the `a⁵`) | ±0.015 | ±0.0015 | [NE] |
| production `a⁴` completeness (`v′−u′`, `R` at depth 5) | −0.0002 | −0.00002 | [DV] |
| RD (non-perturbative) | `e^{−O(Δ)}`; not quantified | same | [D-sketch] |
| **total** `|Δs_d|` | **≤ 0.015** | **≤ 0.0015** | |
| `|ΔFF/N|`, δ = 0.2 (slope −0.1126) | ≤ 0.0017 | ≤ 0.0002 | |
| `|ΔFF/N|`, δ = 0.05 (slope −0.0863) | ≤ 0.0013 | ≤ 0.0001 | |

* With the `K = 5` images, FF/N at the binding state is:

  | δ | Δ_g = 300 (H0 / SB) | Δ_g = 1000 (H0 / SB) |
  |---|---|---|
  | 0.2 | **+0.1427 / +0.1412** | +0.1428 / +0.1411 |
  | 0.05 | +0.0061 / +0.0050 | +0.0062 / +0.0049 |

  These are unchanged from BELOW_FOLD2/3 at the quoted precision.
* HRS_CHECK's direct comparisons (`Δs_d = −0.10 ± 0.35` on [5, 25] at 300) are consistent with 0, as the proposition requires.
* The margin `s_fail − s_d = 1.27` (δ = 0.2) is ≥ 85× the total bound.

## 7. Regrading `H_RS`

| part of `H_RS` | before (HRS_CHECK §6) | now | evidence |
|---|---|---|---|
| scalar, all J | [DV] | [DV]; also an instance of Prop. RS | BELOW_FOLD2 |
| EE at LO | [DV] | [DV] (= order 0 of Prop. RS) | — |
| **EE `k = 2, 4`, J = 0, 2: engine = crystal 1/Δ series, all orders, all real ν** | [NE] ν ≥ 10; [NE, weak] 6–10; [I] 2–6; [I] untested < 2 | **[D]** (Lemmas [P]; hypotheses [DV] §3) | §2, §3 |
| RD for EE (crystal integral Regge-dominated) | (implicit) | [NE] ν ≥ 10; **[D-sketch]** ν < 10; [DV] on the same local system via the scalar | §2.6 |
| **`H_RS` as used by Result Q's binding state `++` J = 0 (300, 1000), every ν window incl. ν < 1** | **[I]** | **[D given RD]**, net **[D-sketch]**; error budget = truncation, `|Δs_d| ≤ 0.015 / 0.0015` | §4–6 |
| `++` J = 2 (1.05Λ) | [NE] above fold / [I] below | same as J = 0 | §3 |
| OO J = 0 | via `H_OO` (reduces to EE `++` + H_RS) | via `Γ_OO = Γ_++` ([P given lemma] order 0, [NE] beyond, `|Δs_d| ≤ 5·10⁻⁶`) + the above. Directly by Prop. RS: [D-sketch] (OO pieces not booked) | PARITY_ODD_CLOSE |
| `+−` J = 4 (2×2 entries), J ≥ 4 | — | [D-sketch]: same construction, bookkeeping not run | — |
| Result D `(336, 28/32)` window | fails entrywise 3–16σ | **not covered**: the issue there is `J²/Δ` non-uniformity and the missing resummed EE engine, not the cycle. Prop. RS implies a fixed-J cycle discrepancy cannot be the cause [D-sketch] | RU issue 2 |

**Recommended wording.**
> `H_RS` (EE beyond LO, real section ≡ crystal cycle) holds to all orders in `1/Δ` at every real ν for the `k = 2, 4` images at J = 0, 2.
> * The Regge-scaled integrand is a product of one-variable twisted local systems times polynomials, so integration by parts makes every
>   order cycle-independent up to one common constant [D]. The hypotheses (same local system as the scalar, no new divisors, no boundary
>   layers, one populated cohomology class) are checked on the engine's pieces [DV].
> * What remains is the premise that the crystal-cycle integral is Regge-dominated ([DV] for the scalar on the same local system, [NE]
>   for EE at ν ≥ 10, [D-sketch] below).
> * The binding `s_d` is then +0.510 ± 0.015 (300) / +0.508 ± 0.002 (1000), the error being the engine's truncation.

For UNIFIED_THEOREM §5, H_RS moves from rank 1 ([I]) to the [D-sketch] group. Its residual is a non-perturbative premise shared with the
already-validated scalar.

## 8. What remains

1. **RD as a proof.** An estimate of the far-region (`|raw x| ~ 1`) weight on the crystal descent cycle for EE4, with the loci
   `m₁₂ = 0` (EE exponent −10.5 vs scalar −6.5) and `m₁₄ = 0` (double pole, as in the scalar), would make the closure [D]. A cheap
   empirical handle: compare the crystal/engine `d(ν)` of HRS_CHECK at two Δ for a non-perturbative drift. No such drift is seen at
   ν ≥ 5.
2. **Independent low-ν numerics.** Still beyond reach (HRS_CHECK §7). Their role is now only a test of RD, not of the reduction.
3. **Bookkeeping not run:**
   * OO (parity-odd `⟨TTO⟩` squared) and `+−` J = 4 (2×2 entries, 75 min per state);
   * the `k = 3` image `Ĉ₃` (`j′ = 1`) used in Q's smeared part;
   * the J/Δ-resummed contraction needed for Result D's spins.

   The argument applies verbatim if B1–B6 hold there.
4. **Production engine hygiene.** Use depth ≥ K + 2 (or deep `v′ − u′`, `R`), so that the last order is complete (§3.5). The effect on
   every quoted number is ≤ 0.0002 in `s_d`.
5. **The HRS_CHECK heavy tails.** The EE integrand is absolutely integrable on `D` at every order (§3.2). The heavy tails of the `k = 4`
   crystal nodes therefore do not come from `R`'s poles, which are in the far region; the likely source is the radial polynomial growth
   (`R ≈ 4a⁴C²`, EE density `∝ X₁^{9.25}`) relative to the proposal [I, not tested]. This is of no consequence for H_RS now.

## 9. Reproduce (`ads4/code/spinfunc/`; every job `ulimit -v 3000000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`)

```
# bookkeeping B1-B5 (90 s per J = 0 state, 11 min for J = 2):
python3 rs2_bookkeep.py 300 0 4          # -> cache/rs2_bookkeep_D300_J0_K4.txt (+ .pkl)
./rs2_queueA.sh                          # Delta = 1000, 2e4, 100 (J = 0)
./rs2_queueB.sh                          # J = 2: Delta = 314.935 (K = 4), 2e4 (K = 2)
# B6 uniformity at the divisors (95 s):
python3 rs2_uniform.py 4                 # -> cache/rs2_uniform_K4.txt
# B4' alpha-odd class on one arc (2 min):
python3 rs2_parity.py 300 0              # -> cache/rs2_parity_D300_J0.txt
# K = 5 images (4.5-5 min each) and truncation tables:
./rs2_queueC.sh                          # rs2_engine5.py 300 0; 1000 0  -> cache/rs2_eeq_D{300,1000}_J0_K5.pkl
./rs2_queueD.sh                          # parity + 'full' variants (deep v'-u', R) -> cache/rs2_eeq_D*_J0_K5full.pkl
python3 rs2_trunc.py 300; python3 rs2_trunc.py 1000; python3 rs2_trunc.py 300 full; python3 rs2_trunc.py 1000 full
                                         # -> cache/rs2_trunc_D{300,1000}[full].txt (+ .pkl)
```
Logs: `cache/rs2_queue*.log`. The first Δ = 300 bookkeeping run (before the round-off-filtered convergence line was added) is kept as
`cache/rs2_bookkeep_D300_J0_K4_v1.*`.

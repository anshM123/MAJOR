# ANGULAR_ROUTE: celestial recurrence, continuous dual Hahn structure, and an all-harmonic positivity mechanism

*25 Nov 2026. Author-supplied research (Ansh, via an independent read-only review of the repository at 1ba3e88). Recorded and partly
verified by the coordinator. Verification code: `code/angular/verify_cdh.py` (exact sympy, runs in about 1 min).*

## 0. Status at a glance

| item | status |
|---|---|
| SP1: three-term recurrence for the celestial chord-Fourier products k_{ℓ,d,m}(ν) | proof supplied (ODE for H_ℓ conjugated by the chord factor); numerically checked by the author (residual 1e-12 for d=1, 5e-6 for d=0) |
| d=1, ℓ=0: k₃/k₂ = (15−4ν²)/64 and k₄/k₂ = (ν⁴ − 35ν²/2 + 569/16)/400 | **[P] verified exactly**; this proves the K_B2A4 k=3 polynomial `1 − (35/2)ν⁻² + (569/16)ν⁻⁴` from the recurrence |
| k₅/k₂ = (39159 − 26956ν² + 2896ν⁴ − 64ν⁶)/921600 | **verified exactly** |
| Recurrence solutions = normalised continuous dual Hahn ₃F₂(−j, a+iν, a−iν; a+b, a+c; 1), with (a,b,c) = (2d+3/2, ½, ℓ+½) and m = j+d+1 | **verified exactly** for j ≤ 8 and (d,ℓ) ∈ {(1,0),(0,0),(1,4),(0,4),(1,2),(2,3)}; the general statement follows from the standard CDH recurrence |
| Cubic light-side polynomial Q₆: P(ν) = R₋Q₆(ν+i) + R₊Q₆(ν−i), with Q₆ = ν⁶ − 245ν⁴/4 + 9619ν²/16 − 51575/64 | **verified exactly** (sympy residual 0); K(L) moment formula reproduces the certified K = 0.99996752314 at L = Λ |
| SP2: ∫e^{−tν²}Q_j dμ > 0 for all j (more generally, for any completely monotone multiplier) | proof checked: interpolation at the j roots of Q_j, sign (−1)^j of the leading coefficient, orthogonality. This is classical-type and correct |
| Uniform angular tail: \|∫_R^∞ e^{−tν²}Q_j dμ\| ≤ e^{−tR²}‖Q_j‖, with ‖Q_j‖² = j!(ℓ+1)_j/((2d+2)_j(2d+ℓ+2)_j) | Cauchy–Schwarz, correct given the norm formula (norm not independently checked) |
| Gamma factorisation of the LARGESPIN_L Riesz measure: A ~ Γ(λ,2), V₀ ~ Γ(λ−½,2), Z ~ N(0,1), M² = AV | verified by Monte Carlo against direct quadrature (0.1%) |
| Boost-weighted rare-event bounds (M<800: 1e-5 including e^{4β}(R/M)⁸, times (1+J/489)^{−(n+d−k)}) | author's exact arithmetic; not re-run |
| Selection rule: the energy-detector transition between distinct conformal families has no harmonics 0, ±1 (momentum Ward identity ∫n^μE = P^μ) | argument correct for momentum eigenstates |
| **Full physical k=2,3 polarisation-matrix positivity** | **open** |

## 1. What remains (the angular route)

1. **Heavy-family transition matrices.** Write each heavy state's contribution to the k=2,3 sum rules as Σ (amplitude bilinears) × harmonic kernels k_{ℓ,d,m}. This needs the first retarded-time moment of E between different families. Translation fixes its ∂P part, but an extra transition matrix remains. C1_FINAL evaluates it only inside the stress-tensor family.
2. **Normalisations and signs.** Fix the channel normalisations k_{ℓ,d,d+1}(ν) (and their signs) that multiply Q_j, so that SP2 applies to the physical combination. It is also needed for the matrix (polarisation) structure, not only the diagonal entries.
3. **Cubic sensitivity must survive.** A functional that is positive harmonic by harmonic on heavy states must still be non-positive on the light stress-tensor contribution (harmonics up to |λ−h| = 4) in the direction of t₄ and α_T. Otherwise the bound is empty. Check this first.
4. If 1–3 hold: Gaussian (or completely monotone) low-frequency smearing gives heavy positivity for every spin and harmonic, with the Gamma tail bounds for the low-mass plane-wave components. That would be unconditional positivity, and the bound would come from the compressed-weight light side (≈ 2.9e-3 at L = 0.05Λ; LOWNU_SCOPE).

## 2. First-moment reduction (author-supplied, 25 Nov)

**(a) Conformal identity.**
* *Claim.* The special conformal generator 𝖪_n along the detector direction gives [𝖪_n, E(n)] = −6E⁽¹⁾(n). This needs the detector density to transform as (1+2εu)^{−5}𝒯(u/(1+2εu)) and the integration-by-parts boundary term to vanish.
* *Consequence.* Between distinct families, E⁽¹⁾_{OT} = −(1/6)(𝖪_{O,n}E_{OT} − E_{OT}𝖪_{T,n}). So the first-moment transition is fixed by the zeroth-moment amplitudes and the conformal representations of O and T; it is not independent data.
* *Coordinator note.* The structure is right: with a weight-h density, δE = −(4h−4)εE⁽¹⁾. The coefficient −6 corresponds to exponent −5 (h = 5/2). That exponent is **not yet checked** against the detector-limit conventions; check it with (b) on the stress-tensor family, where C1_FINAL already has the first moment.

**(b) Monomial rule.**
* *Statement.* With the retarded time centred at the midpoint of x₃ and x₄, I₁^mid(a,b) = (a−b)/(2(a+b−2))·(2z·x)·I₀(a,b), for the (F2) integral I₀ = ∫dα X₁₄^{−a}X₁₃^{−b} of KERNELS_D3.md.
* *Verification* [P]: α − α_mid = −½(X₁₄ + X₁₃), so I₁ = −½[I₀(a−1,b) + I₀(a,b−1)]. The Γ-ratios of (F2) give I₀(a−1,b) = −I₀·(a−1)δ/(a+b−2) and I₀(a,b−1) = I₀·(b−1)δ/(a+b−2), with δ = 2z·x. Combining gives the claimed result.
  * The author also checked 2860 integer-exponent cases and higher moments.
* *Pipeline change* (`code/kernels/pipeline.py`): multiply each expanded monomial by the corresponding factor, with ℓ → ℓ+1 for the extra power of (z·x). The Fourier-transform machinery is reused unchanged.
* *Caveat.* The physical moment is centred at the origin, not the midpoint. The difference is (α_mid)×(zeroth moment), which in momentum space is a momentum derivative of the zeroth-moment amplitude (the translation piece).

## 3. Tridiagonal first-moment operator, explicit weight, telescoping (author-supplied, 25 Nov)

Author claims, with the coordinator's checks (`code/angular/verify_sp3.py`):

| claim | check |
|---|---|
| First-moment matrix M = (i/6)·A_{Δ,J}H, with A tridiagonal: diagonal T_Δ + C_J − 2m², off-diagonal −b_m, where T_Δ = (Δ−2)(Δ−1), C_J = J(J+1). Author checked 46/46 components at J = 2, 4 | not re-checked |
| **SP3**: v†Av = Σ\|√x_m v_m − √y_m v_{m+1}\|² + (T_Δ − C_J)\|v_{−J}\|², hence A ⪰ 0 for twist ≥ 2 | **identity verified algebraically and numerically** (300 random Δ, J) |
| A ⪰ (√T_Δ − √C_J)² I | **verified numerically** (300 random Δ, J) |
| Finite-difference identities for the angular kernels Z_m and Y_m in terms of k_m | not re-checked |
| Explicit weight χ = −ν²(4ν²+9)(4ν²+49)/(32(4ν²+1)cosh πν), with CDH density w = 8ν²/cosh πν (∫w = 1) | ∫w = 1 verified |
| f_m = (½)_{m−1}(1)_{m−1}/(2)_{m−1}² | **verified** against quadrature, m ≤ 6 |
| z_m/f_m = 3(12m²−7m−24)/(2(m+1)(2m−3)); z_m > 0; y_m < 0 | **true up to one overall sign** for 2 ≤ m < 3000: with Z_m as written and f_m > 0, z_m < 0 and y_m > 0. Fixing this global sign against the physical functional is essential. It fails at m = 1. |
| z_{m+1}/z_m > (m−1)/(m+2), m ≥ 2 | **verified exactly**, 2 ≤ m < 3000 (sign-independent) |
| Ladder ratios a_{k+1}/a_k = ((k+4)/(k+1))√(X_k/Y_k) on the b = 0 and b = 4 ladders | repository-verified through J = 12; **open for general J** |
| Telescoping identity, giving a weighted heavy first-moment form ≥ [T_Δ − b(b−1)]z₂\|a₀\|² | not re-checked; depends on the ladder ratios and the global sign |
| Prospective bound √(t₄²+α_T²) ≤ 0.0224 at Δ_g = 1000, from (139/4 − 48 log 2)/(47000 + 303/5120). Note y₄ = 303/5120 in this convention | **not established**; needs physical sum-rule matching, light-side, contact and double-trace terms, and convergence |

The global claim is that the exact rest-frame transition-matrix representation would remove the need for any beyond-Regge bound. That holds only if the physical functional's heavy action equals this rest-frame quadratic form, which is the physical-matching step and is open.

## 4. All-spin helicity normalisation (author-supplied, 25 Nov)

* *Claim.* Gillioz, "Spinors and conformal correlators", eqs. (2.37)–(2.39), gives N_{m+1}/N_m = (J−m)(Δ+m−1)/((J+m+1)(Δ−m−2)). Hence N_m/N_J = C(2J,J+m)·Γ(Δ+m−1)Γ(Δ−m−1)/(Γ(Δ+J−1)Γ(Δ−J−1)).
* *Coordinator check.* The closed form solves the recurrence. It equals the m-dependence of KERNELS_D3 §3.3, C(2J,J+m)(Δ−J−1)_{J−|m|}(Δ−1)_{|m|}, using (Δ−J−1)_{J−|m|}(Δ−1)_{|m|} = Γ(Δ−|m|−1)Γ(Δ+|m|−1)/(Γ(Δ−J−1)Γ(Δ−1)). Identical [P].
* *Effect.* KERNELS_D3 §3.3 (previously "derived for J ≤ 5, used as a formula above") becomes all-J, **once the citation is checked against the paper**: the equation numbers and conventions must match before it goes in the supplement.
* *Remaining.* The ladder problem is now only about the unnormalised energy-detector transition amplitudes.

## 5. Exact smeared coefficients, SP4 (ladder-free) at ν = i/2, and the sign calibration (author-supplied, 25 Nov)

Coordinator checks are in `code/angular/verify_sp4.py`.

| claim | check |
|---|---|
| f_m = (½)_{m−1}/((m−1)!m²) proved analytically (moments of w, plus an arcsine × uniform × uniform product density with h(x) + h(1−x) = 4ρ(x) − 2) | equals §3's formula; values verified by quadrature (m ≤ 6). The proof steps are not re-checked |
| z_{m+1}/z_m − (m−1)/(m+2) = (12m³−17m²+43m−48)/(2(m+1)(m+2)(12m²−7m−24)) > 0 for all m ≥ 2; the numerator at m = n+2 is 12n³+55n²+119n+66 | **verified** (exact, m < 60; the polynomial shift is symbolic). This proves the telescoping weight inequality for all m |
| Sign convention: z_m^raw < 0 and y_m^raw > 0. The positive z used in telescoping is −Z_m | consistent with my §3 finding |
| At ν = i/2: f_m = 1/m and z_m = 12/m (m ≥ 2; z₁ = 13) | **verified** |
| **SP4**: with Z = diag(12/(k+2)) on the helicity sectors b = 0, 4, B = (ZA+AZ)/2 ⪰ (√T−√C)² Z for ANY amplitude vector (no ladder needed). Proof by a positive test vector with nonpositive off-diagonals | **verified numerically**: generalised λ_min(B,Z) ≥ (√T−√C)² on 400 random (b, J ≤ 300, Δ) |
| Real-frequency weights do NOT give arbitrary-amplitude positivity (counterexample at J = 10⁶, Δ = J+1000) | **reproduced**: λ_min(B,Z)/C = −3.7e-3 at J = 10⁶ on the first 25 helicities (+5.9e-3 at J = 10⁴). With real frequencies the physical ladder identities are needed |
| ν = i/2 caveat: the celestial dimension δ = ½+iν = 0, so coincident-direction contact terms are not suppressed. A contact-subtracted, admissible functional is needed, with a finite Einstein light value | **open** |
| R₋+R₊ = −8(4ν²+1)(16ν⁴−728ν²+969)/((4ν²+9)(4ν²+49)) | **verified** (sympy) |
| For χ = −ν²(4ν²+9)(4ν²+49)/(32(4ν²+1)cosh πν): ∫χ(R₋+R₊) = 15/4; Einstein ∫χ·[R₋D(ν+i) + R₊D(ν−i)], D = 1/(ν²+9/4), equals 139/4 − 48 log 2; ∫χP = −13635/8 | **all three verified** to 25 digits (kc_grid V3 convention for the Einstein term) |
| −Ω_light^red = 139/4 − 48 log 2 − (303/5120)\|c\|², \|c\|² = (t₄²+α_T²)/16 | the K_B2A4 normalisation of the \|c\|² coefficient is not re-checked. **The signs are right: cubic coupling reduces the positive budget, and the forward heavy normalisation is positive** |
| Scalar (J = 0) first moment from the Wightman pipeline: M_mid/H = i(Δ²−7)/6; normalisation shift −i(Δ−3)/2 gives M/H = i(Δ−1)(Δ−2)/6 in both parities. This reproduces the J = 0 case of the tridiagonal A. The physical combination is −T_Δ Z₂/6 | not re-run. The algebra (Δ²−7)/6 − (Δ−3)/2·... is consistent: (Δ²−7 − 3Δ + 9)/6 = (Δ−1)(Δ−2)/6 ✓ |
| −2048/((4x+9)(4x+49))·Z₂/k₀ = (R₋+R₊)/6 with Z₂/k₀ = (4x+1)(16x²−728x+969)/1536 | consistent with the verified R₋+R₊ (the identity reduces to it) |

**Where this leaves the route:**
* The sign chain is calibrated in one convention at J = 0: heavy first-moment term, forward Regge kernel and Einstein light term.
* Open:
  * the spinning transition operator for all J (physical matching), or the ladder identities;
  * the ν = i/2 contact subtraction, if the ladder-free SP4 route is used;
  * the K_B2A4 normalisation of the 303/5120 coefficient;
  * convergence, plus double-trace and contact terms for the new weight.

## 6. Generating function, contact regularity, angular light side, and the all-spin conditional bound (author-supplied, 25 Nov)

Coordinator checks are in `code/angular/verify_s6.py`.

| claim | check |
|---|---|
| S(z) = Σf_m z^m = 4[1 − √(1−z) + log((1+√(1−z))/2)], with f₀ = −8(1 − log 2) | **series verified** (m < 12) |
| 𝓕(φ) = 2Re[S(e^{iφ}) − S(1)] has 𝓕(0) = 𝓕′(0) = 0 and 𝓕 ≈ −(4√2/3)\|φ\|^{3/2} | **verified numerically** (ratio → 1 as φ → 0). The reduced Z and Y kernels sum to 0 ⇒ they annihilate δ(φ) and δ′(φ) contact terms (argument not re-checked) |
| Y₀ = 6f₀ + 13f₁ + 2f₂ = −139/4 + 48 log 2 = −(Einstein integral of §5) | **verified exactly**. The angular recurrence reproduces the rational-kernel Einstein term |
| Y₄ = 303/5120, the \|c\|² light coefficient of §5 | **verified exactly** |
| y_m/f_m = 3(48m⁵−176m⁴+303m³−376m²+297m−480)/(16(m−1)(m+1)²(2m−5)(2m−3)), y_m > 0 | **verified** for 3 ≤ m < 80. The numerator at m = n+3 is 48n⁵+544n⁴+2511n³+5807n²+6654n+2616 (all positive), so y_m > 0 for all m ≥ 3; y₂ = 41/64 |
| Heavy form 𝓗 = Σy_{k+2}\|a_k\|² − (1/6)Σz_{k+2}Re[ā_k(Aa)_k] (raw z < 0), so both terms are ≥ 0 given the ladder and the rest-frame reduction | structure consistent with the verified signs; **conditional** |
| b = 4 ladder: 𝓗 ≥ (5/48)[(Δ−1)(Δ−2) − 12]\|a₀\|²; fourth-harmonic W₂⁴ = 100(J−5)(J−4)(J+5)(J+6)/((Δ−7)(Δ−6)(Δ+3)(Δ+4)). Uniform bound > 47266 for twist ≥ 1000, J ≥ 6, via min_x (x+998)⁶/x⁴ = (729/16)·998² | **min and constant verified** (47266.56, argmin x = 1996). The inequality chain from W₂⁴ to y⁶/x⁴ is not re-checked |
| Conditional bound: 47266\|c\|² ≤ 139/4 − 48 log 2 − (303/5120)\|c\|² ⇒ √(t₄²+α_T²) < 0.022375 at twist gap 1000 | **arithmetic verified** (0.0223748). **Conditional** on the all-spin ladder and first-moment identities and on physical matching |
| J = 0: 𝓗 = [41/64 + (5/48)(Δ−1)(Δ−2)]\|a₀\|² ≥ (59/96)\|a₀\|², so light scalars have the safe sign | arithmetic consistent; relies on the §5 scalar first moment (pipeline reproduction pending) |
| ν = i/2 shortcut fails: the finite-part light budget is −14 − (31/15)\|c\|² < 0. Keep the real-frequency route | not re-checked; recommendation adopted |

**Open for an unconditional statement:**
1. the all-J transition ladder (the b = 0 and b = 4 amplitude ratios);
2. the all-J first-moment identity (the tridiagonal A);
3. physical matching of the rest-frame heavy form to the sum rule, with convergence, contact and double-trace terms;
4. the pipeline reproduction of the scalar first moment.

## 7. All-J first-moment operator from the special-conformal Ward identity (author-supplied, 25 Nov)

**Derivation (author):**
* E(n) = 2∫T(Z − αX_∞, −X_∞)dα with α = 2u. The special conformal transformation along n acts by X_∞ → X_∞ − εZ, so the integrand picks up (1+εα)^{−5}, where 5 = Δ_T + J_T. This gives δ_K E = −6E⁽¹⁾; the boundary term is O(\|α\|^{−5}).
* K annihilates the primary T at the origin, so G_{E⁽¹⁾} = (1/6)K_n^𝒪 G, with the explicit position-space generator.
* Fourier transform, then go to the rest frame with spinor polarisations ζ(a,b) and generators B₁, B₂, R. This gives K̂/i = Δ²−3Δ−16 + B₁²+B₂²+R² + (4−2Δ)B₁ − 2RB₂.
* Normalisation of the momentum states and the translation part of E⁽¹⁾ supply +18. So M = (i/6)[T + B₁²+B₂²+R² + (4−2Δ)B₁ − 2RB₂]·H.

**Checks** (`code/angular/verify_ward.py`, exact sympy):

| step | check |
|---|---|
| The rest-frame operator T + B₁²+B₂²+R² + (4−2Δ)B₁ − 2RB₂ on e_m = a^{J+m}b^{J−m} is tridiagonal with diagonal T + C − 2m² | **verified exactly**, J = 0..6. The diagonal holds for all J: B₁²+B₂²−R² is the Casimir C and R² = −m² |
| Unnormalised neighbour coefficients: coef(e_m in We_{m+1}) = −(Δ−m−2)(J+m+1) and coef(e_m in We_{m−1}) = −(Δ+m−2)(J−m+1) | **verified**, J ≤ 3 printed (pattern exact); general J follows from the first-order action of B₁, B₂ on monomials |
| With the all-J norm ratio (§4), the matrix elements between normalised states are exactly −b_m = −√((C−m(m+1))(T−m(m+1))) | **verified exactly**, J = 0..6, at a rational Δ |
| The origin-centred constant Δ²−3Δ−16 | **matches** ANGULAR_GO's independent J=0 origin-centred value i(Δ²−3Δ−16)/6 |
| (1+εα)^{−5} weight, hence the −6 coefficient | consistent with weight Δ_T + J_T = 5 |
| Fourier transform plus rest-frame reduction of K_n^𝒪 to the B/R form, and the +18 from normalisation and translation | **not independently re-derived**; the J = 0 normalised value iT/6 agrees (ANGULAR_GO, with the normalisation sign [D-sketch]). The ladder agent is checking it against the pipeline for J = 1..8 |

**Status.**
* *All-J first-moment identity:* [D] with its algebraic core [P].
* *Still open:*
  * the all-J ladder ratios (b = 0, 4);
  * the Y and Z kernels for general rectangular H (author has it; details to come);
  * interchange and convergence;
  * double-trace and contact terms;
  * parity-odd α_T;
  * recomputing 47266.

## 8. All-J helicity ladder from light-transform conservation (author-supplied, 25 Nov)

**Derivation (author):**
* 𝕃[T](x,z) = ∫dα(−α)^{−5} T(x − z/α, z) has transformed spin −2. Conservation and tracelessness of T give ∂_x·D_z 𝕃[T] = −½∫dα ∂_α[(−α)^{−5}□T(...)] = 0, with no endpoint terms in separated, regulated 3-pt matrix elements.
* Mapping the base point to infinity turns translations into K, giving D_z·K̂^𝒪 H = 0.
* With the §7 generators this becomes, for the source helicity +2 unnormalised coefficients F_m:
  m(m−4)(T+C−2m²)F_m = (m−3)(m−4)(Δ−m−2)(J+m+1)F_{m+1} + m(m−1)(Δ+m−2)(J−m+1)F_{m−1}.
  These are exactly the §7 diagonal and neighbour coefficients, weighted.
* The momentum Ward identity gives F₁ = F₂ = F₃ = 0, and the recurrence decouples at m = 0 and m = 4. With endpoints F_{±(J+1)} = 0:
  F_{−k} = F₀·C(J,k)(4)_k/(Δ−1−k)_k and F_{4+k} = F₄·C(J−4,k)(4)_k/(Δ−5−k)_k.
* The author checked 11,244 exact recurrence identities for J ≤ 250.

**Checks** (`code/angular/verify_ladder.py`):

| step | check |
|---|---|
| Both closed forms satisfy the recurrence at every m, with F₁ = F₂ = F₃ = 0 and F_{±(J+1)} = 0 | **verified exactly**, J = 4..10, symbolic Δ |
| Normalised a = F/√N_m (§4 norms, symmetric in m) gives a_{k+1}/a_k = ((k+4)/(k+1))√((C−s(s+1))/(T−s(s+1))), s = b+k, for b = 0, 4 | **verified**, J ≤ 40, two Δ each (Δ near unitarity and Δ = J+1000). a = F/√N is the same convention as §7's normalised operator |
| Uniqueness by solving from the endpoints (coefficients nonzero for long unitary reps) | argument fine |
| Light-transform conservation ⇒ D_z·K̂H = 0, and that this produces exactly the m(m−4), (m−3)(m−4), m(m−1) weights | **not independently re-derived.** The repository's symbolic ladder through J = 12 (KERNELS_D3) agrees with these closed forms. The ladder agent is extending the pipeline comparison to J ≈ 30 and auditing this step |

**Status.**
* *All-J ladder:* [D] with its solution [P].
* *Remaining for the unconditional 0.0224:*
  1. audit of the conservation/chart step;
  2. the Y and Z kernels for rectangular H;
  3. interchange, convergence, double-trace and contact terms;
  4. parity-odd α_T;
  5. recompute 47266.

## 9. All-spin positivity of the angular form, a direct endpoint sum rule, and the spectral-tail obstruction (author-supplied, 25 Nov)

Coordinator checks are in `code/angular/verify_s9.py`.

| claim | check |
|---|---|
| On ladder vectors (§8): Σw_k Re[ā_k(Aa)_k] = w₀[T−b(b−1)]\|a₀\|² + Σ_{k<J−b}(ρ_k−1)(ρ_k w_{k+1} − w_k)X_k\|a_k\|², with ρ_k = (k+4)/(k+1) | **verified numerically** for random weights (max rel. err. 2e-14, 300 random b, J, Δ) |
| With w_k = −z^raw_{k+2}: each summand ≥ 0 by §5's ratio inequality, so 𝓗 ≥ (5/48)[T−b(b−1)]\|a₀\|² for all J | follows (all pieces verified) — **all-spin positivity of the reduced form** |
| q_m = 6/((m−1)m(m+1)) (m ≥ 2), q₀ = −3, q₁ = 0, q_{−m} = q_m: q_{k+3}/q_{k+2} = (k+1)/(k+4), so the interior telescopes away: Σq_{k+2}Re[ā(Aa)] = [T−b(b−1)]\|a₀\|² **exactly** | **verified** (identity numerically, ratio by inspection) |
| Ỹ_m = q_m + (1−m/2)q_{m−1} + (1+m/2)q_{m+1}: Ỹ₀ = −3, Ỹ₂ = 3/2, Ỹ_m = 0 for m ≥ 3 | **verified exactly** (symbolic m ≥ 3); also Ỹ₁ = 0 |
| Candidate direct sum rule from [E(n₃), E⁽¹⁾(n₁)] = 0: Σ_{O≠T,b}[(T_Δ − b(b−1))/3 + 3/2]\|a_{O,b}\|² = 3 | **derivation not checked**: the physical normalisation (the /3, the RHS 3 from Ỹ₀) is the author's |
| Combined with the fourth-harmonic rule: (T−12)/3/W₂⁴ > 151250 for twist ≥ 1000, J ≥ 6 ⇒ √(t₄²+α_T²) < 0.01782 | **constant and arithmetic verified** (151252.98; 0.017814). **Conditional** |
| Spectral-tail obstruction: truncation at harmonic N+2 leaves B_N = q_{N+2}Σ_{O,b}X_N\|a_{O,b,N}\|²/(N+1). This vanishes per representation, but the infinite sum needs Σ J(J+1)\|a_{O,b,N}\|² = o(N⁴) | **genuine open problem** (author-identified). A positive limit of B_N would enlarge the RHS and break the bound |

**Coordinator's note — the main physics risk.** The gap applies to single traces only. At leading order in 1/c_T, multi-trace states (e.g. [TT]_{n,J}, twist ≈ 2Δ_T + 2n) also enter the heavy sum with **low** T_Δ. They are positive, but they could supply fourth-harmonic weight cheaply, and the bound's W₂⁴ estimate assumes twist ≥ 1000. The original Ω^EE handles them with the [TT] double zeros of the de-shift weight; this direct angular sum rule has no such mechanism visibly. The open list must include either their exclusion at the relevant order or explicit control. Same for contact terms in the mixed E, E⁽¹⁾ commutator and the uniform large-c_T limit.

## 10. Spectral boundary reduced to a two-norm condition; ANEC fourth-harmonic rule; contact kernel (author-supplied, 25 Nov)

Coordinator checks are in `code/angular/verify_s10.py`.

| claim | check |
|---|---|
| Fourth-harmonic rule from [E₄, E₋₄] = 0 in a normalised \|T,+2⟩: Σ_{O,J≥6}\|a_{O,4,2}\|² = \|c\|² + Σ_{O≠T}\|a_{O,0,2}\|² ≥ (t₄²+α_T²)/16 | **argument not re-checked.** Caveat: the Fourier modes of E(n) involve coincident directions. KKSZ commutativity holds at separated n, so [E₄, E₋₄] = 0 also needs the absence of a contact term in [E(n), E(n′)] at n = n′, or its cancellation in this combination |
| Bounded-twist tail lemma: Σ_{J≥L, Δ−J≤g}\|a_{O,b,2}\|² = O(g⁶/L⁶), uniform in c_T, from the ladder plus Σ_O\|a_{O,b,M}\|² ≤ 1 | argument plausible, not re-checked. Multi-trace large-spin escape is controlled. **Fixed-window planar suppression of multi-trace weights is still assumed** |
| On ladders: d₀ = (Aa)₀/a₀ = T − 3C + 2b² + 4b, d₁ = (3/4)T − (3/2)C + (3/4)(b+1)(b+4) | **verified** (rel. err. 6e-16, 2000 random b, J, twist) |
| Coercivity C²\|a₀\|² ≤ (25/4)(\|(Aa)₀\|² + \|(Aa)₁\|²) for twist > 2 | **verified numerically**: the minimum ratio over 2000 random cases is 4.63 (claim ≥ 1) |
| Theorem: if Σ_{O,b}(\|(Aa)₀\|² + \|(Aa)₁\|²) < ∞, then B_N → 0 (twist ≥ 2 part via coercivity; bounded-twist part via the dyadic ladder comparison, O(g²/N²)) | logic follows from the verified pieces; bounded-twist part not re-checked |
| Abstract counterexample (L_n = 2ⁿ, J_n = L_n³, Δ_n = J_n + L_n² + 2) with B_{L_n}⁽ⁿ⁾ → e⁻²/6: ladder plus bounded harmonic norms do not imply convergence | not re-checked; plausible. It shows the two-norm condition is genuinely needed |
| Contact kernel Q(θ) = Σq_m e^{imθ} = (1 − cos θ)[12 log(2\|sin θ/2\|) − 9]; Q(0) = Q′(0) = 0, not C² | **verified** (series vs closed form, three angles, 15 digits). It annihilates δ and δ′ contacts but not higher-derivative ones |
| K(g) = (729/4800)(g−2)²(1 − 49/(2g+9)²); bound √(t₄²+α_T²) ≤ 4√(3/K(g)) → (160/9)/g | **verified**: K(1000) = 151266, bound 0.017814, asymptotic 0.017778. Min of (v+g−2)⁶/v⁴ checked numerically |

**Open (author's list, plus the coordinator's caveat on [E₄, E₋₄]):**
1. **Physical finiteness of the two first-moment norms** Σ(\|(Aa)₀\|² + \|(Aa)₁\|²). Weighted light-ray commutativity is known to fail in some perturbative holographic limits, so this cannot be taken from ANEC alone.
2. **The endpoint sum rule:** normalisation and contact matching (higher-derivative contacts).
3. **Fixed-window planar suppression** of multi-trace transition weights.
4. **The coincident-direction contact term** in [E₄, E₋₄].

## 11. Weaker spectral-convergence condition (author-supplied, 25 Nov)

**Theorem (twist ≥ 2 sector).**
* *Hypotheses:* the ladder holds (§8) and the normalised ANEC harmonics are bounded, Σ_O\|a_{O,b,N}\|² ≤ 1.
* *Statement:* if Σ_{O,b}C^{3/2}\|a₀\|² < ∞, then Σ_{O,b}C\|a_N\|² = o(N⁴), hence B_N → 0 (§10, using X_N ≤ C and q_{N+2}/(N+1) = 1/(D_N(N+1)) ~ 6/N⁴).
* *Corollary:* with §10's coercivity, it suffices that Σ(\|(Aa)₀\|² + \|(Aa)₁\|²)/√(1+C) < ∞.

**Coordinator check of the proof [P]:**
* T ≥ C, so each ladder ratio is ≤ ((k+4)/(k+1))². The product gives \|a_N\|² ≤ D_N²\|a₀\|², with D_N = (N+1)(N+2)(N+3)/6.
* Split at J = εN². Below the split, C ≤ ε²N⁴(1 + O(N⁻²)) and the ANEC bound give ≤ ε²N⁴. Above it, C\|a₀\|² ≤ C^{3/2}\|a₀\|²/(εN²), so the region contributes ≤ (D_N²/(εN²))·(tail of a convergent sum) = o(N⁴).
* Then ε → 0. The corollary follows from C²\|a₀\|² ≤ (25/4)(\|(Aa)₀\|² + \|(Aa)₁\|²). **Correct.**

**Effect.** The unweighted first-moment norm is no longer needed. The open physical input is now a spin-weighted moment of the zeroth-harmonic (endpoint) amplitudes, Σ C^{3/2}\|a₀\|² < ∞.

## 12. Endpoint normalisation, contact terms, multi-traces, and the Regge-convergence obstruction (author-supplied, 25 Nov)

Coordinator checks are in `code/angular/verify_s12.py`; everything else is reasoning, with comments.

| claim | check / comment |
|---|---|
| **Sharpness of the C^{3/2} condition.** Abstract ladder family J_L = L², Δ_L = L²+L+2, \|a_{L,0}\|² = L⁻⁶ (L = 2ⁿ): ΣC^p\|a₀\|² < ∞ for all p < 3/2, but limsup N⁻⁴ΣC\|a_N\|² ≥ e⁻²/36 | not re-checked; plausible. It makes the exponent 3/2 sharp for ladder + ANEC alone |
| **Endpoint rule.** With E⁽¹⁾ = −(i/2){H, n·∂_p} + M and M = (i/6)AH, [H₁, E⁽¹⁾₂] = (i/2){H₂, D₂H₁} + [H₁, M₂]. Boost weight 2 gives D₂H₁ = (1+2cos θ)H₁ + sin θ ∂_θH₁. Per harmonic: \|a_m\|²[(1+2cos θ)cos mθ − m sin θ sin mθ] + (1/3)Re(ā(Aa))cos mθ. With q_m: Y₀ = −3, Y₂ = 3/2, Y_{m≠0,±2} = 0 (verified §9); the interior telescopes (verified §9); the stress-tensor term is −3 and the cubic term is absent (Y₄ = 0). So Σ[(T−b(b−1))/3 + 3/2]\|a₀\|² = 3 | derivation consistent with the verified pieces; the 1/3 and the RHS 3 are now derived, not normalised by hand. **Not independently re-derived** (especially the "reduced first moment vanishes within the T family" step) |
| **[E₄, E₋₄] = 0.** The normalised ANEC is a POVM, 0 ⪯ E(A) ⪯ I. Separated-direction commutativity gives commutativity for disjoint Borel sets (disjoint compacts are separated, then inner regularity). For general A, B, write E(A) = E(A∩B) + E(A∖B) and E(B) = E(A∩B) + E(B∖A); every cross pair is identical or disjoint. So all E(A) commute, and so do their Fourier integrals | **correct** (standard; it needs E to be a genuine bounded POVM on the fixed-energy fibre). This **resolves** the coincidence caveat of §10 |
| d_k = (Aa)_k/a_k = 3[T−(b+k−1)(b+k)]/(k+3) − 3[C−(b+k)(b+k+1)]/(k+1) on ladders | **verified** (rel. err. 6e-16) |
| **Contact terms.** Under the C^{3/2} condition, Σ\|d_k\|\|a_k\|² = o(k²), so the Fourier coefficients of the even mixed commutator are o(m²). A distribution supported at coincidence is Σc_rδ^{(2r)} with polynomial-in-m² coefficients, so only c₀δ survives, and Q(0) = 0 kills it | argument plausible (not re-checked in detail). Higher-derivative contacts are then excluded **by the same spectral condition** — no separate classification needed |
| **Fixed-window multi-traces.** a_{O,m} = √c_T λ_{ttO}K_m(Δ,J) with K_m(Δ_*,J) = 0 at GFF dimensions, γ = O(1/c_T) and λ = O(1), so \|a_{[TT]}\|² = O(1/c_T). Combined with the bounded-twist tail (§10), the whole low-twist fourth-harmonic contribution 𝓛_g(c_T) → 0 | **consistent** with KERNELS_D3 §5.7 double-trace zeros. **Caveat (author):** degenerate multi-trace mixing with γ ~ c_T^{−1/2} would break the Taylor estimate. That needs checking for TT-built multi-traces |
| **Fractional regularity.** The C^{3/2} condition follows from E₂ψ_T ∈ Dom((K⁰)^{3/4}), via H_c = (P⁰+K⁰)/2 ≥ Δ, dilatation conjugation and spectral Jensen; sufficient: ‖(e^{isK⁰} − I)E₂ψ_T‖² = O(s^{3/2+ε}) | formal reduction plausible; **open** — this is now the central physical estimate |
| **Regge convergence (KKSZ 1904.05905 §4.6).** Weighted light rays commute in the Regge limit if J₁+J₂ > J₀+3+n+m (author's reading). For T detectors: EE needs J₀ < 3, EE⁽¹⁾ needs J₀ < 2, E⁽¹⁾E⁽¹⁾ needs J₀ < 1 | **coordinator note:** the citation needs checking. If the table is right, the mixed commutator [E, E⁽¹⁾] that this route uses needs only J₀ < 2. The nonperturbative CFT bound J₀ ≤ 1 satisfies it with room; at tree level in 1/c_T the graviton gives J₀ = 2, exactly marginal. So the issue is the order of limits (c_T → ∞ vs Regge). The k=2 Ω^EE also uses a first light-ray moment (the affine weight), and its order-of-limits question (M2) was resolved in the existing record — **that resolution should be checked for transfer to the angular sum rule** |

**Status.**
* *Derived algebraically:* the endpoint normalisation 3.
* *Resolved:* the ANEC coincidence.
* *Reduced to the spectral condition:* the contact terms.
* *Controlled modulo degenerate mixing:* multi-traces.
* *Single remaining physical estimate:* ΣC^{3/2}\|a₀\|² < ∞ (or the fractional K⁰ regularity), together with the c_T-vs-Regge order of limits for the mixed commutator.

## 13. Contact terms under a weaker condition, the exact defect theorem, the critical tail criterion, and degenerate mixing (author-supplied, 25 Nov)

Coordinator checks are in `code/angular/verify_s13.py`.

| claim | check / comment |
|---|---|
| Contacts need only ΣC\|a₀\|² < ∞. Then S_C(k) = o(k⁶), Σ(T−C)\|a_k\|² = o(k⁵), Σ\|ā_k(Aa)_k\| = o(k⁴), so Ĉ(m) = o(m⁴), leaving c₀δ + c₂δ″. The Ward relation Ĉ(0) = Ĉ(1) (heavy harmonics start at \|m\| = 2; the 2→1 coefficient is 1 − 2/2 = 0) forces c₂ = 0 | growth chain consistent with the verified d_k (§12). **The Ĉ(0) = Ĉ(1) step is not re-checked** |
| Truncated kernel: Q_N(0) = −6/((N+2)(N+3)) | **verified exactly** (N < 60, by telescoping) |
| Truncated first moment: (1/3)Σ_{k≤N} ā_k(Aa)_k/D_k = [(T−b(b−1))/3]\|a₀\|² − B_N, B_N = X_N\|a_N\|²/((N+1)D_N) | **verified numerically** on ladders (rel. err. 7e-10) |
| **Defect theorem.** Given only a δ contact and finite 𝒫: 𝒫 − 3 − B_N = O(N⁻²), so 𝒟 = lim B_N = 𝒫 − 3 ≥ 0 exists and 𝒫 = 3 + 𝒟. Gravitational inequality: t₄² + α_T² ≤ 16(3+𝒟)/K(g) + 16𝓛_g | structure follows from the verified pieces. **This turns the convergence question into an explicit, signed correction** |
| Critical tail criterion: L·F(L) → 0, with F(L) = Σ_{J>L, twist≥2} C\|a₀\|², implies S_C(N) = o(N⁴), hence 𝒟 = 0. Strictly weaker than ΣC^{3/2}\|a₀\|² < ∞ | proof follows §11's split; correct |
| Physical sufficient condition: ‖1_{K⁰>R}v‖² = o(R^{−3/2}), or ‖(e^{isK⁰} − I)v‖² = o(\|s\|^{3/2}), for v = E₂ψ_T | formal reduction plausible; **open: the central physical estimate** |
| **Counterexample.** J_ℓ = 2^ℓ, Δ_ℓ = J_ℓ + ⌈J_ℓ/ℓ⌉ + 2, \|a_{ℓ,0}\|² = 1/(J_ℓ(J_ℓ+1)ℓ³). Its harmonic norms decay faster than any polynomial, yet lim B_N = (1/6)∫x⁻³e^{−2/x}dx = 1/24 | **integral verified** (= 1/24); construction not re-checked. Smoothness of ordinary energy correlators cannot close the gap: CFT dynamics must enter |
| **Degenerate mixing.** a_ε = ε⁻¹K(D_ε)v_ε gives O(1) terms K(Δ_*)v₁ + K′(Δ_*)Vv₀. Basis-independent sufficient conditions: Vv₀ = 0 and K(Δ_*)v₁ = 0, which give \|a\|² = O(1/c_T) even with O(c_T^{−1/2}) splittings. Factorisation route: at first order, the only surviving one-detector transition to a Wick-ordered multi-particle bra is a single 3-pt function (detector, external T, one constituent); the rest vanish by ANEC vacuum annihilation | perturbative algebra correct. Factorisation argument plausible; **analytic justification (differentiable spectral projections, light transform vs projected expansion) open** |
| **Literature (author).** KKSZ give J₀ ≤ 1 nonperturbatively, and their eq. (4.87) puts the mixed first-moment product in the J₀ < 2 range. Weighted light-cone convergence and momentum eigenstates remain to be justified. Córdova–Shao's contact algebra uses an extra closure hypothesis | citations to be checked before use |

**Status of the angular route.**
* *Bound:* √(t₄²+α_T²) ≤ 4√((3+𝒟)/K(g) + 𝓛_g), which is 0.01782 at g = 1000 when 𝒟 = 𝓛_g = 0.
* *𝒟 = 0* follows from the single tail estimate ‖1_{K⁰>R}E₂ψ_T‖² = o(R^{−3/2}).
* *𝓛_g → 0* follows from fixed-window planar suppression, now reduced to Vv₀ = 0 and K(Δ_*)v₁ = 0.
* *Also open:* existence and separated vanishing of the physical mixed commutator, and the order of limits (c_T vs Regge, light cone).

## 14. Subcritical-twist sectors contribute no spectral defect (author-supplied, 25 Nov)

**Theorem.**
* *Hypotheses:* the ladder (§8, audited), normalised ANEC boundedness, and a sector with τ = Δ − J ≤ A(1+J)^α, 0 ≤ α < ½.
* *Statement:* N⁻⁴Σ_sector C\|a_N\|² = O(N^{−(2−4α)/(1−α)}) → 0.
* *Refinement:* limsup N⁻⁴Σ_{τ≤ε√J}C\|a_N\|² ≲ ε⁴. So any sector with τ = o(√J) has vanishing defect; this includes bounded twist and logarithmic trajectories.

**Coordinator check of the proof [P]:**
* **Per shell.** On a dyadic shell L ≤ J < 2L with τ ≤ g_L ≤ J, T − C = (τ−2)(2J+τ−1) ≤ 3g_L J. For k < M_L = ⌊L/(12g_L)⌋ each ladder ratio is ≥ ((k+4)/(k+1))²(1+6g_L/L)^{−1}. So for N ≤ M_L, \|a_{M_L}\|² ≥ (D_{M_L}/D_N)²·exp(−(M_L−N)·6g_L/L)\|a_N\|² ≥ e^{−1/2}(D_{M_L}/D_N)²\|a_N\|².
* **ANEC at the higher harmonic** (Σ\|a_M\|² ≤ 1) gives Σ_shell C\|a_N\|² ≲ L²·min{1, (Ng_L/L)⁶}.
* **Dyadic sum.** With g_L ~ L^α, split at L_* ~ N^{1/(1−α)}. Below L_*: ΣL² ~ N^{2/(1−α)}. Above: N⁶ΣL^{6α−4} ~ N^{2/(1−α)} (convergent since 6α − 4 < 0). Dividing by N⁴ gives the stated rate, negative for α < ½.
* **Refinement.** With g = ε√(2L), L_* ~ ε²N², and both parts are ~ε⁴N⁴.
* **Finitely many shells.** They contribute nothing for large N, since a_N = 0 for N > J − b.

**Correct.**

**Effect.** All states at bounded twist, in particular every heavy state near the gap at any spin, are excluded as a source of 𝒟. The defect question is confined to τ ≳ √J.

## 15. Finite-harmonic obstruction, the comparable-twist sector, the R-tail criterion, the KZ contact, and no-log mixing (author-supplied, 25 Nov)

Checks are in `code/angular/verify_s15.py`.

| claim | check / comment |
|---|---|
| **Finite-harmonic obstruction.** For finitely supported weights w, on large-spin ladders with C/T → x ∈ (0,1): P_w(x) = 3Σw_k D_k²[x^k/(k+3) − x^{k+1}/(k+1)] satisfies ∫₀¹(1−x)P_w = 0. So P_w ≥ 0 forces w = 0. The endpoint weights w_k = 1/D_k give P_w ≡ 1, a telescoping series that converges non-uniformly at x → 1 | **verified** (termwise integrals vanish exactly; telescoping exact). **The slowly decaying endpoint kernel is essential: no smoother certificate exists in this construction** |
| **Comparable twist.** For T > C, \|a_N\|² ≤ D_N²(C/T)^N\|a₀\|², since (C−y)/(T−y) ≤ C/T. So B_N ≤ ((N+2)(N+3)/6)ΣC\|a₀\|²e^{−N/R}, with R = T/(T−C). If τ ≥ δJ then R = O(1): exponential decay, given a finite endpoint budget (which bounds ΣC\|a₀\|²/3) | **verified** (inequality symbolic; the rest is direct) |
| **R-tail criterion.** For μ = ΣC\|a₀\|²δ_R: N²∫e^{−N/R}dμ = ∫e^{−u}N²μ((N/u,∞))du, so μ((L,∞)) = o(L⁻²) implies 𝒟 = 0 (dominated convergence, given sup L²μ((L,∞)) < ∞) | **correct** [P] |
| **Contact.** Korchemsky–Zhiboedov's finite-frequency algebra (eqs. 1.13–1.15) at ω → 0 gives [E₀(n), E_ω(n′)] = ωδ(n,n′)E_ω(n). Hence [E₀(n), E⁽¹⁾(n′)] = iδ(n,n′)E₀(n), which Q(0) = 0 annihilates | consistent with §§12–13. It assumes their algebra holds and that no relevant light scalars appear (an extra restriction). **Citation to be checked** |
| **No-log mixing.** With D_ε = E_*I + εV, r^{D_ε} = r^{E_*}[I + εV log r + …]. First-order factorised correlators of single-trace products have no radial logs if constituent dimensions do not shift at first order, so V = 0. With the kernel zero, \|a_{[TT]}\|² = O(1/c_T) | plausible [D]; trace-changing directions are still to be verified |

**Remaining spectral target, now narrow.** Σ_{τ ≳ √J, T/(T−C) > L} C\|a₀\|² = o(L⁻²). The bounded-twist, subcritical-twist (§14) and comparable-twist sectors are all done.

## 16. The intermediate band as an exact Tauberian problem (author-supplied, 25–26 Nov)

**Setting.**
* *Band:* τ ≥ κ√J and τ/J → 0. Here R = T/(T−C) ≍ J/(2τ), so J ≳ R².
* *Measure:* dμ(R) = ΣC\|a₀\|²δ_R.
* *Boundary:* B_N = ((N+2)(N+3)/6)ΣC\|a₀\|²𝓡_N, with 𝓡_N = (X_N/C)∏_{k<N}(C−s(s+1))/(T−s(s+1)), where s = b + k.

**Claims.**
* 𝓡_N = e^{−N/R}(1 + O(1/N)) uniformly for R/N in compacts.
* B_N → 0 ⇔ μ((L,∞)) = o(L⁻²).
* L²μ((L,∞)) ≲_κ B_{⌊L⌋}.
* If the finite-cutoff commutator gives a finite defect 𝒟, then L²μ((L,∞)) → 3𝒟.

**Coordinator check [P].**
* **Uniform asymptotics.** ∏(C/T) = (1 − 1/R)^N = e^{−N/R}(1 + O(N/R²)). The finite-s corrections total ~Σs²/(RC) ~ N³/(RJ²) ≲ N²/J², which is negligible since J ≳ R² ~ N². X_N/C = 1 − O(N²/J²).
* **Necessity.** For R ≥ N, 𝓡_N ≥ e⁻¹(1 − o(1)). Then B_N ≥ (N²/6)·c·μ((N,∞)).
* **Factor 3.** B_N ≈ (N²/6)∫e^{−N/R}dμ = (N²/6)∫e^{−u}μ((N/u,∞))du. With tail A/L², ∫e^{−u}Au²/N²du = 2A/N², so B_N → A/3. Hence A = 3𝒟 (verified symbolically).
* **Reading.** The last spectral question is now sharp: exclude a critical L⁻² tail, ~3𝒟/L², of the spin-weighted endpoint measure in the intermediate band.

**Order of limits (author).** Establish the spectral identity in each exact finite-c_T CFT (where J₀ ≤ 1 < 2), and only then take c_T → ∞. This is consistent with the coordinator's §12 note; the matching agent is checking the M2 analogue.

## 17. Null-plane Ward route to the spectral estimate; M2 transfer; trace-changing multi-traces; light scalars (author-supplied, 26 Nov)

| claim | check / comment |
|---|---|
| **Double commutator.** From the light-ray algebra [N(g), A(f)] = iA(gf′) (Córdova–Shao) and conservation [H, A(f)] = (i/2)N(f′), with H the positive null-translation generator: [A(f), [H, A(f)]] = ½A((f′)²) | **algebra verified**. The inputs are the formal identities |
| **Quadratic-form identity.** [A,[H,A]] = 2AHA − A²H − HA², hence 2‖H^{1/2}Aψ‖² = 2Re⟨Aψ, AHψ⟩ + ½⟨ψ, A((f′)²)ψ⟩. Domain via a regulated time-average with ρ̂_η ≥ 0 and Fatou, so no circular assumption of H-regularity | expansion verified symbolically; the regulated-limit/Fatou argument is plausible. **Needs the Ward identity as a regulated quadratic-form statement — open** |
| **Harmonic profiles.** f_m = (1+y²)e^{imθ}, y = tan(θ/2): \|f_m′\|² = 4(y²+m²) ≤ 4(1+m²)(1+y²), so the RHS is controlled by total energy | **verified** (sympy). Note A(f_m) is not Hermitian; use the cos/sin parts separately |
| **K⁰ estimate.** Two stereographic charts give ‖(K⁰)^{1/2}A_mψ_m‖² = O((1+\|m\|)²), with ψ_m localised within O(1/\|m\|) of unit rest momentum | plausible. **Caveat:** passing from wavepackets to fixed-momentum normalised amplitudes a_{O,b,N} needs their smoothness on the 1/m momentum scale |
| **Family bound.** ⟨P⁰⟩⟨K⁰⟩ ≥ Δ²‖w‖⁴ (conformal Hamiltonian ≥ Δ, dilatation conjugation), and E commutes with P, so Σ_{O,b}Δ²\|a_N\|² = O(N²) ⇒ ΣC\|a_N\|² = O(N²), B_N = O(N⁻²), μ((L,∞)) = O(L⁻⁴) | logic correct given the K⁰ estimate. **Strictly stronger than the o(L⁻²) needed (§16)** |
| **M2 transfer.** The letter already uses the exact finite-c_T sum rule, expanding only the light side and keeping the heavy side by positivity. For the angular rule: lim_{c_T→∞}[lim_N 𝒫_N(c_T)], with J₀ ≤ 1 < 2 at finite c_T. The spectral estimate also gives absolute convergence: \|(Aa)_N/a_N\| ≲ Δ²/(N+1) and q_{N+2} = O(N⁻³) | consistent; matching agent's independent check pending |
| **Trace-changing multi-traces.** With ε = c_T^{−1/2} and E = ε⁻¹𝕃[t]: at O(ε⁰) the detector sits in a 2-pt contraction (vacuum annihilation). At O(ε¹) the three cases vanish by vacuum annihilation or Wick ordering. So ‖Π_multi E_m\|T⟩‖² = O(1/c_T), basis-independently. With the bounded-twist tail, 𝓛_g → 0 | argument plausible [D], given differentiable fixed-window planar expansion under spectral projection and light transform |
| **Light scalars.** The KZ finite-frequency algebra needs Δ_φ > 1 for every scalar in T×T (d = 3), plus the higher-spin-current restrictions. A higher-spin gap alone does not imply this, so it must be **stated as an assumption** | adopted |

**Result if the regulated null-plane identity holds** (plus physical matching, fixed-window planar regularity and Δ_φ > 1): t₄² + α_T² ≤ 48/K(g), i.e. √(t₄²+α_T²) < 0.01782 at g = 1000.

**Single remaining operator-algebra target:** the regulated quadratic-form version of [A(f), [H, A(f)]] = ½A((f′)²) at finite c_T.

## 18. Final audit: null-plane sign, d=3 contact uniqueness, regulator theorem, wavepacket → fixed momentum (author-supplied, 26 Nov)

| claim | check / comment |
|---|---|
| **Sign and normalisation.** ds² = −du dv + dy². Conservation gives ∂_u𝓔 = ½∂_y𝓝, so [H, A(f)] = (i/2)N(f′). Córdova–Shao (3.24): (1/i)[𝓝(y), 𝓔(z)] = −δ∂_z𝓔 + δ′𝓔, so [N(g), A(f)] = iA(gf′). Hence [A, [H, A]] = +½A((f′)²) | algebra consistent (§17). Citation equation number to be confirmed |
| **d = 3 contact uniqueness.** R = physical − canonical commutator, supported at y = z, R = Σa_r∂^rδ. The H-Jacobi gives ∂_yR(y,z) = ∂_zR(z,y), so r_max is odd. The null-rotation Ward identity gives [4 − (r_max+1)]a_{r_max} = 0, so r_max = 3. Translation gives ∫R dy = 0, leaving R = aδ‴ − 2a′δ″ + a″δ′, i.e. R(g,f) = −⟨a, g′f″⟩. The dilatation and transverse SCT Ward identities (g = y, y²) give a″ = 0 and (ya)″ = 0, so a′ = 0. Locality plus boost charge gives a = 0 | **Step C/D algebra verified** (sympy: the smeared form, and g = 1, y, y² giving a″ = 0, (ya)″ = 0 ⇒ a′ = 0). Steps A/B (Jacobi + null rotation) and the boost-charge step are plausible, not re-derived. **This removes the local-operator closure assumption of Córdova–Shao**, but needs the regulated commutators to exist as local distributions satisfying the Ward/Jacobi identities |
| **Regulator theorem.** I_η = −∫ρ_η′⟨ψ,[A, A_t]ψ⟩/i = 2∫λρ̂(ηλ)dν_v − 2Re⟨v, AHψ⟩ + o(1). With ρ̂ ≥ 0, Fatou gives: liminf I_η < ∞ ⇒ Aψ ∈ Dom(H^{1/2}). An **upper bound** suffices for convergence; exact equality is not needed | correct [P] (spectral theorem + Fatou) |
| **Wavepacket → fixed momentum.** Cutoff f_{m,R}, uniform bound \|f′\|²/(1+y²) ≤ C(1+m²); no atoms in the angular POVM; two charts give K⁰. Packets of width δ_m = ε/(1+\|m\|) give ‖K⁰ψ_m‖ = O(m²). The Doppler/aberration errors are O(\|m\|δ_m) = O(ε), so the target Fourier coefficient stays ≥ ½. Coercivity ⟨P⁰⟩⟨K⁰⟩ ≥ Δ²‖w‖⁴ gives Σ_OΔ²\|a_{O,m}\|² = O(m²) | argument sound. It resolves §17 caveat 2 |
| **Consequences.** ΣC\|a_N\|² = O(N²), B_N = O(N⁻²), 𝒟 = 0, and absolute angular swapping (Σ_k q_{k+2}Σ\|ā(Aa)\| < ∞). Endpoint rule with RHS 3 | follows |
| **Theorem.** limsup_{c_T→∞}√(t₄²+α_T²) ≤ (160/9)/(g−2)·(1 − 49/(2g+9)²)^{−1/2} | **verified**: 0.017814 at g = 1000; same as 4√(3/K(g)) |

**Assumptions of the angular theorem** (final list):
1. *Physical matching:* the ladder, the first-moment matrix and the mixed-commutator normalisation are the physical ones. The algebra is audited (LADDER_ALLJ); it rests on the operator identities [z·K, E] = 6iE⁽¹⁾ and D_z·[K, E] = 0.
2. *The regulated null-plane identity at finite c_T:* specifically liminf_η I_η < ∞ with the profile bound, plus Ward/Jacobi compatibility for the contact uniqueness. **This is the single open nonperturbative statement.**
3. *Fixed-window planar regularity* for multi-traces.
4. *Light scalars:* Δ_φ > 1 for scalars in T×T; conserved higher-spin sectors excluded.

**Order of limits:** finite c_T first, then the planar limit.

## 19. Matching audit (ANGULAR_MATCH.md, 26 Nov) — amendments to §18's assumption list

* *Resolved [P]:*
  * the constant 0.0178135;
  * Ĉ(0) = Ĉ(1);
  * the endpoint 1/3;
  * the first moment vanishing within the T family.
* *Resolved [D]:*
  * the M2 transfer;
  * TT-built multi-traces, in the strong form (no O(c_T^{−1/2}) splittings).
* **Add to assumptions:** separated-direction vanishing of the mixed [E, E⁽¹⁾] at the light-cone edge. This uses KKSZ's light-cone expansion (4.66), which they leave unproved; with Rindler positivity alone the mixed product is exactly marginal. With (4.66) the margin is ≥ ½, or 1 when Δ_φ > 1. KKSZ (4.87) is confirmed verbatim.
* 𝒟(c_T) = 0 at finite c_T remains §18 assumption 2 (the regulated null-plane bound). It does not follow from J₀ ≤ 1.

## 20. Attacks on the finite-c_T null-plane bound (author-supplied, 26 Nov; in progress)

| claim | check / comment |
|---|---|
| **Transverse cancellation.** Q = A(1) = P_v commutes with H and with every 𝓔_t(z), since 𝓔 is a v-integral, invariant under v-translations. So C_t(y,z) = (1/i)⟨ψ,[𝓔₀(y), 𝓔_t(z)]ψ⟩ has vanishing marginals, and ∫∫f(y)f(z)C_t = −½∫∫[f(y)−f(z)]²C_t, exactly, at finite regulator | **correct [P]**: the marginals follow from [P_v, 𝓔_t(z)] = 0, and the identity is f(y)f(z) = ½(f(y)²+f(z)²) − ½(f(y)−f(z))². This buys two powers of \|y−z\| at coincidence. **Note:** ‖f′‖_∞ is infinite for f_m (\|f_m′\|² = 4(y²+m²)), so use the cutoff profiles f_{m,R} of §18, with a uniform bound on \|f′\|²/(1+y²) |
| **All-spin collinear identities.** With s = (Δ+J)/2, h = (Δ−J)/2, p(h) = (h−1)(h−2)(h−3): F_{−2,−2} = p(h)(h−4)(λ₊₋ − 2λ₊₊) and F_{−1,−1} = −16p(h)[2(s−1)λ₊₊ + sλ₊₋], for h ≥ 5, J ≥ 4. Conservation equations were checked for symbolic J, Δ; extreme-helicity detector amplitudes summed by Vandermonde | **not re-checked.** Can be tested against `code/kernels/` structures (KERNELS_D3 conventions) |
| **Crossing estimate.** Reflection positivity plus the crossed-channel OPE give the dyadic bound Σ_{s~S, h~U}4^Δ2^{−J}\|F_{m₁m₂}\|² ≲ S^{11/2+M}U^{11/2−M}, using the lower bound k_{2s}^{(δ)}(1−S⁻²) ≳ 4^sS^{1/2+2δ}. Inverting gives Σ4^Δ2^{−J}(\|λ₊₊\|²+\|λ₊₋\|²) ≲ S^{3/2}U^{3/2}, with constants allowed to depend on the finite-c_T CFT | **not re-checked; in progress.** Still needed: conversion of the OPE-coefficient bound into a bound on the detector amplitudes a_{O,b,0} (kernel normalisation, including the 4^Δ factors), and hence on the intermediate-band tail μ((L,∞)) = o(L⁻²) |

## 21. Claimed proof of the finite-c_T bound via crossing + Casimir (author-supplied, 26 Nov) — UNDER AUDIT

**Claimed chain:**
1. **Collinear all-spin identities** (h = (Δ−J)/2, s = (Δ+J)/2, p(h) = (h−1)(h−2)(h−3); h ≥ 5, J ≥ 4):
   * F_{−2,−2} = p(h)(h−4)(λ₊₋ − 2λ₊₊);
   * F_{−1,−1} = −16p(h)[2(s−1)λ₊₊ + sλ₊₋];
   * F_{−2,−1} = −2(2h−7)(2h−5)(2h−3)λ_odd.

   Checked by the author with symbolic J, Δ.
2. **Reflection-positive crossing (dyadic).** Σ_{s∈[S,2S], h∈[U,2U]} 4^Δ2^{−J}\|F_{m₁m₂}\|² ≲ S^{11/2+M}U^{11/2−M}, with M = m₁+m₂, via the collinear block lower bound k^{(δ)}_{2s}(1−S⁻²) ≳ 4^sS^{1/2+2δ}. Inverting the identities gives Σ4^Δ2^{−J}(\|λ₊₊\|²+\|λ₊₋\|²) ≲ S^{3/2}U^{3/2} and Σ4^Δ2^{−J}\|λ_odd\|² ≲ S^{5/2}U^{5/2}. The constants may depend on the finite-c_T CFT.
3. **Kernel bounds (Stirling).** \|K_even\|² ≲ 4^Δ2^{−J}√(J+1)√(sh)/Δ^{15/2} and \|K_odd\|² ≲ 4^Δ2^{−J}√(J+1)/(√(sh)Δ^{15/2}). Hence Σ_{shell}Σ_b\|a_{O,b,0}\|² ≲ S^{−5}U², and Σ_{τ≥10}Δ²\|a₀\|² < ∞. Low-spin exceptional sectors are handled separately: kernel O(4^ΔΔ^{−29/2}) against a budget Δ^{11}.
4. **Tails.** μ((L,∞)) ≲ Σ_{S≳L}S⁻¹L⁻² ~ L⁻³; in the intermediate band (S ≳ L²), O(L⁻⁴). Together with the bounded-twist ladder estimate: ΣΔ²\|a_N\|² = O(N³), B_N = O(1/N), so 𝒟 = 0.
5. **Casimir domain estimate.** ‖(K⁰)^{1/2}w_O‖² ≲ ‖∇_p w_O‖² + (1+Δ²+J²)‖w_O‖² for wavefunctions on a compact timelike patch. Then E(F)ψ_T ∈ Dom((K⁰)^{1/2}), so A(f)ψ ∈ Dom(H^{1/2}), so lim I_η < ∞ by dominated convergence. **This bypasses the Córdova–Shao closure for finiteness.**

**Coordinator arithmetic check:**
* Item 3: \|a\|² = \|λ\|²\|K\|² with J ~ Δ ~ S (U ≪ S) gives S^{3/2}U^{3/2}·S^{−13/2}U^{1/2} = S^{−5}U² ✓.
* Item 4: the shell sum of C\|a₀\|² is ~S^{−3}U². Summing U ≤ S/L gives S^{−1}L^{−2}, and summing S ≳ L gives L^{−3} ✓; with S ≳ L², L^{−4} ✓.

**Not checked:** items 1, 2, 3 (as derivations), 5, and the O(N³) step. An independent audit has been launched (→ `ANGULAR_AUDIT.md`).

**Retained inputs:** physical matching (all-J ladder and first moment), contact matching (exact BMS coefficient), fixed-window planar regularity, Δ_φ > 1, and KKSZ (4.66) at the light-cone edge.

## 22. Audit of §21 (ANGULAR_AUDIT.md): confirmed with fixes

* **Proved, given standard CFT axioms and two finite-J extrapolations** (verified exactly for J ≤ 30): 𝒟(c_T) = 0 at each finite c_T, via reflection-positive two-scale crossing (mirror ordering) and the kernel bounds.
* **Corrections:**
  * B_N = O(1/N), not O(N⁻²) (§18 overstated). The "absolute angular swapping" sum diverges logarithmically; 𝒟 = 0 is unaffected.
  * The low-spin sectors need the collinear growth of F (DV).
  * The Casimir/K⁰ domain step (§21 item 5) is [D] and unnecessary.
* **To write out:** the symbolic-J collinear identities and the crossing lemma with the ordering stated.
* **Final input list of the angular theorem:**
  1. physical ladder and first-moment matching (LADDER_ALLJ);
  2. the contact coefficient (KZ algebra; Δ_φ > 1);
  3. fixed-window planar regularity;
  4. KKSZ light-cone convergence (4.66) and existence of the mixed commutator.

  **The regulated null-plane bound is no longer an assumption.**

## 23. Free-field test (FREE_NULLPLANE.md)

* **Null-plane algebra.** It is exact in the free scalar and Dirac fermion, with no extra contacts.
* **Correction to §17/§18's quadratic-form statement.** The limit is lim I_η = ½⟨A(f′²)⟩ = 2‖H^{1/2}Aψ‖² − 2Re⟨Aψ, AHψ⟩, i.e. 2‖H^{1/2}Aψ‖² = 2Re⟨Aψ, AHψ⟩ + ½⟨A(f′²)⟩. Consistent with §17.
* **Endpoint rule.** It holds exactly (3/2 + 3/2 = 3), with 𝒟 = 0. The O(N²) family bound is sharp (saturated).
* **Cutoffs.** C¹ cutoffs are needed for f_{m,R}.

## 24. Citation check: all-J helicity norms (Gillioz, arXiv:2112.12540, JHEP 03 (2022) 170)

**Verified from the paper text** (coordinator, 26 Nov):
* Eq. (2.35): ⟨⟨O_{j,−m}(−p)O_{j,m}(p)⟩⟩ = (p²)^{Δ−j−3/2}C_{Δ,j,m}.
* Eq. (2.37): the special conformal Ward identity gives the recursion (j+m+1)(Δ−m−2)C_{Δ,j,m+1} = (j−m)(Δ+m−1)C_{Δ,j,m}.
* Eq. (2.38): C_{Δ,j,m} = (2j)!/((j−m)!(j+m)!)·Γ(Δ+m−1)Γ(Δ−m−1)/(Γ(Δ+j−1)Γ(Δ−j−1))·C_{Δ,j}.
* Eq. (2.39): the equivalent Pochhammer form.
* Eq. (2.40): unitarity Δ ≥ j+1.

This is exactly §4's N_m/N_J and the ratio used in §8 and KERNELS_D3 §3.3, in d = 3 momentum-space helicity components.

**Effect.** The all-m helicity norm moves from [DV] (J ≤ 30) to **[Q]**, quoted with the derivation from the special conformal Ward identity. The all-J ladder normalisation is then [P given the operator identities] with no finite-J extrapolation.

**To do:** update the Supplement's angular section and the S1 legend ("all-m helicity norm [DV] for J ≤ 30" → "[Q], Ref. Gillioz eqs. (2.37)–(2.38)"), add the bib entry, and drop the "(checked for J ≤ 30)" qualifier. This waits until the jargon-cleanup agent releases supplement.tex.

## 25. Full-positivity attempt via the angular route (author-supplied, 26 Nov): flip cancellation, CDH resummation, and a published-weight endpoint obstruction

* **k = 3 flip cancellation (author).** On the physical ladders, Σ_λ H̄_{λ,−2}H_{λ,+2}sin(λ(θ₃−θ₁)) = 0. For even J this follows from the mirror identity of LADDER_ALLJ; for odd J the supports are disjoint. This closes ANGULAR_GO §1.2's open flip item for the physical ladder. Not re-checked.
* **CDH resummation (author).** Σ_n ((2d+2)_n/n!)² rⁿ Q_n^{(d)}(ν²) = P_{−1/2+iν}((1+r)/(1−r))/(1−r)^{2d+2}, with d = 0, 1. Proof: the recurrence becomes an ODE that reduces to the Legendre equation. **Verified numerically** (coordinator): d = 0, 1 at (r,ν) = (0.3,0.7), (0.5,2.0), (0.2,0), 12 digits, sum to n = 400.
* **Consequences (author).**
  * At fixed impact parameter, C/T → tanh²(β/2) and the k = 3 ladder sum tends to P_{−1/2+iν}(cosh β)/(1−r)⁴.
  * The k = 2 first-moment sum Z(r,x) gives 576(1−r)³Z = (16x²−728x+969)((1+r)/(1−r))L + (1896−352x)rL′.
  * Together with R₋P_{ℓ−1} + R₊P_{ℓ+1} = −4608(4x+1)/((4x+9)(4x+49))·(1−r)³Z, the shifted EE kernel of Ω^EE is recovered from the angular ladder.
* **OBSTRUCTION (author; verification in progress).**
  * *Finding.* The published k = 2 continuum weight at Δ_g = 1000 has E_eff(0) = 2Re[χR₋](−i) ∈ (−3.924105, −3.924103)·10⁻⁶ < 0.
  * *Consequence.* The realised leading-Regge image then ~ e^{−β/2}E_eff(0) at large β, below the target floor ~ e^{−β/2}/β. So the published functional's leading-Regge action is negative for β ≳ 3·10⁴ (spins ~ e^{15000} m). This is consistent with T0's realised corollary stopping at β ≲ 200.
  * *Effect on (iii).* (iii) is false for the published weight on those states; T0 for the target weights is unaffected.
  * *Proposed repair.* δx₀ ≈ −1.297683·10⁻⁷ sets E_eff(0) = 0; an adjustment of the first log coefficient matches E_eff′(0⁺) = 0.18978 to the target ψ′(0⁺) = 0.18918.
  * *Assigned.* The ε_W agent (→ EPSW_CERT.md) is verifying, repairing, recertifying T0 for the repaired weight at all β, and recomputing τ₀/K.
* **Exact finite-(Δ,J) form (author).** H_b/|a₀|² = Σ_k W_{b,k}(y_k^Ω + z_k^Ω d_{b,k}/6 + σ_b q_k^Ω), a rational function of T at fixed J. This is a route to certifying every mass at fixed spin by polynomial positivity. Physical matching for Ω^EE at finite (Δ,J) is still open.

## 26. Exact finite-mass k = 2 angular coefficients of Ω^EE (author-supplied, 26 Nov)

**Derivation (author).** In the x₀ = ∞ gauge, with δ = ½+iν, the celestial factor is W_δ = e^{δ/2}(e₁e₃)^{−δ/2}, and K_δ(φ) is its integral over the mean detector angle. Integrating I_j(c) = ∫(cos t − c)^j\|cos t − c\|^{−δ}dt by parts gives the exact insertions. Substituting them into K_B2A4 §4's momentum-space EE expression gives:
* 𝒵 = 12K_δ − 2 sin φ K_δ′;
* 𝒴 = (6 − (ν²+¼)/2)(1+cos φ)K_δ − 7 sin φ K_δ′.

No large-Δ, large-J or large-ν approximation is used.

**Coordinator check [P] (Fourier algebra).**
* *Z.* The Fourier coefficients of 𝒵 are 12k_n − (n−1)k_{n−1} + (n+1)k_{n+1}. This is §3's Z_m up to the overall sign convention (raw z < 0).
* *Y.* 𝒴_n = A k_n + (A/2 − 7(n−1)/2)k_{n−1} + (A/2 + 7(n+1)/2)k_{n+1}, with A = 6 − s/2 and s = ν² + ¼. Using the d = 0, ℓ = 0 recurrence (m+1)²k_{m+1} + (s − 2m²)k_m + (m−1)²k_{m−1} = 0 at m−1, m and m+1, this equals §3's five-term Y_m coefficient by coefficient (k_m: 6−m²; k_{m∓1}: (m∓1)²/2 − s/4 + (13∓7m)/2).

**Significance.** The angular-route coefficients are exactly the k = 2 kernels of the original Ω^EE at finite mass. This is the physical-matching step for the full-positivity route, given the I_j identities and the continuation from −1 < Re δ < 0.

## 27. Finite-mass positivity programme for Ω^EE (author-supplied, 26 Nov)

| claim | check |
|---|---|
| Finite-mass k = 2 action ∝ Σ_k W_{b,k}(y_{k+2} + d_{b,k}z_{k+2}/6). It is diagonal in the ladder basis: no b = 0/4 cross terms (disjoint supports, A tridiagonal), plus the k = 3 flip cancellation | consistent with §§3, 25, 26. Absolute normalisation vs T0's N_s still open |
| Numerical cross-check: (H_ang − H_hyb)/N_s ≤ 0.0026 on 39 states (13 spin/channel pairs, twists 1000/1300/2000; half-step quadrature) | not re-run [NE]; uses the old weight and endpoint normalisation |
| d_{b,k}/3 = T(1/(k+3) − r/(k+1)) + V_{b,k}, with V_{0,k} = 4k/(k+3) and V_{4,k} = 4 + 12/(k+1) | **verified** by hand (u_k/(k+1) − u_{k−1}/(k+3)) |
| Every mass at fixed J reduces to a polynomial 𝒫_{J,b}(w) > 0 on 0 ≤ w ≤ Λ²/((J+998)(J+999)). Bernstein/interval certificate passes for J ≤ 13 (moment uncertainty 1e-8), J ≤ 15 (1e-10), J ≤ 16 (1e-12) | polynomial tests rigorous per author; **moments not certified** (old weight) |
| 0 ≤ r^k − U_{b,k} ≤ r^{k−1}S_{b,k}/T, giving \|H_ang − H_LR^realised\| ≤ 𝒦_b(r₀)/T uniformly in spin for r ≤ r₀ | proof elementary and correct (ratio ≤ r; telescoping); constant not evaluated |
| Large-r limit (repaired weight with χ_re(i) = 0 and the target slope): Z_k ~ z_*k^{−5/2}/log k, Y_k/Z_k → 1/4, Q_k = O(k^{−7/2}); ladder U → e^{−γ atanh(x/γ)} | [D], contour argument not re-checked |
| Step C: ∫₀^γ x^{7/2}U_γ[1/x − 2/x² + 9/(2γ²)]dx = ½∫₀^γ x^{3/2}U_γ dx | **verified numerically** (γ = 3, 10, 50, 1996; 12 digits) |
| 𝔑(γ) = Γ(5/2)⁻¹∫x^{3/2}U_γ, with 1 − 105/(8γ²) ≤ 𝔑 ≤ 1; at γ ≥ 1996, 1 − 𝔑 < 3.3e-6 | **verified**: 𝔑(1996) = 0.999996705605 against a lower bound of 0.999996705586 (tight); 𝔑(50) = 0.99480 ≥ 0.99475 |

**Programme status.**
* *Covered by methods:* fixed J at every mass (polynomial certificate); bounded impact parameter uniformly in spin (comparison lemma); unbounded impact parameter (positive gamma integral).
* *Open:*
  * certified angular moments for the repaired weight;
  * explicit thresholds for the two asymptotic regions;
  * certification of the finite transition region;
  * the absolute finite-mass normalisation against T0.

## 28. Absolute normalisation of the angular action against T0 (author-supplied, 26 Nov)

| claim | check |
|---|---|
| y₀ = 6f₀ + 13f₁ + 2f₂ with f₀ = −8/(4x+1), f₁ = 1, f₂ = (7−4x)/16 gives y₀ = −(16x²−440x+273)/(8(4x+1)); with F(x) = −9216(4x+1)/((4x+9)(4x+49)), F·y₀ = −36ℰ(ν), where ℰ is K_B2A4's exact Einstein light function | **verified exactly** (sympy) |
| The k = 2 physical prefactor is 𝒜₂ = 16πG/36 = 4πG/9; y₄/y₀ = 16Θ(ν) (K_B2A4), consistent with \|c\|² = (t₄²+α_T²)/16 | first part follows; y₄/y₀ not re-checked |
| k = 3 prefactor: 𝒜₃ = 4πG/9 = 𝒜₂, from the published ½Λ⁶∫ψ₃Ĉ₃ and \|c\|² = 3600\|ĝ₃\|² | not re-checked |
| Ω^EE_O = ℬ(Δ,J)Σ_b\|a_{b,0}\|²H_ang,b, with ℬ = (2πGΛ⁴/9)T⁴/(T−C)³ > 0 and H_ang,b = ((1−r)³/T)Σ_k W_{b,k}[2Y_k + (d_{b,k}/3)Z_k + σ_bΛ²Q_k] | **algebra verified** (sympy); same H_ang as §27 |
| OPE basis: 𝔠_h ∝ K_h(Δ,J)λ_h, with the exact kernels of ANGULAR_PROOFS. Ω^EE_O = H_ang,0(\|𝔠₊₊\|² + \|𝔠_odd\|²) + H_ang,4\|𝔠₊₋\|². Double-trace zeros are handled automatically | consistent |
| **Congruence lift of T0.** Define Ω_hyb,O = ℬΣ\|a\|²H_hyb,b and N_s,O^phys = ℬΣ\|a\|²N_s,b. Then the scalar per-channel estimate \|H_ang,b − H_hyb,b\| ≤ (0.1999−ε)N_s,b implies Ω^EE_O ≥ εN_s,O^phys ≥ 0; the normalisation cancels | **logic correct**. What remains to verify is that H_ang,b and T0's H_hyb,b are the same scalar function in the leading-Regge limit, with no relative factor (the matching the FULLPOS agent checks numerically) |
| Impact-parameter correction: r_β = tanh²(β/2) = (C+¼)/(T+¼) ≠ r; cosh β = (1+r)/(1−r) + 1/(2m²); ℬ_β/ℬ = (1 + 1/(4T))⁴, with excess < 1.004e-6 at gap 1000 | consistent (arithmetic) |

**Status of item 4.** Resolved by the congruence argument, given the retained physical identities, once H_ang,b → H_hyb,b with unit relative normalisation is confirmed (FULLPOS_CERT task 2).

## 29. No-stray-factor check and exact finite-T cancellation (author-supplied, 26 Nov)

| claim | check |
|---|---|
| Leading-Regge limit at fixed r: H_ang^{(2)} → 2∫ψ₂P_ℓ(η)dp (factor 2 = 9216/4608) and H_ang^{(3)}/(σ_b u) → ∫ψ₃P_ℓ(η)dp. T0's B₂ + cv₂ and (B₃ + cv₃)/(σ_b u) have the same limits (Mehler–Heine D_b → J₀; cv terms are part of the action). **No stray factor**, with unit relative k = 2/3 normalisation | consistent with §§25, 28; not independently re-derived. FULLPOS agent checking numerically |
| Repaired weight: the limit statement is H_ang,b − H_hyb,b^repaired → 0. Against the ideal target the residual is the certified de-shift error | correct distinction |
| **Summation by parts**, with L_k = (Z_k − (k/(k+3))Z_{k−1})/(k+3): H_ang^{(2)} = w³Σ D_k²U_{b,k}[L_k + (2Y_k − u_{k−1}L_k)/T] exactly | **verified** in exact rationals on 40 random (b, J, Δ) with random moments |
| r^k − U_{b,k} = wΣ_{j≤min(k−1,J−b)} u_jU_{b,j}r^{k−1−j}/(T−u_j), including beyond termination | **verified**, same test |
| Signed finite remainder: H_ang − H_LR = (w³/T)Σ_j U_{b,j}[M_{b,j} − (wu_j/(1−u_j/T))(F_j + σ_bΛ²G_j/T)], with tail generating functions F_j, G_j from the CDH resummation | follows from the two identities |
| Improved ladder bound: for r ≥ ½, 0 ≤ r^k − U_{b,k} ≤ (2(1−r)/C)S_{b,k}r^k; Cw² ≥ 4a²r², Tw² ≥ 4a²r | **verified** on the random heavy-ish cases with r ≥ ½ (no failures). The heavy-state inequalities are not re-checked |
| Tail certificate: if, for k ≥ K, \|P_k − P_k^*\| ≤ P_k^*/4 and \|M_{b,k}\| ≤ (k+b+1)²P_k^*/4 (P_k^* = ((5/2)_k/k!)ℓ_k, ℓ_k nonincreasing), plus three prefix conditions, then the relative k = 2 error is < 80.625/a² ≤ 8.2e-5 at gap 1000 | conditional; logic plausible (negative-binomial moments) |
| Repaired weight needs a two-scale reference: ℓ_k = A + B/log(k+2), A = E_eff(0)/Λ ≈ 9.7e-17, B = 2E_eff′(0⁺)/(πΛ) ≈ 1.21e-4; 𝒮(r) ≥ w^{−5/2}[A + B/log(2 + 5r/(2w))] | candidate; to be certified |
| Full comparison: H_ang − H_hyb^rep = [H_ang − H_LR(r)] + [H_LR(r) − H_LR(r_β)] − [H_flat − H_Bessel]. The last two terms (the Weyl-shift and finite-spin flat corrections) must be kept in the transition region | noted |

## 30. Direct certificate of the finite middle region through J = 26 (author-supplied, 26 Nov; awaiting reproduction)

**Claim.** For every physical channel with J ≤ 26 and every mass (u = Λ²/m² ∈ [0, 999/(999+2J)]) at gap 1000, with the repaired weight:
* \|H_ang,b − H_hyb,b^rep\| ≤ 0.01·H_hyb,b^rep, hence H_ang,b ≥ 0.99·0.1999·N_s,b = 0.1979·N_s,b, keeping the exact Weyl shift and the finite-spin flat polynomials;
* the certificate consists of 74 Bernstein-certified polynomial inequalities after multiplying by 𝒟_{J,b}(u) = (1+Cu/Λ²)⁴∏(1+(C−u_j)u/Λ²);
* curvature integral enclosed by a conical Legendre expansion to degree 110 (truncation < 1.8e-45); de-shift residual < 5.38e-9;
* the cutoff is handled exactly by the Gaussian-convolution recurrence M_{n+1} = cM_n + ns²M_{n−1} + s²∫ν^nG_c(1−θ)′, with an erfc initial value.

**Bug report.** The d = 1 CDH recurrence must divide by (n+4)², not (n+3)², so that Q₁ = (15−4x)/64.
* *Coordinator check:* the committed snapshot 8d6af34 of `code/fullpos/fp_common.py` had (n+3)²; the agent's working copy already has (n+4)² in both files.
* *Follow-up:* the agent is asked to confirm every cache was regenerated after the fix.

**Status.** The author's computation is outside the repository. FULLPOS agent asked to reproduce it independently, extend in J, and find the tail threshold.

## 31. Exact elimination of Y_k (author-supplied, 26 Nov)

| claim | check |
|---|---|
| 2Y_k = ((k+3)²/2)L_k + ((k+4)²/2)L_{k+1}, for an arbitrary harmonic sequence f_m (any smearing weight) | **verified symbolically** with free symbols f₀…f₁₅, k = 1…5 |
| M_{b,k} = D_k²(2Y_k − u_{k−1}L_k) = [(5−2b)k + 5 − b(b−1)]P_k + ((k+1)²/2)(P_{k+1} − P_k); M_{0,k} = 5(k+1)P_k + …, M_{4,k} = −(3k+7)P_k + … | **verified symbolically**, b = 0, 4, k = 1…5 |

**Consequence.** The finite-T correction is a fixed linear combination of P_k and ΔP_k. Only the P_k (the differenced Z moments) need envelopes; no separate Y_k asymptotics are needed. The §29 tail certificate reduces to one coefficient sequence.

## 32. Single-sequence positive representation, forward-strip certificate, cutoff formula (author-supplied, 26 Nov)

| claim | check |
|---|---|
| H_ang,b^{(2)} = (1−r)³Σ_{k=0}^{J−b+1}V_{b,k}P_k, with V_{b,k} = U_{b,k}[1 + ((k+3)²/2 − u_{k−1})/T] + (k²/2T)U_{b,k−1} ≥ 0 | **verified exactly** (40 random exact-rational cases, random harmonic sequences; identity and V ≥ 0) |
| Σ P_k z^k = (2/(Λ(1−z)³))∫E_eff,θ(ν)P_{−1/2+iν}((1+z)/(1−z))dν | consistent with §§25, 29; not re-checked |
| Forward strip b̂ ≤ 1, J ≥ 27, Δ−J ≥ 1000: Cauchy majorants \|P_k\| ≤ 200s^{−k}, \|D_k²Q_k\| ≤ 120s^{−k} (s = Λ⁻²); tight P₀, P₁, P₂; total error < 0.06N_s, giving H_ang,b > 0.1399N_s,b | not re-checked; needs reproduction |
| Tail certificate simplified: \|P_k − P_k^*\| ≤ P_k^*/4 and \|D_k²Q_k\| ≤ (k+1)P_k^* for k ≥ K, plus four prefix constants D₀…D₃; then the relative k = 2 error is < 106/a² < 1.1e-4, and the k = 3 error is ≤ 7Λ²(1−r)/(4a²r) | conditional |
| Asymptotics (repaired weight, A > 0): P_k = ((5/2)_k/k!)[A + B/log(k+2) + O(log log k/log²k)]; D_k²Q_k = ((7/2)_k/k!)[B₃/log(k+2) + …]. So the envelopes hold beyond some finite K | [D]; effective K open |
| Extreme tail: \|H_flat − H_Bessel\| ≤ 10¹³/b̂² (Jacobi ODE, two integrations by parts) and Weyl ≤ 3·10⁶e^{−3β/2} (β ≥ 100); relative < 1e-8 for β ≥ 100 | not re-checked |
| Cutoff: w_* = min{sech²(50), (A/4D₀)^{2/5}, (A/4D₁)^{2/11}, (A/4D₂)^{2/9}, (A/D₃)^{2/7}} and β_* = 2 arcosh(w_*^{−1/2}); for 1 − r ≤ w_* the tail certificate applies | formula given; needs certified D_i and K |

**Remaining.**
* an effective harmonic cutoff K, from certified P_k, Q_k;
* evaluating w_* and β_*;
* certifying the finite middle region (J ≥ 27, b̂ > 1, β < β_*), with error O(1/T) by the Cauchy bounds;
* reproducing the J ≤ 26 and forward-strip certificates.

## 33. Effective P_k remainder (author-supplied, 26 Nov)

**Claim** (pen-and-paper, conservative constants):
* \|k!P_k/(5/2)_k − A − B/log(k+2)\| ≤ 0.02/log²(k+2) for k ≥ ⌈e^{100000}⌉.
* \|D_k²Q_k\| ≤ (k+1)P_k^* for k ≥ ⌈exp(exp(900))⌉.
* **Method.** Effective Darboux transfer with remainder e^{10000}k^{−1/2}, then a contour shift to Im ν = 1/8. On that segment only the logarithmic basis contributes to the imaginary part (exact formula), with \|W_α′\| < 6.
* **Result.** The explicit remainder is (2/Λ)[6/L² + \|f₁\|e^{−L/8}/(πL) + e^{8100−L/8} + e^{10000−L/2}], with L = log k. Not re-derived.

**Coordinator assessment.**
* *The cutoff is intrinsic.* The weight lives at ν ≲ Λ ≈ 10³. The coefficient asymptotics are an expansion in e^{iν log k}, so the 1/log k regime only sets in for log k ≳ Λ, i.e. k ≳ e^{~1000} or beyond. Coefficient-wise envelopes cannot be checked computationally up to K.
* *Proposed pivot.* For the middle region (J ≥ 27, b̂ > 1), work in the impact-parameter variable instead of per harmonic. §29's exact signed remainder writes H_ang − H_LR through the finite ladder weights U_{b,j} and the tail generating functions F_j(r), G_j(r). Those are leading-Regge-type integrals at fixed r (F_j = (w⁻³H_LR^{(2)}(r) − Σ_{k≤j}P_kr^k)/r^{j+1}), which T0's machinery already certifies as functions of r. Combined with U_{b,j} ≤ e^{−j/R} and the improved ladder bound (§29), this avoids enumerating harmonics.

## 34. Bivariate (x, v) certificate for J ≥ 240 and the collective Mellin transform (author-supplied, 26 Nov)

**Context.** The FULLPOS_CERT checkpoint certifies H_ang,b > 0 at every mass for every channel with J ≤ 239 (gap 1000) and J ≤ 285 (gap 300) [certified]. Its finding: separate \|F_j(r)\| tails are useless (majorant ~10⁶² against an action of ~10⁻⁴ at J = 100).

| claim | check |
|---|---|
| With x = b̂²/4, v = (J+½)⁻², λ = Λ²: λ^kU_{b,k} = x^kA_{b,k}(v)/∏_{j<k}[1 + (x/λ)(1−(b+j+½)²v)], where A_{b,k} = ∏_{j<k}[1−(b+j+½)²v]; r_β = x/(λ+x) | **verified exactly** (30 random states; uses m² = T − C) |
| T0's flat partial wave D_b(1−2up²) = Σ_k(−xp²)^kA_{b,k}(v)/(k!)², with u = xv | b = 0 follows from P_J(1−2y) = ₂F₁(−J,J+1;1;y) with (−J)_k(J+1)_k = (−1)^k(J+½)^{2k}A_{0,k}; b = 4 not re-checked |
| J ≥ 240 ⇔ 0 ≤ v ≤ 4/481². A bivariate certificate on (x,v) ∈ [0,121]×[0,4/481²] (b̂ ≤ 22) passes for both channels: \|H_ang − H_hyb^rep\| ≤ 0.05H_hyb^rep (1% for b̂ ≤ 16), so H_ang ≥ 0.1899N_s. Degree-(66,66) Bernstein, 240-digit outward rounding, ideal-density moments via exact NF02 rationals; repaired-weight residual < 0.002537 within a 0.003 allowance | not reproduced (author computation) |
| Beta-mixture identity: Σ_k U_{b,k}p_k = E_B[[z^n](1−z)^{−a−1}F(Bz)]/C(n+a,n), with B ~ Beta(n+2b+1, a), n = J−b, a = Δ−J−2 | plausible (beta moments); not re-checked |
| Mellin multiplier 𝖴_b(s) (Γ-ratio interpolation of U_{b,k}); \|𝖴₀(iξ)\| = 1; \|𝖴_b(iξ) − r_β^{iξ}\| ≤ M_b[(b²+1)\|ξ\| + bξ² + \|ξ\|³/3](1/n² − 1/N²) | not re-checked; boundary term at r = 1 is unresolved (author) |

**Coverage at gap 1000:**
* J ≤ 239: every mass [certified];
* J ≥ 240 with b̂ ≤ 22: author's bivariate certificate (to reproduce);
* J ≥ 240 with b̂ > 22: open (the Mellin route has an unresolved boundary term);
* the extreme large-β limit: analytic (§27).

## 35. Mellin boundary correction, uniform singular-power theorem, and explicit tails (author-supplied, 26 Nov)

| claim | check |
|---|---|
| Correction to §34: the rational Mellin interpolant vanishes only at k = n+1…N, so ℒ(F) = 𝖴_b(𝒟)F(1) − Σ_{k>N}𝖴_b(k)p_k. The extra tail cannot be dropped | correct (the earlier Mellin argument was incomplete) |
| **Uniform singular-power theorem.** For 0 ≤ α ≤ 4, r ≥ ½ and every heavy state at gap 1000: \|𝒱_b((1−z)^{−α})w^α − 1\| ≤ 94/a² < 9.44e-5. It extends to every nonnegative mixture of such powers | **tested numerically**: 40 states with r ≥ ½ (a = 998–3000, J = 3a–30a, α ∈ {0.75, 1.5, 2.5, 4}); the worst \|ratio−1\|·a² is 14.8 ≤ 94 |
| Positive part S₂ of the repaired generating function. With the contour at height ¾, W(y) = −2B(iy)Im f(iy) ≥ 0 on [0, ¾] (sign reversals of B and R₊ at y = ½ cancel). The resulting mixture of powers w^{−5/2+y+j} has nonnegative weights, and its floor is w³S₂ ≥ A√w | sign structure not re-checked |
| Extreme tail 1−r ≤ e^{−5000}: \|H_ang − H_hyb^rep\| < 1e-4·N_s, with no harmonic cutoff (coefficient bounds e^{2600}c_{7/4,k} + e^{7100}) | not re-checked; conservative constants |
| Large T: \|H_ang − H_hyb^rep\| < 0.01N_s for T ≥ e^{30000} and w ≥ e^{−5000}; the flat correction is \|D_b − J₀\| ≤ 5Ju + 2057√u | not re-checked |
| **Remaining region:** 240 ≤ J < e^{15000}, T < e^{30000}, b̂ > 22, 1−r > e^{−5000}. Bounded but astronomically large; exhaustive certification is infeasible (author) | agreed |

**Assessment.**
* *Unbounded directions:* all three now have explicit (if crude) certificates.
* *What is left:* a bounded but astronomically large transition region. The constants e^{~7000} come from the e^{2πν} growth of the Legendre functions against a weight supported up to ν ≈ Λ. Closing the region needs a collective remainder estimate with constants polynomial in Λ, not exponential. That is a research problem beyond this submission's window.

## 36. The complex-power route: exact beta reduction, obstructions, and the remaining signed estimate (author-supplied, 26 Nov)

| claim | check |
|---|---|
| No uniform polynomial bound in \|ξ\|: 𝒱_b((1−z)^{−σ−iξ}) has degree exactly n+1 in α (V_{b,n+1} > 0) | correct (elementary) |
| Exact first-moment reduction: 𝒱_b((1−z)^{−α}) = (1 + B_b/T)S_b(α) + (αA_b/T)[S_b(α+1) − S_b(α)], with A_b = (9+α)/2 − 2b, B_b = (9+α)/2 − b(b−1), S_b(α) = Σ_k U_{b,k}(α)_k/k! | **verified numerically** (30 random states, complex α; rel. err. 6e-41) |
| Beta-integral form of S_b(α) for Re α > 0 | not re-checked |
| Fixed-twist limit 𝔑_a(α) = Γ(α)⁻¹∫₀^{2a}x^{α−1}e^{−2a atanh(x/2a)}dx; turning point at \|ξ\| = a, with exponential growth e^{a g(ξ/a)} above it | not re-checked; consistent with §27 at real α |
| Restricted band: \|𝔑_a(5/2+iξ) − 1\| ≤ 10⁸(1+\|ξ\|)⁴/a² for \|ξ\| ≤ a+60 | not re-checked |
| **Obstruction:** absolute-value frequency integration over 5 ≤ ν ≤ 10 alone costs > 7e-8·√w at a = 998, against an available margin of < 9.5e-9·√w at states with 1−r ~ e^{−4000}. So the absolute route cannot certify the middle region, even with polynomial constants | logic sound (lower bounds from a Taylor-remainder argument); numbers not re-checked |

**Remaining target (author).** A signed estimate of ℰ₂ = (4/Λ)√w Re∫f(ν)B(ν)e^{iν log(1/w)}[w^α𝒱_b(F_α) − 1]dν (α = 5/2+iν), keeping the oscillatory factor, with the origin residue isolated, by contour deformation or integration by parts. **Status: full AdS positivity remains open; this is the precise missing estimate.**

## 37. Finite-state complex beta-kernel bound (author-supplied, 26 Nov; part 1 of a signed estimate)

| claim | check |
|---|---|
| Euler representation K_{n,a,b}(q) = ₂F₁(−n,a;n+a+2b+1;q) = B(a,n+2b+1)⁻¹∫₀^∞x^{a−1}(1+(1−q)x)^n(1+x)^{−2n−a−2b−1}dx | **verified** (n = 50, complex q, 12 digits vs mpmath hyp2f1) |
| Global R(t) = \|1+ct\|/\|1+dt\|² ≤ (1+t/4)⁻¹ (d = 1/(1+q), c = 2d−1, q in the first-quadrant unit disk) via the polynomial (1+pt+st²)² − (1+t/4)²[1+2(p−1)t+(4s−2p+1)t²] having coefficients ≥ 0 on 1 ≤ p ≤ 2, ¼ ≤ s ≤ p/2 | **coefficients computed**: 0, 3/2, (16p²+16p−32s−1)/16, (16ps+7p−16s−3)/8, (2p+16s²−4s−1)/16. Nonnegativity checked on 20,000 samples (no failures) |
| **\|K_{n,a,b}(q)\| ≤ 2\|1+q\|^{−a}** for n ≥ 1024a², a ≥ 998, q in the closed first-quadrant unit disk | proof structure sound; **numerically** max \|K\|·\|1+q\|^a = 1.002 at a = 2, 3, 5 with n = 1024a² (the theorem is stated for a ≥ 998; the small-a test is a sanity check of the mechanism) |

**Status.** This is the finite-state replacement for K ~ (1+q)^{−a}. The author reports an explicit relative error of 1.01e-4 for the leading connection term ℰ₂ in an enlarged large-impact region (details to follow). The remaining connection-function terms and the lower-impact transition region are still open.

## 38. Signed estimate of the leading connection term (author-supplied, 26 Nov)

**Claims (not re-checked; constants conservative).**
* *Second beta-kernel regime.* For n < 1024a² and L = log(1/(1−r)) ≥ 40, a beta-expectation representation K = (1+q)^{−a}E[((1−2eS)/(1−eS)²)^n(1−eS)^{−a−2b−1}] gives \|K\| ≤ 2\|1+q\|^{−a} on a smaller quarter disk. Together with §37 this covers every heavy state with L ≥ 40.
* *Polynomial transform bound.* \|D(t+9i/4)\| ≤ 200(1+t)^{5/4}e^{a g(t/a)} < 5·10¹² on the band, a finite-state polynomial bound, where D = w^α𝒱_b(F_α) − 1.
* *Contour.* The frequency contour runs 0 → ¾i → ¼+¾i → ¼+(9/4)i → ∞+(9/4)i. The local piece uses W ≥ 0 and \|D(iy)\| ≤ 100/a², giving \|ℰ₂,local\| ≤ (100/a²)𝓜(L), with 𝓜(L) = √w[A + (2/Λ)∫₀^{3/4}W e^{−Ly}].
* *Detour.* It costs ≤ (e^{−3L/4}/a²)(1.5·10¹⁰ + 1.2·10¹³/(L−17)) · 4√w/Λ.
* *Horizontal leg.* It costs ≤ 4·10³¹e^{−9L/4} · 4√w/Λ.
* *Coefficient sums.* The Gaussian coefficient sums were certified from the repaired JSON (906 and 1.22·10⁶ for the two widths).
* *Result.* \|ℰ₂\| < 0.001𝓜(L) for L ≥ 40, and < 1.01e-4𝓜 for L ≥ 64. With 𝓜 < 10N_s^rep, \|ℰ₂\| < 0.00101N_s^rep for every spin and mass with 1−r ≤ e^{−64}.

**Coordinator note.** The large-impact threshold improves from 1−r ≤ e^{−5000} (§35) to e^{−64} for the leading connection term.

**Open.**
1. The leading term for L < 64.
2. The remaining hypergeometric connection terms: j ≥ 1 of ₂F₁(a,a;2a;w) and the second branch.
3. The bounded transition region J ≥ 240, b̂ > 22, 1−r > e^{−64}.
4. Reproduction of §§30, 32, 34, 37 and 38.

## 39. Coordinate correction (author, 26 Nov)

L = −log(1−r) = log(1 + (b̂²/4Λ²)[1 − (2J+1)⁻²]) ≥ 0. At gap 1000 the rectangle edge b̂ = 282.84 is **L ≈ 0.0199**; L = 10 ↔ b̂ ≈ 2.96·10⁵ (coordinator recomputed). The coordinator's "L ≈ −4" in the §38 summary was wrong.

**The open band** is therefore L ∈ (≈0.02, L_signed), where L_signed is the signed estimate's threshold (64 now; CONTOUR_BOUNDS pending). Tightening the contour threshold alone cannot close it: a separate moderate-impact estimate is required.

## 40. Full-tail positivity for L ≥ 28, including all connection terms and k = 3 (author-supplied, 26 Nov; to reproduce)

**Claims.**
* *Off-axis bounds.* With the cancellation between Gaussian families kept: \|f(t+¾i)\| < 4.1 (old: 10⁷), and a vertical-segment table rising from 13 to 2.7·10¹³; the horizontal integral is < 10²¹ (old: 4·10³¹). The CONTOUR_BOUNDS agent is computing the same independently.
* *Connection terms, j < 600.* Integration by parts in the beta formula, S_b(α−j) = ((−1)^jΓ(N+1+α−j)/(Γ(N+1)Γ(α)))∫q^{α−1}∂_q^j[(1−q)^aK(q)]dq, plus the differentiated-kernel bound \|∂^j[(1−q)^aK]\| ≤ 2(2a+j)^j\|1−q\|^{a−j}\|1+q\|^{−a}, give a collective bound on terms j = 1…599.
* *Connection terms, j ≥ 600.* A localisation theorem: if \|F(z)\| ≤ K\|1−z\|^m (1 ≤ m ≤ 597), then \|ℒ(F)\| ≤ 50K(4w)^m. So the e^{9400}-size remainder is multiplied by (4w)^{595} < e^{−6400·…} at L ≥ 28.
* *k = 3.* Fully included: \|H_ang^{(3)}\| < 4·10⁵ w√w < 0.529·𝓜(L).
* **Result.** H_ang,b > 0.45·𝓜(L) > 0 for L ≥ 28, i.e. b̂ ≳ 2.4·10⁹, every spin and mass, both channels. The k = 2 error budget at L = 28 is < 0.0035·𝓜.
* *Centred-Taylor bridge.* For the middle region, expand about z_c = ρμ, ρ = n/(n+a), μ = (n+2b+1)/(n+a+2b+1). Monomials transform exactly: ℒ((z−z_c)^ℓ) = Σ_k C(ℓ,k)(−z_c)^{ℓ−k}U_{b,k} (linearity). The remainder is bounded through centred beta moments. A uniform bound M on the Taylor remainder over the region is still open.

**Coverage at gap 1000.**
* *Certified:*
  * J ≤ 239, every mass;
  * J ≥ 240 with b̂ ≤ 100;
  * J ≥ 420 with b̂ ≤ 282.84;
  * J ≥ 651 with b̂ ≤ 400 (L ≈ 0.0393).
* *Derived (to reproduce):* L ≥ 28 (b̂ ≳ 2.4·10⁹).
* ***Open:*** **L ∈ (0.039, 28)**, i.e. 400 ≲ b̂ ≲ 2.4·10⁹, at large spin.

## 41. Local Taylor remainder constants, polynomial in Λ (author-supplied, 26 Nov)

* **Expansion centre.** c = z_c = ρμ = (C − b(b+1))/(T − b(b+1)) = U_{b,1}. **Verified exactly** on 200 random states. So c = r for b = 0, and r − c = 20(1−r)/(T−20) for b = 4. In the band 0.039 < L < 28: c > 0.038 and d = 1 − c > e^{−28}.
* **Local constants.** On \|z − c\| ≤ R_c/2, with R_c = √c·d/1200:
  * M_{2,m} = (116/d³)P_{−1/2}((1+c)/(1−c))(1200/(√c d))^m;
  * M_{3,m} = (68/d⁴)P_{−1/2}((1+c)/(1−c))(1200/(√c d))^m.

  These are of order Λ^m, not e^{O(Λ)}. Not re-checked.
* **Qualification (author).** The bound is local. The beta-contour remainder estimate also needs control of F outside \|z − c\| ≤ R_c/2, and analyticity alone gives no polynomial bound there. **Open.**

## 42. Explicit local remainder, first-moment operator, global obstruction, weighted outer bound (author-supplied, 26 Nov)

* **Inputs (author).** The exponential moments ∫|f|e^{ν/1100} < 28Λ and ∫|ψ₃|e^{Λp/1100} < 33 come from NF02 norms, EPSW_CERT and the Gaussian coefficients. They give sup over |z−c| ≤ R_c of |F₂| ≤ S₂ = 58 d⁻³P_{−1/2}, and of |F₃| ≤ S₃ = 34 d⁻⁴P_{−1/2}. Then M_{i,m} = 2S_i R_c^{−m}.
* **First-moment operator.** 𝒱_b(F) = ℒ(F + 𝒟_bF/T), with 𝒟_b = z(1−z)/2 ∂² + [(7/2−2b)z + 1/2]∂ + [9/2 − b(b−1)]. The constant M̃₂ is valid for \|z−c\| ≤ R_c/4. Since ℒ(z−c) = 0 (because c = U_{b,1}), the bridge starts at the second centred moment.
* **Checked here:**
  * \|QR₋(t−i)\| ≤ 16, exact: it reduces to 0 ≤ 180t² + 621, and the sup is 16 as t → ∞.
  * The β-disk bound: R_c/(√(c−R_c)(d−R_c)) · 1100 ≤ 0.919 < 1 for all c in [0.038, 1).
  * The weighted-exponent maximum N·arctan(N/A) − (A/2)log(1 + N²/A²) with A = a + ½, and its bound N²/(2A). The bound follows from g(x) = x·arctan x − ½log(1+x²) with g′ = arctan x ≤ x.
  * The table matches to all quoted digits: 475.22, 54.50, 0.547, 5.47e-4.
* **Global obstruction (author).** At J = a = 998, b = 0, the contour samples z = −0.49, where \|F₃\| > e^{1137}. Any unweighted global M must exceed e^{1136}. So absolute-value remainders on the whole contour cannot work at low twist.
* **Weighted outer bound (author).** K_a\|F_i\| ≤ A_i δ^{−p_i} exp[N_i·arctan s − ((a+½)/2)·log(1+s²)], with (A,p,N) = (22,3,1045) for F₂ and (13,4,999) for F₃. The peak sits at s = N/(a+½), far from the local disk.
* **Status (coordinator).** Combined local + outer inequality **not closed**:
  * The outer exponent is ≤ 1 only for a ≳ 5.5×10⁵, ≤ 5 for a ≳ 1.1×10⁵, and ≤ 20 for a ≳ 2.7×10⁴.
  * At a ≈ 998 to 10⁴ it is e^{475} to e^{54}. This must be beaten by cancellation, consistent with the obstruction above.
  * The low-twist strip 998 ≤ a ≲ 10⁵, inside 0.039 < L < 28, needs a signed or optimised-contour treatment, or direct certification.

## 43. Contour bounds certified (CONTOUR agent; see CONTOUR_BOUNDS.md, code/contour/)

* **All five §38 inputs certified with Arb.** Each is far sharper than §38 assumed:
  * sup\|f(t+¾i)\| ≤ 3.924, against 10⁷ in §38;
  * the vertical leg weighted by e^{−17(y−¾)} ≤ 227, against 2·10⁹;
  * ∫\|f(t+9i/4)\| ≤ 9.01·10¹³ and ∫\|fB\| ≤ 5.30·10¹¹, against 3·10¹⁷;
  * W = (½−y)²V(y) with V > 0, so W ≥ 0 exactly, and ∫W e^{−Ly} ≥ 0.0545/L at L = 28;
  * \|B\| ≤ 0.104 on the detour and ≤ 0.017 on the leg, against 20.
* **§40 cross-checks.** \|f(t+¾i)\| < 4.1 agrees. All twelve entries of the ⅛-interval table agree.
* **Thresholds for ℰ₂, by the transform bound D used:**
  * **Uniform \|D\| ≤ 5·10¹²** (the §38 number): L ≥ 30.16 for 1%, L ≥ 29.12 for 10%. ρ(28) ≤ 1.19, so **L ≥ 28 is NOT supported on this bound**.
  * **D-profile 200(1+t)^{5/4}e^{a g(t/a)}** (§38 line 636, not re-checked): L ≥ 20.27 for 1%, 19.22 for 10%. §40's ∫\|fBD\| < 10²¹ rests on this profile; the uniform bound gives 2.65·10²⁴.
* **Bottleneck.** About 99% of the bound comes from the horizontal leg at Y = 9/4. The transform bound D is now the limiting input, not f.
* **Open.** Prove the complex-exponent D-profile. That is the complex-α version of the §35 singular-power theorem, D = w^α𝒱_b(F_α) − 1 on α = t + 9i/4. Also still open: the other connection terms and k = 3 at the new threshold.

## 44. FULLPOS results (agent; FULLPOS_CERT.md, FULLPOS_STRIPSCAN.md, code/fullpos/)

* **Certified at gap 1000, H_ang,b > 0 in both channels:**
  * J ≤ 766 at every mass: one Bernstein certificate per J, 1147 polynomials, no failures;
  * J ≥ 240 with b̂ ≤ 100; J ≥ 420 with b̂ ≤ 282.8; J ≥ 651 with b̂ ≤ 400, all via the §34 (x, v) certificate with exact denominators.
* **Gap 300.** J ≤ 285 is certified at every mass.
* **Open.** J ≥ 767 with 400 < b̂ ≤ 44.7√J, i.e. m_gap ≤ m < 5(J+½), including r → 1.
* **Moments.**
  * Y_k and Z_k are certified for k ≤ 1000. There is a certified plateau \|Z_k\| ≈ 10^{672.5} ≈ e^{(π/2)Λ}.
  * The §§29–33 tail envelopes fail for every k ≤ 1000, with an intrinsic threshold. The §29 pivot's triangle bound overshoots by 10^{22} to 10^{127}.
  * §32 is exact only with u₋₁ = (b−1)b.
* **Normalisation.** H_ang → T0's limit with factor exactly 1 (the I2/I3 identities hold to 1e-39). On 75 states, \|H_ang − H_hyb\|/N_s ≤ 4.1e-3.
* **Strip scan (uncertified).** L ≤ 0.29 only, because n ≤ 1000. Min H_ang/N_s = 0.893 over 18 states, and agreement with the hybrid is ≤ 1.3e-4.
* **Rigour patch.** acb.integral was re-run with the holomorphy flag; every certificate stands.
* **Agent's suggested route for J → ∞ at the edge.** Bound the finite-mass frequency kernel per ν; it is O(ν²/T) from Legendre, with no exponential sizes. Alternatively, close the §34 Mellin boundary term at r = 1.

## 45. Exact finite-mass kernel: dual-Hahn convolution and a five-point equation (author-supplied, 26 Nov)

* **Notation.** n = J−b, N = Δ−b−2, a = N−n, d_b = 2b+1. U_k = (−n)_k(n+d_b)_k/((−N)_k(N+d_b)_k). Bare kernel 𝒦_p(ν) = ℒ(F_{p,ν}), with F_{p,ν} = (1−z)^{−p}P_{−1/2+iν}((1+z)/(1−z)).
* **Coefficients.** [z^k]F_{p,ν} = ((p)_k/k!)²·₃F₂(−k, A±iν; p, p; 1), with A = p − ½.
* **Closed form (author).**
  * 𝒦_p(ν) = U_n((p)_n/n!)·Σ_j (−n)_j(½±iν)_j/((p)_j(j!)²)·R_{n−j}, with R_m = ₃F₂(−m, a+1, −2n−a−2b; −2n−2b, 1−p−n; 1).
  * R_m is a dual-Hahn family with M = n+p−1, γ = −2n−2b−1, δ = 1, x = −a−1. This is outside the positive-orthogonality range.
  * Three-term recurrence in m.
  * **Verified here exactly** (Fractions, 40 random states, p ∈ {2,3,4}, b ∈ {0,4}, rational ν²): `code/angular/verify_s45_dualhahn.py`.
* **Not a single Wilson polynomial (author).**
  * The exact residual at (J,Δ) = (5,30) is 1586490242224566444/25787514494155770625.
  * The Wilson difference-equation determinant is nonzero mod 1000003 at four states in the open region.
  * Not re-checked.
* **Five-point equation (author).** [ν² + (p−½)²]K = m²[(∇_N − ∇_n)K + p²K/(N(N+d_b))], with m² = (N−n)(N+n+d_b) and the positive coefficients A_±(n), B_±(N) given by the author.
  * **Verified here exactly** on 40 random states.
  * It is not a maximum principle: ∇_N and ∇_n enter with opposite signs.
* **Dressing (author).**
  * 𝒦₂^{phys} = (1 + (9−b(b−1))/T)𝒦₃ − ((ν²+25/4)/(2T))𝒦₄ + ((7−2b)/T)[(3+ℓ)𝒦₄ − (3+ℓ/2)𝒦₃ − (ℓ/2)𝒦₃(ν+i)], with ℓ = −½+iν.
  * 𝒦₃^{phys} = 𝒦₄.
  * H = w³[(2/Λ)∫f𝒦₂ + (σ_bΛ²/T)∫ψ₃𝒦₃(Λp)dp].
  * Not yet matched against fp_kernel.kfin.
* **Open.** A non-exponential per-frequency envelope, perhaps from a discrete energy estimate for the five-point equation.

## 46. Weight conversion, Darboux expansion with 1/M correction (author, 26 Nov); cancellation size (coordinator)

* **Weight conversion (author).**
  * f(ν) = 2Re[χ_re(ν−i)θ(ν−i)R₋(ν−i)].
  * H^{(2)} = ∫g K_raw dν = (2/Λ)∫f Φ₂ dν, with Φ₂ = w³𝒱_b(F_{3,ν}).
  * K_raw = −((4ν²+9)(4ν²+49)/(4608(4ν²+1)))·Re[R₋(ν)Φ₂(ν+i) + R₊(ν)Φ₂(ν−i)].
* **Coordinator scan** (raw g paired with fp_kernel.kfin; the c_k = 0 Gaussian is not doubled; k = 3 piece omitted), at J = 200, Δ = 1200, b = 0:
  * H^{(2)} = 1.800e-4 and the Legendre-kernel H_LR^{(2)} = 1.762e-4, against N_s = 4.34e-5.
  * But ∫\|g\|·\|K_fin − K_LR\| dν = 2.0·10⁴ N_s. It comes almost entirely from ν = 50…1200: 9.6e2, 1.55e4 and 3.9e3 N_s in the bands [50,200], [200,500] and [500,1200].
  * **Consequence:** the ν-integrand cancels down to about 10⁻⁵ of its absolute size. Any per-frequency envelope must therefore either keep the oscillation (a signed ν-integral) or match the kernel to a relative accuracy of about 10⁻⁵ at ν ~ 10², which a C(ν)/T envelope cannot do. [DV]
* **Darboux expansion (author; diagnostic numerics only).**
  * Variables: M = J + ½ + a/2, α = p − ½ + iν, and 𝔑_a(α) = Γ(α)⁻¹∫₀^{2a}x^{α−1}e^{−2a·atanh(x/2a)}dx.
  * Expansion: 𝒦_p(ν) = 2Re{B(ν)d^{−α}[𝔑_a(α) + ℭ/M]} + O(d^{−Re α}M^{−min(2,Re α)}).
  * First correction: ℭ = a[(2b−α)𝔑_a(α) + (p−2b)𝔑_a(α−1)].
  * The centred index k + p/2 cancels the first-order coefficient correction.
  * Ladder expansion: U_k = ((1−x)/(1+x))^a[1 + (2a/M)(b − (b−p/2)/(1−x²)) + O(M⁻²)], with x = (k + p/2)/M.
  * Not yet checked here.

## 47. Darboux dressing, saddles, turning point, scales (author, 26 Nov)

* **k = 2 dressing.** ℭ̃_{a,b,3} = ℭ_{a,b,3} + (αA_b(α)/(2a))𝔑_a(α+1), with A_b = (9+α)/2 − 2b. The k = 3 piece uses the bare p = 4 expansion.
* **Darboux multiplier.**
  * The ₂F₁ form (2a)^α B(α, a+1)₂F₁(a, α; a+α+1; −1)/Γ(α) equals the atanh integral.
  * Large-a expansion: 𝔑_a = 1 − (α)₃/(12a²) + [(α)₆/288 − (α)₅/80]/a⁴ + O(a⁻⁶). This gives 𝔑_a(5/2) = 1 − 105/(32a²), as in §27.
  * **Verified here** (mpmath, a = 50 and 200, α = 5/2, 5/2 + 3i, 7/2 + 10i): the ₂F₁ form agrees to 1e-31. The expansion residual is 1e-8 to 3e-5 at a = 50 and 3e-12 to 7e-9 at a = 200, consistent with O(a⁻⁶). \|𝔑\| > 1 at α = 7/2 + 10i.
* **Saddle at fixed κ = ν/a (author).**
  * Saddle point: x_* = i(1 − √(1−κ²))/κ.
  * Amplitude: (2/(1+s))^σ s^{−1/2}, with s = √(1−κ²).
  * Phase: a[κ log(2/(1+s)) + κ − arcsin κ].
  * Above κ = 1 the multiplier grows like e^{aG(κ)}, with G = κ·arccos(1/κ) − arcosh κ.
  * Airy variable ζ = 2^{1/3}a^{2/3}(κ−1).
* **Phase scale (author).** log 𝔑 ≈ (σ+1)ν²/(4a²) + iν³/(12a²). Comparison with Legendre fails once ν ~ a^{2/3}. This is consistent with the coordinator's measured cancellation at ν ≳ 50 (§46).
* **Joint saddle (author, formal WKB).**
  * Equation: κ²[(1+τ)² − z²][(1−τ)² − z²] + 4z² = 0, with τ = a/(2M).
  * Turning points at κ = 1 and κ = 1/τ.
  * Dispersion relation: cos ϑ = 1 − ν²/(2y²).
  * Action: S(y) = ∫₀^y log((ℓ²−t²)/(h²−t²))dt + iΘ_ν(y).
* **Consistency checks (author).** At fixed c the result reduces to F_{p,ν}(c) with the first-order profile. At ν = 0 it gives 𝔑_a(p − ½). Both pass.

## 48. Airy-window bound on the Darboux multiplier and weighted edge integrals (author, 26 Nov)

* **Window.** ν₁(a) = a − a^{1/3}/2; at a = 998, ν₁ ≈ 993.003.
* **Bound (author, proof sketched).**
  * Statement: \|B(ν)𝔑_a(5/2+iν)\| ≤ 8a^{−1/3}e^{G_a(ν)} and \|B(ν)𝔑_a(7/2+iν)\| ≤ 16a^{−1/3}e^{G_a(ν)} for a ≥ 998 and ν ≥ ν₁. Here G_a = ν·arccos(a/ν) − a·arcosh(ν/a) for ν > a, and 0 otherwise.
  * Proof outline: deform the contour; integrate by parts on y ≤ 1 − 2a^{−1/3}; use the cubic curvature bound h ≤ G − (a/6)\|v − v_*\|³.
  * **Checked here numerically** (mpmath, ₂F₁ form) at a = 998, 1050 and 1500 over ν from ν₁ to 1.12a. The largest ratio to the bound is 0.19, for both p.
* **Weighted edge integrals (author, stated as certified).** ∫_{993}^∞\|f\|e^{G_{998}} < 5.65e-8; the derivative integral is < 1.62e-7. Not re-checked.
* **Caveat (author and coordinator).** This bounds the limiting multiplier, not the finite-spin kernel. FULLPOS2 §3 finds the j = 0 Darboux reference accurate only when ν·d ≪ 1 (L ≳ 7). So the bound controls the turning-point sliver only in the large-L part of the open region.

## 49. Edge lemma: high-frequency part ν ≥ ν₁ = a − a^{1/3}/2 (author, 26 Nov)

* **NF02 edge (author, exact rationals).** \|ψ₂\| ≤ 9700(1−p)⁵ and \|ψ₃\| ≤ 9500(1−p)⁵ on [0.994, 1]. Not re-checked.
* **Weighted density enclosures (author, outward-rounded Taylor panels).**
  * ∫_{993}^∞\|f\|e^{G₉₉₈} < 5.65e-8;
  * ∫\|f′\|e^{G₉₉₈} < 1.62e-7;
  * sup\|f\|e^{G₉₉₈} < 4.1e-8;
  * these hold for all a ≥ 998, since G decreases in a. Not re-checked.
* **Signed limiting-Darboux edge contribution (author).** \|H^D_edge\| < 0.002·N_s^rep.
* **Exact finite-state tail (author).** For a ≥ 1145, \|(2/Λ)∫_{ν≥ν₁} f·𝒦_phys\| < 0.01·N_s^rep at every spin.
  * Method: a crude coefficient bound e^{2πν}, beaten by the Gaussian decay of f above 1139.
  * **Coordinator check:** the weight's largest Gaussian centre is 1028.7 and its widths are σ ∈ {0.45, 0.9}. So f ≤ 10¹⁰e^{−(ν−1030)²/1.62} has the right centre and width.
* **Remaining for the edge.** For 998 ≤ a < 1145 and ν₁ ≤ ν ≤ 1139, the finite-kernel envelope \|𝒦_phys(ν)\| ≤ 10√w·e^{G_a(ν)} would give \|H_edge\| < 0.0067·N_s^rep.
* **Still open (main gap).** The bulk ν < ν₁ at 0.04 < L < 7 (see FULLPOS2 §§2–3).

## 50. High-twist closure (a ≥ 3·10⁸) and the Jacobi–Racah structure (author, 26 Nov)

* **Signed differential expansion (author).**
  * U_k/c^k = exp[−Σ_ℓ h^ℓ L_ℓ(k)], with h = 1/(T − u₀).
  * H = Σ_{m≤M} h^m P_m(ϑ)G(c) + ℛ_M, where deg P_m ≤ 3m and ϑ = z∂_z, acting on the full integrated generating function G.
  * Majorant: h^mQ_m(k) ≤ [h·k(k+9)(1+dk)/c]^m.
* **Positivity for a ≥ 3·10⁸ (author).**
  * The claim: H_ang,b > 4e-8 on 0.04 < L < 7, uniformly in J and in both channels.
  * Error budget: first-order term 1.035e-7, orders 2 to 1500 1.16e-8, ℛ < e^{−1000}, against a floor of 1.59e-7.
  * It uses LRINF_TIGHT's M₁ floor (> 7e-4 on 0.5 ≤ β ≤ 8.5).
  * **Constants not reproduced here.**
* **Jacobi → Racah (author). Verified here exactly** (`code/angular/verify_s50_racah.py`):
  * U_k − U_{k+1} = (1−c)(−n)_k(n+s)_k/((1−N)_k(N+s+1)_k), with s = 2b+1;
  * ℒ((1−z)Q_j) = (1−c)R_j(n), where Q_j = ₂F₁(−j, j+2; 1; z) and R_j = ₄F₃(−j, j+2, −n, n+s; 1, 1−N, N+s+1; 1), a balanced Racah polynomial;
  * orthogonality over n = 0..N−1 with ω_n = (1−c_n)(2n+s)/s > 0, plus the author's norm formula, for s = 1, 9 and N = 6, 11, 17.

  So for b = 0 and integer N: Σ_n(1−c_n)(2n+1)\|ℒ_n((1−z)P)/(1−c_n)\|² = N(N+1)∫₀¹(1−z)\|P\|²dz, an exact discrete–continuous isometry.
* **Minimum-norm extension (author).** P_min = Σβ_jE_j, with E_k = z^k₂F₁(−ℓ_k, ℓ_k+2k+3; 2k+2; z). Not re-checked.
* **Open.**
  1. The region 998 ≤ a < 3·10⁸ with 0.04 < L < 7 (J from 767 to about 6.6·10¹¹).
  2. A uniform bound on the residual extension's Jacobi norm. Its β_j are built from Taylor coefficients h_k, which are e^{O(Λ)} in size; whether they cancel is unknown.
  3. Physical Δ is real, so N is generically not an integer. The Racah orthogonality as proved needs integer N.

## 51. Racah route resolved structurally; high-twist threshold lowered to a ≥ 2.2·10⁷ (author, 26 Nov)

* **Helicity-adapted extension basis.** E_j = z^j₂F₁(j+1−N, j+N+s+1; 2j+s+1; z), weight z^{s−1}(1−z).
  * ℒ_{n,N}((1−z)E_j) = δ_{nj}d_nŨ_n. **Checked here exactly:** 80 random states (s = 1, 9; integer N), with no off-diagonal nonzeros.
  * Hence ‖P_min‖² = Σ_{j≤n}\|ℒ_{j,N}(R)\|²/A_{j,N,s}. **The minimum-norm certificate is circular:** it re-encodes exactly the lower-spin residuals that we need to bound.
* **Real-mass obstruction (author).** For non-integer N, Saalschütz gives ℒ((1−z)E_j^{M,δ}) ∝ (1+δ−m)_m, which no longer vanishes.
  * Counterexample at b = 0, J = 1000, N = 2000.5, with a polynomial of L²((1−z)dz) norm ≤ 1/√2: \|ℒ((1−z)P)\| ≈ 10^{217.87}. **Checked here:** the formula's log₁₀ is 217.99 + log₁₀ d_n.
  * So no uniform polynomial-space bound can extend the integer isometry to real masses.
* **High-twist closure (author, sharper rerun).**
  * Claim: H_ang,b > 1.6e-7·P_∘ > 0 for a ≥ 2.2·10⁷ throughout 0.04 < L < 7, both channels, all real masses.
  * Inputs:
    * moments of \|ψ_i\| (1024 dyadic Bernstein cells);
    * ∫\|f − ψ₂\| < 2.68e-6;
    * M₁/P_∘ > 0.0064 on 0.5 ≤ β ≤ 8.5;
    * floor 1.50e-6·P_∘.
  * Budget: first-order term 1.116e-6, higher orders 2.175e-7, ℛ₅₀₀ < e^{−650}. **Not reproduced here.**
* **Open.** 998 ≤ a < 2.2·10⁷ with 767 ≤ J < 4.83·10¹⁰ and 0.04 < L < 7. This needs a tame reference specialised to the physical generating function; generic norm arguments are excluded.

## 52. Airy representation of the cubic-resummed ladder; perturbative references disproved (author, 26 Nov)

* **Exact Airy form (author).**
  * Definitions: η = (1−c)/(C−b(b+1)), λ = η^{1/3}, r_* = c·e^{η(b²+b+1/3)}.
  * R₁ = e^{ηb³/3}∫₀^∞[Ai(u)e^{−bλu}F_{3,ν}(r_*e^{−λu}) + 2Re{e^{iπ/6}Ai(e^{iπ/6}u)e^{ibλu}F_{3,ν}(r_*e^{iλu})}]du.
  * **Verified here numerically** (`code/angular/verify_s52_airy.py`):
    * at J = 800, Δ = 1800 it matches the k-sum R₁ to 1e-18 to 1e-20 at ν = 5 and 50, in both channels;
    * at ν = 300 it matches the R₁ value implied by the exponentiated-ladder test.
    * There, ∫\|integrand\|/\|R₁\| ≈ 1.05 for ν ≤ 50.
* **Arc growth (author).** q_* = √r_*·λ/d_* ≤ 1.001/(4a²)^{1/3}. Hence d_*³\|F_{3,ν}\| ≤ e^{νq_*u} on the arc. The Airy envelope's amplification is exp(ν³/(6a²)).
* **No exact short-arc measure (author).** The moments c^k e^{−hL₁(k)} decay superexponentially, which forces an essential singularity in the Cauchy transform.
* **Truncation multipliers (author).** 0 ≤ δ_M(k) ≤ m_M(k)·min(1, A₁x_k^M/(1−x_k)). The correction formula is −(1/2π)∫F_{3,ν}(ce^{it})D_M(t)dt.
* **Counterexample (author; 1100-digit directed arithmetic).** At J = 768, Δ = 1768, b = 0 (inside the open strip), d³R₁ ≈ 1e27 and d³R₂ ≈ 2e17 at ν = 998, against an exact value of 0.054. **Truncating the ladder symbol at cubic or quintic order fails at ν ~ a.** The saddle scale grows like ν/a^{2/3}, and the quintic term becomes O(a).
* **Next reference.** It must keep the full logarithmic ladder symbol at the complex saddle: the §47 joint WKB, S(y) = ∫log((ℓ²−t²)/(h²−t²)) + iΘ_ν.

## 53. Full-ladder saddle: exact two-dimensional contour, phase and Hessian (author, 26 Nov)

* **Contour identity (author; not re-checked).**
  * 𝒦_p(ν) = (𝒞/2πi)∫₀^∞∮_Γ (du/u)(dv/v)·z^{s/2}·exp[ℓ log(u/v) − H log((1+u)/(1+v))]·F_{p,ν}(z), with z = uv/((1+u)(1+v)).
  * Parameters: ℓ = J + ½, H = ℓ + a, 𝒞 = aΓ(n+1)Γ(N+s)/(Γ(n+s)Γ(N+1)).
  * Exact, and valid for real N.
* **Saddle (author). Verified here** at 40 digits on 4 states (both channels, a = 998 to 5000, ν = 300 to 1000):
  * the saddle point u_* = (ℓ+iy)/a, v_* = (ℓ−iy)/a, where the gradient of Ψ₊ vanishes;
  * z_* = (ℓ²+y²)/(H²+y²), y = ν√z_*/(1−z_*), and the saddle equation ν²(ℓ²+y²)(H²+y²) = (H²−ℓ²)²y²;
  * Ψ₊(saddle) = iΦ(ν), with Φ = 2ℓ·arctan(y/ℓ) − 2H·arctan(y/H) + 2ν·atanh√z_*, and Φ′(ν) = β(z_*);
  * det D²Ψ₊ = −a²(ℓ²H² − y⁴)/(H² + y²)²;
  * the saddles coalesce at y² = ℓH, i.e. ν = a.
* **Open.** The global contour deformation with explicit error, uniform through the coalescence at ν = a (Airy/Chester–Friedman–Ursell), and the transport amplitude, all with constants.

## 53a. Correction to the coordinator's exponentiated-ladder table (§51 discussion)

* The Airy evaluation of R₁ gives:
  * J = 800, Δ = 1800: ν = 300 → 0.069456 (consistent with the k-sum); ν = 500 → −5.9e4; ν = 700 → 1.36e15;
  * J = 2000, Δ = 3000, ν = 700 → 8.5e16;
  * J = 800, Δ = 3800, ν = 700 → 0.02994 (consistent).
* The earlier k-sum test truncated the non-terminating sum near k ≈ n + 90. At ν ≥ 500 with c ≈ 0.2–0.45 that cutoff is too early. Its ν = 700 entries at J = 800/Δ = 1800 and J = 2000/Δ = 3000 (relative errors 0.34/0.14 and 1.9/0.5) are therefore **wrong**. R₁ actually blows up there, in agreement with the author's §52 counterexample.
* Airy-envelope amplification ∫\|integrand\|/\|R₁\|: 43 at ν = 300 and 27 at ν = 500 (J = 800, Δ = 1800); 51 at J = 800, Δ = 3800, ν = 700.

## 54. Exact separation of the saddle variables; exact global contour (author, 26 Nov)

* **Angles.** x = e^{iα} = √(B/q), y = e^{iγ} = x(1−q)/(1−B), ξ = √z = sin((γ−α)/2)/sin((γ+α)/2), and e^β = cot(α/2)/cot(γ/2).
* **Separated phase.** Ψ₊ = i[𝓕_S(α) − 𝓕_a(γ)], with 𝓕_m(t) = mt + ν log cot(t/2) and S = H + ℓ = 2J + 1 + a.
  * Saddle points: sin α_* = ν/S (hard, nondegenerate) and sin γ_* = ν/a (soft, coalescing at κ = 1).
* **Checked here (40 digits, 5 states including J = 5·10⁶, a = 2·10⁴):**
  * the saddle phase 𝔉(S,ν) − 𝔉(a,ν), with 𝔉(m,ν) = ν·arcosh(m/ν) + m·arcsin(ν/m), equals the §53 phase Φ(ν) exactly;
  * ξ_*² = z_*;
  * the e^β identity holds at the saddle.
* **Airy variable.** ζ = {(3a/2)[arccos κ − κ·arcosh(1/κ)]}^{2/3}. It depends on a and κ only, not on τ. Checked against 2^{1/3}a^{2/3}(1−κ) at κ = 0.999.
* **Exact global contour (author; diagnostics only).**
  * 𝒦_p = (𝒞/π)Re(I_U − iI_T − iI_L), built from the outgoing Legendre solution.
  * Pieces: a bounded oscillatory arc, a damped vertical tail and a lower ray, with \|ξ\| < 1 and Re ξ ≥ 0 throughout.
  * Not re-checked.
* **Structural consequence (coordinator).** The phase is a hard factor in (S, α), with large parameter S ≥ 2J + 1 + a ≥ 2535, times a soft factor in (a, γ) that depends only on (a, ν).
  * Only the amplitude, through P⁺(β(α, γ)), couples them.
  * So the leading uniform approximation should be: a nondegenerate Gaussian in α with O(1/S) corrections, times an exact one-dimensional soft integral in γ. The latter is a Darboux-multiplier-type special function of (a, ν) that can be evaluated in Arb uniformly through the turning point.

## 55. Exact non-amplifying contour and a certified reference-kernel construction (author, 26 Nov)

* **Global bound (author).**
  * Building blocks: the outgoing Legendre solution P⁺ with \|P⁺\| ≤ e^{−ν Im β}/√(2ν·sinh(Re β)), and Im β ≥ 0 on all three contour pieces.
  * The vertical tail is damped by e^{−(a−ν)t}, the lower ray by e^{−av}.
  * Consequence: \|𝒦_p(μ)\| ≤ (𝒞D_p/π)[3π²/4 + 2π/(a−ν) + 6π/√(aε)] for 1 ≤ ν ≤ a − a^{1/3}/2 and μ ∈ {ν, ν±i}, with D_p = 16e(2/ε)^{p+1}/(ε√ν) and ε = 1/S.
  * **No e^{O(ν)} or e^{O(a)} factors** (the prefactor is polynomial in S).
* **Certified reference (author).** Truncate the contour at δ, T±; use Taylor panels of radius R_an = (ε/100)·min(εδ²/20, 1/(1+ν)); truncate the P⁺ series at K terms.
  * Explicit total error: E_cut + PV·M_an[e^{−Kρ}/(1−e^{−ρ}) + 2(m+3)2^{−m}] + E_round.
  * Efficient form: integrate the carriers e^{iSA−iaC} etc. exactly per panel, with an Euler-integral evaluator for P⁺.
* **Leading saddle (author, diagnostic).** At J = 768, Δ = 1768, b = 0, the absolute error in d³𝒦₃ is 1.4e-5, 5.8e-5, 2.2e-4 and 9.3e-3 at ν = 300, 700, 950 and 993. The transport correction 𝒯_f is given.
* **Transfer to the action (author).** \|H − H_ref\| < 55η when w^p\|𝒦_p − 𝒦_p^ref\| ≤ η.
* **Open (coordinator).**
  1. Numerically verify the contour identity.
  2. An efficient certified evaluator.
  3. **Uniformity over the infinite family** (J up to ~4.8e10, a from 998 to 2.2e7). This needs either Lipschitz control in (a, 1/S) on a grid, or certified asymptotics: the hard α-integral in 1/S, and the soft γ-integral at large a (the 𝔑_a expansion). Each state's reference is a finite computation; the family is not.

## 56. Joint compactification (τ, t) = (a/S, 1/S²); uniform large-twist theorem; parameter Lipschitz/Taylor cells (author, 26 Nov)

* **Contour identity (author, numerical).** Matches the terminating kernel at three states, including non-integer N (N = 7.5, 10.25), to 1e-12 (one case 8e-9). The identities u/v = xy and (1+u)/(1+v) = y/x give the carrier e^{iSα−iaγ}.
* **Exact relations. Checked here by hand:**
  * r = ((1−τ)² − t)/((1+τ)² − t);
  * w = 4τ/D and 1/T = 4t/D, with D = (1+τ)² − t;
  * U_k = ∏_{j<k}((1−τ)² − (s+2j)²t)/((1+τ)² − (s+2j)²t).

  These follow from C = ℓ² − ¼, T = H² − ¼ and C − u_j = ℓ² − (s+2j)²/4. The moderate strip is 0.00022 < τ < 0.68.
* **Uniform theorem (author).** \|H − Σ_{m≤1300} t^m P_m(ϑ; τ)G(c₀)\| < e^{−55} for all states with a ≥ 125000, both channels and real masses.
  * Inputs: sup_{\|z\|=(1+c₀)/2}\|G\| < e^{1710}, and an analytic ladder in complex t.
  * Consequence: H_ref is a degree-1301 polynomial in y = (125000/a)² ∈ [0,1], with τ-dependent coefficients. A Bernstein certificate in y plus enclosures in τ would cover the whole unbounded family.
  * A per-frequency version holds for a ≥ 150000.
* **Compact remainder region (author).** 998 ≤ a < 125000, 767 ≤ J < 2.85·10⁸, 0.04 < L < 7.
  * Exact parameter derivatives: ∂_S^j∂_a^k I = ∫(iα)^j(−iγ)^k𝒲_p.
  * Conservative bounds with e^{G₉₉₇(ν)}; ∫(1+ν²)\|f\|e^{G₉₉₇} < 1e14.
  * Signed parameter-Taylor cells with explicit remainder.
* **Coordinator cautions.**
  1. The retained coefficients P_m(ϑ)G(c₀) involve derivatives of G up to order ~3900 at c₀. They may be individually large and cancel across m; that needs high-precision Taylor models in τ, not naive interval enclosures.
  2. The §50 majorant (e^{3.1m}[B_m³/a²]^m) does not converge at a = 125000. Only the new global-radius argument does; reproduce it before relying on it.

## 57. Local saddle majorants (no positive power of S); tails removed; nonstationary complement open (author, 26 Nov)

* **Saddle widths.** h_A = √ν/(S√cos A_*) and h_C = √ν/(a√cos C_*), with sin A_* = τκ and sin C_* = κ.
* **Local amplitude majorant 𝓜_loc.** Built from ρ(x), θ(x), q±, y± and the bound on P⁺. It replaces the global ε^{−6}(4/ε²)^p.
  * It is sandwiched against the saddle value: (1−δ_*)𝓜_* ≤ \|𝒲\| ≤ 𝓜_*, with δ_* < 0.216/ν.
  * h_Ah_C𝓜_* = O(√τ/(√ν(1−κ²)^{1/4})), with no positive power of S.
  * Diagnostic: integrated local bound 0.0088 at J = 10⁶, against ~1e77 before.
  * q_* = tan(A_*/2)/tan(C_*/2) = τ(1+cos C_*)/(1+cos A_*). Checked.
* **Parameter cells (author).** Centring the carrier at the saddle angles gives \|∂_S^j∂_a^k I_cent\| ≤ 𝓛_loc R_A^j R_C^k. The cell remainder is 𝓛_loc(r_S R_A + r_a R_C)^m/m!, so cells of width r_S = O(S/√ν) and r_a = O(a√cos C_*/√ν) are possible.
* **Tails (author).**
  * Vertical: h(t) ≥ (a−ν)t + νt³/9 (1 − sech t ≥ t²/3 on [0,1]; checked at the endpoints). Uniformly E_{T,all} < 1e-28 for 998 ≤ a ≤ 125000.
  * Lower ray cut at v = 300/a: omitted part < 1e-84.
* **Open.**
  1. The finite nonstationary complement of the contour. The plan is integration by parts with \|B_m\| ≤ M(m/(RD))^m and panel-specific gaps D from 𝓕′_S(α) = S(sin α − sin A_*)/sin α.
  2. The endpoint and diagonal neighbourhoods, including the logarithmic singularity of P⁺ at ξ = 0.
  3. The final signed positivity grid.

## 58. Frequency integration done once: universal rapidity transforms; dressing without shifted frequencies (author, 26 Nov)

* **Universal transforms (author).**
  * 𝒫_h⁺(β) = ∫h(ν)P_ν⁺(β)dν = −(ie^{−β/2}/π)∫₀^∞e^{−u/2}W_h(β+u)/√((1−e^{−u})(1−e^{−2β−u}))du, with W_h(v) = ∫h·coth(πν)·e^{iνv}dν.
  * These are functions of one complex rapidity, independent of (J, a, b). Needed for h = f, ν²f and ψ₃(ν/Λ); 0 ≤ ν ≤ 1 is handled with the full conical function.
  * **Checked here:** the Euler representation of P⁺ against the ₂F₁ form (3e-16 at ν = 2.3, 1.5e-6 at ν = 9.1, the latter limited by oscillatory quadrature), and P⁺ + c.c. = P to 1e-31.
* **First-moment operator without shifts (author). Verified here** to 1e-30 at three (ν, z, b):
  * 𝒟_bF_{3,ν} = (1−z)^{−4}[(6 − b(b−1) + (12 − 6b + b(b−1))z − (ν²+¼)/2)P_ν + (7−2b)ξ∂_βP_ν], with z = tanh²(β/2);
  * it holds for each outgoing solution separately.
* **Consequence.** The dressed contour amplitude needs only 𝒫_f⁺, ∂_β𝒫_f⁺, 𝒫_{ν²f}⁺ and 𝒫_{ψ₃}⁺, precomputed once on the relevant β-region. Each state-cell then costs one two-dimensional contour integral, not ~1000.

## 59. Global localisation of both contour angles; complete cell-certificate recipe (author, 26 Nov)

* **Radial-soft representation (author). Verified here** at n = 10, N = 20.5, b = 0, ν = 1.2, p = 3, against the terminating sum: 1.357403299075 both ways, with imaginary part 3e-24.
  * 𝒦_p = −(i𝒞/2π)∫_{−π/2−iε}^{π/2−iε}(e^{iBα}/sin α)∫₀^∞e^{−av}ξ^{s−1}F_{p,ν}(ξ²)dv dα, where ξ = sinh(v/2)/sinh(v/2+iα) and B = 2J+1.
  * It uses the complete conical function, so there is no ξ = 0 logarithm.
* **Hard localisation (author).** A₀ = min(π/2, 8V/B) and Y = 250/B, with V = 1100. The outer hard intervals move to Im α = Y and the vertical connectors are kept.
  * The \|Re α\| = π/2 boundaries cancel because B is odd.
  * Damping is e^{−0.8By}, with no positive power of spin.
* **Soft localisation (author).** C₀ = min(π/2, max(3A₀, 8V/a)), T = 250/a, lower ray cut at 500/a.
  * Budget: w^p\|𝒦_p − 𝒦_p^trunc\| < e^{−70} for cells of relative half-width ≤ 1%.
* **Retained-contour majorant (author).** < 10⁶⁶, normalised, with no growing power of J.
* **Spectral side (author).**
  * The ν > 1100 tail is < e^{−1200}.
  * Degree-250 Taylor panels for f, with \|f − f̃\|_{L¹} < 1e-110, give a physical error < 1e-43.
  * Exact Fourier integrals of the polynomial panels; a coth expansion to 48 terms.
  * A Mehler form of the frequency-integrated full conical function.
* **Diagonal logarithm (author).** A zero-balanced ₂F₁ expansion in ω = 1 − e^{−2β}, with exact log moments ∫₀^ℓ d^j log d = ℓ^{j+1}[log ℓ/(j+1) − 1/(j+1)²].
* **Parameter cells (author).**
  * The retained variables satisfy \|α\| ≤ 9100/B₀ and \|γ−α\| ≤ 64000/a₀.
  * Carrier Taylor remainder: E ≤ 10⁶⁶e^Z Z^m/m!, with Z = 9100r_B/B₀ + 64000r_a/a₀.
  * Example: 0.1% relative cells need order 512 for 1e-114.
* **Final test.** min_cell H_ref > E_contour + E_spectral + E_parameter + E_quadrature, checked with Bernstein coefficients in (B, a).
* **Status.** The recipe is complete on paper. Not executed: the certified spectral panels, the cell evaluations, the a ≥ 125000 polynomial certificate. No constant in §§55–59 has been independently reproduced. Cost is undetermined; with the conservative 10⁶⁶ majorant the worst-case orders are high.

## 60. FULLPOS3 interim: the signed h-expansion closes high twist down to a ≈ 5000 (first cell) (agent, 26 Nov)

* **Key finding [DV and exact].** The signed finite-mass correction is tiny even at a = 998.
  * Exact Σ_k U_kM_k (with \|M_k\| up to 6e627) equals G(c) to a relative 4.5e-7 to 3.6e-6.
  * The Legendre kernel at c has signed integrated error ≤ 3e-6·N_s at a = 998, J = 800 to 5000, even though its absolute error is ~3e4·N_s.
  * Mechanism [P]: ψ₂ vanishes at p = 0 and has a (1−p)⁵ zero at p = 1. Integrating by parts, the e^{Λ Im ξ} growth sits only in the p = 1 endpoint term, which is ~\|ψ⁽⁵⁾(1)\|/(Λβ)⁶ ~ 1e-11. The smooth part is of leading-Regge size.
  * So per-frequency accuracy is not needed: endpoint-aware sup bounds make the h-expansion terms signed-safe.
* **Rigorous h-expansion certificate** (code/fullpos3/fp3_ht.py, Arb). Ingredients: Lemma A, the EPSW residual by panels, exact P_m for m ≤ 30 plus a majorant up to M ~ 300–1600, stadium c-cells, and a far remainder.
  * First cell (β ∈ [2, 2.05], L ≈ 1.06, b = 0): **a₀ = 5000 passes** (error / F_cert = 0.0019). a₀ = 1500 fails, because the far remainder needs d·n ≳ 2300.
  * Full scan running: 60 β-cells × 2 channels.
* **Other results.**
  * §53 contour identity: proved termwise, and checked to 15 digits at non-integer N.
  * R_M scan: I₄ ≈ 4e-4 N_s at a ≥ 3000.
  * Per-J certificates 767 ≤ J ≤ 1004 (running): all PASS so far.
* **Consequence.** The residual region becomes compact in a: 998 ≤ a < a*(L) ~ few·10³ to 10⁴, with 0.04 < L < 7, J ≥ 1005. It is still unbounded-ish in J (J ≲ a*·2e^L). This is the region for the §§55–59 contour cells, which is now far smaller than 998 ≤ a < 125000.

## 61. Sharper far remainder via the rational ladder function; a cell-endpoint bug in fp3_ht.far_remainder (author, 26 Nov)

* **Lemma (author).** 𝒰_k(h; c) = ∏_{j<k}(c − hu_j)/(1 − hu_j) is analytic for \|h\| < 1/u_{k−1}, a larger region than the log-series radius c/u_{k−1}.
  * Exact sup on \|h\| = ρ: max{(c+x)/(1+x), \|c−x\|/(1−x)}, with x = ρu_j.
  * Cell majorant M_k(ρ) uses c₋ and c₊ jointly.
  * Cauchy remainder M_k(ρ)(h/ρ)^{M+1}/(1 − h/ρ). It stays valid past termination, since the numerator vanishes at j = n.
* **First cell (author).** β_c ∈ [2, 2.05], both channels, a ≥ 1250, M = 850: \|ℛ₈₅₀\| < e^{−51}.
  * Inputs: ρ_k = 3/(4k(k+s)), M_k ≤ e·0.6768^k, \|g_k\| ≤ e^{1650}(k+5)⁵, split at k = 4500.
  * Diagnostic at M = 600: log(far bound) at a = 1500 improves from +220 to −253.
* **Bug (author; confirmed by the coordinator on reading the code).** fp3_ht.py line 219 sets c = cb.upper() and uses φ = log(1 − ru) − log(1 − ru/c) for log(1/E_k). But 1/E_k decreases in c, so the upper endpoint gives a LOWER bound on that factor, not an upper bound. With c^k, the combination c²/(c−x) is not monotone in c. **Cell certificates that used far_remainder() are not valid until this is fixed.** The joint majorant M_k(ρ) fixes it.

## 62. Far remainder at the true floor a = 998 (first cell): vertical-line Cauchy contour plus the termination zero (author, 26 Nov)

* **Claim (author).** \|ℛ₇₀₀\| < e^{−95} for β_c ∈ [2, 2.05], all a ≥ 998, b ∈ {0, 4}. Previously the bound was e^{+238}.
* **Telescoping (author).** For k > n, U_k(v) = (1 − v/h)V_k(v), because c = hu_n. Hence T_M(k, h) = h^M[v^M]V_k(v).
  * **Checked here:** [v^m]U = V_m − V_{m−1}/h, and the h^m-weighted sum telescopes.
* **Vertical-line remainders (author).** ρ_k = θ/u_k with 0 < θ < 1, which lies left of all the poles 1/u_j, j < k.
  * k ≤ n, θ > c: U_k(h) − T_M = (h^{M+1}/2πi)∫_{Re v=ρ_k} U_k(v)/(v^{M+1}(v−h))dv.
  * k > n: −T_M = −(h^M/2πi)∫ V_k/v^{M+1}dv, which needs no ρ_k > h.
  * Closing to the left is justified because U_k is bounded at ∞. Structure checked; the bound along the line (not yet received) and the e^{−95} numerics are not reproduced.

## 63. Far remainder < e^{−74} uniformly over the whole strip at a ≥ 998 (author, 26 Nov)

* **Coefficient envelope (author).** Binning the repaired Gaussians gives ∫\|f\|e^{πν/2} < e^{1624}, hence \|g_k\| ≤ e^{1625}(k+5)⁵.
* **Vertical-line product majorant (author).**
  * v = ρ_k(1+ix), with ρ_k = θ/u_k and θ = (1+c₋)/2 > c₊.
  * Squared modulus of each factor: R(c, z, q) = [q(c−z)² + (1−q)z²]/[q(1−z)² + (1−q)z²], with q = 1/(1+x²) and z_j = θu_j/u_k.
  * **Checked here:** R ≤ 1, since numerator − denominator = q(c−1)(c+1−2z) ≤ 0 for z ≤ (1+c)/2.
  * The maximum over a (c, z, q) box sits at its corners (author's convexity and monotonicity argument). Not re-checked.
* **Bound (author).** \|U_k − T_M\| ≤ (100/(1−θ))(hu_k/θ)^M sup_q q^{(M−1)/2}∏max{0.01, √R_j}.
  * For k ≥ 200 the product is grouped into 128 bins with a lattice-discrepancy correction, giving Σlog ≤ kA + V.
  * The k-dependence is concave, with maximum at k_* = (2M+5)/(−A) − 5.
  * Physical bound h ≤ h_* = min{(1−√c₋)²/998², c₊/1001²}.
* **Cover (author; outward-padded elementary arithmetic, not Arb).** 283 c-cells on [0.038, 0.9991], 662 q-panels, M = 600, both channels.
  * Contributions: k < 200 gives e^{−659}; middle e^{−76.4}; tail e^{−4472}.
  * **Total \|ℛ₆₀₀\| < e^{−74} for all a ≥ 998, J ≥ 1005, 0.04 < L < 7.** M = 650 gives < e^{−120}.
* **Remaining for full positivity.** The local signed terms must pass at a₀ = 998:
  * F_cert > E_residual + E_dressing + Σ_{m≤600}E_m + e^{−74};
  * plus Arb reproduction of this cover and an independent audit.

## 64. Local terms: factor out the endpoint oscillation before differentiating (author, 26 Nov)

* **Problem.** The current local certificate takes absolute values of the coefficients of P_m(ϑ) before applying them to G. That loses the frequency cancellation and the cancellation between orders.
* **Exact split (author).** Ω_ψ(β) = S_ψ(β) + e^{iΛβ}A_{ψ,+}(β) + e^{−iΛβ}A_{ψ,−}(β).
  * S_ψ holds the p = 0 endpoint terms and the damped coth correction.
  * A_{ψ,±} = (±√2/(2πi))∫₀^∞K(β,u)e^{±iΛu}Σ_{j=5}^{14}(−1)^jψ^{(j)}(1)/(±iΛ(β+u))^{j+1}du, with K = [2 sinh(β+u/2) sinh(u/2)]^{−1/2}. The sum starts at j = 5 because of the fifth-order zero at p = 1.
  * Bounds: \|A_{ψ,±}(β)\| ≤ (√2/2π)J_e(β₀−R) on \|β−β₀\| ≤ R, with **no e^{ΛR}**. Derivatives follow by Cauchy, using the J_e tables in fp3_htdata.
* **Use.** Apply P_m(ϑ) to S_ψ and to e^{±iΛβ}A_± separately. On the oscillatory part, the operator acts on the exponential exactly (powers of Λ times the chain rule in β(z)). Combine orders before taking absolute values. Not yet implemented.

## 64b. Local terms, continued (author, 26 Nov)

* **Smooth high orders (author).** (1−c)³Σ_{m=20}^{600}h^m\|P_m(ϑ)G_sm(c)\| < e^{−55} on the whole strip.
  * Method: a closed Faulhaber/Bernoulli majorant with F, H, Q; 283 cells; outward-padded elementary arithmetic.
  * Exact rationals: Σ\|coefficients of ψ₂\| < 9,163,906 and Σ\|coefficients of ψ₃\| < 18,937,299.
  * (1−c)³ sup\|G_sm\| < e^{75}.
  * 𝒟_b in t = log(z/c): ((1−z)/(2z))∂_t² + (4−2b)∂_t + 9/2 − b(b−1).
  * Not reproduced.
* **Endpoint carrier (author).** G_end(ce^t) = e^{iωt}A(t) with ω = Λ√c/(1−c). The remaining phase grows only quadratically in R.
  * Shifted-operator identity: P_m(∂_t)[e^{iωt}A]₀ = P_m(∂_t + iω)A(0).
  * E_m^end ≤ S_A(R)Σ_j\|h^mP_m^{(j)}(iω)\|R^{−j}.
* **Warning (author, diagnostic).** At J = 1996, a = 998, b = 0, the frozen-symbol size \|h^mP_m(1198i)\| reaches 10^{34} around m ~ 80. Individual-order bounds may fail; cancellation across m may be needed.
* **Grouped fallback (author).** 𝒬_M = Σ_{m≤M}h^mP_m.
  * Recurrence A_m(k+1) = A_m(k) + u_k[cA_{m−1}(k+1) − A_{m−1}(k)], with A_m = c^mP_m. **Checked here by hand** from (1−hu_k)E_{k+1} = (1−hu_k/c)E_k.
  * Newton form: 𝒬_M(ϑ)G(c) = Σ_j Δ^j𝒬_M(0)·c^jG^{(j)}(c)/j!, which is exact for polynomial 𝒬. Checked.
  * A Taylor-cell remainder is given explicitly.
* **Residual r = f − ψ₂ (author).** Frequency-bin carriers e^{±iν_Iβ} with an amplitude growth of only e^{Δ_I\|Im β\|}.

## 65. FULLPOS3 progress (2): bug fixed; local terms bind below a ≈ 3000 (agent, 26 Nov)

* **far_remainder() replaced** with the §61 joint-cell circles. The old scan is discarded.
  * Far bounds: β ∈ [2, 2.05] gives 2e-58 at a₀ = 1500 but 8e156 at a₀ = 998; β ∈ [0.39, 0.44] gives 2e-55 at a₀ = 1250.
  * The §62/§63 vertical lines (e^{−74} at a = 998) are **not yet implemented or reproduced**.
* **Local terms bind.** At β ∈ [0.39, 0.44], b = 0:
  * a₀ = 3200 **passes** (err/F_cert = 0.097, exact P_m to m = 70, 81 radii);
  * a₀ = 2000 **fails**: per-order absolute bounds reach 1e158 near m ~ 300, from the ψ₂ p = 1 endpoint and the residual cells at ν ≈ 1000–1060. This matches the author's §64b warning.
* **Scan.** 60 cells × 2 channels, bisection in a₀; about 4 h; expected a* ≈ 3000 at low L.
* **Per-J certificates complete:** every channel with J ≤ 1000, plus b = 4 to J = 1004. **All PASS** at every mass. [certified]
* **Remaining.** 998 ≤ a < a*(L) ≈ 3000 with J ≥ 1005, 0.04 < L < 7. Fix: §64b carrier-shifted or grouped 𝒬_M bounds for the endpoint and high-residual components, plus §63 vertical lines for the far remainder below a ≈ 1500.

## 66. Grouping in h is decisive: frozen-symbol test (agent, 26 Nov) [DV]

* **Test.** At the carrier k = iλ (b = 0, h = h_max): per-order max_m\|h^mP_m(iλ)\| versus grouped \|𝒬_M(iλ) − 1\| → \|e^{−A(iλ)} − 1\|.
  * c = 0.042, λ = 210/222, a₀ = 3200: per-order 16/32, grouped 1.66/1.12.
  * a₀ = 2000: per-order 9.8e3/7.0e4, grouped 1.58/0.064.
  * a₀ = 1400: per-order 9.7e8/6.2e10, grouped 1.72/0.65.
  * a₀ = 998: per-order 8.1e18/4.6e22, grouped 0.34/0.18.
  * c = 0.58, λ = 1810, a₀ = 998: per-order 5.6e34, grouped 1.67.
* **Conclusion.** Per-order absolute bounds lose e^{\|A\|}-type factors. The resummed multiplier stays O(1), at most ≈ 2.7 in modulus (the Darboux-multiplier size).
* **Estimate.** The derivative sum Σ_i\|𝒬^{(i)}(iλ)\|R^{−i} is controlled if λR ≳ 3\|A\|. That holds down to a ≈ 1500–2000 at low c. At a = 998 it needs R ≈ 0.6, where the nonlinear-phase factor is about e^{30}: marginal.
* **Certified (scan).** β ∈ [0.39, 0.44]: a* = 3200 in both channels (err/F = 0.096 and 0.117). The scan continues.

## 67. FULLPOS3 progress (3): certified a*(L) in the first cells; the low-L corner closes (agent)

* **Certified a*, both channels** (fixed far remainder, exact P_m to m = 70, 81 radii):
  * β ∈ [0.39, 0.44]: 3200 (err/F = 0.096/0.117);
  * β ∈ [0.49, 0.54]: 3200 (0.067/0.075);
  * β ∈ [0.59, 0.64]: 3200 (0.063/0.067);
  * β ∈ [8.24, 8.40]: 4000 (0.907/0.907; exact orders 1–30 give 0.90·F_cert, F_cert = 2.7e-7).
  * The remaining cells need about 5 h.
* **Low-L corner closed** (pending the intermediate cells). n ≈ a√c/(1−√c), so every state with L ≤ ~0.058 and a < 3200 has J ≤ 1000, which the per-J certificates cover.
* **Remaining compact family.** 998 ≤ a < a*(L) ≈ 3200–4000, J from 1001 to ~a*√c/(1−√c) (≈ 9e6 at L = 7).
  * Tools: grouped carrier bounds (§§64b, 66) and the vertical-line far remainder (§63). Neither is implemented yet.
  * Caveat (agent): at a ≈ Λ the endpoint carrier sits at 0.5–0.83 of the h-series convergence radius.

## 68. The naive grouped-carrier certificate fails; what is needed instead (agent, 26 Nov) [DV]

* **Test.** β = 0.415, a₀ = 2000, b = 0. Bound Σ_±S_B(R)Σ_i\|q_i\|i!R^{−i}, with q_i the Taylor coefficients of 𝒬₇₀ − 1 at ±iλ (λ = 213).
  * Best case: 2.5e5·F_cert at R = 0.4; per-order: 2.1e6·F.
  * R = 0.6: grouped sum 1.67, but S_B ~ 3e4 (quadratic phase e^{33}).
  * R = 0.1: S_B small, but grouped sum 1e103.
* **Cause.** Grouping fixes the value at the carrier but not the high k-derivatives. A degree-D polynomial approximating e^{−αk³} forces R ≳ α^{1/3}D^{2/3} ≈ 0.4, and there the carrier's nonlinear phase costs e^{ΛR²\|ξ″\|/2}. **Taylor-expanding the amplitude at t = 0 is the lossy step.**
* **What is needed (agent's view).** Treat the ladder multiplier as an integral (Airy-convolution) operator acting on e^{iλt}B(t). It samples B only on \|t\| ≲ few·α^{1/3} ≈ 0.05. That is the §52 Airy form, but with the endpoint carrier factored out and the full ladder symbol, not the cubic truncation. Not attempted.
* **Scan.** 5/60 cells done (a* = 3200 at low L, 4000 at L ≈ 7); continuing.

## 69. §59 contour on the integrated G: full conical function unusable; outgoing split needed (agent, 26 Nov) [DV, Arb]

* **Full-conical form fails.** The near-diagonal contour points (small v, \|x\| ~ ε) land at complex z with Im arccosh X ≈ 0.3–0.7, carrying weight O(1).
  * G at z = 0.0865 + 0.162i is ≈ 5e201 (Σ\|terms\| = 1.4e203), against G(0.3) = 4e-5.
  * ∫\|integrand\|/\|ℒ(G)\| ≳ 1e200.
* **What should work (agent).** The outgoing split G⁺ = ∫ψP⁺ (author §58 𝒫⁺) on the Im β ≥ 0 contour (§55).
  * Its endpoint structure matches Lemma A: the p = 0 terms are floor-sized, and the p = 1 terms e^{iΛβ}/Λ⁶ decay for Im β ≥ 0.
  * So \|G⁺\| should be a few × floor and ∫\|·\| ≈ (log B)·few·F. Not built.
  * Uniformity: parameter cells in B need relative width ~1/8800 before carrier extraction.
* **Scan (restarted after the container restart).** a* = 3200 on β ∈ [0.39, 0.74] (4 cells), 4500 on [7.84, 8.04], 4000 on [8.24, 8.40]. Both channels.

## 70. FULLPOS3 final report (agent hand-back, 26 Nov)

* **Verdict: not closed.** The open band is reduced to a compact family.
* **Certified.**
  * Per-J certificates: every channel with J ≤ 1000, and b = 4 up to J = 1004, at every mass.
  * High-twist theorem (Arb, fp3_ht.py), both channels:
    * a* = 3200 on 0.39 ≤ β_c ≤ 1.04;
    * a* = 4500 on 7.44 ≤ β_c ≤ 8.04;
    * a* = 4000 on 8.24 ≤ β_c ≤ 8.40;
    * 1.04 < β_c < 7.44 still running (nohup; ~13 min per cell; rerun skips finished cells).
  * Low-L corner (L ≲ 0.058) closed.
* **Remaining.** 998 ≤ a < a*(L), J ≥ 1001, up to ~1e7 at L = 7.
* **Not done.**
  * The outgoing-split, frequency-integrated contour (spec in FULLPOS3 §5.2).
  * Reproduction of the §§62–63 vertical-line far remainder.
  * The §54 global contour as literally written gives NaN at a test state; it needs the author's diagonal and branch conventions.

## 71. Lemma A⁺: signed outgoing frequency-integrated amplitude (author, 26 Nov)

* **Representation.** 𝒫⁺_h(β) = −(i/(√2π))∫₀^∞K(β,u)W_h(β+u)du, with K = [2 sinh(β+u/2) sinh(u/2)]^{−1/2} and W_h(v) = ∫h(ν)coth(πν)e^{iνv}dν.
  * Derivative: ∂_β𝒫⁺ = −(i/√2π)∫K[iW_{νh} − ½coth(β+u/2)W_h]du.
* **Ideal weights H_{i,j} = ν^jψ_i(ν/Λ).** Integration by parts in ν terminates exactly (polynomials).
  * Bound: \|W(X+iy)\| ≤ A_{i,j}(X) + 1.004e^{−Λy}E_{i,j}(X) for y ≥ 0, integrating over the full support [0, Λ], since H(0) = 0.
  * The upper-endpoint sum starts at r = 5 (fifth-order zero).
  * **Checked here (exact NF02 polynomials, sympy):** E_{2,0}(0.37) = 5.545e-7, E_{2,1} = 5.63e-4, E_{3,0} = 5.43e-7, E_{3,1} = 5.52e-4; L_{2,0}(0.37) = 1.905, L_{2,1} = 10.84, L_{3,0} = 0.508, L_{3,1} = 2.886. All are at or below the author's bounds. H′_{2,0}(0) = 0.18918 < 0.189275. H^{(r)}(Λ) = 0 for r < 5.
* **Residual.** Per-panel EPSW masses, R_j(y) = Σε_ℓ coth(πl_ℓ)sup ν^je^{−νy}, \|r\|₁ < 2.682e-6. The erf-cutoff difference is < 1e-11.
* **Vertical contour (Im β < 0).** e^{−at−νy} ≤ e^{−t−νθ(A)}e^{G_{a−1}(ν)}, with G₉₉₇(Λ) < 0.055. No e^{πΛ/2} loss.
* **Diagonal.** 𝒫⁺_h = (i/π)[M_h log(2β) − D_h] + E_h, with explicit \|E_h\| and \|E_h′\| bounds for \|β\| ≤ 1/(32·1101). ξ∂_β𝒫⁺ is finite at the diagonal.
* **Closed majorants.** 𝓑_x[g] ≤ e^{−x/2}g(x)/√(1−e^{−2x}) for decreasing g. Leading size (2/Λ)\|𝒫⁺\| = O(e^{−x/2}/(Λ²x)), the leading-Regge floor scale.
* **Not yet evaluated:** the resulting full contour-integrand norm relative to N_s over the cells.

## 72. Lemma B4: helicity reduction; the corner cancellation and how to remove it (author, 26 Nov)

* **Exact.** H₄ = (w³/U_{0,4})ℒ₀(z⁴G₄), with U_{4,k} = U_{0,k+4}/U_{0,4}; coordinator-verified on 50 random states.
  * 𝒞₉ = a/U_{0,4}. **Checked here by hand** from the Γ-ratio with N = N₀ − 4.
  * So the s = 9 and reduced s = 1 contour integrands are identical pointwise. **The shift alone does not remove the cancellation.**
  * r⁴(1 − 20/C) ≤ U_{0,4} ≤ r⁴.
* **A⁺ transfer.** G₄⁺ = A/(1−z)³ − [(49/8)A + B/2 + ξA′ + Λ²D]/(T(1−z)⁴), with the Lemma A⁺ majorants. The helicity factor is \|ξ\|⁸/U_{0,4} ≤ 1.000021(\|ξ\|²/r)⁴.
  * At saddles: ≤ 1.000023/(1−κ²)⁴, which is < 3.161 for ν ≤ a/2.
* **The cancellation.** I_U ≈ iI_T near the shared corner γ = π/2 − iε at small A, where \|ξ\|⁸/U_{0,4} → r^{−4} (up to 4.2e5 for r > 0.0392).
* **The fix.**
  * The diagonal zero is of order 8: ξ⁸G₄⁺ = O(δ⁸ log δ), so derivatives of order ≤ 7 vanish at γ = α.
  * Eight boundary-free integrations by parts along the combined soft path give a factor (ia)^{−8}; the corner terms cancel between the arc and the tail.
  * Practical version: Taylor-expand F_α at the corner and integrate e^{−iaγ}P(γ) exactly with the SAME polynomial on both sides, so the shared endpoint cancels exactly.
* **Not yet certified:** the global contour norm relative to N_s.

## 73. High-twist scan complete: all 60 β-cells, both channels (author's laptop run, 26 Nov)

* **Certified a*(L)** (fp3_ht.py with the fixed far remainder; run_laptop_scan.py, 22 workers; summary reported by the author):
  * β ∈ [0.39, 1.04] (L ≈ 0.04–0.25): a* = 3200;
  * β ∈ [1.04, 3.04]: a* = 3600;
  * β ∈ [3.04, 8.24]: a* = 4500;
  * β ∈ [8.24, 8.40] (L ≈ 7): a* = 4000.
* **Consequence.** H_ang,b > 0 for every heavy state with β_c ∈ [0.39, 8.40] and a ≥ 4500 (≥ 3200/3600 in the lower cells), both channels, all J.
* **Remaining.**
  1. The sliver 998 ≤ a < a*(L), J ≥ 1001 (FULLPOS4).
  2. β > 8.4 (FULLPOS5).
  3. The tail.
* **Note.** The per-cell jsonl records live on the author's laptop (cache/laptop_*.jsonl). They should be pushed for the audit trail.

## 74. Diagonal-endpoint resummation (author, 26 Nov)

* **Context (agent FULLPOS4, at (a, J, b) = (2000, 1500, 0)).** Contour norms in units of N_s ≈ \|H\|:
  * fully absolute 2.2e10;
  * soft-signed 2.6e6;
  * hard-localised ~7e5.

  The remaining cancellation is in the diagonal (z → 0) region: G(0) ≈ 20.8 must cancel to G(c) ≈ 3e-5.
* **Exact hard projection (author).**
  * D_m(v) = (1−e^{−v})^{2m}Σ_{r≤J−m}C(2m+r−1, r)e^{−rv}.
  * a∫e^{−av}D_m dv = U_{0,m} (telescoping beta sum). **Checked here** (mpmath, 4 cases including non-integer a): agrees with U_{0,m} to 15 digits.
  * So the boundary-jet expansion resums to the ladder. But finite truncation is ill-conditioned and the 1/a boundary series diverges when r ≥ a.
* **The construction to use (author).** Change the soft variable to ξ. Then dγ/sin α = 2dξ/((1−ξX)(1−ξ/X)), and the carrier becomes X^B((1−ξX)/(1−ξ/X))^a.
  * Extract the real segment 0 ≤ ξ ≤ ξ₀ = 0.02 (β₀ ≈ 0.04), plus a short connector:
    H_diag = (aw³/(πU_{0,b}))∫₀^{ξ₀}ξ^{2b}G_b(ξ²)𝔥_B(ξ)dξ, with 𝔥_B(ξ) = ∫_{−π/2}^{π/2}X^B(1−ξX)^{a−1}(1−ξ/X)^{−a−1}dA.
  * 𝔥_B is real, so the outgoing logarithm cancels and the full real G is used.
  * 41-term closed form: 𝔥_B ≈ (2(−1)^J e^{Bε}/B)·Re Σ_{j≤40}(−1)^j Z^{(j)}(ε)/B^j, with error < (6/B)e^{−78}. Matches direct GL to 1e-16 (author's diagnostic).
  * Certified prescription: 40 panels × 48-node GL; total error < 3e-22.
* **Diagnostic (author, ideal weights).** The signed diagonal is −1.66e-6 at (1500, 2000, 0), i.e. ≈ −0.10·N_s, with real-integral norm 1.4e-2; b = 4 values are ~1e-17.
* **Connector (author).** Carrier ≤ e^{1.002}; with the outgoing exponential, ≤ e^6; \|β\| > 0.03999.
* **Open.**
  1. The norm of the remaining contour (without the diagonal prefix) relative to N_s.
  2. The certified implementation.

## 75. Cost of a certified sliver contour run: does not fit before submission (agent FULLPOS4, 26 Nov)

* **Status after the §74 diagonal extraction.**
  * (2000, 1500, 0): diagonal ≈ 0.106\|G(c)\|, remaining contour norm ≈ 1e4\|H\|.
  * (998, 20000, 0): ~340\|G(c)\| before extraction.
  * The toy end-to-end check of the full §74 pipeline gives 6e-11 relative.
* **Cost is set by node counts, not by the norm.**
  * Hard angle: ~3000 α-nodes (320κ oscillations after rescaling). Soft path: ~4a nodes per α. Total ~1–2e7 amplitude evaluations per state.
  * Certified balls cost ~10–20 ms each, so ~50–100 core-h per state. With a one-time certified table of 𝒫⁺ (not yet built; ~1–2 days), ~0.3–0.5 core-h.
  * Cells: naive ~2e5; with a (τ, 1/S, a) design perhaps 1e3–1e4 (unproven).
* **Total:** optimistic ~1.5e4 core-h, naive ≥ 1e6 core-h, against ~176 core-h available overnight. Raising ξ₀ lowers the norm but not the node counts.
* **Conclusion.** The sliver 998 ≤ a < a*(L), J ≥ 1001, 0.04 < L ≲ 30 cannot be certified by brute-force contour quadrature before submission. It needs either an asymptotic (saddle) evaluation with certified remainders replacing the quadrature, or weeks of compute.

## 76. Lemma S, part 1: transport expansion in log-tan coordinates (author, 26 Nov)

* **Coordinates.** u = log tan(α/2), v = log tan(γ/2), with β = v − u, dα/sin α = du, dγ = sech v dv. Phase Ψ₊ = i[f_S(u) − f_a(v)], with f_m(x) = 2m·arctan(e^x) − νx.
  * Saddle points: u_* = −arcosh(S/ν), v_* = −arcosh(a/ν).
  * Curvatures: f″ = ν√(1−ν²/m²). **The hard saddle's large parameter is ν, not S.**
  * **Checked here by hand:** f_m(u) = mα + ν log cot(α/2); the saddle condition is m sech x = ν; the curvature follows.
* **Transport machinery.**
  * Amplitude A_{p,b}(β) = tanh^{2b}(β/2)cosh^{2p}(β/2)Q_ν(β), with Q = e^{−iνβ}P⁺ obeying Q″ + (2iν + coth β)Q′ + (¼ + iν coth β)Q = 0.
  * Phase-correction polynomials: E_r = (i/r)Σ_j j c_j X^{j+2}E_{r−j}, evaluated with Gaussian moments.
  * No quadrature anywhere.
* **Order test (author)** at J = 1001, a = 998, b = 0, p = 3, against the exact ladder (1050 digits). Absolute errors after w³ normalisation, at orders 0 / 1 / 2 / 3:
  * ν = 20: 2.5e-6 / 1.4e-3 / 1.2e-4 / 1.8e-5;
  * ν = 50: 4.1e-6 / 1.3e-4 / 4.1e-6 / 2.7e-7;
  * ν = 100: 4.3e-6 / 2.9e-5 / 3.0e-7 / 1.7e-8;
  * ν = 300: 1.6e-5 / 3.0e-7 / 1.1e-8 / 1.8e-11;
  * ν = 700: 6.3e-5 / 1.7e-7 / 4.2e-9 / 9.1e-12;
  * ν = 950: 5.6e-4 / 6.0e-5 / 2.9e-6 / 9.4e-7;
  * ν = 993: 2.5e-2 / 1.3e-2 / 4.3e-2 / 1.4e-1, i.e. fails at the turning point.
* **Assessment.** Excellent in the central band. Problems at low ν (the hard parameter ν is small; the Legendre/T0 kernel is accurate there instead) and near κ → 1 (a uniform Airy form is needed).

## 77. Large-β high-twist certificate (agent FULLPOS5 final report, 26 Nov)

* **New module, frozen.** fp3_ht_large.py and run_large_scan.py (fp3_ht.py untouched).
  * Large-β failures fixed: arc-subdivision floor, precision, J-tables, ρ-grid, iteration cap, overflow. None could cause a wrong PASS.
  * New bounds [P + DV]: Lemma B′ (oscillation-aware \|P_ν\|), Lemma A′ (coth correction by parts), residual rest split by frequency support.
* **Probes [certified], width 0.1, both channels:** a* = 3200 at β = 9, 12, 16, 20; 3600 at β = 30, 45, 63.9. Single runs show 3400 passing at 30 and 64; the [8.24, 8.4] cell now gives 3400 instead of 4000.
* **Laptop run.** run_large_scan.py, 556 cells on [8.4, 64], running since commit 8ded31f. The η-choice change in ea2f379 is rigour-neutral (agent-confirmed; regression check agrees to 7 digits).
* **β ≥ 64 is not reachable by this certificate.** The low-frequency residual mass, known only as an absolute integral, gives G^r/F ≈ 8e-6β², which reaches 1 near β ≈ 350. So β → ∞ needs the contour route (ℰ₂ < 0.01𝓜 for L ≥ 30.16), with its connection terms j ≥ 1 and k = 3 still unchecked. --beyond B works to β ≈ 200.

## 78. Lemma S, part 2: quadrature-free rational reference (author, 26 Nov)

* **Saddle side.** Airy-zone phase 𝓕_a(π/2+t) = aπ/2 + (a−ν)t − νt³/6 − νt⁵/24 + E, with \|E\| ≤ 0.014ν\|t\|⁷. A global CFU remainder is not obtained.
* **Replacement reference.** 𝒦_p^ref = 𝒬_M(ϑ)F_{p,ν}(c), with the full rational ladder Taylor-truncated in h at order M.
  * Claimed uniform remainder w^p\|𝒦_p − 𝒦_p^ref\| < e^{−105} (M = 1000; band-dependent M = 64–850 gives e^{−85}–e^{−129}) for a ≥ 998, J ≥ 1001, 0.04 < L < 30, ν ≤ 1100, both channels.
  * Proof: the §63 vertical-line machinery with normalised coefficient majorant 20e^{πν/2}[1+d(k+5)]⁵; 1215 d-cells; outward-padded elementary arithmetic, **not Arb, not reproduced**.
* **Evaluator.**
  * R_{m,j} recurrence in the scaled falling-factorial basis (author: 588 exact-rational checks); W_j = Σ_m R_{m,j}.
  * Contraction against the scaled jet [t^j]Y, Y(t) = w^pF(c + cdt). The jet comes from a degree-3-coefficient ODE recurrence seeded by P_ν(cosh β_c) and ∂_βP_ν.
  * Cost: O(M²) per state, O(M) per frequency.
  * **Internal cancellation:** coefficients reach ~1e366, so ball arithmetic at ~1300 bits is required.
* **Coordinator caution (to be tested).** Evaluating one state is cheap. A certificate over the CONTINUUM of states (real a and J up to ~1e7) needs parameter cells. With ~1e366 internal cancellation in Σ_j W_j[t^j]Y, a naive interval enclosure over a cell is hopeless unless the cells are tiny or Taylor models in (c, h) track the cancellation. That cost decides feasibility.

## 79. Benchmark of the §78 rational reference (agent FULLPOS4, 26 Nov)

* **(a) Pointwise exact [DV, ball arithmetic, 1600–8000 bits].** 𝒦_ref·w³ against the exact ladder, max\|R_{m,j}\| in brackets:
  * (1001, 998) [1e571]: differences 0 at ν = 20/300/700, 1e-271 at 950, 7e-237 at 993;
  * (1500, 2000) [1e187]: differences of 2e-185 or less;
  * (20000, 998) [≈ 1]: ≤ 2e-276.
  * Cost: W_j takes 5–40 s per state, each ν 0.04–0.4 s, so H_ref takes ~15 min per state.
* **(b) Continuum certificate: does not close tonight.**
  * Naive interval cells in a lose the sign for δ ≳ 1e-550; at J = 20000 the R-recurrence wrapping is ~2^M.
  * Taylor-model cells: δ ≲ 1e-6, so ~1e9 cells per J. J must be continuous (n ≈ 2ae^L, up to ~1e12 at L ≈ 20).
* **Needed instead.** An a priori uniform derivative bound for H in (a, 1/J), of physical size (Lemma U).
* **Contour side, for reference (§74 extraction, soft-signed remaining norm / \|G(c)\|).**
  * (2000, 1500, 0): ξ₀ = 0.02, 0.05, 0.1, 0.2, 0.35 gives 1.6e4, 2.9e3, 790, 198, 61.
  * (998, 20000, 0): ξ₀ = 0.02 gives 97 (16 after localisation).
  * The node counts still rule out an overnight run.

## 80. Lemma U: exact derivative contractions; twist direction as one polynomial; c-cell polynomial construction (author, 26 Nov)

* **Coordinates.** u = β_c = 2atanh√c and x = h/d², with x = 1/(a²(1+c+√(4c+s²d²/a²))).
  * \|∂_{log a} log x\| ≤ 2 and \|∂_u log x\| ≤ 1−√c. The physical β differs negligibly (∂_uβ > 0.9998).
  * **Checked here:** d/h = a(2J+1+a) = aS.
* **Twist direction.** At fixed c, R_{m,j}(c, x) = y^m R_{m,j}(c, x_*), with y = x/x_*. So H_M = P(c, y)/D⁴ with P of degree ≤ M+1 in y, and one polynomial covers the whole twist interval.
  * Certify with the signed Chebyshev conversion P̂_ℓ = (2−δ_{ℓ0})Σ_{m≥ℓ}C(2m, m−ℓ)4^{−m}p_m, and the positivity test P̂₀ − Σ_{ℓ≥1}\|P̂_ℓ\| > D_max⁴E_ref.
* **Exact derivatives.** Y_j = d³(cd)^jG^{(j)}(c)/j! satisfies ∂_uY_j = ([j(1−2c) − 3c]Y_j + (j+1)Y_{j+1})/√c. **Checked here by hand** (dc/du = √c·d). Also x∂_xR_{m,j} = mR_{m,j}.
* **Lipschitz targets.** \|U\| ≤ 4F and \|V\| ≤ 2F give C_β = 6.01 and C_a = 4.01. Not yet evaluated.
* **Impact-parameter cells.** \|c−c₀\| ≤ 0.03d₀ (Δβ ≈ 0.05).
  * Degree-6000 Taylor of G about c₀, with sup bound e^{1770} on \|z−c₀\| ≤ 0.9d₀ and remainder < e^{−600} after contraction (row-sum bound e^{4700}).
  * c^M R_{m,j} is a polynomial in c, so Q(c, y) = c^M P̃ is bivariate of degree (≤ 11105, ≤ 851). A Bernstein positivity check of Q − (E_ref + e^{−600})c^M D⁴ certifies a whole cell.
* **Adjoint recurrence.** A_m = (𝓑^m Y^A)₀, with 𝓑 banded (author: 78 exact checks).
* **Diagnostic (ideal weights, J = 1004, a = 1000).** ∂_βH/H = −0.38 / −0.31 and ∂_{log a}H/H = −0.08 / +0.09 for b = 0 / 4.
* **Open.**
  1. Repaired-density jet evaluation.
  2. Running the signed polynomial checks.
  3. Arb reproduction of the §78 reference error.
  4. Cost of the c-cell construction.

## 81. Lemma R: signed low-frequency residual in the β → ∞ tail (author, 26 Nov)

* **Correction.** The repaired weight has r(0) ≈ +9.69e-14 > 0, not 0; EPSW rejected the slope repair. So a two-sided bound is false asymptotically, since I_r ~ r(0)/(2 sinh(β/2)). It is positive, though.
* **Exact subtraction.** ∫sech(πν)P_{−½+iν}(cosh β)dν = 1/(2cosh(β/2)). **Verified here** (mpmath, β = 0.5, 2, 6, 15 digits).
  * r̃ = r − r(0)sech(πν) has r̃(0) = 0, so ∫W = 0, where W(t) = ∫r̃ cos(νt).
  * Mehler then gives \|I_{r̃}\| ≤ (√2/π)(K₀/β)[B + 8B log 2/β + Eβe^{−β/2}], with E < 2.683e-6 and B < 1.901e-3 (EPSW inputs plus the cutoff correction).
* **Result.** (2/Λ)I_r > −0.027F_*(β) for all β ≥ 30, and > −0.015F_cert for β ≥ 64. The extracted endpoint is positive.
  * Finite mass: the extracted endpoint gives r(0)/Λ·w³𝒱_b((1−z)^{−5/2}) > 0 exactly, since the coefficients and V_{b,k} are nonnegative.
  * The zeroth-order dressing is < 2e-9 of the floor.
* **Open.** The higher-order finite-mass transfer for r̃, i.e. w³\|(ℒ_b − ev_c)(G_{r̃} + T^{−1}𝒟_bG_{r̃})\|, and the higher-frequency residual terms whose FULLPOS5 bounds grow like β.

## 82. Lemma R′: finite-mass signed tail transfer (author, 26 Nov)

* **Dressed multiplier.** D_dr(α) = (1 + (B_b − αA_b)/T)D₀(α) + (αA_b/(Td))D₀(α+1), with D₀ = d^αS_b(α) − 1.
* **Beta representation.** S_b(α) = Γ(N+1+α)/(Γ(N+1)Γ(α))·∫₀¹q^{α−1}(1−q)^aK(q)dq, where K(q) = B(a, n+s)^{−1}∫t^{a−1}(1−t)^{n+s−1}(1−qt)^n dt. **Verified here** (mpmath, two states with non-integer N, both channels): agrees with Σ_kU_k(α)_k/k! to 15 digits.
* **Bounds (author).**
  * Low frequencies: \|D_dr\| < 12/a² on the detour, using the telescoping 0 ≤ c^k − U_k ≤ (2hd/c)c^k[(k)₃/3 + ((s+1)/2)(k)₂]. This proves the multiplier input that the contour certificate had assumed.
  * Kernel: \|K(q)\| ≤ 2\|1+q\|^{−a} (direct proof), with a small-disc variant for a > 1e8.
  * Horizontal line: \|D_dr(5/2 + i(t + 9i/4))\| ≤ 4e4·e^{G₉₉₈(t)}, **replacing the unproven 5e12 / D-profile**.
* **Leading connection term.** \|ℰ₀\| < 6e-4·F_* for β ≥ 30, a ≥ 998.
* **All connection terms (author).**
  * Mellin–Barnes: F_{3,ν} = 2Re Σ_{j≤500}BC_j(1−z)^{−5/2−iν+j} + ℛ_MB, with C_j = (½−iν)_j²/(j!(1−2iν)_j).
  * Derivative beta representation for S_b(α−j).
  * Contributions: Σ_{j≥1} < 2e-6·F_*; \|ℛ_MB\| < e^{−2800}F_*; endpoint matching 3e-8.
* **Result (author).** w³\|(ℒ_b − ev_c)(G_{r̃} + T^{−1}𝒟_bG_{r̃})\| < 0.001F_*(β) for β ≥ 30, all a ≥ 998, both channels. Combined with Lemma R: H_r > −0.029F_*.
* **Not covered.** The ideal-weight finite-mass contribution and the k = 3 piece; the same machinery is said to apply. The CONTOUR_BOUNDS Arb tables are used as inputs; the rest is not reproduced.

## 83. Lemma U twist certificate, first numbers: passes at L = 2 in both channels (agent FULLPOS4) [DV]

* **Setup.** c = 1 − e^{−2}, ideal NF02 weights, M = 1000, 1500-bit midpoint arithmetic, integrated-G jets to j ≤ 3003 (8000-node GL in p).
* **Sanity check.** A₀ = d³G(c) = 2·OmLR(β_c) to 1e-9.
* **A_m behaviour.** A_m is large at moderate m (up to 4e17 at m = 100), then decays (2e-253 at m = 1000). Partial sums are stable from M′ ≈ 300.
* **Twist test** (Chebyshev in y over all a ≥ 998):
  * b = 0: \|P̂₁\| = 1.5e-8, \|P̂₂\| = 2e-14, then noise ~1e-17. (P̂₀ − Σ\|P̂_ℓ\|)/P̂₀ = **0.99757, PASS**.
  * b = 4: **0.99756, PASS**.
  * The finite-mass shift is 0.5% at a = 998.
* **Cost.** ~11 min per c-value on one core.
* **Running.** L = 0.5, 7, 0.1, 20.
* **Caveats.**
  1. DV only: the jets need rigorous enclosures.
  2. The c-direction cells are not yet addressed.

## 84. Large-β scan, laptop progress (author run of run_large_scan.py): 191/556 cells

* **Certified, both channels:**
  * β ∈ [8.4, 11.9]: a* = 3100;
  * β ∈ [11.9, 23.6]: a* = 3200;
  * β ∈ [23.6, 27.5]: a* = 3400.
* At the bisected a*, err/F_cert runs from 0.21 to 0.9989. The pass criterion is err < F_cert, so all cells pass, some only just; a* is the smallest passing rung on the twist ladder.
* β ∈ [27.5, 64]: running. Per-cell jsonl files are on the laptop and should be pushed for the audit trail.

## 85. Lemma C: complex-impact-parameter strip (author, 26 Nov)

* **Proved (author).** On 0.37 ≤ Re u ≤ 31.62, \|Im u\| ≤ 0.02, the reference is analytic and \|D_b\|^{−4} < 1.000081.
* **Infinite-twist face, y = 0 (author; 60-digit outward-rounded intervals over 753 floor cells).** \|H_M(u, 0)\| < 810·F_cert(Re u).
  * Key step: rotate the endpoint ray v = re^{iπ/4}, so ∫\|K\|e^{−Λr/√2} ≤ √(2π/(Λ sinh B)), which recovers a factor Λ^{−1/2}.
  * Components: smooth part 20.7, endpoint 790.4, residual 32.6.
* **Cutoff correction (author).** For f = 2Re g with a complex θ(ν−i), use \|f_θ − f\| ≤ 2\|θ−1\|\|g\|, NOT a bound through \|Re g\|. Corrected: ∫e^{0.02ν}\|f_θ − f\| < 1e-7. Any earlier certificate that used the \|Re g\| form for the complex-shifted cutoff must be rechecked.
* **All twists (author).** \|H_M(u, y)\| ≤ e^{10500}F_cert. Too large to be useful.
* **Open.** sup\|H_M(u, y) − H_M(u, 0)\|/F_cert ≤ 190 on the strip, i.e. the signed finite-twist transfer at complex u. This is the binding step for the sliver's impact-parameter cells. The Taylor-order budget (order 20 gives 0.13%) is ready once K ≤ 1000 holds.

## 86. TAIL CLOSED for L ≥ 20, every twist a ≥ 998 (agent TAIL_CERT final report, 26 Nov)

* **Result.** H_ang,b ≥ (1−ρ)𝓜(L) with ρ ≤ 0.066 for every L ≥ 20 (β_c ≥ 21.4), and ρ ≤ 0.68 for L ≥ 18.
  * Holds for all a ≥ 998, all J, both channels, all connection terms, the dressing and k = 3.
  * Arithmetic is Arb or exact; it rests on [P] lemmas written out in TAIL_CERT.md (unreviewed). CONTOUR and EPSW tables are inputs.
* **Budget at L = 20:** j = 0 gives 8.2e-5; local multiplier 3.3e-6; connection terms j ≥ 1 5.6e-10; remainder < 1e-435; **k = 3: 0.065 (binding)**, falling like e^{−L}.
* **New proofs.**
  * Lemma NB: \|w^βS(β) − 1\| ≤ 1.0003(β)₃/(12a²).
  * Theorem P: oscillatory multiplier bound on the leg, 1.8–20 for t ≤ 45. This replaces the assumed 100/a², 300/a² and 5e12 inputs.
  * A three-term recurrence in β for the connection terms.
  * A localisation lemma with the saddle-chosen radius.
  * k = 3 evaluated with exact Taylor shifts; a ball-overestimation bug was fixed.
* **Audits.**
  * Lemma R reproduced: 0.0249 and 0.0144.
  * The §85 θ-inequality was redone correctly: δE ≤ 1.7e-23 and ∫e^{0.02ν}\|f_θ − f\| ≤ 1.3e-14. No conclusion changes.
* **Consequence.** β → ∞ is closed at every twist. The sliver reduces to 998 ≤ a < a*(L), J ≥ 1001, 0.04 < L < 20. The laptop cells for β > 21.4 are now only overlap.
* **Caveats.** Some paper-bounded steps (terms ≲ 1e-30, a Binet bound, contour-shift justifications) are not ball-evaluated. The shifted-density generating-function identity is checked only on a toy weight.

## 87. Lemma C sublemma (transfer ≤ 190) is false at η = 0.02 (author counterexample, 26 Nov)

* **Counterexample.** b = 0, M = 1000, c = 37/1000 + 37i/10000 (u ≈ 0.39002 + 0.01995i), y = 1/20.
  * \|H_M(u, y) − H_M(u, 0)\| > 0.008899, against F_cert = 1.907e-5: **ratio > 463**.
  * Evaluated at 1100 digits from exact spectral moments, with the signed kernel assembled in a Chebyshev frequency coordinate before integration.
  * Errors: EPSW residual 4.5e-4, cutoff 1.6e-4, ladder replacement e^{−190}.
* **Complex ladder truncation (author).** On \|v\| = 0.04/(k(k+1)), \|U_k\| ≤ 0.055^k, which gives an e^{−190} remainder at this point.
* **Consequence.** The split "face K = 810 plus transfer ≤ 190" fails. The full \|H_M\| ≤ K·F with a moderate K may still hold; it needs a cancellation-preserving proof over the complex strip.
* **Options.**
  1. Certify the full signed H_M on the strip directly.
  2. Narrow η to ~0.015.
  3. Accept a larger K (doubling K costs one more Taylor order).
* **Unsolved in every option:** a rigorous sup over a two-dimensional complex region without cells.

## 88. Lemma U twist test passes at L = 7 as well (agent FULLPOS4) [DV]

* **L = 7** (c = 0.99909), ideal weights, M = 1000:
  * b = 0: margin 0.999986 (\|P̂₁\| = 4.8e-12);
  * b = 4: margin 0.999983;
  * A₀ matches 2·OmLR(8.386).
  * A_m peaks at 1e14 (m ≈ 100); the sum is stable from M′ = 300.
  * Cost: 28 min per c.
* **So far:** L = 2 and L = 7 PASS in both channels. L = 0.1 and 0.5 are queued; complex-u measurements at η = 0.02 and 0.01 (L = 2) are running.

## 89. Lemma Cβ: non-amplifying complex beta measure; continuation obstruction; rotated-line complex truncation (author, 26 Nov)

* **Complex beta measure (author, proof given).** On the saddle contour t(v) = Ae^v/(B+Ae^v), ∫\|dμ_{A,B}\| < 1.015 for \|arg A\|, \|arg B\| ≤ 0.1 and Re ≥ 100, uniformly in size.
  * Key identity: Ψ′(v) = (1−e^v)/(1/A + e^v/B). **Checked here by hand.**
  * The §80 strip parameters satisfy \|arg a\|, \|arg n\| < 0.1.
* **Inner kernel.** \|K(q)\| ≤ 1.05·₂F₁(−0.93Re n, 0.93Re a; 0.93Re(N+s); q) ≤ 1.05 on q ∈ [0, 1].
* **Obstruction.** For non-integer complex n the ladder does not terminate, and U_k can have poles: n = 1001.5, N = 2000 lies in the rectangle. Also S_β(−m) = Σ_k(−1)^k C(m, k)U_k. So the continued beta form cannot be identified naively with H_M; a regularisation is needed.
* **Outer Mellin factor.** The 1/Γ(α) factor is ~e^{1550} at ν = 998. It is cancelled exactly on the contour q = it/(1+it), which gives a factor \|i^α\| = e^{−πν/2}. A joint bound on that curved q-path is still open.
* **Complex truncation without integer spin.** On the rotated line v = e^{iφ}(ρ+it), φ = arg(1−c), Θ = (1−\|c\|²)/(2\|1−c\|):
  * \|1−vp_j\|² − \|c−vp_j\|² = 2\|1−c\|(Θ − ρp_j). **Checked here by hand.** Hence \|U_k\| ≤ 1, with an exponential envelope.
  * Residue identity for U_k − T_{M,k}, and a circle maximum (\|c−x²\| + x\|1−c\|)/(1−x²).
* **Status.** Lemma C (full signed strip bound with K ≲ 1000) is open. Route: a regularised beta transform on the saddle contour plus rotated-line truncation; the joint complex q-path estimate and a uniform truncation budget remain.

## 90. Lemma C: finite-prefix regularisation and uniform complex truncation (author, 26 Nov)

* **Finite prefix.** Keep the exact U_k(h; c) for k with w(k+s/2) ≤ 20000 and p_k Re(e^{−iφ}h) ≤ 0.98Θ, where w = \|1−c\| and Θ = cos φ − w/2. Drop the rest. The rotated-line residue identities then keep every retained term away from its poles. **This resolves the §89 continuation obstruction.**
* **Regularised beta representation (author, proof given).** S_L(α) = Σ_{k≤L}U_k(α)_k/k! = Γ(N+1+α)/(Γ(N+1)Γ(α))∫₀¹q^{α−1}J_L(q)dq, with J_L = (1−q)^N Σ_{k≤L}C(n, k)(n+s)_k/(N+s)_k·(q/(1−q))^k. Exact for complex n, N with Re N > L−1; Re N > L+1 holds on the strip.
* **Uniform truncation (author).** \|H₁₀₀₀ − H_pref\| < e^{−90} on the whole strip, both channels, with the repaired weight.
  * It rests on a product bound \|U_k(v; c)\| ≤ exp(62 − XI/√(ϑT)) and the scalar inequality (1+√(1−w))I > 1.91.
  * That inequality was certified on 6526 interval cells (directed floating-point, not Arb), with smallest margin 1.9100191.
  * On real integer-spin states this also controls the §78 reference error (< e^{−200}).
* **Outer contour.** The curved path q = ir/(1+ir) is ruled out: at n = 1e15, a = ν = 998 the absolute integral exceeds e^{120}, because of the geometry, not a loose estimate.
  * Replacement 0 → 0.85i → quarter-circle → 0.85 → 1: the arc exponent is < 1.6 at a = 998, giving 𝒩_L < 500 at that state.
* **Open.** A uniform joint bound on the replacement contour over all complex parameters, including the signed cancellation of the regularised kernel, and hence Lemma C with a practical K.

## 91. Fixed outer contour disproved; signed low-frequency multipliers proved; regularised recurrence (author, 26 Nov)

* **Counterexample** (b = 0, n = 1001, a = 3000, ν = 998, inside the sliver). On the clipped quarter-circle, 𝒩_L > 2^{305} > e^{211}, from exact integer arithmetic (\|K(i/2)\| > 2^{89}).
  * The SIGNED multiplier is tame: \|d^αS_L(α)\|² ∈ (1.01040, 1.01041). **Checked here independently:** 1500-digit ladder sum gives 1.010403764.
* **Unbounded in general.** Along the family n = a = m → ∞ the absolute norm on the fixed contour is unbounded. No fixed-contour modulus proof can work.
* **Signed low-frequency multiplier theorem (author, proof given).** \|d^αS_L(α) − 1\| ≤ e^{0.02\|τ\|}(Γ(σ)/\|Γ(σ+iτ)\|)(C_σ y + e^{−46}) for σ ∈ {5/2, 7/2} on the whole complex strip, with C_{5/2} = 5.9e-5 and C_{7/2} = 1.1e-4. So the error is < 0.004 for \|τ\| ≤ 5. Proof by telescoping plus negative-binomial factorial moments.
* **Regularised three-term recurrence (author).** [D − α(α−s)]S_L(α−1) = [dD + α(s−2α−1)]S_L(α) + α(α+1)S_L(α+1) + (cD − p_L)U_L(α)_L/L!, with D = 1/h.
* **Open.** A uniform SIGNED high-frequency multiplier bound (τ up to ~1100) on the complex strip, and hence Lemma C.

## 92. High-frequency signed multiplier via an exact spin recurrence (author, 26 Nov)

* **Genuine turning-point growth.** At b = 0, n = 8000, a = 998, α = 5/2 + 1100i: \|d^αS(α)\| ∈ (4220, 4221). **Checked here** with a 1300-digit ladder sum: 4220.389.
  * This is consistent with e^{G_a(τ)}, G₉₉₈(1100) ≈ 30. So the target is a weighted bound with e^{G_a(τ)}, not a constant.
* **Exact spin recurrence at fixed N, s, α.**
  * A_n[S(n+1) − S(n)] + C_n[S(n−1) − S(n)] = αS(n) − (α+L)U_L(α)_L/L!, with an explicit regularisation defect.
  * For b = 0: S(n+1) = S(n−1) + (2(2n+1)/(D − n(n+1)))[αS(n) − defect].
  * **Checked here** (full ladder, defect zero): relative error 1e-61 at N = 57, α = 2.5 + 7.3i, n = 5, 12, 30.
* **Real-state bound (author).** \|d^αS\| < 21000e^{G_a(τ)} for all integer n, real a ≥ 998, σ ≤ 3.5, \|τ\| ≤ 1100, and < 400 for \|τ\| ≤ a.
  * Method: a quadratic energy for the transfer matrix; the Σq_j ≤ −log d cancellation; turning-point growth ≤ G_a.
* **Theorem C-HF, complex strip (author).** \|d^αS_{0,L}\| ≤ e^{12 + 0.043\|τ\| + G₉₉₇.₈(τ)} and \|d^αS_{4,L}\| ≤ e^{29 + …}, uniform in spin and including y → 0.
  * The defect is < e^{−120} (scalar integer-arithmetic enclosures over 410 cells). The complex perturbation sum is < 0.023. b = 4 follows by the helicity shift.
  * Not reproduced in Arb.
* **Open.** The constants are too large for K ≲ 1000. Track the physically selected recurrence branch rather than bounding all solutions; combine with the fifth-order endpoint zero and the signed spectral integration; then an audited uniform truncation.

## 93. FULLPOS4 final report: Lemma U passes at all tested L; Lemma C strip width depends on L (agent, 26 Nov) [DV]

* **Lemma U twist test** (ideal NF02 weights, dressing included, M = 1000). Margins, b = 0 / b = 4:
  * L = 0.04: 0.978 / 0.977;
  * L = 0.1: 0.972 / 0.971;
  * L = 0.5: 0.985 / 0.985;
  * L = 2: 0.9976 / 0.9976;
  * L = 7: 0.99999 / 0.99998.

  All pass. \|R_{m,j}\| reaches 1e1697 at L = 0.04, which needs 6000 bits. Cost 10–28 min per c.
* **Lemma C, measured** (max over y of \|H_M\|/F_cert at u = β_c + iη, b = 0 / b = 4):
  * L = 2, η = 0.01 and 0.02: 1.26 / 1.25 (transfer 0.009);
  * L = 0.5, η = 0.02: 1.27 / 1.40;
  * L = 0.04, η = 0.02: 242 / 5271 (transfer 479, confirming the author's > 463).

  So η = 0.02 is fine for L ≳ 0.5. At very low L it needs ~0.003–0.005, but **L ≲ 0.058 is already closed by the per-J certificates** (a < 3200 ⇒ J ≤ 1000). The relevant low end is L ≈ 0.06–0.5, not yet measured off-axis.
* **Contour route.** It works numerically: the §55 benchmark is reproduced and the §74 pipeline matches to 6e-11. Brute-force certification is infeasible, ≳ 1.5e4 core-h.
* **§78 reference.** Exact at single states; plain interval cells fail through ball wrapping.
* **Missing for a sliver certificate.**
  1. Rigorous jet enclosures of the integrated G, with the realised EPSW-repaired weight.
  2. An Arb reproduction of the §78 E_ref.
  3. The c-direction: ~1.4e3 c-points × 15–30 min, i.e. several nights on 22 cores, before certification overhead.

## 94. Full-kernel recurrence and principal-phase control (author, 26 Nov)

* **Exact annihilator for the complete conical kernel.** [z(𝒟−p)² + 𝒟 − p + λ]F_{p,ν} = 0, with 𝒟 = (1−z)∂_z.
  * Combined with the ladder identity: 𝒫(t) = R₀Q(t−1)Q(t) + R₁QP + R₂P(t+1)P(t), degree 4 after exact cancellation (t⁴ coefficient 16 − 4s − C₀).
  * [𝒫(𝒜_n)K_p]_n = 0, a nine-point lattice equation, with explicit boundary terms for the finite prefix.
  * Author: 60 + 60 + 600 + 3900 modular checks.
* **Principal symbol.** m²t⁴ + ν²(M²+m²−ν²)t² + M²ν⁴, with discriminant ν⁴(a²−ν²)(S²−ν²). **Checked here by hand.**
  * Phases θ_n, θ_N = ±arcsin(ν/a) + arcsin(ν/S), so θ_n + θ_N = 2 arcsin(ν/S).
  * This reproduces the §54 separated phase algebraically.
* **Principal complex-phase bound.** \|e^{±iΦ}\| ≤ exp[G_{a₀}(ν) + (1.15η + 0.00032)ν], i.e. e^{G + 0.0234ν} at η = 0.02.
* **Nonstationarity.** Re ∂_νΦ > 0.31 on the strip, so there are no stationary frequencies. This enables signed integration by parts against the fifth-order NF02 zero.
* **Spectral tail.** ν > 1100 contributes < e^{−1200} to the prefix action.
* **Observation.** At u = 0.39, ν = 998, a single connection term \|C₂₀₀\|d^{200} > e^{366}. So the connection expansion cannot be bounded termwise; the full-kernel recurrence avoids it.
* **Open.** A uniform amplitude estimate for the physical solution of the lattice equation, including the Airy transition near ν = a; then the practical K bound for H_M.

## 95. Full-rectangle K = 1000 is false; restrict Lemma C to the physical domain (author, 26 Nov)

* **Counterexample.** b = 4, c = 0.037 + 0.0037i (u ≈ 0.39 + 0.02i), y = 1, M = 1000: \|H₁₀₀₀\| > 0.13, so \|H\|/F_cert > 6700.
  * Method: signed Chebyshev-in-ν² assembly of the spectral kernel at 1000 digits; exact ideal-weight moment recurrence μ_{j,m}; repaired-residual quadrature; local truncation certified < e^{−120} (128 interval cells).
  * Not reproduced here.
* **Why this does not threaten positivity.** The point is off-axis and in a low-L, low-spin region that is already certified.
* **Restriction.** For J ≥ 1001 and a < 3200: c ≥ 1001·1002/(4201·4202), so **u ≥ 0.486**. More generally y ≤ Y_s(u) = min{1, c/(1001(1001+s)d²x_*)}.
  * So Lemma C only needs complex discs around the restricted real domain.
  * The agent's measurement at L = 0.5 (u ≈ 1.2) gives K ≈ 1.3–1.4 at η = 0.02. A larger K costs one Taylor order per doubling (K = 1e5 needs order 24).
* **Open.** The restricted physical-strip inequality (target K ≤ 1e5), then the cell run.

## 96. Independent audit of TAIL_CERT: tail closed for L ≥ 20 confirmed (TAIL_AUDIT.md, 26 Nov)

* **Verdict: holds** for every L ≥ 20, all a ≥ 998 (including non-integer and a → ∞), every spin, both channels, with ρ ≤ 0.066.
  * Every tailcert script was rerun in a mirror; outputs match digit for digit.
  * About 900 adversarial spot checks: worst ratio to the stated bounds 0.9995 (Lemma NB) and 0.95 (Theorem P).
* **Low-severity gaps; each fix changes no number:**
  1. The Lemma NB step needs a finer argument (excess ≤ 2e-5).
  2. Theorem P at non-integer a needs the general form g = E_V[(1−p/V)₊^a] with a partial deformation.
  3. The recurrence range in the doc should say \|β\| ≤ 6200.
  4. The weight bound beyond ν = 1060 should cite CONTOUR's tail formula.
* **Info.** k = 3 is binding only because Theorem P is loose; the true k = 3 term is 4.8e-5·𝓜 at L = 20, so L_tail could be lowered.
* **Inherited, not re-derived:** the action formula and the generating-function identity E5, CONTOUR's Arb tables (spot-checked), and EPSW's r(0) and sup\|r\|.
* **Open region unchanged:** 998 ≤ a < a*(L), J ≥ 1001, 0.06 ≲ L < 20 (u ≥ 0.486).

## 97. Polynomial-certificate route for the sliver (author, 26 Nov)

* **Restricted domain.** 0.466 ≤ Re u ≤ 19.407, \|Im u\| ≤ 0.02.
  * Independent truncation check: 2260 interval cells, bound 1.900054618.
  * Results: \|H₁₀₀₀ − H_pref\| < e^{−85}, and \|H₁₀₀₀ − H_ang\| < e^{−84} on real physical states.
* **Lemma P6000.** Replace G_A, G_D by degree-6000 Taylor polynomials about an anchor c₀. Then \|H₁₀₀₀ − Ĥ₁₀₀₀\| < e^{−400} on \|u−u₀\| ≤ 0.02, from sup bound e^{1650}, ratio 0.08 and row norm e^{13100}.
* **Certificate.** P_b(c, h) = (1−c)³Σ_{k≤6000}T_{1000,k}(h; c)[(1+q_bh)A_k + h(B_k + σ_bD_k)], bivariate of degree (≤ 6003, ≤ 1001).
  * Positivity of all Bernstein coefficients of P_b − 1e-35(1+q_bh)⁴ on real cells (~473 β-cells of width ≤ 0.04, twist handled by the polynomial in h) implies H_ang,b > 0.
  * No small-K lemma is logically needed.
* **Diagnostics.**
  * 57 off-axis ideal-weight states at Im u = 0.019 in the low-u region: max \|H\|/F_cert = 27.8.
  * One full repaired-weight enclosure near u = 0.5 + 0.019i (b = 4): \|H\|/F_cert < 750.
* **Repaired-moment recurrences.** M_{j+1} = (c+i)M_j + jσ²M_{j−1} + σ²E_j, and J_{j+1} = (z₀+i)J_j + M_j. Not implemented.
* **Coordinator note.**
  * Bernstein coefficients of a high-degree polynomial on a real cell are governed by its size on a complex neighbourhood of the cell. So the Bernstein checks succeed only if \|H\| stays moderate on discs of radius ~ cell width; the diagnostics (≤ 750) suggest it does. The certificate stays rigorous regardless, and it simply fails rather than lying.
  * Cost (not yet estimated): certified degree-6000 jets of the integrated G at ~473 anchors, plus degree-(6003, 1001) Bernstein conversions at ~6000-bit precision. Likely days of compute.
* **Status.** A complete, rigorous route on paper; execution outstanding. The sliver is not yet closed.

## 98. First executed sliver certificates at fixed spin (author, 26 Nov; package not yet in the repository)

* **Method (author).** At fixed physical spin, expand E_k(t) = 15^kU_k(D₀ + Wt) in the mass variable t ∈ [−1, 1] with the exact triangular recurrence e_{k+1,m} = [15(C₀ − k(k+s))e_{k,m} − We_{k+1,m−1}]/(D₀ − k(k+s)).
  * Contract against the Chebyshev-in-ν² conical coefficients to get the signed polynomial 𝒦(t, x) = ΣA_{m,ℓ}t^mT_ℓ(x), then bound by Σ\|A\|.
  * 1400-bit fixed-point integers with explicit rounding errors; residual split into three frequency bands near the gap.
* **Claimed certified intervals.** Lower bounds on the normalised action:
  * J = 1002, b = 0: twist [998, 1002] gives 3.57e-5; [2951, 3049] gives 2.40e-5;
  * J = 1005, b = 4: twist [998, 1002] gives 2.98e-5; [2951, 3049] gives 2.26e-5.

  This uses the terminating ladder, the repaired weight, the dressing and k = 3, and relies on the retained action identity and the EPSW inputs.
* **Status.** Author-executed (C++ package of 18 files, not yet in the repository; must be pushed and rerun independently).
  * Coverage is so far only four (J, twist-interval) pairs.
  * Fixed-spin certificates cannot cover J up to ~1e7+; a spin-uniform version (compactified 1/J) is still needed.

## 99. Fixed-spin sliver certificates: package added and reproduced (coordinator)

* **The package** is in code/sliver/: C++/GMP checkers check_wide, check_gap, check_gap_low, check_gap_mid, verify_gap_bands.py, README and the author's logs.
* **Reproduced here.**
  * check_wide 1002 0 rebuilt from source and rerun (45 s): output identical to the author's log, action_lower = 2.4032e-5, PASS.
  * The hard-coded NF02 weight coefficients (X, Y, k₀) match prl/data/nf02_psi.json exactly (psi2 and psi3, exact rationals).
  * All remaining runs (wide b = 4; gap, gap_low, gap_mid for both spins) rebuilt and rerun; logs identical to the author's. verify_gap_bands.py on the fresh logs: both near-gap certificates PASS (3.5714e-5 and 2.9787e-5).
* **Scope.** Two spins (J = 1002, b = 0; J = 1005, b = 4) × twist intervals [998, 1002] and [2951, 3049]. Not a continuum in spin.

## 100. Analytic sliver attempt: L ≥ 5 closed analytically, pending review (agent, 26–27 Nov)

* **Claim** (ANALYTIC_SLIVER.md §0, §7.3–7.4, §9). For Δ_g = 1000, a ≥ 998 (every twist, not only a < a*(L)), J ≥ 1001, b ∈ {0, 4} and L ≥ 5, d³ℒ(G) ≥ M − (E₂ + E₃ + E_r + E_ρ) with the errors ≤ ρ·F_cert(β_c), ρ < 1. Together with M ≥ F_cert(β_c) (FULLPOS3 §4.1(c)) this gives H_ang,b > 0.
* **Coverage and ρ.**
  * L ∈ [5, 8]: Route T′ with the new Lemma Loc′ (J₀ = 1500), half-unit cells: ρ = 0.947, 0.179, 0.070, 0.042, 0.032, 0.028. The worst cell is [5, 5.5].
  * L ∈ [8, 16]: Route T′ with TAIL's Lemma Loc, ρ ≤ 0.0464.
  * L ≥ 15.5: TAIL+E′, ρ ≤ 0.50 at 15.5 and ≤ 0.171 at 16.
  * L ≥ 20: TAIL_CERT, ρ ≤ 0.066.
* **New ingredient: Lemma Loc′.** A regionwise remainder bound replacing TAIL's single global localisation constant: connection-series tail, B-tiers weighted by (1.02d/b)^998, an angular fraction, dressing factor 41, and the B = E/(2iν) pole handled by Cauchy. log∫|f||𝒱_bℒ(R_1500)| = −894 at L = 5. Spot checks: K₀ = 1.0003 against the bound 1.6; the j-sum growth constant is 1.946 (≤ 2 asserted).
* **Unreviewed.** Theorem E′, Lemma X, the helicity bound, the splitting lemma, Theorem R, Lemma NB with explicit d, Lemma V, Lemma Loc′, and the state relations at L ≥ 4.5 (0.49d ≤ a/N ≤ 0.51d and Td² ≥ 3.9a²; exact algebra, grid spot-checked only).
* **Below L = 5 (open).**
  1. Turning-point residual at a = 998, ν ∈ [1031, 1032]: 0.94·F_cert at L = 5 and 16·F_cert at L = 4.5. This is termwise j-series growth on top of Lemma X; it needs a resummed bound on K_J near ν ≈ a or finer certified residual masses.
  2. Lemma Loc′ fails near L ≈ 4.
  3. Top line (Lemma T) termwise bounds die near L ≈ 2–4, so a full-kernel upper-half-plane bound is needed.
* **Sliver status.** If reviewed, the open sliver shrinks to 998 ≤ a < a*(L), J ≥ 1001 and 0.06 ≲ L < 5 (plus the four fixed-spin intervals of §98–99 inside it).
* **Code.** code/analytic/ (rt_cert.py, whose 5th argument selects Loc′; rt_loc3.py, rt_pp.py, rt_rmass2.py, rt_rem.py, crossing.py, thmE.py, energy*.py, tp/); logs/. cache/ is gitignored.

## 101. Independent review of the L ≥ 5 sliver certificate (referee agent, 27 Nov)

* **Verdict.** No error found that breaks the L ≥ 5 certificate. Full report: REVIEW_SLIVER_L5.md; scripts and logs in code/review_l5/.
* **Reproduced from source.**
  * The Loc′ cells [5, 5.5], [5.5, 6.5] and [6.5, 8] give ρ = 0.94683, 0.17942, 0.069807, 0.041726, 0.031817 and 0.027767.
  * All 16 cells of [8, 16] are identical to the cached output.
  * The Loc′ integrals and rt_rmass2.json are identical.
  * Not rerun: [4.5, 5], rt_rem.py and TAIL+E′.
* **Cell [5, 5.5] survives (5.3% margin).** Every piece is evaluated at the largest d in the cell and is increasing in d.
* **Lemmas.**
  * Loc′, Theorem E′, the helicity bound, Lemma V, the splitting lemma and Theorem R: hold.
  * Imported constants are used in the safe direction.
  * State relations at L ≥ 4.5 are now proven from the exact formula d = a(2N−a+s)/(N(N+s)): 0.5d ≤ a/N ≤ 0.5029d, Dd² ≥ 3.977a², T ≥ D and N ≥ 179a.
* **Rigor gaps to write out (numerically negligible).**
  * G1: E_r omits the θ-cutoff correction ∫\|f_θ − f\|\|K\| (≈ 1e-11·F).
  * G2: the residual on ν ∈ [1150, 3000] uses a crude kernel bound with j ≥ 3 terms uncontrolled, and the integrals add nothing above ν = 1300 (or above 3000 for the residual); the tail mass is e⁻⁹¹⁰³.
  * G3: Lemma NB with explicit d drops lower-order factorial-moment terms (relative effect ≤ 1e-7).
  * G4: Lemma X not re-derived; 134 numerical checks, exact/bound ≤ 2.5e-5.
  * Loc′ also silently needs F/den ≤ 1 (holds: 0.93).
* **Text correction (E1).** The 0.94·F at L = 5 is the whole band ν ∈ [898, ∞), not the single cell [1031, 1032], which contributes 0.22. The bulk is spread over ν ≈ [1000, 1033].
* **Status.** L ≥ 5 is certified, reviewed, and reproduced, modulo writing out G1–G4. The open sliver is 998 ≤ a < a*(L), J ≥ 1001, 0.06 ≲ L < 5.

## 102. Analytic sliver: L ∈ [3.65, 5] closed (agent, 27 Nov; pending review)

* **Claim** (ANALYTIC_SLIVER2.md). H_ang,b > 0 for every sliver state with L ≥ 3.65 (β_c ≥ 5.02): all a ≥ 998, J ≥ 1001, b ∈ {0, 4}. Together with §100–101 this covers L ≥ 3.65.
* **ρ per cell.** In every cell E_ρ ≤ e⁻²⁸ against F ≈ e⁻¹³.

  | L-cell | ρ |
  |---|---|
  | [4.5, 5] | 0.0267 |
  | [4, 4.5] | 0.0894 |
  | [3.75, 4] | 0.223 |
  | [3.65, 3.75] | 0.545 |
* **New and unreviewed ingredients.**
  1. Lemma F (turning point): Theorem E′ stopped K = 200 steps before the end; the last steps are bounded with Riccati frames affine in (τ, a, 1/N, σ) plus a second-order defect bound. It is sharp to e^{2–4}, against e^{10–11} for Lemma X, and lies above the exact value at 30 random states.
  2. Lemma B (β-chain over j): chained similarity transforms, minimised over frame families. It captures the true decay of about 0.52 per step past the turning point (the old absolute-value chain grew by 1.22 per step), and is within e^{0.4–1.8} of exact.
  3. Lemma Loc″: a proved tail ratio \|c_{j+1}/c_j\| < 1 for j ≥ J, with J′ = 2200 for ν ≥ 999.
  4. Per-unit-cell EPSW residual masses on [1, 898); totals match EPSW's band totals.
  5. Top line raised to Y = 6.25 below L = 4; state relations checked for L ≥ 3.75.
* **Referee items handled.** G1 is included. G2 is replaced by a proved full-kernel bound with an e^{πν/2}-weighted CONTOUR tail. F/den ≤ 1 is asserted. E1 is confirmed.
* **Still open: 0.06 ≲ L < 3.65.**
  1. Loc″ remainder now binds: log E_ρ = +26 on [3.58, 3.65]. It uses the worst case e^{πν/2} on the Cauchy circle; going lower needs a phase-keeping (steepest-descent) remainder bound or the remainder-free κ^M representation.
  2. The contour right edge (0.32·F at 3.65) and the top-line j-sum tail (assertion fails near L ≈ 3.58) still use the old code; they should be upgraded with Lemmas F and B.
  3. The contour truncation J = 1500 runs out near L ≈ 3.5 (growth constant 1.983 against 1.99).
  4. Below L ≈ 3, termwise e^{νd/2} growth at large twist needs resummation (₂F₁ form plus corrections).
* **Code.** code/analytic2/:
  * lemmas: tpf3.py, tpframes.py, bchain.py, loc4.py;
  * pipeline: rt2_*.py, rmass3.py;
  * checks: check_F.py, check_B.py, check_state.py;
  * outputs: cache/rt2_cert_*.txt.

## 103. Large-β high-twist scan completed on the author's laptop (27 Nov)

* **Run.** run_large_scan.py (fp3_ht_large.py), 556 cells of width 0.1 covering β_c ∈ [8.4, 64] contiguously, both channels b ∈ {0, 4}. Merged record: code/fullpos3/results/ht_large_done.jsonl; summary in cache/large_summary.txt. The earlier small-β laptop cells (cache/laptop_*.jsonl, 45 shards) are added for the audit trail.
* **Checked here.**
  * All 556 × 2 records have ok = true and a finite a*.
  * ok agrees with err < F_cert in every record.
  * The cells tile [8.4, 64] with no gaps or overlaps.
* **Thresholds.** a* = 3100 on [8.4, 11.9] (35 cells), 3200 on [11.9, 23.6] (117 cells), 3400 on [23.6, 64] (404 cells).
* **Margins.** Maximum err/F = 0.99888, at β ∈ [23.5, 23.6]. Four records exceed 0.99: [11.8, 11.9] and [23.5, 23.6] for both b. At the next cell the threshold steps up (3200, respectively 3400) and err/F drops to 0.39, respectively 0.21. Rerunning these two cells at a* + 100 would widen the margin; this is optional.
* **Status.** High twist a ≥ a*(β) is certified on β ∈ [8.4, 64], at every spin, in both channels. Together with the earlier small-β scan (60 cells, β ∈ [0.39, 8.4]) and TAIL_CERT, this closes the a ≥ a* side.

## 104. Phase-keeping localisation remainder certified on L ∈ [0.8, 3.65] (author, 27 Nov; checker not yet in repo)

* **Claim.** E_ρ,2 + E_ρ,3 < e⁻⁴⁶ for 0.8 ≤ L ≤ 3.65 (t = a/N ∈ (0.012, 0.26)). It holds for both helicity channels, all spins, all real a ≥ 998, the exact dressing and the repaired density, with J₀ = 1500 kept. This removes the Loc″ obstruction of §102.
* **Ingredients.**
  1. **Termwise exact dressing** (no enlargement): \|(1+τ𝒟_b)(1−z)^{μ_j}\| ≤ q^{j−p+½}e^{νA}D_j, with D_j = 1 + τ[rM_j(M_j+1)/(2q) + (4.5r+0.5)M_j/q + 8] and M_j = j+p+ν. The geometric tail is resummed with S₀, S₁, S₂ (ratio ≤ 1 for j ≥ 1500, ν ≤ 1150).
  2. **Phase from the coordinator's Lemma PK.** A = \|arg(1−z)\| = arctan(2rx√(1−x²)/(1−r+2rx²)) and q² = (1−r)² + 4rx², with x = sin(\|θ\|/2), are used on the tail. The full conical function uses \|Im β(z)\| = arctan(√(2(\|z\|−Re z))/(1−\|z\|)), with Cauchy radius δ = t_lo/100 (τ/δ² < 0.013).
  3. **Two-sided beta Chernoff bound**: P(V ≷ v) ≤ e^{−aI(v;m)}, where m = a/(N+s) and I = log(m/v) + (1/m − 1)log((1−m)/(1−v)). It is combined cellwise with the angular mass, giving E_C ≤ P_B·P_W·min{M_tail, M_{F+T}}.
  4. **Coefficient bounds.** The ratio \|c_{j+1}/c_j\| is maximised over ν ∈ [0, V] at an endpoint (quasi-convex in ν²); log\|c_j\| ≤ 547 / 570 / 630 for j ≤ 1500 at V = 999 / 1040 / 1150.
* **Execution.**
  * 190,031 joint (t, v, θ, ν) cells in directed fixed-point integer arithmetic (scale 1e16), with 0 failures.
  * Worst weighted cell bound: −93.0 on t ∈ [0.012, 0.12], −104.8 on [0.12, 0.20], and −177.3 / −129.6 / −2.27 on [0.20, 0.26] for ν ≤ 999 / 1040 / 1150.
  * Kernel remainders on the last t-range: e⁻¹⁶⁸, e⁻¹²⁰ and e⁷. The last band is multiplied by ∫_{1040}^∞\|f_θ\| < e⁻⁵⁵ (Gaussian inventory Σ\|x_j\| < 1.22e6, centres < 1029, widths < 0.901). Total < e⁻⁴⁶.
  * The F₄/G_D piece is covered by the common q⁻⁴ envelope (e₀ = 7/2), with τΛ² ≤ 1 and ∫\|ψ₃\| < 3e7.
* **Checked here** (numerical spot checks, not a reproduction of the certificate):
  * the q², A and \|Im β\| identities (errors 1e-16 and 1e-199);
  * the beta Chernoff bound against scipy's exact tails (400 states × 10 points; worst slack e⁻¹·⁶);
  * the dressing formula and the D_j envelope constants (\|9/2 − b(b−1)\| ≤ 8 and \|(7/2−2b)z + ½\| ≤ 4.5r + 0.5 for b ∈ {0, 4});
  * the S₀, S₁, S₂ resummation;
  * endpoint maximality of the coefficient ratio;
  * the L→t translation: t = 0.2580 at L = 0.8 and 0.01308 at L = 3.65.
* **Status.** Author-executed; the standalone checker still needs to be committed and independently reproduced. This settles the localisation remainder only. **Still open:** the main signed connection sum on [0.8, 3.65] (contour right edge, top line, termwise hump below L ≈ 3), and everything below L = 0.8 down to the certified corner L ≲ 0.06.

## 105. Loc-Joint certificate (§104) reproduced (coordinator, 27 Nov)

* **Package** (author): code/locjoint/ contains check_joint.js (Node, BigInt, scale 1e16), verify_cell_list.py, verify_spectral_inputs.py, README (with the analytic bridge), run.log and summary.json. The 70 MB cell CSVs are not committed; they regenerate in 6 s.
* **Reproduced here.**
  * The source was read in full before running (it writes only to its output directory).
  * `node check_joint.js` on Node v22.23 gives PASS, 190,031 cells. cells_A/B/C.csv are byte-identical to the author's archive, and summary.json is identical apart from elapsed time.
  * verify_cell_list.py on the fresh output: PASS (complete partition, maxima agree).
  * verify_spectral_inputs.py on this checkout: PASS. Repaired-inventory coefficient ceiling 1,212,961; ψ₂/k₀ and ψ₃/k₀ l1 ceilings 9.16e6 and 1.89e7.
  * Log: code/locjoint/rerun_coordinator.log.
* **Source audit** (coordinator reading, not a second implementation):
  * The beta upper and lower tails take safe endpoints (I decreasing in m above the mean and increasing below it; m ∈ [t_lo/(1+9t_hi/998), t_hi]), with a = 998 in the exponent.
  * The angular tail uses x_lo and t_hi in the decay and t_lo in the prefactor; it is monotone in a once k ≥ 1/998.
  * The bounds on q, A and \|Im β\| are outward.
  * The coefficient ratio is bounded at its endpoints and the cumulative log includes c₁₅₀₀. The tail-ratio inequality is checked exactly at X = 1501.
  * The tail budgets e¹⁵ and e²⁰ hold (≈ e^{10.6} and e^{12.9}).
  * log(a+b) ≤ 1 + max is used for F + T_J; the cell sum uses log(max cells per slab) + log 1.6 for K₀.
* **Retained (not reproduced).** The angular-action identity, the EPSW residual and cutoff integral certificates, the small-ν paired-branch constant 208, and the K₀ < 1.6 Stirling argument (both written out in the README, not independently re-derived).
* **Status.** The localisation remainder E_ρ,2 + E_ρ,3 < e⁻⁴⁶ on L ∈ [0.8, 3.65] is reproduced. The main signed connection sum on this range is still open.

## 106. Signed connection-transfer lemma; low-frequency residual theorem (author, 27 Nov; code not yet in repo)

* **Signed transfer.** Let X_j = d^{α−j}S(α−j) and Z_j = c_j d^j X_j, with r_j = c_{j+1}/c_j and D_j = 1 − hβ_j(β_j − s). Then exactly

  D_jZ_{j+1} = r_j[d + hβ_j(s−2β_j−1)]Z_j + r_jr_{j−1}hβ_j(β_j+1)Z_{j−1},

  so Σ_{j<1500}Z_j = P(d,h,ν)X₀ + Q(d,h,ν)X₋₁. The whole varying recurrence is resummed: this is the drift resummation that ANALYTIC_SLIVER3 §6 left open, done by Taylor models instead of a matrix drift calculus.
* **Certificate method.** Taylor models in (d, h, ν) over a box (3 normalised variables), 65-digit integer arithmetic, total degree m, Cauchy remainder on a polydisk of radius R = 32 against an absolute recurrence majorant M_R. The dressing is contracted in, and the omitted indices are bounded separately. The method is independent of spin.
* **Executed boxes.**
  * L ≈ 3.2, a ≈ 1061 (d ∈ [0.04072, 0.04080], ν ∈ [899.75, 900.25]): \|P\| ≤ 1.94–1.97 and \|Q\| ≤ 0.67–0.68 for dressed F₃ and F₄, b = 0 and 4. The termwise transfer is > 2e9. Taylor remainder < 2.2e-8.
  * L ≈ 2.5, a ≈ 1500, ν ≈ 800, b = 4: \|P\| ≤ 1.28–1.29 and \|Q\| ≤ 0.17–0.18. Remainder < 1.2e-7.
* **Low-frequency residual theorem.** d³\|𝒱_b F_{3,ν}\| ≤ 14 exp(0.0065ν + ν²/4200) for ν ≤ 200 and 0.8 ≤ L ≤ 3.65, at all spins and a ≥ 998. It uses no connection expansion: beta–Cauchy with 1−r ≥ 1.555t on v ≥ 0.75t, \|Im β\| ≤ arctan(X/1.555), P(v < 0.75t) < e⁻³⁶, Cauchy radius t/200 for the dressing, and log(1+X²) ≥ (7/8)X². Against the resid3 per-cell masses this gives E_{r,[0,160]} < 0.08·F_cert, cutoff correction included.
* **Remaining for 0.8 < L < 3.65.** E_ψ₂ + E_ψ₃ + E_{r,ν>160} < 0.92·F_cert − e⁻⁴⁶, needing continuum coverage of the signed transfer in (d, h, ν) and assembly of the contour pieces.
* **Status.** Author-executed on two boxes; the code has not yet been received or reproduced. The low-frequency theorem has not been independently checked.

## 107. Signed-transfer checker reproduced; full-twist boxes; majorant scaling in L (coordinator, 27 Nov)

* **Package** (author): code/signed_transfer/ contains check_transfer.js (BigInt Taylor models in (d, h, ν), scale 1e80), check_identity.py, check_lowfreq.py, make_grid.py, benchmarks, logs and README. The fresh reruns are in code/signed_transfer/rerun_coordinator/.
* **Reproduced here (identical output).**
  * check_identity.py: 576 exact recurrence and dressing checks, 48 signed transfer sums.
  * check_lowfreq.py: constant 12.61 + 7e-4 + 2e-39 < 14. E_{r,[0,160]} < 1.451e-7 = 0.0776·F_cert. The successive rows reach 0.30·F_cert and 0.96·F_cert at the next two ceilings, so 160 is the right cut. PASS, conditional on the EPSW masses and F_cert ≥ 1.87e-6.
  * pilot: 6 channel boxes at L ≈ 3.2 and 2.5, e.g. \|P\| < 1.9358, \|Q\| < 0.6709.
  * combo: full twist, 4% in d, ν = 900 ± 10, b = 4: \|P\| < 5.4653, \|Q\| < 2.0611, remainder 8.96e-5 (27 s).
  * lower40: full twist at L ≈ 2.5, m = 40: \|P\| < 2.9384, \|Q\| < 0.7456 (105 s).
  * complex_top: ν = 900 + 6.25i ± 5: \|P\| < 2.6761, \|Q\| < 1.0107.
* **Source audit** (coordinator reading):
  * The cleared-denominator recurrence (A = (d + hβ(s−1−2β))·nr_j·dr_{j−1}, B = hβ(β+1)·nr_j·nr_{j−1}, D = den·dr_j·dr_{j−1}) is correct.
  * The l1 algebra is submultiplicative for complex coefficients. Division uses a Neumann inverse with outward residual bookkeeping; truncated products are exact.
  * The polydisk majorant covers both basis solutions and the dressing. The multivariate Cauchy tail and the R = 1 tail for indices K..1499 are correct.
  * The h-split is an exact binomial shift.
  * h_cap = d_max²/(998²(1+√(1−d_max))²) is rounded up, so h ∈ [0, h_cap] covers every a ≥ 998.
* **Majorant scaling** (coordinator preflight probe, full-twist boxes, b = 4, rerun_coordinator/majorant_probe_boxes.json). The absolute polydisk majorant M_R is the termwise hump and is almost independent of box size. The Taylor order needed for a Cauchy remainder < 1e-3 at the best R ∈ {2, 4, 8}:

  | L | ν = 500 | ν = 900 |
  |---|---|---|
  | 3.2 | 12–14 | 21–23 |
  | 2.5 | 20–23 | 39–44 |
  | 2.0 | 30–36 | 65–73 |
  | 1.6 | 44–53 | 104–118 |
  | 1.3 | 61–75 | 159–176 |
  | 1.0 | 107–119 | 337–410 |
  | 0.8 | 235–252 | 751–836 |

  Cost per box scales roughly as m³·K (27 s at m = 28, 105 s at m = 40), and memory as C(m+3, 3) coefficients. **So the method as built is practical for L ≳ 1.8**. Below about L ≈ 1.5 it needs a resummed (non-termwise) polydisk majorant, or the region must be handed to the small-L Lemma U route.
* **Still needed for 1.8 ≲ L ≤ 3.65:**
  1. Certified bounds M₀, M₁ on X₀ and X₋₁ per box (Theorem E′, Lemma F, helicity).
  2. An adaptive grid in all 4 channels.
  3. Contour assembly for ψ₂/ψ₃ (top line at complex ν, right edge).
  4. The real-axis residual for ν > 160 against the EPSW cell masses.

## 108. No small complex-twist majorant at low L; real-axis Chebyshev norming replaces it (author, 27 Nov; checker not yet in repo)

* **Obstruction (author, rigorous).** Dressed b = 4 transfer at d = 0.449329 (L ≈ 0.8), ν = 900, with h = h_c(1+x) and h_c = 3.34e-8. The full signed 1500-term recurrence at 400 digits with rounding enclosures gives:
  * \|P\| = 1.736 at the real twist midpoint;
  * \|P(x = i)\| > 1e12;
  * \|P(x = 1.2+1.6i)\| > 1e26.

  So no O(1) resummed majorant can hold on a complex twist disk at low L. This is a property of the true function, not of the termwise majorant, and it rules out the coordinator's suggested target in §107.
* **Replacement (author).** In the twist direction use Chebyshev approximation on a Bernstein ellipse E_ρ: ‖P − P_m‖ ≤ 2M_ρρ^{−m} for ρ ≥ 2. Then the cosine-grid norming inequality (Ehlich–Zeller type) gives ‖P‖_∞ ≤ (max_k \|P(cos kπ/N)\| + ε_m)/cos(mπ/(2N)) + ε_m. This needs only scalar recurrence evaluations at real twists, plus a sup bound M_ρ on the ellipse.
* **Executed (author)**, full twist 0 ≤ h ≤ 6.68e-8, L ≈ 0.8, dressed b = 4:
  * ν = 900: \|P\| < 7.745 and \|Q\| < 2.509 from 513 evaluations (ε_m = 7.4e-7);
  * ν = 1031: \|P\| < 6131 and \|Q\| < 2365 from 598 evaluations (ε_m = 4.6e-7).

  These are single (d, ν) cross-sections, continuous in twist only.
* **Open.**
  1. Continuum in (d, ν): a tensor Chebyshev/ellipse version, or Taylor in (d, ν) combined with Chebyshev in h.
  2. Certified M₀ and M₁.
  3. A kernel budget near the turning point, where \|P\| ~ 6e3 at L = 0.8.
  4. Contour assembly.
  5. The checker has not been received or reproduced.

## 109. Norming checker, base-multiplier lemma, first end-to-end cell [3.55, 3.65] (author; reproduced by coordinator, 27 Nov)

* **Package** (author): code/norming/ contains:
  * scalar.cpp: a GMP directed scalar recurrence with exact rational (d, h, ν);
  * check_norm.py: sequential one-variable Chebyshev norming in d, h, ν;
  * base_bounds.py: the new M₀, M₁ lemma;
  * assemble.py and assemble_top.py;
  * identity, input and low-frequency checkers;
  * logs.

  The coordinator's reruns are in code/norming/rerun_coordinator/.
* **Norming theorem.** Each axis in turn is bounded by an ellipse majorant M (uniform over the other real coordinates; partial sum j < K), giving ‖f − p_m‖ ≤ 2Mρ^{−m}/(ρ−1). The node spacing satisfies m(π/2N + 3e-6) < π/3, so the Bernstein–Szegő inequality gives ‖p_m‖ ≤ 2 max over nodes, hence ‖f‖ ≤ 2 max\|f(x_j)\| + 3ε per axis. Rational nodes have error < 1e-12; the real tail j ≥ K is added. No complex polydisk in all three variables is needed, which avoids the §108 obstruction.
* **Base-multiplier lemma** (spin-uniform, crosses the turning point):
  * Use the exact dual-Hahn recurrence S_{k+1} = S_{k−1} + 2αη(k)S_k, with eigenframe T = [[λ₊, λ₋], [1, 1]], λ_± = w ± √(1+w²) and w = (κ+i)y.
  * Exact identity: ‖T⁻¹T′‖_∞ = (√(1+κ²)/2)(1/\|w+i\| + 1/\|w−i\|). Its integral is A(Y, κ) = ½[asinh(((1+κ²)Y+1)/κ) + asinh(((1+κ²)Y−1)/κ)], finite through y = 1 whenever σ > 0.
  * The growth normalised by d^σ uses ∫₀ⁿη = −log d and η′ ≥ η².
  * b = 4 follows from the exact helicity identity.
* **End-to-end cell, L ∈ [3.55, 3.65]**, both b, all twists. Contour pieces are imported from code/analytic3/cache/contour3_3.55_3.65_k4.json for a ≤ 4500; a ≥ 4500 is covered by §73. The real residual (ν ≥ 159) uses closure5: 5 full-twist boxes × 2 channels, 137,082 scalar evaluations, 881 EPSW mass rows in 16 bands. The rest of the budget is low frequency 1.46e-7, cutoff, Gaussian tail and far tail < 1e-100, and localisation < 1.1e-20.

  | | total error | F_lower | ratio | action lower bound |
  |---|---|---|---|---|
  | b = 0 | 2.036e-7 | 1.98e-6 | 0.1028 | > 1.7763e-6 |
  | b = 4 | 2.142e-7 | 1.98e-6 | 0.1082 | > 1.7657e-6 |
* **Reproduced here (identical exact rational bounds).**
  * `JOBS=3 ./reproduce.sh` (126 s): cross-sections 2/2, tensor1 4/4, complex pilot 4/4 and closure5 10/10 give identical P_bound/Q_bound strings. The assembly gives identical totals and ratios, PASS_CONDITIONAL.
  * Identity checks: 576 recurrence and dressing identities plus 48 transfer sums; 132 spin steps plus 48 helicity identities.
  * verify_inputs.py on this checkout: PASS (881 mass rows, 16 band ceilings, contour ceilings, NF02 exact coefficients).
* **Source audit** (coordinator reading):
  * The coefficients of scalar.cpp (2wβ scaling, r_j = rn/rd, dressing c₀ = 8w²(q + A₀), c₁ = 8w²A₁) are exact, with outward ball errors.
  * The majorant terms in check_norm.py (centre plus variation, pole gaps) are correct, and the sequential norming chain is valid in any axis order.
  * The eigenframe identity was checked numerically to 2e-17 relative, and the antiderivative and d^σ normalisation were checked by hand.
* **Retained.** The contour3 pieces (agent-certified, unreviewed), F_cert (LRINF), the EPSW masses and cutoff, the angular-action identity, Loc-Joint (§105) and §73.
* **Cost note.** At L ≈ 1.6 a full-twist box (0.2% in d, ν width 1) needs 209,760 evaluations and 429 s, so the lower part of the band is expensive.

## 110. Sliver positivity certified on L ∈ [2.9, 3.65] (author; recomputed by coordinator, 27 Nov)

* **Package** (author): code/sliver_ext/ (91 files). The coordinator's recomputation is in code/sliver_ext/rerun_coordinator/.
* **New ingredients.**
  1. **Coupled twist coordinate** h = d²x, with 0 ≤ x ≤ 1/(998²(1+√(1−d))²). It enters check_norm.py through the majorant variation \|d²x − d₀²x_c\| ≤ (d₀+d_r)²x_r + x_c(2d₀d_r + d_r²). This excludes non-physical (d, h) pairs; for example, the dressed b = 4 bound on [850, 950.1] drops from 17584 to 25.1.
  2. **major_rect.py**: a rectangle-enclosure majorant with one complex axis at a time. Re β(β+k) is enclosed exactly, using the vertex. The bound on \|r_j\| is quasi-convex in t², and the pole gap is checked through Re D > 0. The code takes the minimum with the disk majorant.
  3. **ellipse_extend.py**: Hadamard three circles on the Joukowski annulus at fixed real (d, x), giving sup_{E_r}\|P\| ≤ U^{1−λ}M^λ with λ = log r/log R. Here U is the certified real-ν bound and M an outer-ellipse majorant of the full 1500-term sum; Q is handled in the same way. This moves the real transfer bounds onto the top and right contour edges, with an exact containing-ellipse corner check.
  4. **contour_direct.py**: bottom and left edges via the exact difference recurrence for Δ_j = X_j − 1. The source term is gβc(βc + 1 + sd) with g = h/d² and c = 1−d, which the coordinator re-derived by hand. The recurrence is seeded by a negative-binomial telescoping bound (`nb_delta`), and the infinite-twist tail uses the ratio bound. Gamma-duplication and Euler-integral envelopes handle B(ν) and \|B H(d)\| on the top line. The right endpoint 998.5 > Λ is handled with an explicit extra segment.
  5. **Sharper Gaussian tail**: 5 groups of centres with exact amplitude sums and prefactor 1062; the mass above 1040 is < 1.26e-37.
* **Result.** PASS_CONDITIONAL in all four cells, both b, every a ≥ 998, with no spin enumeration. The [3.5, 3.65] contour is imported from contour3 for a ≤ 4500 (with §73 for larger a); the contours below 3.5 are new. Error/floor ratios, author's run vs coordinator's recomputation:

  | L-cell | floor | b = 0 author / rerun | b = 4 author / rerun |
  |---|---|---|---|
  | [2.9, 3.2] | 2.5e-6 | 0.07923 / 0.07923 | 0.08840 / 0.08840 |
  | [3.2, 3.35] | 1.87e-6 | 0.24606 / 0.24425 | 0.29755 / 0.29523 |
  | [3.35, 3.5] | 1.87e-6 | 0.16816 / 0.16812 | 0.19135 / 0.19130 |
  | [3.5, 3.65] | 1.98e-6 | 0.13750 / 0.13751 | 0.14886 / 0.14887 |
* **Reproduction.**
  * `JOBS=3 bash reproduce.sh` takes 616 s. The [2.9, 3.2] outputs (residuals, F₄, all four contour files, budget) are identical to the author's apart from timing.
  * The three older cells were recomputed with the current majorants, which choose different grids, as the README warns. All pass, with equal or slightly smaller ratios.
  * check_identity (576 + 48), check_base_identity (132 + 48), check_lowfreq and verify_inputs on this checkout (floor rows, inventory, 881 mass rows, contour ceilings, NF02) all PASS.
* **Not yet independently reviewed.** The `nb_delta` seed lemma (telescoping with factor 2; uses n ≥ 997 and dn ≥ 497), the B_box/inf_BH envelopes, the floors per L-cell taken from lrt3_rcert, the contour3 import on [3.5, 3.65], and the retained inputs (angular-action identity, EPSW masses and cutoff, Loc-Joint).
* **Sliver status.**
  * Certified, conditional on retained inputs: L ≥ 2.9 (L ≥ 5 reviewed; [3.65, 5] under review; [2.9, 3.65] here).
  * Open: 998 ≤ a < a*(L), J ≥ 1001, 0.06 ≲ L < 2.9.

## 111. Positivity on L ∈ [2.6, 2.9]; direct low-L prefix truncation; global moment transform (author; reproduced by coordinator, 27 Nov)

* **Package** (author): code/sliver_next/. The coordinator's rerun is in code/sliver_next/rerun_coordinator/.
* **Cell [2.6, 2.9]** (same machinery as §110: coupled norming, ellipse extension, direct contours). It uses 12 channel boxes and 1,057,106 scalar evaluations. The p = 4 top line is at height 20.25. Contour totals are < 9.80e-10 (b = 0) and < 1.487e-9 (b = 4).
  * b = 0: E/F < 0.091148, action > 2.272e-6.
  * b = 4: E/F < 0.112364, action > 2.219e-6.
  * Floor 2.5e-6. Rerun here (732 s): results, contours, budget, ledger and low_prefix outputs are all **identical** apart from timing.
* **Floors checked** against code/lrtight/cache/lrt3_rcert.json, with β = 2 acosh(e^{L/2}); every cell's β-range is covered by contiguous rows:

  | L-cell | β-range | rows | min F_cert | floor used |
  |---|---|---|---|---|
  | [2.6, 2.9] | [3.948, 4.258] | 14 | 2.994e-6 | 2.5e-6 |
  | [2.9, 3.2] | [4.258, 4.566] | 13 | 2.547e-6 | 2.5e-6 |
  | [3.2, 3.35] | — | — | 2.346e-6 | 1.87e-6 |
  | [3.35, 3.5] | — | — | 2.160e-6 | 1.87e-6 |
  | [3.5, 3.65] | — | — | 1.98662e-6 | 1.98e-6 |
* **Low-L direct prefix theorem** (low_prefix.py, PASS; rerun identical).
  * Statement: for c = 1 − e^{−L} ∈ [0.05, 0.56] (L ≈ 0.051–0.821), all a ≥ 998 and both b, |H_ang − H_prefix| < 2^{−110}. Here H_prefix uses the exact rational ladder Π_{j<k}(c − hj(j+s))/(1 − hj(j+s)), with K = 1193–3248 on six c-cells and pole gap ≥ 0.1135.
  * Proof: the Mehler bound on \|z\| ≤ 0.995 gives \|F_{p,ν}\| ≤ (π/2)·200^p·e^{πν/2}. Cauchy on \|z\| = 0.99 then bounds the global Taylor coefficients of the integrated G (≤ 2^2554), using U_k ≤ c^k; the ν > 1100 Gaussian tail is handled separately.
  * Consequence: the sign problem becomes one bivariate rational numerator in (c, h) with GLOBAL coefficients A_k, B_k, D_k. Its sign is NOT yet certified.
* **Global moment transform** (moment_transform.cpp): 4(k+1)²I_{k+1,j} = (8k² + 8pk + 4p − 1)I_{k,j} − 4Λ²I_{k,j+1} − 4(k−1+p)²I_{k−1,j}, with directed errors. For ideal weights, 3250 entries at 2500 digits give error < 1e-463 (rerun: PASS; 10 pointwise ideal-action checks positive). The realised-weight moments (Gaussian × erf cutoff; the README sketches an exact recurrence) are not done yet.
* **Not independently reviewed.** The low_prefix constants (2^2554, dressing factor < 2, the h-cap formula (1−√c)²/998²), the moment-transform derivation, and the items listed in §110.
* **Sliver status.**
  * Certified, conditional on retained inputs: L ≥ 2.6.
  * Reduced to a rational sign problem, not yet certified: L ∈ [0.051, 0.82].
  * Open: L ∈ [0.82, 2.6], plus the sign problem above.

## 112. Agent ANALYTIC_SLIVER3 final: independent cells [2.6, 2.75] and [2.75, 2.9], cross-check cell [3.4, 3.55], laptop launcher (agent, 27 Nov; unreviewed)

* **Cells run on the author's sliver_ext scripts (called unchanged).** Both are conditional on the same retained inputs as §110–111.
  * [2.75, 2.9]: E/F = 0.122 (b = 0) and 0.176 (b = 4); 523,324 scalar evaluations, 377 s on 3 workers.
  * [2.6, 2.75]: E/F = 0.059 and 0.068; 1,215,108 evaluations, 1028 s.

  These independently re-cover the author's [2.6, 2.9] cell (§111) on a different cell split.
* **Independent cross-check cell [3.4, 3.55].** It uses the agent's own contour: a complex tiling of the whole top line, the right edge by imaginary-axis boxes, and near-real pieces with exact state constants. Result: E/F = 0.142 / 0.160, PASS. It agrees in size with the author's [3.35, 3.5] and [3.5, 3.65] cells.
* **Lemma Loc‴** (cellwise remainder, Chernoff weights, phase-aware): log E_ρ = −2063 on [3.55, 3.65], −1453 at L = 3, −437 at L = 2; it fails near 1.6. This partially reproduces §104 independently for L ≥ 2.
* **The exact state relations matter.** The old 3.9a² and N ≥ 40a relations fail near L ≈ 3.0–3.05. The near-real pieces use g = 1/(998²(1+√(1−d_hi))²).
* **Launcher.** code/analytic3/ext_launch.py (plan / preflight / estimate / run / summary) wraps ext_cell.py, which calls the sliver_ext pipeline (SLIVER_EXT defaults to code/sliver_ext). It is resumable, and results go to cache/ext/<cell>/ledger.json.
  * Measured rate: 9.5e-6 core-s per node·K (1.5e-5 used for margin).
  * Estimates: [1.8, 2.9] ≈ 60 core-h; [1.5, 1.8] ≈ 240 core-h. [1.2, 1.25] alone is ≈ 350 core-h, and K reaches 1000–1240 near L ≈ 1.0.
  * So the launcher is laptop-feasible down to L ≈ 1.5; [1.3, 1.5] takes days; below ≈ 1.2 it is not feasible.
  * PASS margins below 2.6 are untested; a FAIL means split the cell and rerun.

## 113. Low-L realised moments: the rational-pole seed errors reduce to three imaginary-frequency evaluations (author, 27 Nov; code not yet received)

* **Set-up.** For the unshifted Gaussian part of the repaired weight, define F_{j,a} = ∫₀^∞ B(ν)θ(ν)/(4ν²+a)·(ν/1000)^{2j} dν for a ∈ {1, 9, 49}.
* **Recurrence.** Since (ν/1000)² = (4ν²+a−a)/(4·10⁶), the exact recurrence F_{j+1,a} = (S_j − aF_{j,a})/(4·10⁶) holds, where S_j is the plain Gaussian moment. The coordinator checked this.
* **Error propagation.** A seed error ε_a propagates as ΔF_{j,a} = ε_a(−a/(4·10⁶))^j, so its effect on the action is ε_a times the generating series evaluated at ν² = −a/4 (ν = i√a/2). Three correlated errors therefore replace thousands of independent ones.
* **Claim (author).** Over the full 3,564-coefficient Gaussian inventory, \|ΔH_pole\| < 2.4e-20 on 0.05 ≤ c ≤ 0.56, against floors of about 1e-6.
* **Scope.** This controls only the seed errors of the rational-pole moments, including the Gaussian cutoff correction. The plain Gaussian moments S_j (high order) and the full repaired moment vector are not yet certified.

## 114. Independent referee of the author's chain §104–111: no errors (referee agent, 27 Nov)

* **Report:** REVIEW_AUTHOR_CHAIN.md; scripts and logs in code/review_chain/.
* **Verdict.** L ∈ [2.6, 3.65] positivity stands, conditional on retained inputs. The low-L prefix theorem stands as a truncation and pole-exclusion result only; L ∈ [0.051, 0.82] is not proved, because the sign of the rational prefix is uncertified.
* **Per item.**
  1. **`nb_delta` seed, downward steps, Δ-recurrence: hold.** The seed constant was re-derived term by term. Its hypotheses hold on every cell (n ≥ 997; actually dn ≥ 747, so the code's 497 is loose). Checked against exact values: Lemma V to 1e-80; nb_delta ≤ 0.48× exact; near_sum ≤ 0.67× exact on the [2.6, 2.9] bottom/left edges up to y = 20.25.
  2. **Gamma/Euler envelopes and contour normalisation: hold.** The prefactor equals \|coth πν\|/(2π) to 12 digits. The worst ratio is 0.78 against mpmath (300 complex points plus 6 boxes). The 4/Λ and 2τΛ factors and the extra segment [Λ, 998.5] match AS §4. The dressed finite sum versus the undressed infinite-twist main term is consistent with FULLPOS3 §4.1(c). Exact physical transfers against the bounds used: worst 0.96.
  3. **Base-multiplier lemma: holds.** All steps were re-derived, including negative real exponents on the top line. Against exact Arb ladders (9000 bits) the worst ratio is 0.992. A 4000-bit run gives a spurious "violation" from precision loss; the true value is 8.78.
  4. **low_prefix constants: all hold** (dressing factor < 1.082, 2^2554, Cauchy tail, pole gap).
  5. **Coverage: complete,** with no gaps in channels, bands, contour pieces or top-line breakpoints, and J0 = 1500 used consistently. The floors were recomputed from lrt3_rcert and agree with §111. Two documentation gaps:
     * (a) TAIL_CERT's L uses w = e^{−L}, while the code uses d = e^{−L}; for b = 4, d = w(1+20h). The d-slab margins (≥ 2.9e-7 relative) exceed 20·hcap (≤ 2.9e-8), so both conventions are covered.
     * (b) assemble_26.py does not itself assert the F_cert rows for [2.6, 2.9]; they are correct (§111).
  6. **Low-frequency theorem: holds** (principal constant 12.61 < 14).
* **§73 is in the repository,** contrary to the referee's note. code/fullpos3/cache/laptop_*.jsonl plus results/ht_scan_done.jsonl hold the 60 cells contiguously on β ∈ [0.39, 8.4], all ok, both channels. The cell [4.84, 5.04] has a* = 4500 (b = 0, 4), which covers the [3.5, 3.65] β-range. Coordinator check.
* **Still retained (not reviewed):** the angular-action identity; EPSW masses, cutoff and Gaussian inventory; LRINF floors; Loc-Joint internals and K₀; contour3 on [3.5, 3.65] (its state relations hold there); the helicity identity.

## 115. Positivity on L ∈ [2.3, 2.6]; realised-moment machinery for low L (author; reproduced by coordinator, 27 Nov)

* **Package** (author): code/sliver_L23/. The coordinator's rerun (1149 s) is in rerun_coordinator/.
* **Cell [2.3, 2.6].** Coupled norming with 10 p = 3 boxes (1,307,286 evaluations) and 2 p = 4 right-edge boxes (247,104). Contour heights are 15.25 (p = 3) and 30.25 (p = 4), with a normed right edge for p = 4 on Im ν ∈ [0, 25]. The Gaussian bands near 1040 are refined. Floor 3e-6: the lrt3_rcert rows on β ∈ [3.6, 4.0] have minimum 3.46e-6, and assemble_23.py asserts coverage.
  * b = 0: E/F < 0.316041, action > 2.0518e-6.
  * b = 4: E/F < 0.568905, action > 1.2932e-6. This is the thinnest margin so far, 43%.
* **Checker changes (audited by the coordinator).**
  * Norming now uses N ≥ 1.4m nodes, factor 5/2 and +7/2·ε per axis. This is valid because m(π/2N + 3e-6) < 1.13 and cos(1.13) > 2/5, which is checked exactly.
  * scalar_binary.cpp is a binary fixed-point variant. Each coefficient ratio is rounded once (l1 ≤ 2 units), and its error term (\|a\|₁+2)e_u + 2\|u\|₁ + 2 covers both coefficient rounding and propagation.
* **Reproduced.**
  * The ledger, budget, all four contour files, the p = 4 results and pole_bound_results are identical apart from timing.
  * The p = 3 transfer bounds agree to all printed digits; their exact rational tails differ only in roundoff-unit bookkeeping at the 1e-60 level.
  * inventory_cut on this checkout's k2_continuum_D_1000_rep.json gives 1030 groups and cutoff seed < 7.0957e-24. pole_bound gives cutoff pole contributions 2.30e-20 (b = 0 and 4). raw_seed_budget gives a correlated log-seed action error of 7.3e-46, and < 1e-400 for the log cutoff. test_raw_transform passes 33 outputs against 12 exact atoms. gaussian_high_test passes 7 repository centres, 3254 even moments at 2500 digits.
* **Low-L realised-moment route** (README; extends §111 and §113).
  * The raw unshifted weight g = −(R(ν²)B + R̃(ν²)𝓛)θ/Λ has three rational poles at a = 1, 9, 49.
  * Plain Gaussian moments with the erf cutoff come from an exact integration-by-parts recurrence M_{k+1} = (c/V)M_k + k(s²/V²)M_{k−1} + (s²/V)(E⁻_k − E⁺_k), implemented in gaussian_cut_moments.cpp with MPFR intervals.
  * Rational-pole seeds are handled as imaginary-frequency atoms (§113); the log seeds use G/E/K recurrences.
  * The raw triangular transform to I_{m,j} feeds the retained A_k identity.
  * **Not yet run:** the full-inventory realised A_k, and the sign of the low-L bivariate rational numerator.
* **Sliver status.**
  * Certified (conditional on retained inputs): L ≥ 2.3. The chain §104–111 has been refereed with no errors (§114).
  * Open: L ∈ [0.82, 2.3] (the laptop launcher covers [1.5, 2.6]) and the low-L sign problem on [0.051, 0.82].

## 116. First direct low-L positivity cells, c ∈ [0.05, 0.09] (author claim, 27 Nov; package NOT yet received)

* **Claim** (author). The full realised repaired weight (all 3,564 Gaussian centres, binary64 coefficients, directed arithmetic at 1600 digits) gives 1,253 enclosed integrated coefficients A_k, with widths < 1e-375. The widths apply to the seed-substituted vector; the correlated seed discrepancy is charged at the action level (§113, §115).
* **Direct prefix sign certificates** (§111 numerator), continuous in c, all twists and spins, both b:
  * c ∈ [0.05, 0.065]: normalised action > 2.20970e-5 (b = 0) and 2.19548e-5 (b = 4);
  * c ∈ [0.065, 0.09]: > 1.90522e-5 and 1.88525e-5.

  Together these cover L ∈ [0.0513, 0.0943] and overlap the earlier certified corner L ≲ 0.058.
* **Key step** (author). The signed rational Taylor enclosure is converted to a tensor Chebyshev basis before the coefficients are bounded. This takes the first channel's absolute remainder from ≈ 22,240 down to 1.785e-6.
* **Also claimed.** The prefix truncation (< 2^{−110}) is proved on 22 continuous cells covering the rest of the low-L interval.
* **Status.** The uploaded file was byte-identical to the §115 package (same SHA-256), so the direct-certificate package (110 files) still has to be received and reproduced.
* **Open after this claim:** L ∈ [0.0943, 2.3].

## 117. Lemma U certification agent final report: pipeline built, no certificate (agent, 27 Nov)

* **Verdict.** Nothing certified. Write-up: LEMU_CERT.md; code: code/lemU_cert/.
* **Single-c check.** Run with the repaired weight at L = 2 and L = 0.5, both b. It fails with useless ball radii: P̂₀ = ±2.4e439 at L = 2 and ±6.7e115 at L = 0.5.
  * A₀ is correct (6.1e-6 ± 7.6e-8 at L = 2), and the residual and quadrature budgets are as planned.
  * The blow-up is in the ball contraction against R_{m,j}; the cause is not yet located. code/lemU_cert/lu_debug_radius.py locates it.
  * Measured cancellation of the integrated G: e^238 at L = 2, e^70 at L = 0.5.
* **Rigorous pieces built.** Complex-frequency jet bounds via a Laplace integral; Gauss–Legendre on Bernstein ellipses with rigorous nodes; the realised weight via EPSW's code; residual cells.
* **c-direction.** A new exact identity R_{m,j} = x^m d^{2m−j} R′_{m,j}(v), with v = d/c and R′ a polynomial whose coefficients in −v are all nonnegative. On it, a Taylor-model cell in w = v − v₀ that reproduces direct evaluation to 15 digits. Not run end to end.
* **Cost** (untested on a passing cell): ≈ 600–1100 core-h for β ∈ [1.47, 5.1].
* **Low L.** Rebuilding the jets from global A_k is about 50× cheaper than per-c quadrature and matches to 7–11 digits for j ≤ 1000, but needs K ≈ 7000–8000 coefficients for j ≥ 2000 when c ≤ 0.56. The recommendation for L ≲ 0.82 is global realised A_k with K ≈ 8000 (author's tools, §115), then contraction plus a cell model.
* **Assessment (coordinator).** For low L, the author's direct prefix method (§111, §116), which does a Chebyshev sign test on the bivariate numerator, avoids the R-contraction cancellation altogether and has already produced claimed PASS cells. It is the route to pursue. The Lemma U route is parked.

## 118. Direct low-L positivity, c ∈ [0.05, 0.09], with the full realised weight (author; reproduced bit-for-bit by coordinator, 27 Nov)

* **Package** (author): code/sliver_direct/. Large regenerable intermediates are not committed: the shard files g1254_1600_*.txt, the *.in inputs, raw_coefficients.txt and ideal4.txt (SHA256SUMS lists them). The coordinator's logs are in rerun_coordinator/.
* **Inputs.**
  * **Inventory.** weights.json equals ads4/prl/data/k2/k2_continuum_D_1000_rep.json exactly: all 3,564 centres, widths and x_re, plus y_re. The source blob is d7e513e1 (coordinator check).
  * **Moments.** Realised Gaussian moments S_j, j ≤ 1254, at 10^1600 scale, from the exact integration-by-parts recurrence with an erf cutoff (MPFR, directed). Adjacent-centre erfc increments use a Cauchy-bounded Taylor series.
  * **Seeds.** Correlated seeds (F₁, F₉, F₄₉, K₀, LF₉, LF₄₉) have widths < 1e-40. They are charged at action level:
    * Gaussian pole correction < 2.4e-20 and log-seed correction < 7.35e-26 for c ≤ 0.56;
    * < 1.8e-15 and < 1e-20 respectively for c ≤ 0.9.
  * **Coefficients.** Raw triangular transform to A_k, then the retained raw angular identity; the ideal G_D comes from exact NF02 ψ₃. The result is 1,253 coefficient enclosures, widths < 1e-375 for the seed-substituted vector.
* **Prefix theorem** (full_prefix.py): \|H − H_prefix\| < 2^{−110} on 22 contiguous cells tiling c ∈ [0.05, 0.9], with K up to 17,847 and minimum pole gap 0.050.
* **Sign certificate** (direct_sign_cheb.py/.cpp).
  * Spin cap h ≤ c/((1001−b)(1002+b)) covers every J ≥ 1001, with h → 0 at infinite spin.
  * The nested rational sum S_k = F_k + (c−p_k h)/(1−p_k h)S_{k+1} is expanded as an exact Taylor rectangle in (u, v). Its error recurrence is E′ = ((\|a\|+DC+r)E + M·d)/(d−r) + F_r.
  * The analytic tail comes from a Horner majorant on a polydisk with exact pole gaps: M(Rc^{−mc−1} + Rh^{−mh−1})/((1−1/Rc)(1−1/Rh)).
  * The rectangle is converted to tensor Chebyshev. The weights are nonnegative with column sums 1, so the l1 error does not grow.
  * Lower bound = c₀₀ − Σ\|others\| − rounding − tail, normalised by (1−c)³/(1+qh)⁴, minus the 1e-18 seed and truncation allowance.
  * The coordinator audited both files line by line.
* **Results** (PASS_CONDITIONAL, normalised action lower bounds):

  | c-cell | K | (mc, mh) | b = 0 | b = 4 |
  |---|---|---|---|---|
  | [0.05, 0.065] | 684 | (369–371, 241–242) | 2.2097e-5 | 2.1955e-5 |
  | [0.065, 0.09] | 777 | (433–435, 325–326) | 1.9052e-5 | 1.8853e-5 |

  Together these cover L ∈ [0.0513, 0.0943], overlapping the earlier corner L ≲ 0.058.
* **Reproduction** (coordinator, 386 s plus 384 s).
  * Everything was rebuilt from source.
  * With the author's 4-way sharding, S_1254_1600.txt, raw_coefficients.txt and direct_coefficients.txt are byte-identical, and all four sign ledgers are identical apart from timing. (With 3 shards the interval endpoints differ slightly, and the four cells still pass with the same lower bounds to 7 digits.)
  * verify_direct_results.py: PASS (inventory, seed widths, exact moment sum, coefficient radii, Chebyshev transform, seed budget, four ledgers).
  * inventory_cut, pole_bound(_full), raw_seed_budget(_actual/_full): PASS, with the numbers above.
* **Scope.** Coverage for c > 0.09 needs a larger moment table: K ≈ 3,250 at c = 0.56, and the prefix theorem allows up to 17,847 at c = 0.9. The [2.15, 2.3] slab was interrupted and makes no claim.
* **Open after this:** L ∈ [0.0943, 2.3].

## 119. Referee of L ∈ [3.65, 5] (ANALYTIC_SLIVER2): all four cells survive (referee agent, 27 Nov)

* **Report:** REVIEW_SLIVER2.md; scripts in code/review_s2/.
* **Lemmas.**
  * **Lemma F holds.** The frames are affine; the defect is bounded by Taylor with the Hessian enclosed in Arb over the whole box. Against exact values on 148 certificate boxes the bound is never closer than e^{2.49}.
    * Nit S4: in 17% of boxes the float midpoint ± radius misses the box by about 1e-15. Widen the radius; no state is lost.
  * **Lemma B holds.** Each frame family is rigorous over the whole box, so the per-j minimum is valid; it is never closer than e^{1.5} to exact.
  * **Loc″ holds.** The ratio is < 1 on all 2294 cells, and ψ₂ = 0 above Λ, so J′ = 2200 is legitimate.
  * **rmass3** regenerates bit-identically.
* **Gap S1, closed by the referee.** The state relations were needed on [3.65, 3.75] but stated only for L ≥ 3.75. They are now proven for all L ≥ 3.65: a/(Nd) ≤ 0.50337, Dd²/a² ≥ 3.947, N/a ≥ 76. The text's "Dd² ≥ 3.95a²" is false at 3.65, but nothing uses it.
* **Contour on [3.65, 3.75].** The hypotheses hold: (A+B)/den = 1.9745 ≤ 1.99; F/den = 0.930; the j-tail ratio is 0.486 < 0.5 (a narrow pass); the top line stays clear of the B poles. The right edge (0.32·F) rests on **Lemma X (AS §5.5), still unreviewed**; numerically it sits above the truth by at least e^{8.35}.
* **S3, wrong justification, negligible effect.** The claim "\|f_θ\| ≤ \|f\| above Λ" is false, because θ is evaluated at ν − i (θ(1098.5 − i) = 0.49999 + 0.0566i). The correct omitted term is ≤ 7e-11·F.
* **Reruns.**

  | cell | ρ |
  |---|---|
  | [3.65, 3.75] | 0.545 (bit-for-bit) |
  | [3.75, 4] | 0.223 (bit-for-bit) |
  | [4, 4.5] | 0.0894 (bit-for-bit) |
  | [4.5, 5] | cached 0.0267, rerun 0.0177 |

  S2: the cached [4.5, 5] value comes from the older single-family Lemma B call, which reproduces exactly when patched in. Both versions are rigorous.
* **Status.** [3.65, 5] stands. [3.65, 3.75] is conditional on Lemma X, which can be removed by covering [3.65, 3.8] with the author's norming pipeline (§109–110); see §120.

## 120. Direct low-L positivity extended to c ∈ [0.09, 0.258]; conformal prefix theorem (author; reproduced by coordinator, 27 Nov)

* **Package** (author): code/sliver_direct_cont/, a delta on top of code/sliver_direct/. The bundled original zip is byte-identical to the §118 package. The large .in files are not committed. The coordinator's logs are in rerun_coordinator/.
* **New checkers.**
  * direct_sign_triangle.cpp/.py stores A_k, G_k × 2^{−k}. This doubles the Horner ratio, and the dressing uses 2(k+1)²A′_{k+1}. The truncation is a downward-closed Taylor triangle {a_c i + a_h j < Q}, with polydisk radii 2^{a_c}, 2^{a_h} and Cauchy tail v(Q+2)2^{1−Q}.
  * direct_sign_tight.py does the same on [0.223, 0.258]. There the prefix truncation error is relaxed locally to 3e-6: the exact ordinary-power tail at c = 0.258 with K = 1251 is 2.5901e-6, from |A_k|+|G_k| ≤ 2^2417(199/200)^{−k}. Its sign test subtracts 3e-6 + 1e-18.
  * The coordinator audited both: the recurrence and error growth (\|a\|+dc+2r), and the triangle tail count of at most m+1 pairs per weighted degree.
* **Results** (PASS_CONDITIONAL, normalised action lower bounds, both b):

  | c-cell | method | b = 0 | b = 4 |
  |---|---|---|---|
  | [0.09, 0.12] | Chebyshev rectangle | 1.7305e-5 | 1.7064e-5 |
  | [0.12, 0.15] | Chebyshev rectangle | 1.5690e-5 | 1.5385e-5 |
  | [0.15, 0.18] | triangle | 1.4713e-5 | 1.4341e-5 |
  | [0.18, 0.223] | triangle | 1.3289e-5 | 1.2845e-5 |
  | [0.223, 0.258] | tight | 9.9709e-6 | 9.4152e-6 |

  Together with §118 this gives direct coverage of c ∈ [0.05, 0.258], i.e. L ∈ [0.0513, 0.29841].
* **Reproduction** (637 s). All 10 cell .json/.out files and extension_ledger.json are byte-identical (apart from timing). verify_extension.py passes: source transfer, spin caps, pole gaps, Cauchy majorants and tails, final inequalities. prefix_improvements.py passes (direct tail at 0.258: 2.59e-6; conformal tail at 0.9 with J = 5333: 1.33e-34 < 2^{−110}). verify_conformal_recurrence.py passes (symbolic identity plus 1980 finite-spin checks).
* **Conformal prefix theorem** (proved; coefficient table not yet built).
  * Variable: w = (1−√(1−c))/(1+√(1−c)), c = 4w/(1+w)², so w ≤ 0.52 for c ≤ 0.9.
  * The Mehler–Dirichlet bound with \|Im β\| ≤ 4 atan√0.99 < π gives \|Ã_j\| + \|G̃_j\| ≤ 2^4813·0.99^{−j}.
  * Lagrange inversion gives [z^k]w^j = (j/k)4^{−k}C(2k, k−j) ≥ 0, hence T_h[w^j] ≤ w(c)^j (coordinator checked the Lagrange step).
  * Consequence: 5,334 coefficient pairs suffice on the whole band c ≤ 0.9, against 17,849 in the ordinary basis.
  * The five-term transformed ladder recurrence is verified, but a stable continuum summation and the transformed seed-error budget are still open.
* **Sliver status.**
  * Certified (conditional on retained inputs): L ≤ 0.2984 (corner plus direct), and L ≥ 2.3.
  * Open: 0.2984 < L < 2.3, i.e. c ∈ (0.258, 0.8997).

## 121. Lemma X (AS §5.5) refereed: holds (referee agent, 27 Nov)

* **Report:** REVIEW_LEMMA_X.md; scripts and logs in code/review_lemx/.
* **Verdict: holds.** The proof was re-derived step by step; no inequality has the wrong direction and no constant is wrong. There are five imprecisions (N1–N5), each covered by slack already in the code:
  * N1: the first window step starts from y ≤ y_A (about 0.03 nats), covered by the code's +3 steps.
  * N2: phase III drops a boundary term (0.1–0.27 nats).
  * N3: lmax sits about 3% below \|λ₊\| near y_B.
  * N4: some redundant slack.
  * N5: c = 0.5 is used only by the DV driver, not by the certificate.

  An independent bound with N1–N3 fixed stays at least 0.066 nats below rt_pp.cross for every c ∈ {0.7, 1, 1.5}, box [t−1, t] with t ∈ [949, 1150], and σ ≤ 6.5.
* **Hypotheses on the [3.65, 3.75] right edge.** piece_right → PP → cross with τ ∈ [Λ−1, Λ_up] and σ ≤ 5.5 (helicity shift included), for all a ≥ 998 and any spin. Since τ = Λ < 999, every step has y_k < 0.9995: the turning point is never crossed, and only phases I and II are active. The lemma makes no assumption on L, d or N/a.
* **Numerical margins** (exact \|Ŝ\| in Arb at 12000 and 16000 bits; 21 values of a in [998, 4500], L ∈ {3.65, 3.70, 3.75}, σ ∈ [0, 5.5], τ ∈ [Λ−1, Λ]):
  * Lemma X (b = 0): smallest gap e^{7.727}, at a = 998, n = 83372, σ = 0.
  * PP with helicity (b = 4): e^{8.21}.
  * Every intermediate inequality holds on 150 trajectories.
* **Alternative.** base_bounds (§109) would also serve on this edge and is 4–6.5 nats sharper, but would need rewiring.
* **Status.** With §119, the cell [3.65, 3.75] no longer has an unreviewed dependency, so L ∈ [3.65, 5] is fully refereed. The remaining retained inputs are the global ones: the angular-action identity, EPSW masses and cutoff, the LRINF floors, and §73.

## 122. Larger realised coefficient table (K = 3300); direct cells to c = 0.263; driver h-cap limitation (agent, 27 Nov)

* **Table.** run_full_gauss.py with K = 3300 and D = 2600, followed by make_coefficients.py, produces 3299 coefficient pairs. These are enough for need_k up to c = 0.56.
  * Radii: A ≤ 1e-466 and G ≤ 1e-543 at k = 3298. Time: 854 s plus 23 s.
  * Only inputs were changed: K and D in make_coefficients.py, and the log_constants.cpp output scale raised to 10^2600 with 10528 bits. At the old fixed 10^1600 it silently gave A₀ = 26.73 instead of 20.77, so the scale must track D.
  * Cross-checks: all 1253 old coefficient intervals overlap the new ones, and the author's tight cell [0.223, 0.258] reproduces the exact same physical_lower.
  * The tables are saved compressed in code/sliver_direct_ext/tables/.
* **New PASS cells** (author's direct_sign_triangle, unchanged; verified with a script based on verify_extension.py):

  | c-cell | b | K | lower bound |
  |---|---|---|---|
  | [0.258, 0.262] | 0 | 1397 | 1.404e-5 |
  | [0.258, 0.262] | 4 | 1397 | 1.113e-5 |
  | [0.262, 0.264] | 0 | 1405 | 1.094e-5 |
  | [0.262, 0.263] | 4 | 1401 | 1.74e-6 |

  Direct coverage now reaches c = 0.263, i.e. L ≈ 0.305.
* **Limitation.** The drivers bound h by c_hi/((1001−b)(1002+b)), which uses only J ≥ 1001. Above c ≈ 0.251 the physical constraint a ≥ 998, h ≤ (1−√c)²/998², is tighter. The rectangle therefore contains non-physical (c, h) points where the prefix numerator is genuinely negative or huge, whatever K is (b = 0 corner probe: −1.5e-5 at c = 0.266). That is why cells above 0.264 come out inconclusive.
* **Proposed fix (not applied; needs the author's decision).** In the driver, use h ≤ min(c_hi/((1001−b)(1002+b)), (1−√c_lo)²/998²). This is a valid superset of the physical states, since (1−√c)² decreases in c. The C++ is unchanged. Non-rigorous probes with this cap show a smooth positive numerator of about 1.1–1.5e-5 up to c = 0.56 (b = 0). Estimated cost to c ≈ 0.51 is 2–8 CPU-h. Above about 0.51, the smallest h-radius option (2) hits the pole (1.5hK² ≥ 1), so a smaller radius option or a smaller K is needed.

## 123. Direct cell c ∈ [0.268, 0.28] (author claim, 27 Nov; package not yet received)

* **Claim** (author): PASS in both channels on c ∈ [0.268, 0.28]: b = 0 lower bound 1.3553e-5, b = 4 lower bound 1.2991e-5. This is L ∈ [0.3120, 0.3285].
* **Coverage gap to close.** The reproduced coverage ends at c = 0.264 for b = 0 and c = 0.263 for b = 4 (§122). The strip c ∈ [0.263, 0.268] (b = 4) and [0.264, 0.268] (b = 0) is not covered unless the author's package contains cells for it.
* **New idea** (author, not yet received): a partial-fraction representation of the angular ladder with explicit control of its rational poles. It would replace the two-variable Taylor certificates by one-variable polynomial certificates plus a controlled approximation in the twist direction.

## 124. Author: cell [0.258, 0.268]; full moment table; termination-preserving reduction; degree-600 twist approximation (author claims, 27 Nov; package of 116 MB not yet received)

* **Cell.** c ∈ [0.258, 0.268] passes in both channels: b = 0 > 1.3804e-5, b = 4 > 1.3220e-5. Together with the claimed §123 cell, this fills the strip noted in §123 and puts the open band's lower end at c = 0.28 if both are reproduced. The author quotes L = −log(0.732) ≈ 0.312, which corresponds to c = 0.268.
* **Moments.** All 3,564 centres; 6,001 enclosed orders at 3,000 digits; a faster scaled-integer recurrence. Truncation in spectral degree gives an error < 2^{−1910}. New Gaussian seed action error < 1e-45. The coefficient vector has 17,201 pairs; an unsigned-integer error in the ideal transport was found, fixed and regression-tested.
* **Cover.** 21 pole-free continuous rectangles cover the whole remaining band, with maximum ordinary prefix degree 16,989 and minimum denominator gap > 0.1447.
* **Obstruction (rigorous).** The unrestricted rational prefix is NEGATIVE at b = 4, c = 0.35, h = 2.33e-7. That point violates the twist constraint, so it is not a counterexample. It shows that certificates must be restricted to the physical domain; compare the §122 diagnosis and GAP_LEMMAS G1.
* **Termination-preserving coordinate.** Use x = h/c = 1/(n(n+s)), so U_k = c^k Π_{j<k}(1 − x·j(j+s))/(1 − c·x·j(j+s)). The physical domain is x ≤ min(1/((1001−b)(1002+b)), (1−√c)²/(998²c)).
* **Residue checker** (exact arithmetic, PASS on 21/21 rectangles, both b). A degree-600 polynomial in the twist variable approximates the prefix with uniform error ≤ 5.31e-23, continuously in x.
* **Open.** The sign of the degree-600 approximating polynomial on the physical domain. Since 27 Nov the coordinator's agents have been working on it in GAP_LEMMAS.md: the designer on the sign certificate, the prover on G1-capped direct cells from c = 0.28 upward.

## 125. Full-band reduction package: reproduced (fast path) by coordinator (27 Nov)

* **Package** (author): code/sliver_fullband/. Received in 5 split parts, joined (121 MB zip, integrity OK) and unpacked to 254 MB. Large data files are not committed (coeff_long_17200.txt 192 MB, S_6000_3000.txt, coeff_long_1800.txt, old_coefficients.txt, *.in, *.seeds); their SHA-256 sums are in LARGE_FILES_SHA256.txt. Rerun log: rerun_coordinator/reproduce_new.log.
* **Reproduced (fast path, reproduce_new.sh, 86 s, all PASS).**
  * Seeds: Gaussian M0 coherent seed error < 1e-45 for both b, and all 3564 seeds and fixed-point constants (12 shards).
  * Spectral-degree reduction: spectral degree 6000 with error < 2^{−1910}; ordinary K = 17199 with tail < 2^{−140}.
  * Finite reduction of the whole remaining band: 21 contiguous rectangles, K_max = 16989, minimum pole gap 0.14470.
  * Correlated log-seed action error 8.2e-22; inventory cutoff < 7.1e-24.
  * Sign cell c ∈ [0.258, 0.268]: the direct_sign_triangle outputs for both b are **byte-identical** to the archived ones. b = 0 > 1.38044e-5, b = 4 > 1.32204e-5. The executed sign ledger reaches c = 0.268.
  * Transport: 236 exact transport checks (both channels, active spectral truncation); ideal L1 bounds; all 1253 old ideal enclosures overlap the corrected ones (centre difference < 1e-107).
  * Non-physical witness: an exact negative value of the rational extension, −5.91e18, at a point outside the physical domain.
  * Termination-preserving split: 42 channel rectangles and 10214 finite-spin interval checks (distinct n: 4801 for b = 0, 4841 for b = 4). This is a cover only; nothing on it is sign-certified.
* **Not done here.** The full data regeneration (`--regenerate`: 6000-order moments at 3000 digits, long transform); the archived 17,201-pair vector is checked only by verify_long_transport. The [0.268, 0.28] cell (§123) and the degree-600 twist residue checker are not in this package.
* **Status.** Direct coverage of c ≤ 0.268 (L ≤ 0.31197) is reproduced. Open: 0.312 < L < 2.3; the laptop launcher is running [1.5, 2.3], and the gap agents are on the low side.

## 126. Gap agents: designer report; prover cells c ∈ [0.28, 0.56] both b; G1, G2 proved (27 Nov)

* **Prover** (code/gap_prove/ledger.json, GAP_LEMMAS.md; still running):
  * G1 (physical twist cap h ≤ min(c/((1001−b)(1002+b)), (1−√c)²/998²)): PROVED (+ check PASS).
  * G2 (seed/tail allowance 1e-7 + 1e-14 for c ∈ [0.56, 0.9] with coeff_long_17200): PROVED (+ check PASS).
  * G1-capped direct sign cells PASS_CONDITIONAL for c ∈ [0.28, 0.56] in both b (b = 0 to 0.585), i.e. L up to 0.8210. They are verified by verify_gap_cells.py, modelled on the author's verify_extension.py.
* **Designer** (code/gap_design/; finished).
  * On the physical range y = h/hcap ∈ [0, 1] the prefix is affine in y to 1e-8 relative, and it explodes at y ≈ 1.1 for every c.
  * Cell-size rule: the overshoot ω = (1−√c_lo)²/(1−√c_hi)² must be ≤ 1.08; b = 4 is about 5× more sensitive near the threshold. The author's [0.268, 0.35] failure and the c = 0.35 witness both sit at y ≈ 1.4.
  * Above c ≈ 0.5, the author's 2^{−k} storage lets rounding grow by up to 1e427. At D = 100 the output is garbage, though not a false sign. The fix is to choose D from the a-priori growth bound, or use the unscaled checker direct_sign_triangle_u.cpp, which gives identical results.
  * The designer's own PASS cells:

    | cell | b | lower bound |
    |---|---|---|
    | [0.3, 0.3097] | 0 | 1.33e-5 |
    | [0.6, 0.6068] | 0 | 9.88e-6 (two h-pieces) |
    | [0.77, 0.7781] | 0, 4 | 6.95–7.04e-6 (three h-pieces, ≈ 23 min each) |
  * A 102-rectangle plan for [0.56, 0.777] is predicted to pass (L1/centre ≤ 0.053).
  * Route ranking: (A) G1-capped cells — the cheapest; ≈ 15–20 CPU-h for the whole band, ≈ 10 on [0.56, 0.777]. (B) The author's termination-preserving split, a fallback at similar cost. Not recommended: the conformal route, since \|Ã_j\| plateaus at 1e1390; and the monomial twist polynomial or dominant-term idea, which hides 40–90 orders of cancellation, so the right twist basis is Chebyshev in y.
* **Status.** The band c ∈ (0.268, 0.28) is covered by the author's claimed cell (§123), whose package is still pending. Open on the low side: c ∈ (0.56, 0.777), i.e. L ∈ (0.821, 1.5), which is being executed with the designer's plan. The laptop covers [1.5, 2.3].

## 127. FULL SLIVER COVERAGE: c ∈ [0.258, 0.9] both b by coupled-cap norming (author, 27 Nov); fast checks reproduced

* **Package** (author): code/sliver_fullcov/; the proof is in sliver_merge/PROOF.md and the ledger in sliver_merge/FULL_COVERAGE_LEDGER.json. coeff_compact60.txt is committed xz-compressed. The *.in inputs, old_coefficients.txt and the long-table prefix are not committed; their SHA-256 sums are in LARGE_FILES_SHA256.txt.
* **Method.**
  * Coupled coordinate y = h/H(c) ∈ [0, 1], where H(c) is the affine secant of the convex physical twist cap (1−√c)²/998². The secant lies above the cap, so the rectangle contains every physical state, and its overshoot is small: under 2% even on the widest cell. This resolves the non-physical obstruction of §122 and §124.
  * The pole and ratio checks are done at the corners.
  * Approximation error in y comes from the partial-fraction residues of the rational ladder; in c, from a Möbius-disk product majorant.
  * Degree-matched trigonometric norming: Lebesgue constant Λ_m = 1 + 2Σ1/(2k−1), with Markov control of node perturbations. Grouped c-norming uses a self-consistent bound, requiring λ < 0.1.
  * Directed scaled-integer scalar recurrences, with SIMD exact division; a portable GMP twin is included.
  * Shortened group prefixes carry an explicit tail term.
* **Result** (FULL_SLIVER_PASS_CONDITIONAL). 18 cells; the 11 new ones use 5,669,470 scalar evaluations. Each cell is a single rectangle covering all physical J ≥ 1001 and a ≥ 998. The union is c ∈ [129/500, 9/10] for both b, i.e. L ∈ [0.29841, 2.30259], which also covers the w-convention. Common lower bound on the normalised action: **5.14e-6** (b = 0: 5.166e-6, b = 4: 5.148e-6).
* **Reproduced by the coordinator.**
  * check_ledgers.py (111 s): every cell is verified (analytic error sums, rational nodes, pole gaps, coefficient transfers, scaling schedules, node counts, reference Bernstein coefficients, final inequalities). The exact union is c ∈ [0.258, 0.9] for both b.
  * authenticate_original.py against the §125 coeff_long_17200.txt: the original SHA and the exact retained prefix match, and all 17201 compact pairs are enclosed.
* **Still to do.**
  1. The full portable-GMP replay (replay_scalar.py), comparing every archived scalar shard byte for byte. The author's laptop is running it.
  2. An independent referee of PROOF.md's new mathematics: the residue bound behind E_y, the Möbius majorant behind E_c, the norming constants, and the self-consistent grouped norming of §3b.
* **Status.** If items 1 and 2 pass, the sliver is fully certified: together with L ≤ 0.298 (§118, §120, §125) and L ≥ 2.3 (§§100–121), this covers every L. That makes assumption (iii), full AdS positivity at Δ_g = 1000, hold, conditional only on the retained global inputs: the angular-action identity, the EPSW spectral data and cutoff, the LRINF floors, and the §73 high-twist scan.

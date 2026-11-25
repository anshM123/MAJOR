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

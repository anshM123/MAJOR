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

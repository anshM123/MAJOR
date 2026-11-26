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

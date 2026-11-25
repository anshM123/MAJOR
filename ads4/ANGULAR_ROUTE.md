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

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

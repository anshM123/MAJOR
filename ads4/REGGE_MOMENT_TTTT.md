# The Regge-moment formula for heavy blocks: scalar case re-derived and tested, spinning case as far as it goes

Status: written 16 Nov 2026 by a derivation agent; **not refereed**. Code and recorded output:
`ads4/code/reggemoment/` (index in §8). Sources used: LaTeX sources of arXiv:2106.10274 ("CHMRS"), 2311.04271
("CLS-D"), 2201.06602 ("CLPS"), 2102.08160 (Caron-Huot–Li, "CHL"), downloaded to the session scratchpad
`…/scratchpad/reggemoment/src/`; equation numbers are those of the arXiv versions.

Marks: **[Q]** quoted (reference, equation); **[DV]** derived here and verified numerically (script named);
**[D]** derived here, not verified; **[I]** inferred (argument given, no derivation).

## 0. Summary

1. **Scalar case (Stage 0): verified with exact blocks.** The heavy density of the Regge moment `Π̂_{k,ν}` of an
   exact `d = 3` block, computed with no large-`Δ` approximation, equals the CHMRS formula
   `γ_{2Δ_φ+k−1}(ν)²/a_{Δ_φ} · 2P_{−1/2+iν}(η_AdS)/m^{2k}` up to `1/m²`: e.g. `Δ_φ = 3`, `k = 2`, `ν = 0.7`, `β ≈ 0.65`:
   ratio `1.0385, 1.0073, 1.00072` at `m = 17.6, 36.7, 113.1`. Including the paper's `1/m²` and `1/m⁴` terms the ratio
   is `1 + 4·10⁻⁶` and `1 − 1·10⁻⁸` at `m = 113`, with a residual `∝ m⁻⁶`. The same holds for the physical functional
   `C_{k,ν}` through (C.16). **By-product:** the Pochhammer symbol `(Δ_φ−d−1)_3` printed in CHMRS (3.27) must read
   `(2Δ_φ−d−1)_3` (the paper's ancillary file has the latter; the printed form leaves an `O(1/m²)` residual).
2. **Heavy limit of conserved spinning three-point functions (Stage 2a): verified.** In the CHMRS symmetric
   frame, for `l/Δ` fixed, the conserved `⟨JφO_l⟩`, `⟨TφO_l⟩`, `⟨TTO_l⟩` structures (and their Lorentzian-shadow
   partners) become *scalar plane wave × constant polarisation tensor*, and the tensor is **transverse at the
   point of `H²` of the Regge pair and carries helicity `±l_i` only**, up to `O(1/m)`. For `⟨TTO⟩`: exactly two
   parity-even conserved structures; non-helicity harmonics `/` total `= 0.18, 0.017, 0.0017` at `Δ = 60, 600, 6000`.
3. **Transform to the transverse harmonic basis (Stage 2b): derived and verified.** A tensor Regge-moment
   kernel `|y|^{a−3} 𝔻_s P_{−1/2+iν}(−x̂·ŷ)` Laplace-transforms into `C_aγ_a(ν)|p|^{−a}` times an explicit **rational
   function `M_s(ν; a)`** times the harmonic function of the helicity transfer: `P_{−1/2+iν}` (weight 0),
   `P^{(2)}` (currents, flip), `P^{(4)}` (stress tensors, flip). Checked by direct integration to 9 digits.
4. **Answer to the question of Stage 2.** In the helicity-preserving components there is **no helicity-
   dependent function of `β` at leading order**: the `β`-dependence is exactly `P_{−1/2+iν}(cosh β)`. There **is** a
   `ν`-dependent factor, `γ_a(ν)² M(ν;a)²`, but it belongs to the functional (it is the spinning version of the
   `γ²` that CHMRS divide out), it is the same for `++` and `+−` for parity-even kernels, and for the
   `1↔3`-antisymmetric (parity-odd-at-each-pair) kernels it comes with the relative sign `|c^{++}|² − |c^{+−}|²`
   (the `B_3`-type sum rule).
5. **Subleading terms (Stage 3).** With frames and variables built from `(α, ᾱ)`, the `O(1/m)` terms of the left
   and right vertices cancel in the `1↔3`-symmetric part: corrections start at `1/m²` **[DV pointwise for
   currents]**. Their coefficients for spinning operators were **not** derived. The helicity-dependent part
   is fixed only in the flat corner `β → 0`: `+(1/m²)[2λ²(ν/β)J_1(νβ) + λν²J_0(νβ)]`, `λ = 4` for `+−` **[DV in flat
   space; AdS form I]**.
6. **Not achieved:** an end-to-end test on a spinning block (no spinning block code); the overall
   normalisation relative to mean-field theory in the helicity basis (inferred from CHL); the identification
   of the CLS-D functionals with definite Regge-moment kernels in `d = 3`. Exact blocking steps in §7.

## 1. Stage 0 — the scalar case

### 1.1 Definitions (all [Q] CHMRS)

* Cross-ratios `u′ = ww̄`, `v′ = (1−w)(1−w̄)`; `G = Σ p_{Δ,J} G^s_{Δ,J}`, with `⟨φφφφ⟩ = G/((−x_13²)^{Δ_φ}(−x_24²)^{Δ_φ})` (3.3), i.e.
  `G^s_{Δ,J} = u′^{−Δ_φ} g_{Δ,J}(w,w̄)`, `g ∼ w^{(Δ−J)/2} w̄^{(Δ+J)/2}` (the normalisation in which (2.33) is the mean-field
  coefficient; checked at `n = 0`).
* (2.36)–(2.37) `B_{k,v}[F] = Γ(k/2)²/(π²Γ(k−1)) ∫_v^∞ dv′ ∫_0^{(√v′−√v)²} du′ (v′−u′) dDisc_s[F(u′,v′)] / ((u′v′)^{k/2}[v² − 2(u′+v′)v + (u′−v′)²]^{(3−k)/2})`,
  `k = 2, 4, …`, `v > 0`; sum rule (2.38) `Σ p_{Δ,J} B_{k,v}[G^s_{Δ,J}] = (−1)^{k/2−1}`.
* (2.39)–(2.41) Regge moments `Π_{k,η}[G] = ∫_0^{r_max(η)} dr r^{k−2} dDisc_s G(r,η)`, `r = √(ρρ̄)`, `η = (ρ+ρ̄)/(2√(ρρ̄))`,
  `ρ = 1/(√(1−w)+√(−w))²` (`u`-channel radial coordinates; Regge limit `r → 0`).
* (2.45)–(2.47): `B_{k,v} ∼ (2^{3k−2}Γ(k/2)²/(π²Γ(k−1))) ∫_{√v}^∞ η dη (η²−v)^{(k−3)/2} Π_{k,η} + O(Π_{k+2})`, inverted by
  `((−1)^{k/2} πΓ(k)/(2^{3k−2}Γ(k/2)²)) ∫_{η²}^∞ dv (v−η²)_+^{−(k+1)/2} B_{k,v} ∼ Π_{k,η} + O(Π_{k+2})`.
* (2.50) `Π̂_{k,ν} = ∫_1^∞ [dη] P_{(2−d)/2+iν}(η) Π_{k,η}`, `[dη] = 2^{d−2}(η²−1)^{(d−3)/2} dη`; `P_J(η) = ₂F₁(−J, J+d−2; (d−1)/2; (1−η)/2)`.
* (2.55)–(2.57) `γ_a(ν) = Γ((1+a−d/2−iν)/2)Γ((1+a−d/2+iν)/2)`, `a_{Δ_φ} = π^{(d−3)/2}Γ(Δ_φ)²Γ(Δ_φ−(d−2)/2)²/(2^{d−5}Γ((d−1)/2))`,
  `C_{k,ν} ≡ (a_{Δ_φ}/γ_{2Δ_φ+k−1}(ν)²) · ((−1)^{k/2}πΓ(k)/(2^{3k−2}Γ(k/2)²)) ∫_1^∞ [dη] P_{(2−d)/2+iν}(η) ∫_{η²}^∞ dv (v−η²)_+^{−(k+1)/2} B_{k,v}`,
  so that (2.58) `C_{k,ν} ∼ (a_{Δ_φ}/γ²) Π̂_{k,ν} + O(Π_{k+2})`; exact first correction (C.16).
* Kinematic map (2.28)–(2.30): `α = Δ−J−d+1`, `ᾱ = Δ+J−1`, `m² = αᾱ`, `η_AdS = cosh β = (α²+ᾱ²)/(2αᾱ)`
  (so `tanh(β/2) = (ᾱ−α)/(ᾱ+α) ≈ J/Δ`). The constants are fixed by the Weyl group and remove odd powers of `1/m`.
* Heavy density (2.31)–(2.34): `ω[Δ,J] = (−1)^J ω[G^s_{Δ,J}] / (q_{Δ,J} · 2sin²(π(τ−2Δ_φ)/2))`,
  `q_{Δ,J} = (1/(π p^{MFT}_{Δ,J})) (2n_J/m^{2d−4}) Γ(Δ−1)(2Δ−d)/Γ(Δ−d+2)`, `n_J = 2^{d+1}π^{(d−1)/2}(J+1)_{d−3}(2J+d−2)/Γ((d−1)/2)`;
  `⟨f⟩ = Σ_{τ>Δ_gap} 2sin²(…) q_{Δ,J} p_{Δ,J} f(Δ,J)`.
* Heavy-block formulas: (3.19)–(3.20) `lim Π̂_{k,ν}[Δ,J] = (γ_{2Δ_φ+k−1}(ν)²/a_{Δ_φ}) · 2P_{(2−d)/2+iν}(η_AdS)/m^{2k}`;
  (2.59) `C_{k,ν}[Δ,J] = (2P_{(2−d)/2+iν}(η_AdS)/m^{2k})(1 + O(ν²/m², 1/m²))`, regime `m → ∞` at fixed `η_AdS`, `ν ≪ m`;
  (3.27), with `a = 2Δ_φ+k−d/2`, `P = P_{(2−d)/2+iν}(η_AdS)`:
  `(m^{2k}/2) C_{k,ν}[Δ,J] = P + (1/6m²){P′[(η²−1)(ν²−3a²+6a(2Δ_φ−d)−2) − (d−4)(d−2)] + ηP[(ν²−a(a+d−2))(7d/2+3a−1−6Δ_φ) − 2X − 2(a+(d−2)/2)(1−a²−a(d−2)/2)]} + O(m⁻⁴)`,
  printed with `X = (Δ_φ−d−1)_3`; **the correct value is `X = (2Δ_φ−d−1)_3` [DV]** (§1.3). Structure: the coefficient of
  `P^{(n)}/m^{2j}` grows at most like `η^{j+n}`, and like `ν^{2j}` near `η = 1`; so the expansion parameter is `η_AdS/m²` and `ν²/m²`.
* Origin of the harmonic function (3.4)–(3.18): Lorentzian shadow representation of the block, gauge-fixed to
  the frame `(x_1,x_2,x_3,x_4;x_5) = (−y, x/x², y, −x/x²; e)` (3.23), in which `r = |x||y|`, `η = −x·y/(|x||y|)`; the product of
  three-point functions becomes `𝒯(x,y) → e^{2p·(x+ȳ)}/(x²y²)^{d−Δ_φ}`, `p = (Δ, J e_z)` (C.8); and (3.18)/(B.2)
  `∫_{x>0} d^dx |x|^{a−d} P_{(2−d)/2+iν}(−x·y/|x||y|) e^{x·p} = 2^{a−1}π^{(d−2)/2} γ_a(ν) |p|^{−a} P_{(2−d)/2+iν}(−p·y/|p||y|)`.
  **This identity is the core analytic step** (the harmonic function reproduces itself; all `ν`-dependence is `γ_a`).
* Smearing (5.26)–(5.29), (5.41)–(5.44): `ω_AdS = ∫_0^M (dν/M) f(ν/M) C^{imp}_{2,ν}`; for `m ≫ M`, `β ≫ 1/M` the `ν`-integral is
  dominated by its endpoints; with `f ∼ f_0 p^n` at `0` and `f ∼ (1−p)^l`, `l > n+1−d/2`, at `1`:
  `m⁴ω_AdS → (2f_0/M^{n+1}) Q_n(β)`, `Q_n(β) = ∫_0^∞ dν ν^n P_{(2−d)/2+iν}(cosh β)` (analytically continued, (B.15)), plus an oscillatory
  term from `ν = M`. For `f = Σ a_n p^n`: `Σ a_n Q_n/M^{n+1}`; signs in `d = 3`: `POSITIVITY_LEMMA.md`.

### 1.2 Independent computation: what was implemented

* **Blocks [DV]** (`blocks3d.py`, `s0_block_checks.py`). `g_{Δ,J} = (4r_s)^Δ Σ_{n even} r_s^n Σ_j B_{n,j} P_j(η_s)` (`s`-channel radial
  coordinates), `B_{0,J} = J!/(1/2)_J`. The Casimir operator in `(r,η)`, `d = 3`, multiplied by `(1−r²)((1+r²)²−4r²η²)`,
  gives a recursion of depth three in `r²` and band `j → j, j±2` (formulas in the header of `blocks3d.py`; derived
  by hand). Checks: the block satisfies the Casimir equation in `(z,z̄)` to `10⁻¹⁵`–`10⁻¹⁸`; float64 and an independent
  mpmath implementation agree to `10⁻¹⁴`; normalisation `g/(z^{(Δ−J)/2}z̄^{(Δ+J)/2}) → 1`; convergence for `|ρ_s| → 1`.
* **Change of variables [D].** `dr dη = |ρ−ρ̄| dρ dρ̄/(4r²)`, each `(r,η)` covered twice by `(ρ,ρ̄) ∈ (0,1)²`, hence
  `Π̂_{k,ν}[G^s]/(2sin²) = ¼ ∫_0^1∫_0^1 dρ dρ̄ |ρ−ρ̄|^{d−2}(ρρ̄)^{(k−d−1)/2} P_{(2−d)/2+iν}(η) G^s_{Δ,J}`
  (agrees with the kernel `Ω` of (3.16) after `dw dw̄ |w−w̄|^{d−2}`, `w = −(1−ρ)²/(4ρ)`). For `w, w̄ < 0` the block is needed at
  `ρ_s = −((1−√ρ)/(1+√ρ))²`: the Regge limit `ρ → 0` is `ρ_s → −1`, the boundary of convergence of the radial series, and the
  integrand peaks at `1−|ρ_s| ∼ (2Δ_φ+k)/Δ`. Up to `N = 5400` radial levels were used (`Δ = 120`).
* **Consistency of the quoted normalisations [DV]** (`chmrs_formulas.py`): `4^{2Δ_φ−1}π^{d−2}/(b_{Δ,J} q_{Δ,J} m^{2(2Δ_φ−1)}) → 2/a_{Δ_φ}` with
  relative error `1.5·10⁻⁴, 1.3·10⁻⁶, 1.3·10⁻⁸` at `Δ = 40, 400, 4000` (`d = 3`, `Δ_φ = 1.2`): `1/m²`, as it must be. (3.45) to 15
  digits. (B.2) in `d = 3` by direct two-dimensional integration with `p` boosted: 12 digits (`s2_tensor_laplace.py`, part 0).

### 1.3 Result (`s0_regge_moment.py`, `out_s0_regge_moment.txt`)

`R_n :=` (exact heavy density of `Π̂_{k,ν}`) `/` ((3.20) with the corrections of the paper's ancillary file through `1/m^{2n}`);
`R_1′`: the same with the `1/m²` term of (3.27) as printed; `S_2`: exact heavy density of `C_{k,ν}` via (C.16) over the
ancillary series through `1/m⁴`. `d = 3`, `Δ_φ = 3`, `k = 2`:

| `Δ, J` | `m` | `β` | `ν` | `R_0` | `R_1` | `R_2` | `m⁶(R_2−1)` | `R_1′` | `S_2` |
|---|---|---|---|---|---|---|---|---|---|
| 20.3, 6 | 17.64 | 0.721 | 0.7 | 1.03848 | 1.00711 | 0.99951 | −1.5e4 | 0.97616 | 0.99867 |
| 40.3, 12 | 36.73 | 0.668 | 0.7 | 1.00730 | 1.00038 | 0.999992 | −1.9e4 | 0.99318 | 0.999980 |
| 80.3, 24 | 74.90 | 0.643 | 0.7 | 1.00166 | 1.000022 | 0.99999989 | −1.9e4 | 0.99830 | 0.9999997 |
| 120.3, 36 | 113.05 | 0.635 | 0.7 | 1.00072 | 1.0000042 | 0.99999999 | −1.9e4 | 0.99925 | 0.99999998 |
| 120.3, 36 | 113.05 | 0.635 | 4.0 | 1.00511 | 1.000021 | 1.0000001 | 2.2e5 | 0.99927 | 1.0000005 |
| 40.3, 2 | 38.72 | 0.129 | 0.7 | 1.00511 | 1.00020 | 0.999998 | −8.3e3 | 0.99487 | 0.999993 |
| 80.3, 60 | 50.49 | 2.030 | 0.7 | 1.01440 | 1.00091 | 0.99995 | −7.7e5 | 0.98906 | 0.99996 |

(also `k = 4`, `ν = 0, 2`, and `Δ_φ = 2`, in the output file: same pattern; `P` changes sign between `ν = 2` and `4`, and the
formula tracks it). Conclusions: (2.59)/(3.20) **verified** with its normalisation; the corrections are `1/m²` with the
ancillary-file coefficient; the printed (3.27) is off at `O(1/m²)` by the Pochhammer misprint (at `Δ_φ = 2`, `d = 3` both
versions vanish and `R_1′ = R_1`); corrections grow with `β` and `ν` as stated.

## 2. Stage 1 — heavy spinning blocks in `d = 3`

* **[Q]** (Kravchuk–Simmons-Duffin, arXiv:1612.08987; Erramilli et al., arXiv:1907.11247, 2011.01959; from memory, not
  re-read in this session.) Three-point structures in the collinear frame are labelled by the SO(2) charges `[q_1 q_2 q_3]`,
  `q_1+q_2+q_3 = 0`; four-point structures by `[q_1q_2q_3q_4]`; blocks are functions `g^{ab}_{[q_1q_2q_3q_4]}(r,θ)`.
* **[D]** (representation theory of the leading OPE term.) The primary contributes `r^Δ` times the matrix element
  of a rotation by `θ` in the spin-`J` representation between the states selected by the two three-point structures:
  `g^{ab}_{[q_1q_2q_3q_4]} = (4r)^Δ [n^a_{q_1q_2} n^b_{q_3q_4} d^J_{−q_1−q_2, q_3+q_4}(θ) + O(r)]`,
  the Wigner function replacing `P_J(cos θ)` (signs and phases of the labels are convention dependent; `n^a` are
  structure normalisations). The labels of `d^J` are the **`s`-channel** charge sums `q_1+q_2` and `q_3+q_4`.
* **[Q]+[DV]** scalar large-`Δ` form at fixed `r < 1`: `g_{Δ,J} ≈ (4r)^Δ N_J P_J(η)/[(1−r²)^{1/2}((1+r²)²−4r²η²)^{1/2}]`, `N_J = J!/(1/2)_J`
  (`s1_large_delta.py`): exact/asymptotic `= 0.9906, 0.9974, 0.99932, 0.99983` at `Δ = 10, 40, 160, 640` (`ρ,ρ̄ = 0.35, 0.2`), also for
  `J = 0.4Δ`. It **fails** in the region that controls dispersive functionals, `1−|ρ_s| ∼ 1/Δ`: there the ratio tends to
  `≈ 0.88 ≠ 1` (`0.914, 0.898, 0.890, 0.886` at `Δ = 20, 40, 80, 160`; `1−|ρ| = 6/Δ, 3/Δ`, `J = 0.3Δ`). For spinning blocks the analogous prefactor (a matrix in the `q` labels) was
  not reconstructed here; it is not needed below, for the same reason.
* **[Q]** CLS-D §4.1 ("Conservation at large `Δ`"): for a conserved operator and `Δ_O → ∞` (fixed `J`) conservation becomes
  algebraic, `v^μ_{1,23}⟨O_{1;μ…} O_2 O⟩ = 0`, `v_{1,23} = ½∂_1 log(x_12²/x_13²)`.
  **[DV]** (`s1_large_delta.py`, part b): in the collinear frame, `l = 6` fixed, the conserved `⟨JφO_l⟩` and `⟨TφO_l⟩` structures have
  (components with a time index)`/`(charge `±l_1` part) `= 0.38, 0.035, 0.0035, 0.00035` resp. `0.48, 0.040, 0.0040, 0.00040` at
  `Δ = 20, 200, 2000, 20000`. So for conserved operators at large `Δ` only `q_i = ±l_i` survive: **`q_i` is the helicity**, and the
  leading block components are `d^J_{λ,μ}(θ)` with `λ = −(q_1+q_2)`, `μ = q_3+q_4 ∈ {0, ±4}` for `⟨TTTT⟩`: `d^J_{0,0}` (`++→++`), `d^J_{4,4}`
  (`+−→+−`), `d^J_{4,−4}`, `d^J_{4,0}`, `d^J_{0,0}` (`++→−−`) — the list of CLPS (2.12)–(2.19).
* **[Q]** CHL (equation numbers not extracted; LaTeX labels given): in `d = 3` the helicity of a conserved current is
  conformally invariant (`h J^μ(p) = −iε^{μνσ}p_σ J_ν(p)/|p|`, §2.2, label `h in p`); helicity three-point structures `T^{±±}_{123}`
  (§2.1, boxed, label `Tpm`); the pairing (label `pairing general`) and the shadow maps (labels `Seasy`, `hard shadow`) are
  diagonal in helicity; mean-field OPE data of conserved currents and stress tensors are diagonal (§3.2, label
  `evenOPEcoeff`; the factor 2 from spurious poles is called there a "tentative guess", justified by the Lorentzian
  computation of their §3.3). Their "same"/"opp" structures exist for `J ≥ |J_1−J_2|` resp. `J ≥ J_1+J_2`.
* **Not done:** the requested check of a spinning block through differential operators acting on scalar blocks
  (needs scalar blocks with unequal external dimensions and the operators `D_{11}, D_{12}` in a conformal frame).

## 3. Stage 2a — heavy limit of spinning three-point functions in the Regge regime

### 3.1 Claim

Frame (3.23); `n = (1, e_z)`; `p = (Δ, J e_z)`, `p̄ = (Δ, −J e_z)`; `p̂, \hat{p̄} ∈ H²` (unit timelike), `−p̂·\hat{p̄} = cosh β`; at each point the
tangent frame `n_∥` (along the geodesic, parallel transported: `n_∥(\hat{p̄}) = (−J, Δ, 0)/m`, `n_∥(p̂) = (J, Δ, 0)/m`), `e_⊥ = e_y`,
`ε_± = (n_∥ ± i e_⊥)/√2`. Operators 1, 3 sit at `∓y` and are paired with `p̄` (`𝒯 ∼ e^{2p·x + 2p̄·y}`); operators 2, 4 at `±x/x²` with `p`;
their polarisations are referred to the origin by the inversion, `z_2 = I(x_2) z_2′`.

**Lemma V [D + DV].** Let `O_1` be conserved of spin `l_1` (`Δ_1 = l_1+1`). For `Δ, l → ∞` with `l/Δ` fixed and `x = X/m`, `y = Y/m`,
`⟨O_1(x_1,z_1) O_2(x_2,z_2) O(x_5,n)⟩ = (scalar structure) × [τ(z_1,z_2′) + O(1/m)]`, where `τ` is independent of `X, Y` and, as a
polynomial in `z_1`, is a combination of `(z_1·ε_+)^{l_1}` and `(z_1·ε_−)^{l_1}` at `\hat{p̄}`: transverse, helicity `±l_1`. The same holds for
`O_2` at `p̂`, and for the structures of the Lorentzian shadow `(3−Δ, −1−l)` at `(x_3,x_4)`.
*Argument.* The scalar factor is `e^{S}` with `S = O(m)`; `V_i`, `H_{ij}`, `H_{i5}/V_5` are `O(1)` and vary by `O(1/m)` over the region
`|x|,|y| ∼ 1/m`. Conservation `∂_{x_1}·D_{z_1}(e^S τ) = 0` gives at order `m`: `(∂_{x_1}S)·D_{z_1}τ = 0`, and `∂_{x_1}S ∝ p̄ + O(1)` (the
exponent is `p̄·y`-linear). A symmetric traceless tensor transverse to the timelike `\hat{p̄}` in three dimensions has only the
two components `ε_±^{⊗l_1}`. This extends the statement of CLS-D §4.1 from fixed `J` to the Regge regime, where `V_5^J`
contributes to `∂S`. ∎

### 3.2 Numerical verification (`threept.py`, `s2_three_point_heavy.py`, `s2_TTO_heavy.py`)

Embedding-space structures `V_1^{m_1}V_2^{m_2}V_5^{m_5}H_{12}^{n_{12}}H_{15}^{n_{15}}H_{25}^{n_{25}}` (and `ε_{15}`-structures), conservation imposed
numerically (divergence by high-precision differentiation, Todorov operator for spin 2; null spaces by SVD at 50–60
digits), `l/Δ = 0.4`, `X = (1.3, 0.4, −0.5)`, `Y = (0.9, −0.3, 0.2)`:

| vertex | quantity | `Δ = 50` | `500` | `5000` |
|---|---|---|---|---|
| `⟨JφO⟩` even | `t·\hat{p̄}` (unit `t`); direction | 0.040; `n_∥` | 0.0038 | 0.00038 |
| `⟨JφO⟩` odd | `t·\hat{p̄}`; direction | 6·10⁻⁴; `e_⊥` | 5·10⁻⁵ | 5·10⁻⁶ |
| `⟨TφO⟩` even | form mod `z²`: `q_nn : q_ee`; `q_pn` | 0.978 : −1; 0.064 | 0.9972 : −1; 0.0062 | 0.99971 : −1; 0.00062 |
| `⟨TφO⟩` odd | `q_ne`; `q_pe` | 1; 0.065 | 1; 0.0062 | 1; 0.00062 |
| shadow `⟨O^S Jφ⟩`, `⟨O^S Tφ⟩` | same pattern | `t·\hat{p̄} = 0.015`, `q_pn = 0.0096` | | `−1.5·10⁻⁴`, `−9·10⁻⁵` |

`⟨TTO_l⟩`, parity even, conservation at both points, no Bose symmetry imposed (`Δ = 60, 600, 6000`): the conservation
system on the 14 monomials has exactly three null vectors (singular values `10⁻⁵⁹`), one of which is the `d = 3`
linear relation among monomials: **two conserved structures**, as in arXiv:1708.05718. With `z_1(θ_1) = \hat{p̄} + cos θ_1 n_∥ +
sin θ_1 e_⊥`, `z_2′(θ_2)` likewise at `p̂`, the vertex contains only the harmonics `e^{±2iθ_1}e^{±2iθ_2}`: all others `/` total
`= 0.176, 0.0173, 0.00172`. The coefficients satisfy `(2,2) = \overline{(−2,−2)}`, `(2,−2) = \overline{(−2,2)}` (parity) and the two
structures span both orbits (determinant `0.84`).

### 3.3 Dictionary of labels [D]

`e^{∓2iθ}` at a point `↔` tensor `ε_±²`. Total charge under a common rotation of both transverse frames (which becomes
the collinear-frame SO(2) when `β → 0`): harmonics `(2,−2), (−2,2)` have charge `0` — the structure exists for all even `J`
— this is `c^{++}` (`λ = h_1−h_2 = 0`); harmonics `(2,2), (−2,−2)` have charge `±4` — this is `c^{+−}` (`J ≥ 4`). So

`τ_L = c^{++}[ε_+²(\hat{p̄})ε_−²(p̂) + c.c.] + c^{+−}[ε_+²(\hat{p̄})ε_+²(p̂) + c.c.]` (parity-even theory),

and with the right vertex its conjugate (helicity-diagonal pairing, **[Q]** CHL), the block tensor at the pair `(1,3)` is
`ε_{h_1}² ⊗ \bar ε_{h_3}²`: SO(2) weight `h_1 − h_3 ∈ {0, ±4}`, i.e. the **helicity transfer of the Regge pair**; likewise at `(2,4)`.
Helicity-preserving terms `|c^{++}|²`, `|c^{+−}|²`: weights `(0,0)`. `(c^{++})²`: `ε_+⁴` at `\hat{p̄}`, `ε_−⁴` at `p̂` — opposite, i.e. `Φ^ν_{4,4}` in the
convention of Theorem R of `REGGE_FROM_BLOCKS.md`. `|c^{+−}|²` with both flipped: `ε_+⁴`, `ε_+⁴`, i.e. `Φ^ν_{4,−4}`. `c^{++}c^{+−}`: `(0,4)`.
This derives the table of `REGGE_FROM_BLOCKS.md` §5 from the three-point structures (its status there was [I]).

## 4. Stage 2b — the transform to the transverse harmonic basis

### 4.1 Set-up

A functional acts on a block through `∫d³x d³y K_I(x,y) 𝒯^I(x,y)` (frame (3.23), `I` = tensor indices of the four
operators, 2 and 4 inverted). Its **Regge-moment expansion** is the expansion of `K` in powers of `|x||y|`; the leading
term is homogeneous, `K_I = (|x||y|)^{a−3} 𝒦_I(x̂,ŷ)`, and `𝒦` is decomposed in tensor harmonics on `H² × H²`:
`𝒦 = Σ κ_{ss′}(ν) 𝔻^{(13)}_s(ŷ) 𝔻^{(24)}_{s′}(x̂) P_{−1/2+iν}(−x̂·ŷ)`. For scalars `𝔻 = 1` and `a = 2Δ_φ+k−1` (from (3.25) and (C.8)). By Lemma V
`𝒯^I → e^{2p·x+2p̄·y} τ^I`, so the leading heavy action is the Fourier–Laplace transform of `𝒦_I` contracted with a
**constant, transverse** tensor.

Structures at one pair (polarisations `z_1, z_3` at `ŷ`): CHP (2.41) `D_1 = z_1·z_3 + (z_1·ŷ)(z_3·ŷ)`, `D_2 = −(z_1·ŷ)(z_3·ŷ)`,
`D_3 = (z_1·ŷ)(z_3·∇) + (z_3·ŷ)(z_1·∇)`, `D_4 = −(z_1·∇)(z_3·∇) + ½D_1∇²`; parity odd `D_5 = ε(z_1,z_3,ŷ)`. For stress tensors: `D_1²`, `D_1D_2`,
`D_2²`, … and `𝔄 := D_5 D_1`.

### 4.2 Lemma L (tensor Fourier–Laplace transform) [DV]

Write `Φ_a(y) = |y|^{a−3}P_{−1/2+iν}(−ŷ·x̂)`, `C_a = 2^{a−1}√π`, `γ_a(ν) = |Γ((a−½+iν)/2)|²`, `P^{(n)} := sinh^nξ d^nP/d(cosh ξ)^n`, `ξ = dist(x̂,p̂)`,
`ε_±` at `p̂` built on the geodesic from `x̂`. Then

`∫_{y>0} d³y e^{y·p} [𝔻_s Φ_a](ε_+, ε_−) = C_a γ_a(ν) |p|^{−a} M⁰_s(ν;a) P_{−1/2+iν}(cosh ξ)` (helicity preserved),
`∫_{y>0} d³y e^{y·p} [𝔻_s Φ_a](ε_+, ε_+) = C_a γ_a(ν) |p|^{−a} M^{flip}_s(ν;a) P^{(2)}` (currents) resp. `P^{(4)}` (stress tensors),

with (`𝖽_b := ν² + (a−b)²`, `g_1 := γ_{a−1}(ν)/γ_a(ν)`):

| structure | `M⁰` (helicity preserved) | `M^{flip}` |
|---|---|---|
| `D_1` | `(ν² + 2(a−2)² + ¼)/(2𝖽_{5/2})` | `1/(2𝖽_{5/2})` |
| `D_2` | `(ν² − 2a + 17/4)/(2𝖽_{5/2})` | `−1/(2𝖽_{5/2})` |
| `D_3` | `−(a−1)(ν²+¼)/𝖽_{5/2}` | `(a−3)/𝖽_{5/2}` |
| `D_4` | `−(ν²+¼)(ν²+9/4)/(4𝖽_{5/2})` | `−(ν² + 2(a−3)² + ¼)/(4𝖽_{5/2})` |
| `D_5` (odd) | `∓ i (a−1) g_1/2` for helicity `±` | `0` |
| `D_1²` | `N_1/(64 𝖽_{5/2}𝖽_{9/2})` | `1/(4𝖽_{5/2}𝖽_{9/2})` |
| `D_1D_2` | `−N_2/(64 𝖽_{5/2}𝖽_{9/2})` | `−1/(4𝖽_{5/2}𝖽_{9/2})` |
| `D_2²` | `N_3/(64 𝖽_{5/2}𝖽_{9/2})` | `1/(4𝖽_{5/2}𝖽_{9/2})` |
| `𝔄 = D_5D_1` (odd) | `∓ i g_1 (a−1)(ν² + 2(a−3)² + ¼)/(4𝖽_{7/2})` | `0` |

`N_1 = 64a⁴ − 768a³ + 64a²ν² + 3280a² − 320aν² − 5840a + 16ν⁴ + 488ν² + 3705`,
`N_2 = 64a³ − 32a²ν² − 584a² + 224aν² + 1720a − 16ν⁴ − 424ν² − 1641`, `N_3 = 128a² − 128aν² − 800a + 16ν⁴ + 360ν² + 1113`.
(`M⁰_{D_1} + M⁰_{D_2} = 1`, as it must for `z_1·z_3`.)

*Proof.* Each `𝔻_sΦ_a` is a polynomial in `y` times flat derivatives of `Φ_b`, e.g. `D_1Φ_a = (z_1·z_3)Φ_a + (z_1·y)(z_3·y)Φ_{a−2}`,
`D_3Φ_a = [(z_1·y)(z_3·∂) + (z_3·y)(z_1·∂)]Φ_a + 2(a−3)(z_1·y)(z_3·y)Φ_{a−2}`,
`|y|^{a−3}(z_1·∇)(z_3·∇)P = (z_1·∂)(z_3·∂)Φ_{a+2} + (a−2)[(z_1·y)(z_3·∂)+(z_3·y)(z_1·∂)]Φ_a + (a−1)(z_1·z_3)Φ_a + (a−3)²(z_1·y)(z_3·y)Φ_{a−2}`
(checked against the definition with tangential projections, exactly, at a random point). Under the transform `y^μ → ∂/∂p_μ`,
`∂_y → −p`, `Φ_b → C_bγ_b|p|^{−b}P(−p̂·x̂)` by (B.2); the rest is differentiation and Legendre's equation (sympy, rational in
`tanh(ξ/2)`). ∎ Numerical check by direct integration of the tensor kernels (`a = 4.5`, `ν = 1.1`, `ξ = 0.9`): `D_1, D_2, D_3, D_5` agree
to 9–10 digits in both columns, `D_4` to `10⁻⁵` (quadrature-limited). The stress-tensor rows are obtained by the same
code path and were **not** separately integrated numerically. (`s2_tensor_laplace.py`, `s2_tensor_laplace_odd.py`.)

Remarks. (i) The transform is SO(2,1)-covariant but **not local on `H²`**: a kernel that is pure weight 0 at `ŷ` (`D_1`) feeds the
flip component at `p̂` and vice versa (`D_4` feeds weight 0). Only the weight **at `p̂`**, fixed by the block's transverse tensor,
is meaningful; the kernel's own weight at `ŷ` is not. (ii) All factors are rational in `ν` apart from `γ_a`, `γ_{a−1}`, with poles at `ν = ±i(a−5/2), ±i(a−9/2)` (`±i(a−7/2)` for `𝔄`);
`M⁰_{D_1} > 0`, `→ ½` as `ν → ∞`; `M⁰_{D_1²} → ¼` and is positive for real `ν` in the cases evaluated (`a = 6, 7, 9`). (iii) For two pairs the transforms in `y` and `x` commute with the derivative operators of
the other point, so the factors multiply and, by Theorem R of `REGGE_FROM_BLOCKS.md`, the harmonic function of the
component with helicity transfers `(m_L, m_R)` is `Π_{k<|m_L|}…Π_{k<|m_R|}… × Φ^ν_{m_L,±m_R}(β)`.

### 4.3 Leading heavy action

**[D]** (from Lemma V, Lemma L, and the scalar reduction (3.25)–(C.8)) For a kernel with leading Regge moment
`(|x||y|)^{a−3} Σ κ_{ss′}𝔻_s𝔻_{s′}P_{−1/2+iν}` acting on the block with helicity structures `(h_1h_2)`, `(h_3h_4)`:

`ω[G^{(h_1h_2),(h_3h_4)}_{Δ,J}] = 𝒩_{Δ,J} · (π/4) γ_a(ν)² m^{−2a} · Σ_{ss′} κ_{ss′}(ν) M_s^{(h_1,h_3)}(ν;a) M_{s′}^{(h_2,h_4)}(ν;a) · 𝖧_{m_L,m_R}(β) · (1 + O(1/m²))`,

`𝖧_{0,0} = P_{−1/2+iν}(cosh β)`; one pair flipped: `𝖧 ∝ P^{(4)}`; both flipped: `𝖧 ∝ Π_{k<4}(ν²+(k+½)²) Φ^ν_{4,±4}(β)` (Theorem R); `𝒩_{Δ,J}` is the `ν`-independent
normalisation of the shadow representation of the spinning block (for scalars `4^{2Δ_φ}(−1)^J/b_{Δ,J}`). In particular:

* `++ → ++` and `+− → +−`: `𝖧 = P_{−1/2+iν}(cosh β)`, no other dependence on `β`. For a parity-even kernel `𝔖⊗𝔖`, `𝔖 = D_1²` (or any
  combination of `D_1², D_1D_2, D_2², …`) the factor `M_𝔖(ν;a)²` is common to both: the functional sees `ρ^{++} + ρ^{+−}`.
* For `𝔄⊗𝔄` the factors are `(∓i…)(∓i…)`: sign `−sgn(h_1)sgn(h_2)` in the parallel frame, i.e. `+` for `++`, `−` for `+−`: the
  functional sees `ρ^{++} − ρ^{+−}` (the pattern of CLPS (2.33) versus (2.35a), and of CLS-D's even/odd `j′`).
* `𝔖⊗𝔄` sees `Σ sgn(h_2)|c^{h_1h_2}|²`, which vanishes in a parity-invariant theory.

## 5. Stage 4 — statement

**Proposition (leading Regge-moment formula for `⟨TTTT⟩` in `d = 3`; conditional).** Let `O` be a primary in `T × T` with
`(Δ, J)`, `m² = αᾱ`, `cosh β = η_AdS` as in (2.28)–(2.30) at `d = 3`. Assume:

* (H1) the Lorentzian shadow representation of spinning blocks, with the right three-point structure the shadow
  continuation of the left one and helicity-diagonal normalisation **[Q for scalars, CHMRS (3.4); I for `⟨TTTT⟩`, from
  Costa et al. differential operators and CHL's helicity-diagonal pairing and shadow maps (labels `pairing general`,
  `Seasy`, `hard shadow`)]**;
* (H2) `ω` is a swappable dispersive functional whose kernel, in the frame (3.23), has a Regge-moment expansion
  with leading term `(|x||y|)^{a−3}[κ_+(ν) 𝔖⊗𝔖 + κ_−(ν) 𝔄⊗𝔄] P_{−1/2+iν}(−x̂·ŷ)` plus structures of non-zero helicity transfer
  **[I: no spinning functional has been expanded in this form]**;
* (H3) `m → ∞` at fixed `β > 0` and fixed `ν`.

Then, with `ρ^{h}_{Δ,J} ≥ 0` the helicity-diagonal block coefficients (`|λ^{h}_{TTO}|²` in the normalisation of (H1), each summed
over the two parity images; the overall sign of the `κ_−` term is a convention for `𝔄`),

`ω|_{hel.-pres.}[Δ,J] = 𝒩_{Δ,J} (π/4) γ_a(ν)² m^{−2a} P_{−1/2+iν}(cosh β) { κ_+ M_𝔖(ν;a)² (ρ^{++}+ρ^{+−}) + κ_− |M_𝔄(ν;a)|² (ρ^{++}−ρ^{+−}) } (1 + O(1/m²))`,

`M_𝔖 = N_1/(64𝖽_{5/2}𝖽_{9/2})`, `|M_𝔄| = g_1(a−1)(ν²+2(a−3)²+¼)/(4𝖽_{7/2})` (table of §4.2).

Status of the ingredients: plane-wave form of the scalar factor **[Q, (C.8); DV through §1.3]**; transversality and
helicity content of the vertices **[D; DV for `⟨JφO⟩`, `⟨TφO⟩` and their shadows, and for `⟨TTO⟩` (left vertex only)]**; `P_{−1/2+iν}` and the factors `M` **[DV]**; relative
sign in the `𝔄⊗𝔄` term **[D]**; `O(1/m²)` rather than `O(1/m)` **[D; DV pointwise for currents, §6.1]**; `𝒩_{Δ,J}` and the pairing **[I]**.

**Normalised form (proposal) [I for the normalisation].** Define, by analogy with (2.56),
`C^±_{k,ν} := a_{Δ_T}/(γ_a(ν)² M_±(ν;a)²) × (the functional with pure 𝔖⊗𝔖, resp. 𝔄⊗𝔄, leading Regge moment)`, `a = 2Δ_T+k−1 = k+5`,
and heavy densities by dividing by the mean-field data in the helicity basis (CHL §3.2, label `evenOPEcoeff`) as in
(2.31)–(2.32). Then

`C^±_{k,ν}[Δ,J] = (2P_{−1/2+iν}(η_AdS)/m^{2k}) (|c^{++}_{Δ,J}|² ± |c^{+−}_{Δ,J}|²) (1 + O(ν²/m², 1/m²))`, `|c^{h}|² := ρ^{h}/ρ^{h}_{MFT}`,

i.e. exactly the scalar formula (2.59) with `p/p^{MFT}` replaced by `|c^{++}|² ± |c^{+−}|²` (`k` even for `+`, odd `≥ 3` for `−`: [Q] CLS-D
§3.5.2, subtraction factors `f^±_{k;J_e}`). In terms of CFT data: `ρ^{h_1h_2} = (λ^{h_1h_2}_{TTO})²/c_T²`, with `λ^{h_1h_2}` the coefficient of CHL's
helicity structure `T^{h_1h_2}_{TTO}` (a fixed linear combination of the parity-even/odd `λ^{(a)}_{TTO}`: `λ^{±±} ∝ λ^{E,same} ± λ^{O,same}`,
`λ^{±∓} ∝ λ^{E,opp} ± λ^{O,opp}`) for unit-normalised `O` and `⟨TT⟩ ∝ c_T`; `ρ ∼ 1/c_T` for single-trace `O`; `++ ≡` "same", `+− ≡` "opp".
The relation to the numbers `a^{++}_O`, `a^{+−}_O` of `KERNELS_D3.md` was not worked out.
**What is and is not needed for positivity:** `q`-type normalisations are positive and drop out; what matters is the
shape `P_{−1/2+iν}(cosh β)`, the factor `γ_a² M²` (known exactly once the kernel is known), and the relative weight of
different functionals on the same block (fixed by their kernels, not by the three-point normalisation).

**Consequence for the uplift of CLPS (3.3).** If the CFT functionals uplifting `B^{(1)}_{2}`, `B^{(1)}_3` are of the types `C^+_2`, `C^−_3`
(**[I]**; this is the identification in item 2 of §7), then for `m → ∞` at fixed `β ≫ 1/M`
`F^{±} → (2/m⁴)[a_1Q_1(β)/M² + a_2Q_2(β)/M³](|c^{++}|²+|c^{+−}|²) + (1/m⁶)(…)(|c^{++}|²−|c^{+−}|²) + …`, with the **scalar** moments `Q_n`: the
positivity input is Lemma Q, as inferred in `REGGE_FROM_BLOCKS.md` §7.3; the kernel statement behind it is now
derived, its normalisation is not.

## 6. Stage 3 — subleading terms

### 6.1 Order of the first correction [D + DV pointwise]

Let `F(x,y;Δ,l)` be the left vertex of an analytic family of conserved structures. The right vertex (operators at
`(x_3,x_4) = (y,−x/x²)`, shadow quantum numbers) is `F(−x,−y; 3−Δ, −1−l)`. Hence under
`(x,y;α,ᾱ) → (−x,−y;−α,−ᾱ)` the tensor `τ_L⊗τ_R` goes to its transpose (`1↔3`, `2↔4`), while the scalar factor and the
measure are invariant (CHMRS, below (3.26)). With `x = X/m`: the `1↔3`-symmetric part of the integrand is even in `m`
at fixed `(X, Y, η_AdS)`. So **for kernels symmetric under transposition (`𝔖⊗𝔖`, `𝔄⊗𝔄`) the expansion is in `1/m²`**, in the
variables `(m², η_AdS)` and with frames built from `(α,ᾱ)`: `tanh(β/2) = (2J+1)/(2Δ−3)`.
Check (`s3_weyl_even.py`, `⟨JφO⟩ ⊗ ⟨O^S Jφ⟩`, symmetric part `W` of `t_L⊗t_R`, `l/Δ ≈ 0.4`): each vertex has non-transverse
components of order `1/m` (§3.2), but in the `(α,ᾱ)` frame `m²W_{pn} = 0.317, 0.326, 0.331, 0.333, 0.335` and
`m²W_{ne} = −0.4244, −0.4241, −0.4240, −0.4240, −0.4241` at `m = 44, 90, 182, 365, 732`; in the frame built from `(Δ,−l)` one finds
instead `m W_{pn} → −1.20`. Not checked for `⟨TTO⟩`.

### 6.2 Structure and what is known of the coefficients

* **[I]** By the mechanism of CHMRS §3.3 (polynomials in `x, y` times the plane wave, integrated by `∂_p`), the `1/m²` term
  is `(1/m²)[𝒜(η,ν) P + ℬ(η,ν) P′]` with `𝒜, ℬ` polynomial in `η` (degree `≤ 1, 2`) and in `ν²` (degree `≤ 1`), now depending on the
  kernel structure and on the helicity channel. For scalars it is (3.27) (with the correction of §1). **Not derived for
  stress tensors**: it needs the second-order expansion of the conserved `⟨TTO⟩` structures in the frame (3.23).
* **Helicity-dependent part, flat corner [DV in flat space, I in AdS].** For the flat partial waves
  `d̃^J_{λ,λ}(1−2p²/m²) = ₂F₁(λ−J, J+λ+1; 1; p²/m²)`, at fixed `b`, `J(J+1) = b²m²/4`:
  `d̃^J_{λ,λ} − d̃^J_{0,0} = (1/m²)[2λ²(p/b)J_1(bp) + λp²J_0(bp)] + O(m⁻⁴)` (derived from the series; `s3_flat_helicity_term.py`:
  `m²(d̃_{4,4}−d̃_{0,0}) = 4.5699062` at `m = 427`, `b = 3`, `p = 0.7`, prediction `4.5697976`; four `(b,p)` checked). By the bulk-point dictionary
  (CLS-D) and the overlap of regimes, the `+−` density of `C^+_{k,ν}` differs from the `++` one, for `1/m ≪ β ≪ 1`, by
  `(2/m^{2k})(1/m²)[32(ν/β)J_1(νβ) + 4ν²J_0(νβ)]` (the second term is present only if the functional uplifts the *stripped* partial
  wave `d̃`, as in CLPS (2.33); for the unstripped `d^J_{λ,λ}` only the `λ²` term remains). The natural completion at `β ∼ 1` is
  `(1/m²)[−2λ²P′_{−1/2+iν}(η_AdS) + λν²P]` up to terms that vanish as `β → 0`; **the `β`-dependence at `β ∼ 1` is not derived.** After smearing with `f = a_1p + …` the first piece
  gives `−(32a_1/(M²m²)) ∂_βQ_1/sinh β`, positive since `Q_1` decreases (tables of `REGGE_FROM_BLOCKS.md` §7.1): if the completion
  is right it favours positivity in the `+−` channel.
* **The term of opposite sign (`B_3`) [D at leading order].** It is not a correction to `C^+_2` but a different component:
  the `𝔄⊗𝔄` kernel, whose heavy action is `∝ P_{−1/2+iν}(cosh β)(|c^{++}|² − |c^{+−}|²)` with the same harmonic function and the
  factor `γ_a² g_1²(a−1)²(ν²+2(a−3)²+¼)²/(16𝖽_{7/2}²)`. In a composite functional it enters at relative order `M²/m²`, as in flat space.
* Parametric sizes: `1/m²`, `ν²/m²`, `η_AdS/m²`; in the flat corner the helicity-dependent term is `λ²/(mβ)² ∼ λ²/J²`.

## 7. What is missing for a theorem; exact blocking steps

1. **No end-to-end test on a spinning block.** The chain here is: scalar reduction tested on exact blocks (§1) +
   vertices tested pointwise (§3) + transforms tested (§4). The link not tested is (H1): the spinning block as the
   shadow integral with a helicity-diagonal pairing and its normalisation `𝒩_{Δ,J}`. *Blocking step:* an exact spinning
   block in the region `1−|ρ_s| ∼ 1/Δ`, e.g. `⟨JφJφ⟩` from the operators `D_{11}, D_{12}` acting on scalar blocks with unequal
   external dimensions (the radial recursion of `blocks3d.py` would have to be extended to `Δ_{12}, Δ_{34} ≠ 0`), and a
   tensor Regge moment of it.
2. **The functionals.** CLS-D's `Ψ^{(a),(c)}_{k,ν,j′}` are defined in general `d` with light-transformed operators 1, 3; their
   Regge-moment expansion (the analogue of (2.45)–(2.47), (C.13)), hence `κ_{ss′}(ν)` and `a`, was not computed, and their
   `d = 3` reduction (where `j′` is only a parity label) is not in the literature. Until then (H2) is an assumption, and
   "which CFT functional uplifts `B^{(1)}_2`, `B^{(1)}_3`" is inferred. The sum rules beyond the lowest-subtracted ones are
   missing altogether (CLS-D §5, Discussion).
3. **Normalisation relative to mean-field theory** in the helicity basis: needs the Lorentzian `b`-coefficients for
   CHL's structures and the large-`Δ` limit of CHL's mean-field data (their spurious-pole factor is itself a stated guess).
4. **`1/m²` coefficients** for `⟨TTTT⟩`: second-order expansion of the conserved structures (mechanical; not done), and
   the check of §6.1 for `⟨TTO⟩`.
5. **Uniformity.** In `β`: the scalar corrections grow like `(η_AdS/m²)^j`, and at `β = 2` the `1/m⁴`-accurate formula is already
   off by `4·10⁻³` at `m = 50`, `ν = 4`, `Δ_φ = 2` (output file); large `β` (large-spin, fixed-twist) needs a separate treatment. In `J`: `β → 0`
   overlaps the bulk-point regime only if the spinning analogue of (2.60) holds for all `ν < m` (CLS-D derive it for
   `ν, m ≫ 1` and "expect" the rest). In `ν`: the factors `M(ν;a)` have their poles off the real axis, but zeros of `M`
   on the real axis were excluded only in examples, and `C^±` inherit the `e^{πν}` growth of `1/γ²`, as in the scalar
   case (CHMRS App. A).
6. **Threshold region `m ∼ M`** (`ν` up to `M` comparable to `m`): the Regge formula does not apply; the scalar treatment
   (CHMRS App. F.1, the oscillatory endpoint) has no spinning counterpart.
7. **Convergence and swapping** for spinning functionals: Regge boundedness for `⟨TTTT⟩` and the subtraction factors
   are discussed in CLS-D §3.5.2; nothing was examined here. **Discrete series**: for flip components the transform
   produces `Φ^ν_{m,m′}` times polynomials vanishing at `ν = ±i(k+½)`; whether functionals see those points was not examined.
8. **The light side**: graviton exchange in `AdS_4` in each helicity component, double traces `[TT]`, contact terms.
9. Parity-odd theories (`|c^{++}| ≠ |c^{−−}|`): the `𝔖⊗𝔄` components do not vanish; not treated.

## 8. Code index (`ads4/code/reggemoment/`)

| script | content | output |
|---|---|---|
| `chmrs_formulas.py` | quoted formulas; `(3.19) → (3.20)`; (3.45) | `out_chmrs_formulas.txt` |
| `chmrs_series.py`, `data/chmrs_regge_moments_order8.m` | parser of the paper's ancillary series (unchanged file); misprint of (3.27) | `out_chmrs_series.txt` |
| `blocks3d.py`, `s0_block_checks.py` | exact `d = 3` scalar blocks (Casimir recursion in radial coordinates); validation | `out_s0_block_checks.txt` |
| `s0_regge_moment.py` | exact Regge moments of blocks vs (3.20), (3.27), (C.16) | `out_s0_regge_moment.txt` |
| `s1_large_delta.py` | scalar `h_∞`; conserved vertices in the collinear frame | `out_s1_large_delta.txt` |
| `threept.py`, `s2_three_point_heavy.py`, `s2_TTO_heavy.py` | embedding three-point structures, numerical conservation, heavy limit in the frame (3.23) | `out_s2_three_point_heavy.txt`, `out_s2_TTO_heavy.txt` |
| `s2_tensor_laplace.py`, `s2_tensor_laplace_odd.py` | Lemma L: symbolic derivation and numerical check | `out_s2_tensor_laplace*.txt` |
| `s3_weyl_even.py` | cancellation of `1/m` terms in the symmetric vertex tensor | `out_s3_weyl_even.txt` |
| `s3_flat_helicity_term.py` | flat-space helicity-dependent `1/m²` term | `out_s3_flat_helicity_term.txt` |

Run time: a few minutes in total (`s0_regge_moment.py` about 3 min, `s2_tensor_laplace.py` about 2 min).

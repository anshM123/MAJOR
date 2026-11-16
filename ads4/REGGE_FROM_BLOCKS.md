# Spinning Regge kernels in CFT3 from conformal Regge theory: derivation in the smallest sectors

Status: written 16 Nov 2026 by a derivation agent; **not refereed**. Code and recorded output:
`ads4/code/regge/` (`01`–`09_*.py`, `outputs.txt`, `out08.txt`). Sources (PDF text) are in the session scratchpad
`…/scratchpad/regge/` and `…/scratchpad/chmrs/p.txt`.

Every formula carries one of three marks:

* **[Q]** quoted from the literature, with reference and equation number;
* **[D]** derived here and verified symbolically (sympy) or numerically (mpmath, 12–28 digits) by the named script;
* **[I]** inferred (argument given), **not verified** against a block computation.

## 0. Summary of what is and is not established

1. **[D]** Starting from the literature's Regge-limit formula (differential operators acting on the harmonic
   function `Ω_{iν}(L)` of `H_{d−1}`), specialised to `d = 3` and decomposed into SO(2) weights about the
   transverse point: for `⟨JJOO⟩` the weight-0 (helicity-preserving) component is `β_1(ν) Ω_{iν}` and the weight-`±2`
   (helicity-flip) component is `−½ β_4(ν) (ν tanh πν/4π) P²_{−1/2+iν}(cosh L)`; for `⟨TTOO⟩` the weight-0 component is
   `β_1(ν) Ω_{iν}`, the weight-`±4` component is `¼ β_3(ν) (ν tanh πν/4π) P⁴_{−1/2+iν}(cosh L)`, **the operator `D̂_2`
   vanishes identically and `D̂_3` has no weight-0 component**. `P^m = (−1)^m Π_{k<m}(ν²+(k+½)²) P^{−m}`.
2. **[D]** Both pairs spinning: `m` helicity raisings at each end of `P_{−1/2+iν}(cosh L(x,y))` give exactly
   `2^{−m} Π_{k=0}^{m−1}(ν²+(k+½)²) · Φ^ν_{m,±m}(β)` (Theorem R, proved; checked for `m ≤ 4` to 24 digits). The factor is a
   polynomial, `|Γ(½+m+iν)/Γ(½+iν)|²`; it is neither `1` nor an irreducible ratio of Gamma functions. One-sided:
   `2^{−m/2}(½+iν)_m Φ^{−ν}_{m,0} = (−1)^m 2^{−m/2} P^m`.
3. **[D] + [Q] + [I]** The SO(2) weights are the helicity transfers of the *Regge pairs* (the two operators that
   are the same particle before and after the collision), not the `s`-channel labels of the Wigner functions.
   The corrected dictionary of `SPINNING_KERNELS.md` Part B (version after review) is **confirmed** in the
   smallest sectors; see §5 for the table and for which coupling combination multiplies which kernel.
4. **[Q] + [I]** Measure: in the scalar case the kernel is the unit-normalised `P_{−1/2+iν}(cosh β)` against `dν f(ν/M)`
   by construction of `C_{k,ν}` (2106.10274 (2.56), (2.59), (5.27), (5.41)–(5.43)). For spinning correlators the
   same construction is **not in the literature**; the normalisation that makes the Regge formula agree with
   the flat sum rule in the overlap region `β → 0` at *every* fixed `ν` is: unit-normalised `Φ_{0,0}` and `Φ_{4,4}`;
   `P^{−4}` for one flip; `Φ_{4,−4}/Π_{k<4}(ν²+(k+½)²)` for `+− → −+` (two flips, total `|ΔJ_z| = 8`) (§6).
5. **Consequences (§7).** The explicit functionals `F_{g3}, F_{g4}` of 2201.06602 (3.3) use only the fixed-`u` MHV sum
   rules `B^{(1)}`, whose heavy densities `|c^{++}|²` and `|c^{+−}|²` both have pair weights `(0,0)`. Their leading
   curved-regime action is therefore `∝ a_1 Q_1/M² + a_2 Q_2/M³` with the **scalar** moments, `a_1 = 65, a_2 = −40`
   (resp. `99, −15`), in both channels. The original Part C conclusion ("the first moment `1/24 − m²/2` is negative,
   so positivity does not uplift in the `+−` channel") does **not** stand; the retraction in the rewritten
   `SPINNING_KERNELS.md` is correct. `Φ_{4,4}` multiplies the interference term `(c^{++})²` of the all-plus sum
   rule, where only `A_{00} ≥ |A_{44}|` is required.
6. **Not derived:** the spinning heavy-block Regge-moment formula itself (the analogue of 2106.10274
   (3.17)–(3.20) and (2.59) for stress tensors in `d = 3`); everything marked [I]. List in §8.

## 1. Conventions

* `d = 3`, `h = d/2 = 3/2`. Transverse space `H² = H_{d−1}`, curvature `−1`; `L` (or `β`) the geodesic distance;
  `s := ½ + iν`, `λ := ν² + ¼`, `(∇² + λ) F = 0`, i.e. `F″ + coth L F′ + λF = 0` for radial `F`.
* `P^μ_{−1/2+iν}(x)`, `x > 1`: associated conical functions "of type 3" (mpmath `legenp(…, type=3)`), so that for integer
  `m ≥ 0`: `P^m = (x²−1)^{m/2} d^m P/dx^m = (sinh L · d/d cosh L)^m P`. **[D]** (`02`, `03`):
  `P^m = (−1)^m (½+iν)_m (½−iν)_m P^{−m} = (−1)^m Π_{k=0}^{m−1}(ν²+(k+½)²) P^{−m}`; `P^{−m}(cosh β) → (β/2)^m/m!` as `β → 0`.
* `Φ^ν_{m,n}(β)`: matrix coefficients of the unitary spherical principal series, as defined in
  `SPINNING_KERNELS.md` Part A. `Φ_{0,0} = P_{−1/2+iν}`, `Φ_{m,m} = cosh^{2m}(β/2) ₂F₁(½+m+iν, ½+m−iν; 1; −sinh²(β/2))`.
  **[D]** (`03`, numerically): `Φ^ν_{m,0} = (½+iν)_m P^{−m}`, `Φ^ν_{0,m} = (−1)^m (½−iν)_m P^{−m}`; `Φ_{m,±m}` are real and even in `ν`;
  flat limit (`β → 0`, `νβ` fixed) `Φ_{m,n} → J_{|m−n|}(νβ)` (checked for `(4,4), (4,−4), (4,0)`).
* Polarisations. A bulk massless particle moving along a null geodesic has its transverse polarisations
  tangent to `H²`. At a point with unit tangent `n` to the geodesic towards/away from the other impact point
  and `n_⊥` its rotation by `+π/2`: `ε_± := (n ± i n_⊥)/√2`, `ε_+·ε_− = 1`, `ε_±·ε_± = 0`. In `d = 3` every null complex
  tangent vector is a multiple of `ε_+` or `ε_−`. Frames at two points are identified by parallel transport
  along the geodesic. A structure has **weight `m`** at a point if it is multiplied by `e^{imφ}` (up to the
  orientation convention) under rotation of the frame at that point.
* Regge pair. Costa–Hansen–Penedones label the operators so that `O_1, O_2` (positions `x_1, x_2` in
  different Poincaré patches, `x ≈ x_1 − x_2` small) are the *same* bulk particle before and after the
  interaction, likewise `O_3, O_4`. The `t`-channel OPE is `(12)(34)`. In the all-incoming labelling of
  2201.06602 at fixed `u = −(p_1+p_3)²` the Regge pairs are `(1,3)` and `(2,4)`.

## 2. Literature: formulas quoted

### 2.1 Costa, Hansen, Penedones, arXiv:1707.07689 ("CHP")

* **[Q]** (2.7) `σ² = x² x̄²`, `cosh ρ = −x·x̄/(|x||x̄|)`, Regge limit `σ → 0` at fixed `ρ`.
* **[Q]** (2.12)–(2.16) `A(z,z̄) = Σ_J ∫dν b_J(ν²) F_{ν,J}`, and after continuation, to leading order,
  `F_{ν,J} ≈ −iπ^h 4^J σ^{1−J} γ(ν)γ(−ν) Ω_{iν}(ρ)`; (A.16) `γ(ν) = Γ((2Δ_1+J+iν−h)/2) Γ((2Δ_3+J+iν−h)/2)`.
* **[Q]** (2.22)–(2.23) Regge pole: `A(σ,ρ) = ∫dν α(ν) σ^{1−j(ν)} Ω_{iν}(ρ)`,
  `α(ν) = −(i cot(πj(ν)/2) − 1) π^{h+1} 4^{j(ν)} γ(ν)γ(−ν) j′(ν) r(ν)/(4ν)`, (2.21) `r(ν) = C_{12j(ν)} C_{34j(ν)} K_{h±iν,j(ν)}`.
* **[Q]** (2.25)–(2.27), (2.31)–(2.32) impact-parameter transform `A(x,x̄) = (−1)^{−Δ_1−Δ_3} ∫dp dp̄ e^{−2ip·x−2ip̄·x̄} B(p,p̄)`,
  `B = 𝓑(S,L)/((−p²)^{h−Δ_1}(−p̄²)^{h−Δ_3})`, `S = 4|p||p̄|`, `cosh L = −p·p̄/(|p||p̄|)`,
  `𝓑(S,L) = ∫dν β(ν) S^{j(ν)−1} Ω_{iν}(L)`, `α(ν) = γ(ν)γ(−ν) β(ν)`.
* **[Q]** (A.13), (A.17), (A.18):
  `Π_{Δ−h}(ρ) = π^{1−h} Γ(Δ−1)/(2Γ(Δ−h+1)) · e^{(1−Δ)ρ} ₂F₁(h−1, Δ−1; Δ−h+1; e^{−2ρ})`,
  `Ω_{iν} = (iν/2π)(Π_{iν} − Π_{−iν})`, `(∇²_{H_{d−1}} + ν² + (h−1)²) Ω_{iν} = 0`.
* **[Q]** (2.40) spinning Regge pole: `A(x,x̄,z_1,z_2) ≈ ∫dν Σ_k α_k(ν) D_k Ω_{iν}(ρ) / ((−x²+iε_x)^{Δ_1+(j(ν)−1)/2} (−x̄²+iε_x̄)^{Δ_3+(j(ν)−1)/2})`.
  (2.41) currents: `D_1 = z_1·z_2 + (z_1·x̂)(z_2·x̂)`, `D_2 = −(z_1·x̂)(z_2·x̂)`, `D_3 = (z_1·x̂)(z_2·∇) + (z_2·x̂)(z_1·∇)`,
  `D_4 = −(z_1·∇)(z_2·∇) + D_1∇²/(d−1)`; (2.42) ten operators for stress tensors.
* **[Q]** (B.9)–(B.10) Casimir equation in the Regge limit,
  `[½(J^{(x̂)}+J^{(z_1)}+J^{(z_2)})² + ν² + (h−1)²] f(x̂, x̄̂, z_1, z_2) = 0`, solved by `f = 𝒟 Ω_{iν}(ρ)` for any `𝒟` built from
  `(z_i·x̂)`, `(z_i·∇)`, `(z_1·z_2)`, which commute with the total rotation generators.
* **[Q]** (2.47) `𝓑(p,p̄,z_1,z_2) ≈ ∫dν S^{j(ν)−1} Σ_k β_k(ν) D̂_k Ω_{iν}(L)`; (2.48)–(2.49) conservation `p·∂_{z_i} B = 0`
  gives `β_2 = β_3 = 0` for `⟨JJOO⟩` and `β_4 = … = β_{10} = 0` for `⟨TTOO⟩`.
* **[Q]** (2.54) on `H_{d−1}`: `D̂_1{}^μ{}_ν = δ^μ_ν`, `D̂_4{}^μ{}_ν = −∇^μ∇_ν − (ν² + ((d−2)/2)²)/(d−1) · δ^μ_ν`.
* **[Q]** (2.55): `D̂_1 = δ^{\{μ}_{\{ρ} δ^{ν\}}_{σ\}}`, `D̂_2 = δ^{\{μ}_{\{ρ} ∇^{ν\}}∇_{σ\}} + (ν² + ((d−2)/2)²)/(d−1) · δ^{\{μ}_{\{ρ}δ^{ν\}}_{σ\}}`,
  `D̂_3 = ½(∇^{\{μ}∇_{\{ρ}∇^{ν\}}∇_{σ\}} + ∇_{\{ρ}∇^{\{μ}∇_{σ\}}∇^{ν\}}) + (4ν²+d²−2d−2)/(2(d+1)) · δ^{\{μ}_{\{ρ}∇^{ν\}}∇_{σ\}}`
  (curly brackets: symmetrise and remove traces).
* **[Q]** (2.56)–(2.58) `𝒟(ν) = ε^{*ν} ε_μ Σ_{k=1,4} β_k D̂_k` (currents), `ε^{*ρ}ε^{*σ}ε_μ ε_ν Σ_{k=1,2,3} β_k D̂_k` (stress tensors), "polarizations are
  chosen to be complex conjugate `z_1 = z_2^* ≡ ε` because it is only in this configuration that we can relate the
  correlator to the norm of a state"; `β_k/β_j = r_k/r_j = C^{(k)}_{12j(ν)}/C^{(j)}_{12j(ν)}`; (4.7)
  `β_k(ν) = (i cot(πj/2) − 1) π^{h+1} 4^{j−1} (−ij′(ν)/(iν)) r_k(ν)`.
* **[Q]** (4.11)–(4.14), (4.16)–(4.17): with `n` the unit tangent of the geodesic,
  `𝒟Π_{iν} = Π |ε|² β_1 [1 + (β_4/β_1) a(ν,L)(|n·ε|²/|ε|² − 1/(d−1))]`, `a(ν,L)` given in (4.14);
  `𝒟Π_{iν} = Π|ε|⁴β_1[1 + t_2(|n·ε|²/|ε|² − 1/(d−1)) + t_4(|n·ε|⁴/|ε|⁴ − 2/((d+1)(d−1)))]`, `t_2 = −(r_2/r_1)a + (r_3/r_1) f_2`, `t_4 = (r_3/r_1) f_4`,
  `f_2, f_4` in (E.1). Below (4.20): "the three-point function of three stress-tensors … has only two parity
  even tensor structures in `d = 3`".
* **[Q]** (D.9), (D.13), (D.15): `β_k = (ζ(ν,1)/4) Σ M_{kl}(ν, j, Δ_1) α_l` with
  `ζ(ν,n) = 4π^{2−d}/[Γ((2Δ_1+j−h+iν)/2+n) Γ((2Δ_1+j−h−iν)/2+n) Γ((2Δ_3+j−h+iν)/2) Γ((2Δ_3+j−h−iν)/2)]`; for stress tensors the
  matrix is only in the paper's Mathematica notebook (not reproduced here).

### 2.2 Kulaxizi, Parnachev, Zhiboedov, arXiv:1705.02934 ("KPZ")

* **[Q]** (A.1) `Ω_{iν}(ρ) = ν sinh(πν) Γ(h−1+iν)Γ(h−1−iν) / (2^{d−1} π^{(d+1)/2} Γ((d−1)/2)) · ₂F₁(h−1+iν, h−1−iν; h−½; −sinh²(ρ/2))`
  (printed with `sinh(ν)`; the argument must be `πν`, checked in `01`).
* **[Q]** (4.25)–(4.26) `⟨JψψJ⟩`: `δ^{mn}(S,L) = −4c̃_0 Σ_{k=1}^4 ∫dν (π^{d/2} e^{−iπj(ν)/2}/sin(πj(ν)/2)) β_k(ν) S^{j(ν)−1} D̂_k^{mn} Ω_{iν}(L)`,
  `D̂_1 = η^{mn} − p^mp^n/p²`, `D̂_2 = p^mp^n/p²`, `D̂_3 = p^m∂^n + p^n∂^m`, `D̂_4 = p²∂^m∂^n + (p^m∂^n+p^n∂^m) − (η^{mn} − p^mp^n/p²) p²∂²/(d−1)`.
* **[Q]** (5.45) **both pairs spinning**, `⟨J_m J_a J_b J_n⟩`:
  `δ^{mn,ab}(S,L) = 16c̃_0 Σ_{k,ℓ=1}^{4} ∫dν (π^{d/2} e^{−iπj(ν)/2}/sin(πj(ν)/2)) S^{j(ν)−1} β_k(ν) β_ℓ(ν) D̂_k^{mn} D̂_ℓ^{ab} Ω_{iν}(L)`.
  This is the only both-ends formula found; it is written for currents. For stress tensors the same
  structure follows from CHP (B.9)–(B.10) applied to each pair, and from factorisation of the residue. **[I]**

### 2.3 Other sources

* **[Q]** Costa, Gonçalves, Penedones, arXiv:1209.4355, (56)–(58): `A(σ,ρ) ≈ 2πi ∫dν α(ν) σ^{1−j(ν)} Ω_{iν}(ρ)`,
  `α(ν) = −π^{h−1} 2^{j(ν)−1} e^{iπj(ν)/2}/sin(πj(ν)/2) · γ(ν)γ(−ν)β(ν)`; (188) Mellin–Barnes and `₂F₁` forms of `Ω_{iν}`
  (same as KPZ (A.1)). Scalars only.
* **[Q]** Li, Meltzer, Poland, arXiv:1705.03453: (4.5) `Ω_{iν}(ρ) = ν sin(ρν)/(4π² sinh ρ)` in `d = 4`; (C.11) completeness
  `∫_{H_{d−1}} dw′ Ω_{iν}(w,w′) Ω_{iν̄}(w′,w″) = ½(δ(ν−ν̄)+δ(ν+ν̄)) Ω_{iν}(w,w″)`; around (3.36): for spinning correlators the
  polarisations must be chosen "in a reflection positive way" and one must project on Lorentz representations;
  "we have not yet understood a rigorous argument" for part of the spinning statement.
* **[Q]** Afkhami-Jeddi, Hartman, Kundu, Tajdini, arXiv:1610.09378, (3.1)–(3.4), (3.16): `⟨ψ ε_1·T ε_2·T ψ⟩` with polarisations
  related by Rindler reflection, chaos bound `Im δG ≤ 0`; arXiv:1805.07393 §7 (`d = 3`): (7.3)
  `E_{TT}(ρ) ∼ 8π(−iπ ñ_odd AB + 12(4ñ_f − ñ_s)(A²+B²))/(15(1−ρ)⁴)`, (7.5)–(7.7) `ñ_s = 4ñ_f`, `ñ_odd = 0` at infinite gap and
  `|ñ_s − 4ñ_f|/c_T, |ñ_odd|/c_T ≲ ln Δ_gap/Δ_gap⁴`; (7.14) for `⟨TTψ⟩`. In `d = 3` there is **one parity-even and one
  parity-odd** non-Einstein structure: this is the real and imaginary part of the single complex weight-`±4`
  coefficient found in §3.3 below. **[I]** for the identification.
* **[Q]** Chang, Landau, Simmons-Duffin, arXiv:2311.04271: spinning dispersive CFT functionals and their
  bulk-point dictionary in general `d`; "We leave an exploration of the Regge regime for four-point functions
  of spinning operators to future work"; their sum rules reproduce only the flat sum rules that "can be
  expressed in terms of shock amplitudes" (lowest-subtracted), "we have to construct the additional CFT sum
  rules that give rise to the other flat space sum rules".
* **[Q]** Caron-Huot, Li, Parra-Martinez, Simmons-Duffin, arXiv:2201.06602 ("CLPS"):
  (2.1)–(2.3) `M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴ f(s,u)`, `M(1⁺2⁺3⁺4⁻) = ([12][13]⟨14⟩)⁴ g`, `M(1⁺2⁺3⁺4⁺) = [12]²[34]²/(⟨12⟩²⟨34⟩²) h`;
  (2.12) `M = 16π Σ_J (2J+1) a_J^{\{h\}}(s) d^J_{h_{12},h_{34}}(1+2t/s)`; (2.14) `Im a_J^{+−+−} = |c^{+−}_{J,s}|²`, `Im a_J^{++−−} = |c^{++}_{J,s}|²`;
  (2.16) `Im f(m²,−p²) = (16π/m⁸) Σ_{J≥4}(2J+1)|c^{+−}|² d̃^J_{4,4}(1−2p²/m²)` (`s`-cut),
  `Im f(−p², p²−m²) = (16π/m⁸) Σ_{J even}(2J+1)|c^{++}|² d̃^J_{0,0}(1−2p²/m²)` (`t`-cut);
  (2.17)–(2.19) `g ↔ c^{++}c^{+−} d̃^J_{4,0}`, `h ↔ (c^{++})² d̃^J_{0,0}`;
  (2.26)–(2.27) sum rules `B_k^{(1)}` (fixed `u`, `f`), `B_k^{(2)}, B_k^{(3)}, B_k^{(4)}` (fixed `t`: `f(s,p²−s)`, `g`, `h`);
  (2.33) `B_k^{(1)}|_high = ⟨(2m²−p²)(|c^{++}|² d̃^J_{0,0}(x) + |c^{+−}|² d̃^J_{4,4}(x))/(m^{k+4}(m²−p²)^{(k−2)/2})⟩`, `k` even;
  (2.35a) `B_k^{(1)}|_high = ⟨(|c^{++}|² d̃^J_{0,0} − |c^{+−}|² d̃^J_{4,4})/(m^{k+3}(m²−p²)^{(k−3)/2})⟩`, `k` odd;
  (2.35b) `{B^{(2)},B^{(3)},B^{(4)}}|_high = (2m²−p²)⟨{(−1)^J|c^{+−}|² d̃^J_{4,−4}/(m^{k+8}(m²−p²)^{(k+2)/2}), c^{++}c^{+−} d̃^J_{4,0}/(m^{k+8}(m²−p²)^{(k−2)/2}), (c^{++})² d̃^J_{0,0}/(m^k(m²−p²)^{(k+2)/2})}⟩`;
  (A.1)–(A.3) `d^J_{h,h′}(x) = N^J_{h,h′}((1+x)/2)^{(h+h′)/2}((1−x)/2)^{(h−h′)/2} ₂F₁(h−J, J+h+1; h−h′+1; (1−x)/2)`, `d̃` = the same without the two
  power prefactors; (A.11) `lim_{J,m→∞} d̃^J_{h,h′}(1−2p²/m²) = J_{h−h′}(bp)/(p/m)^{h−h′}`, `b = 2J/m`;
  (3.3a) `F_{g3} = ∫_{m_IR}^M dp (1−p)³[p(65+155p+47p²)B_2^{(1)} + 2(81−1927p+4601p²)B_3^{(1)}] − 5∂_{p²}B_4^{(1)}|_{p=0}` and (3.3b)
  `F_{g4}` with `p(99+282p+49p²)`, `(21−474p+1216p²)`, `−(7+∂_{p²})B_4^{(1)}|_0`; "These use only `B^{(1)}`";
  (3.4a) `|ĝ_3|²M⁸ ≤ 37.8 log(M/m_IR) − 45.4 − F_matter`; (3.5) `lim_{m→∞} F[m², bm/2] = (4/m⁴)∫_0^1 dp ψ(p)(|c^{++}|²+|c^{+−}|²)J_0(pb) + …`,
  "In this limit, `F^{++}` and `F^{+−}` coincide"; (3.6) `≈ 160/(m⁴b³) − 260 m_IR² J(m_IR b)/(m⁴ m_IR b)`;
  (3.16) minimal exponents of the improved spin-2 wavefunctions `(n^{(1)}, n^{(2)}, n^{(3)}, n^{(4)})_min = (1, 10, 6, 1)`;
  (3.19)–(3.25) the eikonal matrix of sum rules, diagonal `½B_2^{(1)}` with `J_0`, off-diagonal `p̂⁴B_2^{(3)}` with `J_4`, heavy side
  a positive matrix times `δ(2J/m − b)`; (4.4) optimal `|ĝ_3|²M⁸ ≤ 24.9 log(M/m_IR) − 27.6`, obtained with "all improved sum
  rules `B_2^{imp}` and `B_3^{imp}`", heavy densities being `n×n` matrices (their App. C).
  *Reading of (3.3a).* The extracted text does not show the bracketing. **[D]** (`06`): only the reading in
  which `(1−p)³` multiplies both polynomials and `p` only the `B_2` one is non-negative on all sampled `(m, J)`
  (`m ∈ [1,16]`, `J ≤ 100`, both channels; minimum `m^{10}F = 5.06` at `m = 1, J = 4` in `+−` and `7.39` at `m = 1, J = 0` in `++`); the
  other two readings are negative somewhere. The factor `4/m⁴` in (3.5) is fixed by consistency with (3.6).
* **[Q]** Caron-Huot, Mazáč, Rastelli, Simmons-Duffin, arXiv:2106.10274 ("CHMRS"): see §6.1.

## 3. One spinning pair in `d = 3` (Task 2)

### 3.1 The harmonic function

**[D]** (`01`, 30 digits; CHP (A.13)+(A.17) and KPZ (A.1) agree):

`Ω_{iν}(L) = (ν tanh πν / 4π) · P_{−1/2+iν}(cosh L)` in `d = 3`.

### 3.2 Reading of the operators, cross-checked against the source

**[D]** (`05`). With `D̂_4 = −∇∇ − λ_d δ/(d−1)` acting on any radial solution, `ε^{*}·D̂_4·ε = −(F″ − coth L F′)(|n·ε|² − |ε|²/(d−1))`. Hence
CHP's `a(ν,L)` must equal `−(Π″ − coth L Π′)/Π`: verified against the printed (4.14) to 12 digits in `d = 3, 4, 6`.
Likewise the coefficient of `|n·ε|⁴` produced by `D̂_3` is `(sinh L ∂_{cosh L})⁴ F`, so CHP's `f_4` must equal
`(sinh L ∂_{cosh L})⁴Π/Π`: verified against the printed (E.1) to 12 digits in `d = 3, 4, 6` (this also fixed the
transcription of (E.1) from the PDF text). So the reading of (2.54)–(2.55) used below is the authors'.

### 3.3 Results

Geodesic polar coordinates centred at the impact point of the scalar pair; `ε_±` at the impact point of
the spinning pair; `z_1` contracts the indices of operator 1, `z_2` those of operator 2. All identities below
hold modulo the radial equation and were obtained by explicit covariant differentiation (`02`, exact).

**Currents.** **[D]**

| `(z_1, z_2)` | weight | `D̂_1 F` | `D̂_4 F` |
|---|---|---|---|
| `(ε_+, ε_−)` (CHP's `z_2 = z_1^*`, helicity preserved) | 0 | `F` | `0` |
| `(ε_±, ε_±)` (helicity flipped) | `±2` | `0` | `−½ P²` , `P² = F″ − coth L F′` |

For general tangent `z_i = z_i^+ ε_+ + z_i^− ε_−`:

`𝓑^{JJOO} ≈ ∫dν S^{j(ν)−1} (ν tanh πν/4π) [ β_1(ν) (z_1^+z_2^− + z_1^−z_2^+) P_{−1/2+iν}(cosh L) − ½ β_4(ν) (z_1^+z_2^+ + z_1^−z_2^−) P²_{−1/2+iν}(cosh L) ]`,

`P²_{−1/2+iν} = (ν²+¼)(ν²+9/4) P^{−2}_{−1/2+iν}`. Equivalently `D̂_4F = −P²[(z_1·n)(z_2·n) − ½ z_1·z_2]`, which is CHP (4.13) at `d = 3`.

**Stress tensors.** **[D]**

| `(z_1, z_2)` | weight | `D̂_1 F` | `D̂_2 F` | `D̂_3 F` |
|---|---|---|---|---|
| `(ε_+, ε_−)` | 0 | `F` | `0` | `0` |
| `(ε_±, ε_±)` | `±4` | `0` | `0` | `¼ P⁴` , `P⁴ = (sinh L ∂_{cosh L})⁴ F` |

So in `d = 3`, for polarisation tensors `e_i = a_i ε_+ε_+ + b_i ε_−ε_−`,

`𝓑^{TTOO} ≈ ∫dν S^{j(ν)−1} (ν tanh πν/4π) [ β_1(ν) (a_1b_2 + b_1a_2) P_{−1/2+iν}(cosh L) + ¼ β_3(ν) (a_1a_2 + b_1b_2) P⁴_{−1/2+iν}(cosh L) ]`,

`P⁴_{−1/2+iν} = (ν²+¼)(ν²+9/4)(ν²+25/4)(ν²+49/4) P^{−4}_{−1/2+iν}`.

Remarks. (i) `D̂_2 ≡ 0` and the weight-0 part of `D̂_3` vanishes: in `d = 3` there are exactly two parity-even
structures (weight 0 and weight `±4`), in agreement with CHP's remark and with the single non-Einstein
parity-even coefficient of 1805.07393 (7.3). (ii) **In the helicity-preserving sector the derivative operators
produce no `ν`-dependent factor at all**: the kernel is the scalar harmonic function times the coupling `β_1(ν)`.
(iii) In the flip sector the operators produce `P^{+m}`, i.e. `Π_{k<m}(ν²+(k+½)²)` times the function `P^{−m}` that is
`ν`-independent at small `L`; relative to the unitary coefficient, `P^m = (−1)^m(½−iν)_m Φ^ν_{m,0}`, a factor of modulus
`[Π_{k<m}(ν²+(k+½)²)]^{1/2}`. (iv) General statement **[D]** (`02`, `m ≤ 4`, exact): `(ε_±·∇)^m F = 2^{−m/2} P^m_{−1/2+iν}(cosh L)`
(derivative taken at the point at distance `L` from the centre, `n` pointing away from the centre).
(v) Parity-odd structures of `d = 3` are not in CHP's basis; in the helicity basis they are the imaginary part
of the weight-`±m` coefficient (different coefficients for `+m` and `−m`) and an antisymmetric weight-0
structure `∝ i(z_1^+z_2^− − z_1^−z_2^+)`, with the same kernels. **[I]**

## 4. Both pairs spinning (Task 3)

### 4.1 Theorem R (raising operators on `H²`)

Disk model, `ds² = e^{2σ}dw dw̄`, `e^{2σ} = 4/(1−|w|²)²`; `ε_− = √2 e^{−σ}∂_w`, `ε_+ = √2 e^{−σ}∂_{w̄}`; on the `(w…w)` component of a
rank-`k` tensor `∇_w = ∂_w − 2k ∂_wσ`, `2∂_wσ = 2w̄/(1−|w|²)`. Plane waves `E_s(w;b) := ((1−|w|²)/|1−w b̄|²)^s`, `|b| = 1`.

**Lemma.** `∇_w^m E_s = (s)_m A^m E_s`, `A := (b̄ − w̄)/((1−|w|²)(1−w b̄))`.
*Proof.* `∂_w E_s = sAE_s` and `∂_wA = A² + 2(∂_wσ)A`, because `A + w̄/(1−|w|²) = b̄/(1−w b̄)`. Then
`(∂_w − 2k∂_wσ)(A^kE_s) = k A^{k−1}(A² + 2∂_wσ A)E_s + sA^{k+1}E_s − 2k∂_wσ A^kE_s = (s+k)A^{k+1}E_s`. ∎ (Also checked by sympy for `m ≤ 5`, `03`.)

**Theorem R.** Let `F(x,y) = P_{−1/2+iν}(cosh L(x,y))`, frames at `x` and `y` parallel along the geodesic, `β = L(x,y)`. Then

1. `(ε_−·∇_x)^m (ε_+·∇_y)^m F = 2^{−m} Π_{k=0}^{m−1}(ν²+(k+½)²) · Φ^ν_{m,m}(β)`;
2. `(ε_−·∇_x)^m (ε_−·∇_y)^m F = 2^{−m} Π_{k=0}^{m−1}(ν²+(k+½)²) · Φ^ν_{m,−m}(β)`;
3. more generally `(ε_−·∇_x)^m (ε_+·∇_y)^n F = 2^{−(m+n)/2} (s)_m (1−s)_n Φ^{−ν}_{m,n}(β)`; in particular
   `(ε_−·∇_x)^m F = 2^{−m/2}(½+iν)_m Φ^{−ν}_{m,0}(β) = (−1)^m 2^{−m/2} P^m_{−1/2+iν}(cosh β)`.

*Proof.* `F(x,y) = (1/2π)∫_0^{2π} E_s(x;e^{iθ}) E_{1−s}(y;e^{iθ}) dθ` (both sides are invariant under the isometry group,
since `E_s(gx;gb)|g′(b)|^s = E_s(x;b)`; at `x = 0` it is the integral defining `Φ_{0,0}`). Differentiate under the
integral with the Lemma (and its conjugate at `y`), then put `x = 0`, `y = tanh(β/2)`: `√2e^{−σ}A|_{x=0} = 2^{−1/2}e^{−iθ}`,
`√2e^{−σ}Ā|_y = 2^{−1/2}e^{iθ′}` with `e^{iθ′}` the boundary Möbius map of the definition of `Φ`, and
`E_{1−s}(y) = (cosh β − sinh β cos θ)^{−1/2+iν}`. This gives `2^{−(m+n)/2}(s)_m(1−s)_nΦ^{−ν}_{m,n}`; `(s)_m(1−s)_m = Π_{k<m}(ν²+(k+½)²)`, and
`Φ_{m,±m}` are even in `ν`. ∎

**[D]** Numerical verification (`03`): symbolic two-point derivatives in the disk, evaluated with exact
derivatives of `P_{−1/2+iν}`, against the integral definition of `Φ^ν_{m,n}`; `m = 1, 2, 3, 4`, three `(ν, β)` each,
all three statements; maximal relative deviation `6·10^{−24}`. The repository's earlier check
`−(F″+F′/sinh β)/(ν²+¼) = Φ_{1,1}` is the case `m = 1` of statement 1.

**Answer to "which factor".** For `(m, ±m)` it is the **polynomial** `Π_{k=0}^{m−1}(ν²+(k+½)²) = |Γ(½+m+iν)/Γ(½+iν)|²`
(times `2^{−m}` from unit-normalised helicity vectors): `(ν²+¼)(ν²+9/4)` for currents, `Π_{k=0}^{3}(ν²+(k+½)²) = ν⁸ + 21ν⁶ + (987/8)ν⁴ + (3229/16)ν² + 11025/256`
for stress tensors. For `(m,0)` relative to the unitary coefficient it is the Pochhammer symbol `(½+iν)_m` (complex);
relative to `P^{−m}` it is again the polynomial. Its zeros `ν = ±i(k+½)`, `k < m`, are the points of the discrete
series containing the `K`-type `m`: the derivative-generated kernels span only the principal-series part
of the weight-`m` harmonic analysis on `H²` (cf. `SPINNING_KERNELS.md`, Theorem K(5)).

### 4.2 The Regge amplitude in the helicity basis

Inserting §3.3 and Theorem R in KPZ (5.45) (currents **[Q]+[D]**; stress tensors with the same structure **[I]**),
with the residue written as `β_k^Lβ_ℓ^R` and all common factors (`S^{j(ν)−1}`, signature factor, `ν tanh πν/4π`) dropped:

| Regge-pair weights | currents (`m = 2`) | stress tensors (`m = 4`) |
|---|---|---|
| `(0,0)` | `β_1^Lβ_1^R · P_{−1/2+iν}` | `β_1^Lβ_1^R · P_{−1/2+iν}` |
| `(m,0)` | `−½ β_4^Lβ_1^R · P² = −½β_4β_1 (ν²+¼)(ν²+9/4) P^{−2}` | `¼ β_3^Lβ_1^R · P⁴ = ¼β_3β_1 Π_{k<4}(ν²+(k+½)²) P^{−4}` |
| `(m,m)` (flat limit `J_0`) | `¼ β_4^Lβ_4^R (ν²+¼)(ν²+9/4) Φ_{2,2}` | `(1/16) β_3^Lβ_3^R Π_{k<4}(ν²+(k+½)²) Φ_{4,4}` |
| `(m,−m)` (flat limit `J_{2m}`) | `¼ β_4^Lβ_4^R (ν²+¼)(ν²+9/4) Φ_{2,−2}` | `(1/16) β_3^Lβ_3^R Π_{k<4}(ν²+(k+½)²) Φ_{4,−4}` |

So the *`t`-channel* (Regge-exchange) description of the weight-`(m,±m)` component is "coupling × polynomial ×
unitary matrix coefficient". The polynomial belongs with the three-point coupling of the exchanged
trajectory to the flipped pair (e.g. it vanishes at `ν² = −9/4`, the stress-tensor point, for `m ≥ 2`: the
*harmonic* function has no component there, while the *propagator* `Π_{iν}` with the same operators does not
vanish — graviton exchange with the non-minimal coupling does flip helicity, CLPS (3.25)).

## 5. Which helicity process has which kernel (dictionary; confirms the corrected Part B)

**Claim.** The SO(2) weights are the helicity transfers of the Regge pairs; in the all-incoming labelling of
CLPS at fixed `u`: `m_L = h_1 − h_3^{out}`, `m_R = h_2 − h_4^{out}` (i.e. `h_1 + h_3`, `h_2 + h_4` in all-incoming helicities, up to the
orientation sign of the second particle). They are **not** the labels `λ = h_1 − h_2`, `μ = h_3 − h_4` of `d^J_{λμ}`.

Evidence:

1. **[Q]+[D]** CHP (2.56): the configuration computing a norm (same particle in and out, same helicity) is
   `z_2 = z_1^*`; by §3.3 this is weight 0 and its kernel is `Ω_{iν}` with no derivative. Flipping the helicity of
   the outgoing operator gives weight `±m` and `P^{±m}`.
2. **[Q]+[D]** Flat benchmark. CLPS (A.11): `d̃^J_{h,h′} → J_{h−h′}(bp)/(p/m)^{h−h′}`; AdS: `Φ_{m,n} → J_{|m−n|}(νβ)` (`03`). The Bessel order
   is the total change of angular momentum about the collision axis, `|ΔJ_z^{(a)} + ΔJ_z^{(b)}|`, i.e. it is fixed by the
   helicity transfers. For `B^{(1)}` both transfers vanish (`J_0`), although the `s`-channel labels are `(4,4)` for `|c^{+−}|²`
   and `(0,0)` for `|c^{++}|²`. `d̃^J_{4,4}` and `d̃^J_{0,0}` differ only at relative order `1/J² ∼ 1/(bm)²` at fixed `b`
   (`04`(a): `m²(d̃_{4,4} − d̃_{0,0}) → 4.55` at `b = 3, p = 0.7`), an effect that vanishes in the regime `m → ∞`, `β` fixed.
3. **[Q]** CLPS (3.20)–(3.25): in the eikonal matrix (helicity of particle 1 → 4 against a target) the diagonal
   entries are `½B_2^{(1)}` with `J_0` and heavy weight `|c^{++}|² + |c^{+−}|²`, the off-diagonal ones `p̂⁴B_2^{(3)}` with `J_4`.
4. **[Q]** CHMRS (3.17)–(3.19): a heavy `s`-channel block is, in the Regge limit, a plane wave `e^{x·p+y·p̄}`, `p = (Δ, J e_z)`,
   i.e. a point in CHP's impact-parameter space with `cosh L = η_AdS`; the functional's `ν` is the harmonic label
   on the transverse `H_{d−1}`, on which the stabiliser of the *Regge* limit acts. **[I]** For spinning operators the
   same plane wave carries the polarisation tensors of the two pairs, and its component in the weight-`(m_L,m_R)`
   sector is the tensor harmonic of §4 evaluated at `β`.

The label confusion of the first version ("continue `J → −½+iν`, `θ → iβ` in `d^J_{4,4}`") exchanges label and
argument: physically `J ↔ mβ/2` and `θ ↔ p ↔ ν`.

**Dictionary** (stress tensors; flat sum rules of CLPS; kernels at leading order `m → ∞`, `β` fixed; normalisations
from §6.2):

| flat sum rule | process (physical in → out) | pair weights | heavy coupling, flat partial wave | type | AdS kernel |
|---|---|---|---|---|---|
| `B_k^{(1)}` (fixed `u`, `f`), `t`-cut | `++ → ++` | `(0,0)` | `|c^{++}|² d̃^J_{0,0}` | diagonal, `≥ 0` | `P_{−1/2+iν}(cosh β)` **[D]** given the structure |
| `B_k^{(1)}`, `s`-cut | `+− → +−` (MHV, helicities preserved) | `(0,0)` | `|c^{+−}|² d̃^J_{4,4}` (sign `(−1)^k` relative to `++`) | diagonal, `≥ 0` | `P_{−1/2+iν}(cosh β)` **[D]** given the structure |
| `B_k^{(2)}` (fixed `t`, `f`) | `+− → −+` (both flip; angular-momentum transfers add, `|ΔJ_z| = 8`, flat `J_8`) | `(4,−4)` | `(−1)^J |c^{+−}|² d̃^J_{4,−4}` | diagonal coupling, sign alternates with `J` | `m⁸ Φ_{4,−4}(β)/Π_{k<4}(ν²+(k+½)²)` **[I]** |
| `B_k^{(3)}` (`g`) | `++ → +−` (one flip) | `(4,0)` | `c^{++}c^{+−} d̃^J_{4,0}` | interference | `m⁴ P^{−4}_{−1/2+iν}(cosh β)` (up to a phase convention) **[I]** |
| `B_k^{(4)}` (`h`) | `++ → −−` (both flip `+ → −`; angular-momentum transfers cancel, `ΔJ_z = 0`, flat `J_0`) | `(4,4)` | `(c^{++})² d̃^J_{0,0}` | interference (phase free) | `Φ_{4,4}(β)`, unit normalised **[I]** |

Two clarifications of the corrected Part B table. (a) The `(−1)^J Φ_{4,−4}` term is not an addition to the `+− → +−`
kernel: it is a different polarisation component (a different sum rule, `B^{(2)}`), which contains the same
coupling `|c^{+−}|²`. (b) For currents the same table holds with `4 → 2` and `P^{−4} → P^{−2}`.

## 6. The measure and the normalisation (Task 4)

### 6.1 Scalar case, quoted

* **[Q]** CHMRS (2.49)–(2.50), (2.52)–(2.58): one seeks functionals with `ω[Δ,J] ∼ P_{(2−d)/2+iν}(η_AdS)`; the harmonic transform
  of the Regge moment `Π̂_{k,ν} = ∫_1^∞[dη] P_{(2−d)/2+iν}(η) Π_{k,η}` acts on heavy blocks as (2.54)
  `Π̂_{k,ν}[Δ,J] ∝ γ_{2Δ_φ+k−1}(ν)² P_{(2−d)/2+iν}(η_AdS)`, `γ_a(ν) = Γ((1+a−d/2−iν)/2)Γ((1+a−d/2+iν)/2)`, and `C_{k,ν}` is **defined** (2.56) by dividing
  by `γ²` ("the CFT angle-of-approach `η` is a smeared version of the bulk impact parameter `η_AdS`"), so that
* **[Q]** (2.59) `C_{k,ν}[Δ,J] = (2P_{(2−d)/2+iν}(η_AdS)/m^{2k})(1 + O(ν²/m², 1/m²))` in the Regge regime, and (2.60)
  `C_{k,ν}[Δ,J] = (2m²−ν²)P_J(1−2ν²/m²)/((m²−ν²)(m²(m²−ν²))^{k/2}) · (1 + O(J²/m²))` in the bulk-point regime, "valid for all `ν ∈ [0, m)`";
  "the `β → 0` limit of the Regge limit overlaps with the bulk point limit" (below (5.19)).
* **[Q]** (3.17)–(3.20): the origin of `γ²`: the heavy block is `∝ ∫d^dx d^dy (|x||y|)^{k−d−1+2Δ_φ} P_{(2−d)/2+iν}(−x·y/|x||y|) e^{x·p+y·p̄}`, and the
  Fourier–Laplace transform of a Gegenbauer function is a Gegenbauer function times `γ_a(ν)` (3.18). This is the
  same transform as CHP (D.6)–(D.9) (`α = γγβ`).
* **[Q]** (5.26)–(5.27) `ω_flat = ∫_0^M (dp/M) f(p/M) C^{imp}_{2,−p²}`, `ω_AdS = ∫_0^M (dν/M) f(ν/M) C^{imp}_{2,ν}`; (5.28) "`β ≫ 1/M`, `m` fixed: improvement
  dominates; `β` fixed, `m ≫ M`: unimproved dominates"; (5.29) `f ∼ (1−p)^l`, `l > n+1−d/2`;
  (5.36) `lim m⁴ω_flat = 2∫_0^1 dp f(p) J̃(Mβp)`; (5.41)–(5.43) `lim_{m→∞} m⁴ω_AdS → 2∫_0^∞ (dν/M) f_0 (ν/M)^n P_{(2−d)/2+iν}(cosh β) + osc. = (2f_0/M^{n+1}) Q_n(β)`,
  `Q_n(β) := ∫_0^∞ dν ν^n P_{(2−d)/2+iν}(cosh β)`; (5.44) asymptotics; "`Q_n(β)` is everywhere positive … for `0 ≤ n ≤ 1` in `d = 3`".

So in the scalar case the reduction is: (i) `C_{2,ν}` normalised so that the heavy action is the unit-normalised
harmonic function; (ii) plain measure `dν` (the image of `dp`), **not** the Plancherel measure `ν tanh πν dν`;
(iii) endpoint analysis in `ν`: the region `ν ∼ 1` selects the Taylor coefficients of `f` at `p = 0`, giving
`Σ_n a_n Q_n(β)/M^{n+1}` for `f = Σ a_n p^n`; the endpoint `ν = M` is oscillatory and subleading under (5.29).

### 6.2 Spinning case: what replaces `ν^n P`

There is no derivation in the literature (2311.04271: "future work"). What can be said:

1. **[D] given the structure of §2 / [I] overall.** In the sector of pair weights `(m_L, m_R)` the `β`-dependence at fixed `ν`
   is `Φ^ν_{m_L,m_R}(β)` up to a factor `c(ν)` (and up to a positive `ν`-independent function of `β`, which can be
   absorbed in the heavy density). This is what §3–§4 derive from the Regge-theory operators.
2. **The factor `c(ν)` is a convention for the family of functionals, not physics** — exactly as the division by
   `γ²` in CHMRS (2.56). The physical statement is about the functional as a whole. Unitarity does *not* fix
   `c(ν)`: `g ↦ c(ν)Φ^ν_{m,m}(g)` is a positive-definite function on the group for every `c(ν) > 0`. Positivity of
   the heavy density fixes the *shape* (a diagonal matrix coefficient of a unitary representation) and
   the sign, nothing more.
3. **[I]** What fixes `c(ν)` for "the uplift of a given flat functional `∫dp ψ(p)B(p²)`" is the requirement, built
   into CHMRS's definition, that the CFT functional labelled by `ν` reproduce the flat sum rule at `u = −ν²`
   in the bulk-point regime **for all `ν < m`**, and that the two regimes overlap at `β → 0`. At `p ≪ m/J` the flat
   stripped partial waves are `p`-independent, `d̃^J_{h,h′}(1) = N^J_{h,h′}`. Hence the AdS kernel must have a
   `ν`-independent limit as `β → 0`, equal to the large-`J` form of `N^J_{h,h′}`. **[D]** (`07`) this singles out:
   * `(0,0)`: `P_{−1/2+iν}(cosh β) → 1`;
   * `(4,4)`: `Φ^ν_{4,4}(β) → 1` (unit normalisation; **not** `Π_{k<4}(ν²+(k+½)²)Φ_{4,4}`);
   * `(4,0)`: `m⁴P^{−4}_{−1/2+iν}(cosh β) → (mβ/2)⁴/4! ≈ J⁴/4! ≈ N^J_{4,0}`;
   * `(4,−4)`: `m⁸Φ_{4,−4}(β)/Π_{k<4}(ν²+(k+½)²) → (mβ/2)⁸/8! ≈ N^J_{4,−4}` (with `ν⁸` in place of the polynomial the limit
     depends on `ν`: ratio `390` at `ν = 1`, `5.1` at `ν = 3`).
   In the flip channels the flat "stripping" of `(p/m)^{h−h′}` is thus replaced in AdS by division by the
   polynomial of Theorem R — the polynomial is *not* absorbed in the heavy coupling there; in the `(4,4)`
   channel it is absent.
4. **Consequently**, for the diagonal `(0,0)` and `(4,4)` channels the combination replacing `ν^nP` is `ν^nΦ^ν_{m,m}(β)`
   with unit normalisation and measure `dν`; in the flat limit all reduce to `ν^nJ_0(νβ)`. Level of rigour: the
   kernel shape is derived from the Regge-theory operators; the normalisation rests on the overlap
   argument, which presupposes spinning analogues of CHMRS (2.59)–(2.60) that nobody has derived.
5. A residual freedom that is genuinely a choice of functional: `f(ν/M) → f(ν/M) + M^{−2}f_1(ν/M) + …`. It changes the
   curved-regime action at the order of the `a_2Q_2/M³` term (through `Q_0/M³`) but not the `a_1Q_1/M²` term.

## 7. Consequences (Task 5)

### 7.1 Moments of the kernels

**[D]** (`04`). Regulated moments `∫_0^∞ e^{−κν}ν^nΦ^ν_{m,m}dν` at `κ = 0.4` by direct integration of the `₂F₁` form agree with
the Chebyshev–Mehler `t`-integral to 8–10 digits (`m = 2, 4`; `n = 0, 1, 2`), and the `κ → 0` limits agree with the
closed forms of `SPINNING_KERNELS.md` Theorem K to 7 digits (`m = 0, 2, 4`; `β = 0.6, 1, 2`); `Q_1^{(m)}(0⁺) = 1/24 − m²/2`
(`0.04166, −1.9579, −7.9516` at `β = 0.02`). The mathematics of Theorem K(1)–(4) stands.

| `β` | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.8 | 1.0 | 1.5 | 2.0 | 3.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| `Q_1^{(0)}` | 0.04166 | 0.04162 | 0.04149 | 0.04127 | 0.04058 | 0.03897 | 0.03756 | 0.03319 | 0.02820 | 0.01858 |
| `Q_2^{(0)}` | −7997.5 | −998.8 | −124.4 | −36.62 | −7.754 | −1.804 | −0.8834 | −0.2248 | −0.0770 | −0.0129 |
| `Q_1^{(4)}` | −7.916 | −7.791 | −7.307 | −6.551 | −4.519 | −1.406 | 0.0015 | 0.7590 | 0.0483 | −0.1657 |
| `Q_2^{(4)}` | −8155.9 | −1075.6 | −158.2 | −54.47 | −11.30 | 1.663 | 3.051 | 0.8773 | −0.4958 | −0.0679 |
| `|Q_1^{(4)}|/Q_1^{(0)}` | 190.0 | 187.2 | 176.1 | 158.7 | 111.4 | 36.1 | 0.04 | 22.9 | 1.71 | 8.92 |

### 7.2 The weights of the CLPS functionals in the `++` and `+−` channels

**[Q]+[D]** From (3.3a), (2.33), (2.35a), in units `M = 1`:

* `B_2^{(1)}`: `ψ_2(p) = (1−p)³p(65+155p+47p²) = 65p − 40p² − 223p³ + 259p⁴ − 14p⁵ − 47p⁶`, **the same for `|c^{++}|²` and `|c^{+−}|²`**, heavy weight
  `(2m²−p²)/m⁶`. It starts at order `p` with `a_1 = +65`, `a_2 = −40`. (For `F_{g4}`: `99p − 15p² − 500p³ + …`.)
* `B_3^{(1)}`: `ψ_3(p) = 2(1−p)³(81−1927p+4601p²) = 162 − 4340p + 21250p² − …`, entering with sign `+` for `|c^{++}|²` and `−` for `|c^{+−}|²`,
  heavy weight `1/m⁶` (suppressed by `M²/m²` relative to `B_2`); it starts at order `p⁰`.
* `−5∂_{p²}B_4^{(1)}|_{p=0}`: the same sign in both channels, `5(4d̃′(1) − 1)/m^{10}` with `d̃′_{0,0}(1) = J(J+1)/2`, `d̃′_{4,4}(1) = (J−4)(J+5)/2`.
* In the scaling limit `m → ∞`, `b` fixed, only `B_2` survives and the two channels coincide (their (3.5));
  `160/(m⁴b³)` in (3.6) is `(4/m⁴)·|a_2|/b³`, and the `a_1` term contributes only through the IR cutoff,
  `−260 m_IR J_1(m_IR b)/(m⁴b)`, with the *wrong* sign.
* The optimal functionals (their (4.4)) use all four improved spin-2 sum rules with `ψ_i = Σ_n c_{in}(1−p)²p^n`,
  `n ≥ (1, 10, 6, 1)`: order `p` for `B^{(1)}` and `B^{(4)}`, `p^{10}` for `B^{(2)}`, `p⁶` for `B^{(3)}`; individual coefficients are not printed.

### 7.3 Curved-regime action

1. **Diagonal MHV-type components (`B^{(1)}`).** By §5 both `|c^{++}|²` and `|c^{+−}|²` have pair weights `(0,0)`; by §3.3 there is
   no `ν`-dependent operator factor in that sector; by §6.2 the normalisation is the scalar one. Hence, for
   `m → ∞` at fixed `β ≫ 1/M` (CHMRS's "unimproved" regime),

   `F^{++}, F^{+−} → (c/m⁴)[a_1 Q_1(β)/M² + a_2 Q_2(β)/M³ + O(M^{−4})]`, `(a_1,a_2) = (65, −40)`, `Q_n` the **scalar** moments **[I]**

   (`c = 2` from the density (2.33), `c = 4` in the normalisation of CLPS (3.5); `m` in units of the gap).

   The positivity input is Lemma Q (`Q_1 > 0`, `Q_2 < 0`), in both channels. Where flat space has the IR-cutoff
   term with the wrong sign, AdS4 has `+a_1Q_1/M²`, `Q_1(0⁺) = 1/24`.
2. **At finite `m/M` the other two terms of `F_{g3}` compete**, as in the scalar case (CHMRS (5.28)), and one of them
   distinguishes the channels: in units of the gap for `m`,
   `F^{±} ≈ (2/(m⁴M))[65Q_1/M + 40|Q_2|/M²] ± 162 Q_0(β)/(m⁶M) + 10M²(−∂_{ν²}P_{−1/2+iν})_{ν=0}/m⁸ + …` (`+` for `++`, `−` for `+−`). **[I]**
   The `B_3` term (negative for `+−` when `∫ψ_3P > 0`) is dominated by the first term for `m² ≫ 1.2 M Q_0/Q_1` and by the
   third for `m² ≪ M³(−∂_{ν²}P)_0/(16Q_0)`; these ranges overlap for large `M` (crude estimate: `Q_0` stands for
   `M∫ψ_3P/162`, which is accurate only for `Mβ ≫ 1`). So the additional positivity inputs are `Q_0 > 0`,
   `(−∂_{ν²})P|_{ν=0} > 0` and dominance estimates — the same list as in the scalar case (Lemma M), now with a
   channel-dependent sign. A direct evaluation of this model (`08`, `09`; conditional, corrections of relative
   order `1/m²` unknown):

   | `M` | `β` | `m = 3` | `10` | `30` | `100` | `1000` |
   |---|---|---|---|---|---|---|
   | 20 | 0.3 | 275 / 2192 | 705 / 877 | 777 / 796 | 786 / 788 | 787 / 787 |
   | 20 | 1.0 | 4836 / 4380 | 67.8 / 26.8 | 12.8 / 8.20 | 10.2 / 9.81 | 10.0 / 10.0 |
   | 20 | 2.0 | 15314 / 15053 | 139 / 115 | 6.80 / 4.19 | 4.11 / 3.87 | 3.98 / 3.98 |
   | 100 | 0.3 | 2.81e5 / 2.71e5 | 2729 / 1810 | 113 / 10.7 | 38.9 / 29.7 | 34.1 / 34.1 |
   | 100 | 1.0 | 2.88e6 / 2.87e6 | 23447 / 23144 | 310 / 276 | 9.45 / 6.41 | 5.62 / 5.59 |

   | 100 | 2.0 | 9.49e6 / 9.49e6 | 76929 / 76794 | 960 / 945 | 12.1 / 10.7 | 3.73 / 3.72 |

   (entries: `m⁴M²F^{++}` / `m⁴M²F^{+−}`). All sampled values are positive; the `+−` channel has a dip at intermediate
   `m` (e.g. `10.7` against `34` at `M = 100`, `β = 0.3`, `m = 30`). A scan over `m ∈ [1.5, 6600]` (110 values) for
   `M ∈ {30, 100, 300}`, `β ∈ {0.1, 0.2, 0.3, 0.5, 1}` (`09`) finds the minimum of `m⁴M²F^{+−}` positive in all 15 cases
   (smallest `5.09`, at `M = 300`, `β = 1`; `m → ∞` value there `5.12 ≈ 2·65·Q_1(1)`). This is an illustration of the
   conditional model, **not** a positivity proof.
3. **Does Part C stand?** The *original* Part C (first moment `1/24 − m²/2 < 0` in the `+−` channel, so positivity
   of a weight `a_1p + a_2p²`, `a_1 > 0 > a_2`, does not uplift) does **not** stand: its kernel was assigned to the wrong
   process. The value `Q_1^{(4)}(0⁺) = −191/24` is correct for the unit-normalised `Φ_{4,4}` — and by §6.2 unit
   normalisation is the right one for the `(4,4)` channel if the overlap argument holds — but that kernel
   multiplies the interference term `(c^{++})²`. The retraction and items 1–2 of the rewritten Part C agree
   with the derivation here. One amendment to its remark (ii): in the scalar/`(0,0)` sector the measure is
   fixed (`dν`, unit-normalised kernel, CHMRS (5.27), (5.41)); the Plancherel alternative is not open.
4. **Matrix positivity when helicity-changing components are used.** For even `J` the heavy data are
   `(c^{++}, c̄^{++}, c^{+−})` and a functional `Σ_i∫ψ_iB^{(i)} + c.c.` acts as the quadratic form
   `A_{00}^{++}|c^{++}|² + [A_{00}^{+−} + (−1)^JA_{4,−4}]|c^{+−}|² + Re(A_{44}(c^{++})²) + Re(A_{40}c^{++}c^{+−})`,
   which must be positive semi-definite for every `(m, J)`: necessary are `A_{00}^{++} ≥ |A_{44}|` and `A_{00}^{+−} ≥ |A_{4,−4}|`, and the full
   condition is the `3×3` determinant condition including `A_{40}`; for odd `J` only `c^{+−}` exists. **[Q]** for flat space
   (CLPS (2.13), App. C); **[I]** for AdS with the kernels of §5. In the curved regime, with `ψ^{(4)} = a_1^{(4)}p + a_2^{(4)}p² + …`,

   `a_1^{(1)}Q_1^{(0)}(β) + a_2^{(1)}Q_2^{(0)}(β)/M ≥ |a_1^{(4)}Q_1^{(4)}(β) + a_2^{(4)}Q_2^{(4)}(β)/M|` (same normalisation of `A_{44}` as of `A_{00}`),

   and in the window `M^{−1/3} ≪ β ≪ 1`, where the first-moment terms dominate, this requires
   `|a_1^{(4)}| ≤ a_1^{(1)} Q_1^{(0)}/|Q_1^{(4)}|`, i.e. `≈ a_1^{(1)}/191` for asymptotically large `M` (the ratio is `176` at `β = 0.2`, `111` at
   `β = 0.5`; the window needs `β ≫ 2.5 M^{−1/3}` for `(a_1,a_2) = (65,−40)`) — a constraint with no flat-space counterpart (in flat space
   `∫dp p J_0(pb) = 0` for `b > 0`, so the order-`p` coefficients are unconstrained at large `b`). `A_{44}` oscillating in
   sign is harmless. Which functionals are affected: `F_{g3}`, `F_{g4}` of CLPS (3.3) use **only** the diagonal
   components (`B^{(1)}`); the optimal bound (4.4) uses all four, so its uplift needs the matrix condition and
   the off-diagonal kernels with the normalisations of §6.2.

## 8. What remains underived

1. **The spinning Regge-moment formula itself**: the analogue of CHMRS (3.17)–(3.20), (2.56), (2.59) for
   `⟨TTTT⟩` (or `⟨JJJJ⟩`) in `d = 3` — the heavy `s`-channel spinning block at fixed `β` in the helicity basis, the
   definition of functionals `C^{(m_L,m_R)}_{k,ν}`, and their swappability. Everything in §5–§7 marked [I] is
   conditional on it. The cleanest test remains: one heavy spinning block (currents) at fixed `β`.
2. The both-ends structure for stress tensors (KPZ (5.45) is for currents); the statement that a heavy block
   is a point in impact-parameter space *with a polarisation tensor whose helicity matrix is positive*.
3. The normalisations of §6.2 for the flip channels (overlap argument only); the phase conventions of
   `c^{++}c^{+−}` against `Φ_{4,0}`; the sign conventions `(4,4)` versus `(4,−4)` beyond the flat-limit matching.
4. Corrections of relative order `1/m²` (they distinguish `d̃^J_{4,4}` from `d̃^J_{0,0}` in flat space and are needed
   where the `B_3`/`B_4` terms compete, §7.3 item 2); uniform remainders; the region `m ∼ M`; the endpoint `ν ∼ M`.
5. Discrete series: a weight-`(m,m)` function on `H²` localised at `L = β` has, besides the real-`ν` spectrum,
   components at `ν = ±i(k+½)`, `k < m` (the zeros of the polynomial of Theorem R). Whether a dispersive functional
   in a flip channel sees them was not examined.
6. Parity-odd structures in `d = 3` (not in CHP's basis); the relation of `β_1, β_3` to `⟨TTO⟩` OPE data and to
   `(t_4, α_T)` in `d = 3` (CHP's (C.18)–(C.23), (D.15) are in a Mathematica notebook that was not processed).
7. The light side (graviton pole with the non-minimal coupling in the flip channels, double traces), and
   the positivity of the flat functionals themselves, known only numerically.

## 9. Code index (`ads4/code/regge/`)

| script | content | result |
|---|---|---|
| `01_omega_d3.py` | `Ω_{iν}` of CHP and KPZ in `d = 3` vs `ν tanh πν/(4π) P_{−1/2+iν}` | agree to `10^{−30}` |
| `02_mixed_sector.py` | covariant derivatives on `H²`, CHP (2.54)–(2.55) in the helicity basis | tables of §3.3, exact |
| `03_bitensor.py` | Lemma (sympy, `m ≤ 5`); Theorem R, `m ≤ 4` | max rel. dev. `6·10^{−24}` |
| `04_moments_and_channels.py` | flat `d̃_{4,4}` vs `d̃_{0,0}`; moments `Q_n^{(m)}`; tables | §5, §7.1 |
| `05_chp_crosscheck.py` | reading of `D̂_4`, `D̂_3` against CHP (4.14), (E.1), `d = 3, 4, 6` | 12 digits |
| `06_fg3_reading.py` | bracketing of CLPS (3.3a) by its positivity | reading C |
| `07_small_beta_normalisation.py` | `β → 0` limits fixing the normalisations | §6.2 |
| `08_fg3_curved_prediction.py`, `09_fg3_curved_scan.py` | conditional model of `F_{g3}` in the curved regime (output in `out08.txt`) | §7.3 table and scan |

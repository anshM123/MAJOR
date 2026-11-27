# A2′ (detector decoupling) tested on tree-level Einstein gravity in AdS₄

27 Nov 2026. Unrefereed. Nothing outside this directory was edited, and nothing was committed.
Inputs read (read-only): `prl/supplement.tex` (row A2′ of Table `tab:inputs`, l. 339; "Order of limits", l. 688–704; Sec.
`sec:angular`/`sec:angcross`, l. 2093–2190), `code/no_eft/NOTE.md` (§3: (ii′_E), (ii′), Lemma 4, Corollary A′, Theorem A),
`code/kernels/closed_forms.py`, `code/angproofs/{ap3_collinear,ap6_envelopes}.py` and `angproofs/cache/symbolic_structures.pkl`.

Grades: [computed] = done here; [cited] = taken from the literature; [estimated] = power counting, not computed.

## Verdict

**A2′ holds for tree-level Einstein gravity, with a large margin in spin.** The referee's worry has a definite answer. With
γ_{n,J} ~ 1/J, the weighted detector sum converges absolutely. Its sum is O(1/c_T), so 𝓔_light → 0, and c_T·sup_J γ² → 0
uniformly.

Two defects turned up in the *supporting lemmas* of the paper (not in A2′ for Einstein). Both are fixable. Details are in §4.

1. Lemma 4's extension to twist below 10 (h < 5), graded D-sketch, does not work with the collinear components it uses. The
   identities for F_{m₁m₃} with m ≤ −1 vanish identically exactly on the light double-twist families h = 1, 2, 3, 4 (even)
   and h = 3/2, 5/2, 7/2 (odd). The components with m ≥ +1 repair this for every family except n = 0 parity-even.
2. At twist 2 (h = 1) the two canonical parity-even structures coincide (rank 1). A second conserved twist-2 structure N(1)
   exists. Its detector amplitude is **not** at a zero. Its O(1/c_T) coefficient is a datum that (ii′), the OPE-data form,
   does not mention. The detector form (ii′_E), which is the actual A2′ of Table `tab:inputs`, is unaffected. For Einstein this
   datum is harmless, even in the worst case with no decay.

## 1. Anomalous dimensions at large J

| statement | grade |
|---|---|
| At J > 2 the tree-level [TT]_{n,J} data come from dDisc_t of the single T block. The t-channel double traces enter at O(c_T⁻²), and Einstein contact terms have finite spin support. | [cited] Caron-Huot 1703.00278; Kravchuk–Simmons-Duffin (spinning) |
| For every pair of collinear components (T_{m₁}, T_{m₂}), and so for every helicity/tensor structure: c_T γ_{n,J} = κ_{n,X}/J + O(J⁻²). The exponent is τ_T = 1. When the coefficient vanishes (component weight ½), the decay is faster. | exponent [computed] `s1_gamma_largeJ.py` (exact SL(2) inversion, slope −1.000004 at h̄ = 10⁵); coefficients [cited] Li–Meltzer–Poland 1511.08025 (3d, [TT]_{0,J}, negative iff collider bounds) |
| Impact-parameter form in AdS₄: γ ≈ −c·G h³/h̄ for h̄ ≫ h ≫ 1, i.e. κ_n ~ n³. Helicity flip is subleading. | [cited] Cornalba–Costa–Penedones–Schiappa hep-th/0611122, hep-th/0611123; KKSDZ 1904.05905 |
| Uniformity: at tree level c_Tγ_{n,J} does not depend on c_T and is finite at every J. Hence sup_J c_Tγ² = sup_J κ²/(c_T J²) = O(Δ_g^{p}/c_T) → 0. The first clause of (ii′) holds. | [computed]/[cited] |

## 2. Large-J decay of the detector amplitudes

Per light state, in the normalisation of Eq. (3.1) of NOTE.md:
|a_{O,b,0}|² = C_T·λ̂_X²·|K_X|². Here λ̂_X is the unit-normalised OPE coefficient on the canonical structure X, and K_X is
the closed-form kernel, which carries sin²(πh) (even) or cos²(πh) (odd).

| ingredient | result at fixed light twist, J → ∞ | grade |
|---|---|---|
| stripped kernel 4^{−Δ}2^{J}\|K_X\|²/(zero factor) | J^{−13/2} (pp, pm); J^{−15/2} (odd, oddJ), on every light family including h = 1, 3/2 | [computed] `s2_kernel_decay.py` |
| GFF OPE coefficients 4^{Δ}2^{−J}\|λ̂_X\|² | J^{1/2} (pp, pm); J^{3/2} (odd). From the 1d GFF coefficients c_n = (2a)_n(2b)_n/(n!(2a+2b+n−1)_n) and c_n 4^H ~ H^{2a+2b−3/2}, together with the collinear identities. Every component gives the same power. | [computed exponents] `s3_gff_envelope.py`; constants not tracked |
| zero factor | sin²(πγ/2) ≈ π²γ²/4 (likewise cos² on the odd families) | [cited] L1; [P] closed form |
| **product** | **\|a_{[TT]_{n,J},b,0}\|² ≈ C_T·A_n·J^{−6}·γ_{n,J}² ≈ A_nκ_n²·J^{−8}/c_T** | [computed exponents] |

This agrees with the Supplement's shell bound Σ_shell Σ_b|a|² ≲ S⁻⁵U², which is J⁻⁶ per state, now checked directly on the
light families.

## 3. The weighted sums

The summand of 𝓔_light is (1+J)^{5/2}log(2+J)·w_{O,b}·|a|², with w ≈ T_Δ/3 ≈ J²/3. It goes like:

* J^{5/2}·J²·C_T J^{−6}·γ² = c_Tγ²·J^{−3/2}log J;
* with γ = κ/(c_T J), that is **κ²J^{−7/2}log J/c_T**.

| quantity | result | grade |
|---|---|---|
| summand slope, q = 5/2 (pp h=1, pm h=1, odd h=3/2) | −3.39 including the log, i.e. J^{−7/2}log J | [computed] `s4_assemble.py` (i) |
| c_T·𝓔_light, families n ≤ 3, both ladders, odd n ≤ 1, J_max = 2000/8000, c_T = 10⁴, 10⁶, 10⁸ | 3755.31 … 3755.40: converged in J_max and independent of c_T. So **𝓔_light = O(1/c_T) → 0** | [computed]; A_n = κ_n = 1 [estimated constants] |
| with κ_n ∝ h_n³ (eikonal) | finite (7.3×10⁴ in the same units). The n-sum is finite at fixed Δ_g, so 𝓔 = O(poly(Δ_g)/c_T) | [estimated] |
| twist-2 N-structure added: B_J = 1/J / B_J = 1 (no decay) | c_T𝓔 = 4996.09 / 2.55×10⁴, both convergent (worst-case tail J^{−3/2}log J) | [computed given the estimate] |
| spin weight available, (1+J)^q | γ ~ 1/J: converges iff **q < 5**. Only sup_J c_Tγ² → 0, with no decay of γ: iff **q < 3**. The paper's q = 5/2 passes both. | [computed] (iv) |
| Theorem A consequence: c_T\|R\| ≤ (2πc_TGΛ⁴/3)·max(C_h, C′)·𝓔_light | O(C_h·poly(Δ_g)/c_T) → 0. Numerically meaningful only for c_T ≳ 10⁶⁹⁰ at Δ_g = 1000, as NOTE §4 already states | [D-sketch via L0] |
| clause 2 (other light states) | Tree-level ⟨TTTT⟩ contains only 1, T, [TT]_{n,J} (all structures, both parities), so it holds trivially at O(1/c_T). Multi-traces enter ⟨TTTT⟩ at O(c_T⁻²) (one loop); their uniformity in spin is outside what a tree-level check can test. | [cited] standard large-N counting |

No modification of the spin weight is needed. If one wants a weight that survives even without the 1/J decay of γ (i.e. only
the sup form of (ii′)), any q < 3 works, and 5/2 is inside it. The proofs (Theorem A via L0) need q = 5/2 with a log, fixed by
the functional side (L0: |ĥ| ≲ T_Δ(C_h + C′J^{5/2}log J)), not by the spectrum.

## 4. The two defects (in the lemmas, not in A2′ for Einstein)

**(a) Lemma 4 / Corollary A′ below h = 5** [computed, `s3_gff_envelope.py` (b)].

* The Supplement's collinear components are F_{−2,−2} ∝ p(h)(h−4), F_{−1,−1} ∝ p(h) and F_{−2,−1} ∝ (2h−7)(2h−5)(2h−3).
  They vanish at h ∈ {1, 2, 3, 4} and h ∈ {3/2, 5/2, 7/2}. These are exactly the parity-even [TT]_{n≤3} and parity-odd
  [TT]_{n≤2} twists.
* This is physical, not a normalisation artefact. T_m has collinear weights ((3−m)/2, (3+m)/2), so T_{m≤−1}×T_{m≤−1} cannot
  reach h < 5 in GFF.
* Inverting F → λ there gives |λ|² ≲ bound/γ², which cancels the double zero. So "for 1 ≤ h < 5 the same inequalities hold"
  is false as stated.
* Repair: use F_{2,2} (no h-zero), F_{1,1} (zero only at h = 1) and F_{1,2}/F_{2,1} (odd, no h-zero). They give the same per-state
  powers (s^{1/2} even, s^{3/2} odd), and hence the same S⁻⁵ shell envelope, for every light family except n = 0 parity-even.

**(b) Twist 2, parity even** [computed, `s3_gff_envelope.py` (c)].

* At h = 1 the canonical S^{++} and S^{+−} become proportional (rank 1 for J = 6, 12, 40):
  S^{+−} = c(J)·S^{++}.
* N(1) := lim_{h→1}(S^{+−} − cS^{++})/(h−1) is a second conserved twist-2 structure (13 non-zero monomials). Its b = 4 detector
  amplitude lim M₄^{+−}(h)/(h−1) is non-zero.
* In GFF its coefficient vanishes, because ⟨T|ℰℰ|T⟩_GFF = 0 by Wick's theorem with L[T]|Ω⟩ = 0 [D-sketch]. So the
  amplitude is still O(c_T^{-1/2}), and L1's conclusion at O(1/c_T) is unaffected.
* But |a_{[TT]_{0,J},4,0}|² ≈ C_T|B_J/c_T|²·J⁻⁶ is set by the O(1/c_T) coefficient B_J of N, not by γ_{0,J}. So (ii′) ⇒ (ii′_E)
  (Corollary A′) needs one more clause.
* For Einstein, B_J ~ 1/J relative to GFF [estimated, the same T-exchange mechanism as γ]. Even with B_J = O(1) the sum
  converges (§3).
* Also note: at h = 1 the parity-even double zero of L1 holds only for the S^{++}-direction. The N-direction has a simple
  O(1/c_T) amplitude with no γ factor.

## 5. Proposed Supplement wording

Row A2′, after "…uniformly in spin shells.", add:

> On the light $[TT]_{0,J}$ family the canonical even structures coincide at twist 2; the OPE-data form then also requires the $O(1/c_T)$ coefficient $B_J$ of the second conserved twist-2 structure to obey $c_T|B_J/c_T|^2\to0$ uniformly in $J$. Tree-level Einstein gravity satisfies A2$'$: $c_T\gamma_{n,J}=\kappa_{n}/J+O(J^{-2})$ for every structure ($\tau_T=1$), the stripped detector amplitudes of the light double twists fall as $J^{-6}$, so the summand of $\mathcal E_{\rm light}$ is $O(J^{-7/2}\log J/c_T)$ and $\mathcal E_{\rm light}=O(1/c_T)$ [computed exponents; constants not tracked]. The sum converges for any weight $(1+J)^q$ with $q<5$, and for $q<3$ using only $\sup_J c_T\gamma^2\to0$.

Sec. `sec:angcross`, "Consequences", after "…$\mathcal D=0$.", add:

> For $h<5$ the components $F_{-2,-2},F_{-1,-1},F_{-2,-1}$ vanish identically at the double-twist values $h=1,2,3,4$ (even) and $\tfrac32,\tfrac52,\tfrac72$ (odd), because $T_{m\le-1}$ has collinear weight $(3-m)/2\ge2$. There one uses $F_{2,2}$, $F_{1,1}$ and $F_{1,2}$ (crossing powers $S^{19/2}$, $S^{15/2}$, $S^{17/2}$), whose identities have no zero for $h\ge2$ (even) or $h\ge\tfrac32$ (odd). They give the same per-shell bounds. At $h=1$ the two even structures degenerate and only $\lambda_{+-}-2\lambda_{++}$ is bounded.

NOTE.md Lemma 4, "Light twist": replace "For 1 ≤ h < 5 the same inequalities hold with explicit Gamma-function constants
[D-sketch]" with the m ≥ +1 statement above, and add the B_J clause to (ii′).

## Files

| file | content | runtime |
|---|---|---|
| `common_a2.py` | read-only path setup; writes only to `out/` | — |
| `s1_gamma_largeJ.py` → `out/s1_gamma_largeJ.txt` | τ_T = 1 exponent of the single-block inversion for every component weight | < 1 s |
| `s2_kernel_decay.py` → `out/s2_kernel_decay.txt` | stripped kernel exponents; factorised vs closed form, agreement 10⁻⁵⁰ | ~ 1 s |
| `s3_gff_envelope.py` → `out/s3_gff_envelope.txt` | 1d GFF coefficients (closed form vs forward solve, 10⁻³³); h-zeros of all 50 collinear components; per-state λ² powers; twist-2 degeneracy and N(1) amplitude | ~ 5 s |
| `s4_assemble.py` → `out/s4_assemble.txt` | weighted sums, c_T scaling, N-structure, available spin weight | ~ 100 s |

Run each with `python3 <script>` from this directory.

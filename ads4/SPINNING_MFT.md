# Mean-field ("unitarity") normalisation of the spinning heavy formula: exact MFT OPE coefficients of `[Jφ]`, `[Tφ]` in `d = 3`

Status: written 16 Nov 2026 by a derivation/numerics agent; **not refereed**. Answers item 4 of `SPINNING_E2E.md` §11 (and part of
item 2); section numbers of `SPINNING_E2E.md` refer to its revised version (after `REFEREE_SPINNING_E2E.md`). Code and recorded outputs: `ads4/code/mft/` (`sh run_all.sh`; index in §9). Imports (does not modify) `ads4/code/e2e/`.
References: Caron-Huot–Li, *Helicity basis for three-dimensional conformal field theory*, arXiv:2102.08160 ("CHL"; LaTeX in the
session scratchpad, labels quoted); Karateev–Kravchuk–Simmons-Duffin, *Harmonic analysis and mean field theory*, arXiv:1809.05111
("KKSD", via CHL §3.1 only); CHMRS = arXiv:2106.10274; Fitzpatrick–Kaplan arXiv:1111.6972 ("FK"). Chang–Landau–Simmons-Duffin
(2311.04271) is not needed beyond what `SPINNING_E2E.md` §9 already says.

Marks: **[D]** derived; **[DV]** derived and verified numerically; **[NE]** numerical evidence (read off numbers, no derivation);
**[I]** inferred; **[Q]** quoted. "Exact" numerical evidence = equality of exact rationals.

## 0. Summary

1. **Exact MFT OPE coefficients, first principles [DV for the method; NE (exact) for the closed forms].** The coefficients `p^MFT` of
   `⟨OφOφ⟩_MFT = ⟨OO⟩⟨φφ⟩` (`O = J, T`) were computed **without conformal blocks**, as norms of two-particle primary states in the Fock
   space of MFT (§2). Spectrum: exactly one primary for each `(N, J)`, `l ≤ J ≤ N + l` (`Δ = Δ_O + Δ_φ + N`), parity even for
   `τ − τ_0 ∈ 2ℕ`, odd for `τ − τ_0 ∈ 2ℕ + 1` (`τ_0 = Δ_φ + 1`). All 145 coefficients computed (`N ≤ 8` / `6`, `Δ_φ = 7/4, 5/2, 13/10`) equal,
   as exact rationals, the closed forms (with `P_FK(Δ_1,Δ_2;Δ,J)` the Fitzpatrick–Kaplan coefficient continued in `n`,
   `H_l := (−J)_l(Δ−1)_l/((J+1)_l(2−Δ)_l)`):

   | | structure | `p^MFT` |
   |---|---|---|
   | `E`, `l = 1` | `𝔖^cons` of `SPINNING_E2E.md` §1.2 | `2^J P_FK(1, Δ_φ) · H_1 · (Δ−Δ_φ−J−1)²/4` |
   | `E`, `l = 2` | `𝔖^cons` | `2^J P_FK(1, Δ_φ) · H_2 · [(Δ−Δ_φ−J−1)(Δ−Δ_φ−J−3)]²/96` |
   | `O`, `l = 1` | `𝔒 = ε(Z_1,P_1,Z_5,P_2,P_5) v_5^{J−1}/(…)` (§3.2) | `2^{3J−3} P_FK(2, Δ_φ) · H_1` |
   | `O`, `l = 2` | §3.2 | `2^J P_FK(2, Δ_φ) · H_2/12` (in the frame normalisation of §3.2) |

   i.e. scalar MFT data of an **effective scalar of dimension `τ_O = 1` (even) or `τ_O + 1 = 2` (odd)** times CHL's helicity factor
   `H_l`, as CHL found for two conserved currents [Q, CHL eq. `free data`]. Not derived in closed form here (§3.3).
2. **Heavy limit [D from the closed forms + DV].** `P_FK(Δ_1,Δ_2)/b^{(Δ_12)}` is an exact ratio of four Gamma functions (§4.1), hence
   `p^MFT` is asymptotically the inverse Lorentzian shadow coefficient times explicit powers. All `L_0` factors and the "normalisation
   polynomials" `Π[(J+½)²−(i−½)²][(Δ−3/2)²−(i−½)²]` of `SPINNING_E2E.md` §6.3/§7.2 cancel **exactly** against the MFT density.
3. **MFT-normalised formula [D; constants DV].** With `q^{(P)}_{Δ,J} := (2n_J(2Δ−3))/(π m² p^{MFT}_P(Δ,J))` (CHMRS (2.32) with the MFT
   coefficient of the same parity sector `P = E, O`) and the functional normalised as CHMRS do, with `Δ_φ → (Δ_O, Δ_φ)`,
   `C^{(q_1q_3)}_{k,ν} := a_{Δ_O,Δ_φ}/(γ_{2Δ_φ+k−1}(ν)γ_{2Δ_O+k−1}(ν)) · (−1)^{J+l}Π̂^{(q_1q_3)}_{k,ν}/(q^{(P)}·2sin²)`, `a_{Δ_1,Δ_2} := 4Γ(Δ_1)Γ(Δ_1−½)Γ(Δ_2)Γ(Δ_2−½)`:

   **`C^{(q_1q_3)}_{k,ν}[Δ,J] = (1/(4m^{2k})) Σ_{h_1,h_3=±} c^{h_1} c̄^{−h_3} M^{(q_1q_3)}_{h_1h_3}(ν) F^ν_{l|h_1+h_3|}(β) · (1 + O(1/m²))`**,

   `F^ν_0 = P_{−½+iν}(cosh β)` (helicity preserved, `h_3 = −h_1`), `F^ν_{2l}` (flipped), `M` the table of `SPINNING_E2E.md` §4.4, and
   MFT-normalised helicity amplitudes `c^± = ĉ_E` (even), `c^± = ±ĉ_O` (odd), `|ĉ_P|² = p/p^MFT_P` (`= 1` in MFT). Equivalently
   `(1/(4m^{2k}))[(|c^+|²+|c^−|²)·½M_pres P + 2Re(c^+c̄^−)·½M_flip F^ν_{2l}]` for parity-definite exchanges (the `S`-parts of `M` cancel for
   both parities, §4.4).
   Status: even parity `l = 1, 2` **[D from closed form; DV]** (the normalised leading formula equals this × an analytic factor `Q_E → 1`,
   checked to 15 digits in 72 cases, and end-to-end on the exact moments of `SPINNING_E2E.md`, §5); odd parity `l = 1` **[D, vertex
   DV pointwise; overall sign fixed by positivity, I]**; odd `l = 2` **[ND]** (the MFT data exist; the heavy vertex was not computed).
4. **Answer to "what replaces `|c^{++}|² ± |c^{+−}|²`".** For one spinning pair: the constant is **`1/4` per helicity pair `(h_1,h_3)`
   per operator**, relative to the scalar normalisation with `a_{Δ_O,Δ_φ}`; the preserved component carries `|c^+|² + |c^−|²`, the flipped
   one `2Re(c^+c̄^−)` (= `|ĉ_E|²` resp. `−|ĉ_O|²`): MFT is helicity diagonal, flip contributions of `E` and `O` cancel on the MFT
   density. Spin-traced (kernel `(z_1·z_3)^l`, `M_{h,−h} = 1`, `M_{h,h} = 0`) and summed over parities, the MFT heavy density of `[Oφ]` is
   **exactly that of scalar double twists `[χφ]` with `Δ_χ = Δ_O`** (`P/m^{2k}` per twist interval 2) [D for `l = 1`; for `l = 2` the
   even half D, odd half I]. In the normalisation `C^{pres} := (8/M_pres)C`, each MFT operator gives exactly `2P/m^{2k}` — the CHMRS
   form; the "2" is then the helicity count, not the identical-particle factor. For `⟨TTTT⟩` this suggests `C^± = (2P/m^{2k})(…)` with
   per-pair factors `1/4` absorbed in the kernel normalisation [I; not computed].
5. **`O(1/m²)` (Goal 2) [D + NE].** MFT normalisation removes the "pol" term of `SPINNING_E2E.md` §6.3 identically [D]. After it, the
   helicity-preserving `c_2` of `l = 2` minus that of `l = 1` equals `−6P′/P` (`= −2(2²−1²)P′/P`) to `≤ 0.07` for `β = 0.4` and `β = 1`
   separately (three `ν` each, no fitted parameter) [NE]. Equivalently, per-`β` fits of the `l = 2` spin term with the `l = 1` data as
   reference give `d_β = −7.996` (`β = 0.4`), `−8.017` (`β = 1`) for the coefficient of `P′/P` (prediction `−8`) — this **resolves the referee's
   point M1** on `SPINNING_E2E.md` §6.3 (`d_β = −8.34, −10.79` there): the `β`-dependence came from using the scalar skeleton `(3, Δ_φ)` as
   reference. Moreover `c_2^{MFT}(l) + 2l²P′/P` equals the MFT-normalised `c_2` of the skeleton `(2, Δ_φ)` to `≤ 0.15` (`β ≤ 1`) for both
   `l` [NE]; the fitted "extra `l = 2` piece" is the difference of the skeletons `(2,Δ_φ)` and `(3,Δ_φ)`, not a spin effect. The term
   `−2l²P′/P` itself remains **not derived** (inferred from flat space; now confirmed at fixed `β` for both `l`).


## 1. Conventions

As `SPINNING_E2E.md` §1 (d = 3, `Δ_φ = 7/4` unless stated, `A_+ = Δ−3/2`, `A_− = J+½`, `α = A_+−A_−`, `ᾱ = A_++A_−`, `m² = αᾱ`,
`cosh β = (α²+ᾱ²)/(2m²)`, `γ_a(ν) = |Γ((a−½+iν)/2)|²`, `b^{(Δ_34)}_{Δ,J}` the Lorentzian shadow coefficient of §4.1 there, `𝔖^cons`, `L_n`,
`L_0` of §1.2 there, structures evaluated with `|s_ij|^κ`). Unit normalisations `⟨φφ⟩ = 1/x^{2Δ_φ}`, `⟨O(x,z)O(0,z′)⟩ = (z·I(x)·z′)^l/x^{2Δ_O}`;
the MFT correlator is `⟨O_1O_3⟩⟨φ_2φ_4⟩` and the channel is `(12) → (34)`.
`p_P(Δ,J) := λ²` where `⟨O_1φ_2X⟩ = λ·(structure)` for unit-normalised `X`; by `SPINNING_E2E.md` §1.4 [DVP there] the even exchange
contributes `λ²(−½)^J W_diff/L_0²`. `P_FK(Δ_1,Δ_2;Δ,J)` is FK's coefficient (as in `SPINNING_E2E.md` §2.1) with `n = (Δ−Δ_1−Δ_2−J)/2`
continued through Gamma functions. `n_J = 16π(2J+1)` (CHMRS (2.32) in d = 3).

**Frame at infinity.** `X(∞,ζ) := lim_{L→∞} L^{2Δ} X(Le, I(e)ζ)`, so `⟨X(∞,ζ′)X(0,ζ)⟩ = (ζ′·ζ)^J`. In this limit (`x := x_1 − x_2`)
`S_n → (−1)^{l−n}(−½)^n (z_1·x)^{l−n}(z_1·ζ)^n(ζ·x)^{J−n}(x²)^{(N−J−l+2n)/2}` [D; from `V_1 → −z_1·x`, `H_15 → −½L²z_1·ζ`,
`V_5 → L²ζ·x/x²`]; the structure is polynomial in `x` iff `N − J − l + 2n ∈ 2ℕ` for every `n` present.


## 2. Method: two-particle Fock space (`fock.py`) [D; DV]

* One-particle states of `O` (resp. `φ`) are represented by polynomial wavefunctions `f(x,z) = ⟨f|O(x,z)|0⟩`; their inner product is the
  one whose reproducing kernel is the radially-quantised two-point function
  `K_O(y,w;x,z) = ⟨0|O(y,w)^†O(x,z)|0⟩ = (w·I(y)I(y/y²−x)·z)^l/D^{Δ_O} = A^l/D^{2l+1}`, `D = 1−2x·y+x²y²`, `K_φ = D^{−Δ_φ}`
  (`O(y,w)^† = |y|^{−2Δ}O(y/y², I(y)w̄)`; `A` is a polynomial of bidegree `≤ (2,2)`, computed symbolically). If `K = Σ y^α M_{αβ} x^β`
  (polarisations reduced mod `z² = 0`), then `‖ψ‖² = ψ^T M^+ ψ̄` on the range of `M` (the range is the conserved subspace for `O`; the code
  asserts that every wavefunction used lies in it) [D: `K = f(y)^†G^{−1}f(x)` for any basis of states].
* A primary `X` of dimension `Δ_O + Δ_φ + N`, spin `J` has wavefunction `q_X = ⟨0|X(∞,ζ)O(x_1,z_1)φ(x_2)|0⟩`, which is translation invariant
  (`⟨X|P = 0`), i.e. a polynomial three-point structure `T(x_1−x_2)` conserved in `x_1`; conversely each such structure is a primary.
  With `‖X,ζ̄‖² = N_X(ζ·ζ̄)^J`: **`p = (ζ·ζ̄)^J/‖𝒮(∞)‖²_Fock`** for the structure `𝒮`. Everything is exact (python-flint rationals; `ζ = (1,i,0)`).
* Checks: (i) scalar × scalar (`l = 0`, `Δ_1 = 5/4`, `Δ_2 = 7/4`, `n ≤ 2`, `J ≤ 3`): `p = 2^J P_FK` exactly (12 cases), which also confirms the
  `(−½)^J` of `SPINNING_E2E.md` §1.4; (ii) `N = 0`: `‖ζ·z‖² = 2`, `‖(ζ·z)²‖² = 4` (unit-normalised `Oφ`); (iii) the null space of the
  divergence on the even polynomial structures is one-dimensional at every `(N,J)` and proportional to `Σ_n L_n S_n(∞)` with the `L_n` of
  `structures.py` (an independent check of the conserved structures, at infinity); (iv) the counts below agree with one primary per `(N,J)`.

## 3. Exact MFT data and closed forms (`mft_exact.py`, `out_mft_exact.txt`)

### 3.1 Spectrum [NE (exact), all `N ≤ 5`]
Exactly one conserved polynomial structure, hence one primary, for each `(N,J)` with `l ≤ J ≤ N + l`, none for `J < l` or `J > N + l`.
It is parity even iff `N − J − l` is even, i.e. `τ − τ_0 ∈ {0, 2, 4, …}`, and odd for `τ − τ_0 ∈ {1, 3, …}`. At `τ = τ_0` (and `τ_0 + 2` for
`l = 2`) the even structure is the pure-`H` one (`L_0 = 0` there; `𝔖^cons` is singular, `𝔏 := Σ_n L_n S_n = L_0𝔖^cons` is not).

### 3.2 Closed forms [NE: exact rational agreement, 145/145 primaries]

`l = 1`: 25 + 12 + 12 even and 20 + 9 + 9 odd primaries (`Δ_φ = 7/4, 5/2, 13/10`); `l = 2`: 16 + 9 + 9 even, 12 + 6 + 6 odd.

* even: **`p[𝔏] = 2^J P_FK(1,Δ_φ;Δ,J) / (c_l (J−l+1)_{2l} (Δ−l−1)_{2l})`, `c_1 = 4`, `c_2 = 24`** (`(J−l+1)_{2l} = Π_{i≤l}[(J+½)²−(i−½)²]`,
  `(Δ−l−1)_{2l} = Π_{i≤l}[(Δ−3/2)²−(i−½)²]`); `p[𝔖^cons] = L_0² p[𝔏]`, giving the table of §0.
  Sample (`l = 1`, `Δ_φ = 7/4`, `p[𝔏]`): `(N,J) = (0,1): 4/121`, `(2,1): 88/70395`, `(4,3): 33440/687664971`; (`l = 2`): `(0,2): 16/731025`.
* odd, `l = 1`: structure at infinity `ε(z_1,ζ,x)(ζ·x)^{J−1}(x²)^{(N−J)/2}` (unique, automatically conserved):
  **`p = 2^J P_FK(2,Δ_φ) J(Δ−1)/(2(J+1)(Δ−2))`**. Covariantly `𝔒 := ε(Z_1,P_1,Z_5,P_2,P_5)v_5^{J−1}/(|s_12|^α|s_15|^β|s_25|^γ)` with the
  exponents of `S3(Δ_1=2,Δ_φ)`; `𝔒(∞) = −(−½)^{J−1}ε(z_1,ζ,x)(ζ·x)^{J−1}|x|^{N−J}` [DV, `odd_checks.py` (ii)], hence `p[𝔒] = 2^{3J−3}P_FK(2,Δ_φ)H_1`.
  (`ε(A,…,E)` := determinant of the upper components `(P⁺,P⁻,x)`.)
* odd, `l = 2`: the conserved structure is `(Δ−Δ_φ)ε(z_1,ζ,x)(z_1·ζ)(ζ·x)^{J−2}|x|^{p+2} − (Δ−Δ_φ−J−2)ε(z_1,ζ,x)(z_1·x)(ζ·x)^{J−1}|x|^p`,
  `p = N−J−1` (coefficient read off and checked exactly on all 24 cases), with **`p = 2^J P_FK(2,Δ_φ) H_2/12`**.

The dependence on `Δ_φ` is entirely in `P_FK` (`p/P_FK` is `Δ_φ`-independent at fixed `(N,J)`, checked at three `Δ_φ`).

### 3.3 Status of the closed forms
They were **identified** from the exact data (rational functions recognised on the data, no continuous parameter fitted) and then
confirmed exactly on all computed points; they are **not derived** here. Supporting evidence of a structural kind: CHL obtain, by
Euclidean inversion in the helicity basis, MFT data of two conserved currents equal to scalar MFT data with external dimensions equal to
the twists (`c(1,1)` even, `½c(1,2)` odd) times `diag[(−J)_{J_1+J_2}(Δ−1)_{|J_1−J_2|}/…]` [Q, CHL eqs. `free data`, `evenOPEcoeff`];
the forms above are the one-scalar analogue (`(−J)_l(Δ−1)_l/((J+1)_l(2−Δ)_l)` is their matrix entry with `J_2 = 0`), and the effective
dimensions `1` (even) and `2` (odd) are exactly those of the seeds whose shadow coefficients govern the spinning blocks: `b^{(1−Δ_φ)}`
for the even conserved structure (`SPINNING_E2E.md` §2.3, §4.2) and `b^{(2−Δ_φ)}` for the odd one (§4.4 below). A proof would follow from
KKSD's pairing formula (three shadows and a pairing, CHL eq. `cMFT`) with a hard-shadow coefficient for a generic scalar leg, which CHL
only "observe" for two currents. CHL's factor-2 issue (spurious poles at `n ≥ 1`) does not arise in the Fock computation.


## 4. Heavy limit (`mft_heavy.py`, `out_mft_heavy.txt`)

### 4.1 Scalar identity and the mixed-scalar MFT normalisation [DV]
With `S := Δ_1+Δ_2`, **`P_FK(Δ_1,Δ_2;Δ,J)/b^{(Δ_12)}_{Δ,J} = (2/π) A_+A_− m² G(α)G(ᾱ)/(Γ(Δ_1)Γ(Δ_1−½)Γ(Δ_2)Γ(Δ_2−½))`**, `G(x) = Γ((x+S−2)/2)/Γ((x+4−S)/2)`,
exactly (the four Gamma functions with `±Δ_12/2` cancel between `P_FK` and `b`; check (1): `10^{−29}` at random points). Since
`G(x) = (x/2)^{S−3}(1 − (4/3)B_3((S−2)/2)/x² + O(x^{−4}))` (no `1/x` term because the two shifts add to 1), for the mixed scalar `[φ_1φ_2]`
with `q := 2n_J(2Δ−3)/(πm² P_FK)`:
`(−1)^JΠ̂_{k,ν}[block]/(q·2sin²)` per unit `P_FK`-normalised coefficient `→ (γ_{2Δ_2+k−1}γ_{2Δ_1+k−1}/a_{Δ_1,Δ_2}) · P_{−½+iν}(η)/m^{2k}`,
with the exact finite-`m` factor `(4/m²)^{S−3}G(α)G(ᾱ) = 1 − (8/3)B_3((S−2)/2) cosh β/m² + O(m^{−4})` (check (2): five `(Δ_1,Δ_2)`, `m` up to
4500). For `Δ_1 = Δ_2` this is CHMRS (2.59) (`a_{Δ_φ,Δ_φ} = a_{Δ_φ}`; their `p^MFT = 2P_FK` gives the `2`); the mixed case has `1·P`.

### 4.2 Asymptotics of `p^MFT` [D from §3.2]
`p^MFT_E[𝔖^cons] = b^{(1−Δ_φ)}_{Δ,J} · 2^{J+1} L_0² A_+A_−m² G(α)G(ᾱ) / (π c_l Π_l Γ(½)Γ(Δ_φ)Γ(Δ_φ−½))` (`S = 1+Δ_φ`, `Π_l` the product of §3.2),
`→ b^{(1−Δ_φ)} 2^{J+1} 4^{2−Δ_φ} L_0² m^{2Δ_φ−2}/(π^{3/2}c_l (A_+A_−)^{2l−1} Γ(Δ_φ)Γ(Δ_φ−½)) (1+O(m^{−2}))`;
`p^MFT_O[𝔒] = b^{(2−Δ_φ)} 2^{3J−2} H_1 A_+A_−m²G(α)G(ᾱ)/(π Γ(2)Γ(3/2)Γ(Δ_φ)Γ(Δ_φ−½))` (`S = 2+Δ_φ`, `l = 1`); `H_l → 1`; the `b`'s are exactly the shadow coefficients that appear in the heavy formulas (§4.3, §4.4), which is why everything cancels.

### 4.3 Parity even: MFT-normalised heavy formula [D; DV]
Multiply the leading formula of `SPINNING_E2E.md` §4.5 (per `λ²`: `|T_l|/(2^J L_0²)`, prefactor `(π/4)4^{l+1+Δ_φ}/b^{(1−Δ_φ)}`,
`m^{−a_x−a_y}`, `a_y = 2Δ_O+k−1`) by `1/q^{(E)}`. `2^J`, `L_0²` and `b` cancel; with §4.1 the result is
`(γ_{a_x}γ_{a_y}/m^{2k}) · Y_l · [M_pres P + M_flip F^ν_{2l}] · (4/m²)^{Δ_φ−2}G(α)G(ᾱ)/a_{1,Δ_φ}`,
`Y_l = 4^l|T_l|/(m^{2l}c_l Π_l) = w_l (A_+A_−)^{2l}/Π_l`, `w_1 = ½`, `w_2 = 1/6` (using `𝖭_1 = A_+A_−m`, `𝖭_2 = ½(A_+A_−)²m²` of
`SPINNING_E2E.md` §4.3). Since `a_{Δ_O,Δ_φ}/a_{1,Δ_φ} = l!(½)_l` and **`w_l·l!(½)_l = ¼` for `l = 1, 2`**:

`C^{(q_1q_3)}_{k,ν}[Δ,J] = Q_E · ¼[M^{(q_1q_3)}_pres P_{−½+iν}(cosh β) + M^{(q_1q_3)}_flip F^ν_{2l}(β)] · |ĉ_E|²`, exactly at leading order, with
`Q_E = (4/m²)^{Δ_φ−2}G(α)G(ᾱ) · (A_+A_−)^{2l}/Π_l = 1 + [κ_l(1/sinh²(β/2)+1/cosh²(β/2)) − (8/3)B_3((Δ_φ−1)/2)cosh β]/m² + …`, `κ_1 = ¼`, `κ_2 = 5/2`.

Check (3): `predict.prediction` (the code verified end-to-end in `SPINNING_E2E.md`) × `1/q^{(E)}` × `a/(γγ)` × `m^{2k}`, divided by
`¼[M_pres P + M_flip F]` and by the analytic `Q_E`: `1.000000000000000` in all 72 cases (`l = 1, 2`; `β = 0.4, 1, 2`; `m ≈ 20, 100, 1000`;
`k = 4, 6`; `ν = 0, 3`); `m²(Q_E−1)` converges to the stated coefficient (e.g. `l = 2`, `β = 0.4`: `72.12, 64.34, 63.995 → 63.991`).
Consequences: (i) the factor `(A_+A_−)^{2l}/Π_l` in `Q_E` is the inverse of the exact normalisation-polynomial ratio that
`SPINNING_E2E.md` found in the exact moments (§6.3 "pol", §7.2 fixed-`J` factor), so **the MFT-normalised density has neither** — it is
uniform down to `J = l` at leading order [D, given those findings]; (ii) the `ν`-, `k`-, `β`-dependence is unchanged.

### 4.4 Parity odd, `l = 1` [D; pieces DV; sign I]
* Differential basis [DV]: with `D̃ := ε(Z_1,P_1,∂_{P_1},P_2,∂_{P_2})` (spin `+1` at point 1, dimensions unchanged),
  `D̃ S3(Δ_1,Δ_2;Δ,J) = −J(Δ−1) ε(Z_1,P_1,Z_5,P_2,P_5) S3/v_5`, i.e. `D̃S3(2,Δ_φ) = −J(Δ−1)𝔒` (hand derivation; `odd_checks.py` (i) to 20
  digits at two points).
* Shadow coefficient [D]: applying `D̃^{(12)}D̃^{(34)}` to the scalar Lorentzian shadow representation of the seed `(2,Δ_φ,2,Δ_φ)`
  (`SPINNING_E2E.md` §4.1, DVE there) gives the odd partial wave `W^O := D̃D̃[K g]` with coefficient `1/b^{(2−Δ_φ)}` and right structure
  `D̃S3^S = −(J+1)(Δ−2)𝔒^S` (a single word, so no conservation-consistency condition is needed). By the OPE-limit argument of
  `SPINNING_E2E.md` §1.4 the contribution is `λ_𝔒²(−1)^J 8^{−J}W^O/(J(Δ−1))²` [D, analogous; not tested numerically].
* Heavy vertex [DV pointwise]: in the frame of `SPINNING_E2E.md` §4.3, `ε/v_5 → −2 z·e_⊥` at both vertices (`odd_checks.py` (iii):
  `−1.99644, −1.99996, −2.0000000` at `m = 10, 100, 10⁴`); `z·e_⊥ = (z·ε_+ − z·ε_−)/(i√2)`: odd, transverse, helicity `±1` (Lemma V).
  So the tensor is `4J(J+1)(Δ−1)(Δ−2)(z_1·e_⊥)(z_3·e_⊥) = −2J(J+1)(Δ−1)(Δ−2)Σ_h sgn(h_1)sgn(h_3)(z_1·ε_{h_1})(z_3·ε_{h_3})`.
* `M` with weights `sgn(h_1)sgn(h_3)`: `S^{q_1q_3}_{−h_1,−h_3} = −S^{q_1q_3}_{h_1h_3}` exactly (check (5), all modes, `l = 1, 2`), so the `S`-part
  cancels and the result is `−M_pres P + M_flip F` up to the overall factor.
* Product with `p^MFT_O[𝔒] = 2^{3J−3}P_FK(2,Δ_φ)H_1`: `[2J(J+1)(Δ−1)(Δ−2)/(8^J J²(Δ−1)²)] × 2^{3J−3}P_FK(2,Δ_φ)J(Δ−1)/((J+1)(Δ−2)) = ¼P_FK(2,Δ_φ)`:
  all `J, Δ` factors cancel exactly; with §4.1 (`S = 2+Δ_φ`, `a_{2,Δ_φ} = a_{Δ_O,Δ_φ}`): **`C^{(q_1q_3)} = Q_O·¼[M_pres P − M_flip F^ν_2]|ĉ_O|²`**,
  `Q_O = (4/m²)^{Δ_φ−1}G(α)G(ᾱ) = 1 − (8/3)B_3(Δ_φ/2) cosh β/m² + …` (check (4)); no normalisation polynomial at all.
* Sign [I]: the Euclidean MFT computation uses the real Euclidean `ε`; the real-section `ε` is its continuation times `±i`, so `λ²` changes
  sign between the two, and `(z_1·e_⊥)(z_3·e_⊥)` carries `−sgn·sgn`. The overall sign above is the one for which the helicity-preserving
  component is positive (reflection positivity of the preserved, `1↔3`-symmetric kernel); the relative sign preserved/flip (`−`) is derived.

### 4.5 What is not done for parity odd, `l = 2`
The MFT data and closed form exist (§3.2); the differential-basis odd structure (`D̃` composed with `D_11`, `D_12` on seeds `(4,Δ_φ)`,
`(3,Δ_φ+1)`), its conserved combination, shadow coefficient (expected `b^{(2−Δ_φ)}`) and heavy vertex were not computed. Expected by
analogy: the same `¼[M_pres P − M_flip F^ν_4]` [I].


## 5. End-to-end on the exact moments of `SPINNING_E2E.md` (check (6))

Since both sides of `SPINNING_E2E.md` §5 use the same `M`'s, the MFT-normalised ratio is `R^{MFT} = R_{e2e} × Q_E` exactly; for the
helicity components `R^{MFT}_{pres,flip} = R_{pres,flip} × Q_E`. Sample (`k = 6`, largest `m`; full table in `out_mft_heavy.txt`):

| `l` | `β` | `m` | `ν` | `R_pres` (e2e) | `R_pres^MFT` | `R_flip^MFT` | `m²(R_pres−1)` | `m²(R_pres^MFT−1)` |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.4 | 101.8 | 0 / 3 | 0.99072 / 0.99753 | 0.99133 / 0.99813 | 0.98535 / 0.98990 | −96.2 / −25.6 | −89.9 / −19.3 |
| 2 | 0.4 | 101.8 | 0 / 3 | 0.98529 / 0.99513 | 0.99140 / 1.00131 | 0.97446 / 0.97964 | −152.5 / −50.5 | −89.1 / +13.5 |
| 1 | 1.0 | 100.8 | 0 / 3 | 0.98760 / 0.97992 | 0.98770 / 0.98002 | 0.98193 / 0.99359 | −125.9 / −203.8 | −124.9 / −202.9 |
| 2 | 1.0 | 100.8 | 0 / 3 | 0.98669 / 0.97702 | 0.98777 / 0.97808 | 0.97581 / 0.98498 | −135.1 / −233.3 | −124.2 / −222.5 |

So the MFT-normalised formula `¼[M_pres P + M_flip F]` (no fitted constant) is approached as `1 + c_2/m²` by the exact `⟨JφJφ⟩`, `⟨TφTφ⟩`
moments, for both helicity components [DV, inherits the numerical accuracy of `SPINNING_E2E.md` §5.1]. After normalisation the
preserved `c_2` of `l = 1` and `l = 2` nearly coincide at `ν = 0` (−89.9 vs −89.1; −124.9 vs −124.2), the difference at `ν = 3` being the
helicity term (§6).


## 6. `O(1/m²)` of the helicity-preserving component (Goal 2) (`goal2.py`, `out_goal2.txt`)

Input: the measured `c_2` of `SPINNING_E2E.md` §6.2 (`k = 6`, their fits (B)) and the scalar-skeleton `c_2` of the same table. MFT
normalisation adds the derived `q_2 = κ_l(1/sinh²(β/2)+1/cosh²(β/2)) − (8/3)B_3((Δ_φ−1)/2)cosh β` (§4.3), i.e. it removes "pol".

| `β` | `ν` | `−6P′/P` | `c_2^{MFT}(l=2) − c_2^{MFT}(l=1)` | residual `c_2^{MFT}(l)+2l²P′/P − c_2^{MFT,skel(2,Δ_φ)}`, `l = 1 / 2` |
|---|---|---|---|---|
| 0.4 | 0 / 1 / 3 | 0.724 / 3.696 / 33.128 | 0.758 / 3.698 / 33.128 | −0.03, 0.003, 0.01 / 0.007, 0.005, 0.01 |
| 1.0 | 0 / 1 / 3 | 0.608 / 3.479 / −19.737 | 0.616 / 3.466 / −19.804 | −0.15, 0.04, −0.07 / −0.14, 0.03, −0.14 |
| 2.0 | 0 / 1 / 3 | 0.349 / 3.962 / −9.197 | 0.574 / 3.574 / −8.426 | −1.9, 1.3, −4.1 / −1.7, 0.9, −3.3 |

(`c_2^{MFT,skel}` = skeleton `c_2 − (8/3)B_3(Δ_φ/2)cosh β`; quoted fit errors of the inputs: `≲ 0.3` for `β ≤ 1`, a few units at `β = 2`.)

Per-`β` fits (the "direct test" of `SPINNING_E2E.md` §6.3, but with the `l = 1` data as reference): fit
`c_2^{MFT}(l=2) − [c_2^{MFT}(l=1) + 2P′/P] = c_β + d_β P′/P` over `ν = 0, 1, 3`:

| `β` | `d_β` (prediction `−8`) | `c_β` | max residual with `d = −8` |
|---|---|---|---|
| 0.4 | −7.996 | 0.02 | 0.022 |
| 1.0 | −8.017 | −0.01 | 0.043 |
| 2.0 | −7.51 | 0.07 | 0.59 |

Findings [NE, no parameter fitted]: (a) after MFT normalisation the `l`-dependence of the preserved `c_2` is `−2l²P′/P` to `≤ 0.07`
(`β = 0.4` and `1` separately), so the curved form `−2l²P′/P` holds for `l = 2` at fixed `β` once the `l = 1` (or skeleton `(2,Δ_φ)`)
reference is used — the `β`-dependent `P′/P` coefficient found for `l = 2` in `SPINNING_E2E.md` §6.3 (referee point M1) is an artefact of the
skeleton `(3,Δ_φ)` reference; (b) `c_2^{MFT}(l) = c_2^{MFT}[scalar skeleton (2,Δ_φ)] − 2l²P′/P` within the input errors, for both `l`
(`β ≤ 1`; consistent at `β = 2`). Hence the fitted "extra `l = 2` piece" `−2sinh²βP′/P − 17.3 cosh β` equals `c_2^{skel(2,Δ_φ)} − c_2^{skel(3,Δ_φ)}`
(+ the `l = 1` constants `0.07 + 0.10 cosh β`), e.g. `−18.5` vs `−18.49` at `β = 0.4`, `ν = 0`. The "pol" term is exactly the MFT normalisation [D].
Why the reference for both `l` is the skeleton with `Δ_1 = 2` (and not `Δ_O`) is not understood [observation only].
**Not derived**: the helicity term `−2l²P′/P` and the skeleton `c_2` itself (the analogue of CHMRS (3.27) for unequal dimensions and
tensor kernels); the flip `c_2` was not re-analysed. The mechanical route is still item 2 of `SPINNING_E2E.md` §11.


## 7. Answer to the question, stated precisely

For `⟨OφOφ⟩` (one spinning Regge pair, `O = J, T`) and the tensor Regge moments `Π̂^{(q_1q_3)}_{k,ν}`, with heavy densities defined as in
CHMRS (2.31)–(2.32) but with `p^MFT` of the exchanged operator's parity sector, and the functional normalised by `a_{Δ_O,Δ_φ}/(γ_{2Δ_φ+k−1}γ_{2Δ_O+k−1})`:

`C^{(q_1q_3)}_{k,ν}[Δ,J] → (1/(4m^{2k})) [ M^{(q_1q_3)}_pres(ν) P_{−½+iν}(η_AdS) (|ĉ_E|² + |ĉ_O|²) + M^{(q_1q_3)}_flip(ν) F^ν_{2l}(β) (|ĉ_E|² − |ĉ_O|²) ]`,

each term from operators of the indicated parity (`E` at `τ−τ_0` even, `O` odd). In helicity language the weight of the pair
`(h_1,h_3)` is `¼c^{h_1}c̄^{−h_3}`, so "`|c^{++}|² ± |c^{+−}|²`" is replaced, for one pair, by `|c^+|²+|c^−|²` (preserved) and `2Re(c^+c̄^−)` (flip),
each with the exact constant `1/8` times the summed transform factor. Constants: `¼` exact (`w_l l!(½)_l = ¼`, `l = 1, 2`); `a_{Δ_O,Δ_φ}`
(not `a_{1,Δ_φ}` and not `a_{Δ_φ}`); "`2P/m^{2k}`" is recovered per MFT operator in the normalisation `C^{pres} = (8/M_pres)C` and, spin-traced
and summed over parities, the density equals that of a scalar of dimension `Δ_O`. Status: `E`, `l = 1, 2` [D from NE-exact closed forms; DV];
`O`, `l = 1` [D; sign I]; `O`, `l = 2` [ND]; `⟨TTTT⟩` (two pairs: expected `1/16` per helicity quadruple and the identical-particle `2`) [I, not computed].


## 8. Status list

| # | statement | mark |
|---|---|---|
| 1 | Fock-space method for exact MFT OPE coefficients of `[Oφ]` with `O` conserved | [D]; validated on scalars (FK exact) and normalisations |
| 2 | Spectrum: one primary per `(N,J)`, `l ≤ J ≤ N+l`, parity `(−1)^{N−J−l}` | [NE, exact, `N ≤ 5`] |
| 3 | Closed forms of `p^MFT` (even `l = 1,2`; odd `l = 1,2`) | [NE, exact on 145/145]; not derived (CHL analogue [Q]) |
| 4 | `P_FK/b` Gamma identity; mixed-scalar MFT normalisation `1/a_{Δ_1,Δ_2}`, exact `1/m²` factor | [DV] |
| 5 | MFT-normalised even formula `¼[M_pres P + M_flip F]`, `Q_E` exact | [D from 3; DV to 15 digits; end-to-end DV] |
| 6 | MFT normalisation cancels the normalisation polynomials ("pol", fixed-`J` factor) | [D] |
| 7 | `D̃S3 = −J(Δ−1)εS3/v_5`; odd vertex `→ −2z·e_⊥`; `S`-cancellation | [DV] |
| 8 | Odd `l = 1` formula `¼[M_pres P − M_flip F]`, `Q_O` | [D]; shadow coefficient and OPE normalisation by analogy [D, not tested]; sign [I] |
| 9 | Odd `l = 2` heavy formula | [ND] |
| 10 | Spin-traced, parity-summed MFT density = scalar `(Δ_O,Δ_φ)` | [D] `l = 1`; [I] `l = 2` (odd half) |
| 11 | After MFT normalisation, `l`-dependence of preserved `c_2` is `−2l²P′/P` at fixed `β` (`d_β = −8.00, −8.02` for `l = 2`); extra `l = 2` piece = skeleton difference | [NE] |
| 12 | Analytic `−2l²P′/P` and skeleton `c_2` | [ND] |
| 13 | `⟨TTTT⟩` normalised form | [I] |

**What remains.** (i) A derivation of the closed forms (KKSD/CHL pairing formula with one generic scalar leg, or a closed evaluation of
the Fock norm); (ii) the odd `l = 2` heavy vertex and shadow coefficient (mechanical with `D̃` and the `e2e` word machinery), and a direct
numerical test of the odd shadow coefficient (adapt `t5_shadow_direct.py`); (iii) parity-odd blocks, to test §4.4 end-to-end and to
decompose the MFT correlator into blocks directly (not needed for §3, which is block-free); (iv) Goal 2 analytically (`SPINNING_E2E.md`
§11 item 2); (v) two spinning pairs.

## 9. Code index (`ads4/code/mft/`; `sh run_all.sh`)

| file | content | output |
|---|---|---|
| `fock.py` | Fock-space kernels, exact norms, structures at infinity, conservation | — |
| `mft_exact.py` | exact `p^MFT` for `l = 1, 2`, three `Δ_φ`, vs closed forms (cache in `cache/`, git-ignored) | `out_mft_exact.txt` |
| `odd_checks.py` | `D̃` identity; odd structure at infinity and in the heavy frame | `out_odd_checks.txt` |
| `mft_heavy.py` | checks (1)–(6) of §4–5 | `out_mft_heavy.txt` |
| `goal2.py` | §6 | `out_goal2.txt` |

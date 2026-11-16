# Mean-field ("unitarity") normalisation of the spinning heavy formula: exact MFT OPE coefficients of `[Jφ]`, `[Tφ]` in `d = 3`

Status: written 16 Nov 2026 by a derivation/numerics agent; **not refereed**. Answers item 4 of `SPINNING_E2E.md` §11 (and part of
item 2); section numbers of `SPINNING_E2E.md` refer to its revised version (after `REFEREE_SPINNING_E2E.md`). Revised after the
referee report `ads4/REFEREE_SPINNING_MFT.md` (no errors found; four label/attribution issues; changes listed in §10). Code and recorded outputs: `ads4/code/mft/` (`sh run_all.sh`; index in §9). Imports (does not modify) `ads4/code/e2e/`.
References: Caron-Huot–Li, *Helicity basis for three-dimensional conformal field theory*, arXiv:2102.08160 ("CHL"; LaTeX in the
session scratchpad, labels quoted); Karateev–Kravchuk–Simmons-Duffin, *Harmonic analysis and mean field theory*, arXiv:1809.05111
("KKSD", via CHL §3.1 only); CHMRS = arXiv:2106.10274; Fitzpatrick–Kaplan arXiv:1111.6972 ("FK"). Chang–Landau–Simmons-Duffin
(2311.04271) is not needed beyond what `SPINNING_E2E.md` §9 already says.

Marks: **[D]** derived; **[DV]** derived and verified numerically; **[NE]** numerical evidence (read off numbers, no derivation);
**[I]** inferred; **[Q]** quoted. "Exact" numerical evidence = equality of exact rationals.

## 0. Summary

1. **Exact MFT OPE coefficients, first principles [DV for the method; NE (exact) for the closed forms].** The coefficients `p^MFT` of
   `⟨OφOφ⟩_MFT = ⟨OO⟩⟨φφ⟩` (`O = J, T`) were computed **without conformal blocks**, as norms of two-particle primary states in the Fock
   space of MFT (§2). Spectrum (agrees with the character count, `REFEREE_SPINNING_MFT.md`): exactly one primary for each `(N, J)`, `l ≤ J ≤ N + l` (`Δ = Δ_O + Δ_φ + N`), parity even for
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
   `H_l`, as CHL found for two conserved currents [Q, CHL eq. `free data`]. Not derived in closed form here (§3.3). Literature (§3.4):
   the leading-twist (`τ_0`) even coefficients agree, up to the structure normalisation `4/L_1²` (exactly, all `J`, three `Δ_φ`), with
   Sleight–Taronna arXiv:1804.09334 eq. (`opejphi`); the all-`n` closed forms for conserved `O` in `d = 3`, in particular the parity-odd
   ones, and the heavy-limit cancellations of §4 were not found in the papers checked.
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
   Status: even parity `l = 1, 2` **[D from the closed forms]**; the bridge "Fock-space `p^MFT` ↔ block normalisation" on which it rests
   is **tested against the exact e2e blocks** by decomposing the MFT correlator itself (leading-twist family, all `J`, §5.1:
   ratio `1 + O(u)`, e.g. `1 + 1.7·10⁻⁷`); the remaining checks of the heavy algebra (§4.3 check (3), §5.3) are algebraic identities, not
   independent tests. Odd parity `l = 1` **[D]**, with its normalisation **and sign** fixed by a block-level OPE test (§5.2:
   `σ = −1.0002, −1.0003, −0.9987` at `δ = 10⁻⁴` for three polarisation pairs, prediction `|σ| = 1`), no longer by positivity; its shadow
   coefficient and heavy vertex are derived (vertex verified pointwise), not tested in moments. Odd `l = 2` **[ND]**.
4. **Answer to "what replaces `|c^{++}|² ± |c^{+−}|²`".** For one spinning pair: the constant is **`1/4` per helicity pair `(h_1,h_3)`
   per operator**, relative to the scalar normalisation with `a_{Δ_O,Δ_φ}`; the preserved component carries `|c^+|² + |c^−|²`, the flipped
   one `2Re(c^+c̄^−)` (= `|ĉ_E|²` resp. `−|ĉ_O|²`): MFT is helicity diagonal, flip contributions of `E` and `O` cancel on the MFT
   density. Spin-traced (kernel `(z_1·z_3)^l`, `M_{h,−h} = 1`, `M_{h,h} = 0`) and summed over parities, the MFT heavy density of `[Oφ]`
   equals that of scalar double twists `[χφ]` with `Δ_χ = Δ_O` (`P/m^{2k}` per twist interval 2), consistently with the exact identity
   "polarisation trace of `⟨O_1O_3⟩` = `tr I^{(l)} ⟨χχ⟩`, `tr I^{(l)} = 1`" [D for `l = 1`, now that the odd sign is fixed by §5.2; for `l = 2`
   the even half D, the odd half I]. "In the normalisation `C^{pres} := (8/M_pres)C` each MFT operator gives `2P/m^{2k}`" is **true by
   definition of that normalisation** (it is `¼ × 8 = 2`), not a result; reading the "2" as the helicity count is interpretation. For
   `⟨TTTT⟩` nothing here is tested: the relative weights and signs of the `±` terms, the two structures per parity, the identical-particle
   factor and positivity are all untouched (CHL give the MFT `[TT]` data in the helicity basis; not used) [ND].
5. **`O(1/m²)` (Goal 2) [D + NE].** MFT normalisation removes the "pol" term of `SPINNING_E2E.md` §6.3 identically [D]; it adds only
   `ν`-independent terms (`κ_l·pol`, `B_3 cosh β`), so it **cannot change any coefficient of `P′/P`**. The genuine zero-parameter result is
   the **difference test** `c_2(l=2) − c_2(l=1) + (κ_2−κ_1)·pol = −6P′/P` (`κ_2−κ_1 = 9/4`, pol derived), satisfied to `≤ 0.07` at `β = 0.4`
   and at `β = 1` separately (three `ν` each) [NE]: the `l`-dependence `−2(l_2²−l_1²)P′/P` is confirmed with no free parameter. The
   absolute per-`l` form `−2l²P′/P` at fixed `β` (`d_β = −7.996, −8.017` for `l = 2`) holds **relative to the scalar skeleton `(2,Δ_φ)`**
   (or the `l = 1` data), a single discrete, post-hoc choice of reference that is **not explained** [NE, reference choice unexplained];
   with the skeleton `(3,Δ_φ)` used in `SPINNING_E2E.md` §6.3 one gets the referee's `−8.34, −10.79` (point M1). So M1 is resolved only in
   the sense that the `l`-dependence is clean; the "extra `l = 2` piece" equals the difference of the two skeletons. `−2l²P′/P` itself is
   **not derived**.


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

### 3.1 Spectrum [NE (exact); agrees with the character count]
Exactly one conserved polynomial structure, hence one primary, for each `(N,J)` with `l ≤ J ≤ N + l` (all `N ≤ 8` for `l = 1`, `N ≤ 6` for
`l = 2`); none for `J < l` or `J > N + l` (scanned for `N ≤ 5`; in general this follows from the character of the conserved current,
`q^{l+1}(χ_l − qχ_{l−1})/(1−q)³`, i.e. level-`N` primaries `spin l ⊗ Sym^N ⊖ spin (l−1) ⊗ Sym^{N−1}`, as checked by the referee).
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

### 3.4 Literature check (novelty)
Checked (abstracts and, for the first two, the LaTeX sources): Sleight–Taronna arXiv:1804.09334 and 1807.05941; Karateev–Kravchuk–
Simmons-Duffin arXiv:1809.05111; Hofman–Li–Meltzer–Poland–Rejon-Barrera arXiv:1603.03771; Antunes–Costa–Gonçalves–Vilas Boas
arXiv:2111.05453; Albayrak–Meltzer–Poland arXiv:2006.07374.
* **Sleight–Taronna 1804.09334** give, in general `d` and for generic twists `τ_1` (spin-1) and `τ_2` (scalar), the MFT coefficients of the
  **leading-twist** `[JΦ]_{0,l}` (eq. `opejphi`), `[TΦ]_{0,l}` (eq. `2020mftope`) and `[O_JΦ]_{0,l}`, in their Mellin CPW basis. At `τ_1 = 1`,
  `d = 3` our `τ_0` family satisfies `p[S_1] = 4·a^{ST}_{(J,0)}` exactly for `J ≤ 9` and `Δ_φ = 7/4, 5/2, 13/10` (`st_compare.py`;
  `S_1 = 𝔏/L_1` is the pure-`H` structure on `τ_0`): a constant structure normalisation, i.e. **agreement**; the leading-twist even data
  are therefore known. Subleading twists (`n ≥ 1`) are not given there for `[JΦ]`.
* **KKSD 1809.05111** give a general MFT formula (pairings, shadow coefficients, Plancherel measure) and work out the 3d examples
  `⟨JJJJ⟩` and `⟨TTTT⟩` (identical currents / stress tensors), not `⟨JφJφ⟩`/`⟨TφTφ⟩`; their formalism would produce our coefficients.
* 1603.03771 (collider bounds; leading-twist `[Tφ]`-type data in the lightcone), 2111.05453 (higher-point lightcone) and 2006.07374 (3d
  fermion MFT): no overlap beyond leading twist.
Assessment [I]: likely new are the closed forms for all `n` for conserved `O` in `d = 3`, in particular the **parity-odd** ones, their
block-free derivation route, and the **exact cancellations** of §4 (MFT density against the spinning shadow coefficients and the
normalisation polynomials, the constant `¼`). The structural form "scalar MFT of twist dimension × helicity factor" is CHL's.


## 4. Heavy limit (`mft_heavy.py`, `out_mft_heavy.txt`)

### 4.1 Scalar identity and the mixed-scalar MFT normalisation [DV]
With `S := Δ_1+Δ_2`, **`P_FK(Δ_1,Δ_2;Δ,J)/b^{(Δ_12)}_{Δ,J} = (2/π) A_+A_− m² G(α)G(ᾱ)/(Γ(Δ_1)Γ(Δ_1−½)Γ(Δ_2)Γ(Δ_2−½))`**, `G(x) = Γ((x+S−2)/2)/Γ((x+4−S)/2)`,
exactly (the four Gamma functions with `±Δ_12/2` cancel between `P_FK` and `b`; check (1): `10^{−29}` at random points). Since
`G(x) = (x/2)^{S−3}(1 − (4/3)B_3((S−2)/2)/x² + O(x^{−4}))` (no `1/x` term because the two shifts add to 1), for the mixed scalar `[φ_1φ_2]`
with `q := 2n_J(2Δ−3)/(πm² P_FK)`:
`(−1)^JΠ̂_{k,ν}[block]/(q·2sin²)` per unit `P_FK`-normalised coefficient `→ (γ_{2Δ_2+k−1}γ_{2Δ_1+k−1}/a_{Δ_1,Δ_2}) · P_{−½+iν}(η)/m^{2k}`,
with the exact finite-`m` factor `(4/m²)^{S−3}G(α)G(ᾱ) = 1 − (8/3)B_3((S−2)/2) cosh β/m² + O(m^{−4})` (check (2): five `(Δ_1,Δ_2)`, `m` up to
4500; a check of the Gamma-function identity and its expansion, not a block check — the heavy scalar formula itself is DVE in
`SPINNING_E2E.md` §4.1). For `Δ_1 = Δ_2` this is CHMRS (2.59) (`a_{Δ_φ,Δ_φ} = a_{Δ_φ}`; their `p^MFT = 2P_FK` gives the `2`); the mixed case has `1·P`.

### 4.2 Asymptotics of `p^MFT` [D from §3.2]
`p^MFT_E[𝔖^cons] = b^{(1−Δ_φ)}_{Δ,J} · 2^{J+1} L_0² A_+A_−m² G(α)G(ᾱ) / (π c_l Π_l Γ(½)Γ(Δ_φ)Γ(Δ_φ−½))` (`S = 1+Δ_φ`, `Π_l` the product of §3.2),
`→ b^{(1−Δ_φ)} 2^{J+1} 4^{2−Δ_φ} L_0² m^{2Δ_φ−2}/(π^{3/2}c_l (A_+A_−)^{2l−1} Γ(Δ_φ)Γ(Δ_φ−½)) (1+O(m^{−2}))`;
`p^MFT_O[𝔒] = b^{(2−Δ_φ)} 2^{3J−2} H_1 A_+A_−m²G(α)G(ᾱ)/(π Γ(2)Γ(3/2)Γ(Δ_φ)Γ(Δ_φ−½))` (`S = 2+Δ_φ`, `l = 1`); `H_l → 1`; the `b`'s are exactly the shadow coefficients that appear in the heavy formulas (§4.3, §4.4), which is why everything cancels.

### 4.3 Parity even: MFT-normalised heavy formula [D; the Fock ↔ block bridge DV in §5.1]
Multiply the leading formula of `SPINNING_E2E.md` §4.5 (per `λ²`: `|T_l|/(2^J L_0²)`, prefactor `(π/4)4^{l+1+Δ_φ}/b^{(1−Δ_φ)}`,
`m^{−a_x−a_y}`, `a_y = 2Δ_O+k−1`) by `1/q^{(E)}`. `2^J`, `L_0²` and `b` cancel; with §4.1 the result is
`(γ_{a_x}γ_{a_y}/m^{2k}) · Y_l · [M_pres P + M_flip F^ν_{2l}] · (4/m²)^{Δ_φ−2}G(α)G(ᾱ)/a_{1,Δ_φ}`,
`Y_l = 4^l|T_l|/(m^{2l}c_l Π_l) = w_l (A_+A_−)^{2l}/Π_l`, `w_1 = ½`, `w_2 = 1/6` (using `𝖭_1 = A_+A_−m`, `𝖭_2 = ½(A_+A_−)²m²` of
`SPINNING_E2E.md` §4.3). Since `a_{Δ_O,Δ_φ}/a_{1,Δ_φ} = l!(½)_l` and **`w_l·l!(½)_l = ¼` for `l = 1, 2`**:

`C^{(q_1q_3)}_{k,ν}[Δ,J] = Q_E · ¼[M^{(q_1q_3)}_pres P_{−½+iν}(cosh β) + M^{(q_1q_3)}_flip F^ν_{2l}(β)] · |ĉ_E|²`, exactly at leading order, with
`Q_E = (4/m²)^{Δ_φ−2}G(α)G(ᾱ) · (A_+A_−)^{2l}/Π_l = 1 + [κ_l(1/sinh²(β/2)+1/cosh²(β/2)) − (8/3)B_3((Δ_φ−1)/2)cosh β]/m² + …`, `κ_1 = ¼`, `κ_2 = 5/2`.

Check (3) (an algebraic consistency check of the above simplification, not an independent test): `predict.prediction` (the code verified end-to-end in `SPINNING_E2E.md`) × `1/q^{(E)}` × `a/(γγ)` × `m^{2k}`, divided by
`¼[M_pres P + M_flip F]` and by the analytic `Q_E`: `1.000000000000000` in all 72 cases (`l = 1, 2`; `β = 0.4, 1, 2`; `m ≈ 20, 100, 1000`;
`k = 4, 6`; `ν = 0, 3`); `m²(Q_E−1)` converges to the stated coefficient (e.g. `l = 2`, `β = 0.4`: `72.12, 64.34, 63.995 → 63.991`).
Consequences: (i) the factor `(A_+A_−)^{2l}/Π_l` in `Q_E` is the inverse of the exact normalisation-polynomial ratio that
`SPINNING_E2E.md` found in the exact moments (§6.3 "pol", §7.2 fixed-`J` factor), so **the MFT-normalised density has neither** — it is
uniform down to `J = l` at leading order [D, given those findings]; (ii) the `ν`-, `k`-, `β`-dependence is unchanged.

### 4.4 Parity odd, `l = 1` [D; normalisation and sign DV at block level (§5.2); vertex DV pointwise]
* Differential basis [DV]: with `D̃ := ε(Z_1,P_1,∂_{P_1},P_2,∂_{P_2})` (spin `+1` at point 1, dimensions unchanged),
  `D̃ S3(Δ_1,Δ_2;Δ,J) = −J(Δ−1) ε(Z_1,P_1,Z_5,P_2,P_5) S3/v_5`, i.e. `D̃S3(2,Δ_φ) = −J(Δ−1)𝔒` (hand derivation; `odd_checks.py` (i) to 20
  digits at two points).
* Shadow coefficient [D]: applying `D̃^{(12)}D̃^{(34)}` to the scalar Lorentzian shadow representation of the seed `(2,Δ_φ,2,Δ_φ)`
  (`SPINNING_E2E.md` §4.1, DVE there) gives the odd partial wave `W^O := D̃D̃[K g]` with coefficient `1/b^{(2−Δ_φ)}` and right structure
  `D̃S3^S = −(J+1)(Δ−2)𝔒^S` (a single word, so no conservation-consistency condition is needed) [D; not tested numerically]. By the
  OPE-limit argument of `SPINNING_E2E.md` §1.4 (D-operators commute with the OPE limit of `Kg`) `W^O → 8^J (J(Δ−1))²⟨𝔒,𝔒⟩` at leading
  order, and the real-section contribution is `σ |λ_𝔒|² (−1)^J ⟨𝔒,𝔒⟩ = σ p[𝔒](−1)^J 8^{−J}W^O/(J(Δ−1))²` with **`σ = −1` measured** in §5.2.
* Heavy vertex [DV pointwise]: in the frame of `SPINNING_E2E.md` §4.3, `ε/v_5 → −2 z·e_⊥` at both vertices (`odd_checks.py` (iii):
  `−1.99644, −1.99996, −2.0000000` at `m = 10, 100, 10⁴`); `z·e_⊥ = (z·ε_+ − z·ε_−)/(i√2)`: odd, transverse, helicity `±1` (Lemma V).
  So the tensor is `4J(J+1)(Δ−1)(Δ−2)(z_1·e_⊥)(z_3·e_⊥) = −2J(J+1)(Δ−1)(Δ−2)Σ_h sgn(h_1)sgn(h_3)(z_1·ε_{h_1})(z_3·ε_{h_3})`.
* `M` with weights `sgn(h_1)sgn(h_3)`: `S^{q_1q_3}_{−h_1,−h_3} = −S^{q_1q_3}_{h_1h_3}` exactly (check (5), all modes, `l = 1, 2`), so the `S`-part
  cancels and the result is `−M_pres P + M_flip F` up to the overall factor.
* Product with `p^MFT_O[𝔒] = 2^{3J−3}P_FK(2,Δ_φ)H_1`: `[2J(J+1)(Δ−1)(Δ−2)/(8^J J²(Δ−1)²)] × 2^{3J−3}P_FK(2,Δ_φ)J(Δ−1)/((J+1)(Δ−2)) = ¼P_FK(2,Δ_φ)`:
  all `J, Δ` factors cancel exactly; with §4.1 (`S = 2+Δ_φ`, `a_{2,Δ_φ} = a_{Δ_O,Δ_φ}`): **`C^{(q_1q_3)} = Q_O·¼[M_pres P − M_flip F^ν_2]|ĉ_O|²`**,
  `Q_O = (4/m²)^{Δ_φ−1}G(α)G(ᾱ) = 1 − (8/3)B_3(Δ_φ/2) cosh β/m² + …` (check (4)); no normalisation polynomial at all.
* Sign [D, given `σ = −1` from §5.2]: on the real section an even exchange contributes `p(−½)^J W_diff` (verified in §5.1), whose
  heavy moment carries `(−1)^{J+l}` (`SPINNING_E2E.md` §4.5), i.e. sign `(−1)^l` per unit `p` in the preserved component. The odd exchange
  contributes `σ p[𝔒](−1)^J 8^{−J}W^O/κ²`, and the heavy tensor of `W^O` is `(−1)^J·κκ′·4(z_1·e_⊥)(z_3·e_⊥)` = `(−1)^J·(−2κκ′)Σ sgn(h_1)sgn(h_3)…`;
  for the preserved pairs (`sgn·sgn = −1`) this gives sign `σ·(+1) = −1 = (−1)^l` (`l = 1`): **the same sign as the even exchange**. Hence,
  in the sign convention of `C` (the positive quantity `(−1)^{J+l}Π̂` of `SPINNING_E2E.md`), the odd exchange gives `+¼[M_pres P − M_flip F^ν_2]`.
  This is the sign that reflection positivity requires; it is now an output, not an input. (`σ = −1` is what one expects from
  `λ_E² = −|λ|²` for a real Euclidean `ε`-structure — inversion is orientation reversing — combined with the Wick factor `i²` of the
  two `ε`'s [I; the measured value is what is used].)

### 4.5 What is not done for parity odd, `l = 2`
The MFT data and closed form exist (§3.2); the differential-basis odd structure (`D̃` composed with `D_11`, `D_12` on seeds `(4,Δ_φ)`,
`(3,Δ_φ+1)`), its conserved combination, shadow coefficient (expected `b^{(2−Δ_φ)}`) and heavy vertex were not computed. Expected by
analogy: the same `¼[M_pres P − M_flip F^ν_4]` [I].


## 5. Tests against blocks

### 5.1 Even parity: block decomposition of the MFT correlator, leading twist (`lightcone_check.py`, `out_lightcone_check.txt`) [DV]

On the real section of `SPINNING_E2E.md`, take `x_12 = T(1, (1−ε)n̂)` nearly null (`u → 0` at fixed `z̄ ≈ −1.33` (point A) or `−1.11`
(point B); all other `x_ij` timelike) and compare `F_MFT = ⟨O_1O_3⟩⟨φ_2φ_4⟩` (with `|x_ij²|` powers) with the sum over the `τ_0` family of
the exact e2e partial waves, `Σ_{J=l}^{J_max} p[𝔏](J−l, J) (−½)^J W_diff^{(τ_0+J, J)}` (`p` from the Fock closed form, exact on `N ≤ 8/6`;
on `τ_0`, `𝔏 = L_1S_1` resp. `L_2S_2` is finite while `L_0 = 0`; `J_max = 34/30`, last term `< 10⁻¹⁴`). In the component `z_1 ∥ n_{12}`
(the null direction of `x_12`), where `z_1·x_12 = 0` suppresses the structures of higher twist, the ratio is
`1 + c·ε` with **no free constant and no sign choice**:

| `l` | point | `z_3` | ratio−1 at `ε = 10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶` |
|---|---|---|---|
| 1 | A | `n_{12}` / `(1,cos2.6,sin2.6)` / `(1,−1,0)` | `1.7e−4, 1.7e−5, 1.7e−6, 1.7e−7` / `−2.0e−3 … −2.0e−6` / `−4.4e−2 … −4.6e−5` |
| 1 | B | same | `1.8e−4 … 1.8e−7` / `2.0e−5 … 2.0e−8` / `−3.1e−4 … −3.1e−7` |
| 2 | A | same | `1.7e−4, 1.7e−5, 1.7e−6` / `−6.3e−3 … −6.3e−5` / `−1.2e−1 … −1.4e−3` |
| 2 | B | same | `1.8e−4 … 1.8e−6` / `−2.7e−4 … −2.7e−6` / `−1.2e−3 … −1.2e−5` |

The control with `(+½)^J` instead of `(−½)^J` gives `−6.19, −5.05` (`l = 1`), `7.60, 6.12` (`l = 2`). This tests, for every `J` of the
leading-twist family at once, the full chain Fock-space `p` → `⟨Oφ X⟩ = λ𝔏` → "contribution `λ²(−½)^J W_diff`" of `SPINNING_E2E.md` §1.4, with
the actual blocks. (For generic `z_1` the ratio tends to an `O(1)` constant `≠ 1`: there the odd `τ_0+1` family and higher even twists
contribute at the same power of `u`, as expected; not a failure.) Not tested: families with `τ − τ_0 ≥ 2` (where `𝔖^cons` is regular)
other than through the closed form, which is exact on the Fock data.

### 5.2 Odd parity, `l = 1`: OPE test of normalisation and sign (`ope_odd_check.py`, `out_ope_odd_check.txt`) [DV]

`x_12 = δv → 0` on the real section. The primaries with `N ≤ 1` are `(0,1)` even, `(1,1)` odd, `(1,2)` even. Subtracting the two even
exchanges with their exact e2e blocks and Fock coefficients from `F_MFT` leaves `R = (odd (1,1) leading OPE term) + O(δ²)`, which is
compared with `p[𝔒](−1)^J⟨𝔒,𝔒⟩` (`p[𝔒] = 7/11`, the spin-`J` contraction of `t2_spinning_blocks.ope_check`, `ε` the Lorentzian determinant):

| polarisations `(θ_1, θ_3)` | `σ = R/(p(−1)^J⟨𝔒,𝔒⟩)` at `δ = 10⁻², 10⁻³, 10⁻⁴` | control `(F − (0,1) only)/F` |
|---|---|---|
| (0.7, 2.1) | `−1.0241, −1.0024, −1.00024` | `4.7·10⁻³, 4.7·10⁻⁴, 4.7·10⁻⁵` (order `δ`: the `(1,2)` term is needed) |
| (1.9, 0.4) | `−1.0355, −1.0035, −1.00035` | `1.3·10⁻³ …` |
| (3.0, 4.4) | `−0.8637, −0.9871, −0.99871` | `9.4·10⁻³ …` |

So `|σ| → 1` linearly in `δ`: the odd Fock coefficient, the odd structure `𝔒` and the contraction rule are consistent at block level, and
the real-section sign is `σ = −1` (used in §4.4). This also re-tests `p[𝔏](0,1)` and `p[𝔏](1,2)` (the subtraction must leave an `O(δ)`
remainder). Not tested: the odd shadow coefficient `b^{(2−Δ_φ)}` and the odd heavy formula in moments (no odd blocks in the Regge region).

### 5.3 Consequence for the exact moments of `SPINNING_E2E.md` (check (6)) — not an independent test

Since both sides of `SPINNING_E2E.md` §5 use the same `M`'s, the MFT-normalised ratio is `R^{MFT} = R_{e2e} × Q_E` exactly; for the
helicity components `R^{MFT}_{pres,flip} = R_{pres,flip} × Q_E`. Sample (`k = 6`, largest `m`; full table in `out_mft_heavy.txt`):

| `l` | `β` | `m` | `ν` | `R_pres` (e2e) | `R_pres^MFT` | `R_flip^MFT` | `m²(R_pres−1)` | `m²(R_pres^MFT−1)` |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.4 | 101.8 | 0 / 3 | 0.99072 / 0.99753 | 0.99133 / 0.99813 | 0.98535 / 0.98990 | −96.2 / −25.6 | −89.9 / −19.3 |
| 2 | 0.4 | 101.8 | 0 / 3 | 0.98529 / 0.99513 | 0.99140 / 1.00131 | 0.97446 / 0.97964 | −152.5 / −50.5 | −89.1 / +13.5 |
| 1 | 1.0 | 100.8 | 0 / 3 | 0.98760 / 0.97992 | 0.98770 / 0.98002 | 0.98193 / 0.99359 | −125.9 / −203.8 | −124.9 / −202.9 |
| 2 | 1.0 | 100.8 | 0 / 3 | 0.98669 / 0.97702 | 0.98777 / 0.97808 | 0.97581 / 0.98498 | −135.1 / −233.3 | −124.2 / −222.5 |

By construction this is `R_e2e × Q_E` (exact), so it adds nothing to the evidence beyond §5.1 and `SPINNING_E2E.md` §5; it shows
what the MFT-normalised ratio looks like: `1 + c_2/m²` for both helicity components, with the "pol" part removed. After normalisation the
preserved `c_2` of `l = 1` and `l = 2` nearly coincide at `ν = 0` (−89.9 vs −89.1; −124.9 vs −124.2), the difference at `ν = 3` being the
helicity term (§6).


## 6. `O(1/m²)` of the helicity-preserving component (Goal 2) (`goal2.py`, `out_goal2.txt`)

Input: the measured `c_2` of `SPINNING_E2E.md` §6.2 (`k = 6`, their fits (B)) and the scalar-skeleton `c_2` of the same table. MFT
normalisation adds the derived `q_2 = κ_l(1/sinh²(β/2)+1/cosh²(β/2)) − (8/3)B_3((Δ_φ−1)/2)cosh β` (§4.3), i.e. it removes "pol". **`q_2` is
`ν`-independent, so MFT normalisation cannot change any coefficient of `P′/P`**; what changes the per-`β` fits below relative to
`SPINNING_E2E.md` §6.3 is the reference.

| `β` | `ν` | `−6P′/P` | `c_2^{MFT}(l=2) − c_2^{MFT}(l=1)` | residual `c_2^{MFT}(l)+2l²P′/P − c_2^{MFT,skel(2,Δ_φ)}`, `l = 1 / 2` |
|---|---|---|---|---|
| 0.4 | 0 / 1 / 3 | 0.724 / 3.696 / 33.128 | 0.758 / 3.698 / 33.128 | −0.03, 0.003, 0.01 / 0.007, 0.005, 0.01 |
| 1.0 | 0 / 1 / 3 | 0.608 / 3.479 / −19.737 | 0.616 / 3.466 / −19.804 | −0.15, 0.04, −0.07 / −0.14, 0.03, −0.14 |
| 2.0 | 0 / 1 / 3 | 0.349 / 3.962 / −9.197 | 0.574 / 3.574 / −8.426 | −1.9, 1.3, −4.1 / −1.7, 0.9, −3.3 |

(`c_2^{MFT,skel}` = skeleton `c_2 − (8/3)B_3(Δ_φ/2)cosh β`; quoted fit errors of the inputs: `≲ 0.3` for `β ≤ 1`, a few units at `β = 2`.)

Per-`β` fits (the "direct test" of `SPINNING_E2E.md` §6.3, but with the `l = 1` data as reference — a post-hoc choice of reference): fit
`c_2^{MFT}(l=2) − [c_2^{MFT}(l=1) + 2P′/P] = c_β + d_β P′/P` over `ν = 0, 1, 3`:

| `β` | `d_β` (prediction `−8`) | `c_β` | max residual with `d = −8` |
|---|---|---|---|
| 0.4 | −7.996 | 0.02 | 0.022 |
| 1.0 | −8.017 | −0.01 | 0.043 |
| 2.0 | −7.51 | 0.07 | 0.59 |

Findings: (a) **zero-parameter difference test** [NE]: `c_2(l=2) − c_2(l=1) + (9/4)·pol = −6P′/P` (equivalently the third column
minus the second) to `≤ 0.07` at `β = 0.4` and at `β = 1` separately; the referee reproduced it from the full-precision fits
(residuals `≤ 0.06`). This confirms the `l`-dependence `−2(l_2² − l_1²)P′/P`, with nothing chosen. (b) **Reference-dependent statement**
[NE, reference choice unexplained]: relative to the scalar skeleton `(2,Δ_φ)` (or the `l = 1` data), the absolute per-`l` term is
`−2l²P′/P` at fixed `β` for both `l` (`d_β = −7.996, −8.017` for `l = 2`), and `c_2^{MFT}(l) = c_2^{MFT}[skeleton (2,Δ_φ)] − 2l²P′/P` within the
input errors. Relative to the skeleton `(3,Δ_φ) = (Δ_O, Δ_φ)` of `SPINNING_E2E.md` §6.3 the per-`β` coefficients are `−8.34, −10.79` (M1).
Choosing `(2,Δ_φ)` is a single discrete, post-hoc choice; it is *not* derived. (Natural candidates from this paper would be the
effective scalar `(1,Δ_φ)` of the MFT data and of `b^{(1−Δ_φ)}`, or `(Δ_O,Δ_φ)` of `a_{Δ_O,Δ_φ}`; neither is `(2,Δ_φ)` for `l = 2`, and the
skeleton `c_2` for `(1,Δ_φ)` was not computed.) So M1 is resolved only for the `l`-dependence; the fitted "extra `l = 2` piece"
`−2sinh²βP′/P − 17.3 cosh β` equals `c_2^{skel(2,Δ_φ)} − c_2^{skel(3,Δ_φ)}` (+ the `l = 1` constants `0.07 + 0.10 cosh β`), e.g. `−18.5` vs `−18.49`
at `β = 0.4`, `ν = 0`. The "pol" term is exactly the MFT normalisation [D].
**Not derived**: the helicity term `−2l²P′/P` and the skeleton `c_2` itself (the analogue of CHMRS (3.27) for unequal dimensions and
tensor kernels); the flip `c_2` was not re-analysed. The mechanical route is still item 2 of `SPINNING_E2E.md` §11.


## 7. Answer to the question, stated precisely

For `⟨OφOφ⟩` (one spinning Regge pair, `O = J, T`) and the tensor Regge moments `Π̂^{(q_1q_3)}_{k,ν}`, with heavy densities defined as in
CHMRS (2.31)–(2.32) but with `p^MFT` of the exchanged operator's parity sector, and the functional normalised by `a_{Δ_O,Δ_φ}/(γ_{2Δ_φ+k−1}γ_{2Δ_O+k−1})`:

`C^{(q_1q_3)}_{k,ν}[Δ,J] → (1/(4m^{2k})) [ M^{(q_1q_3)}_pres(ν) P_{−½+iν}(η_AdS) (|ĉ_E|² + |ĉ_O|²) + M^{(q_1q_3)}_flip(ν) F^ν_{2l}(β) (|ĉ_E|² − |ĉ_O|²) ]`,

each term from operators of the indicated parity (`E` at `τ−τ_0` even, `O` odd). In helicity language the weight of the pair
`(h_1,h_3)` is `¼c^{h_1}c̄^{−h_3}`, so "`|c^{++}|² ± |c^{+−}|²`" is replaced, for one pair, by `|c^+|²+|c^−|²` (preserved) and `2Re(c^+c̄^−)` (flip),
each with the exact constant `1/8` times the summed transform factor. Constants: `¼` exact (`w_l l!(½)_l = ¼`, `l = 1, 2`); `a_{Δ_O,Δ_φ}`
(not `a_{1,Δ_φ}` and not `a_{Δ_φ}`). "`2P/m^{2k}`" per MFT operator in the normalisation `C^{pres} = (8/M_pres)C` holds by definition of
that normalisation. Spin-traced and summed over parities, the density equals that of a scalar of dimension `Δ_O`. Status: `E`, `l = 1, 2`
[D from NE-exact closed forms; Fock ↔ block bridge DV (§5.1)]; `O`, `l = 1` [D; normalisation and sign DV at block level (§5.2); shadow
coefficient not tested]; `O`, `l = 2` [ND]; `⟨TTTT⟩` [ND: no tested content here].


## 8. Status list

| # | statement | mark |
|---|---|---|
| 1 | Fock-space method for exact MFT OPE coefficients of `[Oφ]` with `O` conserved | [D]; validated on scalars (FK exact) and normalisations |
| 2 | Spectrum: one primary per `(N,J)`, `l ≤ J ≤ N+l`, parity `(−1)^{N−J−l}` | [NE, exact, `N ≤ 8/6`]; character count (referee) |
| 3 | Closed forms of `p^MFT` (even `l = 1,2`; odd `l = 1,2`) | [NE, exact on 145/145]; not derived (CHL analogue [Q]); leading twist agrees with Sleight–Taronna [NE, exact] |
| 4 | `P_FK/b` Gamma identity; mixed-scalar MFT normalisation `1/a_{Δ_1,Δ_2}`, exact `1/m²` factor | [DV] (identity checks, not block checks) |
| 5 | MFT-normalised even formula `¼[M_pres P + M_flip F]`, `Q_E` exact | [D from 3]; Fock ↔ block normalisation DV on the `τ_0` family (§5.1); check (3) and §5.3 algebraic |
| 6 | MFT normalisation cancels the normalisation polynomials ("pol", fixed-`J` factor) | [D] |
| 7 | `D̃S3 = −J(Δ−1)εS3/v_5`; odd vertex `→ −2z·e_⊥`; `S`-cancellation | [DV] |
| 8 | Odd `l = 1` formula `¼[M_pres P − M_flip F]`, `Q_O` | [D]; normalisation and sign `σ = −1` DV (OPE test §5.2); shadow coefficient [D, not tested] |
| 9 | Odd `l = 2` heavy formula | [ND] |
| 10 | Spin-traced, parity-summed MFT density = scalar `(Δ_O,Δ_φ)` | [D] `l = 1`; [I] `l = 2` (odd half) |
| 11 | `c_2(2) − c_2(1) + (9/4)pol = −6P′/P` (zero-parameter); per-`l` `−2l²P′/P` relative to skeleton `(2,Δ_φ)`; extra `l = 2` piece = skeleton difference | [NE]; reference choice post hoc, unexplained |
| 12 | Analytic `−2l²P′/P` and skeleton `c_2` | [ND] |
| 13 | `⟨TTTT⟩` normalised form | [ND] (no tested content) |

**What remains.** (i) A derivation of the closed forms (KKSD/CHL pairing formula with one generic scalar leg, or a closed evaluation of
the Fock norm); (ii) the odd `l = 2` heavy vertex and shadow coefficient (mechanical with `D̃` and the `e2e` word machinery), and a direct
numerical test of the odd shadow coefficient (adapt `t5_shadow_direct.py`); (iii) parity-odd blocks, to test §4.4 in moments and to
extend §5.1 beyond leading twist; (iv) Goal 2 analytically (`SPINNING_E2E.md` §11 item 2), including why the skeleton `(2,Δ_φ)` is the
right reference; (v) two spinning pairs (`⟨TTTT⟩`: two structures per parity, `±` weights, identical-particle factor).

## 9. Code index (`ads4/code/mft/`; `sh run_all.sh`)

| file | content | output |
|---|---|---|
| `fock.py` | Fock-space kernels, exact norms, structures at infinity, conservation | — |
| `mft_exact.py` | exact `p^MFT` for `l = 1, 2`, three `Δ_φ`, vs closed forms (cache in `cache/`, git-ignored) | `out_mft_exact.txt` |
| `odd_checks.py` | `D̃` identity; odd structure at infinity and in the heavy frame | `out_odd_checks.txt` |
| `mft_heavy.py` | checks (1)–(6) of §4–5 | `out_mft_heavy.txt` |
| `goal2.py` | §6 | `out_goal2.txt` |
| `lightcone_check.py` | §5.1: leading-twist block decomposition of the MFT correlator | `out_lightcone_check.txt` |
| `ope_odd_check.py` | §5.2: odd normalisation and sign by the OPE expansion | `out_ope_odd_check.txt` |
| `st_compare.py` | §3.4: comparison with Sleight–Taronna | `out_st_compare.txt` |

## 10. Revisions after `REFEREE_SPINNING_MFT.md`

* Issue 1: block-level test of the Fock ↔ block normalisation added (§5.1, leading-twist family, `l = 1, 2`); check (3) and the old §5 relabelled as
  algebraic consequences; "end-to-end DV" removed.
* Issue 2: the odd sign is no longer fixed by positivity: measured by an OPE-level block test (§5.2, `σ = −1`, which also tests the odd normalisation)
  and propagated to the heavy formula (§4.4); the spin-traced statement for `l = 1` now rests on it.
* Issue 3: §0.5 and §6 rewritten — MFT normalisation cannot change `d_β`; the zero-parameter content is the difference test; the per-`l` form relies on
  the post-hoc, unexplained skeleton `(2,Δ_φ)`.
* Issue 4: "`2P/m^{2k}` per operator" stated as true by definition; `⟨TTTT⟩` extension marked [ND].
* Novelty: §3.4. Minor: run time, spectrum range, check (2) label.

# The light side of the AdS4 functional `Ω_{g3}` at leading order in `1/c_T`

Status: written 16 Nov 2026 by a derivation/numerics agent; **revised after `REFEREE_ADS4_LIGHT.md`** (M1: input (E) separated and
its leverage quantified; M2: positivity-threshold dependence of the bound; minors m1–m6). **Not a referee-approved result. Not a theorem; no AdS4 bound on `ĝ₃` follows from it
without the heavy-side positivity of `ADS4_FUNCTIONAL.md`, which is itself open.** Code: `ads4/code/adslight/` (`bash run_all.sh`;
logs in `adslight/cache/`, git-ignored). Imports nothing from, and modifies nothing in, other code directories.

Marks: **[D]** derived here (elementary or with the inputs named in the line); **[DV]** derived and checked numerically; **[Q]**
quoted (CHMRS = arXiv:2106.10274, LaTeX source read; MP = Meltzer–Perlmutter arXiv:1712.04861, LaTeX source read; CDP =
Chowdhury–David–Prakash 1707.03007 via `DICTIONARY_T4.md` §1); **[I]** inferred, not computed; **[A]** assumed.
Notation as in `ADS4_FUNCTIONAL.md`: AdS units `R = 1`, `d = 3`, `Λ = √((Δ_g−2)(Δ_g−1))`, `ν = Λp`, weights `ψ₂, ψ₃`, `μ₂_j = ∫₀¹ψ₂p^j`,
`μ₃_j = ∫₀¹ψ₃p^j`, `x² := |ĝ₃|²Λ⁸` (ĝ₃ in AdS units), `t₄² + α_T² = 57600|ĝ₃|²` (`DICTIONARY_T4.md`), and

`Ω_{g3} = ∫₀¹dp[ψ₂Λ⁴C⁺_{2,Λp} + ½ψ₃Λ⁶C⁻_{3,Λp}] − 5Λ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (5/2)Λ¹⁰C⁻_{5,0}`,  sum rule `Ω|_heavy = −Ω|_light`.

## 0. Results

1. **Scalar warm-up [DV].** CHMRS's `−C_{2,ν}|_grav = 8πG/(ν² + 9/4)` (d = 3, normalisation in which a heavy state gives `2P/m⁴`) is
   reproduced numerically from the stress-tensor single-trace Mellin poles alone (primary + descendants, CHMRS (MGravityRegge)),
   relative accuracy `≤ 5·10⁻⁹` at `Δφ = 1.2, 1.75, 2.6`, `ν = 0 … 12`; the primary pole alone is *not* of this form (ratio 0.69–0.99).
   The implementation is calibrated on two contact terms with known exact answers (`−C₂ = 2`, `ν² − (2Δφ+½)(2Δφ−½)`; agreement `10⁻¹⁴`).
   Dictionary: **`8πG = 256/c_T`** (`G = 32/(πc_T)`), identical to `DICTIONARY_T4.md`'s `c_T = 32L²/(πG)` [Q + D].
2. **Pole [D, DV].** The transverse space is `H₂` (= `H_{d−1}`), not `H₃`. A spin-`J` exchange of dimension `Δ` gives the `H₂` propagator
   of a field of dimension `Δ − 1`, i.e. the pole `1/(ν² + (Δ − 3/2)²)`; for `T` (`Δ = 3`): **`ν² + 9/4`**. `ν² + 1/4` would be a `Δ = 2`
   exchange. This is a statement about the exchanged operator only, so it holds for every external spin and helicity component.
3. **Graviton term of `⟨TTTT⟩` [D given (F), (R)].** In the preserved–preserved moments (the only ones in `Ω_{g3}`):
   **`−C⁺_{2,ν}|_{T-exch} = 16πG/(ν² + 9/4) + 16πG c₁ t₄ + G t₄²·(polynomial of degree ≤ 2 in ν², feeding b₂, b₁, b₀)`**.
   **That there is no graviton-exchange term at all in `C⁻₃`, `C⁺₄`, `C⁻₅` is [I]: input (E), §2.1**. (F) alone excludes only polynomial
   terms there, not decaying residues such as `G/(ν²+9/4)³`, and such residues would dominate `Ω` by factors `10⁴` (`C⁻₃`) to `10¹³` (forward
   terms) at `Δ_g = 100` (`t4`).
   The residue at the pole is independent of `t₄` and `α_T`: it is the diagonal, in the helicity basis, of the Hofman–Maldacena/CDP flux
   matrix, `1 + t₄(¼ − ¼) = 1` [D from Q]. The Einstein part has **no** polynomial remainder (forced by the flat limit, §2.2). The only
   allowed linear term is a `ν`-independent `16πG c₁ t₄` (an `L = 0` contact in impact parameter; `c₁` not computed); it moves the bound by
   relative `O(x⁻¹Λ⁻²)`. Parity forbids a linear `α_T` term. **The conjecture `16πG/(ν² + 9/4)` of `ADS4_FUNCTIONAL.md` is therefore right,
   up to these two inputs.**
4. **`|ĝ₃|²` terms [D given (F), (P), (E)].** They come from the `t`-channel diagram with two `W³` vertices, which in the `u`-channel Regge
   limit acts like a contact: `−C⁺_{2,ν} ⊃ 2πG|ĝ₃|²(ν⁶ + b₂ν⁴ + b₁ν² + b₀)`, `−C⁻_{3,ν} ⊃ −4πG|ĝ₃|²(ν⁴ + f₁ν² + f₀)`,
   `−∂_{ν²}C⁺_{4,ν}|₀ ⊃ 4πG|ĝ₃|²` **exactly**, `−C⁻_{5,0} ⊃ 0` **exactly**. Hence
   `k_AdS(Λ) = k(0) + (−9677b₂/13860 + 225f₁/14)/Λ² + (−20b₁/9 + 7f₀/3)/Λ⁴ − 2591b₀/(210Λ⁶)`, `k(0) = 330151/12012` (`t2`);
   `b_i`, `f_i` are `O(1)` numbers that are not computed. "Exactly" means: no non-polynomial residue, which is (E) again.
5. **Contacts [D given (F), (P), (E)].** `g₄` drops out exactly; **`g₅` cancels exactly between `∂_{ν²}C⁺₄` and `C⁻₅`**, so the `g₅`-free
   completion survives the uplift with no AdS residue; `g₆, g₆′` and every higher contact leave residues suppressed by `Λ⁻²` relative to
   their natural size. They are `O(G)` in `Ω` against the leading `O(GΛ² log Λ)`, if the contacts obey EFT power counting [A].
6. **Double traces [Q for scalars, I for `⟨TTTT⟩`].** The claim that they "vanish at leading order by double zeros" is **only partly true**.
   Families `[TT]_{n,ℓ}` with `n ≥ k/2` sit at double zeros of the `C_k` functional and enter at `O(γ²) = O(1/c_T²)`. The first `k/2`
   trajectories (and the leading-trajectory OPE corrections) enter at `O(1/c_T)`. They are not an independent unknown, because their
   `O(1/c_T)` data is fixed by the same bulk diagrams. CHMRS's contour formula already includes them: the `1/γ²` normalisation of `C_{k,ν}`
   cancels their poles, which is how the clean `1/(ν²+9/4)` arises (item 1). The first effect not captured by the tree-level light side is
   `O(1/c_T²)`, i.e. one loop, of relative size `~ GΛ² ≈ 10Λ²/c_T`.
7. **Conditional AdS4 statement.** If `Ω|_heavy ≥ 0` **for all heavy states down to `m ≈ Λ`** and the inputs of §4 hold:
   **`|ĝ₃|²Λ⁸ ≤ A log Λ + B_AdS + C₁/Λ + O(log Λ/Λ²)`**, with `A = 1040/k(0) = 37.8387` [D given (F)], `B_AdS = −60.7875` and
   `C₁ = 5765760π/330151 = 54.865` [both D given (F), (R)]. Equivalently, **`t₄² + α_T² ≤ 57600(A log Λ + B_AdS + C₁/Λ)/Λ⁸`**,
   e.g. `√(t₄²+α_T²) ≤ 2.7·10⁻⁵` at `Δ_g = 100` (`t3`). **Caveat (M2):** positivity down to `m ≈ Λ` is open and is *contradicted at
   leading Regge order* near `m ≈ 6.3Λ` (`ADS4_FUNCTIONAL.md` §3(a); rescued only in the [M] hybrid model). With positivity only for
   `m ≥ nΛ_eff`, one must use `Λ_eff = Λ/n`, and the bound weakens by `≈ n⁴`: with the Regge-side threshold `n = 7.4`, it is `4.9·10⁻²` at
   `Δ_g = 100` and `8.1·10⁻⁶` at `Δ_g = 1000` (§4, `t4`). `B_AdS = −45.4452 − 1040 log(3/2)/k(0)`: the flat logarithm is cut off at
   `ν = 3/2`, i.e. `m_IR → 3/(2Λ)` in flat units [D]. No eikonal or IR hypothesis is needed on the light side.

## 1. Scalar warm-up (`t1_scalar_mellin.py`)

**1.1 Formula [Q + D].** CHMRS (CfromBMellin) with (BMellinDefinition) for `k = 2`: if `M(S,U) = A(U)S² + O(S)` at fixed `U`,
the residue at `S′ = ∞` picks only `A`, and `B̂_{2,T} = −4A(2Δφ − T)` [D; the orientation is fixed by the `M_{1,0}` check]. With
`T = 2z`, `c = Δφ + ¼ − iν/2`:

`−C_{2,ν}[A] = (π^{3/2}a_{Δφ}/2)|Γ(c)|⁻⁴ ∫_{Re z = z₀} dz/(2πi) Γ(z)²Γ(c−z)Γ(c̄−z) A(2Δφ − 2z)`.

For a constant `A` the Barnes lemma gives a `ν`-independent result, and with CHMRS's `c₂` it is exactly `2`. This is `−C₂[M_{1,0}] = 2` [D].

**1.2 Checks [DV, `cache/t1.log`].**
* (a) `M_{1,0}`: `−C₂ = 2.0` and (b) `M_{0,1} = c₃STU`: `−C₂ = ν² − (2Δφ+½)(2Δφ−½)` (CHMRS (C2C4Light) at d = 3), both to `≤ 6·10⁻¹⁵`.
  The tested values are `Δφ ∈ {1.2, 1.75, 2.6}` and `ν ∈ {0, 0.7, 3, 9}`.
* (c) Graviton: `A(U) = Σ_m R_m/(U − 1 − 2m)`, where `R_m` is CHMRS (MGravityRegge). These are the `T` primary (`m = 0`) and its descendants: the
  `T` block in Mellin space. The contour runs between the `T` poles and the `Γ(c−z)` poles. The terms `m < 12` are integrated directly; the tail is summed
  to `m = 2·10⁴` via the moments `∫Γ²ΓΓ z^j`. Result: **`Σ_m = 8πG/(ν²+9/4)` to `≤ 4·10⁻¹⁴`** for `ν ≤ 4`, and to `≤ 6·10⁻⁹` at
  `ν = 12`, where tail truncation limits the accuracy. `(MGravityRegge)` therefore includes the `8πG` once. The `m = 0` term alone gives
  0.69–0.99 of the target, so the exact single pole is a property of the full `T` contribution divided by `γ_{2Δφ+1}(ν)²`.
* (d) `H₂` harmonic analysis: `∫₀^∞ dν/(2π) πν tanh(πν) P_{−½−iν}(g)/(ν² + a²) = h_{D₀}(g − 1)`, with `D₀ = ½ + a` (CHMRS's `h`). Checked
  for `a² = 9/4 ↔ D₀ = 2` and `a² = 1/4 ↔ D₀ = 1` to `≤ 2·10⁻³` at `g ≤ 3` and `1.2·10⁻²` at `g = 10` (oscillatory tail cut at `ν = 150`). The constant: we find ratio 1, while our
  reading of CHMRS's prefactor `vol S^{d−2}/2^{d−2}` gives `π`. Unresolved: it may be our reading of their conventions, and it is not claimed to be an
  error in CHMRS. It is immaterial, because only the pole ↔ dimension identification is used.
* (e) `8πG·dΓ(d/2)³/(4π^{d/2}Γ(d+2)) = 1/c_T` [Q CHMRS] ⇒ `8πG = 256/c_T` and `f²_{φφT} = 3Δφ²/(8c_T)`. This agrees with the standard
  `dΔφ²/(4(d−1))·c_T^free/c_T`.

**1.3 What this fixes.** In the normalisation in which a heavy state contributes `2P_{−½+iν}(cosh β)/m⁴` with CHMRS's positive measure, the
light side is `(256/c_T)/(ν² + 9/4)` [DV]. For `⟨TTTT⟩`, the corresponding absolute constant needs the `l = 2` heavy normalisation (`SPINNING_TTTT.md` item 5, open)
[I]. **It does not enter the bound**: the bound uses only the ratio of the graviton and `|ĝ₃|²` terms, both computed in the same correlator,
and the heavy side enters only through its sign.

## 2. `⟨TTTT⟩`: light side of the pres–pres moments

**2.1 Structure [Q + I].** CHMRS (eftformula)/(moregeneralformula): for a tree-level EFT correlator, the light side of a pure `Π_k`
moment is a Laurent coefficient of the Regge expansion. Only terms of Regge spin `n ≥ k` contribute, and the light side is not the
(divergent) sum of moment actions over light operators. We apply this to the tensor moments of `SPINNING_E2E.md` §3.1 [I: the proof
in CHMRS App. (omega light) uses only the contour structure of dDisc-based functionals, so it should carry over]. Regge spins: the
graviton exchange in any channel, with any vertices, has spin 2. So **only `C⁺₂` sees the `u`-channel pole**. `C⁻₃`, `C⁺₄`, `C⁻₅`
receive only contact-like (polynomial) terms [D given this].

**Input (E) [I]: the spinning carry-over of CHMRS (eftformula).** The light side of a pure tensor `Π_k` moment is *exactly* the Regge
Laurent coefficient, so a term of Regge spin `n < k` gives exactly zero, not merely a decaying function of `ν`. (P), the polynomial form of
the contact-like terms, is the same statement for `n ≥ k`. Evidence:
* (i) Scalars [D given Q]: in Mellin space `B̂_{k,T}` is a residue at `S′ = ∞` (CHMRS (BMellinDefinition)). An exchange amplitude bounded by `S²`
  at fixed `U` gives an integrand `O(S′^{−3})` for `k ≥ 4`, hence `B̂_{k≥4} ≡ 0` and `C_{k≥4,ν}|_grav ≡ 0` for all `ν`.
* (ii) Flat space [D, `t2`]: the per-channel Einstein pieces of `−B₃` are `±8πG/p⁴`, which are non-zero. Their sum vanishes exactly because
  `f ~ s⁻²` makes the residue at infinity zero. This is the same mechanism.

For `⟨TTTT⟩` neither a spinning Mellin representation of the pres–pres component nor the contour argument for the tensor kernels was
written down here. (F) does **not** imply (E): (F) fixes only the fixed-`p` scaling and allows `O(G)` residues decaying like `ν⁻⁶` (`C⁻₃`) or
`ν⁻⁸` (`C⁺₄`), and it says nothing at `ν = 0`. **Leverage** (`t4`): consider a test residue with the same coefficient `16πG` as the pole. The table gives its contribution to
`−Ω` divided by the leading `16πGΛ²L_AdS`:

| test residue | `Δ_g = 20` | `100` | `1000` |
|---|---|---|---|
| `C⁺₂`: `r = (ν²+9/4)⁻²` (this is (R), not (E)) | 0.14 | 0.073 | 0.042 |
| `C⁻₃`: `(ν²+9/4)⁻³` (`ψ₃(0) = 162`, scales `~Λ³`) | 1.4·10² | 2.6·10⁴ | 1.8·10⁷ |
| `C⁺₄` forward: `(ν²+9/4)⁻⁴` (`~Λ⁸`) | 5·10⁷ | 1.6·10¹³ | 10²¹ |
| `C⁻₅` forward: `(ν²+9/4)⁻⁵` (`~Λ⁸`) | 6·10⁶ | 2·10¹² | 10²⁰ |

So only an exact zero is usable: a `C⁻₃` residue must be below `~10⁻³Λ⁻³`, and a forward residue below `~Λ⁻⁸`, relative to the pole coefficient.
The bound therefore stands or falls with (E).

**2.2 Flat-limit lemma [D given (F)].** Hypothesis (F) [Q for scalars: CHMRS (BknuFlatSpaceLimit) and §(light position space); I for
`⟨TTTT⟩`]: at fixed `p = ν/Λ`, with flat-unit couplings `G_f = GΛ²`, `ĝ_f = ĝ₃Λ⁴`, `g_f` fixed, `ε_kΛ^{2k}(−C^±_{k,Λp})` tends to the
flat `−B_k(p)` as `Λ → ∞`. A term `G^aĝ₃^b ν^{2j}` in `−C^±_k` scales as `Λ^{E}` with `E = 2k − 2a − 4b + 2j` in flat units. It must have
`E < 0` unless it is the flat term (`E = 0`). Consequences (`a = 1`):

| moment | Einstein (`b = 0`) | linear `ĝ₃`/`t₄` (`b = 1`) | `|ĝ₃|²` (`b = 2`) |
|---|---|---|---|
| `C⁺₂` | no polynomial at all (only the pole) | constant only (`j = 0`), suppressed `Λ⁻²` | degree ≤ 3, top = flat `2πG ν⁶` |
| `C⁻₃` | no polynomial ⁽*⁾ | no polynomial ⁽*⁾ | degree ≤ 2, top = flat `−4πG ν⁴` |
| `C⁺₄` | no polynomial ⁽*⁾ | no polynomial ⁽*⁾ | degree ≤ 1, top = flat `4πG ν²` ⇒ `∂_{ν²}` exact ⁽*⁾ |
| `C⁻₅` | no polynomial ⁽*⁾ | no polynomial ⁽*⁾ | degree 0 would be a flat term, absent ⇒ **0** ⁽*⁾ |

⁽*⁾ Non-polynomial decaying residues are not excluded by (F). They are excluded only by (E), §2.1.

For contacts with flat term `g_f p^{2j_f}`, the AdS polynomial has degree `≤ j_f` with the top coefficient equal to the flat one. A contact absent from
flat `B_k` contributes only with degree `< j_f`, i.e. suppressed by `Λ⁻²` [D given (F) and (P): the light side of contact-like terms is
a polynomial in `ν²`. This is Q for scalar contacts (CHMRS (BknuContactsLeading)) and I for spinning].

**2.3 The `u`-channel `T` exchange.**
* Pole at `ν² = −9/4` [D, §0 item 2; DV for the pole ↔ dimension map].
* The leading Regge term of a single exchange is (polynomial vertex factors in `ν`) × `1/(ν² + 9/4)` × (ratio `R(ν)` of double-twist
  `Γ`-factors to the `γ²` divided out in `C^±`). For scalars `R ≡ 1` [DV, §1]. **(R) [I]:** `R ≡ 1` also for the spinning `C^±`
  normalisation of `SPINNING_4PT.md` (`γ_a²` factors).
* Residue [D from Q]: at the pole the pres–pres vertex is MP (TTdiagonal) evaluated at `ν₀ = −3i/2`, i.e. the HM structure
  `1 + t₄(|n·ε·n|²/|ε|² − ¼)`. For circular polarisation `ε = e_±e_±` one has `|n·ε·n|² = ¼|ε|²`, so the `t₄` term vanishes. The parity-odd term
  enters the CDP matrix `[[1−t₄/4, α_T/4],[α_T/4, 1+t₄/4]]` only off-diagonally in the helicity basis. **Diagonal = 1**, which equals the Ward-identity
  normalisation fixed by `c_T`. Helicity flip: `|⟨∓|E|±⟩| = |c| = √(t₄²+α_T²)/4` with a spin-4 harmonic on `H₂`. This enters only
  flip moments (`𝒦_{±4,·}`), which are not in `Ω_{g3}`.
* Einstein part: pole + polynomial, residue `16πG` (flat matching `16πG/p²`, `FLAT_GRAVITON_BOUND.md` (2.31a)), polynomial ≡ 0 by the
  lemma ⇒ **`16πG/(ν² + 9/4)` exactly** [D given (F), (R)]. Independently [Q]: CEMZ's AdS shock analysis (`DICTIONARY_T4.md` §3.1) has
  Einstein phase shifts helicity-independent at every impact parameter.
* `t₄`-dependence: `V(ν)/(ν²+9/4)` with `V(−3i/2)` `t₄`-free ⇒ `t₄[V(ν)−V(−3i/2)]/(ν²+9/4)` is polynomial, and by the lemma it is the
  constant `16πG c₁t₄`. Its contribution `−Ω ⊃ (10364/105)πG c₁t₄Λ⁴ = −23690πG c₁ x cos φ` (`t₄ = −240 Re ĝ₃`, `t2`) is `O(G x)` against
  `O(GΛ²x²)`. It shifts `x` by `≈ 431|c₁|/Λ²` (`= 23690/(2k(0))`). At `Δ_g = 100` this shifts the bound on `x²` by `0.95|c₁|`, more than the
  retained `C₁/Λ = 0.557`. So quoting `C₁/Λ` as the next term is meaningful only for `|c₁| ≲ 0.5`, or asymptotically: at `Δ_g = 1000`
  the shift is `0.012|c₁|` against `0.055` (`t4`).
* `C⁺_{2,0} = −64πG/9` is finite. AdS regulates the forward limit (CHMRS footnote, sec. holographic sum rules) [Q].

**2.4 `|ĝ₃|²` terms (`t`-channel `W³⊗W³`) and the forward terms.** Flat (`t2`, exact): `−B₂ ⊃ 2πG|ĝ₃|²p⁶`, `−B₃ ⊃ −2πG|ĝ₃|²p⁴`,
`−B₄ ⊃ 4πG|ĝ₃|²p²`. All of these come from `t`-channel poles ("`++`" channel, heavy analogue `|c⁺⁺|²`, `FLAT_GRAVITON_BOUND.md` §1.2). The Einstein graviton
splits as `8πG/p² (+−) + 8πG/p² (++)` in `−B₂` and `±8πG/p⁴` (cancelling) in `−B₃`. With `B_k ↔ ε_kΛ^{2k}C^{(−)^k}_k` and the lemma: item 4 of §0.
Coefficient at leading order [D]: `|ĝ₃|² = |α₃ + iα̃₃|² = (t₄² + α_T²)/57600`. AdS corrections may treat `t₄²` and `α_T²` differently only
inside the unknown `b_i`, `f_i` [I].

**2.5 Contacts [D given (F), (P)].** Flat `−B₄ = 2g₄ + g₅p² + (g₆+g₆′)p⁴ + …`, `−B₅ = g₅ + (g₆−g₆′)p² + …` (`t2`). AdS:
`−C⁺₄ ⊃ Λ⁻⁸[2g₄] + Λ⁻¹⁰g₅(ν² + e₅) + (g₆+g₆′)Λ⁻¹²(ν⁴ + a₁ν² + a₀)` and `−C⁻₅ ⊃ 2Λ⁻¹⁰g₅ + 2Λ⁻¹²(g₆−g₆′)(ν² + a₀′)` (flat-unit couplings).
So (ĝ₃ in AdS units, contact couplings in flat units) `−5Λ¹⁰∂_{ν²}(−C⁺₄)|₀ + (5/2)Λ¹⁰(−C⁻₅)|₀ = −20πG|ĝ₃|²Λ¹⁰ + 5[a₀′(g₆−g₆′) − a₁(g₆+g₆′)]/Λ²`. In it, `g₄` is gone, `g₅` cancels exactly,
`g₆, g₆′` leave `Λ⁻²` residues. Contacts in `C⁺₂, C⁻₃` (absent in flat space) likewise contribute `O(Λ⁻²)` relative. With EFT power counting
`g_f = O(G_f)` [A, CHMRS's technical assumption] all residues are `O(G)` in `Ω`, i.e. `O(Λ⁻²)` relative to the leading terms. With the
flat two-sided bounds `g_f ≲ G_f log Λ`, they are `O(log Λ/Λ²)`.

**2.6 Assembled light side** (tree level, `πG = 32/c_T` in the scalar-normalised convention, up to the overall `l = 2` constant):

`−Ω_light = πGΛ²[16 L_AdS(Λ) − k_AdS(Λ) x²] − 23690πG c₁ x cos φ + O(G) + O(G²Λ⁴)`,
`L_AdS(Λ) = ∫₀¹ψ₂/(p² + 9/(4Λ²)) = 65 log(2Λ/3) − 1171/15 + 30π/Λ + (2007/4)log(2Λ/3)/Λ² − 4389/(8Λ²) + 6993π/(16Λ³) + …`.

The closed form is exact (`t3`): with `a = 3/(2Λ)`, `ψ₂ = (p²+a²)Q_a + αp + β`, `α = 65 + 223a² − 14a⁴`, `β = a²(40 + 259a² + 47a⁴)`, one has
`L_AdS = ∫₀¹Q_a + (α/2)log(1 + a⁻²) + (β/a)arctan(1/a)`. It agrees with quadrature to 15 digits at `Λ = 3 … 1000`.

## 3. Double traces `[TT]_{n,ℓ}` (§0 item 6)

[Q CHMRS §(Mellin), after (BknuContactsLeading)]: for scalars, `C_{k,ν}` receives contributions only from the anomalous dimensions of the
first `k/2` double-twist trajectories and from the anomalous OPE coefficients of the leading one. All other families sit at double zeros
(`2sin²(π(τ−2Δφ)/2)` in the measure, with no compensating pole) and enter at `O(γ²) = O(G²)`. [I] for `⟨TTTT⟩`: the same holds for each
helicity structure, with the double-twist points of that structure, and the number of `O(1/c_T)` families is `⌈k/2⌉` for the
odd kernels (not derived). Caveat: the `[TT]` component twists are `2 + n` for **all** integers `n ≥ 0`, both parities, while a single
`sin²(π(τ−τ₀)/2)` vanishes on one parity class only. So it must be shown that each spinning kernel (even `D₁⊗D₁`, odd `D₅⊗D₅`, per helicity
structure) has double zeros at *its own* families. **This is not shown**; if it fails, some MFT families enter at `O(1/c_T)` or even `O(1)`. The `O(1/c_T)` data is that of the tree-level EFT, so it is already contained in §2. The explicit evidence is §1.2(c): the clean
`1/(ν²+9/4)` requires dividing by `γ²`, whose double zeros remove exactly the double-trace poles. At `O(1/c_T²)` (one loop) light double
traces contribute with positive measure `≈ (π²γ²/2)p`, but the functional's action on light states has no sign. This is an error term of
relative size `GΛ² = 32Λ²/(πc_T)` (log-enhanced), so the tree-level statement needs **`c_T ≫ Δ_g²`** [A]. So the naive claim is refuted
as stated, but the correct version (no independent light unknown at `O(1/c_T)`) holds.

## 4. Conditional AdS4 statement

Insert §2.6 in `Ω|_heavy = −Ω|_light ≥ 0`. Then `k_AdS x² ≤ 16L_AdS + O(x/Λ²) + O(log Λ/Λ²)`, which gives (`t3`):

**`|ĝ₃|²Λ⁸ ≤ 37.8387 log Λ − 60.7875 + 54.865/Λ + O(log Λ/Λ²)`**,  `t₄² + α_T² ≤ 57600 × (same)/Λ⁸`.

| `Δ_g` | 10 | 20 | 50 | 100 | 1000 |
|---|---|---|---|---|---|
| bound on `√(t₄²+α_T²)` | 0.24 | 1.5·10⁻² | 4.1·10⁻⁴ | 2.7·10⁻⁵ | 3.4·10⁻⁹ |

The value at `Δ_g ≤ 20` is meaningless: there the `O(Λ⁻²)` terms are uncontrolled.

**Dependence on the heavy positivity threshold (M2).** The table assumes `Ω|_heavy ≥ 0` down to `m ≈ Λ`. `ADS4_FUNCTIONAL.md` gives the following:
* The leading Regge action is negative near `m ≈ 6.3Λ`, `Λβ ≈ 8.7` (`Δ_g ≥ 50`).
* The Regge-side statement R gives positivity for `m ≥ 7.4Λ`. With envelopes the threshold is `≈ 8.1Λ`, and `≈ 22Λ` at large `Δ_g` (referee's reading).
* Positivity below that rests on the [M] hybrid model.
* The `+−` near-zero needs `Λ ≳ 444` if the sign of the `β²` term is not trusted.

If positivity is available only for `m ≥ nΛ_eff`, take the smearing scale `Λ_eff = Λ/n` [I: this assumes the threshold scales with the smearing
scale]. The bound becomes `√(t₄²+α_T²) ≤ 240√(A log Λ_eff + B_AdS + C₁/Λ_eff)/Λ_eff⁴`, about `n⁴` weaker (`t4`):

| `Δ_g` | `n = 1` | `6.3` | `7.4` | `8.1` | `22` |
|---|---|---|---|---|---|
| 50 | 4.1·10⁻⁴ | 0.33 | 0.56 | 0.75 | — |
| 100 | 2.7·10⁻⁵ | 2.8·10⁻² | 4.9·10⁻² | 6.8·10⁻² | 1.7 |
| 300 | 3.8·10⁻⁷ | 4.4·10⁻⁴ | 8.1·10⁻⁴ | 1.1·10⁻³ | 4.6·10⁻² |
| 1000 | 3.4·10⁻⁹ | 4.4·10⁻⁶ | 8.1·10⁻⁶ | 1.2·10⁻⁵ | 5.2·10⁻⁴ |

`Λ_eff ≳ 20` is also needed for the `O(Λ_eff⁻²)` light-side terms. The sign-agnostic `+−` near-zero condition `Λ_eff ≳ 444` needs `Δ_g ≳ 3300` at
`n = 7.4`. **Honest headline:** `2.7·10⁻⁵` at `Δ_g = 100` is conditional on positivity down to `m ≈ Λ`, which is open and contradicted at leading Regge
order. With the currently available Regge-side threshold (`n = 7.4`, itself conditional on H_env), the bound is `≈ 5·10⁻²` at `Δ_g = 100`.

Inputs:
| input | status |
|---|---|
| pole at `ν² = −9/4`; transverse space `H₂` | [D], pole ↔ dimension [DV] |
| pres–pres residue `t₄`, `α_T`-independent | [D] from [Q] MP (TTdiagonal), CDP |
| `A = 1040/k(0) = 37.8387` | [D] given (F) |
| Einstein term exactly `16πG/(ν²+9/4)` ⇒ `B_AdS = −60.7875`, `C₁ = 54.865` | [D] given (F) and **(R) [I]**. If (R) fails, with `−C⁺₂ = 16πG[1/(ν²+9/4) + r]` and `r = O(ν⁻⁴)`: `B_AdS` shifts by `(1040/k(0))∫₀^∞ν r dν`, and `C₁` by `−(640/k(0))∫₀^∞ν² r dν` (from `ψ₂ = 65p − 40p² + …`) |
| **(E)**: no Einstein term in `C⁻₃, C⁺₄, C⁻₅`; no non-polynomial residue in any `|ĝ₃|²`/contact term (spinning carry-over of CHMRS (eftformula)) | **[I]** (scalar analogue [D given Q]); leverage `~Λ³` (`C⁻₃`) and `~Λ⁸` (forward terms), §2.1 table: **load-bearing** |
| scalar analogue of all of the above, incl. normalisation `8πG = 256/c_T` | [DV] |
| `|ĝ₃|²` forward terms exact; `g₅` cancellation exact | [D] given (F), (P), (E) |
| `k_AdS − k(0) = O(Λ⁻²)` with `b₂, f₁, …` | [D] structure; **coefficients not computed** |
| linear `c₁t₄` term | [D] structure; `c₁` not computed; effect `O(x/Λ²)` |
| (F) spinning flat limit of the light side; (P) polynomial light side of spinning contacts and `t`-channel diagrams | [I] (Q for scalars) |
| `[TT]` at `O(1/c_T)` contained in the EFT light side; one-loop error `O(Λ²/c_T)` | [Q] scalars / [I] spinning |
| `Ω|_heavy ≥ 0` for all heavy states down to `m ≈ Λ` (else use `Λ_eff = Λ/n`, table above; heavy `⟨TTTT⟩` formula, spinning bulk point, uniform remainder, `k = 2, 3` swappability, parity-odd heavy exchanges) | [A] (`ADS4_FUNCTIONAL.md` §4 gaps 1–4, 8: open) |
| no light single-trace operators other than `T` below `Δ_g` (or their AdS contribution `≥ 0`) | [A] |
| EFT power counting of contacts; tree level `c_T ≫ Δ_g²` | [A] |

Not done: the spinning Witten diagrams themselves (Einstein `u`-channel in `C^±` normalisation for (R), `W³⊗W³` `t`-channel for `b_i, f_i`,
mixed for `c₁`). A spinning Mellin or conformal-Regge computation of these three diagrams would turn every [I] above except (A)-type inputs into [D].

## 5. Code (`ads4/code/adslight/`, `bash run_all.sh`, ≈ 5–6 min; `set -o pipefail`)
| file | content | log |
|---|---|---|
| `t1_scalar_mellin.py` | CHMRS Mellin `C_{2,ν}`; contact calibrations; graviton = `T` poles ⇒ `8πG/(ν²+9/4)`; `H₂` propagator; `G ↔ c_T` | `cache/t1.log` |
| `t2_flat_split.py` | flat `−B₂…−B₅` by channel (sympy, exact); moments `μ₂_j, μ₃_j`; `k_AdS(Λ)` in terms of the AdS coefficients; linear-`t₄` term | `cache/t2.log` |
| `t3_light_omega.py` | exact `L_AdS(Λ)`, expansion, check vs quadrature; constants `A, B_AdS, C₁`; CFT3 form | `cache/t3.log` |
| `t4_sensitivity.py` | leverage of residues not excluded by (F) (M1); bound vs positivity threshold `n` (M2); linear-`t₄` shift vs `C₁/Λ` (m1) | `cache/t4.log` |

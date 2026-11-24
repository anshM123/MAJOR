# GAP_UNIFORM: can the near-gap beyond-leading-Regge corrections of Ω^EE be bounded uniformly for all Δ_g ≥ 1000?

*24 Nov 2026, gap-uniformity agent. Unreviewed, not committed. No existing file was edited. Nothing in `gapnear/` was touched (its
caches were read only).*
* *New code:* `ads4/code/gapuni/` (`run_all.sh`; ≈ 0.6 CPU-h in total; at most 2 processes; `ulimit -v 3000000`; `OMP_NUM_THREADS=1`).
  Every script installs `gu_guard.py`, which refuses writes outside `gapuni/cache/`. That directory's `.gitignore` ignores everything
  but itself.
* *Read-only imports:* `poserr/large_m_map` (→ `spinfunc/ee1000_common.hyb`, the chain's hybrid states and node weights),
  `agproof/ag_exact` (Arb exact scalar), `p4win/p4w_common` (Mehler–Dirichlet P), `prl/data/nf02_psi.json`.
* *Read-only data:* `gapnear/cache/gn_hmono_*.json` and `gn_bind_base1000.json` (the other agent's K₂ rows at Δ ≈ 1000+J, 1600 and
  2000+J; I only read them).
* *Literature:* LaTeX sources and arXiv PDFs of 2205.01495, 2106.10274, 2203.11170, 2311.04271, 2210.15683, 2312.02273 (equation
  numbers below are those of the arXiv PDFs). INSPIRE citation lists read by title: 2106.10274 (141 records), 2311.04271 (15),
  2203.11170 (22), 2205.01495 (103); arXiv keyword searches.

Grades: **[P]** proved; **[certified]**; **[D]** derived; **[D-sketch]**; **[DV]** derived + float check; **[Q]** quoted; **[NE]**
numerical evidence; **[M]** model/convention; **[A]** assumption; **[I]** open. Notation as GAP_SCALING: `Λ = √((Δ_g−2)(Δ_g−1))`,
`Λ₀ = Λ(1000) = 998.50`, `μ = m/Λ`, `x = ν/m`, `p = ν/Λ = μx`, `ξ = ν/Δ`, `Γ = m² × (beyond-hybrid correction)/(flat envelope)`
("c₂ units"). **Chart A** = same `(J, μ, h)` at a different gap.

---

> **Update (GAP_C3.md, unreviewed): the Laplace remainder constant C₃ is explicit and certified (Arb Taylor-model evaluator on a complex polydisc around the crystal; Lemma L3 [P]), but it is ~10⁸–10¹² too large to use at Δ_g = 1000. It gives a certified threshold Δ₁ ≈ 5·10⁵–10⁷ above which K₂ truncation is accurate enough. The premise B_F is still needed on 1000 ≤ Δ_g ≲ 10⁶.**
> * **Certified:** the scalar k=2 integrand with Q_J at J ∈ {0, 2, 4–8} × ξ ∈ {0.1 … 0.85} (35 nodes), giving Δ₁ = 4.7·10⁵ … 10⁷; graviton EE and C3H at J = 0, 2, giving Δ₁ = 1.0–1.8·10⁶.
> * **Eventual positivity:** GAP_UNIFORM's Λ\* is now explicit, Δ_g\* ≲ Δ₁, at the certified nodes. This holds given:
>   * Lemma R (contour) [D-sketch];
>   * the chain propagation factor Π ≤ 3 [M];
>   * the K₂ part of Γ for all m [NE];
>   * continuity in ξ [D-sketch].
> * **The obstruction is the method:** the genuine cubic size of Φ_ξ (S₃ ≈ 10–12) caps the radius at ~0.04 against a Gaussian width of 0.032 at Δ = 1000. Ideal sup norms give the same Δ₁. The true next Laplace term is ~2·10⁻⁵ relative at Δ = 1000, so K₂ is accurate and the certificate is pessimistic by ~10¹⁰.
> * **Next lever:** a certified holomorphic Morse chart (rough estimate Δ₁ ≈ 3·10³–1.3·10⁴), ξ as a Taylor-model variable, and the J ≥ 4 two-structure / odd-J graviton amplitudes.

## 0. Verdict

**No unconditional uniform bound was obtained. What was obtained is a sharp reduction, one proved lemma, a structural derivation that
makes H_mono-A a standard Laplace statement, and an effective statement for J ≤ 7 that rests on one explicit "no-growth" premise
with a ×12–×3·10⁴ safety factor.** Nothing in the literature gives an explicit (constant-bearing) bound on heavy-block corrections to
the flat limit; every source is asymptotic and non-effective (§1).

| # | statement | grade | coverage |
|---|---|---|---|
| **T** | **Chart-A tolerance lemma.** With the chain's own definitions, the beyond-hybrid correction enters a state as `dF = Σ_n a_n Γ(x_n)/m²` with node coefficients `a_n` (and the flat part, N_f, the node positions `x_n = p_n/μ`) **exactly Λ-independent**; the Lemma-C bound decreases in Λ. Hence the sign-free tolerance on `sup_x|Γ|` satisfies `G_cert(Λ) ≥ G_cert(Λ₀)(Λ/Λ₀)²`. Positivity at gap Λ follows from `sup_x|Γ_Λ(x)| ≤ G_cert(Λ₀)(Λ/Λ₀)²` (and the window terms). | **[P]** (structure) + [DV] (numbers, 366 states J ≤ 60 at 1000/2000/5000: `G_cert(Λ)/(Λ/Λ₀)² ≥ G_cert(Λ₀)` everywhere; hybrid-floor form within 0.976–1.000 of the Λ² law) | every near-gap state, every Λ ≥ Λ₀ |
| **S** | **Fixed-phase structure of chart A.** At fixed `ξ = ν/Δ` the crystal integral is `∫_{D′} e^{Δ·Φ_ξ(z)} g_J(z; Δ) dz` with `Φ_ξ = (iξ/2) log A + ½ log B` **exactly Δ-independent** and `g_J` independent of Δ (scalar) or rational in Δ (EE ⟨TTO⟩ coefficients). So "Γ in c₂ has a finite m → ∞ limit at fixed (J, x)" (H_mono-A) is the ordinary Laplace asymptotics of a fixed-phase integral, with remainder `O(1/m)` in c₂, **provided** the crystal is the only contributing critical point of `Φ_ξ` at or above its height. | **[D]** (code-level identity) + [Q] (Laplace theorem) | x ∈ compacts of (0, 1) |
| **R** | **Cycle rigidity in chart A.** The descent flow of `Re Φ_ξ` is Δ-independent, so the thimble coefficients of the physical chain D′ (= PCYC's identity) are Δ-independent integers (up to `1 + O(e^{−πξΔ})` Hankel factors at the ν-dependent atom exponents). A direct-vs-K₂ agreement at **one** Δ therefore excludes a higher or near-degenerate non-crystal thimble at **every** Δ, except for the non-generic case of an accidentally tiny amplitude. The 1000 direct checks (J = 5, 6, 8, 12, 16, 29) are such tests. | **[D-sketch]** + [NE] | ξ at the checked nodes |
| **E** | **Eventual positivity (non-effective).** Given R and the window behaviour (§4.5), every fixed near-gap state (J, μ > 1, h) is positive for all Λ ≥ Λ\*(J, μ) < ∞. Λ\* is not computed. | **[D-sketch]** | every (J, μ > 1) |
| **U7** | **All Δ_g ≥ 1000 for J ≤ 7** (states J = 0, 2, 4, 5, 6, 7 at the 1000 twist edge and m/Λ = 1.02, 1.05, 1.1, 1.3, plus their seam images): positive at every Δ_g ≥ 1000 **if** `sup_{x∈[0.05,0.9]}|Γ_{Δ_g}(x)| ≤ F_J·S_J` (premise **B_F**; Γ the true correction, `S_J` the measured K₂ `sup|Γ|` at 1000, J = 7: at 2007) with the window conventions of NEARGAP2 [M]. Seam states (§3.3) have tolerances ≥ 0.97 × the scaled ones for J ≤ 7. **F_J = 11.8 (J = 7), 14.9 (6), 43 (4), 743 (5), 799 (2), 3.4·10⁴ (0)** against the certified flat floor minus Lemma C. | **[P given B_F + windows]**; B_F **[NE]** | J ≤ 7 near gap |
| B_F | Evidence: graviton K₂ rows 1000 → 2000: `sup|Γ|` ratio 0.97–1.06 (J ≤ 8); K₂ step (last Laplace term) m-independent in c₂ (44/39, 593/602, 350/358, 1232/1245, 2029/2048); exact scalar (Arb, no asymptotics) 1000 → 10⁴ in chart A: `sup|Γ|` drift +1.2 % (J = 5), −0.2 % (12), −2.6 % (26), 0 (8 at 1.3Λ), worst single node +3.1 %; `Λ²S` converges; approach ∝ 1/m² in the bulk (§5.2). Measured direct − K₂ at 1000 (≤ 13 c₂ at J = 5, ≤ 1.9σ at J = 6) is ×50–×2000 below what B_F tolerates. | [NE] | — |
| ✗ | **J = 8:** sign-free certified version fails (F = 0.85 at 1.02Λ); hybrid-floor version F = 2.7 [NE]. **J ≥ 12:** sign-free fails at Δ_g = 1000 already (F_hyb = 0.46 at J = 12, 0.12 at 16, 0.06 at 21); only signed budgets work, and at J ≥ 16 the K₂ rows cannot be error-bounded at 1000 (L < 1), so there is nothing to transfer from K₂ alone. **Prop. D′ (RD_CLOSE)** cannot be used in the bulk of chart A: its far-region bound is absolute, while the signed integral is `e^{−πν}` smaller, so it controls only `ξ < 2 artanh(r/2)/π = 0.107` (§4.6). | [D] | — |

**What remains for a proof (in order of value).**
1. **The explicit Laplace remainder constant.** By S, every hypothesis of a quantitative Laplace theorem refers to the
   Δ-independent phase Φ_ξ and to the amplitude on a **fixed** complex polydisc around the crystal. So one sup-norm evaluation per node
   (interval arithmetic on an explicit elementary function in 7 complex variables, ~36 nodes per spin) would give
   `|Γ_m(x) − Γ_{K₂,m}(x)| ≤ C₃(ξ, J)/m` with an explicit C₃, decreasing in m. This is "uniformity for free": verifying at Δ₀ covers
   every Δ ≥ Δ₀. With it, U7 becomes [D + certified] modulo R (and B_F is no longer a premise). Not attempted here: it needs a
   pointwise complex evaluator of the spinning EE/C3H amplitudes, which the project does not have (only TPS around the saddle and the
   J = 0 pointwise integrand of `repro_engine`).
2. `sup_{m≥m₀} |Γ_{K₂,m}(x)|`: K₂ is an explicit function of 1/m; rows at a third gap (the other agent's 5000 rows) plus
   analyticity in 1/m would bound it.
3. Lemma R's genericity residual and the windows (ν < 50: Regge engine + CHMRS's matching property; x > 0.9: fold distance
   `1 − 1/μ²`, harmless at fixed μ > 1, Airy-type as μ → 1 but killed by ψ_k ~ (1−p)⁵, §4.5).
4. J = 8–16: either certify the hybrid floor (GAP_SCALING: "≈ 1 day") and use the incremental form U_Δ (§2, §3.4), or compute C₃.

---

## 1. Literature: what is (and is not) available

### 1.1 Caron-Huot, Li, Parra-Martinez, Simmons-Duffin, arXiv:2205.01495 (PRL letter)
* The only AdS statement is one paragraph before eq. (30): "As explained in [CHMRS], since the scattering processes under
  consideration take place in a region of small size ∼ 1/M, flat-space dispersive bounds uplift in AdS to rigorous bounds on
  holographic CFTs, up to corrections suppressed by 1/(M R_AdS) = 1/Δ_gap."
* Eq. (30): `|a − c|/c ≤ 23/Δ_gap² + O(1/Δ_gap⁴)` (AdS₅/CFT₄), with the remark that it "could potentially be improved at the ∼ 5 %
  level".
* **Assessment.** The uplift is asserted by citation; the quantity is the final coupling bound; the correction is an unquantified
  `O(1/Δ_gap⁴)` (there is no constant and no statement of uniformity over heavy states). Nothing in the letter can be made explicit: it
  contains no heavy-state analysis. [Q]

### 1.2 Caron-Huot, Mazáč, Rastelli, Simmons-Duffin (CHMRS), arXiv:2106.10274, §3–5 and App. C, F
* **Exact representation.** (3.10): the action of any kernel on a block is an exact spacetime integral with elementary factors; (3.13):
  its heavy limit is a Laplace transform. In the symmetric frame (3.23) the Regge moment is exact, (3.25)–(3.26) (Appendix C, (C.5)–(C.8)).
* **Large-m expansion at fixed (ν, β).** (3.27): first 1/m² correction; a series to `1/m⁸` is in the ancillary file
  `anc/regge_moments_order8.m`. Stated properties (text after (3.27)): (i) η_AdS enters only through `P^{(n)}(η)` with polynomial
  coefficients; the coefficient of `P^{(n)}/m^{2j}` grows no faster than `η^{j+n}`, "so the 1/m series remains applicable in the Regge
  limit no matter the spin"; (ii) the ν-dependence at each order is entire, growing no faster than `ν^{2j}` near η → 1, "so the
  regimes ν → 0 and ν ∼ m are smoothly connected by a series in ν/m". (3.28)–(3.30): at fixed (ν, J), `(m^{2k}/2)C_{k,ν} = 1 +
  Q_{k,2}/m² + Q_{k,4}/m⁴ + …` with `Q_{k,n}` polynomials of degree n in ν²; the leading large-ν terms (3.30) match flat space.
* **Bulk-point saddle.** (3.31)–(3.38): four saddles (3.33); for ν < m the contour is homologous to the thimble of `p_{−,−}` (argued
  from 1D steepest-descent plots, Fig. 8, exact in d = 2 where the integrand factorises); result (3.35), (3.38). "Since we start with
  an integral over real spacetime, we know *exactly* when the saddle point is valid, namely for |ν| < Δ_gap" (§3.4) — a statement about
  which saddle dominates, not an error bound. App. F, (F.1)–(F.3): the large-(Δ, J, ν) saddle and formula `S_ν`.
* **The uplift argument (§5).** (5.7): `C_{k,0}[Δ, J] = 2P(cosh β)/m^{2k} × (1 + O(1/M²))`, "assuming … M is sufficiently large";
  (5.11): two-sided bounds "× (1 + O(1/Δ_gap²))", with "breathing room to deal with 1/m² effects in (5.7), for example by slightly
  increasing the coefficient of C_{2,0}". §5.3, (5.26)–(5.29): `ω_AdS = ∫ f(ν/M) C^{imp}_{2,ν}`; "Suppose that ω_flat is heavy-positive.
  We will argue that, **for sufficiently large M**, ω_AdS, or perhaps a small correction to it, is heavy-positive." "By construction,
  ω_flat and ω_AdS agree up to 1/M² corrections in the bulk point regime of fixed J; hence ω_AdS is automatically positive in that
  regime for sufficiently large M." The remaining regimes are handled by two asymptotic dominance arguments: improvement terms for
  `β ≫ 1/M` ((5.30)–(5.35), using monotonicity of the ratio (5.35) in β), the unimproved term for `m ≫ M` ((5.36)–(5.44), endpoint
  asymptotics and the functions `Q_n(β)`), under the endpoint condition (5.29) `f(p) ~ (1−p)^l`, `l > n + 1 − d/2`. App. F.1 checks
  that the oscillatory endpoint term is not enhanced for m ∼ M.
* **Truncation errors (§5.2):** the technical assumption (5.3) (`Σ O(g/R²) = O(1/R²)`) is supported by a heuristic `n^j(1−ε)^n`
  estimate; "A rigorous proof … is left to the future."
* **Assessment.** The uplift is an **existence** statement ("for sufficiently large M"), with all corrections `O(1/M²)` at fixed J
  and no constants. Nothing bounds the 1/m² remainder uniformly in (J, ν) or gives its constant. What *is* usable: (a) the exact
  integral representations (3.10)/(3.25) (the project's crystal integral is the CLS-D spinning analogue); (b) the structural
  properties after (3.27), which are exactly the chart-A/chart-B structure (fixed-x series with Δ-independent coefficients); (c) the
  endpoint condition (5.29): NF02's ψ_k vanish like `(1−p)^5` at p = 1 (checked here), which is what makes the light-cone edge
  harmless (§4.5). [Q]

### 1.3 Knop & Mazáč, arXiv:2203.11170 (AdS₂)
* **Mechanism of uniformity.** Every functional is written as (3.21) and, by contour deformation, as **(3.25)**
  `ω[G] = ∫ dz f(z) dDisc G(z)`, i.e. (3.26) a positive-kernel integral against the double discontinuity of the block. The kernels
  (4.6)–(4.7) of `C_k^{B,F}` are explicit ₃F̃₂'s. **(4.8)** `f_k(z) > 0` on z < 0 is *not proved* ("we were able to verify it by
  numerically computing the kernels in many cases, and also by observing positivity of coefficients in their Taylor expansion around
  z = ∞"). Given (4.8), **(4.9)** `C_k[G_Δ] ≥ 0` for **all** `Δ ≥ 2Δ_φ + k` — positivity uniform in Δ (hence in Δ_gap) because it is
  pointwise in z: there is no smearing over an oscillating impact-parameter kernel.
* **Heavy asymptotics** are only asymptotic: (3.36) `ω[Δ] = Σ_{n<N} b_n Π_{k+2n}[Δ] + O(Δ^{1−2k−4N})`; (4.19)
  `C_k[Δ] = 8/(π[Δ(Δ−1)]^{k−½})·[1 + (48Δ_φ(Δ_φ−1)k − 4k³ + 4k + 3)/(24Δ(Δ−1)) + O(Δ⁻⁴)]`; bounds (4.25)–(4.27) and (4.30) carry
  `O(M⁻⁴)`; "assuming Δ_gap is sufficiently large, so that the leading-order behaviour (4.19) can be trusted".
* **Assessment.** The uniform-in-Δ_gap control comes from manifest (pointwise) kernel positivity, which d = 1 allows. For Ω^EE the
  ν-smearing against oscillating `P_{−½+iν}(cosh β)` (or its spinning analogue) destroys pointwise positivity; that is the reason the
  project needs heavy-state numerics at all. The idea transfers only to a *different* functional with a manifestly positive spacetime
  kernel; Ω^EE is fixed. Nothing quantitative transfers. [Q]

### 1.4 Chang, Landau, Simmons-Duffin (CLS-D), arXiv:2311.04271
* §3.3: `f_{ν,m}` (3.27); integral (3.30); "sixteen different saddles" (3.31); dominant locus for ν < m (3.32) ("the analysis is
  essentially identical to [CHMRS]"); expansion `f = −πν + iν log 2 − ½ wᵀMw + …` (3.33) with `det M` (3.34); leading result
  (3.35)–(3.37); shadow coefficient at large Δ by a second saddle (3.38)–(3.43); the "scattering crystal" (3.50); final heavy action
  (3.51)–(3.52) and spinning analogues (§3.5).
* **Only the leading order is computed; there are no corrections and no error estimates.** Footnote in §1: "The 'uplifting' procedure
  only guarantees positivity of the functional in the 'bulk point' regime … One must additionally check other regimes to derive a
  rigorous CFT bound." §5 (Discussion), remaining tasks: (i) the CFT sum rules not expressible through shock amplitudes; (ii) "check
  the positivity of the sum rules in the Regge limit, where Δ and J are both large"; (iii) the light-side dictionary.
* **Usable here:** (3.33) gives `Re f = −πν` at the crystal — the `e^{−πν}` cancellation that limits any absolute-value bound on the
  real section (§4.6). The crystal formula is what the project's K₂ series expands around. [Q]

### 1.5 Other literature
* **van Rees & Zhao, arXiv:2210.15683.** Rigorous flat-space statements (dispersion relation, unitarity, Froissart–Martin) under the
  *assumption* that the limit exists in a Euclidean region (eq. (9)–(10)); tools: positivity of the spectral density (16) plus Montel
  compactness (footnote below (17)) — existence of limits, no rates. Their footnote to (15) [the block-sum representation] is directly
  on point: "This would follow if the subleading corrections to … (14) are uniformly small for all h and h̄. We believe this can be
  proved rigorously, but it is not immediate: the large Δ limit of the conformal block as given in equation (13) is not uniform in r."
* **van Rees & Zhao, arXiv:2312.02273** (OPE density → partial waves): evidence-level, explicitly delineates divergence regions.
  **Komatsu, Paulos, van Rees, Zhao, arXiv:2007.13745**: 'bad regions'/AdS Landau diagrams where the flat limit diverges; no rates.
  **Gadde & Sharma, arXiv:2204.06462**: formal 1/R perturbation theory. **Y.-Z. Li, arXiv:2106.04606**: frameworks, no errors.
  **Caron-Huot, Chakravarty, Namjou, arXiv:2502.14963**: bulk-point factorisation in general geometries, leading order.
  **Chang, Dommes, Kravchuk, Poland, Simmons-Duffin, arXiv:2509.14307**: error norms for block approximations in numerical bootstrap
  (positivity for all Δ by construction of approximants) — a different setting (fixed external data, not the flat limit).
* **Mathematics.** Explicit-constant Laplace/steepest-descent bounds exist in 1D (Olver, *Asymptotics and Special Functions*, ch. 3–4)
  and, non-explicitly, in n dimensions (Fedoryuk, *The saddle-point method*; Pham's Picard–Lefschetz theory; Witten, arXiv:1001.2933,
  for the multivalued/twisted Picard–Lefschetz setting). These are the inputs used in §4.
* **Searches** (24 Nov 2026): the INSPIRE citation lists above, and arXiv abstract searches for flat-space-limit error bounds, uniform
  heavy-block asymptotics, finite-gap dispersive bounds. **Nothing found gives an explicit, uniform-in-Δ_gap bound on heavy-state
  corrections to a flat-space functional.** The closest statements are CHMRS §5 (non-effective) and Knop–Mazáč (4.9) (uniform only by
  pointwise positivity in d = 1, itself numerical).

---

## 2. The inequality that is needed

Per state the chain value is (ee1000_common.hyb / dF_state / lam2)
```
λ(Λ)·N_s = B₂ + B₃ + cv₂(Λ) + cv₃(Λ) + dF(Λ),      dF = Σ_n [p2_n env(h,x_n,2) Γ_C2(x_n) + u p3_n env(h,x_n,3) Γ_C3H(x_n)] / m²
```
with `p_k,n = ψ_k(p_n) w_n` on fixed nodes `p_n ∈ (0,1)`, `x_n = p_n/μ`, `u = μ⁻²`, `m = μΛ` (+ off-diagonal analogue for even J).
Define `A = Σ_n (|p2_n| env₂ + u|p3_n| env₃)` and the per-N_s cost of a uniform |Γ| = 1 c₂, `W = A/(m² N_s)`.

**Needed (sign-free form, U_sf).** For all Λ ≥ Λ₀:
`sup_x |Γ_{Λ}(x; J, h)| · A/(μΛ)² + (window terms) < N_f/5 − LC(Λ)` ⇐ `sup_x|Γ_Λ| ≤ G_cert(Λ₀)·(Λ/Λ₀)²` with
`G_cert = (N_f/5 − LC)/(W N_s)` (minus the window share, §3.4).

**Needed (incremental/signed form, U_Δ).** With the 1000 state value and budget `λ_K(Λ₀), E₀` known:
`λ(Λ) ≥ λ_hyb,min − (Λ₀/Λ)²[|δλ₀| + E₀ + W₀·sup_x|Γ_Λ(x) − Γ_{Λ₀}(x)|]`, so positivity needs only that the **change** of Γ grows
slower than Λ².

Either form says the same thing: **the beyond-hybrid relative correction `Γ/m²` must not grow** (it may even stay O(1) relative
to its 1000 size times F); a 1/m² decay (H_mono-A) is far more than needed.

## 3. Lemma T and the tolerances

### 3.1 Lemma T [P]
*Statement.* At fixed (J, μ, h): (i) `B₂, B₃, N_f, A` and the nodes are Λ-independent; (ii) `LC(Λ)` (Lemma C) is non-increasing in Λ
(GAP_SCALING §2.1 [P]); (iii) hence `G_cert(Λ)·(Λ₀/Λ)² ≥ G_cert(Λ₀)` whenever `G_cert(Λ₀) > 0`; (iv) if `sup_n|Γ_Λ(x_n)| < G_cert(Λ)` (1×1)
then λ(Λ) > 0; for the even-J 2×2 matrix (the hybrid base matrix is diagonal, `ee1000_common.lam2`), by Weyl, `λ_min ≥ min_h rig_h − sup|Γ|(max_h W_h + W_off)`, giving `G_W`.
*Proof.* (i): `ψ_k(p)` are fixed functions of p (the functional's definition); `env(h, x, k)` and the flat parts depend on
(J, u, h, x) only; the factor `1/m² = 1/(μΛ)²` is explicit. (ii) is GAP_SCALING's monotonicity. (iii): `G_cert = μ²Λ²(N_f/5 − LC)/A`.
(iv): `|dF| ≤ A sup|Γ|/m²` and `B₂ + B₃ + cv ≥ N_f/5 − LC` (CERT_CLOSE flat certificate + Lemma C). ∎
*Caveat.* The functional's min-norm de-shift weight differs from ψ₂ by `S(p)/Λ² + O(Λ⁻⁴)` on p ≥ 0.1 (GAP_SCALING §1.2) and by a
Λ-independent inner profile at ν = O(1); that changes A by a factor `1 + O(Λ⁻²)` in the bulk and enters only the ν < 50 window
(gauge factor, cost ∝ Λ⁻⁴). [DV]

### 3.2 Numerical check (`gu1_tol.py`, `gu1b_weyl.py` → `cache/gu1_tol.txt`, `gu1_summary.txt`, `gu1b_weyl.txt`) [DV]
* 366 states, J = 0–60, m/Λ ∈ {1000 twist edge, 1.02, 1.05, 1.1, 1.3}, both helicities: `W·Λ²` at 2000 / 5000 equals the 1000 value to
  ≤ 10⁻⁴ for J ≤ 10 and ≤ 1.3·10⁻² up to J = 60 (all of it the drift of N_s through the LR curvature; A itself is exactly invariant). `G_cert(Λ)/(Λ/Λ₀)²` is ≥ 1 wherever `G_cert(1000) > 0` (1.0000–1.15 for J ≤ 13; it
  grows because LC shrinks). The hybrid-floor tolerance `G_hyb = λ_hyb/W` scales as Λ² within 0.976–1.000 (λ_hyb drifts down as the
  positive LR curvature decays).
* `G_cert(1000) > 0` up to J = 13 at the twist edge (J ≤ 11 at every tabulated mass), as GAP_SCALING §2.2 found; at 5000 for J ≤ 23.

| J | 0 | 2 | 4 | 5 | 6 | 7 | 8 | 10 | 12 | 16 | 21 | 26 | 30 | 40 | 50 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G_cert, min over masses (c₂) | 1.99e5 | 6.34e4 | 1.10e4 (W) | 3.01e4 | 3.60e3 (W) | 2.34e3 | 316 (W) | 153 (W) | 30.5 (W) | < 0 | < 0 | < 0 | < 0 | < 0 | < 0 |
| G_hyb, min over masses (c₂) | 9.97e5 | 3.17e5 | 1.12e4 | 3.59e4 | 4.17e3 | 6.28e3 | 994 | 636 | 458 | 232 | 166 | 92 | 62 | 29 | 16 |

(W = 2×2 Weyl form; odd J ≥ 5 and J = 0, 2 are 1×1. The binding state J = 5 at its twist edge: 3.57·10⁴, matching NEARGAP2's 3.6·10⁴.)

### 3.3 Seam states (`gu3_seam.py` → `cache/gu3_seam.txt`) [DV]
Γ is a property of the **block** (x, J, m); it does not know Λ. A seam state (J, μ_s, Λ′) with μ_s between the Δ′ and 1000 twist
edges therefore needs Γ(x; J, μ_sΛ′) on x ≤ 1/μ_s, i.e. the same function as its chart-A neighbours plus a thin strip
x ∈ (1/μ_edge(1000), 1/μ_s]. Computed at Δ′ = 2000, 5000, 10⁵ for J ≤ 16: the seam tolerances are 0.65–1.0 (hybrid) / ≥ 0.75
(certified, J ≤ 12) × the scaled 1000 twist-edge tolerance, and the strip carries ≤ 10⁻³ of the state's sign-free weight (≤ 0.9/900 in
Λ²W units at J = 16, 10⁵). **The seam is not a separate problem in the formulation of §2**: it is the x → 1 edge window of the
chart-A states (§4.5).

### 3.4 Growth allowances (`gu5_cover.py`, `gu2_growth.py` → `cache/gu5_cover.txt`, `gu2_growth.txt`)
`sup|Γ|` over all diagonal kinds/helicities on x ∈ [0.05, 0.9] (gapnear K₂ rows; J = 7: the 2007 row for both):

| J | sup\|Γ\| ~1000 | ~2000 | ratio | G_cert | G_hyb | **F_cert = G_cert/sup** | F_hyb |
|---|---|---|---|---|---|---|---|
| 0 | 5.9 | 5.8 | 0.98 | 1.99e5 | 9.97e5 | **3.4·10⁴** | 1.7·10⁵ |
| 2 | 79.3 | 77.0 | 0.97 | 6.34e4 | 3.17e5 | **799** | 4.0·10³ |
| 4 | 254.5 | 254.6 | 1.00 | 1.10e4 | 1.12e4 | **43** | 44 |
| 5 | 38.4 | 40.5 | 1.06 | 3.01e4 | 3.59e4 | **743** | 885 |
| 6 | 241.8 | 234.2 | 0.97 | 3.60e3 | 4.17e3 | **14.9** | 17.3 |
| 7 | (197.5) | 197.5 | — | 2.34e3 | 6.28e3 | **11.8** | 31.8 |
| 8 | 359.5 | 370.2 | 1.03 | 316 | 994 | 0.85 | 2.7 |
| 12 | 978.9 | 984.9 | 1.01 | 30.5 | 458 | 0.03 | 0.46 |
| 16 | 1783 | 1998 | 1.12 | < 0 | 232 | — | 0.12 |
| 21 | 1758 | 2786 | 1.58 | < 0 | 166 | — | 0.06 |

* **Windows.** At 1000 the NEARGAP2 window terms (ν < 50 at G_low = 10⁴ c₂ with the gauge factor; x > 0.9 at G_edge) use ≤ 4 % of the
  certified margin for J ≤ 6 and ≤ 10 % for J = 7; in chart A they scale as Λ⁻⁴ and Λ⁻² (gapnear R = 0.065 and 0.25 at 2000). They
  reduce F_J by ≤ 10 %.
* **Reading.** For J ≤ 7 the 1000 certified margin is ≥ 12× the measured `sup|Γ|`, and Lemma T makes the allowance grow like Λ². For
  J = 8 only the (uncertified) hybrid floor gives a margin (×2.7). From J ≈ 12 on, the sign-free form fails already at 1000.
* **Incremental form (U_Δ), `gu2_growth.py inc`.** Under the two-term law `Γ_m = g₀ + g₁/m` fitted on the ~1000/~2000 rows
  (max |g₁|/m₀ = 1.1 / 9.6 / 31 / 17 / 61 / 102 / 216 c₂ for J = 0/2/4/5/6/8/12; roughly 1.6·10³ J²/m, the J²/m Laplace parameter),
  `min_{t≥1}[λ_min t² − |δλ₀| − E₀ − W₀|g₁|/m₀ (1 − 1/t)]/λ_min` is ≥ +0.07 (sign-free E₀) for every 1000 state with J ≤ 8, and ≥ +0.38
  with the λ-level budget for J ≤ 12. For J = 16, 21 it is negative, but there the 1000 K₂ budget alone already exceeds λ (L < 1):
  K₂ cannot carry these spins at any gap; HIGHSPIN's direct rows do at 1000, and nothing comparable exists at 2000. [NE/M]
* **Caveat (data).** The published 1600 rows are not on the same protocol as gapnear's own rows for even-J `++` (their `sup|Γ|` is
  ×3–7 smaller, e.g. J = 4 `C2 ++`: 219.7 / 33.2 / 215.2 at 1004 / 1600 / 2004), so the 1600 residual column of `gu2_growth.txt` is
  meaningful only for odd J and `+−`.

## 4. The structure that makes uniformity plausible, and what a proof needs

### 4.1 Proposition S (fixed phase in chart A) [D]
In the crystal engine (`spinfunc/l1_crystal.build`, and the EE/parity-odd builders that reuse its geometry) the integrand is
`exp[(iν/2) log A(z) + (Δ/2) log B(z)] · exp G(z)`, where A and B are cross-ratio-type functions of the seven integration variables
only, and G collects measure, `(−x²)^{−2d}`, FP, dual, celestial, subtraction and three-point factors (J enters through exponents and
the spin contraction Q_J, a polynomial). For the EE kinds G also carries the conserved ⟨TTO⟩ coefficients, rational in Δ (closed form
at J = 0 in REPRO_ENGINE §1). Hence, with `ξ = ν/Δ`,
```
I(Δ, ξ) = c(Δ, ν) ∫_{D′} e^{Δ Φ_ξ(z)} g_J(z; Δ) dz,      Φ_ξ = (iξ/2) log A + ½ log B   (Δ-independent),
```
c a z-independent prefactor (`1/β_{Δ,J}`, `A_{k,ν}`, …) that cancels between state and its flat normalisation or is explicit. The
crystal (the code's `saddle(nu, Delta)`: `a = (Δ − √(Δ²−ν²))/ν = ξ/(1 + √(1−ξ²))`) depends on ξ only. Chart A's nodes sit at fixed
`x = ν/m`, i.e. `ξ = x(1 − 3/(2m) + O(J²/m²))`; this is a smooth O(1/m) reparametrisation (`Δ = 3/2 + √(m² + (J+½)²)`, analytic in 1/m
for |1/m| < 1/(J+½)).
*Consequence.* For a Δ-independent phase on a Δ-independent contour, with amplitude analytic in a fixed complex neighbourhood of a
non-degenerate saddle and analytic in 1/Δ, the Laplace expansion is asymptotic **uniformly** in Δ ≥ Δ₀, with remainder after K orders
bounded by `C_K Δ^{−K−1}`, C_K determined by sup norms on that fixed neighbourhood and by the height gap to the rest of the contour
[Q: Fedoryuk, *The saddle-point method*; R. Wong, *Asymptotic Approximations of Integrals*, chapter on multidimensional Laplace integrals]. In c₂ units the K₂ truncation (relative O(1/m²) kept) leaves
`Γ_m(x) = Γ_{K₂,m}(x) + O(1/m)` (O(1/m²) if the odd orders cancel in c₂, as they do for the exact scalar, §5.2). **H_mono-A at fixed (J, x) is therefore not an independent hypothesis**: it holds whenever the crystal
is the only contributing critical point at or above its own height.
*Check in the data.* The K₂ step (the last Laplace term, `|Γ_{K=2} − Γ_{K=1}|`) is m-independent in c₂, as the structure requires:
44 / 39 (J = 0), 593 / 602 (4), 350 / 358 (5), 1232 / 1245 (6), 2029 / 2048 (8), 3761 / 3790 (12) at ~1000 / ~2000.

### 4.2 Lemma R (cycle rigidity) [D-sketch]
* The descent (Lefschetz) flow of `Re(ΔΦ_ξ)` has Δ-independent flow lines; on the cover on which Φ_ξ is single-valued, the
  decomposition `D′ = Σ_σ n_σ 𝒥_σ` (relative to `Re Φ_ξ → −∞`, i.e. the B-zero boundary of D′ and the spiral ends at the A-divisors)
  has **Δ-independent integer** coefficients. Δ enters only through Hankel regularisations at atom divisors whose exponents contain
  `iν/2 = iξΔ/2` (e.g. `m₂₀^{−¼ − iν/2}`), giving factors `1 + O(e^{−πξΔ})`.
* The amplitude's own branch divisors are Δ-independent (RD_CLOSE §1.2; EE = a·Q·scalar with no new divisor at J = 0).
* Therefore: a non-crystal thimble σ with `Re Φ_ξ(σ) ≥ Re Φ_ξ(σ*)` and `n_σ ≠ 0` would contribute ≥ `e^{0}` × its amplitude ratio at
  **every** Δ. The direct D′ RQMC rows at Δ ≈ 1000 (J = 5: 9 nodes, J = 6: 5, J = 8: 4, J = 12, 16, 29: full rows; NEARGAP2, HIGHSPIN)
  agree with K₂ to ~10⁻⁵ relative; that excludes such σ at those ξ **for all Δ**. A lower thimble (δ_σ = Re Φ(σ*) − Re Φ(σ) > 0)
  contributes `e^{−δ_σΔ}·poly(Δ)`, decreasing for Δ > deg/δ_σ.
* **Residual (non-generic):** a thimble with δ_σ ≲ |log ε|/Δ₀ and an amplitude accidentally ≤ ε at Δ₀ is not excluded; its
  polynomial growth could outrun `e^{−δ_σΔ}` for a while. This is RD_CLOSE §4.2's "sharp failure mode" in chart-A form.
* RD_CLOSE's census (critical points 30–94 units *above* the crystal at fixed ν, Δ = 300) is consistent with this: such points must
  have n_σ = 0, which the scalar periods and the direct rows show.
* This is what removes PCYC_CLOSE's remaining "[I] at finite Δ (no Stokes jump for a ∈ (0, 1/Δ])" **in chart A**: at fixed ξ there is
  no Δ-dependent Stokes phenomenon at all. (It does not remove it in the fixed-ν Regge chart, where the phase depends on ν/Δ.)

### 4.3 Lemma L (quantitative Laplace remainder) — the missing constants
In Morse coordinates `w` at the crystal (`z = z* + Lw`, `Lᵀ(−∂²Φ)L = 1`, the code's L), a standard Cauchy-estimate argument bounds the
K-th remainder by `C_K(n, ρ, θ, σ_Φ, η_g, δ) Δ^{−K−1}` given: (a) a polydisc `|w_i| ≤ ρ` on which the cubic-and-higher part of Φ_ξ is
bounded by σ_Φ and `g_J(·; Δ)` by η_g for all Δ ≥ Δ₀; (b) `Re[Φ_ξ(z* + Lw) − Φ_ξ(z*)] ≤ −(1−θ)|w|²/2` on the real tangent plane inside
the polydisc; (c) a height gap δ to the rest of the (Lemma-R) cycle. [Q: Fedoryuk; Olver for 1D explicit forms]
**All of (a)–(c) are statements about Δ-independent objects** (g's Δ-dependence is rational and monotone), so checking them once at
each node gives a bound valid for every Δ ≥ Δ₀, decreasing in Δ. This is the concrete route to an effective theorem; it needs a
pointwise complex evaluator of the spinning EE/C3H amplitudes (not available in the project) and interval sup bounds on a 7-complex-
dimensional polydisc per node. Not done.

### 4.4 Theorem E (eventual positivity) [D-sketch]
Given R at all nodes and the window statements of §4.5: for every near-gap state (J, μ > 1, h),
`sup_{x≤1/μ}|Γ_Λ(x)| ≤ K(J, μ) < ∞` for all Λ ≥ Λ₀, so `λ(Λ) ≥ λ_hyb(Λ) − K·A/(μΛ)²/N_s → flat ≥ N_f/5 > 0`. Hence λ(Λ) > 0 for
Λ ≥ Λ\*(J, μ). The constant K (and Λ\*) is not computed — this is CHMRS's "sufficiently large M" (§1.2), now at the level of
individual near-gap states of the spinning EE functional, with the reason for uniformity identified (Prop. S, Lemma R).

### 4.5 Windows
* **Small ν (ν ≲ 50, x → 0).** At fixed ν the engine's 1/Δ series applies (RD_CLOSE Prop. D′: engine = ∫_{D′} + O(Δ^N e^{−0.336Δ}),
  [D] for the J = 0 kinds; the spinning near-gap engines `po_ee_odd`, `l9` are [NE]). The intermediate range 50 ≲ ν ≪ m is CHMRS's
  matching property (ii) (§1.2) [D-sketch]. In chart A this window's weight ∝ (50/Λ)²·(1/m²): its cost falls like Λ⁻⁴ at fixed |Γ|
  (gapnear: E_low ratio 0.065 at 2000 = 1/15.4).
* **Light-cone edge (x → 1/μ).** At fixed μ > 1 the crystal stays a distance ∝ √(1−1/μ²) from the fold (saddle coalescence at
  ξ = 1), so S/R/L apply with constants depending on that distance. As μ → 1 (twist edges at large Δ_g) the expansion parameter is
  `1/(m(1−x²)^{3/2})` and an Airy region of width m^{−2/3} appears; there Γ can be O(m²). Since NF02's ψ₂, ψ₃ vanish like (1−p)^{5.0}
  at p = 1 (checked: local exponents 4.92 → 4.99), the edge cost is ≲ ∫(1−x)^5 (1−x)^{−γ}dx with γ ≈ 3–4.3 (exact scalar,
  GAP_SCALING §3.3), finite uniformly in μ ≥ 1; the Airy cap costs ~m^{−4}·m² relative. [D-sketch]

### 4.6 Why Prop. D′ (RD_CLOSE/PCYC_CLOSE) does not do the job in chart A [D]
* Prop. D′ is an expansion at **fixed ν** with an **absolute** far-region bound `|B|^{Δ/2} ≤ e^{−2Δ artanh(r/2)}` on D′.
* On D′, `|A^{iν/2}| = 1`, while the signed integral is `e^{−πν}` times the near-region scale (CLS-D (3.33): `Re f = −πν` at the
  crystal; the exact-skeleton cancellation `e^{−πν}`, EXACT_SKELETON §1). So the far-region error relative to the value is
  `≲ Δ^N e^{−(2 artanh(r/2) − πξ)Δ}`, useful only for `ξ < 2 artanh(r/2)/π`: **ξ < 0.107** at the guaranteed r = 1/3 (0.163 with the
  measured c(½) = 0.511).
* Its near region (raw radius ≤ 1/3) contains the crystal only for ξ < 0.6 (`a(ξ) = ξ/(1+√(1−ξ²)) = 1/3`), and its "engine series"
  is the fixed-ν a-series, which must be resummed in ξ at ν ∝ Δ.
* So Prop. D′ covers at most the low window ξ ≲ 0.1, where ψ₂ ~ 188p carries ~1 % of a near-gap state's weight. The decisive bulk
  x ∈ [0.1, 0.9] requires the crystal (deformed-contour) analysis of §4.1–4.3.

## 5. Numerical checks of the structure

### 5.1 Graviton K₂ rows (gapnear, read only) [NE]
Covered in §3.4 and §4.1: `sup|Γ|` ratios 2000/1000 = 0.97–1.06 for J ≤ 8; per-node changes ≤ |g₁|/m₀ ≈ 1.6·10³ J²/m₀ c₂; K₂ steps
m-independent. At J = 5 the variation is concentrated between 1005 and 1600 (x = 0.4: 22.3 → 26.3 → 26.8; x = 0.5: 19.8 → 24.9 →
25.3), i.e. faster than 1/m afterwards.

### 5.2 Exact scalar heavy moment, chart A to Δ_g = 10⁴ (`gu4_scalar.py` → `cache/gu4_scalar*.json`, `gu4_scalar.txt`) [DV]
Arb Laplace/Barnes series (no asymptotics), k = 2, Δ_φ = 3, fixed (J, μ), Δ_g = 1000, 1500, 2000, 3000, 5000, 10000
(2·10⁴ exceeds the 3 GB limit):

| (J, m/Λ) | Λ²S (k = 2 cost × Λ²), 1000 → 10⁴ | sup\|Γ\| on x ≤ 0.9, 1000 → 10⁴ | sup\|Γ\| on x ≤ x₁ | worst node `max_m\|Γ_m\|/\|Γ_{1000}\|` |
|---|---|---|---|---|
| (5, 1.05) | 48.608 → 48.621 → 48.626 → 48.631 → 48.636 → 48.639 | 2484 → 2514 (+1.2 %) | 4717 → 4806 (+1.9 %) | 1.021 (x = 0.93) |
| (12, 1.05) | −2.9508 → −2.9377 → −2.9318 → −2.9284 → −2.9267 → −2.9255 | 2400 → 2395 (−0.2 %) | 13006 → 13110 (+0.8 %) | 1.031 (x = 0.93) |
| (26, 1.05) | 0.0520 → 0.0447 → 0.0420 → 0.0401 → 0.0392 → 0.0388 | 7852 → 7649 (−2.6 %) | 61130 → 59343 (−2.9 %) | 1.024 (x = 0.40, 0.95) |
| (8, 1.3) | 10.304 → 10.303 → 10.303 → 10.302 → 10.302 → 10.306 | 300.2 → 300.2 | 300.2 | 1.001 |

* Γ(x) in c₂ converges at every node; the node-wise fit `g₀ + g₁/m + g₂/m²` has residuals ≤ 5·10⁻³ c₂ on x ≤ 0.85 (≤ 10⁻² at J = 26), rising to 0.1–52 c₂ at x = 0.9–0.95 where the fold sets in. In the bulk the approach is **1/m², not 1/m** (e.g. J = 5, x = 0.85: successive differences 2.19, 0.78 for 1000→1500→2000, ratio 2.8 against 2.9 for 1/m² and 2.0 for 1/m): the odd Laplace orders cancel in c₂, as CHMRS's symmetric-frame statement (only even powers, §1.2) suggests; the graviton K₂ rows are consistent with this at J = 5 (§5.1). (The 10⁴ entries carry ~2·10⁻³ c₂ of quadrature noise at small x.)
* `sup_{m≥m₀}|Γ_m|/|Γ_{m₀}|` (worst node) ≤ 1.024 on x ≤ 0.9 and ≤ 1.031 on x ≤ 0.95: the 1000 value is *not* literally an upper bound
  for larger gaps (the approach to g₀ can be from below), but the excess is ≤ 3 % — one (J = 6, 7) to four (J = 0) orders of magnitude
  inside F_J.
* This is the fixed-phase mechanism of Prop. S on an object with an exact, non-asymptotic evaluation (the scalar shares Φ_ξ with the
  graviton kinds; RD_CLOSE §1).

## 6. Reproduce
```
cd ads4/code/gapuni && bash run_all.sh        # ~0.6 CPU-h, <= 2 processes; outputs in cache/
```

| file | content | grade |
|---|---|---|
| `gu_guard.py` | write guard (only `gapuni/cache`) | — |
| `gu1_tol.py` | Lemma T check; 1×1 tolerances J = 0–60 at 1000/2000/5000 (`gu1_tol.txt/.json`, summary `gu1_summary.txt`) | [DV] |
| `gu1b_weyl.py` | 2×2 Weyl tolerances, even J ≤ 16 (`gu1b_weyl.txt`) | [DV] |
| `gu2_growth.py` | Γ rows 1000/1600/2000: sup, g₁, K₂ steps; `inc`: incremental criterion (`gu2_growth.txt/.json`) | [NE/M] |
| `gu3_seam.py` | seam-state tolerances and edge-strip weights at 2000 / 5000 / 10⁵ (`gu3_seam.txt`) | [DV] |
| `gu4_scalar.py` | exact scalar chart-A scan to Δ_g = 10⁴ (`gu4_scalar*.json`, `gu4_scalar.txt`) | [DV] |
| `gu5_cover.py` | growth-allowance table F_J (`gu5_cover.txt`) | [NE] |

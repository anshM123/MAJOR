# NU_TAIL_CHECK: do the published weights reach ν above a near-gap heavy state's mass, and does the exponential growth of the exact per-ν image matter?

*24 Nov 2026, checker. Unreviewed, not committed. No existing file was edited. New code: `ads4/code/nutail/` (`probe.py`,
`arbprobe.py`, `tail.py`, `altcut.py`); outputs in `ads4/code/nutail/cache/` (git-ignored by its own `.gitignore`). Read-only
imports: `exactskel/exactmom.py` (double-precision contour), `exactskel/skelexact3.py` (contour-free Arb Laplace/Barnes series),
`uniform/skel.py`. Data read: `prl/data/d_abel_X_{300,1000}.txt`, `prl/data/k2/xq_abel_X_{300,1000}.txt`,
`prl/data/k2/k2_continuum_{D,Q}_{300,1000}.json` (+ tables), `prl/data/kernels.json`, `prl/data/q_deshift_xp_*.txt`.
Compute ≈ 2 CPU-h, ≤ 2 processes, 1 thread, `ulimit -v 3000000`. The caution checked is DSKETCH_UP §0 (iii) / UT3 "Caution (unchecked)".*

Grades: [P] proved, [D] derived, [DV] derived and checked numerically, [NE] numerical evidence, [M] model, [I] open.

---

## 0. Verdict

**Split verdict. Nothing computed is affected; one published *definition* is.**

| what | verdict | size |
|---|---|---|
| **Every published per-state number** (Ω^EE, Ω^Q_up, Ω^Q_low evidence: crystal/D′ integrals, K₂ rows, direct n1/n3 points, GK15 edge nodes, LR/hybrid floors, EDGE-R, Theorems W/L/L′) | **(a) no effect** [P] (by inspection of the definitions) | all are evaluated in the ψ-form on `ν ≤ Λ ≤ m` (`x ≤ x₁ = Λ/m ≤ 1`); the exact image is never sampled above the mass |
| **Published light-side constants** (`R₀`, `L_A`, `2L_AdS`, `K`, the (2a)/(2b) numbers) | **(a) no effect** | light images are polynomial/pole-type; the heavy growth does not enter |
| **Q's C4/C5 pieces** (Ω^Q_low, Ω^Q_up), **all k = 3 pieces** | **(a) no effect** [P] | supports `ν ≤ 123` (C4), `≤ 82` (C5), `≤ Λ` (ψ₃); heavy `m ≥ Λ ≥ 298.5`, so `x ≤ 0.41` |
| **k = 2 de-shift pieces realised by the continuum weights** (`k2_continuum_{D,Q}_*.json`) **or by the exact min-norm X** | **(a) no effect** [DV, scalar proxy] / [M] for EE | tail beyond `m` ≤ **2·10⁻¹⁰** (D) / **1.2·10⁻⁶** (Q) of the state's own action; margin ≥ 10⁵ |
| **k = 2 de-shift pieces realised by the stored grid files at Δ_g = 1000** (`d_abel_X_1000.txt`, `k2/xq_abel_X_1000.txt`, which run to ν = 1310.72 = Λ + 312) **with the published cutoff θ(Λ+100, 10)** | **(c) REAL PROBLEM (definitional)** [DV, scalar proxy] / [M] for EE | the double-precision floor of the grid (`|X| ≈ 2–13·10⁻¹⁵` D, `2–10·10⁻¹⁶` Q) times the exact image (`10¹⁶–10²⁶` at ν = Λ+100…Λ+150) gives a tail ≈ **1.0× (D) / 0.5× (Q) the whole k = 2 action of the J = 0 gap state**, and ≥ 2× (J = 2) to ≥ 50× (J = 10) the scalar-proxy action of higher-J near-gap states (see §3.3) |
| same, **Δ_g = 300** grids (end at ν = 327.68 = Λ + 29) | **(a) no effect** | tail ≤ 1.6·10⁻⁹ (D), ≤ 1.1·10⁻⁵ (Q) of the state's action |

**What the real problem is, precisely.** No number in the project is wrong. But if `Ω^EE_θ` / `Ω^Q_up` at Δ_g = 1000 are *defined* by the
stored grid weights (as `prl/data/README.md` §4 "Functional. `Re Σ_k dν X_k Ψ̂_{EE,ν_k}/(h₋+h₊)(ν_k)`" and UT3 §2.2(a) read), then
the functional's action on near-gap heavy states with `m ≲ Λ + 120` contains a term that no computation evaluated, of the same size
as (or larger than) the action itself, whose sign is set by double-precision round-off in the grid near ν ≈ Λ + 120. (Pos) for that
object is then untested on those states, and HK_CLOSE §5's "every computed per-state value is unchanged to ≲ 10⁻⁴¹ (D)" and §6's "D
cutoff (Λ+100, 10): none" are false for it (HK §6's flat model has no growth at J = 0; the AdS growth is J-independent).

**The fix is free and changes no number:** state that the k = 2 de-shift weights are the continuum weights (`k2_continuum_*.json`,
already "the JSON formula is the definition"), *or* keep the grids but cut them at `(ν_c, a) = (Λ + 60, 6)` or tighter (tail ≤ 5·10⁻⁹
(D) / 4.6·10⁻⁷ (Q) of the action; the cutoff's leakage `½erfc(10) ≈ 10⁻⁴⁵` at ν = Λ is the same as HK's), or set the grid to zero beyond
`Λ + 10`. With any of these, the answer to the caution is **(a) no effect**, with ≥ 10⁵ margin.

---

## 1. Item 1: the ν-support of the published weights

`Λ = 298.4996 / 998.4999`. Every heavy state (`Δ − J ≥ Δ_g`) has `m² = (Δ−J−2)(Δ+J−1) ≥ (Δ_g−2)(Δ_g+2J−1) ≥ Λ²`, so **`m ≥ Λ`, with
equality only for J = 0 at τ = Δ_g**; near the gap `m − Λ ≈ J`.

| piece | weight (file) | support in ν | reaches ν > m for some heavy state? |
|---|---|---|---|
| Ω^EE k = 3; Ω^Q smeared k = 3 | `ψ₃` (`nf02_psi.json`; Q3\*) on `p = ν/Λ ∈ [0, 1]`, ∝ `(1−p)⁵` (NF02) / `(1−p)³` (Q3\*) | `ν ≤ Λ` | **no** (`ν ≤ Λ ≤ m`) |
| ψ-form of the k = 2 pieces (`∫₀¹ψ₂Ĉ′`, the form every evidence computation uses, §2) | `ψ₂` on `[0, 1]` | `ν ≤ Λ` | **no** |
| Q C4 value/derivative (Ω^Q_low, Ω^Q_up) | XP grid `q_deshift_xp_C4_*.txt` (ν ≤ 122.86, last value `6·10⁻⁶⁷`); continuum `q_deshift_continuum_c4{v,d}.json` (Gaussians to 60, `< 10⁻²⁵` of max beyond); θ(150, 6) | `ν ≲ 123` | **no** (`x ≤ 123/298.5 = 0.41`) |
| Q C5 | XP grid to ν = 81.9; continuum to ≈ 60 | `ν ≲ 82` | **no** |
| **Ω^EE k = 2 (D), Ω^Q smeared k = 2 (X^Q)**, X-form `Re∫X θ Γ` | min-norm Abel X, θ(Λ+100, 10) | see below | **yes**, for `m < ν` in the tail |

**The k = 2 de-shift weights beyond Λ** (`tail.py`, header of `cache/tail_{300,1000}.txt`; `|X|`, max |X| ≈ 21.2 (D), 1.77 (Q)):

| weight | (Λ, Λ+10] | (Λ+10, Λ+30] | (Λ+30, Λ+45] | beyond | grid end |
|---|---|---|---|---|---|
| D grid, 300 | exact tail `2.5·10⁻⁸ e^{−1.52(ν−Λ)}` | floor `≤ 6.8·10⁻¹⁴` | — | — | Λ + 29.2 |
| D grid, **1000** | exact tail `1.2·10⁻¹⁰ e^{−1.36(ν−Λ)}` | floor `≤ 8·10⁻¹⁵` | floor `≤ 1·10⁻¹⁴` | **floor `≤ 1.3·10⁻¹⁴` (median 3·10⁻¹⁵) to Λ+312** | **Λ + 312.2** |
| Q grid, 300 | exact tail `2.8·10⁻⁵ e^{−1.53(ν−Λ)}` | `≤ 5.7·10⁻¹²` | — | — | Λ + 29.2 |
| Q grid, **1000** | exact tail `3.8·10⁻⁷ e^{−1.56(ν−Λ)}` | `≤ 1.2·10⁻¹³` | floor `≤ 9·10⁻¹⁶` | **floor `≤ 1.0·10⁻¹⁵` (median 2·10⁻¹⁶) to Λ+312** | **Λ + 312.2** |
| D / Q continuum, 300 | exact tail | `≤ 3.9·10⁻¹² / 2.6·10⁻¹¹` | last Gaussians `≤ 2.9·10⁻¹² / 1.7·10⁻¹¹` | super-Gaussian (`σ = 0.9`, last centre Λ + 30.2) | formula |
| D / Q continuum, 1000 | exact tail | `≤ 8·10⁻¹⁴ / 1.6·10⁻¹³` | `≤ 1.3·10⁻¹² / 2.3·10⁻¹³` | super-Gaussian | formula |

(Continuum values: the JSON formula evaluated in 40-digit mpmath on [Λ−5, Λ+45], `X = χ(R₋+R₊)`; agrees with the published table to the
table's double-precision noise, `≈ 10⁻¹³` absolute.)

* **The exact weights decay like `e^{−(π/2)(ν−Λ)}` beyond Λ** (fitted rates 1.36–1.56; π/2 = 1.571). This is the same rate as the
  asymptotic growth of a block's per-ν image (§3.1), so the exact weights are safe only because the growth rate `c(x)` stays well below
  π/2 at the relevant `x = ν/m ≤ 1.3` (§3.1).
* **The grids at Δ_g = 1000 keep a non-decaying double-precision floor out to Λ + 312.** The published θ(Λ+100, 10) leaves it fully
  active up to ν ≈ Λ + 100 and Gaussian-suppressed only beyond. The Δ_g = 300 grids stop at Λ + 29, before any harm.
* "De-shift tails up to Λ+100" is therefore true only for these floors and for θ's nominal reach. The exact X's own tail is negligible
  beyond Λ + 10.

---

## 2. Item 2: how the heavy-state action is actually computed

* **All per-state evidence uses the ψ-form on `ν ≤ Λ`, at fixed `x = ν/m ≤ x₁ = Λ/m ≤ 1`** [P, by inspection]:
  * rows (K₂, w2/l9 chains, crystal/D′ = the CLS-D physical cycle per PCYC_CLOSE) are evaluated at fixed `x` on `[0, x₁]`, and
    assembled with ψ₂/ψ₃ Gauss–Legendre weights ("the chain's local form uses psi2 node weights", N2_SUBREGGE §2, `w2_chain.py:163`;
    `dF = Σ p₂ env Γ/m²`);
  * direct n1/n3 points: windows `ν < 50` and bulk rows, plus GK15 edge nodes on `[0.8, min(x₁, 0.945)]`; nodes beyond `x = 0.945` are
    skipped and covered sign-free (HIGHSPIN, EXT_FOLD); n1 itself fails for `ν/Δ ≳ 0.94`;
  * the one place where X (not ψ₂) multiplies rows is the window `ν < 50` (N2_SUBREGGE), far below any mass;
  * AG_PROOF's Lemma P4W-χ / (EDGE-R) is `E_χ = ∫_{p_w}^{1}ψ₂χG dp`, i.e. `ν ≤ Λ`; P4_WINDOW's EDGE term is the amplitude at `p = 1`.
    These handle the **light-cone edge below the mass** (flat `E_k(x)P_J ∝ (1−x²)^{−1−k/2}`, the exact `(1−x²)^{−3}` growth of G),
    not anything above it.
  * Leading-Regge images `P_{−½+iν}(cosh β)` (HK, CERT, Lemma LR∞, Theorems W/L) are bounded for all real ν, so their use of the
    full X grid is harmless.
* **So the exponential growth above m never enters a computed number.** It enters only the *definition*: the functional is the X-form
  `Re∫Xθ Γ dν` over all ν (UT3 §2.1–2.3, README §4), and on a near-gap state its value is
  (computed ψ-form value) + (X-form − ψ-form on `ν ≤ Λ`: the N2/PW_IMAGES non-local and monopole terms, already budgeted elsewhere)
  + **the tail `∫_{ν>Λ} Xθ Γ`** — the subject of §3.
* **Does the per-ν heavy image grow above m?** In the scalar CHMRS representation: **yes**, smoothly through `x = 1` and then like
  `e^{m G(x)}` (§3.1). The EE raw rows are not available above `x ≈ 0.95`; the EE integrand is the scalar one times a rational factor
  with no new divisor (RD_CLOSE §1) and the growth comes from the common `N_com⁻¹ ~ γ_a⁻² ~ e^{πν}` normalisation, so the scalar profile
  is the natural proxy ([M] for EE).
* **Ill-defined or edge-treatment?** Neither the light-cone-edge treatment (it lives at `x ≤ x₁`) nor AG_PROOF's lemma applies above m.
  The total is **well defined** for every weight: θ (or the continuum weights' Gaussians, or the exact `e^{−(π/2)(ν−Λ)}` tail) beats any
  `e^{c(x)ν}` with `c < π/2`. The issue is **size**, not existence: harmless for the exact/continuum weights, O(1) for the grid floors.

---

## 3. Item 3: numerics

### 3.1 The exact per-ν image above the mass (scalar CHMRS, Δφ = 3, k = 2; `Γ := m⁴C_{2,ν}/2 → 1` at ν → 0)

**Method and validation [DV].** `skelexact3`'s contour-free series (Barnes' first lemma; positive `g_m`, no `e^{−πν}` cancellation) in
Arb, N0 = 240, nq = 320:
* reproduces DSKETCH's contour values above the mass: Δ = 40 / 20 / 10 / (12, J=2) at ν = 60: `7.2652·10⁶ / 5.1210·10¹⁶ /
  8.99·10²³ / 3.448·10²⁴` (DSKETCH `7.265·10⁶ / 5.122·10¹⁶ / 9.2·10²³ / 3.46·10²⁴`; the 2 % at Δ = 10 is the contour's error);
* agrees with the double-precision contour (`exactmom`) where the latter is stable (Δ = 100: 2·10⁻⁸ at ν = 90; the contour drifts by
  2·10⁻⁵ / 3·10⁻⁴ / 8 % / ×5 at ν = 100 / 110 / 130 / 150 as its cancellation reaches 10¹⁵);
* converged in `(N0, nq, M, prec)`: `(120, 160)` is **not** enough above the mass at Δ = 1000 (4.6 % at ν = 1310); `(240, 320)`,
  `(400, 400)`, `(600, 500)` agree to `6·10⁻¹⁰`; Arb radius `≤ 10⁻¹¹⁴`.

**Values** (`cache/tail_{300,1000}.txt`, `arbprobe.py`):

| state (Δ, J), m | Γ(0.9m) | Γ(m) | Γ(Λ+30) | Γ(Λ+100) | Γ(Λ+150) | Γ(grid end) |
|---|---|---|---|---|---|---|
| (300, 0), m = Λ | 33.6·½… = 7.2ᵃ | 347 | 3.6·10⁶ | 2.2·10²⁴ | 1.4·10⁴¹ | 1.8·10⁶ (Λ+29) |
| (1000, 0), m = Λ | 16.8 | 1.80·10³ | 5.1·10⁵ | 1.59·10¹⁶ | 2.8·10²⁶ | 5.8·10⁶⁹ (Λ+312) |
| (1002, 2), m = Λ+2 | — | 1.63·10³ | 3.3·10⁵ | 8.6·10¹⁵ | 1.4·10²⁶ | 2.6·10⁶⁹ |
| (1010, 10), m = Λ+10 | — | −480 | 6.6·10⁴ | 3.4·10¹⁵ | 6.0·10²⁵ | 1.2·10⁶⁹ |

ᵃ Γ = m⁴C/2; the flat bulk-point form at 0.9m is 16.9 (1000).

**Growth rate `d ln Γ/dν` as a function of `x = ν/m`** (J = 0): 0.185 / 0.129 at `x = 1` (300 / 1000), 0.235 / 0.208 at 1.02,
0.31 / 0.31 at 1.05, 0.42 / 0.43 at 1.1, 0.58 at 1.2, 0.69 at 1.3 (1000); J = 2, 10 agree to ±0.03 beyond `x = 1.05`. So
**`ln Γ ≈ m G(x)`** with `G′(1) ≈ 0.13–0.19`, rising slowly to π/2 only as `x → ∞` (DSKETCH: 1.38 at `x ≈ 17`) [NE]. Through `x = 1` the
exact image is smooth (entire in ν), finite at the mass (`Γ(m) = 347 / 1802`, against the flat `E₂ ∝ (1−x²)⁻²` divergence).

### 3.2 Tails of the k = 2 functional on near-gap states (`tail.py`)

Per state and weight: `A_ψ = ∫₀^Λ ψ Γ dν` (the ψ-form action: what the evidence computes), `T_m = Re Σ_{ν>m} dν X θ Γ` (signed tail
beyond the mass), `|T_m| = Σ_{ν>m} dν |X| θ |Γ|`; `S_ψ = ∫|ψ|dν` (action of a state with Γ ≡ 1: 3.10·10³ D / 303 Q at 300;
1.04·10⁴ D / 1.01·10³ Q at 1000).

**Δ_g = 1000** (`cache/tail_1000.txt`, `tail_1000_hi.txt`):

| state (Δ, J) | m − Λ | A_ψ (D) | **D grid**: T_m (\|T_m\|/\|A_ψ\|; /S_ψ) | **Q grid**: T_m (/\|A_ψ\|) | D cont | Q cont | D exact-tail |
|---|---|---|---|---|---|---|---|
| (1000, 0) | 0 | 1.42·10⁴ | **+1.41·10⁴ (1.2; 1.7)** | **+1.03·10³ (0.68)** | 1.8·10⁻¹⁰ | 4.6·10⁻⁷ | 1.3·10⁻¹¹ |
| (1002, 2) | 2.0 | 3.16·10³ | **+7.39·10³ (2.9; 0.87)** | **+5.40·10² (22)** | 5·10⁻¹⁰ | 1.2·10⁻⁶ | 3.5·10⁻¹² |
| (1010, 10) | 9.9 | 55.2 | **+3.07·10³ (68; 0.36)** | **+2.24·10² (94)** | 6.5·10⁻⁹ | 2·10⁻⁸ | 10⁻¹⁵ |
| (1030, 30) | ⟦J30⟧ | | | | | | |
| (1060, 60) | ⟦J60⟧ | | | | | | |
| (1100, 100) | ⟦J100⟧ | | | | | | |
| (1150, 0) | ⟦1150⟧ | | | | | | |
| (1250, 0) | ⟦1250⟧ | | | | | | |

(cont / exact-tail columns: `|T_m|/|A_ψ|`.) The grid tail peaks at **ν ≈ Λ + 121–124** (θ's Gaussian against the growth) and is
positive in every case computed (the floor's local pattern there is fixed; the sign is not controlled).

**Δ_g = 300** (`cache/tail_300.txt`; grids end at Λ + 29.2):

| state (Δ, J) | m − Λ | A_ψ (D) | D grid \|T_m\|/\|A_ψ\| | Q grid | D cont | Q cont | D exact-tail |
|---|---|---|---|---|---|---|---|
| (300, 0) | 0 | 4.25·10³ | 1.6·10⁻⁹ | 1.1·10⁻⁵ | 1.5·10⁻⁸ | 1.2·10⁻⁵ | 1.6·10⁻⁹ |
| (302, 2) | 2.0 | 967 | 3.2·10⁻¹⁰ | 8.6·10⁻⁵ | 3.2·10⁻⁸ | 1.6·10⁻⁴ | 2.8·10⁻¹⁰ |
| (310, 10) | 9.8 | 17.6 | 4.6·10⁻¹⁰ | 1.7·10⁻⁹ | 4.4·10⁻⁷ | 3.9·10⁻⁵ | 9·10⁻¹⁴ |
| (330, 30) | 28.6 | 0.93 | 5·10⁻¹¹ | 1.5·10⁻¹⁰ | 2.4·10⁻⁷ | 8.5·10⁻⁵ | 7·10⁻²⁴ |
| (350, 0) | 50.0 | 3.80·10³ | 0 (no ν > m) | 0 | 0 | 0 | 0 |

(Against `S_ψ` every Δ_g = 300 entry is ≤ 1.6·10⁻⁶. The Q entries are dominated by the exact X^Q tail just past Λ, whose relative size
at 300 is 10⁻⁵, not by growth.) The scalar proxy's `A_ψ` is small at J ≥ 10 only because `P_J(1−2x²)` oscillates; relative to `S_ψ`
the numbers are smaller still.

### 3.3 Alternative cutoffs for the grid weights (`altcut.py`, J = 0 gap state, Δ_g = 1000)

| cutoff (ν_c, a) | D grid \|T_m\|/A_ψ | Q grid \|T_m\|/A_ψ | leakage ½erfc((ν_c−Λ)/a) |
|---|---|---|---|
| (Λ+100, 10) (published) | **1.2** | **0.68** | 10⁻⁴⁵ |
| (Λ+60, 6) | 5.3·10⁻⁹ | 4.6·10⁻⁷ | 10⁻⁴⁵ |
| (Λ+40, 4) | 2.1·10⁻¹¹ | 4.6·10⁻⁷ | 10⁻⁴⁵ |
| (Λ+20, 2) | 1.5·10⁻¹¹ | 4.6·10⁻⁷ | 10⁻⁴⁵ |

(The Q residue 4.6·10⁻⁷ is X^Q's exact tail at `ν ∈ (Λ, Λ+10]`, not the floor.) The published `ν_c = Λ+100` was chosen for the
leakage `erfc(10)` (HK §4.4); `(Λ+60, 6)` has the same leakage and removes the problem.

---

## 4. Item 4: verdict, what is affected, wording

**(a) No effect** on: every computed per-state value, every published constant, Q's C4/C5 pieces, all k = 3 pieces, and the k = 2 pieces
whenever they are realised by the continuum weights or the exact X. Reason: the ψ-form weights vanish for `ν ≥ Λ ≤ m`; the exact and
continuum de-shift weights decay beyond Λ at ≈ π/2 (or super-Gaussianly beyond Λ + 30), faster than the exact image grows at
`x ≤ 1.3` (rate ≤ 0.7); the resulting tails are ≤ 2·10⁻¹⁰ (D) / 1.2·10⁻⁶ (Q) of the state's action [DV scalar; M for EE].

**(c) Real problem, definitional, Δ_g = 1000 only:** the stored grids `d_abel_X_1000.txt` and `k2/xq_abel_X_1000.txt` as the definition of
the k = 2 pieces of `Ω^EE_θ` and `Ω^Q_up` (with θ(Λ+100, 10)). Affected statements:
1. UT3 §2.2(a) / README §4 ("Functional. `Re Σ_k dν X_k Ψ̂…`" on the stored grid) and UT3 §2.3 (X^Q "now published" as a grid): at
   Δ_g = 1000 these define a functional whose near-gap actions (states with `m ≲ Λ + 120`: J ≲ 120 at the twist edge, and heavy J = 0, 2
   states up to that mass) differ from the computed ones by a round-off-controlled term of the size of the action [DV scalar; M EE].
2. HK_CLOSE §5 "Every computed per-state value is therefore unchanged to ≲ e^{−180} (Q) and ≲ 10⁻⁴¹ (D) relative" — true for the computed
   values, false as a statement about the stored-grid functional on near-gap states; HK §6's table row "D cutoff (Λ+100, 10): none; worst
   −19.6 at J = 0" — its flat model has no growth at J = 0; with the AdS growth the stored 1000 grid gives `+ln(1.7·10⁴) ≈ +9.8`.
3. UT2/UT3 §2.1's claim "Relative effect on any computed number `≲ 10⁻⁴¹` (k = 2)" stands for computed numbers; it should not be read as
   "the stored grid and the cutoff weight are the same functional".
4. Not affected: (Pos) as *assumed* in the theorem is unaffected if the functional is the continuum one; the theorem's constants are
   light-side and weight-normalisation numbers (`K`, `dk₀`, `L_A`), which do not see the heavy tail.

**Proposed wording** (UT3 §2.1/§2.2, README §4):
> The k = 2 de-shift weights are the continuum weights `k2_continuum_{D,Q}_Δ_g.json` (the JSON formula is the definition); the stored
> grids `d_abel_X_*`, `xq_abel_X_*` tabulate them on `ν ≤ Λ` (beyond Λ the grids carry a double-precision floor and must be read as 0).
> Above a heavy state's mass the exact per-ν image grows like `e^{mG(ν/m)}` (`G′ ≈ 0.13` at `x = 1`, 0.43 at `x = 1.1`); the exact and
> continuum weights decay at ≈ π/2 beyond Λ, so the tail is ≤ 10⁻⁶ of any state's action. (Equivalently: use the cutoff `(Λ+60, 6)`
> for the k = 2 pieces; its leakage is the same `½erfc(10)`.)

**Open [I], not needed:** whether the *uncut* exact functional is per-state finite (exact weight decays at π/2, image rate → π/2 from below
as `x → ∞`, so the far tail is marginal); the cutoff makes this moot. The EE raw rows above `x ≈ 0.95` (n1 fails there), so the EE growth
profile itself is [M].

---

## 5. Grades

| statement | grade |
|---|---|
| ψ-form weights (ψ₂, ψ₃), Q C4/C5 weights never reach ν > m for any heavy state | **[P]** (supports + `m ≥ Λ`) |
| every published per-state computation evaluates images only at `x ≤ x₁ ≤ 1` | [P] by inspection of the cited definitions/code |
| exact scalar per-ν image above the mass: smooth through `x = 1`, `ln Γ ≈ mG(x)`, rates 0.13–0.69 on `x ∈ [1, 1.3]` | **[DV]** (Arb series, converged; cross-checked against DSKETCH and the contour) |
| same growth for the EE rows | [M] (shared normalisation; no EE data above `x ≈ 0.95`) |
| exact X / X^Q tails decay at ≈ π/2 beyond Λ | [NE] (fits on the grids before the floor) |
| continuum / exact-tail weights: tail ≤ 2·10⁻¹⁰ (D), 1.2·10⁻⁶ (Q) of the action | [DV] scalar proxy; [M] EE |
| **stored Δ_g = 1000 grids × θ(Λ+100, 10): tail ≈ 1× the action (J = 0) and more at J > 0** | **[DV] scalar proxy; [M] EE** |
| cutoff (Λ+60, 6) or tighter removes it | [DV] scalar proxy |

---

## 6. Reproduce

```
cd ads4/code/nutail; ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 -B arbprobe.py 1000 0 2 998.5,1098.5,1310 400 240 320    # exact image above the mass (seconds)
python3 -B probe.py 100 0 2 90,100,110,130,150 0,-15              # double-precision contour, for comparison
python3 -B tail.py 300          # ~1.5 min   -> cache/tail_300.txt
python3 -B tail.py 1000         # ~45 min    -> cache/tail_1000.txt (J = 0, 2, 10, 30)
python3 -B tail.py 1000 hi      # ~90 min    -> cache/tail_1000_hi.txt (J = 60, 100; m = Λ+150, Λ+250)
python3 -B altcut.py 1000       # ~1 min     -> cache/altcut_1000.txt
```

# QBIND_COROLLARY: error budgets for the binding states of Ω^Q_up at Δ_g = 1000, and an energy-flux corollary for the Letter

*25 Nov 2026, research agent. Unreviewed, not committed by this agent (another agent's WIP commits 715e470 / a9e5d8d swept in
intermediate snapshots of `code/qbind/`; the working tree holds the final versions). No existing file was edited; `C1_FINAL.md` and
`code/c1final/` were not touched.*
* *New code:* `ads4/code/qbind/`:
  * `qb_guard.py`: write guard; every write outside `qbind/cache/` raises.
  * `qb_rows.py` = `neargap2/ng_rows.py`: own-mass K₂ rows, cache redirected.
  * `qb_points.py` = `ng_points.py`: direct n1 RQMC points.
  * `qb_chain.py` = `ng_chain.py`: the unchanged EE_NEARGAP chain. Its row loaders also read `neargap2/cache` (read-only).
  * `qb_val.py` = `ng_val.py`: direct vs K₂ node by node, with the state evaluated with Q's weights.
  * `qb_bind.py`: new. The per-state error budget of Ω^Q_up, built on `ng_bind`'s protocol.
  * `qb_lane.sh`: sequential job lane.
* *Caches:* `ads4/code/qbind/cache/`. Its `.gitignore` ignores everything but itself.
* *Read-only imports:*
  * `spinfunc/` (`ee1000_common`, `w2_chain`, `q3_analysis`, `w1_rows`, `hc_points`, `hc_nodes`, …);
  * `resultq/` (`stage6_model`, `common_q`; tables `stage6_tab_Q3-full-delta0.2_998_2300.pkl`; Q's min-norm Abel weight
    `rv_k2eng_Y_1000.pkl`);
  * NEARGAP2's own rows J = 4, 5, 6 in `neargap2/cache`.
* *Compute:* ≈ 1.3 CPU-h. At most 2 processes at a time, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.
* *Literature (Task 2), read as data:*
  * Chowdhury–David–Prakash arXiv:1707.03007 (CDP), eqs. (1.5)–(1.6), (4.39)–(4.41);
  * Buchel–Escobedo–Myers–Paulos–Sinha–Smolkin arXiv:0911.4257, eqs. (3.6)–(3.8);
  * Kats–Petrov arXiv:0712.0743, eq. (1.3), via its abstract and search summary.

**Grades:** [P], [D], [DV], [NE±] (numerical evidence with an explicit error budget: every model input carries a stated safety factor or
convention, and the sign survives it; NEARGAP2's definition), [NE], [NE + M], [M], [Q], [I].

---

## 0. Task 1 result: the binding states of Ω^Q_up at Δ_g = 1000 (δ = 0.2 functional, H0 gauge)

`λ` = `FF/N` (min eigenvalue of the normalised matrix for even J ≥ 4). The budget is sign-free and uses NEARGAP2's conventions (§1.2):
×3 on truncation, interpolation and the measured forward terms; `G = 10⁴ c₂` where nothing is computed. `L = λ_central/E_tot`. In the
`L` column, the first number uses the published C5 gauge (XP, `A₅ = 1.80`) and the number in brackets the worst Stage-6 C5 gauge
(`A₅ = 7.66`).

| state (role) | own row | `λ` central (published) | budget: trunc / interp / low ν / edge / `s_d` / C5 / other | `E_tot` | lower bound | **L** | grade |
|---|---|---|---|---|---|---|---|
| **`++` J = 0, m = Λ (binding)** | Δ = 1000 (new) | **+0.14282** (+0.1428) | 6.6e-5 / 1e-8 / 5.7e-4 / 1.1e-4 / 5.1e-4 / 1.02e-3 / `s_v` 1.3e-5 | 2.29e-3 | **≥ +0.1405** | **62** (25) | **[NE±]** (was [NE + M] given H_RS) |
| **OO J = 0, m = Λ (co-binding)** | same row (Γ_OO = Γ_++) | **+0.14256** (+0.1426) | as above; `s_d` 6.8e-4; OO−++ 7e-7 | 2.46e-3 | **≥ +0.1401** | **58** (25) | **[NE±]** |
| `++` J = 2, 1.118Λ (channel min.) | Δ = 1117.83 (new) | +0.16652 (+0.1665) | 1.4e-4 / 3e-8 / 4.3e-4 / 0 / 4.0e-4 / 4.8e-4 | 1.46e-3 | ≥ +0.1651 | 114 (55) | [NE±] |
| OO J = 2, 1.118Λ | same | +0.16645 (+0.1664) | same | 1.46e-3 | ≥ +0.1650 | 114 (55) | [NE±] |
| `++` J = 2: 1.002 (gap) / 1.10 / 1.135Λ | 1002 (new) / 1117.83 | +0.21289 / +0.16712 / +0.16804 | — | 2.0 / 1.5 / 1.4 e-3 | — | 105 / 110 / 120 | [NE±] |
| `+−` J = 5, 1.187Λ (channel min.) | Δ = 1186.73 (new) | +0.21748 (+0.2175) | 7.0e-4 / 3e-8 / 2.6e-4 / 0 / 1.2e-3 / 2.2e-4 | 2.40e-3 | ≥ +0.2151 | 90 (70) | [NE±]; `s_d` transferred from the gap [M-transfer, inside ×3] |
| `+−` J = 5: 1.16 / 1.21Λ | 1186.73 | +0.21909 / +0.21922 | — | 2.4e-3 | — | 90 / 92 | [NE±] |
| `+−` J = 6, gap 1.006Λ (channel min.; 2×2) | Δ = 1006 (NG2) + K₂ edge nodes to 0.9935 (new) | +0.21309 (+0.2135) | 1.1e-3 / 1e-6 / 2.3e-4 / 4.7e-5 / 2.8e-4 / 2.8e-4 | 1.97e-3 | ≥ +0.2111 | 108 (74) | [NE±] **given K₂ on x ∈ (0.955, 0.994]**; strict sign-free edge: L = **1.1** (thin) |
| `+−` J = 5, gap 1.005Λ | 1005 (NG2) + K₂ edge to 0.994 (new) | +0.29139 (+0.2912) | 7.5e-4 / 4e-7 / 2.2e-4 / 1.2e-4 / 5.1e-4 / 3.8e-4 | 1.98e-3 | ≥ +0.2894 | 148 (91) | [NE±] given K₂ edge; strict: L = 1.2 |
| **`+−` J = 4, gap 1.004Λ (channel min.; 2×2)** | 1004 (NG2) + K₂ edge to 0.9955 (new) | +0.25940 (+0.2594) | 2.7e-4 / 7e-8 / 7.3e-4 / 2.4e-5 / 1.6e-3 / 9.7e-4 | 3.63e-3 | ≥ +0.2558 | 71 (38) | [NE±] given K₂ edge; **strict sign-free edge: L = 0.4, does not close → [NE + M]** |

**Verdict.**
* **The binding pair `++`/OO J = 0 at the gap is now error-bounded [NE±]: `λ ≥ +0.1405` (`++`) and `≥ +0.1401` (OO), L = 62 / 58.**
  This is comparable to Ω^EE's binding state (`L = 78`). With the worst admissible C5 gauge, L = 25.
* The near-binding channel minima `++`/OO J = 2 (1.118Λ) and `+−` J = 5 (1.187Λ) are [NE±] with L ≥ 90.
* No state came out negative.
* **Caution, reported prominently.** The `+−` gap states (J = 4, 5, 6) carry 70–80 % of their sign-free edge weight at
  x = ν/m ∈ (0.955, √u ≈ 0.995], where the `+−` envelope `(1 − y)⁻⁴` peaks.
  * There, direct n1 RQMC fails (replica rel. err. 0.6–3; the contour failure HIGHSPIN noted for ν → Δ).
  * The new K₂ edge nodes give `|Γ_{+−}| ≤ 0.05 c₂` and falling, with `|K₂ − K₁| → 0`. K₂ agrees with direct n1 node by node within
    1σ on x ≤ 0.94–0.955 (§1.4). The K₂ edge rows close all three states with L ≥ 71.
  * If x > 0.955 is instead charged at the protocol's 10⁴ c₂ floor, the result is L = 0.4 (J = 4), 1.2 (J = 5) and 1.1 (J = 6).
  * **So `+−` J = 4 at the gap is [NE±] only if K₂ is trusted on the last 4 % of the ν range, and [NE + M] under the strict protocol.**
    It is a channel minimum (+0.259), not the binding state (+0.143).

---

## 1. Task 1: method and details

### 1.1 Which states bind (sources)
* **Binding pair.** RESULT_Q_STATEMENT §0 and RESULT_Q_FWD_CLOSE §4: `++` J = 0 and OO J = 0 at the gap, FF/N +0.1427 / +0.1426
  (H0), equal to 10⁻⁵.
* **Near-binding channels.** These are the five with `s_fail` within 4× the measured `|s_d|` (FWD_CLOSE §1): `++`/OO J = 0, 2 and
  `+−` J = 4, 5, 6. Their Δ_g = 1000 channel minima, from `spinfunc/cache/qf_ffscan.txt`, are:
  * `++`/OO J = 2: +0.1665 at 1.118Λ;
  * `+−` J = 4: +0.2594 at the gap;
  * `+−` J = 5: +0.2175 at 1.187Λ;
  * `+−` J = 6: +0.2135 at the gap.
* Every other measured channel has a minimum ≥ +0.295.
* **The lowest HYB-only state** is `+−` J = 26 at 6.27Λ, with +0.197 [M]. It lies in POS_ERRORS §2.5's `m > 3Λ` region ([NE + M],
  L ≥ 1.65) and is **not** budgeted here.
* "Q's binding state" in HRS_CLOSE, DG_SCAN, RD_CLOSE, REPRO and REPRO_ENGINE refers to `++` J = 0 at the gap. Its `s_d` and H_RS
  status are taken from there.

### 1.2 Per-state value and budget (`qb_bind.py`)
The value of a state is
`FF_h = HYB_h − a·2u³·s_d(h, J) + δB₂,h + δB₃,h`.

**The pieces:**
* **HYB.** `stage6_model.hyb`: the flat action with exact Wigner functions, plus the leading-Regge AdS curvature.
* **`s_d`.** The measured analytic-radial-engine k = 4 derivative correction:
  * `++` J = 0: +0.5077 (K = 5, H0; HRS_CLOSE §4);
  * OO J = 0: +0.510;
  * `++`/OO J = 2: +0.509 / +0.510;
  * `+−` J = 4, 5, 6: −0.510, −0.512, −0.515. For the 2×2 entries, `++` J = 4 is +0.510, `++` J = 6 is +0.511, and the off-diagonal
    entries are `|s_off| ≤ 0.015` (J = 4) and `≤ 0.76` (J = 6), the latter being the full basis range of FWD_CLOSE §3.3, taken at the
    worse sign.
* **`δB₂,₃`.** The smeared k = 2, 3 graviton corrections. They are Result D's K₂ objects projected on **Q's** `ψ₂`, `ψ₃`, exactly as
  in RESULT_Q_REVISION §1. They are now computed on **own-mass rows** through the unchanged chain.
* **Normalisation.** `N` is the Stage-0 normalisation.

Error budget (sign-free):

| item | definition | source of the convention |
|---|---|---|
| `E_trunc` | 3 × Σ_nodes `|w_Q|·|Γ_K2 − Γ_K1|` on the own row | `ng_bind` |
| `E_interp` | 3 × max(spline vs linear-in-x, every-other-node) | `ng_bind` |
| `E_low` | 10⁴ c₂ × W_Q(ν < max(50, x₀m)) × max(1, `|Re X^Q_minnorm/ψ₂^Q|`) (Q's own min-norm Abel weight; factor ≤ 6.89) | `ng_bind` (gauge factor recomputed for Q) |
| `E_edge` | max(10⁴, 10 × max\|Γ\| of the last 4 own-row nodes) c₂ × W_Q(x > x_max) | `ng_bind` |
| `E_sd` | 3 × (engine truncation + m-transfer) × `a·2u³/N`. The truncations are the a⁴-term sizes: 0.0015 (J = 0, H0), 0.002 (OO J = 0; J = 2 incl. transfer 1.05 → 1.118Λ), 0.005 (`+−` J = 4), 0.012 (`+−` J = 5 at 1.187Λ: 0.004 + the 300-measured gap → 1.19Λ drift 0.008), 0.003 (J = 6) | HRS_CLOSE §6, FWD_CLOSE §2–3, PARITY_ODD_CLOSE §2 |
| `E_sv` | 3 × `c·(2u³/Λ²)·13.4/N` (k = 4 value piece, Λ⁻²-suppressed) | BELOW_FOLD2 §2 |
| `E_C5` | `a u³ A₅ G₅/(m²N)`, `G₅ = 10⁴ c₂`, `A₅ = 1.80` (XP, published; worst gauge 7.66) | DESHIFT_CLOSE §1.3 bound; no k = 5 engine exists |
| `E_OO` | 0.12 c₂ × all node weights + `a·2u³·5·10⁻⁶/N` | PARITY_ODD_CLOSE §1(F); FWD_CLOSE §2.3 |
| `E_stat` | 3 × jackknife (direct-node variant only) | `ng_val` |
| base | flat value (exact polynomial quadrature) and LR curvature (tables); Lemma C bound printed | POS_ERRORS §2.1 |

### 1.3 The binding state in detail (`++` J = 0, Δ = 1000, `cache/qb_bind_even1000.txt`)
* **Chain:**
  * HYB: +0.20000 (flat +0.200000; LR curvature −1.0·10⁻⁸ N, Lemma C bound 2.3·10⁻⁷ N [P]);
  * with `s_d`: +0.14282;
  * with `δB` from the own row: +0.14282 (`δF = +3.5·10⁻⁶`; published rows give the same).
  * The K₁ chain gives +0.14283, so the λ-level truncation is 1.9·10⁻⁵ (×3).
* **The budget is dominated by three model or measured inputs, none of which is the K₂ chain:**
  * the C5 convention: 1.02·10⁻³, 44 %. The tolerance is `G₅* = 1.4·10⁶ c₂`; the measured k = 4 beyond-LR corrections are
    10–330 c₂;
  * the ν < 50 window at 10⁴ c₂ × gauge factor: 5.7·10⁻⁴. The engine's k = 2 image there is ≤ 0.19 c₂ (RESULT_Q_REVISION §1.3), so
    this is ×5·10⁴ conservative;
  * the `s_d` truncation ×3: 5.1·10⁻⁴.
* **Direct n1 check** (5 nodes on the own row, x = 0.1–0.85; `cache/qb_val_J0_1000.txt`):
  * every node agrees with K₂ within 1.1σ;
  * the state shift is `+2·10⁻⁶ ± 3·10⁻⁶`.
* **Gauge.** In the SB gauge (`s_d` = +0.5228) the value is +0.1411 with the same budget (L ≈ 62). With the published C4 weights (XP,
  `s_d` = +0.5132) it is +0.1420.
* **Remaining non-numerical inputs:**
  * H_RS for `s_d` is [D] (PCYC_CLOSE; unreviewed);
  * the C5 and ν < 50 conventions above;
  * the OO identification `Γ_OO = Γ_++` (k = 2, 3: 10⁻⁸ relative on K₂ rows [NE]; k = 4: measured, `|Δs_d| ≤ 5·10⁻⁶`).

### 1.4 The `+−` gap edge (`cache/qb_bind_pmgap1000*.txt`, `qb_val_J4_1004.txt`, `qb_val_J6_1006.txt`)
* **Why the gap states differ.** On the `+−` entries, `W_Q ∝ |ψ|(2 − y)/(1 − y)⁴ m⁻²` with `1 − y ≥ 1 − u_gap ≈ 2J/Δ_g`.
  * About 75 % of the x > 0.9 weight sits in x ∈ [0.95, 0.995]; 0.06 % sits beyond 0.995.
  * The uniform `Γ` that would zero the state is 2.7·10³ c₂ (J = 4), 1.3·10⁴ (J = 6) and 3·10³ (J = 5). This matches
    RESULT_Q_REVISION §1.4.
* **NEARGAP2's rows stop at x = 0.9** (0.875 at J = 6 after despiking). With the 10⁴ floor on x > 0.9 the states get L = 0.3 (J = 4),
  0.8 (J = 5) and 1.3 (J = 6) (`qb_bind_even1000`, `odd1000`).
* **Direct n1 points.** 14 per state, x = 0.85–0.9955.
  * Good to x = 0.94 (J = 4 C3H) or 0.955 (J = 6). Beyond that the RQMC fails (max rel. err. 0.6–3, masked samples).
  * Where they work, the `+−` entries agree with K₂ within 1σ (|diff| ≤ 6 c₂, mostly ≤ 0.4 c₂).
  * `++` J = 4 differs by −50 … −340 c₂ (≈ 2σ), inside its K₂−K₁ term (+130 … +280, ×3).
  * The state values with the good direct nodes injected are +0.25940 and +0.21310 (jackknife 4·10⁻⁶ / 3·10⁻⁵), unchanged.
* **K₂ edge rows** (new, x = 0.92–0.9955).
  * `Γ_{+−}`: J = 4 ≤ 0.01 c₂ beyond 0.9; J = 6 −2.4 → 0 (C2) and +4.6 → 0 (C3H) over x = 0.9 → 0.99. In both cases
    `|K₂ − K₁| → 0` monotonically.
  * `Γ_{++}` grows to −780/−1450 c₂ (J = 4) and −5300/−7100 c₂ (J = 6), with a K₂−K₁ of 30–40 %. These enter the `++` entries
    (λ = 0.75 / 0.85) with the mild `(2 − y)` envelope, so they cost ≤ 1.4·10⁻³ ×3.
* **What is not established.** K₂ on x ∈ (0.955, 0.995] is validated by nothing but its own convergence. The n1 contour fails there,
  and no other method reaches ν → m. Hence the two-level grade in §0.

### 1.5 Not done
* Δ_g = 300 and 2000 (the referee's "Ω^Q_up at 2000"). The same code runs at 2000 if K₂ rows at Δ ≈ 2000 + J are added, and
  `stage6` tables for Λ(2000) would need to be computed. That is a write outside `qbind`, so it is not done here.
* The HYB-only large-mass states (`m > 3Λ`, e.g. `+−` J = 26 at 6.27Λ, +0.197 [M]).
* A k = 5 engine (would replace `E_C5`'s convention).
* A method for ν → m on `+−` gap states (would settle §1.4).

### 1.6 Reproduction (≈ 1.3 CPU-h; at most 2 lanes)
```
cd ads4/code/qbind   # every python call: ulimit -v 3000000; OMP_NUM_THREADS=1
python3 qb_rows.py even EE,C3H 0 1000 dense                    # 5 min
python3 qb_rows.py even EE,C3H 2 1117.8256594 dense            # 8 min
python3 qb_rows.py even EE,C3H 2 1002 dense                    # 8 min
python3 qb_rows.py odd  EE,C3H 5 1186.7321127 dense            # 18 min
python3 qb_rows.py even EE,C3H 4 1004 0.92,0.94,0.955,0.965,0.972,0.978,0.983,0.987,0.99,0.992,0.994,0.9955
python3 qb_rows.py even EE,C3H 6 1006 0.92,0.94,0.955,0.965,0.972,0.978,0.983,0.987,0.99,0.992,0.9935
python3 qb_rows.py odd  EE,C3H 5 1005 0.92,0.94,0.955,0.965,0.972,0.978,0.983,0.987,0.99,0.992,0.994
python3 qb_points.py 4 1004 0.85,0.9,0.92,0.94,0.955,0.965,0.972,0.978,0.983,0.987,0.99,0.992,0.994,0.9955
python3 qb_points.py 6 1006 0.85,0.875,0.9,0.92,0.94,0.955,0.965,0.972,0.978,0.983,0.987,0.99,0.992,0.9935
python3 qb_points.py 0 1000 0.1,0.3,0.5,0.7,0.85
python3 qb_bind.py even1000 "0:++:1.0:0.5077:0.0015" "0:OO:1.0:0.510:0.002" "2:++:1.002000002,1.10,1.118,1.135:0.509:0.002" \
    "2:OO:1.118:0.510:0.002"        # (J = 4, 6 entries in this run predate the edge rows; superseded by pmgap1000)
python3 qb_bind.py odd1000 "5:+-:1.0049925423,1.16,1.187,1.21:-0.512:0.012"
python3 qb_bind.py pmgap1000 "4:+-:1.003996020:-0.510:0.005:0.015:0.510" "5:+-:1.0049925423:-0.512:0.004" \
    "6:+-:1.005988077:-0.515:0.003:0.76:0.511"
QB_XCAP=0.955 python3 qb_bind.py pmgap1000_cap955 (same states)       # strict sign-free edge variant
QB_INJ="4:1004;6:1006" python3 qb_bind.py pmgap1000_inj (J = 4, 6)    # direct nodes injected
python3 qb_val.py 1000 0 1000 1.0 0.5077 ; python3 qb_val.py 1000 4 1004 1.003996020 0.510,-0.510 ; python3 qb_val.py 1000 6 1006 1.005988077 0.511,-0.515
```

---

## 2. Task 2 corollary: the stress-tensor energy flux is isotropic to 7·10⁻⁹

### 2.1 Statement
**Corollary (energy-flux isotropy).** Assume the Letter's (i)–(iii) at `Δ_g = 1000` for the cubic functional `Ω^EE`:
* a higher-spin gap in the single traces of `T×T`;
* tree-level weak coupling below the gap, at leading order in `1/c_T`;
* positivity of `Ω^EE` on heavy states.

Create a state with `ε_ij T^{ij}`, with any polarisation (linear, circular or mixed), and measure the energy flux at null infinity as
a function of the detector angle `θ`. Then

    ⟨𝓔(θ)⟩ = (E/2π) [ 1 + A_ε cos(4θ − φ_ε) ],    A_ε ≤ τ/4 ≤ 7.23·10⁻⁹ ,

so `(𝓔_max − 𝓔_min)/(𝓔_max + 𝓔_min) ≤ 7.23·10⁻⁹` for every polarisation. Helicity eigenstates are exactly isotropic.

**Free-field form.** Write `⟨TTT⟩ = n_s⟨TTT⟩_real boson + n_f⟨TTT⟩_Majorana + p_T⟨TTT⟩_odd`, with `c_T = n_s + n_f`. Then

    |n_s − n_f| ≤ 7.23·10⁻⁹ c_T ,     |p_T| ≤ 1.11·10⁻¹⁰ c_T   (CDP normalisation of the odd structure).

**Comparison with ANEC.** The conformal-collider (ANEC) bound allows `A ≤ 1`.
* A free scalar (`t₄ = +4`) or free Majorana fermion (`t₄ = −4`) has `A = 1`: its flux vanishes in four directions.
* Large-N Chern–Simons–matter theories lie on the circle `A = 1` at every 't Hooft coupling (CDP §5.2).
* A holographic CFT₃ with `Δ_g = 1000` is isotropic 1.4·10⁸ times more tightly than ANEC requires. Its `⟨TTT⟩` is that of equal
  numbers of free real bosons and free Majorana fermions, to 7 parts in 10⁹.

**Grade.**
* "Main bound ⇒ corollary" is [D]: exact kinematics, using quoted formulas [Q] checked here [DV].
* The corollary inherits every grade of the cubic bound:
  * "(Pos) ⇒ bound" is [D];
  * (Pos) itself is proved in flat space, [NE±] on the most constraining AdS states, and conditional or modelled elsewhere.
* It says nothing beyond the main bound; it is the main bound read as an observable.

### 2.2 Assumptions, exactly
1. **The Letter's (i)–(iii) for `Ω^EE` at `Δ_g = 1000`.** The corollary uses **no** quartic input: no `κ` and no light-scalar
   condition.
2. **`c₁ = 0`.** This is now [D] (C1_FINAL.md, unreviewed). For general `c₁`, `τ ≤ R₀(x + √(1+x²))` with `R₀ = 2.89126·10⁻⁸` and
   `x = 1.670·10⁻⁴|c₁|` (C1_COMPUTE). This gives `A ≤ 7.35·10⁻⁹` at `|c₁| = 100`.
3. **Leading order in `1/c_T`.** The statement bounds the tree-level `t₄` and `α_T`. `O(1/c_T)` corrections are outside the Letter's
   scope.
4. **Larger gaps.** For `Δ_g > 1000` the constant improves along `τ ≤ 1.3194·10⁴ √(ln(Λ/8.7415) + 6.384/Λ)/Λ⁴` (DG_SCAN, [DV] for
   the implication). Positivity there rests on the Letter's no-growth assumption. Nothing new is claimed below 1000; at 300 the
   number would be `A ≤ 7.8·10⁻⁷`, with H_LJ(300) open.

### 2.3 Derivation
* **The energy matrix [Q].** In `d = 3` the transverse polarisation of a `T` insertion is a symmetric traceless 2×2 tensor, i.e. a
  2-vector in the basis `(ε₊, ε×) = ((xx − yy)/√2, (xy + yx)/√2)`. CDP (4.39) gives the energy at a calorimeter along `y` as

      Ê(T) = (E/4π) [[1 − t₄/4, α_T/4], [α_T/4, 1 + t₄/4]].

  ANEC positivity of the eigenvalues `(E/4π)(1 ± τ/4)` gives `t₄² + α_T² ≤ 16` (CDP (4.41)).
* **Angular dependence [D].** Rotating the detector by `θ` rotates a spin-2 polarisation by `2θ` in the `(ε₊, ε×)` plane, so the
  traceless part of `Ê` rotates by `4θ`.
  * For a normalised `ε` the flux is `∝ 1 + ε†M(θ)ε`, with `M(θ) = (τ/4)·R(4θ)·diag(−1, 1)·R(−4θ)` up to a fixed rotation.
  * Its amplitude is `A_ε = (τ/4)·|ε₊² + ε×²|/(|ε₊|² + |ε×|²) ≤ τ/4`. This identity was checked algebraically:
    `(|ε₁|² − |ε₂|²)² + 4(Re ε₁*ε₂)² = |ε₁² + ε₂²|²`.
  * Equality holds for linear polarisation. The amplitude is zero for helicity eigenstates `ε ∝ (1, ±i)`.
* **Parity-even part [Q + DV].** It reproduces Hofman–Maldacena / Buchel et al.:
  `⟨𝓔(n)⟩ = (E/2π)[1 + t₄(|ε_ij n_i n_j|²/ε*ε − 1/4)]`. For `ε = ε₊` the bracket is `cos²2θ/2 − 1/4 = cos4θ/4`.
  * **Literature note.** Buchel et al. (0911.4257) eq. (3.8) prints the `d = 3` constant as **1/8**.
  * The angular average of `|ε_ij n_i n_j|²/ε*ε` over the circle is **1/4**, which is `2/(d² − 1)` of their own general eq. (3.6) at
    `d = 3`. Zero net flux therefore requires 1/4. With 1/8 the positivity bounds would come out asymmetric (`−4 ≤ t₄ ≤ 8`) instead
    of `|t₄| ≤ 4` (CDP; free fields at `t₄ = ±4`).
  * We read (3.8) as a misprint and cite CDP for the constant.
* **Free-field basis [Q].** From CDP (1.6)/(4.40): `t₄ = −4(n_f − n_s)/(n_f + n_s)` and `α_T = 8π⁴p_T/(3(n_f + n_s))`, with
  `n_s`, `n_f` counting real bosons and Majorana fermions. In the Letter's normalisation `c_T = n_s + n_f` (DICTIONARY_T4; a Majorana
  fermion has the `c_T` of a real scalar in `d = 3`).
* **Numbers.** `τ ≤ 2.89126·10⁻⁸` (the Letter's headline; reproduced in REPRO) gives:
  * `A ≤ 7.228·10⁻⁹`;
  * `|n_s − n_f|/c_T ≤ 7.23·10⁻⁹`;
  * `|p_T|/c_T = 3|α_T|/(8π⁴) ≤ 1.113·10⁻¹⁰`;
  * a max/min flux ratio `≤ 1 + 1.45·10⁻⁸`.

### 2.4 Considered and not used
* **(b) Transport.** There is no clean `d = 3` route:
  * four-derivative terms do not correct `η/s` in AdS₄ (the Weyl² correction carries `(D − 4)`; Kats–Petrov (1.3));
  * `W³` (∝ `ĝ₃`) and `R⁴` both enter, and the `R⁴` contraction that controls `η/s` is not the single `g₄` that `Ω^Q_up` bounds,
    which is bounded from above only;
  * no verified AdS₄ `W³ → η/s` coefficient was found in the literature.

  Skipped.
* **(c) Ruling out a class.** The contrapositive (`τ > 2.89·10⁻⁸` ⇒ a gap below 1000, or a failure of (ii)/(iii)) is immediate. But
  every example we found with `τ = O(1)` (vector models, CS–matter) already has twist-1 higher-spin currents. It is therefore a
  remark, not a new exclusion.
* **Flavour currents.** `a₂, α_J` (`⟨JJT⟩`) are **not** bounded by this work; they are not `c_T` data.

### 2.5 Draft text
**Letter** (2–3 sentences, e.g. after "Why it matters"):
> In physical terms, a state created by the stress tensor deposits energy at null infinity with angular profile
> $(E/2\pi)[1+A\cos(4\theta-\varphi)]$, $A\le\tau/4$. Causality only requires $A\le1$, which a free scalar or free fermion saturates;
> under (i)–(iii), every holographic CFT$_3$ with $\Delta_g=1000$ has $A\le7.2\times10^{-9}$. Its energy flux is isotropic to one
> part in $10^8$, and its $\langle TTT\rangle$ equals that of equal numbers of free real bosons and free Majorana fermions,
> $|n_s-n_f|\le7.2\times10^{-9}c_T$.

**Supplement paragraph:**
> *Energy-flux corollary.* In $d=3$ the energy measured along $y$ after a $T_{ij}$ insertion is, in the polarisation basis
> $((xx-yy)/\sqrt2,(xy+yx)/\sqrt2)$, the matrix $(E/4\pi)\,[[1-t_4/4,\ \alpha_T/4],[\alpha_T/4,\ 1+t_4/4]]$ [CDPrakash, Eq. (4.39)].
> Rotating the detector by $\theta$ rotates its traceless part by $4\theta$. So for any normalised polarisation
> $\langle\mathcal E(\theta)\rangle=(E/2\pi)[1+A_\varepsilon\cos(4\theta-\varphi_\varepsilon)]$, with
> $A_\varepsilon=(\tau/4)\,|\varepsilon_+^2+\varepsilon_\times^2|\le\tau/4$; helicity eigenstates are isotropic. ANEC gives $\tau\le4$.
> The cubic bound $\tau\le2.89\times10^{-8}$ at $\Delta_g=1000$ ($c_1=0$) gives $A\le7.23\times10^{-9}$. With
> $t_4=-4(n_f-n_s)/(n_f+n_s)$, $\alpha_T=8\pi^4p_T/[3(n_f+n_s)]$ and $c_T=n_s+n_f$ [CDPrakash], this reads
> $|n_s-n_f|\le7.23\times10^{-9}c_T$ and $|p_T|\le1.11\times10^{-10}c_T$. The corollary uses only the cubic functional, holds at
> leading order in $1/c_T$ under the same assumptions, and has the same status as the cubic bound. Ref. [Buchel] Eq. (3.8) prints the
> $d=3$ constant as $1/8$; zero net flux requires $1/4$, consistent with $|t_4|\le4$.

**Bib.** Nothing new is needed. CDP is `CDPrakash` in `prl/refs_verified.bib` (JHEP 11 (2017) 171, arXiv:1707.03007) and is
already cited in the Letter.

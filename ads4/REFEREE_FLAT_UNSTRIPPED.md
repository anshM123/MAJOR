# Referee report: `FLAT_UNSTRIPPED.md` and `ads4/code/flatunstripped/`

Referee: independent, hostile; 16 Nov 2026. I did not edit the files under review. My scripts are in `ads4/code/referee_unstripped/`
(`check.py`, ≈ 23 s; `matter.py`; logs in `cache/`, which is git-ignored). They import nothing from `flatunstripped/`. Sources checked:
the arXiv:2201.06602 LaTeX source ((2.1), (2.7a), (2.12), (2.16), (2.23) and the footnote "J₀ < 2", (2.26), App. A), and the CHMRS
LaTeX source (arXiv:2106.10274, `eq:bulkpointbknu`). `bash run_all.sh` gives PASS (19 s).

## Verdicts
* **Correctness: correct.** I found no errors in any formula, identity or number. Everything I recomputed agrees exactly.
* **Soundness:** sound for the flat statements (items 1–4). The AdS statements (item 5) hold only under the [I] dictionary, and two
  of them are worded more strongly than that allows (M2).
* **Labels:** mostly right. Some headings overstate (m1, M2).
* **Applicability:** the flat bound is unchanged. This is a change of basis, so no new constraint arises. The note also contradicts
  `ADS4_FUNCTIONAL.md` without saying so in the place where it matters (M1).
* **Novelty:** low. The identities are elementary algebra. The value of the note is that it corrects the project's own AdS uplift.

## Independent verification
1. **Unstripped amplitude.** On the `s`-cut, `t⁴ = m⁸((1+x)/2)⁴`. So `Im 𝓜 = 16π d₄₄` (full Wigner function) on the `s`-cut and
   `16π P_J` on the `t`-cut.
   * Checked at exact rational `(w, q, J)`: `B′₂` density `= (2w−q)Z²(P+d₄₄)`, `B′₃` density `= Z²(P−d₄₄)`.
   * Regge: the paper bounds the smeared unstripped `M` by `|s|` (footnote: only `J₀ < 2` is used). The kernel is `s^{−k−1}`, so
     `B′_k` converges for `k ≥ 2`.
   * The leading term is `2/m^{2k}` in the measure `dm²/m²`.
2. **Low sides** (my own Laurent-coefficient code; sign convention fixed against the paper's (2.31a)). `−B′₂ … −B′₇` agree term by term
   with the §3 table. The matter low sides of `B′₂` and `B′₃` agree as well.
3. **Identities.** I checked them as integrand identities for a generic non-symmetric test `f`, using reflection antisymmetry
   (`∮g(s) = −∮g(q−s)`). Verified:
   * forward identities, `k = 2, 3`;
   * inverse identities for `B′₂`, `B′₃`;
   * the same identities on the EFT low sides;
   * the heavy-density identity with explicit `P_J` and `₂F₁` at `J = 4, 6, 9`.
   * Algebra: `a⁴ ± b⁴` and the `√Z` series `r⁴ − 2qr⁵ + 2q²r⁶ − …` are correct.
4. **`φ_j` and forward terms.** I rebuilt the `φ_j` from the forward identities and they match §0 item 5. The forward map is
   `−a∂B′₄|₀ + 5aB′₅|₀`, so `25B′₅ → (25/2)Λ¹⁰C⁻_{5,0}` (since `ε₅ = ½`) is right, and `5/2` was wrong. With these, `F_{g3}` in the `B′`
   basis has `g₄, g₅, g₆, g₆′` equal to 0 and `k(0) = 330151/12012`.
5. **`F_U`.**
   * Low-side contact coefficients are exactly `60526/3465`, `−2746189/360360`, `729983/90090`, `−358681/180180`.
   * By mpmath at 30 digits, starting from the definitional densities (with a numerical forward derivative, not the normal form):
     `m⁴F_U(++, J=2, 213/200) = −3.3338483442325` and `m⁴F_U(++, J=2, 6/5) = −1.1277435944154`. Both agree with the document's values.
6. **Item 7.**
   * For even `k`, the `B′_k` density in `dm²/m²` is literally CHMRS `eq:bulkpointbknu`: `(2m²−ν²)P/(m^k(m²−ν²)^{(k+2)/2})`.
   * At leading order it agrees with `C^± = 2P/m^{2k}(|c_same|² ± |c_opp|²)`.
   * It is also the reading favoured by `SKELETON_TEST.md` §8, where spinning moments ≈ the exact scalar skeleton × `D/P` with full `D`.
   * The [I] label is correct.

## Issues (no fatal ones)
**M1 (major, consistency).** Under the `B′` dictionary the `++` density also changes, by `(1−q/w)^{−2}` relative to the stripped one.
That change is what makes `F_U` negative. This contradicts statements in `ADS4_FUNCTIONAL.md`:
* §0 item 3: the bulk-point action is the stripped `F`;
* §0 item 9 and §7.2–7.4: the flat limit is `F′`, with `++` unchanged, and the consequence is "harmless or favourable".

`FLAT_UNSTRIPPED.md` only says that `F′` "is not a statement about any sum rule". It should say explicitly that those items are
superseded, and that on the available evidence (CHMRS, SKELETON §8) `F′` is the less supported reading.

**M2 (major, labels and applicability).** Calling `Ω^U_{g3}` "the correct uplift" overstates the result. Only its flat (bulk-point)
limit is established, and only given the dictionary. Three things are missing:
* its Regge-regime (`β ~ 1`) heavy action, which is needed for positivity in AdS;
* the `C⁺₆`, `C⁻₇`, `C⁻₅` moments, which require odd-`k` `C⁻` moments with two spinning pairs, never computed;
* the AdS light sides of `C₄…C₇`.

Also, "§7.2 … is false" is stated without the condition "given the dictionary".

**Minor.**
* **m1.** The heading of item 5 is labelled "[D + certified]", but the AdS conclusion is [I]. Only the flat facts are [D]/certified.
* **m2.** Item 4, completeness argument. `++` allows only even `J`, and even `P_J` are not complete on `[1−2/m², 1]` when `m < √2`. The
  argument must use `m ≥ √2` or the `+−` channel. Also, §5's "every sum-rule combination" covers only fixed-`u` `f` rules, not `B^{(2)}`
  or the `g`, `h` sum rules.
* **m3.** "The LP optima are equal [D]" holds only for compatible truncations: weight degrees shift by up to +10 in `p`, and `k` by up
  to +5.
* **m4.** "`F_U` bounds nothing" ignores the fact that `g₄`, `g₅`, … are bounded by other sum rules. The point is moot, because `F_U` is
  not positive.
* **m5.** Small corrections to the text:
  * the dropped factor is `16(2J+1)` (with `dm²`), not `16(2J+1)/π`;
  * "`E, O ∝ 1/w⁴`" should read `E ∝ (2w−q)/w⁴`, `O ∝ 1/w⁴`;
  * the header says the run takes ≈ 1 min; it takes ≈ 20 s.
* **m6.** "Exact" bulk-point limit:
  * Even for scalars, CHMRS carries `(1+O(J²/m²))` corrections, i.e. `O(β²)`. These vanish at fixed `b̂` as `Λ → ∞`, so the dictionary
    should be read as a flat-limit statement.
  * The odd-`k` prefactor `Z^{(k+1)/2}` and the even-type prefactor `(2w−q)Z^{(k+2)/2}/2` differ only at `O(q²/w²)`, below the resolution
    of SKELETON §8 (`x ≤ 0.2`). So `ε_odd` and the odd prefactor are untested [I].

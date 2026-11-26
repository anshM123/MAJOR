# Cover letter: draft for the authors to approve

*Draft of 27 Nov 2026 for Letter v17 (positivity of the cubic functional established on every heavy state at Δ_g = 1000). A plain-text copy ready to paste into the APS
submission form is `COVER_LETTER.txt`; the two say the same thing. Fill in the bracketed items before sending.*

---

To the Editors of *Physical Review Letters*

**Submission: "How Close to Einstein Must Holographic Gravity Be?", by Ansh Mishra and Aryan Senthilkumar**

Dear Editors,

We submit the enclosed Letter, with its Supplemental Material, for publication in *Physical Review Letters*.

Holography predicts that a conformal field theory with a large central charge and a large gap in its higher-spin spectrum is
described by Einstein gravity, but so far only as a scaling law: causality forces the non-Einstein graviton couplings to vanish
as powers of the gap, with no explicit constants. Our Letter supplies the numbers for three-dimensional theories. Assuming only a
higher-spin gap Δ_g ≥ 1000 and weak coupling below it, the
stress-tensor three-point function matches Einstein gravity to one part in 10⁸ at leading order in 1/c_T:
√(t₄² + α_T²) ≤ 2.89 × 10⁻⁸, or |ĝ₃|Λ⁴ ≤ 120 in units of the gap, including the parity-odd coupling that exists only in three
dimensions. Most of the factor 10⁸ relative to the conformal-collider bound is the Λ⁻⁴ scaling already implied by causality;
what is new is the coefficient. A direct physical consequence: the energy flux created by the stress tensor is isotropic to one part
in 10⁸ (anisotropy at most 7.2 × 10⁻⁹, where the conformal-collider bound allows 1). We also bound the first stringy (R⁴) correction from above, under a further positivity assumption and two stated conditions; ABJM theory, whose
R⁴ coefficient is known from localisation, sits at 2.6% of that bound, a consistency check.

The method uses causality in its sharpest CFT form: energy detectors (light-transformed stress tensors) commute at spacelike
separation, which gives sum rules on the four-point function ⟨TTTT⟩. The bounds follow once fixed linear functionals are
positive on every heavy state, which has long been the soft spot of bootstrap arguments in anti-de Sitter space. For the cubic
functional we close it. We prove positivity in the flat-space limit; we certify that the leading Regge part keeps one fifth of
the flat-space margin on every heavy state at every gap above 300; and at a gap of 1000 we establish positivity of the full
action on every heavy state, at every spin and mass, by a computer-assisted proof in interval and directed rounding arithmetic.
The proof covers the heavy spectrum with certificates organised by spin, twist and impact parameter; it has to control a
cancellation in which terms of size 10⁶⁸⁵ sum to a number of order one. The certificates for large spin near the twist edge,
where the difficulty lies, were re-run from source on separate hardware with identical output, and the mathematics of each new
lemma was re-derived independently. The cubic bound therefore needs no positivity assumption at a gap of 1000. Two new ingredients make the AdS problem tractable: a dual tensor structure in which curvature causes no helicity
mixing at leading Regge order, and "de-shifted" weights for which every heavy state contributes its flat-space value at that
order. In bulk language the cubic bound reads |ĝ₃|Λ⁴ ≤ 120: defining a scale by |ĝ₃| = M_*⁻⁴, whatever generates the R³
coupling cannot set in below about 0.3 of the gap. The Letter states its assumptions once and does not hide which parts are
proved, which are certified, and which are numerical or conditional (positivity of the quartic functional remains an assumption); the Supplemental Material grades every statement against a
stated legend.

Every number is reproducible. The functional weights are published as exact rationals and tables, separate code reproduces every
constant from them, an independent check rebuilt the key reductions from the written statements alone, and the positivity
certificates are published with their inputs, outputs and replay scripts. Data and code are available at https://doi.org/10.5281/zenodo.24571131.

**Why *Physical Review Letters*.** That a large gap forces Einstein gravity is one of the load-bearing beams of holography, and until now
it has been a statement about exponents. This Letter turns it into a number. In any
three-dimensional holographic CFT with a higher-spin gap of 1000 or more and weak coupling below it, the graviton's cubic couplings deviate from Einstein gravity by at
most 2.9 × 10⁻⁸ at leading order in 1/c_T; equivalently, the energy flux of the stress tensor is isotropic to one part in 10⁸, where causality alone allows
order one. The first stringy correction is capped too, with ABJM theory inside the cap. The positivity that drives such bounds has
long been the soft spot of bootstrap arguments in anti-de Sitter space. Here, for the cubic functional, it is proved on the entire
heavy spectrum at a gap of 1000, by a computer-assisted proof; to our knowledge this is the first such positivity proof at finite
gap in AdS. Separately, without any positivity assumption, the algebra of energy detectors (conservation, conformal symmetry and causality) together
with crossing symmetry bounds the deviation by 0.018 at a gap of 1000, falling as the inverse gap, with an explicit list of inputs. The results answer, with numbers, how much of Einstein gravity
quantum mechanics and causality force on a holographic theory. That question matters across holography and quantum gravity, the
conformal and S-matrix bootstraps, string theory, and the fast-growing study of energy correlators. The main text is written for that
broad readership; the technical detail is in the Supplemental Material.

The manuscript is not published or under consideration for publication elsewhere. It has been posted to the arXiv as
arXiv:[XXXX.XXXXX]. Suggested referees are listed in the
submission form.

Thank you for considering our manuscript.

Sincerely,

Ansh Mishra (independent researcher), ansh.mishra2025@gmail.com
Aryan Senthilkumar (independent researcher), aryan.senthilk@gmail.com

---

## Notes for the authors (not part of the letter)

* **[AUTHORS] arXiv number.** The data URL is filled in (Zenodo DOI). Fill in the arXiv number, or delete the sentence "It has been posted to the arXiv …" if you
  submit to PRL first.
* **AI-use sentence.** Removed from the cover letter (it was optional): APS asks for disclosure in the manuscript, and the
  acknowledgments contain it unchanged. Re-add it here if you prefer.
* **Every number in the letter is a number in the Letter:** 1000, 10⁸, 2.89 × 10⁻⁸, 120, 7.2 × 10⁻⁹, 2.6%. The conditions on the R⁴ bound ("two further stated conditions") are κ > 0 or τ = 0, and no scalar lighter than 0.46Λ
  coupling to TT; the Letter states both.
* **Wording checked against the Letter (v17).** "Prove" for the flat-space limit; "certify" for the leading-Regge floor (T0) and for
  the full action at Δ_g = 1000 (Theorem FP, Supplemental Sec. S11), both computer-assisted. "10⁶⁸⁵" and "re-run ... re-derived"
  are the Letter's wording. Theorem FP is conditional on the same identities the Letter already uses (the angular-route operator
  identities, the certified spectral data, the certified leading-Regge floor); the cover letter does not say "unconditional". The first-of-
  its-kind claim is phrased as the Letter phrases it ("to our knowledge neither has been derived before from ⟨TTTT⟩"); the
  cover letter keeps both qualifiers ("to our knowledge", "derived from ⟨TTTT⟩ in AdS₄/CFT₃"), because CHLPS record flat-space
  and higher-dimensional precedents.

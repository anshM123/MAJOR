# Cover letter: draft for the authors to approve

*Draft of 26 Nov 2026 for Letter v19. A plain-text copy ready to paste into the APS submission form is `COVER_LETTER.txt`;
the two say the same thing. Fill in the bracketed items before sending.*

---

To the Editors of *Physical Review Letters*

**Submission: "How Close to Einstein Must Holographic Gravity Be?", by Ansh Mishra and Aryan Senthilkumar**

Dear Editors,

We submit the enclosed Letter, with its Supplemental Material, for publication in Physical Review Letters.

Holography predicts that a conformal field theory with a large central charge and a large gap in its higher-spin spectrum is described by Einstein gravity. So far this is known only as a scaling law: causality forces the non-Einstein graviton couplings to vanish as powers of the gap, without explicit constants. For three-dimensional theories we compute the constants. Assuming a higher-spin gap Δ_g ≥ 1000, weak coupling below it and no light single trace in T × T, the stress-tensor three-point function agrees with Einstein gravity to one part in 10⁸ at leading order in 1/c_T. In formulas, √(t₄² + α_T²) ≤ 2.89 × 10⁻⁸, or |ĝ₃|Λ⁴ ≤ 120 in units of the gap, and this includes the parity-odd coupling that exists only in three dimensions. Most of the factor 10⁸ relative to the conformal-collider bound is the Λ⁻⁴ scaling already implied by causality; our result is its coefficient. Physically, the energy flux created by the stress tensor is isotropic to one part in 10⁸ (anisotropy at most 7.2 × 10⁻⁹, where the conformal-collider bound allows 1). Under a further positivity assumption and two stated conditions we also bound the first stringy (R⁴) correction from above. The ABJM value, known from localization, is 2.6% of that bound, which is a consistency check.

The method uses the commutativity of energy detectors (light-transformed stress tensors) at spacelike separation, which gives sum rules on the four-point function ⟨TTTT⟩. The bounds follow once fixed linear functionals are positive on every heavy state. Establishing that positivity has been the main open step in bootstrap bounds in anti-de Sitter space, and for the cubic functional we establish it. We certify positivity in the flat-space limit, and certify that in AdS the leading Regge part loses less than 0.1% of the flat-space margin on every heavy state at a gap of 1000. At a gap of 1000 we then establish positivity of the full action on every heavy state, at every spin and mass, by a computer-assisted proof in interval and directed-rounding arithmetic, given the detector identities that also underlie the positivity-free bound. The proof covers the heavy spectrum with certificates organized by spin, twist and impact parameter, and it has to control a cancellation in which terms of size 10⁶⁸⁵ sum to a number of order one. The certificates for large spin near the twist edge were rerun from source on separate hardware with identical output, and the main new lemmas behind them were rederived independently. The cubic bound therefore needs no positivity assumption at a gap of 1000. For the cubic bound, two new ingredients make the AdS problem tractable: a dual tensor structure in which curvature causes no helicity mixing at leading Regge order, and weights chosen so that every heavy state contributes its flat-space value at that order. In bulk terms, if one defines a scale by |ĝ₃| = M_*⁻⁴, whatever generates the R³ coupling cannot set in below 0.30 of the gap if its coefficient is of order one.

The Letter states its assumptions in one place and marks each result as proved, certified, numerical or conditional; positivity of the quartic functional, for instance, remains an assumption. The Supplemental Material grades every statement against a stated legend. The functional weights are published as exact rationals and tables, separate code reproduces every constant from them, key reductions were rederived independently from the written statements alone, and the positivity certificates are published with their inputs, outputs and replay scripts. Data and code are available at https://doi.org/10.5281/zenodo.24571131.

We believe the Letter meets the criteria of Physical Review Letters. That a large gap forces Einstein gravity is a central expectation of holography, but it has so far been established only as a scaling law; this Letter gives the constant. The positivity that drives such bounds is established here, for the cubic functional, on the entire heavy spectrum at a gap of 1000 by a computer-assisted proof. To our knowledge this is the first such positivity proof at finite gap in AdS. Separately, without any positivity assumption, the algebra of energy detectors together with crossing symmetry bounds the deviation by 0.018 at a gap of 1000, falling as the inverse gap, with an explicit list of inputs. Read in reverse, this holds at any gap and turns a measured deviation τ into an upper bound, about 2 + 17.8/τ, on the twist of the lightest higher-spin state, that is, on the string scale of the bulk. The Letter states these results as one Einstein-rigidity theorem valid at every gap, and formulates the positivity proof as a framework that applies at any finite gap: the exact rewriting does not depend on the gap, so positivity at another gap reduces to a weight with the same analytic properties and a new certificate cover. The results quantify how closely unitarity and causality force a holographic theory to follow Einstein gravity. The question is relevant to holography and quantum gravity, the conformal and S-matrix bootstraps, string theory, and the study of energy correlators. We have written the main text for that readership and put the technical detail in the Supplemental Material.

The manuscript is not published or under consideration for publication elsewhere. It has been posted to the arXiv as arXiv:[XXXX.XXXXX]. Suggested referees are listed in the submission form.

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
* **Wording checked against the Letter (v19).** "Prove" for the flat-space limit; "certify" for the leading-Regge floor (T0) and for
  the full action at Δ_g = 1000 (Theorem FP, Supplemental Sec. S11), both computer-assisted. "10⁶⁸⁵" and "rerun ... rederived"
  are the Letter's wording. Theorem FP is conditional on the same identities the Letter already uses (the angular-route operator
  identities, the certified spectral data, the certified leading-Regge floor); the cover letter does not say "unconditional". The first-of-
  its-kind claim is phrased as the Letter phrases it ("to our knowledge neither has been derived before from ⟨TTTT⟩"); the
  cover letter keeps both qualifiers ("to our knowledge", "derived from ⟨TTTT⟩ in AdS₄/CFT₃"), because CHLPS record flat-space
  and higher-dimensional precedents.

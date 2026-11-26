# Executed full sliver coverage

**FULL_SLIVER_PASS_CONDITIONAL.** The exact merged ledger covers every c in [129/500, 9/10] for both b=0 and b=4, at Δ_g=1000, a≥998 and all physical integer J≥1001. This proves the requested 0.3<L<2.3 in both repository conventions, conditional on the retained angular-action identity, authenticated spectral coefficient enclosures, and retained coherent/ordinary-tail bounds. No remaining interval is untested in the continuum certificate.

The common normalized action lower bound is **5.140000E-6**. Channel-specific minima, rounded downward for display, are 5.166280E-6 (b=0) and 5.148037E-6 (b=4). The 11 new cells contain 5,669,470 directed scalar evaluations across both channels, and all 22 new channel-cell inequalities passed the separate analytic/ledger checker. Seven earlier local cells were rechecked and merged.

The last inspected repository ledger (`ads4/code/gap_prove/ledger.json`, branch `claude/eloquent-cori-nftfjy`, blob SHA `a4fa33817e148bb246fa35edf2da6947d10448cd`) records common coverage from c=.28 to c=.56. The present local calculation extends beyond that endpoint; no changes were made to the GitHub repository.

| c interval (exact decimal endpoints) | b=0 lower bound | b=4 lower bound |
|---|---:|---:|
| [0.2580, 0.2680] | 1.380437E-5 | 1.322035E-5 |
| [0.2680, 0.2800] | 1.355361E-5 | 1.299169E-5 |
| [0.2800, 0.2900] | 1.351826E-5 | 1.297521E-5 |
| [0.2900, 0.3000] | 1.339302E-5 | 1.286763E-5 |
| [0.3000, 0.3500] | 1.050459E-5 | 1.026366E-5 |
| [0.3500, 0.4462] | 8.358087E-6 | 8.362481E-6 |
| [0.4462, 0.5305] | 7.981095E-6 | 7.972298E-6 |
| [0.5305, 0.6041] | 1.043479E-5 | 1.026382E-5 |
| [0.6041, 0.6677] | 9.575659E-6 | 9.449694E-6 |
| [0.6677, 0.7187] | 8.819155E-6 | 8.723942E-6 |
| [0.7187, 0.7596] | 8.152688E-6 | 8.079011E-6 |
| [0.7596, 0.7929] | 7.543729E-6 | 7.485170E-6 |
| [0.7929, 0.8204] | 7.010982E-6 | 6.964830E-6 |
| [0.8204, 0.8434] | 6.529089E-6 | 6.491804E-6 |
| [0.8434, 0.8628] | 6.092176E-6 | 6.062101E-6 |
| [0.8628, 0.8794] | 5.690831E-6 | 5.666348E-6 |
| [0.8794, 0.8936] | 5.338072E-6 | 5.318014E-6 |
| [0.8936, 0.9000] | 5.166280E-6 | 5.148037E-6 |

The table uses downward-rounded lower bounds; the full rational numbers are in `sliver_merge/FULL_COVERAGE_LEDGER.json`. The coverage proof uses exact rational endpoints, not the displayed decimals. Since e^(0.3−20/998²)>1/(1−.258) and e^2.3<10 (checked by rational series bounds), the closed c interval covers the target even if L is defined with the physical w correction log(1+qh).

## Reproducibility and limits

Start with `REPRODUCTION.md`. The mathematical norming and directed-arithmetic proof is in `sliver_merge/PROOF.md`. `python check_ledgers.py` recomputes the analytic bounds and exact sign/coverage ledger, and `python replay_scalar.py --jobs 4` independently reruns every scalar shard using portable GMP. The executed scans use corrected exact integer SIMD division; seven representative node subsets have been compared against GMP as exact integers. The entire slower portable replay has not been rerun here. The original 3,000-digit source table was authenticated and all compact enclosures verified; raw moment regeneration and the angular-action identity remain retained dependencies. See `authenticate_original.py` to check provenance against `coeff_long_17200.txt` in the prior full-band reduction archive. Historical subdirectory READMEs describe their own earlier checkpoints, not the status of this assembled certificate.

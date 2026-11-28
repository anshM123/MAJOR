ANGULAR RACAH DUAL: UNIFORM HEAVY-STATE FEASIBILITY
Research proof package — November 27, 2026

MAIN RESULT
For real g>=4, n=floor(g/2), and every real Delta>=J+g, the dual
lambda_k D_k^2 = -W_B (-n)_k (n)_k / (k!)^2
satisfies both heavy-state inequalities implemented in angular_opt:
  helicity 0: J=0,2,4,...;
  helicity 4: every integer J>=4, including odd J.
Here g is the TWIST gap used by the source code.

CONTENTS
  angular_racah_feasibility.pdf   Complete mathematical proof and LP consequence.
  angular_racah_feasibility.tex   Editable LaTeX source.
  verify_angular_racah_feasibility.py
                                Unified exact symbolic verifier.
  coefficient_ledger.json        All 338 positive integer coefficients.
  verification_summary.txt       Verification and source provenance.

REPRODUCTION
  python verify_angular_racah_feasibility.py
Requires Python 3 and SymPy. The verifier recreates coefficient_ledger.json
next to itself. The proof assertions use exact symbolic and integer arithmetic;
there are no numerical scans or floating-point tolerance checks.
The PDF can be rebuilt with pdflatex, run twice for cross-references.

WHAT IS CLOSED
The two-ladder heavy-state feasibility step, uniformly in spin and continuous
real dimension. Its proof uses exact recurrences and decreasing energies, with
one finite polynomial positivity certificate and explicit low-spin certificates.

LP CONSEQUENCE
Under the parent LP's moment-balance and endpoint-budget assumptions,
L=W_B n^2(n^2-1)/400 and tau^2<=16(3+D)/L.
For D=0: tau<=240/[n sqrt(n^2-1) sqrt(g^2+5g-3/2)]
             =960 g^(-3) [1+O(g^(-1))].
This result does not establish an independent viscosity theorem or g^(-4)
scaling. Those are different claims from the feasibility step proved here.

SOURCE CONVENTIONS
Repository: https://github.com/anshM123/MAJOR
Commit: 01b01e6811825e446632d23ee8b33aab15b52b59
Files: ads4/code/angular_opt/lp_dual.py, lp_primal.py, racah_explore.py.
The source repository was read; this package is a standalone addition.

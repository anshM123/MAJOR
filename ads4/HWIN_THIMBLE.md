# H_win by per-wave thimble integration (Result D, Ω^EE fixed: NF02, min-norm gauge)

*20–21 Nov 2026, H_win-thimble agent. Unrefereed, not committed. No existing document or code file edited. New code
`ads4/code/spinfunc/pw_*.py`, caches `spinfunc/cache/pw_*`. ≤ 2 cores (`OMP_NUM_THREADS=1`), < 3 GB per process.*

Grades: [P] proved; [D] derived; [NE] numerical evidence; [M] model; [I] open.

**STATUS: in progress (written early for hand-over; sections are filled as results arrive).**

## Hand-over (kept current)
* Code: `pw_wave.py` (exact waves `W± = ½[P_J ∓ (2i/π)Q_J]`, derivatives; self-test `python3 pw_wave.py`), `pw_core.py`
  (per-wave evaluator on the n3_hybrid contour; wave saddles by homotopy Newton).
* Nothing queued yet.

# Replay of the uniform-gap package (26 Nov 2026)

The seven reproduction commands of RESEARCH_NOTE.md Sec. 6 were re-run unchanged in a separate copy of this directory
(Python 3.13, SymPy, NumPy). Total wall time 57 s.

| step | result |
|---|---|
| `check_uniform.py 32` | 46/46 spin-channel cases PASS |
| `check_physical.py 48 33` | 24/24 PASS |
| `check_physical.py 64 49` | 24/24 PASS; minimum normalised lower bound 2.7185352228409318e-8 (J=63, b=4) |
| `impact_tensor.py --X 400 --K 110 --mq 14 --mv 40 --Jmin 65 --r 1/10` | both channels PASS; H lower bound 0.015195461585320941 (b=4) |
| `check_M1.py` | PASS, compact M1 lower bound 9.083e-5 |
| `check_cutoff.py` | PASS (94 channels > 1e-8; cutoff errors < 2^-500; bhat<=40 transfer) |
| `check_residual.py` | PASS; K(t) exact coefficients as in the note; K at Delta_g=1000 = 0.999967523127696; 16 L_A <= 3022.34 log(2Lambda/3) - 5307.87 |

All JSON outputs are identical to the archived ones apart from timing fields. This is a re-run (reproducibility), not an
independent check; the independent review is in `../review_uniform_gap/`.

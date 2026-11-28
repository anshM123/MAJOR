# Light single-trace repair (collaborator, 2026-11-28)

Replace W_B in the Racah dual by A = 3(sqrt2 - 1) n^2 (n = floor(g/2)).

* Heavy inequalities are preserved: A <= W_B, and both heavy conditions are monotone in the overall scale of
  the dual (W_0 + A R >= W_0 + W_B R when R < 0, and W_4 - A R >= W_4 - W_B R when R > 0; the other signs are trivially positive).
* Light scalars: R_{n,0}(0,Delta) = 1, so W_0 + A > 0.
* Light spin 2, Delta in (3, g+2): W_0 + A R_{n,0}(2,Delta) >= 0.
  Numerical check (coordinator, mpmath, g = 4..10^4, 2e4 points per g): min (W_0 + A R)/W_0 = 1.3e-7 >= 0,
  approached at large g near Delta ~ 0.66 g (tight).
* Cubic coefficient: L = A n^2 (n^2-1)/400, tau <= 80 sqrt(3)/(sqrt(A) n sqrt(n^2-1)) -> 994.415.../g^3.

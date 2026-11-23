"""deshift/c5_binding.py -- DESHIFT_CLOSE.md §1.3: how much can ANY admissible C5 weight move Q's binding per-state values?

At leading Regge every solution c of the C5 de-shift identity has the same image sigma_h P_{-1/2} (Abel injectivity), so only the
beyond-LR part depends on c.  Write the per-harmonic C5 image as P(nu - ij) + Gamma_j(nu)/m^2 (Gamma in c2 units, as Stage 7 §2.1 for C4).
Then delta B5 = (u^3/m^2) Re sum_j int c hhat_j Gamma_j, |delta B5| <= (u^3/m^2) A5 G,  A5 = int |c| sum_j |hhat_j|,  G = sup |Gamma|,
and since FW = -a (dB4 - B5):  |delta(FF/N)| <= a u^3 A5 G/(m^2 N).  At the binding ++ J = 0 gap state u = 1, m^2 = Lam^2.
Inputs: a, N, FF from repro/q_binding.hyb_J0 (read-only import), the quoted FF/N = 0.1427 / 0.1428 (with s_d), A5 and the XP -> new
gauge change dA from cache/ab_check_c5_p015.txt.  Output: cache/c5_binding.txt
"""
import os, re, sys
import kern
sys.path.insert(0, os.path.join(kern.CODE, 'repro')); sys.dont_write_bytecode = True
import q_binding as QB

out = []
txt = open(os.path.join(kern.CACHE, 'ab_check_c5_p015.txt')).read()
A_new = float(re.search(r'A0: new ([0-9.]+)', txt).group(1)); A_xp = float(re.search(r'XP ([0-9.]+);', txt).group(1))
dA = float(re.search(r'dA = int \|w - w_XP\| sum\|R_j\| = ([0-9.e+-]+)', txt).group(1))
for Dg, ffn in ((300, 0.1427), (1000, 0.1428)):
    H = QB.hyb_J0(Dg); Lam2 = H['Lam']**2; N = H['N']; a = QB.a
    s = (f"Dg = {Dg}: a = {a:.5f}, N = {N:.5f}, Lam^2 = {Lam2:.4e}, quoted FF/N = {ffn}.  Tolerance G* = (FF/N) N Lam^2/(a A5): "
         f"A5 = {A_new:.3f} (new) -> {ffn*N*Lam2/(a*A_new):.2e} c2;  A5 = {A_xp:.3f} (XP) -> {ffn*N*Lam2/(a*A_xp):.2e} c2;  "
         f"A5 = 7.66 (worst Stage-6 gauge G3) -> {ffn*N*Lam2/(a*7.66):.2e} c2.  XP -> new gauge change: |d(FF/N)| <= "
         f"{a*dA/(Lam2*N):.2e} G")
    print(s); out.append(s)
open(os.path.join(kern.CACHE, 'c5_binding.txt'), 'w').write("\n".join(out) + "\n")

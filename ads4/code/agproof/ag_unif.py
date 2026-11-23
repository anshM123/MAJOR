"""agproof/ag_unif.py (AG_PROOF.md §1.3): Lemma F's flat waves satisfy the UNIFORM-rate premise (AG_U) on Lemma P4W's whole rectangle.
For theta(nu) = 2 arcsin(nu/m): Im theta(nu_r + i y) >= 2 asinh(y/m) for |nu_r| < m (alpha = (|z+1| + |z-1|)/2 increases with |Re z|), so
the rigorous bound of Lemma F gives |W+(theta(nu_r+iy))| <= (2/pi) C_n e^{-(n+1) 2 asinh(y/m)} (1 - e^{-4 asinh(y/m)})^{-1/2}.  Checked
here: min over y in (0, T], T = 25/(0.97 beta) (nu units), of (n+1) 2 asinh(y/m) / (0.97 beta y), at every sampled state (W47, L').
Also: the claimed monotonicity of Im arcsin(x + i eta) in x on [0, 1) (grid check).  -> cache/unif.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ag_flat as AF
L = AF.L
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))
xs = np.linspace(0, 0.999, 1000)
mono = all(np.all(np.diff(np.arcsin(xs + 1j*e).imag) >= -1e-15) for e in (1e-4, 1e-2, 0.1, 0.5, 1.0, 3.0))
pr(__doc__.split('->')[0])
pr(f'monotonicity of Im arcsin(x + i eta) in x on [0, 0.999] (eta = 1e-4 .. 3): {mono}')
worst = []
for name, sts in (('W47', AF.states_W()), ("L'", AF.states_L())):
    for J, mfl in sts:
        m = mfl*L; beta = 2*np.arcsinh((J + 0.5)/m); T = 25/(0.97*beta)
        y = np.linspace(T*1e-6, T, 4000)
        r = (J + 1)*2*np.arcsinh(y/m)/(0.97*beta*y)
        worst.append((float(r.min()), name, J, mfl))
w = min(worst)
pr(f'min over states and y of (n+1) 2 asinh(y/m) / (0.97 beta y) = {w[0]:.4f} (at {w[1]} J = {w[2]}, m/L = {w[3]:.3f}); '
   f'>= 1 means the flat waves decay at least at the uniform rate 0.97 beta on the whole rectangle')
open(os.path.join(AF.CACHE, 'unif.txt'), 'w').write('\n'.join(OUT) + '\n')

# Abel-summation amplitude bound: |F_n(e^{2 i th})| <= 2/|sin Re th| (Im th >= 0), hence
# |W+(th)| <= (2/pi) C_n e^{-(n+1) Im th} min{2/|sin Re th|, (1 - e^{-2 Im th})^{-1/2}}.
import mpmath as mp
worst = 0.0
for n in (26, 40, 100, 344):
    for tr in (0.2, 0.8, 1.5, 2.3, 2.9):
        for ti in (0.0, 1e-3, 0.05, 0.4):
            th = mp.mpc(tr, ti)
            lhs = abs(AF.Wplus(n, th))
            rhs = 2/mp.pi*AF.Cn(n)*mp.exp(-(n + 1)*ti)*min(2/abs(mp.sin(tr)), (1/mp.sqrt(1 - mp.exp(-2*ti)) if ti > 0 else mp.inf))
            worst = max(worst, float(lhs/rhs))
pr(f'Abel amplitude bound: max |W+| / bound = {worst:.3f} (n = 26..344, Re th in [0.2, 2.9], Im th in [0, 0.4]; must be <= 1)')
open(os.path.join(AF.CACHE, 'unif.txt'), 'w').write('\n'.join(OUT) + '\n')

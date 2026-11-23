"""agproof/ag_pwmeasure.py (AG_PROOF.md §6): how much of the heavy factor's positive plane-wave measure (LARGESPIN_L Corollary B) lies
BEYOND the twist edge nu = Lam, i.e. at |p| < Lam, for the near-gap states (Monte Carlo, exact sampling of the measure):
    p = s n_phi + q,  s ~ Gamma(J+1),  phi uniform,  q ~ (-q^2)^{lambda - 3/2} e^{-q^0} d^3q on V_+  (lambda = (tau - 1)/2),
i.e. q = r(cosh chi, sinh chi cos a, sinh chi sin a) with chi ~ sinh chi / cosh^{2 lambda} chi, r | chi ~ Gamma(2 lambda, cosh chi);
|p|^2 = r^2 + 2 s r (cosh chi - sinh chi cos(a - phi)).   Also with CHMRS's moment weight w = |p|^{-2(k + 2 Delta_phi - 1)} (k = 2,
Delta_phi = 3).  Output: P(|p| < Lam), P(|p| < Lam (1 + 0.02)), the relative width of |p|, against 1 - x1 = 1 - Lam/m.
-> cache/pwmeasure.txt
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(12345)
Dg = 1000; L = np.sqrt((Dg - 2)*(Dg - 1))
OUT = []
pr = lambda s: (print(s, flush=True), OUT.append(s))


def sample(J, tau, N=400000):
    lam = (tau - 1)/2
    # chi: density ~ sinh(chi) cosh(chi)^{-2 lam}; with c = cosh chi: density in c ~ c^{-2 lam} on [1, inf): c = U^{-1/(2 lam - 1)}
    c = rng.random(N)**(-1/(2*lam - 1))
    sh = np.sqrt(c*c - 1)
    r = rng.gamma(2*lam, 1.0, N)/c
    s = rng.gamma(J + 1, 1.0, N)
    psi = rng.random(N)*2*np.pi
    p2 = r*r + 2*s*r*(c - sh*np.cos(psi))
    return np.sqrt(p2)


def main():
    pr(__doc__.split('->')[0])
    pr('  J   tau   m/L    1-x1   | sd|p|/m | P(|p|<Lam)  P(|p|<1.02 Lam) | w-weighted: P(|p|<Lam)  P(|p|<1.02 Lam)')
    for J in (26, 30, 38, 42, 50, 64, 100, 150, 344):
        tau = Dg
        Delta = Dg + J; al = Delta - J - 2.0; ab = al + 2*J + 1; m = np.sqrt(al*ab)
        P = sample(J, tau)
        w = (P/m)**(-2*(2 + 2*3 - 1)); w = w/w.sum()
        f = lambda a: float(np.mean(P < a*L)); fw = lambda a: float(np.sum(w*(P < a*L)))
        pr(f"  {J:3d} {tau:5d} {m/L:.4f} {1 - L/m:.4f} | {np.std(P)/m:.4f}  | {f(1.0):.3f}       {f(1.02):.3f}           | {fw(1.0):.3f}                 {fw(1.02):.3f}")
    open(os.path.join(HERE, 'cache', 'pwmeasure.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()

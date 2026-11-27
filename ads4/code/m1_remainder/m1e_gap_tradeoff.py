r"""
m1_remainder/m1e_gap_tradeoff.py  [estimated below Delta' = 300; computed points at 300 and 1000]

Route 3 cost estimate.  The light action is ~ e^{pi Lam'/2} for a functional of this design built at gap Delta'
(Lam'^2 = (Delta'-1)(Delta'-2)): the e^{pi Lam/2} is the value of the target heavy image continued to the back-to-back point
(imaginary impact parameter i pi/2), so it is set by how far in nu the weight reaches.  Capping the light action at c_T
therefore means using a functional whose reach is Lam' ~ (2/pi) ln(c_T): any theory with gap Delta_g >= Delta' has gap Delta'
(Remark of Sec. allgap.A), states of twist in [Delta', Delta_g) are then heavy and positive.

  log10 c_T*(Lam') = (pi/(2 ln 10)) Lam' + delta(Lam'),   delta from m1c at Delta_g = 300, 1000 (linear in log Lam),
  tau_0(Lam')      = (240/Lam'^4) (16 L_A/K)^{1/2},  16 L_A = A ln(2Lam'/3) + B_inf + b1/Lam' + (b2L ln(2Lam'/3) + b2)/Lam'^2,
                     K = 1 - 68.1/Lam'^2   (Supplement, Sec. on the expansion; table values at 300, 500, 1000 reproduced).
Compared with the angular bound tau_*(Delta_g) = (160/9)/Delta_g, whose light remainder has no exponential constant.
Positivity of the Delta'-functional on heavy states is certified only at 1000 (Theorem FP); at 300 it is conditional.
Usage: python3 m1e_gap_tradeoff.py
"""
import math, os

A, Binf, b1, b2L, b2 = 3022.3377, -5327.2096, 19293.4, -13549.0, -14705.0
HERE = os.path.dirname(os.path.abspath(__file__))


def tau0(L):
    x = math.log(2*L/3)
    l16 = A*x + Binf + b1/L + (b2L*x + b2)/L**2
    K = 1 - 68.1/L**2
    return 240/L**4*math.sqrt(l16/K)


def main():
    out = []
    pr = lambda s: (print(s), out.append(s))
    s = math.pi/(2*math.log(10))
    L1000, L300 = math.sqrt(998*999), math.sqrt(298*299)
    pr(f'check tau_0: 300 -> {tau0(L300):.4e} (table 3.1321e-6), 1000 -> {tau0(L1000):.4e} (table 2.8913e-8)')
    # delta from m1c (signed J <= Jmax sums, n_0 from m1_n0): log10 c_T* - s Lam
    pts = {'signed': {1000: 685.31, 300: 209.08}, 'termwise |.|': {1000: 688.29, 300: 211.33}}
    for lab, d in pts.items():
        d1 = d[1000] - s*L1000; d3 = d[300] - s*L300
        slope = (d1 - d3)/math.log10(L1000/L300)
        pr(f'\n[{lab}] delta(1000) = {d1:.2f}, delta(300) = {d3:.2f}')
        pr('  Delta\'   Lam\'    log10 c_T*   tau_0(Lam\')     tau_*(1000) = 0.0178')
        for Dp in (20, 30, 40, 50, 75, 100, 150, 200, 300, 500, 1000):
            L = math.sqrt((Dp - 1)*(Dp - 2))
            dl = d3 + slope*math.log10(L/L300)
            pr(f'  {Dp:6d}  {L:7.1f}  {s*L + dl:10.1f}   {tau0(L):.3e}   {"better" if tau0(L) < 0.0178 else "worse"} than tau_*(1000)')
    open(os.path.join(HERE, 'out', 'm1e_gap_tradeoff.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()

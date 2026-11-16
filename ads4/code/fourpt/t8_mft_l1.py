"""
Mean-field (unitarity) normalisation of the <JJJJ> heavy formula in the exact helicity basis, and its 1/m structure.
Inputs:
  [Q] CHL (arXiv:2102.08160, eq. 'free data'): MFT data of two conserved currents in their parity-even basis T^{E,opp}, T^{E,same}
      are diagonal:  p[T^{E,h}] = c^E(Delta,J)^s * H_h * C_1^2 (x 2 for n >= 1),  H_opp = (-J)_2/(J+1)_2,  H_same = (Delta-1)_2/(2-Delta)_2,
      c^E^s = scalar MFT data for Delta_1 = Delta_2 = 1 (CHL normalisation).
  [D + exact identification, t1b_helicity_l1.py]  T^{E,h} = Kc * (2^{(J+1)/2} alpha/2) * a_h * hhat_h,  Kc common (frame/convention),
      a_opp = -2/(J(J-1) alpha Y),  a_same = -2/(Delta(Delta-1) alpha Y),  Y = Delta^2 - 3 Delta + J^2 + J.
  [identified + verified, t3/predict4] leading heavy vertex of hhat_h on its own helicity:
      V_opp = A_-^2 m^2 (A_+^2 + A_-^2)/2,   V_same = A_+^2 m^2 (A_+^2 + A_-^2)/2   (plus O(m^-2)-relative other-helicity leakage).
Output: Q_h := p[hhat_h] V_h^2 / [common factor], exact rational; check (i) Q_opp/Q_same -> 1 (equal MFT weights of the two
helicities), (ii) Q_h are Weyl-even (invariant under (A_+, A_-) -> (-A_+, -A_-)), hence expand in 1/m^2 at fixed beta;
(iii) the scalar factor c^E^s/b^(0): for Delta_1 = Delta_2 = 1, P_FK(1,1)/b^(0) = (8/pi^2) A_+ A_- exactly (SPINNING_MFT sect. 4.1 with S = 2),
which is Weyl-even.
"""
import sympy as sp
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

Ap, Am, m, b, eps = sp.symbols('A_p A_m m beta epsilon', positive=True)
D = Ap + sp.Rational(3, 2); Jv = Am - sp.Rational(1, 2)
al = Ap - Am
Y = D**2 - 3*D + Jv**2 + Jv
H = {'opp': (-Jv)*(-Jv + 1)/((Jv + 1)*(Jv + 2)), 'same': (D - 1)*D/((2 - D)*(3 - D))}
a = {'opp': -2/(Jv*(Jv - 1)*al*Y), 'same': -2/(D*(D - 1)*al*Y)}
V = {'opp': Am**2*(Ap**2 - Am**2)*(Ap**2 + Am**2)/2, 'same': Ap**2*(Ap**2 - Am**2)*(Ap**2 + Am**2)/2}
pr("# MFT normalisation of the helicity components of <JJJJ> (fourpt/t8_mft_l1.py)")
Q = {}
for h in ('opp', 'same'):
    # p[hhat] = p[T^E] (Kc 2^{(J+1)/2} alpha/2 a_h)^2 ;  drop the h-independent factors Kc^2 2^{J+1} c^E^s C_1^2
    Q[h] = sp.factor(sp.simplify(H[h]*(al/2*a[h])**2*V[h]**2))
    pr(f"  Q_{h} = p[hhat_{h}] V_{h}^2 / (h-independent factor) = {Q[h]}")
    refl = sp.simplify(Q[h].subs({Ap: -Ap, Am: -Am}, simultaneous=True) - Q[h])
    pr(f"     Weyl reflection (A_+,A_-) -> (-A_+,-A_-): Q - Q_reflected = {refl}")
r = sp.factor(sp.simplify(Q['opp']/Q['same']))
pr(f"  Q_opp/Q_same = {r}")
# expansion at fixed beta: A_+ = m cosh(beta/2), A_- = m sinh(beta/2), m = 1/eps
sub = {Ap: sp.cosh(b/2)/eps, Am: sp.sinh(b/2)/eps}
for h in ('opp', 'same'):
    lead = sp.simplify(sp.limit(Q[h].subs(sub)*eps**4, eps, 0))
    ser = sp.series(sp.simplify(Q[h].subs(sub)*eps**4/lead), eps, 0, 3).removeO()
    pr(f"  Q_{h} = {sp.simplify(lead)} m^4 * (1 + {sp.simplify(ser.coeff(eps, 1))}/m + ({sp.simplify(sp.expand_trig(ser.coeff(eps, 2)))})/m^2 + ...)")
ser = sp.series(sp.simplify(r.subs(sub)), eps, 0, 3).removeO()
pr(f"  Q_opp/Q_same = 1 + {sp.simplify(ser.coeff(eps, 1))}/m + ({sp.simplify(ser.coeff(eps, 2))})/m^2 + ...")
pr("  scalar factor (Delta_1 = Delta_2 = 1):  P_FK(1,1)/b^(0) = (8/pi^2) A_+ A_-  [SPINNING_MFT sect. 4.1, S = 2: G(x) = 2/x]: Weyl-even, no 1/m term")
import os
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t8_mft_l1.txt'), 'w').write("\n".join(out) + "\n")

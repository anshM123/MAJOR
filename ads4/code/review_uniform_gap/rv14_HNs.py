"""rv14_HNs.py -- normalised margins: H_ang/N_s for the NEW family (ideal, uncut; H_ang = w^3 * closed form) at the twist edge,
all physical J<=64, with N_s from the independent xcheck_ladder/xl_norm.hybrid (mpmath, read-only import).
Usage: python3 rv14_HNs.py Dg"""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'xcheck_ladder')))
import mpmath as mp
import xl_norm
from fractions import Fraction as Fr
from rv_closed import closed_form
Dg = int(sys.argv[1]); lam2 = (Dg-1)*(Dg-2); ag = Dg-2
rows = []
for J in range(0, 65):
    for b in (0, 4):
        if (b == 0 and J % 2) or J < b:
            continue
        T = Fr((J+Dg-1)*(J+Dg-2)); m2 = T-J*(J+1)
        S = closed_form(J, b, m2, lam2)
        w3 = (m2/T)**3
        H = w3*S
        mfl = mp.nstr(mp.sqrt(mp.mpf(m2.numerator)/lam2), 30)
        hy = xl_norm.hybrid(Dg, J, mfl, '++' if b == 0 else '+-', npn=200, dps=30)
        Hf = float(H)
        rows.append(dict(J=J, b=b, H_over_Ns=Hf/float(hy['Ns']), Hhyb_over_Ns=float(hy['Hhyb']/hy['Ns']), S=float(S)))
        print(J, b, f"H/Ns = {rows[-1]['H_over_Ns']:.6f}  Hhyb/Ns = {rows[-1]['Hhyb_over_Ns']:.6f}", flush=True)
m = min(rows, key=lambda r: r['H_over_Ns'])
print("MIN H/Ns at twist edge, Dg =", Dg, ":", m)
json.dump(rows, open(f'out/rv14_HNs_{Dg}.json', 'w'), indent=1)

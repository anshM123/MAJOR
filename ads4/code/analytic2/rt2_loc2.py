"""rt2_loc2.py La Lb [J] [Jp]: E_rho by Lemma Loc'' with (i) per-cell f-mass bounds on [898, 1150) and (ii) a different truncation
J' (Jp) for nu >= Lam (where only the residual r = f pairs with K, so the j-truncation index may differ from the contour's J).
  int_cell |f| <= width * sup_cell |psi2(nu/Lam)| + m_r(cell)   (nu < Lam; psi2 = NF02 psi2/k0, exact Taylor shift: rt_cert.psi_abs)
  int_cell |f| <= m(cell)                                          (Lam <= nu < 1059: rt_rmass2 panel masses, f = r there)
  int_T^{T+1} |f| <= tail(T)                                       (1059 <= nu < 1150: CONTOUR)
  nu < 898: sup|f| <= 21.34 times width (as rt_loc3).  nu >= 1150 is not split (rt2_resid far term)."""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common2
import loc4 as L4
import rt_loc3 as R3
import rt_pp as PPm
import rt_cert as RC
from flint import acb
La, Lb = float(sys.argv[1]), float(sys.argv[2]); J = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
Jp = int(sys.argv[4]) if len(sys.argv) > 4 else J
d = math.exp(-La); t0 = time.time()
for JJ in (J, Jp):
    A, B, den, Fi = PPm.rec_step(math.hypot(JJ + 3, 1151) + 2, 9, d)
    print('REC at J =', JJ, ': (A+B)/den =', (A + B)/den, ' F/den =', Fi/den, flush=True)
A, B, den, Fi = PPm.rec_step(math.hypot(J + 3, 1151) + 2, 9, d)
assert (A + B)/den <= 1.99 and Fi/den <= 1.0      # contour J (near-real tails of rt_cert.jsum_near use growth 2 and F/den <= 1)
FSUP = R3.FSUP; LAMf = RC.LAMf
rm = {int(c[0]): c[2] for c in json.load(open(os.path.join(common2.AN, 'cache', 'rt_rmass2.json'))) if c[1] < 1e8}
ls2 = []; worst = []
def add(lm, b, lo):
    ls2.append(lm + b); worst.append((lm + b, lo))
b = L4.cell_bound(3, J, La, 0.0, 0.4, True, Lb); add(math.log(0.4*FSUP), b[0], 0.0)
edges = [0.4, 1.0] + [float(k) for k in range(2, 898)]
for n1, n2 in zip(edges[:-1], edges[1:]):
    b = L4.cell_bound(3, J, La, n1, n2, True, Lb); add(math.log(FSUP*(n2 - n1)), b[0], n1)
for T in range(898, 1150):
    subs = [(float(T), min(float(T + 1), LAMf)), (max(float(T), LAMf), float(T + 1))] if T < LAMf < T + 1 else [(float(T), float(T + 1))]
    for (lo, hi) in subs:
        if hi - lo <= 0: continue
        JJ = J if T < 999 else Jp          # J(nu) = J below 999, J' from 999 on (EPSW panel edges: cell T >= 999 is exactly [T, T+1])
        b = L4.cell_bound(3, JJ, La, lo, hi, True, Lb)
        if lo < LAMf:
            mid = (lo + hi)/2/LAMf; rad = (hi - lo)/2/LAMf + 1e-12
            ps = RC.psi_abs(RC.Q2, acb(mid), rad)
            mass = (hi - lo)*ps + rm[T]
        elif T < 1059:
            mass = rm[T]
        else:
            mass = math.exp(R3.tail_mass(T))
        lm = math.log(mass) if mass > 0 else R3.tail_mass(T)
        add(lm, b[0], lo)
lI2 = R3.logsum(ls2)
ls3 = []
b = L4.cell_bound(4, J, La, 0.0, 0.4, False, Lb); ls3.append(math.log(0.4/LAMf) + b[0])
lo = 0.4
while lo < LAMf:
    hi = min(LAMf, math.floor(lo) + 1.0)
    b = L4.cell_bound(4, J, La, lo, hi, False, Lb)
    ls3.append(math.log((hi - lo)/LAMf) + b[0]); lo = hi
lI3 = R3.logsum(ls3)
s3 = sum(abs(float(c)) for c in RC.P3)
lE2 = 3*math.log(d) + math.log(2/RC.LAMf) + math.log(2) + lI2 + math.log(2)
lE3 = 3*math.log(d) + math.log(d*d/RC.TAU_COEF) + 2*math.log(RC.LAMup) + math.log(s3) + math.log(2) + lI3 + math.log(2)
worst.sort(reverse=True)
out = dict(La=La, Lb=Lb, J=J, Jp=Jp, lI2=lI2, lI3=lI3, lE2=lE2, lE3=lE3, worst2=worst[:6])
json.dump(out, open(os.path.join(common2.HERE, 'cache', f'loc_{La}_{Lb}_{J}.json'), 'w'))
print(out, f'[{time.time()-t0:.0f}s]')

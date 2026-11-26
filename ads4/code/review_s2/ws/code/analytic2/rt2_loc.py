"""rt2_loc.py La Lb [J]: E_rho by Lemma Loc'' (loc4.py) for nu < 1150 (k = 2, p = 3, dressed) and the k = 3 piece (p = 4, nu <= Lam);
nu >= 1150 is NOT split (the residual uses the full kernel there, see rt2_cert.far_term).  Outputs log integrals and log E_rho2/3."""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common2
import loc4 as L4
import rt_loc3 as R3
import rt_pp as PPm
import rt_cert as RC
La, Lb = float(sys.argv[1]), float(sys.argv[2]); J = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
d = math.exp(-La); t0 = time.time()
# j-sum growth hypotheses used by the geometric tails (REC, M'_j x2 per step) up to |beta| = hypot(J+3, 1151) + 2
A, B, den, Fi = PPm.rec_step(math.hypot(J + 3, 1151) + 2, 9, d)
assert (A + B)/den <= 1.99 and Fi/den <= 1.0, (A, B, den, Fi)
FSUP = R3.FSUP; LAM = math.sqrt(998*999)
ls2 = []; worst = (-1e9, None)
b = L4.cell_bound(3, J, La, 0.0, 0.4, True, Lb); ls2.append(math.log(0.4*FSUP) + b[0])
edges = [0.4, 1.0] + [float(k) for k in range(2, 1031)]
for n1, n2 in zip(edges[:-1], edges[1:]):
    b = L4.cell_bound(3, J, La, n1, n2, True, Lb)
    ls2.append(math.log(FSUP*(n2 - n1)) + b[0])
    if b[0] > worst[0]: worst = (b[0], n1)
for T in range(1030, 1150):
    b = L4.cell_bound(3, J, La, float(T), float(T + 1), True, Lb)
    ls2.append(R3.tail_mass(T) + b[0])
    if b[0] + R3.tail_mass(T) > -50: print('  warn cell', T, b[0], R3.tail_mass(T))
lI2 = R3.logsum(ls2)
ls3 = []
b = L4.cell_bound(4, J, La, 0.0, 0.4, False, Lb); ls3.append(math.log(0.4/LAM) + b[0])
lo = 0.4
while lo < LAM:
    hi = min(LAM, math.floor(lo) + 1.0)
    b = L4.cell_bound(4, J, La, lo, hi, False, Lb)
    ls3.append(math.log((hi - lo)/LAM) + b[0]); lo = hi
lI3 = R3.logsum(ls3)
s3 = sum(abs(float(c)) for c in RC.P3)
lE2 = 3*math.log(d) + math.log(2/RC.LAMf) + math.log(2) + lI2 + math.log(2)
lE3 = 3*math.log(d) + math.log(d*d/RC.TAU_COEF) + 2*math.log(RC.LAMup) + math.log(s3) + math.log(2) + lI3 + math.log(2)
out = dict(La=La, Lb=Lb, J=J, lI2=lI2, lI3=lI3, lE2=lE2, lE3=lE3, worst2=worst)
json.dump(out, open(os.path.join(common2.HERE, 'cache', f'loc_{La}_{Lb}_{J}.json'), 'w'))
print(out, f'[{time.time()-t0:.0f}s]')

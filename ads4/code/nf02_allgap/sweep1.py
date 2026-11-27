"""sweep1.py -- numerical sweep (NOT a certificate) of S=H/w^3 and H/H_LR in the open region J>=65, bhat>40,
at several finite gaps, including large impact parameter / near-gap twist at large J (never sampled before
at gaps other than 1000).  Usage: python3 sweep1.py DELTA_G K  -> out/sweep1_<DG>.json"""
import sys, json, math, time, os
from fractions import Fraction as Fr
from hev import Coeffs, S_state, H_LR, state_params
DG = Fr(sys.argv[1]); K = int(sys.argv[2])
lam2 = (DG-2)*(DG-1); ag = DG-2
t0 = time.time()
cf = Coeffs(lam2, K)
print('coeffs', time.time()-t0, flush=True)
Js = sorted(set([65, 66, 69, 70, 80, 100, 130, 131, 160, 200, 260, 330, 420, 530, 670, 850, 1070, 1350, 1700, 2150, 2700,
           3400, 4300, 5400, 6800, 8600, 10800, 13600, 17200, 21700, 27300, 34400, 43400, 54700, 69000, 87000, 110000, 140000, 200000]))
xs = [0, Fr(1, 1000), Fr(1, 100), Fr(3, 100), Fr(1, 10), Fr(3, 10), 1, 2, 4, 10, 30]
res = []
for J in Js:
    for b in (0, 4):
        if b == 0 and J % 2: pass  # both parities kept: H_ang,0 is evaluated for every J (stronger than needed)
        for x in xs:
            a = ag*(1+x)
            sp = state_params(lam2, J, a)
            if sp['bhat'] <= 40: continue
            S, w, info = S_state(cf, J, b, a)
            lr = H_LR(cf, J, b, a)
            H = (w**3*S)
            rec = dict(J=J, b=b, a=float(a), bhat=sp['bhat'], ell=sp['ell'], r=sp['r'], S=float(S.mid()), Srad=float(S.rad()),
                       H=float(H.mid()), HLR=float(lr.mid()), ratio=float((H/lr).mid()) if lr.mid() != 0 else None,
                       trunc=info['truncated'], last=max(abs(t) for t in info['last_terms']) if info['last_terms'] else 0)
            res.append(rec)
            flag = '' if (S > 0) else ' <<<< NONPOS'
            print(J, b, float(a), '%.3g %.4g S=%.6g H/HLR=%s last=%.2g%s' % (sp['bhat'], sp['ell'], rec['S'], rec['ratio'], rec['last'], flag), flush=True)
os.makedirs('out', exist_ok=True)
json.dump(res, open('out/sweep1_%s.json' % str(float(DG)), 'w'), indent=0)
print('done', time.time()-t0)

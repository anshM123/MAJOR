"""sweep2.py -- larger-K numerical sweep (NOT a certificate) of the moving near-gap region at large spin.
Usage: python3 sweep2.py DELTA_G K -> out/sweep2_<DG>.json"""
import sys, json, math, time, os
from fractions import Fraction as Fr
from hev import Coeffs, S_state, H_LR, state_params
DG = Fr(sys.argv[1]); K = int(sys.argv[2])
lam2 = (DG-2)*(DG-1); ag = DG-2
t0 = time.time(); cf = Coeffs(lam2, K); print('coeffs', time.time()-t0, flush=True)
Js = [int(round(1000*1.4**i)) for i in range(0, 25)]
xs = [0, Fr(1, 200), Fr(1, 20), Fr(1, 5), Fr(1, 2), 1, 3]
res = []
for J in Js:
    for b in (0, 4):
        for x in xs:
            a = ag*(1+x); sp = state_params(lam2, J, a)
            if sp['bhat'] <= 40: continue
            S, w, info = S_state(cf, J, b, a)
            valid = (not info['truncated']) or max(abs(t) for t in info['last_terms']) < 1e-25*abs(float(S.mid()) or 1e-300)
            rec = dict(J=J, b=b, a=float(a), bhat=sp['bhat'], ell=sp['ell'], S=float(S.mid()), Srad=float(S.rad()), valid=bool(valid and S.rad() < 1e-10*abs(S.mid())),
                       last=max(abs(t) for t in info['last_terms']), trunc=info['truncated'])
            res.append(rec)
            if not rec['valid']: continue
            print(J, b, float(a), 'bhat=%.4g ell=%.4g S=%.6g w^3S=%.4g' % (sp['bhat'], sp['ell'], rec['S'], float((w**3*S).mid())), '' if S > 0 else '<<<< NONPOS', flush=True)
json.dump(res, open('out/sweep2_%s.json' % str(float(DG)), 'w'), indent=0)
print('done', time.time()-t0, 'valid', sum(r['valid'] for r in res), 'of', len(res))

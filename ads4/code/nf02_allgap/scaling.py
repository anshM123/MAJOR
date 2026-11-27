"""scaling.py -- numerical test of the scaling structure (J,a,Lam) -> (lam J, lam a, lam Lam):
the leading-Regge part is invariant up to O(1/a); does H/H_LR - 1 (finite-mass correction) shrink like a power of 1/Lam?
Numerical evidence only.  Usage: python3 scaling.py DG K"""
import sys, json
from fractions import Fraction as Fr
from hev import Coeffs, S_state, H_LR, state_params
DG = Fr(sys.argv[1]); K = int(sys.argv[2]); lam = DG/1000
lam2 = (DG-2)*(DG-1)
cf = Coeffs(lam2, K)
base = [(65, 0, 998), (66, 4, 998), (100, 0, 998), (101, 4, 998), (200, 0, 998), (201, 4, 998), (400, 0, 998), (401, 4, 998),
        (100, 0, 1996), (101, 4, 1996), (300, 0, 2994), (301, 4, 2994), (700, 0, 1497), (701, 4, 1497), (1000, 0, 4990), (1001, 4, 4990)]
out = []
for (J0, b, a0) in base:
    J = int(round(J0*lam)); J += (J - b) % 2 if b == 0 and J % 2 else 0
    a = Fr(a0)*lam
    sp = state_params(lam2, J, a)
    if sp['r'] > 0.4: continue
    S, w, info = S_state(cf, J, b, a); lr = H_LR(cf, J, b, a); H = w**3*S
    rec = dict(J0=J0, b=b, a0=a0, J=J, a=float(a), bhat_over_lam=sp['bhat']/float(lam), ell=sp['ell'], H=float(H.mid()), HLR=float(lr.mid()),
               dev=float((H/lr-1).mid()), last=max(abs(t) for t in info['last_terms']), Srad=float(S.rad()))
    out.append(rec); print(rec, flush=True)
json.dump(out, open('out/scaling_%s.json' % sys.argv[1].replace('/', '_'), 'w'), indent=0)

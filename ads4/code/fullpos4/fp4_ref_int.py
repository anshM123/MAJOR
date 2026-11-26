"""
fullpos4/fp4_ref_int.py (FULLPOS4.md §5.1) [DV]: integrated §78 reference vs the direct ladder at a physical state, weight psi2:
  H_ref = 2 int_0^1 psi2(p) w^3 K_ref(Lam p) dp,  H_dir = 2 int psi2 w^3 sum_{k<=n} U_k g_k(Lam p) dp   (same Gauss-Legendre p-nodes),
  and the leading-Regge value w^3 G(c).   Midpoint values (DV).   Usage: fp4_ref_int.py J a b [npts] [M] [prec]
"""
import sys, os, json, time, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'largespin'))
import numpy as np
from flint import arb, ctx
import fp4_ref as R, fp4_int as I
J = int(sys.argv[1]); a = float(sys.argv[2]); b = int(sys.argv[3])
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 1600; M = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
ctx.prec = int(sys.argv[6]) if len(sys.argv) > 6 else 4000
t0 = time.time()
w2 = I.Weight(2)
x, wx = np.polynomial.legendre.leggauss(npts); ps = (x + 1)/2; wp = wx/2
P = R.params(J, a, b); W, big = R.Wcoef(P['c'], P['h'], P['d'], P['s'], M)
Hr = 0.0; Hd = 0.0; Ha = 0.0; mx = 0.0
for p_, w_ in zip(ps, wp):
    nu = I.LAM*p_; psi = float(w2.psi(p_))
    r = float(R.kref(3, nu, P, W).mid()); e = float(R.kexact(3, nu, J, a, b).mid())
    Hr += 2*psi*w_*r; Hd += 2*psi*w_*e; Ha += 2*abs(psi)*w_*abs(e); mx = max(mx, abs(r - e))
import ls_core as LSC
c = float(P['c'].mid()); w3 = float((P['w']**3).mid())
Gc = (1 - c)**(-3)*2*LSC.OmLR(2, [math.acosh((1 + c)/(1 - c))])[0]
out = dict(J=J, a=a, b=b, M=M, prec=ctx.prec, npts=npts, H_ref=Hr, H_dir=Hd, rel=abs(Hr - Hd)/abs(Hd), max_nu_diff=mx, abs_int=Ha,
           w3Gc=w3*Gc, H_dir_over_w3Gc=Hd/(w3*Gc), log10maxR=big, secs=time.time() - t0)
print(json.dumps(out), flush=True)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'fp4_ref_int.jsonl'), 'a').write(json.dumps(out) + '\n')

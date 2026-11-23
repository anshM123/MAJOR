"""t_deshift.py: the de-shift identity with the PUBLISHED grid weight X (d_abel_X_1000.txt) and kernels R_-+ (kernels.json):
   Re sum dnu X A[P_beta](nu) / Lam  ==  sum dnu psi2 P_beta(nu) / Lam      (A = local de-shift operator of the chains)
 and the 'raw-read-as-stripped' mismatch Re sum dnu (X - psi2) P_beta / Lam.  -> cache/t_deshift.txt"""
import os, sys
sys.dont_write_bytecode = True
import numpy as np
import n2_common as N
Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
D = N.xdata(Dg); L = N.lam(Dg); nu, X, ps, dn = D['nu'], D['X'], D['psi2'], D['dnu']
gp, gm = N.gpm(nu)
out = []
for beta in (0.0, 0.02, 0.05, 0.1, 0.1836, 0.3, 0.5, 1.0):
    P0 = N.Pnu(beta, nu).real
    Pp = N.Pnu(beta, nu + 1j); Pm = N.Pnu(beta, nu - 1j)
    lhs = (X*(gp*Pp + gm*Pm)).real.sum()*dn/L
    rhs = (ps*P0).sum()*dn/L
    raw = ((X - ps)*P0).real.sum()*dn/L
    w25 = ((X - ps)*P0)[nu < 25].real.sum()*dn/L
    out.append('beta=%.4f  Re sum X A[P] = %.10e   sum psi2 P = %.10e   rel.diff %.1e   | Re sum (X-psi2) P = %+.3e (nu<25: %+.3e)  ratio to psi2-action %+.2f'
               % (beta, lhs, rhs, abs(lhs/rhs - 1), raw, w25, raw/rhs))
txt = '\n'.join(out); print(txt)
open(os.path.join(N.CACHE, 't_deshift_%d.txt' % Dg), 'w').write(txt + '\n')

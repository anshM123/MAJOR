"""
pwimg/pw5_theoremL.py (PW_IMAGES.md §6): Theorem L (J >= J*(m), b_hat >= 1.0421 Lam) with premise P_L replaced by the exact
de-shift identity (pw1) and the exact cost formula (pw3).  For gap states (and 1.3 m_gap) J = 857 ... 10^6:
  * F_LR(beta, u, sigma) (lsl_core.LR_action, read-only), min over sigma = +-1, u = Lam^2/m^2;
  * monopole tolerance  Q_tol := min_h F_LR / (4 |D|)   (the monopole uses half of Theorem L's per-entry budget 1/2 min F_LR; the
    k = 2 entry carries ring weight 2, so the monopole cost of an entry is 2 |D| |Q|);
  * the monopole of a scalar-sized O(1/m^2) correction: |Q| <= c* Lam_amp (|kappa_A| + |kappa_B|/4) eta/m^2, kappa^sc = (6.46, 11/6),
    c* = max_s |c_s(i/2)| over the O(1/m^2) radial shifts (pw2: 0.50, 0.914; 1.25 incl. (3,3));
  * exact per-unit costs (c2 units) of the ring classes and index-shifted harmonics in EE de-shift form R1 and in the c-ratio form R1'
    (pw4 definition), relative to 1/2 min F_LR: K_tol.  P_{-1/2+i nu}(cosh beta) by Mehler with nphi = 600 + 600 beta; the R1 cost of
    eta P is cross-checked against lsl_core.Omega (closed-form Phi).
Output: cache/pw5_theoremL.txt, cache/pw5_theoremL.json
"""
import os, sys, io, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'lsl')); sys.path.insert(0, os.path.join(HERE, '..', 'n2sub'))
import numpy as np
import lsl_core as C           # read-only
import n2_common as N          # read-only
from pw_common import cfac
import pw4_theoremW as P4      # CS (exact c_s), SC

OUT = io.StringIO()
pr = lambda s: (print(s, flush=True), OUT.write(s + '\n'))
D = json.load(open(os.path.join(HERE, 'cache', 'pw3_cost.json')))['1000']['D']
XD = N.xdata(1000); L = N.lam(1000)
sel = XD['nu'] < L
NU, X, PS, DN = XD['nu'][sel], XD['X'][sel], XD['psi2'][sel], XD['dnu']
KA, KB = 86.75/6 - 8.0, 11/6


def Pd(beta, s=0.0):
    nph = int(600 + 600*beta)
    return N.Pnu(beta, NU - 1j*s, nphi=nph)


def main():
    pr(__doc__.split('Output')[0])
    pr(f'D = {D:+.6e}')
    cst = {s: P4.CS[s](0.5j).real for s in P4.CS}
    c2 = max(abs(cst[(1, 1)]), abs(cst[(0, 2)])); c4 = max(abs(v) for v in cst.values())
    rows = []
    for J in (857, 1000, 2000, 3000, 10000, 100000, 1000000):
        for fm in (1.0, 1.3):
            m = C.mgap(1000, J)*fm
            u = L**2/m**2 if m > L else 1.0
            beta = C.beta_of(J, m); eta = np.cosh(beta)
            if 2*L*np.sinh(beta/2) < 1.0421*L:
                continue
            FLR = min(C.LR_action(beta, u, s) for s in (1.0, -1.0))
            Qtol = FLR/(4*abs(D))
            Qest2 = c2*(abs(KA) + abs(KB)/4)*eta/m**2      # Lam_amp = 1, O(1/m^2) shifts
            P0 = Pd(beta).real
            dP = (Pd(beta*(1 + 1e-4)).real - Pd(beta*(1 - 1e-4)).real)/(2e-4*beta)
            sh1 = Pd(beta, 1.0).real; sh2 = Pd(beta, 2.0).real
            S = {'A': eta*P0, 'B': eta*NU**2*P0, 'C1': dP/np.sinh(beta), 'C2': np.sinh(beta)*NU**2*dP, 'C3': np.sinh(beta)*dP,
                 'sh1': sh1, 'sh2': sh2}
            Sh = {'A': eta, 'B': -eta/4, 'C1': 0, 'C2': 0, 'C3': 0, 'sh1': (1 + eta)/2, 'sh2': (eta + (3*eta**2 - 1)/2)/2}
            om_chk = (PS*P0).sum()*DN/L, C.Omega(2, beta) if beta < 2 else C.Omega_form(2, beta)
            K = {}
            for k in S:
                r1 = abs(2*(PS*S[k]).sum()*DN/L)/m**2
                r1p = max(abs(2*((PS*(P4.CS[s](NU)/cfac(NU))*S[k]).sum()*DN/L + D*P4.CS[s](0.5j).real*Sh[k])) for s in ((1, 1), (0, 2)))/m**2
                K[k] = (0.5*FLR/r1 if r1 > 0 else np.inf, 0.5*FLR/r1p)
            row = dict(J=J, fm=fm, m=m, beta=beta, eta_m2=eta/m**2, FLR=FLR, Qtol=Qtol, Qest=Qest2, Lam_tol=Qtol/Qest2, K=K,
                       omchk=om_chk)
            rows.append(row)
            pr(f'J = {J:7d} m = {fm}*m_gap (m/Lam = {m/L:.3f}) beta = {beta:.3f}: min F_LR = {FLR:.3e}; Q_tol = {Qtol:.2e}; '
               f'eta/m^2 = {eta/m**2:.2e} (x(tau-2)^2 = {eta/m**2*998**2:.3f}); scalar-sized |Q| <= {Qest2:.2e} (Lam_amp = 1) '
               f'=> tolerated Lam_amp = {Qtol/Qest2:.0f}   [Omega2 check grid/lsl: {om_chk[0]:.4e}/{om_chk[1]:.4e}]')
            pr('      K_tol (c2 units, half budget) R1 / R1\': ' + '  '.join(f'{k}: {v[0]:.3g}/{v[1]:.3g}' for k, v in K.items()))
            sc = {'A': abs(KA), 'B': abs(KB), 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}
            vt = [1/sum(sc[k]/K[k][i] for k in sc) for i in (0, 1)]
            row['vec_tol'] = vt
            shl = K['sh2'][1]*m**2/eta
            row['sh2_lc'] = shl
            pr(f'      scalar-shaped O(1/m^2) vector (kappa^sc of CHMRS (3.20) + residual): tolerated Lam_amp R1 {vt[0]:.3g} / R1\' {vt[1]:.3g};'
               f'  sh2 with lightcone-consistent coefficient kappa eta/m^4: K_tol(R1\') = {shl:.3g}')
    json.dump(dict(D=D, cs={str(k): v for k, v in cst.items()}, rows=rows), open(os.path.join(HERE, 'cache', 'pw5_theoremL.json'), 'w'), default=float)
    open(os.path.join(HERE, 'cache', 'pw5_theoremL.txt'), 'w').write(OUT.getvalue())


if __name__ == '__main__':
    main()

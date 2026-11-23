"""
lrtight/lrt5_tolerances.py (LRINF_TIGHT.md §4-§5): Theorem L (PW_IMAGES Theorem L') premise tolerances at the new J*, the
certified floor against the numerical one, and the hybrid-vs-F_LR relation in the newly covered band.  [DV] (floating point;
the floor used for the tolerances is either the numerical min_sigma F_LR, as in pw5_theoremL.py, or the CERTIFIED lower bound
F_cert of lrt3 (valid for every u <= 1, sigma = +-1 on the beta-cell containing the state).

Read-only imports: pwimg/pw5_theoremL.py (and through it lsl_core, n2_common, pw4_theoremW, pw_common); largespin/ls_core.py
(hyb_fast).  The per-state formulas are those of pw5_theoremL.main(), re-run at new (J, m); nothing of theirs is written.
  * Q_tol = F/(4|D|), D = -9.983e-4 (monopole cost per unit F(i/2));
  * scalar-sized monopole |Q| <= c* (|kA| + |kB|/4) eta/m^2 (Lam_amp = 1) => tolerated Lam_amp (monopole only) = Q_tol/|Q|;
  * full scalar-shaped O(1/m^2) vector in the exact c-ratio form R1' (worst radial shift (1,1)/(0,2)): tolerated Lam_amp;
  * K_tol for index shifts sh+-1 (R1') and the ring class C2.
Also: (i) sanity [DV]: F_cert <= numerical min_{u<=1,sigma} F_LR and |E_k| (two-term form vs direct Mehler) <= certified bound;
(ii) hybrid (ls_core.hyb_fast) vs F_LR at J = 139 ... 857, gap and 1.3 m_gap.
Output: cache/lrt5_tolerances.txt, cache/lrt5_tolerances.json
"""
import os, sys, json, io
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PW = os.path.join(HERE, '..', 'pwimg')
sys.path.insert(0, os.path.join(HERE, '..', 'largespin')); sys.path.insert(0, PW)
from fractions import Fraction as Fr
import numpy as np
import pw5_theoremL as P5          # read-only
import lsl_core as C               # read-only (via pw5's path)
import ls_core as LS               # read-only
from lrt_common import CACHE

OUT = io.StringIO()
def pr(s=''):
    print(s, flush=True); OUT.write(s + '\n')

L = C.lam(1000)
R3 = json.load(open(os.path.join(CACHE, 'lrt3_rcert.json')))
ROWS = [(Fr(r['ba']), Fr(r['bb']), r) for r in R3['rows']]


def cell_of(beta):
    for a, b, r in ROWS:
        if float(a) <= beta <= float(b):
            return r
    return None


def Fcert(beta):
    r = cell_of(beta)
    return None if r is None else r['Fcert']


def tolerances(J, m, F):
    """pw5_theoremL.main() per-state block, with the floor F supplied."""
    P4 = P5.P4
    NU, PS, DN = P5.NU, P5.PS, P5.DN
    u = L**2/m**2 if m > L else 1.0
    beta = C.beta_of(J, m); eta = np.cosh(beta)
    cst = {s: P4.CS[s](0.5j).real for s in P4.CS}
    c2 = max(abs(cst[(1, 1)]), abs(cst[(0, 2)]))
    Qtol = F/(4*abs(P5.D))
    Qest = c2*(abs(P5.KA) + abs(P5.KB)/4)*eta/m**2
    P0 = P5.Pd(beta).real
    dP = (P5.Pd(beta*(1 + 1e-4)).real - P5.Pd(beta*(1 - 1e-4)).real)/(2e-4*beta)
    sh1 = P5.Pd(beta, 1.0).real
    S = {'A': eta*P0, 'B': eta*NU**2*P0, 'C1': dP/np.sinh(beta), 'C2': np.sinh(beta)*NU**2*dP, 'C3': np.sinh(beta)*dP, 'sh1': sh1}
    Sh = {'A': eta, 'B': -eta/4, 'C1': 0, 'C2': 0, 'C3': 0, 'sh1': (1 + eta)/2}
    K = {}
    for k in S:
        r1 = abs(2*(PS*S[k]).sum()*DN/L)/m**2
        r1p = max(abs(2*((PS*(P4.CS[s](NU)/P5.cfac(NU))*S[k]).sum()*DN/L + P5.D*P4.CS[s](0.5j).real*Sh[k]))
                  for s in ((1, 1), (0, 2)))/m**2
        K[k] = (0.5*F/r1 if r1 > 0 else np.inf, 0.5*F/r1p)
    sc = {'A': abs(P5.KA), 'B': abs(P5.KB), 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}
    vt = [1/sum(sc[k]/K[k][i] for k in sc) for i in (0, 1)]
    return dict(beta=beta, u=u, Qtol=Qtol, eta_m2=eta/m**2, Lam_mono=Qtol/Qest, Lam_vec_R1p=vt[1], Lam_vec_R1=vt[0],
                Ksh1=K['sh1'][1], KC2=K['C2'][1])


def main():
    pr(__doc__.split('Output')[0])
    J4 = json.load(open(os.path.join(CACHE, 'lrt4_jstar.json')))
    b0r = Fr(R3['beta0_r_1.0']); b0t = Fr(R3['beta0_rt_1.0'])
    pr(f'certified thresholds: beta0(r<1) = {float(b0r)}, beta0(rtilde<1) = {float(b0t)}')

    pr('\n(i) sanity [DV]: certified floor vs numerical min_{u<=1, sigma} F_LR = 2 Omega_2 - |Omega_3| (lsl_core, direct Mehler)')
    for beta in (0.25, 0.3, 0.37, 0.45, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0):
        o2 = C.Omega(2, beta, nphi=12000); o3 = C.Omega(3, beta, nphi=12000)
        Fnum = 2*o2 - abs(o3)
        r = cell_of(beta)
        f2 = C.Omega_form(2, beta); f3 = C.Omega_form(3, beta)
        pr(f'   beta = {beta:4.2f}: F_num = {Fnum:.5e}  F_cert = {r["Fcert"]:.5e}  (ratio {r["Fcert"]/Fnum:.3f});  '
           f'|E2|_num = {abs(o2 - f2):.2e} <= {r["E2"]:.2e};  |E3|_num = {abs(o3 - f3):.2e} <= {r["E3"]:.2e}:  '
           f'{abs(o2 - f2) <= r["E2"] and abs(o3 - f3) <= r["E3"] and r["Fcert"] <= Fnum}')

    pr('\n(ii) Theorem L\' tolerances (pw5 definitions) at new J*; floor = numerical min_sigma F_LR at the state\'s u / certified F_cert')
    pr('| state | J | m/Lam | beta | min_s F_LR (num) | F_cert | Q_tol num / cert | Lam_amp tol (monopole) num / cert | '
       'Lam_amp tol (scalar vector, R1\') num / cert | K_tol sh+-1 (R1\') cert | K_tol C2 cert |')
    res = []
    states = []
    for J in (139, 175, 206, 223, 250, 276, 300, 324, 344, 400, 500, 600, 700, 857, 3000, 10000, 100000, 1000000):
        states.append(('gap', J, C.mgap(1000, J)))
    for J in (250, 307, 400, 600):
        states.append(('1.1 m_gap' if J == 250 else '1.3 m_gap', J, C.mgap(1000, J)*(1.1 if J == 250 else 1.3)))
    states += [('1.5 Lam', 279, 1.5*L), ('2 Lam', 372, 2*L), ('5 Lam', 929, 5*L)]
    for lab, J, m in states:
        beta = C.beta_of(J, m); u = L**2/m**2 if m > L else 1.0
        Fn = min(C.LR_action(beta, u, s) for s in (1.0, -1.0))
        Fc = Fcert(beta)
        if Fc is None or Fc <= 0:
            pr(f'| {lab} | {J} | {m/L:.3f} | {beta:.3f} | {Fn:.3e} | (beta below certified range) |'); continue
        tn = tolerances(J, m, Fn); tc = tolerances(J, m, Fc)
        res.append(dict(state=lab, J=J, m_over_L=m/L, beta=beta, Fnum=Fn, Fcert=Fc, num=tn, cert=tc))
        pr(f'| {lab} | {J} | {m/L:.3f} | {beta:.3f} | {Fn:.3e} | {Fc:.3e} | {tn["Qtol"]:.2e} / {tc["Qtol"]:.2e} | '
           f'{tn["Lam_mono"]:.0f} / {tc["Lam_mono"]:.0f} | {tn["Lam_vec_R1p"]:.0f} / {tc["Lam_vec_R1p"]:.0f} | '
           f'{tc["Ksh1"]:.2e} | {tc["KC2"]:.2e} |')

    pr('\n(iii) hybrid (flat NF02 partial wave + curvature, ls_core.hyb_fast) vs F_LR, (hybrid - F_LR)/F_LR  [NE]')
    hy = []
    for mf in (1.0, 1.3):
        for J in (139, 175, 223, 276, 344, 500, 700):
            m = C.mgap(1000, J)*mf; mfl = m/L; u = 1/mfl**2
            row = f'   J = {J:4d}, m = {mf} m_gap (m/Lam = {mfl:.3f}):'
            for h in ('++', '+-'):
                hb = LS.hyb_fast(1000, J, mfl, h)
                F = C.LR_action(hb['beta'], u, C.SIGMA[h])
                d = (hb['base'] - F)/F
                hy.append(dict(J=J, mf=mf, h=h, rel=d, beta=hb['beta']))
                row += f'  {h}: beta = {hb["beta"]:.4f}, rel = {d:+.3e}'
            pr(row)
    json.dump(dict(tol=res, hyb=hy), open(os.path.join(CACHE, 'lrt5_tolerances.json'), 'w'), default=float, indent=0)
    open(os.path.join(CACHE, 'lrt5_tolerances.txt'), 'w').write(OUT.getvalue())


if __name__ == '__main__':
    main()

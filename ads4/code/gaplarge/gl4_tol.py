"""
gaplarge/gl4_tol.py (GAP_LARGESPIN.md §2.1): Theorem L' per-state tolerances (PW_IMAGES / LRINF_TIGHT §5 definitions, re-typed from
lrtight/lrt5_tolerances.tolerances and pwimg/pw5_theoremL) as functions of Delta_g, at states with the SAME gap-unit coordinates
(beta, mu = m/Lam).  [DV] (floating point).

Per gap Dg in {1000, 2000, 5000}:
  * Lam = sqrt((Dg-2)(Dg-1)); grid nu_k = (k + 1/2) 0.04 < Lam; PS = psi2(nu/Lam) (exact NF02 polynomial; Lam-independent in p);
  * D(Dg) = Re sum X A[1/c] dnu/Lam - sum psi2/c dnu/Lam  (pw3 definition at t = 0, P_0 = 1), X the k = 2 Abel weight re-solved here
    with deshift/d_abel.Abel (read-only import, its main() not called) on dnu = 0.04, nu_max = Lam + 40;  checked against PW4's D(1000);
  * floors: 'num' = min_sigma F_LR(beta, u, sigma) (lsl_core.LR_action, read-only, at this Dg);
            'unif' = f_U(beta-cell)/Lam^2, the CERTIFIED Dg-uniform gap-unit floor of gl2 (valid for every Dg >= 1000);
  * P_{-1/2+i nu}(cosh beta) by Mehler with nphi = max(600 + 600 beta, 1.5 Lam beta + 200) (lrt5's 600 + 600 beta under-resolves
    nu beta at Lam = 5000); states with beta > 3.1 (J1000 = 1e5, 1e6) only at Dg = 1000 (cost);
  * Q_tol = F/(4|D|);  Lam_amp (monopole) = Q_tol/[c2 (|kA| + |kB|/4) eta/m^2];  Lam_amp (scalar-shaped vector, c-ratio form R1');
    K_tol for sh+-1 and C2 (R1'), exactly as lrt5.
States: the LRINF_TIGHT §5 rows mapped to the same (beta, mu): m = mu Lam, J = round(mu Lam sinh(beta/2) - 1/2) (beta recomputed),
plus the twist-edge J*(Dg) of the Dg-uniform threshold beta0 = 0.37.
Output: cache/gl4_tol.{txt,json}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
CODE = os.path.abspath(os.path.join(HERE, '..'))
for sub in ('deshift', 'n2sub', 'lsl'):
    sys.path.insert(0, os.path.join(CODE, sub))
from fractions import Fraction as Fr
import numpy as np
import d_abel as DA          # read-only
import n2_common as N        # read-only
import lsl_core as C         # read-only

OUT = []
def pr(s=''):
    print(s, flush=True); OUT.append(s)


def cfac(mu):
    m2 = 4*mu*mu
    return 32*(m2 + 1)*(m2 + 9)/((m2 + 25)*(m2 + 81))


CS = {(1, 1): lambda v: 32*(16*v**4 + 120*v**2 + 119)/((4*v**2 + 49)*(4*v**2 + 121)),
      (0, 2): lambda v: 32*(4*v**2 + 9)*(4*v**2 + 49)/((4*v**2 + 81)*(4*v**2 + 169))}      # pw4_theoremW.CS (copied)
KA, KB = 86.75/6 - 8.0, 11/6                                                                   # pw5_theoremL (copied)
DNU = 0.04
G2 = json.load(open(os.path.join(HERE, 'cache', 'gl2_allgaps.json')))
FL = [(Fr(f['a']), Fr(f['b']), f['fU']) for f in G2['floor']]
LRT3 = json.load(open(os.path.join(CODE, 'lrtight', 'cache', 'lrt3_rcert.json')))['rows']


def fU(beta):
    for a, b, v in FL:
        if float(a) <= beta <= float(b):
            return v
    return None


def Fcert1000(beta):
    for r in LRT3:
        if float(Fr(r['ba'])) <= beta <= float(Fr(r['bb'])):
            return r['Fcert']
    return None


class Gap:
    def __init__(self, Dg):
        self.Dg = Dg; self.L = N.lam(Dg)
        t0 = time.time()
        A = DA.Abel(Dg, DNU, self.L + 40.0)
        Y, res, it = A.solve()
        X = A.a*Y
        sel = A.nu < self.L
        nu = A.nu[sel]; X = X[sel]; ps = A.phi[sel]
        raw = (X*N.Aop(lambda v: 1/cfac(v), nu)).real.sum()*DNU/self.L
        ds = (ps/cfac(nu)).real.sum()*DNU/self.L
        self.D = raw - ds
        self.NU, self.PS = nu, ps
        self.info = f'Abel solve N = {A.N}, CGLS its {it}, residual {res:.1e}, {time.time() - t0:.0f}s'

    def Pd(self, beta, s=0.0):
        return N.Pnu(beta, self.NU - 1j*s, nphi=int(max(600 + 600*beta, 1.5*self.L*beta + 200)))   # Mehler nodes resolve nu*beta up to Lam*beta

    def tol(self, J, m, F):
        L, NU, PS, D = self.L, self.NU, self.PS, self.D
        beta = C.beta_of(J, m); eta = np.cosh(beta)
        cst = {s: CS[s](0.5j).real for s in CS}
        c2 = max(abs(cst[(1, 1)]), abs(cst[(0, 2)]))
        Qtol = F/(4*abs(D))
        Qest = c2*(abs(KA) + abs(KB)/4)*eta/m**2
        P0 = self.Pd(beta).real
        dP = (self.Pd(beta*(1 + 1e-4)).real - self.Pd(beta*(1 - 1e-4)).real)/(2e-4*beta)
        sh1 = self.Pd(beta, 1.0).real
        S = {'A': eta*P0, 'B': eta*NU**2*P0, 'C1': dP/np.sinh(beta), 'C2': np.sinh(beta)*NU**2*dP, 'C3': np.sinh(beta)*dP, 'sh1': sh1}
        Sh = {'A': eta, 'B': -eta/4, 'C1': 0, 'C2': 0, 'C3': 0, 'sh1': (1 + eta)/2}
        K = {}
        for k in S:
            r1p = max(abs(2*((PS*(CS[s](NU)/cfac(NU))*S[k]).sum()*DNU/L + D*CS[s](0.5j).real*Sh[k])) for s in CS)/m**2
            K[k] = 0.5*F/r1p
        sc = {'A': abs(KA), 'B': abs(KB), 'C1': 1/6, 'C2': 1/6, 'C3': 11.75/6}
        vt = 1/sum(sc[k]/K[k] for k in sc)
        return dict(Qtol=Qtol, Lam_mono=Qtol/Qest, Lam_vec=vt, Ksh1=K['sh1'], KC2=K['C2'])


def jstar_edge(Dg, s0=0.37211):
    J = 0
    while (2*J + 1)**2 < s0**2*(Dg - 2)*(Dg + 2*J - 1):
        J += 1
    return J


ROWS = [('gap', 223), ('gap', 300), ('gap', 500), ('gap', 857), ('gap', 10**4), ('gap', 10**5), ('gap', 10**6),
        ('1.3 m_gap', 307), ('1.5 Lam', 279), ('2 Lam', 372), ('5 Lam', 929)]


def main():
    pr(__doc__.split('Output')[0])
    L0 = C.lam(1000)
    base = []
    for lab, J in ROWS:
        if lab == 'gap':
            m = C.mgap(1000, J)
        elif lab == '1.3 m_gap':
            m = 1.3*C.mgap(1000, J)
        else:
            m = float(lab.split()[0])*L0
        base.append((lab, J, C.beta_of(J, m), m/L0))
    res = {}
    for Dg in (1000, 2000, 5000):
        g = Gap(Dg)
        pr(f'\n## Dg = {Dg}: Lam = {g.L:.4f}; D = {g.D:+.7e}; Lam^2 D = {g.L**2*g.D:+.3f}  ({g.info})')
        if Dg == 1000:
            pr(f'   check: PW4 D(1000) = -9.9833211e-4 (published 0.04 grid); difference {g.D + 9.9833211e-4:+.1e}')
        pr('   state (1000 label, J1000) | J | m/Lam | beta | F num | F unif-cert (F_cert@1000) | Q_tol num/unif | Lam_amp mono num/unif | '
           'Lam_amp vec R1\' num/unif | K_tol sh1 unif | K_tol C2 unif')
        rows = []
        states = [s for s in base if Dg == 1000 or s[2] < 3.1] + [('twist edge J*(Dg)', None, None, None)]   # beta > 3.1: 1000 only (cost)
        for lab, J0, beta0, mu in states:
            if J0 is None:
                J = jstar_edge(Dg); m = np.sqrt((Dg - 2.0)*(Dg - 1.0 + 2*J))
            else:
                m = mu*g.L
                J = int(round(m*np.sinh(beta0/2) - 0.5))
            beta = C.beta_of(J, m); u = min(1.0, g.L**2/m**2)
            Fn = min(C.LR_action(beta, u, s, Dg) for s in (1.0, -1.0))
            fu = fU(beta)
            Fu = fu/g.L**2 if fu is not None else None
            tn = g.tol(J, m, Fn); tu = g.tol(J, m, Fu)
            fc = Fcert1000(beta) if Dg == 1000 else None
            rows.append(dict(label=lab, J0=J0, J=J, mu=m/g.L, beta=beta, Fn=Fn, Fu=Fu, num=tn, unif=tu))
            pr(f"   {lab:9s} {str(J0):>8s} | {J:8d} | {m/g.L:6.3f} | {beta:6.3f} | {Fn:.3e} | {Fu:.3e}"
               + (f" ({fc:.3e})" if fc else '') + f" | {tn['Qtol']:.2e}/{tu['Qtol']:.2e} | {tn['Lam_mono']:.0f}/{tu['Lam_mono']:.0f} | "
               f"{tn['Lam_vec']:.0f}/{tu['Lam_vec']:.0f} | {tu['Ksh1']:.2e} | {tu['KC2']:.2e}")
        res[Dg] = dict(Lam=g.L, D=g.D, rows=rows)
    pr('\n## scaling: tolerance(Dg)/tolerance(1000) at the same (beta, mu), unif-cert floor; (Lam/Lam0)^2 = '
       + ', '.join(f'{Dg}: {(res[Dg]["Lam"]/res[1000]["Lam"])**2:.3f}' for Dg in (2000, 5000)))
    for (lab, J0, b0_, _) in base:
        if b0_ >= 3.1:
            continue
        get = lambda Dg: next(r['unif'] for r in res[Dg]['rows'] if r['label'] == lab and r['J0'] == J0)
        r0 = get(1000)
        pr(f'   {lab:9s} {J0:>8d}: ' + '; '.join(
            f"Dg={Dg}: Q_tol x{get(Dg)['Qtol']/r0['Qtol']:.3f}, mono x{get(Dg)['Lam_mono']/r0['Lam_mono']:.2f}, "
            f"vec x{get(Dg)['Lam_vec']/r0['Lam_vec']:.2f}, sh1 x{get(Dg)['Ksh1']/r0['Ksh1']:.2f}, "
            f"C2 x{get(Dg)['KC2']/r0['KC2']:.2f}" for Dg in (2000, 5000)))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE, 'cache', 'gl4_tol.json'), 'w'), default=float, indent=0)
    open(os.path.join(HERE, 'cache', 'gl4_tol.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()

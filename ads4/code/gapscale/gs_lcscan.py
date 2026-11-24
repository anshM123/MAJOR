"""
gapscale/gs_lcscan.py (GAP_SCALING.md §2.2): where does the SIGN-FREE leading-Regge floor "certified flat value minus the Lemma-C
curvature bound" stay positive, as a function of Delta_g?
At fixed gap-unit state (J, u = Lam^2/m^2, h) the flat parts (B2, B3) do not depend on Lam (exact) and the Lemma-C bound LC decreases
monotonically in Lam [P] (gs_state.py docstring).  Hence lb(J, u; Lam) = (B2 + B3 - LC)/(N_f + LC) is increasing in Lam, and the set
{(J, u): lb > 0} at Delta_g = 1000 is contained in the same set at every larger gap.  This script tabulates, per m/Lam and helicity,
J_LC(Lam) = the largest J such that lb > 0 for EVERY J' <= J (spins allowed for the helicity), at several gaps, using the flat parts of
omegau/nf_common.flat_parts_nf through spinfunc/q3_analysis (read-only), as ee1000_common.hyb does.
(The flat values themselves are quadrature numbers here; the certificate CERT_CLOSE gives flat >= N_f/5 for every state, so replacing
B2 + B3 by N_f/5 gives the fully rigorous variant lb_cert = (N_f/5 - LC)/(N_f + LC), also tabulated.)
Usage: python3 gs_lcscan.py -> cache/lcscan.txt
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import ee1000_common as E
Q = E.Q
JMAX = 600
DGS = (300, 1000, 2000, 5000, 20000, 100000)
MFS = (1.0, 1.02, 1.05, 1.1, 1.3, 1.5, 2.0, 3.0, 5.0)
OUT = []


def pr(s=''):
    print(s, flush=True); OUT.append(s)


def main():
    t0 = time.time()
    F = Q.NearGap(1000).F                 # NF02 weights (Lam-independent functions of p)
    pf, wf = Q.gl01(2*JMAX + 200)
    pr("# gs_lcscan: J_LC = largest J with sign-free (flat - LemmaC) > 0 for all J' <= J  [quadrature flat | certified N_f/5 floor]")
    pr("# m/Lam   h  | " + " | ".join(f"Dg={D}" for D in DGS))
    for mfl in MFS:
        u = 1/mfl**2
        for h in ('++', '+-'):
            Js, B2, B3, _ = Q.flat_parts_nf(F, h, u, JMAX, pf, wf)
            p2 = F.psi2(pf)*wf; p3 = F.psi3(pf)*wf
            row = []
            for Dg in DGS:
                Lm = E.lam_of(Dg); m = mfl*Lm
                jl_q = jl_c = None
                for i, J in enumerate(Js):
                    J = int(J)
                    if h == '+-' and J < 4:
                        continue
                    if h == '++' and J % 2 == 1:
                        continue
                    B = (2*J + 1)/m; beta = 2*np.arcsinh(B/2)
                    lc = (1 - 1/np.cosh(beta/2)) + np.minimum(2.0, pf*Lm*(B - beta))
                    LC = (np.abs(p2)*2*lc).sum() + (np.abs(p3)*u*lc).sum()
                    Nf = abs(B2[i]) + abs(B3[i])
                    lbq = (B2[i] + B3[i] - LC)/(Nf + LC); lbc = (Nf/5 - LC)/(Nf + LC)
                    if jl_q is None and lbq <= 0:
                        jl_q = J - 1
                    if jl_c is None and lbc <= 0:
                        jl_c = J - 1
                    if jl_q is not None and jl_c is not None:
                        break
                row.append(f"{'>=%d' % JMAX if jl_q is None else jl_q:>6} {'>=%d' % JMAX if jl_c is None else jl_c:>6}")
            pr(f"{mfl:5.2f}   {h} | " + " | ".join(row))
    pr(f"# fit: J_LC ~ c Lam^(1/3) expected (LC/N_f ~ J^6/Lam^2 near the gap).  [{time.time()-t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'lcscan.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()

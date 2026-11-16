"""
End-to-end test at k = 2, 3 (and k = 4, 6 as controls): tensor Regge moments of the exact spinning partial waves of <O phi O phi>
(O = J: l = 1, O = T: l = 2; Dphi = 1.75) with the corner completed (corner.TensorMomentsLC), against the derived leading heavy
formula (ads4/code/e2e/predict.py; SPINNING_E2E.md sect. 4.5).  Same (Delta, J) ladder as e2e (cases.cases_fixed_beta).
Helicity components from the two transverse modes (e2e/t7_corrections.measured_components).

Error of every quoted R (column 'err'): propagated through the 2x2 component solve from the relative errors of the modes (l,-l), (l,l):
   err(mode) = err_lc (|A - B| of the two LC settings, summed over the LC points) + err_drop (power-law extrapolation of the
               points still missing, counted entirely) + err_rad (radial-series truncation at the converged points: the e2e mask
               criterion X - 10 log X > 25 bounds it by e^{-25} X^p ~ 1e-9 per point; taken as 1e-9)
               + err_grid (quadrature: from q1_grid.py, see SPINNING_LOWK.md; not included here).
Usage: python3 t1_lowk_e2e.py l  [build beta]       Output: out_t1_lowk_e2e_l{l}.txt, data_t1_l{l}.pkl (in cache/)
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
from corner import TensorMomentsLC, CACHE, e2e_file, refined_base
from moments import TensorMoments
from predict import components, scalar_prefactor
from formulas import Fnu_mp
from t4_spinning_e2e import main_cases, moments_of, DPHI
from t7_corrections import measured_components
from cases import info
mp.mp.dps = 25
KS = (2, 3, 4, 6); NUS = (0.0, 1.0, 3.0)
ERR_RAD = 1e-9

def build_ref(l, Delta, J):
    """corner-completed moments on the corner-refined grid (corner.refined_base): all nodes with u < 0.03 from the LC
    representation, the e2e radial-series data above (the e2e grid under-resolves the log-oscillating kernel P_{-1/2+i nu}(cosh xi)
    near u = 0 at k = 2; out_q1_grid.txt, out_v1_scalar.txt)."""
    fn = os.path.join(CACHE, f"lcmodes_ref_l{l}_D{Delta:.6f}_J{J}.npz")
    if os.path.exists(fn):
        return TensorMomentsLC.load(fn)
    return TensorMomentsLC.build(l, Delta, J, base=refined_base(l, Delta, J), fn=fn)

def components_with_error(l, Delta, J, k, nu, res, e1, e2):
    Rp, Rf, worst, beta = measured_components(l, Delta, J, DPHI, k, nu, res)
    pref, _ = scalar_prefactor(l, Delta, J, DPHI, k, nu)
    X0 = Fnu_mp(nu, 0, beta); Y0 = Fnu_mp(nu, 2*l, beta)
    A1, B1 = components(l, k, nu, l, -l); A2, B2 = components(l, k, nu, l, l)
    E1 = abs(res[(k, nu, l, -l)]/pref); E2 = abs(res[(k, nu, l, l)]/pref)
    det = abs(A1*B2 - A2*B1)
    dX = (e1*E1*abs(B2) + e2*E2*abs(B1))/det; dY = (e1*E1*abs(A2) + e2*E2*abs(A1))/det
    return float(mp.re(Rp)), float(mp.re(Rf)), float(dX/abs(X0)), float(dY/abs(Y0)), float(worst)

if __name__ == "__main__":
    l = int(sys.argv[1])
    if len(sys.argv) > 3 and sys.argv[2] == 'build':
        # build the corner-completed caches for one beta only (to run several betas in parallel); no output file
        beta = float(sys.argv[3])
        for c in main_cases():
            if abs(float(info(*c)[1]) - beta) < 0.15 and os.path.exists(e2e_file(l, *c)):
                build_ref(l, *c)
        sys.exit(0)
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    pr(f"# l = {l}, Dphi = {DPHI}: R_pres, R_flip = measured/predicted (helicity preserved / flipped at the spinning pair), corner completed.")
    pr("# err: propagated numerical error (see docstring); old: e2e value without corner completion (radial series, masked); lc: share of the")
    pr("# (l,-l) / (l,l) mode coming from the LC points; other: max deviation of the remaining modes from the two-component form.")
    data = {}
    for c in main_cases():
        Delta, J = c
        if not os.path.exists(e2e_file(l, Delta, J)):
            continue
        t0 = time.time()
        mom = build_ref(l, Delta, J)
        old = TensorMoments.load(e2e_file(l, Delta, J))
        m, b_, eta = info(Delta, J)
        res = moments_of(mom, KS, NUS); res_old = moments_of(old, KS, NUS)
        pr(f"# (Delta, J) = ({Delta:.4f}, {J}), m = {float(m):.2f}, beta = {float(b_):.3f}   [{time.time()-t0:.0f} s]   LC info: {mom.info if isinstance(mom.info, str) else repr(mom.info)}"[:400])
        for k in KS:
            for nu in NUS:
                v1, lc1, dr1, sh1 = mom.moment_err(k, nu, l, -l)
                v2, lc2, dr2, sh2 = mom.moment_err(k, nu, l, l)
                e1 = lc1 + dr1 + ERR_RAD; e2 = lc2 + dr2 + ERR_RAD
                Rp, Rf, dRp, dRf, worst = components_with_error(l, Delta, J, k, nu, res, e1, e2)
                Rpo, Rfo, _, _ = measured_components(l, Delta, J, DPHI, k, nu, res_old)
                data[(round(float(b_), 1), k, nu, c)] = dict(m=float(m), beta=float(b_), Rp=Rp, Rf=Rf, dRp=dRp, dRf=dRf, Rp_old=float(mp.re(Rpo)),
                                                             Rf_old=float(mp.re(Rfo)), lc=(sh1, sh2), err_modes=(e1, e2), other=worst,
                                                             drop=(dr1, dr2), errlc=(lc1, lc2))
                pr(f"   k={k} nu={nu:3.1f}:  R_pres {Rp:.9f} +- {dRp:.1e}   R_flip {Rf:.9f} +- {dRf:.1e}   | old {float(mp.re(Rpo)):.6f} / {float(mp.re(Rfo)):.6f}"
                   f"   | lc {sh1:.1e}/{sh2:.1e}  drop {dr1:.0e}/{dr2:.0e}  other {worst:.1e}")
    os.makedirs(CACHE, exist_ok=True)
    pickle.dump(data, open(os.path.join(CACHE, f'data_t1_l{l}.pkl'), 'wb'))
    open(os.path.join(HERE, f'out_t1_lowk_e2e_l{l}.txt'), 'w').write("\n".join(out) + "\n")

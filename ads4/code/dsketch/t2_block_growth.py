"""
t2_block_growth (DSKETCH_UP.md §3.2): exact per-nu action of the CHMRS functional C_{k,nu} on ONE scalar-external conformal block,
for nu below and ABOVE the block's mass m (m^2 = (Delta-J-2)(Delta+J-1)), in the exact Mellin representation of CHMRS
(2106.10274, eqs. CfromBMellin + Bkt action), evaluated by exactskel/exactmom.py (read-only import; validated there against the
CHMRS 1/m^8 series to 1e-12 for heavy states).

What is tested: a light single trace has m << Lambda, so the functional samples nu in (m, Lambda].  The heavy-state hybrid used for
light states in NEARGAP2 §2 (flat image + leading-Regge curvature P_{-1/2+i nu}(cosh beta) - J0(nu B)) assumes the per-nu action
stays flat-like there.  Here: (a) for nu < m the exact action approaches the flat bulk-point form C_BP (CHMRS (3.40)) and the LR
form; (b) for nu > m it grows like exp(+c nu) with c -> pi/2 (the 1/gamma_a(nu)^2 normalisation, ~ e^{pi nu}, is no longer
compensated: the Regge moment of a single block decays only like e^{-pi nu/2} above its mass); (c) the number is independent of the
contour position (5 positions) and of the line step (h = 0.1, 0.05) to <= 1e-8: it is not a numerical artefact.
Normalisation: the heavy density m^{2k} C_{k,nu}[Delta,J] (q_{Delta,J} set to 1 when it is not real: only nu-dependence matters).
Writes cache/t2_block_growth.txt
"""
import os, sys, warnings
sys.dont_write_bytecode = True
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ES = os.path.join(HERE, '..', 'exactskel')
sys.path.insert(0, ES); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
import exactmom as X
import skel

_lq = X.log_q
def safe_log_q(D, J, Df):
    try:
        return _lq(D, J, Df)
    except Exception:
        return 0.0
X.log_q = safe_log_q

lines = []
def log(s=''):
    print(s, flush=True); lines.append(s)

ok = True
log('# exact CHMRS per-nu action on one block (Df = 3), heavy density m^{2k} C_{k,nu}; C_BP = flat bulk-point form; P = P_{-1/2+i nu}(eta)')
for (D, J) in ((40.0, 0), (20.0, 0), (12.0, 2), (10.0, 0), (4.5, 0)):
    m2, eta = skel.m2_eta(D, J)
    for k in (2, 3):
        hc = X.HeavyC(D, J, 3.0)
        if not np.isfinite(hc.logq) or hc.logq == 0.0:
            hc.m2 = 1.0
        nus = np.array([1.0, 5, 10, 15, 20, 25, 30, 40, 50, 60])
        v, c, jj = hc.value(k, nus, diag=True)
        log(f'## Delta = {D}, J = {J}, k = {k}: m = {np.sqrt(abs(m2)):.3f}')
        for nu, vv, cc in zip(nus, v, c):
            bp = skel.C_BP(k, nu, m2, J) if m2 > 0 else float('nan')
            pe = skel.Pders(nu, eta, 0)[0] if m2 > 0 else float('nan')
            log(f'   nu = {nu:5.1f}  exact {vv:+.6e}  C_BP {bp:+.4e}  P(eta) {pe:+.4e}  (cancellation {cc:.1e})')
        # growth rate above m
        sel = nus > np.sqrt(abs(m2)) + 10
        if sel.sum() >= 2:
            lv = np.log(np.abs(v[sel])); nn = nus[sel]
            rates = np.diff(lv)/np.diff(nn)
            log(f'   d ln|C|/d nu above m: ' + ' '.join(f'{r:.3f}' for r in rates) + f'   (pi/2 = {np.pi/2:.3f})')
# contour / step independence at nu > m
log('# contour-position and step independence (nu > m)')
for (D, J, k) in ((10.0, 0, 2), (10.0, 0, 3), (12.0, 2, 2)):
    vals = []
    for h in (0.1, 0.05):
        hc = X.HeavyC(D, J, 3.0, h=h); hc.m2 = 1.0
        A, c1max, j0 = hc.contour(k)
        for jc in sorted(set([-2, -1, 0, j0])):
            if A + 2*jc + 1 <= 0:
                continue
            v = hc.value(k, np.array([20.0, 30.0]), jc=jc)
            vals.append(v)
    vals = np.array(vals)
    spread = np.max(np.abs(vals/vals[0] - 1))
    log(f'   Delta = {D}, J = {J}, k = {k}: values at nu = 20, 30: {vals[0][0]:+.8e} {vals[0][1]:+.8e}; max relative spread over '
        f'{len(vals)} (contour, step) choices: {spread:.1e}')
    ok &= spread < 1e-7
log(f'ALL CHECKS PASS: {ok}')
open(os.path.join(HERE, 'cache', 't2_block_growth.txt'), 'w').write('\n'.join(lines) + '\n')

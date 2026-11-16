# Referee mutation tests on cached mode data: does the c1 fit detect O(1/m) errors in the prediction?
import sys, os, numpy as np, mpmath as mp
sys.path.insert(0, sys.argv[1]); os.chdir(sys.argv[1])
from t4_spinning_e2e import main_cases, get_moments, moments_of, fname, DPHI
from t7_corrections import measured_components
from formulas import log_b
from cases import info
mp.mp.dps = 25
l = int(sys.argv[2])
muts = {
 'none':      lambda D, J: 1,
 'b(0)':      lambda D, J: mp.exp(log_b(D, J, 0.0) - log_b(D, J, 1 - DPHI)),        # R_mut = R * pred/pred_mut
 'b(l+1-Dphi)': lambda D, J: mp.exp(log_b(D, J, l + 1 - DPHI) - log_b(D, J, 1 - DPHI)),
 'b(1-Dphi)*(1+1/m)': lambda D, J: 1 + 1/info(D, J)[0],
 'x(1+3/m^2)': lambda D, J: 1 + 3/info(D, J)[0]**2,
 'shadow 2x': lambda D, J: 2,
}
for beta in (0.4, 1.0):
    cs = [c for c in main_cases() if abs(float(info(*c)[1]) - beta) < 0.05 and os.path.exists(fname(l, *c))]
    data = {}
    for c in cs:
        res = moments_of(get_moments(l, *c), (6,), (0.0, 1.0, 3.0))
        for nu in (0.0, 1.0, 3.0):
            Rp, Rf, w, bb = measured_components(l, c[0], c[1], DPHI, 6, nu, res)
            data[(c, nu)] = (float(info(*c)[0]), float(mp.re(Rp)), float(mp.re(Rf)))
    for nu in (0.0, 1.0, 3.0):
        for name, f in muts.items():
            ms = np.array([data[(c, nu)][0] for c in cs])
            out = []
            for idx in (1, 2):
                R = np.array([data[(c, nu)][idx]*float(f(*c)) for c in cs]) - 1
                cA = np.linalg.solve(np.stack([1/ms[-3:], 1/ms[-3:]**2, 1/ms[-3:]**4], 1), R[-3:])
                # 4-parameter fit with c3 included (four largest m)
                cC = np.linalg.solve(np.stack([1/ms[-4:], 1/ms[-4:]**2, 1/ms[-4:]**3, 1/ms[-4:]**4], 1), R[-4:])
                out.append(f"R(m~100)={R[-1]+1:.5f} c1(A)={cA[0]:+.4f} c2(A)={cA[1]:8.2f} c1(c1..c4)={cC[0]:+.4f}")
            print(f"l={l} beta={beta} nu={nu} {name:18s} PRES {out[0]} | FLIP {out[1]}", flush=True)

"""
Spin-dependent part of the O(1/m^2) correction (task 4): c2(spinning, helicity preserved / flipped) - c2(scalar skeleton with the same
dimensions (l+1, Dphi), the same k, nu, beta), compared with
  (i)  the exact normalisation polynomials of the conserved structures,   pol = -kappa_l [1/sinh^2(beta/2) + 1/cosh^2(beta/2)],
       kappa_1 = 1/4, kappa_2 = 5/2    [(J+1/2)^{2l} -> J(J+1) resp. (J-1)J(J+1)(J+2), and the same for Delta - 3/2];
  (ii) the helicity-dependent term of flat-space partial waves in its natural AdS form,   hel = -2 l^2 P'_{-1/2+i nu}(eta)/P_{-1/2+i nu}(eta)
       (-> 2 l^2 (nu/beta) J_1(nu beta)/J_0(nu beta) for beta -> 0;  REGGE_MOMENT_TTTT.md sect. 6.2 with lambda = l).
Usage: python3 t7c_spin_dependence.py     (needs data_t7_l*.pkl, data_t7b.pkl)
"""
import pickle, os, numpy as np, mpmath as mp
from formulas import Pnu
mp.mp.dps = 20
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
sb = pickle.load(open('data_t7b.pkl', 'rb'))
for l in (1, 2):
    if not os.path.exists(f'data_t7_l{l}.pkl'): continue
    s7 = pickle.load(open(f'data_t7_l{l}.pkl', 'rb'))
    kap = 0.25 if l == 1 else 2.5
    pr(f"# l = {l}:   D_pres = c2_pres - c2_scalar,  D_flip = c2_flip - c2_scalar")
    pr("#  beta   k   nu    c2_pres    c2_flip   c2_scalar |   D_pres   pol+hel  (pol, hel) |   D_flip")
    for (beta, k, nu), v in sorted(s7.items()):
        if k == 2: continue
        cs = sb[(l, beta, k, nu)]
        eta = mp.cosh(beta)
        P = Pnu(nu, eta); dP = mp.diff(lambda x: Pnu(nu, x), eta)
        pol = -kap*(1/np.sinh(beta/2)**2 + 1/np.cosh(beta/2)**2)
        hel = float(mp.re(-2*l*l*dP/P))
        c2p = v[5] if v[5] == v[5] else v[0]; c2f = v[7] if v[7] == v[7] else v[1]     # three-parameter fit (B) where available
        pr(f"  {beta:4.1f}   {k}  {nu:4.1f}  {c2p:9.2f}  {c2f:9.2f}  {cs:9.2f}  | {c2p-cs:8.2f}  {pol+hel:8.2f}  ({pol:7.2f}, {hel:7.2f}) | {c2f-cs:8.2f}")
# Empirical functional form of the helicity-preserving spin-dependent correction (a fit of the FORM of the O(1/m^2) correction;
# no normalisation is involved):   D_pres - pol  =  a0 + a1 cosh(beta) + a2 P'/P + a3 sinh^2(beta) P'/P,   P' = dP/d eta.
pr("# least-squares fit  D_pres - pol = a0 + a1 cosh(beta) + a2 P'/P + a3 sinh^2(beta) P'/P   (all available (beta, nu) at fixed k);")
pr("#   c2_pres from the three-parameter fit (B) of t7 where available")
for l in (1, 2):
    if not os.path.exists(f'data_t7_l{l}.pkl'): continue
    s7 = pickle.load(open(f'data_t7_l{l}.pkl', 'rb'))
    kap = 0.25 if l == 1 else 2.5
    for k in (6, 4):
        rows = []
        for (beta, kk, nu), v in sorted(s7.items()):
            if kk != k: continue
            c2p = v[5] if v[5] == v[5] else v[0]
            eta = mp.cosh(beta)
            dl = float(mp.re(mp.diff(lambda x: Pnu(nu, x), eta)/Pnu(nu, eta)))
            pol = -kap*(1/np.sinh(beta/2)**2 + 1/np.cosh(beta/2)**2)
            rows.append((beta, nu, c2p - sb[(l, beta, k, nu)] - pol, dl))
        rows = np.array(rows)
        if len(rows) < 6:
            pr(f"  l = {l}, k = {k}: only {len(rows)} points, no fit"); continue
        b = rows[:, 0]
        A = np.stack([np.ones(len(b)), np.cosh(b), rows[:, 3], np.sinh(b)**2*rows[:, 3]], 1)
        c = np.linalg.lstsq(A, rows[:, 2], rcond=None)[0]
        pr(f"  l = {l}, k = {k} ({len(rows)} points):  a0 = {c[0]:.2f}, a1 = {c[1]:.2f}, a2 = {c[2]:.2f}, a3 = {c[3]:.2f};   residuals: " + " ".join(f"{x:+.2f}" for x in (A @ c - rows[:, 2])))
        # the same with the coefficient of pol left free (prediction: 1) and with beta <= 1
        polv = -kap*(1/np.sinh(b/2)**2 + 1/np.cosh(b/2)**2)
        A5 = np.concatenate([polv[:, None], A], 1)
        c5 = np.linalg.lstsq(A5, rows[:, 2] + polv, rcond=None)[0]
        pr(f"       with a free coefficient p of pol (prediction p = 1):  p = {c5[0]:.3f}, a0 = {c5[1]:.2f}, a1 = {c5[2]:.2f}, a2 = {c5[3]:.2f}, a3 = {c5[4]:.2f}")
        sel = b < 1.5
        c = np.linalg.lstsq(A[sel], rows[sel, 2], rcond=None)[0]
        pr(f"       beta <= 1 only ({sel.sum()} points):  a0 = {c[0]:.2f}, a1 = {c[1]:.2f}, a2 = {c[2]:.2f}, a3 = {c[3]:.2f};   residuals: " + " ".join(f"{x:+.2f}" for x in (A[sel] @ c - rows[sel, 2])))
# Per-beta fits (the four-function fit above on two beta values is equivalent to a separate fit c_beta + d_beta P'/P at each beta,
# extrapolated linearly in sinh^2 beta; the per-beta coefficients d_beta are the direct test of -2 l^2 P'/P).
pr("# per-beta fits at k = 6:  D_pres - pol = c_beta + d_beta P'/P  over nu = 0, 1, 3   (prediction of the flat-space term: d_beta = -2 l^2)")
for l in (1, 2):
    if not os.path.exists(f'data_t7_l{l}.pkl'): continue
    s7 = pickle.load(open(f'data_t7_l{l}.pkl', 'rb'))
    kap = 0.25 if l == 1 else 2.5
    for beta in (0.4, 1.0, 2.0):
        rows = []
        for nu in (0.0, 1.0, 3.0):
            v = s7.get((beta, 6, nu))
            if v is None: continue
            c2p = v[5] if v[5] == v[5] else v[0]
            eta = mp.cosh(beta); dl = float(mp.re(mp.diff(lambda x: Pnu(nu, x), eta)/Pnu(nu, eta)))
            pol = -kap*(1/np.sinh(beta/2)**2 + 1/np.cosh(beta/2)**2)
            rows.append((nu, c2p - sb[(l, beta, 6, nu)] - pol, dl))
        r = np.array(rows)
        A = np.stack([np.ones(len(r)), r[:, 2]], 1); c = np.linalg.lstsq(A, r[:, 1], rcond=None)[0]
        cc = np.mean(r[:, 1] + 2*l*l*r[:, 2])
        pr(f"  l = {l}, beta = {beta}:  d_beta = {c[1]:7.2f}  (prediction {-2*l*l}),  c_beta = {c[0]:7.2f};   residuals with d fixed to {-2*l*l}: " + " ".join(f"{x:+.2f}" for x in (cc - 2*l*l*r[:, 2] - r[:, 1])))
open('out_t7c_spin_dependence.txt', 'w').write("\n".join(out) + "\n")

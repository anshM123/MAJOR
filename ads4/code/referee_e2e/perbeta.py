import pickle, numpy as np, mpmath as mp, sys, os
sys.path.insert(0, os.getcwd())
from formulas import Pnu
sb = pickle.load(open('data_t7b.pkl', 'rb'))
for l in (1, 2):
    s7 = pickle.load(open(f'data_t7_l{l}.pkl', 'rb')); kap = 0.25 if l == 1 else 2.5
    for beta in (0.4, 1.0, 2.0):
        rows = []
        for nu in (0.0, 1.0, 3.0):
            v = s7.get((beta, 6, nu));
            if v is None: continue
            c2p = v[5] if v[5] == v[5] else v[0]
            eta = mp.cosh(beta); dl = float(mp.re(mp.diff(lambda x: Pnu(nu, x), eta)/Pnu(nu, eta)))
            pol = -kap*(1/np.sinh(beta/2)**2 + 1/np.cosh(beta/2)**2)
            rows.append((nu, c2p - sb[(l, beta, 6, nu)] - pol, dl))
        r = np.array(rows)
        A = np.stack([np.ones(len(r)), r[:, 2]], 1); c, res, *_ = np.linalg.lstsq(A, r[:, 1], rcond=None)
        print(f"l={l} beta={beta}: per-beta fit D-pol = c + d P'/P :  d = {c[1]:.2f} (pred -2l^2 = {-2*l*l}),  c = {c[0]:.2f}, resid = {np.round(A@c - r[:,1],2)}, P'/P = {np.round(r[:,2],3)}")
        # with d fixed to -2l^2
        cc = np.mean(r[:, 1] + 2*l*l*r[:, 2]); print(f"      d fixed to {-2*l*l}: residuals {np.round(cc - 2*l*l*r[:,2] - r[:,1],2)};  d = 0: residuals {np.round(np.mean(r[:,1]) - r[:,1],2)}")

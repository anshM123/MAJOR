"""
altcut.py (NU_TAIL_CHECK.md §4): the grid weights' tail beyond m on the J = 0 gap state (Delta = Dg, m = Lam) under alternative analytic
cutoffs theta(nu; nu_c, a), and its scaling with the floor; scalar proxy Gamma = m^4 C_{2,nu}/2 (Arb, as tail.py).
usage: python3 -B altcut.py Dg   -> cache/altcut_<Dg>.txt
"""
import os, sys, warnings
sys.dont_write_bytecode = True
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel')); sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import numpy as np
from scipy.special import erf
import skelexact3 as S3
DATA = os.path.join(HERE, '..', '..', 'prl', 'data')
Dg = int(sys.argv[1]); Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
out = []
def log(s=''):
    print(s, flush=True); out.append(s)
grids = {}
for tag, fn in (('Dgrid', f'd_abel_X_{Dg}.txt'), ('Qgrid', os.path.join('k2', f'xq_abel_X_{Dg}.txt'))):
    a = np.loadtxt(os.path.join(DATA, fn), comments='#')
    grids[tag] = (a[:, 0], a[:, 1] + 1j*a[:, 2], a[:, 5])
nu = grids['Dgrid'][0]; dn = nu[1] - nu[0]
D, J = float(Dg), 0
m = np.sqrt((D - 2)*(D - 1))
numax = nu[-1]
L = S3.Laplace3(D, J, 3.0, ks=(2,), prec=400, N0=240, nq=320, M=int(3*numax**2/(2*D) + 6*np.sqrt(numax + 1) + 60))
nb = np.linspace(0.02, Lam, 400); nt = np.arange(np.floor(Lam) - 2, numax + 0.5, 0.5)
Gb = np.array([float(L.value(2, v).mid())/2 for v in nb]); Gt = np.array([float(L.value(2, v).mid())/2 for v in nt])
g = np.where(nu <= nt[0], np.interp(nu, nb, Gb), np.exp(np.interp(nu, nt, np.log(np.abs(Gt)))))
log(f'# Dg = {Dg}, J = 0 gap state, m = Lam = {Lam:.4f}; grid end Lam{nu[-1]-Lam:+.1f}; tails = Re sum_(nu>m) dnu X theta Gamma')
for tag, (nu_, X, psi) in grids.items():
    inL = nu_ <= Lam
    Apsi = np.sum(psi[inL]*g[inL])*dn
    for (nc, a) in ((Lam + 100, 10), (Lam + 60, 6), (Lam + 40, 4), (Lam + 30, 3), (Lam + 20, 2), (Lam + 15, 1.5)):
        th = 0.5*(erf((nu_ + nc)/a) - erf((nu_ - nc)/a))
        T = np.real(np.sum(X[~inL]*th[~inL]*g[~inL]))*dn
        Ta = np.sum(np.abs(X[~inL])*th[~inL]*g[~inL])*dn
        leak = 0.5*(1 - erf((nc - Lam)/a))          # theta's leakage at the edge of X's support (HK's eps_c criterion)
        log(f'  {tag}: cutoff (Lam{nc-Lam:+.0f}, {a}): tail {T:+.3e}, |tail| {Ta:.3e} = {Ta/abs(Apsi):.1e} A_psi '
            f'(A_psi = {Apsi:.4e}); 1 - theta(Lam) = {leak:.1e}')
open(os.path.join(HERE, 'cache', f'altcut_{Dg}.txt'), 'w').write('\n'.join(out) + '\n')

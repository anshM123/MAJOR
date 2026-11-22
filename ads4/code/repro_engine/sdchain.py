"""repro_engine/sdchain.py -- from a state/reference ratio rho(nu) to Gamma_K(nu), s_d and FF/N (own code; definitions as read from
spinfunc/bg_qan.py, bf_qan.py and resultq/stage6_model.py, identical to REPRO.md's chain):
  T(nu) = sum_j R_j(nu) P_{-1/2 + i(nu - i j)}(cosh beta),  cosh beta = 1 + (2J+1)^2/(2 m^2),  m^2 = (D-J-2)(D+J-1);
  Gamma_exact = m^2 (c T_ref rho - T_st),  c: Gamma_exact(0.02) = 0;   Gamma_K = Gamma_exact - 0.75 sum_j R_j (nu - ij)^2 - const N^_4,
  Gamma_K(0.02) = 0;   s_d = 0.04 sum_nu chi_d(nu) Re Gamma_K(nu)   (H0 min-norm gauge, Stage 6 cache, read as data);
  FF/N = (FF - 2 a u^3 s_d - c (2u^3/Lam^2) s_v)/N  with the HYB FF, N of the binding state (REPRO q_binding numbers, re-derived there).
Data read (not code): resultq/cache/stage4_hj.pkl ['EE3'] (R_j as sympy strings), resultq/cache/stage6_gauge.pkl (chi), the engine caches
spinfunc/cache/bg_qan_D*_J0.pkl, rs2_eeq_D*_J0_K5*.pkl (for comparison only).
"""
import os, pickle
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
NU0 = 0.02


def legP(mu, eta, nmax=400):
    """P_{-1/2 + i mu}(eta) = 2F1(1/2 + i mu, 1/2 - i mu; 1; (1 - eta)/2) by its series (|1 - eta| << 1)."""
    mu = np.asarray(mu, complex); w = (1 - eta)/2
    term = np.ones_like(mu); tot = term.copy()
    for n in range(1, nmax):
        term = term*((n - 0.5)**2 + mu*mu)*w/(n*n)
        tot = tot + term
        if np.max(np.abs(term)) < 1e-19*np.max(np.abs(tot)):
            break
    return tot


_RJ = {}


def Rj(nus):
    if not _RJ:
        D = pickle.load(open(os.path.join(ROOT, 'resultq', 'cache', 'stage4_hj.pkl'), 'rb'))['EE3']
        nu = sp.Symbol('nu')
        for j, v in D.items():
            _RJ[j] = sp.lambdify(nu, sp.sympify(v), 'numpy')
    return {j: np.array(f(np.asarray(nus, complex)), complex)*np.ones(len(nus)) for j, f in _RJ.items()}


def mass(D, J=0):
    al = D - J - 2.0; m2 = al*(al + 2*J + 1)
    return np.sqrt(m2), 1 + (2*J + 1)**2/(2*m2)


def chi(gauge='H0'):
    G = pickle.load(open(os.path.join(ROOT, 'resultq', 'cache', 'stage6_gauge.pkl'), 'rb'))
    key = {'H0': 'H0 min-norm', 'H2': 'H2 (1+(nu/5)^2)', 'H4': 'H4 amp-L1'}[gauge]
    return G[('C4', 0.04, 'nu')], G[('C4v', 0.04)][key]['chi'], G[('C4d', 0.04)][key]['chi']


def gammaK(nus, rho, D, DREF=2e4):
    """nus must start at 0.02 (normalisation point)."""
    assert abs(nus[0] - NU0) < 1e-12
    m, eta = mass(D); _, etar = mass(DREF)
    R = Rj(nus)
    Tst = sum(R[j]*legP(nus - 1j*j, eta) for j in R)
    Tref = sum(R[j]*legP(nus - 1j*j, etar) for j in R)
    NH = sum(R.values()); sabs = sum(np.abs(R[j]) for j in R)
    tot = Tref*rho; c = Tst[0]/tot[0]
    Gex = m*m*(c*tot - Tst)
    GK = Gex - 0.75*sum(R[j]*(nus - 1j*j)**2 for j in R)
    GK = GK - (GK[0]/NH[0])*NH
    return dict(GK=GK, Gex=Gex, sabs=sabs, m=m, Tst=Tst, Tref=Tref)


def sd(nus, GK, gauge='H0', lo=0.0, hi=41.0):
    nuw, cv, cd = chi(gauge)
    k = min(len(nus), len(cd))
    assert np.allclose(nus[:k], nuw[:k])
    sel = (nus[:k] >= lo) & (nus[:k] < hi)
    return 0.04*np.sum(cd[:k][sel]*GK[:k][sel].real), 0.04*np.sum(cv[:k][sel]*GK[:k][sel].real)


def engine_ratio(D, K):
    """engine rho(nu) on the 0.04 grid: K <= 4 from the production bg_qan cache, K = 5 from the rs2 K = 5 images (own evaluation of
    I = exp(l0 + i nu l1) sum_{k<=K, j} D^-k (i nu/2)^j M[(k, j)])."""
    if K <= 4:
        d = pickle.load(open(os.path.join(ROOT, 'spinfunc', 'cache', 'bg_qan_D%g_J0.pkl' % D), 'rb'))
        r = d['rs'][K]
        return r['nus'], r['ratio']
    db = pickle.load(open(os.path.join(ROOT, 'spinfunc', 'cache', 'rs2_eeq_D%g_J0_K5%s.pkl' % (D, '' if K == 5 else 'full')), 'rb'))
    def val(rec, Kt):
        r = rec['res'][('EE4', 0, 0)]; nus = rec['nus']; Dl = rec['Delta']
        pre = np.exp(r['lconst'][0] + 1j*nus*r['lconst'][1])
        tot = np.zeros(len(nus), complex)
        for (k, j), v in r['M'].items():
            if k <= Kt:
                tot += Dl**(-k)*(0.5j*nus)**j*v
        return nus, pre*tot
    nus, Ist = val(db['st'], 5); _, Iref = val(db['ref'], db['ref']['K'])
    n = 1024
    return nus[:n], (Ist/Iref)[:n]


# HYB action of the binding state (++ J = 0, u = 1), Q3-full-delta0.2 -- values re-derived by REPRO (q_binding.py: FF, N, a, c)
def hyb_numbers(Dg):
    import json, importlib.util
    spec = importlib.util.spec_from_file_location('qb', os.path.join(ROOT, 'repro', 'q_binding.py'))
    qb = importlib.util.module_from_spec(spec); spec.loader.exec_module(qb)
    H = qb.hyb_J0(Dg)
    return dict(FF=H['FF'], N=H['N'], Lam=H['Lam'], a=qb.a, c=qb.c)


def ffn(Dg, sd_, sv_, H=None):
    H = H or hyb_numbers(Dg)
    return (H['FF'] - 2*H['a']*sd_ - H['c']*(2/H['Lam']**2)*sv_)/H['N'], H

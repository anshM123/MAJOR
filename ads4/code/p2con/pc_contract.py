"""
p2con/pc_contract.py (P2_CONTRACT.md §1.3): contraction of block tensors (pc_sym) with the exact leading-Regge EE engine (pc_engine).

Side tables (beta-INDEPENDENT, computed once per nu set): A_mode[nu, k2, k4, aX, a, k], B_mode[nu, k1, k3, aY, b, k'] in the
polarisation-mode basis of each leg (rest frames of p0, pb0) and rest-frame insertion monomials.  For a block tensor
{(k1, k2, k3, k4): {(aX, aY): c}} the raw (eta-kernel) image is
    raw(nu) = sum c sum_{a b} G^{ab}(beta) sum_{k k'} A_mode[.., a, k] B_mode[.., b, k'] M_{k k'}(nu, beta).
Normalisation: the engine's raw image of the helicity configuration (z2', z4', z1, z3) = (eps-, eps+, eps+, eps-) equals
2 pi (R_- P_{l-1} + R_+ P_{l+1}) with the PUBLISHED R_-+ (kernels.json) exactly (checked, pc_checks), so the published normalised image
is raw / (2 pi n_h (R_- + R_+)) with n_h the class normalisation of the leading tensor.
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'n2sub'))
import numpy as np
import pc_engine as E
import n2_common as N

MODES = E.MODES
WM = np.array([E.mode_weights(k) for k in MODES])           # (mode, grid)


def tables(nus, alX, alY, spin=True, chunk=40):
    nus = np.atleast_1d(np.asarray(nus, float))
    TA = []; TB = []; tails = []
    for i0 in range(0, len(nus), chunk):
        nn = nus[i0:i0 + chunk]
        a, ta = E.side_table('A', spin, nn, E.POLPAIRS, alX, [0, 1, 2])
        b, tb = E.side_table('B', spin, nn, E.POLPAIRS, alY, [0, 1, 2])
        TA.append(a); TB.append(b); tails.append(max(ta, tb))
    TA = np.concatenate(TA); TB = np.concatenate(TB)
    n = len(nus)
    TA = TA.reshape(n, 5, 5, len(alX), 3, -1); TB = TB.reshape(n, 5, 5, len(alY), 3, -1)
    TA = np.einsum('ka,lb,nab...->nkl...', WM, WM, TA)
    TB = np.einsum('ka,lb,nab...->nkl...', WM, WM, TB)
    return dict(nus=nus, A=TA, B=TB, alX=list(alX), alY=list(alY), tail=max(tails))


def nnz_list(tensor, tab):
    iX = {a: i for i, a in enumerate(tab['alX'])}; iY = {a: i for i, a in enumerate(tab['alY'])}
    I, J, C = [], [], []
    nX, nY = len(tab['alX']), len(tab['alY'])
    for (k1, k2, k3, k4), d in tensor.items():
        for (aX, aY), c in d.items():
            I.append(((k2 + 2)*5 + (k4 + 2))*nX + iX[aX]); J.append(((k1 + 2)*5 + (k3 + 2))*nY + iY[aY]); C.append(c)
    return np.array(I, int), np.array(J, int), np.array(C, complex)


def TG(tensor, tab, beta, sel=None):
    """T_G[nu, k, k'] = sum_ab G^ab sum_nnz c A[.., a, k] B[.., b, k']."""
    G = E.gradient_metric(beta)
    A = tab['A'] if sel is None else tab['A'][sel]
    B = tab['B'] if sel is None else tab['B'][sel]
    n = A.shape[0]
    A = A.reshape(n, -1, 3, A.shape[-1]); B = B.reshape(n, -1, 3, B.shape[-1])
    Bt = np.einsum('ab,nibk->niak', G, B)
    I, J, C = nnz_list(tensor, tab)
    out = np.zeros((n, A.shape[-1], B.shape[-1]), complex)
    if len(C) == 0:
        return out
    step = max(1, int(4e7/(n*3*A.shape[-1])))
    for s0 in range(0, len(C), step):
        ii, jj, cc = I[s0:s0 + step], J[s0:s0 + step], C[s0:s0 + step]
        out += np.einsum('niak,i,nial->nkl', A[:, ii], cc, Bt[:, jj], optimize=True)
    return out


def Mpair(nus, beta, thr=12.0):
    """M[nu, k, k'] as numbers: quadrature for nu beta < thr, closed form (P, P' reduced) above."""
    nus = np.atleast_1d(np.asarray(nus, float))
    out = np.zeros((len(nus), 2*E.KH + 1, 2*E.KH + 1), complex)
    lo = nus*beta < thr
    if np.any(lo):
        idx = np.where(lo)[0]
        for i0 in range(0, len(idx), 200):
            ii = idx[i0:i0 + 200]
            out[ii] = E.M_quad(nus[ii], beta)
    if np.any(~lo):
        idx = np.where(~lo)[0]
        for i0 in range(0, len(idx), 400):
            ii = idx[i0:i0 + 400]
            cP, cD = E.M_closed(nus[ii], beta)
            P, dP = E.legendre_PdP(nus[ii], beta)
            out[ii] = cP*P[:, None, None] + cD*dP[:, None, None]
    return out


def Mpair_PD(nus, beta):
    """closed-form (coefficient of P_l, coefficient of P_l') pairing (for asymptotics; stable for nu beta >> 1)."""
    return E.M_closed(np.atleast_1d(np.asarray(nus, float)), beta)


def raw_image(tensor, tab, beta, M=None, sel=None):
    T = TG(tensor, tab, beta, sel)
    nus = tab['nus'] if sel is None else tab['nus'][sel]
    if M is None:
        M = Mpair(nus, beta)
    return np.einsum('nkl,nkl->n', T, M)


def norm_hpm(nus):
    """2 pi (R_- + R_+)(nu) (published R, kernels.json): the engine's forward normalisation of the unit helicity configuration."""
    Rm, Rp = N.Rk()
    nus = np.asarray(nus, complex)
    return 2*np.pi*(Rm(nus) + Rp(nus))


def lead_shape(nus, beta):
    """2 pi (R_- P_{l-1} + R_+ P_{l+1}) on nus."""
    Rm, Rp = N.Rk()
    nus = np.asarray(nus, float)
    return 2*np.pi*(Rm(nus + 0j)*N.Pnu(beta, nus + 1j) + Rp(nus + 0j)*N.Pnu(beta, nus - 1j))

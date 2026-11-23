"""
dc_fit (DIRCV.md): N2_SUBREGGE's raw-form fits of the window data with JACKKNIFE-propagated errors (the stored fits quote the
residual-scaled lstsq error, which with 1 dof is not an error estimate), and paired (same-CRN) estimator differences.
  F1   : y - (A-1)kappa = c0 A[D] + s A[x D'/(2J+1)] + r A[1]            (A + C ring class + constant r)
  F1B  : + b A[nu^2 D/100]                                              (A + B + C)
Jackknife: the fit is redone on each leave-one-replica-out record set (a_h, bases, cv factors all recomputed, dc_pipe).
"""
import dc_guard, os, sys
import numpy as np
import n2_common as N

_COL = {}


def Acol(fx, nu, m, h):
    e = lambda x: N.env(h, x, 2)
    return (N.A_on_x(lambda x: e(x)*fx(x), nu, m)/e(nu/m)).real


def dprime(h, J, x, hh=1e-6):
    return (N.Dsh(h, J, x + hh) - N.Dsh(h, J, x - hh))/(2*hh)


def columns(J, m, nu, h):
    key = (J, round(m, 6), tuple(np.round(nu, 9)), h)
    if key not in _COL:
        x = nu/m
        AD = Acol(lambda z: N.Dsh(h, J, z), nu, m, h)
        ADp = Acol(lambda z: dprime(h, J, z)*z/(2*J + 1), nu, m, h)
        AB = Acol(lambda z: (z*m)**2*N.Dsh(h, J, z)/100.0, nu, m, h)
        A1 = Acol(lambda z: 0*z + 1.0, nu, m, h)
        dk = Acol(lambda z: N.kappa('C2', h, J, m, z), nu, m, h) - N.kappa('C2', h, J, m, x).real
        _COL[key] = dict(AD=AD, ADp=ADp, AB=AB, A1=A1, dk=dk)
    return _COL[key]


def select(o, nus_sel):
    nus = np.asarray(o['nus'])
    idx = [int(np.argmin(np.abs(nus - v))) for v in nus_sel]
    assert all(abs(nus[i] - v) < 1e-6 for i, v in zip(idx, nus_sel)), (nus, nus_sel)
    return np.array(idx)


def fit1(cols, y):
    A = np.array(cols).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c, y - A @ c


def fit(o, h, nus_sel, model='F1', kind='C2'):
    """returns dict(c, r, r_jk_err, rms, res, chi2 (diag jk errors), ...)."""
    idx = select(o, nus_sel); J, m = o['J'], o['m']; nu = np.asarray(o['nus'])[idx]
    C = columns(J, m, nu, h)
    cols = [C['AD'], C['ADp']] + ([C['AB']] if model == 'F1B' else []) + [C['A1']]
    y = o[(kind, h)]['y'][idx] - C['dk']
    yj = o[(kind, h)]['yj'][:, idx] - C['dk'][None, :]
    R = yj.shape[0]
    c, res = fit1(cols, y)
    cj = np.array([fit1(cols, yy)[0] for yy in yj])
    ej = np.sqrt((R - 1)/R*np.sum((cj - cj.mean(0))**2, axis=0))
    ey = np.sqrt((R - 1)/R*np.sum((yj - yj.mean(0))**2, axis=0))
    chi2 = float(np.sum((res/np.maximum(ey, 1.0))**2))
    dof = len(y) - len(cols)
    return dict(c=c, cj=cj, ec=ej, r=c[-1], er=ej[-1], rms=float(np.sqrt(np.mean(res**2))), res=res, ey=ey, y=y, yj=yj,
                chi2=chi2, dof=dof, nu=nu, model=model)


def paired(oa, ob, h, nus_sel, model='F1', kind='C2'):
    """difference of the fitted r (and of y per node) between two estimator configurations on the SAME common-seed replicas:
    jackknife of the difference (correlations included)."""
    fa, fb = fit(oa, h, nus_sel, model, kind), fit(ob, h, nus_sel, model, kind)
    R = fa['cj'].shape[0]
    dr = fa['cj'][:, -1] - fb['cj'][:, -1]
    edr = np.sqrt((R - 1)/R*np.sum((dr - dr.mean())**2))
    dy = fa['yj'] - fb['yj']
    edy = np.sqrt((R - 1)/R*np.sum((dy - dy.mean(0))**2, axis=0))
    return dict(dr=fa['r'] - fb['r'], edr=edr, dy=fa['y'] - fb['y'], edy=edy, fa=fa, fb=fb)

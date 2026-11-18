"""
l8_analysis (M6 round 4, doc 22): helicity-resolved graviton C^_2 action relative to the scalar k = 2 action, with
 (a) the LOG (cumulant) form of the Laplace series:  log E[p e^A] = log p0 + (p2 + p4)/p0 - (p2/p0)^2/2 + O(eps^6)  -- the nu-only
     factor exp(i a1/nu + a2/nu^2 + ...) sits INSIDE the Laplace terms (Im K1 m = a1/x, a1 ~ 21.7) and is resummed only in this form;
 (b) a REAL x-independent helicity basis T0(Delta) per state, from the order-0 (saddle) action (rank-1 decomposition over x) -- the basis
     defined by the saddle evaluation of the exact-Delta <TTO> structures, common to C^_2 and C^_3;
 (c) EXACT removal of every nu-only factor (incl. any helicity-dependent one) by the heavy reference state at the same nu:
     D_h(nu; Delta) = L_h(nu; Delta) - L_h(nu; Delta_ref),  L_h = log(M~_hh / Y_scal / W_h(x)).
     D_h = log(C_h(Delta)/C_h(Delta_ref)) [per-state constant, real] + I_h(x, m),  I_h -> 0 at fixed nu as m -> oo.
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l7_dc2_analysis as An
import l8_points as P

HEL = ('++', '+-')


def from_l7(J):
    """l7 fixed-x records (cumulative orders) -> l8 format, keyed (Delta, nu)."""
    fn = os.path.join(HERE, 'cache', 'l7_J%d.pkl' % J)
    if not os.path.exists(fn):
        return {}
    out = {}
    for (D, x), r in pickle.load(open(fn, 'rb')).items():
        M = r['M']; s = r['sS']
        out[(D, round(x*r['m'], 9))] = dict(M0=M[:, :, 0], M2=M[:, :, 1] - M[:, :, 0], M4=M[:, :, 2] - M[:, :, 1], logC=r['logC'],
                                             logK=r['logK'], slog=r['sl'], sS=np.array([s[0], s[1] - s[0], s[2] - s[1]]), m=r['m'])
    return out


def records(J):
    R = from_l7(J); R.update(P.load(J)); return R


def realvec(v):
    v = np.asarray(v, complex); ph = v[np.argmax(np.abs(v))]/abs(v[np.argmax(np.abs(v))])
    w = v/ph
    return w.real, np.abs(w.imag).max()/np.abs(w).max()


def basis(J, D, R, xs_min=0.15):
    """real helicity basis T0 (columns ++, +-) at Delta from the order-0 action at all cached nu with x >= xs_min."""
    m = P.mphys(D, J)
    if J < 4:
        return dict(T=np.eye(1), Ti=np.eye(1), sv3=0.0, imag=0.0)
    pts = sorted((k[1]/m, R[k]) for k in R if k[0] == D and k[1]/m >= xs_min)
    if len(pts) < 4:
        return None
    xs = np.array([p[0] for p in pts])
    sv, vs = An.rank1([p[1]['M0'] for p in pts])
    vv = [realvec(v) for v in vs]
    T = np.array([v[0] for v in vv]).T
    Ti = np.linalg.inv(T)
    # identify ++ by the shape d/W_++ flattest (relative to scalar), as l7
    flat = []
    for j in range(2):
        q = np.array([(Ti @ p[1]['M0'] @ Ti.T)[j, j]*np.exp(p[1]['logC'] - p[1]['logK'] - p[1]['slog'])/p[1]['sS'][0]
                      for p in pts])/An.W('++', J, xs)
        flat.append(np.std(np.abs(q/q[0])))
    if flat[0] > flat[1]:
        T = T[:, ::-1]; Ti = np.linalg.inv(T)
    return dict(T=T, Ti=Ti, sv3=sv[2], imag=max(v[1] for v in vv))


def logser(p0, p2, p4):
    r = p2/p0
    return np.log(p0) + r + p4/p0 - 0.5*r*r


def point(J, D, nu, R, B):
    """L_h (h = ++, +-), mixing kappa = M~_{+-}/sqrt(M~_++ M~_+-) (resummed diag), and raw K1 m, at (Delta, nu)."""
    r = R[(float(D), round(float(nu), 9))]
    x = nu/r['m']
    Mt = [B['Ti'] @ r[k] @ B['Ti'].T for k in ('M0', 'M2', 'M4')]
    ls = r['slog'] + logser(*r['sS'])
    out = dict(x=x, m=r['m'], ls=ls, lc=r['logC'] - r['logK'])
    for j, h in enumerate(HEL[:Mt[0].shape[0]]):
        lh = logser(Mt[0][j, j], Mt[1][j, j], Mt[2][j, j])
        out[h] = lh + r['logC'] - r['logK'] - ls - np.log(An.W(h, J, x) + 0j)
        out['K1m' + h] = r['m']*Mt[1][j, j]/Mt[0][j, j]
        out['l' + h] = lh
    if Mt[0].shape[0] == 1:
        return out
    out['r01'] = sum(0.5*(Mk[0, 1] + Mk[1, 0]) for Mk in Mt)/np.exp(out['l++'])      # M~_{+-}/M~_{++} (resummed ++)
    off = 0.5*sum(Mk[0, 1] + Mk[1, 0] for Mk in Mt)
    out['kappa'] = off/np.sqrt(np.exp(out['l++'] + out['l+-']))
    out['asym'] = abs(sum(Mk[0, 1] - Mk[1, 0] for Mk in Mt))/abs(np.sqrt(np.exp(out['l++'] + out['l+-'])))
    return out

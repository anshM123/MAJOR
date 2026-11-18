"""
n1_pw (doc 24): POINT-WISE backend for the crystal integrands of l1_crystal / l5_spin / l6_spin_action / l9_points, so that the exact
7-variable CLS-D integrand can be evaluated at many points z = z* + L xi (xi real) instead of as a truncated power series.

Mechanism: the integrand modules call l1_tps (T.mul, T.log, T.powr, ...) and l5_spin.bm.  Fresh copies of l1_crystal, l5_spin,
l6_spin_action and l9_points are exec'd with sys.modules['l1_tps'] -> this module, so every "TPS" becomes an array of point values
(last axis = point index).  Index 0 of every batch is the saddle (xi = 0): log and fractional powers are continued from the saddle,
log(a) = log(a_0) + Log(a/a_0) (principal Log of the ratio), which is the analytic continuation used by the TPS engine (principal
log at the saddle + series) as long as a/a_0 does not wind around 0 along the ray from the saddle.  `WIND` records the largest
|arg(a/a_0)| met (per batch), as a diagnostic of contour problems.

Exports: crystal (Cr), spin (S), A6, P9 point-wise module copies; set_batch(n).
"""
import os, sys, types
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

import l1_tps as _TT                      # the real TPS engine (DEG from L1_DEG; used for saddle/Hessian and shadow matrices)
import l1_crystal as _Cr_real             # noqa: F401  (cached in sys.modules before the swap)
import l5_spin as _S_real                 # noqa: F401
import l5_shadow_matrix as _SM_real       # noqa: F401
import l6_spin_action as _A6_real         # noqa: F401

NV = 7
DEG = 6
N = 1
WIND = [0.0]
WINDP = np.zeros(1)                       # per-point max |arg(a/a_0)| over all logs / powers of the current batch
MONS = _TT.MONS[:NV + 1]
_I = np.zeros(1, int); _J = np.zeros(1, int); _K = np.zeros(1, int)


def set_batch(n):
    global N
    global WINDP
    N = int(n)
    S.N = N
    WIND[0] = 0.0
    WINDP = np.zeros(N)


def const(c):
    return np.full(N, c, dtype=complex)


def mul(a, b):
    return a*b


def recip(a):
    return 1.0/a


def div(a, b):
    return a/b


MODE = ['plain']                          # 'plain' | 'record' | 'replay' | 'replay_record'  (path continuation, see evaluate_path)
REC = [[], [], 0]                         # [recorded logs of the previous level, logs being recorded, replay index]


def path_mode(mode):
    MODE[0] = mode
    if mode in ('record', 'replay_record'):
        if mode == 'replay_record':
            REC[0] = REC[1]
        REC[1] = []
    if mode == 'replay':
        REC[0] = REC[1]; REC[1] = []
    REC[2] = 0


def _rlog(a):
    a = np.asarray(a, complex)
    if MODE[0] in ('replay', 'replay_record'):
        ref = REC[0][REC[2]]; REC[2] += 1
        lr = np.log(a*np.exp(-ref))
        out = ref + lr
    else:
        lr = np.log(a/a[..., :1])
        out = np.log(a[..., :1]) + lr
    if MODE[0] in ('record', 'replay_record'):
        REC[1].append(out)
    global WINDP
    ai = np.abs(lr.imag)
    if ai.ndim > 1:
        ai = ai.reshape(-1, ai.shape[-1]).max(axis=0)
    if ai.shape[-1] == WINDP.shape[0]:
        WINDP = np.maximum(WINDP, ai)
    w = ai.max()
    if w > WIND[0]:
        WIND[0] = w
    return out


def log(a, branch0=None):
    out = _rlog(a)
    if branch0 is not None:
        out = out - out[..., :1] + branch0
    return out


def powr(a, p):
    return np.exp(p*_rlog(a))


def _load(name, swaps):
    saved = {k: sys.modules.get(k) for k in swaps}
    sys.modules.update(swaps)
    try:
        M = types.ModuleType(name + '_pw'); M.__dict__['__file__'] = os.path.join(HERE, name + '.py')
        src = open(os.path.join(HERE, name + '.py')).read()
        exec(compile(src, name + '_pw', 'exec'), M.__dict__)
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    return M


_self = sys.modules[__name__]
Cr = _load('l1_crystal', {'l1_tps': _self})
S = _load('l5_spin', {'l1_tps': _self})


def _bm(a, b):
    return np.asarray(a, complex)*np.asarray(b, complex)


S.bm = _bm
S.conserved = _S_real.conserved            # coefficient vectors: computed by the real (reduced-TPS) engine, cached there
A6 = _load('l6_spin_action', {'l1_tps': _self, 'l1_crystal': Cr, 'l5_spin': S})
_env = os.environ.get('L1_DEG')
os.environ['L1_DEG'] = '6'
P9 = _load('l9_points', {'l1_tps': _self, 'l1_crystal': Cr, 'l5_spin': S, 'l6_spin_action': A6})
if _env is None:
    os.environ.pop('L1_DEG', None)
else:
    os.environ['L1_DEG'] = _env


def points(zs, L, xi):
    """list of 7 point arrays z_i = zs_i + (L xi)_i; xi: (npts, 7) real; the saddle is prepended (index 0)."""
    X = np.vstack([np.zeros((1, NV)), xi])
    Z = zs[None, :] + X @ L.T
    set_batch(Z.shape[0])
    return [Z[:, i].copy() for i in range(NV)]

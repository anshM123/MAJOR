"""Vectorised ladders / per-harmonic coefficients (float64) for arrays of states.  Same conventions as core.py."""
import numpy as np


def rmat(b, J, D, K):
    b = np.asarray(b, float); J = np.asarray(J, float); D = np.asarray(D, float)
    C = J * (J + 1); T = (D - 1) * (D - 2)
    S = len(J); r = np.zeros((S, K + 1)); r[:, 0] = 1.0
    for k in range(K):
        s = b + k
        f = ((k + 4.0) / (k + 1.0)) ** 2 * (C - s * (s + 1)) / (T - s * (s + 1))
        f = np.where(k + 1 <= J - b, f, 0.0)
        r[:, k + 1] = r[:, k] * f
    return r


def dmat(b, J, D, K):
    b = np.asarray(b, float)[:, None]; J = np.asarray(J, float)[:, None]; D = np.asarray(D, float)[:, None]
    C = J * (J + 1); T = (D - 1) * (D - 2); j = np.arange(K + 1)[None, :]; s = b + j
    return 3 * (T - (s - 1) * s) / (j + 3) - 3 * (C - s * (s + 1)) / (j + 1)


def hmat(b, J, D, N, r=None):
    if r is None or r.shape[1] < N + 2:
        r = rmat(b, J, D, N + 1)
    d = dmat(b, J, D, N + 1)
    S = r.shape[0]; h = np.zeros((S, N + 1))
    for n in range(2, N + 1):
        h[:, n] = r[:, n - 2] * (1 + d[:, n - 2] / 3) / 2 + (n + 1) / 4.0 * (r[:, n - 3] if n >= 3 else 0) - (n - 1) / 4.0 * r[:, n - 1]
    return h


def Wv(b, D):
    b = np.asarray(b, float); D = np.asarray(D, float)
    return ((D - 1) * (D - 2) - b * (b - 1)) / 3.0 + 1.5


def states(g, Jlin, Jmax, nJlog, ts, light_h=None):
    """heavy: twist = g + t for t in ts; light: twist in light_h (2<h<g)."""
    Js = np.array(sorted(set(list(range(0, Jlin + 1)) + [int(x) for x in np.geomspace(Jlin + 1, Jmax, nJlog)])))
    B, JJ, DD, L = [], [], [], []
    for J in Js:
        for b in (0, 4):
            if (b == 4 and J < 4) or (b == 0 and J % 2):
                continue
            for t in ts:
                B.append(b); JJ.append(J); DD.append(J + g + t); L.append(0)
            if light_h is not None:
                for h in light_h:
                    B.append(b); JJ.append(J); DD.append(J + h); L.append(1)
    return np.array(B), np.array(JJ, float), np.array(DD, float), np.array(L)

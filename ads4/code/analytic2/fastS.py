"""fastS: vectorised dual-Hahn recurrence (complex128 / longdouble-free) for many exponents at once; DV only."""
import numpy as np, math

def S_vec(n, a, s, alphas):
    al = np.asarray(alphas, dtype=complex)
    N = n + a
    if s == 1:
        Sm = np.ones_like(al); S = np.ones_like(al)
        # track log scale to avoid overflow
        lg = np.zeros(al.shape)
        for m in range(n):
            eps = 2*(2*m + 1)/((N - m)*(N + m + 1))
            Sm, S = S, Sm + eps*al*S
            if m % 64 == 0:
                sc = np.maximum(np.abs(S), np.abs(Sm)); sc[sc == 0] = 1
                S = S/sc; Sm = Sm/sc; lg += np.log(sc)
        return S, lg
    Sm = np.zeros_like(al); S = np.ones_like(al); lg = np.zeros(al.shape)
    for m in range(n):
        Bm = (N + m + s)*(m + s)*(N - m)/((2*m + s)*(2*m + s + 1))
        Dm = -m*(N + m + s)*(N - m)/((2*m + s - 1)*(2*m + s)) if m > 0 else 0.0
        Sm, S = S, (al*S + (Bm + Dm)*S - Dm*Sm)/Bm
        if m % 64 == 0:
            sc = np.maximum(np.abs(S), np.abs(Sm)); sc[sc == 0] = 1
            S = S/sc; Sm = Sm/sc; lg += np.log(sc)
    return S, lg

def d_of(n, a, s):
    N = n + a
    return a*(2*n + a + s)/(N*(N + s))

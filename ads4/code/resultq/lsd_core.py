"""
lsd_core (LIGHTSIDE_DIST.md): the swapped-form functional Psi_k on EIN's product correlator T' = sigma x T_free, with the sigma part
in its Mellin representation (lsd_mellin), evaluated three ways (Omega_line, Omega_OPE, and the residue at infinity).

  H[c](t1, t3) = w1 w3 t1^5 t3^5 f_k G_free[c] (X14 X23)^-2 / 2        (the 1/2: EIN's I_Delta = MB/2, checked in lsd_run 'val')
  K[c](d12, d13) = sum_{t1,t3} H[c] exp(-d12 LA + d13 LB),  LA = A(t1) + A(t3),  LB = B(t1) + B(t3) - log 8y     (separable)
  E[p x c] = (1/(2 pi i)^2) int int dd12 dd13 Gamma(d12)^2 Gamma(d13)^2 Gamma(d14)^2 M_p(d) K[c](d)       (Omega_line)
  Omega_OPE[p x c] = (1/(2 pi i)) int dd13 Gamma(d13)^2 sum_{left poles in d12} Res[Gamma(d12)^2 Gamma(d14)^2 M_p K[c]]
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import loggamma, digamma
import og_free as F
import lsd_mellin as LM

CLS = ('1234', '1243', '1324', 's', 't', 'u')
PIECES = ('c0', 'c2', 'c3', 'c4', 'xs', 'xt', 'xu')
# Mellin contours (Re d12, Re d13) per piece: every s-, t-, u-channel pole on its own side (validated pointwise, lsd_run 'val')
CONT = {'c0': (2/3, 2/3), 'c2': (2/3, 2/3), 'c3': (2/3, 2/3), 'c4': (2/3, 2/3), 'xs': (1.3, 0.35), 'xt': (0.35, 0.35),
        'xu': (0.35, 1.3), 'xu4': (0.35, 1.3), 'r4': (2/3, 2/3), 'r4u': (0.35, 1.3)}


class Setup:
    """grid, H matrices per class, column functions A, B for one (k, phi1, phi3, p2, p4)."""

    def __init__(self, k, phi1, phi3, p2, p4, s0=1e-6, m=12, q=3.0, weight=5, classes=CLS):
        self.k = k
        t, w, th = LM.contour_C(s0, m, q)
        self.t, self.w = t, w
        y = 2*np.sin((phi1 - phi3)/2)**2
        self.y = y
        z1 = F.zvec(phi1); z3 = F.zvec(phi3)
        nt = len(t)
        self.H = {c: np.empty((nt, nt), complex) for c in classes}
        ch = max(1, 40000//nt)
        for i0 in range(0, nt, ch):
            T1 = t[i0:i0 + ch, None]; T3 = t[None, :]
            Lk = F.links(T1, T3, z1, z3, p2, p4, y13=y)
            G = F.correlator(Lk, classes)
            T1b, T3b = np.broadcast_arrays(T1, T3)
            h = (T1b*T3b)**weight*F.fk(k, T1b, T3b, z1, z3, y13=y)*(w[i0:i0 + ch, None]*w[None, :])
            h = h*((2*T1b + 1)*(2*T3b + 1))**-2.0/2
            for c in classes:
                self.H[c][i0:i0 + ch] = h*G[c]
        self.A = LM.Acol(t); self.B = LM.Bcol(t)
        self.l8y = np.log(8*y)

    def K(self, d12, d13, deriv=False):
        """K[c] for arrays d12 (shape (n,)) at scalar d13; with deriv also dK/dd12.  Returns dict c -> (n,) arrays."""
        d12 = np.atleast_1d(np.asarray(d12, complex))
        u = np.exp(-d12[:, None]*self.A[None, :] + d13*self.B[None, :])            # (n, N)
        pref = np.exp(-d13*self.l8y)
        out = {}
        for c, H in self.H.items():
            W = u @ H                                                                 # (n, N): sum over t1 index
            K = np.einsum('ij,ij->i', W, u)*pref
            if deriv:
                WA = (u*self.A[None, :]) @ H
                dK = -(np.einsum('ij,ij->i', WA, u) + np.einsum('ij,ij->i', W, u*self.A[None, :]))*pref
                out[c] = (K, dK)
            else:
                out[c] = K
        return out


def lmeas(d12, d13):
    return 2*(loggamma(d12) + loggamma(d13) + loggamma(2 - d12 - d13))


def omega_line(S, pieces=PIECES, Y12=30.0, h12=0.05, Y13=20.0, h13=0.05, cont=CONT, cutoff=None):
    """Omega_line[p x c] by trapezoid on both Mellin lines.  cutoff: optional |Im d12| truncation (symmetric) for convergence studies.
    Returns {(p, c): value}; also the partial values at a list of symmetric Im d12 cutoffs if cutoff is a list."""
    out = {}
    cuts = cutoff if cutoff is not None else [Y12]
    for (c1, c2) in set(cont[p] for p in pieces):
        ps = [p for p in pieces if cont[p] == (c1, c2)]
        y12 = np.arange(-Y12, Y12 + h12/2, h12); y13 = np.arange(-Y13, Y13 + h13/2, h13)
        acc = {(p, c, Yc): 0j for p in ps for c in S.H for Yc in cuts}
        d12 = c1 + 1j*y12
        for yy in y13:
            d13 = c2 + 1j*yy
            Kd = S.K(d12, d13)
            lm = lmeas(d12, d13)
            for p in ps:
                fac = np.exp(lm)*LM.M_amp(p, d12, d13)
                for c in S.H:
                    v = fac*Kd[c]
                    for Yc in cuts:
                        acc[(p, c, Yc)] += np.sum(v[np.abs(y12) <= Yc + 1e-12])
        for key, v in acc.items():
            out[key] = v*h12*h13/(4*np.pi**2)
    if cutoff is None:
        return {(p, c): v for (p, c, _), v in out.items()}
    return out


def residues_d12(S, p, d13, N):
    """left-pole residues in d12 of Gamma(d12)^2 Gamma(d14)^2 M_p K[c] at fixed d13: double poles d12 = -n (n = 0..N), plus the
    s-exchange pole d12 = 1 for xs.  Returns {c: array (N+1,)} (and the xs single-trace residue under key (c, 'st'))."""
    n = np.arange(N + 1)
    d12 = -n.astype(complex)
    d14 = 2 - d12 - d13
    Kd = S.K(d12, d13, deriv=True)
    lg = 2*loggamma(d14) - 2*loggamma(n + 1.0)
    g = np.exp(lg)
    M = LM.M_amp(p, d12, d13); dM = LM.dM_amp(p, d12, d13)
    psi14 = digamma(d14); psin = digamma(n + 1.0)
    out = {}
    for c, (K, dK) in Kd.items():
        G = M*K
        dG = -2*psi14*M*K + dM*K + M*dK
        out[c] = g*(dG + 2*psin*G)
    if p == 'xs':
        K1 = S.K(np.array([1.0 + 0j]), d13)
        Q = (d13 - (1 - d13))**2/4 + 1/12
        for c in S.H:
            out[(c, 'st')] = np.exp(2*loggamma(1 - d13))*Q*K1[c][0]
    return out


def residues_right(S, p, d13, M_):
    """right-pole residues in d12 (t-channel: Gamma(d14)^2 double poles d14 = -m, i.e. d12 = 2 - d13 + m; plus the t-exchange pole
    d14 = 1 for xt).  Returned as Res (counter-clockwise convention), {c: array (M_+1,)} and {(c, 'st'): value} for xt."""
    m = np.arange(M_ + 1)
    d12 = 2 - d13 + m.astype(complex)
    Kd = S.K(d12, d13, deriv=True)
    g = np.exp(2*loggamma(d12) - 2*loggamma(m + 1.0))
    Mv = LM.M_amp(p, d12, d13); dMv = LM.dM_amp(p, d12, d13)
    psi12 = digamma(d12); psim = digamma(m + 1.0)
    out = {}
    for c, (K, dK) in Kd.items():
        G = Mv*K
        dG = 2*psi12*Mv*K + dMv*K + Mv*dK
        out[c] = g*(dG - 2*psim*G)
    if p == 'xt':
        # M = Q_t/(d14 - 1), d14 - 1 = 1 - d12 - d13: simple pole at d12 = 1 - d13 with residue -Q_t Gamma(d12)^2 Gamma(1)^2 K
        d12s = 1 - d13
        K1 = S.K(np.array([d12s]), d13)
        Q = (d13 - d12s)**2/4 + 1/12
        for c in S.H:
            out[(c, 'st')] = -np.exp(2*loggamma(d12s))*Q*K1[c][0]
    return out


def omega_closed(S, pieces, N=60, Y13=15.0, h13=0.1, cont=CONT):
    """(1/(2 pi)) int dy13 Gamma(d13)^2 [sum_left Res + sum_right Res] = the d12 contour at infinity (all poles).
    Returns {(p, c): (left, right)}."""
    y13 = np.arange(-Y13, Y13 + h13/2, h13)
    out = {}
    for p in pieces:
        c1, c2 = cont[p]
        accL = {c: 0j for c in S.H}; accR = {c: 0j for c in S.H}
        for yy in y13:
            d13 = c2 + 1j*yy
            w13 = np.exp(2*loggamma(d13))
            RL = residues_d12(S, p, d13, N); RR = residues_right(S, p, d13, N)
            for c in S.H:
                accL[c] += w13*(np.sum(RL[c]) + RL.get((c, 'st'), 0.0))
                accR[c] += w13*(np.sum(RR[c]) + RR.get((c, 'st'), 0.0))
        for c in S.H:
            out[(p, c)] = (accL[c]*h13/(2*np.pi), accR[c]*h13/(2*np.pi))
    return out


def omega_ope(S, pieces=PIECES, N=200, Y13=20.0, h13=0.05, cont=CONT, Nlist=None):
    """Omega_OPE[p x c] = (1/(2 pi)) int dy13 Gamma(d13)^2 sum_res; partial sums at the n-cutoffs Nlist.  Returns {(p, c, Nc): value}."""
    Nlist = Nlist or [N]
    y13 = np.arange(-Y13, Y13 + h13/2, h13)
    out = {}
    for p in pieces:
        c1, c2 = cont[p]
        acc = {(c, Nc): 0j for c in S.H for Nc in Nlist}
        for yy in y13:
            d13 = c2 + 1j*yy
            R = residues_d12(S, p, d13, N)
            w13 = np.exp(2*loggamma(d13))
            for c in S.H:
                cs = np.cumsum(R[c])
                st = R.get((c, 'st'), 0.0)
                for Nc in Nlist:
                    acc[(c, Nc)] += w13*(cs[Nc] + st)
        for (c, Nc), v in acc.items():
            out[(p, c, Nc)] = v*h13/(2*np.pi)
    return out


def closed_fast(S, pieces, N=60, Y13=12.0, h13=0.1, cont=CONT, Ncuts=(20, 40)):
    """as omega_closed, but K, dK are computed once per (d13 line, node) and shared by all pieces on that line.
    Returns {(p, c): {'L': left sum, 'R': right sum, 'L<Nc>': partial left sums}} (all with the 1/(2 pi) int dy13 Gamma(d13)^2)."""
    y13 = np.arange(-Y13, Y13 + h13/2, h13)
    n = np.arange(N + 1); nf = n.astype(float)
    out = {(p, c): {'L': 0j, 'R': 0j, **{f'L{Nc}': 0j for Nc in Ncuts}} for p in pieces for c in S.H}
    lines = sorted(set(cont[p][1] for p in pieces))
    for c2 in lines:
        ps = [p for p in pieces if cont[p][1] == c2]
        for yy in y13:
            d13 = c2 + 1j*yy
            w13 = np.exp(2*loggamma(d13))*h13/(2*np.pi)
            dL = -nf.astype(complex); dR = 2 - d13 + nf
            pts = np.concatenate([dL, dR, [1.0 + 0j, 1 - d13]])
            Kd = S.K(pts, d13, deriv=True)
            gL = np.exp(2*loggamma(2 - d13 + nf) - 2*loggamma(nf + 1)); psiL14 = digamma(2 - d13 + nf); psin = digamma(nf + 1)
            gR = np.exp(2*loggamma(dR) - 2*loggamma(nf + 1)); psiR12 = digamma(dR)
            for p in ps:
                ML = LM.M_amp(p, dL, d13); dML = LM.dM_amp(p, dL, d13)
                MR = LM.M_amp(p, dR, d13); dMR = LM.dM_amp(p, dR, d13)
                for c, (K, dK) in Kd.items():
                    KL, dKL = K[:N + 1], dK[:N + 1]; KR, dKR = K[N + 1:2*N + 2], dK[N + 1:2*N + 2]
                    rL = gL*((-2*psiL14*ML + dML)*KL + ML*dKL + 2*psin*ML*KL)
                    rR = gR*((2*psiR12*MR + dMR)*KR + MR*dKR - 2*psin*MR*KR)
                    st_L = st_R = 0j
                    if p == 'xs':
                        st_L = np.exp(2*loggamma(1 - d13))*((d13 - (1 - d13))**2/4 + 1/12)*K[2*N + 2]
                    if p == 'xt':
                        d12s = 1 - d13
                        st_R = -np.exp(2*loggamma(d12s))*((d13 - d12s)**2/4 + 1/12)*K[2*N + 3]
                    o = out[(p, c)]
                    cs = np.cumsum(rL)
                    o['L'] += w13*(cs[-1] + st_L); o['R'] += w13*(np.sum(rR) + st_R)
                    for Nc in Ncuts:
                        o[f'L{Nc}'] += w13*(cs[Nc] + st_L)
    return out

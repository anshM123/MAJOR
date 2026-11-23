r"""
rh_t2_heft (R4ARC_HEFT.md sec. 3): H_EFT(a) for heavy exchanges on spinning externals, AT FIXED CELESTIAL ANGLES, split by the Wick
class of the spin-2 factor (product family T' = sigma x T of LIGHTSIDE_DIST / F_NORM sec. 3.4).

Heavy sigma-exchange (s + t + u), dimension Delta_h = tau_h + J_h, spin J_h in {0, 2, 4} (even), Delta_sigma = 2, d = 3:
    W_h = sum_m w_m [ rho_m(d13 - d14)/(d12 - dm) + rho_m(d13 - d12)/(d14 - dm) + rho_m(d12 - d14)/(d13 - dm) ],
    dm = 2 - tau_h/2 - m,   w_m rho_m(y) = K_m Qpoly^m(y)  (exact d = 3 Mack polynomials and Mack weights, exactskel/exactmom.BlockMellin
    with Df = 2; y = d13 - d14 is the Mack variable y = T' - 2Df + tau/2 + m on the s-pole; parity (-1)^J makes the three channels
    1<->3 symmetric).  For J_h = 0 this is fn_spin's exchange.
EFT (1/m) expansion: 1/(x - dm) = -sum_j x^j/dm^{j+1}; the exact remainder after order R is
    X^(R+1) = sum_m w_m rho_m(y) x^{R+1}/(dm^{R+1} (x - dm))     (per channel, x = d12, d14, d13),
so  TP_R - DT = -V[X^(R+1)]  with no cancellation (V = symmetrised term-wise value; DT = V[X^(0)] = double-trace part of W_h).
Single-trace part S_h: residues at d12 = dm (m = 0..MS).
Classes: s, t, 1324 (dispersive: only finitely many d12 = -n residues are non-zero, n <= 5; summed to n <= NDISP) and 1243
(non-dispersive: n <= N160 = 160, as in F_NORM B3).
Usage: python3 rh_t2_heft.py run <k> <ab|ba> <tag>   ;   python3 rh_t2_heft.py report <k> <tag>
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
import numpy as np
from scipy.special import loggamma, digamma
import rh_common as R
import exactmom as EM

CLASSES = ('s', 't', '1324', '1243')
DISP = ('s', 't', '1324')
NDISP, NFULL, MS = 12, 160, 160
TAUS = (41.0, 61.0, 81.0)
JHS = (0, 2, 4)
MODE = os.environ.get('RH_MODE', 'deg')     # 'deg': EFT organised by total Mellin degree (= powers of 1/m); 'j': by propagator order
J0S = tuple([0] + list(range(1, 41)))          # remainders X^(j0), j0 = 0 (full exchange) ... 40


class Exch:
    def __init__(self, tau, J):
        B = EM.BlockMellin(tau + J, J, 2.0, N0=400)
        self.tau, self.J = tau, J
        lr = B.logK0 + 2*B.g1
        self.Lr = lr.max()
        self.w = np.exp(lr - self.Lr)*B.wn                      # all nodes (direct + Euler-Maclaurin tail)
        self.ns = B.ns
        self.dm = 2 - tau/2 - B.ns
        self.G = B.G if not B.cheb else None                     # [nn, J+1] coefficients in y
        self.lr_direct = (lr - self.Lr)[:MS + 1]                 # integer m = 0..MS (direct nodes)

    def rho(self, y, nodes=slice(None)):
        """rho_m(y) and d rho/dy, y: (P,) -> (P, nn)."""
        G = self.G[nodes]
        y = y[:, None]
        Q = np.zeros(y.shape[:1] + (G.shape[0],), complex) + G[None, :, -1]
        dQ = np.zeros_like(Q)
        for i in range(self.J - 1, -1, -1):
            dQ = dQ*y + Q
            Q = Q*y + G[None, :, i]
        return Q, dQ

    def remainders(self, d12, d13, j0s):
        """M = X^(j0) and dM/dd12 at fixed d13 for every j0 in j0s; d12: (P,).  Returns arrays (len(j0s), P).
        X^(j0) = sum over channels and i of y^i x^e A_{i,e}(x),  A_{i,e}(x) = sum_m w_m g_i(m) dm^-e/(x - dm),
        e = max(j0 - i, 0) (MODE 'deg': remainder after total Mellin degree j0 - 1) or e = j0 (MODE 'j': after propagator order j0 - 1)."""
        d14 = 2 - d12 - d13
        dm = self.dm
        J = self.J
        Emax = max(j0s)
        epow = dm[:, None]**(-np.arange(Emax + 1)[None, :])                      # (nn, E+1)
        Wm = (self.w[:, None, None]*self.G[:, :, None]*epow[:, None, :]).reshape(len(dm), -1)   # (nn, (J+1)(E+1))
        out = np.zeros((len(j0s), len(d12)), complex); dout = np.zeros_like(out)
        chans = [(d12, 1.0, 2*d13 + d12 - 2, 1.0), (d14, -1.0, d13 - d12, -1.0), (np.full(d12.shape, d13, complex), 0.0, 2*d12 + d13 - 2, 2.0)]
        for x, dx, y, dy in chans:
            Inv = 1.0/(x[:, None] - dm[None, :])
            A = (Inv @ Wm).reshape(len(x), J + 1, Emax + 1)
            B = ((Inv*Inv) @ Wm).reshape(len(x), J + 1, Emax + 1)
            for a, j0 in enumerate(j0s):
                for i in range(J + 1):
                    e = max(j0 - i, 0) if MODE == 'deg' else j0
                    yi = y**i; dyi = i*y**(i - 1) if i > 0 else 0.0
                    xe = x**e; dxe = e*x**(e - 1) if e > 0 else 0.0
                    out[a] += yi*xe*A[:, i, e]
                    dout[a] += dyi*dy*xe*A[:, i, e] + yi*(dxe*A[:, i, e] - xe*B[:, i, e])*dx
        return out, dout


def run(k, order, tag):
    a, b, p2, p4 = R.P.config(tag) if tag in ('P1', 'P2', 'P3') else SPECIAL[tag]
    aa, bb = (a, b) if order == 'ab' else (b, a)
    t0 = time.time()
    S = R.setup(k, aa, bb, p2, p4, classes=CLASSES)
    print(f'setup {len(S.t)} [{time.time()-t0:.0f}s]', flush=True)
    EX = {(tau, J): Exch(tau, J) for tau in TAUS for J in JHS}
    d13s, h = R.d13_line(Y=8.0, h=0.2)
    nf = np.arange(NFULL + 1).astype(float); d12 = -nf.astype(complex)
    acc = {(e, j0, c): np.zeros(NFULL + 1, complex) for e in EX for j0 in J0S for c in CLASSES}
    sh = {(e, c): np.zeros(MS + 1, complex) for e in EX for c in CLASSES}
    for i, d13 in enumerate(d13s):
        w13 = np.exp(2*loggamma(d13))*h/(2*np.pi)
        Kd = S.K(d12, d13, deriv=True)
        d14 = 2 - d12 - d13
        g = np.exp(2*loggamma(d14) - 2*loggamma(nf + 1)); psi14 = digamma(d14); psin = digamma(nf + 1)
        for e, E in EX.items():
            Mv, dMv = E.remainders(d12, d13, J0S)
            for c, (K, dK) in Kd.items():
                r = w13*g[None, :]*((-2*psi14[None, :]*Mv + dMv)*K[None, :] + Mv*dK[None, :] + 2*psin[None, :]*Mv*K[None, :])
                for aj, j0 in enumerate(J0S):
                    acc[(e, j0, c)] += r[aj]
        # single-trace residues at d12 = dm, m = 0..MS (direct integer nodes)
        for tau in TAUS:
            dmv = (2 - tau/2 - np.arange(MS + 1)).astype(complex)
            Ks = S.K(dmv, d13)
            lg = 2*loggamma(dmv) + 2*loggamma(2 - dmv - d13)
            for J in JHS:
                E = EX[(tau, J)]
                Q, _ = E.rho(2*d13 + dmv - 2, nodes=slice(0, MS + 1))
                rq = np.diagonal(Q) if Q.shape[0] == Q.shape[1] else None
                rq = np.array([Q[m, m] for m in range(MS + 1)])
                wt = np.exp(lg + E.lr_direct)*rq
                for c, K in Ks.items():
                    sh[((tau, J), c)] += w13*wt*K
        if i % 10 == 0:
            print(f'  node {i}/{len(d13s)} [{time.time()-t0:.0f}s]', flush=True)
    pickle.dump({'acc': acc, 'sh': sh, 'k': k, 'order': order, 'tag': tag}, open(f'{R.CACHE}/t2{"" if MODE == "j" else "d"}_{tag}_k{k}_{order}.pkl', 'wb'))
    print(f'done [{time.time()-t0:.0f}s]', flush=True)


SPECIAL = {}


def report(k, tag):
    D = {o: pickle.load(open(f'{R.CACHE}/t2_{tag}_k{k}_{o}.pkl', 'rb')) for o in ('ab', 'ba')}
    sgn = (-1.0)**k
    L = []
    pr = lambda s: (print(s, flush=True), L.append(s))
    pr(f'H_EFT(a) at fixed angles, tag {tag}, k = {k}; V = E(a,b) + (-1)^k E(b,a); dispersive classes summed to n <= {NDISP}, 1243 to n <= {NFULL}')
    for e in sorted(set(x[0] for x in D['ab']['acc'])):
        tau, J = e
        for c in CLASSES:
            nmax = NDISP if c in DISP else NFULL
            V = {j0: (np.sum(D['ab']['acc'][(e, j0, c)][:nmax + 1]) + sgn*np.sum(D['ba']['acc'][(e, j0, c)][:nmax + 1])) for j0 in J0S}
            Sh = np.sum(D['ab']['sh'][(e, c)]) + sgn*np.sum(D['ba']['sh'][(e, c)])
            DT = V[0]
            # first EFT order with non-negligible value
            eft = {j: V[j] - V[j + 1] for j in J0S[:-1]}
            scale = max(abs(DT), abs(Sh), 1e-300)
            jl = next((j for j in sorted(eft) if abs(eft[j]) > 1e-6*scale), None)
            errs = ' '.join(f'{abs(V[j0])/scale:.1e}' for j0 in (1, 5, 9, 10, 12, 15, 20, 25, 30, 40))
            lead = (DT/eft[jl] - 1).real if jl is not None else np.nan
            pr(f'  tau {tau:4.0f} J {J}  {c:5s}  S_h {Sh.real:+.4e}  DT {DT.real:+.4e}  (S_h+DT)/|S_h| {abs(Sh + DT)/max(abs(Sh),1e-300):.1e}'
               f'  | |TP_R - DT|/scale at R+1 = 1,5,9,10,12,15,20,25,30,40: {errs}  | lead order j = {jl}, DT/EFT_j - 1 = {lead:+.3e}')
    open(f'{R.CACHE}/t2_{tag}_k{k}.txt', 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]), sys.argv[3], sys.argv[4])
    else:
        report(int(sys.argv[2]), sys.argv[3])

"""
p5 (PARITY_ODD_ROUTEA.md sect. 4): the normalisation N of the parity-odd exchange, from the crystal action at the SAME state.

PARITY_ODD_HEAVY.md sect. 4 found W[S,S]/W[X^E_opp, X^E_opp] = N = 3.77 / 3.30 (J = 19 / 35) for CHL-normalised structures (X^O = i S), where
CHL's construction (equal helicity moduli of X^O and X^E) suggests |N| = 1.  Here both blocks are built as <T1T2O>_a beta^{-1} <O~T3T4>_b with
the shadow coefficient computed (p2), so no seed normalisation enters:
  odd J : X^E_opp lives in the 1<->2 ANTISYMMETRIC parity-even conserved space (2-dim, with X^E_same); X^O_opp = (1/kappa_O) x p1's odd
          structure (p1b).  Crystal structure set [anti-even (2), odd (1)] (P3_ANTI=1).
  even J: X^E_same in the symmetric even space (l5_spin.conserved), X^O_same = (1/kappa_O) p1.
  N_kind(x) := [M_OO/kappa_O^2] / [lambda^T M_EE lambda],  lambda = coefficients of X^E_kind (fit to chl_xspace at matched points);
  both entries at the same (nu, Delta, J), so every common factor (fast exponent, nu-only factors, scalar reference) cancels.
Expectation if the helicity heavy formula holds with CHL's normalisation: N -> +1 at leading order (sign +: 34-side coupling = conjugate).
Usage: P3_ANTI=1 python3 p5_N.py fit Delta J      (CHL fits)
       P3_ANTI=1 python3 p5_N.py run J Delta x1 ... (crystal points via p3, both kinds)
       P3_ANTI=1 python3 p5_N.py report
"""
import os, sys, functools, itertools
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
TT = os.path.join(HERE, '..', 'fourpt', 'tttt'); sys.path.insert(0, TT)
import numpy as np
import l4_tto_even as L4
import p1_tto_odd as P1


@functools.lru_cache(maxsize=None)
def conserved_even_anti(Delta, J):
    """(mons, sol) of the conserved parity-even structures with f(1<->2) = -f (numerical null space, l4 code), real basis."""
    import l5_spin as S
    L = S._l4r()
    rng = np.random.default_rng(3)
    mons = L.monomials(J); K = len(mons); npts = 3*K
    pts = [L.rand_pt(rng) for _ in range(npts)]
    Fv = np.array([L.structures(*p, Delta, J) for p in pts])
    Cv = np.array([L.conservation_row(p, Delta, J, 1) for p in pts] + [L.conservation_row(p, Delta, J, 2) for p in pts])
    Sv = np.array([L.structures(p[2], p[3], p[0], p[1], p[4], p[5], Delta, J) for p in pts]) + Fv
    rn = lambda A_: A_/np.maximum(np.abs(A_).max(axis=1, keepdims=True), 1e-300)
    sc = np.maximum(np.abs(Fv).max(axis=1, keepdims=True), 1e-300)
    Fv, Cv, Sv = Fv/sc, rn(Cv), Sv/sc
    U, sF, Vh = np.linalg.svd(Fv); rF = int(np.sum(sF > 1e-9*sF[0]))
    B = Vh[:rF].conj().T
    A = np.vstack([Cv @ B, Sv @ B])
    sA = np.linalg.svd(A, compute_uv=False); n0 = int(np.sum(sA < 1e-8*sA[0]))
    _, _, VhA = np.linalg.svd(A)
    sol = B @ VhA[len(VhA) - n0:].conj().T
    U2, sv_, _ = np.linalg.svd(np.hstack([sol.real, sol.imag]), full_matrices=False)
    nreal = int(np.sum(sv_ > 1e-8*sv_[0])); assert nreal == n0, (Delta, J, sv_)
    sol = U2[:, :nreal]; sol = sol/np.abs(sol).max(axis=0, keepdims=True)
    return tuple(mons), sol


def even_emb(Dv, J, X1, Z1, X2, Z2, X5, Z5, mons):
    """l4 parity-even basis structures at embedding vectors (numbers; exact powers)."""
    L = L4
    t1 = t2 = 5.0; t5 = Dv + J
    V1, V2, V5 = P1.eV(L, X1, Z1, X2, X5), P1.eV(L, X2, Z2, X5, X1), P1.eV(L, X5, Z5, X1, X2)
    H12, H15, H25 = P1.eH(L, X1, Z1, X2, Z2), P1.eH(L, X1, Z1, X5, Z5), P1.eH(L, X2, Z2, X5, Z5)
    X12, X15, X25 = P1.eXX(L, X1, X2), P1.eXX(L, X1, X5), P1.eXX(L, X2, X5)
    den = L.P(X12, -(t1 + t2 - t5)/2)*L.P(X15, -(t1 + t5 - t2)/2)*L.P(X25, -(t2 + t5 - t1)/2)
    return [den*V1**a1*V2**a2*V5**a5*H12**b12*H15**b15*H25**b25 for (a1, a2, a5, b12, b15, b25) in mons]


def fit_chl(Dv, J, kind, ncfg=12):
    """lambda (real-basis coefficients) with  sum_a lambda_a E_a = X^E_kind (CHL, chl_xspace), and residuals; also kappa_O (p1/X^O_kind)."""
    import sympy as sp, mpmath as mp
    import odd3 as O, chl_xspace as CX
    from ft_helicity import AX, BX, WX, XX, AB, AW, BW
    import p1b_validate_odd3 as V
    import l5_spin as S
    mp.mp.dps = 30
    if J % 2 == 0:
        mons, sol = S.conserved(Dv, J)[:2]
    else:
        mons, sol = conserved_even_anti(Dv, J)
    T = CX.xspace_structure(2, sp.nsimplify(Dv), sp.Integer(J), kind)
    syms = [AX, BX, WX, XX, AB, AW, BW]
    f = sp.lambdify(syms, T, 'mpmath')
    cfgs = O.configs(ncfg, seed=202)
    cfgs = [([c/mp.sqrt(sum(cc*cc for cc in x)) for c in x], a, b, z) for (x, a, b, z) in cfgs]   # |x| = 1: no Delta-power overflow
    A, b = [], []
    for (x, a, bb, z) in cfgs:
        I = O._inv(x, a, bb, z); b.append(complex(f(*[I[s] for s in syms])))
        y = [-c for c in x]; yl, al, bl, zl = V.lc(y), V.lc(a), V.lc(bb), V.lc(z)
        X1 = P1.emb_X(L4, yl); Z1 = P1.emb_Z(L4, yl, al)
        X2 = [1.0, 0, 0, 0, 0]; Z2 = [0.0, 0.0] + bl; X5 = [0.0, 1.0, 0, 0, 0]; Z5 = [0.0, 0.0] + zl
        A.append(np.array(even_emb(Dv, J, X1, Z1, X2, Z2, X5, Z5, list(mons))) @ sol)
    A, b = np.array(A), np.array(b)
    lam, *_ = np.linalg.lstsq(A, b, rcond=None)
    res = np.linalg.norm(A @ lam - b)/np.linalg.norm(b)
    okind = 'same' if J % 2 == 0 else 'opp'
    mo, so, _ = P1.conserved_odd(Dv, J)
    ro = np.array(V.chl(Dv, J, okind, cfgs[:6]))
    oo = np.array([V.ours(Dv, J, c, mo, so) for c in cfgs[:6]])
    kap = oo/ro
    return dict(lam=lam, res=res, imag=np.abs(lam.imag).max()/np.abs(lam).max(), kappa=kap[0].real,
                kspread=np.abs(kap/kap[0] - 1).max())


def ratio(J, D, x, kind='C2', lam=None, kap=None, order=0):
    import p3_odd_action as P3
    m = P3.mphys(D, J)
    r = P3.get(kind, J, D, [x*m], verbose=False)[0]
    M = r['M0'] + (r['M2'] if order >= 1 else 0) + (r['M4'] if order >= 2 else 0)
    nE = M.shape[0] - 1
    num = M[nE, nE]/kap**2
    den = lam @ M[:nE, :nE] @ lam
    return num/den, M


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'fit':
        Dv, J = float(sys.argv[2]), int(sys.argv[3])
        kind = 'same' if J % 2 == 0 else 'opp'
        r = fit_chl(Dv, J, kind)
        print(f"Delta={Dv} J={J} X^E_{kind}: lambda={r['lam']} res={r['res']:.1e} imag={r['imag']:.1e}; kappa_O={r['kappa']:.10e} "
              f"(spread {r['kspread']:.1e})")


def report_state(J, D, xs=(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), pr=print, kinds=('C2',)):
    kind_chl = 'same' if J % 2 == 0 else 'opp'
    f = fit_chl(D, J, kind_chl)
    lam, kap = f['lam'].real, f['kappa']
    out = {}
    for kd in kinds:
        for order in (0, 2):
            Ns = []
            for x in xs:
                try:
                    n, M = ratio(J, D, x, kd, lam, kap, order)
                except Exception:
                    continue
                Ns.append((x, n))
            if Ns:
                v = np.array([n for _, n in Ns])
                pr(f"  {kd} J={J} Delta={D:.0f} order<={order}: N(x) = " + " ".join(f"{x:.1f}:{n.real:+.5f}{n.imag:+.1e}i" for x, n in Ns) +
                   f"   (CHL fit res {f['res']:.0e}, kappa spread {f['kspread']:.0e})")
                out[(kd, order)] = v
    return out

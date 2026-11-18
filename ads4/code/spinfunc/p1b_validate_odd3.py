"""
p1b: validation of the covariant parity-odd conserved <T T O_{Delta,J}> (p1_tto_odd) against odd3.py (fourpt/tttt; CHL's parity-odd helicity
structures X^O_kind built in position space from their momentum-space definition, independently of any conservation solve).

Matched points: Euclidean configurations (x, a, b, zeta) of odd3.configs (x real, a, b, zeta complex null).  odd3 compares X^O_kind(x) (CHL
frame O1(0), O2(x), O3(inf)) with embedding words at P1 = (1, y^2, y), Z1 = (0, 2y.a, a), P2 = (1,0,0), Z2 = (0,0,b), P5 = (0,1,0),
Z5 = (0,0,zeta), y = -x.  Our structure is evaluated at the SAME embedding vectors, with the Euclidean vectors mapped to our lightcone
components by u = e1 + i e2, v = -(e1 - i e2), y = e3 (so that our Lorentzian dot product equals the Euclidean one); then our eps (det in
(+,-,u,v,y) components) = 2i x the Euclidean det.  Claim checked: p1(odd, conserved, symmetric) = kappa(Delta, J) X^O_kind, configuration-
independent kappa, with kind = 'same' for even J and 'opp' for odd J (CHL: the identical-T odd structure feeds a^{++} at even J, a^{+-}
at odd J), and kappa REAL: our Lorentzian-real structure is a real multiple of CHL's X^O = i S (S Euclidean-real).
Usage: python3 p1b_validate_odd3.py [Delta:J ...]  -> cache/p1b_validate_odd3.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
TT = os.path.join(HERE, '..', 'fourpt', 'tttt'); sys.path.insert(0, TT)
import numpy as np
import mpmath as mp
import sympy as sp
import l4_tto_even as L4
import p1_tto_odd as P


def lc(e):
    e = [complex(c) for c in e]
    return [e[0] + 1j*e[1], -(e[0] - 1j*e[1]), e[2]]


def ours(Dv, J, cfg, mons, sol):
    x, a, b, z = cfg
    y = [-c for c in x]
    yl, al, bl, zl = lc(y), lc(a), lc(b), lc(z)
    X1 = P.emb_X(L4, yl); Z1 = P.emb_Z(L4, yl, al)
    X2 = [1.0, 0.0, 0.0, 0.0, 0.0]; Z2 = [0.0, 0.0] + bl
    X5 = [0.0, 1.0, 0.0, 0.0, 0.0]; Z5 = [0.0, 0.0] + zl
    v = np.array(P.structures_odd_emb(L4, X1, Z1, X2, Z2, X5, Z5, Dv, J, list(mons)))
    return complex(v @ sol[:, 0])


def chl(Dv, J, kind, cfgs):
    import odd3 as O
    from ft_helicity import AX, BX, WX, XX, AB, AW, BW
    from odd3 import EX1, EX2, EP, xspace_odd, _inv, _eps3
    mp.mp.dps = 30
    T = xspace_odd(2, Dv, J, kind, 30)
    syms = [AX, BX, WX, XX, AB, AW, BW]
    f = sp.lambdify(syms + [EX1, EX2, EP], T, 'mpmath')
    out = []
    for (x, a, b, z) in cfgs:
        I = _inv(x, a, b, z)
        out.append(complex(f(*[I[s] for s in syms], _eps3(a, x, z), _eps3(b, x, z), _eps3(a, b, z))))
    return out


def run(Dv, J, ncfg=8, pr=print):
    import odd3 as O
    mons, sol, sv = P.conserved_odd(Dv, J)
    if sol is None:
        pr(f"  Delta={Dv} J={J}: no conserved symmetric odd structure"); return None
    cfgs = O.configs(ncfg, seed=101)
    res = {}
    for kind in ('same', 'opp'):
        c = chl(Dv, J, kind, cfgs)
        o = [ours(Dv, J, cfg, mons, sol) for cfg in cfgs]
        r = np.array(o)/np.array(c)
        res[kind] = r
        pr(f"  Delta={Dv} J={J} vs X^O_{kind:4s}: ratio ours/CHL = {r[0]:.10e}; spread max|r/r0 - 1| = {np.abs(r/r[0] - 1).max():.1e}; "
           f"|Im r0/Re r0| = {abs(r[0].imag/r[0].real) if r[0].real != 0 else np.inf:.1e}")
    return res


if __name__ == '__main__':
    args = sys.argv[1:] or ['23.7:2', '23.7:5', '43.697092:5', '101.1:4', '101.1:7', '43.697092:19']
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    pr("# p1b: covariant odd <TTO> (p1) vs odd3 CHL X^O_same / X^O_opp at 8 matched Euclidean points (expect: constant real ratio with")
    pr("#      'same' for even J, 'opp' for odd J; the other kind must NOT be proportional)")
    for a in args:
        Dv, J = a.split(':'); run(float(Dv), int(J), pr=pr)
    open(os.path.join(HERE, 'cache', 'p1b_validate_odd3.txt'), 'w').write("\n".join(lines) + "\n")

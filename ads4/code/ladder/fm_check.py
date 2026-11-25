"""LADDER_ALLJ (A): pipeline check of the ANGULAR_ROUTE §7 first-moment operator for J = 0,2,4,5,6,7,8 (J = 1,3 have no
conserved <TTO> structure), both parities, symbolic Delta, source helicity +2, detector n = (1,0).
Tested statement:   M^orig_m = (i/6) [ (W - 18) H ]_m ,   W = T + B1^2 + B2^2 + R^2 + (4-2Delta) B1 - 2 R B2   (spinor generators
on e_m = a^{J+m} b^{J-m}; (W H)_m' = sum_m coef(e_m' in W e_m) H_m),  M^orig = u-moment about the T insertion point (x3 = 0),
H = zeroth moment, both unnormalised (fixed rest-frame polarisations; pipeline component m = coeff of t^{J-m}).
Also tested: midpoint moment M^mid = M^orig + (i/2)[(Delta+3) - c B1] H, i.e. the translation piece (i/2) z.dF with
z.dF = (Delta+1 scaling) + (2, detector weight) - (boost of the fixed bra polarisation); Lorentz covariance predicts c = 1.
Usage: python3 -B fm_check.py [J ...]  -> cache/fm_check_J*.txt"""
import sys
from common import *

def W_apply(Hvec, J, Dv):
    """Hvec: dict m -> value; returns dict m' -> (W H)_{m'} with Delta = Dv (sympy)."""
    T = (Dv - 1)*(Dv - 2)
    B1 = lambda f: sp.Rational(1, 2)*(b_*sp.diff(f, a_) + a_*sp.diff(f, b_))
    B2 = lambda f: sp.I/2*(-b_*sp.diff(f, a_) + a_*sp.diff(f, b_))
    R = lambda f: -sp.I/2*(a_*sp.diff(f, a_) - b_*sp.diff(f, b_))
    f = sum(Hvec.get(m, 0)*a_**(J + m)*b_**(J - m) for m in range(-J, J + 1))
    g = sp.expand(T*f + B1(B1(f)) + B2(B2(f)) + R(R(f)) + (4 - 2*Dv)*B1(f) - 2*R(B2(f)))
    P = sp.Poly(g, a_, b_)
    return {m: P.coeff_monomial(a_**(J + m)*b_**(J - m)) for m in range(-J, J + 1)}, \
        sp.Poly(sp.expand(B1(f)), a_, b_)

def run(J):
    out = []
    pr = lambda x: (print(x, flush=True), out.append(str(x)))
    for parity in ('even', 'odd'):
        monos, basis = kl.conserved(J, parity)
        if not basis:
            continue
        kin = pl.Kin('plus')
        e0 = [pl.scri_mono(m, kin) for m in monos]
        emid = [moment_expr(e, parity, 'mid') for e in e0]
        eorg = [moment_expr(e, parity, 'orig') for e in e0]
        K = pl.Kbounds(e0 + emid + eorg, parity)
        a0 = [pl.amplitude(e, kin, parity, K=K)[0] for e in e0]
        am = [pl.amplitude(e, kin, parity, K=K)[0] for e in emid]
        ao = [pl.amplitude(e, kin, parity, K=K)[0] for e in eorg]
        for bi, coeffs in enumerate(basis):
            get = lambda polys: {m: kl.poly_to_sym(v).subs(S, (Dl + J)/2)
                                 for m, v in pl.tcoeffs(kl.combine(polys, coeffs), J).items()}
            H, Mm, Mo = get(a0), get(am), get(ao)
            Morig = {m: sp.expand(Mm.get(m, 0) + Mo.get(m, 0)) for m in range(-J - 1, J + 2)}
            WH, B1H = W_apply(H, J, Dl)
            pred = {m: sp.I/6*(WH[m] - 18*H.get(m, 0)) for m in range(-J, J + 1)}
            bad = [m for m in range(-J, J + 1) if sp.simplify(Morig[m] - pred[m]) != 0]
            extra = [m for m in (-J - 1, J + 1) if Morig[m] != 0]
            pr(f'J={J} {parity} structure {bi}: nonzero H_m at m = {sorted(m for m in H if H[m] != 0)}')
            pr(f'   origin-centred  M = (i/6)(W-18)H : {"EXACT for all m" if not bad and not extra else "FAILS at m=" + str(bad + extra)}')
            # midpoint: M^mid - M^orig = -(translation piece) ; test = (i/2)[(Delta+3) H - 2 B1 H] ?
            B1d = {m: B1H.coeff_monomial(a_**(J + m)*b_**(J - m)) for m in range(-J, J + 1)}
            for cB in (1, -1, 2, -2, 0):
                ok = all(sp.simplify(Mm.get(m, 0) - Morig[m] - sp.I/2*((Dl + 3)*H.get(m, 0) - cB*B1d[m])) == 0
                         for m in range(-J, J + 1))
                if ok:
                    pr(f'   midpoint: M^mid = M^orig + (i/2)[(Delta+3) - {cB} B1] H  EXACT')
                    break
            else:
                pr('   midpoint: no fit of the form M^orig + (i/2)[(Delta+3) - c B1]H, c in {1,-1,2,-2,0}')
    return out

if __name__ == '__main__':
    Js = [int(x) for x in sys.argv[1:]] or [0, 2, 4, 5, 6, 7, 8]
    for J in Js:
        write(f'fm_check_J{J}.txt', run(J))

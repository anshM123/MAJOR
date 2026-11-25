"""LADDER_ALLJ (B): pipeline zeroth-moment amplitudes <O,m|E(n)|T,h> (rest frame, n = (1,0)), h = +2 and h = -2, all conserved
<TTO> structures, at rational Delta, compared with the closed forms of ANGULAR_ROUTE §8 / KERNELS_D3 (R0),(R4):
   h=+2:  F_{-k} = F_0 C(J,k)(4)_k/(Delta-1-k)_k ,  F_{4+k} = F_4 C(J-4,k)(4)_k/(Delta-5-k)_k ,  F_1=F_2=F_3=0,
          and the §8 recurrence  m(m-4)(T+C-2m^2)F_m = (m-3)(m-4)(D-m-2)(J+m+1)F_{m+1} + m(m-1)(D+m-2)(J-m+1)F_{m-1}  at every m;
   h=-2:  F^{(-2)}_m = sigma F^{(+2)}_{-m}  (sigma = +1 parity-even, -1 parity-odd structures), i.e. the mirrored ladders.
Usage: python3 -B ladder_data.py Jmin Jmax  -> cache/ladder_J*.txt (one line summary per J, parity, structure, Delta)"""
import sys, time
from fractions import Fraction as Fr
from common import *

def closed(J, D):
    F = {m: Fr(0) for m in range(-J, J + 1)}
    t0 = {-k: Fr(sp.binomial(J, k)*sp.rf(4, k))/Fr(str(sp.rf(sp.Rational(D.numerator, D.denominator) - 1 - k, k))) for k in range(J + 1)}
    t4 = {4 + k: Fr(sp.binomial(J - 4, k)*sp.rf(4, k))/Fr(str(sp.rf(sp.Rational(D.numerator, D.denominator) - 5 - k, k))) for k in range(max(J - 3, 0))}
    return t0, t4

def comps(J, parity, coeffs, polys):
    tot = kl.combine(polys, coeffs)
    out = {}
    for m, v in pl.tcoeffs(tot, J).items():
        d = v.to_dict()
        re = sum((Fr(int(c.p), int(c.q)) for (et, ed, es, eI), c in d.items() if eI == 0), Fr(0))
        im = sum((Fr(int(c.p), int(c.q)) for (et, ed, es, eI), c in d.items() if eI == 1), Fr(0))
        out[m] = (re, im)
    return out

def run(J, Dvals):
    lines = []
    for parity in ('even', 'odd'):
        for D in Dvals:
            sval = (D + J)/2
            t1 = time.time()
            monos, basis = kl.conserved(J, parity, sval=str(sval))
            if not basis:
                continue
            amps = {}
            for src in ('plus', 'minus'):
                kin = pl.Kin(src)
                ex = [pl.scri_mono(m, kin) for m in monos]
                K = pl.Kbounds(ex, parity)
                amps[src] = [pl.amplitude(e, kin, parity, sig=sval, K=K)[0] for e in ex]
            t0, t4 = closed(J, D)
            T = (D - 1)*(D - 2); C = J*(J + 1)
            for bi, coeffs in enumerate(basis):
                coeffs = [sp.Rational(c) for c in coeffs]
                Fp = comps(J, parity, coeffs, amps['plus']); Fm = comps(J, parity, coeffs, amps['minus'])
                z = (Fr(0), Fr(0))
                g = lambda m: Fp.get(m, z)
                # sector content and closed forms (complex ratios)
                sec0 = g(0) != z; sec4 = g(4) != z
                ok = all(g(m) == z for m in (1, 2, 3)) and all(g(m) == z for m in (-J - 1, J + 1))
                cplx = lambda u: complex(float(u[0]), float(u[1]))
                def ratio_ok(m, ref, r):        # F_m == r * F_ref exactly (componentwise)
                    return g(m) == (g(ref)[0]*r, g(ref)[1]*r)
                if sec0:
                    ok &= all(ratio_ok(m, 0, t0[m]) for m in t0)
                else:
                    ok &= all(g(m) == z for m in t0)
                if J >= 4 and sec4:
                    ok &= all(ratio_ok(m, 4, t4[m]) for m in t4)
                else:
                    ok &= all(g(m) == z for m in range(4, J + 1))
                # recurrence at every m in [-J, J]  (real and imaginary parts separately)
                rec = True
                for m in range(-J, J + 1):
                    for q in (0, 1):
                        lhs = m*(m - 4)*(T + C - 2*m*m)*g(m)[q]
                        rhs = (m - 3)*(m - 4)*(D - m - 2)*(J + m + 1)*g(m + 1)[q] + m*(m - 1)*(D + m - 2)*(J - m + 1)*g(m - 1)[q]
                        rec &= (lhs == rhs)
                # h = -2 mirror
                sig = 1 if parity == 'even' else -1
                mir = all(Fm.get(m, z) == (sig*g(-m)[0], sig*g(-m)[1]) for m in range(-J - 1, J + 2))
                mir_other = all(Fm.get(m, z) == (-sig*g(-m)[0], -sig*g(-m)[1]) for m in range(-J - 1, J + 2))
                lines.append(f'J={J:2d} {parity:4s} str{bi} Delta={str(D):>10s} sectors(0,4)=({int(sec0)},{int(sec4)}) '
                             f'closed-forms+zeros: {ok}  recurrence(all m): {rec}  h=-2 mirror F^-_m = {sig:+d}F^+_-m: {mir}'
                             + ('' if mir else f' (opposite sign: {mir_other})') + f'  [{time.time()-t1:.0f}s]')
                print(lines[-1], flush=True)
    return lines

if __name__ == '__main__':
    J0, J1 = int(sys.argv[1]), int(sys.argv[2])
    for J in range(J0, J1 + 1):
        Dvals = [Fr(J + 1) + Fr(7, 3), Fr(J + 1) + Fr(1000) + Fr(1, 7)]
        write(f'ladder_J{J}.txt', run(J, Dvals))

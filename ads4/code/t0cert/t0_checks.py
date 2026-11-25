"""
t0cert/t0_checks.py (T0_CERT.md §5): float / point-Arb sanity checks [DV] of the identities used by the certificate (not part of the proof).
  1. g1(x) := -(x^2/8) r(x) - (x^3/24) r'(x)  ==  (1/24) int chi'' cos(p x) dp   (Lemma 2(iii)), Arb at points;
  2. Lam^2 curv_k(bhat, Lam) vs a1 M1(beta) + F_k(bhat)/24 at Delta_g = 300, 1000, 10^4 (largespin.ls_core Mehler-Dirichlet quadrature);
  3. the T0 ratio lambda = H/N_s on a sample of states (ee1000_common.hyb pieces via largespin.ls_core.hyb_fast) at 300 / 1000 / 10^4,
     including the near-sharp flat states (+-, J = 5 and ++, J = 4 at m = Lam).
Output: cache/t0_checks.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import t0_common as T
import ls_core as LS

lines = []; pr = T.pr_to(lines)


def M1f(beta):
    f = lambda ph: (lambda sig: 1/sig**2 - 1/((2*mp.asinh(sig/2))**2*mp.sqrt(1 + sig**2/4)))(2*mp.sinh(beta/2)*mp.sin(ph))
    return float(2/mp.pi*mp.quad(f, [0, mp.pi/4, mp.pi/2]))


def main():
    pr('1. g1 identity (Arb, points):')
    for k in (2, 3):
        for x in (Fr(1, 2), Fr(3), Fr(11), Fr(40)):
            d = T.r_derivs_point(k, x, 1)
            g1 = -(T.A(x)**2/8)*d[0] - (T.A(x)**3/24)*d[1]
            # (1/24) int chi'' cos(px) dp via the cosine series of chi''
            M = int(2*float(x)) + 40
            mu = [sum(ci/(i + j + 1) for i, ci in enumerate(T.CHI2[k])) for j in range(2*M + 2)]
            s = sum(((-1)**m)*T.A(mu[2*m])*T.A(x)**(2*m)/mp.factorial(2*m) for m in range(M)) if False else None
            from flint import arb, ctx
            ctx.prec = 400
            acc = arb(0); X = T.A(x)
            fac = arb(1)
            for m in range(M):
                if m > 0:
                    fac = fac*(2*m - 1)*(2*m)
                t = T.A(mu[2*m])*X**(2*m)/fac
                acc += t if m % 2 == 0 else -t
            ctx.prec = 128
            pr(f'   k={k} x={float(x):5.1f}: g1 = {g1.str(10)}   (1/24)int chi\'\' cos = {(acc/24).str(10)}')
    pr('2. Lam^2 curv_k vs a1 M1(beta) + F_k/24:')
    for Dg in (300, 1000, 10000):
        L = LS.lam(Dg)
        for b in (3, 12, 40, 150, 0.9*2*L):
            beta = 2*np.arcsinh(b/(2*L))
            row = []
            for k in (2, 3):
                cv = L**2*(LS.OmLR(k, [beta], Dg)[0] - LS.OmFlat(k, [b], Dg)[0])
                F, _ = T.F_point(k, Fr(b).limit_denominator(10**6))
                pred = float(T.A1[k])*M1f(beta) + float((F/24).mid())
                row.append(f'k={k}: {cv:.6f} / {pred:.6f}')
            pr(f'   Dg={Dg} bhat={b:8.2f}: ' + ';  '.join(row))
    pr('3. T0 ratio H/N_s (float hybrid, ls_core.hyb_fast):')
    mn = (9, None)
    for Dg in (300, 1000, 10000):
        L = LS.lam(Dg)
        for h, Js in (('++', [0, 2, 4, 6, 10, 20, 50, 120, 400]), ('+-', [4, 5, 6, 7, 9, 15, 33, 101, 400])):
            for J in Js:
                for mfl in (LS.mgap(Dg, J), 1.0, 1.02, 1.3, 2.0, 5.0, 20.0):
                    if mfl < LS.mgap(Dg, J) and mfl != 1.0:
                        continue
                    d = LS.hyb_fast(Dg, J, mfl, h)
                    lam = d['base']/d['Ns']
                    if lam < mn[0]:
                        mn = (lam, (Dg, h, J, mfl, d['bhat']))
    pr(f'   min over the sample: {mn[0]:.10f} at (Dg, h, J, m/Lam, bhat) = {mn[1]}  (certified c = {T.C} = {float(T.C)})')
    open(os.path.join(T.CACHE, 't0_checks.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()

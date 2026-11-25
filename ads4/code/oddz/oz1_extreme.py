"""PARITY_ODD_DZ §2: the parity-odd extreme detector amplitudes for SYMBOLIC J and their zeros.
(a) re-runs ANGULAR_PROOFS ap5 (read-only; output captured here) for S^odd (M_{-J}) and S^oddJ (M_{+J}): the symbolic-J
    Chu-Vandermonde amplitude divided by the closed form is exactly 1;
(b) rewrites 1/Gamma(3/2-h) = -cos(pi h) Gamma(h-1/2)/pi (reflection) and shows: for h > 1/2 the amplitude vanishes
    exactly at h = 3/2 + n (n = 0, 1, ...), simply, with an explicit nonzero derivative; all other factors are finite and
    nonzero there;
(c) ladder factors (all m) at Delta0 = J+3+2n: finite and positive for every k, symbolic J and n;
(d) control: the even structures (S^++, S^+-) carry 1/Gamma(1-h) = sin(pi h) Gamma(h)/pi: zeros at h = 1 + n.
Output: cache/oz1_out.txt"""
import subprocess
from common import *
out = []
pr = lambda x: (print(x, flush=True), out.append(str(x)))

# (a)
r = subprocess.run(['python3', '-B', 'ap5_kernels.py', 'odd:minusJ', 'oddJ:plusJ', 'odd:plusJ', 'oddJ:minusJ'],
                   cwd=ANGP, capture_output=True, text=True, env={**os.environ, 'KERNELS_CACHE': os.path.join(ANGP, 'cache')})
pr('(a) ANGULAR_PROOFS ap5 re-run (read-only):')
for line in r.stdout.strip().splitlines():
    pr('    ' + line)

# (b)
n = sp.Symbol('n', integer=True, nonnegative=True)
Jv = sp.Symbol('J', integer=True, positive=True)
h, s_, D = sp.symbols('h s Delta', positive=True)
refl = sp.simplify(sp.gamma(sp.Rational(3, 2) - h)*(-sp.cos(sp.pi*h)*sp.gamma(h - sp.Rational(1, 2))/sp.pi))
pr(f'(b) Gamma(3/2-h) * [-cos(pi h) Gamma(h-1/2)/pi] = {sp.simplify(sp.gammasimp(refl))}  (reflection formula)')
Modd = lambda J_, s, Dl, hh: 8*(-1)**(J_/2)*sp.exp(sp.I*sp.pi*s)*sp.pi**3*2**(-J_)*sp.factorial(J_ + 3) \
    / (sp.gamma(s + sp.Rational(1, 2))*sp.gamma(Dl + 3))*(-sp.cos(sp.pi*hh)*sp.gamma(hh - sp.Rational(1, 2))/sp.pi)
ModdJ = lambda J_, s, Dl, hh: -8*sp.I*(-1)**((J_ - 5)/2)*sp.exp(sp.I*sp.pi*s)*sp.pi**3*2**(-J_)*sp.factorial(J_ - 1) \
    / (sp.gamma(s + sp.Rational(1, 2))*sp.gamma(Dl - 1))*(-sp.cos(sp.pi*hh)*sp.gamma(hh - sp.Rational(1, 2))/sp.pi)
pr('    M_{-J}[S^odd]  = -8 (-1)^{J/2} e^{i pi s} pi^2 2^{-J} (J+3)!  cos(pi h) Gamma(h-1/2) / (Gamma(s+1/2) Gamma(Delta+3))')
pr('    M_{+J}[S^oddJ] =  8i (-1)^{(J-5)/2} e^{i pi s} pi^2 2^{-J} (J-1)! cos(pi h) Gamma(h-1/2) / (Gamma(s+1/2) Gamma(Delta-1))')
pr('    at h = 3/2+n: s = J+3/2+n, Delta = J+3+2n; Gamma(h-1/2) = n!, Gamma(s+1/2) = (J+1+n)!, Gamma(Delta+3) = (J+5+2n)!,')
pr('    Gamma(Delta-1) = (J+1+2n)!: finite and nonzero; cos(pi h) has a simple zero (d/dh cos(pi h) = -pi sin(pi h) = (-1)^n pi).')
# derivative with respect to Delta at Delta0 (h = (Delta-J)/2): d/dDelta cos(pi h) = -(pi/2) sin(pi h)
for nm, M, Jt in (('S^odd', Modd, 'even'), ('S^oddJ', ModdJ, 'odd')):
    for Jn in ((4, 6, 10, 30) if Jt == 'even' else (5, 7, 11, 31)):
        for nn in (0, 1, 2, 7):
            h0 = sp.Rational(3, 2) + nn; s0 = Jn + h0; D0 = Jn + 2*h0
            eps = sp.Symbol('eps')
            expr = M(Jn, s0 + eps/2, D0 + eps, h0 + eps/2)
            val0 = sp.simplify(expr.subs(eps, 0))
            der = sp.simplify(sp.diff(expr, eps).subs(eps, 0))
            pred = sp.simplify((M(Jn, s0, D0, h0)/(-sp.cos(sp.pi*h0)) if False else 0))
            if nn == 0 or nn == 7:
                pr(f'    {nm:6s} J={Jn:2d} n={nn}: M(Delta0) = {val0},  dM/dDelta(Delta0) = {sp.nsimplify(der) if der.is_number else der}'
                   .replace('\n', ' ')[:220])
            assert val0 == 0 and der != 0
pr('    (asserted for J in {4,6,10,30}/{5,7,11,31}, n in {0,1,2,7}: value 0, derivative nonzero)')
# symbolic-J statement of the derivative
hh0 = sp.Rational(3, 2) + n
pr('    general: dM_{-J}[S^odd]/dDelta |_{Delta0} = -8(-1)^{J/2} e^{i pi s0} pi^2 2^{-J}(J+3)! * (-(pi/2) sin(pi h0)) * n!/((J+1+n)!(J+5+2n)!),'
   f'  sin(pi h0) = {sp.simplify(sp.sin(sp.pi*hh0))}  -> nonzero for every J, n')

# (c) ladder factors at Delta0
k = sp.Symbol('k', integer=True, nonnegative=True)
pr('(c) ladder (KERNELS_D3 (R0),(R4); ANGULAR_ROUTE §8):  F_{-k}/F_{-J} = C(J,k)(4)_k (Delta-1-J)_J / ((Delta-1-k)_k (4)_J),')
pr('    F_{4+k}/F_{J} = C(J-4,k)(4)_k (Delta-1-J)_{J-4} / ((Delta-5-k)_k (4)_{J-4}).  At Delta0 = J+3+2n:')
pr('    (Delta0-1-k)_k = (J+2+2n-k)...(J+1+2n) >= (2+2n)_k > 0 for 0 <= k <= J; (Delta0-1-J)_J = (2+2n)_J > 0;')
pr('    (Delta0-5-k)_k = (J-2+2n-k)...(J-3+2n) > 0 for 0 <= k <= J-4 (J >= 5).  Every ratio is finite and > 0.')
for Jn in (4, 5, 6, 9, 12, 25, 40):
    for nn in (0, 3):
        D0 = Jn + 3 + 2*nn
        r0 = [sp.binomial(Jn, kk)*sp.rf(4, kk)*sp.rf(D0 - 1 - Jn, Jn)/(sp.rf(D0 - 1 - kk, kk)*sp.rf(4, Jn)) for kk in range(Jn + 1)]
        r4 = [sp.binomial(Jn - 4, kk)*sp.rf(4, kk)*sp.rf(D0 - 1 - Jn, Jn - 4)/(sp.rf(D0 - 5 - kk, kk)*sp.rf(4, Jn - 4))
              for kk in range(Jn - 3)] if Jn >= 5 else [1]
        assert all(x > 0 for x in r0 + r4)
pr('    asserted exactly for J in {4,5,6,9,12,25,40}, n in {0,3}.')

# (d) control
pr('(d) control: M_{-J}[S^++] and M_{+J}[S^+-] carry 1/Gamma(1-h) = sin(pi h) Gamma(h)/pi: simple zeros at h = 1+n')
pr('    (Delta = J+2+2n, the parity-even [TT]); at odd twist they are nonzero (|sin| = 1). Conversely the odd structures')
pr('    are nonzero at even twist (|cos| = 1).  The two zero sets are disjoint and interlaced.')
write('oz1_out.txt', out)

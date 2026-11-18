"""
Exact Mack polynomial P_{Delta,J}(s,t) (same polynomial and output layout as mack.mack_P_fast / mack_flint.mack_P_flint) with an
integer-matrix organisation of the Costa--Goncalves--Penedones sum, for J ~ 100-150:

  P = sum_r pre_r (l1)_r (l2)_r  sum_{k13+k14+k23+k24 = K} c * (s1)_{k13} (s1)_{k24} (sig)_{k14} (sig)_{k23},  s1 = -s/2, sig = (s+t)/2,
  c = (-1)^{k13+k24}/(k13! k14! k23! k24!) (a)_{L-k13-k14} (a)_{L-k23-k24} (b)_{L-k13-k23} (b)_{L-k14-k24},  K = J-2r, L = J-r.
With Delta = N/64: (a)_i = RA_i/128^i, (b)_i = RB_i/128^i with integers RA, RB, and the four Pochhammer indices always sum to 2J, so
  K! 128^{2J} 2^K * inner(s, sig) = sum_{pairs (k13,k24)} [2^{K-k13-k24} A_{k13} A_{k24}](s) [sum_{k14} (+-) multinom RA RA RB RB pb_{k14} pb_{k23}](sig)
(A_k(s) = prod_{q<k}(2q - s), pb_k = (sig)_k): an integer matrix product Amat^T Bmat over the pairs (fmpz_mat, multimodular),
then sig -> (s+t)/2 by Horner in fmpq_mpoly.  Self-test: python3 mack_mat.py (identity with mack_flint.mack_P_flint).
"""
from fractions import Fraction as Fr
from math import factorial, comb
import functools
import flint

CTX = flint.fmpq_mpoly_ctx.get(('s', 't'), 'lex')
S, T = CTX.gens()


def _rf(x, k):
    r = Fr(1)
    for q in range(k):
        r *= (x + q)
    return r


def _poch(lin, k):
    out = CTX.constant(1)
    for q in range(k):
        out = out*(lin + q)
    return out


@functools.lru_cache(None)
def mack_P_mat(Delta, J):
    D = Fr(Delta)
    assert (D*64).denominator == 1, "Delta must be a multiple of 1/64"
    h = Fr(3, 2)
    hp = D; hm = 2*h - D
    a = 1 - (hp + J)/2; b = 1 - (hm + J)/2
    Na, Nb = int(a*128), int(b*128)
    assert Fr(Na, 128) == a and Fr(Nb, 128) == b
    RA = [1]; RB = [1]
    for q in range(J + 1):
        RA.append(RA[-1]*(Na + 128*q)); RB.append(RB[-1]*(Nb + 128*q))
    den = _rf(hp - 1, J)*_rf(hm - 1, J)
    xs = flint.fmpz_poly([0, 1])
    Apol = [flint.fmpz_poly([1])]
    for q in range(J):
        Apol.append(Apol[-1]*(2*q - xs))
    pb = [flint.fmpz_poly([1])]
    for q in range(J):
        pb.append(pb[-1]*(xs + q))
    sig = (S + T)/2
    P = CTX.constant(0)
    fact = [factorial(i) for i in range(J + 1)]
    for r in range(J//2 + 1):
        aJr = Fr((-1)**r*factorial(J), 2**(2*r)*factorial(r)*factorial(J - 2*r))/_rf(h + J - 1 - r, r)
        pre = aJr*Fr(2)**(J + 2*r)*factorial(J - 2*r)/den
        l1 = flint.fmpq((hp - J).numerator, 2*(hp - J).denominator) - T/2
        l2 = flint.fmpq((hm - J).numerator, 2*(hm - J).denominator) - T/2
        pr = _poch(l1, r)*_poch(l2, r)
        K = J - 2*r; L = J - r
        pairs = [(k13, k24) for k13 in range(K + 1) for k24 in range(K + 1 - k13)]
        Arows, Brows = [], []
        pbb = {}
        for k13, k24 in pairs:
            rest = K - k13 - k24
            combp = flint.fmpz_poly([])
            for k14 in range(rest + 1):
                k23 = rest - k14
                c = fact[K]//(fact[k13]*fact[k14]*fact[k23]*fact[k24])
                c *= RA[L - k13 - k14]*RA[L - k23 - k24]*RB[L - k13 - k23]*RB[L - k14 - k24]
                if c == 0:
                    continue
                key = (min(k14, k23), max(k14, k23))
                if key not in pbb:
                    pbb[key] = pb[k14]*pb[k23]
                combp += pbb[key]*c
            if (k13 + k24) % 2:
                combp = -combp
            if combp.is_zero():
                continue
            ap = Apol[k13]*Apol[k24]*(2**rest)
            Arows.append([int(v) for v in ap.coeffs()] + [0]*(K + 1 - (ap.degree() + 1)))
            Brows.append([int(v) for v in combp.coeffs()] + [0]*(K + 1 - (combp.degree() + 1)))
        if not Arows:
            continue
        Am = flint.fmpz_mat(Arows); Bm = flint.fmpz_mat(Brows)
        C = Am.transpose()*Bm                       # C[i, j]: coefficient of s^i sig^j  (times K! 128^{2J} 2^K)
        scale = Fr(1, fact[K]*128**(2*J)*2**K)
        inner = CTX.constant(0)
        for j in range(K, -1, -1):
            row = CTX.from_dict({(i, 0): int(C[i, j]) for i in range(K + 1) if C[i, j] != 0}) if any(C[i, j] != 0 for i in range(K + 1)) else CTX.constant(0)
            inner = inner*sig + row
        P = P + pr*inner*flint.fmpq((pre*scale).numerator, (pre*scale).denominator)
    out = {}
    for mon, c in zip(P.monoms(), P.coeffs()):
        out[(int(mon[0]), int(mon[1]))] = Fr(int(c.p), int(c.q))
    return out


if __name__ == '__main__':
    import sys, time, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import mack_flint
    for J in (0, 1, 2, 3, 6, 9, 14, 21):
        for D in (900.25, 301.125, 3000.015625):
            A = mack_flint.mack_P_flint(D, J); B = mack_P_mat(D, J)
            ok = (set(A) == set(B)) and all(A[k] == B[k] for k in A)
            print('J', J, 'Delta', D, 'identical to mack_P_flint:', ok, flush=True)
    for J in [int(a) for a in sys.argv[1:]]:
        t = time.time(); mack_P_mat(900.25, J); print('J', J, 'mack_mat %.1fs' % (time.time() - t), flush=True)

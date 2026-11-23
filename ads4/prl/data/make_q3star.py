"""make_q3star.py -- publish Q3*, the flat target of Omega^Q_up, as exact factored coefficients.

Q3* (CERT_CLOSE.md Sec. 2.1) is Q3-full-delta0.2 (q3_full_delta0.2.json) with the factor (1-p)^3 of its weights made exact:
    psi2(p) = p (1-p)^3 X(p),    psi3(p) = (1-p)^3 Y(p),
where X, Y are the exact quotients of the stored IEEE doubles (read as exact binary rationals) by p(1-p)^3 and (1-p)^3
(synthetic division by (p-1) three times; the remainders, ~1e-11 .. 1e-8, are dropped).  This is definitional: the stored
doubles carry (1-p)^3 only to rounding (psi2(1) = 1.8e-11), which makes the flat action of the stored weights negative at +-,
m = M, J >= 1990; Q3* has minimum +0.2000 there.  A, kappa, a, c are unchanged.

Standalone (Python 3 standard library only); reproduces the fingerprint printed by ads4/code/cert/q_cert.py (load_q):
    sha256("psi2=<coeffs>;psi3=<coeffs>;a=<a>;c=<c>") = 80e0834a6b64da575b730cc6f6d3e8102aa27435e42e9fc1d422c094ba3833f2
Usage (from ads4/prl/data):  python3 make_q3star.py   -> q3star_factored.json
"""
import json, hashlib, os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'q3_full_delta0.2.json')
OUT = os.path.join(HERE, 'q3star_factored.json')
FINGERPRINT = '80e0834a6b64da575b730cc6f6d3e8102aa27435e42e9fc1d422c094ba3833f2'


def pmul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def factor(P):
    """P = p^e (p-1)^3 q + rem (e = 1 iff P[0] == 0). Returns (quotient Z with P ~ p^e (1-p)^3 Z, remainders)."""
    e = 1 if P[0] == 0 else 0
    q = list(P[e:]); rems = []
    for _ in range(3):
        out = [Fr(0)] * (len(q) - 1); acc = Fr(0)
        for i in reversed(range(1, len(q))):
            acc += q[i]; out[i - 1] = acc
        rems.append(acc + q[0]); q = out
    return e, [-x for x in q], rems          # (p-1)^3 = -(1-p)^3


d = json.load(open(SRC))
P2 = [Fr(float(x)) for x in d['psi2']]; P3 = [Fr(float(x)) for x in d['psi3']]
a = Fr(float(d['a'])); c = Fr(float(d['c']))
e2, X, r2 = factor(P2); e3, Y, r3 = factor(P3)
assert (e2, e3) == (1, 0)
omp3 = [Fr(1), Fr(-3), Fr(3), Fr(-1)]
psi2 = [Fr(0)] + pmul(omp3, X); psi3 = pmul(omp3, Y)
canon = 'psi2=' + ','.join(map(str, psi2)) + ';psi3=' + ','.join(map(str, psi3)) + ';a=%s;c=%s' % (a, c)
fp = hashlib.sha256(canon.encode()).hexdigest()
assert fp == FINGERPRINT, fp
A = 2 * psi2[1]
dev2 = max(abs(float(x - y)) for x, y in zip(psi2, P2 + [Fr(0)] * (len(psi2) - len(P2))))
dev3 = max(abs(float(x - y)) for x, y in zip(psi3, P3 + [Fr(0)] * (len(psi3) - len(P3))))
out = {
    'name': 'Q3* (Q3-full-delta0.2 with exact (1-p)^3 factor)',
    'source': 'q3_full_delta0.2.json (sha256[:16] %s), by make_q3star.py' % hashlib.sha256(open(SRC, 'rb').read()).hexdigest()[:16],
    'definition': 'psi2(p) = p (1-p)^3 X(p), psi3(p) = (1-p)^3 Y(p), X(p) = sum_n X[n] p^n, Y(p) = sum_n Y[n] p^n (lowest power '
                  'first, exact rationals as strings); flat action as in q3_full_delta0.2.json with a, c unchanged',
    'X_exact': [str(x) for x in X], 'Y_exact': [str(y) for y in Y],
    'X_float': [repr(float(x)) for x in X], 'Y_float': [repr(float(y)) for y in Y],
    'psi2_expanded_exact': [str(x) for x in psi2], 'psi3_expanded_exact': [str(x) for x in psi3],
    'a_exact': str(a), 'c_exact': str(c),
    'A': repr(float(A)), 'A_stored': d['A'], 'kappa_stored': d['kappa'],
    'dropped_remainders_psi2': ['%.3e' % float(x) for x in r2], 'dropped_remainders_psi3': ['%.3e' % float(x) for x in r3],
    'max_coefficient_change_psi2': '%.3e' % dev2, 'max_coefficient_change_psi3': '%.3e' % dev3,
    'sha256_canonical': fp,
    'canonical_string_rule': 'psi2=<psi2_expanded_exact joined by ,>;psi3=<psi3_expanded_exact joined by ,>;a=<a_exact>;c=<c_exact>',
}
json.dump(out, open(OUT, 'w'), indent=1)
print('Q3*: fingerprint OK', fp)
print('A = 2 psi2\'(0) =', float(A), '(stored', d['A'], ')')
print('remainders psi2', out['dropped_remainders_psi2'], 'psi3', out['dropped_remainders_psi3'])
print('max coefficient change psi2 %s, psi3 %s' % (out['max_coefficient_change_psi2'], out['max_coefficient_change_psi3']))
print('deg X =', len(X) - 1, ' deg Y =', len(Y) - 1, '; wrote', OUT)

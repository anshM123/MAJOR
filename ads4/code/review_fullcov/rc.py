"""Independent review helpers (referee of sliver_fullcov PROOF.md).

Independent fixed-point evaluator of
    F_K(c,h) = (1-c)^3 sum_{k<=K} (A_k + h T_{k,b}) U_k(c,h)
from the compact coefficient centres (radius <= 2e-60 ignored: it is
far below every quantity tested here).  c,h are exact Fractions; the
evaluation is absolute fixed point with W fractional bits, error
<= 4(K+1) 2^-W whenever all |ratios| <= 1.
"""
import hashlib, json, os
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / 'sliver_fullcov'
MERGE = ROOT / 'sliver_merge'
COEFF = os.environ.get('SLIVER_COEFF',
    '/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/fullcov/sliver_merge/coeff_compact60.txt')
SHA = '811e15cdd95319309edf64f5cf862bb674f43fbc73f551a002577c4b56b1ee06'

_cache = {}

def coeffs():
    if 'a' in _cache:
        return _cache['a'], _cache['g']
    raw = open(COEFF, 'rb').read()
    assert hashlib.sha256(raw).hexdigest() == SHA, 'coefficient SHA mismatch'
    lines = raw.decode().split('\n')
    N, D = map(int, lines[0].split()); assert D == 60
    a = []; g = []
    for ln in lines[1:N + 2]:
        x, r, y, t = map(int, ln.split()); a.append(x); g.append(y)
    _cache['a'], _cache['g'] = a, g
    return a, g

def tcoef_num(k, b):
    """2*T_k*10^60 as an integer (exact from centres)."""
    a, g = coeffs()
    q = b * (b + 1); sg = 1 if b == 0 else -1
    k1 = (k + 1) ** 2; k2 = -k * k + (8 - 4 * b) * k + 9 - 2 * b * (b - 1) + 2 * q
    return k1 * a[k + 1] + k2 * a[k] + 2 * sg * g[k]

def fixed_tables(K, b, W):
    key = (K, b, W)
    if key in _cache: return _cache[key]
    a, g = coeffs(); S = 10 ** 60
    Af = [(a[k] << W) // S for k in range(K + 1)]
    Tf = [(tcoef_num(k, b) << W) // (2 * S) for k in range(K + 1)]
    _cache[key] = (Af, Tf)
    return Af, Tf

def F_K(c, h, K, b, W=256, prefactor=True):
    """Return F_K(c,h) as a Fraction (fixed point, abs error ~ 4(K+1)2^-W)."""
    c = Fr(c); h = Fr(h)
    s = 2 * b + 1
    Af, Tf = fixed_tables(K, b, W)
    cn, cd = c.numerator, c.denominator; hn, hd = h.numerator, h.denominator
    v = 0
    for k in range(K, -1, -1):
        if k < K:
            p = k * (k + s)
            num = cn * hd - p * hn * cd       # (c - p h) * cd*hd
            den = cd * (hd - p * hn)          # (1 - p h) * cd*hd
            assert den != 0
            v = (v * num) // den   # floor: error < 1 unit
        v += Af[k] + (Tf[k] * hn) // hd
    val = Fr(v, 1 << W)
    if prefactor:
        val *= (1 - c) ** 3
    return val

def cell(stem):
    z = json.load(open(MERGE / (stem + '.json')))
    l, u, Hl, Hu = map(Fr, (z['l'], z['u'], z['Hl'], z['Hu']))
    return z, l, u, Hl, Hu

def Hsec(c, l, u, Hl, Hu):
    return Hl + (Hu - Hl) * (Fr(c) - l) / (u - l)

def ledger(stem):
    for suf in ('_ref_ledger.json', '_group_ledger.json', '_short_group_ledger.json'):
        p = MERGE / (stem + suf)
        if p.exists():
            return json.load(open(p)), suf
    raise FileNotFoundError(stem)

def Rval(P, x, y):
    return sum(Fr(P[a][d]) * Fr(x) ** a * Fr(y) ** d for a in range(3) for d in range(3)) / 10 ** 20

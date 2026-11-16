"""
Parser for the ancillary file of arXiv:2106.10274 (data/chmrs_regge_moments_order8.m, copied unchanged from the
arXiv source).  It defines
   seriesPi : large-m expansion of  a_{Dphi} m^{2k} / gamma_{2 Dphi+k-1}(nu)^2 * Pihat_{k,nu}[Delta,J]
   seriesC  : large-m expansion of  m^{2k} C_{k,nu}[Delta,J]
both as  2 P + c_1/m^2 + c_2/m^4 + c_3/m^6 + c_4/m^8,   P = P_{(2-d)/2+i nu}(eta_AdS), with derivatives pj[.,.,n].
In the file  alpha := 2 Dphi + k - 1  (established by comparison with (3.27), see __main__).
"""
import re, os, sys
import mpmath as mp

_HERE = os.path.dirname(os.path.abspath(__file__))
_cache = {}

def _load(name):
    if name in _cache:
        return _cache[name]
    txt = open(os.path.join(_HERE, 'data', 'chmrs_regge_moments_order8.m')).read()
    start = txt.index({'Pi': 'series\\[CapitalPi] = ', 'C': 'seriesC = '}[name])
    end = txt.index({'Pi': 'seriesC = ', 'C': 'seriesCfixedJ = '}[name])
    body = txt[start:end]
    body = body[body.index('=')+1:].strip().rstrip(';').strip()
    body = body.replace('\\[Eta]AdS', 'eta').replace('\\[Alpha]', 'al').replace('\\[CapitalDelta]\\[Phi]', 'Dphi').replace('\\[Nu]', 'nu')
    body = re.sub(r'\s+', ' ', body)
    body = re.sub(r'pj\[\(2 - d\)/2 \+ I\*nu,\s*eta,\s*(\d)\]', r'pj[\1]', body)
    # split into orders: "... )/m2 + ( ... )/m2^2 + ..."
    parts = re.split(r'/m2(?:\^(\d))?\s*\+?', body)
    # parts = [chunk0, None, chunk1, '2', chunk2, ...]
    chunks = parts[0::2]; orders = parts[1::2]
    first = chunks[0]
    lead, rest = first.split('+', 1)          # "2*pj[0] " + " (....)"
    out = {0: lead.strip()}
    cur = rest.strip()
    for i, o in enumerate(orders):
        out[1 if o is None else int(o)] = cur
        cur = chunks[i+1].strip() if i+1 < len(chunks) else ''
    for key in out:
        out[key] = out[key].replace('^', '**')
    _cache[name] = out
    return out

def series(name, k, nu, eta, Dphi, d=3, max_order=2):
    """Returns the list [c_0, c_1, ..., c_max_order] of the coefficients of 1/m2^n (mpmath numbers)."""
    src = _load(name)
    j = mp.mpf(2-d)/2 + 1j*nu
    f = lambda x: mp.hyp2f1(-j, j+d-2, mp.mpf(d-1)/2, (1-x)/2)
    pj = [f(eta)] + [mp.diff(f, eta, n) for n in range(1, 2*max_order+1)]
    env = dict(d=mp.mpf(d), k=mp.mpf(k), al=2*Dphi + k - 1, Dphi=mp.mpf(Dphi), nu=mp.mpf(nu), eta=mp.mpf(eta), pj=pj)
    sys.setrecursionlimit(100000)
    res = []
    for n in range(max_order+1):
        code = re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'mpf(\1)', src[n])     # exact arithmetic for integers
        code = re.sub(r'pj\[mpf\((\d)\)\]', r'pj[\1]', code)
        res.append(eval(code, dict(mpf=mp.mpf, I=mp.mpc(0, 1)), env))
    return res

if __name__ == "__main__":
    from chmrs_formulas import C_series_first, m2_eta
    mp.mp.dps = 25
    print("# seriesC first order vs CHMRS (3.27)  [alpha = 2 Dphi + k - 1; the printed (Dphi-d-1)_3 must read (2Dphi-d-1)_3]; and seriesPi vs seriesC at orders 1, 2")
    for (k, nu, Delta, J, Dphi, d) in ((2, mp.mpf('1.3'), mp.mpf('41.7'), 14, mp.mpf(3), 3), (4, mp.mpf('0.4'), mp.mpf('30.2'), 6, mp.mpf('1.6'), 3),
                                       (2, mp.mpf('2.1'), mp.mpf('55.5'), 20, mp.mpf('2.2'), 4)):
        m2, eta = m2_eta(Delta, J, d)
        cC = series('C', k, nu, eta, Dphi, d, 2); cP = series('Pi', k, nu, eta, Dphi, d, 2)
        P0, P01 = C_series_first(k, nu, Delta, J, Dphi, d)
        P0, P01p = C_series_first(k, nu, Delta, J, Dphi, d, printed=True)
        print(f"k={k} nu={nu} d={d}: (3.27) as printed {mp.nstr(((P01p-P0)*m2*2).real, 15)}; with (2Dphi-d-1)_3 {mp.nstr(((P01-P0)*m2*2).real, 15)}  file seriesC c1 {mp.nstr(cC[1], 15)}  seriesPi c1 {mp.nstr(cP[1], 15)}")
        print(f"      c2: seriesC {mp.nstr(cC[2], 15)}   seriesPi {mp.nstr(cP[2], 15)}")

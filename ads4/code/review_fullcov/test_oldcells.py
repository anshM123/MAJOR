"""Random-point sanity for the three earlier norming cells merged into the ledger
(sliver_close/cell_30_35, sliver_fast/cell_3500_4462, cell_4462_5305; Chebyshev-Lobatto
grid, constant midpoint reference, C=5/2): every sampled F_K must lie within
[R - (R - prefix_lower), R + (R - prefix_lower)] i.e. |F_K - R| <= certified residual."""
import json, random
from fractions import Fraction as Fr
import rc
random.seed(99)
for d, stem in [('sliver_close', 'cell_30_35'), ('sliver_fast', 'cell_3500_4462'), ('sliver_fast', 'cell_4462_5305')]:
    z = json.load(open(rc.ROOT / d / (stem + '.json'))); led = json.load(open(rc.ROOT / d / (stem + '_ledger.json')))
    l, u, Hl, Hu = map(Fr, (z['l'], z['u'], z['Hl'], z['Hu'])); K = z['K']
    for r in led['results']:
        b = r['b']; R = Fr(r['reference']); L = Fr(r['prefix_lower']); M = R - L; mx = 0
        pts = [(Fr(i), Fr(j)) for i in (0, 1) for j in (0, 1)] + [(Fr(random.random()), Fr(random.random())) for _ in range(16)]
        for x, y in pts:
            c = l + (u - l) * x; f = rc.F_K(c, rc.Hsec(c, l, u, Hl, Hu) * y, K, b); mx = max(mx, abs(f - R))
        print(f'{stem} b={b} K={K}: max|F_K-R|={float(mx):.3e} certified residual={float(M):.3e} prefix_lower={float(L):.4e}', 'OK' if mx <= M else 'FAIL')

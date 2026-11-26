"""lemU_cert/lu_cell_check.py [validation, DV]: the cell Taylor model of lu_cell.py against direct point evaluations.
G is replaced by a single-frequency F_{3,nu} (Y^A = its jets, Y^D = 0), M small; the w-polynomial P(w, y) of lu_cell is evaluated
at w = +-delta and compared with the direct (lu_core) polynomial at c(v0 + w), on h = hbar y.
Usage: lu_cell_check.py u0 du nu M"""
import sys, os, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, arb_poly, ctx
import lu_core as C, lu_cert1 as C1
u0, du, nu, Mt = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
C1.M = Mt; C1.K = 3*Mt + 3
import lu_cell as LC
LC.M = Mt
captured = {}
def fake_anchor(c0, prec, pr, Kj):
    ctx.prec = prec
    y = C.jets(3, arb(nu), c0, Kj, mid=True)
    import lu_bounds as B
    env = C1.env_floats(B.Envelope(c0))
    captured['c0'] = c0
    return [arb(v) for v in y], [arb(0)]*(Kj + 1), dict(F=1.0, env=env)
LC.anchor = fake_anchor
# capture the per-k columns by wrapping C.cheb
orig_cheb = C.cheb
cols_seen = []
def cheb_spy(pm, mid=False):
    cols_seen.append(pm)
    return orig_cheb(pm, mid)
C.cheb = cheb_spy
os.environ['LU_CHANNELS'] = '0'
res = LC.cell(u0 - du, u0 + du, 600, 3*Mt + 5, 3*Mt + 5, 40, None, 6.0, 600, KC=200, Jextra=60)
# reconstruct P(w, y) at y = 1 from captured columns (k = 0.. ) and compare with direct evaluation
ctx.prec = 600
c0 = captured['c0']; d0 = 1 - c0; v0 = d0/c0
cols = cols_seen[:len(cols_seen)]
s, q_b, sig, kb, e = C.chan(0)
c_lo = (arb(u0 - du)/2).tanh()**2
hbar = arb((C.xstar(c_lo, s)*(1 - c_lo)**2).upper())
for w in (arb(0), arb(res[0]['delta_w']), -arb(res[0]['delta_w'])):
    Pw = sum((sum(col, arb(0))*w**k for k, col in enumerate(cols)), arb(0))     # y = 1
    v = v0 + w; c = 1/(1 + v); d = 1 - c
    Kd = 3*Mt + 3
    YA = C.jets(3, arb(nu), c, Kd, mid=True); YD = [arb(0)]*(Kd + 1)
    YB = C.dress(YA, YD, c, 0)
    xs = hbar/(d*d)                       # x such that h = x d^2 = hbar
    (Am, _), (Bm, _) = C.contract(c, xs, s, [YA, YB], Mt)
    pm = C.poly_y(Am, Bm, c, xs, 0)
    direct = sum(pm, arb(0))
    print('w =', float(w.mid()), ' cell model P(w, y=1) =', Pw.str(15), ' direct =', direct.str(15))

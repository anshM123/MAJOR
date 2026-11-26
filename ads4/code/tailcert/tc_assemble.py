"""
tailcert/tc_assemble.py (TAIL_CERT.md §3): the k = 2 contour tail, all connection terms j < J0, certified for L >= L0.

H^(2) = (2/Lam) w^3 sum_{j<J0} I_j + (remainder, tc_remainder.py),   I_j = 2 Re int_0^oo f B c_j V_b(E_{alpha-j}) dnu,
deformed (each j separately) to: quarter circle at 0 (residue f(0) X_j(0)/2), segment nu = iy (0 <= y <= 3/4) with integrand W(y) X_j(iy),
short detour t + 3i/4 (0 <= t <= 1/4), vertical 1/4 + iy (3/4 <= y <= 9/4), leg t + 9i/4 (t >= 1/4).
X_j(nu) := w^{alpha-j} V_b(E_{alpha-j}),  |w^3 f B c_j V_b(E_{alpha-j})| = sqrt(w) e^{-L Im nu} |f B| |c_j| w^j |X_j|.
Main term (j = 0, segment):  M(L) >= (2/Lam) sqrt(w) int W e^{-Ly}(1 - delta_seg) dy   (W >= 0, CONTOUR cb_w cells; f(0) X_0(0) >= 0 dropped).
Bounds on X_j:  X_0 via tc_power (Lemma NB on segment/detour, Theorem P on the leg, Ŝ(alpha) and Ŝ(alpha+1)), dressing exact,
                j >= 1 via the exact recurrence (tc_exact E3): |Ŝ(b-1)| <= (1+1e-6)|Ŝ(b)| + |b(b+1)|/(3.9 a^2) |Ŝ(b+1)|.
f-integrals: CONTOUR certified rows (cb_short, cb_vert, cb_horiz_2.25: per-cell upper bounds of int|fB|), leg tail t >= 1054 recomputed.
Output: cache/tc_assemble_L<L0>.txt / .json
Usage: python3 tc_assemble.py L0 [J0]
"""
import os, sys, json, math, time
from flint import arb, acb, ctx
ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import tc_power as P
CB = os.path.join(ADS4, 'code', 'contour', 'cache')
A_MIN = 998
LAM = math.sqrt(998*999)


def up(x):
    return float(arb(x).upper())


def cj_bounds(nu_box, J0):
    """upper bounds of |c_j(nu)| over the box, j = 0..J0-1."""
    out = [1.0]
    c = acb(1)
    a = acb(0.5) - acb(0, 1)*nu_box
    b = acb(1) - 2*acb(0, 1)*nu_box
    for j in range(J0 - 1):
        c = c*(a + j)*(a + j)/((b + j)*(j + 1))
        out.append(up(abs(c)))
    return out


def jsum(m0, mm1, absalpha, cbs, w, J0, real=False):
    """sum_{j>=1} |c_j| w^j |X_j| given m0 >= |Ŝ(alpha)|, mm1 >= |Ŝ(alpha+1)|, |alpha| <= absalpha.  Includes dressing.
    Returns (sum_j>=1, |X_0| bound)."""
    a2 = 3.9*A_MIN**2
    m = [mm1, m0]          # m[-2], m[-1]
    def Xb(mj, mjm1, bmod):
        A1 = bmod*(bmod + 1)/2 + 4*bmod
        A0 = bmod*bmod/2 + 5*bmod + 8
        return mj*(1 + A0*w*w/a2) + (w/a2)*A1*mjm1
    X0 = Xb(m0, mm1, absalpha)
    tot = 0.0
    for j in range(1, J0):
        bm = absalpha + (j - 1)            # |beta| for beta = alpha - (j-1)
        mj = (1 + 1e-6)*m[-1] + bm*(bm + 1)/a2*m[-2]
        m.append(mj)
        term = cbs[j]*w**j*Xb(mj, m[-2], absalpha + j)
        tot += term
        if term < 1e-300 and j > 5:
            break
    return tot*(1 + 1e-9), X0


def main():
    L0 = float(sys.argv[1]); J0 = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    w = math.exp(-L0)
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    t00 = time.time()
    # ---------------- segment: W cells
    Wd = json.load(open(os.path.join(CB, 'cb_w.json')))
    cells = Wd['cells']                  # [y0, y1, minW, maxW] on [0, 3/4]
    nb25 = P.NB_bound(2.5, 0); 
    seg = []                             # (y0, y1, minW, maxW, delta_seg(y), jsum(y))
    for y0, y1, mn, mx in cells:
        sig = 2.5 - y1                   # beta in [5/2 - y1, 5/2 - y0], real
        nbb = max(P.NB_bound(2.5 - y0, 0), P.NB_bound(2.5 - y1, 0))
        nbb1 = max(P.NB_bound(3.5 - y0, 0), P.NB_bound(3.5 - y1, 0))
        m0 = 1 + nbb; mm1 = 1 + nbb1
        cbs = cj_bounds(acb(0, arb((y0 + y1)/2, (y1 - y0)/2)), J0)
        js, X0 = jsum(m0, mm1, 2.5, cbs, w, J0)
        # X_0 >= 1 - nbb - dressing(|A1| w/(3.9a^2) mm1)  ; |A1| <= 2.5*3.5/2 + 4*2.5
        dseg = nbb + (w/(3.9*A_MIN**2))*(2.5*3.5/2 + 10)*mm1 + (2.5**2/2 + 12.5 + 8)*w*w/(3.9*A_MIN**2)*m0
        seg.append((y0, y1, mn, mx, dseg, js))
    pr(f'segment: {len(seg)} W-cells; max delta_seg = {max(s[4] for s in seg):.4e}; max j>=1 factor = {max(s[5] for s in seg):.3e}')
    # residue j>=1 part: f(0) = r(0) (EPSW) times sum_j c_j(0) w^j |X_j| / 2
    ew = json.load(open(os.path.join(ADS4, 'code', 'epswcert', 'cache', 'ew_cert_1000_re_rep.json')))
    f0 = ew['r0_hi']*(1 + 1e-12)
    js0, _ = jsum(1 + nb25, 1 + P.NB_bound(3.5, 0), 2.5, cj_bounds(acb(0), J0), w, J0)
    res_err = f0/2*js0
    # ---------------- detour short: nu = t + 3i/4, beta = 7/4 + i t (t in [0, 1/4]); vertical: nu = 1/4 + iy, beta = 5/2 - y + i/4
    def load(fn):
        return json.load(open(os.path.join(CB, fn)))
    sh = load('cb_short.json')['rows']
    short = []
    for r in sh:
        t0, t1 = r['t0'], r['t1']
        m0 = 1 + P.NB_bound(1.75, 0.25); mm1 = 1 + P.NB_bound(2.75, 0.25)
        cbs = cj_bounds(acb(arb((t0 + t1)/2, (t1 - t0)/2), 0.75), J0)
        js, X0 = jsum(m0, mm1, abs(complex(1.75, 0.25)), cbs, w, J0)
        short.append((r['intfB'], X0 + js))
    ve = load('cb_vert.json')['rows']
    vert = []
    for r in ve:
        y0, y1 = r['y0'], r['y1']
        sl = 2.5 - y1
        m0 = 1 + max(P.NB_bound(2.5 - y0, 0.25), P.NB_bound(sl, 0.25))
        mm1 = 1 + max(P.NB_bound(3.5 - y0, 0.25), P.NB_bound(sl + 1, 0.25))
        cbs = cj_bounds(acb(0.25, arb((y0 + y1)/2, (y1 - y0)/2)), J0)
        js, X0 = jsum(m0, mm1, abs(complex(2.5 - y0, 0.25)), cbs, w, J0)
        vert.append((y0, r['intfB'], X0 + js))
    pr(f'detour: short {len(short)} cells (max |X|-factor {max(s[1] for s in short):.6f}), vertical {len(vert)} cells '
       f'(max {max(v[2] for v in vert):.6f})  [{time.time()-t00:.0f}s]')
    # ---------------- leg Y = 9/4: beta = 1/4 + i t
    hz = load('cb_horiz_2.25.json')
    leg = []
    for i, r in enumerate(hz['rows']):
        t0, t1 = r['t0'], r['t1']
        m0, _ = P.K_osc(0.25, t0, t1, L0)
        mm1, _ = P.K_osc(1.25, t0, t1, L0)
        cbs = cj_bounds(acb(arb((t0 + t1)/2, (t1 - t0)/2), 2.25), J0)
        js, X0 = jsum(m0, mm1, abs(complex(0.25, t1)), cbs, w, J0)
        leg.append((t0, t1, r['intfB'], X0, js))
        if i % 2000 == 0:
            pr(f'   leg cell {i}/{len(hz["rows"])} t = {t0:.2f}: |X0| <= {X0:.4g}, j-sum {js:.3e} [{time.time()-t00:.0f}s]')
    # leg tail t >= T_END: per unit cell, int_{t1}^oo |f| * sup|B| (CONTOUR's tail formula, monotone in T) times K
    import importlib.util
    sys.path.insert(0, os.path.join(ADS4, 'code', 'contour'))
    import cb_common as K
    import cb_lines as CL
    W = K.load()
    TEND = hz['summary']['T_END']
    tail = 0.0
    T = TEND
    while T < 3000:
        T2 = T + 1
        intf, Bsup = CL.tail(W, 2.25, T)
        l0 = P.logK_osc(0.25, T, T2, L0); l1 = P.logK_osc(1.25, T, T2, L0)
        cbs = cj_bounds(acb(arb(T + 0.5, 0.5), 2.25), J0)
        js, X0 = jsum(1.0, math.exp(min(l1 - l0, 600)) if l1 > l0 else 1.0, abs(complex(0.25, T2)), cbs, w, J0)
        lt = float((intf*Bsup).log().upper()) + l0 + math.log(X0 + js)
        tail += math.exp(lt) if lt > -700 else 0.0
        T = T2
    # t >= 3000: |Ŝ(sigma+it)| <= Gamma(sigma) |1/Gamma(sigma+it)| Ŝ(sigma) <= Gamma(1/4) Gamma(3/4) cosh(pi t)/pi * 1.01 (reflection);
    # sigma = 5/4 adds 1/|1/4+it| <= 1.  int_3000^oo |f| e^{pi t} <= via cb_lines.tail with weight: bounded crudely by e^{pi*T}*tail(T)
    # summed over unit cells up to 6000 and e^{-t^2/4}-dominated beyond (Gaussian erfc factors ~ e^{-(T-1029)^2/1.62}).
    T = 3000; tail2 = 0.0
    g14 = math.gamma(0.25)*math.gamma(0.75)/math.pi*1.01
    while T < 6000:
        intf, Bsup = CL.tail(W, 2.25, T)
        tail2 += up(intf*Bsup*arb(g14)*(arb(math.pi)*(T + 1)).exp()*(1 + w*(T + 2)**J0))
        T += 1
    pr(f'leg tail: [T_END={TEND}, 3000): {tail:.3e};  [3000, 6000): {tail2:.3e}  (beyond 6000: e^{{-1.6e7}})')
    # ---------------- assemble rho(L) on a grid of L >= L0 (each term decreasing in L; Main decreasing): rho(L) <= Err(Lk)/Main(L_{k+1})
    def main_lb(L):
        tot = 0.0
        for y0, y1, mn, mx, dseg, js in seg:
            # int_{y0}^{y1} e^{-Ly} = (e^{-L y0} - e^{-L y1})/L
            tot += mn*(1 - dseg)*(math.exp(-L*y0) - math.exp(-L*y1))/L
        return tot*(1 - 1e-12)
    def err(L):
        e_seg = sum(mx*js*(math.exp(-L*y0) - math.exp(-L*y1))/L for y0, y1, mn, mx, dseg, js in seg)
        e_res = res_err
        e_short = 2*math.exp(-0.75*L)*sum(i*x for i, x in short)
        e_vert = 2*sum(i*x*math.exp(-L*y0) for y0, i, x in vert)
        e_leg = 2*math.exp(-2.25*L)*(sum(i*(X0 + js) for t0, t1, i, X0, js in leg) + tail + tail2)
        return e_seg, e_res, e_short, e_vert, e_leg
    pr(f'leg sum int|fB|(|X0| + j-sum) over CONTOUR cells: {sum(i*(X0+js) for t0,t1,i,X0,js in leg):.4e} '
       f'(vs int|fB| = {sum(i for t0,t1,i,X0,js in leg):.4e}; t<45 part {sum(i*(X0+js) for t0,t1,i,X0,js in leg if t1<=45.0001):.4e})')
    grid = [L0 + 0.05*k for k in range(int((80 - L0)/0.05) + 1)]
    rows = []
    for La, Lb in zip(grid[:-1], grid[1:]):
        e = err(La); m = main_lb(Lb)
        rows.append((La, Lb, sum(e)/m, [x/m for x in e]))
    worst = max(r[2] for r in rows)
    pr(f'L in [{L0}, 80] (cells 0.05): max rho = {worst:.4e}')
    for r in rows[::40]:
        pr(f'  L in [{r[0]:.2f},{r[1]:.2f}]: rho <= {r[2]:.3e}   (seg_j {r[3][0]:.1e}, res {r[3][1]:.1e}, short {r[3][2]:.1e}, '
           f'vert {r[3][3]:.1e}, leg {r[3][4]:.1e})')
    # beyond L = 80: Main >= (int_0^{3/4} minW e^{-Ly}) >= 0.0586 (1 - e^{-0.01 L})/L * (1-dseg); errors <= C e^{-3L/4} (all
    # e^{-L Im nu} with Im nu >= 3/4 off the segment) + seg_j (w^1 factor) + res (w factor):  ratio decreasing -> checked at 80.
    json.dump(dict(L0=L0, J0=J0, rows=rows, worst=worst), open(os.path.join(HERE, 'cache', f'tc_assemble_L{L0}.json'), 'w'))
    open(os.path.join(HERE, 'cache', f'tc_assemble_L{L0}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()

"""
rt2_resid.py La Lb [part nparts]: residual term E_r of Route T'' for the L-cell [La, Lb] (ANALYTIC_SLIVER2 §6), as a ratio to F:
  E_r = d^3 (2/Lam) [ 2 sum_bands m_band sup_band |K_J| + 2 sup_{nu<Lam}|K_J| * 1.66e-23 (theta-cutoff, G1) ]  +  E_far,
  (for nu >= Lam psi2 = 0, r_theta = f_theta and |f_theta| <= |f| = |r| pointwise, so the EPSW masses already bound int|r_theta|)
  E_far = d^3 (2/Lam) int_{1150}^oo |f_theta| |V_b K(nu)| dnu,  |V_b K| <= 1.6 * 41 * (pi/2) e^{pi nu/2} (0.235 d_lo)^{-3}
          (no j-splitting above 1150: Loc' representation, K0 <= 1.6, dressing 41, FULLPOS3 Lemma B), with the e^{pi nu/2}-weighted
          CONTOUR tail (wtail below).
Bands: [0, 0.05] and [0.05, 0.4] as rt_cert.residual (max principle / old per-power bounds); [0.4, 40], [40, 898] (EPSW range masses)
by rt2_kernel.K_cell on sub-cells of width <= 2; [898, 1059) unit cells (rt_rmass2 masses) by K_cell; [1059, 1150] (EPSW Gaussian tail
mass) by K_cell without frames.  Part p of nparts processes the cell list round-robin; 'assemble' merges.
"""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common2
import rt_cert as RC
import rt2_kernel as RK
from flint import arb, acb
HERE = common2.HERE
JP = int(os.environ.get('RT2_JP', '1500'))


def wtail(T, w=math.pi/2):
    """upper bound of int_T^oo |f_theta(t)| e^{w t} dt (real axis), CONTOUR cb_lines.tail with the weight folded into each Gaussian."""
    sys.path.insert(0, os.path.join(common2.AN, '..', 'contour'))
    import cb_common as KC
    W = KC.load()
    ymax = arb(1); wa = arb(w)
    th = (ymax*ymax/(W.acut*W.acut)).exp()
    R0 = arb(T)
    qbar = 32*R0*(2*R0 + 3)*(2*R0 + 5)*(2*R0 + 11)/((2*R0 - 1)*(2*R0 - 7)*(2*R0 - 3)*(2*R0 - 7))
    tot = arb(0)
    for (c, s, x), (ca, sa, xa, inv) in zip(W.sc, W.sc_arb):
        z = (R0 - ca - wa*sa*sa)/(sa*arb(2).sqrt())
        assert z > 0
        tot += abs(xa)*(ymax*ymax*inv).exp()*(wa*ca + wa*wa*sa*sa/2).exp()*sa*(arb.pi()/2).sqrt()*z.erfc()
    ysum = sum((abs(y) for y in W.y), arb(0))
    m2 = R0 + ymax + 2
    Pb = ysum*m2**4*(m2*m2 + 1)*(2*m2.log() + arb.pi() + 1)*2*(2*m2 + 11)**6/(R0*(2*R0 - 7)**3)
    assert Pb*(ymax*ymax/2).exp() < (R0*R0/4).exp()
    ytail = (wa*wa).exp()*2/(R0 - 2*wa)*(-((R0 - 2*wa)**2)/4).exp()
    return 2*th*(qbar*tot + ytail)


FINE = os.environ.get('RT2_FINE', '0') == '1'


def cell_list():
    if FINE:
        # per-unit-cell EPSW panel masses on [1, 898) (rmass3.py, same machinery as rt_rmass2); [0.4, 1] with EPSW's mass of [0, 1]
        cells = [('band', 0.4, 1.0, 1.1063e-6, 0.6)]
        for fn in ('rmass3_1_40.json', 'rmass3_40_898.json'):
            for (T0, T1, m) in json.load(open(os.path.join(HERE, 'cache', fn))):
                cells.append(('cell', float(T0), float(T1), m, 2.0))
    else:
        cells = [('band', 0.4, 40.0, 1.1063e-6 + 6.0850e-7, 2.0), ('band', 40.0, 898.0, 1.7690e-7, 2.0)]
    rm = json.load(open(os.path.join(common2.AN, 'cache', 'rt_rmass2.json')))
    for (T0, T1, m) in rm:
        if T1 > 1e8:
            cells.append(('tail', 1059.0, 1150.0, m, 13.0))
        else:
            cells.append(('cell', float(T0), float(T1), m, 2.0))
    return cells


def run_part(La, Lb, part, nparts):
    cells = cell_list()
    out = []
    t0 = time.time()
    for idx, (kind, lo, hi, m, w) in enumerate(cells):
        nsub = max(1, int(math.ceil((hi - lo)/w - 1e-9)))
        subs = [(lo + (hi - lo)*i/nsub, lo + (hi - lo)*(i + 1)/nsub) for i in range(nsub)]
        for (a_, b_) in subs:
            key = (idx, a_)
            if hash(key) % nparts != part and nparts > 1:
                pass
        mine = [sb for k, sb in enumerate(subs) if (idx + k) % nparts == part]
        if not mine:
            continue
        best = 0.0
        for (a_, b_) in mine:
            jm = RK.JMAX if La >= 3.9 else 110
            Jx = JP if a_ >= 999.0 else RK.J0            # J(nu) = J' for nu >= 999 (rt2_loc2), the contour's J below
            v, ar = RK.K_cell(a_, b_, La, Lb, use_frames=(kind != 'tail'), jmax=jm, J0x=Jx)
            best = max(best, v)
        out.append(dict(idx=idx, kind=kind, lo=lo, hi=hi, m=m, K=best))
        print(f'{kind} [{lo},{hi}] m={m:.3e}: sup|K_J| <= {best:.4e}   [{time.time()-t0:.0f}s]', flush=True)
    json.dump(out, open(os.path.join(HERE, 'cache', f'resid_{La}_{Lb}_p{part}of{nparts}.json'), 'w'))


def small_nu(La):
    """[0, 0.05] and [0.05, 0.4] exactly as rt_cert.residual (returns the sum of 2 K m)."""
    d = math.exp(-La)
    R0 = 0.05; nb = 64; best = 0.0
    for k in range(nb):
        th = 2*math.pi*(k + 0.5)/nb
        mid = complex(R0*math.cos(th), R0*math.sin(th)); r = R0*math.pi/nb*1.5
        Bv = RC.B_abs_cell(mid, r)
        box = acb(arb(mid.real, r), arb(mid.imag, r))
        re_hi = 2.5 + abs(mid.imag) + r; re_lo = 2.5 - abs(mid.imag) - r
        js = RC.jsum_near(box, (re_lo, re_hi), abs(mid.real) + r, 600, d, 3, True)
        cs = RC.cj_list(box, 60)
        js1 = sum(cs[j]*d**(j - re_hi) for j in range(60))*1.01
        best = max(best, Bv*(js + js1))
    t1 = 2*best*7.8938e-7
    K = RC.K_sup_real(0.05, 0.4, La, 14)
    t2 = 2*K*1.1063e-6
    return t1, t2


def assemble(La, Lb, nparts):
    d = math.exp(-La); d_lo = math.exp(-Lb)
    rows = []
    for p in range(nparts):
        rows += json.load(open(os.path.join(HERE, 'cache', f'resid_{La}_{Lb}_p{p}of{nparts}.json')))
    byidx = {}
    for r in rows:
        if r['idx'] not in byidx or r['K'] > byidx[r['idx']]['K']:
            byidx[r['idx']] = r
    cells = cell_list()
    assert sorted(byidx) == list(range(len(cells))), 'missing cells'
    F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    pref = d**3*(2/RC.LAMf)
    t1, t2 = small_nu(La)
    tot = t1 + t2; Ksup = 0.0
    parts = [('[0,0.05]', t1), ('[0.05,0.4]', t2)]
    worst = []
    for i in sorted(byidx):
        r = byidx[i]
        tot += 2*r['K']*r['m']
        if r['lo'] < RC.LAMf:
            Ksup = max(Ksup, r['K'])     # G1 only where psi2 lives: for nu >= Lam, r_theta = f_theta and |f_theta| <= |f| (0 <= theta <= 1)
        worst.append((2*r['K']*r['m']*pref/F, r['lo'], r['hi']))
    worst.sort(reverse=True)
    g1 = 2*Ksup*1.66e-23
    wt = wtail(1150.0)
    far_log = float((arb(1.6*41*math.pi/2*1.000001).log() - 3*arb(0.235*d_lo).log() + wt.log()).upper())
    far = math.exp(max(far_log, -690.0))      # floor at e^-690 so the float stays an upper bound
    Er = pref*(tot + g1) + pref*far
    res = dict(La=La, Lb=Lb, F=F, Er_over_F=Er/F*(1 + 1e-9), small=[t1*pref/F, t2*pref/F], G1=g1*pref/F, far=pref*far/F, far_log=far_log,
               worst=worst[:8], Ksup=Ksup)
    json.dump(res, open(os.path.join(HERE, 'cache', f'resid_{La}_{Lb}.json'), 'w'))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    La, Lb = float(sys.argv[1]), float(sys.argv[2])
    if len(sys.argv) > 3 and sys.argv[3] == 'assemble':
        assemble(La, Lb, int(sys.argv[4]))
    else:
        part = int(sys.argv[3]) if len(sys.argv) > 3 else 0
        nparts = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        run_part(La, Lb, part, nparts)

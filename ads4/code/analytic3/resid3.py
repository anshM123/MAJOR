"""
resid3.py La Lb part nparts | La Lb assemble nparts:  residual term E_r of Route T''' (ANALYTIC_SLIVER3 §4) for the L-cell [La, Lb].
Identical bookkeeping to rt2_resid.py with RT2_FINE=1 (EPSW per-unit-cell masses on [1, 898), rt_rmass2 cells on [898, 1059), the EPSW
Gaussian tail mass on [1059, 1150], G1 theta-cutoff term, G2 far term with the full kernel above 1150, small-nu [0, 0.4] as rt_cert),
except that sup|K_J| per cell is kernel3.K_cell (a restricted to [998, 4500], finer a-cover).  J = 1500 everywhere (no J' split:
Lemma Loc''' needs no REC hypothesis), so the contour's kappa_J and the real-axis K_J are the same truncation.
"""
import sys, os, math, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kernel3 as K3
os.environ['RT2_FINE'] = '1'
import rt2_resid as R2
import rt_cert as RC
from flint import arb

RATIO = float(os.environ.get('R3_RATIO', '1.03'))
TURN = tuple(float(x) for x in os.environ.get('R3_TURN', '995,1059').split(','))   # nu-band with refined boxes
TURN_NU = int(os.environ.get('R3_TURN_NU', '2')); TURN_L = int(os.environ.get('R3_TURN_L', '1')); TURN_W = float(os.environ.get('R3_TURN_W', '0.25'))


def run_part(La, Lb, part, nparts):
    cells = R2.cell_list()
    out = []
    t0 = time.time()
    for idx, (kind, lo, hi, m, w) in enumerate(cells):
        nsub = max(1, int(math.ceil((hi - lo)/w - 1e-9)))
        subs = [(lo + (hi - lo)*i/nsub, lo + (hi - lo)*(i + 1)/nsub) for i in range(nsub)]
        mine = [sb for k, sb in enumerate(subs) if (idx + k) % nparts == part]
        if not mine:
            continue
        best = 0.0
        for (a_, b_) in mine:
            if TURN[0] <= a_ < TURN[1] and kind != 'tail':
                v, ar = K3.K_cell_refined(a_, b_, La, Lb, ratio=RATIO, nsub_nu=TURN_NU, nsub_L=TURN_L, w_turn=TURN_W, jmax=110, J0x=1500)
            else:
                v, ar = K3.K_cell(a_, b_, La, Lb, ratio=RATIO, use_frames=(kind != 'tail'), jmax=110, J0x=1500)
            best = max(best, v)
        out.append(dict(idx=idx, kind=kind, lo=lo, hi=hi, m=m, K=best))
        print(f'{kind} [{lo},{hi}] m={m:.3e}: sup|K_J| <= {best:.4e}   [{time.time()-t0:.0f}s]', flush=True)
    json.dump(out, open(os.path.join(HERE, 'cache', f'resid3_{La}_{Lb}_p{part}of{nparts}.json'), 'w'))


def assemble(La, Lb, nparts):
    d = math.exp(-La); d_lo = math.exp(-Lb)
    rows = []
    for p in range(nparts):
        rows += json.load(open(os.path.join(HERE, 'cache', f'resid3_{La}_{Lb}_p{p}of{nparts}.json')))
    byidx = {}
    for r in rows:
        if r['idx'] not in byidx or r['K'] > byidx[r['idx']]['K']:
            byidx[r['idx']] = r
    cells = R2.cell_list()
    assert sorted(byidx) == list(range(len(cells))), 'missing cells'
    F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    pref = d**3*(2/RC.LAMf)
    t1, t2 = R2.small_nu(La)
    tot = t1 + t2; Ksup = 0.0
    worst = []
    for i in sorted(byidx):
        r = byidx[i]
        tot += 2*r['K']*r['m']
        if r['lo'] < RC.LAMf:
            Ksup = max(Ksup, r['K'])
        worst.append((2*r['K']*r['m']*pref/F, r['lo'], r['hi']))
    worst.sort(reverse=True)
    g1 = 2*Ksup*1.66e-23
    wt = R2.wtail(1150.0)
    far_log = float((arb(1.6*41*math.pi/2*1.000001).log() - 3*arb(0.235*d_lo).log() + wt.log()).upper())
    far = math.exp(max(far_log, -690.0))
    Er = pref*(tot + g1) + pref*far
    bands = {}
    for (v, lo, hi) in worst:
        key = '[0.4,40)' if lo < 40 else ('[40,800)' if lo < 800 else ('[800,998)' if lo < 998 else ('[998,1059)' if lo < 1059 else '[1059,1150]')))
        bands[key] = bands.get(key, 0.0) + v
    res = dict(La=La, Lb=Lb, F=F, Er_over_F=Er/F*(1 + 1e-9), small=[t1*pref/F, t2*pref/F], G1=g1*pref/F, far=pref*far/F, far_log=far_log,
               worst=worst[:8], bands=bands, Ksup=Ksup)
    json.dump(res, open(os.path.join(HERE, 'cache', f'resid3_{La}_{Lb}.json'), 'w'))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    La, Lb = float(sys.argv[1]), float(sys.argv[2])
    if sys.argv[3] == 'assemble':
        assemble(La, Lb, int(sys.argv[4]))
    else:
        run_part(La, Lb, int(sys.argv[3]), int(sys.argv[4]))

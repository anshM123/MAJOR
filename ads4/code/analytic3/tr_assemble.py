"""
tr_assemble.py La Lb: Route T''' with the signed transfer (ANALYTIC_SLIVER3 §9) for the L-cell [La, Lb].
Inputs: cache/jobs_La_Lb.jsonl + cache/results_La_Lb.jsonl (tr_jobs.py: transfer boxes 'res' and 'top2'), cache/contour3_La_Lb_k4.json
(termwise bottom/left/right of psi2 and all psi3 pieces), cache/loc3_La_Lb_1500.json (Loc''').
  * residual: cells with lo in [NU_RES): sup|kappa_J| <= max over d-slabs of |B| d_slab^{-5/2} (P M0 + Q M1) (transfer; M0, M1 by
    transfer_kernel.M_uniform), else (or if the box failed) the termwise kernel3.K_cell; bookkeeping as resid3 (G1, far, small nu);
  * psi2 top line: t >= TOP_T0 by the transfer, |kappa_J - kappa_inf| <= |B||d^{-alpha}|[(P M0 + Q M1) + (P + Q) + 1e-300]
    (kappa_inf truncated = the h = 0 face of the same box with X = 1; the j >= 1500 tail of kappa_inf is < e^{-3000}); M0 >= |S^(alpha)|,
    M1 >= |S^(alpha+1)| at Re alpha = 1/4 - 4 by rt_pp.PP + the REC chain (as contour3.piece_top); t < TOP_T0 termwise (contour3).
"""
import sys, os, json, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import transfer_kernel as TK
import tr_jobs as TJ
import kernel3 as K3
import contour3 as C3
os.environ['RT2_FINE'] = '1'
import rt2_resid as R2
import rt_cert as RC
import rt_pp as PPm
from flint import arb, acb


def load(La, Lb):
    jobs = [json.loads(l) for l in open(os.path.join(HERE, 'cache', f'jobs_{La}_{Lb}.jsonl')) if l.strip()]
    res = TJ.load_results(os.path.join(HERE, 'cache', f'results_{La}_{Lb}.jsonl'))
    return jobs, res


def box_for(jobs, res, kind, slab, lo, hi):
    for j in jobs:
        if j['kind'] == kind and j['slab'] == slab and j['f'][0] <= lo + 1e-9 and hi <= j['f'][1] + 1e-9:
            r = res.get(j['id'])
            if r and r['status'] == 'PASS':
                return j, r['P'], r['Q']
            return j, None, None
    return None, None, None


def residual(La, Lb, jobs, res, nslab):
    d_hi = math.exp(-La)
    sl = TJ.slabs(La, Lb, nslab)
    cells = R2.cell_list()
    rows = []; t0 = time.time(); nfall = 0
    for idx, (kind, lo, hi, m, w) in enumerate(cells):
        nsub = max(1, int(math.ceil((hi - lo)/w - 1e-9)))
        best = 0.0; how = 'termwise'
        if kind == 'cell' and TJ.NU_RES[0] <= lo and hi <= TJ.NU_RES[1] + 0.2:
            vals = []
            for si, (dlo, dhi) in enumerate(sl):
                j, P, Q = box_for(jobs, res, 'res', si, lo, hi)
                if P is None:
                    vals = None; break
                k, M0, M1 = TK.K_cell_transfer(lo, hi, dlo, dhi, P, Q)      # |B| d_slab^{1/2} (P M0 + Q M1)
                vals.append(k/d_hi**3)                                      # resid3 normalisation (caller multiplies by d_hi^3)
            if vals is not None:
                best = max(vals); how = 'transfer'
        if how == 'termwise':
            if kind == 'cell' and TJ.NU_RES[0] <= lo < TJ.NU_RES[1]: nfall += 1
            subs = [(lo + (hi - lo)*i/nsub, lo + (hi - lo)*(i + 1)/nsub) for i in range(nsub)]
            for (a_, b_) in subs:
                v, _ = K3.K_cell(a_, b_, La, Lb, use_frames=(kind != 'tail'), jmax=110, J0x=1500)
                best = max(best, v)
        rows.append(dict(idx=idx, kind=kind, lo=lo, hi=hi, m=m, K=best, how=how))
    print(f'  residual cells done ({nfall} transfer fallbacks) [{time.time()-t0:.0f}s]', flush=True)
    F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    pref = d_hi**3*(2/RC.LAMf)
    t1, t2 = R2.small_nu(La)
    tot = t1 + t2; Ksup = 0.0; bands = {}; worst = []
    for r in rows:
        v = 2*r['K']*r['m']; tot += v
        if r['lo'] < RC.LAMf: Ksup = max(Ksup, r['K'])
        key = r['how'] + (' <500' if r['lo'] < 500 else (' [500,998)' if r['lo'] < 998 else ' >=998'))
        bands[key] = bands.get(key, 0.0) + v*pref/F
        worst.append((v*pref/F, r['lo'], r['how']))
    worst.sort(reverse=True)
    g1 = 2*Ksup*1.66e-23
    wt = R2.wtail(1150.0)
    d_lo = math.exp(-Lb)
    far_log = float((arb(1.6*41*math.pi/2*1.000001).log() - 3*arb(0.235*d_lo).log() + wt.log()).upper())
    far = math.exp(max(far_log, -690.0))
    Er = pref*(tot + g1) + pref*far
    return dict(Er_over_F=Er/F*(1 + 1e-9), small=[t1*pref/F, t2*pref/F], G1=g1*pref/F, far=pref*far/F, bands=bands, worst=worst[:8],
                fallbacks=nfall)


def top_psi2(La, Lb, jobs, res, nslab, k0=4):
    """psi2 top line (dressed p = 3, Y = 6.25): termwise for t < TOP_T0, transfer above; returns the raw sum (before fac)."""
    d = math.exp(-La); d0 = d*1.0001
    p = 3; e0 = 2.5; Y = e0 - 0.25 + k0
    dt = d*d/RC.TAU_COEF
    sl = TJ.slabs(La, Lb, nslab)
    tot_tw = 0.0; tot_tr = 0.0; nfall = 0
    t = RC.DELTA
    while t < RC.LAMf:
        step = 0.05 if t < 3 else (0.25 if t < 30 else 1.0)
        t1 = min(t + step, RC.LAMup)
        box = acb(arb((t + t1)/2, (t1 - t)/2), arb(Y))
        w = RC.psi_abs(RC.Q2, acb((t + t1)/2, Y)/RC.LAM, (t1 - t)/2/RC.LAMf + 1e-12)*RC.B_abs_cell(complex((t + t1)/2, Y), (t1 - t)/2*1.0001)
        b0 = PPm.PP(0.25, t, t1, d0); b1 = PPm.PP(1.25, t, t1, d0)
        ch = C3.chain_abs(b0, b1, 0.25, t1, k0 + 162, d)
        xs = ch[k0:]                     # xs[0] >= |S^(alpha+1)|, xs[1] >= |S^(alpha)|
        done = False
        if t >= TJ.TOP_T0:
            vals = []
            for si in range(nslab):
                j, P, Q = box_for(jobs, res, 'top2', si, t, t1)
                if P is None:
                    vals = None; break
                vals.append(P*xs[1] + Q*xs[0] + P + Q + 1e-300)
            if vals is not None:
                tot_tr += w*max(vals)*d**(k0 - 0.25)*(t1 - t); done = True
            else:
                nfall += 1
        if not done:
            js = C3.jsum_general(box, xs, d, dt, 0.25 - k0, 160, 1500, abs(complex(0.25 - k0, t1)))
            tot_tw += w*js*(t1 - t)
        t = t1
    return tot_tw, tot_tr, nfall


def main(La, Lb, nslab=4):
    t0 = time.time()
    jobs, res = load(La, Lb)
    F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    d = math.exp(-La)
    ct = json.load(open(os.path.join(HERE, 'cache', f'contour3_{La}_{Lb}_k4.json')))
    lc = json.load(open(os.path.join(HERE, 'cache', f'loc3_{La}_{Lb}_1500.json')))
    tw, tr, nf = top_psi2(La, Lb, jobs, res, nslab)
    fac = d**3*(2/RC.LAMf)*2
    top = fac*(tw + tr)
    print(f'  psi2 top: termwise part {fac*tw/F:.3e}, transfer part {fac*tr/F:.3e} (F units), fallbacks {nf} [{time.time()-t0:.0f}s]', flush=True)
    rs = residual(La, Lb, jobs, res, nslab)
    psi2 = dict(ct['psi2']); psi2['top'] = top
    E2 = sum(v for k, v in psi2.items() if k != 'Y'); E3 = sum(v for k, v in ct['psi3'].items() if k != 'Y')
    Er = rs['Er_over_F']*F
    Erho = math.exp(max(lc['lE2'], -700)) + math.exp(max(lc['lE3'], -700))
    rho = (E2 + E3 + Er + Erho)/F*(1 + 1e-9)
    cpu = sum(r['secs'] for r in res.values())/3600
    lines = [f"Route T''' + signed transfer, L in [{La}, {Lb}] (beta_c in [{RC.beta_of_L(La):.3f}, {RC.beta_of_L(Lb):.3f}]), a in [998, 4500]: F_cert >= {F:.4e}",
             f"   psi2: " + '  '.join(f"{k} {psi2[k]/F:.2e}" for k in ('bottom', 'left', 'top', 'right')) + f"  (top: termwise t<{TJ.TOP_T0:.0f} {fac*tw/F:.2e}, transfer {fac*tr/F:.2e})",
             f"   psi3: " + '  '.join(f"{k} {ct['psi3'][k]/F:.2e}" for k in ('bottom', 'left', 'top', 'right')),
             f"   residual E_r/F = {rs['Er_over_F']:.4e} (small-nu {rs['small'][0]:.1e}+{rs['small'][1]:.1e}, G1 {rs['G1']:.1e}, far {rs['far']:.1e}); bands: "
             + ', '.join(f"{k} {v:.2e}" for k, v in sorted(rs['bands'].items())),
             f"      worst cells: " + ', '.join(f"[{w[1]}] {w[0]:.2e} ({w[2]})" for w in rs['worst'][:5]),
             f"   Loc''' remainders: log E_rho2 = {lc['lE2']:.1f}, log E_rho3 = {lc['lE3']:.1f}",
             f"   transfer jobs: {len(res)} run, {sum(1 for r in res.values() if r['status'] == 'PASS')} PASS, checker CPU {cpu:.2f} h",
             f"   rho <= {rho:.4e}   {'PASS' if rho < 1 else 'FAIL'}   [assembly {time.time()-t0:.0f}s]"]
    print('\n'.join(lines), flush=True)
    open(os.path.join(HERE, 'cache', f'cert3tr_{La}_{Lb}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 4)

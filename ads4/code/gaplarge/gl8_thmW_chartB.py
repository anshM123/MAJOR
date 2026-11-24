"""
gaplarge/gl8_thmW_chartB.py (GAP_LARGESPIN.md §4): Theorem W (Lemma P4W-chi form of AG_PROOF §5.2: premises P0, P1', W, AG_U, EDGE-R)
re-run at Delta_g = 2000 and 5000 on the CHART-B images (same beta, same m/Lam up to O(1/J)) of its 47 sample states.  [DV]

Delta_g as a parameter (no existing file is edited; GAP_SCALING §6.1 item 4): the Delta_g-dependent module constants of
  lsw/lsw_common (Dg, L, NU, X, PS, DN), p4win/p4w_common (L, D1000), p2con/pc_run (D1000), agproof/ag_flat (L), agproof/ag_chi (L)
are re-set IN THIS PROCESS by set_gap(Dg) (read-only imports; every function of those modules reads these globals at call time).
  * grid nu_k = (k + 1/2) 0.04 < Lam, PS = psi2(nu/Lam) (the Lam-independent NF02 weight), X = the k = 2 Abel weight re-solved with
    deshift/d_abel.Abel (dnu = 0.04, nu_max = Lam + 40), D(Dg) = PW3's monopole defect computed from that X (check: D(1000) = -9.98332e-4);
  * ag_exact is not needed here (it enters only the exact-scalar edge measurements, quoted).
Chart B: J' + 1/2 = (J + 1/2) Lam'/Lam0 rounded to the nearest J' of the SAME parity, m' = (J' + 1/2)/sinh(beta/2) (beta kept exact),
m'/Lam' = m/Lam0 (1 + O(1/J)).  The computed O(1/m^2) corrections (P2_CONTRACT: block side pc_run + kernel O(r) pc_ksh_run) are, per
unit leading image and in c2 units, functions of (beta, nu) ONLY (their fits x depend on q = e^{beta/2} alone) [D: code inspection,
pc_run.run_state / pc_ksh_run.main take only beta and the parity of J]; they are therefore reused at the image state, and only their
COSTS (psi2(nu/Lam'), D(Lam'), 1/m'^2) and the floors H' (ls_core.hyb_fast at Dg') are recomputed.
Per state: d_chi/H (everything but the edge integral E_chi, Lam_AG = 1, plain / +- transferred), the tolerated |E_chi| in units of S_E
(plain / transferred), the transfer factor (1 - x1^2)^-2, and the Lemma P4W ('local' AG') tolerances of p4w_states.analyse.
Env GL8_TEXT=1: ag_chi's rectangle-height grid T/T0 in {0.5..1} is extended by {1.15, 1.3, 1.6} (premise AG_U then needed up to 1.6 T0).
Usage: python3 gl8_thmW_chartB.py Dg [J1,J2,...]   -> cache/gl8_W_<Dg>.{txt,json}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for sub in ('deshift', 'n2sub', 'lsw', 'p2con', 'p4win', 'agproof'):
    sys.path.insert(0, os.path.join(CODE, sub))
import numpy as np
import d_abel as DA                 # read-only
import n2_common as N               # read-only
import p4w_common as C              # read-only (globals re-set)
import p4w_contour as K             # read-only
import p4w_states as PS             # read-only
import ag_flat as AF                # read-only (L re-set)
import ag_chi as AGC                # read-only (L re-set)
PJ = PS.PJ; Rn = PJ.Rn; W = C.W
assert Rn.W is W

OUT = []
def pr(s=''):
    print(s, flush=True); OUT.append(s)

L0 = N.lam(1000)
cfac = Rn.cfac


def set_gap(Dg):
    L = N.lam(Dg)
    A = DA.Abel(Dg, 0.04, L + 40.0)
    Y, res, it = A.solve()
    X = A.a*Y
    sel = A.nu < L
    nu, X, ps = A.nu[sel], X[sel], A.phi[sel]
    D = (X*N.Aop(lambda v: 1/cfac(v), nu)).real.sum()*0.04/L - (ps/cfac(nu)).real.sum()*0.04/L
    W.Dg = Dg; W.L = L; W.NU = nu; W.X = X; W.PS = ps; W.DN = 0.04
    C.L = L; C.D1000 = D; Rn.D1000 = D; AF.L = L; AGC.L = L
    return L, D, f'Abel N = {A.N}, its {it}, residual {res:.1e}'


def H_of(fl, J):
    f = {h: v['floor'] for h, v in fl.items()}
    return f['++']*f['+-']/(f['++'] + f['+-']) if J % 2 == 0 else f['+-']


def image(s, L):
    J, beta = s['J'], s['beta']
    t = (J + 0.5)*L/L0 - 0.5
    cand = [j for j in (int(np.floor(t)) - 1, int(np.floor(t)), int(np.floor(t)) + 1, int(np.floor(t)) + 2) if j % 2 == J % 2]
    Jp = min(cand, key=lambda j: abs(j - t))
    m = (Jp + 0.5)/np.sinh(beta/2)
    return Jp, m, m/L


def main():
    Dg = int(sys.argv[1])
    if os.environ.get('GL8_TEXT'):
        AGC.GRID_T = (1.6, 1.3, 1.15) + tuple(AGC.GRID_T)   # taller rectangles allowed (top side e^{-25 T/T0}; see GAP_LARGESPIN §4)
    only = [int(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else None
    t0 = time.time()
    pr(__doc__.split('Usage')[0])
    L, D, info = set_gap(Dg)
    pr(f'Dg = {Dg}: Lam = {L:.4f}, D = {D:+.7e} (Lam^2 D = {L*L*D:+.3f}); {info}')
    runs = PJ.load_runs(); ksh = PJ.load_runs('ksh_*_*.json')
    states = json.load(open(os.path.join(CODE, 'p2con', 'cache', 'jstar_bh.json')))['states']
    states = sorted(states, key=lambda z: (z['J'] % 2, z['J'], z['mfl']))
    if only:
        states = [s for s in states if s['J'] in only]
    ref = {(r['J'], round(r['mfl'], 6)): r for r in json.load(open(os.path.join(CODE, 'agproof', 'cache', 'chi.json')))}
    pr(' J1000 m/L1000 -> J\'  m\'/L\'  beta | H\' (H\'/H x (L\'/L0)^2) | x1  (1-x1^2)^-2 | d_chi/H plain / tr | eps_E tol/S_E plain / tr '
       '| Lam_AG tol (P4W local) plain / tr | [Dg = 1000 reference: d_tr, tol_tr]')
    rows = []
    for s in states:
        key0 = (s['J'], round(s['mfl'], 6))
        if Dg == 1000:
            Jp, m, mfl = s['J'], s['m'], s['mfl']
        else:
            Jp, m, mfl = image(s, L)
        beta = runs[key0]['beta']
        fl = W.floors(Jp, mfl)
        H = H_of(fl, Jp)
        sp = dict(J=Jp, mfl=mfl, m=m, beta=float(2*np.arcsinh((Jp + 0.5)/m)), H=H)
        key = (Jp, round(mfl, 6))
        runs[key] = runs[key0]; ksh[key] = ksh[key0]
        b = AGC.best_for_state(sp, runs, ksh)
        SEp, _ = AGC.scale_SE(Jp, mfl, b['u'], 0.0)
        SEt, _ = AGC.scale_SE(Jp, mfl, b['u'], -2.0)
        tolp = (b['H'] - b['d_plain']*b['H'])/SEp; tolt = (b['H'] - b['d_tr']*b['H'])/SEt
        o = PS.analyse(sp, runs, ksh)
        x1 = 1/mfl; tf = (1 - x1*x1)**-2
        r0 = ref.get(key0, {})
        row = dict(J1000=s['J'], mfl1000=s['mfl'], J=Jp, mfl=mfl, beta=sp['beta'], H=H, Hratio=H/s['H']*(L/L0)**2, x1=x1, transfer=tf,
                   d_plain=b['d_plain'], d_tr=b['d_tr'], tolE_plain=tolp, tolE_tr=tolt, SEp=SEp, SEt=SEt, fixed=b['fixed'],
                   LamAG_pl=o['local']['LamAG_tol_plain'], LamAG_tr=o['local']['LamAG_tol'], K4w_tol=o['local']['K4w_tol'],
                   win=o['local']['win']/H, bd=o['local']['bd']/H)
        rows.append(row)
        pr(f" {s['J']:4d} {s['mfl']:.4f} -> {Jp:5d} {mfl:.4f} {sp['beta']:.4f} | {H:.3e} ({row['Hratio']:.3f}) | {x1:.4f} {tf:8.1f} | "
           f"{b['d_plain']:.3f} / {b['d_tr']:.3f} | {tolp:9.3g} / {tolt:9.3g} | {row['LamAG_pl']:8.3g} / {row['LamAG_tr']:8.3g} | "
           + (f"[{r0.get('d_tr', float('nan')):.3f}, {r0.get('tolE_tr', float('nan')):.3g}]" if r0 else '')
           + f"  ({time.time() - t0:.0f}s)")
        json.dump(rows, open(os.path.join(HERE, 'cache', 'gl8_W_%d%s.json' % (Dg, '_sub' if only else '')), 'w'), default=float)
    pr('summary Dg = %d: max d_chi/H plain %.3f, transferred %.3f; min eps_E tol/S_E plain %.3g, transferred %.3g; '
       'min Lam_AG tol (local) plain %.3g, transferred %.3g' % (Dg, max(r['d_plain'] for r in rows), max(r['d_tr'] for r in rows),
                                                              min(r['tolE_plain'] for r in rows), min(r['tolE_tr'] for r in rows),
                                                              min(r['LamAG_pl'] for r in rows), min(r['LamAG_tr'] for r in rows)))
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(HERE, 'cache', 'gl8_W_%d%s.txt' % (Dg, '_sub' if only else '')), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()

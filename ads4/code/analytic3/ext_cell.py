#!/usr/bin/env python3
"""
ext_cell.py: one L-cell end to end with the author's sliver_extension package (ANALYTIC_SLIVER3 §11), reusing its scripts and config
format unchanged: make_resid.py -> check_norm.py (coupled coordinate for the high bands, as resid_coupled_hi/p4_coupled of the
[2.9, 3.2] run) -> resid_budget.py -> ellipse_extend.py (p, b) -> the ledger comparison of assemble_full.py.
Nothing is written into the package; all outputs go to cache/ext/La_Lb/.  The package's scalar.cpp is compiled unchanged into
analytic3/bin/scalar (SHA256 identical to code/norming/scalar.cpp) and check_norm.HERE is pointed there.

  python3 ext_cell.py La Lb [--jobs 3] [--floor F] [--preflight]

Steps are resumable (a step is skipped when its output exists).  Records wall time per step and the scalar-evaluation count
(sum of norming nodes over all boxes and channels) in cache/ext/La_Lb/ledger.json.
"""
import os, sys, json, time, math, argparse, subprocess
from fractions import Fraction as F
from pathlib import Path
HERE = Path(__file__).resolve().parent
EXT = Path(os.environ.get('SLIVER_EXT', str(HERE / '..' / 'sliver_ext')))
if not (EXT / 'ellipse_extend.py').exists():
    EXT = Path('/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/ext/sliver_extension')
EXT = EXT.resolve()
PY = sys.executable
SHIM = HERE / 'ext_shim.py'


def sh(args, cwd, log=None):
    env = dict(os.environ, PYTHONPATH=str(EXT), EXT_DIR=str(EXT), BIN_DIR=str(HERE / 'bin'), PYTHONDONTWRITEBYTECODE='1')
    t0 = time.time()
    r = subprocess.run([PY, str(SHIM)] + args, cwd=cwd, env=env, capture_output=True, text=True)
    if log:
        open(log, 'a').write(r.stdout + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f'{args[0]} failed: {r.stderr[-1500:]}')
    return r.stdout, time.time() - t0


def floor_for(La, Lb):
    sys.path.insert(0, str(HERE / '..' / 'analytic'))
    import rt_cert as RC
    Fc = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    return F(int(Fc*1e8), 10**8)*F(9999, 10000)          # rounded down


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('La'); ap.add_argument('Lb'); ap.add_argument('--jobs', type=int, default=3)
    ap.add_argument('--floor'); ap.add_argument('--preflight', action='store_true'); ap.add_argument('--max-nodes', type=int, default=400000)
    A = ap.parse_args()
    La, Lb = A.La, A.Lb; tag = f'{float(La)}_{float(Lb)}'
    od = HERE / 'cache' / 'ext' / tag; od.mkdir(parents=True, exist_ok=True)
    log = od / 'run.log'
    led_fn = od / 'ledger.json'
    led = json.load(open(led_fn)) if led_fn.exists() else dict(L=[La, Lb], steps={})
    def save(): json.dump(led, open(led_fn, 'w'), indent=1)
    floor = F(A.floor) if A.floor else floor_for(float(La), float(Lb))
    led['floor'] = str(floor); save()
    # 1. configs (make_resid.py; then the coupled variant for the bands >= 700 and p = 4, exactly as the [2.9, 3.2] run)
    if not (od / f'resid_{tag}.json').exists():
        out, w = sh(['make_resid.py', La, Lb], od, log); led['steps']['make_resid'] = dict(wall=w, out=out[-2000:]); save()
    lowc = od / f'resid_low_{tag}.json'; hic = od / f'resid_coupled_hi_{tag}.json'; p4c = od / f'p4_coupled_{tag}.json'
    if not hic.exists():
        rs = json.load(open(od / f'resid_{tag}.json')); p4 = json.load(open(od / f'p4_{tag}.json'))
        d0, dd = F(rs[0]['d0']), F(rs[0]['dd'])
        sys.path.insert(0, str(EXT)); from contour_direct import ds
        root = F(math.isqrt(int((1 - d0 - dd)*10**40)), 10**20)              # <= sqrt(1 - d_hi)
        xcap = 1/(998**2*(1 + root)**2)                                       # >= h/d^2 for every state of the box
        xh = F(math.ceil(xcap/2*10**15) + 1, 10**15)
        json.dump(rs[:2], open(lowc, 'w'), indent=2)
        hi_ = []
        p4all = [dict(c, name=c['name'].replace('real_3', 'real_4'), channels=[[4, 0], [4, 4]]) for c in rs]          # p = 4 real rows on every band (for the psi3 top-line extension)
        for c in rs[2:] + p4all[2:]:
            c = dict(c); c.update(h0=ds(xh), dh=ds(xh), coupled=True); hi_.append(c)
        json.dump(hi_[:len(rs) - 2], open(hic, 'w'), indent=2); json.dump(hi_[len(rs) - 2:], open(p4c, 'w'), indent=2)
        json.dump(json.load(open(lowc)) + [dict(c) for c in p4all[:2]], open(lowc, 'w'), indent=2)
    if A.preflight:
        for fn in (lowc, hic, p4c):
            out, w = sh(['check_norm.py', str(fn), '--preflight'], od)
            for l in out.splitlines():
                r = json.loads(l); print(fn.name, r['name'], r['p'], r['b'], r['nodes'], [a['m'] for a in r['axes']], flush=True)
        return
    # 2. norming runs
    for fn in (lowc, hic, p4c):
        res = od / (fn.stem + '_results.json')
        if not res.exists():
            out, w = sh(['check_norm.py', str(fn), '--jobs', str(A.jobs), '--max-nodes', str(A.max_nodes), '--output', str(res)], od, log)
            rows = json.load(open(res))
            led['steps'][fn.stem] = dict(wall=w, nodes=sum(r.get('nodes', 0) for r in rows),
                                         status=[(r['name'], r['b'], r['status']) for r in rows if r['status'] not in ('ENCLOSED', 'PASS')]); save()
    combo = od / f'resid_combo_{tag}_results.json'; comb = od / f'combined_{tag}.json'
    if not comb.exists():
        r = json.load(open(od / (lowc.stem + '_results.json'))) + json.load(open(od / (hic.stem + '_results.json')))
        json.dump(r, open(combo, 'w'), indent=2)
        json.dump(r + json.load(open(od / (p4c.stem + '_results.json'))), open(comb, 'w'), indent=2)
    # 3. residual budget
    bud = od / f'resid_combo_{tag}_results_budget.json'
    if not bud.exists():
        out, w = sh(['resid_budget.py', str(combo), str(floor)], od, log); led['steps']['resid_budget'] = dict(wall=w); save()
    # 4. contour (ellipse extension), all (p, b)
    for p in (3, 4):
        for b in (0, 4):
            ct = od / f'ct_{tag}_{p}_{b}.json'
            if not ct.exists():
                out, w = sh(['ellipse_extend.py', str(comb), '10.25', str(ct), str(p), str(b)], od, log)
                led['steps'][f'ct_{p}_{b}'] = dict(wall=w); save()
    # 5. ledger (assemble_full.py's comparison for this cell)
    sys.path.insert(0, str(EXT))
    from check_lowfreq import exp_u
    from check_norm import parameters
    from base_bounds import hi
    rows = json.load(open(combo))
    d, dd, *_ = parameters(rows[0]['config']); dl, du = d - dd, d + dd
    assert dl <= 1/exp_u(F(Lb)) and du >= exp_u(-F(La))
    assert all(parameters(r['config'])[:2] == [d, dd] for r in rows)
    res = []
    for r in json.load(open(bud)):
        b = r['b']; assert F(r['dlo']) == dl and F(r['dhi']) == du
        contour = F(0); parts = {}
        for p in (3, 4):
            c = json.load(open(od / f'ct_{tag}_{p}_{b}.json'))[0]
            assert c['p'] == p and c['b'] == b and F(c['dlo']) == dl and F(c['dhi']) == du
            contour += F(c['total']); parts[p] = {k: float(F(c[k])) for k in ('bottom', 'left', 'top', 'right', 'extra')}
        E = F(r['total_exact']) + contour
        res.append(dict(b=b, total_exact=str(E), floor_exact=str(floor), residual_total=float(F(r['total_exact'])), contour=float(contour),
                        contour_parts=parts, error_ratio_upper=float(hi(E/floor)), status='PASS_CONDITIONAL' if E < floor else 'FAIL'))
        print(json.dumps(res[-1]), flush=True)
    led['result'] = res
    led['scalar_evaluations'] = sum(v.get('nodes', 0) for v in led['steps'].values())
    led['wall_total'] = sum(v.get('wall', 0) for v in led['steps'].values())
    save()
    print(json.dumps(dict(cell=tag, scalar_evaluations=led['scalar_evaluations'], wall_total_s=led['wall_total'], jobs=A.jobs)))


if __name__ == '__main__':
    main()

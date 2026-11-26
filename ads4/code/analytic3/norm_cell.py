#!/usr/bin/env python3
"""
norm_cell.py: one L-cell of Route T''' with the author's norming package (code/norming/, ANGULAR_ROUTE §109; imported read-only,
its scalar.cpp compiled unchanged into analytic3/bin/scalar).  Generalises norming/assemble.py (which is hard-wired to [3.55, 3.65])
to any cell [La, Lb] and adds the complete top line.  ANALYTIC_SLIVER3 §10.

  python3 norm_cell.py inputs   La Lb            write cache/norm/La_Lb/{resid,top}.json (check_norm configs)
  python3 norm_cell.py run      La Lb [--jobs 3] run check_norm on both inputs (resumable per box: skips boxes already in the JSONL)
  python3 norm_cell.py assemble La Lb            total error / F_lower per channel b in {0, 4}; PASS iff < 1

Pieces (per channel b):
  low frequency [0, 160]: 1.46e-7 (author's theorem, check_lowfreq.py, uniform on 0.8 <= L <= 3.65);
  residual [159, 1040.1]: 16 EPSW sub-bands (norming/band_masses.csv) x enclosing real norming box over the whole d-range of the cell,
      full twist; kernel 2 sqrt(d_hi)|B| (P M0 + Q M1), M0, M1 = base_bounds.base_bound on the sub-band (as assemble.py);
  cutoff: 2/998 * max kernel * 1.34e-14;  [1040, 1150]: recurrence majorant x Gaussian envelope;  > 1150: unsplit full kernel;
  localisation: < e^{-46} (Loc-Joint, §104-105; Loc''' gives the same independently for L >= 2);
  contour: psi2 and psi3 top line tiled completely by complex norming boxes (assemble_top.bound_segment, x in [1/4, 999]);
      bottom/left/right of psi2 and psi3 from contour3.py (termwise, Lemma NB/V and Lemma F/B; valid where its asserts pass, L >= 3.05).
  F_lower: LRINF_TIGHT F_cert minimum over the beta-range of the cell (rt_cert.fcert_min), rounded down to 4 digits.
"""
import sys, os, json, math, argparse, subprocess, time
from fractions import Fraction as F
from pathlib import Path
HERE = Path(__file__).resolve().parent
NORM = (HERE / '..' / 'norming').resolve()
sys.path.insert(0, str(NORM)); sys.path.insert(1, str(HERE))
import check_norm as CN
CN.HERE = HERE / 'bin'                       # compiled copy of norming/scalar.cpp (identical source, SHA256 checked)
from base_bounds import base_bound, sqrtup, hi, hcap as HCAP
from check_lowfreq import exp_u
import assemble_top as AT
import assemble as AS
OUT = HERE / 'cache' / 'norm'

RES_BANDS = [(159, 400.1), (400, 700.1), (700, 850.1), (850, 950.1), (950, 1040.1)]
RHOS = ["2", "4", "8", "16", "32", "64", "128", "512", "2048", "8192", "32768", "131072", "1048576"]


def _dstr(x, k):
    """exact decimal string of a Fraction with denominator dividing 10^k."""
    v = x*10**k; assert v.denominator == 1
    v = int(v); return f'{v//10**k}.{v % 10**k:0{k}d}'


def dbox(La, Lb):
    """decimal d-box [d0 - dd, d0 + dd] containing [e^{-Lb}, e^{-La}] (7 decimals), and the full-twist h (3 significant digits, up)."""
    lo = 1/exp_u(F(Lb)); hi_ = exp_u(F(-La))        # exp_u >= exp, so lo <= e^{-Lb}, hi_ >= e^{-La}
    s = 10**7
    dl = F(math.floor(lo*s) - 1, s); dh = F(math.ceil(hi_*s) + 1, s)
    if (dh - dl)*s % 2: dh += F(1, s)
    d0 = (dl + dh)/2; dd = (dh - dl)/2
    req = (d0 + dd)**2/(998**2*(1 + F(math.isqrt(int((1 - d0 - dd)*10**30)), 10**15))**2)
    e = -int(math.floor(math.log10(float(req/2)))) + 2
    h = F(math.ceil(req/2*10**e) + 1, 10**e)
    box = dict(d0=_dstr(d0, 7), dd=_dstr(dd, 7), h0=_dstr(h, e), dh=_dstr(h, e))
    assert F(box['d0']) - F(box['dd']) == dl and F(box['d0']) + F(box['dd']) == dh and 2*F(box['h0']) >= req
    return box, dl, dh


def inputs(La, Lb, top_width=(20, 10, 5)):
    box, dl, dh = dbox(La, Lb)
    od = OUT / f'{La}_{Lb}'; od.mkdir(parents=True, exist_ok=True)
    res = []
    for (lo, hi_) in RES_BANDS:
        c = dict(name=f'res_{lo}-{hi_}', **box, nu0=f'{(lo + hi_)/2:.2f}', dnu=f'{(hi_ - lo)/2:.2f}', fullTwist=True, K=120, digits=65,
                 epsilon='0.001', rhos=RHOS[:9], channels=[[3, 0], [3, 4]])
        CN.cap_check(c); res.append(c)
    top = []
    for (p, Y) in ((3, '6.25'), (4, '7.25')):
        x = F(1, 4)
        while x < 999:
            w = 1 if x < 20 else (top_width[0] if x < 400 else (top_width[1] if x < 900 else top_width[2]))
            x1 = min(F(999), x + w)
            c = dict(name=f'top{p}_{float(x):.2f}', **box, nu0=f'{float((x + x1)/2):.3f}', dnu=f'{float((x1 - x)/2):.3f}', fullTwist=True,
                     digits=65, K=90, channels=[[p, 0], [p, 4]], epsilon='0.00001', rhos=RHOS, nuImag=Y)
            CN.cap_check(c); top.append(c); x = x1
    import norm_extra as NE
    right = NE.right_inputs(box, F('6.25'), 3) + NE.right_inputs(box, F('7.25'), 4)
    for c in right: CN.cap_check(c)
    json.dump(res, open(od / 'resid.json', 'w'), indent=1); json.dump(top, open(od / 'top.json', 'w'), indent=1)
    json.dump(right, open(od / 'right.json', 'w'), indent=1)
    print(f'{len(res)} residual boxes, {len(top)} top boxes, {len(right)} right-edge boxes -> {od}')


def run(La, Lb, jobs=3, which=('resid', 'top', 'right'), max_nodes=400000):
    """run check_norm box by box (resumable: results appended to {kind}_results.jsonl, finished box names skipped)."""
    od = OUT / f'{La}_{Lb}'
    for kind in which:
        cfg = json.load(open(od / f'{kind}.json'))
        rfn = od / f'{kind}_results.jsonl'
        done = set()
        if rfn.exists():
            for l in open(rfn):
                if l.strip():
                    r = json.loads(l); done.add((r['name'], r['p'], r['b']))
        for c in cfg:
            if all((c['name'], p, b) in done for p, b in c['channels']): continue
            tmp = od / f'_box.json'; json.dump([c], open(tmp, 'w'))
            t0 = time.time()
            r = subprocess.run([sys.executable, str(HERE / 'norm_run.py'), str(tmp), '--jobs', str(jobs), '--max-nodes', str(max_nodes),
                                '--output', str(od / '_out.json')], capture_output=True, text=True)
            if r.returncode != 0:
                print(c['name'], 'ERROR', r.stderr[-300:], flush=True); continue
            rows = json.load(open(od / '_out.json'))
            with open(rfn, 'a') as fo:
                for row in rows:
                    row['wall'] = time.time() - t0; fo.write(json.dumps(row) + '\n')
            print(c['name'], [(row['b'], row['status'], row.get('nodes'), round(row.get('P_approx', -1), 2), round(row.get('Q_approx', -1), 2)) for row in rows],
                  f'{time.time()-t0:.0f}s', flush=True)


def load(od, kind):
    rows = []
    fn = od / f'{kind}_results.jsonl'
    if fn.exists():
        rows = [json.loads(l) for l in open(fn) if l.strip()]
    return rows


def assemble(La, Lb):
    od = OUT / f'{La}_{Lb}'
    box, DLO, DHI = dbox(La, Lb)
    sys.path.insert(0, str(HERE / '..' / 'analytic'))
    import rt_cert as RC
    Fc = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
    floor = F(int(Fc*1e10), 10**10)*F(9999, 10000)
    AS.DLO, AS.DHI = DLO, DHI                       # assemble.py's kernel normalisation uses its module constants
    res = load(od, 'resid'); top = load(od, 'top')
    import norm_extra as NE
    nr_fn = HERE / 'cache' / f'nearreal3_{La}_{Lb}_k4.json'; c3_fn = HERE / 'cache' / f'contour3_{La}_{Lb}_k4.json'
    ct = json.load(open(nr_fn)) if nr_fn.exists() else json.load(open(c3_fn))
    right_rows = load(od, 'right')
    import csv
    masses = list(csv.DictReader(open(NORM / 'band_masses.csv')))
    out = []
    for b in (0, 4):
        residual = F(0); maxk = F(0)
        for r in masses:
            l, u, m = map(F, [r['lo'], r['hi'], r['mass']])
            cand = [s for s in res if s['b'] == b and s['p'] == 3 and s['status'] in ('PASS', 'ENCLOSED')
                    and F(s['config']['nu0']) - F(s['config']['dnu']) <= l and F(s['config']['nu0']) + F(s['config']['dnu']) >= u]
            assert cand, ('missing residual band', b, str(l), str(u))
            K, M0, M1 = min((AS.kernel_from_result(s, l, u) for s in cand), key=lambda t: t[0])
            residual += F(2, 998)*K*m; maxk = max(maxk, K)
        cutoff = F(2, 998)*maxk*F('1.34e-14')
        gauss = F(0)
        for l, u in zip([1040, 1050, 1070, 1100], [1050, 1070, 1100, 1150]):
            c = dict(**box, nu0=str(float((l + u)/2)), dnu=str(float((u - l)/2)), fullTwist=True); CN.cap_check(c)
            ma, _, _ = CN.major(c, 3, b, [F(1)]*3)
            P = F(ma, CN.S); Q = P/DLO; M0, M1 = base_bound(3, b, DLO, DHI, F(l), F(u))
            gauss += F(2, 998)*AS.normalized(P, Q, M0, M1, F(l))*AS.gaussian_mass(l)
        far = F(2, 998)*200*AS.gaussian_mass(1150, F(11, 7)); assert far < F('1e-100')
        low = F('1.46e-7'); loc = F('1.1e-20')
        # top line, psi2 (p = 3) and psi3 (p = 4), complete tiling of x in [1/4, 999]
        tops = {}
        for (p, Y) in ((3, F('6.25')), (4, F('7.25'))):
            cands = []
            for r in top:
                c = r['config']; d, dd, h, hh, nu, dnu, im = CN.parameters(c)
                if r['p'] == p and r['b'] == b and r['status'] in ('PASS', 'ENCLOSED') and d - dd <= DLO and d + dd >= DHI and im == Y:
                    cands.append((nu - dnu, nu + dnu, r))
            pos = F(1, 4); tot = F(0)
            while pos < 999:
                el = [(l, u, r) for l, u, r in cands if l <= pos < u]
                if not el:
                    tot = None; break
                l, u, r = max(el, key=lambda z: z[1]); end = min(F(999), u)
                while pos < end:
                    q = min(end, pos + 5); tot += AT.bound_segment(r, p, b, DLO, DHI, pos, q, Y); pos = q
            tops[p] = tot
        if tops[3] is None or tops[4] is None:
            out.append(dict(b=b, status='INCOMPLETE_TOP')); continue
        near = sum(F(repr(ct[pp][k])) for pp in ('psi2', 'psi3') for k in ('bottom', 'left'))
        rights = {p: NE.right_bound(right_rows, p, b, DLO, DHI, Y) for (p, Y) in ((3, F('6.25')), (4, F('7.25')))}
        if any(v is None for v in rights.values()):
            if all('right' in ct[pp] for pp in ('psi2', 'psi3')):
                rights = {3: F(repr(ct['psi2']['right'])), 4: F(repr(ct['psi3']['right']))}       # termwise (contour3)
            else:
                out.append(dict(b=b, status='INCOMPLETE_RIGHT')); continue
        # real segment [LAM_R, Lam] (right edge placed at LAM_R): length x sup|psi2| x (|K_J| + |K_inf|), K_J from the [950, 1040.1] box
        seg = F(0)
        cand = [s_ for s_ in res if s_['b'] == b and s_['p'] == 3 and F(s_['config']['nu0']) + F(s_['config']['dnu']) >= 999]
        Kseg, _, _ = min((AS.kernel_from_result(s_, F('998.4998'), F(999)) for s_ in cand), key=lambda t: t[0])
        Kinf = 2*sqrtup(DHI)*F(1)/sqrtup(1 - DHI)            # 2 d^{1/2} |B H(d)| <= 2 d^{1/2} coth(pi nu)/(2 sqrt(1-d)) <= d^{1/2}/sqrt(1-d)
        seg = F(2, 998)*F(13, 100000)*AT.psi_sup(3, F('998.4998'), AT.LAMU, F(0))*(Kseg + Kinf)
        contour = tops[3] + tops[4] + near + rights[3] + rights[4] + seg
        total = low + residual + cutoff + gauss + far + loc + contour
        r = dict(b=b, low=float(low), residual=float(residual), cutoff=float(cutoff), gaussian=float(gauss), loc=float(loc),
                 top_psi2=float(tops[3]), top_psi3=float(tops[4]), right_psi2=float(rights[3]), right_psi3=float(rights[4]), near_real=float(near), segment=float(seg),
                 contour_other=float(contour - tops[3] - tops[4]), total=float(total), floor=float(floor), ratio=float(total/floor),
                 status='PASS' if total < floor else 'FAIL')
        out.append(r); print(json.dumps(r), flush=True)
    wall = sum(r.get('wall', 0) for r in res + top)/2   # wall recorded once per box for both channel rows
    summ = dict(La=La, Lb=Lb, d=[str(DLO), str(DHI)], rows=out, norming_wall_s=wall, boxes=len(res) + len(top))
    json.dump(summ, open(od / 'assembly.json', 'w'), indent=1)
    return summ


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('La', type=float); ap.add_argument('Lb', type=float)
    ap.add_argument('--jobs', type=int, default=3); ap.add_argument('--which', default='resid,top,right')
    A = ap.parse_args()
    if A.cmd == 'inputs': inputs(A.La, A.Lb)
    elif A.cmd == 'run': run(A.La, A.Lb, A.jobs, tuple(A.which.split(',')))
    elif A.cmd == 'assemble': assemble(A.La, A.Lb)

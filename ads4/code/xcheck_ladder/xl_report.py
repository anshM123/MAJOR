"""
xcheck_ladder/xl_report.py: merge out/xcheck_<Dg>.csv into out/results.csv with the comparison columns, and print the summary.
  rounding  = half a unit of the last printed digit of the direct value
  diff      = exact - direct;  diff/budget;  diff/stat (jackknife / Monte Carlo error, where stated)
  d_exact   = exact - hybrid,  d_direct = direct - hybrid  (beyond-leading-Regge correction, both in units of N_s)
"""
import csv, os, sys
sys.dont_write_bytecode = True
import xl_common as X
import xl_states as ST


def ndig(v):
    s = repr(v)
    return len(s.split('.')[1]) if '.' in s else 0


def main():
    rows = []
    for Dg in (1000, 2000, 5000):
        fn = os.path.join(X.OUT, f'xcheck_{Dg}.csv')
        if os.path.exists(fn):
            rows += list(csv.DictReader(open(fn)))
    out = []
    for r in rows:
        ex, di, hy = float(r['exact']), float(r['direct']), float(r['hyb'])
        bud = float(r['budget']) if r['budget'] not in ('', 'None') else float('nan')
        st = float(r['stat']) if r['stat'] not in ('', 'None') else float('nan')
        rnd = 0.5*10**(-ndig(di))
        o = dict(Dg=r['Dg'], J=r['J'], mass=r['mass'], m_over_Lam=r['m_over_Lam'], block=r['block'], exact_channel=r['exact_b'],
                 exact=f'{ex:.6f}', exact_pp=r['ratio_pp'][:9], exact_pm=r['ratio_pm'][:9], hybrid=f'{hy:.6f}',
                 direct=r['direct'], budget=f'{bud:.4g}', stat=r['stat'] if r['stat'] not in ('', 'None') else '',
                 rounding=f'{rnd:.1g}', diff=f'{ex-di:+.5f}', diff_over_budget=f'{(ex-di)/bud:+.3f}',
                 diff_over_stat=f'{(ex-di)/st:+.2f}' if st == st else '',
                 within_rounding='yes' if abs(ex-di) <= rnd else 'no', d_exact=f'{ex-hy:+.5f}', d_direct=f'{di-hy:+.5f}',
                 exact_numerical_spread=r['num_spread'], nu_gt_Lam_contribution=r['tail_nu_gt_Lam'],
                 exact_cert=r.get('exact_cert', ''), exact_quadGL24=r.get('exact_q24', ''), exact_grid=r.get('exact_grid', ''),
                 note=r['note'])
        out.append(o)
    with open(os.path.join(X.OUT, 'results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        for o in out:
            w.writerow(o)
    # summary
    print(f'{len(out)} comparisons')
    print('max |diff|/budget:', max(abs(float(o['diff_over_budget'])) for o in out))
    neg = [o for o in out if float(o['exact']) <= 0]
    print('exact values <= 0:', len(neg), ' min exact:', min(float(o['exact']) for o in out))
    for Dg in ('1000', '2000', '5000'):
        sub = [o for o in out if o['Dg'] == Dg]
        if not sub:
            continue
        big = [o for o in sub if o['diff_over_stat'] and abs(float(o['diff_over_stat'])) > 2]
        print(f'Dg={Dg}: |diff| > 2 stat (where a stat error is quoted):',
              ', '.join(f"J={o['J']} {o['mass']} {o['block']} {o['diff']} ({o['diff_over_stat']} stat)" for o in big))
        top = sorted(sub, key=lambda o: -abs(float(o['diff_over_budget'])))[:5]
        print(f'   largest |diff|/budget:', ', '.join(f"J={o['J']} {o['mass']} {o['block']} {o['diff_over_budget']}" for o in top))
        print(f'   max |diff| for J<=13: {max(abs(float(o["diff"])) for o in sub if int(o["J"]) <= 13):.5f};'
              f' J>=14: {max(abs(float(o["diff"])) for o in sub if int(o["J"]) >= 14):.5f}' if any(int(o['J'])>=14 for o in sub) else '')
    print('artefacts (exact at the same state):')
    for Dg, J, mass, ch, v, what in ST.ARTEFACTS:
        e = next(o for o in out if o['Dg'] == str(Dg) and o['J'] == str(J) and o['mass'] == mass and o['block'] == ch)
        print(f'  Dg={Dg} J={J} {mass}: chain/variant {v:+.3f} ({what})  vs exact {float(e["exact"]):+.5f}, direct {e["direct"]}')


if __name__ == '__main__':
    main()

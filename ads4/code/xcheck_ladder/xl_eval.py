"""
xcheck_ladder/xl_eval.py: evaluate the exact finite-mass action (Eq. fp-H) at every directly computed state of xl_states.py,
normalise by the state's N_s (xl_norm.py, the convention of the direct rows), and compare with the direct value.

Normalisation map (see README.md):  lambda_min/N_s  =  min_b  H_ang,b(Delta,J) / N_{s,b}(m),   b over the structures of the block:
   even J 2x2 (++,+-): b = 0 (++) and b = 4 (+-);  odd J (parity-odd, 1x1 '+-'): b = 4;  OO (parity-odd even J): b = 0 with N_s(++).
The factor B(Delta,J)|a_{O,b,0}|^2 of Eq. (fp-target) is common to M_Omega and to the flat normalisation N^phys_s (Eq. fp-lift), so it
cancels in the ratio.

Usage: python3 xl_eval.py [Dg ...]     -> out/xcheck_<Dg>.csv, printed table
"""
import sys, os, json, math, time, csv
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp
from flint import arb, ctx
import xl_common as X
import xl_norm as N
import xl_states as ST

PREC = 1400


def load_pmom(tag, K):
    d = json.load(open(os.path.join(X.OUT, f'pmom_{tag}.json')))
    ctx.prec = PREC
    Sr = [[arb(s) for s in S] for S in d['S']]
    tot = [sum((S[j] for S in Sr), arb(0)) for j in range(len(Sr[0]))]
    inside = Sr[0]                                       # nu <= Lam
    Zt, Yt = X.moments_from_power(tot, K)
    Zi, Yi = X.moments_from_power(inside, K)
    return dict(Z=Zt, Y=Yt, Zin=Zi, Yin=Yi, ranges=d['ranges'])


def T_and_mfl(Dg, J, mass):
    lam = X.lam_int(Dg); C = J*(J + 1)
    if mass == 'tw':
        D = J + Dg; T = Fr((D - 1)*(D - 2)); label = 'tw'
    elif mass == 'img':
        T = Fr(lam*(998 + 2*J + 1), 999) + C; label = 'img'
    else:
        T = Fr(mass)**2*lam + C; label = mass
    with mp.workdps(40):
        mfl = mp.sqrt(mp.mpf(T.numerator)/T.denominator - C)/mp.sqrt(lam)
        Delta = (3 + mp.sqrt(1 + 4*mp.mpf(T.numerator)/T.denominator))/2
        return T, mp.nstr(mfl, 30), float(Delta), label


def main():
    Dgs = [int(a) for a in sys.argv[1:]] or [1000, 2000, 5000]
    K = 40
    for Dg in Dgs:
        lam = X.lam_int(Dg)
        ctx.prec = PREC
        Q = X.Q_exact(K, lam)
        src = {}
        if Dg == 1000:
            src['cert'] = dict(zip(('Z', 'Y'), X.load_cert_moments(1000, 1000, kmax=K)))
            src['q24'] = load_pmom('cont_1000_n24_p1400', K)
            src['q32'] = load_pmom('cont_1000_n32_p1400', K)
            src['q24p900'] = load_pmom('cont_1000_n24_p900', K)
            src['grid'] = load_pmom('grid_1000_p1400', K)
            main_src = 'cert'
        else:
            src['grid'] = load_pmom(f'grid_{Dg}_p1400', K)
            main_src = 'grid'
        hf = None
        if Dg == 1000:
            sys.path.insert(0, os.path.join(X.ADS4, 'code', 'largespin'))
            import ls_core as hf
        rows = []
        for s in [r for r in ST.S if r['Dg'] == Dg]:
            J = s['J']
            T, mfl, Delta, label = T_and_mfl(Dg, J, s['mass'])
            Ta = X.to_arb(T)
            bs = {'min': [(0, '++'), (4, '+-')], '+-': [(4, '+-')], 'OO': [(0, '++')]}[s['ch']]
            per_b = []
            for b, h in bs:
                Hs = {}
                for nm, M in src.items():
                    Hs[nm] = X.H_ang(J, b, Ta, M['Z'], M['Y'], Q, lam)
                h2tot, _ = X.H_ang(J, b, Ta, src[main_src if main_src != 'cert' else 'q24']['Z'],
                                   src[main_src if main_src != 'cert' else 'q24']['Y'], Q, lam, parts=True)
                Mi = src['q24'] if Dg == 1000 else src['grid']
                h2in, _ = X.H_ang(J, b, Ta, Mi['Zin'], Mi['Yin'], Q, lam, parts=True)
                tail = h2tot - h2in                          # contribution of nu > Lam
                hy = N.hybrid(Dg, J, mfl, h)
                Ns = float(hy['Ns'])
                Nsf = float('nan')
                if hf is not None:
                    Nsf = hf.hyb_fast(Dg, J, float(mfl), h)['Ns']
                Hm = float(Hs[main_src].mid())
                spread = max(abs(float((v - Hs[main_src]).mid())) for v in Hs.values()
                             if v is not Hs[main_src]) if len(Hs) > 1 else float('nan')
                rec = dict(b=b, h=h, H=Hm, Ns=Ns, Ns_hybfast=Nsf, Hhyb=float(hy['Hhyb']), ratio=Hm/Ns,
                           tail_over_Ns=float(tail.mid())/Ns, spread_over_Ns=spread/Ns)
                for nm, v in Hs.items():
                    rec['ratio_' + nm] = float(v.mid())/Ns
                per_b.append(rec)
            best = min(per_b, key=lambda r: r['ratio'])
            ex = best['ratio']
            diff = ex - s['val']
            row = dict(Dg=Dg, J=J, mass=label, m_over_Lam=mfl[:12], Delta=f'{Delta:.6f}', block=s['ch'],
                       exact=ex, exact_b=best['h'], direct=s['val'], budget=s['budget'], L=s['L'], stat=s['stat'],
                       diff=diff, diff_over_budget=diff/s['budget'] if s['budget'] else float('nan'),
                       diff_over_stat=diff/s['stat'] if s['stat'] else float('nan'),
                       hyb=min(r['Hhyb']/r['Ns'] for r in per_b),
                       ratio_pp=next((r['ratio'] for r in per_b if r['b'] == 0), float('nan')),
                       ratio_pm=next((r['ratio'] for r in per_b if r['b'] == 4), float('nan')),
                       num_spread=max(r['spread_over_Ns'] for r in per_b),
                       tail_nu_gt_Lam=max(abs(r['tail_over_Ns']) for r in per_b),
                       Ns_rel_check=max(abs(r['Ns'] - r['Ns_hybfast'])/r['Ns'] for r in per_b) if hf else float('nan'),
                       note=s['note'])
            for nm in src:
                row['exact_' + nm] = min(r['ratio_' + nm] for r in per_b)
            rows.append(row)
            print(f"Dg={Dg} J={J:2d} {label:>6s} {s['ch']:>3s} m/L={mfl[:9]}  exact={ex:+.5f} ({best['h']})  direct={s['val']:+.5f}"
                  f"  diff={diff:+.5f}  diff/budget={row['diff_over_budget']:+.3f}  hyb={row['hyb']:+.5f}"
                  f"  spread={row['num_spread']:.1e} tail={row['tail_nu_gt_Lam']:.1e} Nschk={row['Ns_rel_check']:.1e}", flush=True)
        fn = os.path.join(X.OUT, f'xcheck_{Dg}.csv')
        with open(fn, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print('written', fn, flush=True)


if __name__ == '__main__':
    main()

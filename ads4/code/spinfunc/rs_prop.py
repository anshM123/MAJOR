"""
rs_prop (HRS_CHECK.md): propagation of the crystal-minus-engine node differences dG(nu) (rs_compare) into s_d and FF/N of the
binding state (++ / OO, J = 0 at the gap) through the Result Q chain (Stage 7 §2: s_d = 0.04 sum_nu chi_d(nu) Re Gamma_K(nu);
FF -> FF - a 2u^3 s_d, bg_qan.ffn, read-only).
  dG(nu) per harmonic (c2 units of Gamma/sum|R_j|), measured at the nodes; between nodes linear interpolation; outside the measured range
  NOT TESTED (set to 0, and its weight reported separately).  Gamma = sum|R_j| dG; the chain's a_h subtraction -(Gamma(nu_min)/NH(nu_min)) NH
  needs dGamma at nu = 0.02 (not measured; = 0 by the engine's normalisation; its projection 0.04 sum chi_d NH is printed).
  Errors: Gaussian resampling of the node values (independent nodes; the kappa-fit correlation is included by refitting per sample).
Also: per window, the uniform per-harmonic |dG| that would move s_d by the margins s_fail - s_d (delta = 0.2: 1.78 - 0.51; delta = 0.05:
0.58 - 0.51), and the 2 sigma bound on |dG| measured there.
Usage: python3 rs_prop.py Dg tags [est]      (state = J = 0 gap, Delta = Dg)   -> stdout + cache/rs_prop_D<Dg>.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import rs_compare as RC

SFAIL = {0.2: 1.78, 0.05: 0.58}


def main(Dg, tags, est='dir', J=0, exclude=()):
    import bf_qan as BQ
    import bg_qan as BG
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    rows = RC.main(J, float(Dg), tags)
    rows = [o for o in rows if not any(abs(o['nu'] - x) < 1e-6 for x in exclude)]
    db, r, _ = RC.engine(float(Dg), J)
    nus = r['nus']; R = BQ.Rj(nus); SR = sum(np.abs(v) for v in R.values()); NH = r['NH']
    W = BQ.weights()
    nn = np.array([o['nu'] for o in rows])
    order = np.argsort(nn); rows = [rows[i] for i in order]; nn = nn[order]
    infl = float(os.environ.get('RS_INFL', '1'))    # error inflation (independent-seed repeats at nu = 12.02: jackknife low by x2-3)
    d = np.array([o['d_' + est] for o in rows]); e = infl*np.array([o['e_' + est] for o in rows])
    g = np.array([o['m']**2*o['Tst']/o['SR'] for o in rows])
    def dG_nodes(dv):
        w = 1.0/np.maximum(np.abs(e)**2, 1e-30)
        ref = dv[int(np.argmax(w))]; dv = np.array([ref + RC.wrap(x - ref) for x in dv])
        k = np.sum(w*dv)/np.sum(w)
        return (g*(dv - k)).real
    def proj(dGn):
        prof = np.where((nus >= nn[0]) & (nus <= nn[-1]), np.interp(nus, nn, dGn), 0.0)
        out = {}
        for gname, (_, cv, cd) in W.items():
            k = min(len(nus), len(cd))
            out[gname] = 0.04*np.sum(np.asarray(cd[:k])*SR[:k]*prof[:k])
        return out
    dG0 = dG_nodes(d)
    P0 = proj(dG0)
    rng = np.random.default_rng(7)
    samp = []
    for _ in range(2000):
        dv = d + rng.normal(size=d.shape)*e.real + 1j*rng.normal(size=d.shape)*e.imag
        samp.append(proj(dG_nodes(dv)))
    pr(f'# rs_prop: Delta_g={Dg} J={J} est={est} tags={tags}; node errors x {infl}; nodes {[round(float(v), 2) for v in nn]} (measured range {nn[0]:.2f}-{nn[-1]:.2f}; outside: not tested)')
    m = r['m']; u = (np.sqrt((Dg - 2.0)*(Dg - 1.0))/m)**2
    pr('  0.04 sum chi_d NH (projection of the a_h class; ~0 means the unmeasured dGamma(0.02) does not matter): ' +
       ', '.join(f'{gn} {0.04*np.sum(np.asarray(cd[:min(len(nus), len(cd))])*NH[:min(len(nus), len(cd))]).real:+.2e}' for gn, (_, _, cd) in W.items()))
    for gname in W:
        sd_eng = 0.04*np.sum(np.asarray(W[gname][2][:len(nus)])*r['GK'][:len(W[gname][2])].real)
        s = np.array([x[gname] for x in samp])
        f = BG.ffn('++', J, u, int(Dg), sd_eng + P0[gname], 0.0)
        f0 = BG.ffn('++', J, u, int(Dg), sd_eng, 0.0)
        pr(f'  {gname}: engine s_d = {sd_eng:+.4f}; crystal-engine on the measured range: Delta s_d = {P0[gname]:+.4f} +- {s.std():.4f}'
           f'  -> FF/N (s_v = 0 here) d=0.2: {f0[0.2][1]:+.4f} -> {f[0.2][1]:+.4f} (+- {f[0.2][2]*s.std():.4f}); d=0.05: {f0[0.05][1]:+.4f} -> '
           f'{f[0.05][1]:+.4f} (+- {f[0.05][2]*s.std():.4f}); slope a2u^3/N = {f[0.2][2]:.4f} / {f[0.05][2]:.4f}')
    # window tolerances (H0)
    cd = np.asarray(W['H0'][2][:len(nus)])
    sd_eng = 0.04*np.sum(cd*r['GK'][:len(cd)].real)
    pr('  windows (H0): weight 0.04 sum|chi_d| sum|R_j|; uniform |dG| breaking delta=0.2 / 0.05 (margin s_fail - s_d); measured |dG| + 2 sigma there')
    for a, b in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 6), (6, 8), (8, 10), (10, 13), (13, 18), (18, 41)):
        sel = (nus[:len(cd)] >= a) & (nus[:len(cd)] < b)
        wgt = 0.04*np.sum(np.abs(cd[sel]*SR[:len(cd)][sel]))
        th = [(SFAIL[dl] - sd_eng)/wgt for dl in (0.2, 0.05)]
        inw = [(dG0[i], np.hypot(abs(g[i].real)*e[i].real, abs(g[i].imag)*e[i].imag)) for i in range(len(nn)) if a <= nn[i] < b]
        meas = ', '.join(f'{x:+.1f}+-{s_:.1f}' for x, s_ in inw) if inw else 'not measured'
        pr(f'    nu [{a:2d},{b:2d}): weight {wgt:.4f}; breaking |dG| >= {th[0]:8.1f} (0.2) / {th[1]:7.2f} (0.05); measured: {meas}')
    open(os.path.join(RC.C, 'rs_prop_D%g_%s_x%g%s.txt' % (Dg, est, infl, '_ex' + '_'.join('%g' % x for x in exclude) if exclude else '')), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    Dg = float(sys.argv[1]); tags = sys.argv[2].split(',')
    est = sys.argv[3] if len(sys.argv) > 3 else 'dir'
    ex = [float(v) for v in sys.argv[4].split(',')] if len(sys.argv) > 4 else ()
    main(Dg, tags, est, exclude=ex)

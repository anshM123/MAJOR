"""
pwimg/pw3_cost.py (PW_IMAGES.md §4): the exact k = 2 cost of a raw correction in SCALAR de-shift form, and the monopole defect.

By pw1 every raw EE image is  raw = A_g[F]  (scalar de-shift, F even, analytic in |Im mu| <= 3/2).  In the project's normalisation
(n2sub/n2_common: published X, R_-+ of EE, A[f] = (R_- f(nu+i) + R_+ f(nu-i))/(R_- + R_+)), R_-+^EE = R_-+^sc c(nu +- i), so
A_g[F] = A[F/c] with c(mu) = 32(4mu^2+1)(4mu^2+9)/((4mu^2+25)(4mu^2+81)).  EE de-shift form (R1) <=> f = F/c strip-analytic
<=> F(+-i/2) = 0.  Measured here:
  (a) defect(t) := Re sum X A[P_t/c] dnu/Lam - sum psi2 P_t/c dnu/Lam  for t in [0, 12]  (claim: constant = D, since P_t(+-i/2) = 1);
      Omega2(t) = sum psi2 P_t, Omega_c(t) = sum psi2 P_t/c, Omega~(t) = Omega_c + D = the exact cost of the scalar-de-shift atom A_g[P_t];
      Delta_g = 1000 and 300; grid check (every second node).
  (b) controls: F = c P_t (EE de-shift: defect 0), F = P_t/(nu^2 + 9/4) (pole at 3i/2 only, outside the shift strip: defect 0?).
  (c) => cost(A_g[F]) = sum psi2 F/c + D F(i/2)   for every even F = int sigma_F(dt) P_t.
Output: cache/pw3_cost.txt, cache/pw3_cost.json
"""
import os, sys, io, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'n2sub'))
import numpy as np
import n2_common as N          # read-only
from pw_common import cfac

OUT = io.StringIO()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.write(s + '\n')


def grid(Dg, step=1):
    XD = N.xdata(Dg); L = N.lam(Dg)
    nu, X, ps = XD['nu'], XD['X'], XD['psi2']
    sel = nu < L
    return nu[sel][::step], X[sel][::step], ps[sel][::step], XD['dnu']*step, L


def costs(Dg, F, step=1):
    nu, X, ps, dn, L = grid(Dg, step)
    raw = (X*N.Aop(F, nu)).real.sum()*dn/L
    ds = (ps*F(nu)).real.sum()*dn/L
    return raw, ds


def P(t):
    nph = int(600 + 80*t)
    return lambda v: N.Pnu(t, v, nphi=nph)


def main():
    pr(__doc__.split('Output')[0])
    ts = [0.0, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
    res = {}
    for Dg in (1000, 300):
        pr(f'\n(a) Delta_g = {Dg}:   t | Omega2 = sum psi2 P_t | Omega_c = sum psi2 P_t/c | raw cost of A_g[P_t] = Omega~ | defect | '
           f'defect (every 2nd node) | control F = cP_t: raw - ds | F = P_t/(nu^2+9/4): raw - ds')
        rows = []
        for t in ts:
            Pt = P(t)
            Om2 = costs(Dg, Pt)[1]
            raw, dsc = costs(Dg, lambda v: Pt(v)/cfac(v))
            raw2, dsc2 = costs(Dg, lambda v: Pt(v)/cfac(v), step=2)
            r3, d3 = costs(Dg, lambda v: cfac(v)*Pt(v))
            r4, d4 = costs(Dg, lambda v: Pt(v)/(v*v + 2.25))
            rows.append(dict(t=t, Om2=Om2, Omc=dsc, Omt=raw, D=raw - dsc, D2=raw2 - dsc2, c1=r3 - d3, c1rel=(r3 - d3)/abs(d3), c2=r4 - d4))
            pr(f'   {t:5.2f} | {Om2:+.5e} | {dsc:+.5e} | {raw:+.5e} | {raw - dsc:+.7e} | {raw2 - dsc2:+.7e} | {r3 - d3:+.1e} ({(r3 - d3)/abs(d3):.0e} rel) | {r4 - d4:+.3e}')
        Ds = np.array([r['D'] for r in rows])
        pr(f'   => D({Dg}) = {Ds.mean():+.7e}, spread over t in [0, 12]: {Ds.max() - Ds.min():.1e} '
           f'(grid halving changes D by {abs(np.mean([r["D2"] for r in rows]) - Ds.mean()):.1e});  sup_t |Omega~(t)| = {max(abs(r["Omt"]) for r in rows):.4e}'
           f' at t = {rows[int(np.argmax([abs(r["Omt"]) for r in rows]))]["t"]}')
        res[Dg] = dict(D=float(Ds.mean()), rows=rows)
    # (c) generic even F: a mixture, check the formula
    pr('\n(c) formula check on mixtures F = sum_k w_k q_k(nu) P_{t_k}, q_k even polynomials / c-ratios (Delta_g = 1000):')
    D = res[1000]['D']
    rng = np.random.default_rng(3)
    for trial in range(4):
        tk = rng.uniform(0.1, 4, 3); wk = rng.normal(size=3); ak = rng.uniform(-0.3, 0.3, 3)
        F = lambda v: sum(w*(1 + a*v*v)*P(t)(v) for w, a, t in zip(wk, ak, tk))
        Fi2 = sum(w*(1 + a*(-0.25)) for w, a in zip(wk, ak))
        raw, dsc = costs(1000, lambda v: F(v)/cfac(v))
        pr(f'   trial {trial}: raw cost {raw:+.6e}; sum psi2 F/c + D F(i/2) = {dsc + D*Fi2:+.6e}; rel diff {abs(raw - dsc - D*Fi2)/abs(raw):.1e}')
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE, 'cache', 'pw3_cost.json'), 'w'))
    open(os.path.join(HERE, 'cache', 'pw3_cost.txt'), 'w').write(OUT.getvalue())


if __name__ == '__main__':
    main()

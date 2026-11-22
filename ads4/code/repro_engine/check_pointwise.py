"""repro_engine/check_pointwise.py -- CROSS-CHECK ONLY (runs project code read-only): own point-wise EE k = 4 integrand (integrand.py)
vs the project's exact point-wise integrand bf_core.grav_parts (n1_pw / l1_crystal / w0_ee_common / l5_spin.conserved), and own
conserved <TTO> coefficients vs l5_spin.conserved.  Nothing of the project is used by the reproduction itself.
Expected: log|own| - log|project| = const (a-powers, Jacobian constants) at fixed Delta; equal signs up to a global sign; equal log A.
Usage: python3 check_pointwise.py  -> cache/check_pointwise.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
sys.dont_write_bytecode = True
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import integrand as IG

out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


def main():
    import bf_core as B
    import l5_spin as S5
    rng = np.random.default_rng(11)
    n = 24
    X1, X2, Y1, Y2 = (rng.gamma(k, size=n) for k in (6, 4, 2.5, 2.0))
    s = rng.uniform(-0.9, 0.9, n); phi = rng.uniform(-1.3, 1.3, n)
    psi = rng.uniform(0.15, np.pi - 0.15, n); al = psi + phi
    for Dl in (300.0, 1000.0, 2e4):
        mons, sol, _ = S5.conserved(Dl, 0)
        mine = IG.coef(Dl)
        idx = {tuple(m): i for i, m in enumerate(mons)}
        theirs = np.array([sol[idx[(2, 2, 0, 0, 0, 0)], 0], sol[idx[(1, 1, 0, 1, 0, 0)], 0], sol[idx[(0, 0, 0, 2, 0, 0)], 0]]).real
        pr(f'Delta={Dl:g}: conserved coef own/project (normalised to H12^2): {mine/mine[2]} vs {theirs/theirs[2]}; '
           f'max rel diff {np.max(np.abs(mine/mine[2] - theirs/theirs[2])/np.abs(mine/mine[2])):.1e}  (other project monomials: '
           f'{[m for m in mons if tuple(m) not in idx or 0]} ; |sol| outside basis = '
           f'{np.max(np.abs(np.delete(sol[:, 0], [idx[(2,2,0,0,0,0)], idx[(1,1,0,1,0,0)], idx[(0,0,0,2,0,0)]]))) if len(mons) > 3 else 0:.1e})')
        S = B.Sample.__new__(B.Sample)
        S.X1, S.X2, S.Y1, S.Y2, S.s, S.phi, S.al = X1, X2, Y1, Y2, s, phi, np.mod(al, np.pi)
        S.theta = 0.0; S.n = n
        P = B.grav_parts(S, Dl, 0, 'EE')
        lp = (P['G'] + 0.5*Dl*P['lB']).real + np.log(np.abs(P['pref'][0, 0]*P['R']))
        sp_ = np.sign((np.exp(1j*(P['G'] + 0.5*Dl*P['lB']).imag)*P['pref'][0, 0]*P['R']).real)
        M = IG.point(X1, X2, Y1, Y2, s, phi, al, Dl)
        d = M['lG'] - lp
        sg = M['sG']*sp_
        pr(f'   log|own| - log|project|: mean {d.mean():+.6f}, spread (max-min) {d.max() - d.min():.2e};  sign ratio own/project: '
           f'{np.unique(sg)};  max |Im part of project phase| {np.max(np.abs(np.sin((P["G"] + 0.5*Dl*P["lB"]).imag))):.1e}')
        pr(f'   log A own - project: max |diff| {np.max(np.abs(M["lA"] - P["lA"].real)):.2e}; max |Im lA project| {np.max(np.abs(P["lA"].imag)):.1e}')
    open(os.path.join(HERE, 'cache', 'check_pointwise.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()

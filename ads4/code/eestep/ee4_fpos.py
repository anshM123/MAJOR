r"""
ee4_fpos (EESTEP_CLOSE.md §3): f_pos(+-3i/2) for Result Q's actual cut-off de-shift weights.

LRC_CLOSE §3.1: the k = 4 residue at nu = 3i/2 is proportional to f_pos(3i/2) = (position weight x c_l)(3i/2); since c_l has a simple zero
there, f_pos(3i/2) != 0 iff the position weight has a pole at 3i/2.  The position weight is chi_c(nu) x (heavy -> position factor), and
the latter (const/N_com for the scalar-type map; the same for the EE kernel by ee2, whose dressing is nu-independent) is regular and
non-zero near +-3i/2.  chi is even, so the question is whether chi has a pole at nu = -3i/2.

Exact algebra [P]: R_{-3}^{EE3}(-3i/2) = -524288 and R_{-1}^{EE3}(-3i/2) = 6193152/25 (both non-zero), whereas R_{-3}^{sc3}(-3i/2) = 0.
If chi had a simple pole at -3i/2 with residue rho_chi, then g_{-1} := chi R_{-1} would too, and for t -> +infinity
    F_{-1}(t) := int_R chi R_{-1} e^{-i nu t} dnu  ~  -2 pi i rho_chi R_{-1}(-3i/2) e^{-3t/2}.
So  A(t) := F_{-1}(t) e^{3t/2} / (-2 pi i R_{-1}(-3i/2))  -> rho_chi.   (Lemma S4 [D] already makes chi R_{-3} analytic at -3i/2, and with
R_{-3}^{EE3}(-3i/2) != 0 that gives rho_chi = 0 for the exact de-shift solution; this measures it on the stored weights.)
Weights: Q's XP C4_value / C4_deriv (functionals_q7_xprec.json, 64 digits, EE3 R_j), times the HK cutoff c(nu) (nu_c = 150, a = 6);
control: the double-precision sc3 weights of functionals_q4_deshift.json (scalar-type R_j, for which R_{-3}(-3i/2) = 0).
Output: cache/ee4_fpos.txt.
"""
import os, sys, json, time
from fractions import Fraction
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
RQ = os.path.join(HERE, '..', 'resultq')
import numpy as np
import mpmath as mp
import sympy as sp

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def cutoff(nu, nuc=150, a=6):
    return (mp.erf((nu + nuc)/a) - mp.erf((nu - nuc)/a))/2


def load(key_file, key, chikey):
    nus = sp.Symbol('nu')
    if key_file == 'q7':
        x = json.load(open(os.path.join(RQ, 'functionals_q7_xprec.json')))[chikey]
        Rj = json.load(open(os.path.join(RQ, 'functionals_q4_deshift.json')))['EE3']['R_j']
        fr = Fraction(x['dnu']); dnu = mp.mpf(fr.numerator)/fr.denominator
        chi = [mp.mpf(v) for v in x['chi']]
    else:
        q = json.load(open(os.path.join(RQ, 'functionals_q4_deshift.json')))[key]
        Rj = q['R_j']
        dnu = mp.mpf(q['dnu'])
        chi = [mp.mpf(v) for v in json.loads(q[chikey])] if isinstance(q[chikey], str) else [mp.mpf(v) for v in q[chikey]]
    Rs = {int(j): sp.sympify(v, locals={'nu': nus, 'I': sp.I}) for j, v in Rj.items()}
    Rf = {j: sp.lambdify(nus, e, 'mpmath') for j, e in Rs.items()}
    nu = [(k + mp.mpf(1)/2)*dnu for k in range(len(chi))]
    chi = [c*cutoff(v) for c, v in zip(chi, nu)]
    return dict(dnu=dnu, nu=nu, chi=chi, R={j: [Rf[j](v) for v in nu] for j in Rf}, Rs=Rs, nus=nus)


def F(W, j, t):
    """F_j(t) = int_R chi R_j e^{-i nu t} dnu, chi even, R_j(-nu) = R_{-j}(nu)."""
    s = mp.mpc(0)
    for c, v, r, rm in zip(W['chi'], W['nu'], W['R'][j], W['R'][-j]):
        e = mp.expjpi(-v*t/mp.pi)
        s += c*(r*e + rm/e)
    return W['dnu']*s


def run(label, W, ts, dps):
    mp.mp.dps = dps
    nus = W['nus']
    Rm1 = complex(W['Rs'][-1].subs(nus, -3*sp.I/2))
    Rm3 = complex(sp.limit(W['Rs'][-3], nus, -3*sp.I/2)) if sp.denom(sp.together(W['Rs'][-3])).subs(nus, -3*sp.I/2) == 0 else complex(W['Rs'][-3].subs(nus, -3*sp.I/2))
    pr(f'\n== {label}: {len(W["chi"])} nodes, nu_max = {mp.nstr(W["nu"][-1], 6)};  R_-1(-3i/2) = {Rm1:.6g},  R_-3(-3i/2) = {Rm3:.6g}')
    pr('      t     |F_-1| e^{3t/2}    |F_-1| e^{2.5t}    local rate -dlog|F_-1|/dt    |A(t)| = |F_-1 e^{3t/2}/(2 pi R_-1(-3i/2))|     |F_-3| e^{3t}')
    prev = None
    rows = []
    for t in ts:
        f1 = F(W, -1, mp.mpf(t)); f3 = F(W, -3, mp.mpf(t))
        a1 = abs(f1)
        rate = '' if prev is None else f'{-(float(mp.log(a1)) - float(mp.log(prev[1])))/(t - prev[0]):8.3f}'
        A = a1*mp.e**(mp.mpf(1.5)*t)/(2*mp.pi*abs(Rm1))
        rows.append((t, float(A)))
        pr(f'    {t:5.1f}   {mp.nstr(a1*mp.e**(1.5*t), 6):>14}   {mp.nstr(a1*mp.e**(2.5*t), 6):>14}   {rate:>10}                  '
           f'{mp.nstr(A, 6):>14}                              {mp.nstr(abs(f3)*mp.e**(3*t), 6)}')
        prev = (t, a1)
    return rows


if __name__ == '__main__':
    pr('# ee4_fpos: does Q\'s cut-off de-shift weight have a pole at nu = -+3i/2 (f_pos(+-3i/2) != 0)?  (EESTEP_CLOSE.md §3)')
    ts_xp = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0]
    out = {}
    for chikey in ('C4_value', 'C4_deriv'):
        mp.mp.dps = 90                      # BEFORE parsing the 64-digit weights
        W = load('q7', 'EE3', chikey)
        out[chikey] = run(f'Q {chikey} (XP, EE3, x erf cutoff 150/6)', W, ts_xp, 90)
        pr(f'  [{time.time() - T0:.0f}s]')
    ts_dp = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    for chikey in ('chi_value', 'chi_deriv'):
        mp.mp.dps = 40
        W = load('q4', 'sc3', chikey)
        out['sc3_' + chikey] = run(f'CONTROL sc3 {chikey} (double precision, scalar-type R_j)', W, ts_dp, 40)
        pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'ee4_fpos.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')

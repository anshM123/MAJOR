"""
transfer_kernel.py: the signed connection transfer (ANGULAR_ROUTE §106–107, code/signed_transfer/check_transfer.js, imported as a
black box and run unchanged) turned into per-(d, nu)-box kernel bounds.  ANALYTIC_SLIVER3 §9.

Exact identity (signed_transfer README): with X_j = d^{alpha-j} S(alpha - j), Z_j = c_j d^j X_j,
    sum_{j<1500} [dressed Z_j]  = (P X_0 + Q X_{-1})/(1 + q h)   (p = 3, dressed);   = P X_0 + Q X_{-1}   (p = 4),
and the checker certifies |P| <= P_bound, |Q| <= Q_bound over a box d in [d0 +- dd], h in [0, h_cap] (fullTwist: every a >= 998),
nu in [nu0 +- dnu] (real, or complex with fixed nuImag / along nuAxis = imag).  Hence for one branch
    |kappa_J(nu)| = |B(nu)| |d^{-alpha}| |sum_j dressed Z_j| <= |B(nu)| d^{-Re alpha} (P_bound M0 + Q_bound M1),
M0 >= |X_0| = |S^(alpha)|, M1 >= |X_{-1}| = |S^(alpha + 1)| for every state in the box (both channels).

M0, M1 here: max over the a-cover kernel3.a_ranges (998 <= a <= A* = 4500, a >= A* imported) of rt2_kernel.PP_box (Theorem E',
Lemma X, Lemma F, helicity bound for b = 4) -- valid uniformly over the d-slab and the nu-cell.
"""
import os, sys, json, math, subprocess, time, hashlib
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kernel3 as K3
import rt2_kernel as RK
import rt_cert as RC
CHECKER = os.path.abspath(os.path.join(HERE, '..', 'signed_transfer', 'check_transfer.js'))
CACHE = os.path.join(HERE, 'cache', 'transfer')
os.makedirs(CACHE, exist_ok=True)
DIG = 14
DIG_D, DIG_H, DIG_NU = 10, 14, 2     # decimal places (the checker's sparse products must stay exact at scale 1e80)


def _dec_up(x, dig=DIG):
    """decimal string >= x (Fraction) with dig decimal places."""
    s = 10**dig
    v = -((-x.numerator*s)//x.denominator)
    return _fmt(v, dig)


def _fmt(v, dig):
    s = 10**dig
    neg = v < 0; v = abs(v)
    ip, fp = divmod(v, s)
    return ('-' if neg else '') + (f'{ip}.{fp:0{dig}d}' if dig else f'{ip}')


def _dec_near(x, dig=DIG):
    s = 10**dig
    v = round(x*s)
    return _fmt(int(v), dig), Fr(int(v), s)


def hcap_of(d_max):
    """h_cap >= d_max^2/(998^2 (1 + sqrt(1 - d_max))^2) (rounded up; sqrt rounded down)."""
    dm = Fr(d_max)
    sq = Fr(math.isqrt(int((1 - dm)*10**40)), 10**20)       # <= sqrt(1 - d_max)
    return dm*dm/(Fr(998)**2*(1 + sq)**2)


def make_box(name, d_lo, d_hi, nu_lo, nu_hi, m, K=200, R=8, hsplits=16, channels=((3, 0), (3, 4)), nu_imag=None, axis='real',
             max_rem='0.001'):
    dlo, dhi = Fr(d_lo), Fr(d_hi)
    d0s, d0 = _dec_near((dlo + dhi)/2, DIG_D)
    dd = max(dhi - d0, d0 - dlo) + Fr(1, 10**DIG_D)
    dds = _dec_up(dd, DIG_D)
    assert d0 - Fr(dds) <= dlo and d0 + Fr(dds) >= dhi
    hc = hcap_of(float(d0 + Fr(dds)) * (1 + 1e-12))
    hh = hc/2*(1 + Fr(1, 10**6))
    hs = _dec_up(hh, DIG_H)
    nlo, nhi = Fr(nu_lo), Fr(nu_hi)
    n0s, n0 = _dec_near((nlo + nhi)/2, DIG_NU)
    dn = max(nhi - n0, n0 - nlo) + Fr(1, 10**DIG_NU)
    box = dict(name=name, d0=d0s, dd=dds, h0=hs, dh=hs, nu0=n0s, dnu=_dec_up(dn, DIG_NU), m=m, K=K, R=R, hSplits=hsplits,
               channels=[list(c) for c in channels], fullTwist=True, maxRemainder=max_rem, maxTransfer='1000000')
    if nu_imag is not None:
        box['nuImag'] = _dec_up(Fr(nu_imag), DIG_NU)
    if axis == 'imag':
        box['nuAxis'] = 'imag'
    return box


def _key(box):
    return hashlib.sha1(json.dumps(box, sort_keys=True).encode()).hexdigest()[:16]


def run_box(box, verbose=False):
    """run the checker on one box (cached); returns {(p,b): row}."""
    k = _key(box)
    fn = os.path.join(CACHE, f'{box["name"]}_{k}.json')
    if os.path.exists(fn):
        rows = json.load(open(fn))
    else:
        inp = fn + '.in'
        json.dump([box], open(inp, 'w'))
        t0 = time.time()
        r = subprocess.run(['node', CHECKER, inp, fn + '.tmp'], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(r.stderr[-2000:])
        os.replace(fn + '.tmp', fn); os.remove(inp)
        rows = json.load(open(fn))
        if verbose:
            print(f'   checker {box["name"]} m={box["m"]}: {time.time()-t0:.0f}s', flush=True)
    return {(r['p'], r['b']): r for r in rows}


def PQ(d_lo, d_hi, nu_lo, nu_hi, m, channels=((3, 0), (3, 4)), name='box', tries=3, verbose=False, **kw):
    """certified (P, Q) maxima over the channels, adaptively raising m; None if every try fails."""
    for t in range(tries):
        box = make_box(name, d_lo, d_hi, nu_lo, nu_hi, m + 6*t, channels=channels, **kw)
        res = run_box(box, verbose)
        if all(r['status'] == 'PASS' for r in res.values()):
            P = max(float(Fr(r['P_bound']))*(1 + 1e-12) for r in res.values())
            Q = max(float(Fr(r['Q_bound']))*(1 + 1e-12) for r in res.values())
            return P, Q, box['m'], res
        if verbose:
            print('   ', name, 'm', box['m'], {k: (v['status'], v.get('remainder', v.get('error', ''))[:12]) for k, v in res.items()}, flush=True)
    return None


def m_for(nu, L):
    """Taylor order heuristic from the coordinator's preflight table (ANGULAR_ROUTE §107), plus margin."""
    if L >= 3.0:
        return 12 if nu < 400 else (16 if nu < 600 else (22 if nu < 800 else 28))
    if L >= 2.5:
        return 18 if nu < 400 else (24 if nu < 600 else (34 if nu < 800 else 44))
    return 30 if nu < 400 else (40 if nu < 600 else (60 if nu < 800 else 76))


_MC = {}
def M_uniform(sig, t0, t1, d_lo, d_hi):
    """max over 998 <= a <= 4500 (kernel3.a_ranges) of the certified |S^(sigma' + i tau)| bound, sigma' <= sig, tau in [t0, t1]."""
    key = (sig, t0, t1, d_lo, d_hi)
    if key not in _MC:
        _MC[key] = max(RK.PP_box(sig, t0, t1, a_lo, a_hi, d_lo, d_hi) for (a_lo, a_hi) in K3.a_ranges(t1))
    return _MC[key]


def K_cell_transfer(T0, T1, d_lo, d_hi, P, Q):
    """|kappa_J(nu)| d^{1/2} bound for nu in [T0, T1] (one branch, real axis, dressed p = 3): |B| d^{1/2} (P M0 + Q M1), d <= d_hi."""
    M0 = M_uniform(2.5, T0, T1, d_lo, d_hi); M1 = M_uniform(3.5, T0, T1, d_lo, d_hi)
    return RC.B_abs_real(T0)*math.sqrt(d_hi)*(P*M0 + Q*M1), M0, M1


if __name__ == '__main__':
    d_lo, d_hi, n_lo, n_hi, m = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
    t0 = time.time()
    r = PQ(d_lo, d_hi, n_lo, n_hi, m, name='test', verbose=True)
    print('P, Q, m =', r[:3], f'{time.time()-t0:.0f}s')
    for T in range(int(n_lo), int(n_hi)):
        k, M0, M1 = K_cell_transfer(float(T), T + 1.1 if T + 1.1 <= n_hi else n_hi, d_lo, d_hi, r[0], r[1])
        print(f'  nu {T}: M0 {M0:.3g} M1 {M1:.3g}  d^1/2 |kappa| <= {k:.4g}')

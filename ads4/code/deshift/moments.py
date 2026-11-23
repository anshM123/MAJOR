"""deshift/moments.py -- DESHIFT_CLOSE.md §3: the "5x" moment disagreement (RESULT_Q_STATEMENT App. A, last paragraph).

Computes, for every stored k = 4 de-shift weight, the exact (signed) moments  m_n = sum_j int chi R_j (nu - i j)^{2n} dnu  and the
absolute moments  |m|_n = int |chi| sum_j |R_j| |nu - i j|^{2n} dnu,  n = 0..9, by midpoint sums on the weight's own grid, in 50-digit
mpmath for the XP weights (64-digit data) and numpy double for the double-precision weights:
  (a) EE3, Stage 4 (functionals_q4_deshift.json['EE3'], dnu = 1/50, N = 2048, nu_max = 40.95)
  (b) EE3, Stage 7 XP (functionals_q7_xprec.json, dnu = 1/25, N = 3072, nu_max = 122.86)   <- the theorem's Omega^Q_up
  (c) sc3, Stage 4 (functionals_q4_deshift.json['sc3'], the SCALAR-proxy kernel; dnu = 1/50, N = 2048)
and, for (b), the truncation diagnostic  m_n^{tail}(nu_c) = sum over nu > nu_c  (how the signed high moments are built up).
Output: cache/moments.txt
"""
import os, json
import numpy as np
import mpmath as mp
import sympy as sp
import kern

mp.mp.dps = 50
out = []


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)


def mom_np(nu, chi, R, dnu, nmax=10):
    Rv = {j: R[j].np(nu) for j in R}
    m = [sum(np.sum(dnu*chi*Rv[j]*(nu - 1j*j)**(2*n)) for j in R) for n in range(nmax)]
    a = [sum(np.sum(dnu*np.abs(chi*Rv[j])*np.abs(nu - 1j*j)**(2*n)) for j in R) for n in range(nmax)]
    return m, a


def mom_mp(nu, chi_s, R, dnu, nmax=10, nu_c=None):
    m = [mp.mpc(0)]*nmax; a = [mp.mpf(0)]*nmax
    cf = {j: ([mp.mpc(mp.mpf(sp.re(q).p)/sp.re(q).q, mp.mpf(sp.im(q).p)/sp.im(q).q) for q in R[j].num],
              [mp.mpc(mp.mpf(sp.re(q).p)/sp.re(q).q, mp.mpf(sp.im(q).p)/sp.im(q).q) for q in R[j].den]) for j in R}
    for x, cs in zip(nu, chi_s):
        if nu_c is not None and x <= nu_c:
            continue
        c = mp.mpf(cs); z = mp.mpf(repr(float(x)))
        for j in R:
            num = mp.polyval(cf[j][0], z); den = mp.polyval(cf[j][1], z)
            rv = num/den; w = z - 1j*j
            for n in range(nmax):
                p = w**(2*n)
                m[n] += c*rv*p; a[n] += abs(c*rv)*abs(p)
    return [x*dnu for x in m], [x*dnu for x in a]


if __name__ == '__main__' and not os.environ.get('MOM_CONT'):
    R4 = kern.load_C4()
    st = json.load(open(os.path.join(kern.CODE, 'resultq', 'functionals_q4_deshift.json')))
    Rsc = {int(j): kern.Rat(v) for j, v in st['sc3']['R_j'].items()}
    pr("# moments: m_n = sum_j int chi R_j (nu - ij)^(2n) (signed, Re) and |m|_n = int |chi| sum_j |R_j| |nu - ij|^(2n); n = 0..9")
    for lab, key, R in (('(a) EE3 Stage-4 stored (nu_max 40.95)', 'EE3', R4), ('(c) sc3 Stage-4 stored (nu_max 40.95)', 'sc3', Rsc)):
        e = st[key]; dnu = e['dnu']; nu = (np.arange(e['N']) + 0.5)*dnu
        for w in ('chi_value', 'chi_deriv'):
            m, a = mom_np(nu, np.array(e[w]), R, dnu)
            pr(f"{lab} {w}:\n    m_n  = " + " ".join(f"{x.real:+.2e}" for x in m) + "\n    |m|_n = " + " ".join(f"{x:.3g}" for x in a))
    X = kern.load_xp()
    for key in ('C4_value', 'C4_deriv'):
        nu = X[key]['nu']; dnu = mp.mpf(1)/25
        m, a = mom_mp(nu, X[key]['w'], R4, dnu)
        pr(f"(b) EE3 Stage-7 XP (nu_max 122.86) {key}:\n    m_n  = " + " ".join(mp.nstr(x.real, 3) for x in m) +
           "\n    |m|_n = " + " ".join(mp.nstr(x, 4) for x in a))
        for nc in (30.0, 60.0):
            mt, at = mom_mp(nu, X[key]['w'], R4, dnu, nu_c=nc)
            pr(f"    tail nu > {nc:.0f}: m_n = " + " ".join(mp.nstr(x.real, 2) for x in mt))
    open(os.path.join(kern.CACHE, 'moments.txt'), 'w').write("\n".join(out) + "\n")


def continuum():
    """(d) the continuum (cut-off) weights q_deshift_continuum_c4v/c4d (tables at nu = k/50, trapezoid on [0, 60]; the integrand
    chi sum_j R_j (nu - ij)^{2n} is even and analytic in a strip, so the trapezoid rule is exponentially accurate)."""
    R4 = kern.load_C4()
    for piece in ('c4v', 'c4d'):
        rows = [l.split() for l in open(os.path.join(kern.ADS4, 'prl', 'data', f'q_deshift_continuum_{piece}_table.txt')) if l[0] != '#']
        nu = [mp.mpf(r[0]) for r in rows]; chi = [mp.mpf(r[1]) for r in rows]
        cf = {j: ([mp.mpc(mp.mpf(sp.re(q).p)/sp.re(q).q, mp.mpf(sp.im(q).p)/sp.im(q).q) for q in R4[j].num],
                  [mp.mpc(mp.mpf(sp.re(q).p)/sp.re(q).q, mp.mpf(sp.im(q).p)/sp.im(q).q) for q in R4[j].den]) for j in R4}
        m = [mp.mpc(0)]*10; a = [mp.mpf(0)]*10
        for i, (z, c) in enumerate(zip(nu, chi)):
            wt = mp.mpf(1)/50*(mp.mpf(1)/2 if i == 0 else 1)
            if z == 0:
                continue
            for j in R4:
                rv = mp.polyval(cf[j][0], z)/mp.polyval(cf[j][1], z); w = z - 1j*j
                for n in range(10):
                    p = w**(2*n); m[n] += wt*c*rv*p; a[n] += wt*abs(c*rv)*abs(p)
        pr(f"(d) continuum {piece} (exact definition in prl/data; tail beyond nu = 60 < 1e-25 relative):\n    m_n  = " +
           " ".join(mp.nstr(x.real, 3) for x in m) + "\n    |m|_n = " + " ".join(mp.nstr(x, 4) for x in a) +
           "\n    ratio |m|_(n+1)/|m|_n = " + " ".join(mp.nstr(a[n + 1]/a[n], 3) for n in range(9)))


if __name__ == '__main__' and os.environ.get('MOM_CONT'):
    out.clear()
    continuum()
    open(os.path.join(kern.CACHE, 'moments.txt'), 'a').write("\n".join(out) + "\n")

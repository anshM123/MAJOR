"""
h13_exact: exact scalar-skeleton remainder rho_k and the exact D-safe action F = F_hyb + S on large-impact-parameter states (H13').
Method: h13_laplace.Laplace3H (= skelexact3.Laplace3 with the n-sum tail quadrature in Arb; Arb Laplace/Barnes, no contour; even k = 2,4,6,8,10 exact, odd k = 3,5,7 by quartic interpolation in k of
m^{2k}C_k/E_k [M], exactly as e8_states_arb.py), with the Mack polynomial from mack_mat / mack_flint (FLINT; identical output to mack.mack_P_fast,
~200-1000x faster), which makes J ~ 100 reachable.  States: Delta on a 1/64 grid (h13_common.exact_state).  Both skeletons Df = 3, 2 share
the Mack polynomial.  rho_k is the e8 definition; S and F_hyb, N_s as e7_nf.py / nf_common.action_nf (64-node p rule, K* kernel:
normalised response ratio = 1, evenkernel.py / e10).
Usage: python3 h13_exact.py Dg LIST   LIST = "mfl:J,mfl:J,..."   -> cache/h13_ex_<Dg>.pkl (incremental, resumable; key (mfl_target, J))
n-sum error estimate: several (N0, nq, L) configurations; 'n-sum err' = max |S - S_alt|.  NCONF history (stored per record):
2 configs (200,320,1),(320,320,.5) for the first records; then (120,320,1),(200,320,1),(320,320,.5); final (480,480,.5),(320,480,.25),
(320,320,.5); follow-up (H13 v2): (200,960,1),(320,1280,.35),(320,960,.5) -- the Dg = 1000 large-J failures were the tail Gauss
order nq (nq = 320 unresolved at nu -> Lam; nq >= 800 converged to 1e-11), while N0 >~ 600 hurts near the twist gap (Dg = 300).
Env H13_REDO=1 recomputes.
The Euler-Maclaurin tail is reliable only while J << tau = Delta - J (it fails near the twist gap at Dg = 100, J >~ 45): records
with n-sum err/N_s > 1e-3 are flagged as unreliable in h13_summary.py.
Env: H13_DFS (default "3,2"), H13_PREC (bits, default 400; raised automatically if the Arb radius is > 1e-25 relative).
"""
import os, sys, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0, eval_legendre
import h13_common as HC
from h13_common import F, KS, SIGMA, Lam_of, exact_state, flat_rows, hybrid
from common import gl01
import mack, mack_flint, mack_mat
from fractions import Fraction as _Fr


_MCACHE = {}


def _mackP(Delta, J):
    # same exact polynomial as mack.mack_P_fast (self-tests in mack_flint.py, mack_mat.py); mack_mat needs Delta in (1/64)Z.
    # Disk cache (exact rationals, pickled) for J >= 60: cache/mack/<Delta*64>_<J>.pkl
    key = (float(Delta), int(J))
    if key in _MCACHE:
        return _MCACHE[key]
    fn = os.path.join(HC.CACHE, 'mack', '%d_%d.pkl' % (round(Delta*64), J)) if (_Fr(Delta)*64).denominator == 1 else None
    if fn and os.path.exists(fn):
        P = pickle.load(open(fn, 'rb'))
    elif (_Fr(Delta)*64).denominator == 1:
        P = mack_mat.mack_P_mat(Delta, J)
        if J >= 60:
            os.makedirs(os.path.dirname(fn), exist_ok=True)
            pickle.dump(P, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
    else:
        P = mack_flint.mack_P_flint(Delta, J)
    _MCACHE[key] = P
    return P


mack.mack_P_fast = _mackP
import h13_laplace
h13_laplace.mack = mack
from h13_laplace import Laplace3H as Laplace3
# n-sum configurations (N0, nq, L-factor): the result is taken from the last one, the spread between them is the n-sum error estimate
# (t_h13_nsum: at J >~ 40 the Euler-Maclaurin tail error grows with N0 beyond ~400 and is large for small N0, so two moderate N0).
NCONF = ((200, 960, 1.0), (320, 1280, 0.35), (320, 960, 0.5))
from nf_common import envk
import remk

NP = 64
pn, pw = gl01(NP)


def nodes_for(bh):
    # Gauss-Legendre in p must resolve oscillations of frequency up to ~1.35 bhat: 64 nodes suffice for bhat <= 100 (error
    # ~ (e w/(8n))^{2n} < 1e-18); beyond, n = 0.75 bhat + 32.  (h13_exact NP column; checked against 64 nodes at bhat <= 100)
    n = 64 if bh <= 100 else int(0.75*bh) + 32
    return n, gl01(n)


KE = (2, 4, 6, 8, 10)


def Ek(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2)


def rho_state(Dg, Df, st, prec, pn=pn):
    Lam = Lam_of(Dg)
    J = st['J']; u = 1/st['mfl']**2
    x = pn*np.sqrt(u); nus = Lam*pn
    exs = []
    for N0, nq, Lf in NCONF:
        while True:
            Ld = max(4.0*N0, 0.5*st['Delta']**2/max(2*Df - 1.5, 0.5))*Lf
            L3 = Laplace3(st['Delta'], J, Df, ks=KE, numax=float(nus.max()), prec=prec, N0=N0, nq=nq, L=Ld)
            exE, rad = {}, 0.0
            for k in KE:
                vals = [L3.value_fast(k, nu) for nu in nus]      # = L3.value (t_h13_fast.py)
                exE[k] = np.array([float(v.mid()) for v in vals])
                rad = max(rad, max(float(v.rad())/max(abs(float(v.mid())), 1e-300) for v in vals))
            if rad < 1e-25 or prec >= 1600:
                break
            prec *= 2
        exs.append(exE)
    out = [_rho_from(exE, J, x, nus, st) for exE in exs]
    r = out[-1]
    r.update(rho_alt=out[0]['rho'], rho_alts=[o['rho'] for o in out[:-1]], rad=rad, prec=prec, M=L3.M)
    return r


def _rho_from(exE, J, x, nus, st):
    ex = {k: exE[k] for k in (2, 4, 6)}
    for ko in (3, 5, 7):
        r = sum(exE[ka]/Ek(ka, x)*np.prod([(ko - kb)/(ka - kb) for kb in KE if kb != ka]) for ka in KE)
        ex[ko] = r*Ek(ko, x)
    # odd-k interpolation sensitivity: cubic from k = 2..8 (drop k = 10)
    ex3 = {}
    for ko in (3, 5, 7):
        K4 = (2, 4, 6, 8)
        r = sum(exE[ka]/Ek(ka, x)*np.prod([(ko - kb)/(ka - kb) for kb in K4 if kb != ka]) for ka in K4)
        ex3[ko] = r*Ek(ko, x)
    Pd = remk.Pders_vec(nus, st['eta'], nmax=0)[0]
    rho, rho3 = {}, {}
    base = lambda k: Ek(k, x)*eval_legendre(J, 1 - 2*x*x) + 2*(Pd - j0((2*J + 1)*x))
    for k in KS:
        E = Ek(k, x)
        rho[k] = (ex[k] - base(k))/E
        rho3[k] = (ex3[k] - base(k))/E if k in ex3 else rho[k]
    return dict(rho=rho, rho3=rho3, ex=ex, Pd=Pd)


def S_of(h, mfl, rho, pn=pn, pw=pw):
    s = SIGMA[h]; u = 1/mfl**2; x = pn*np.sqrt(u)
    S = Sabs = 0.0; Sk = {}
    for k in KS:
        wk = u**(k - 2)*F.phi(k, pn)*pw*envk(k, x)
        Sk[k] = float(np.sum((s**k)*wk*rho[k])); S += Sk[k]; Sabs += float(np.sum(np.abs(wk*rho[k])))
    return S, Sabs, Sk


def main():
    Dg = int(sys.argv[1])
    todo = [(float(a.split(':')[0]), int(a.split(':')[1])) for a in sys.argv[2].split(',') if a]
    dfs = [float(a) for a in os.environ.get('H13_DFS', '3,2').split(',')]
    prec0 = int(os.environ.get('H13_PREC', '400'))
    pk = os.path.join(HC.CACHE, 'h13_ex_%d.pkl' % Dg)
    db = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    for mflt, J in todo:
        st = exact_state(Dg, mflt, J)
        if st is None:
            print(f"Dg={Dg} m/Lam={mflt} J={J}: unphysical, skipped", flush=True); continue
        for Df in dfs:
            key = (round(mflt, 6), J, Df)
            if key in db and not os.environ.get('H13_REDO'):
                continue
            t0 = time.time()
            tm = time.time(); mack.mack_P_fast(st['Delta'], J); tmack = time.time() - tm
            npn, (qn, qw) = nodes_for((2*J + 1)/st['mfl'])
            r = rho_state(Dg, Df, st, prec0, pn=qn)
            rec = dict(Dg=Dg, Df=Df, mfl=st['mfl'], mflt=mflt, J=J, Delta=st['Delta'], beta=st['beta'], eta=st['eta'], m=st['m'],
                       bh=(2*J + 1)/st['mfl'], rho=r['rho'], rho3=r['rho3'], rho_alt=r['rho_alt'], NCONF=NCONF, NP=npn, rad=r['rad'], prec=r['prec'], M=r['M'], hel={})
            for h in ('++', '+-'):
                if (h == '++' and J % 2) or (h == '+-' and J < 4):
                    continue
                B2, B3 = flat_rows(h, 1/st['mfl']**2, [J])[J]
                hy = hybrid(h, Dg, st['mfl'], J, B2, B3, beta=st['beta'])
                S, Sabs, Sk = S_of(h, st['mfl'], r['rho'], qn, qw)
                S3, _, _ = S_of(h, st['mfl'], r['rho3'], qn, qw)
                Salts = [S_of(h, st['mfl'], ra, qn, qw)[0] for ra in r['rho_alts']]
                Salt = max(Salts, key=lambda v: abs(v - S))
                rec['hel'][h] = dict(Fh=hy['Fh'], Ns=hy['Ns'], FR=hy['FR'], Fflat=hy['Fflat'], S=S, Sabs=Sabs, Sk=Sk, S3=S3, Salt=Salt, Salts=Salts,
                                     Fex=hy['Fh'] + S)
            rec['time'] = time.time() - t0; rec['tmack'] = tmack
            db[key] = rec
            msg = '  '.join(f"{h}: F/Ns {v['Fex']/v['Ns']:+.5f} S/Ns {v['S']/v['Ns']:+.2e} (n-sum err {abs(v['S'] - v['Salt'])/v['Ns']:.0e}, odd-k {abs(v['S'] - v['S3'])/v['Ns']:.0e}) |S|abs/Ns {v['Sabs']/v['Ns']:.2e}"
                            for h, v in rec['hel'].items())
            print(f"Dg={Dg} Df={Df:g} m/Lam={st['mfl']:.4f} J={J:3d} bhat={rec['bh']:6.1f} beta={st['beta']:.3f}  {msg}  "
                  f"[NP {npn}, rad {r['rad']:.0e}, prec {r['prec']}, M {r['M']}, mack {tmack:.0f}s, {rec['time']:.0f}s]", flush=True)
            pickle.dump(db, open(pk + '.tmp', 'wb')); os.replace(pk + '.tmp', pk)
    print('# done', flush=True)


if __name__ == '__main__':
    main()

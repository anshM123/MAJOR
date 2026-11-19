"""
u4 (EE_FUNCTIONAL.md §4): (a) real zeros of the forward normalisation of S(lam) = EE + lam q^2 (exact rational R_-+, u2);
(b) NF02 per-state check and NF02-family LP (s2_lp machinery, model M': diagonal = stripped flat B2+B3 + leading AdS curvature, identical
for every crystal-flip-free structure with the q^2 stripped density) with the exact leading-Regge off-diagonal of S(lam):
    o_EE = 0 identically (OFFDIAG §3.4)  =>  Y_S(nu, beta) = lam o_q/(n_EE + lam n_q) = Y_q(nu, beta) * lam/(r(nu) + lam),  r = n_EE/n_q
(r exact from u2; Y_q from the s1 tables; OFF=psi, i.e. no de-shift on the off-diagonal, as s2/s3).  lam = 0 is EE alone (off-diagonal 0).
Modes:  zeros                     -> cache/u4_zeros.txt
        check|lp  (env U4_DG, U4_LAM (comma list), U4_BMAX)   -> cache/u4_<mode>_<Dg>.txt (+ .pkl for lp)
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import s2_lp as S2
import u2_hclosed_struct as U2

CACHE = os.path.join(HERE, 'cache')
F = U2.numfuncs()
nE = lambda v: F['EE'][0](np.asarray(v, complex)) + F['EE'][1](np.asarray(v, complex))
nQ = lambda v: F['q2'][0](np.asarray(v, complex)) + F['q2'][1](np.asarray(v, complex))
rfun = lambda v: nE(v)/nQ(v)


class Ylam:
    def __init__(self, Yq, lam):
        self.Yq, self.lam, self.B = Yq, lam, Yq.B

    def __call__(self, nu, be):
        nu = np.asarray(nu, float)
        if self.lam == 0.0:
            return np.zeros(len(nu), complex)
        return self.Yq(nu, be)*self.lam/(rfun(nu) + self.lam)


def zeros_mode():
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    nu = np.concatenate([np.linspace(1e-4, 20, 2000001), np.geomspace(20, 1e5, 20001)[1:]])
    r = rfun(nu)
    pr("# u4 zeros: S(lam) = EE + lam q^2 has a real zero of h_- + h_+ at nu iff r(nu) := n_EE/n_q = -lam (r complex in general)")
    s = np.sign(r.imag); idx = np.where(s[:-1]*s[1:] < 0)[0]
    pr(f"Im r = 0 crossings on (0, 1e5): {len(idx)}")
    forb = []
    for i in idx:
        t = r.imag[i]/(r.imag[i] - r.imag[i + 1]); v = nu[i] + t*(nu[i + 1] - nu[i]); rr = r.real[i] + t*(r.real[i + 1] - r.real[i])
        forb.append(-rr)
        pr(f"   nu = {v:.6f}: Re r = {rr:+.6f}  -> forbidden lam = {-rr:+.6f}")
    pr(f"small nu: r(1e-4) = {r[0]:.6f}; large nu: r(1e5) = {r[-1]:.8f}")
    pr("\nrobustness: m(lam) := min_nu |n_EE + lam n_q|/(|n_EE| + |lam n_q|)  (1 = no cancellation; 0 = real zero), and where:")
    for lam in (-0.9, -0.5, -0.2, -0.1, -0.05, -0.02, 0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0):
        q = np.abs(r + lam)/(np.abs(r) + abs(lam))
        j = np.argmin(q[nu < 200])
        pr(f"   lam = {lam:+6.2f}: m = {q[j]:.4f} at nu = {nu[j]:.3f};  large-nu off-diagonal factor lam/(1+lam) = {lam/(1 + lam):+.4f}")
    open(os.path.join(CACHE, 'u4_zeros.txt'), 'w').write(buf.getvalue())


def lp_mode(mode):
    Dg = int(os.environ.get('U4_DG', '300'))
    bmax = float(os.environ.get('U4_BMAX', '2.2'))
    lams = [float(v) for v in os.environ.get('U4_LAM', '0').split(',')]
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    Lm = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    Yq = S2.Yint(S2.load_table(('C2a', 'C2b')))
    mfls = [float(v) for v in os.environ.get('U4_MFLS', '1.02,1.05,1.1,1.2,1.3,1.5,2.0,3.0,5.0').split(',')]
    sts = S2.grid_states(Dg, bmax, mfls)
    v_nf = S2.R3.nf_vector()
    res = {}
    for lam in lams:
        Yf = Ylam(Yq, lam)
        blocks = [S2.state_block(Dg, mfl, J, Yf, 'psi' if lam != 0 else 'none') for mfl, J in sts]
        pr(f"# u4 [{mode}] Delta_g = {Dg} (Lam = {Lm:.1f}); S = EE + {lam} q^2 (off-diagonal = Y_q lam/(r + lam), OFF=psi); "
           f"{len(blocks)} states, beta <= {bmax}  [{time.time() - t0:.0f}s]")
        if mode == 'check':
            allneg = []
            for mfl in mfls:
                bl = sorted([b for b in blocks if b['mfl'] == mfl], key=lambda b: b['bh'])
                if not bl:
                    continue
                ls = [(b['bh'], S2.R3.lam(b, v_nf)[0], b['beta']) for b in bl]
                allneg += [l for _, l, _ in ls]
                cross = next((bh for bh, l, be in ls if l < 0), None)
                pr(f"  m/Lam = {mfl}: first negative at b-hat = {cross}; min lam/N_s = {min(l for _, l, _ in ls):+.4f}; (b-hat, beta, lam/N_s): " +
                   " ".join(f"({bh:.0f},{be:.2f},{l:+.3f})" for bh, l, be in ls[::4]))
            pr(f"  overall min lam/N_s = {min(allneg):+.4f}")
            res[lam] = min(allneg)
        else:
            G = S2.NL.NFGrid(S2.D, S2.RR)
            S2.R3.summary(blocks, v_nf, pr, f'NF02, S = EE + {lam} q2')
            st, v, info = S2.R3.solve(G, blocks, 0.2, 0.0, Lm, 'A', v0=v_nf)
            if st != 0:
                pr(f"  LP + AdS rows: status {st}: {info.get('msg')}"); res[lam] = None
            else:
                pr(f"  LP + AdS rows: A = {info['A']:.2f}, B_AdS = {info['B_AdS']:.2f}, bound 100 / 300 / 1000: "
                   f"{S2.R3.bound(info['A'], info['B_AdS'], 100):.3e} / {S2.R3.bound(info['A'], info['B_AdS'], 300):.3e} / "
                   f"{S2.R3.bound(info['A'], info['B_AdS'], 1000):.3e}")
                res[lam] = dict(v=v, **info)
        pr(f"  [{time.time() - t0:.0f}s]")
    tag = f"{Dg}_b{bmax}_" + '_'.join(f"{l:g}" for l in lams)
    open(os.path.join(CACHE, f'u4_{mode}_{tag}.txt'), 'w').write(buf.getvalue())
    if mode == 'lp':
        pickle.dump(res, open(os.path.join(CACHE, f'u4_lp_{tag}.pkl'), 'wb'))


if __name__ == '__main__':
    m = sys.argv[1]
    zeros_mode() if m == 'zeros' else lp_mode(m)

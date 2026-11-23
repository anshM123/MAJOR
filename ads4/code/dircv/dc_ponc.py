"""
dc_ponc (DIRCV.md §8): "both waves on one contour" twin: the state integrand with P_J = W+ + W- (scalar and graviton, mode 0) on the
per-wave contour C+ (flow and saddle of the W+ wave flow function, mode 1).  P_J is entire, so there is no Q cut to sweep: unbiased
like the n3 P_J twin, with (hopefully) the wave contour's variance.  Reference: P_J on n3 (as production), same CRN.
Saved as pw_points records (modes [0], 'ponc' flag) to dircv/cache/pw_pts_J<J>_D<D><tag>P.pkl (state + its reference).
Usage: python3 dc_ponc.py J Delta nu R seed tag [noref]
"""
import dc_guard, os, sys, pickle, time
os.environ.setdefault('L1_DEG', '6'); os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
import pw_points as PP
import pw_core as Co
import pw_wave as WV
import n3_flow as F
import n3_hybrid as Hy
import z3_nodes as Z3


class PonC(Co.WEv):
    """flow (phi) of mode +1, integrand (pint, grav) of P_J."""
    def pint(self, xi):
        zl, P, t = self._common(xi, self.nu)
        return P, WV.wave_derivs(self.J, t, 0, 0)[0]


def make_ponc(nu, D, J, kinds=('EE', 'C3H')):
    nuf = Co.nuf_of(nu)
    c, wc, L, gm = Co.saddle(nuf, D, J, 1)
    ev = PonC(nu, D, J, c, wc, L, 1, 'S', nuf=nuf)
    F.set_branch(ev, fns=('phi', 'pint'))
    gk = {}
    for k in kinds:
        ge = Co.WEv(nu, D, J, c, wc, L, 0, k, nuf=nuf, lam=1.0)      # mode 0 = P_J contraction, on the C+ points
        F.set_branch(ge, fns=('gravlog',)); ge.saddle_refs['grav'] = ge.saddle_refs['gravlog']
        gk[k] = ge
    return ev, gk


if __name__ == '__main__':
    J, Ds, nu, R, seed, tag = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    noref = len(sys.argv) > 7 and sys.argv[7] == 'noref'
    PP.setup(); t0 = time.time()
    prop = Hy.adapt(Co.make_ev(nu, Ds, J, 1), 3.5)
    out = {}
    ev, gk = make_ponc(nu, Ds, J)
    res = Co.block([ev], prop, 9, R, seed, gkinds={id(ev): gk})
    out[(Ds, round(nu, 9))] = dict(J=J, Delta=Ds, nu=nu, m=Z3.mphys(Ds, J), modes=[0], R=R, mq=9, nuf=ev.nuf, adapt=True, scale=1.0, ponc=True,
                                   per_mode=[{k: v for k, v in o.items() if k != 'g'} for o in res],
                                   g={k: [o['g'][k] for o in res] for k in gk}, time=time.time() - t0, wbranch=[None])
    print(f"  PonC state nu={nu}: scalar rel err {abs(res[0]['err']/res[0]['E']):.4f}; graviton EE rel err "
          f"{np.round(np.abs(res[0]['g']['EE']['err']/res[0]['g']['EE']['E']), 4).tolist()}; diag {res[0]['diag']}; {time.time()-t0:.0f}s", flush=True)
    if not noref:
        evr = Co.make_ev(nu, Z3.DREF, J, 0)
        gr = {k: Co.make_gev(evr, k) for k in ('EE', 'C3H')}
        rr = Co.block([evr], prop, 9, R, seed, gkinds={id(evr): gr})
        out[(Z3.DREF, round(nu, 9))] = dict(J=J, Delta=Z3.DREF, nu=nu, m=Z3.mphys(Z3.DREF, J), modes=[0], R=R, mq=9, nuf=evr.nuf, adapt=True,
                                            scale=1.0, per_mode=[{k: v for k, v in o.items() if k != 'g'} for o in rr],
                                            g={k: [o['g'][k] for o in rr] for k in gr}, time=time.time() - t0, wbranch=[None])
    pickle.dump(out, open(os.path.join(dc_guard.CACHE, 'pw_pts_J%d_D%g%sP.pkl' % (J, Ds, tag)), 'wb'))
    print(f"## done {tag} in {time.time()-t0:.0f}s", flush=True)

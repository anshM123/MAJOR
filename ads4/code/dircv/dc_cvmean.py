"""
dc_cvmean (DIRCV.md §2.1): is the control variate's mean (the exact scalar state/reference ratio) exactly known?
At (1100, 100) nodes: lx_exact(state) - lx_exact(ref) with (i) production (rho_2 from sm_exact, NLOC=10 local barycentric continuation
to nu +- i), (ii) rho_2 from the even / odd halves of the nodes, (iii) NLOC = 8, 12, (iv) rho_2 := 0 (total size of the non-hybrid part);
reference: rho_2 := 0 (production) vs a rho_2(ref) scaled from the state as (m_st/m_ref)^2 (bound on the reference's omission).
Errors expressed as scalar Delta c2 (= m^2 Re(d) ex2/E2) and as relative log error (multiply by the graviton factor F_g, dc_gfac).
"""
import dc_guard, os, sys
import numpy as np
import z3_nodes as Z3, fw_scalar as FS, n3_exact as NX
J, D = 100, 1100.0
FS.install(J, D)
rs = Z3.exact_rec(D, J); rr = Z3.exact_rec(Z3.DREF, J)
lq = Z3.log_q(D, J) - Z3.log_q(Z3.DREF, J)
m = Z3.mphys(D, J)
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
pr(f"# CV mean at (1100,100); state rho_2 src {rs['src']}; {len(rs['nus'])} nodes nu in [{rs['nus'].min():.2f}, {rs['nus'].max():.2f}]")
pr(f"   max |rho_2| on nu<40: {np.max(np.abs(rs['rho'][2][rs['nus'] < 40])):.3e};  m^2 * that = {m*m*np.max(np.abs(rs['rho'][2][rs['nus'] < 40])):.1f}")
ev = np.arange(0, len(rs['nus']), 2); od = np.arange(1, len(rs['nus']), 2)
def lx(nu, r, nodes=None, zero=False):
    rr_ = dict(r); 
    if zero:
        rr_['rho'] = {2: 0*r['rho'][2]}
    return np.log(complex(NX.X(nu, rr_, nodes=nodes)[0]))
for nu in (10.0, 13.3, 16.3, 16.7, 20.8, 25.0, 33.4):
    x = nu/m; E2 = (2 - x*x)/(1 - x*x)**2
    ex2 = float(np.interp(nu, rs['nus'], rs['ex'][2])); fac = m*m*ex2/E2
    l0 = lx(nu, rs)
    le, lo_ = lx(nu, rs, ev), lx(nu, rs, od)
    old = NX.NLOC
    NX.NLOC = 8; l8 = lx(nu, rs); NX.NLOC = 12; l12 = lx(nu, rs); NX.NLOC = old
    lz = lx(nu, rs, zero=True)
    lref = lx(nu, rr)
    # reference: production rho_ref = 0; bound = state rho_2 scaled by (m_st/m_ref)^2 (1/m^2 law of rho_2)
    rr2 = dict(rr); mr = rr['m']
    rho_s = np.interp(rr['nus'], rs['nus'], rs['rho'][2])
    rr2['rho'] = {2: rho_s*(m/mr)**2}
    lref2 = lx(nu, rr2)
    d = lambda a: (a - l0).real
    pr(f"nu={nu:5.1f}: scalar c2-factor {fac:+.3g};  rel log errors: even {d(le):+.2e} odd {d(lo_):+.2e} NLOC8 {d(l8):+.2e} NLOC12 {d(l12):+.2e} "
       f"| rho2:=0 {d(lz):+.2e} | ref rho scaled {(lref2 - lref).real:+.2e};  in scalar c2: interp <= "
       f"{abs(fac)*max(abs(d(le)), abs(d(lo_)), abs(d(l8)), abs(d(l12))):.2f}, rho2 total {fac*d(lz):+.1f}, ref {fac*(lref2-lref).real:+.3f}")
open(os.path.join(dc_guard.CACHE, 'dc_cvmean.txt'), 'w').write('\n'.join(lines) + '\n')

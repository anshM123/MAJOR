"""
largespin/ls_struct.py (LARGESPIN_ANALYTIC.md §3-§4): per-state impact-parameter tolerance profile and the cost of candidate
beyond-hybrid correction STRUCTURES, over the H_LJ(1000) and H_R' regions.

Exact ingredients (ls_core): flat NF02 parts (exact recursion), leading-Regge curvature (Mehler-Dirichlet), OmLR_k(t).
Per state (worse helicity): floor = N_flat/5 + min(0, curv)   [certified flat floor + exact curvature].
(1) Impact-parameter tolerance profile (k = 2 piece, sign-free, exact):  a correction whose H2 density puts weight w at radius t costs
    <= 2 w |OmLR_2(t)|;  tol(t) := floor/(2|OmLR_2(t)|)  = the weight (in units of the leading-Regge ring weight 2) tolerated at radius t.
    Printed at t = beta (ring), 0.8 beta, beta/3, beta/10, 1/Lam (bulk point), and the sup-form tau(eps) for supp in [(1-eps)beta, inf).
(2) Structures (k = 2 image, correction dI(nu) = coefficient x structure; tolerated |coefficient| a_tol = floor/|2 int psi2 S|):
    A  (a/m)  P(cosh beta)                       amplitude (ring)                       natural a ~ 21 (engine nu-flat, J = 0..100)
    B  (a/m^2) nu^2 P(cosh beta)                 nu-dependent amplitude (ring derivative: -(Lap+1/4) OmLR)
    C  (a/m)  beta d/dbeta P(cosh beta)          displacement / phase of the ring
    D  (a/m)  Re P_{1/2+i nu}(cosh beta)          index-shifted Regge function (H2 density inside the disc t <= beta)
    E  (a/m)  P(cosh beta) (1/4)/(nu^2 + 1/4)    H2-resolvent smearing of the ring (density spread over t ~ O(1))
    F  (a/m^2) * 1                               smooth, nu-flat absolute (density at t <~ 1/Lam): the 'smooth uniform c2' class
    G  (a/m^2) exp(-(nu beta/3)^2)               window-smooth bump (density smooth on scale beta/3)
    (a for F, G is in c2 units; for A, C, D, E it multiplies 1/m relative to the ring.)
Output: cache/struct_<Dg>.txt, cache/struct_<Dg>.json.   Usage: python3 ls_struct.py [Dg]
"""
import os, sys, json, time
sys.dont_write_bytecode = True
import numpy as np
import ls_core as C

Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
L = C.lam(Dg)
F = C.NFw(Dg)
pg, wg = C.gl01(3000)
ps2 = F.psi2(pg)*wg
ph, wph = C.gl01(3000); ph = ph*np.pi/2; wph = wph*np.pi/2


def Pmat(beta, shift=False):
    """P_{-1/2+i Lam p}(cosh beta) (or Re P_{1/2 + i Lam p}) on the p-nodes, Mehler-Dirichlet in the Lemma-C variables."""
    s = 2*np.arcsinh(np.sinh(beta/2)*np.sin(ph))
    wt = wph/np.cosh(s/2)*(np.cosh(s) if shift else 1.0)
    out = np.empty_like(pg)
    for i0 in range(0, len(pg), 500):
        out[i0:i0 + 500] = (2/np.pi)*(np.cos(np.outer(L*pg[i0:i0 + 500], s))*wt).sum(1)
    return out


def states():
    S = []
    for J in [26, 30, 34, 38, 42, 46, 50, 64, 100, 150, 200, 250, 300, 344, 33, 51, 101, 201, 301]:
        mg = C.mgap(Dg, J)
        for f in (1.0, 1.1, 1.2):
            if mg*f <= 1.3:
                S.append(('LJ', J, mg*f))
    for mfl in (1.3, 1.5, 2, 3, 5, 10, 20):
        Jmax = int(((mfl*L)**2/(Dg - 2) - (Dg - 2) - 1)/2)
        for b in (5, 7, 10, 15, 20, 30, 50, 80, 120, 200, 300, 500, 800, 1200, 2000):
            J = int(round((b*mfl - 1)/2))
            if 4 <= J <= min(Jmax, 3000):
                S.append(('R', J, mfl))
    return S


rows = []
t0 = time.time()
for reg, J, mfl in states():
    m = mfl*L
    best = None
    for h in ('++', '+-'):
        if (h == '++' and J % 2) or (h == '+-' and J < 4):
            continue
        st = C.hyb_fast(Dg, J, mfl, h)
        Nf = abs(st['B2']) + abs(st['B3']); cv = st['cv2'] + st['cv3']
        fl = Nf/5 + min(0.0, cv)
        r = dict(h=h, floor=fl, Ns=st['Ns'], lam=st['base']/st['Ns'], beta=st['beta'])
        if best is None or fl/st['Ns'] < best['floor']/best['Ns']:
            best = r
    if best is None:
        continue
    beta = best['beta']; fl = best['floor']; Ns = best['Ns']
    Om = lambda t: C.OmLR(2, [t], Dg)[0]
    prof = {k: fl/(2*abs(Om(t))) for k, t in [('ring', beta), ('0.8b', 0.8*beta), ('b/3', beta/3), ('b/10', beta/10), ('1/L', 1.0/L)]}
    # sup-form over [(1-eps) beta, inf): OmLR_2 is monotone non-increasing for t >= 2e-3 (checked on the grid in ls_barrier), so sup = value
    tau = {e: fl/(2*abs(Om((1 - e)*beta))) for e in (0.0, 0.2, 0.5)}
    P0 = Pmat(beta); P1 = Pmat(beta, shift=True)
    nu = L*pg
    h1 = 1e-3*beta
    Op, O0, Om_ = Om(beta + h1), Om(beta), Om(beta - h1)
    d1 = (Op - Om_)/(2*h1); d2 = (Op - 2*O0 + Om_)/h1**2
    lapO = d2 + d1/np.tanh(beta)
    cost = dict(A=2*O0/m, B=2*(-(lapO + 0.25*O0))/m**2, C=2*beta*d1/m, D=2*(ps2*P1).sum()/m,
                E=2*(ps2*P0*0.25/(nu**2 + 0.25)).sum()/m, F=2*ps2.sum()/m**2, G=2*(ps2*np.exp(-(nu*beta/3)**2)).sum()/m**2)
    chkB = 2*(ps2*nu**2*P0).sum()/m**2
    atol = {k: fl/abs(v) for k, v in cost.items()}
    row = dict(reg=reg, J=J, mfl=mfl, m=m, bhat=(2*J + 1)/mfl, beta=beta, Ns=Ns, lam=best['lam'], h=best['h'], floor_rel=fl/Ns,
               prof=prof, tau=tau, atol=atol, chkB=chkB/cost['B'] - 1)
    rows.append(row)
    print('%3s J=%4d m/L=%6.3f b=%8.2f beta=%.4f lam=%.4f fl/Ns=%.4f | ring %.3f b/3 %.2e 1/L %.2e | aA %.3g aB %.3g aC %.3g aD %.3g aE %.3g aF %.3g aG %.3g | chkB %.1e  (%.0fs)'
          % (reg, J, mfl, row['bhat'], beta, row['lam'], row['floor_rel'], prof['ring'], prof['b/3'], prof['1/L'], atol['A'], atol['B'],
             atol['C'], atol['D'], atol['E'], atol['F'], atol['G'], row['chkB'], time.time() - t0), flush=True)
json.dump(rows, open(os.path.join(C.CACHE, 'struct_%d.json' % Dg), 'w'))
out = []
w = out.append
w('Delta_g = %d. Impact-parameter tolerance profile (weight tolerated at radius t, units of the ring weight) and structure tolerances' % Dg)
for reg in ('LJ', 'R'):
    rr = [r for r in rows if r['reg'] == reg]
    w('\nregion %s (%d states)' % (reg, len(rr)))
    for k in ('ring', '0.8b', 'b/3', 'b/10', '1/L'):
        v = [r['prof'][k] for r in rr]; w('  tol(t = %-5s): min %.3e  max %.3e' % (k, min(v), max(v)))
    for e in (0.0, 0.2, 0.5):
        v = [r['tau'][e] for r in rr]; w('  tau(eps=%.1f): min %.4f' % (e, min(v)))
    for k in 'ABCDEFG':
        v = [r['atol'][k] for r in rr]
        i = int(np.argmin(v)); r = rr[i]
        w('  a_tol[%s]: min %.4g (J=%d, m/Lam=%.3f, bhat=%.1f)   max %.4g' % (k, v[i], r['J'], r['mfl'], r['bhat'], max(v)))
txt = '\n'.join(out)
print(txt)
open(os.path.join(C.CACHE, 'struct_%d.txt' % Dg), 'w').write(txt + '\n')

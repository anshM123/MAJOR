"""
largespin/ls_curvbound.py (LARGESPIN_ANALYTIC.md §2): Lemma C' -- a sign-free bound on the leading-Regge AdS curvature of the hybrid
that uses the decay of Phi_k(x) = int psi_k cos(px) dp, so that it survives to large impact parameter.

Exact representation (Mehler-Dirichlet, Lemma-C substitution):  with sigma(phi) = 2 sinh(beta/2) sin(phi) (so Lam sigma = bhat sin phi)
and s = 2 asinh(sigma/2),
    curv_k = w_k (2/pi) int_0^{pi/2} [ Phi_k(Lam s)/cosh(s/2) - Phi_k(Lam sigma) ] dphi ,      w_2 = 2,  w_3 = u sigma_h.
Integration by parts (psi_k(0) = 0, psi_k has the factor (1-p)^5):
    |Phi_k(x)|  <= min(M0, C2/x^2),   M0 = int|psi|,   C2 = |psi'(0)| + |psi'(1)| + int|psi''|      [P]
    |Phi_k'(x)| <= min(M1, C3/x^2, C3b/x^3),   M1 = int p|psi|, C3 = int |g''|, C3b = |g''(0)| + |g''(1)| + int |g'''|, g = p psi  [P]
and 1 - sech(s/2) <= s^2/8, 0 <= sigma - s <= sigma^3/24, s >= sigma(1 - sigma^2/24). Hence, pointwise in phi,
    |...| <= min(M0, C2/(Lam s)^2) s^2/8 + Lam sigma^3/24 * min(M1, C3/(Lam s)^2)
which is non-decreasing in phi, so the right-endpoint Riemann sum is an upper bound B_k(beta).  [P]; constants [DV]
Statement: lam_hyb >= rig' := (N_flat/5 - B_2 - u B_3)/N_s  for every state, N_flat/5 certified (CERT_CLOSE).
Output: cache/curvbound_<Dg>.txt.  Usage: python3 ls_curvbound.py [Dg]
"""
import os, sys
sys.dont_write_bytecode = True
import numpy as np
import ls_core as C

Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
L = C.lam(Dg); F = C.NFw(Dg)


def consts(k):
    f = F.psi2 if k == 2 else F.psi3
    # polynomial coefficients: psi = (1-p)^r p X(p) -> exact polynomial via numpy poly
    P = np.polynomial.Polynomial
    base = P([1, -1])**F.r*P([0, 1])*P(F.X if k == 2 else F.Y)
    d1, d2 = base.deriv(1), base.deriv(2); q2 = (P([0, 1])*base).deriv(2)
    # piecewise Gauss on 400 equal panels (integrand |poly| has finitely many kinks; error << 1e-10 relative)
    edges = np.linspace(0, 1, 401); xs, ws = np.polynomial.legendre.leggauss(40)
    def absint(poly):
        tot = 0.0
        for a, b in zip(edges[:-1], edges[1:]):
            t = (xs + 1)/2*(b - a) + a; tot += (np.abs(poly(t))*ws).sum()*(b - a)/2
        return tot
    M0 = absint(base); M1 = absint(P([0, 1])*base)
    C2 = abs(d1(0.0)) + abs(d1(1.0)) + absint(d2); C3 = absint(q2)
    q3 = q2.deriv(1); C3b = abs(q2(0.0)) + abs(q2(1.0)) + absint(q3)
    chk = np.max(np.abs(base(np.linspace(0, 1, 101)) - f(np.linspace(0, 1, 101))))
    pad = 1 + 1e-6
    return dict(M0=M0*pad, M1=M1*pad, C2=C2*pad, C3=C3*pad, C3b=C3b*pad, chk=chk)


K = {2: consts(2), 3: consts(3)}

# numerical envelopes E0(x) >= sup_{y >= x}|Phi(y)|, E1(x) >= sup_{y >= x}|Phi'(y)| on a grid (step H) with a Lipschitz pad from the
# analytic decay bounds, and the analytic bounds beyond XMAX.  [DV: float evaluation of Phi on the grid; pads rigorous]
H, XMAX = 0.01, 6000.0
XG = np.arange(0.0, XMAX + H, H)


def envelopes(k):
    c = K[k]; F_ = F.psi2 if k == 2 else F.psi3
    p, w = C.gl01(4000); ps = F_(p)*w
    ph0 = np.empty_like(XG); ph1 = np.empty_like(XG)
    for i0 in range(0, len(XG), 2000):
        xx = XG[i0:i0 + 2000]; cs = np.outer(xx, p)
        ph0[i0:i0 + 2000] = np.cos(cs) @ ps; ph1[i0:i0 + 2000] = -(np.sin(cs) @ (p*ps))
    # Lipschitz pads on [x_i, x_i + H]: |Phi'| <= min(M1, C3/x^2, C3b/x^3), |Phi''| <= int p^2|psi| =: M2 (crude, fine for the pad)
    xs = np.maximum(XG, 1e-9)
    L1 = np.minimum(np.minimum(c['M1'], c['C3']/xs**2), c['C3b']/xs**3)
    M2 = float((np.abs(F_(p))*p*p*w).sum())*1.000001
    Pp = np.polynomial.Polynomial
    base = Pp([1, -1])**F.r*Pp([0, 1])*Pp(F.X if k == 2 else F.Y)
    g2 = (Pp([0, 0, 1])*base).deriv(2)
    e = np.linspace(0, 1, 401); xs_, ws_ = np.polynomial.legendre.leggauss(40)
    C4 = sum((np.abs(g2((xs_ + 1)/2*(b - a) + a))*ws_).sum()*(b - a)/2 for a, b in zip(e[:-1], e[1:]))*1.000001
    L2 = np.minimum(M2, C4/xs**2)                                      # |Phi''| <= int|(p^2 psi)''|/x^2  [P]
    a0 = np.abs(ph0) + H*L1; a1 = np.abs(ph1) + H*L2
    tail0 = c['C2']/XMAX**2; tail1 = min(c['C3']/XMAX**2, c['C3b']/XMAX**3)
    E0 = np.maximum.accumulate(np.append(a0, tail0)[::-1])[::-1]       # E0[i] = sup over grid cells j >= i and the tail
    E1 = np.maximum.accumulate(np.append(a1, tail1)[::-1])[::-1]
    return E0, E1, tail0, tail1


ENV = {k: envelopes(k) for k in (2, 3)}


def Ek(k, which, x):
    """sup_{y >= x} |Phi^{(which)}(y)| (upper bound)."""
    E0, E1, t0, t1 = ENV[k]; E = E0 if which == 0 else E1
    i = np.floor(np.asarray(x)/H).astype(int)                         # cell containing x: its left node <= x
    out = np.where(i < len(XG), E[np.minimum(i, len(XG))], (t0 if which == 0 else t1))
    return out


def Bk_env(k, beta, n=20000):
    """cell-wise upper sum: on [phi_a, phi_b], |...| <= E0(Lam s_a) s_b^2/8 + Lam sigma_b^3/24 E1(Lam s_a)."""
    ph = np.arange(0, n + 1)*(np.pi/2)/n
    sg = 2*np.sinh(beta/2)*np.sin(ph); s = sg*(1 - sg**2/24)
    sa, sb, sgb = s[:-1], 2*np.arcsinh(sg[1:]/2), sg[1:]
    t = Ek(k, 0, L*sa)*sb**2/8 + L*sgb**3/24*Ek(k, 1, L*sa)
    return (2/np.pi)*t.sum()*(np.pi/2)/n


def Bk(k, beta, n=20000):
    c = K[k]
    phi = np.arange(1, n + 1)*(np.pi/2)/n                       # right endpoints (integrand non-decreasing in phi)
    sg = 2*np.sinh(beta/2)*np.sin(phi)
    s = np.maximum(sg*(1 - sg**2/24), 1e-300)
    t1 = np.minimum(c['M0'], c['C2']/(L*s)**2)*s**2/8
    t2 = L*sg**3/24*np.minimum(np.minimum(c['M1'], c['C3']/(L*s)**2), c['C3b']/(L*s)**3)
    return (2/np.pi)*((t1 + t2).sum()*(np.pi/2)/n)


out = []
w = out.append
w('Delta_g = %d. Constants (padded 1e-6): %s' % (Dg, {k: {q: '%.6g' % v for q, v in K[k].items()} for k in K}))
w('rig\' = (N_flat/5 - 2 B_2(beta) - u B_3(beta))/N_s, worse helicity; computed exact curvature for comparison')
w('%6s %6s %8s %8s | %10s %10s %10s | %8s %8s' % ('m/Lam', 'J', 'bhat', 'beta', 'N_s', '|bound|', 'curv', 'rig\'', 'lam_hyb'))
res = {}
for mfl in (1.3, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0):
    Jmax = int(((mfl*L)**2/(Dg - 2) - (Dg - 2) - 1)/2)
    res[mfl] = []
    for b in (5, 10, 20, 30, 40, 50, 60, 70, 80, 91.2, 100, 120, 150, 200, 250, 300, 350, 400, 450, 500, 530, 600, 800, 1000, 1500, 1997):
        J = int(np.floor((b*mfl - 1)/2))
        if J < 4 or J > min(Jmax, 4000):
            continue
        u = 1/mfl**2
        worst = None
        for h in ('++', '+-'):
            if h == '++' and J % 2:
                continue
            st = C.hyb_fast(Dg, J, mfl, h)
            Nf = abs(st['B2']) + abs(st['B3'])
            bnd_a = 2*Bk(2, st['beta']) + u*Bk(3, st['beta'])
            bnd = min(bnd_a, 2*Bk_env(2, st['beta']) + u*Bk_env(3, st['beta']))
            rig = (Nf/5 - bnd)/st['Ns']
            r = (rig, st['Ns'], bnd, st['cv2'] + st['cv3'], st['base']/st['Ns'], st['bhat'], st['beta'], (Nf/5 - bnd_a)/st['Ns'])
            if worst is None or rig < worst[0]:
                worst = r
        rig, Ns, bnd, cv, lamh, bh, beta, riga = worst
        res[mfl].append((bh, rig))
        w('%6.2f %6d %8.2f %8.5f | %10.3e %10.3e %10.3e | %8.4f %8.4f | analytic-constant version %8.4f' % (mfl, J, bh, beta, Ns, bnd, cv, rig, lamh, riga))
        print(out[-1], flush=True)
w('\nlargest sampled bhat per m/Lam up to which rig\' > 0 on every sampled state (sampled at the listed bhat, floor(J)):')
for mfl, v in res.items():
    bm = 0
    for bh, rg in v:
        if rg <= 0:
            break
        bm = bh
    w('  m/Lam = %g: rig\' > 0 for sampled bhat <= %.1f' % (mfl, bm))
txt = '\n'.join(out)
print(txt)
open(os.path.join(C.CACHE, 'curvbound_%d.txt' % Dg), 'w').write(txt + '\n')

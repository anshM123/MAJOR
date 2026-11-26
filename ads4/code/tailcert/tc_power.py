"""
tailcert/tc_power.py (TAIL_CERT.md §2): certified upper bounds for the ladder transform of a pure power,
    S^(beta) := w^beta S(beta),  S(beta) = L((1-z)^{-beta}) = sum_k U_k (beta)_k/k!,
uniformly over every tail state: a >= 998, b in {0,4}, w = 1 - r <= e^{-L0}.

(I) Oscillatory bound (Theorem P, TAIL_CERT §2.2) for Re beta = sigma > 0, tau = Im beta in [tau1, tau2] (tau >= 0; conjugate otherwise):
  S(beta) = kappa(beta) [ E_V J_{lambda_V}(beta) - int_0^{q1} q^{beta-1} T + int_{q1}^oo q^{beta-1} g ],
  kappa = Gamma(N+s)/(Gamma(beta)Gamma(N+s-beta)), J_lam = int_0^{q1} q^{beta-1}((1-lam q)/(1+q))^a dq, V ~ Beta(N+1, n+s), lam = 1/V - 1.
  Path for J: q = i y (0 <= y <= rho), arc |q| = rho down to q1 = rho.  In x = 2 a q: Y1 = 2 a rho.
  Regime A (u2 := tau2/a <= 0.52):  rho = u2, Y1 = 2 tau2; arc exponent <= 0 analytically (f <= u(pi/2 - phi) - 2 rho cos(phi)/(1+rho^2) <= 0).
  Regime B (u2 > 0.52, hence a <= tau2/0.52):  rho = 0.92, arc exponent max certified on phi-boxes for the a-cell.
  |S^| <= e^{eps_G} { G Y1^sigma [ e^{c1}/sigma + (pi/2) e^{a Ehat + c2} ] + |1/Gamma(beta)| [Gamma(sigma, 0.995 Y1)/0.995^sigma + tiny] + NC }
  with G = sup_cell |1/Gamma(beta)| e^{-pi tau/2}, c1, c2 the lambda-corrections (delta_lam <= 7 sqrt(w)), NC <= 1e-30 (see doc).
(II) Near-real bound (Lemma NB, TAIL_CERT §2.1) for Re beta = sigma >= 1/4, |tau| <= 1/4:
  |S^(beta) - 1| <= C_tau (1 + 1e-6) [(sigma)_3/3 + (1+s/2)(sigma)_2 d + ... ] * 0.2501/a^2 + 1e-12,
  C_tau = exp(tau^2 (1/sigma^2 + 1/sigma)/2), and S^(beta) - 1 <= 0 for real beta (E_k <= 1).
"""
import os, sys, math, json
from flint import arb, acb, ctx
ctx.prec = 128
PI = arb.pi()
A_MIN = 998


def up(x):
    return float(arb(x).upper())


def G_cell(sigma, t1, t2):
    """sup over tau in [t1,t2] of |1/Gamma(sigma + i tau)| e^{-pi tau/2}; Arb enclosure on sub-boxes."""
    nb = max(1, int(math.ceil((t2 - t1)/0.05)))
    best = 0.0
    for k in range(nb):
        a_ = t1 + (t2 - t1)*k/nb; b_ = t1 + (t2 - t1)*(k + 1)/nb
        tau = arb((a_ + b_)/2, (b_ - a_)/2*1.0000001)
        v = abs(acb(arb(sigma), tau).rgamma())*(-PI*tau/2).exp()
        best = max(best, up(v))
    return best


def rg_cell(sigma, t1, t2):
    """sup |1/Gamma(sigma + i tau)| over the cell."""
    nb = max(1, int(math.ceil((t2 - t1)/0.05)))
    best = 0.0
    for k in range(nb):
        a_ = t1 + (t2 - t1)*k/nb; b_ = t1 + (t2 - t1)*(k + 1)/nb
        tau = arb((a_ + b_)/2, (b_ - a_)/2*1.0000001)
        v = abs(acb(arb(sigma), tau).rgamma())
        best = v if best == 0.0 or v > best else best
    return arb(best.upper())


def Ehat(u, rho, nphi=4000):
    """certified upper bound of max_{phi in [0,pi/2]} u(pi/2 - phi) + (1/2) log((1-2 rho c + rho^2)/(1+2 rho c + rho^2))."""
    u = arb(u); rho = arb(rho)
    best = -1e300
    for k in range(nphi):
        ph = PI/2*(arb(k) + arb(0.5, 0.5))/nphi
        c = ph.cos()
        f = u*(PI/2 - ph) + ((1 - 2*rho*c + rho*rho)/(1 + 2*rho*c + rho*rho)).log()/2
        best = max(best, up(f))
    return best


def Ehat_d(u, rho, delta, n1=300, n2=1500):
    """certified upper bound of max_phi u(pi/2-phi) + log((|1-q| + delta)/|1+q|), q = rho e^{i phi}, phi in [0, pi/2]
    (the lambda-perturbed arc exponent: |1 - lam q| <= |1-q| + |lam-1| rho, delta >= |lam-1| rho).  Two-level phi grid."""
    u = arb(u); rho = arb(rho); de = arb(delta)
    best = -1e300
    cut = PI/2 - arb('0.05')
    for (lo, hi, n) in ((arb(0), cut, n1), (cut, PI/2, n2)):
        for k in range(n):
            ph = lo + (hi - lo)*(arb(k) + arb(0.5, 0.5))/n
            c = ph.cos()
            f = u*(PI/2 - ph) + ((1 - 2*rho*c + rho*rho).sqrt() + de).log() - (1 + 2*rho*c + rho*rho).log()/2
            best = max(best, up(f))
    return best


def inc_gamma_upper(sigma, X):
    """Gamma(sigma, X) for sigma <= 1, X > 0: <= X^{sigma-1} e^{-X} (since t^{sigma-1} <= X^{sigma-1} on [X, oo)); for sigma > 1 use
    Arb's exact incomplete gamma (upper)."""
    s = arb(sigma); X = arb(X)
    if sigma <= 1:
        return arb((X**(s - 1)*(-X).exp()).upper())
    return arb(acb(X).gamma_upper(acb(s)).real.upper())


_ETAB = None
ETAB_U0, ETAB_DU, ETAB_N = 0.5, 0.001, 700            # u in [0.5, 1.2]  (tau <= 1146, a >= 998)
ETAB_L0 = float(os.environ.get('TC_L0TAB', '20'))
ETAB_DELTA = 7*math.exp(-ETAB_L0/2)*0.92               # >= |lam - 1| rho for every tail state with L >= ETAB_L0


def Ehat_table():
    """certified Ehat(u, 0.92) at u = u_i (grid); Ehat is nondecreasing in u (d/du = pi/2 - phi >= 0), so Ehat(u) <= table[ceil]."""
    global _ETAB
    if _ETAB is None:
        fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', f'ehat_d_table_L{ETAB_L0:g}.json')
        if os.path.exists(fn):
            _ETAB = json.load(open(fn))
        else:
            _ETAB = [Ehat_d(ETAB_U0 + ETAB_DU*i, 0.92, ETAB_DELTA) for i in range(ETAB_N)]
            json.dump(_ETAB, open(fn, 'w'))
    return _ETAB


def Ehat_up(u):
    T = Ehat_table()
    if u <= ETAB_U0:
        return T[0]
    i = int(math.ceil((u - ETAB_U0)/ETAB_DU))
    assert i < len(T), u
    return T[i]


def K_osc(sigma, t1, t2, L0, a_cells=None):
    """certified upper bound of |S^(sigma + i tau)| for tau in [t1, t2] (t2 >= t1 >= 0), all tail states (a >= 998, w <= e^{-L0}).
    Returns (bound, info)."""
    assert sigma > 0 and t2 >= t1 >= 0
    assert t2 <= 1146.0          # validity range of the NC / q >= 0.92 estimates (pi tau/2 <= 1800)
    sw = math.exp(-L0/2)
    dlam = 7*sw                                  # delta_lambda <= w + 6.6 sqrt(w) <= 7 sqrt(w)
    G = G_cell(sigma, t1, t2)
    RG = rg_cell(sigma, t1, t2)
    eG = math.exp(2*(t2 + 5)**2/(998*math.exp(L0)))     # e^{eps_Gamma}, x = N + s >= 998 e^{L0}
    tau2 = max(t2, 1e-9)
    res = []
    # ---- regime A: a >= tau2/0.52 (rho = tau2/a <= 0.52, Y1 = 2 tau2), uniform in a
    if True:
        Y1 = 2*tau2
        c1 = 3*dlam*Y1*Y1/(8*max(A_MIN, tau2/0.52))
        c2 = dlam*tau2/0.48
        X = 0.995*Y1
        bA = eG*(G*Y1**sigma*(math.exp(c1)/sigma + (math.pi/2)*math.exp(c2))
                 + up(RG*inc_gamma_upper(sigma, X))/0.995**sigma) + 1e-30
        res.append(('A', bA))
    # ---- regime B: 998 <= a < tau2/0.52, rho = 0.92
    amax = tau2/0.52
    if amax > A_MIN:
        rho = 0.92
        # a-cells, geometric
        edges = [A_MIN]
        while edges[-1] < amax:
            edges.append(min(amax, edges[-1]*1.02))
        for a1, a2 in zip(edges[:-1], edges[1:]):
            uhi = t2/a1
            E = Ehat_up(uhi)
            Y1 = 2*a2*rho
            c1 = 3*dlam*Y1*Y1/(8*a1)
            assert L0 >= ETAB_L0
            c2 = 0.0                      # the lambda-correction is inside the delta-table (Ehat_d)
            loss = a2*max(E, 0.0) if E > 0 else a1*E
            X = 0.995*2*a1*rho
            lg = math.log(eG) + math.log(G*Y1**sigma) + max(c1 - math.log(sigma), math.log(math.pi/2) + loss + c2) + math.log(2)
            lrt = math.log(eG) + float((RG*inc_gamma_upper(sigma, X)).log().upper()) - sigma*math.log(0.995)
            lB = max(lg, lrt) + math.log(2)
            bB = math.exp(lB) if lB < 700 else float('inf')
            res.append(('B', bB, a1, a2, E, lB))
    b = max(r[1] for r in res)
    return b, res


def logK_osc(sigma, t1, t2, L0):
    """log of the K_osc bound (no overflow)."""
    b, res = K_osc(sigma, t1, t2, L0)
    l = math.log(res[0][1])
    for r in res[1:]:
        l = max(l, r[5])
    return l


def NB_bound(sigma, tau, s=9):
    """|S^(beta) - 1| for beta = sigma + i tau, |tau| <= 1/4, sigma >= 1/4: Lemma NB (TAIL_CERT §2.1)."""
    sg = arb(sigma); ta = arb(abs(tau))
    Ct = (ta*ta*(1/(sg*sg) + 1/sg)/2).exp()
    p3 = sg*(sg + 1)*(sg + 2); p2 = sg*(sg + 1)
    # E[K^3/3 + s K^2/2] (c/d)^{-3} terms: (sigma)_3/3 + d(...) ; d <= 1e-8 for L0 >= 20, a-factor 0.2501/a^2
    main = (p3/3 + arb('1e-6')*(3*p2 + sg + s*(p2 + sg)))*arb('0.2501')/arb(A_MIN)**2
    return up(Ct*main*arb('1.0001') + arb('1e-12'))


if __name__ == '__main__':
    for (sg, t1, t2) in ((0.25, 0.25, 0.3), (0.25, 1, 1.05), (0.25, 10, 10.05), (0.25, 44.95, 45), (0.25, 300, 300.05),
                         (0.25, 600, 600.05), (0.25, 998, 998.05), (0.25, 1053.95, 1054), (0.25, 1099.95, 1100), (1.25, 10, 10.05),
                         (3.5, 10, 10.05), (3.5, 998.95, 999)):
        b, res = K_osc(sg, t1, t2, 20)
        print(f'sigma={sg} tau in [{t1},{t2}]: |S^| <= {b:.4g}   regimes: ' + ', '.join(f'{r[0]}:{r[1]:.3g}' for r in res[:1]) +
              (f' ... B max {max(r[1] for r in res if r[0]=="B"):.3g}' if any(r[0] == 'B' for r in res) else ''), flush=True)
    for sg, ta in ((2.5, 0), (1.75, 0), (1.75, 0.25), (0.25, 0.25)):
        print(f'NB: sigma={sg}, tau={ta}: |S^ - 1| <= {NB_bound(sg, ta):.4e}')

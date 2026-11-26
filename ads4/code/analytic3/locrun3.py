"""
locrun3.py La Lb [J]: E_rho (k = 2, dressed, p = 3; k = 3, p = 4) for the L-cell [La, Lb] by Lemma Loc''' (loc5.py), nu < 1150.
  E_rho2 = d^3 (2/Lam) 2 int_0^{1150} |f| |V_b L(R_J)| dnu * 2,   |f| := |psi2/k0| + |r_theta| <= FSUP = 21.34 on [0, 1300] (TAIL_CERT §4,
           EPSW: sup|psi| + sup|r| + sup|f_theta - f|);
  E_rho3 = d^3 tau_T Lam^2 int_0^1 |psi3(p)| |L(R'_J)(Lam p)| dp * 2 * 2,  tau_T <= d^2/(3.9 a^2) (D d^2 >= 3.9 a^2 for L >= 3.01),
           sup|psi3| <= sum|coef|  (same prefactors as rt2_loc.py).
nu >= 1150: not split (rt2_resid far term, full kernel).  Output cache/loc3_La_Lb_J.json."""
import sys, os, math, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import loc5 as L5
import rt_cert as RC
FSUP = 21.33815


def logsum(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


def run(La, Lb, J):
    t0 = time.time()
    d = math.exp(-La)
    th = L5.t_range(La, Lb)[1]
    cD = min(3.9, (2 - th)**2/(1 + 9*th/998)*(1 - 1e-9))     # D d^2/a^2 = (2 - t + e)^2/(1 + e) >= (2 - t_hi)^2/(1 + e_max)
    G2 = L5.geometry(La, Lb, True); G3 = L5.geometry(La, Lb, False)
    ls2 = []; worst = (-1e9, None)
    b = L5.cell_bound(G2, 3, J, 0.0, 0.4, True); ls2.append(math.log(0.4*FSUP) + b[0])
    edges = [0.4, 1.0] + [float(k) for k in range(2, 1151)]
    for n1, n2 in zip(edges[:-1], edges[1:]):
        b = L5.cell_bound(G2, 3, J, n1, n2, True)
        v = math.log(FSUP*(n2 - n1)) + b[0]; ls2.append(v)
        if v > worst[0]: worst = (v, n1)
    lI2 = logsum(ls2)
    ls3 = []
    b = L5.cell_bound(G3, 4, J, 0.0, 0.4, False); ls3.append(math.log(0.4/RC.LAMf) + b[0])
    lo = 0.4
    while lo < RC.LAMf:
        hi = min(RC.LAMup, math.floor(lo) + 1.0)
        b = L5.cell_bound(G3, 4, J, lo, hi, False)
        ls3.append(math.log((hi - lo)/RC.LAMf) + b[0]); lo = hi
    lI3 = logsum(ls3)
    s3 = sum(abs(float(c)) for c in RC.P3)
    lE2 = 3*math.log(d) + math.log(2/RC.LAMf) + math.log(2) + lI2 + math.log(2)
    lE3 = 3*math.log(d) + math.log(d*d/(cD*998**2)) + 2*math.log(RC.LAMup) + math.log(s3) + math.log(2) + lI3 + math.log(2)
    out = dict(La=La, Lb=Lb, J=J, lI2=lI2, lI3=lI3, lE2=lE2 + 1e-6, lE3=lE3 + 1e-6, worst2=worst, t=[G2['t_lo'], G2['t_hi']],
               K0=L5.K0_bound(), secs=time.time() - t0)
    json.dump(out, open(os.path.join(HERE, 'cache', f'loc3_{La}_{Lb}_{J}.json'), 'w'))
    print(out, flush=True)
    return out


if __name__ == '__main__':
    La, Lb = float(sys.argv[1]), float(sys.argv[2]); J = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
    run(La, Lb, J)

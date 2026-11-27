"""P1: worst-case |F(0)| from Matsubara data |F(n)| <= eps, n = 1..Nmax, in three a-priori classes.

All three are linear programs over a discretised measure (scipy/HiGHS). Units: 2*pi*T = 1, so s = n are
the Matsubara points and F(s) is the Laplace transform of the relaxation-function deviation.

Class A ("H_gamma / weighted-TV Laplace"):  F(s) = int x^s dtau(x), x = e^{-t} in (0,1],
          int x^{-gamma} d|tau| <= 1.   Then |F| <= 1 on Re s >= -gamma (a sub-class of Theorem C's (S)).
          Prediction (Theorem C + Muntz): W(eps) ~ (log 1/eps)^{-2 gamma}.
Class B ("completely monotone relaxation / Stieltjes, signed, bounded TV"):
          F(s) = int dsigma(lam)/(lam + s), supp sigma in [gamma, inf), int d|sigma| <= 1.
          This is E - E_ref with both relaxation functions CM (positive mixtures of decaying exponentials)
          and total masses Phi(0), Phi_ref(0) bounded (a sum-rule quantity).
          Prediction (this README, Sec. P1.4): W(eps) ~ exp(-c sqrt(gamma log 1/eps)).
Class C ("one-signed deviation"): as A but tau >= 0.  Lyapunov: W = eps^{gamma/(1+gamma)} exactly
          (single atom), already from n = 1.

Output: out_p1_rates.txt
"""
import numpy as np
from scipy.optimize import linprog

NMAX = 60
EPS = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]


def solve(phi_data, phi_target, weight, eps, signed=True):
    """max sum target*tau  s.t. |sum data_n*tau| <= eps, sum weight*|tau| <= 1."""
    K = len(phi_target)
    if signed:
        # tau = p - q
        c = -np.concatenate([phi_target, -phi_target])
        A_tv = np.concatenate([weight, weight])[None, :]
        A_d = np.hstack([phi_data, -phi_data])
    else:
        c = -phi_target
        A_tv = weight[None, :]
        A_d = phi_data
    A = np.vstack([A_tv, A_d, -A_d])
    b = np.concatenate([[1.0], eps * np.ones(2 * phi_data.shape[0])])
    res = linprog(c, A_ub=A, b_ub=b, bounds=(0, None), method="highs",
                  options=dict(primal_feasibility_tolerance=1e-10, dual_feasibility_tolerance=1e-10))
    assert res.status == 0, res.message
    return -res.fun


def class_A(gamma, eps, signed=True):
    x = np.concatenate([np.logspace(-25, -1, 2500), np.linspace(0.1, 1, 500)])
    n = np.arange(1, NMAX + 1)[:, None]
    # variable u = x^{-gamma} tau (unit TV weight) for conditioning
    data = x[None, :] ** (n + gamma)
    return solve(data, x ** gamma, np.ones_like(x), eps, signed)


def class_B(gamma, eps):
    lam = gamma * np.logspace(0, 7, 3000)
    n = np.arange(1, NMAX + 1)[:, None]
    data = 1.0 / (lam[None, :] + n)
    return solve(data, 1.0 / lam, np.ones_like(lam), eps, True)


if __name__ == "__main__":
    out = []
    for gamma in (0.5, 1.0):
        out.append(f"gamma = {gamma}")
        out.append(f"{'eps':>8} {'L=ln1/eps':>9} | {'A: W':>10} {'W*L^2g':>8} | {'B: W':>10} "
                   f"{'-lnW/sqrt(gL)':>13} | {'C: W':>10} {'eps^(g/(1+g))':>13}")
        for eps in EPS:
            L = np.log(1 / eps)
            wa = class_A(gamma, eps)
            wb = class_B(gamma, eps)
            wc = class_A(gamma, eps, signed=False)
            out.append(f"{eps:8.0e} {L:9.2f} | {wa:10.3e} {wa * L ** (2 * gamma):8.3f} | {wb:10.3e} "
                       f"{-np.log(wb) / np.sqrt(gamma * L):13.3f} | {wc:10.3e} {eps ** (gamma / (1 + gamma)):13.3e}")
        out.append("")
    txt = "\n".join(out)
    print(txt)
    open(__file__.replace("p1_rates_lp.py", "out_p1_rates.txt"), "w").write(txt + "\n")

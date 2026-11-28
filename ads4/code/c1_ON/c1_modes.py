"""Method B, n != 0 modes (double precision, fastflat).  Output per n:
  d0_n = log(8 w Pi0_n(0)),  I_n = int_0^inf D_n du,  B_n = d0_n/3 + I_n,  A_n = int_0^inf delta0_n du.
"""
import sys, time, numpy as np
import fastflat as ff
x, wgl = np.polynomial.legendre.leggauss(48)

C_T = -128*0.229575957618702205  # eps -> C_T (w^2 - p^2/3)/p^5  (OPE, check_flat_deriv.py)

def asym(u, w2):
    p2 = u + w2
    e = C_T*((2/3)*w2*p2**-2.5 - (1/3)*p2**-1.5)
    e1 = C_T*((2/3)*w2*(-2.5)*p2**-3.5 - (1/3)*(-1.5)*p2**-2.5)
    e2 = C_T*((2/3)*w2*(8.75)*p2**-4.5 - (1/3)*(3.75)*p2**-3.5)
    D = -(u/3)*(-e1/p2 + e2)
    return e, D

def tails(Um, w2):
    # integrate the asymptotic forms from Um to infinity: u = Um/t^2, t in (0,1]
    t = 0.5*(x+1); wt = 0.5*wgl
    u = Um/t**2; jac = 2*Um/t**3
    e, D = asym(u, w2)
    return (e*jac) @ wt, (D*jac) @ wt

def nodes(n, Umax=None):
    w2 = (2*np.pi*n)**2
    Umax = 1e4*max(w2, 1e4) if Umax is None else Umax
    ed = [0.0, w2/64, w2/16, w2/4]
    while ed[-1] < Umax:
        ed.append(ed[-1]*2)
    ed = np.array(ed)
    U = np.concatenate([0.5*(b-a)*x+0.5*(a+b) for a, b in zip(ed[:-1], ed[1:])])
    W = np.concatenate([0.5*(b-a)*wgl for a, b in zip(ed[:-1], ed[1:])])
    return U, W, ed[-1]

def mode(n):
    w2 = (2*np.pi*n)**2
    U, W, Um = nodes(n)
    P, P1, P2 = ff.Pi0(n, U)
    p2 = U + w2
    D = -(U/3)*(P2/P - 3/(4*p2**2))
    th = ff.pi_th(n, U)[0]
    eps = -(2/np.pi)*np.arctan(2*ff.M0/np.sqrt(p2)) + 8*np.sqrt(p2)*th   # = 8 p Pi0 - 1 without cancellation
    dl = np.log1p(eps)
    # tails: D ~ a u^{-5/2}, dl ~ b u^{-3/2} (fit from the largest nodes)
    eA, eD = asym(U[-1], w2)
    a = D[-1]/eD; b = dl[-1]/eA      # ratios to the asymptotic forms at the last node (should be ~1)
    tA, tD = tails(Um, w2)
    I = D @ W + tD
    A = dl @ W + tA
    P00 = ff.Pi0(n, np.array([0.0]))[0][0]
    d0 = np.log(8*np.sqrt(w2)*P00)
    return d0, I, d0/3 + I, A, a, b

if __name__ == '__main__':
    n0, n1 = int(sys.argv[1]), int(sys.argv[2])
    for n in range(n0, n1+1):
        t = time.time()
        r = mode(n)
        print(n, ' '.join(f'{v:.15e}' for v in r), f'{time.time()-t:.1f}s', flush=True)

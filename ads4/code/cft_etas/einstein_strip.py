"""
Einstein strip data for the Matsubara viscosity function
   E(s) = [dG(0) - dG(i*1.5*s)]/(1.5*s)      (units r_+=1: 2 pi T = 1.5, p = 1, eta_KSS = s_ent/4pi = 1),
so E(0+) = 4 pi eta/s = 1 and E(n), n = 1,2,..., are Euclidean (Matsubara) data.  omega = i*1.5*s.
(i) lowest zero-momentum shear QNM (zero of A(omega)) -> relaxation gap gamma_E = -Im(omega_QNM)/(2 pi T);
(ii) sup_t |E(-g + i t)| on lines Re s = -g (by the maximum principle this is sup over Re s >= -g).
"""
import numpy as np
import mpmath as mp
from einstein_ref import G

def A_of(w):
    return G(w)[1]

def E(s):
    s = complex(s)
    w = 1j*1.5*s
    g, _ = G(w)
    dG = g - 1j*w**3
    return -dG/(1.5*s)        # dG(0) = 0 in this contact scheme (checked in out_einstein_ref.txt)

if __name__ == '__main__':
    # (i) QNM: Newton/secant on A(omega)=0 starting from literature-like guess
    f = lambda w: complex(A_of(complex(w)))
    w0 = mp.findroot(lambda w: f(complex(w)), mp.mpc(1.85, -2.66), solver='secant', tol=1e-20)
    print('lowest k=0 shear QNM omega (r_+=1) =', complex(w0), '  -> omega/(2 pi T) =', complex(w0)/1.5)
    gE = -complex(w0).imag/1.5
    print('relaxation gap gamma_E = -Im omega/(2 pi T) = %.6f' % gE)
    w1 = mp.findroot(lambda w: f(complex(w)), mp.mpc(3.3, -4.9), solver='secant', tol=1e-14, verify=False)
    print('second QNM omega/(2piT) =', complex(w1)/1.5)
    # check E on the Euclidean axis
    print('E(n) n=1..4:', [round(E(n).real, 10) for n in range(1, 5)], ' E(1e-6)=', E(1e-6))
    # (ii) sup on lines
    for g in (0.25, 0.5, 0.75, 0.9, 1.1, 1.25, 1.5):
        ts = np.concatenate([np.linspace(0, 6, 121), np.linspace(6.5, 30, 48)])
        vals = [abs(E(-g + 1j*t)) for t in ts]
        i = int(np.argmax(vals))
        print('gamma=%.2f  sup_t |E(-gamma+it)| ~ %.5f at t=%.2f ;  |E| at t=30: %.2e' % (g, vals[i], ts[i], vals[-1]))

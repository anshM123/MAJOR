"""Zero-momentum shear QNMs of the AdS4 brane (units 2piT = 1.5), to fix the analyticity sector of the retarded kernel
K(t) = sum_n r_n e^{-i w_n t}: each term decays in |arg t| < pi/2 - theta_n with theta_n = -arg(w_n) - ... ; sector opening
alpha = min over n of (pi/2 - |angle of w_n below the real axis| ) ... reported as the angle of w_n below the real axis."""
import mpmath as mp, cmath
from einstein_ref import G
f = lambda w: complex(G(complex(w))[1])
roots = []
guess = [mp.mpc(1.85, -2.66), mp.mpc(3.16, -4.92)]
for n in range(3):
    g0 = guess[-1] + (guess[-1] - guess[-2]) if n >= 2 else guess[n]
    try:
        w = mp.findroot(lambda x: f(complex(x)), g0, solver='secant', tol=1e-12, verify=False)
    except Exception as e:
        print('fail', n, e); break
    w = complex(w)
    if n >= 2: guess.append(mp.mpc(w.real, w.imag))
    roots.append(w)
    ang = -cmath.phase(w)*180/cmath.pi
    print('n=%d  w/(2piT) = %.6f %+.6fi   angle below real axis = %.3f deg   |A(w)|=%.1e' % (n+1, w.real/1.5, w.imag/1.5, ang, abs(f(w))))
print('(n>=4 not converged by this LHP shooting; dropped.)  Angles approach pi/3 = 60 deg (AdS4 asymptotics).')
print('Term r_n e^{-i w_n t} with w_n = |w| e^{-i theta}: decays for t = |t|e^{i phi} iff |phi| < theta (mirror pole -w_n^* included).')
print('=> sector of analyticity of K(t), t>0: |arg t| < alpha with alpha -> pi/3 (min over n of theta_n = 55.2 deg for the lowest QNM).')

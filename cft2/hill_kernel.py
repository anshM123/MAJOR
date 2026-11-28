"""t_p from the stress tensor alone. For a chiral half with uniformiser f and L = -(1/2){f,x} (so psi'' = L psi,
f = psi1/psi2), the frame-horocyclic chiral length is l(x,y) = 2 log(K(x,y)/eps), where
    K(x,y) = (f(y)-f(x))/sqrt(f'(x) f'(y)) = psi(y),  psi'' = L psi,  psi(x)=0, psi'(x)=1.
So t_p^chiral = log[K(p-,p) K(p,p+)/(eps K(p-,p+))] needs only <T> along the arc p- -> p+ (on the universal cover).
Check numerically for a boundary-graviton-dressed BTZ uniformiser."""
import numpy as np
from scipy.integrate import solve_ivp
a = 1.7
f   = lambda x: np.exp(a*(x + 0.3*np.sin(x)))
fp  = lambda x: a*(1 + 0.3*np.cos(x))*f(x)
def L(x):   # -(1/2) Schwarzian of f, computed analytically from g = a(x + 0.3 sin x): {e^g} = {g} - g'^2/2
    g1, g2, g3 = a*(1 + 0.3*np.cos(x)), -0.3*a*np.sin(x), -0.3*a*np.cos(x)
    schw = g3/g1 - 1.5*(g2/g1)**2 - 0.5*g1**2
    return -0.5*schw
def K_hill(x, y):
    sol = solve_ivp(lambda t, u: [u[1], L(t)*u[0]], (x, y), [0.0, 1.0], rtol=1e-12, atol=1e-14)
    return sol.y[0, -1]
K_unif = lambda x, y: (f(y) - f(x))/np.sqrt(fp(x)*fp(y))
worst = 0
for x, y in [(0.1, 0.9), (-2.0, 1.5), (0.3, 0.3 + 2*np.pi), (1.0, 5.0)]:
    worst = max(worst, abs(K_hill(x, y)/K_unif(x, y) - 1))
    print(f"K_hill({x:.2f},{y:.2f}) = {K_hill(x, y):.12e}   K_unif = {K_unif(x, y):.12e}")
print("max relative difference:", f"{worst:.1e}")

"""Min over real w of the Einstein total normalised viscosity 4 pi eta_E(w)/s = deta(w) + w^2 (units p=1, eta(0)=1, vacuum rho_0 = w^3)."""
import numpy as np
from einstein_ref import deta
ws = np.concatenate([np.linspace(1e-3, 3, 61), np.linspace(3.2, 12, 23)])
vals = np.array([deta(w) + w**2 for w in ws])
i = int(np.argmin(vals))
print('min_w [deta+w^2] on grid = %.6f at w = %.3f (w in units r+=1; 2piT=1.5)' % (vals[i], ws[i]))
print('values at w = 0.001, 0.5, 1, 2, 3:', [round(float(deta(w) + w**2), 6) for w in (1e-3, 0.5, 1, 2, 3)])

"""rv_theta.py: the theta-cutoff bookkeeping above Lam (ANALYTIC_SLIVER2 §6, G1).
The text argues |f_theta| <= |f| pointwise for nu >= Lam because "0 <= theta <= 1 on the real axis".  But f_theta = 2 Re[theta(nu - i) g(nu)]
(TAIL_CERT R.2, tc_cutoff.py): theta is evaluated at nu - i (complex), cutoff centre nu_c = Lam + 100 ~ 1098.5.  We (i) evaluate
theta(nu - i) on [Lam, 1150] to show it is not real / can exceed 1 in modulus, (ii) bound the omitted term
  2 int_{nu >= Lam} |f_theta - f| |K_J|  <=  2 [ Ksup[Lam,1059) * dE  +  Ksup[1059,1150] * (tail_theta(1059) + m_EPSW(1059)) ]
with dE = int|f_theta - f| <= 1.66e-23 (tc_cutoff) and tail_theta(1059) = CONTOUR's bound of int_1059^oo |f_theta| (cb_lines.tail)."""
import sys, os, math, json
sys.path.insert(0, '/home/user/MAJOR/ads4/code/contour')
import cb_common as KC, cb_lines as CL
from flint import arb, acb
W = KC.load()
Lam = W.Lam; nuc = Lam + 100
print('nu_c =', nuc)
# theta as used by CONTOUR: inspect the definition
import inspect
src = inspect.getsource(KC)
i = src.find('def theta')
print(src[i:i + 600] if i >= 0 else 'theta def not found in cb_common')
# (i) theta(nu - i) on [Lam, 1150]
for nu in (1000.0, 1030.0, 1059.0, 1080.0, 1098.5, 1110.0, 1150.0):
    z = acb(nu, -1)
    th = KC.theta_series(W, z)
    print(f'nu = {nu}: theta(nu - i) = {complex(th.real.mid(), th.imag.mid()):.6g}, |theta| = {float(abs(th).mid()):.6f}')
# (ii) omitted term with the cached sup|K_J| of the review rerun (or the original cache)
CL_t, _ = CL.tail(W, 0.0, 1059)
print('CONTOUR bound int_1059^oo |f_theta| <=', CL_t)
for cell in ('3.65_3.75', '3.75_4.0', '4.0_4.5', '4.5_5.0'):
    for base in ('/home/user/MAJOR/ads4/code/review_s2/ws/code/analytic2/cache', '/home/user/MAJOR/ads4/code/analytic2/cache'):
        rows = []
        for fn in os.listdir(base):
            if fn.startswith(f'resid_{cell}_p') and 'of' in fn:
                rows += json.load(open(os.path.join(base, fn)))
        if rows: break
    if not rows:
        continue
    res = json.load(open(os.path.join(base, f'resid_{cell}.json'))) if os.path.exists(os.path.join(base, f'resid_{cell}.json')) else None
    Lam_f = float(Lam.mid())
    K_mid = max([r['K'] for r in rows if Lam_f - 1 <= r['lo'] < 1059.0] + [0.0])
    K_tail = max([r['K'] for r in rows if r['kind'] == 'tail'] + [0.0])
    La = float(cell.split('_')[0]); d = math.exp(-La)
    F = res['F'] if res else float('nan')
    pref = d**3*(2/Lam_f)
    om = 2*(K_mid*1.66e-23 + K_tail*(float(CL_t.upper()) + 5.896e-247))
    print(f'[{cell}] ({base.split("/")[-3]}): sup K on [Lam,1059) = {K_mid:.3e}, on [1059,1150] = {K_tail:.3e};'
          f' omitted 2 int_(nu>=Lam)|f_theta - f||K_J| <= {om:.3e} -> /F = {pref*om/F:.3e}')

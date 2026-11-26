"""cert3.py La Lb [k0] [J]: assemble Route T''' for the L-cell [La, Lb] (ANALYTIC_SLIVER3 §5):
   rho = (E2 + E3 + E_r + E_rho2 + E_rho3)/F_cert(beta-cell) * (1 + 1e-9);  PASS iff rho < 1.
Inputs: cache/contour3_La_Lb_k{k0}.json (contour3.py), cache/resid3_La_Lb.json (resid3.py ... assemble), cache/loc3_La_Lb_J.json (locrun3.py)."""
import sys, os, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kernel3
import rt_cert as RC
La, Lb = float(sys.argv[1]), float(sys.argv[2]); k0 = int(sys.argv[3]) if len(sys.argv) > 3 else 4; J = int(sys.argv[4]) if len(sys.argv) > 4 else 1500
C = lambda f: json.load(open(os.path.join(HERE, 'cache', f)))
ct = C(f'contour3_{La}_{Lb}_k{k0}.json'); rs = C(f'resid3_{La}_{Lb}.json'); lc = C(f'loc3_{La}_{Lb}_{J}.json')
F = RC.fcert_min(RC.beta_of_L(La), RC.beta_of_L(Lb))
assert abs(F - rs['F']) < 1e-12*F and abs(F - ct['F']) < 1e-12*F
E2 = sum(v for k, v in ct['psi2'].items() if k != 'Y'); E3 = sum(v for k, v in ct['psi3'].items() if k != 'Y')
Er = rs['Er_over_F']*F
Erho = math.exp(max(lc['lE2'], -700)) + math.exp(max(lc['lE3'], -700))
rho = (E2 + E3 + Er + Erho)/F*(1 + 1e-9)
lines = [f"Route T''' certificate, L in [{La}, {Lb}] (beta_c in [{RC.beta_of_L(La):.3f}, {RC.beta_of_L(Lb):.3f}]), a in [998, {kernel3.A_STAR:.0f}]: F_cert >= {F:.4e}",
         f"   psi2: " + '  '.join(f"{k} {ct['psi2'][k]/F:.2e}" for k in ('bottom', 'left', 'top', 'right')) + f"  (Y = {ct['psi2']['Y']})",
         f"   psi3: " + '  '.join(f"{k} {ct['psi3'][k]/F:.2e}" for k in ('bottom', 'left', 'top', 'right')),
         f"   residual E_r/F = {rs['Er_over_F']:.4e}  (small-nu {rs['small'][0]:.1e}+{rs['small'][1]:.1e}, G1 {rs['G1']:.1e}, far nu>=1150 {rs['far']:.1e});",
         f"      bands: " + ', '.join(f"{k} {v:.3e}" for k, v in sorted(rs['bands'].items())),
         f"      worst cells: " + ', '.join(f"[{w[1]},{w[2]}] {w[0]:.2e}" for w in rs['worst'][:5]),
         f"   Loc''' remainders (J = {J}): log E_rho2 = {lc['lE2']:.1f}, log E_rho3 = {lc['lE3']:.1f}",
         f"   rho <= {rho:.4e}   {'PASS' if rho < 1 else 'FAIL'}"]
print('\n'.join(lines))
open(os.path.join(HERE, 'cache', f'cert3_{La}_{Lb}.txt'), 'w').write('\n'.join(lines) + '\n')

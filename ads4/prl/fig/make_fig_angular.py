"""make_fig_angular.py -- Fig. 1 of letter_v3: bounds on tau = sqrt(t4^2 + alpha_T^2) against the higher-spin gap g (log-log).

Curves (all at leading order in 1/c_T; see the caption of letter_v3):
  * collider bound tau <= 4 (every unitary CFT_3);
  * tau_*(g) = 160/(9(g-2)) [1 - 49/(2g+9)^2]^(-1/2): the earlier positivity-free angular bound (no large constants);
  * the new positivity-free bound, exact finite-g formula
        tau <= 80 sqrt(3) / ( sqrt(A) n sqrt(n^2-1) ),  n = floor(g/2),  A = 3(sqrt2-1) n^2,
    (Racah dual, heavy feasibility proved, light single traces included; code/angular_racah_proof);
  * the LP limit of the same rules, 4 sqrt(21600)/g^3 = 587.9/g^3 (three-state no-go, dashed);
  * the dispersive functional at its own gap, tau0(Lam) = (240/Lam^4) sqrt(3023 log(2 Lam/3) - 5307), Lam^2 = (g-2)(g-1)
    (dotted; positivity certified only at g = 1000, where it gives the filled circle 2.89e-8; elsewhere it is the expected
    sqrt(log)/g^4 scaling with positivity assumed, shown for orientation only).
Usage (from ads4/prl/fig): python3 make_fig_angular.py -> fig_angular.pdf, fig_angular_preview.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
                     'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False, 'font.size': 8,
                     'axes.linewidth': 0.6, 'pdf.fonttype': 42})
INK, INK2 = '#1a1a1a', '#555555'
C_NEW, C_OLD, C_DISP, C_LP = '#2a78d6', '#3a9b5c', '#c0392b', '#7f7f7f'


def tau_new(g):
    n = np.floor(g / 2.0)
    A = 3 * (np.sqrt(2) - 1) * n ** 2
    return 80 * np.sqrt(3) / (np.sqrt(A) * n * np.sqrt(n ** 2 - 1))


def tau_star(g):
    return 160.0 / (9 * (g - 2)) / np.sqrt(1 - 49.0 / (2 * g + 9) ** 2)


def tau_disp(g):
    lam = np.sqrt((g - 2) * (g - 1))
    return 240.0 / lam ** 4 * np.sqrt(3023 * np.log(2 * lam / 3) - 5307)


fig, ax = plt.subplots(figsize=(3.375, 2.55))
fig.subplots_adjust(left=0.16, right=0.97, bottom=0.15, top=0.97)
g = np.logspace(np.log10(9.0), 5, 6000)
ax.axhline(4, color=INK2, lw=0.7, ls=(0, (5, 2, 1, 2)))
ax.plot(g, tau_star(g), color=C_OLD, lw=1.1, ls=(0, (1, 1.4)))
ax.plot(g, tau_new(g), color=C_NEW, lw=1.5)
ax.plot(g, 587.88 / g ** 3, color=C_LP, lw=0.9, ls=(0, (4, 2)))
gd = np.logspace(np.log10(120.0), 5, 400)
ax.plot(gd, tau_disp(gd), color=C_DISP, lw=0.9, ls=(0, (2, 1.5)), alpha=0.75)
ax.plot([1000], [2.8914e-8], 'o', ms=4, color=C_DISP, zorder=5)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(9, 1e5); ax.set_ylim(3e-14, 30)
ax.set_xlabel(r'higher-spin gap $\Delta_g$', labelpad=1)
ax.set_ylabel(r'bound on $\tau=(t_4^2+\alpha_T^2)^{1/2}$', labelpad=1)
ax.text(1.2e4, 5.6, r'collider, $\tau\leq4$', fontsize=7, color=INK2, ha='center')
ax.text(2.2e3, 2.3e-2, r'$\tau_*\propto\Delta_g^{-1}$', fontsize=7, color=C_OLD, ha='left')
ax.text(1.6e3, 2.5e-6, r'this work: $994/\Delta_g^{3}$', fontsize=7, color=C_NEW, ha='left')
ax.annotate(r'LP limit $588/\Delta_g^{3}$', xy=(2.0e4, 587.88/2.0e4**3), xytext=(1.5e4, 8e-9), fontsize=7, color=C_LP,
            ha='center', arrowprops=dict(arrowstyle='-', lw=0.5, color=C_LP))
ax.annotate(r'$\propto\sqrt{\log\Delta_g}/\Delta_g^{4}$' + '\n(positivity assumed)', xy=(5e3, tau_disp(5e3)), xytext=(4.5e2, 2e-12),
            fontsize=6.5, color=C_DISP, ha='center', arrowprops=dict(arrowstyle='-', lw=0.5, color=C_DISP))
ax.annotate('certified\n(positivity)', xy=(1000, 2.8914e-8), xytext=(140, 4e-10), fontsize=6.5, color=C_DISP,
            ha='center', arrowprops=dict(arrowstyle='-', lw=0.5, color=C_DISP))
ax.tick_params(which='both', direction='in', labelsize=7)
fig.savefig(os.path.join(HERE, 'fig_angular.pdf'))
fig.savefig(os.path.join(HERE, 'fig_angular_preview.png'), dpi=220)
print('tau_new(1000)=%.4e  tau_star(1000)=%.4e  tau_disp(1000)=%.4e' % (tau_new(1000.0), tau_star(1000.0), tau_disp(1000.0)))

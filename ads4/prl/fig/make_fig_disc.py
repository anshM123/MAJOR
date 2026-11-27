#!/usr/bin/env python3
"""Fig. 1 of letter_v2: the universal disc of thermal spectral weights.

Axes: the shear sum rule S_E (horizontal) and the Hall sum rule S_H (vertical), both normalized by the energy density.
  S_E - 3/8 = t4/32 (Chowdhury-David-Prakash; checked on free fields, the AdS4 brane and large-N CS-fermion),
  |S_H| = |alpha_T|/32 (this work; checked on the AdS4 brane and on the CS-fermion OPE side).
Every unitary CFT3 obeying the sum-rule conditions has t4^2 + alpha_T^2 <= 16, hence lies in the disc of radius 1/8 about
(3/8, 0). Free scalar and free fermion sit at the ends of the horizontal diameter; large-N Chern-Simons-matter theories trace
the circle; Einstein gravity is the centre. A higher-spin gap Delta_g shrinks the disc to radius tau_*(Delta_g)/32
(no positivity input), and to 9.0e-10 at Delta_g >= 1000 with the certified bound tau <= 2.8914e-8.
Usage (from ads4/prl/fig): python3 make_fig_disc.py -> fig_disc.pdf, fig_disc_preview.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2 = '#1a1a1a', '#555555'
C1, C2 = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
                     'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False, 'font.size': 7.5,
                     'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'axes.linewidth': 0.6, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'text.color': INK, 'axes.labelcolor': INK, 'pdf.fonttype': 42})
FS = 7.0


def tau_star(dg):
    return 160.0 / (9.0 * (dg - 2.0)) / np.sqrt(1.0 - 49.0 / (2.0 * dg + 9.0) ** 2)


fig, ax = plt.subplots(figsize=(3.375, 2.75))
fig.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.97)
c0 = 3.0 / 8.0
R = 1.0 / 8.0
ax.set_aspect('equal')
ax.set_xlim(c0 - 0.162, c0 + 0.162)
ax.set_ylim(-0.142, 0.142)
ax.set_xlabel(r'shear sum rule $\mathcal{S}_E$', labelpad=2)
ax.set_ylabel(r'Hall sum rule $\mathcal{S}_H$', labelpad=1)
ax.set_xticks([0.25, 0.3125, 0.375, 0.4375, 0.5])
ax.set_xticklabels([r'$\frac{1}{4}$', '', r'$\frac{3}{8}$', '', r'$\frac{1}{2}$'])
ax.set_yticks([-0.125, 0, 0.125])
ax.set_yticklabels([r'$-\frac{1}{8}$', r'$0$', r"$\frac{1}{8}$"])
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)

# universal disc (all unitary CFT3)
ax.add_patch(Circle((c0, 0), R, fc='#f1f1f1', ec='none', zorder=0))
# its boundary: large-N Chern-Simons-matter theories (t4 = -4 cos pi lam, alpha_T = 4 sin pi lam)
th = np.linspace(0, 2 * np.pi, 721)
ax.plot(c0 + R * np.cos(th), R * np.sin(th), color=C2, lw=1.1, zorder=2)
# mixed free theories fill the diameter
ax.plot([c0 - R, c0 + R], [0, 0], color=INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)

# gap discs, radius tau_*(Delta_g)/32
for dg, alpha in ((10, 0.30), (30, 0.50), (100, 0.85)):
    r = tau_star(dg) / 32.0
    ax.add_patch(Circle((c0, 0), r, fc=C1, alpha=alpha * 0.45, ec=C1, lw=0.6, zorder=3))
labs = {10: (0.058, 0.090), 30: (0.072, 0.062), 100: (0.082, 0.034)}
for dg, (dx, dy) in labs.items():
    r = tau_star(dg) / 32.0
    a = np.radians(50 if dg == 10 else (35 if dg == 30 else 20))
    ax.annotate(r'$\Delta_g=%d$' % dg, xy=(c0 + r * np.cos(a), r * np.sin(a)), xytext=(c0 + dx, dy),
                fontsize=6.3, color=C1, ha='left', va='center',
                arrowprops=dict(arrowstyle='-', color=C1, lw=0.5, shrinkA=0.5, shrinkB=0.5))

# Einstein gravity at the centre
ax.plot([c0], [0], 'o', ms=2.6, color=INK, zorder=5)
ax.annotate('Einstein', xy=(c0, 0), xytext=(c0 - 0.004, -0.030), fontsize=FS, ha='center', va='top',
            arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5, shrinkA=0.5, shrinkB=2))
# free fields on the boundary
ax.plot([c0 + R], [0], 's', ms=3.2, color=INK, zorder=5)
ax.plot([c0 - R], [0], 'D', ms=3.0, color=INK, zorder=5)
ax.text(c0 + R - 0.006, -0.012, 'free\nscalar', ha='right', va='top', fontsize=6.5, linespacing=0.95)
ax.text(c0 - R + 0.006, -0.012, 'free\nfermion', ha='left', va='top', fontsize=6.5, linespacing=0.95)
# Chern-Simons-matter label on the circle
a = np.radians(128)
ax.annotate('large-$N$ Chern-Simons-matter', xy=(c0 + R * np.cos(a), R * np.sin(a)),
            xytext=(c0 - 0.158, 0.137), fontsize=6.5, color=C2, ha='left', va='top',
            arrowprops=dict(arrowstyle='-', color=C2, lw=0.5, shrinkA=0.5, shrinkB=0.5))
ax.text(c0 + 0.158, -0.137, 'all CFT$_3$', ha='right', va='bottom', fontsize=6.5, color=INK2)
ax.text(c0 + 0.158, 0.137, r'$\Delta_g\geq1000$: radius $9\times10^{-10}$', ha='right', va='top', fontsize=6.3,
        color=C1)

out = os.path.join(HERE, 'fig_disc.pdf')
fig.savefig(out)
fig.savefig(os.path.join(HERE, 'fig_disc_preview.png'), dpi=220)
print('wrote', out)
for dg in (10, 30, 100, 1000):
    print(dg, tau_star(dg) / 32.0)

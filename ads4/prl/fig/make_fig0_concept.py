#!/usr/bin/env python3
"""Fig. 1 of the Letter (v19): concept figure.
(a) Energy flux at null infinity created by the stress tensor, E(theta) ~ 1 + A cos(4 theta - phi), measured by energy
    detectors on the celestial circle. Dashed: A = 1, allowed in every unitary CFT_3 (free scalar / free fermion).
    Solid: the isotropic pattern; at Delta_g >= 1000 the bound gives A <= 7.2e-9, invisible at this scale.
(b) Schematic spectrum of T x T in the (spin, twist) plane. Below the gap: T, its double traces and light single traces
    of spin 0 and 2 (the light side of the sum rule, fixed by the bulk couplings). Above the gap: heavy states, on which
    the functional must be positive.
Usage (from ads4/prl/fig): python3 make_fig0_concept.py -> fig0_concept.pdf, fig0_concept_preview.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID = '#1a1a1a', '#555555', '#e6e6e6'
C1, C2, HEAVY = '#2a78d6', '#eb6834', '#dbe9f6'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
                     'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False, 'font.size': 7.5,
                     'axes.linewidth': 0.6, 'axes.edgecolor': INK2, 'text.color': INK, 'axes.labelcolor': INK,
                     'pdf.fonttype': 42})
FS = 7.5
fig = plt.figure(figsize=(3.375, 1.95))
axA = fig.add_axes([0.0, 0.0, 0.47, 1.0])
axB = fig.add_axes([0.575, 0.15, 0.41, 0.80])

# ---------------------------------------------------------------- (a) energy flux at null infinity
axA.set_aspect('equal'); axA.axis('off')
axA.set_xlim(-1.55, 1.55); axA.set_ylim(-1.62, 1.48)
R = 1.18                                   # celestial circle (null infinity)
axA.add_patch(Circle((0, 0), R, fill=False, ec=INK2, lw=0.6, ls=(0, (2, 2))))
th = np.linspace(0, 2 * np.pi, 721)
r0 = 0.62
axA.plot(r0 * (1 + np.cos(4 * th)) / 2 * 1.45 * np.cos(th), r0 * (1 + np.cos(4 * th)) / 2 * 1.45 * np.sin(th),
         color=C2, lw=0.9, ls=(0, (3, 1.6)))
axA.fill(r0 * np.cos(th), r0 * np.sin(th), color=C1, alpha=0.15, lw=0)
axA.plot(r0 * np.cos(th), r0 * np.sin(th), color=C1, lw=1.3)
axA.plot([0], [0], 'o', ms=3.2, color=INK)
axA.text(0.17, -0.15, r'$T$', ha='left', va='top', fontsize=FS)
# two energy detectors on the celestial circle
for ang, lab, dx, ha in ((62, r'$\mathcal{E}(\vec n_1)$', 0.06, 'left'), (118, r'$\mathcal{E}(\vec n_2)$', -0.06, 'right')):
    a = np.radians(ang)
    x, y = R * np.cos(a), R * np.sin(a)
    axA.add_patch(FancyBboxPatch((x - 0.07, y - 0.07), 0.14, 0.14, boxstyle='round,pad=0.01,rounding_size=0.03',
                                 fc=INK, ec='none', zorder=4))
    axA.text(x + dx, y + 0.12, lab, ha=ha, va='bottom', fontsize=FS)
axA.text(-1.5, 1.45, '(a)', ha='left', va='top', fontsize=FS)
# key below the circle
axA.plot([-1.30, -1.05], [-1.50, -1.50], color=C2, lw=0.9, ls=(0, (3, 1.6)))
axA.text(-0.98, -1.50, r'any CFT$_3$', ha='left', va='center', fontsize=6.8)
axA.plot([0.12, 0.37], [-1.50, -1.50], color=C1, lw=1.3)
axA.text(0.44, -1.50, r'$\Delta_g\geq1000$', ha='left', va='center', fontsize=6.8)

# ---------------------------------------------------------------- (b) light and heavy states
axB.set_xlim(-0.6, 12.6); axB.set_ylim(0, 10)
for sp in ('top', 'right'):
    axB.spines[sp].set_visible(False)
axB.set_xticks([]); axB.set_yticks([])
axB.set_xlabel(r'spin $J$', labelpad=2)
axB.set_ylabel(r'twist $\Delta-J$', labelpad=2)
GAP = 5.6
axB.fill_between([-0.6, 12.6], GAP, 10, color=HEAVY, lw=0, zorder=0)
axB.axhline(GAP, color=INK2, lw=0.8, zorder=1)
axB.text(12.4, GAP + 0.25, r'$\Delta_g$', ha='right', va='bottom', fontsize=FS)
axB.text(6.0, 8.55, 'heavy states', ha='center', va='center', fontsize=FS)
axB.text(6.0, 7.15, r'need $\Omega\geq0$', ha='center', va='center', fontsize=FS)
# light side: T (twist 1, spin 2), double traces [TT]_{n,J} (twist 2 + 2n), light single traces of spin 0 and 2
Js = np.arange(0, 9, 2)
axB.plot(Js, 1.9 + 0 * Js, 'o', ms=2.4, color=INK2, zorder=3)
axB.plot(Js, 3.1 + 0 * Js, 'o', ms=2.4, color=INK2, zorder=3)
axB.text(9.3, 2.5, r'$[TT]$', ha='left', va='center', fontsize=6.8, color=INK2)
axB.plot([2], [0.85], '*', ms=7, color=C1, zorder=4)
axB.text(2.6, 0.82, r'$T$', ha='left', va='center', fontsize=FS)
axB.plot([0, 2], [4.2, 4.2], 'o', ms=3.0, mfc='white', mec=C2, mew=0.9, zorder=3)
axB.text(12.4, 4.85, r'light: $G$, $\hat g_3$, $g_4$', ha='right', va='center', fontsize=6.8)
axB.text(-0.2, 9.85, '(b)', ha='left', va='top', fontsize=FS)

out = os.path.join(HERE, 'fig0_concept.pdf')
fig.savefig(out)
fig.savefig(os.path.join(HERE, 'fig0_concept_preview.png'), dpi=220)
print('wrote', out)

"""make_fig.py -- Fig. 1 of the Letter: how close to Einstein (a) and the R^4 upper end (b), versus the higher-spin gap.

Reads only the published functional data in ../data/ (nf02_psi.json: Omega^EE weights, exact rationals; q3star_factored.json:
Omega^Q_up flat target Q3*), and evaluates the exact light-side integrals
    16 L_A(Lam)  = 16 int_0^1 psi2^D(p) dp / (p^2 + 9/(4 Lam^2))        (Result (1); K = 1; c1 = 0 exactly, C1_FINAL)
    2 L_AdS(Lam) =  2 int_0^1 psi2^Q(p) dp / (p^2 + 9/(4 Lam^2))        (Result (2a), published constants)
term by term in mpmath (recursion I_n = 1/(n-1) - a^2 I_{n-2}, a = 3/(2 Lam)), as ads4/code/repro/bounds.py does.
Plotted (v18: Computer Modern fonts, no overlapping labels, Delta_g >= 1000 only; same numbers as v14):
  (a) log scale: the bound tau0 = 240 sqrt(16 L_A)/Lam^4 on tau = sqrt(t4^2 + alpha_T^2), against the conformal-collider bound
      tau <= 4 that holds in every unitary CFT_3 (Supplement Sec. S5.F; saturated by a free scalar or free fermion), with the
      factor 4/tau0 at Delta_g = 1000 drawn as an arrow. Comparison curve (not an AdS bound, no status markers): the flat D = 4
      optimum of Caron-Huot-Li-Parra-Martinez-Simmons-Duffin, |g3|^2 M^8 <= 24.9 ln(M/m_IR) - 27.6, transplanted with M -> Lam,
      m_IR -> 3/(2L), converted with tau = 240 |g3|. The headline label also gives the gap-unit constant |g3| Lam^4 = sqrt(16 L_A).
  (b) log scale: the upper end of g4~ Lam^6/8piG: (2a) 2 L_AdS; (2b) 2 L_AdS + 2 eps L_A, eps = 0.040; type II / ABJM-IIA value
      2 zeta(3), which in these units does not depend on the gap (Lam <-> M L), drawn as a horizontal line.
Status code (unchanged from v11): filled dot = Delta_g = 1000, the headline; open diamond = 300 (near-gap large-spin positivity
an open hypothesis); faded lines = values of the implication only (positivity examined only at the markers); panel (a) solid for
Delta_g >= 1000 = carried over from 1000 under the no-growth premise (Supplement Sec. S8), with small open circles at 2000 and
5000 where near-gap states were recomputed. Panel (b) has no solid part (Omega^Q_up was examined only at 300 and 1000). Nothing
below 300 is claimed (grey band).
History: v6 (referee 3, m14) diamonds/key; v10 solid part above 1000; v11 (referee 5 m9) comparison curve without status
markers, (2b) label "needs both functionals"; v14 annotation without c1; v15 log scales, collider reference line, 10^8 arrow.
Usage (from ads4/prl/fig):  python3 make_fig.py  -> fig1_bounds.pdf, fig1_bounds_preview.png and a printed check of the headline
constants (tau0, |g3|Lam^4, 2L_AdS, (2b), ratio to the transplanted optimum, ABJM fraction, collider ratio).
Colours: first three slots of the dataviz reference palette (all-pairs CVD-validated); every series also has its own line style,
and every series is direct-labelled, so identity never rests on colour.
"""
import json, os
from fractions import Fraction as Fr
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
s2 = [Fr(x) for x in nf['psi2_exact']]; s3 = [Fr(x) for x in nf['psi3_exact']]
mom = lambda c, j: sum(cn / (n + j + 1) for n, cn in enumerate(c))
k0 = -2 * mom(s2, 6) + 2 * mom(s3, 4)                     # normalisation k(0) = 1
psiD = [x / k0 for x in s2]
q = json.load(open(os.path.join(DATA, 'q3star_factored.json')))
psiQ = [Fr(x) for x in q['psi2_expanded_exact']]


def Lpole(c, Lam):
    a = mp.mpf(3) / (2 * mp.mpf(Lam)); a2 = a * a
    I = [mp.atan(1 / a) / a, mp.log(1 + 1 / a2) / 2]
    for n in range(2, len(c)):
        I.append(mp.mpf(1) / (n - 1) - a2 * I[n - 2])
    return sum(mp.mpf(cn.numerator) / cn.denominator * I[n] for n, cn in enumerate(c))


lam = lambda Dg: mp.sqrt((mp.mpf(Dg) - 2) * (mp.mpf(Dg) - 1))
EPS = mp.mpf('0.040')


def vals(Dg):
    L = lam(Dg)
    L16 = 16 * Lpole(psiD, L); two = 2 * Lpole(psiQ, L)
    return dict(L=float(L), g3=float(mp.sqrt(L16)), tau=float(240 * mp.sqrt(L16) / L**4),
                q2a=float(two), q2b=float(two + 2 * EPS * L16 / 16),
                chlps=float(mp.sqrt(mp.mpf('24.9') * mp.log(2 * L / 3) - mp.mpf('27.6'))),
                chlps_tau=float(240 * mp.sqrt(mp.mpf('24.9') * mp.log(2 * L / 3) - mp.mpf('27.6')) / L**4))


for Dg in (300, 1000, 2000, 5000, 10000, 100000):
    v = vals(Dg)
    print('Delta_g = %4d: tau0 = %.5e  |g3|Lam^4 <= %.2f  2L_AdS = %.4f  (2b) %.2f  CHLPS-transplant %.2f  ratio %.1f  '
          'ABJM 2zeta(3)/2L_AdS = %.2f %%' % (Dg, v['tau'], v['g3'], v['q2a'], v['q2b'], v['chlps'], v['g3'] / v['chlps'],
                                              100 * 2 * float(mp.zeta(3)) / v['q2a']))
print('collider bound tau <= 4 over tau0(1000) = %.3e;  A = tau0/4 <= %.3e at 1000'
      % (4 / vals(1000)['tau'], vals(1000)['tau'] / 4))

grid = sorted(set([int(round(10 ** (3.0 + 0.02 * i))) for i in range(36)] + [1000, 2000, 5000]))
grid = [g for g in grid if 1000 <= g <= 5000]
V = {Dg: vals(Dg) for Dg in grid}
xs = grid
ZETA3x2 = 2 * float(mp.zeta(3))
TAU_COLL = 4.0          # conformal-collider bound tau <= 4 (every unitary CFT_3; Supplement Sec. S5.F)

# ------------------------------------------------------------------ style: single PRL column, Computer Modern like REVTeX
INK, INK2, GRID = '#1a1a1a', '#555555', '#e6e6e6'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
FS = 7.5
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
                     'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False,
                     'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'axes.edgecolor': INK2, 'axes.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
                     'xtick.direction': 'in', 'ytick.direction': 'in',
                     'text.color': INK, 'axes.labelcolor': INK, 'pdf.fonttype': 42})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.375, 3.5), sharex=True,
                               gridspec_kw=dict(hspace=0.08, height_ratios=[1.25, 1]))
XMIN, XMAX = 900, 5600
for ax in (ax1, ax2):
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(XMIN, XMAX)
    ax.grid(True, which='major', color=GRID, lw=0.5, zorder=0)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
LAB = dict(fontsize=FS, color=INK)

# ------------------------------------------------------------------ (a) how close to Einstein
ax1.axhline(TAU_COLL, color=INK2, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=2)
ax1.text(1500, TAU_COLL * 1.8, r'collider bound, every CFT$_3$', ha='left', va='bottom', **LAB)
ax1.plot(xs, [V[x]['chlps_tau'] for x in xs], color=C2, ls=(0, (4, 2)), lw=1.0, zorder=3)
ax1.plot(xs, [V[x]['tau'] for x in xs], color=C1, lw=1.4, alpha=0.45, zorder=3)   # beyond 1000: assumes positivity at own gap
ax1.plot([2000, 5000], [V[2000]['tau'], V[5000]['tau']], 'o', ms=3.2, mfc='white', mec=C1, mew=0.9, zorder=4)
t1000 = V[1000]['tau']
ax1.plot([1000], [t1000], 'o', ms=4.8, mfc=C1, mec='white', mew=0.8, zorder=5, clip_on=False)
ang = [160.0 / 9.0 / (x - 2) * (1 - 49.0 / (2 * x + 9) ** 2) ** -0.5 for x in xs]
ax1.plot(xs, ang, color=C3, ls=(0, (1, 1.4)), lw=1.4, zorder=3)
# direct labels, each in empty space next to its own line
v26 = vals(2600)
ang26 = 160.0 / 9.0 / (2600 - 2) * (1 - 49.0 / (2 * 2600 + 9) ** 2) ** -0.5
ax1.text(2600, ang26 * 2.2, 'without positivity', ha='left', va='bottom', **LAB)
ax1.text(1080, t1000 * 2.0, r'this work, $\tau\leq2.89\times10^{-8}$', ha='left', va='bottom', **LAB)
ax1.text(1080, vals(2300)['chlps_tau'] / 1.8, 'flat-space optimum', ha='left', va='top', color=INK2, fontsize=FS)
# the factor between the collider bound and this work, as an arrow at Delta_g = 1000
ax1.annotate('', xy=(1000, t1000 * 2.4), xytext=(1000, TAU_COLL / 1.6),
             arrowprops=dict(arrowstyle='-|>', color=INK2, lw=0.7, mutation_scale=6, shrinkA=0, shrinkB=0), zorder=4)
ax1.text(1080, 1e-4, r'$1.4\times10^{8}$ times' + '\nsmaller', ha='left', va='center', linespacing=1.15, **LAB)
ax1.set_ylim(1e-12, 1e3)
ax1.set_yticks([1e-12, 1e-9, 1e-6, 1e-3, 1e0])
ax1.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax1.set_ylabel(r'bound on $\tau$')
ax1.text(0.97, 0.97, '(a)', transform=ax1.transAxes, ha='right', va='top', **LAB)

# ------------------------------------------------------------------ (b) the R^4 upper end (Omega^Q_up examined at 1000 only)
for key, col, ls in (('q2b', C3, (0, (1.2, 1.4))), ('q2a', C1, '-')):
    ax2.plot(xs, [V[x][key] for x in xs], color=col, ls=ls, lw=1.2, alpha=0.45, zorder=3)
    ax2.plot([1000], [V[1000][key]], 'o', ms=4.8, mfc=col, mec='white', mew=0.8, zorder=5, clip_on=False)
ax2.axhline(ZETA3x2, color=INK2, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=2)
ax2.text(1500, ZETA3x2 * 1.35, r'ABJM, $2\zeta(3)$', ha='left', va='bottom', **LAB)
ax2.text(1500, vals(1500)['q2b'] * 1.45, 'light scalars allowed', ha='left', va='bottom', **LAB)
ax2.text(1500, vals(1500)['q2a'] / 1.5, 'no light scalar', ha='left', va='top', **LAB)
ax2.set_ylim(1, 3000)
ax2.set_yticks([1, 10, 100, 1000]); ax2.set_yticklabels(['1', '10', '100', '1000'])
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.set_ylabel(r'bound on $\tilde{g}_4\Lambda^6/8\pi G$')
ax2.set_xlabel(r'higher-spin gap $\Delta_g$')
ax2.set_xticks([1000, 2000, 3000, 5000]); ax2.set_xticklabels(['1000', '2000', '3000', '5000'])
ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.text(0.97, 0.97, '(b)', transform=ax2.transAxes, ha='right', va='top', **LAB)

fig.align_ylabels((ax1, ax2))
fig.subplots_adjust(left=0.15, right=0.975, top=0.985, bottom=0.10)
out = os.path.join(HERE, 'fig1_bounds.pdf')
fig.savefig(out)
fig.savefig(os.path.join(HERE, 'fig1_bounds_preview.png'), dpi=220)
print('wrote', out)

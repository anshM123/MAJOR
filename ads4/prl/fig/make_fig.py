"""make_fig.py -- Fig. 1 of the Letter: how close to Einstein (a) and the R^4 upper end (b), versus the higher-spin gap.

Reads only the published functional data in ../data/ (nf02_psi.json: Omega^EE weights, exact rationals; q3star_factored.json:
Omega^Q_up flat target Q3*), and evaluates the exact light-side integrals
    16 L_A(Lam)  = 16 int_0^1 psi2^D(p) dp / (p^2 + 9/(4 Lam^2))        (Result (1); K = 1; c1 = 0 exactly, C1_FINAL)
    2 L_AdS(Lam) =  2 int_0^1 psi2^Q(p) dp / (p^2 + 9/(4 Lam^2))        (Result (2a), published constants)
term by term in mpmath (recursion I_n = 1/(n-1) - a^2 I_{n-2}, a = 3/(2 Lam)), as ads4/code/repro/bounds.py does.
Plotted (v15 redesign, for a broad reader; same numbers as v14):
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

grid = [int(round(10 ** (2.48 + 0.02 * i))) for i in range(int((3.74 - 2.48) / 0.02) + 1)]
grid = sorted(set([g for g in grid if g >= 300] + [300, 1000, 2000, 5000]))
V = {Dg: vals(Dg) for Dg in grid}
xs = grid
ZETA3x2 = 2 * float(mp.zeta(3))
TAU_COLL = 4.0          # conformal-collider bound tau <= 4 (every unitary CFT_3; Supplement Sec. S5.F)

# ------------------------------------------------------------------ style (single-column PRL width 3.375 in; size as v14)
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
FS = 6.6
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'mathtext.fontset': 'dejavuserif',
                     'font.size': 7.5, 'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'axes.edgecolor': INK2, 'axes.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
                     'text.color': INK, 'axes.labelcolor': INK, 'pdf.fonttype': 42})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.375, 3.6), sharex=True,
                               gridspec_kw=dict(hspace=0.10, height_ratios=[1.22, 1]))
XMIN, XMAX = 230, 5600
for ax in (ax1, ax2):
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(XMIN, XMAX)
    ax.axvspan(XMIN, 300, color='#f0efec', lw=0, zorder=0)
    ax.grid(True, which='major', color=GRID, lw=0.5, zorder=0)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def series(ax, key, color, ls, solid_from=None, status=True):
    ys = [V[x][key] for x in xs]
    if not status:   # comparison curve: full strength, no status markers (referee 5, m9)
        ax.plot(xs, ys, color=color, ls=ls, lw=1.0, zorder=3)
        return
    ax.plot(xs, ys, color=color, ls=ls, lw=1.0, alpha=0.45, zorder=3)   # implication only between markers
    if solid_from is not None:   # Delta_g >= 1000 for Omega^EE, given the no-growth premise (Supplement Sec. S8)
        xx = [x for x in xs if x >= solid_from]
        ax.plot(xx, [V[x][key] for x in xx], color=color, ls=ls, lw=1.4, zorder=3)
        ax.plot([2000, 5000], [V[2000][key], V[5000][key]], 'o', ms=3.0, mfc='white', mec=color, mew=0.9, zorder=4)
    ax.plot([300], [V[300][key]], 'D', ms=3.6, mfc='white', mec=color, mew=1.1, zorder=4)    # one hypothesis open
    ax.plot([1000], [V[1000][key]], 'o', ms=4.6, mfc=color, mec='white', mew=0.8, zorder=5)  # headline


# ------------------------------------------------------------------ (a) how close to Einstein
ax1.axhline(TAU_COLL, color=INK2, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=2)
ax1.text(XMAX, TAU_COLL * 2.2, r'any unitary CFT$_3$: $\tau\leq4$ (conformal collider)', fontsize=FS, color=INK,
         ha='right', va='bottom')
series(ax1, 'chlps_tau', C2, (0, (4, 2)), status=False)
series(ax1, 'tau', C1, '-', solid_from=1000)
t1000 = V[1000]['tau']
ax1.annotate('', xy=(1000, t1000 * 1.9), xytext=(1000, TAU_COLL / 1.6),
             arrowprops=dict(arrowstyle='-|>', color=INK2, lw=0.8, mutation_scale=7, shrinkA=0, shrinkB=0), zorder=4)
ax1.text(870, 10 ** ((mp.log10(TAU_COLL) + mp.log10(t1000)) / 2 + 0.3), r'$\times\,1.4\cdot10^{8}$' + '\ntighter',
         fontsize=FS + 0.4, color=INK, ha='right', va='center', linespacing=1.1)
ax1.text(1130, t1000 * 2.5, r'$\tau\leq2.89\cdot10^{-8}$' + '\n' + r'$|\hat g_3|\Lambda^4\leq120$', fontsize=FS,
         color=INK, ha='left', va='bottom', linespacing=1.15)
ax1.text(330, V[330]['tau'] * 0.62 if 330 in V else 1e-6, 'this work', fontsize=FS, color=INK, ha='left', va='top',
         rotation=-33)
ax1.text(1500, V[1514]['chlps_tau'] * 0.42 if 1514 in V else 1e-10,
         'flat-space optimum, transplanted\n(comparison, not an AdS bound)', fontsize=FS - 0.6, color=INK2,
         ha='left', va='top', linespacing=1.1)
ax1.set_ylim(1e-11, 1e3)
ax1.set_yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0, 1e2])
ax1.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax1.set_ylabel(r'bound on $\tau=\sqrt{t_4^2+\alpha_T^2}$')
ax1.text(0.015, 0.985, '(a)', transform=ax1.transAxes, va='top', fontsize=7.5, color=INK)
from matplotlib.lines import Line2D
key = [Line2D([], [], marker='o', ls='', ms=4.0, mfc=INK2, mec='white', label=r'$\Delta_g=1000$: headline'),
       Line2D([], [], marker='D', ls='', ms=3.3, mfc='white', mec=INK2, mew=1.0, label=r'$\Delta_g=300$: one hypothesis open'),
       Line2D([], [], color=INK2, lw=1.4, marker='o', ms=3.0, mfc='white', mec=INK2, mew=0.9,
              label=r'solid: $\Delta_g\geq1000$, given premises'),
       Line2D([], [], color=INK2, alpha=0.45, lw=1.0, label='faded: implication only')]
ax1.legend(handles=key, loc='upper right', bbox_to_anchor=(1.0, 0.80), fontsize=5.7, frameon=False, handlelength=1.4,
           borderaxespad=0.2, labelspacing=0.22)
ax1.text(262, 3e-5, 'nothing claimed', rotation=90, fontsize=6.0, color=INK2, ha='center', va='center')

# ------------------------------------------------------------------ (b) the R^4 upper end
series(ax2, 'q2b', C3, (0, (1.2, 1.4)))
series(ax2, 'q2a', C1, '-')
ax2.axhline(ZETA3x2, color=INK2, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=2)
ax2.text(XMAX, ZETA3x2 * 1.25, r'type II strings / ABJM: $2\zeta(3)$', fontsize=FS, color=INK, ha='right', va='bottom')
ax2.text(1080, ZETA3x2 * 4.2, r'$2.6\%$ of the bound' + '\n' + r'at $\Delta_g=1000$', fontsize=FS - 0.3, color=INK2,
         ha='left', va='bottom', linespacing=1.1)
ax2.annotate('', xy=(1000, ZETA3x2 * 1.15), xytext=(1000, V[1000]['q2a'] / 1.3),
             arrowprops=dict(arrowstyle='-|>', color=MUTED, lw=0.7, mutation_scale=6, shrinkA=0, shrinkB=0), zorder=2)
ax2.text(320, V[320]['q2b'] * 1.45 if 320 in V else 200, r'(2b) light scalars allowed; needs both functionals',
         fontsize=FS - 0.3, color=INK, ha='left', va='bottom')
ax2.text(1300, V[1300]['q2a'] / 1.5 if 1300 in V else 50, r'(2a) no scalar with $m<0.46\Lambda$',
         fontsize=FS - 0.3, color=INK, ha='left', va='top')
ax2.set_ylim(1, 1000)
ax2.set_yticks([1, 10, 100, 1000]); ax2.set_yticklabels(['1', '10', '100', '1000'])
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.set_ylabel(r'upper end of $\tilde{g}_4\Lambda^6/8\pi G$')
ax2.set_xlabel(r'higher-spin gap $\Delta_g$')
ax2.set_xticks([300, 500, 1000, 2000, 5000]); ax2.set_xticklabels(['300', '500', '1000', '2000', '5000'])
ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.text(0.015, 0.97, '(b)', transform=ax2.transAxes, va='top', fontsize=7.5, color=INK)
ax2.text(262, 31, 'nothing claimed', rotation=90, fontsize=6.0, color=INK2, ha='center', va='center')

fig.align_ylabels((ax1, ax2))
fig.subplots_adjust(left=0.165, right=0.98, top=0.985, bottom=0.10)
out = os.path.join(HERE, 'fig1_bounds.pdf')
fig.savefig(out)
fig.savefig(os.path.join(HERE, 'fig1_bounds_preview.png'), dpi=220)
print('wrote', out)

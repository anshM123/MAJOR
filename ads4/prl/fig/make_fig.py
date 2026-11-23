"""make_fig.py -- Fig. 1 of the Letter: the bounds of Results D and Q as functions of the higher-spin gap.

Reads only the published functional data in ../data/ (nf02_psi.json: Omega^EE weights, exact rationals; q3star_factored.json:
Omega^Q_up flat target Q3*), and evaluates the exact light-side integrals
    16 L_A(Lam)  = 16 int_0^1 psi2^D(p) dp / (p^2 + 9/(4 Lam^2))        (Result D; K = 1, c1 = 0)
    2 L_AdS(Lam) =  2 int_0^1 psi2^Q(p) dp / (p^2 + 9/(4 Lam^2))        (Result Q, published constants (2a))
term by term in mpmath (recursion I_n = 1/(n-1) - a^2 I_{n-2}, a = 3/(2 Lam)), as ads4/code/repro/bounds.py does.
Plotted:
  (a) |g3| Lam^4 <= sqrt(16 L_A)  [equivalently sqrt(t4^2 + alpha_T^2) <= 240 sqrt(16 L_A)/Lam^4], and for comparison the
      flat D = 4 optimum of Caron-Huot-Li-Parra-Martinez-Simmons-Duffin, |g3|^2 M^8 <= 24.9 ln(M/m_IR) - 27.6, transplanted with
      M -> Lam, m_IR -> 3/(2L) (not an AdS bound).
  (b) the upper end of g4~ Lam^6/8piG: (2a) 2 L_AdS; (2b) 2 L_AdS + 2 eps L_A, eps = 0.040; type II / ABJM-IIA value 2 zeta(3).
Values between and beyond Delta_g = 300, 1000 are values of the implication with the functionals fixed; positivity was
examined only at 300 and 1000. Nothing below 300 is claimed.
Usage (from ads4/prl/fig):  python3 make_fig.py  -> fig1_bounds.pdf (and a printed check of the headline constants).
Colours: first three slots of the dataviz reference palette (all-pairs CVD-validated); every series also has its own line style.
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
                chlps=float(mp.sqrt(mp.mpf('24.9') * mp.log(2 * L / 3) - mp.mpf('27.6'))))


for Dg in (300, 1000):
    v = vals(Dg)
    print('Delta_g = %4d: tau0 = %.5e  |g3|Lam^4 <= %.2f  2L_AdS = %.4f  (2b) %.2f  CHLPS-transplant %.2f  ratio %.1f  '
          'ABJM 2zeta(3)/2L_AdS = %.2f %%' % (Dg, v['tau'], v['g3'], v['q2a'], v['q2b'], v['chlps'], v['g3'] / v['chlps'],
                                              100 * 2 * float(mp.zeta(3)) / v['q2a']))

grid = [int(round(10 ** (2.48 + 0.02 * i))) for i in range(int((3.74 - 2.48) / 0.02) + 1)]
grid = sorted(set([g for g in grid if g >= 300] + [300, 1000]))
V = {Dg: vals(Dg) for Dg in grid}
xs = grid

# ------------------------------------------------------------------ style (single-column PRL width 3.375 in)
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'mathtext.fontset': 'dejavuserif',
                     'font.size': 7.5, 'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'axes.edgecolor': INK2, 'axes.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
                     'text.color': INK, 'axes.labelcolor': INK, 'pdf.fonttype': 42})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.375, 3.6), sharex=True, gridspec_kw=dict(hspace=0.12))
XMIN, XMAX = 230, 5600
for ax in (ax1, ax2):
    ax.set_xscale('log'); ax.set_xlim(XMIN, XMAX)
    ax.axvspan(XMIN, 300, color='#f0efec', lw=0, zorder=0)
    ax.grid(True, which='major', color=GRID, lw=0.5, zorder=0)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def series(ax, key, color, ls, label, lab_xy, lab_ha='left'):
    ax.plot(xs, [V[x][key] for x in xs], color=color, ls=ls, lw=1.0, alpha=0.45, zorder=3)  # positivity not examined between markers
    ax.plot([300], [V[300][key]], 'o', ms=4.5, mfc='white', mec=color, mew=1.2, zorder=4)
    ax.plot([1000], [V[1000][key]], 'o', ms=4.5, mfc=color, mec='white', mew=0.8, zorder=4)
    ax.annotate(label, xy=lab_xy, color=INK, fontsize=6.8, ha=lab_ha, va='center')


# (a) cubic
series(ax1, 'g3', C1, '-', r'this work, item (1): $\sqrt{16L_A}$', (1150, 92))
series(ax1, 'chlps', C2, (0, (4, 2)), r'flat $D{=}4$ optimum [CHLPS], transplanted', (330, 23))
ax1.set_ylim(0, 160); ax1.set_yticks([0, 50, 100, 150])
ax1.set_ylabel(r'upper bound on $|\hat{g}_3|\,\Lambda^4$')
ax1.annotate(r'$\sqrt{t_4^2+\alpha_T^2}\leq 2.89\times10^{-8}$', xy=(1000, V[1000]['g3']), xytext=(1150, 145),
             fontsize=6.8, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.5))
ax1.text(0.015, 0.97, '(a)', transform=ax1.transAxes, va='top', fontsize=7.5)
ax1.text(262, 75, 'not claimed', rotation=90, fontsize=6.3, color=INK2, ha='center', va='center')

# (b) quartic
series(ax2, 'q2b', C3, (0, (1.2, 1.4)), r'(2b): any light scalars, $\epsilon=0.04$', (330, 205))
series(ax2, 'q2a', C1, '-', r'(2a): no scalar with $m<0.46\Lambda$ in $T\times T$', (560, 50))
ax2.axhline(2 * float(mp.zeta(3)), color=INK2, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=2)
ax2.annotate(r'type II / ABJM-IIA: $2\zeta(3)$', xy=(1500, 2.404), xytext=(1500, 20), fontsize=6.8, color=INK, ha='left',
             arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.5))
ax2.set_ylim(0, 230); ax2.set_yticks([0, 50, 100, 150, 200])
ax2.set_ylabel(r'upper end of $\tilde{g}_4\Lambda^6/8\pi G$')
ax2.set_xlabel(r'higher-spin gap $\Delta_g$')
ax2.set_xticks([300, 500, 1000, 2000, 5000]); ax2.set_xticklabels(['300', '500', '1000', '2000', '5000'])
ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.text(0.015, 0.97, '(b)', transform=ax2.transAxes, va='top', fontsize=7.5)
ax2.text(262, 115, 'not claimed', rotation=90, fontsize=6.3, color=INK2, ha='center', va='center')

fig.subplots_adjust(left=0.155, right=0.98, top=0.985, bottom=0.10)
out = os.path.join(HERE, 'fig1_bounds.pdf')
fig.savefig(out)
fig.savefig(os.path.join(HERE, 'fig1_bounds_preview.png'), dpi=220)
print('wrote', out)

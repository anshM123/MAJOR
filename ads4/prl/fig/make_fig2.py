"""make_fig2.py -- Fig. 2 of the Letter (v6): coverage of (Pos) for Omega^EE at Delta_g = 1000 in the (J, m/Lambda) plane,
coloured by evidence grade. No computation of physics: every boundary below is transcribed from the deliverables named here.

Plane.  x = J + 1/2 (log scale; J = 0 sits at x = 1/2), y = m/Lambda, m^2 = (Delta-J-2)(Delta+J-1),
Lambda = sqrt((Dg-2)(Dg-1)), Dg = 1000.  Physical states lie on or above the twist edge
    m_gap(J)/Lambda = sqrt(1 + 2J/(Dg-1))         (supplement Sec. S2.A; crosses 1.3 between J = 344 and 345).
Impact parameter b = (2J+1) Lambda / m, i.e. b < B  <=>  x < (B/2) y.

Grades (ordinal, strongest first), with sources:
  certified   flat limit m -> infinity at fixed b: Omega^EE >= N/5 on every heavy state [certified]          CERT_CLOSE; Table S5
  [NE+-]      error budget with x3 on every model term (not a bound):
                near gap J = 0, 2 (to 1.2 Lambda; POS_ERRORS Sec. 3.1), J = 4, 5, 6, 8 (to 1.05; NEARGAP2, J = 8 on the
                lambda-level budget), J = 12, 16 (to 1.05; HIGHSPIN Sec. 4, SMALLNU_WINDOW Sec. 3.2), J = 21, 29 (to 1.10;
                SMALLNU_WINDOW Sec. 0, 3.1, 3.3; thin at the twist edge);
                strip m >= 1.3 Lambda, b < 5 (POS_ERRORS Sec. 2.4; one bulk model input tolerated x76; computed to 20 Lambda)
  theorem     conditional theorems (implications with measured/modelled premises):
                W: near gap, m <= 1.3 Lambda, 50 <= J <= 344 (40 sampled states, even J 26-344, odd 33-301; PW_IMAGES TW',
                   LARGESPIN_W Sec. 4, 9; P2_SPINNING Sec. 7); its premises are listed in the Letter/Supplement S7.C
                L: b >= 1.0421 Lambda (J + 1/2 >= 0.52105 m), J <= 10^6 (PW_IMAGES Sec. 6; LARGESPIN_L Sec. 3)
  [NE+M]/[M]  computed inside models or interpolated: near-gap J <= 32 (chain spins, masses to 5 Lambda; interpolated
                even J = 26, 30 included), even J = 36, 40, 48 (HIGHSPIN, POS_ERRORS Sec. 3.2, Table S5/S6);
                m >= 1.3 Lambda, 5 <= b < b_A(m) (POS_ERRORS Sec. 2.4: b_A = 8.5, 8.7, 10.5, 13, 15, 15.3, 13.8
                at m/Lambda = 1.3, 1.5, 2, 3, 5, 10, 20)
  hypothesis  no computation and no theorem: covered only by H_LJ(1000) / H_R' [D-sketch + M] -- near-gap odd J = 33-49 and
                even J = 34, 38, 42-46 (interpolated), m >= 1.3 Lambda with b_A(m) <= b < 1.0421 Lambda (this contains the
                near-gap band 345 <= J < 857, m = 1.3-1.65 Lambda, which neither theorem covers), and J > 10^6
                (LARGESPIN_ANALYTIC Sec. 5.2; LARGESPIN_L Sec. 8 inventory; Table S11)
  The star marks (Delta, J) = (1100, 100), the one state at which the spinning O(1/m^2) magnitudes (Lambda_amp ~ 2-3) were
  measured (DIRCV Sec. 8; LARGESPIN_W C2).
Usage (from ads4/prl/fig):  python3 make_fig2.py  -> fig2_coverage.pdf, fig2_coverage_preview.png
Colours: one-hue ordinal blue ramp (dataviz reference sequential steps 700/550/400/250; validated with
validate_palette.js --ordinal), neutral grey + hatch for "hypothesis only"; every region also carries a direct label.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
DG = 1000
LAM = np.sqrt((DG - 2) * (DG - 1))

# ---------------------------------------------------------------- geometry
def y_gap(x):                      # twist edge, x = J + 1/2
    return np.sqrt(1 + 2 * (np.asarray(x) - 0.5) / (DG - 1))


def x_edge(y):                     # inverse of the twist edge: largest physical x at mass y
    return (np.asarray(y) ** 2 - 1) * (DG - 1) / 2 + 0.5


BA_M = np.array([1.3, 1.5, 2, 3, 5, 10, 20]); BA_V = np.array([8.5, 8.7, 10.5, 13, 15, 15.3, 13.8])
def b_A(y):
    return np.interp(np.log(y), np.log(BA_M), BA_V)


XL_COEF = 0.52105 * LAM
XMIN, XMAX = 0.4, 4.0e6
JL_MAX = 1.0e6 + 0.5

NEPM = [(0, 1.2), (2, 1.2), (4, 1.05), (5, 1.05), (6, 1.05), (8, 1.05), (12, 1.05), (16, 1.05), (21, 1.10), (29, 1.10)]
NEM_TICKS = [36, 40, 48]

# ---------------------------------------------------------------- style
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e4e3df'
C_CERT, C_NEPM, C_THM, C_NEM = '#0d366b', '#1c5cab', '#3987e5', '#86b6ef'
C_HYP, C_HATCH, SURF = '#ecebe7', '#a9a8a2', '#ffffff'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'mathtext.fontset': 'dejavuserif',
                     'font.size': 7.0, 'axes.labelsize': 7.2, 'xtick.labelsize': 6.6, 'ytick.labelsize': 6.6,
                     'axes.edgecolor': INK2, 'axes.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2,
                     'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.minor.width': 0.35,
                     'ytick.minor.width': 0.35, 'text.color': INK, 'axes.labelcolor': INK, 'pdf.fonttype': 42,
                     'hatch.linewidth': 0.45})

W_IN, H_IN = 3.375, 3.25
fig = plt.figure(figsize=(W_IN, H_IN))
gs = fig.add_gridspec(3, 1, height_ratios=[0.11, 1.0, 0.82], hspace=0.0, left=0.13, right=0.985, top=0.985, bottom=0.235)
axS = fig.add_subplot(gs[0]); axU = fig.add_subplot(gs[1], sharex=axS); axD = fig.add_subplot(gs[2], sharex=axS)
for ax in (axS, axU, axD):
    ax.set_xscale('log'); ax.set_xlim(XMIN, XMAX); ax.set_facecolor(SURF)
    ax.tick_params(which='both', direction='out')

# ---------------------------------------------------------------- top strip: flat limit
axS.axvspan(XMIN, XMAX, color=C_CERT, lw=0)
axS.set_ylim(0, 1); axS.set_yticks([])
axS.text(np.sqrt(XMIN * XMAX), 0.5, r'flat limit $m\to\infty$: certified', color='white', ha='center', va='center',
         fontsize=6.6)
axS.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
for s in ('top', 'right', 'left', 'bottom'):
    axS.spines[s].set_visible(False)


def fbx(ax, y, xl, xr, color, hatch=None, z=1):
    xl = np.broadcast_to(xl, y.shape).astype(float); xr = np.broadcast_to(xr, y.shape).astype(float)
    ok = xr > xl
    if not ok.any():
        return
    kw = dict(color=color, lw=0, zorder=z, where=ok, interpolate=False)
    if hatch:
        ax.fill_betweenx(y, xl, np.where(ok, xr, xl), facecolor=color, edgecolor=C_HATCH, hatch=hatch, lw=0, zorder=z)
    else:
        ax.fill_betweenx(y, xl, np.where(ok, xr, xl), **kw)


# ---------------------------------------------------------------- upper panel: 1.3 <= m/Lambda <= 20 (log)
YU0, YU1 = 1.3, 20.0
yu = np.geomspace(YU0, YU1, 1500)
xe = np.minimum(x_edge(yu), XMAX)
fbx(axU, yu, XMIN, xe, C_HYP, hatch='//////', z=1)                                   # hypothesis only (H_R')
fbx(axU, yu, XMIN, np.minimum(xe, np.where(yu <= 5.0, 33.0, XMIN)), C_NEM, z=2)    # chain spins J <= 32 to 5 Lambda
fbx(axU, yu, 2.5 * yu, np.minimum(xe, b_A(yu) / 2 * yu), C_NEM, z=2)               # 5 <= b < b_A(m)
fbx(axU, yu, XMIN, np.minimum(xe, 2.5 * yu), C_NEPM, z=3)                          # b < 5 strip
fbx(axU, yu, XL_COEF * yu, np.minimum(xe, JL_MAX), C_THM, z=3)                     # Theorem L
axU.set_yscale('log'); axU.set_ylim(YU0, YU1)
axU.set_yticks([1.3, 2, 5, 10, 20]); axU.set_yticklabels(['1.3', '2', '5', '10', '20'])
axU.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

# ---------------------------------------------------------------- lower panel: 1 <= m/Lambda <= 1.3 (linear)
YD0, YD1 = 1.0, 1.3
yd = np.linspace(YD0, YD1, 1200)
xed = np.minimum(x_edge(yd), XMAX)
fbx(axD, yd, XMIN, np.minimum(xed, 33.0), C_NEM, z=2)                               # near-gap J <= 32
fbx(axD, yd, 33.0, np.minimum(xed, 50.0), C_HYP, hatch='//////', z=1)               # H_LJ: odd 33-49, interpolated even
fbx(axD, yd, 50.0, np.minimum(xed, 345.0), C_THM, z=2)                              # Theorem W (sampled 50-344)
axD.set_ylim(YD0, YD1); axD.set_yticks([1.0, 1.1, 1.2]); axD.set_yticklabels(['1.0', '1.1', '1.2'])
for J in NEM_TICKS:                                                                # same spins, chain masses to 5 Lambda
    axU.plot([J + 0.5] * 2, [YU0, 5.0], color=C_NEM, lw=1.1, solid_capstyle='butt', zorder=3)
for J in NEM_TICKS:                                                                # computed even J = 36, 40, 48
    axD.plot([J + 0.5] * 2, [y_gap(J + 0.5), YD1], color=C_NEM, lw=1.1, solid_capstyle='butt', zorder=3)
for J, ytop in NEPM:                                                               # [NE+-] near-gap budgets
    axD.plot([J + 0.5] * 2, [max(y_gap(J + 0.5), YD0), ytop], color=C_NEPM, lw=1.9, solid_capstyle='butt', zorder=4)

# twist edge, Theorem-L end, panel seam
for ax, y in ((axU, yu), (axD, yd)):
    ax.plot(x_edge(y), y, color=INK, lw=0.8, zorder=6)
axU.plot([JL_MAX, JL_MAX], [y_gap(JL_MAX), YU1], color=INK2, lw=0.5, ls=(0, (2, 1.5)), zorder=6)
for s in ('right',):
    axU.spines[s].set_visible(False); axD.spines[s].set_visible(False)
axU.spines['top'].set_visible(False)
axU.spines['bottom'].set_color(INK2); axU.spines['bottom'].set_linewidth(0.6); axU.spines['bottom'].set_linestyle((0, (1, 1.5)))
axD.spines['top'].set_visible(False)
axU.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

# x ticks at J + 1/2
JT = [0, 10, 100, 1000, 10**4, 10**5, 10**6]
axD.set_xticks([j + 0.5 for j in JT])
axD.set_xticklabels(['0', '10', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$'])
axD.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axD.set_xlabel(r'spin $J$ (log scale)', labelpad=1.5)
fig.text(0.018, 0.60, r'mass $m/\Lambda$', rotation=90, va='center', ha='center', fontsize=7.2)

# star: the one state where Lambda_amp was measured
m_star = np.sqrt((1100 - 100 - 2) * (1100 + 100 - 1)) / LAM
axD.plot([100.5], [m_star], marker='*', ms=6.5, mfc='white', mec=INK, mew=0.6, zorder=7)

# ---------------------------------------------------------------- direct labels
kw = dict(fontsize=6.2, ha='center', va='center', zorder=8)
axD.text(112, 1.25, 'Thm W', color='white', **kw)
axD.annotate(r'$(1100,100)$: premise' '\n' r'measured here only', xy=(100.5, m_star), xytext=(2.5e3, 1.055),
             fontsize=5.9, color=INK, ha='left', va='center', zorder=8,
             arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5, shrinkA=0, shrinkB=3))
axD.text(9.0e4, 1.20, r'no states ($\Delta-J<\Delta_g$)', color=INK2, fontsize=6.0, ha='center', va='center', zorder=8)
axU.text(2.0, 7.0, r'$b<5$', color='white', rotation=0, **kw)
axU.text(4.0e4, 14.0, 'Thm L', color='white', **kw)
axU.text(5.0e5, 1.75, 'no states', color=INK2, fontsize=6.0, ha='center', va='center', zorder=8)
axU.text(3.0e2, 9.0, r'hypothesis only' '\n' r'($H_R^\prime$, $H_{\rm LJ}$)', color=INK, fontsize=6.0, ha='center',
         va='center', zorder=8, bbox=dict(boxstyle='round,pad=0.15', fc=C_HYP, ec='none'))
# the near-gap band 345 <= J < 857 that neither theorem covers
axU.annotate(r'$345\leq J<857$: no theorem', xy=(560, 1.47), xytext=(9.0, 1.62), fontsize=5.9, color=INK,
             ha='center', va='center', zorder=8, bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none'),
             arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5, shrinkA=0, shrinkB=1))
axU.plot([345.5, 857.5], [1.36, 1.36], color=INK, lw=0.6, zorder=7)
for xx in (345.5, 857.5):
    axU.plot([xx, xx], [1.33, 1.39], color=INK, lw=0.6, zorder=7)

# ---------------------------------------------------------------- legend (ordinal, strongest first)
handles = [Patch(fc=C_CERT, ec='none', label='certified (flat limit)'),
           Patch(fc=C_NEPM, ec='none', label=r'error-budgeted [NE$\pm$]'),
           Patch(fc=C_THM, ec='none', label='conditional theorem (W, L)'),
           Patch(fc=C_NEM, ec='none', label='model / interpolated [NE+M], [M]'),
           Patch(fc=C_HYP, ec=C_HATCH, hatch='//////', lw=0.0, label=r'hypothesis only (not covered)')]
fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=6.0, frameon=False, bbox_to_anchor=(0.53, 0.0),
           handlelength=1.6, handleheight=0.9, columnspacing=0.9, labelspacing=0.3, borderaxespad=0.15)

out = os.path.join(HERE, 'fig2_coverage.pdf')
fig.savefig(out)
fig.savefig(os.path.join(HERE, 'fig2_coverage_preview.png'), dpi=240)
print('wrote', out, '  aspect (width/height) = %.4f' % (W_IN / H_IN),
      '  twist edge 1.3 at J = %.1f; b = 1.0421 Lambda at the gap: J = %.1f' % (x_edge(1.3) - 0.5,
      (lambda y: XL_COEF * y)(1.6480) - 0.5), '  star m/Lambda = %.4f' % m_star)

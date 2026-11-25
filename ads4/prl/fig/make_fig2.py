"""make_fig2.py -- Fig. 2 of the Letter (v14: BATCH5_FOLD Sec. 0 -- J = 17, 18, 19, 20 (twist edge to 1.05; 17, 19 held above the
twist edge, 20 at 1.05 on the lambda-level budget) and J = 25 (twist edge to 1.05) error-controlled; J = 31 twist edge on the signed
window only (pale, like J = 29). The even J = 22-32 twist edges stay [NE+M]; the partial [NE+-] mass ranges of J = 22, 24, 26, 28, 30
(held or lambda-level, Supplement Table S11) are deliberately NOT drawn, so the dense bars at J >= 17 do not suggest coverage of
those twist edges. v12: J = 9, 10, 11, 13, 14, 15 error-controlled from the twist edge to 1.05 Lambda on
own direct rows, J = 23 from its twist edge (thin, L = 1.18), J = 31 from 1.04 Lambda; NEARGAP_FILL Sec. 0, 2.4, 3. v9; v7 layout; v9: Theorem W drawn as its 47 sampled states over a hatched
"interpolated" band, J = 29 twist edge pale (signed window only), REFEREE_REPORT_4 M4.4, m8): coverage of positivity for Omega^EE at Delta_g = 1000 in the (J, m/Lambda) plane,
shaded by strength of evidence. No computation of physics: every boundary below is transcribed from the deliverables named here.

Plane.  x = J + 1/2 (log scale; J = 0 sits at x = 1/2), y = m/Lambda, m^2 = (Delta-J-2)(Delta+J-1),
Lambda = sqrt((Dg-2)(Dg-1)), Dg = 1000.  Physical states lie on or above the twist edge
    m_gap(J)/Lambda = sqrt(1 + 2J/(Dg-1))         (supplement Sec. S2.A; crosses 1.3 between J = 344 and 345).
Impact parameter b = (2J+1) Lambda / m = 2 Lambda sinh(beta/2), i.e. b < B  <=>  x < (B/2) y.

Grades (ordinal, strongest first), with sources:
  certified   flat limit m -> infinity at fixed b: Omega^EE >= N/5 on every heavy state [certified]          CERT_CLOSE; Table S5
  [NE+-]      error budget with x3 on every model term (not a bound):
                near gap J = 0, 2 (to 1.2 Lambda; POS_ERRORS Sec. 3.1), J = 4, 5, 6, 8 (to 1.05; NEARGAP2, J = 8 on the
                lambda-level budget), J = 7 (twist edge to 1.05; EXT_FOLD Sec. 0), J = 12, 16 (to 1.05; HIGHSPIN Sec. 4,
                SMALLNU_WINDOW Sec. 3.2), J = 21, 29 (to 1.10; SMALLNU_WINDOW Sec. 0, 3.1, 3.3; thin at the twist edge),
                J = 27 (twist edge to 1.10; EXT_FOLD), J = 23 (twist edge to 1.10; the twist edge thin, L = 1.18, after the
                representation term is evaluated on the windowed lambda_c, NEARGAP_FILL Sec. 3) and J = 31 (1.04-1.10 only;
                EXT_FOLD incl. its batch-2 addendum, NEARGAP_FILL Sec. 3: its twist edge stays [NE+M]);
                v12: J = 9, 10, 11, 13, 14, 15 (twist edge to 1.05; NEARGAP_FILL Sec. 2.4, own direct rows); strip m >= 1.3 Lambda, b < 5 (POS_ERRORS Sec. 2.4; one bulk model input tolerated x76;
                computed to 20 Lambda)
  theorem     conditional theorems (implications whose premises include the uniform-rate assumption (AG_U) [A] and, for W,
              the graviton edge integral (EDGE-R) [M]; AG_PROOF Sec. 5, 8; EDGE_R Sec. 0):
                W (chi-form, Lemma P4W-chi): near gap, m <= 1.3 Lambda, sampled states with 26 <= J <= 344 (47 sampled (J, m)
                   states, even J = 26-50, 64, 80, 100, 120, 150, 200, 250, 300, 344, odd 33, 51, 101, 201, 301; P2_CONTRACT Sec. 4,
                   AG_PROOF Sec. 5.2: sampled J* = 26 given (EDGE-R). v8: the v7 dashed line at J = 42 (Omega^U +- transfer under
                   the v7 premise (AG')) is removed: that alternative is superseded)
                L' : b >= 0.37211 Lambda (J + 1/2 >= 0.186055 m), J <= 10^6 (LRINF_TIGHT Sec. 3, 5; AG_PROOF Sec. 6.2: uniform rate);
                   J*(m) = ceil((0.37211 m - 1)/2) = 223 at the twist edge; overlaps W at the gap on 223 <= J <= 344
  [NE+M]/[M]  computed inside models or interpolated: near-gap J <= 25 not error-budgeted (chain spins, masses to 5 Lambda;
                v14: the twist edges of even J = 22-32, J = 28 at 1.03, J = 30 at 1.03-1.04 and J = 32 above its twist edge)
                and J = 26-32 at masses outside W's samples (drawn under W), even J = 36, 40, 48 above 1.3 Lambda;
                m >= 1.3 Lambda, 5 <= b < b_A(m) (POS_ERRORS Sec. 2.4: b_A = 8.5, 8.7, 10.5, 13, 15, 15.3, 13.8
                at m/Lambda = 1.3, 1.5, 2, 3, 5, 10, 20)
  hypothesis  no computation and no theorem: covered only by H_R' [D-sketch + M] -- m >= 1.3 Lambda with
                b_A(m) <= b < 0.37211 Lambda, and J > 10^6 (off scale) (LARGESPIN_ANALYTIC Sec. 5.2; Table S11)
  v7: the v6 bracket "345 <= J < 857: no theorem" is removed (that band lies inside Theorem L', LRINF_TIGHT Sec. 6), and the
  v6 star at (1100, 100) is removed (the O(1/m^2) premise is now computed at all 47 W states, P2_CONTRACT).
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


XL_COEF = 0.37211 / 2 * LAM      # Theorem L': 2J + 1 >= 0.37211 m  (LRINF_TIGHT Sec. 3)
XMIN, XMAX = 0.4, 4.0e6
JL_MAX = 1.0e6 + 0.5

NEPM = [(0, None, 1.2), (2, None, 1.2), (4, None, 1.05), (5, None, 1.05), (6, None, 1.05), (7, None, 1.05), (8, None, 1.05),
        (9, None, 1.05), (10, None, 1.05), (11, None, 1.05), (12, None, 1.05), (13, None, 1.05), (14, None, 1.05),
        (15, None, 1.05), (16, None, 1.05), (17, None, 1.05), (18, None, 1.05), (19, None, 1.05), (20, None, 1.05),
        (21, None, 1.10), (23, None, 1.10), (25, None, 1.05), (27, None, 1.10), (29, 1.04, 1.10),
        (31, 1.04, 1.10)]  # (J, lowest budgeted mass or None = twist edge, highest); v14: BATCH5_FOLD adds J = 17-20, 25 and the J = 31 twist edge (signed)
NE_SIGNED = [(29, 1.04), (31, 1.04)]  # J = 29, 31: twist edge to 1.04 only with the signed window (sign-free L = 0.94 / 0.90; SMALLNU_WINDOW, BATCH5_FOLD)
# Theorem W's 47 sampled states (J, m/Lambda), from code/agproof/cache/chi.txt (AG_PROOF Sec. 5.2)
W_STATES = [(J, m) for J, ms in [(26, (1.026, 1.128, 1.231)), (30, (1.030, 1.133, 1.236)), (34, (1.033, 1.137, 1.240)),
            (38, (1.037, 1.141, 1.245)), (42, (1.041, 1.145, 1.249)), (46, (1.045, 1.150, 1.254)),
            (50, (1.049, 1.154, 1.259)), (64, (1.062, 1.168, 1.275)), (80, (1.077, 1.185, 1.293)), (100, (1.096, 1.205)),
            (120, (1.114, 1.225)), (150, (1.140, 1.254)), (200, (1.183,)), (250, (1.225,)), (300, (1.265,)), (344, (1.299,)),
            (33, (1.033, 1.136, 1.239)), (51, (1.050, 1.155, 1.260)), (101, (1.096, 1.206)), (201, (1.184,)), (301, (1.266,))]
            for m in ms]
assert len(W_STATES) == 47
JW_LO = 26                         # Theorem W: sampled J* given (EDGE-R) (AG_PROOF Sec. 5.2; v7's alternative 42 superseded)
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


def fbx(ax, y, xl, xr, color, hatch=None, z=1, hc=C_HATCH):
    xl = np.broadcast_to(xl, y.shape).astype(float); xr = np.broadcast_to(xr, y.shape).astype(float)
    ok = xr > xl
    if not ok.any():
        return
    kw = dict(color=color, lw=0, zorder=z, where=ok, interpolate=False)
    if hatch:
        ax.fill_betweenx(y, xl, np.where(ok, xr, xl), facecolor=color, edgecolor=hc, hatch=hatch, lw=0, zorder=z)
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
fbx(axD, yd, JW_LO, np.minimum(xed, 345.0), SURF, hatch='\\\\\\', z=2.4, hc=C_NEM)                  # Theorem W: between samples, interpolated
fbx(axD, yd, XL_COEF * yd, xed, C_THM, z=2.5)                                        # Theorem L' (J >= 223 at the gap)
axD.set_ylim(YD0, YD1); axD.set_yticks([1.0, 1.1, 1.2]); axD.set_yticklabels(['1.0', '1.1', '1.2'])
for J in NEM_TICKS:                                                                # same spins, chain masses to 5 Lambda
    axU.plot([J + 0.5] * 2, [YU0, 5.0], color=C_NEM, lw=1.1, solid_capstyle='butt', zorder=3)
for J, ybot, ytop in NEPM:                                                         # [NE+-] near-gap budgets
    y0 = max(y_gap(J + 0.5), YD0) if ybot is None else ybot
    axD.plot([J + 0.5] * 2, [y0, ytop], color=C_NEPM, lw=1.6, solid_capstyle='butt', zorder=4)
for J, ytop in NE_SIGNED:                                                          # signed window only
    axD.plot([J + 0.5] * 2, [y_gap(J + 0.5), ytop], color=C_NEM, lw=1.6, solid_capstyle='butt', zorder=4)
axD.scatter([J + 0.5 for J, m in W_STATES], [m for J, m in W_STATES], s=3.0, color=C_THM, edgecolors='none', zorder=4.5)
# Theorem L' boundary in the near-gap panel
yl = np.linspace(y_gap(223.5), YD1, 50)
axD.plot(XL_COEF * yl, yl, color='white', lw=0.7, ls=(0, (0.8, 1.0)), zorder=5)


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


# ---------------------------------------------------------------- direct labels
kw = dict(fontsize=6.2, ha='center', va='center', zorder=8)
axD.text(76, 1.2, 'theorem 1\n(sampled\nstates)', color=INK, fontsize=5.4, ha='center', va='center', zorder=8,
         bbox=dict(boxstyle='round,pad=0.12', fc=SURF, ec='none'))
axD.text(9.0e4, 1.20, r'no states ($\Delta-J<\Delta_g$)', color=INK2, fontsize=6.0, ha='center', va='center', zorder=8)
axU.text(2.0, 7.0, r'$b<5$', color='white', rotation=0, **kw)
axU.text(3.0e4, 14.0, 'theorem 2', color='white', **kw)
axU.text(5.0e5, 1.75, 'no states', color=INK2, fontsize=6.0, ha='center', va='center', zorder=8)
axU.text(330.0, 9.0, r'hypothesis' '\n' r'only', color=INK, fontsize=6.0, ha='center',
         va='center', zorder=8, bbox=dict(boxstyle='round,pad=0.15', fc=C_HYP, ec='none'))
# ---------------------------------------------------------------- legend (ordinal, strongest first)
handles = [Patch(fc=C_CERT, ec='none', label='certified (flat limit)'),
           Patch(fc=C_NEPM, ec='none', label='numerics with error control'),
           Patch(fc=C_THM, ec='none', label='conditional theorem'),
           Patch(fc=SURF, ec=C_NEM, hatch='\\\\\\', lw=0.0, label='between theorem-1 samples'),
           Patch(fc=C_NEM, ec='none', label='model / interpolated'),
           Patch(fc=C_HYP, ec=C_HATCH, hatch='//////', lw=0.0, label='hypothesis only')]
fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=6.0, frameon=False, bbox_to_anchor=(0.53, 0.0),
           handlelength=1.6, handleheight=0.9, columnspacing=0.9, labelspacing=0.3, borderaxespad=0.15)

out = os.path.join(HERE, 'fig2_coverage.pdf')
fig.savefig(out)
fig.savefig(os.path.join(HERE, 'fig2_coverage_preview.png'), dpi=240)
jl = next(J for J in range(2000) if 2 * J + 1 >= 0.37211 * LAM * y_gap(J + 0.5))
print('wrote', out, '  aspect (width/height) = %.4f' % (W_IN / H_IN),
      '  twist edge 1.3 at J = %.1f; Theorem L\' J* at the twist edge = %d (expect 223)' % (x_edge(1.3) - 0.5, jl))

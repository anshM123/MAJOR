"""Fig. 1: the two chiral quotient surfaces of a rotating three-boundary wormhole and one cut class.

Schematic. Each panel shows a hyperbolic pair of pants with three funnels: the outer annulus and the two
inner annuli are the funnels (asymptotic regions), dashed circles are the horizons (closed geodesics
bounding the convex core), thick circles are the conformal boundaries. The region A consists of an
interval of boundary 1 (between the marked points p, p') and the whole of boundary 2. The same cut class c
(an arc from p to p' going around funnel 2) is drawn on X_L and on X_R, whose geometries differ.
Usage: python3 fig1_wormhole.py <output.pdf>
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from scipy.interpolate import splprep, splev
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
style.apply()
out = sys.argv[1] if len(sys.argv) > 1 else "fig1_wormhole.pdf"

FUNNEL = "#e9eef5"
CORE = "#ffffff"

def smooth(pts, n=300):
    pts = np.array(pts).T
    tck, _ = splprep(pts, s=0, k=3)
    return splev(np.linspace(0, 1, n), tck)

def panel(ax, holes, cut_pts, label, ell, theta=(58, 122)):
    R_bdy, R_hor = 1.0, 0.80
    ax.add_patch(Circle((0, 0), R_bdy, color=FUNNEL, zorder=0, lw=0))
    ax.add_patch(Circle((0, 0), R_hor, color=CORE, zorder=1, lw=0))
    # region A on the outer boundary: interval between the marked points
    t0, t1 = theta
    ax.add_patch(Wedge((0, 0), R_bdy, t0, t1, width=R_bdy - R_hor, color=style.BLUE, alpha=0.13, lw=0, zorder=1))
    ax.add_patch(Circle((0, 0), R_bdy, fill=False, ec=style.INK, lw=1.4, zorder=4))
    ax.plot(*[R_bdy * f(np.radians(np.linspace(t0, t1, 60))) for f in (np.cos, np.sin)], color=style.BLUE, lw=2.6, zorder=5,
            solid_capstyle="butt")
    ax.add_patch(Circle((0, 0), R_hor, fill=False, ec=style.INK2, lw=0.8, ls=(0, (3, 2)), zorder=3))
    for k, (cx, cy, rh, rb) in enumerate(holes):
        inA = (k == 0)
        ax.add_patch(Circle((cx, cy), rh, color=FUNNEL, zorder=2, lw=0))
        if inA:
            ax.add_patch(Circle((cx, cy), rh, color=style.BLUE, alpha=0.13, zorder=2, lw=0))
        ax.add_patch(Circle((cx, cy), rb, color="white", zorder=2, lw=0))
        ax.add_patch(Circle((cx, cy), rh, fill=False, ec=style.INK2, lw=0.8, ls=(0, (3, 2)), zorder=3))
        ax.add_patch(Circle((cx, cy), rb, fill=False, ec=style.BLUE if inA else style.INK, lw=2.2 if inA else 1.4, zorder=4))
        ax.text(cx, cy, str(k + 2), ha="center", va="center", fontsize=7, color=style.INK2)
    # marked points
    for t, nm, dx in ((t0, r"$p'$", 0.06), (t1, r"$p$", -0.06)):
        x, y = np.cos(np.radians(t)), np.sin(np.radians(t))
        ax.plot([x], [y], "o", ms=3.2, color=style.INK, zorder=6)
        ax.text(1.12 * x + dx, 1.10 * y, nm, ha="center", va="center", fontsize=8)
    # the cut
    xs, ys = smooth(cut_pts)
    # homology region: between the cut and the boundary interval (funnel 2 lies inside it)
    tt = np.radians(np.linspace(theta[0], theta[1], 80))
    ax.fill(np.concatenate([xs, np.cos(tt)]), np.concatenate([ys, np.sin(tt)]), color=style.BLUE, alpha=0.10,
            lw=0, zorder=1.5)
    ax.plot(xs, ys, color=style.RED, lw=1.5, zorder=5)
    ax.text(*ell[0], ell[1], color=style.RED, fontsize=8, ha="center", va="center")
    ax.text(-1.04, -0.98, label, fontsize=9, ha="left", va="bottom")
    ax.text(0.0, -0.9, "1", ha="center", va="center", fontsize=7, color=style.INK2)
    ax.set_xlim(-1.08, 1.08); ax.set_ylim(-1.08, 1.12)
    ax.set_aspect("equal"); ax.axis("off")

fig, axs = plt.subplots(1, 2, figsize=(3.4, 1.85))
p = (np.cos(np.radians(122)), np.sin(np.radians(122))); pp = (np.cos(np.radians(58)), np.sin(np.radians(58)))
holesL = [(-0.30, -0.05, 0.25, 0.13), (0.36, -0.12, 0.20, 0.10)]
cutL = [p, (-0.42, 0.55), (-0.62, 0.05), (-0.40, -0.40), (-0.05, -0.35), (0.06, 0.15), (0.30, 0.55), pp]
holesR = [(-0.18, 0.08, 0.17, 0.09), (0.34, -0.22, 0.27, 0.14)]
cutR = [p, (-0.34, 0.58), (-0.43, 0.10), (-0.24, -0.16), (0.02, -0.05), (0.12, 0.30), (0.33, 0.58), pp]
panel(axs[0], holesL, cutL, r"$X_L$", ((-0.53, -0.45), r"$\ell_L(c)$"))
panel(axs[1], holesR, cutR, r"$X_R$", ((-0.30, -0.36), r"$\ell_R(c)$"))
fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.04)
fig.savefig(out)

"""Fig. 2: mechanism.
(a) Cut-and-paste on one arrangement of minimal cuts (strong subadditivity, local picture on the quotient
    surface, Klein model so that geodesics are straight).  Cells carry bit strings (x_AB, x_BC); the outputs
    delta_B and delta_ABC are boundaries of unions of cells, resolved at the vertex.  A test geodesic g
    crosses the inputs at least as often as the outputs.
(b) The endpoint class tau_p in the half-plane chart with the marked point p at infinity.  Every geodesic of
    tau_p has its feet on both sides of [q_-, q_+], which contains everything else on the ideal boundary; it
    crosses the cut arc ending at p exactly once.
(c) The resulting graph: a finite cell graph realising G, plus one edge of weight t_p(eps) per entangling
    point p joining the parties adjacent at p.
Usage: python3 fig2_mechanism.py <output.pdf>
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
style.apply()
out = sys.argv[1] if len(sys.argv) > 1 else "fig2_mechanism.pdf"

PARTY = {"A": style.BLUE, "B": style.VIOLET, "C": style.GREEN, "O": style.GRAY}
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.35), gridspec_kw=dict(width_ratios=[1.0, 1.25, 1.0]))

# ---------------- (a) ----------------
ax = axs[0]
def pt(deg, r=1.0): return np.array([r * np.cos(np.radians(deg)), r * np.sin(np.radians(deg))])
ax.add_patch(Circle((0, 0), 1, color="#f6f8fb", lw=0, zorder=0))
arcs = {"A": (135, 225), "B": (45, 135), "C": (-45, 45), "O": (225, 315)}
for nm, (t0, t1) in arcs.items():
    tt = np.radians(np.linspace(t0, t1, 60))
    ax.plot(np.cos(tt), np.sin(tt), color=PARTY[nm], lw=2.6, solid_capstyle="butt", zorder=3)
    m = pt((t0 + t1) / 2, 1.13)
    ax.text(*m, nm, ha="center", va="center", fontsize=8, color=style.INK)
for d in (45, 135, 225, 315):
    ax.plot(*pt(d), "o", ms=2.6, color=style.INK, zorder=4)
# a funnel (hole) of the quotient surface in the O-cell
ax.add_patch(Circle((0.0, -0.58), 0.13, color="white", ec=style.INK2, lw=0.8, ls=(0, (3, 2)), zorder=2))
ax.text(0.0, -0.58, r"$\cdot$", ha="center", va="center", color=style.INK2)
# inputs: gamma_AB from 225 to 45, gamma_BC from 135 to 315 (straight in the Klein model)
for (d0, d1), lab, off in (((225, 45), r"$\gamma_{AB}$", (0.42, 0.24)), ((135, 315), r"$\gamma_{BC}$", (0.44, -0.30))):
    P, Q = pt(d0), pt(d1)
    ax.plot([P[0], Q[0]], [P[1], Q[1]], color=style.INK, lw=1.0, zorder=3)
    ax.text(*off, lab, fontsize=7.5, ha="left", va="center")
# outputs: resolve at the vertex (origin)
def resolved(d0, d1, corner, k=60):
    """boundary of a union of cells, pushed off the vertex: P -> corner -> Q with a rounded corner"""
    P, Q, C = pt(d0, 0.995), pt(d1, 0.995), np.array(corner)
    s = np.linspace(0, 1, k)
    a = np.outer(1 - s, P) + np.outer(s, C)
    b = np.outer(1 - s, C) + np.outer(s, Q)
    pts = np.vstack([a[:-8], b[8:]])
    from scipy.interpolate import splprep, splev
    tck, _ = splprep(pts.T, s=0.0008, k=3)
    return np.array(splev(np.linspace(0, 1, 300), tck)).T
dB = resolved(135, 45, (0.0, 0.16))
dABC = resolved(225, 315, (0.0, -0.16))
ax.plot(dB[:, 0], dB[:, 1], color=style.RED, lw=1.4, ls=(0, (3, 1.3)), zorder=4)
ax.plot(dABC[:, 0], dABC[:, 1], color=style.ORANGE, lw=1.4, ls=(0, (3, 1.3)), zorder=4)
ax.text(0.27, 0.47, r"$\delta_{B}$", color=style.INK, fontsize=7.5, ha="center")
ax.text(0.0, -0.40, r"$\delta_{ABC}$", color=style.INK, fontsize=7.5, ha="center", va="center")
for xy, b in (((-0.62, 0.0), "10"), ((0.0, 0.62), "11"), ((0.62, 0.0), "01"), ((0.33, -0.80), "00")):
    ax.text(*xy, b, fontsize=6.5, color=style.INK2, ha="center", va="center", family="monospace")
# a test geodesic
P, Q = pt(255), pt(100)
ax.plot([P[0], Q[0]], [P[1], Q[1]], color=style.INK2, lw=0.9, ls=(0, (1, 1.2)), zorder=3)
ax.text(-0.33, 0.80, r"$g$", fontsize=7.5, color=style.INK2)
ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(a)", loc="left", fontsize=8, pad=1)

# ---------------- (b) ----------------
ax = axs[1]
qm, qp, x0, Y = -1.0, 1.3, 0.35, 3.9
ax.axhspan(Y, 5.0, color="#eeeeee", lw=0, zorder=0)
ax.text(-5.45, Y + 0.62, r"$B_p(\epsilon)=\{y>1/\epsilon\}$", fontsize=6.5, color=style.INK2, va="center")
ax.axhline(Y, color=style.INK2, lw=0.6, zorder=1)
ax.plot([-5.6, qm], [0, 0], color=PARTY["A"], lw=2.6, solid_capstyle="butt", zorder=3)
ax.plot([qp, 5.6], [0, 0], color=PARTY["B"], lw=2.6, solid_capstyle="butt", zorder=3)
ax.plot([qm, qp], [0, 0], color=style.INK, lw=0.8, zorder=3)
rng = np.random.default_rng(3)
dust = np.sort(rng.uniform(qm, qp, 40))
ax.plot(dust, np.zeros_like(dust), "|", ms=4, color=style.INK, mew=0.6, zorder=4)
ax.text(-3.3, -0.42, r"party $i$", fontsize=7, ha="center", va="center")
ax.text(3.45, -0.42, r"party $j$", fontsize=7, ha="center", va="center")
ax.annotate("", xy=(qm, -0.55), xytext=(qp, -0.55), arrowprops=dict(arrowstyle="<->", lw=0.6, color=style.INK2))
ax.text((qm + qp) / 2, -0.85, r"$d_p=q_+-q_-$", fontsize=6.5, ha="center", va="center", color=style.INK2)
for (a, b) in ((-1.6, 1.9), (-2.6, 2.2), (-1.3, 3.6), (-3.9, 3.3), (-2.2, 4.9), (-4.6, 4.6)):
    c, r = (a + b) / 2, (b - a) / 2
    th = np.linspace(0, np.pi, 200)
    xx, yy = c + r * np.cos(th), r * np.sin(th)
    ax.plot(xx, np.minimum(yy, 10), color=style.INK2, lw=0.7, alpha=0.85, zorder=2)
ax.text(-4.2, 2.65, r"$\tau_p$", fontsize=8, color=style.INK2)
ax.plot([x0, x0], [0, Y], color=style.RED, lw=1.5, zorder=5)
ax.plot([x0, x0], [Y, 5.0], color=style.RED, lw=1.0, ls=(0, (1.5, 1.5)), zorder=5)
ax.text(x0 + 0.15, 1.05, r"$c$", color=style.INK, fontsize=8)
ax.text(x0 + 0.15, 4.65, r"$\to p$", color=style.INK, fontsize=7)
ax.text(1.9, 2.95, r"$m\{g\in\tau_p:\,g\cap c^{\epsilon}\neq\emptyset\}$" "\n" r"$=\log\frac{1}{\epsilon\, d_p}+o(1)$",
        fontsize=6.5, ha="left", va="center", color=style.INK,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
ax.set_xlim(-5.6, 5.6); ax.set_ylim(-1.05, 5.0); ax.axis("off")
ax.set_title("(b)", loc="left", fontsize=8, pad=1)

# ---------------- (c) ----------------
ax = axs[2]
term = {"A": (-0.95, 0.0), "B": (0.0, 0.95), "C": (0.95, 0.0), "O": (0.0, -0.95)}
internal = {"u": (-0.30, 0.22), "v": (0.32, 0.18), "w": (0.0, -0.35)}
edges = [("A", "u", 2.2), ("B", "u", 1.2), ("B", "v", 1.6), ("C", "v", 2.0), ("u", "v", 0.8), ("u", "w", 1.0),
         ("v", "w", 1.0), ("O", "w", 2.4), ("A", "w", 0.6)]
pos = {**term, **internal}
for a, b, w in edges:
    ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color=style.INK2, lw=0.5 + 0.7 * w, zorder=2,
            solid_capstyle="round")
# Bell-pair edges: one per entangling point, joining the two parties adjacent there
bell = [("A", "B", 0.55), ("B", "C", 0.55), ("C", "O", 0.55), ("O", "A", 0.55)]
for a, b, bend in bell:
    P, Q = np.array(pos[a]), np.array(pos[b])
    mid = (P + Q) / 2; n = mid / np.linalg.norm(mid)
    s = np.linspace(0, 1, 100)
    curve = np.outer(1 - s, P) + np.outer(s, Q) + np.outer(bend * np.sin(np.pi * s), n)
    ax.plot(curve[:, 0], curve[:, 1], color=style.ORANGE, lw=1.3, ls=(0, (3, 1.5)), zorder=1)
    lab = mid + (bend + 0.12) * n
    ax.text(*lab, r"$t_p$", fontsize=7, ha="center", va="center", color=style.INK)
for nm, (x, y) in term.items():
    ax.add_patch(FancyBboxPatch((x - 0.11, y - 0.11), 0.22, 0.22, boxstyle="round,pad=0.0,rounding_size=0.04",
                                fc=PARTY[nm], ec="white", lw=1.0, zorder=4))
    ax.text(x, y, nm, fontsize=7, color="white", ha="center", va="center", zorder=5)
for nm, (x, y) in internal.items():
    ax.add_patch(Circle((x, y), 0.055, fc="white", ec=style.INK, lw=0.8, zorder=4))
ax.text(0.0, -1.42, r"$4G_N S^{\epsilon}=G+\sum_p t_p(\epsilon)\,b_p$", fontsize=7.5, ha="center", va="center")
ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.6, 1.5); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(c)", loc="left", fontsize=8, pad=1)

fig.subplots_adjust(left=0.005, right=0.995, top=0.93, bottom=0.01, wspace=0.06)
fig.savefig(out)

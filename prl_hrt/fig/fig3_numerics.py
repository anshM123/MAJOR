"""Fig. 3: minimum facet value of the holographic entropy cone over random configurations.

Reads the logs written by the Schottky-surface test (two hyperbolic metrics with a common marking,
interval and whole-boundary parties) and plots, for every configuration, the smallest value of
Q(S)/max|S| over all facets Q of C_N (N = 3, 4, 5: 7, 20, 372 facets) for
  G        the cut-off independent part of the decomposition F = G + sum_p t_p b_p,
  F^eps    the regulated functional itself at eps = 1e-5,
  ctl 1    G - (W/2) sum_p b_p   (wrong subtraction of the endpoint terms),
  ctl 2    min_c max_j l_j(c)    (not a positive combination of lengths).
Usage: python3 fig3_numerics.py <log directory> <output.pdf>
"""
import re, sys, os, glob
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style
style.apply()

logdir = sys.argv[1] if len(sys.argv) > 1 else "../../fullcone"
out = sys.argv[2] if len(sys.argv) > 2 else "fig3_numerics.pdf"

runs = [("pants", 3, "run_pants_N3.log"), ("pants", 4, "run_pants_N4.log"), ("pants", 5, "run_pants_N5.log"),
        ("pants", 4, "run_hard1_pants_N4.log"),
        ("torus", 3, "run_torus_N3.log"), ("torus", 4, "run_torus_N4.log"), ("torus", 5, "run_torus_N5.log"),
        ("torus", 4, "run_hard2_torus_N4.log")]
pat = re.compile(r"facet\(F\)=([-+0-9.e]+) facet\(G\)=([-+0-9.e]+) ctl=([-+0-9.e]+) max=([-+0-9.e]+)")

fig, ax = plt.subplots(figsize=(3.4, 2.55))
x0 = 0
ticks, labels, bounds = [], [], []
total = 0
allG, allF = [], []
for k, (top, N, fn) in enumerate(runs):
    rows = [tuple(map(float, m.groups())) for m in map(pat.search, open(os.path.join(logdir, fn))) if m]
    rows = np.array(rows)
    n = len(rows); total += n
    xs = x0 + np.arange(n)
    F, G, C1, C2 = rows.T
    allG += list(G); allF += list(F)
    kw = dict(s=7, linewidths=0, zorder=3)
    ax.scatter(xs, C2, marker="s", color=style.AQUA, label=r"control: $\min_c\max_j\ell_j$" if k == 0 else None, **kw)
    ax.scatter(xs, C1, marker="v", color=style.ORANGE, label=r"control: $G-\frac{W}{2}\sum_p b_p$" if k == 0 else None, **kw)
    ax.scatter(xs, F, marker="o", facecolors="none", edgecolors=style.GRAY, s=9, linewidths=0.5, zorder=3,
               label=r"$F^{\epsilon}$, $\epsilon=10^{-5}$" if k == 0 else None)
    ax.scatter(xs, G, marker="o", color=style.BLUE, s=5, linewidths=0, zorder=4,
               label=r"$G$ (cut-off independent part)" if k == 0 else None)
    bad = np.where(F < -1e-9)[0]
    for i in bad:
        ax.scatter([xs[i]], [F[i]], marker="o", facecolors="none", edgecolors=style.INK, s=22, linewidths=0.8, zorder=5)
        ax.annotate(r"$F^{\epsilon}$ at $\epsilon>\epsilon_1$ ($t_p<0$): SA" "\n" r"fails; holds at $\epsilon=10^{-12}$",
                    xy=(xs[i], F[i]), xytext=(xs[i] - 40, -3e-7), fontsize=6, color=style.INK, ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=style.INK, shrinkB=3))
    ticks.append(x0 + n / 2 - 0.5)
    labels.append(("P" if top == "pants" else "T") + f"{N}" + ("*" if "hard" in fn else ""))
    x0 += n + 6
    bounds.append(x0 - 3)

ax.set_yscale("symlog", linthresh=1e-14, linscale=0.35)
ax.axhspan(-1e-13, 1e-13, color="#eeeeee", zorder=0, lw=0)
ax.axhline(0, color=style.INK2, lw=0.5, zorder=1)
for b in bounds[:-1]:
    ax.axvline(b, color=style.GRID, lw=0.5, zorder=0)
ax.set_xticks(ticks); ax.set_xticklabels(labels)
ax.set_xlim(-3, x0 - 3)
ax.set_ylim(-3, 1e-13)
ax.set_yticks([-1, -1e-2, -1e-4, -1e-6, -1e-8, -1e-10, -1e-12, 0])
ax.set_yticklabels([r"$-1$", r"$-10^{-2}$", r"$-10^{-4}$", r"$-10^{-6}$", r"$-10^{-8}$", r"$-10^{-10}$", r"$-10^{-12}$", r"$0$"])
ax.set_ylabel(r"$\min_{Q}\;Q(\cdot)/\max_A|S_A|$")
ax.text(x0 - 5, -2.5e-13, r"$|\cdot|<10^{-13}$", fontsize=6, color=style.INK2, ha="right", va="top")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="x", length=0)
leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
                handletextpad=0.2, columnspacing=0.8, borderpad=0.1, labelspacing=0.2, markerscale=1.6)
fig.tight_layout(pad=0.2)
fig.savefig(out)
print(f"configurations: {total}; min G facet {min(allG):.2e}; F violations {sum(f < -1e-9 for f in allF)}")

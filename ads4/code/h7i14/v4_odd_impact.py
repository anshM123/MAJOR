"""v4 (H7_I14_LIGHTSIDE.md §23.3): parity-odd sector (odd J, +-) near-gap test at Delta_g = 300, 1000 with C^_2' (spinning de-shift), and the
combined even + odd n.  Read-only imports of spinfunc/p4_odd_analysis, p6_odd_neargap (nothing written there).

Exactly the PARITY_ODD_ROUTEA.md §6 protocol (p4 grid: K2 linear chain, Delta_ref = 2e4, a_{+-} from C^_2 on x <= 0.25, nu >= 15, de-shift;
p6: l9_neargap conventions, J = 5, 7, 9 covered, odd J >= 11 base), with ONE change, as v3_heavy_impact does for the even sector: in the
k = 2 de-shift the scalar coefficients g/(g + g~), g~/(g + g~) -> (1 + u_T)/2, (1 - u_T)/2, u_T = (h_- - h_+)/(h_- + h_+) the EXACT rational
function of v4_hclosed (§23.1; agrees with the v3 fits to 2e-11).  k = 3 (C^_3) keeps the scalar g (its shift relation is open).
No shadow coefficient is touched (the §12 branch rule is irrelevant here: the p3 caches already strip the phase on the even-integer rows used).
Runs: baseline (must reproduce 0.2072 / 0.2023), C^_2' everywhere, and C^_2' only at nu >= 15 (diagnostic).
Combined n: the even-sector C^_2' numbers are those of v3_heavy_impact (cache/v3_heavy_impact.txt), recomputed here with the exact u_T
(l9_neargap through v3_heavy_impact.run with the exact shift) for consistency.
Usage: python3 v4_odd_impact.py [even]   -> cache/v4_odd_impact.txt
"""
import os, sys, io, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
sys.dont_write_bytecode = True
import numpy as np
import pickle
import p4_odd_analysis as A
import p6_odd_neargap as P6
from v4_chi import u_spin

_orig_shift = A.shift


def make_shift(nu_min=0.0):
    def shift_spin(J, h, k, xs, m):
        if k != 2:
            return _orig_shift(J, h, k, xs, m)
        base = _orig_shift(J, h, k, xs, m)
        f = lambda t: A.env(h, t, k)*A.Dc(h, J, t)
        out = []
        for x, b0 in zip(xs, base):
            if x*m < nu_min:
                out.append(b0); continue
            u = complex(u_spin(x*m))
            out.append((0.5*((1 + u)*f(x + 1j/m) + (1 - u)*f(x - 1j/m)) - f(x)).real)
        return np.array(out)
    return shift_spin


def odd_run(label, shift, pr, Dgs=(300, 1000)):
    A.shift = shift
    try:
        grids = {J: A.grid(J, K=2, kinds=('C2', 'C3H')) for J in (5, 7, 9)}
    finally:
        A.shift = _orig_shift
    P6.GRID_OVERRIDE.clear(); P6.GRID_OVERRIDE.update(grids)
    res = P6.main(list(Dgs), write=False, quiet=True)
    P6.GRID_OVERRIDE.clear()
    pr(f"## odd sector, {label}")
    for Dg in Dgs:
        r = res[Dg]
        pr(f"  Delta_g = {Dg}: n = 1: {r['ok']}; min F/N_s = {r['min']:+.4f} at (m/Lam, J, Df) = {r['wmin']}; t* = ({r['tstar'][0]:+.1f}, "
           f"{r['tstar'][1]:+.1f}); per J: " + ", ".join(f"J={J}: {r['perJ'][J]:+.4f} (m/Lam {r['perJ_arg'][J]})" for J in sorted(r['perJ'])))
    return grids, res


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(str(s) + "\n")
    t0 = time.time()
    gb, rb = odd_run("baseline (scalar g; PARITY_ODD_ROUTEA §6)", _orig_shift, pr)
    # consistency with the stored p4 grid
    old = pickle.load(open(os.path.join(SF, 'cache', 'p4_grid.pkl'), 'rb'))
    dmax = max(np.max(np.abs(gb[J][D][kd] - old[J][D][kd])) for J in gb for D in gb[J] if D in old.get(J, {}) for kd in ('C2', 'C3H'))
    pr(f"  baseline grid vs stored cache/p4_grid.pkl: max |dGamma| = {dmax:.1e}")
    gs, rs = odd_run("C^_2' (spinning de-shift, exact u_T) at all nodes", make_shift(0.0), pr)
    g15, r15 = odd_run("C^_2' only at nu >= 15", make_shift(15.0), pr)
    pr("## change of the odd-sector grid Gamma_{2,+-} (C^_2' - baseline), max over x, per row (c2 units):")
    for J in (5, 7, 9):
        pr(f"  J={J}: " + "; ".join(f"Delta={D:.0f} (m={gb[J][D]['m']:.0f}, x_min nu={gb[J][D]['xs'][0]*gb[J][D]['m']:.1f}): "
                                    f"C2 {np.max(np.abs(gs[J][D]['C2'] - gb[J][D]['C2'])):.2g}, C3H {np.max(np.abs(gs[J][D]['C3H'] - gb[J][D]['C3H'])):.2g}, "
                                    f"a {gs[J][D]['a']/gb[J][D]['a'] - 1:+.1e}" for D in sorted(gb[J])))
    if 'even' in sys.argv[1:]:
        import v3_heavy_impact as H3
        H3.pr = pr
        def shift_even(J, h, k, xs, m):
            if k != 2:
                return H3._orig_shift(J, h, k, xs, m)
            f = lambda t: H3.A.env(h, t, k)*H3.A.Dc(h, J, t)
            out = []
            for x in xs:
                u = complex(u_spin(x*m))
                out.append((0.5*((1 + u)*f(x + 1j/m) + (1 - u)*f(x - 1j/m)) - f(x)).real)
            return np.array(out)
        H3.A.shift = shift_even
        H3.run("EVEN sector, C^_2' with the exact u_T (v4); l9 Laplace rows", (300, 1000), direct=False)
        H3.A.shift = H3._orig_shift
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v4_odd_impact.txt'), 'w').write(buf.getvalue())

"""
Validation of the exact <JJJJ> partial waves W[v_a, v_b] (fourpt/spinblock2):
 (a) conservation at x_1 and at x_2 (and, by the (12)<->(34) structure, the same code path at x_3, x_4 is tested at x_3):
     |d^mu W_mu| / max_mu |d_mu W_nu| by 4th-order central differences (double precision), for the four basis combinations,
     and a control (a single word, not conserved);
 (b) pair-exchange symmetry  W[v_a, v_b](1,2,3,4) = W[v_b, v_a](3,4,1,2)  (exact symmetry of the s-channel partial wave);
 (c) the scalar limit of the machinery: the word pair (H, H) is H_12 H_34 times the scalar partial wave with external
     dimensions (3,3,3,3) -- checked against an independent evaluation of the scalar block (e2e/blocks.CBlock).
"""
import os, sys, numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from spinblock import block_coords_from_uv      # noqa: E402
from blocks import CBlock                        # noqa: E402
from spinblock2 import SpinBlock4, invariants4, monomial_values, mdot
from tmoments4 import basis_coeffs
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def W_eval(sb, combos, X, Z):
    """X: dict i -> [npts, 3]; Z: dict i -> [npts, npol, 3].  Returns {name: [npts, npol]} = W (unstripped)."""
    s, z = invariants4(X, Z)
    r, chi, A0, E0 = block_coords_from_uv(s)
    ls, C = sb.eval_coeffs(s, r, chi, A0, E0, combos)
    mv = monomial_values(sb.phi[1], z)
    strip = 2*np.log(s['s13']) + 2*np.log(s['s24'])
    return {k: np.exp(ls - strip)[:, None]*np.einsum('mp,mpk->pk', v, mv) for k, v in C.items()}

def config(eps=0.0):
    x = 0.35*np.array([np.cosh(0.7), np.sinh(0.7)*np.cos(0.3), np.sinh(0.7)*np.sin(0.3)])
    y = 0.45*np.array([np.cosh(0.2), np.sinh(0.2), 0.0])
    x2n = mdot(x, x)
    return {1: -y, 2: x/x2n, 3: y, 4: -x/x2n}

def divergence(sb, combos, pt, h=2e-3, seed=1):
    """divergence at operator pt of W (other polarisations fixed random vectors)"""
    rng = np.random.default_rng(seed)
    base = config()
    others = {i: rng.normal(size=3) for i in (1, 2, 3, 4)}
    E = np.eye(3)
    steps = [-2, -1, 1, 2]; wts = [1/12, -2/3, 2/3, -1/12]
    # positions: [mu (derivative direction), step] ; polarisation of pt: the three unit vectors
    Xs = {i: [] for i in (1, 2, 3, 4)}
    for mu in range(3):
        for st in steps:
            for i in (1, 2, 3, 4):
                Xs[i].append(base[i] + (st*h*E[mu] if i == pt else 0))
    Xd = {i: np.array(v) for i, v in Xs.items()}
    npts = Xd[1].shape[0]
    Zd = {}
    for i in (1, 2, 3, 4):
        if i == pt:
            Zd[i] = np.broadcast_to(E[None, :, :], (npts, 3, 3)).copy()
        else:
            Zd[i] = np.broadcast_to(others[i][None, None, :], (npts, 3, 3)).copy()
    Wv = W_eval(sb, combos, Xd, Zd)
    res = {}
    for k, v in Wv.items():
        v = v.reshape(3, 4, 3)          # [mu, step, nu]
        dW = np.einsum('s,msn->mn', wts, v)/h        # d_mu W_nu
        div = -dW[0, 0] + dW[1, 1] + dW[2, 2]
        res[k] = abs(div)/np.max(np.abs(dW))
    return res

if __name__ == "__main__":
    pr("# validation of the <JJJJ> partial waves (fourpt/t2_validate_l1.py)")
    for (Delta, J) in ((9.75, 3), (14.5, 4), (21.25, 6)):
        sb = SpinBlock4(1, Delta, J, N=600)
        ws, vecs = basis_coeffs(1, Delta, J)
        combos = {(a+1, b+1): (vecs[a], vecs[b]) for a in range(2) for b in range(2)}
        one = {w: (1.0 if w == 'AC' else 0.0) for w in ws}
        combos['control: word AC x v_1'] = (one, vecs[0])
        for pt in (1, 2, 3):
            d = divergence(sb, combos, pt)
            pr(f"  (Delta,J)=({Delta},{J}) conservation at x_{pt}: " + ", ".join(f"{k}: {v:.1e}" for k, v in d.items()))
        # (b) pair exchange
        rng = np.random.default_rng(5)
        base = config()
        Z = {i: rng.normal(size=(1, 2, 3)) for i in (1, 2, 3, 4)}
        X1 = {i: base[i][None, :] for i in (1, 2, 3, 4)}
        Xs = {1: X1[3], 2: X1[4], 3: X1[1], 4: X1[2]}; Zs = {1: Z[3], 2: Z[4], 3: Z[1], 4: Z[2]}
        Wa = W_eval(sb, {'12': (vecs[0], vecs[1])}, X1, Z)['12']
        Wb = W_eval(sb, {'21': (vecs[1], vecs[0])}, Xs, Zs)['21']
        pr(f"  (Delta,J)=({Delta},{J}) pair exchange W[v1,v2](1234) / W[v2,v1](3412) - 1 = {np.abs(Wa/Wb - 1).max():.1e}")
        # (c) (H,H) word pair vs H12 H34 x scalar partial wave with dimensions (3,3,3,3)
        Hc = {w: (1.0 if w == 'H' else 0.0) for w in ws}
        Wh = W_eval(sb, {'HH': (Hc, Hc)}, X1, Z)['HH']
        s, z = invariants4(X1, Z)
        r, chi, A0, E0 = block_coords_from_uv(s)
        u = s['s12']*s['s34']/(s['s13']*s['s24'])
        blk = CBlock(Delta, J, 0.0, 0.0, 600)
        # g at rho, rhobar from (r, chi) of the reflected configuration: g^{0,0}(u,v) = (-1)^J g^{0,0}(u/v, 1/v)  [e2e convention]
        S = blk.series(r, chi, sig=+1)[:, 0, 0]
        gt = (4*r)**Delta*S
        H12 = z['z1z2']*s['s12']/2 - z['z12']*z['z21']; H34 = z['z3z4']*s['s34']/2 - z['z34']*z['z43']
        K = 1/(s['s12']**3*s['s34']**3)
        Wscal = (-1)**J*(K*gt)[:, None]*H12*H34
        pr(f"  (Delta,J)=({Delta},{J}) (H,H) word pair / [H12 H34 K g_(3333)] - 1 = {np.abs(Wh/Wscal - 1).max():.1e}")
    open(os.path.join(HERE, 'out_t2_validate_l1.txt'), 'w').write("\n".join(out) + "\n")

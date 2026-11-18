"""
Referee checks (f) bulk-point density ratios of OBJECT_IDENTITY.md, and coverage of the D-safe (NF02) exact-skeleton scan (H13').
 (f): densities from CHMRS bulkpointresult (Pi-hat), bulkpointbknuagain (C: E_k), and the scalar bare commutator with weight r^{k+1};
      ratio to C expressed through the saddle r* = (1-s)/(1+s), s = sqrt(1-x^2); compared with the closed forms (1-r*^2)^k/(1+r*^2) etc.
 coverage: reads exactskel/cache/e7_<Dg>_recs.pkl (read-only) and prints the scanned (m/Lam, J, bhat) envelope and the minimum margins,
      and which states set the minimum; also the twist-gap check m >= Lam  <=>  tau >= Delta_g.
Output: cache/rf4_misc.txt
"""
import os, sys, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ES = os.path.join(HERE, '..', 'exactskel', 'cache')
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

pr("## (f) bulk-point density ratios R(k,x) relative to C_k (OBJECT_IDENTITY 2.5 table, k = 4)")
for k in (4, 2, 6):
    row1, row2 = [], []
    for x in (0.3, 0.6, 0.9, 0.99):
        s = np.sqrt(1 - x*x); r = (1 - s)/(1 + s)
        Ek = (2 - x*x)/s**(k + 2)
        Pi = 2**(2*k - 1)/(s*s*(1 + s)**(2*k - 2))
        Ps = 2**(2*k + 3)/(1 + s)**(2*k + 2)
        a1 = Pi/Ek; a2 = (1 - r*r)**k/(1 + r*r)
        b1 = Ps/Ek; b2 = (1 - r*r)**(k + 2)/(1 + r*r)
        assert abs(a1 - a2) < 1e-12*max(1, a1) and abs(b1 - b2) < 1e-12*max(1, b1)
        row1.append(f"{a1:.4g}"); row2.append(f"{b1:.4g}")
    pr(f"  k={k}: Pi-hat/C at x=0.3/0.6/0.9/0.99: {', '.join(row1)};   Psi/C: {', '.join(row2)}   (closed forms agree to 1e-12)")
# small-x expansion
for k in (4,):
    x = 1e-2; s = np.sqrt(1 - x*x); r = (1 - s)/(1 + s)
    pr(f"  small x: (1 - Pi/C)/x^4 at x=0.01, k=4: {(1-(1-r*r)**k/(1+r*r))/x**4:.5f} (predicted (k+1)/16 = {5/16:.5f})")
# normalisation of densities at x = 0
pr(f"  x=0 normalisation: E_k = 2, Pi-hat = 2^(2k-1)/2^(2k-2) = 2, Psi = 2^(2k+3)/2^(2k+2) = 2  -> all equal at x=0 (OK)")

pr("## coverage of the NF02 (D-safe) exact-skeleton scan (e7 records)")
for Dg in (100, 300, 1000):
    recs = pickle.load(open(os.path.join(ES, 'e7_%d_recs.pkl' % Dg), 'rb'))
    mf = np.array([r['mfl'] for r in recs]); J = np.array([r['J'] for r in recs]); bh = np.array([r['bh'] for r in recs])
    F = np.array([r['Fex']/r['Ns'] for r in recs])
    pr(f"  Dg={Dg}: {len(recs)} records; m/Lam in [{mf.min():.2f}, {mf.max():.2f}], J <= {J.max()}, bhat <= {bh.max():.1f}")
    for lo, hi in ((1.0, 1.1), (1.1, 1.5), (1.5, 2.0), (2.0, 3.01), (3.01, 1e9)):
        sel = (mf >= lo) & (mf < hi)
        if sel.any():
            pr(f"     m/Lam in [{lo},{hi}): {sel.sum():4d} records, J <= {J[sel].max():2d}, bhat <= {bh[sel].max():5.1f}, min F/Ns = {F[sel].min():+.4f}")
        else:
            pr(f"     m/Lam in [{lo},{hi}): NO records")
    i = F.argmin(); r = recs[i]
    pr(f"     global min F/Ns = {F[i]:+.4f} at {r['h']} J={r['J']} m/Lam={r['mfl']} Df={r['Df']}  (odd J +- state = parity-odd exchange: {r['h']=='+-' and r['J']%2==1})")
    # minimum restricted to parity-even exchanges only
    ev = np.array([not (r['h'] == '+-' and r['J'] % 2 == 1) for r in recs])
    j = np.where(ev)[0][F[ev].argmin()]
    pr(f"     min over parity-even exchanges only: {F[j]:+.4f} at {recs[j]['h']} J={recs[j]['J']} m/Lam={recs[j]['mfl']}")
    # number of states in OMEGA_NF's 'series-failing' region bhat in (20,40], m >= 1.5 Lam covered exactly
    sel = (bh > 20) & (mf >= 1.5)
    pr(f"     exact records with bhat > 20 and m >= 1.5 Lam: {sel.sum()} (OMEGA_NF series failures concentrated at bhat >= 10..40)")

pr("## twist gap: m^2 = (Delta-J-2)(Delta+J-1) >= Lam^2 = (Dg-2)(Dg-1)  <=>  tau = Delta - J >= Dg ;  a DIMENSION gap Delta >= Dg allows m < Lam:")
Dg = 100; Lam = np.sqrt((Dg - 2)*(Dg - 1))
for J in (4, 10, 20, 50, 90):
    D = max(Dg, J + 1)
    m = np.sqrt((D - J - 2)*(D + J - 1)) if D - J - 2 > 0 else 0.0
    pr(f"   Dg=100, J={J}, Delta={D}: m/Lam = {m/Lam:.4f}")
open(os.path.join(HERE, 'cache', 'rf4_misc.txt'), 'w').write("\n".join(out) + "\n")

"""S1 [computed exponent; coefficients cited]. Large-spin law of the tree-level [TT]_{n,J} anomalous dimensions.

At tree level (O(1/c_T)) and J > J_0 = 2 (Einstein contact terms have finite spin support), Lorentzian inversion
(Caron-Huot 1703.00278; spinning: Kravchuk-Simmons-Duffin) gives the [TT] data from dDisc_t of the single t-channel
T block (the t-channel double traces enter dDisc only at O(gamma^2) = O(c_T^-2)).  Collinearly every pair of external
components (T_{m1}, T_{m2}) is a pair of sl(2) primaries, and T's lowest-twist component has twist tau_T = 1.
The exact SL(2) inversion of a t-channel power ((1-zb)/zb)^eta  (Caron-Huot eq. (4.7), a=b=0) is
      I_eta(hb) = Gamma(hb-eta-1)/Gamma(hb+eta+1) * Gamma(hb)^2/Gamma(2hb-1) / Gamma(-eta)^2 .
The log z coefficient of a twist-tau' exchange has eta = tau'/2 - p_ext, the identity eta = -p_ext (p_ext = the external
collinear weight); gamma(hb) is proportional to their ratio.  This script measures the exponent of that ratio for
tau' = 1 and every collinear weight p_ext in {1/2, 1, 3/2, 2, 5/2} of the T components: it is -tau' = -1 in all cases.
So for every helicity / tensor structure  c_T gamma_{n,J} = kappa_{n,struct}/J + O(J^-2)  (kappa may vanish for a
structure, in which case the decay is faster).  Coefficients: Li-Meltzer-Poland 1511.08025 (3d, [TT]_{0,J}; negative iff
the collider bounds hold); n-growth from the AdS eikonal (Cornalba-Costa-Penedones-Schiappa hep-th/0611122, -0611123):
gamma ~ -c G h^3/hb for hb >> h >> 1 in AdS4, i.e. kappa_n ~ n^3."""
import mpmath as mp
from common_a2 import Tee
mp.mp.dps = 30
log = Tee('s1_gamma_largeJ.txt')

def I(eta, hb):
    return mp.gamma(hb - eta - 1)/mp.gamma(hb + eta + 1)*mp.gamma(hb)**2/mp.gamma(2*hb - 1)*mp.rgamma(-eta)**2

def ratio(tau, p, hb):
    return I(tau/2 - p, hb)/I(-p, hb)

log("S1: exponent d log(gamma)/d log(hb) of the single-block inversion, tau' = tau_T = 1")
log(" p_ext   hb=1e2     hb=1e3     hb=1e4     hb=1e5")
for p in (mp.mpf(1)/2, 1, mp.mpf(3)/2, 2, mp.mpf(5)/2):
    row = []
    if p == mp.mpf(1)/2:
        log(f" {float(p):4.1f}   leading coefficient = 0 (Gamma(-eta) = Gamma(0)): faster decay"); continue
    for hb in (100, 1000, 10**4, 10**5):
        hb = mp.mpf(hb)
        row.append(mp.log(ratio(1, p, 2*hb)/ratio(1, p, hb))/mp.log(2))
    log(f" {float(p):4.1f}  " + "  ".join(f"{float(x):+.6f}" for x in row))
log("=> gamma_{n,J} ~ kappa/(c_T J^{tau_T}) = kappa/(c_T J) for every component pair (tau_T = 1 in d = 3).")
log("   p_ext = 1/2 (= tau_T/2) has Gamma(-eta) = Gamma(0): that component pair's leading coefficient vanishes")
log("   and its decay is faster (only helps).")
log("Consequence used below: sup_J c_T gamma_{n,J}^2 = sup_J kappa_{n,J}^2/(c_T J^2) = O(1/c_T) -> 0 uniformly in J,")
log("since at tree level c_T gamma_{n,J} is c_T-independent and finite at every J.")

"""
v4c_wu (RESULT_Q_V4PRIME.md §3): the u-channel-identity Wick term W_u = <T4 T2><T3 T1> at the edge-corner of the swapped k = 4 form.

For light-ray polarisations the numerator N = (z1.z3) x13^2 - 2 (z1.x13)(z3.x13) vanishes IDENTICALLY on the cones (any t, phi), so W_u = 0
for phi1 != phi3.  At the collinear point phi1 = phi3 (the corner u' -> 0, r -> 0, rho -> 1) x13^2 = 0 as well: the restriction of the
Wightman 2-pt function to the product of the two light cones is not defined there without a prescription (its wave-front covector
(k, -k), k ~ z, is conormal to the cone product).  We compare
  (A) a timelike-tube regularisation  x1 -> x1 + i eps e  (Im x13 = eps e in the forward/backward tube; e timelike),
  (B) a cone-preserving one  t -> t + i eps (Im x13 along the cone generators; degenerates to a null direction at phi1 = phi3),
and measure  G(eps) = int d(dphi) A13^{(2+delta)/2} K_nu(dphi) U_eps(dphi),  U = (-2 pi i)^2 Res_{t1,t3=1/2}[t1^3 t3^3 f W_u],
for f = f4 (k = 4 edge residue) and, as a control, the same restriction without f's 1/x13^2 pole.
K_nu(dphi) = int dphibar [A10 A30]^{-delta/2} (rotation-invariant part of the celestial structure, z0 fixed).
Analytic expectation: U_eps ~ eps^4/(x13^2 + i eps)^6 with x13^2 ~ -dphi^2/4  =>  G ~ eps^{-1/4} (k = 4, f's pole) vs eps^{+3/4} (control).
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v4c_core as V

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
lines = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)


T42 = None


def t42():
    """<T4^{mu nu} T2_{mu nu}> for x24 = -2e, C_T = 1: [1/2 (I_mr I_ns + I_ms I_nr) - 1/3 eta_mn eta_rs] eta^mr eta^ns / (x24^2)^3."""
    X = -2*V.E0; X2 = V.dot(X, X)
    Il = V.ETA - 2*np.outer(V.lower(X), V.lower(X))/X2      # I_{mu nu}
    Iu = np.linalg.inv(V.ETA)
    A = np.einsum('mr,ns,mr,ns->', Il, Il, Iu, Iu); B = np.einsum('ms,nr,mr,ns->', Il, Il, Iu, Iu)
    Cc = np.einsum('mn,rs,mr,ns->', V.ETA, V.ETA, Iu, Iu)
    return (0.5*(A + B) - Cc/3)/X2**3


def U_eps(dphi, eps, scheme='A', withf=True, r=0.2, M=24, fk='f4', pw=5.0):
    """edge residue of t1^3 t3^3 f T42 N^2/(X13^2)^pw with the (13) pair regularised; phi1 = dphi/2, phi3 = -dphi/2.
    fk = 'f4' or 'f5' (the subtraction factor whose x12^2, x34^2 poles give the edge); pw = 5 is the T 2-pt function;
    pw = 5 - tau/2 models the u-channel lightcone suppression (x13^2)^{tau/2} of a twist-tau u-channel exchange (O(G): tau_T = 1)."""
    dphi = np.atleast_1d(dphi)
    z1 = V.zvec(dphi/2); z3 = V.zvec(-dphi/2)
    th = 2*np.pi*np.arange(M)/M; w = r*np.exp(1j*th)
    T1 = (0.5 + w)[:, None, None]; T3 = (0.5 + w)[None, :, None]
    if scheme == 'A':
        x1 = -T1[..., None]*z1[None, None] + 1j*eps*V.E0
        x3 = T3[..., None]*z3[None, None]
    else:                                     # cone-preserving: complex t along the generators
        x1 = -(T1 + 1j*eps)[..., None]*z1[None, None]
        x3 = (T3 + 1j*eps)[..., None]*z3[None, None]
    x13 = x1 - x3
    X13 = V.dot(x13, x13)
    N = V.dot(z3[None, None], z1[None, None])*X13 - 2*V.dot(z3[None, None], x13)*V.dot(z1[None, None], x13)
    # f4 with the regularised x13^2; x12^2, x34^2 etc. from the UNshifted ray points (the edge residue is in t)
    x1r = -T1[..., None]*z1[None, None]; x3r = T3[..., None]*z3[None, None]
    s12 = V.dot(x1r + V.E0, x1r + V.E0); s34 = V.dot(x3r - V.E0, x3r - V.E0)
    s14 = V.dot(x1r - V.E0, x1r - V.E0); s23 = V.dot(-V.E0 - x3r, -V.E0 - x3r); s24 = -4.0
    up = s12*s34/(X13*s24); vp = s14*s23/(X13*s24)
    f = ((vp - up)**3 - (vp - up))/(4*up*vp) if fk == 'f4' else -((up - vp)**2 - 1)/(4*up*vp)
    g = (T1**3*T3**3)*T42*N**2/X13**pw
    if withf:
        g = g*f*(2*T1 - 1)*(2*T3 - 1)
        co = np.fft.fft2(g, axes=(0, 1))/(M*M)
        return (-2j*np.pi)**2*co[0, 0]/4.0
    # control: restriction at t = 1/2 (no residue), no f
    co = np.fft.fft2(g, axes=(0, 1))/(M*M)
    return co[0, 0]


def U_comm(dphi, eps, t=0.3, fk='f4'):
    """commutator-form analogue at an INTERIOR point of the large-eta corner (t1 = t3 = t, rho = 2t < 1, no edge residue):
    Delta_13 [t1^3 t3^3 f T42 N^2/(X13^2)^5] = g(+eps) - g(-eps), x1 -> x1 +- i eps e (the (13) i-eps flips in [T3, T1]; f's 1/x13^2
    pole carries the same i-eps, Lemma C)."""
    dphi = np.atleast_1d(dphi)
    z1 = V.zvec(dphi/2); z3 = V.zvec(-dphi/2)
    out = 0
    for sgn in (+1, -1):
        x1 = -t*z1 + 1j*sgn*eps*V.E0; x3 = t*z3; x13 = x1 - x3
        X13 = V.dot(x13, x13)
        N = V.dot(z3, z1)*X13 - 2*V.dot(z3, x13)*V.dot(z1, x13)
        x1r = -t*z1; x3r = t*z3
        s12 = V.dot(x1r + V.E0, x1r + V.E0); s34 = V.dot(x3r - V.E0, x3r - V.E0)
        s14 = V.dot(x1r - V.E0, x1r - V.E0); s23 = V.dot(-V.E0 - x3r, -V.E0 - x3r)
        up = s12*s34/(X13*(-4.0)); vp = s14*s23/(X13*(-4.0))
        f = ((vp - up)**3 - (vp - up))/(4*up*vp) if fk == 'f4' else (vp - up if fk == 'f2' else 1.0)
        out = out + sgn*t**6*f*T42*N**2/X13**5
    return out


def Kfun(dphi, nu):
    """K_nu(dphi) = int_0^{2pi} dpb [A10 A30]^{-delta/2}, A_i0 = 2(1 - cos(pb -+ dphi/2)) (phi0 = 0 by rotation invariance).
    Singular (|x|^{-1/2 - i nu}) at pb = +-dphi/2: mpmath tanh-sinh on the pieces between the singular points."""
    import mpmath as mp
    d = mp.mpf(0.5) + 1j*mp.mpf(nu)
    h = mp.mpf(abs(dphi))/2
    f = lambda pb: (16*mp.sin((pb + h)/2)**2*mp.sin((pb - h)/2)**2)**(-d/2)
    return complex(mp.quad(f, [-h, h]) + mp.quad(f, [h, 2*mp.pi - h]))


_KC = {}


def Kcache(x, nu):
    key = (round(float(np.log(x)), 6), nu)
    if key not in _KC:
        _KC[key] = Kfun(x, nu)
    return _KC[key]


def G(eps, nu, scheme='A', withf=True, ns=1601, smax=60.0, fk='f4', pw=5.0, dshift=0.0):
    """int over dphi in (-pi, pi) of A13^{(2+delta)/2} K U; dphi = sqrt(eps) s on the core, plus the outer region.
    dshift: delta -> delta + dshift in the A13 power only (the contour shift nu -> nu - i dshift that the de-shift weight performs on
    the Q^{i nu/2} factor, Q = A13/(A10 A30); K is kept at real nu, i.e. only the A13 power that controls the corner is shifted)."""
    d = 0.5 + 1j*nu + dshift
    se = np.sqrt(eps)
    # core: |dphi| <= smax sqrt(eps) (log grid), outer up to pi
    s = np.concatenate([-np.logspace(np.log10(smax), -4, ns)[::-1][::-1], np.logspace(-4, np.log10(smax), ns)])
    s = np.sort(s)
    dph = se*s
    outer = np.concatenate([-np.logspace(np.log10(np.pi*0.999), np.log10(smax*se), 400), np.logspace(np.log10(smax*se), np.log10(np.pi*0.999), 400)])
    dph = np.unique(np.concatenate([dph, outer]))
    dph = dph[np.abs(dph) < np.pi]
    kg = np.logspace(np.log10(np.abs(dph).min()), np.log10(np.pi), 40)
    Kv = np.array([Kcache(x, nu) for x in kg])
    lk = np.log(kg)
    Ki = np.interp(np.log(np.abs(dph)), lk, Kv.real) + 1j*np.interp(np.log(np.abs(dph)), lk, Kv.imag)
    A13 = 2*(1 - np.cos(dph))
    if scheme == 'comm':
        Uv = U_comm(dph, eps, fk=fk)
    else:
        Uv = np.concatenate([U_eps(dph[i:i + 400], eps, scheme, withf, fk=fk, pw=pw) for i in range(0, len(dph), 400)])
    integ = A13**((2 + d)/2)*Ki*Uv
    return np.trapezoid(integ, dph), np.max(np.abs(Uv))


if __name__ == '__main__':
    T42 = t42()
    pr('# v4c_wu: u-channel identity (light-ray polarised) at the edge-corner, k = 4;  T42 =', T42)
    # pointwise: W_u numerator vanishes identically off the collinear point (unregularised)
    for dph in (0.3, 0.03, 0.003):
        z1 = V.zvec(np.array([dph/2])); z3 = V.zvec(np.array([-dph/2]))
        x13 = -0.5*z1 - 0.5*z3
        N = V.dot(z3, z1)*V.dot(x13, x13) - 2*V.dot(z3, x13)*V.dot(z1, x13)
        pr(f'  unregularised N(dphi={dph}) = {N[0]:.2e}  (identically 0: W_u = 0 off phi1 = phi3)')
    pr('  (B) cone-preserving prescription (t -> t + i eps): x13 = -a z1 - b z3 with complex a, b gives X13^2 = 2ab z1.z3, z1.x13 = -b z1.z3,'
       ' z3.x13 = -a z1.z3, hence N = 2ab (z1.z3)^2 - 2ab (z1.z3)^2 = 0 IDENTICALLY (exact algebra): G_B = 0 for every eps.'
       ' But Im x13 = eps (z1 + z3) becomes NULL at phi1 = phi3, so (B) is not an interior tube direction at the corner.')
    EPS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)

    def scan(label, expect, nu, **kw):
        pr(f'  nu = {nu}, {label}:')
        rows = []
        sch = kw.pop('scheme', 'A')
        for eps in EPS:
            gA, uA = G(eps, nu, sch, **kw)
            rows.append((eps, gA))
            pr(f'    eps={eps:.0e}:  (A) timelike tube G = {gA.real:+.4e}{gA.imag:+.4e}i  |G| = {abs(gA):.4e}  max|U_A| = {uA:.2e}')
        e = np.array([r[0] for r in rows]); a = np.array([abs(r[1]) for r in rows])
        sl = np.polyfit(np.log(e[-4:]), np.log(a[-4:]), 1)[0]
        pr(f'    fitted |G_A| ~ eps^{sl:.3f}  (power counting: {expect:+.2f}; log corrections at nu = 0)')
        return sl
    for nu in (0.0, 1.0):
        scan('k=4 edge residue (f4 has a 1/x13^2 pole at the corner)', -0.25, nu)
        scan('control: restriction at t = 1/2 without f', 0.75, nu, withf=False)
    pr('  k = 5 (f5 = -((u\'-v\')^2-1)/(4u\'v\'): its edge residue -(v\'^2-1)/(4v\') x x13^2 has NO 1/x13^2 pole):')
    scan('k=5 edge residue, j\'=0-type kernel magnitude', 0.75, 0.0, fk='f5')
    # the j'=1 celestial structure is odd under phi1 <-> phi3 while U_eps (symmetric regularisation) is even:
    xs = np.array([1e-3, 1e-2, 0.1])
    ev = U_eps(xs, 1e-4, fk='f5'); od = U_eps(-xs, 1e-4, fk='f5')
    pr(f'    U_eps(f5) even in dphi: max |U(dphi) - U(-dphi)|/|U| = {np.max(np.abs(ev - od)/np.abs(ev)):.1e}  -> its pairing with the'
       ' j\'=1 (1<->3-odd) celestial structure vanishes identically, in any 1<->3-symmetric regularisation')
    pr('  de-shifted C^4 (Q\'s functional): the 4-shift weight moves the nu-contour by up to 3i on the factor Q^{i nu/2}, Q = A13/(A10 A30)'
       ' ~ dphi^2 at the corner (the same shift that turns the raw eta^{5/2} kernel growth into eta^{-1/2} log eta); test with delta -> delta + 3:')
    scan('k=4 edge residue, contour-shifted kernel (delta + 3)', 1.25, 0.0, dshift=3.0)
    pr('  commutator (unswapped) form at an interior large-eta-corner point t1 = t3 = 0.3 (rho = 0.6), Delta_13 of f W_u:')
    scan('commutator form, f4 (raw k = 4)', -0.25, 0.0, scheme='comm', fk='f4')
    scan('commutator form, f2 (raw k = 2: f2 = v\'-u\' also ~ 1/x13^2 at fixed rho)', -0.25, 0.0, scheme='comm', fk='f2')
    scan('commutator form, f = 1 (unsubtracted control)', 0.75, 0.0, scheme='comm', fk='one')
    pr('  O(G) power counting: a twist-tau u-channel exchange carries an extra (x13^2)^{tau/2} relative to the u-channel identity in the'
       ' (13) lightcone limit; tau_T = 1 (graviton/T exchange), same light-ray numerator structure (natural-size model, [D-sketch]):')
    scan('k=4 edge residue, u-channel twist-1 (T-exchange) natural size, raw kernel', 0.25, 0.0, pw=4.5)
    scan('k=4 edge residue, u-channel twist-1, de-shifted kernel (delta + 3)', 1.75, 0.0, pw=4.5, dshift=3.0)
    open(os.path.join(CACHE, 'v4c_wu.txt'), 'w').write('\n'.join(lines) + '\n')

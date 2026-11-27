"""Item 1/2: Hall response for theta(r)-weighted P and Q with GENERIC theta (horizon locality test).  [computed]"""
import io, contextlib
import sympy as sp
with contextlib.redirect_stdout(io.StringIO()):
    import ef_engine as E
r = E.r
th = E.th          # theta(r) Function defined in the engine


def ibp_reduce(sig):
    """Write sig = c(r) theta + d/dr B[theta]; return (c, B)."""
    sig = sp.expand(sig)
    B = 0
    for k in range(6, 0, -1):
        dk = sp.Derivative(th, (r, k))
        sk = sp.expand(sig).coeff(dk)
        if sk == 0:
            continue
        # sk theta^(k) = d(sk theta^(k-1)) - sk' theta^(k-1)
        lower = sp.Derivative(th, (r, k - 1)) if k > 1 else th
        B += sk * lower
        sig = sp.expand(sig - sk * dk - sp.diff(sk, r) * lower + 0) + 0
        sig = sp.expand(sig + 0)
        sig = sp.expand(sig) 
        sig = sig  # (sig now: original minus d/dr(sk theta^(k-1)) )
    return sp.simplify(sig), sp.simplify(B)


for name, Op in (('P', E.P_), ('Q', E.Q_)):
    with contextlib.redirect_stdout(io.StringIO()):
        sig = E.hall(E.cross_part(sp.expand(r ** 2 * th * Op[2])), name, integrate=False)
    c, B = ibp_reduce(sig)
    # check: sig == c theta + dB/dr
    chk = sp.simplify(sp.expand(sig - c * th - sp.diff(B, r)))
    T = sp.symbols('t0:4')
    Bh = sp.simplify(B.subs({sp.Derivative(th, (r, 2)): T[2], sp.Derivative(th, r): T[1]}).subs(th, T[0]).subs(r, 1))
    print(f"[theta {name}] sigma = {sp.simplify(sig)}")
    print(f"[theta {name}]   = c(r) theta + B'  with c = {c},  B = {B};  check residual = {chk}")
    print(f"[theta {name}]   => 4 pi eta_H/s = int_1^oo sigma = -B(1) (theta, theta' -> 0 fast enough at oo) = {sp.expand(-Bh)}"
          f"   (t0 = theta(1), t1 = theta'(1), t2 = theta''(1))")

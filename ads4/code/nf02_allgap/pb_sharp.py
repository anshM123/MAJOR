"""numerical (Arb balls, exact inputs) scan of d -> L[(1-z)^d] at large-n twist-edge states: where does the positive
basis stop being positive?  (diagnostic for the sharpness of Lemma PB+; numerical evidence)"""
from flint import arb, ctx
ctx.prec = 6000
def Lval(J, b, a, d):
    T = arb((a+J+1)*(a+J)); C = J*(J+1); n = J-b
    s = arb(0); U = arb(1); bin_ = arb(1)
    for i in range(min(d, n)+1):
        s += (-1)**i*bin_*U
        uj = (b+i)*(b+i+1); U = U*(C-uj)/(T-uj); bin_ = bin_*(d-i)/(i+1)
    return s
for (J, b, a) in [(2000, 0, 998), (20000, 0, 998), (2000, 4, 998)]:
    out = []
    for d in [999, 1100, 1500, 2000, 3000, 5000, 10000, 20000]:
        v = Lval(J, b, a, d)
        out.append((d, v.str(3, radius=False)))
    print(J, b, a, out, flush=True)

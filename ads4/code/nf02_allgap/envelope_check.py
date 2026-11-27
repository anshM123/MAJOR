"""numerical check (Arb balls at sample points; NOT a proof) of the RESEARCH_NOTE Sec.4 envelope
(k!)^2 t^k |g[k,3]((nu+-i)^2)| <= (1+y^2)^k,  y=(nu+2)/L, t=1/L^2, for k up to 201 (note proves it for k<=64)."""
from flint import arb, acb, ctx
ctx.prec = 3000
L2 = arb(997002); L = L2.sqrt(); t = 1/L2
worst = (-1e9, None)
for nu in [0, 1, 3, 10, 30, 100, 300, 600, 900, 998, 1200, 1500, 2000, 3000, 5000, 1e4, 1e5, 1e6]:
    for sgn in (1, -1):
        z = acb(nu, sgn)**2
        y = (arb(nu)+2)/L
        gm, g0 = acb(0), acb(1)
        fact2 = arb(1); tk = arb(1)
        for k in range(0, 202):
            if k >= 1:
                val = (fact2*tk*abs(g0)).log() - k*(1+y*y).log()
                v = float(val.mid())
                if v > worst[0]: worst = (v, (nu, sgn, k))
            g1 = ((8*k*k+24*k+11 - 4*z)*g0 - 4*(k+2)**2*gm)/(4*(k+1)**2)
            gm, g0 = g0, g1
            fact2 *= (k+1)**2; tk *= t
print('max over samples of log[(k!)^2 t^k |g_k| / (1+y^2)^k], k=1..201:', worst)

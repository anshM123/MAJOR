"""pernu.py -- numerical study of the gap-independent per-frequency kernel of Lemma DIL,
K_st(nu) = (1-r)^3 sum_k U_k [g_k + (gam_k g_k + (k+1)^2 g_{k+1}/2)/T],  g_k = g[k,3](nu^2),
against its leading-Regge limit P_{-1/2+i nu}(eta), eta=(1+r)/(1-r).  Numerical evidence only.
Prints max over nu-grid of |K|, |P|, |K-P| and the 'relative envelope' |K-P|/env, env=max(|P|) over a local window."""
import sys, math
from flint import arb, acb, ctx
def kernel(J, b, a, nus, prec=4000):
    ctx.prec = prec
    T = arb((a+J+1)*(a+J)); C = J*(J+1); n = J-b; r = arb(C)/T
    out = []
    for nu in nus:
        y = arb(nu)**2
        g = [arb(1)]; gm = arb(0)
        for k in range(n+1):
            g1 = ((8*k*k+24*k+11 - 4*y)*g[-1] - 4*(k+2)**2*gm)/(4*(k+1)**2)
            gm = g[-1]; g.append(g1)
        S = arb(0); U = arb(1)
        for k in range(n+1):
            gam = arb(-k*k+(8-4*b)*k+9-2*b*(b-1))/2
            S += U*(g[k] + (gam*g[k] + arb((k+1)**2)/2*g[k+1])/T)
            uj = (b+k)*(b+k+1); U = U*(C-uj)/(T-uj)
        K = (1-r)**3*S
        eta = (1+r)/(1-r)
        P = acb((1-eta)/2).hypgeom_2f1(acb(0.5, -nu), acb(0.5, nu), acb(1)).real
        out.append((nu, K, P))
    return out
if __name__ == '__main__':
    for (J, b, a) in [(420, 4, 998), (420, 0, 998), (2000, 0, 998), (1000, 4, 1500)]:
        nus = [i*5 for i in range(0, 201)]
        res = kernel(J, b, a, nus)
        T = (a+J+1)*(a+J); r = J*(J+1)/T; beta = math.acosh((1+r)/(1-r))
        mK = max(abs(float(K.mid())) for _, K, P in res); mP = max(abs(float(P.mid())) for _, K, P in res)
        mD = max(abs(float((K-P).mid())) for _, K, P in res)
        worst = max(res, key=lambda t: abs(float((t[1]-t[2]).mid())))
        print('J=%d b=%d a=%d beta=%.3f  max|K|=%.3e max|P|=%.3e max|K-P|=%.3e at nu=%g; |K-P|/|P| at nu=0: %.3e; maxrad %.1e' % (
            J, b, a, beta, mK, mP, mD, worst[0], abs(float(((res[0][1]-res[0][2])/res[0][2]).mid())), max(float(K.rad()) for _, K, P in res)), flush=True)
        for (nu, K, P) in res[::20]:
            print('   nu=%5g K=% .4e P=% .4e K-P=% .3e' % (nu, float(K.mid()), float(P.mid()), float((K-P).mid())))

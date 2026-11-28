# Independent check of the endpoint integral in Prop. S4 (renormalised Crofton):
#   m{u>u0, v>v0, uv<Y^2}, dm = du dv/(u+v)^2  ->  log(Y/d) + o(1),  d=u0+v0
# and Remark S4(d): m{u>u0,v>v0,u+v<2Y} = log(2Y/d) - 1 + o(1).
import numpy as np
from scipy import integrate
def m_cross(u0,v0,Y):
    # integrate over u from u0 to Y^2/v0 of [ int_{v0}^{Y^2/u} dv/(u+v)^2 ]
    f=lambda u: 1/(u+v0) - 1/(u+Y**2/u)
    val,err=integrate.quad(f,u0,Y**2/v0,limit=500,points=[Y])
    return val
def m_disj(u0,v0,Y):
    d=u0+v0
    f=lambda s:(s-d)/s**2
    return integrate.quad(f,d,2*Y)[0]
for (u0,v0) in [(0.3,0.7),(0.05,1.95),(1.0,1.0)]:
    d=u0+v0
    for Y in [1e2,1e3,1e4]:
        a=m_cross(u0,v0,Y); b=m_disj(u0,v0,Y)
        print(f"u0={u0} v0={v0} Y={Y:.0e}: cross-mass - log(Y/d) = {a-np.log(Y/d):+.2e};  disjoint-mass - (log(2Y/d)-1) = {b-(np.log(2*Y/d)-1):+.2e}")

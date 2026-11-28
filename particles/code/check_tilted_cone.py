#!/usr/bin/env python3
"""Check the static positive-cone failure of the proposed chiral formula.

This is a counterexample to a length identification, not an entropy-cone
violation. The geometry and acausal homology proof are in
TILTED_CONE_COUNTEREXAMPLE.md. Regulator: subtract 2 log(rho_cut).
"""
from pathlib import Path
import json
import math
import mpmath as mp

mp.mp.dps = 80
theta=3*mp.pi/2
dphi=11*mp.pi/16
dt=9*mp.pi/16
N=2*(mp.cos(dt)-mp.cos(dphi))
L0=mp.log(N)
Lm=mp.log(2*(mp.cos(dt)-mp.cos(dphi-theta)))
F0=mp.log(4*mp.sin(mp.pi/16))
Fm=mp.log(4*mp.sin(mp.pi/8))
gap=F0-L0
assert dt < min(dphi,theta-dphi)
assert dphi < min(mp.pi,theta/2)
assert N>0 and L0<Lm and F0<Fm
assert abs(gap+mp.log(mp.sin(5*mp.pi/8)))<mp.mpf('1e-75')
assert abs(L0-mp.log(4*mp.sin(5*mp.pi/8)*mp.sin(mp.pi/16)))<mp.mpf('1e-75')
p0=mp.matrix([0,1,0,1])
p1=mp.matrix([mp.sin(dt),mp.cos(dphi),mp.sin(dphi),mp.cos(dt)])
def dot(a,b):
    return -a[0]*b[0]+a[1]*b[1]+a[2]*b[2]-a[3]*b[3]
res=mp.mpf(0)
for s in map(mp.mpf,[-10,-4,-1,0,1,4,10]):
    X=(mp.exp(s)*p0+mp.exp(-s)*p1)/mp.sqrt(N)
    V=(mp.exp(s)*p0-mp.exp(-s)*p1)/mp.sqrt(N)
    res=max(res,abs(dot(X,X)+1),abs(dot(V,V)-1),abs(dot(X,V)))
assert res<mp.mpf('1e-65')
report={
 'status':'PASS; mismatch of formulas, not an entropy-cone violation',
 'decimal_precision':mp.mp.dps,
 'theta':'3*pi/2','delta_phi':'11*pi/16','delta_t':'9*pi/16',
 'chiral_sweeps':['5*pi/4','pi/8'],
 'boundary_arc_slopes':['9/11','-9/13'],
 'physical_direct_length':str(L0),'physical_other_length':str(Lm),
 'relaxed_direct_half_sum':str(F0),'relaxed_other_half_sum':str(Fm),
 'relaxed_minimum_minus_physical_minimum':str(gap),
 'different_branch_mismatches':[str(F0-L0),str(Fm-Lm)],
 'minimum_circumferential_radius_squared':str((1+mp.cos(dphi))/(mp.cos(dt)-mp.cos(dphi))),
 'centered_plane_slope_bound':str(mp.sin(dt/2)/mp.sin(dphi/2)),
 'ambient_unit_geodesic_check_max_residual':str(res),
 'cutoff':'Common subtraction 2 log(rho_cut); differences are finite.'
}
ROOT=Path(__file__).resolve().parents[1]
(ROOT/'data').mkdir(exist_ok=True)
(ROOT/'figures').mkdir(exist_ok=True)
(ROOT/'data/tilted_cone_validation.json').write_text(json.dumps(report,indent=2)+'\n')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
xx=np.linspace(0,2.,500)
phi=float(dphi)
ang=float(theta)
def hyp(delta):
    delta=abs(delta)
    return 2*np.log(2*np.sin(delta/2)) if delta<math.pi else 2*np.log(2.)
lor=np.log(2*(np.cos(xx)-np.cos(phi)))
rel=np.array([min((hyp(phi+t)+hyp(phi-t))/2,
                  (hyp(phi+t-ang)+hyp(phi-t-ang))/2) for t in xx])
fig,ax=plt.subplots(1,2,figsize=(9.8,3.9),layout='constrained')
ax[0].plot(xx,lor,lw=2.3,color='#16697a',label='Regular spacetime HRT')
ax[0].plot(xx,rel,lw=2.1,color='#ae5d32',ls='--',label='Relaxed chiral-cone half-sum')
ax[0].set(xlabel='Endpoint time separation $\\Delta t$',ylabel='Renormalized length',title='The proposed formula separates from HRT')
ax[0].legend(frameon=False,fontsize=8.7)
ax[1].plot(xx,rel-lor,lw=2.3,color='#805a92')
ax[1].axvline(float(mp.pi-dphi),ls=':',lw=1,color='.5')
ax[1].scatter([float(dt)],[float(gap)],color='#805a92',s=35,zorder=4)
ax[1].annotate('$\\Delta t=9\\pi/16$\nGap $=0.07917359\\ldots$',xy=(float(dt),float(gap)),
               xytext=(.13,.75),textcoords='axes fraction',fontsize=9,
               arrowprops={'arrowstyle':'->','color':'.4'})
ax[1].set(xlabel='Endpoint time separation $\\Delta t$',ylabel='Half-sum minus HRT',title='Failure starts when a chiral sweep exceeds $\\pi$')
fig.suptitle('Static spinless cone: $\\theta=3\\pi/2$, $\\Delta\\phi=11\\pi/16$ — a length-identity counterexample',fontsize=11)
fig.savefig(ROOT/'figures/tilted_cone_counterexample.png',dpi=190)
fig.savefig(ROOT/'figures/tilted_cone_counterexample.pdf')
print(json.dumps(report,indent=2))

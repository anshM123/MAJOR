"""Verify the explicit all-angle filled-topology Crofton obstruction.

The analytic proof is in CONE_INTERSECTION_AUDIT.md. This script evaluates
the two intersections at 80 decimal digits and makes a topological drawing.
It is a check of closed formulas, not an interval-arithmetic certificate.
"""
from pathlib import Path
import json
import mpmath as mp

mp.mp.dps = 80
ROOT = Path(__file__).resolve().parent


def bisect(f, a, b):
    fa, fb = f(a), f(b)
    assert fa * fb < 0
    for _ in range(280):
        mid = (a + b) / 2
        fm = f(mid)
        if fa * fm <= 0:
            b, fb = mid, fm
        else:
            a, fa = mid, fm
    return (a + b) / 2


def example(theta, beta, eta):
    delta = theta - beta + 2 * eta
    assert 0 < theta < 2 * mp.pi
    assert max(0, theta-mp.pi) < beta < theta/2
    assert 0 < eta < beta/2
    assert delta < min(mp.pi, theta)
    kc, kg = mp.cos(beta/2), mp.cos(delta/2)
    muc, mug = beta/2, (theta+beta)/2
    f = lambda p: kc*mp.cos(p-mug) - kg*mp.cos(p-muc)
    p1 = bisect(f, beta-eta, beta)
    p2 = theta+beta-p1
    r1 = mp.atanh(kc/mp.cos(p1-muc))
    r2 = mp.atanh(kc/mp.cos(p2-(theta+muc)))
    rg1 = mp.atanh(kg/mp.cos(p1-mug))
    rg2 = mp.atanh(kg/mp.cos(p2-mug))
    assert beta-eta < p1 < beta
    assert theta < p2 < theta+eta
    assert abs(r1-rg1) < mp.mpf('1e-70')
    assert abs(r2-rg2) < mp.mpf('1e-70')
    assert abs(r1-r2) < mp.mpf('1e-70')
    return dict(theta=str(theta), beta=str(beta), eta=str(eta),
                theta_over_pi=str(theta/mp.pi), beta_over_pi=str(beta/mp.pi),
                eta_over_pi=str(eta/mp.pi), test_sweep_over_pi=str(delta/mp.pi),
                intersection_angles=[str(p1),str(p2-theta)],
                intersection_radius=str(r1), actual_intersections=2,
                filled_intersection_number=0, punctured_intersection_number=2)


def main():
    cases=[example(mp.pi/2,mp.pi/6,mp.pi/48),
           example(3*mp.pi/2,5*mp.pi/8,mp.pi/32)]
    for n in range(1,40):
        theta=mp.pi*n/20
        lower=max(mp.mpf(0),theta-mp.pi)
        beta=(lower+theta/2)/2
        eta=min(beta/8,(mp.pi-theta+beta)/8,(2*mp.pi-theta)/8)
        example(theta,beta,eta)
    result={'analytic_statement':'All angles 0<theta<2pi; see proof note',
            'arithmetic':'mpmath, 80 decimal digits; formula check, not interval proof',
            'angle_grid_checked':39,'explicit_cases':cases,'status':'PASS'}
    (ROOT/'cone_bigon_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(cases)
    print(json.dumps(result,indent=2))


def plot(cases):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(9,4.7),layout='constrained')
    def curve(a,b,theta):
        pp=np.linspace(a+1e-6,b-1e-6,1500)
        rr=np.arctanh(np.cos((b-a)/2)/np.cos(pp-(a+b)/2))
        rho=np.tanh(rr/2)
        return rho*np.cos(2*np.pi*pp/theta),rho*np.sin(2*np.pi*pp/theta)
    for ax,d,title in zip(axs,cases,[r'$\theta=\pi/2<\pi$',r'$\theta=3\pi/2>\pi$']):
        theta,beta,eta=[float(d[k]) for k in ['theta','beta','eta']]
        cc=np.linspace(0,2*np.pi,300)
        ax.plot(np.cos(cc),np.sin(cc),color='.5',lw=1)
        ax.plot(*curve(0,beta,theta),lw=2.6,color='#2166ac',label='Minimizing cut $c$')
        ax.plot(*curve(beta-eta,theta+eta,theta),lw=2.2,color='#d6604d',label='Regular test geodesic $g$')
        rr=float(d['intersection_radius']);rho=np.tanh(rr/2)
        for phi in map(float,d['intersection_angles']):
            ax.scatter([rho*np.cos(2*np.pi*phi/theta)],[rho*np.sin(2*np.pi*phi/theta)],
                       s=42,facecolor='white',edgecolor='black',zorder=6)
        ax.scatter([0],[0],marker='*',s=80,color='black',zorder=6)
        ax.text(0,-.12,'cone',ha='center',fontsize=9)
        ax.set_aspect('equal');ax.set_xlim(-1.08,1.08);ax.set_ylim(-1.08,1.08)
        ax.set_axis_off();ax.set_title(title,fontsize=13)
        ax.text(0,-1.22,r'$\#(g\cap c)=2,\quad i_{\rm filled}(g,c)=0$',ha='center',fontsize=11)
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=2,frameon=False)
    fig.suptitle('Cone curvature obstructs filled-class Crofton counting\nTopological disk rendering; angles and distances are distorted',fontsize=12)
    fig.savefig(ROOT/'cone_bigon.png',dpi=200)
    fig.savefig(ROOT/'cone_bigon.pdf')
    plt.close(fig)


if __name__=='__main__':
    main()

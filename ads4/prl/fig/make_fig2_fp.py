#!/usr/bin/env python3
"""Fig. 2 (v17): the certificate cover of the heavy spectrum for Omega^EE at Delta_g = 1000 (Theorem FP).
Plane: impact parameter L = -ln(1-c) (ladder variable) vs twist a = Delta-J-2. Shown for the b = 0 ladder (s = 1).
a*(beta) from code/fullpos3 high-twist records."""
import json, glob, math, os
from fractions import Fraction as F
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__))
FP3=os.path.join(HERE,'..','..','code','fullpos3')
rows=[]
for f in glob.glob(os.path.join(FP3,'cache','laptop_*.jsonl'))+[os.path.join(FP3,'results','ht_scan_done.jsonl'),os.path.join(FP3,'results','ht_large_done.jsonl')]:
    if os.path.exists(f):
        for l in open(f):
            if l.strip(): rows.append(json.loads(l))
cells=sorted({(float(F(r['ba'])),float(F(r['bb']))):max(r['a0_b0'],r['a0_b4']) for r in rows}.items())
def L_of_beta(b): return 2*math.log(math.cosh(b/2))
amin,amax=998,2.0e5
Lmin,Lmax=0.02,80
plt.rcParams.update({'font.size':7.5,'axes.linewidth':0.6})
fig,ax=plt.subplots(figsize=(3.4,2.55))
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(Lmin,Lmax); ax.set_ylim(amin,amax)
bands=[(0.0513,0.2984,'#deebf7','direct prefix'),(0.2984,2.3026,'#c6dbef','coupled-cap norming'),
       (2.3,3.65,'#9ecae1','signed transfer'),(3.65,5,'#6baed6',r'route T$^{\prime\prime}$'),(5,20,'#4292c6',r'route T$^{\prime}$'),(20,Lmax,'#2171b5','tail theorem')]
for a,b,col,lab in bands:
    ax.axvspan(a,b,color=col,lw=0)
    ax.axvline(a,color='w',lw=0.5)
    xm=math.sqrt(a*b)
    ax.text(xm,1.08e3,lab,ha='center',va='bottom',fontsize=5.6,rotation=90,color='k' if col in ('#deebf7','#c6dbef','#9ecae1') else 'w')
# high twist region a >= a*(beta): hatched overlay plus its boundary
xs=[];ys=[]
for (ba,bb),ast in cells:
    La,Lb=max(L_of_beta(ba),Lmin),min(L_of_beta(bb),Lmax)
    if Lb<=La: continue
    ax.fill_between([La,Lb],[ast,ast],[amax,amax],facecolor='#fd8d3c',alpha=0.35,lw=0)
    xs+= [La,Lb]; ys+=[ast,ast]
ax.plot(xs,ys,color='#a63603',lw=0.8)
ax.text(1.0,4.0e4,'high twist, $a\\geq a^*(\\beta)$',ha='center',fontsize=6.5,color='#7f2704')
# J <= 1000 region (b=0): a <= 1000 (1-sqrt c)/sqrt c
L=np.logspace(math.log10(Lmin),math.log10(Lmax),400); c=1-np.exp(-L); aJ=1000*(1-np.sqrt(c))/np.sqrt(c)
ax.fill_between(L,amin,np.clip(aJ,amin,amax),where=aJ>amin,color='#969696',alpha=0.85,lw=0)
ax.text(0.024,1.12e3,'$J\\leq1000$',fontsize=6.0,color='w')
# b-hat <= 400 strip (ell <= 0.039)
ax.axvspan(Lmin,0.039,facecolor='#d9d9d9',edgecolor='k',hatch='////',lw=0.0,zorder=0)
ax.text(0.0275,1.2e4,'$\\hat b\\leq400$',fontsize=6.0,rotation=90,va='center',ha='center',bbox=dict(fc='w',ec='none',pad=0.5))
ax.set_xlabel('impact parameter $L=-\\ln(1-c)$')
ax.set_ylabel('twist $a=\\Delta-J-2$')
fig.tight_layout(pad=0.3)
out=os.path.join(HERE,'fig2_certificate.pdf'); fig.savefig(out); fig.savefig(os.path.join(HERE,'fig2_certificate_preview.png'),dpi=200)
print('wrote',out,len(cells),'high-twist cells')

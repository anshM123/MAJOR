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
from matplotlib.patches import Patch
amin,amax=998,3.0e4
Lmin,Lmax=0.012,50
INK,INK2='#1a1a1a','#555555'
plt.rcParams.update({'font.family':'serif','font.serif':['cmr10'],'mathtext.fontset':'cm','axes.formatter.use_mathtext':True,
                     'axes.unicode_minus':False,'font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,
                     'axes.linewidth':0.6,'axes.edgecolor':INK2,'xtick.color':INK2,'ytick.color':INK2,
                     'xtick.direction':'out','ytick.direction':'out','text.color':INK,'axes.labelcolor':INK,'pdf.fonttype':42,
                     'hatch.linewidth':0.5,'hatch.color':'#8c8c8c'})
fig,ax=plt.subplots(figsize=(3.375,2.85))
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(Lmin,Lmax); ax.set_ylim(amin,amax)
# a*(L): step function from the high-twist records
steps=[]
for (ba,bb),ast in cells:
    La,Lb=max(L_of_beta(ba),Lmin),min(L_of_beta(bb),Lmax)
    if Lb>La: steps.append((La,Lb,ast))
def astar(L):
    for La,Lb,a in steps:
        if La<=L<=Lb: return a
    return None
# 1. high twist: one opaque colour above a*(L)
HT='#fdd9b5'
_x=[];_y=[]
for La,Lb,a in steps: _x+=[La,Lb]; _y+=[a,a]
ax.fill_between(_x,_y,amax,facecolor=HT,lw=0,zorder=1)
# 2. routes for large spin near the twist edge: opaque sequential blues, only below a*(L), white 0.6pt gaps
bands=[(0.0513,0.2984,'#deebf7','direct prefix'),(0.2984,2.3026,'#c6dbef','coupled-cap norming'),
       (2.3026,3.65,'#9ecae1','signed transfer'),(3.65,5,'#6baed6',r'route T$^{\prime\prime}$'),
       (5,20,'#3787c0',r'route T$^{\prime}$'),(20,Lmax,'#1c5a99','tail theorem')]
for a,b,col,lab in bands:
    Ls=np.logspace(math.log10(a),math.log10(b),300)
    top=[astar(L) or amax for L in Ls]
    ax.fill_between(Ls,amin,top,facecolor=col,lw=0,zorder=1)
    ax.plot([a,a],[amin,astar(a*1.0001) or amax],color='w',lw=0.6,zorder=1.5)
xs=[];ys=[]
for La,Lb,a in steps: xs+=[La,Lb]; ys+=[a,a]
ax.plot(xs,ys,color='#a63603',lw=0.9,zorder=4)
# 3. b-hat <= 400 strip (ell <= 0.039), hatched, beneath the spin-by-spin region
ax.axvspan(Lmin,0.039,facecolor='#f2f2f2',edgecolor='#8c8c8c',hatch='/////',lw=0,zorder=2)
# 4. J <= 1000, certified spin by spin (b = 0): a <= 1000 (1 - sqrt c)/sqrt c, opaque grey on top
L=np.logspace(math.log10(Lmin),math.log10(Lmax),500); c=1-np.exp(-L); aJ=1000*(1-np.sqrt(c))/np.sqrt(c)
ax.fill_between(L,amin,np.clip(aJ,amin,amax),where=aJ>amin,facecolor='#a6a6a6',lw=0,zorder=3)
ax.plot(xs,ys,color='#a63603',lw=0.9,zorder=4)
# in-plot text: only the two regions with room for it
ax.text(1.2,1.25e4,r'high twist, $a\geq a^*$',ha='center',va='center',fontsize=7.5,zorder=5)
ax.text(0.016,1.15e3,r'$J\leq1000$',ha='left',va='bottom',fontsize=7.5,color='w',zorder=5)
ax.set_xlabel(r'impact parameter $L$')
ax.set_ylabel(r'twist $a=\Delta-J-2$')
handles=[Patch(facecolor=col,edgecolor='none',label=lab) for _,_,col,lab in bands]
handles+=[Patch(facecolor='#a6a6a6',edgecolor='none',label=r'spin by spin'),
          Patch(facecolor='#f2f2f2',edgecolor='#8c8c8c',hatch='/////',lw=0,label=r'$\hat b\leq400$'),
          Patch(facecolor=HT,edgecolor='none',label='high twist')]
fig.legend(handles=handles,loc='upper center',ncol=3,frameon=False,fontsize=7,handlelength=1.4,handleheight=0.9,
           columnspacing=1.0,handletextpad=0.5,borderaxespad=0.2,bbox_to_anchor=(0.55,1.0))
for sp in ('top','right'): ax.spines[sp].set_visible(False)
fig.subplots_adjust(left=0.15,right=0.97,bottom=0.13,top=0.78)
out=os.path.join(HERE,'fig2_certificate.pdf'); fig.savefig(out); fig.savefig(os.path.join(HERE,'fig2_certificate_preview.png'),dpi=200)
print('wrote',out,len(cells),'high-twist cells')

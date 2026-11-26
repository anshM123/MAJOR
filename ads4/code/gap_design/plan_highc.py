# Route-A (G4) cell plan for c in [c0,c1]: cells on the 1/2000 grid with twist overshoot omega<=OM, h-pieces m,
# K = least integer with the author's exact ordinary tail < 1e-7 (tight_spectral_tail.tail, envelope 2^2337),
# estimated time from the measured rate (1.3e-7 s * K * terms; terms from the cost model at anchor points).
import sys,math,json
from fractions import Fraction as F
sys.path.insert(0,sys.argv[4])          # directory with the author's tight_spectral_tail.py
from tight_spectral_tail import tail
c0,c1=F(sys.argv[1]),F(sys.argv[2]);OM=float(sys.argv[3]);G=F(1,2000)
def tw(c):return (1-math.sqrt(c))**2
def Kof(u):
    lo,hi=1000,17199
    while hi-lo>1:
        m=(lo+hi)//2
        if tail(u,m)<F(1,10**7):hi=m
        else:lo=m
    return hi
# per-channel CPU-s per unit c (cost model, omega=1.08, scaled checker) at anchors; linear interpolation in log
anch=[(0.56,6.5),(0.6,9.8),(0.65,15),(0.7,29),(0.74,43),(0.77,66),(0.80,100)]
def rate(c):
    for (a,ra),(b,rb) in zip(anch,anch[1:]):
        if a<=c<=b:return math.exp(math.log(ra)+(math.log(rb)-math.log(ra))*(c-a)/(b-a))*3600
    return anch[-1][1]*3600
cells=[];c=c0;tot=0
while c<c1:
    u=c+G
    while u+G<=c1 and tw(float(c))/tw(float(u+G))<=OM:u+=G
    K=Kof(u);m=2 if float(c)<0.7 else 3
    est=rate(float(c))*float(u-c)
    cells.append(dict(c_lo=str(c),c_hi=str(u),K=K,m=m,omega=round(tw(float(c))/tw(float(u)),4),est_s_per_channel=round(est)));tot+=est;c=u
json.dump(dict(rule=f'omega<={OM}, grid 1/2000, K: tail<1e-7 (2^2337 envelope), m=2 (c<0.7) / 3, allowance 1e-7+1e-14',cells=cells,est_total_cpu_h_both_channels=round(2*tot/3600,1)),open(f'plan_{sys.argv[1]}_{sys.argv[2]}.json','w'),indent=1)
print(len(cells),'cells; est. total CPU-h both channels',round(2*tot/3600,1))
for x in cells[:3]+cells[-3:]:print(x)

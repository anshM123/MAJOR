# cross-check: K=3299 table (sliver_direct_ext) vs author's 17,201 vector, coefficient midpoints and prefix values
import sys,math
import common,common2
from flint import arb
for c,b,K in [(0.3,0,1403),(0.5,4,2438),(0.55,0,2860)]:
    h=common2.hcap_tw(c)*0.7
    v1=common.S0(arb(c),arb(h),b,K,norm=False);v2=common2.Sh(c,h,b,K)
    print('c',c,'b',b,'K',K,'S(3299 table) %.12e  S(17200) %.12e  rel.diff %.2e'%(float(v1.mid()),float(v2.mid()),abs(float((v1-v2).mid())/float(v2.mid()))),flush=True)
A1,G1=common.coeffs();common2.load(3297);A2=common2._A;G2=common2._G
mx=0
for k in range(0,3297,7):
    d=abs(float(((A1[k]-A2[k])/A2[k]).mid()));mx=max(mx,d)
print('max rel diff A_k (sampled k):%.2e'%mx)

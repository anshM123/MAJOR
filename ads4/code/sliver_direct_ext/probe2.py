import sys,mpmath as mp
from probe import S0,caps
cl,cu,b,K=mp.mpf(sys.argv[1]),mp.mpf(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
c1,c2=caps(cl,cu,b);hc=min(c1,c2);print('cap1 %.4e cap2 %.4e  hp_K %.3f'%(c1,c2,hc*K*K))
for c in (cl,(cl+cu)/2,cu):
 print('c',float(c),' '.join('%.3e'%float(S0(c,hc*t,b,K)) for t in [0,.5,.9,1.0]))

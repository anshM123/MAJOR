from mixture_test import *
worst=1
for J in [5,10,20,30,40,60,80,100,150]:
    for b in (0,4):
        for h in [10.5, 12.3, 20.7, 30.5, 50.5, 100.5]:
            if J-b<0: continue
            P=Ydist(m_seq(b,J,mp.mpf(J+h)))
            mn=min(P)
            print(J,b,h,mp.nstr(mn,5))

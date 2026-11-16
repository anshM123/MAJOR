import sys, pickle, numpy as np
U='/home/user/MAJOR/ads4/code/uniform'; sys.path.insert(0,U); sys.path.insert(0,U+'/../robust')
from rcommon import gl01
pn,_=gl01(48)
for Dg in (20,100,1000):
    C=pickle.load(open(f'{U}/cache/skel_t3_{Dg}.pkl','rb'))
    st={}
    for (mfl,J),sk in C.items():
        if sk is None: continue
        b='<1.5' if mfl<1.5 else ('<3' if mfl<3 else '>=3')
        x=pn/mfl
        unc2=~sk[2][2]; unc3=~sk[3][2]; fw=not sk['fw'][4]
        d=st.setdefault(b,dict(n=0,any=0,fw=0,node=0,node_lowx=0,Jmax_node=0))
        d['n']+=1; d['any']+= bool(unc2.any() or unc3.any() or fw); d['fw']+=fw
        if unc2.any() or unc3.any():
            d['node']+=1; d['Jmax_node']=max(d['Jmax_node'],J)
            if (x[unc2|unc3]<0.6).any(): d['node_lowx']+=1
    print(Dg, st)

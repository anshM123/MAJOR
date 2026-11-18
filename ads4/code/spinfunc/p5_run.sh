#!/bin/bash
# p5: crystal points for the N question (antisymmetric even + odd at odd J; symmetric even + odd at even J)
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1 L1_DEG=6 P3_ANTI=1
for J in 5 4 7 9 6; do for D in 1600.37 3200.37 400.37; do
python3 -c "
import p3_odd_action as P
m=P.mphys($D,$J)
P.get('C2',$J,$D,[x*m for x in (0.2,0.3,0.4,0.5,0.6,0.7)])
"
done; done
P3_XS=0.3,0.5 P3_DS=1600 python3 -c "
import p3_odd_action as P
for J in (5,4):
    m=P.mphys(1600.0,J); P.get('C3H',J,1600.0,[x*m for x in (0.3,0.5)])
"
echo done

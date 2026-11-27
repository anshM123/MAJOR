#!/bin/bash
# wait for the 106-130 run to finish, then launch 161-180 and 181-192 (193-200 after sweep2 finishes)
cd /home/user/MAJOR/ads4/code/nf02_allgap/cert_copy
while pgrep -f "check_physical.py 130 106" > /dev/null; do sleep 20; done
nohup python3 check_physical.py 180 161 > phys_161_180.log 2>&1 &
while pgrep -f "sweep2.py 1000 6000" > /dev/null; do sleep 20; done
nohup python3 check_physical.py 192 181 > phys_181_192.log 2>&1 &
while pgrep -f "check_physical.py 160 131" > /dev/null; do sleep 20; done
nohup python3 check_physical.py 200 193 > phys_193_200.log 2>&1 &

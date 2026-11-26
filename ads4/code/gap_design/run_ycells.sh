#!/bin/bash
# curvilinear smoothness scan over the author's 21-rectangle c-grid up to 0.8051
cd /home/user/MAJOR/ads4/code/gap_design
for cell in "0.4638 0.5407" "0.5407 0.6003" "0.6003 0.6478" "0.6478 0.6864" "0.6864 0.7184" "0.7184 0.7453" "0.7453 0.7682" "0.7682 0.7879"; do
 for b in 0 4; do python3 box_cheb_y.py $cell $b 24 8; done
done

#!/usr/bin/env python3
"""Compare bundled low-frequency masses against logs in a local checkout."""
from pathlib import Path
from fractions import Fraction as F
import csv,re,sys
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path('.')
rows={}
for part in (0,1):
    p=root/f'ads4/code/analytic3/logs/resid3_3.55_3.65_p{part}.log'
    for lo,hi,m in re.findall(r'cell \[([\d.]+),([\d.]+)\] m=([\deE.+-]+):',p.read_text()):
        lo=F(lo);hi=F(hi);m=F(m)
        if 1<=lo<=199 and (lo not in rows or m>rows[lo][1]):rows[lo]=(hi,m)
for r in csv.DictReader((Path(__file__).parent/'residual_masses.csv').open()):
    assert rows[F(r['lo'])]==(F(r['frequency_ceiling']),F(r['logged_mass'])),r
assert len(rows)==199
print('PASS: all 199 bundled logged masses match source logs exactly.')

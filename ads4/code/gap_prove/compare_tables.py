#!/usr/bin/env python3
"""Exact overlap check of two direct coefficient tables (format: 'N D' header, rows 'A Ar G Gr' at scale 10^D).
Usage: compare_tables.py small.txt long.txt  -> asserts every ball pair intersects for k<=N_small; reports max |mid diff|."""
import sys
sys.set_int_max_str_digits(0)
f1,f2=open(sys.argv[1]),open(sys.argv[2])
N1,D1=map(int,f1.readline().split());N2,D2=map(int,f2.readline().split());assert N2>=N1 and D2>=D1
sc=10**(D2-D1);worst=None;wk=0
for k in range(N1+1):
 a,ar,g,gr=map(int,f1.readline().split());A,AR,G,GR=map(int,f2.readline().split())
 for m1,r1,m2,r2 in ((a*sc,ar*sc,A,AR),(g*sc,gr*sc,G,GR)):
  d=abs(m1-m2);assert d<=r1+r2,(k,'NO OVERLAP')
  rel=len(str(d)) if d else 0
  if worst is None or rel>worst:worst,wk=rel,k
print('PASS overlap all',N1+1,'pairs; max |mid diff| ~ 10^%d at k=%d (scale 10^-%d)'%(worst-D2,wk,0),flush=True)

#!/bin/bash
# EE_DG100.md: chain variants after the restart (main, can, ex_*, jhi* already done before it) -> cache/y2_chain100_<tag>.txt
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Y2_XROWS=1 Y2_XEVEN=0 Y2_MC=40 Y2_TOL=1 python3 y2_chain100.py xrows_odd > /dev/null 2>&1
Y2_XROWS=1 Y2_MC=40 Y2_TOL=1 python3 y2_chain100.py xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_GAUGE=can Y2_XROWS=1 Y2_XEVEN=0 Y2_MODES=nl python3 y2_chain100.py can_xrows > /dev/null 2>&1
for e in hold zero; do
  Y2_KINDS=EE Y2_EXTRAP=$e Y2_XROWS=1 Y2_XEVEN=0 Y2_MODES=nl python3 y2_chain100.py ex_${e}_xrows > /dev/null 2>&1
done
Y2_KINDS=EE Y2_NUIM=0 Y2_XROWS=1 Y2_XEVEN=0 Y2_MODES=nl python3 y2_chain100.py imK2_xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_AFITXR=p8 Y2_XROWS=1 Y2_XEVEN=0 Y2_MODES=nl,loc python3 y2_chain100.py afitp8_xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_AFITNU=0 Y2_MODES=nl,loc python3 y2_chain100.py afit0 > /dev/null 2>&1
Y2_KINDS=EE Y2_GERR=300 Y2_XROWS=1 Y2_XEVEN=0 Y2_MODES=nl python3 y2_chain100.py gerr300_xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_DIRECT=0 Y2_XROWS=1 Y2_XEVEN=0 Y2_MODES=loc,nl python3 y2_chain100.py k2only_xrows > /dev/null 2>&1
echo ALLDONE

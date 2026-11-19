#!/bin/bash
# EE_DG100.md: all chain variants (y2_chain100.py) -> cache/y2_chain100_<tag>.txt  (run sequentially, one process)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Y2_MC=40 Y2_TOL=1 python3 y2_chain100.py main > /dev/null 2>&1
Y2_XROWS=1 Y2_MC=40 Y2_TOL=1 python3 y2_chain100.py xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_GAUGE=can Y2_MODES=nl Y2_MC=30 python3 y2_chain100.py can > /dev/null 2>&1
Y2_KINDS=EE Y2_GAUGE=can Y2_XROWS=1 Y2_MODES=nl python3 y2_chain100.py can_xrows > /dev/null 2>&1
for e in hold zero; do
  Y2_KINDS=EE Y2_EXTRAP=$e Y2_MODES=nl python3 y2_chain100.py ex_$e > /dev/null 2>&1
  Y2_KINDS=EE Y2_EXTRAP=$e Y2_XROWS=1 Y2_MODES=nl python3 y2_chain100.py ex_${e}_xrows > /dev/null 2>&1
done
Y2_KINDS=EE Y2_NUIM=0 Y2_XROWS=1 Y2_MODES=nl python3 y2_chain100.py imK2_xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_AFITNU=0 Y2_XROWS=1 Y2_MODES=nl,loc python3 y2_chain100.py afit0_xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_GERR=300 Y2_XROWS=1 Y2_MODES=nl python3 y2_chain100.py gerr300_xrows > /dev/null 2>&1
Y2_KINDS=EE Y2_DIRECT=0 Y2_XROWS=1 Y2_MODES=loc,nl python3 y2_chain100.py k2only_xrows > /dev/null 2>&1
Y2_JHI=1 Y2_KINDS=EE Y2_MODES=nl,loc python3 y2_chain100.py jhi > /dev/null 2>&1
Y2_JHI=1 Y2_BELOW=scale Y2_KINDS=EE Y2_MODES=nl,loc python3 y2_chain100.py jhi_scale > /dev/null 2>&1
echo ALLDONE

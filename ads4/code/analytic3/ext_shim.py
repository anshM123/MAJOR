#!/usr/bin/env python3
"""ext_shim.py SCRIPT [args...]: run one of the sliver_extension package's scripts unchanged (EXT_DIR), with
check_norm.HERE pointed at BIN_DIR (where the identical scalar.cpp is compiled), so nothing is written into the package."""
import os, sys, runpy
from pathlib import Path
EXT = Path(os.environ['EXT_DIR']); BIN = Path(os.environ['BIN_DIR'])
sys.path.insert(0, str(EXT))
import check_norm as CN
CN.HERE = BIN
script = sys.argv[1]
sys.argv = [str(EXT / script)] + sys.argv[2:]
if script == 'check_norm.py':
    CN.main()
else:
    runpy.run_path(str(EXT / script), run_name='__main__')

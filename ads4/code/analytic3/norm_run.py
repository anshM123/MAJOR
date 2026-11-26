#!/usr/bin/env python3
"""norm_run.py: run code/norming/check_norm.py unchanged, but with its scalar evaluator taken from analytic3/bin/scalar
(compiled from the identical norming/scalar.cpp), so that nothing is written into code/norming/."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str((HERE / '..' / 'norming').resolve()))
import check_norm as CN
CN.HERE = HERE / 'bin'
sys.argv = ['check_norm.py'] + sys.argv[1:]
CN.main()
